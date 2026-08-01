//! Online (streaming) inference — the first phase of salmon's dual-phase
//! estimator (`processMiniBatch` in `SalmonQuantify.cpp`).
//!
//! # Why there are two phases
//!
//! salmon's final abundance estimate comes from an offline EM over equivalence
//! classes. But several models the EM depends on — the fragment-length
//! distribution and the sequence/GC/positional bias models — have to be *learned
//! from the data first*, and learning them correctly needs a weight per
//! observation: if a fragment could have come from three transcripts, how much
//! should each of them count towards the bias model?
//!
//! A uniform split is wrong, because the transcripts are not equally abundant.
//! The right weight is the posterior probability, which needs abundances, which
//! is what we are trying to estimate. The online phase breaks that circularity:
//! as reads stream past it maintains a running, ever-improving abundance estimate
//! and uses it to weight each fragment's contribution.
//!
//! # How the running estimate works
//!
//! While reads are mapped, each transcript carries a running **log-mass**
//! `logAdd(priorMass, accumulated)` (prior `log(α·len)`, α=0.05). For each
//! fragment, the per-transcript posterior is `softmax(mass_t + auxWeight_t)`
//! over its candidate transcripts; the masses are then incremented by
//! `logForgettingMass + log(posterior)` (an online-EM update with a decaying
//! forgetting-mass schedule, `forgettingFactor` 0.65).
//!
//! The "forgetting mass" is what makes this an *online* EM: early fragments were
//! weighted by a poor abundance estimate, so their contributions are gradually
//! discounted relative to later, better-informed ones.
//!
//! Crucially the posterior is **abundance-aware**, so it is the correct weight
//! for collecting the observed bias models — which is exactly what the
//! offline-only path lacked (see the `salmon-online-phase` analysis). The offline
//! EM still produces the final point estimate; the online phase exists to develop
//! these weights and the fragment-length distribution.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// `log(exp(a)+exp(b))` with `-inf` as log-zero.
///
/// A local copy of [`salmon_core::math::log_add`]; see that module for why
/// probabilities live in log space and how the stable formulation works.
#[inline]
fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

/// An `f64` accumulator supporting lock-free log-space addition.
///
/// Same bit-pattern-in-an-integer trick as [`salmon_core::AtomicF64`]: CPUs have
/// atomic integer instructions but no atomic floats, and mapping threads update
/// these masses continuously.
struct AtomicF64(AtomicU64);
impl AtomicF64 {
    fn new(v: f64) -> Self {
        Self(AtomicU64::new(v.to_bits()))
    }
    #[inline]
    fn load(&self) -> f64 {
        f64::from_bits(self.0.load(Ordering::Relaxed))
    }
    /// `self = logAdd(self, x)` (CAS retry loop).
    ///
    /// "CAS" is compare-and-swap: store the new value only if nobody else changed
    /// the slot since we read it, otherwise retry on top of their value.
    #[inline]
    fn log_add_assign(&self, x: f64) {
        // Adding log-zero is a no-op; skipping it avoids a pointless CAS.
        if x == f64::NEG_INFINITY {
            return;
        }
        let mut cur = self.0.load(Ordering::Relaxed);
        loop {
            let new = log_add(f64::from_bits(cur), x).to_bits();
            match self
                .0
                .compare_exchange_weak(cur, new, Ordering::Relaxed, Ordering::Relaxed)
            {
                Ok(_) => break,
                Err(e) => cur = e,
            }
        }
    }
}

/// Shared online-inference state, updated concurrently by the mapping workers.
pub struct OnlineInference {
    /// per-transcript prior log-mass `log(α·len)`
    ///
    /// Length-proportional because, absent any evidence, a longer transcript is
    /// expected to yield proportionally more fragments. α scales how much this
    /// prior matters relative to observed data.
    prior_mass: Vec<f64>,
    /// per-transcript accumulated log-mass (starts at log 0 = `-inf`)
    mass: Vec<AtomicF64>,
    /// forgetting factor (salmon `ffactor`, default 0.65)
    forgetting_factor: f64,
    /// lazily grown forgetting-mass schedule, indexed by minibatch timestep
    ///
    /// Behind a `Mutex` because it is extended once per *batch* (rare) rather
    /// than per fragment, so the lock is never contended in practice.
    log_fm: Mutex<Vec<f64>>,
    /// fragments assigned so far (for the burn-in cutoff)
    num_assigned: AtomicU64,
    /// stop updating models after this many assigned fragments (salmon 5,000,000)
    ///
    /// Past this point the models have converged and further updates only cost
    /// time — and, since threads race, would add nondeterminism for no gain.
    burnin_frags: u64,
}

impl OnlineInference {
    /// Build online state for `ref_lens` transcript lengths. `alpha` is the
    /// per-base prior (salmon 0.05), `forgetting_factor` the decay (0.65),
    /// `burnin_frags` the model-update cutoff (5,000,000).
    pub fn new(ref_lens: &[u64], alpha: f64, forgetting_factor: f64, burnin_frags: u64) -> Self {
        let prior_mass: Vec<f64> = ref_lens
            .iter()
            // `.max(1)` keeps `ln` defined for a zero-length reference.
            .map(|&l| (alpha * (l.max(1)) as f64).ln())
            .collect();
        // Starts at log-zero: no fragments observed yet, so all the initial
        // signal comes from the prior.
        let mass = ref_lens
            .iter()
            .map(|_| AtomicF64::new(f64::NEG_INFINITY))
            .collect();
        Self {
            prior_mass,
            mass,
            forgetting_factor,
            log_fm: Mutex::new(vec![0.0]), // lfm[0] = log 1 = 0
            num_assigned: AtomicU64::new(0),
            burnin_frags,
        }
    }

    /// Current online log-mass (abundance estimate) of transcript `tid`:
    /// `logAdd(priorMass, accumulated)`.
    ///
    /// Prior plus evidence, so a transcript with no fragments yet still has a
    /// non-zero (prior) mass and can therefore still win posterior weight.
    #[inline]
    pub fn mass_log(&self, tid: usize) -> f64 {
        log_add(self.prior_mass[tid], self.mass[tid].load())
    }

    /// Grab the next minibatch's `logForgettingMass` (one per mapping batch),
    /// extending the schedule as needed:
    /// `lfm[k] = lfm[k-1] + ff·ln(k) - ln((k+1)^ff - 1)`, `lfm[0] = 0`.
    ///
    /// The schedule decreases with `k`, so each successive batch contributes less
    /// new mass than the last; that is what "forgets" the early, badly-weighted
    /// fragments. The exact recurrence is salmon's, chosen so the implied weights
    /// follow a `k^-ff` decay.
    pub fn next_log_fm(&self) -> f64 {
        let ff = self.forgetting_factor;
        let mut v = self.log_fm.lock().unwrap();
        let k = v.len();
        let fm = v[k - 1] + ff * (k as f64).ln() - ((k as f64 + 1.0).powf(ff) - 1.0).ln();
        // Memoized so the recurrence is never recomputed from scratch.
        v.push(fm);
        fm
    }

    /// Whether model collection is still in the (pre-burn-in) online window.
    #[inline]
    pub fn collecting(&self) -> bool {
        self.num_assigned.load(Ordering::Relaxed) < self.burnin_frags
    }

    /// Fragments assigned so far (salmon's `numAssignedFragments`); the caller
    /// uses this against `numPreBurninFrags` to decide whether to fold the
    /// fragment-length / auxiliary terms into the per-fragment probability.
    #[inline]
    pub fn num_assigned(&self) -> u64 {
        self.num_assigned.load(Ordering::Relaxed)
    }

    /// Process one fragment's compatible `(tid, log_aux)` mappings, where
    /// `log_aux` is the abundance-*independent* log auxiliary probability
    /// (`logFragCov + startPosProb + logFragProb` — alignment score, length
    /// normalization, and fragment-length probability). Returns the
    /// abundance-aware posterior over the same transcripts (aligned to input
    /// order) and updates the online masses by `logForgettingMass + log(posterior)`.
    /// `log_fm` is the batch's forgetting mass. The posteriors are the correct
    /// weights for bias collection.
    ///
    /// This is the whole online EM in one function: an E-step (compute the
    /// posterior from the current masses) immediately followed by an M-step
    /// (fold that posterior back into the masses).
    pub fn assign_fragment(&self, maps: &[(u32, f64)], log_fm: f64) -> Vec<f64> {
        let n = maps.len();
        if n == 0 {
            return Vec::new();
        }
        // unnormalized log-posterior = mass(t) + log_aux
        //
        // Adding logs multiplies probabilities: "how abundant is this transcript"
        // times "how well does this fragment fit it". `denom` accumulates the
        // normalizer in log space as we go, so the loop is single-pass.
        let mut unnorm = Vec::with_capacity(n);
        let mut denom = f64::NEG_INFINITY;
        for &(tid, log_aux) in maps {
            let u = self.mass_log(tid as usize) + log_aux;
            unnorm.push(u);
            denom = log_add(denom, u);
        }
        let mut post = Vec::with_capacity(n);
        if denom == f64::NEG_INFINITY {
            // degenerate; fall back to uniform
            //
            // Every candidate had zero probability (e.g. an auxiliary term
            // vanished). Dividing by zero would give NaN, so split the fragment
            // evenly rather than discarding it.
            let u = 1.0 / n as f64;
            for &(tid, _) in maps {
                self.mass[tid as usize].log_add_assign(log_fm + u.ln());
                post.push(u);
            }
        } else {
            for (i, &(tid, _)) in maps.iter().enumerate() {
                // Dividing by the normalizer = subtracting it in log space; this
                // is exactly a softmax over the candidates.
                let lp = unnorm[i] - denom; // log posterior
                self.mass[tid as usize].log_add_assign(log_fm + lp);
                post.push(lp.exp());
            }
        }
        self.num_assigned.fetch_add(1, Ordering::Relaxed);
        post
    }
}
