//! CUDA backend of the GPU warm start (Linux), written with cuTile Rust
//! (NVIDIA CUDA Rust's Tile track; kernels are JIT-compiled through Tile IR on
//! first launch, which needs CUDA 13.3 and compute capability >= 8.0).
//!
//! Every kernel works on tiles of `B` lanes (classes or transcripts) and only
//! writes its own dense, partitioned output: gathers, no scatters, no atomics.
//! Variable-length rows are walked with masked pointer gathers up to the
//! tile's longest row.
//!
//! Untested: this has not run on NVIDIA hardware; CI (`gpu-cuda-check`) only
//! compiles it against the CUDA 13 headers.
#![allow(unsafe_code)]

use cutile::prelude::*;

use super::Problem;

/// Lanes per tile.
const B: usize = 128;

#[cutile::module]
mod kernels {
    use cutile::core::*;

    /// VBEM expectations `exp(ψ(alpha + prior))`, 0 at or below `DIGAMMA_MIN`:
    /// the f32 port of `components.rs` `exp_digamma`, branch-free on tiles.
    #[cutile::entry()]
    fn et_kernel<const B: i32>(
        et: &mut Tensor<f32, { [B] }>,
        alpha: &Tensor<f32, { [-1] }>,
        prior: &Tensor<f32, { [-1] }>,
    ) {
        let a: Tile<f32, { [B] }> = alpha.load_like(et);
        let p: Tile<f32, { [B] }> = prior.load_like(et);
        let x: Tile<f32, { [B] }> = a + p;
        let zero: Tile<f32, { [B] }> = constant(0.0f32, et.shape());
        let one: Tile<f32, { [B] }> = constant(1.0f32, et.shape());
        let two: Tile<f32, { [B] }> = constant(2.0f32, et.shape());
        let ten: Tile<f32, { [B] }> = constant(10.0f32, et.shape());
        let mut num: Tile<f32, { [B] }> = zero;
        let mut den: Tile<f32, { [B] }> = one;
        let mut j: Tile<f32, { [B] }> = zero;
        for _i in 0i32..5i32 {
            let a: Tile<f32, { [B] }> = x + j;
            let b: Tile<f32, { [B] }> = a + one;
            num = num * (a * b) + (a + b) * den;
            den = den * (a * b);
            j = j + two;
        }
        let small: Tile<bool, { [B] }> = lt_tile(x, ten);
        let y: Tile<f32, { [B] }> = select(small, x + ten, x);
        let s: Tile<f32, { [B] }> = select(small, num / den, zero);
        let r: Tile<f32, { [B] }> = one / y;
        let r2: Tile<f32, { [B] }> = r * r;
        let c1: Tile<f32, { [B] }> = constant(0.5f32, et.shape());
        let c2: Tile<f32, { [B] }> = constant(1.0f32 / 12.0f32, et.shape());
        let c3: Tile<f32, { [B] }> = constant(1.0f32 / 120.0f32, et.shape());
        let c4: Tile<f32, { [B] }> = constant(1.0f32 / 252.0f32, et.shape());
        let c5: Tile<f32, { [B] }> = constant(1.0f32 / 240.0f32, et.shape());
        let c6: Tile<f32, { [B] }> = constant(1.0f32 / 132.0f32, et.shape());
        let c7: Tile<f32, { [B] }> = constant(691.0f32 / 32760.0f32, et.shape());
        let u: Tile<f32, { [B] }> =
            zero - c1 * r - r2 * (c2 - r2 * (c3 - r2 * (c4 - r2 * (c5 - r2 * (c6 - r2 * c7)))));
        let e: Tile<f32, { [B] }> = y * exp(u - s);
        let live: Tile<bool, { [B] }> = gt_tile(x, constant(1e-10f32, et.shape()));
        et.store(select(live, e, zero));
    }

    /// `scale[c] = count[c] / Σ basis[label]·w` over class `c`'s members,
    /// 0 when the sum is not positive.
    #[cutile::entry()]
    fn class_scale<const B: i32>(
        scale: &mut Tensor<f32, { [B] }>,
        lo: &Tensor<i32, { [-1] }>,
        len: &Tensor<i32, { [-1] }>,
        count: &Tensor<f32, { [-1] }>,
        labels: &Tensor<i32, { [-1] }>,
        w: &Tensor<f32, { [-1] }>,
        basis: &Tensor<f32, { [-1] }>,
    ) {
        let lo: Tile<i32, { [B] }> = lo.load_like(scale);
        let len: Tile<i32, { [B] }> = len.load_like(scale);
        let longest: Tile<i32, { [] }> = reduce_max(len, 0i32);
        let n: i32 = tile_to_scalar(longest);
        let p: PointerTile<*const i32, { [] }> = pointer_to_tile(labels.as_ptr());
        let p: PointerTile<*const i32, { [1] }> = p.reshape(shape![1]);
        let lab_at: PointerTile<*const i32, { [B] }> = p.broadcast(scale.shape());
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(w.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let w_at: PointerTile<*const f32, { [B] }> = p.broadcast(scale.shape());
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(basis.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let basis_at: PointerTile<*const f32, { [B] }> = p.broadcast(scale.shape());
        let zero: Tile<f32, { [B] }> = constant(0.0f32, scale.shape());
        let mut denom: Tile<f32, { [B] }> = zero;
        for k in 0i32..n {
            let kk: Tile<i32, { [B] }> = broadcast_scalar(k, scale.shape());
            let live: Tile<bool, { [B] }> = lt_tile(kk, len);
            let at: Tile<i32, { [B] }> = lo + kk;
            // SAFETY: unmasked lanes are inside class rows of the uploaded
            // CSR, whose labels index `basis`.
            let (lab, _t0): (Tile<i32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(lab_at, at),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0i32),
                    None,
                    Latency::<0>,
                )
            };
            let (wt, _t1): (Tile<f32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(w_at, at),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0.0f32),
                    None,
                    Latency::<0>,
                )
            };
            let (b, _t2): (Tile<f32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(basis_at, lab),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0.0f32),
                    None,
                    Latency::<0>,
                )
            };
            // `basis >= 0`, so this is also VBEM's "positive expectations only"
            denom = denom + b * wt;
        }
        let cnt: Tile<f32, { [B] }> = count.load_like(scale);
        let ok: Tile<bool, { [B] }> = gt_tile(denom, constant(1.17549435e-38f32, scale.shape()));
        scale.store(select(ok, cnt / denom, zero));
    }

    /// `out[t] = single[t] + basis[t] · Σ w·scale[cls]` over transcript `t`'s classes.
    #[cutile::entry()]
    fn txp_gather<const B: i32>(
        out: &mut Tensor<f32, { [B] }>,
        lo: &Tensor<i32, { [-1] }>,
        len: &Tensor<i32, { [-1] }>,
        cls: &Tensor<i32, { [-1] }>,
        w: &Tensor<f32, { [-1] }>,
        scale: &Tensor<f32, { [-1] }>,
        basis: &Tensor<f32, { [-1] }>,
        single: &Tensor<f32, { [-1] }>,
    ) {
        let lo: Tile<i32, { [B] }> = lo.load_like(out);
        let len: Tile<i32, { [B] }> = len.load_like(out);
        let longest: Tile<i32, { [] }> = reduce_max(len, 0i32);
        let n: i32 = tile_to_scalar(longest);
        let p: PointerTile<*const i32, { [] }> = pointer_to_tile(cls.as_ptr());
        let p: PointerTile<*const i32, { [1] }> = p.reshape(shape![1]);
        let cls_at: PointerTile<*const i32, { [B] }> = p.broadcast(out.shape());
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(w.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let w_at: PointerTile<*const f32, { [B] }> = p.broadcast(out.shape());
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(scale.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let scale_at: PointerTile<*const f32, { [B] }> = p.broadcast(out.shape());
        let mut acc: Tile<f32, { [B] }> = constant(0.0f32, out.shape());
        for k in 0i32..n {
            let kk: Tile<i32, { [B] }> = broadcast_scalar(k, out.shape());
            let live: Tile<bool, { [B] }> = lt_tile(kk, len);
            let at: Tile<i32, { [B] }> = lo + kk;
            // SAFETY: unmasked lanes are inside transcript rows of the uploaded
            // transpose, whose class ids index `scale`.
            let (c, _t0): (Tile<i32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(cls_at, at),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0i32),
                    None,
                    Latency::<0>,
                )
            };
            let (wt, _t1): (Tile<f32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(w_at, at),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0.0f32),
                    None,
                    Latency::<0>,
                )
            };
            let (sc, _t2): (Tile<f32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(scale_at, c),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0.0f32),
                    None,
                    Latency::<0>,
                )
            };
            acc = acc + wt * sc;
        }
        let b: Tile<f32, { [B] }> = basis.load_like(out);
        let s: Tile<f32, { [B] }> = single.load_like(out);
        out.store(s + b * acc);
    }
}

fn err(e: impl std::fmt::Debug) -> String {
    format!("{e:?}")
}

/// `v` zero-padded to a whole, non-zero number of tiles, on the device.
fn upload<T: DType + Default>(mut v: Vec<T>) -> Result<Tensor<T>, String> {
    v.resize(v.len().div_ceil(B).max(1) * B, T::default());
    api::copy_host_vec_to_device(&Arc::new(v))
        .sync()
        .map_err(err)
}

/// CSR offsets as per-row `(start, length)` lanes.
fn rows(start: &[u32]) -> (Vec<i32>, Vec<i32>) {
    start
        .windows(2)
        .map(|s| (s[0] as i32, (s[1] - s[0]) as i32))
        .unzip()
}

/// Whether a CUDA device answers at all (tests skip otherwise).
#[cfg(test)]
pub(super) fn available() -> bool {
    api::zeros::<f32>(&[B]).sync().is_ok()
}

pub(super) struct Device {
    vb: bool,
    c_lo: Tensor<i32>,
    c_len: Tensor<i32>,
    count: Tensor<f32>,
    labels: Tensor<i32>,
    w: Tensor<f32>,
    scale: Tensor<f32>,
    t_lo: Tensor<i32>,
    t_len: Tensor<i32>,
    t_cls: Tensor<i32>,
    t_w: Tensor<f32>,
    single: Tensor<f32>,
    prior: Tensor<f32>,
    et: Tensor<f32>,
    /// ping-pong iterates; `alpha[cur]` is the latest
    alpha: [Tensor<f32>; 2],
    cur: usize,
    nt: usize,
}

impl Device {
    pub(super) fn new(
        prob: &Problem,
        init: &[f64],
        prior: &[f64],
        vb: bool,
    ) -> Result<Self, String> {
        let nt = prob.num_txps();
        let f32s = |v: &[f64]| -> Vec<f32> { v.iter().map(|&x| x as f32).collect() };
        // `Problem::new` bounds every index by `i32::MAX`.
        let i32s = |v: &[u32]| -> Vec<i32> { v.iter().map(|&x| x as i32).collect() };
        let (c_lo, c_len) = rows(&prob.c_start);
        let (t_lo, t_len) = rows(&prob.t_start);
        Ok(Self {
            vb,
            c_lo: upload(c_lo)?,
            c_len: upload(c_len)?,
            count: upload(prob.count.clone())?,
            labels: upload(i32s(&prob.labels))?,
            w: upload(prob.w.clone())?,
            scale: upload(vec![0f32; prob.num_classes()])?,
            t_lo: upload(t_lo)?,
            t_len: upload(t_len)?,
            t_cls: upload(i32s(&prob.t_cls))?,
            t_w: upload(prob.t_w.clone())?,
            single: upload(prob.single.clone())?,
            prior: upload(f32s(prior))?,
            et: upload(vec![0f32; nt])?,
            alpha: [upload(f32s(init))?, upload(vec![0f32; nt])?],
            cur: 0,
            nt,
        })
    }

    /// Run `n` M-steps.
    pub(super) fn run(&mut self, n: u32) -> Result<(), String> {
        for _ in 0..n {
            let (lo, hi) = self.alpha.split_at_mut(1);
            let (src, dst) = if self.cur == 0 {
                (&lo[0], &mut hi[0])
            } else {
                (&hi[0], &mut lo[0])
            };
            if self.vb {
                kernels::et_kernel((&mut self.et).partition([B]), src, &self.prior)
                    .sync()
                    .map_err(err)?;
            }
            let basis = if self.vb { &self.et } else { src };
            kernels::class_scale(
                (&mut self.scale).partition([B]),
                &self.c_lo,
                &self.c_len,
                &self.count,
                &self.labels,
                &self.w,
                basis,
            )
            .sync()
            .map_err(err)?;
            kernels::txp_gather(
                dst.partition([B]),
                &self.t_lo,
                &self.t_len,
                &self.t_cls,
                &self.t_w,
                &self.scale,
                basis,
                &self.single,
            )
            .sync()
            .map_err(err)?;
            self.cur = 1 - self.cur;
        }
        Ok(())
    }

    /// `(previous, latest)` iterates.
    pub(super) fn read(&self) -> Result<(Vec<f32>, Vec<f32>), String> {
        let host = |t: &Tensor<f32>| -> Result<Vec<f32>, String> {
            let mut v = t
                .dup()
                .sync()
                .map_err(err)?
                .to_host_vec()
                .sync()
                .map_err(err)?;
            v.truncate(self.nt); // drop the tile padding
            Ok(v)
        };
        Ok((
            host(&self.alpha[1 - self.cur])?,
            host(&self.alpha[self.cur])?,
        ))
    }
}
