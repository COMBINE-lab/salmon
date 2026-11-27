//! CUDA backend of the GPU warm start (Linux), written with cuTile Rust
//! (NVIDIA CUDA Rust's Tile track; kernels are JIT-compiled through Tile IR on
//! first launch, which needs CUDA 13.3 and compute capability >= 8.0).
//!
//! Every kernel works on tiles of `B` lanes and only writes its own dense,
//! partitioned output: gathers, no scatters, no atomics. `chunk_dot` has one
//! lane per chunk of at most `CHUNK_LEN` incidences, so a tile's lanes are
//! about equally busy; `class_finish`/`txp_finish` then sum each row's chunks.
//! Rows are walked with masked pointer gathers up to the tile's longest row.
//!
//! Untested: this has not run on NVIDIA hardware; CI (`gpu-cuda-check`) only
//! compiles it against the CUDA 13 headers.
#![allow(unsafe_code)]

use cutile::cuda_core;
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
            let xa: Tile<f32, { [B] }> = x + j;
            let xb: Tile<f32, { [B] }> = xa + one;
            num = num * (xa * xb) + (xa + xb) * den;
            den = den * (xa * xb);
            j = j + two;
        }
        let small: Tile<bool, { [B] }> = lt_tile(x, ten);
        let y: Tile<f32, { [B] }> = select(small, x + ten, x);
        let s: Tile<f32, { [B] }> = select(small, num / den, zero);
        let r: Tile<f32, { [B] }> = one / y;
        let r2: Tile<f32, { [B] }> = r * r;
        // literals, not `1.0 / 12.0`: `constant` lowers to a Tile IR dense attribute
        let c1: Tile<f32, { [B] }> = constant(0.5f32, et.shape());
        let c2: Tile<f32, { [B] }> = constant(0.083_333_336f32, et.shape());
        let c3: Tile<f32, { [B] }> = constant(0.008_333_334f32, et.shape());
        let c4: Tile<f32, { [B] }> = constant(0.003_968_254f32, et.shape());
        let c5: Tile<f32, { [B] }> = constant(0.004_166_667f32, et.shape());
        let c6: Tile<f32, { [B] }> = constant(0.007_575_758f32, et.shape());
        let c7: Tile<f32, { [B] }> = constant(0.021_092_797f32, et.shape());
        let u: Tile<f32, { [B] }> =
            zero - c1 * r - r2 * (c2 - r2 * (c3 - r2 * (c4 - r2 * (c5 - r2 * (c6 - r2 * c7)))));
        let e: Tile<f32, { [B] }> = y * exp(u - s);
        let live: Tile<bool, { [B] }> = gt_tile(x, constant(1e-10f32, et.shape()));
        et.store(select(live, e, zero));
    }

    /// `partial[k] = Σ x[idx]·w` over chunk `k` (at most `CHUNK_LEN`
    /// incidences), for the class pass and the transcript pass alike.
    #[cutile::entry()]
    fn chunk_dot<const B: i32>(
        partial: &mut Tensor<f32, { [B] }>,
        lo: &Tensor<i32, { [-1] }>,
        len: &Tensor<i32, { [-1] }>,
        idx: &Tensor<i32, { [-1] }>,
        w: &Tensor<f32, { [-1] }>,
        x: &Tensor<f32, { [-1] }>,
    ) {
        let lo: Tile<i32, { [B] }> = lo.load_like(partial);
        let len: Tile<i32, { [B] }> = len.load_like(partial);
        let longest: Tile<i32, { [] }> = reduce_max(len, 0i32);
        let n: i32 = tile_to_scalar(longest);
        let p: PointerTile<*const i32, { [] }> = pointer_to_tile(idx.as_ptr());
        let p: PointerTile<*const i32, { [1] }> = p.reshape(shape![1]);
        let idx_at: PointerTile<*const i32, { [B] }> = p.broadcast(partial.shape());
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(w.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let w_at: PointerTile<*const f32, { [B] }> = p.broadcast(partial.shape());
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(x.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let x_at: PointerTile<*const f32, { [B] }> = p.broadcast(partial.shape());
        let mut acc: Tile<f32, { [B] }> = constant(0.0f32, partial.shape());
        for k in 0i32..n {
            let kk: Tile<i32, { [B] }> = broadcast_scalar(k, partial.shape());
            let live: Tile<bool, { [B] }> = lt_tile(kk, len);
            let at: Tile<i32, { [B] }> = lo + kk;
            // SAFETY: unmasked lanes are inside chunks of the uploaded CSR,
            // whose indices address `x`.
            let (i, _t0): (Tile<i32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(idx_at, at),
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
            let (xv, _t2): (Tile<f32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(x_at, i),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0.0f32),
                    None,
                    Latency::<0>,
                )
            };
            // `x >= 0` (basis or class scale), so no positivity mask is needed
            acc = acc + xv * wt;
        }
        partial.store(acc);
    }

    /// `Σ partial` over each lane's chunks `lo..lo + len`.
    #[cutile::entry()]
    fn class_finish<const B: i32>(
        scale: &mut Tensor<f32, { [B] }>,
        lo: &Tensor<i32, { [-1] }>,
        len: &Tensor<i32, { [-1] }>,
        partial: &Tensor<f32, { [-1] }>,
        count: &Tensor<f32, { [-1] }>,
    ) {
        let lo: Tile<i32, { [B] }> = lo.load_like(scale);
        let len: Tile<i32, { [B] }> = len.load_like(scale);
        let longest: Tile<i32, { [] }> = reduce_max(len, 0i32);
        let n: i32 = tile_to_scalar(longest);
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(partial.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let p_at: PointerTile<*const f32, { [B] }> = p.broadcast(scale.shape());
        let zero: Tile<f32, { [B] }> = constant(0.0f32, scale.shape());
        let mut denom: Tile<f32, { [B] }> = zero;
        for k in 0i32..n {
            let kk: Tile<i32, { [B] }> = broadcast_scalar(k, scale.shape());
            let live: Tile<bool, { [B] }> = lt_tile(kk, len);
            let at: Tile<i32, { [B] }> = lo + kk;
            // SAFETY: unmasked lanes index this class's chunks in `partial`.
            let (v, _t0): (Tile<f32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(p_at, at),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0.0f32),
                    None,
                    Latency::<0>,
                )
            };
            denom = denom + v;
        }
        let cnt: Tile<f32, { [B] }> = count.load_like(scale);
        let ok: Tile<bool, { [B] }> = gt_tile(denom, constant(1.17549435e-38f32, scale.shape()));
        scale.store(select(ok, cnt / denom, zero));
    }

    /// `out[t] = single[t] + basis[t] · Σ partial` over transcript `t`'s chunks.
    #[cutile::entry()]
    fn txp_finish<const B: i32>(
        out: &mut Tensor<f32, { [B] }>,
        lo: &Tensor<i32, { [-1] }>,
        len: &Tensor<i32, { [-1] }>,
        partial: &Tensor<f32, { [-1] }>,
        basis: &Tensor<f32, { [-1] }>,
        single: &Tensor<f32, { [-1] }>,
    ) {
        let lo: Tile<i32, { [B] }> = lo.load_like(out);
        let len: Tile<i32, { [B] }> = len.load_like(out);
        let longest: Tile<i32, { [] }> = reduce_max(len, 0i32);
        let n: i32 = tile_to_scalar(longest);
        let p: PointerTile<*const f32, { [] }> = pointer_to_tile(partial.as_ptr());
        let p: PointerTile<*const f32, { [1] }> = p.reshape(shape![1]);
        let p_at: PointerTile<*const f32, { [B] }> = p.broadcast(out.shape());
        let mut acc: Tile<f32, { [B] }> = constant(0.0f32, out.shape());
        for k in 0i32..n {
            let kk: Tile<i32, { [B] }> = broadcast_scalar(k, out.shape());
            let live: Tile<bool, { [B] }> = lt_tile(kk, len);
            let at: Tile<i32, { [B] }> = lo + kk;
            // SAFETY: unmasked lanes index this transcript's chunks in `partial`.
            let (v, _t0): (Tile<f32, { [B] }>, Token) = unsafe {
                load_ptr_tko(
                    addptr_tile(p_at, at),
                    ordering::Weak,
                    None::<scope::TileBlock>,
                    Some(live),
                    Some(0.0f32),
                    None,
                    Latency::<0>,
                )
            };
            acc = acc + v;
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

/// Whether a CUDA device answers at all (else the CPU runs; tests skip).
pub(crate) fn available() -> bool {
    api::zeros::<f32>(&[B]).sync().is_ok()
}

pub(super) struct Device {
    vb: bool,
    /// class side: chunks of `labels`/`w`, then each class's chunks
    c_lo: Tensor<i32>,
    c_len: Tensor<i32>,
    c_first: Tensor<i32>,
    c_nchunk: Tensor<i32>,
    c_partial: Tensor<f32>,
    count: Tensor<f32>,
    labels: Tensor<i32>,
    w: Tensor<f32>,
    scale: Tensor<f32>,
    /// transcript side, likewise over `t_cls`/`t_w`
    t_lo: Tensor<i32>,
    t_len: Tensor<i32>,
    t_first: Tensor<i32>,
    t_nchunk: Tensor<i32>,
    t_partial: Tensor<f32>,
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
        let (c_lo, c_len) = rows(&prob.c_chunk);
        let (c_first, c_nchunk) = rows(&prob.c_first);
        let (t_lo, t_len) = rows(&prob.t_chunk);
        let (t_first, t_nchunk) = rows(&prob.t_first);
        let (ncc, ntc) = (c_lo.len(), t_lo.len());
        Ok(Self {
            vb,
            c_lo: upload(c_lo)?,
            c_len: upload(c_len)?,
            c_first: upload(c_first)?,
            c_nchunk: upload(c_nchunk)?,
            c_partial: upload(vec![0f32; ncc])?,
            count: upload(prob.count.clone())?,
            labels: upload(i32s(&prob.labels))?,
            w: upload(prob.w.clone())?,
            scale: upload(vec![0f32; prob.num_classes()])?,
            t_lo: upload(t_lo)?,
            t_len: upload(t_len)?,
            t_first: upload(t_first)?,
            t_nchunk: upload(t_nchunk)?,
            t_partial: upload(vec![0f32; ntc])?,
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
            kernels::chunk_dot(
                (&mut self.c_partial).partition([B]),
                &self.c_lo,
                &self.c_len,
                &self.labels,
                &self.w,
                basis,
            )
            .sync()
            .map_err(err)?;
            kernels::class_finish(
                (&mut self.scale).partition([B]),
                &self.c_first,
                &self.c_nchunk,
                &self.c_partial,
                &self.count,
            )
            .sync()
            .map_err(err)?;
            kernels::chunk_dot(
                (&mut self.t_partial).partition([B]),
                &self.t_lo,
                &self.t_len,
                &self.t_cls,
                &self.t_w,
                &self.scale,
            )
            .sync()
            .map_err(err)?;
            kernels::txp_finish(
                dst.partition([B]),
                &self.t_first,
                &self.t_nchunk,
                &self.t_partial,
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

/// The k-mer table on the device (`kmers.rs`), with its lookup kernel in CUDA
/// C++ compiled by NVRTC (cuTile's tile kernels have no per-thread loops).
pub(crate) struct KmerDevice {
    ctx: std::sync::Arc<cuda_core::Context>,
    f: cuda_core::Function,
    keys: cuda_core::DeviceBuffer<u64>,
    vals: cuda_core::DeviceBuffer<u32>,
    kstart: cuda_core::DeviceBuffer<u32>,
    useq: cuda_core::DeviceBuffer<u64>,
    prm: cuda_core::DeviceBuffer<u32>,
}

impl KmerDevice {
    pub(crate) fn new(
        keys: &[u64],
        vals: &[u32],
        kstart: &[u32],
        useq: &[u64],
        prm: &[u32],
    ) -> Result<Self, String> {
        if keys.is_empty() {
            return Err("this backend takes the minimizer table built on the CPU".into());
        }
        let nvrtc = super::dl::open(&[
            "libnvrtc.so",
            "libnvrtc.so.13",
            "libnvrtc.so.12",
            "/usr/local/cuda/lib64/libnvrtc.so",
        ])?;
        // PTX for the lowest architecture cuTile supports; the driver JITs it up.
        let arch = [c"--gpu-architecture=compute_80"];
        let ptx = super::dl::rtc_compile(&nvrtc, "nvrtc", "PTX", &super::dl::kmer_src(), &arch)?;
        let ptx = String::from_utf8_lossy(&ptx);
        let ctx = cuda_core::Context::new(0).map_err(err)?;
        let module = ctx
            .load_module_from_ptx_src(ptx.trim_end_matches('\0'))
            .map_err(err)?;
        let stream = ctx.new_stream().map_err(err)?;
        let up = |v: &[u32]| cuda_core::DeviceBuffer::from_host(&stream, v).map_err(err);
        Ok(Self {
            f: module.load_function("kmer_hits").map_err(err)?,
            keys: cuda_core::DeviceBuffer::from_host(&stream, keys).map_err(err)?,
            vals: up(vals)?,
            kstart: up(kstart)?,
            useq: cuda_core::DeviceBuffer::from_host(&stream, useq).map_err(err)?,
            prm: up(prm)?,
            ctx,
        })
    }

    /// See `KmerGpu::device_bytes`.
    // shortcut: not reported on this backend; sum the uploads if it matters
    pub(crate) fn device_bytes(&self) -> Option<u64> {
        None
    }

    /// See `KmerGpu::set_occurrences`.
    // shortcut: the selection kernel is Metal only (it sums over SIMD
    // groups); here every read replays its selection on the CPU. Port it once
    // a device can test it.
    pub(crate) fn set_occurrences(
        &mut self,
        _off: &[u32],
        _ent: &[u64],
        _prm: &[u32; 3],
    ) -> Result<(), String> {
        Err("the GPU target selection is not wired on this backend".into())
    }

    /// See `KmerGpu::hits`; `sel` is never set here (no `set_occurrences`).
    pub(crate) fn run(
        &self,
        seq: &[u8],
        off: &[u32],
        hits: &mut [u32],
        counts: &mut [u32],
        sel: Option<(&mut [u32], &mut [u32])>,
    ) -> Result<(), String> {
        if sel.is_some() {
            return Err("no occurrence lists on the device".into());
        }
        self.hits(seq, off, hits, counts)
    }

    /// See `KmerGpu::hits`; `hits` and `counts` are already sized. Each call
    /// has its own stream and buffers, so mapping threads overlap.
    pub(crate) fn hits(
        &self,
        seq: &[u8],
        off: &[u32],
        hits: &mut [u32],
        counts: &mut [u32],
    ) -> Result<(), String> {
        let stream = self.ctx.new_stream().map_err(err)?;
        // a zero-length upload is no allocation; the kernel never reads it
        let seq = if seq.is_empty() { &[0u8][..] } else { seq };
        let seq_d = cuda_core::DeviceBuffer::from_host(&stream, seq).map_err(err)?;
        let off_d = cuda_core::DeviceBuffer::from_host(&stream, off).map_err(err)?;
        let hits_d = cuda_core::DeviceBuffer::<u32>::zeroed(&stream, hits.len()).map_err(err)?;
        let counts_d =
            cuda_core::DeviceBuffer::<u32>::zeroed(&stream, counts.len()).map_err(err)?;
        let mut ptrs = [
            seq_d.cu_deviceptr(),
            off_d.cu_deviceptr(),
            self.keys.cu_deviceptr(),
            self.vals.cu_deviceptr(),
            self.kstart.cu_deviceptr(),
            self.useq.cu_deviceptr(),
            self.prm.cu_deviceptr(),
            hits_d.cu_deviceptr(),
            counts_d.cu_deviceptr(),
        ];
        let mut n = counts.len() as u32;
        let mut params: Vec<*mut std::ffi::c_void> = ptrs
            .iter_mut()
            .map(|p| std::ptr::from_mut(p).cast())
            .collect();
        params.push((&mut n as *mut u32).cast());
        // SAFETY: one pointer per `kmer_hits` parameter, in order and of its
        // type (eight device pointers, then the count); the buffers outlive
        // the stream's work, which `copy_to_host` waits for.
        unsafe {
            cuda_core::launch_kernel_on_stream(
                &self.f,
                (n.div_ceil(256), 1, 1),
                (256, 1, 1),
                0,
                &stream,
                &mut params,
            )
        }
        .map_err(err)?;
        hits_d.copy_to_host(&stream, hits).map_err(err)?;
        counts_d.copy_to_host(&stream, counts).map_err(err)
    }
}
