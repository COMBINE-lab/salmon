//! wgpu/WGSL backend: runs the banded-DP kernel on the GPU (Metal on Apple,
//! Vulkan on Linux/NVIDIA) and returns scores bit-identical to the CPU
//! [`reference`](crate::reference).
//!
//! One GPU invocation scores one alignment; the kernel ([`banded_dp.wgsl`]) is a
//! transliteration of [`crate::reference::run_dp`]. A whole mini-batch of
//! candidate alignments is packed into a few storage buffers and dispatched at
//! once, which is what makes the GPU worthwhile despite each individual banded
//! DP being tiny.
//!
//! [`banded_dp.wgsl`]: ../../src/banded_dp.wgsl

use wgpu::util::DeviceExt;

use crate::reference::{banded_extz_score_dna5, BandedParams};

/// Bandwidths above this are scored on the CPU reference: the kernel stores each
/// DP row's band in fixed-size private arrays sized for this maximum. Must match
/// `MAX_W` in the shader. (salmon's default `--dpBandwidth` is 15.)
const MAX_W: i32 = 40;

const WORKGROUP_SIZE: u32 = 64;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct ParamsU {
    match_score: i32,
    mismatch_pen: i32,
    gap_open_pen: i32,
    gap_extend_pen: i32,
    bandwidth: i32,
    n_tasks: u32,
    _pad0: u32,
    _pad1: u32,
}

/// One alignment to score: the query and target as DNA5 codes (0..=4).
#[derive(Clone, Copy)]
pub struct GpuTask<'a> {
    pub query5: &'a [u8],
    pub target5: &'a [u8],
}

/// A ready-to-use GPU device, queue, and compiled banded-DP pipeline. Construct
/// once (it picks an adapter) and reuse across mini-batches.
pub struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}

impl GpuContext {
    /// Acquire a GPU and build the pipeline, or `None` if no adapter is available
    /// (in which case the caller falls back to the CPU aligner).
    pub fn new() -> Option<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))?;

        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("salmon-gpu device"),
                required_features: wgpu::Features::empty(),
                required_limits: adapter.limits(),
                memory_hints: wgpu::MemoryHints::Performance,
            },
            None,
        ))
        .ok()?;

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("banded_dp"),
            source: wgpu::ShaderSource::Wgsl(include_str!("banded_dp.wgsl").into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("banded_dp bgl"),
            entries: &[
                uniform_entry(0),
                storage_entry(1, true),
                storage_entry(2, true),
                storage_entry(3, true),
                storage_entry(4, false),
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("banded_dp layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("banded_dp pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "main",
            compilation_options: Default::default(),
            cache: None,
        });

        Some(Self {
            device,
            queue,
            pipeline,
            bind_group_layout,
        })
    }

    /// Score a batch of alignments on the GPU. Result `i` corresponds to task
    /// `i`. When the bandwidth exceeds [`MAX_W`] the whole batch is scored on the
    /// CPU reference (same algorithm), so the output is always complete and
    /// bit-identical to the reference.
    pub fn batch_align(&self, tasks: &[GpuTask], p: &BandedParams) -> Vec<i32> {
        let n = tasks.len();
        if n == 0 {
            return Vec::new();
        }
        // The kernel's band storage is sized for bandwidths up to MAX_W; wider
        // bands (rare; salmon's default is 15) are scored on the CPU reference.
        if p.bandwidth > MAX_W {
            return tasks
                .iter()
                .map(|t| banded_extz_score_dna5(t.query5, t.target5, p))
                .collect();
        }

        // Pack queries and targets into flat one-base-per-u32 buffers, recording
        // each task's offsets and lengths.
        let mut qbuf: Vec<u32> = Vec::new();
        let mut tbuf: Vec<u32> = Vec::new();
        let mut meta: Vec<[u32; 4]> = Vec::with_capacity(n);
        for t in tasks {
            let q_off = qbuf.len() as u32;
            let t_off = tbuf.len() as u32;
            qbuf.extend(t.query5.iter().map(|&b| b as u32));
            tbuf.extend(t.target5.iter().map(|&b| b as u32));
            meta.push([q_off, t.query5.len() as u32, t_off, t.target5.len() as u32]);
        }
        // Storage buffers must be non-empty.
        if qbuf.is_empty() {
            qbuf.push(0);
        }
        if tbuf.is_empty() {
            tbuf.push(0);
        }

        let params = ParamsU {
            match_score: p.match_score,
            mismatch_pen: p.mismatch_pen,
            gap_open_pen: p.gap_open_pen,
            gap_extend_pen: p.gap_extend_pen,
            bandwidth: p.bandwidth,
            n_tasks: n as u32,
            _pad0: 0,
            _pad1: 0,
        };

        let params_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("params"),
                contents: bytemuck::bytes_of(&params),
                usage: wgpu::BufferUsages::UNIFORM,
            });
        let meta_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("tasks"),
                contents: bytemuck::cast_slice(&meta),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let qbuf_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("qbuf"),
                contents: bytemuck::cast_slice(&qbuf),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let tbuf_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("tbuf"),
                contents: bytemuck::cast_slice(&tbuf),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let out_size = (n * std::mem::size_of::<i32>()) as u64;
        let out_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("scores"),
            size: out_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: out_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("banded_dp bg"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: meta_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: qbuf_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: tbuf_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: out_buf.as_entire_binding(),
                },
            ],
        });

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("banded_dp encoder"),
            });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("banded_dp pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            let groups = (n as u32).div_ceil(WORKGROUP_SIZE);
            pass.dispatch_workgroups(groups, 1, 1);
        }
        encoder.copy_buffer_to_buffer(&out_buf, 0, &readback, 0, out_size);
        self.queue.submit(Some(encoder.finish()));

        // Map and read back the scores.
        let slice = readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv().expect("map_async channel").expect("map failed");
        let scores: Vec<i32> = bytemuck::cast_slice(&slice.get_mapped_range()).to_vec();
        readback.unmap();
        scores
    }
}

fn uniform_entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn storage_entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}
