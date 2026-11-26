# Graph Report - elegant-stonebraker-f1024b  (2026-11-26)

## Corpus Check
- 165 files · ~339,507 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 3293 nodes · 7249 edges · 189 communities (171 shown, 18 thin omitted)
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 430 edges (avg confidence: 0.83)
- Token cost: 101,408 input · 0 output

## Community Hubs (Navigation)
- EM Component Packing
- Atomic Math & FLD CDF
- Adaptive Decode Broker
- Read Processor & RAD Records
- RAD Baked Header Tags
- Genome BAM Projection
- Alignment Scoring Core
- Mate Pairing Candidates
- BAM-to-RAD Grouping
- Sequence Bias Model
- Alignment Error Model
- Deterministic BAM RAD Writer
- RAD Input Quantification
- GC Bias Binning
- MEM Chaining
- Hit Searcher & Unimems
- Docs Workflow CI
- CLI Main & Args
- RAD Reader Concurrency
- MEM Extension & Revcomp
- Uncertainty & Ambiguity
- Debug Trace Tools
- Concurrent Chunk Writer
- Input Validation Tests
- 2-bit Ref Decoding
- Decoy Scoring & RefProvider
- RefView Chain Alignment
- CI Workflow
- Library Type Compatibility
- EqClass Builder
- Packed EqClass Plans
- Windowed Chain Alignment
- Positional Bias Expected
- RAD Hit Encoding
- Read Strandedness Detection
- Collapsed EqClasses Infer
- Read Mapper Assembly
- RAD Roundtrip Tests
- Fragment Source Results
- Reference Index Loading
- Dependabot Config
- map_to_rad Quant API
- DAAREM Accelerator
- Reference Sequences Build
- Diagnostics & Disk Checks
- Bias Config Merge
- CIGAR Mapping Records
- Release 2.12.0 Notes
- Intermediate RAD Cleanup
- map_to_rad Quant API (2)
- Astro Website Config
- Gene Map Aggregation
- Quantmerge
- Gap Score Cache
- Mapping Modes & Citations Docs
- Index cDBG Build
- Bias Placements Grouping
- CLI Subcommands
- Release 2.0/2.1 Features
- Release 2.5/2.11 IO Features
- Orphan & Fragment Probs
- Hash Feed Batching
- MEM Collection
- BAM Compression Workers
- Alignment Scoring Core (2)
- Align Output Contract Tests
- GC Fragment Model
- Spline Band Matrix
- Determinism & Decoy Releases
- EM Step Benchmark
- Alignment Stdin Tests
- RAD Profile Schema
- AlignQuant FLD Policy
- Index Build & Load
- Release 2.9 VBEM & cf1
- Docker Build Cache
- Metadata Truth Tests
- EqClass Builder (2)
- Packed EqClass Plans (2)
- BAM Output Lifecycle
- Release 2.0 Mapping Features
- Baked FLD & Fixed-Point
- Mapping Modes & Citations Docs (2)
- Genome Projection Pipes
- Library Type Compatibility (2)
- EM Determinism Audit
- Gibbs & Bootstrap Replicates
- RAD Reader Concurrency (2)
- Reference Sequences Build (2)
- Reference Digests
- Collapsed EqClasses Infer (2)
- 2.12 API Changes
- Release Workflow Jobs
- Decoy & Readiness Tests
- EM Stop Rule
- SAM Record Writing
- Migration & 2.0 Index
- Salmon 2.0 Rust Rewrite
- CLI Main & Args (2)
- Compressed Input Open
- Bias Dumps
- Release Testing Checklist
- Docker Workflow Env
- Mapping Modes & Citations Docs (3)
- Align Benchmarks
- Global Allocator
- Dump EqClass Tests
- DP Aligner Kernel
- Bias Config Merge (2)
- CLI Flag Status
- Audit Open Issues
- Release Workflow Steps
- Mapping Modes & Citations Docs (4)
- Mapping Modes & Citations Docs (5)
- Skip-Quant Tests
- Inferential Replicates Output
- Release 2.11 Audit
- Workspace Crates
- Progress Reporting
- Bias Correction Apply
- Alignment Input Readers
- Flag Matrix Tests
- CLI Output Contract Tests
- RAD Train Workers
- Component Summary
- Release 2.4 BAM Threads
- DAAREM & Allocator Releases
- SSHash Index v2
- Release Workflow Matrix
- Orphan & Fragment Probs (2)
- Library Type Compatibility (3)
- Index Build & Load (2)
- EmAccel Enum
- Bias Convolution Bench
- Atomic Math & FLD CDF (2)
- Reference Sequences Build (3)
- Index Provenance Types
- Mapping Parity vs C++
- Decode Deadlock & Seeding
- Dovetail Pairing 2.12
- Genome Projection Mode
- Docker Workflow Steps
- Bump & Publish Script
- DP Aligner Kernel (2)
- FLD Policy Arg
- SalmonError
- GTF/GFF Attribute Parsing
- Posterior Method
- EM Thread Scaling Chart
- Release Workflow Env
- Release Host Outputs
- Release Artifact Upload
- EM RAD Bench Script
- Installation Methods
- Library Types & Quick Start
- Mapping Modes & Citations Docs (6)
- Core Placements
- Piscem Bulk Records
- Progress-Aware Logging
- Phase Timing Bench
- Audit Reconciliation 1140
- Release 2.6/2.10 Fixes
- SSHash Bucket Cache
- Docker Jobs
- Docker Platform Matrix
- Docker Secrets Env
- Release Workflow Env (2)
- Mapping Modes & Citations Docs (7)
- Mapping Modes & Citations Docs (8)
- Refactor Bookmark
- EffLens Newtype
- BGZF Worker Count
- RAD Input Profile
- Release 2.1.2 Piscem Fix
- Release Announcement Env
- Website TSConfig
- AlignQuant FLD Policy (2)
- Bias Config Merge (3)
- Gene Quant Summary
- Reference Sequences Build (4)
- Packed EqClass Plans (3)
- Community Guidelines
- Progress Reporting (2)
- BAM Compression Workers (2)
- CigarKind
- Reference Sequences Build (5)
- Reference Sequences Build (6)
- Decode Crossover Script
- Profile Script

## God Nodes (most connected - your core abstractions)
1. `quantify_rad()` - 65 edges
2. `LibraryFormat` - 55 edges
3. `write_alignment_rad()` - 42 edges
4. `MateStatus` - 40 edges
5. `SalmonIndex` - 40 edges
6. `MappingCandidate` - 32 edges
7. `AlignConfig` - 32 edges
8. `Mem` - 31 edges
9. `RadOutputWriter` - 31 edges
10. `ScoredMapping` - 31 edges

## Surprising Connections (you probably didn't know these)
- `Rejected: shared AtomicF64 accumulator / cached alpha*w scratch / parallel VBEM sum` --semantically_similar_to--> `rapidgzip DecoderPool::release lost wakeup (fixed rapidgzip-core 0.3.1)`  [AMBIGUOUS] [semantically similar]
  docs/em-architecture-and-optimization.md → notes/pool-lost-wakeup-deadlock.md
- `Flag status legend (Implemented/Default/Accept+warn/Rejected/Covered)` --semantically_similar_to--> `Accepted-but-ignored options (parse + warn)`  [INFERRED] [semantically similar]
  docs/cli-flag-status.md → MIGRATION.md
- `Selective-alignment chain pruning defaults off (orphan/postMergeChainSubThresh = 0.0)` --semantically_similar_to--> `Chain sub-optimality threshold default difference (~80% of SA gap)`  [INFERRED] [semantically similar]
  MIGRATION.md → docs/mapping-parity-differences.md
- `Coordinate-sorted BAM rejected up front` --semantically_similar_to--> `Name-grouped BAM requirement`  [INFERRED] [semantically similar]
  docs/release-testing-checklist.md → website/src/content/docs/guides/genome-projection.mdx
- `collect_bias_fragment()` --references--> `groups`  [INFERRED]
  crates/salmon-align/src/rad.rs → .github/dependabot.yml

## Import Cycles
- 1-file cycle: `crates/salmon-align/src/lib.rs -> crates/salmon-align/src/lib.rs`
- 2-file cycle: `crates/salmon-map/src/chain.rs -> crates/salmon-map/src/mem.rs -> crates/salmon-map/src/chain.rs`
- 2-file cycle: `crates/salmon-align/src/lib.rs -> crates/salmon-align/src/rad.rs -> crates/salmon-align/src/lib.rs`
- 2-file cycle: `crates/salmon-align/src/genome_project.rs -> crates/salmon-align/src/lib.rs -> crates/salmon-align/src/genome_project.rs`
- 2-file cycle: `crates/salmon-align/src/bam_rad.rs -> crates/salmon-align/src/lib.rs -> crates/salmon-align/src/bam_rad.rs`
- 2-file cycle: `crates/salmon-infer/src/lib.rs -> crates/salmon-infer/src/uncertainty.rs -> crates/salmon-infer/src/lib.rs`
- 2-file cycle: `crates/salmon-infer/src/components.rs -> crates/salmon-infer/src/lib.rs -> crates/salmon-infer/src/components.rs`
- 2-file cycle: `crates/salmon-rad/src/lib.rs -> crates/salmon-rad/src/writer.rs -> crates/salmon-rad/src/lib.rs`
- 2-file cycle: `crates/salmon-rad/src/lib.rs -> crates/salmon-rad/src/record.rs -> crates/salmon-rad/src/lib.rs`
- 3-file cycle: `crates/salmon-rad/src/lib.rs -> crates/salmon-rad/src/writer.rs -> crates/salmon-rad/src/schema.rs -> crates/salmon-rad/src/lib.rs`
- 3-file cycle: `crates/salmon-rad/src/lib.rs -> crates/salmon-rad/src/writer.rs -> crates/salmon-rad/src/record.rs -> crates/salmon-rad/src/lib.rs`

## Hyperedges (group relationships)
- **2.1.0 decoy correctness fixes** — docs_release_notes_2_1_0_decoy_boundary_fix, docs_release_notes_2_1_0_n_aware_decoy_indexing, docs_release_notes_2_1_0_sketch_decoy_handling, docs_release_notes_2_1_0_allow_decoy_orphans, docs_release_notes_2_1_0_index_version [EXTRACTED 1.00]
- **EM thread-scaling benchmark: baseline vs candidate across two hosts** — docs_benchmarks_em_thread_scaling_scaling_baseline_1169, docs_benchmarks_em_thread_scaling_scaling_selected_128_shard_candidate, docs_benchmarks_em_thread_scaling_scaling_local_epyc_9575f_host, docs_benchmarks_em_thread_scaling_scaling_newton_xeon_e5_2699_v4_host, docs_benchmarks_em_thread_scaling_scaling_speedup_t1_tp [EXTRACTED 1.00]
- **Jointly required deadlock conditions (adaptive broker, paired locks, lost wakeup)** — notes_pool_lost_wakeup_deadlock_adaptive_broker, notes_pool_lost_wakeup_deadlock_paraseq_paired_reader, notes_pool_lost_wakeup_deadlock_decoderpool_lost_wakeup [EXTRACTED 1.00]
- **RAD provenance carried into meta_info completeness** — website_src_content_docs_reference_output_formats_rad_output, website_src_content_docs_reference_output_formats_baked_header_tags, website_src_content_docs_reference_output_formats_alignment_provenance, website_src_content_docs_reference_output_formats_metadata_completeness, docs_release_notes_2_6_0_meta_info_json [EXTRACTED 1.00]
- **Causes of SA mapping-rate gap vs C++ salmon** — docs_mapping_parity_differences_sshash_orientation_bug, docs_mapping_parity_differences_chain_sub_optimality_default, docs_mapping_parity_differences_tie_break_residual, docs_mapping_parity_differences_polya_clipping_mismatch, docs_mapping_parity_differences_seeding_sensitivity_tail [EXTRACTED 1.00]
- **Small dovetail pairing change (#1356/#1357)** — docs_release_notes_2_12_0_small_dovetail_pairing, docs_release_notes_2_12_0_five_prime_end_fragment_rule, docs_release_notes_2_12_0_insert_only_scoring, docs_release_notes_2_12_0_dovetail_tolerance_flag, docs_release_notes_2_12_0_allow_dovetail_inward_only [EXTRACTED 1.00]
- **Adaptive threading stack** — docs_release_notes_2_5_0_thread_broker, docs_release_notes_2_5_0_paraseq_temp, docs_release_notes_2_5_1_bounded_worker_pool, docs_release_notes_2_5_0_rapidgzip_core [EXTRACTED 1.00]
- **Fragment-length distribution modeling** — docs_release_notes_2_0_1_logfragprob_eq_weights, docs_release_notes_2_0_1_abundance_aware_fld_training, docs_release_notes_2_1_0_duplicate_symmetry_fld_snapshot, docs_release_notes_2_1_0_ambiguous_frag_length_prob, docs_release_notes_2_4_0_fld_policy [INFERRED 0.75]
- **sshash lookup speedups inherited via piscem** — docs_sshash_bucket_cache_bucket_bounds_cache, docs_sshash_v6_port_canonical_mmer_minimizer, docs_sshash_bucket_cache_piscem_streaming_query, docs_release_notes_2_12_0_sshash_lib_0_8_2 [INFERRED 0.75]
- **Cross-thread-count EM determinism (data-derived shards, class ordering, invariants, V1 fix)** — docs_em_architecture_and_optimization_shardplan, docs_em_architecture_and_optimization_determinism_invariants, docs_em_architecture_and_optimization_range_factorized_class_ordering, docs_audit_1140_reconciliation_v1_thread_count_reproducibility, docs_em_architecture_and_optimization_max_rel_diff [INFERRED 0.85]
- **EM optimizers and accelerators** — docs_release_notes_2_8_0_vbem, docs_release_notes_2_8_0_squarem, docs_release_notes_2_8_0_daarem_rework, docs_release_notes_2_9_0_component_em_engine [INFERRED 0.85]
- **Index format evolution v2 to v3** — docs_release_notes_2_7_0_index_format_v2, docs_release_notes_2_7_0_sshash_v6, docs_release_notes_2_9_0_index_format_v3, docs_release_notes_2_7_0_sshash_rs [INFERRED 0.85]
- **Inferential replicate generation, storage and downstream use** — website_src_content_docs_reference_cli_numbootstraps, website_src_content_docs_reference_cli_numgibbssamples, website_src_content_docs_guides_inferential_replicates_bootstrap_replicates, website_src_content_docs_guides_inferential_replicates_gibbs_sampling, website_src_content_docs_reference_output_formats_bootstrap_dir, docs_release_notes_2_6_0_meta_info_json, website_src_content_docs_guides_inferential_replicates_tximport, website_src_content_docs_guides_inferential_replicates_fishpond_swish [INFERRED 0.85]
- **RAD-based decoupled and deterministic quantification** — docs_release_notes_2_2_0_rad_decoupled_quant, docs_release_notes_2_2_0_deterministic_mode, docs_release_notes_2_2_0_baked_fld_header, docs_release_notes_2_2_0_fixed_point_u128_accumulation, docs_release_notes_2_2_0_rad_chunk_compression, docs_release_notes_2_4_0_fld_policy [INFERRED 0.85]
- **Deterministic two-phase RAD quantification pipeline** — website_src_content_docs_guides_rad_and_determinism_writerad, website_src_content_docs_guides_rad_and_determinism_rad_format, website_src_content_docs_guides_rad_and_determinism_fldpolicy, website_src_content_docs_guides_rad_and_determinism_deterministic_quant, website_src_content_docs_guides_rad_and_determinism_fixed_point_accumulation [INFERRED 0.85]

## Communities (189 total, 18 thin omitted)

### Community 0 - "EM Component Packing"
Cohesion: 0.07
Nodes (61): a_capped_global_run_blames_only_the_unsettled_component(), a_gathered_step_of_a_large_component_has_the_same_l1_bits(), a_global_run_that_stops_at_its_minimum_while_drifting_is_flagged(), BATCH_INCIDENCES, bits(), Buffers, CompEnd, ComponentPlan (+53 more)

### Community 1 - "Atomic Math & FLD CDF"
Cohesion: 0.05
Nodes (32): ExplicitFldArgs, warn_baked_fld_overrides(), AtomicF64, add_identities(), add_matches_naive(), approx(), LOG_0, LOG_1 (+24 more)

### Community 2 - "Adaptive Decode Broker"
Cohesion: 0.06
Nodes (26): default_engagement_policy(), engagement_policy(), Input, MappingConsumer, MappingPlan, plan(), sketch_uses_the_value_measured_for_salmon_not_piscem(), t() (+18 more)

### Community 3 - "Read Processor & RAD Records"
Cohesion: 0.08
Nodes (24): is_compatible(), observed_paired_format(), accumulate_vm_stats(), build_rad_record(), Counters, flush_bam(), flush_processed(), flush_rad() (+16 more)

### Community 4 - "RAD Baked Header Tags"
Cohesion: 0.05
Nodes (44): BAKED_FLAGS_TAG, BAKED_INDEX_PROV, BAKED_MAP_COUNTERS, clip_pair_start(), clip_pair_start_keeps_the_end(), FIRST_DECOY_INDEX_TAG, for_placement_frag_len_convention(), FRAG_LEN_UNPAIRED (+36 more)

### Community 5 - "Genome BAM Projection"
Cohesion: 0.09
Nodes (26): are_mates(), best_scoring_true_pair_wins(), build_projected_record(), CHUNK_FLUSH_BYTES, ga(), genome_alignment_score(), GenomeProjectOptions, kind_to_sam_code() (+18 more)

### Community 6 - "Alignment Scoring Core"
Cohesion: 0.07
Nodes (35): a_wrong_length_is_a_mismatch_but_absent_or_a_clipped_tail_is_not(), alignment_score_value_decoding(), chain_fits(), coordinate_sorted_unusable(), empty_reference_fasta_is_an_error(), escape_header_value(), fill_missing_bases(), fld_policy_defaults_to_baked() (+27 more)

### Community 7 - "Mate Pairing Candidates"
Cohesion: 0.13
Nodes (40): MappingCandidate, allowed_dovetail_fragment_len_is_five_prime_span(), cand(), DEFAULT_DOVETAIL_TOLERANCE, dovetail_filtered_by_default_but_kept_when_enabled(), dovetail_rejected(), emit_best_orphans(), emits_all_coverage_tied_pairs_per_target() (+32 more)

### Community 8 - "BAM-to-RAD Grouping"
Cohesion: 0.09
Nodes (22): AlignRadSummary, analyze_group(), CHUNK_FLUSH_BYTES, ERROR_MODEL_ALPHA, GroupData, GroupWorker, MINIBATCH, reduce_summary() (+14 more)

### Community 9 - "Sequence Bias Model"
Cohesion: 0.10
Nodes (25): base2bit(), bench_factor_build(), build_expected(), build_expected_respects_num_targets_bound(), complement_bit(), CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT (+17 more)

### Community 10 - "Alignment Error Model"
Cohesion: 0.08
Nodes (27): ScoreWriteWorker<'a>, AlignmentModel, ALN_DASH, ALN_HARD_CLIP, ALN_PAD, ALN_REF_SKIP, ALN_SOFT_CLIP, AlnOp (+19 more)

### Community 11 - "Deterministic BAM RAD Writer"
Cohesion: 0.13
Nodes (39): write_alignment_rad(), a_constant_as_is_reported_too(), a_partially_tagged_bam_is_reported(), a_placement_on_a_reference_missing_from_t_does_not_win(), a_secondary_without_bases_is_scored_from_its_primary(), a_wrong_explicit_lib_type_is_diagnosed(), an_as_less_bam_is_reported(), an_orphan_beside_a_proper_pair_takes_no_mass() (+31 more)

### Community 12 - "RAD Input Quantification"
Cohesion: 0.13
Nodes (31): quantify_rad(), all_bias_rad_is_thread_invariant(), auto_orientation_derives_fld_from_majority_bucket(), baked_fld_is_reported_as_the_frag_length_source(), baked_fld_is_used_not_derived(), baked_library_format_is_authoritative_under_auto(), bias_ref_seq_fixture(), compressed_rad_quantifies_identically() (+23 more)

### Community 13 - "GC Bias Binning"
Cohesion: 0.10
Nodes (30): ctx_frac_table(), DEFAULT_COND_BINS, DEFAULT_GC_BINS, frag_frac_table(), gc_context_matches_gc_desc(), gc_context_tables_match_desc(), gc_desc(), gc_desc_context_window_geometry() (+22 more)

### Community 14 - "MEM Chaining"
Cohesion: 0.11
Nodes (22): clip_mems(), anti_colinear_mems_split(), cfg(), chain_mems(), chain_two_mems(), ChainConfig, chains_sorted_by_descending_score(), ChainScratch (+14 more)

### Community 15 - "Hit Searcher & Unimems"
Cohesion: 0.15
Nodes (30): consensus_filter(), collect_candidates(), debug_best_mapping(), DebugMapping, map_read_pair_into(), map_single_read(), map_single_read_into(), MapConfig (+22 more)

### Community 16 - "Docs Workflow CI"
Cohesion: 0.06
Nodes (37): defaults, runs-on, steps, concurrency, cancel-in-progress, group, run, environment (+29 more)

### Community 17 - "CLI Main & Args"
Cohesion: 0.07
Nodes (19): DEFAULT_FLD_MAX, DEFAULT_FLD_MEAN, DEFAULT_FLD_SD, deterministic_requant_forwards_every_supported_knob(), deterministic_requant_leaves_unset_knobs_at_their_defaults(), GLOBAL_RULE_MIN_ITER, inert_flag_reason_agrees_with_how_many_were_passed(), inert_in_mode_message() (+11 more)

### Community 18 - "RAD Reader Concurrency"
Cohesion: 0.09
Nodes (20): arrival_order_within_a_group_is_preserved(), BIAS_MIN_ALPHA, degenerate_inputs(), duplicate_tids_combine_by_logsumexp(), grouped(), load_baked_abundances(), load_baked_fld(), load_baked_library_format() (+12 more)

### Community 19 - "MEM Extension & Revcomp"
Cohesion: 0.12
Nodes (23): revcomp_into(), base_eq(), candidates_from_raw_hits_true_unimems(), candidates_from_raw_hits_unimems(), collect_read_true_unimems(), collect_read_unimems(), extend_grows_seed_to_full_exact_match(), extend_mem() (+15 more)

### Community 20 - "Uncertainty & Ambiguity"
Cohesion: 0.13
Nodes (24): ambiguity_counts(), ambiguity_counts_split(), block_rng(), bootstrap(), bootstrap_mean_near_point_estimate(), bootstrap_replicates_follow_the_per_nucleotide_prior(), bootstrap_variance_grows_with_ambiguity(), draw_class() (+16 more)

### Community 21 - "Debug Trace Tools"
Cohesion: 0.14
Nodes (10): dump_kmers(), main(), MateTrace, trace_mate(), ScoredMapping, map_read_pair_sketch(), map_read_pair_sketch_into(), map_single_read_sketch() (+2 more)

### Community 22 - "Concurrent Chunk Writer"
Cohesion: 0.08
Nodes (10): BAKED_ABUND, BAKED_FLD, BAKED_LIBFMT, BAKED_SCORE_KIND, MapCounters, FragmentChunkBuf, LOW_SPACE_BYTES, partial_path_for() (+2 more)

### Community 23 - "Input Validation Tests"
Cohesion: 0.19
Nodes (32): a_failed_rerun_does_not_leave_the_previous_quant_sf(), a_fifo_input_quantifies_like_the_file_it_carries(), a_gzip_stream_cut_short_is_refused(), a_plain_stream_cut_mid_record_is_refused_and_named(), a_run_that_maps_nothing_fails_but_keeps_its_output(), a_truncated_fastq_is_refused_rather_than_partially_quantified(), ambig_info_has_one_row_per_quant_row_with_decoys(), an_impossible_k_or_m_is_refused_before_writing_anything() (+24 more)

### Community 24 - "2-bit Ref Decoding"
Cohesion: 0.12
Nodes (24): contested_sam(), counts(), discard_orphans_takes_effect_on_rad_input(), hard_filter_takes_effect_on_alignment_input(), orphan_sam(), os(), quant(), reads_mapping_flags_warn_in_alignment_and_rad_modes() (+16 more)

### Community 25 - "Decoy Scoring & RefProvider"
Cohesion: 0.14
Nodes (23): RefProvider, both_dedup_paths_agree_including_on_ties(), decoy_below_threshold_is_ignored(), decoy_domination_drops_read(), DecoyAt, dedup_by_table(), DEDUP_HASH_THRESHOLD, dedup_hasher() (+15 more)

### Community 26 - "RefView Chain Alignment"
Cohesion: 0.12
Nodes (27): RefView, Bytes, Packed, align_chain_above(), align_chain_clipped_view(), align_chain_exact(), align_chain_view(), AlignConfig (+19 more)

### Community 27 - "CI Workflow"
Cohesion: 0.06
Nodes (32): runs-on, steps, strategy, timeout-minutes, env, CARGO_TERM_COLOR, RUST_TOOLCHAIN, RUSTDOCFLAGS (+24 more)

### Community 28 - "Library Type Compatibility"
Cohesion: 0.11
Nodes (15): compatible_paired(), compatible_single(), format_ids_are_unique_and_dense(), lib_format_summary_matches_cpp_derivation(), LibFormatCountsArray, LibFormatCountsFile, LibFormatSummary, NUM_LIB_FORMATS (+7 more)

### Community 29 - "EqClass Builder"
Cohesion: 0.12
Nodes (16): bins_split_classes_with_same_transcripts(), builder_aggregates_counts_and_weights(), bytemuck_cast(), combined_weights_divide_by_eff_length_and_normalize(), eq_class_dump_collapses_only_without_weights(), finish_orders_the_complete_range_factorized_label(), group_normalizes_and_hashes_consistently(), NAIVE_NO_FMT (+8 more)

### Community 30 - "Packed EqClass Plans"
Cohesion: 0.15
Nodes (21): class_log_term(), compressed_plan_indices_and_ranges_are_consistent(), compressed_steps_match_dense_reference_for_the_same_plan(), Contributor, dense_em_reference(), dense_vbem_reference(), DisjointOut, DisjointOut<'a> (+13 more)

### Community 31 - "Windowed Chain Alignment"
Cohesion: 0.20
Nodes (28): align_chain(), align_in_window(), anchored_and_full_agree_on_clean_read(), anchored_scores_two_mem_chain_with_exact_gap(), BASE_CODE, clipped_read_through_scores_only_insert_bases(), comp(), DNA5_LUT (+20 more)

### Community 32 - "Positional Bias Expected"
Cohesion: 0.13
Nodes (23): BiasInputs, build_expected_pos(), build_expected_pos_respects_num_targets_bound(), corrected_effective_length_full(), EPSILON, mass_diff(), POS_SMOOTH_FRAC, positional_factor() (+15 more)

### Community 33 - "RAD Hit Encoding"
Cohesion: 0.10
Nodes (5): from_mate_status(), ori_ref_roundtrip(), RadHit, SalmonBulkContext, SalmonBulkRecord

### Community 34 - "Read Strandedness Detection"
Cohesion: 0.15
Nodes (15): ReadType, PairedEnd, SingleEnd, DEFAULT_SAMPLES_NEEDED, final_format_without_lockin_infers_from_partial(), infer_format_from_counts(), LibraryTypeDetector, paired_end_detects_isr() (+7 more)

### Community 35 - "Collapsed EqClasses Infer"
Cohesion: 0.19
Nodes (24): ambiguous(), build(), conserves_total_count(), daarem_conserves_total_count(), daarem_matches_plain_em_fixpoint(), daarem_objective_is_monotone_under_plain_steps(), daarem_reduces_iters_on_slow_mixer(), daarem_with_vbem_is_rejected() (+16 more)

### Community 36 - "Read Mapper Assembly"
Cohesion: 0.23
Nodes (22): best_per_target(), map_read_pair(), assembled_paired_read_mapper_concordant(), assembled_single_read_mapper(), break_kmers(), build_index_for(), collected_candidate_validates_by_alignment(), dovetail_pair_rejected_by_default_accepted_with_flag() (+14 more)

### Community 37 - "RAD Roundtrip Tests"
Cohesion: 0.18
Nodes (19): abandoned_write_never_claims_the_final_path(), check_roundtrip(), chunk_codec_tag_advertised(), concurrent_writer_counts_chunks(), finalize_renames_and_marks_complete(), mapping_pass_provenance_round_trips(), mixed_orphan_sides_are_marked_in_the_encoded_hit(), read_all() (+11 more)

### Community 38 - "Fragment Source Results"
Cohesion: 0.14
Nodes (14): AlignQuantOptions, AlignQuantResult, all_bias_alignment_mode_smoke(), bias_ignores_fasta_polya_tail_beyond_header_length(), quantify_two_phase(), write_fixture(), align_auto_library_type_detects_from_alignments(), quantify_two_phase() (+6 more)

### Community 39 - "Reference Index Loading"
Cohesion: 0.09
Nodes (9): IndexInfo, PackedRefseq, read_info(), RefseqLoad, Bytes, None, Packed, SalmonIndex (+1 more)

### Community 40 - "Dependabot Config"
Cohesion: 0.08
Nodes (25): patterns, update-types, exclude-patterns, patterns, update-types, patterns, update-types, prefix (+17 more)

### Community 41 - "map_to_rad Quant API"
Cohesion: 0.10
Nodes (3): map_to_rad(), QuantOptions, deterministic_mapping_pass_reports_incompatible_fragments()

### Community 42 - "DAAREM Accelerator"
Cohesion: 0.12
Nodes (23): A1, chunked_mut(), CYCL_MON_TOL, daarem_em_stops_as_close_as_plain_em(), daarem_iterates_stay_non_negative(), daarem_loop(), damping_find(), gram_and_rhs() (+15 more)

### Community 43 - "Reference Sequences Build"
Cohesion: 0.22
Nodes (20): bam_content_is_independent_of_thread_count(), bam_write_failure_is_reported(), counts_by_name(), gen_seq(), mates_orphaning_to_different_transcripts_are_one_orphan_fragment(), pseudoalignment_quantification_tracks_truth(), pseudoalignment_strand_filter_respects_library_type(), quantify_two_phase() (+12 more)

### Community 44 - "Diagnostics & Disk Checks"
Cohesion: 0.12
Nodes (13): Diagnostic, EmComponentsMeta, free_disk_bytes(), input_diagnostics(), MetaFieldSource, Bam, Index, Rad (+5 more)

### Community 45 - "Bias Config Merge"
Cohesion: 0.14
Nodes (12): BiasCfg, BiasLocal, derive_fld(), LocalTallies, RadTallies, run_bias_pass(), run_rad_eq_and_bias_pass(), run_rad_pass() (+4 more)

### Community 46 - "CIGAR Mapping Records"
Cohesion: 0.15
Nodes (20): Cigar, cigar_string(), CigarOp, emit_fragment_records(), IS_RC, mapping_quality(), mapping_quality_follows_the_star_convention(), MATE_RC (+12 more)

### Community 47 - "Release 2.12.0 Notes"
Cohesion: 0.10
Nodes (21): -a reads SAM/BAM from a pipe, Outputs written via temp name and rename, Audit of 2.11.0 (#1314), bramble-rs 0.2.3, cf1-rs 0.6.1, salmon 2.12.0 release notes, --errorModel with partial -t fix, FASTQ read in batches of 4096 (+13 more)

### Community 48 - "Intermediate RAD Cleanup"
Cohesion: 0.18
Nodes (18): cleanup_intermediate_rad(), cleanup_intermediate_rad_removes_file_and_partials(), gene_map_is_validated_before_the_run(), GeneMapOpts, intermediate_rad_path(), intermediate_rad_path_placement(), load_gene_map(), long_read_redirect() (+10 more)

### Community 49 - "map_to_rad Quant API (2)"
Cohesion: 0.18
Nodes (17): MapResult, clamps_out_of_range_decoy_block(), CmdInfo, decoys_as_clean_suffix_no_shorts(), decoys_then_shorts_skips_decoys_keeps_shorts(), legacy_index_without_num_decoys_drops_suffix(), no_decoys_emits_all(), paths_to_strings() (+9 more)

### Community 50 - "Astro Website Config"
Cohesion: 0.11
Nodes (18): astro, @astrojs/starlight, sharp, dependencies, astro@^7.4.7, @astrojs/starlight@^0.42.9, sharp@^0.35.5, name (+10 more)

### Community 51 - "Gene Map Aggregation"
Cohesion: 0.20
Nodes (17): a_fallback_row_colliding_with_a_real_gene_name_merges_into_it(), a_stale_unmatched_list_is_removed_when_everything_matches(), aggregation_sums_and_weights(), full_match_does_not_trip_the_low_rate_warning(), ignore_tx_version_matches_versioned_names_to_an_unversioned_map(), ignore_tx_version_resolves_a_collision_independently_of_hash_order(), ignore_tx_version_strips_the_gene_map_side_as_well(), LOW_MATCH_RATE (+9 more)

### Community 52 - "Quantmerge"
Cohesion: 0.13
Nodes (13): a_short_row_is_refused_and_named(), a_target_repeated_in_one_sample_is_refused(), duplicate_or_empty_sample_names_are_refused(), HEADER, merge_in(), merge_quants(), MergeColumn, Elen (+5 more)

### Community 53 - "Gap Score Cache"
Cohesion: 0.18
Nodes (14): cached_flank_score(), entries_survive_rehashing_on_growth(), every_dp_kind_is_distinct_for_the_same_bytes(), gap_hash(), GapCache, GapEntry, GapKind, FlankLeftAnchored (+6 more)

### Community 54 - "Mapping Modes & Citations Docs"
Cohesion: 0.12
Nodes (15): --allowDovetail, Selective alignment & sketch mode guide, Bray et al. 2016 - Near-optimal probabilistic RNA-seq quantification (kallisto), Suzuki & Kasahara 2018 - difference recurrence relations (ksw2), ksw2rs (native-Rust port of ksw2), MEM finding and minimap2-style chaining, Li 2018 - Minimap2 (Bioinformatics), --minScoreFraction (default 0.65) (+7 more)

### Community 55 - "Index cDBG Build"
Cohesion: 0.11
Nodes (13): CDBG_PREFIX, INDEX_FORMAT_VERSION, INDEX_PREFIX, INFO_FILE, MAX_DECODE_THREADS, MIN_READABLE_INDEX_VERSION, POLYA_CLIP_LEN, PreprocessResult (+5 more)

### Community 56 - "Bias Placements Grouping"
Cohesion: 0.21
Nodes (10): bias_placements_group_and_keep_their_geometry(), BiasPlacement, core(), grouping_an_ordered_run_does_not_reorder_it(), ordered(), push_and_push_sorted_agree_on_ordered_input(), push_sorted_rejects_an_out_of_order_placement(), the_sorted_flag_is_derived_from_what_was_pushed() (+2 more)

### Community 57 - "CLI Subcommands"
Cohesion: 0.11
Nodes (19): AlevinArgs, Cli, Command, Alevin, DebugMap, Index, Quant, Quantmerge (+11 more)

### Community 58 - "Release 2.0/2.1 Features"
Cohesion: 0.13
Nodes (20): 2.0 index format (incompatible with pufferfish), --sketch alignment-free mode, --allowDecoyOrphans, cf1-rs (cDBG builder), Concordant-pairing correctness (joinReadsAndFilter), Decoy boundary fix (first_decoy_index, num_decoys, #1019), duplicate_clusters.tsv under --keepDuplicates (#1015), Duplicate-transcript symmetry via per-fragment FLD snapshot (+12 more)

### Community 59 - "Release 2.5/2.11 IO Features"
Cohesion: 0.14
Nodes (18): FIFO/process-substitution input (#1269), In-reader truncated FASTQ detection, Content-based compression detection (gzip/bzip2/xz/zstd), --decoder {auto,serial,parallel,parallel=N}, --geneMap gzip GTF/GFF and early validation, Packed CSR equivalence-class layout (u64 offsets), paraseq-temp (resizable worker pool fork of paraseq), Mapping/quant performance campaign (per-worker tallies, scratch reuse) (+10 more)

### Community 60 - "Orphan & Fragment Probs"
Cohesion: 0.25
Nodes (19): best_kept_score(), collect_bias_fragment(), hit(), is_orphan(), no_single_frag_prob_penalizes_only_orphans(), piscem_unpaired_hit_has_no_read_length(), PlacementOrder, Arrival (+11 more)

### Community 61 - "Hash Feed Batching"
Cohesion: 0.16
Nodes (10): Batch, HashFeed, HashItem, Kept, Own, is_pure_acgt(), longest_acgt_run(), open_fasta() (+2 more)

### Community 62 - "MEM Collection"
Cohesion: 0.22
Nodes (10): best_per_target_is_the_same_above_the_hash_threshold(), candidates_from_raw_hits(), collect_read_mems(), forward_kmers_project_and_chain(), fw_table(), MemCollectorConfig, project_raw_hits(), rc_kmers_use_reflected_read_coordinate() (+2 more)

### Community 63 - "BAM Compression Workers"
Cohesion: 0.15
Nodes (16): base_code(), CHUNK_TARGET, compression_worker_count(), compression_workers_are_bounded(), compression_workers_track_the_measured_balance_point(), DEFLATE_MIB_PER_SEC, encode_header(), encode_record() (+8 more)

### Community 64 - "Alignment Scoring Core (2)"
Cohesion: 0.16
Nodes (8): canonical_name(), GroupingCheck, read_errors_name_the_file(), ReadGrouper, ReadGrouper<R>, warn_mass_on_unreported_refs(), write_outputs(), asctime_now()

### Community 65 - "Align Output Contract Tests"
Cohesion: 0.21
Nodes (11): rad_without_decoy_boundary_emits_every_reference(), requant_outputs_honor_the_contract(), write_decoy_rad(), lib_counts(), rad_requant_auto_detects_and_filters(), rad_requant_baked_library_format_is_authoritative(), rad_requant_lib_format_counts_report_incompatible_fragments(), requant_reports_run_spanning_time_and_inference_path() (+3 more)

### Community 66 - "GC Fragment Model"
Cohesion: 0.16
Nodes (9): bin_frac(), bin_frac_coarse(), bin_frac_identity_at_101(), build_expected_gc(), gc_ratio(), idx_lut_matches_bin_frac(), normalize_makes_rows_sum_to_one(), ratio_clamped() (+1 more)

### Community 67 - "Spline Band Matrix"
Cohesion: 0.22
Nodes (5): BandMatrix, extrapolates_without_panic(), interpolates_knots_exactly(), monotone_line_is_reproduced(), Spline

### Community 68 - "Determinism & Decoy Releases"
Cohesion: 0.15
Nodes (15): Decoys excluded from quant.genes.sf (#1276), Alignment mode scores by BAM AS tag by default, Pre-release behavior audit (#1140, 111 findings), Decoy references excluded from quant.sf on RAD paths, Warnings for accepted-but-inert flag x mode combinations, meta_info.json inference_path field, Intermediate lz4-compressed RAD, lib_format_counts.json C++ semantics restored (+7 more)

### Community 69 - "EM Step Benchmark"
Cohesion: 0.21
Nodes (6): bench_em(), bounded(), pick(), Rng, single_step(), synthetic()

### Community 70 - "Alignment Stdin Tests"
Cohesion: 0.15
Nodes (6): dev_stdin_from_a_pipe_quantifies_like_the_file(), quant(), SALMON, sam_to_bam(), write_genome_sam(), write_txome_sam()

### Community 71 - "RAD Profile Schema"
Cohesion: 0.17
Nodes (5): RadProfile, SelectiveAlignment, Sketch, build_prelude(), int_desc()

### Community 72 - "AlignQuant FLD Policy"
Cohesion: 0.21
Nodes (13): align_bootstrap_output_aligns_with_quant_and_conserves_mass(), align_gibbs_output_aligns_with_quant_and_conserves_mass(), check(), quantify_two_phase(), read_quant(), read_replicates(), run(), write_sam() (+5 more)

### Community 73 - "Index Build & Load"
Cohesion: 0.24
Nodes (15): a_repeated_transcript_name_is_refused(), a_truncated_reference_store_is_refused(), build(), default_minimizer_len(), effective_threads(), failed_rebuild_keeps_the_existing_index(), IndexBuildOptions, rejects_future_index() (+7 more)

### Community 74 - "Release 2.9 VBEM & cf1"
Cohesion: 0.15
Nodes (16): VBEM (default optimizer), cf1-rs (cDBG construction, PHastR MPHF), salmon-infer components module (component EM engine), kallistors (Denis Demin), --noEffectiveLengthCorrection, --noLengthCorrection fix (#1243), salmon 2.9.0 release notes, Reproducible index builds (+8 more)

### Community 75 - "Docker Build Cache"
Cohesion: 0.12
Nodes (17): with, cache-from, cache-to, context, if-no-files-found, merge-multiple, name, outputs (+9 more)

### Community 76 - "Metadata Truth Tests"
Cohesion: 0.26
Nodes (13): a_fastq_parse_error_names_the_file(), bias_model_shapes_are_recorded(), decoy_targets_are_counted(), fixture(), frag_length_source_names_the_users_input(), library_types_reports_the_applied_format_and_detected_the_observed_one(), meta(), SALMON (+5 more)

### Community 77 - "EqClass Builder (2)"
Cohesion: 0.19
Nodes (3): EquivalenceClassBuilder, TGValue, weight_to_fixed()

### Community 78 - "Packed EqClass Plans (2)"
Cohesion: 0.21
Nodes (13): assert_same_result(), DIGAMMA_MIN, em_result_is_independent_of_thread_count(), empty_plan_has_one_empty_shard(), fill_exp_theta_seq(), fixture(), MIN_EQ_CLASS_WEIGHT, packed() (+5 more)

### Community 79 - "BAM Output Lifecycle"
Cohesion: 0.22
Nodes (3): BamOutput, BamScratch, WriterFailure

### Community 80 - "Release 2.0 Mapping Features"
Cohesion: 0.15
Nodes (16): sshash-rs, LibraryTypeDetector / compatibleHit port, Library-type auto-detection applied as strand filter (-l A), uni-MEM seeding default (--sparseSeeds), num_mapped counted after strand filter (#1025/#1026), homebrew-core distribution, Always-on library-type mismatch detection, Richer meta_info.json fields and diagnostics (+8 more)

### Community 81 - "Baked FLD & Fixed-Point"
Cohesion: 0.22
Nodes (16): Baked FLD in RAD header (1-vs-2-pass dispatch), Deterministic and faster bias models (FFT, obs-exp table), --deterministic quantification, Fixed-point u128 accumulation, libradicl 0.14.0, piscem map-bulk RAD compatibility, RAD chunk compression (--radCompress lz4/zstd), Decoupled mapping/quantification via RAD (--writeRad, --rad) (+8 more)

### Community 82 - "Mapping Modes & Citations Docs (2)"
Cohesion: 0.22
Nodes (9): RAD I/O & deterministic quantification guide, Intermediate RAD disk use and free-space checks, --online one-pass path (removed in 2.12.0), piscem map-bulk, RAD (Reduced Alignment Data) file, --rad quantification mode, --radScratchDir intermediate placement, --writeRad / --skipQuant (+1 more)

### Community 83 - "Genome Projection Pipes"
Cohesion: 0.32
Nodes (13): a_fifo_alignment_input_reads_like_its_file(), a_fifo_genome_alignment_input_projects_like_its_file(), alignment_format_is_read_from_content_not_name(), bounded(), deterministic(), fifo_of(), missing_or_empty_alignment_input_is_reported(), opts() (+5 more)

### Community 84 - "Library Type Compatibility (2)"
Cohesion: 0.13
Nodes (12): format_id_matches_canonical_name(), ReadOrientation, Away, None, Same, Toward, ReadStrandedness, A (+4 more)

### Community 85 - "EM Determinism Audit"
Cohesion: 0.14
Nodes (10): V1: quant.sf not reproducible across thread counts (fixed #1167), EmAccel (None / SQUAREM / DAAREM), finalize_truncate_redistribute (mass-preserving min_alpha truncation), max_rel_diff convergence check (fixed 2048-element chunks), M-step kernels (em_step_seq/par, vbem_step_seq/par), optimize_packed_with_init (production EM entry), PackedEqClasses CSR structure, Complete (txps, bins) range-factorized class ordering (2477cd28) (+2 more)

### Community 86 - "Gibbs & Bootstrap Replicates"
Cohesion: 0.20
Nodes (12): Inferential replicates guide, edgeR catchSalmon, Simulated FDR evaluation (swish, edgeR catchSalmon), fishpond / swish, Nearly flat likelihood across near-identical isoforms, Inferential replicates, --legacyBootstraps (hidden flag), VBEM sparse prior zeroing isoforms (+4 more)

### Community 87 - "RAD Reader Concurrency (2)"
Cohesion: 0.15
Nodes (4): GatedChunks, GatedChunks<'_, R>, QueueGate, read_rad_parallel()

### Community 88 - "Reference Sequences Build (2)"
Cohesion: 0.16
Nodes (3): &'a RefSeqs, RefSeqs, load_ref_seqs()

### Community 89 - "Reference Digests"
Cohesion: 0.16
Nodes (6): compute_ref_hashes(), is_printable(), processed_name(), read_decoy_names(), RefDigests, RefHashes

### Community 90 - "Collapsed EqClasses Infer (2)"
Cohesion: 0.29
Nodes (11): BootstrapMode, Legacy, RandomStartEm, EffLens, EmOptions, finalize_truncate_redistribute(), InitAlphas, optimize_packed_with_init() (+3 more)

### Community 91 - "2.12 API Changes"
Cohesion: 0.14
Nodes (12): salmon_quant::map_to_rad (replaces quantify), MapResult (replaces QuantResult), salmon_align::quantify_rad, alevin-fry redirect, Introduction, Love et al. 2016 alpine sequence bias, Patro et al. 2017 Salmon (Nature Methods), Decoupled mapping/quantification via RAD (+4 more)

### Community 92 - "Release Workflow Jobs"
Cohesion: 0.14
Nodes (13): env, if, needs, runs-on, jobs, announce, name, on (+5 more)

### Community 93 - "Decoy & Readiness Tests"
Cohesion: 0.24
Nodes (10): decoys_are_excluded_and_mass_is_conserved_with_bias(), SALMON, an_unsatisfiable_thread_count_is_clamped_not_a_panic(), fixture(), gene_map_is_honoured_in_the_default_reads_mode(), meta_number(), rc(), recovered_pairs_keep_their_strand_when_mate_two_anchors() (+2 more)

### Community 94 - "EM Stop Rule"
Cohesion: 0.18
Nodes (11): EmAccel, EmAccelArg, Daarem, None, Squarem, EmStop, EmStopArg, Global (+3 more)

### Community 95 - "SAM Record Writing"
Cohesion: 0.24
Nodes (7): record_base(), AlignmentRecord, complement(), SamWriter, write_fragment(), write_record(), write_sequence()

### Community 96 - "Migration & 2.0 Index"
Cohesion: 0.19
Nodes (12): Sketch-mode degenerate AS:i:0, Gene-map unmatched transcript reporting (genemap_unmatched_txps.json), --ignoreTxVersion gene-map version-suffix matching, Migration Guide: C++ salmon 1.12 to salmon 2.0 (Rust), New 2.0 index format (C++ pufferfish indices rejected), quant.sf / bootstrap output compatibility with tximport/fishpond, Removed quant/index options (error + migration note), --sketch alignment-free pseudoalignment mode (+4 more)

### Community 97 - "Salmon 2.0 Rust Rewrite"
Cohesion: 0.19
Nodes (12): salmon alevin removed (alevin-fry), BSD-3-Clause license, C++ salmon 1.12 (final C++ release), MIGRATION.md, Patro et al. 2017, Nature Methods (Salmon), salmon 2.0.0 (Rust rewrite), Selective-alignment quantification, Abundance-aware fragment-length distribution training (+4 more)

### Community 98 - "CLI Main & Args (2)"
Cohesion: 0.18
Nodes (10): apply_stop_rule(), bam_compress_threads_is_optional_and_positive(), bam_compress_threads_requires_write_bam(), bootstrap_mode(), fld_policy_defaults_to_baked_and_parses_all_variants(), fld_prior_args_distinguish_supplied_from_default(), global_rule_uses_salmon_1x_minimum_iterations(), legacy_bootstraps_flag_selects_the_old_replicates() (+2 more)

### Community 99 - "Compressed Input Open"
Cohesion: 0.29
Nodes (9): maybe_decompress_stream(), open_maybe_compressed(), open_maybe_compressed_raw(), reads_all_members_of_a_concatenated_gzip(), sniffs_every_supported_format_by_content(), strip_compression_extension(), strips_one_compression_extension(), very_short_file_opens_as_plain_text() (+1 more)

### Community 100 - "Bias Dumps"
Cohesion: 0.32
Nodes (10): BiasDump, dump_bias_models_to_file(), gz_segments_round_trip(), gz_write(), write_aux_bias_dumps(), write_f64_gz(), write_fld_dump(), write_gz_segments() (+2 more)

### Community 101 - "Release Testing Checklist"
Cohesion: 0.18
Nodes (10): Parity vs C++ salmon, Determinism / robustness checks, salmon pre-release testing checklist, MIN_READABLE_INDEX_VERSION index gate, #1019 abundance-phase hang regression guard, Mass conservation check (sum NumReads ~ num_mapped), Release test datasets (sample_data, simulated human, GRCh38 decoy, SRR1039508), sshash v6 port (index format v2) (+2 more)

### Community 102 - "Docker Workflow Env"
Cohesion: 0.17
Nodes (11): env, DOCKERHUB_IMAGE, GHCR_IMAGE, name, on, push, workflow_dispatch, permissions (+3 more)

### Community 103 - "Mapping Modes & Citations Docs (3)"
Cohesion: 0.17
Nodes (11): Model & bias options (--seqBias, --gcBias, --posBias), aux_info/ambig_info.tsv, cmd_info.json, Output format specification, aux_info/eq_classes.txt.gz, aux_info/fld.gz, libParams/flenDist.txt, Index directory (piscem-rs format) (+3 more)

### Community 104 - "Align Benchmarks"
Cohesion: 0.29
Nodes (8): anchored(), batch(), CASES, full_window(), gen(), MISMATCHES, READ_LENS, run()

### Community 105 - "Global Allocator"
Cohesion: 0.18
Nodes (8): GLOBAL, NAME, Global allocator choice (mimalloc), scripts/bench_em_rad.sh benchmark harness, examples/cache_bench.rs single-thread allocation microbenchmark, Cross-thread free pattern in mapping phase, mimalloc global allocator, sysalloc diagnostic build feature

### Community 106 - "Dump EqClass Tests"
Cohesion: 0.31
Nodes (9): Dump, dump_eq_describes_the_classes_in_both_modes(), dump_eq_with_skip_quant_says_it_has_nothing_to_dump(), fixture(), num_mapped(), quant(), read_dump(), SALMON (+1 more)

### Community 107 - "DP Aligner Kernel"
Cohesion: 0.24
Nodes (11): dna5_encode_rev_into(), flank_shortcut(), flank_shortcut_matches_the_dp(), flank_shortcut_stays_out_of_soft_clipping(), flank_softclip_recovers_read_overhang(), FlankShortcut, Below, Score (+3 more)

### Community 108 - "Bias Config Merge (2)"
Cohesion: 0.33
Nodes (3): enriched_5p_has_higher_early_weight(), SimplePosBias, uniform_observations_project_near_flat()

### Community 109 - "CLI Flag Status"
Cohesion: 0.20
Nodes (9): Warn-and-ignore mode x flag applicability table (#1156), --bamCompressThreads BGZF compression pool, CLI flag status table: C++ salmon vs Rust port, Flag status legend (Implemented/Default/Accept+warn/Rejected/Covered), --fldPolicy (baked/derive/prior) RAD requant FLD source, --noLengthCorrection / --noEffectiveLengthCorrection semantics (#1243), One-pass --online path removed in 2.12.0, --reduceGCMemory rank-bitvector GC representation (now default) (+1 more)

### Community 110 - "Audit Open Issues"
Cohesion: 0.20
Nodes (11): #1180 right-orphan AS:i:0 fix, Tracking issue #1179 (unresolved audit items), U1: EM parallel scaling ceiling past ~16 threads (deferred), U2 #1141: mapping-record self-consistency (CIGAR vs AS), U3 #1142: mapping-output completeness (QUAL, M5, MC/MQ, @RG), Deterministic EM thread-scaling raw TSV results, SRR21186103 / GENCODE v49 benchmark workload (26.1M fragments), EM/VBEM phase deterministic thread-scaling optimization record (+3 more)

### Community 111 - "Release Workflow Steps"
Cohesion: 0.40
Nodes (11): steps, steps, steps, steps, steps, id, if, name (+3 more)

### Community 112 - "Mapping Modes & Citations Docs (4)"
Cohesion: 0.24
Nodes (7): --decoder gzip decoder selection, Command-line interface reference, Mapping options (--sketch, --allowDovetail, --minScoreFraction), salmon alevin (removed; redirects to alevin-fry), salmon debug-map, salmon quant, salmon quantmerge

### Community 113 - "Mapping Modes & Citations Docs (5)"
Cohesion: 0.20
Nodes (8): tximport average transcript length offset, aux_info/bootstrap output (names.tsv.gz, bootstraps.gz), salmon issue #1229 (whole-blood technical replicates), tximport, --geneMap transcript-to-gene map, aux_info/genemap_unmatched_txps.json, aux_info/bootstrap/ (names.tsv.gz, bootstraps.gz), quant.sf

### Community 114 - "Skip-Quant Tests"
Cohesion: 0.31
Nodes (8): fixture(), SALMON, sequence(), skip_quant_does_not_count_wrong_strand_fragments_as_mapped(), skip_quant_in_default_reads_mode_says_the_gene_map_is_skipped(), skip_quant_under_deterministic_alignment_keeps_the_rad_in_the_output_dir(), skip_quant_under_deterministic_maps_and_stops(), write_sam()

### Community 115 - "Inferential Replicates Output"
Cohesion: 0.27
Nodes (10): Inferential replicates (bootstrap and Gibbs), quant.sf output, tximport / tximeta / fishpond / swish, Alignment mode (-a) bootstrap/Gibbs support, Bootstrap output aligned row-for-row with quant.sf, Min-alpha redistribution and inference_truncated_mass, salmon 2.1.1, DAAREM (Henderson & Varadhan 2019) (+2 more)

### Community 116 - "Release 2.11 Audit"
Cohesion: 0.29
Nodes (9): Audit of 2.10.0 (#1283), -a refuses BAM with non-adjacent mates (#1279), Non-zero exit when no fragments processed/mapped (#1277), Orphan mates judged on own strand (#1282), salmon_core RefView / ref_view, salmon 2.11.0 release notes, Sketch short-mate stale hits fix (#1271), salmon_core::write_via_rename (+1 more)

### Community 117 - "Workspace Crates"
Cohesion: 0.53
Nodes (10): salmon-align, salmon-cli, salmon-core, salmon-eqclass, salmon-index, salmon-infer, salmon-map, salmon-model (+2 more)

### Community 120 - "Alignment Input Readers"
Cohesion: 0.31
Nodes (7): AlignmentRecords, Bam, Sam, BamReader, open_alignment_input(), rad_library_is_paired(), SamReader

### Community 121 - "Flag Matrix Tests"
Cohesion: 0.36
Nodes (6): fixture(), mode_flag_matrix_end_to_end(), os(), salmon(), SALMON, sequence()

### Community 122 - "CLI Output Contract Tests"
Cohesion: 0.42
Nodes (8): an_orphan_record_carries_its_own_strand(), deterministic_flow_honors_the_output_contract_with_decoys(), fixture(), quant_names(), revcomp(), SALMON, sequence(), single_end_frag_length_source_is_the_prior()

### Community 124 - "Component Summary"
Cohesion: 0.25
Nodes (5): ComponentSummary, EmStop, Componentwise, Global, Projected

### Community 125 - "Release 2.4 BAM Threads"
Cohesion: 0.28
Nodes (9): --bamCompressThreads, BGZF dedicated rayon pool, Dependency refresh (#1068/#1070), Docker images include procps (Nextflow), --fldPolicy (baked/derive/prior), frag_length_source in meta_info.json, salmon 2.4.0, --writeBam BAM mapping output (#1065) (+1 more)

### Community 126 - "DAAREM & Allocator Releases"
Cohesion: 0.28
Nodes (8): libradicl, docs/allocator-choice.md, CRAN daarem package daarem_base_objfn, salmon-infer EmOptions::validate, salmon-cli global_alloc module (mimalloc), RAD v2 versioned format (piscem >= 0.24), salmon 2.8.0 release notes, SQUAREM EM accelerator

### Community 127 - "SSHash Index v2"
Cohesion: 0.25
Nodes (8): Cologni & Pibiri, Proposition 24 (centre-closest tie-break), Index format v2, salmon 2.7.0 release notes, info.json sshash_format_version, sshash-rs 0.7, C++ SSHash v6.0.0, docs/sshash-v6-port.md, Index format v3 (sshash-lib 0.8 dictionary 6.0)

### Community 128 - "Release Workflow Matrix"
Cohesion: 0.22
Nodes (9): container, if, name, needs, runs-on, strategy, build-local-artifacts, fail-fast (+1 more)

### Community 129 - "Orphan & Fragment Probs (2)"
Cohesion: 0.25
Nodes (6): PlacementInfo, MateStatus, PairedEndLeft, PairedEndPaired, PairedEndRight, SingleEnd

### Community 130 - "Library Type Compatibility (3)"
Cohesion: 0.29
Nodes (3): FragCfg, LibraryFormat, .MAX_FORMAT_ID

### Community 131 - "Index Build & Load (2)"
Cohesion: 0.43
Nodes (5): build_and_load_roundtrip(), build_composes_short_transcript_short_and_n_decoys(), finish_refseq_store(), header_name(), write_refseq_store()

### Community 132 - "EmAccel Enum"
Cohesion: 0.25
Nodes (6): EmAccel, Daarem, None, Squarem, InitAlphas<'a>, .NONE

### Community 133 - "Bias Convolution Bench"
Cohesion: 0.39
Nodes (3): Lcg, main(), GC_SAMP_STRIDE

### Community 134 - "Atomic Math & FLD CDF (2)"
Cohesion: 0.25
Nodes (5): LengthCorrection, Effective, None, Raw, smoothed_effective_length()

### Community 135 - "Reference Sequences Build (3)"
Cohesion: 0.39
Nodes (6): bam_records_for_a_fragment_are_collated(), quant_to_sam_and_bam(), read_bam_records(), read_sam_records(), record_sort_key(), sam_and_bam_mapping_outputs_are_semantically_equal()

### Community 136 - "Index Provenance Types"
Cohesion: 0.25
Nodes (6): IndexProvenance, MappingType, Alignment, Mapping, Pseudo, WriterProvenance

### Community 137 - "Mapping Parity vs C++"
Cohesion: 0.32
Nodes (6): Chain sub-optimality threshold default difference (~80% of SA gap), Mapping parity vs C++ salmon investigation, PairingConfig::default (crates/salmon-map/src/pair.rs), Poly-A clipping mismatch artifact (simulated-data residual cause 1), Seed-enumeration sensitivity tail (147 reads, cause 2), Benign symmetric tie-break/chaining residual (~20%)

### Community 138 - "Decode Deadlock & Seeding"
Cohesion: 0.25
Nodes (8): mismatchSeedSkip false start (partial symptom), piscem-rs mapper (derives k-mer orientation from query), Pufferfish SSHash streaming k-mer orientation bug (hitFW), Adaptive decode broker (DecodeAllocation::Adaptive vs PinnedPerFile), Deadlock: broker-driven resizing + paired parallel decode, rapidgzip DecoderPool::release lost wakeup (fixed rapidgzip-core 0.3.1), paraseq PairedReader separate R1/R2 locks, tests/pool_resize_flush_safety.rs regression test

### Community 139 - "Dovetail Pairing 2.12"
Cohesion: 0.25
Nodes (6): --dovetailTolerance flag, 5' end fragment rule, Insert-only scoring of dovetailed pairs, salmon_map::PairingConfig.dovetail_tolerance, piscem-rs 0.13.0, Sketch mode one-base fragment length fix

### Community 140 - "Genome Projection Mode"
Cohesion: 0.32
Nodes (5): Coordinate-sorted BAM rejected up front, bramble projection, Genome-alignment quantification, Genome projection mode (-a + --annotation), Name-grouped BAM requirement

### Community 141 - "Docker Workflow Steps"
Cohesion: 0.36
Nodes (8): steps, steps, id, if, name, run, uses, working-directory

### Community 142 - "Bump & Publish Script"
Cohesion: 0.46
Nodes (6): bump_manifest(), die(), print_cmd(), run(), bump_and_publish.sh script, usage()

### Community 143 - "DP Aligner Kernel (2)"
Cohesion: 0.33
Nodes (5): gap_pass(), gap_pass_matches_scalar(), gap_pass_scalar(), gap_scan4(), KswScratch

### Community 144 - "FLD Policy Arg"
Cohesion: 0.33
Nodes (6): FldPolicy, FldPolicyArg, Baked, Derive, Prior, warn_unsupported_fld_policy()

### Community 145 - "SalmonError"
Cohesion: 0.29
Nodes (6): Result, SalmonError, InvalidValue, Io, Other, UnknownLibraryFormat

### Community 146 - "GTF/GFF Attribute Parsing"
Cohesion: 0.33
Nodes (6): annotate_read_error(), ensembl_gff3_resolves_gene_via_parent(), extract_attr(), gtf_and_gff_attr_extraction(), open_maybe_compressed(), read_transcript_gene_map()

### Community 147 - "Posterior Method"
Cohesion: 0.33
Nodes (4): PosteriorMethod, Bootstrap, Gibbs, None

### Community 148 - "EM Thread Scaling Chart"
Cohesion: 0.38
Nodes (7): #1169 baseline EM, Deterministic EM speedup by thread count (chart), Deterministic EM, Local host 2x EPYC 9575F, Newton host 2x Xeon E5-2699 v4, Selected 128-shard EM candidate, Speedup metric T1/Tp (five-run medians)

### Community 149 - "Release Workflow Env"
Cohesion: 0.33
Nodes (7): env, needs, runs-on, env, BUILD_MANIFEST_NAME, GH_TOKEN, build-global-artifacts

### Community 150 - "Release Host Outputs"
Cohesion: 0.29
Nodes (7): env, if, needs, outputs, runs-on, host, val

### Community 151 - "Release Artifact Upload"
Cohesion: 0.29
Nodes (7): with, merge-multiple, name, path, pattern, persist-credentials, submodules

### Community 152 - "EM RAD Bench Script"
Cohesion: 0.52
Nodes (6): affinity_for(), json_scalar(), phase_seconds(), run_one(), bench_em_rad.sh script, usage()

### Community 153 - "Installation Methods"
Cohesion: 0.33
Nodes (7): cargo install salmon-cli, conda (bioconda salmon), Installation, Docker image combinelab/salmon, Build from source, Install script, Nextflow entrypointOverride

### Community 154 - "Library Types & Quick Start"
Cohesion: 0.29
Nodes (7): Decoy-aware indexing (-d decoys.txt), Quick start, --ignoreTxVersion for Ensembl cDNA + GTF, -l A auto-detection, Fragment library types, Library-type string (orientation/strandedness/direction), TopHat library type correspondence

### Community 155 - "Mapping Modes & Citations Docs (6)"
Cohesion: 0.29
Nodes (4): em_components per-component convergence, lib_format_counts.json, -l A library type resolution, salmon-rs meta_info extensions

### Community 157 - "Piscem Bulk Records"
Cohesion: 0.33
Nodes (4): mapping_type_to_status(), PiscemBulkReadRecord, RadFragment, SalmonBulkRecord

### Community 160 - "Audit Reconciliation 1140"
Cohesion: 0.33
Nodes (6): Audit reconciliation #1140 (2.6.0 deterministic-default readiness), Deterministic two-phase quantification mode (RAD-based), Disk story: free-space preflight, partial RAD cleanup, --radScratchDir (#1146), Release-blocking output contract fixes (#1153: decoy leak, cmd_info clobber, meta_info fields), V2: EM 2.14x under-parallelization (fixed #1172), Decoy-aware indexing (-d/--decoys)

### Community 161 - "Release 2.6/2.10 Fixes"
Cohesion: 0.33
Nodes (5): ksw2rs KSW_EZ_NO_MAX, salmon 2.10.0 release notes, Selective-alignment flank scoring bounds and early stop (#1261), VBEM negligible-transcript pruning, SA positions from projected read ends (#1231)

### Community 162 - "SSHash Bucket Cache"
Cohesion: 0.47
Nodes (4): sshash-lib / tiny-dict 0.8.2, sshash bucket-bounds cache notes, Temporary [patch.crates-io] for sshash-lib, PiscemStreamingQuery

### Community 163 - "Docker Jobs"
Cohesion: 0.33
Nodes (6): runs-on, jobs, build, merge, needs, runs-on

### Community 164 - "Docker Platform Matrix"
Cohesion: 0.33
Nodes (6): strategy, platform, runner, include, fail-fast, matrix

### Community 165 - "Docker Secrets Env"
Cohesion: 0.33
Nodes (6): DOCKERHUB, DOCKERHUB_TOKEN, MINOR, platform, VER, env

### Community 166 - "Release Workflow Env (2)"
Cohesion: 0.33
Nodes (6): plan, publishing, tag-flag, env, outputs, runs-on

### Community 168 - "Mapping Modes & Citations Docs (8)"
Cohesion: 0.33
Nodes (5): RAD I/O & determinism options (--writeRad, --deterministic, --radCompress), Alignment provenance (@PG chain), chunk_codec RAD chunk compression, libradicl RAD format, RAD output (--writeRad)

### Community 169 - "Refactor Bookmark"
Cohesion: 0.40
Nodes (4): CMake dependency floor checks (spdlog/fmt, parallel_hashmap, nlohmann/json, Eigen>=3.4), FQFeeder modernization (salmon_fqfeeder::ParserConfig), Quant hot-path optimization pass (read-based mode), Refactor Checkpoint 2026-03-08 (codex/develop-refactor)

### Community 172 - "RAD Input Profile"
Cohesion: 0.50
Nodes (4): detect_input_profile(), RadInputProfile, PiscemBulk, Salmon

### Community 173 - "Release 2.1.2 Piscem Fix"
Cohesion: 0.40
Nodes (5): piscem-rs, align_chain out-of-bounds guard, Mapping panic fix in seed collector (#1038), piscem-rs 0.6.1, salmon 2.1.2

### Community 174 - "Release Announcement Env"
Cohesion: 0.40
Nodes (5): ANNOUNCEMENT_BODY, ANNOUNCEMENT_TITLE, PRERELEASE_FLAG, RELEASE_COMMIT, env

### Community 175 - "Website TSConfig"
Cohesion: 0.40
Nodes (4): astro/tsconfigs/strict, exclude, extends, include

### Community 176 - "AlignQuant FLD Policy (2)"
Cohesion: 0.50
Nodes (4): FldPolicy, Baked, Derive, Prior

### Community 177 - "Bias Config Merge (3)"
Cohesion: 0.50
Nodes (4): FragmentOutcome, Mapped, MappedOrphan, Unmapped

### Community 181 - "Community Guidelines"
Cohesion: 0.67
Nodes (3): GitHub Bug Report Issue Template, Contributor Covenant Code of Conduct v1.4, Contributing Guidelines (PRs to develop, conventional commits, BSD-3 inbound=outbound)

### Community 184 - "CigarKind"
Cohesion: 0.67
Nodes (3): CigarKind, Match, SoftClip

## Ambiguous Edges - Review These
- `rapidgzip DecoderPool::release lost wakeup (fixed rapidgzip-core 0.3.1)` → `Rejected: shared AtomicF64 accumulator / cached alpha*w scratch / parallel VBEM sum`  [AMBIGUOUS]
  notes/pool-lost-wakeup-deadlock.md · relation: semantically_similar_to

## Knowledge Gaps
- **253 isolated node(s):** `Effective`, `None`, `Raw`, `ScoreWriteWorker<'a>`, `TrainWorker<'a>` (+248 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 922 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **18 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `salmon quant` connect `Mapping Modes & Citations Docs (4)` to `Mapping Modes & Citations Docs (7)`, `Mapping Modes & Citations Docs (8)`, `Mapping Modes & Citations Docs (3)`, `Mapping Modes & Citations Docs (5)`, `Gibbs & Bootstrap Replicates`?**
  _High betweenness centrality (0.096) - this node is a cross-community bridge._
- **Are the 45 inferred relationships involving `quantify_rad()` (e.g. with `quantify_two_phase()` and `deterministic()`) actually correct?**
  _`quantify_rad()` has 45 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Effective`, `None`, `Raw` to the rest of the system?**
  _253 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `EM Component Packing` be split into smaller, more focused modules?**
  _Cohesion score 0.06741573033707865 - nodes in this community are weakly interconnected._
- **Why does `Inference options` connect `Gibbs & Bootstrap Replicates` to `Mapping Modes & Citations Docs (4)`, `Mapping Modes & Citations Docs (6)`?**
  _High betweenness centrality (0.047) - this node is a cross-community bridge._
- **Are the 25 inferred relationships involving `write_alignment_rad()` (e.g. with `coordinate_sorted_unusable()` and `load_ref_bytes()`) actually correct?**
  _`write_alignment_rad()` has 25 INFERRED edges - model-reasoned connections that need verification._

### Low-confidence Hints
_AMBIGUOUS edges — the extractor was unsure. Verify before acting on these._

- **What is the exact relationship between `rapidgzip DecoderPool::release lost wakeup (fixed rapidgzip-core 0.3.1)` and `Rejected: shared AtomicF64 accumulator / cached alpha*w scratch / parallel VBEM sum`?**
  _Edge tagged AMBIGUOUS (relation: semantically_similar_to) - confidence is low._