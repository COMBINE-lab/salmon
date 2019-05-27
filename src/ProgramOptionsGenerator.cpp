#include <vector>
#include <iostream>

#include "ProgramOptionsGenerator.hpp"
#include "SalmonDefaults.hpp"

namespace salmon {

  po::options_description ProgramOptionsGenerator::getBasicOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;
    po::options_description basic("\n"
                                  "basic options");
    basic.add_options()("version,v", "print version string")
      ("help,h", "produce help message")
      ("output,o", po::value<std::string>()->required(), "Output quantification directory.")
      ("seqBias", po::bool_switch(&(sopt.biasCorrect))->default_value(salmon::defaults::seqBiasCorrect),
       "Perform sequence-specific bias correction.")
      ("gcBias", po::bool_switch(&(sopt.gcBiasCorrect))->default_value(salmon::defaults::gcBiasCorrect),
       "[beta for single-end reads] Perform fragment GC bias correction")
      ("threads,p",
       po::value<uint32_t>(&(sopt.numThreads))->default_value(sopt.numThreads),
       "The number of threads to use concurrently.")
      ("incompatPrior",
       po::value<double>(&(sopt.incompatPrior))->default_value(salmon::defaults::incompatPrior),
       "This option "
       "sets the prior probability that an alignment that disagrees with the "
       "specified "
       "library type (--libType) results from the true fragment origin.  "
       "Setting this to 0 "
       "specifies that alignments that disagree with the library type should be "
       "\"impossible\", "
       "while setting it to 1 says that alignments that disagree with the "
       "library type are no "
       "less likely than those that do")
      ("geneMap,g", po::value<string>(),
       "File containing a mapping of transcripts to genes.  If this file is "
       "provided "
       "salmon will output both quant.sf and quant.genes.sf files, where the "
       "latter "
       "contains aggregated gene-level abundance estimates.  The transcript to "
       "gene mapping "
       "should be provided as either a GTF file, or a in a simple tab-delimited "
       "format "
       "where each line contains the name of a transcript and the gene to which "
       "it belongs "
       "separated by a tab.  The extension of the file is used to determine how "
       "the file "
       "should be parsed.  Files ending in \'.gtf\', \'.gff\' or \'.gff3\' are "
       "assumed to "
       "be in GTF "
       "format; files with any other extension are assumed to be in the simple "
       "format. In GTF / GFF format, the \"transcript_id\" is assumed to "
       "contain the "
       "transcript identifier and the \"gene_id\" is assumed to contain the "
       "corresponding "
       "gene identifier.")
      ("meta", po::bool_switch(&(sopt.meta))->default_value(salmon::defaults::metaMode),
       "If you're using Salmon on a metagenomic dataset, consider setting this "
       "flag to disable parts of the "
       "abundance estimation model that make less sense for metagenomic data.");

    // use rich eq classes by default
    sopt.noRichEqClasses = salmon::defaults::noRichEqClasses;
    // mapping cache has been deprecated
    sopt.disableMappingCache = salmon::defaults::disableMappingCache;

    return basic;
  }

  po::options_description ProgramOptionsGenerator::getMappingSpecificOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    po::options_description mapspec("\n"
                                    "options specific to mapping mode");
    mapspec.add_options()
      ("discardOrphansQuasi",
       po::bool_switch(&(sopt.discardOrphansQuasi))->default_value(salmon::defaults::discardOrphansQuasi),
       "[Quasi-mapping mode only] : Discard orphan mappings in quasi-mapping "
       "mode.  If this flag is passed "
       "then only paired mappings will be considered toward quantification "
       "estimates.  The default behavior is "
       "to consider orphan mappings if no valid paired mappings exist.  This "
       "flag is independent of the option to "
       "write the orphaned mappings to file (--writeOrphanLinks).")
      ("validateMappings",
       po::bool_switch(&(sopt.validateMappings))->default_value(salmon::defaults::validateMappings),
       "[Quasi-mapping mode only] : Validate mappings using alignment-based verifcation. "
       "If this flag is passed, quasi-mappings will be validated to ensure that they could give "
       "rise to a reasonable alignment before they are further used for quantification.")
      ("consensusSlack",
       po::value<float>(&(sopt.consensusSlack))->default_value(salmon::defaults::consensusSlack),
       "[Quasi-mapping mode only] : The amount of slack allowed in the quasi-mapping consensus "
       "mechanism.  Normally, a transcript must cover all hits to be considered for mapping.  "
       "If this is set to a fraction, X, greater than 0 (and in [0,1)), then a transcript can fail to cover up to "
       "(100 * X)% of the hits before it is discounted as a mapping candidate.  The default value of this option "
       "is 0.2 if --validateMappings is given and 0 otherwise."
       )
      ("minScoreFraction",
       po::value<double>(&sopt.minScoreFraction),
       "[Quasi-mapping mode (w / mapping validation) only] : The fraction of the optimal possible alignment score that a "
       "mapping must achieve in order to be considered \"valid\" --- should be in (0,1].\n"
       "Salmon Default 0.65 and Alevin Default 0.87"
       )
      ("maxMMPExtension",
       po::value<int32_t>(&sopt.maxMMPExtension)->default_value(salmon::defaults::maxMMPExtension),
       "[Quasi-mapping mode (w / mapping validation) only] : Sets the maximum allowable MMP extension when collecting "
       "suffix array intervals to be used in chaining.  This prevents MMPs from becoming too long, and potentially masking "
       "intervals that would uncover other good quasi-mappings for the read.  This heuristic mimics the idea of the "
       "maximum mappable safe prefix (MMSP) in selective alignment.  Setting a smaller value will potentially allow "
       "for more sensitive, but slower, mapping."
       )
      ("ma",
       po::value<int16_t>(&sopt.matchScore)->default_value(salmon::defaults::matchScore),
       "[Quasi-mapping mode (w / mapping validation) only] : The value given to a match between read and reference nucleotides "
       "in an alignment."
       )
      ("mp",
       po::value<int16_t>(&sopt.mismatchPenalty)->default_value(salmon::defaults::mismatchPenalty),
       "[Quasi-mapping mode (w / mapping validation) only] : The value given to a mis-match between read and reference nucleotides "
       "in an alignment."
       )
      ("go",
       po::value<int16_t>(&sopt.gapOpenPenalty)->default_value(salmon::defaults::gapOpenPenalty),
       "[Quasi-mapping mode (w / mapping validation) only] : The value given to a gap opening in an alignment."
       )
      ("ge",
       po::value<int16_t>(&sopt.gapExtendPenalty)->default_value(salmon::defaults::gapExtendPenalty),
       "[Quasi-mapping mode (w / mapping validation) only] : The value given to a gap extension in an alignment."
       )
      ("bandwidth",
       po::value<int32_t>(&sopt.dpBandwidth)->default_value(salmon::defaults::dpBandwidth),
       "[Quasi-mapping mode (w / mapping validation) only] : The value used for the bandwidth passed to ksw2.  A smaller "
       "bandwidth can make the alignment verification run more quickly, but could possibly miss valid alignments."
       )
      ("allowDovetail",
       po::bool_switch(&(sopt.allowDovetail))->default_value(salmon::defaults::allowDovetail),
       "[Quasi-mapping mode only] : allow dovetailing mappings."
       )
      ("recoverOrphans",
       po::bool_switch(&(sopt.recoverOrphans))->default_value(salmon::defaults::recoverOrphans),
       "[Quasi-mapping mode only] : Attempt to recover the mates of orphaned reads. This uses edlib for "
       "orphan recovery, and so introduces some computational overhead, but it can improve sensitivity."
       )
      ("mimicBT2", // horrible flag name, think of something better
       po::bool_switch(&(sopt.mimicBT2))->default_value(salmon::defaults::mimicBT2),
       "[Quasi-mapping mode (w / mapping validation) only] : Set flags to mimic parameters similar to "
       "Bowtie2 with --no-discordant and --no-mixed flags.  This increases disallows dovetailing reads, and "
       "discards orphans. Note, this does not impose the very strict parameters assumed by RSEM+Bowtie2, "
       "like gapless alignments.  For that behavior, use the --mimiStrictBT2 flag below."
       )
      ("mimicStrictBT2", // horrible flag name, think of something better
       po::bool_switch(&(sopt.mimicStrictBT2))->default_value(salmon::defaults::mimicStrictBT2),
       "[Quasi-mapping mode (w / mapping validation) only] : Set flags to mimic the very strict parameters used by "
       "RSEM+Bowtie2.  This increases --minScoreFraction to 0.8, disallows dovetailing reads, "
       "discards orphans, and disallows gaps in alignments."
       )
      ("hardFilter",
       po::bool_switch(&(sopt.hardFilter))->default_value(salmon::defaults::hardFilter),
       "[Quasi-mapping mode (w / mapping validation) only] : Instead of weighting mappings by their alignment score, "
       "this flag will discard any mappings with sub-optimal alignment score.  The default option of soft-filtering "
       "(i.e. weighting mappings by their alignment score) usually yields slightly more accurate abundance estimates "
       "but this flag may be desirable if you want more accurate 'naive' equivalence classes, rather "
       "than range factorized equivalence classes."
       )
      ("allowOrphansFMD",
       po::bool_switch(&(sopt.allowOrphans))->default_value(salmon::defaults::allowOrphansFMD),
       "[FMD-mapping mode only] : Consider orphaned reads as valid hits when "
       "performing lightweight-alignment.  This option will increase "
       "sensitivity (allow more reads to map and "
       "more transcripts to be detected), but may decrease specificity as "
       "orphaned alignments are more likely "
       "to be spurious.")
      ("writeMappings,z", po::value<string>(&sopt.qmFileName)
       ->default_value(salmon::defaults::quasiMappingDefaultFile)
       ->implicit_value(salmon::defaults::quasiMappingImplicitFile),
       "If this option is provided, then the quasi-mapping "
       "results will be written out in SAM-compatible "
       "format.  By default, output will be directed to "
       "stdout, but an alternative file name can be "
       "provided instead.")
      ("consistentHits,c",
       po::bool_switch(&(sopt.consistentHits))->default_value(salmon::defaults::consistentHits),
       "Force hits gathered during "
       "quasi-mapping to be \"consistent\" (i.e. co-linear and "
       "approximately the right distance apart).")
      ("fasterMapping",
       po::bool_switch(&(sopt.fasterMapping))->default_value(salmon::defaults::fasterMapping),
       "[Developer]: Disables some extra checks during quasi-mapping. This "
       "may make mapping a "
       "little bit faster at the potential cost of returning too many "
       "mappings (i.e. some sub-optimal mappings) "
       "for certain reads.")
      ("quasiCoverage,x",
       po::value<double>(&(sopt.quasiCoverage))->default_value(salmon::defaults::quasiCoverage),
       "[Experimental]: The fraction of the read that must be covered by "
       "MMPs (of length >= 31) if "
       "this read is to be considered as \"mapped\".  This may help to "
       "avoid \"spurious\" mappings. "
       "A value of 0 (the default) denotes no coverage threshold (a single "
       "31-mer can yield a mapping).  "
       "Since coverage by exact matching, large, MMPs is a rather strict "
       "condition, this value should likely "
       "be set to something low, if used.");
      return mapspec;
  }

  po::options_description ProgramOptionsGenerator::getMappingInputOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    po::options_description mapin("\n"
                                    "mapping input options");
    mapin.add_options()
      ("libType,l", po::value<std::string>()->required(), "Format string describing the library type")
      ("index,i", po::value<string>()->required(), "salmon index")
      ("unmatedReads,r",
       po::value<vector<string>>(&(sopt.unmatedReadFiles))->multitoken(),
       "List of files containing unmated reads of (e.g. single-end reads)")
      ("mates1,1", po::value<vector<string>>(&(sopt.mate1ReadFiles))->multitoken(),
       "File containing the #1 mates")
      ("mates2,2", po::value<vector<string>>(&(sopt.mate2ReadFiles))->multitoken(),
       "File containing the #2 mates");
    return mapin;
  }



  po::options_description ProgramOptionsGenerator::getAlignmentInputOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    po::options_description alignin("\n"
                                      "alignment input options");
    alignin.add_options()
    ("discardOrphans",
       po::bool_switch(&(sopt.discardOrphansAln))->default_value(salmon::defaults::discardOrphansAln),
       "[Alignment-based mode only] : Discard orphan alignments in the input "
       ".  If this flag is passed, then only paired alignments will be "
       "considered toward quantification estimates.  The default behavior is "
       "to consider orphan alignments if no valid paired mappings exist.")
    ("libType,l", po::value<std::string>()->required(), "Format string describing the library type")
    ("alignments,a", po::value<vector<string>>()->multitoken()->required(),
     "input alignment (BAM) file(s).")
    ("targets,t", po::value<std::string>()->required(),
     "FASTA format file containing target transcripts.");

    return alignin;
  }

  po::options_description ProgramOptionsGenerator::getAlevinDevsOptions() {
    po::options_description alevindevs("\n"
                                       "alevin-developer Options");
    alevindevs.add_options()
      (
       "indrop", po::bool_switch()->default_value(alevin::defaults::isInDrop),
       "Use inDrop (not extensively tested) Single Cell protocol for the library. must specify w1 too.")
      (
       "w1", po::value<std::string>(),
       "Must be used in conjunction with inDrop;")
      (
       "dumpBarcodeEq", po::bool_switch()->default_value(alevin::defaults::dumpBarcodeEq),
       "Dump JointEqClas with umi-barcode count.")
      (
       "iupac,u",po::value<std::string>(),
       "<Deprecated>iupac code for cell-level barcodes.")
      (
       "end",po::value<uint32_t>(),
       "Cell-Barcodes end (5 or 3) location in the read sequence from where barcode has to"
       "be extracted. (end, umiLength, barcodeLength)"
       " should all be provided if using this option")
      (
       "umiLength",po::value<uint32_t>(),
       "umi length Parameter for unknown protocol. (end, umiLength, barcodeLength)"
       " should all be provided if using this option")
      (
       "barcodeLength",po::value<uint32_t>(),
       "umi length Parameter for unknown protocol. (end, umiLength, barcodeLength)"
       " should all be provided if using this option")
      (
       "noem",po::bool_switch()->default_value(alevin::defaults::noEM),
       "do not run em")
      (
       "trimRight",po::value<uint32_t>()->default_value(alevin::defaults::trimRight),
       "The number of bases to trim off the 5' (right) end of the read seequence.")
      (
       "naiveEqclass", po::bool_switch()->default_value(alevin::defaults::naiveEqclass),
       "Run naive per equivalence class deduplication, generating only total number of UMIs")
      (
       "noDedup", po::bool_switch()->default_value(alevin::defaults::noDedup),
       "Stops the pipeline after CB sequence correction and quasi-mapping reads.")
      (
       "freqThreshold", po::value<uint32_t>()->default_value(alevin::defaults::freqThreshold),
       "threshold for the frequency of the barcodes");
    return alevindevs;
  }

  po::options_description ProgramOptionsGenerator::getAlevinBasicOptions(SalmonOpts& sopt) {
    po::options_description alevinspec("\n"
                                       "alevin-specific Options");
    alevinspec.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ("output,o", po::value<std::string>()->required(), "Output quantification directory.")
      ("threads,p",
       po::value<uint32_t>(&(sopt.numThreads))->default_value(sopt.numThreads),
       "The number of threads to use concurrently.")
      (
       "tgMap", po::value<std::string>(), "transcript to gene map tsv file")
      (
       "hash", po::value<std::string>(), "Secondary input point for Alevin "
       "using Big freaking Hash (bfh.txt) file. Works Only with --chromium")
      (
       "dropseq", po::bool_switch()->default_value(alevin::defaults::isDropseq),
       "Use DropSeq Single Cell protocol for the library")
      (
       "chromiumV3", po::bool_switch()->default_value(alevin::defaults::isChromiumV3),
       "Use 10x chromium v3 Single Cell protocol for the library.")
      (
       "chromium", po::bool_switch()->default_value(alevin::defaults::isChromium),
       "Use 10x chromium v2 Single Cell protocol for the library.")
      (
       "gemcode", po::bool_switch()->default_value(alevin::defaults::isGemcode),
       "Use 10x gemcode v1 Single Cell protocol for the library.")
      (
       "celseq", po::bool_switch()->default_value(alevin::defaults::isCELSeq),
       "Use CEL-Seq Single Cell protocol for the library.")
      (
       "celseq2", po::bool_switch()->default_value(alevin::defaults::isCELSeq2),
       "Use CEL-Seq2 Single Cell protocol for the library.")
      (
       "whitelist", po::value<std::string>(),
       "File containing white-list barcodes")
      (
       "noQuant", po::bool_switch()->default_value(alevin::defaults::noQuant),
       "Don't run downstream barcode-salmon model.")
      (
       "numCellBootstraps",po::value<uint32_t>()->default_value(alevin::defaults::numBootstraps),
       "Generate mean and variance for cell x gene matrix quantification"
       " estimates.")
      (
       "forceCells",po::value<uint32_t>()->default_value(alevin::defaults::forceCells),
       "Explicitly specify the number of cells.")
      (
       "expectCells",po::value<uint32_t>()->default_value(alevin::defaults::expectCells),
       "define a close upper bound on expected number of cells")
      (
       "mrna", po::value<std::string>(),
       "path to a file containing mito-RNA gene, one per line")
      (
       "rrna", po::value<std::string>(),
       "path to a file containing ribosomal RNA, one per line")
      (
       "keepCBFraction", po::value<double>(),
       "fraction of CB to keep, value must be in range (0,1], use 1 to quantify all CB."
       )
      (
       "dumpfq", po::bool_switch()->default_value(alevin::defaults::dumpFQ),
       "Dump barcode modified fastq file for downstream analysis by"
       " using coin toss for multi-mapping.")
      (
       "dumpBfh", po::bool_switch()->default_value(alevin::defaults::dumpBFH),
       "dump the big hash with all the barcodes and the UMI sequence.")
      (
       "dumpUmiGraph", po::bool_switch()->default_value(alevin::defaults::dumpUmiGraph),
       "dump the per cell level Umi Graph.")
      (
       "dumpFeatures", po::bool_switch()->default_value(alevin::defaults::dumpFeatures),
       "Dump features for whitelist and downstream analysis.")
      (
       "dumpMtx", po::bool_switch()->default_value(alevin::defaults::dumpMtx),
       "Dump cell v transcripts count matrix in sparse mtx format.")
      (
       "lowRegionMinNumBarcodes", po::value<uint32_t>()->default_value(alevin::defaults::lowRegionMinNumBarcodes),
       "Minimum Number of CB to use for learning Low confidence region (Default: 200).")
      (
       "maxNumBarcodes", po::value<uint32_t>()->default_value(alevin::defaults::maxNumBarcodes),
       "Maximum allowable limit to process the cell barcodes. (Default: 100000)");
    return alevinspec;
  }

  po::options_description ProgramOptionsGenerator::getAlignmentSpecificOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    sopt.useMassBanking = salmon::defaults::useMassBanking;
    sopt.noRichEqClasses = salmon::defaults::noRichEqClasses;

    auto onErrorModel = [&sopt](bool noErrorModel) -> void {
      sopt.useErrorModel = !noErrorModel;
    };

    po::options_description alnspec("\n"
                                      "alignment-specific options");
    alnspec.add_options()
      ("noErrorModel",
       po::bool_switch()->default_value(salmon::defaults::noErrorModel)->notifier(onErrorModel),
       "Turn off the alignment error model, which takes into "
       "account the the observed frequency of different types of mismatches / "
       "indels when computing the likelihood of a given alignment. Turning this off can "
       "speed up alignment-based salmon, but can harm quantification accuracy.")
    ("numErrorBins",
     po::value<uint32_t>(&(sopt.numErrorBins))->default_value(salmon::defaults::numErrorBins),
     "The number of bins into which to divide "
     "each read when learning and applying the error model.  For example, a "
     "value of 10 would mean that "
     "effectively, a separate error model is leared and applied to each 10th "
     "of the read, while a value of "
     "3 would mean that a separate error model is applied to the read "
     "beginning (first third), middle (second third) "
     "and end (final third).")
      (
       "sampleOut,s",
       po::bool_switch(&(sopt.sampleOutput))->default_value(salmon::defaults::sampleOutput),
       "Write a \"postSample.bam\" file in the output directory "
       "that will sample the input alignments according to the estimated "
       "transcript abundances. If you're "
       "going to perform downstream analysis of the alignments with tools "
       "which don't, themselves, take "
       "fragment assignment ambiguity into account, you should use this "
       "output.")
      (
       "sampleUnaligned,u",
       po::bool_switch(&(sopt.sampleUnaligned))->default_value(salmon::defaults::sampleUnaligned),
       "In addition to sampling the aligned reads, also write "
       "the un-aligned reads to \"postSample.bam\".")
      ("gencode", po::bool_switch(&(sopt.gencodeRef))->default_value(salmon::defaults::gencodeRef),
       "This flag will expect the input transcript fasta to be "
       "in GENCODE format, and will split the transcript name at the first "
       "\'|\' character.  These reduced names will be used in "
       "the output and when looking for these transcripts in a gene to "
       "transcript GTF.")
      (
       "mappingCacheMemoryLimit",
       po::value<uint32_t>(&(sopt.mappingCacheMemoryLimit))
       ->default_value(salmon::defaults::mappingCacheMemoryLimit),
       "If the file contained fewer than this "
       "many mapped reads, then just keep the data in memory for subsequent "
       "rounds of inference. Obviously, this value should "
       "not be too large if you wish to keep a low memory usage, but setting it "
       "large enough to accommodate all of the mapped "
       "read can substantially speed up inference on \"small\" files that "
       "contain only a few million reads.");

    return alnspec;
  }

  po::options_description ProgramOptionsGenerator::getAdvancedOptions(int32_t& numBiasSamples, SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    po::options_description advanced("\n"
                                     "advanced options");
    advanced.add_options()
      ("alternativeInitMode",
       po::bool_switch(&(sopt.alternativeInitMode))->default_value(salmon::defaults::alternativeInitMode),
       "[Experimental]: Use an alternative strategy (rather than simple "
       "interpolation between) the "
       "online and uniform abundance estimates to initalize the EM / VBEM "
       "algorithm.")
      ("auxDir",
       po::value<std::string>(&(sopt.auxDir))->default_value(salmon::defaults::auxDir),
       "The sub-directory of the quantification directory where auxiliary "
       "information "
       "e.g. bootstraps, bias parameters, etc. will be written.")
      ("skipQuant", po::bool_switch(&(sopt.skipQuant))->default_value(salmon::defaults::skipQuant),
       "Skip performing the actual transcript quantification (including any Gibbs sampling or bootstrapping)."
       )
      ("dumpEq", po::bool_switch(&(sopt.dumpEq))->default_value(salmon::defaults::dumpEq),
       "Dump the equivalence class counts "
       "that were computed during quasi-mapping")
      ("dumpEqWeights,d",
       po::bool_switch(&(sopt.dumpEqWeights))->default_value(salmon::defaults::dumpEqWeights),
       "Includes \"rich\" equivlance class weights in the output when "
       "equivalence "
       "class information is being dumped to file.")
      ("minAssignedFrags",
       po::value<std::uint64_t>(&(sopt.minRequiredFrags))->default_value(salmon::defaults::minAssignedFrags),
       "The minimum number of fragments that must be assigned to the "
       "transcriptome for "
       "quantification to proceed.")
      ("reduceGCMemory",
       po::bool_switch(&(sopt.reduceGCMemory))->default_value(salmon::defaults::reduceGCMemory),
       "If this option is selected, a more memory efficient (but slightly "
       "slower) representation is "
       "used to compute fragment GC content. Enabling this will reduce "
       "memory usage, but can also reduce "
       "speed.  However, the results themselves will remain the same.")
      ("biasSpeedSamp",
       po::value<std::uint32_t>(&(sopt.pdfSampFactor))->default_value(salmon::defaults::biasSpeedSamp),
       "The value at which the fragment length PMF is down-sampled "
       "when evaluating sequence-specific & GC fragment bias.  Larger "
       "values speed up effective "
       "length correction, but may decrease the fidelity of bias modeling "
       "results.")
      ("fldMax",
       po::value<size_t>(&(sopt.fragLenDistMax))->default_value(salmon::defaults::maxFragLength),
       "The maximum fragment length to consider when building the empirical "
       "distribution")
      ("fldMean",
       po::value<size_t>(&(sopt.fragLenDistPriorMean))->default_value(salmon::defaults::fragLenPriorMean),
       "The mean used in the fragment length distribution prior")
      ("fldSD",
       po::value<size_t>(&(sopt.fragLenDistPriorSD))->default_value(salmon::defaults::fragLenPriorSD),
       "The standard deviation used in the fragment length distribution "
       "prior")
      ("forgettingFactor,f",
       po::value<double>(&(sopt.forgettingFactor))->default_value(salmon::defaults::ffactor),
       "The forgetting factor used "
       "in the online learning schedule.  A smaller value results in "
       "quicker learning, but higher variance "
       "and may be unstable.  A larger value results in slower learning but "
       "may be more stable.  Value should "
       "be in the interval (0.5, 1.0].")
      ("initUniform",
       po::bool_switch(&(sopt.initUniform))->default_value(salmon::defaults::initUniform),
       "initialize the offline inference with uniform parameters, rather "
       "than seeding with online parameters.")
      ("maxReadOcc,w",
       po::value<uint32_t>(&(sopt.maxReadOccs))->default_value(salmon::defaults::maxReadOccs),
       "Reads \"mapping\" to more than this many places won't be "
       "considered.")
      ("noLengthCorrection",
       po::bool_switch(&(sopt.noLengthCorrection))->default_value(salmon::defaults::noLengthCorrection),
       "[experimental] : Entirely disables length correction when "
       "estimating "
       "the abundance of transcripts.  This option can be used with "
       "protocols "
       "where one expects that fragments derive from their underlying "
       "targets "
       "without regard to that target's length (e.g. QuantSeq)")
      (
       "noEffectiveLengthCorrection",
       po::bool_switch(&(sopt.noEffectiveLengthCorrection))
       ->default_value(salmon::defaults::noEffectiveLengthCorrection),
       "Disables "
       "effective length correction when computing the "
       "probability that a fragment was generated "
       "from a transcript.  If this flag is passed in, the "
       "fragment length distribution is not taken "
       "into account when computing this probability.")
      ("noFragLengthDist",
       po::bool_switch(&(sopt.noFragLengthDist))->default_value(salmon::defaults::noFragLengthDist),
       "[experimental] : "
       "Don't consider concordance with the learned fragment length "
       "distribution when trying to determine "
       "the probability that a fragment has originated from a specified "
       "location.  Normally, Fragments with "
       "unlikely lengths will be assigned a smaller relative probability "
       "than those with more likely "
       "lengths.  When this flag is passed in, the observed fragment length "
       "has no effect on that fragment's "
       "a priori probability.")
      ("noBiasLengthThreshold",
       po::bool_switch(&(sopt.noBiasLengthThreshold))->default_value(salmon::defaults::noBiasLengthThreshold),
       "[experimental] : "
       "If this option is enabled, then no (lower) threshold will be set on "
       "how short bias correction can make effective lengths. This can "
       "increase the precision "
       "of bias correction, but harm robustness.  The default correction "
       "applies a threshold.")
      ("numBiasSamples",
       po::value<int32_t>(&numBiasSamples)->default_value(salmon::defaults::numBiasSamples),
       "Number of fragment mappings to use when learning the "
       "sequence-specific bias model.")
      ("numAuxModelSamples",
       po::value<uint32_t>(&(sopt.numBurninFrags))->default_value(salmon::defaults::numBurninFrags),
       "The first <numAuxModelSamples> are used to train the "
       "auxiliary model parameters (e.g. fragment length distribution, "
       "bias, etc.).  After ther first <numAuxModelSamples> observations "
       "the auxiliary model parameters will be assumed to have converged "
       "and will be fixed.")
      ("numPreAuxModelSamples",
       po::value<uint32_t>(&(sopt.numPreBurninFrags))
       ->default_value(salmon::defaults::numPreBurninFrags),
       "The first <numPreAuxModelSamples> will have their "
       "assignment likelihoods and contributions to the transcript "
       "abundances computed without applying any auxiliary models.  The "
       "purpose "
       "of ignoring the auxiliary models for the first "
       "<numPreAuxModelSamples> observations is to avoid applying these "
       "models before thier "
       "parameters have been learned sufficiently well.")
      ("useEM", po::bool_switch(&(sopt.useEM))->default_value(salmon::defaults::useEM),
       "Use the traditional EM algorithm for optimization in the batch passes.")
      ("useVBOpt", po::bool_switch(&(sopt.useVBOpt))->default_value(salmon::defaults::useVBOpt),
       "Use the Variational Bayesian EM [default]")
      ("rangeFactorizationBins",
       po::value<uint32_t>(&(sopt.rangeFactorizationBins))->default_value(salmon::defaults::rangeFactorizationBins),
       "Factorizes the likelihood used in quantification by adopting a new "
       "notion of equivalence classes based on "
       "the conditional probabilities with which fragments are generated "
       "from different transcripts.  This is a more "
       "fine-grained factorization than the normal rich equivalence "
       "classes.  The default value (0) corresponds to "
       "the standard rich equivalence classes, and larger values imply a "
       "more fine-grained factorization.  If range factorization "
       "is enabled, a common value to select for this parameter is 4.")
      ("numGibbsSamples",
       po::value<uint32_t>(&(sopt.numGibbsSamples))->default_value(salmon::defaults::numGibbsSamples),
       "Number of Gibbs sampling rounds to "
       "perform.")
      ("noGammaDraw",
       po::bool_switch(&(sopt.noGammaDraw))->default_value(salmon::defaults::noGammaDraw),
       "This switch will disable drawing transcript fractions from a Gamma distribution during Gibbs sampling.  In this case "
       "the sampler does not account for shot-noise, but only assignment ambiguity")
      ("numBootstraps",
       po::value<uint32_t>(&(sopt.numBootstraps))->default_value(salmon::defaults::numBootstraps),
       "Number of bootstrap samples to generate. Note: "
       "This is mutually exclusive with Gibbs sampling.")
      ("bootstrapReproject",
       po::bool_switch(&(sopt.bootstrapReproject))->default_value(salmon::defaults::noGammaDraw),
       "This switch will learn the parameter distribution from the bootstrapped counts for each sample, but "
        "will reproject those parameters onto the original equivalence class counts.")
      ("thinningFactor",
       po::value<uint32_t>(&(sopt.thinningFactor))->default_value(salmon::defaults::thinningFactor),
       "Number of steps to discard for every sample kept from the Gibbs "
       "chain. "
       "The larger this number, the less chance that subsequent samples are "
       "auto-correlated, but the slower sampling becomes.")
      ("quiet,q", po::bool_switch(&(sopt.quiet))->default_value(salmon::defaults::quiet),
       "Be quiet while doing quantification (don't write informative "
       "output to the console unless something goes wrong).")
      ("perTranscriptPrior", po::bool_switch(&(sopt.perTranscriptPrior))->default_value(salmon::defaults::perTranscriptPrior),
       "The "
       "prior (either the default or the argument provided via --vbPrior) "
       "will "
       "be interpreted as a transcript-level prior (i.e. each transcript "
       "will "
       "be given a prior read count of this value)")
      ("sigDigits", po::value<uint32_t>(&(sopt.sigDigits))->default_value(salmon::defaults::sigDigits),
       "The number of significant digits to write when outputting the EffectiveLength and NumReads columns")
      ("vbPrior", po::value<double>(&(sopt.vbPrior))->default_value(salmon::defaults::vbPrior),
       "The prior that will be used in the VBEM algorithm.  This is "
       "interpreted "
       "as a per-nucleotide prior, unless the --perTranscriptPrior flag "
       "is also given, in which case this is used as a transcript-level "
       "prior")
      ("writeOrphanLinks",
       po::bool_switch(&(sopt.writeOrphanLinks))->default_value(salmon::defaults::writeOrphanLinks),
       "Write the transcripts that are linked by orphaned reads.")
      ("writeUnmappedNames",
       po::bool_switch(&(sopt.writeUnmappedNames))->default_value(salmon::defaults::writeUnmappedNames),
       "Write the names of un-mapped reads to the file unmapped_names.txt "
       "in the auxiliary directory.");
    return advanced;
  }

  po::options_description ProgramOptionsGenerator::getHiddenOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    po::options_description hidden("\nhidden options");
    hidden.add_options()
      ("numGCBins", po::value<size_t>(&(sopt.numFragGCBins))->default_value(salmon::defaults::numFragGCBins),
       "Number of bins to use when modeling fragment GC bias")
      ("conditionalGCBins",
       po::value<size_t>(&(sopt.numConditionalGCBins))->default_value(salmon::defaults::numConditionalGCBins),
       "Number of different fragment GC models to learn based on read start/end "
       "context")
      ("numRequiredObs,n",
       po::value(&(sopt.numRequiredFragments))->default_value(salmon::defaults::numRequiredFrags),
       "[Deprecated]: The minimum number of observations (mapped reads) "
       "that must be observed before "
       "the inference procedure will terminate.")
      ("maxHashResizeThreads", po::value<uint32_t>(&(sopt.maxHashResizeThreads))->default_value(salmon::defaults::maxHashResizeThreads),
       "Maximum number of threads to allow cuckoo hash map to use when / if it resizes");
    return hidden;
  }

  po::options_description ProgramOptionsGenerator::getTestingOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    po::options_description testing("\n"
                                    "testing options");
    testing.add_options()
      ("posBias", po::bool_switch(&(sopt.posBiasCorrect))->default_value(salmon::defaults::posBiasCorrect),
                          "[experimental] Perform positional bias correction")
      ("noRichEqClasses",
      po::bool_switch(&(sopt.noRichEqClasses))->default_value(salmon::defaults::noRichEqClasses),
        "[TESTING OPTION]: Disable \"rich\" equivalent classes.  If this flag is "
        "passed, then "
        "all information about the relative weights for each transcript in the "
        "label of an equivalence class will be ignored, and only the relative "
        "abundance and effective length of each transcript will be considered.")
      ("noFragLenFactor",
      po::bool_switch(&(sopt.noFragLenFactor))->default_value(salmon::defaults::noFragLengthFactor),
        "[TESTING OPTION]: Disable the factor in the likelihood that takes into "
        "account the "
        "goodness-of-fit of an alignment with the empirical fragment length "
        "distribution")
      ("rankEqClasses",
      po::bool_switch(&(sopt.rankEqClasses))->default_value(salmon::defaults::rankEqClasses),
        "[TESTING OPTION]: Keep separate equivalence classes for each distinct "
        "ordering of transcripts in the label.")
      ("noExtrapolateCounts",
      po::bool_switch(&(sopt.dontExtrapolateCounts))->default_value(salmon::defaults::dontExtrapolateCounts),
        "[TESTING OPTION]: When generating posterior counts for Gibbs sampling, "
        "use the directly re-allocated counts in each iteration, rather than "
        "extrapolating "
        "from transcript fractions.");
    return testing;
  }

  po::options_description ProgramOptionsGenerator::getDeprecatedOptions(SalmonOpts& sopt) {
    using std::string;
    using std::vector;

    po::options_description deprecated(
        "\ndeprecated options about which to inform the user");
    return deprecated;
  }

}
