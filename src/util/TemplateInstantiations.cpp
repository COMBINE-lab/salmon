#include "salmon/internal/alignment/AlignmentLibrary.hpp"
#include "salmon/internal/alignment/AlignmentModel.hpp"
#include "salmon/internal/alignment/ONTAlignmentModel.hpp"
#include "salmon/internal/alignment/ReadPair.hpp"
#include "salmon/internal/alignment/UnpairedRead.hpp"
#include "salmon/internal/quant/ReadExperiment.hpp"

template class ReadExperiment<EquivalenceClassBuilder<TGValue>>;
template class ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;

template class AlignmentLibrary<ReadPair, EquivalenceClassBuilder<TGValue>,
                                AlignmentModel>;
template class AlignmentLibrary<UnpairedRead, EquivalenceClassBuilder<TGValue>,
                                AlignmentModel>;
template class AlignmentLibrary<UnpairedRead, EquivalenceClassBuilder<TGValue>,
                                ONTAlignmentModel>;
