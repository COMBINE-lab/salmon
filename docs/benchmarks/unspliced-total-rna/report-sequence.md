### Per library

| mode     | run        |   num_processed |   mapping_rate |   decoy_fraction |   unspliced_target_fraction |   retained_intron_share |
|:---------|:-----------|----------------:|---------------:|-----------------:|----------------------------:|------------------------:|
| gentrome | SRR5693352 |        13657343 |         0.6412 |           0.3051 |                      0.0096 |                  0.0470 |
| gentrome | SRR5693353 |        15958033 |         0.6562 |           0.2880 |                      0.0097 |                  0.0461 |
| gentrome | SRR5693399 |        10298793 |         0.6961 |           0.2389 |                      0.0099 |                  0.0425 |
| gentrome | SRR5693400 |        10266590 |         0.6693 |           0.2635 |                      0.0095 |                  0.0455 |
| intron   | SRR5693352 |        13657343 |         0.9126 |           0.0159 |                      0.3067 |                  0.0324 |
| intron   | SRR5693353 |        15958033 |         0.9119 |           0.0152 |                      0.2901 |                  0.0321 |
| intron   | SRR5693399 |        10298793 |         0.9161 |           0.0145 |                      0.2521 |                  0.0305 |
| intron   | SRR5693400 |        10266590 |         0.9126 |           0.0156 |                      0.2769 |                  0.0319 |
| premrna  | SRR5693352 |        13657343 |         0.9139 |           0.0147 |                      0.3163 |                  0.0334 |
| premrna  | SRR5693353 |        15958033 |         0.9132 |           0.0140 |                      0.2994 |                  0.0330 |
| premrna  | SRR5693399 |        10298793 |         0.9176 |           0.0131 |                      0.2604 |                  0.0311 |
| premrna  | SRR5693400 |        10266590 |         0.9141 |           0.0143 |                      0.2856 |                  0.0329 |

### Retained-intron share vs decoy share, across libraries

- gentrome: Spearman rho = 1.00 (n = 4)
- intron: Spearman rho = 0.80 (n = 4)
- premrna: Spearman rho = 0.80 (n = 4)

### Gene avgTxLength between technical replicate libraries

Genes with > 250 spliced reads in every library under every mode: 5042.

| mode | genes | SD log2 ratio, SRR5693352/SRR5693353 | SD log2 ratio, SRR5693399/SRR5693400 | r of site difference, replicate a vs b |
|---|---|---|---|---|
| gentrome | 5042 | 0.1392 | 0.1937 | 0.201 |
| intron | 5042 | 0.1577 | 0.2079 | 0.191 |
| premrna | 5042 | 0.1588 | 0.2105 | 0.196 |

### Time and memory

| step                      |   wall_s |   max_rss_GiB |
|:--------------------------|---------:|--------------:|
| index_gentrome            |      422 |          24.7 |
| index_intron              |      638 |          35.6 |
| index_premrna             |      550 |          38.4 |
| quant_gentrome_SRR5693352 |      203 |          13.2 |
| quant_intron_SRR5693352   |      299 |          18.8 |
| quant_premrna_SRR5693352  |      346 |          19.0 |
| quant_gentrome_SRR5693353 |      382 |          13.4 |
| quant_intron_SRR5693353   |      350 |          18.9 |
| quant_premrna_SRR5693353  |      440 |          18.9 |
| quant_gentrome_SRR5693399 |      226 |          13.2 |
| quant_intron_SRR5693399   |      467 |          18.7 |
| quant_premrna_SRR5693399  |      341 |          18.8 |
| quant_gentrome_SRR5693400 |      172 |          13.2 |
| quant_intron_SRR5693400   |      183 |          18.8 |
| quant_premrna_SRR5693400  |      389 |          18.8 |
- index_gentrome: 9.8 GiB on disk
- index_intron: 14.2 GiB on disk
- index_premrna: 14.3 GiB on disk
