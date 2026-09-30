### Per library

| mode     | run        |   num_processed |   mapping_rate |   decoy_fraction |   unspliced_target_fraction |   retained_intron_share |
|:---------|:-----------|----------------:|---------------:|-----------------:|----------------------------:|------------------------:|
| gentrome | SRR5693352 |        13657343 |         0.6412 |           0.3051 |                      0.0096 |                  0.0470 |
| gentrome | SRR5693353 |        15958033 |         0.6562 |           0.2880 |                      0.0097 |                  0.0461 |
| gentrome | SRR5693399 |        10298793 |         0.6961 |           0.2389 |                      0.0099 |                  0.0425 |
| gentrome | SRR5693400 |        10266590 |         0.6693 |           0.2635 |                      0.0095 |                  0.0455 |
| intron   | SRR5693352 |        13657343 |         0.9125 |           0.0160 |                      0.3123 |                  0.0311 |
| intron   | SRR5693353 |        15958033 |         0.9118 |           0.0152 |                      0.2951 |                  0.0311 |
| intron   | SRR5693399 |        10298793 |         0.9160 |           0.0145 |                      0.2548 |                  0.0294 |
| intron   | SRR5693400 |        10266590 |         0.9126 |           0.0156 |                      0.2808 |                  0.0304 |
| premrna  | SRR5693352 |        13657343 |         0.9139 |           0.0148 |                      0.3190 |                  0.0320 |
| premrna  | SRR5693353 |        15958033 |         0.9132 |           0.0140 |                      0.3015 |                  0.0318 |
| premrna  | SRR5693399 |        10298793 |         0.9176 |           0.0131 |                      0.2626 |                  0.0301 |
| premrna  | SRR5693400 |        10266590 |         0.9141 |           0.0143 |                      0.2889 |                  0.0313 |

### Retained-intron share vs decoy share, across libraries

- gentrome: Spearman rho = 1.00 (n = 4)
- intron: Spearman rho = 0.80 (n = 4)
- premrna: Spearman rho = 0.80 (n = 4)

### Gene avgTxLength between technical replicate libraries

Genes with > 250 spliced reads in every library under every mode: 5024.

| mode | genes | SD log2 ratio, SRR5693352/SRR5693353 | SD log2 ratio, SRR5693399/SRR5693400 | r of site difference, replicate a vs b |
|---|---|---|---|---|
| gentrome | 5024 | 0.1392 | 0.1946 | 0.202 |
| intron | 5024 | 0.1679 | 0.2237 | 0.191 |
| premrna | 5024 | 0.1692 | 0.2263 | 0.188 |

### Time and memory

| step                      |   wall_s |   max_rss_GiB |
|:--------------------------|---------:|--------------:|
| index_intron              |      580 |          31.5 |
| index_premrna             |      569 |          33.9 |
| quant_gentrome_SRR5693352 |      361 |          13.2 |
| quant_intron_SRR5693352   |      385 |          18.2 |
| quant_premrna_SRR5693352  |      343 |          18.3 |
| quant_gentrome_SRR5693353 |      356 |          13.4 |
| quant_intron_SRR5693353   |      316 |          18.4 |
| quant_premrna_SRR5693353  |      296 |          18.5 |
| quant_gentrome_SRR5693399 |      275 |          13.1 |
| quant_intron_SRR5693399   |      269 |          18.1 |
| quant_premrna_SRR5693399  |      290 |          18.3 |
| quant_gentrome_SRR5693400 |      267 |          13.2 |
| quant_intron_SRR5693400   |      258 |          18.1 |
| quant_premrna_SRR5693400  |      266 |          18.2 |
- index_gentrome: 9.8 GiB on disk
- index_intron: 12.1 GiB on disk
- index_premrna: 12.2 GiB on disk
