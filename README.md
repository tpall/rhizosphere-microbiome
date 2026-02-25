# rhizosphere-microbiome

Analysis of microbial community composition in *Salix alba* cuttings supplemented with *Salicornia europaea* rhizosphere microbiome.

## Study design

Soil samples were collected from three zones (Z1, Z2, Z3) with or without *S. europaea* rhizosphere supplement (+rhiz / -rhiz). 16S rRNA amplicon sequencing was processed through mothur and exported as BIOM files.

**Dataset summary (before rarefaction):**

| Group | OTUs |
|---|---|
| All samples | 10,137 |
| Without supplement (-rhiz) | 6,697 |
| With supplement (+rhiz) | 8,693 |
| Total sequences | 1,439,201 |

**After rarefaction (`rngseed = 11`):**

| Group | OTUs |
|---|---|
| All samples | 8,255 |
| Without supplement (-rhiz) | 3,231 |
| With supplement (+rhiz) | 4,380 |

## Analysis (`scripts/main.R`)

1. **Data import & preprocessing** — imports BIOM file via `phyloseq`, validates, fixes taxonomy, subsets to zone samples, and adds sample metadata (zone, rhizo treatment).

2. **OTU overview** — compositional barplot of top 10 genera across all samples.

3. **Graph-based permutation test** — tests community composition differences by zone/treatment using Jaccard distance kNN graph (`phyloseqGraphTest`).

4. **Venn diagrams** — OTU overlap across zones and between treatment groups:
   - Common OTUs across all three unsupplemented zones: **321**
   - Common OTUs across all three supplemented zones: **346**
   - OTUs common to all six groups: **198**
   - Bayesian binomial models (`brms`) test whether supplementation changes the proportion of shared OTUs within and across zones.

5. **Alpha diversity** — Observed richness, Shannon index, and Inverse Simpson index modelled with Bayesian regression (`brms`, negative binomial / Student-t families), accounting for library size and zone × rhizo interactions.

6. **Beta diversity** — NMDS ordination (Bray-Curtis dissimilarities) and distances to group centroids modelled with Bayesian location-scale regression. NMDS shows a gradient across zones; no clear within-zone separation by rhizosphere treatment.

7. **Taxonomic composition** — CLR-transformed heatmap of top 20 genera (highlighting halophilic taxa: *Halomonadaceae*, *Idiomarina*, *Halomonas*, *Marinobacteria*), and compositional barplots at phylum level and within Proteobacteria.

## Outputs

| File | Description |
|---|---|
| `plots/venn_diagrams.tiff` | Venn diagram panel (A–C + model estimates) |
| `plots/alpha_diversity.tiff` | Observed, Shannon, InvSimpson panel |
| `plots/beta_diversity.tiff` | NMDS + distance-to-centroid panel |
| `plots/heatmap_top20.tiff` | CLR heatmap of top 20 genera |
| `plots/comp_barplot.png` | Phylum and Proteobacteria barplots |
| `models/` | Cached `brms` model fits |

## Dependencies

```r
tidyverse, phyloseq, microbiome, microViz, phyloseqGraphTest,
brms, tidybayes, modelr,
ggvenn, ggVennDiagram, patchwork, cowplot, ggsci,
centroidr, here
```

## Usage

```r
source("scripts/main.R")
```

Data are downloaded automatically from the Galaxy HPC server on first run and cached to `results/`.
