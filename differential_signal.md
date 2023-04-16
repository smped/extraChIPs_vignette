extraChIPs: Differential Signal Analysis
================
true

# Introduction

The [GRAVI](https://github.com/smped/GRAVI) workflow, for which this
package is designed, uses sliding windows for differential signal
analysis in a manner similar to the package `csaw`, but also
incorporating `macs2` peaks. The workflow itself extends to integrating
multiple ChIP targets and external data sources, and as such, this
package introduces a handful of functions to simplify and enable these
analyses. Whilst many existing approaches refer to this type of analysis
as Differential Binding analysis, we prefer the term *Differential
Signal Analysis* as this more accurately captures the range of ChIP
targets which are likely to be investigated.

The majority of examples below use heavily reduced datasets to provide
general guidance on using the functions. Some results may appear trivial
as a result, but will hopefully prove far more useful in a true
experimental context. All data, along with this vignette are available
[here](https://github.com/smped/extraChIPs_vignette). Please place all
contents of the data directory in a directory named data in your own
working directory.

# Setup

## Installation

In order to use the package `extraChIPs` and follow this vignette, we
recommend using the package `BiocManager` hosted on CRAN. Once this is
installed, the additional packages required for this vignette
(`tidyverse`, `Rsamtools`, `csaw`, `BiocParallel` and `rtracklayer`) can
also be installed.

``` r
if (!"BiocManager" %in% rownames(installed.packages()))
  install.packages("BiocManager")
pkg <- c(
  "tidyverse", "Rsamtools", "csaw", "BiocParallel", "rtracklayer", "edgeR", 
  "patchwork", "extraChIPs", "plyranges", "scales", "ggside"
)
BiocManager::install(pkg, update = FALSE)
```

Once these packages are installed, we can load them easily

``` r
library(tidyverse)
library(Rsamtools)
library(csaw)
library(BiocParallel)
library(rtracklayer)
library(edgeR)
library(patchwork)
library(extraChIPs)
library(plyranges)
library(scales)
library(ggside)
```

## Data

All data for this vignette is expected to be in a sub-directory of the
working directory named “data”, and all paths will be predicated on
this. Please ensure you have all data in this location, obtained from
[here](https://github.com/smped/extraChIPs_vignette).

The data itself is ChIP-Seq data targeting the histone mark H3K27ac, and
is taken from the cell-line MDA-MB-453 under Vehicle and DHT-stimulated
conditions. Using CRCh37 as the reference genome, a subset of regions
found on chromosome 10 are included in this dataset for simplicity.

`Seqinfo` objects are the foundation of working with GRanges, so let’s
define a suitable object for consistency throughout th analysis.

``` r
hg19 <- GenomeInfoDb::getChromInfoFromUCSC("hg19")
grch37 <- hg19 %>% 
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y"))) %>% 
  mutate(genome = "GRCh37") %>% 
  dplyr::select(
    seqnames = chrom, seqlengths = size, isCircular = circular, genome
  ) %>% 
  as("Seqinfo")
```

# Working With Peaks

The provided dataset includes six files produced by `macs2 callpeak`
(Zhang et al. 2008) in the `narrowPeak` format, and these are able to be
easily parsed using `extraChIPs`.

``` r
peakFiles <- list.files("data", pattern = "narrowPeak", full.names = TRUE)
peaks <- importPeaks(peakFiles, seqinfo = grch37)
```

This will import the peaks from all files as a single `GRangesList`
object, adding the file-name to each element by default. We can easily
modify these names if we so wish.

``` r
names(peaks) <- str_remove_all(names(peaks), "_peaks.narrowPeak")
```

Once loaded, we can easily check how similar our replicates are.

``` r
plotOverlaps(peaks, min_size = 10, .sort_sets = FALSE)
```

![*UpSet plot showing overlapping peaks across all
replicates*](differential_signal_files/figure-gfm/plot-overlaps-1.png)

Optionally, specifying a column and a suitable function will produce an
additional panel summarising that value. In the following, we’ll show
the maximum score obtained, highlighting that for peaks identified in
only one or two replicates, the overall signal intensity is generally
lower.

``` r
plotOverlaps(peaks, min_size = 10, .sort_sets = FALSE, var = "score", f = "max")
```

![*UpSet plot showing overlapping peaks across all replicates, with the
maximum score across all replicates shown in the upper
panel.*](differential_signal_files/figure-gfm/plot-overlaps-score-1.png)

A common task at this point may be to define consensus peaks within each
treatment group, by retaining only the peaks found in 2 of the 3
replicates. The default approach is to take the union of all ranges,
with the returned object containing logical values for each sample, as
well as the number of samples where an overlapping peak was found.

If we wish to retain any of the original columns, such as the
`macs2 callpeak` score, we can simply pass the column names to
`makeConsensus()`

``` r
consensus_veh <- peaks %>% 
  .[str_detect(names(.), "Veh")] %>% 
  makeConsensus(p = 2/3, var = "score")
consensus_dht <- peaks %>% 
  .[str_detect(names(.), "DHT")] %>% 
  makeConsensus(p = 2 / 3, var = "score")
```

Alternatively, we could find the centre of the peaks as part of this
process, by averaging across the estimated peak centres for each sample.
Whilst this is very common for *transcription factor* peaks, this may be
less informative for other types of ChIP targets, such as the histone
marks we have.

``` r
consensus_veh <- peaks %>% 
  .[str_detect(names(.), "Veh")] %>% 
  endoapply(mutate, centre = start + peak) %>% 
  makeConsensus(p = 2/3, var = "centre") %>% 
  mutate(centre = vapply(centre, mean, numeric(1)))
consensus_dht <- peaks %>% 
  .[str_detect(names(.), "DHT")] %>% 
  endoapply(mutate, centre = start + peak) %>% 
  makeConsensus(p = 2/3, var = "centre") %>% 
  mutate(centre = vapply(centre, mean, numeric(1)))
```

We can also inspect these using `plotOverlaps()` provided we use a
`GRangesList` for the input.

``` r
GRangesList(
  Veh = granges(consensus_veh), DHT = granges(consensus_dht)
) %>% 
  plotOverlaps(set_col = c("grey70", "red"))
```

![*Overlap between consensus peaks identified in a treatment-specific
manner*](differential_signal_files/figure-gfm/venn-consensus-overlap-1.png)

We could go one step further and define the set of peaks found in either
treatment. Given we’re being inclusive here, we can leave p = 0 so any
peak found in either treatment is included.

``` r
all_consensus <- GRangesList(
  Veh = select(consensus_veh, centre), DHT = select(consensus_dht, centre)
) %>% 
  makeConsensus(var = "centre") %>% 
  mutate(centre = vapply(centre, mean, numeric(1)))
```

# Differential Signal Analysis

## Sliding Windows

The standard approach of tools such as DiffBind (Ross-Innes et al. 2012)
is to take a set of peaks, re-centre them, then set all regions to be
the same width. From there, data is passed to `edgeR` (Chen, Lun, and
Smyth 2016) or `DESeq2` (Love, Huber, and Anders 2014) for analysis.
Whilst this approach can be replicated using the defined peaks above,
the suggested approach for `extraChIPs` is to us sliding windows, as per
`csaw`. The resultant *variable width regions* can be particularly
advantageous for ChIP targets such as H3K27ac where regions marked by
histone-acetylation can vary greatly in size.

The starting point for differential signal analyses using `extraChIPs`
is to define a set of sliding windows across the genome, then count
reads from a set of bam files, defined as a `BamFileList.` Commonly one
or more IP input/control samples is also produced during a ChIP-Seq
experiment, and these should be included at this stage of the analysis.
The example files provided here contain a small subset of reads from
chromosome 10 across two experimental conditions and one input sample,
and we will define them all as a `BamFileList`.

``` r
bfl <- list.files("data", pattern = "bam$", full.names = TRUE) %>% 
  BamFileList()
names(bfl) <- str_remove_all(names(bfl), ".bam")
```

**NB:** It should also be noted that counting all reads across a
`BamFileList` using sliding windows, **will require a significant amount
of RAM** and will be beyond the capacity of most laptops as of the time
of writing. When working with complete datasets, this step is best
performed on an HPC or a similar interactive server.

The approach taken below is to first define a set of sliding windows
across the genome, using the capabilities of `csaw`. After counting
reads across all windows, a set of pre-defined regions is then used
guide the function `dualFilter()` which will discard low-signal windows,
retaining only those a) above a minimum signal level and b) with signal
notably above that of any input samples. These regions can be obtained
from any external resource, or can even be taken from `macs2`-defined
peaks from the same samples.

First we can define our windows and count the alignments using the
existing capabilities and functions provided in the `csaw` package (Lun
and Smyth 2016). In the following, we’ll use a sliding window of 120bp
and a step size of 40bp, meaning each nucleotide is covered by 3
windows. In addition, we’ll exclude blacklisted and greylisted regions
as provided in the dataset. These can be obtained easily by using the
`GreyListChIP` package, which is beyond the scope of this vignette.

``` r
greylist <- import.bed("data/chr10_greylist_subset.bed", seqinfo = grch37)
blacklist <- import.bed("data/chr10_blacklist_subset.bed", seqinfo = grch37)
rp <- readParam(
  pe = "none",
  dedup = TRUE,
  restrict = "chr10",
  discard = c(greylist, blacklist)
)
wincounts <- windowCounts(
  bam.files = bfl,
  spacing = 40,
  width = 120,
  ext = 200,
  filter = length(bfl),
  param = rp
)
```

This produces a `RangesSummarizedExperiment` with windows included which
passed the minimum threshold of 7 total reads. We can check which
windows passed this threshold using `rowRanges()`

``` r
rowRanges(wincounts)
```

    ## GRanges object with 262151 ranges and 0 metadata columns:
    ##            seqnames            ranges strand
    ##               <Rle>         <IRanges>  <Rle>
    ##        [1]    chr10 42820161-42820280      *
    ##        [2]    chr10 42820201-42820320      *
    ##        [3]    chr10 42820241-42820360      *
    ##        [4]    chr10 42820281-42820400      *
    ##        [5]    chr10 42821001-42821120      *
    ##        ...      ...               ...    ...
    ##   [262147]    chr10 99998721-99998840      *
    ##   [262148]    chr10 99998761-99998880      *
    ##   [262149]    chr10 99998801-99998920      *
    ##   [262150]    chr10 99998841-99998960      *
    ##   [262151]    chr10 99998881-99999000      *
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome

We can also add some key information to the `colData` element of this
object, which will also be propagated to all downstream objects.

``` r
wincounts$sample <- names(bfl)
wincounts$treat <- str_extract(names(bfl), "(Veh|DHT)") %>% 
  fct(levels = c("Veh", "DHT"))
colData(wincounts)
```

    ## DataFrame with 7 rows and 6 columns
    ##            bam.files    totals       ext      rlen      sample    treat
    ##          <character> <integer> <integer> <integer> <character> <factor>
    ## DHT_1 data/DHT_1.bam    252973       200        74       DHT_1      DHT
    ## DHT_2 data/DHT_2.bam    295028       200        74       DHT_2      DHT
    ## DHT_3 data/DHT_3.bam    295105       200        74       DHT_3      DHT
    ## Input data/Input.bam     82276       200        72       Input      NA 
    ## Veh_1 data/Veh_1.bam    250970       200        74       Veh_1      Veh
    ## Veh_2 data/Veh_2.bam    280898       200        74       Veh_2      Veh
    ## Veh_3 data/Veh_3.bam    294356       200        74       Veh_3      Veh

A density plot can be simply drawn of these counts, with the vast
majority of windows receiving very low counts, due to the nature of
transcription factor binding, where long stretches are unbound. The
windows with higher counts tend to be associated with the samples
targeting a transcription factor (TF), as seen in the two treatment
group samples.

``` r
plotAssayDensities(wincounts, colour = "treat", trans = "log1p") +
  theme_bw()
```

![*Read Densities for all returned windows across all
samples*](differential_signal_files/figure-gfm/plot-densities-1.png)

## Filtering of Sliding Windows

After counting all reads in the sliding genomic windows, the next step
is to discard windows for which counts are unlikely to represent true
signal from our ChIP target. The strategy employed in `extraChIPs` uses
a set of consensus peaks to automatically set thresholds based on 1)
counts strongly above the counts from the input sample, and 2) the
windows with the overall highest signal. Thresholds are determined such
that a proportion (e.g. `q = 0.5`) of the windows which overlap one of
the supplied consensus peaks will be returned. Higher values for `q`
will return more windows, however many of these will tend to only
marginally overlap a peak in one of the tail regions, and these will
most likely be covered by neighbouring windows. Experience has shown
that values such as `q = 0.5` tend to return a considerable proportion
of windows containing true signal from the ChIP target.

The we can pass these to the function `dualFilter()` which utilises the
strategy described above. On large datasets, this can be quite
time-consuming, as can the initial counting step. Multiple alternative
filtering strategies are also provided by the package `csaw` and these
can be accessed using `?csaw::filterWindows`

``` r
filtcounts <- dualFilter(
  x = wincounts[, !is.na(wincounts$treat)],
  bg = wincounts[, is.na(wincounts$treat)], 
  ref = all_consensus,
  q = 0.6
)
```

Thus we have reduced our initial set of 262,151 sliding windows to the
16,141 windows most likely to contain true signal from our ChIP target.
The returned object will by default contain `counts` and `logCPM`
assays, with the complete library sizes used for the calculation of
`logCPM` values. Similarly, *the input sample is no longer included* in
the data object, although additional columns can easily be added to the
returned object using any number of strategies.

``` r
dim(wincounts)
```

    ## [1] 262151      7

``` r
dim(filtcounts)
```

    ## [1] 16141     6

``` r
assays(filtcounts)
```

    ## List of length 2
    ## names(2): counts logCPM

We can once again check our signal distributions, this time on the
logCPM values.

``` r
plotAssayDensities(filtcounts, assay = "logCPM", colour = "treat") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw()
```

![*Densities for logCPM values across all samples after discarding
windows less likely to contain H3K27ac
signal*](differential_signal_files/figure-gfm/plotcpm-1.png)

The `rowData` element of the returned object will contain a logical
column indicating where each specific retained window overlapped one of
the supplied consensus peaks.

``` r
rowRanges(filtcounts)
```

    ## GRanges object with 16141 ranges and 1 metadata column:
    ##           seqnames            ranges strand | overlaps_ref
    ##              <Rle>         <IRanges>  <Rle> |    <logical>
    ##       [1]    chr10 43047561-43047680      * |         TRUE
    ##       [2]    chr10 43047601-43047720      * |         TRUE
    ##       [3]    chr10 43047641-43047760      * |         TRUE
    ##       [4]    chr10 43047681-43047800      * |         TRUE
    ##       [5]    chr10 43047721-43047840      * |         TRUE
    ##       ...      ...               ...    ... .          ...
    ##   [16137]    chr10 99894961-99895080      * |         TRUE
    ##   [16138]    chr10 99895001-99895120      * |         TRUE
    ##   [16139]    chr10 99895041-99895160      * |         TRUE
    ##   [16140]    chr10 99895081-99895200      * |         TRUE
    ##   [16141]    chr10 99895121-99895240      * |         TRUE
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome

``` r
mean(rowRanges(filtcounts)$overlaps_ref)
```

    ## [1] 0.9997522

## Initial Visualisation

Inspecting your data is a common first step, and a common QC step is
Relative Log-Expression (RLE) (Gandolfo and Speed 2018). In the
following, we first inspect the RLE across the entire dataset, followed
by RLE grouping *within treatments*. This can be particularly useful
when distributions vary significantly between treatment groups, such as
may occur with a cytoplasmic to nuclear shift by a given ChIP target.
Here, however, there is minimal difference between the two approaches as
H3K27ac signal tends to be broadly consistent between these treatment
groups.

``` r
a <- plotAssayRle(filtcounts, assay = "logCPM", fill = "treat") +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("RLE: Across All Samples") +
  theme_bw()
b <- plotAssayRle(
  filtcounts, assay = "logCPM", fill = "treat", rle_group = "treat"
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("RLE: Within Treatment Groups") +
  theme_bw()
a + b + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A")
```

![*RLE plots across all samples (A) and with values calculated within
treatment groups
(B).*](differential_signal_files/figure-gfm/plot-assay-rle-1.png)

We can also check the samples using a PCA plot, again colouring the
points by treatment group and adding labels, which will repel by default
if the points are shown.

``` r
plotAssayPCA(filtcounts, "logCPM", colour = "treat", label = "sample") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw()
```

![*PCA plot based on the logCPM
assay*](differential_signal_files/figure-gfm/plot-pca-1.png)

## Statistical Testing

Multiple methods are enabled in the package `extraChIPs` via the
function `fitAssayDiff()`, with the possibility of incorporating any
additional normalisation strategies from external packages. The two
basic strategies are 1) Quasi-Likelihood Fits (Lund et al. 2012) and 2)
`limma-trend` (Law et al. 2014). This first (method = “qlf”) uses counts
with any of the provided normalisation strategies from
`edgeR::calcNormFactors()`, and setting the normalisation method to
“none” is the equivalent of library-size normalisation, which replicates
the default normalisation strategy from DiffBind (Ross-Innes et al.
2012). If choosing to normalise within treatment groups, a factor can be
provided via the groups argument, essentially adding this as an option
for all methods provided in `edgeR::calcNormFactors()`. The second
method (method = “lt”) is specifically for logCPM values and these can
be provided as output by `dualFilter()` or may be normalised using any
number of additional methods. In addition to the above methods, a
range-based $H_0$ (McCarthy and Smyth 2009) can be specified by
providing a value to the `fc` or `lfc` arguments.

Here, we’ll fit our data using Quasi-Likelihood Fits, library-size
normalisation and setting a change in signal beyond the range of $\pm$
20% as being of interest. By default, the returned object, will contain
the results from model fitting in the `rowData()` element as these are
result associated with each row element in the `SummarizedExperiment`
object. If the object is a `RangedSummarizedExperiment` object, setting
`asRanges = TRUE` will simply return the set of GRanges along with the
testing results.

``` r
X <- model.matrix(~treat, data = colData(filtcounts))
fit_gr <- fitAssayDiff(filtcounts, design = X, fc = 1.2, asRanges = TRUE)
```

## Merging Windows

After an analysis has been performed, common values contained in the
output may be estimated signal (`logCPM`), estimated change (`logFC`)
with both raw and adjusted p-values. Given the dependency of
neighbouring windows, any adjusted p-values will not be appropriate and
a merging of overlapping and/or neighbouring windows should be
performed. Multiple `csaw` methods are wrapped using `mergeByCol()`,
`mergeBySig()` with minor changes to the returned object, such as the
inclusion of the representative range in the column `keyval_range`.

For this vignette, we’ll merge using the asymptotically exact harmonic
mean p-value, which can also be used for merging dependent p-values
(Wilson 2019). When merging windows using the harmonic mean p-values,
instead of values from a representative window, weighted averages for
the expression and logFC estimates are returned using the weights
$w_i = \frac{1}{p_i}$. A representative window, corresponding to the
original window with the lowest p-value is returned.

``` r
results_gr <- mergeByHMP(fit_gr, inc_cols = "overlaps_ref", merge_within = 120)
results_gr$status <- case_when(
  results_gr$hmp_fdr > 0.05 ~ "Unchanged",
  results_gr$logFC > 0 ~ "Increased",
  results_gr$logFC < 0 ~ "Decreased"
)
arrange(results_gr, hmp)[1:5]
```

    ## GRanges object with 5 ranges and 10 metadata columns:
    ##       seqnames            ranges strand | n_windows      n_up    n_down
    ##          <Rle>         <IRanges>  <Rle> | <integer> <integer> <integer>
    ##   [1]    chr10 79266481-79268400      * |        38         3         0
    ##   [2]    chr10 43689161-43690240      * |        25         4         0
    ##   [3]    chr10 58717481-58717840      * |         7         3         0
    ##   [4]    chr10 67671561-67671800      * |         4         3         0
    ##   [5]    chr10 67670801-67671160      * |         7         1         0
    ##       overlaps_ref            keyval_range    logCPM     logFC         hmp
    ##          <logical>               <GRanges> <numeric> <numeric>   <numeric>
    ##   [1]         TRUE chr10:79266721-79266840   7.60592   2.33531 1.35064e-13
    ##   [2]         TRUE chr10:43689401-43689520   7.43019   2.28678 4.79764e-13
    ##   [3]         TRUE chr10:58717521-58717640   6.75830   2.52756 4.78855e-11
    ##   [4]         TRUE chr10:67671601-67671720   6.64579   2.54747 1.10884e-10
    ##   [5]         TRUE chr10:67670801-67670920   6.56879   2.62167 3.73516e-10
    ##           hmp_fdr      status
    ##         <numeric> <character>
    ##   [1] 7.50956e-11   Increased
    ##   [2] 1.33375e-10   Increased
    ##   [3] 8.87478e-09   Increased
    ##   [4] 1.54129e-08   Increased
    ##   [5] 4.15350e-08   Increased
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome

In the above, we returned 50 ranges which we might consider using the
significance threshold $\alpha$ = 0.05. A particularly beneficial
feature of this approach is that the final ranges will be of highly
variable width, with this select region of chromosome 10 producing
merged windows ranging from 120 to 17360bp, as may be expected for
H3K27ac signal.

We can also quickly check out results using a heatmap across the
retained windows which correspond to our results for differential
signal, although this will not provide any specific genomic context.

``` r
filtcounts %>%
  subsetByOverlaps(subset(results_gr, hmp == min(hmp))) %>%
  plotAssayHeatmap(assay = "logCPM", ysideline = TRUE, yside_col = "treat") +
  scale_fill_viridis_c() +
  scale_colour_brewer(palette = "Set1") +
  scale_ysidex_continuous(expand = expansion(0.1), minor_breaks = NULL) +
  theme_bw()
```

![*Heatmap showing retained sliding windows whih correspond to the most
highly-ranked region for differential
signal.*](differential_signal_files/figure-gfm/plot-assay-heatmap-1.png)

## Mapping of Windows To Genes

Once the changes in signal for our given ChIP target have been
determined, a common next step is to assess which genes are likely to be
impacted. Whilst no definitive, single methodology exists for this
process, the function `mapByFeature()` offers an intuitive approach,
taking into account any previously defined regulatory features. These
regulatory features may be defined by simple proximity to TSS regions,
by histone marks, downloaded from external repositories or any other
possibility. Whilst these features can improve the precision of mapping,
even without these this function can still enable a useful assignment of
target gene to binding event.

The process undertaken inside `mapByFeature()` is a sequential checking
of each range’s association with regulatory features and the most likely
target as a result. These steps are:

1.  **Check for any HiC interactions**

- All genes which directly overlap an interaction anchor are considered
  part of the regulatory network for that interaction, and as such, all
  genes associated with both anchors are assigned to a peak which
  overlaps a HiC Interaction

2.  **Check for any overlaps with a promoter**

- All genes regulated by that promoter are assigned as regulatory
  targets. By default, this is by direct promoter/gene overlap
  (`prom2gene = 0`)

3.  **Check for any overlaps with an enhancer**

- Peaks which overlap an enhancer are assigned to *all* genes within the
  distance specified by `enh2gene` (default = 100kb)

4.  **Check for genes with no previous mappings**

- Peaks *with no previous mappings* are assigned to all directly
  overlapping genes, or the nearest gene within a specified distance
  (default `gr2gene` = 100kb)

As a result, if no promoters, enhancers or long-range interactions are
supplied, all genes will be mapped to peaks using step 4.

A set of annotated regions has been provided for our subset of
chromosome 10, as has the set of genes within this region. These can be
loaded as follows.

``` r
regions <- read_rds("data/chr10_region_subset.rds")
genes <- import.gff("data/chr10_gene_subset.gtf")
```

For our mapping steps, we’ll simply use the `promoters` element as these
are the regions directly overlapping a TSS for all transcripts within
this region of the genome.

``` r
results_gr <- mapByFeature(results_gr, genes = genes, prom = regions$promoters)
```

Now we have our regions showing changed signal along with the likely
regulatory targets. Our top-ranked region is as follows, and this
appears to be associated with the gene *KCNMA1*.

``` r
arrange(results_gr, hmp)[1]
```

    ## GRanges object with 1 range and 12 metadata columns:
    ##       seqnames            ranges strand | n_windows      n_up    n_down
    ##          <Rle>         <IRanges>  <Rle> | <integer> <integer> <integer>
    ##   [1]    chr10 79266481-79268400      * |        38         3         0
    ##       overlaps_ref            keyval_range    logCPM     logFC         hmp
    ##          <logical>               <GRanges> <numeric> <numeric>   <numeric>
    ##   [1]         TRUE chr10:79266721-79266840   7.60592   2.33531 1.35064e-13
    ##           hmp_fdr      status         gene_id       gene_name
    ##         <numeric> <character> <CharacterList> <CharacterList>
    ##   [1] 7.50956e-11   Increased ENSG00000156113          KCNMA1
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome

### Mapping of Windows to Regions

If we also have a set of annotated regions, we can easily map our
regions to the region that it has the greatest proportion of overlap.
First, we’ll arrange our regions to be a single `GRanges` object with
the column “region”. From these we cn simply add the column `bestRegion`

``` r
regions_gr <- regions %>% 
  lapply(select, region) %>% 
  GRangesList() %>% 
  unlist() %>% 
  sort() %>% 
  setNames(NULL)
region_levels <- vapply(regions, function(x) x$region[1], character(1))
results_gr$bestRegion <- results_gr %>% 
  bestOverlap(regions_gr, var = "region") %>% 
  fct(levels = region_levels)
```

# Visualisation of Results

## Association with Annotated Features

The association of windows or peaks with defined features, such as
histone marks or regulatory elements can be important for describing the
binding characteristics of any given transcription factor. We have
already defined the association of the merged windows with annotated
regions, and we can easily visualise these using `plotPie()`.

``` r
results_gr %>% 
  plotPie(fill = "bestRegion")
```

However, given the default plots can often be slightly unsatisfactory,
`plotPie()` is heavily customisable, in particular taking advantage of
the `glue` syntax to customise labels.

``` r
results_gr %>% 
  plotPie(
    fill = "bestRegion", min_p = 0.05, total_size = 5,
    cat_alpha = 0.5,
    cat_glue = "{str_wrap(bestRegion, 10)}\n{n}\n{percent(p, 0.1)}"
  )
```

![*Pie chart showing the best overlapping region for each merged
windows. The best overlap is determined simply by the region with the
lergest amount of
overlap.*](differential_signal_files/figure-gfm/plot-pie-1.png)

We can also scale by additional columns. Here, we’ll find the proportion
of each window which overlaps each type of region, providing a complete
summary of the proportion of our ranges which overlaps each feature.

``` r
regions %>% 
  lapply(function(x) propOverlap(results_gr, x)) %>% 
  as_tibble() %>% 
  mutate(range = as.character(results_gr)) %>% 
  pivot_longer(cols = names(regions), names_to = "region", values_to = "p") %>% 
  dplyr::filter(p > 0) %>% 
  mutate(
    region = fct(region_levels[region], levels = region_levels),
    w = p * width(GRanges(range)) / 1e3
  ) %>% 
  plotPie(
    scale_by = "w", fill = "region", min_p = 0.05,
    total_glue = "{round(N, 1)}kb", total_size = 5,
    cat_glue = "{str_wrap(region, 10)}\n{percent(p, 0.1)}",
    cat_alpha = 0.5, cat_size = 4
  )
```

![*Pie chart with regions scaled by width to show the proportions of the
total regions which overlap each type of annotated
region.*](differential_signal_files/figure-gfm/plot-pie-width-1.png)

These results can be extended further using `plotSplitDonut()` to show
more complex results.

``` r
results_gr %>% 
  plotSplitDonut(inner = "status", outer = "bestRegion")
```

Again, this plot is heavily customisable, and is able to utilise
separate palettes for the inner and outer rings if preferred. Specific
slices can also be *exploded* for emphasis.

``` r
region_col <- hcl.colors(length(region_levels), "Viridis", rev = TRUE)
results_gr %>% 
  subset(status != "Unchanged") %>% 
  plotSplitDonut(
    inner = "status", outer = "bestRegion", min_p = 0.01,
    inner_palette = c("#3333E6", "#E6334D"),
    outer_palette = region_col,
    inner_glue = "H3K27ac\n{status}\n{n}", inner_label_alpha = 0.8,
    outer_glue = "{str_wrap(bestRegion, 10)}\n{n}", outer_label = "text",
    explode_outer = "Promoter", explode_r = 0.2
  ) 
```

![*Summary of regions with changed signal and the annotated region with
the largest
overlap*](differential_signal_files/figure-gfm/plot-split-donut-1.png)

## Profile Heatmaps

A very common approach to visualising the results of altered TF binding
is to plot *profile heatmaps* centred around the window (or peak), and
extending out a given number of of bases. The data required for this is
referred to in `extraChIPs` as profile data, and given that extracting
this from a set of `BigWigFile`s can be time consuming, this step is
performed prior to the actual plotting, so that ranges can be added or
excluded as desired.

First we need to define a `BigWigFileList` as these are conventionally
very large files which shouldn’t be retained in memory, but are just
accessed to import the key regions for a particular process. Typically,
these can be generated during peak calling but they can also be created
directly from `bam` files if needed. Once again, our reduced dataset
makes this trivial, but the process may require significant computation
resources if being performed on a complete dataset. In the following,
we’ll create the BigWig files an immediately create our
`BigWigFileList`.

``` r
sbp <- ScanBamParam(which = GRanges("chr10:40000000-100000000"))
## Sum the coverage across all Vehicle samples then create a merged BigWig
cov_veh <- names(bfl) %>% 
  str_subset("Veh") %>% 
  lapply(function(x) coverage(bfl[[x]], param = sbp)$chr10) %>% 
  RleList() %>% 
  unlist()
cov_veh <- GRanges(RleList(chr10 = cov_veh))
export(cov_veh, "data/Veh_merged.bw")
## Repeat for DHT
cov_dht <- names(bfl) %>% 
  str_subset("DHT") %>% 
    lapply(function(x) coverage(bfl[[x]], param = sbp)$chr10) %>% 
  RleList() %>% 
  unlist()
cov_dht <- GRanges(RleList(chr10 = cov_dht))
export(cov_dht, "data/DHT_merged.bw")
bwfl <- BigWigFileList(
  list(Veh = "data/Veh_merged.bw", DHT = "data/DHT_merged.bw")
)
```

Now we have our `BigWigFileList` we can define the ranges to plot along
with the profile data. In order to find the region within each range
with the maximal signal, we’ll compare back to the original results
before windows were merged, then return the original window with the
maximal signal (i.e. logCPM). Note that we already had the windows in
“keyval_range” which corresponded to the ranges with the most
statistically significant change. Realistically we can choose any
strategy for centring ranges that we wish.

``` r
gr_increase <- subset(results_gr, status == "Increased")
max_peaks <- gr_increase %>%
  granges() %>%
  join_overlap_left(
    mutate(fit_gr, window = as.character(fit_gr))
  ) %>% 
    arrange(1/logCPM) %>% 
    distinctMC(.keep_all = TRUE) %>% 
    colToRanges("window") %>% 
    granges()
```

Now we’ve found the original windows with the strongest signal, we can
define the profile data using the surrounding regions. Our widest region
with significant change in signal is 5200bp, so let’s set our profiles
to be created stretching 10kb either side of the peak centre.

``` r
pd <- getProfileData(bwfl, max_peaks, upstream = 1e4, log = FALSE)
pd
```

    ## GRangesList object of length 2:
    ## $Veh
    ## GRanges object with 33 ranges and 1 metadata column:
    ##        seqnames            ranges strand |
    ##           <Rle>         <IRanges>  <Rle> |
    ##    [1]    chr10 79246860-79266859      * |
    ##    [2]    chr10 88335820-88355819      * |
    ##    [3]    chr10 79297340-79317339      * |
    ##    [4]    chr10 76771460-76791459      * |
    ##    [5]    chr10 72997140-73017139      * |
    ##    ...      ...               ...    ... .
    ##   [29]    chr10 67661660-67681659      * |
    ##   [30]    chr10 63687060-63707059      * |
    ##   [31]    chr10 79255940-79275939      * |
    ##   [32]    chr10 44966940-44986939      * |
    ##   [33]    chr10 79254900-79274899      * |
    ##                                                 profile_data
    ##                                         <SplitDataFrameList>
    ##    [1]       1.1206:u1:-9900,0.7100:u2:-9700,1.3750:u3:-9500
    ##    [2] 0.879397:u1:-9900,0.665000:u2:-9700,0.075000:u3:-9500
    ##    [3]    1.23618:u1:-9900,0.99000:u2:-9700,0.37000:u3:-9500
    ##    [4] 0.266332:u1:-9900,0.100000:u2:-9700,0.030000:u3:-9500
    ##    [5]          0.375:u1:-9900,0.000:u2:-9700,0.365:u3:-9500
    ##    ...                                                   ...
    ##   [29]          0.000:u1:-9900,0.290:u2:-9700,0.075:u3:-9500
    ##   [30] 0.296482:u1:-9900,0.135000:u2:-9700,0.235000:u3:-9500
    ##   [31] 19.89500:u1:-9900,11.28500:u2:-9700, 5.71642:u3:-9500
    ##   [32] 0.738693:u1:-9900,1.260000:u2:-9700,2.790000:u3:-9500
    ##   [33]        8.590:u1:-9900,13.750:u2:-9700, 9.565:u3:-9500
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
    ## 
    ## $DHT
    ## GRanges object with 33 ranges and 1 metadata column:
    ##        seqnames            ranges strand |
    ##           <Rle>         <IRanges>  <Rle> |
    ##    [1]    chr10 79246860-79266859      * |
    ##    [2]    chr10 88335820-88355819      * |
    ##    [3]    chr10 79297340-79317339      * |
    ##    [4]    chr10 76771460-76791459      * |
    ##    [5]    chr10 72997140-73017139      * |
    ##    ...      ...               ...    ... .
    ##   [29]    chr10 67661660-67681659      * |
    ##   [30]    chr10 63687060-63707059      * |
    ##   [31]    chr10 79255940-79275939      * |
    ##   [32]    chr10 44966940-44986939      * |
    ##   [33]    chr10 79254900-79274899      * |
    ##                                                 profile_data
    ##                                         <SplitDataFrameList>
    ##    [1]    5.36181:u1:-9900,2.61500:u2:-9700,2.20000:u3:-9500
    ##    [2]             0.00:u1:-9900,0.95:u2:-9700,0.50:u3:-9500
    ##    [3] 0.708543:u1:-9900,1.500000:u2:-9700,0.740000:u3:-9500
    ##    [4]             0.00:u1:-9900,0.28:u2:-9700,0.09:u3:-9500
    ##    [5]             0.00:u1:-9900,1.05:u2:-9700,0.04:u3:-9500
    ##    ...                                                   ...
    ##   [29]             0.37:u1:-9900,0.00:u2:-9700,0.03:u3:-9500
    ##   [30]             0.00:u1:-9900,0.74:u2:-9700,0.37:u3:-9500
    ##   [31]       28.210:u1:-9900,16.725:u2:-9700, 7.270:u3:-9500
    ##   [32]    2.88945:u1:-9900,0.56000:u2:-9700,1.73500:u3:-9500
    ##   [33]       10.950:u1:-9900,15.230:u2:-9700,16.265:u3:-9500
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

This produces a `GRangesList` with a `GRanges` element for every file in
the `BigWigFileList`, which has the profile data stored in the final
column. Each element of these columns is a `DataFrame` with the region
broken into a defined number of bins, and an average coverage value
calculated. We can then simply plot this data by specifying this column
in the function `plotProfileHeatmap()`, which produces a `ggplot2`
object able to be customised in the conventional manner. Here, we’ll add
a colour scale and `theme_bw()`

``` r
plotProfileHeatmap(pd, "profile_data") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(fill = "Counts") +
  theme_bw()
```

![*Profile Heatmap of the regions on chromosome 10 showing evidence for
increased H3K27ac signal. Regions are centred at the point of maximal
signal.*](differential_signal_files/figure-gfm/profile-heatmap-1.png)

Note that H3K27ac often contains regions with significantly lower signal
where the transcriptional machinery is instead occupying the DNA and as
such, H3K27ac peaks are often centred using peaks from an additional
transcription factor. Here, we’ve also just plotted raw counts. `macs2`
is able to produce BigWigFiles containing logCPM values, or enrichment
over input and these are often used in preference to raw counts.

## Inspection of Ranges

Another important step in the analysis of ChIP-Seq data is to look at
the binding patterns using coverage, and inspect these in reference to
genes and any other feature of interest. The function `plotHFGC()`
provides a simple, standardised layout using the visualisation tools
from `Gviz`. If supplied, tracks will be drawn in the order 1) HiC; 2)
Features; 3) Genes, and 4) Coverage. All tracks are optional, and if
none are provided a simply cytogenetic plot will be produced. Whilst a
simple and intuitive function to use, it also provides a great deal of
flexibility for advanced customisation. All plots require a `GRanges`
object to define the plotting region along with a set of cytogenetic
bands. These are provided in the package for GRCh37 and GRCh38.

Let’s start by plotting the first promoter showing changed signal in
`results_gr` using the minimal data possible. This is a `GRanges` object
and some cytogenetic bands.

``` r
data("grch37.cytobands")
gr <- results_gr %>%
  subset(hmp_fdr < 0.05) %>% 
  subset(vapply(gene_name, function(x) "PRXL2A" %in% x, logical(1))) %>% 
  subset(logCPM == max(logCPM))
plotHFGC(gr, cytobands = grch37.cytobands)
```

![*Basic output from plotHFGC without any of the optional
tracks*](differential_signal_files/figure-gfm/plot-empty-hfgc-1.png)

### Including Coverage

In order to show our changed signal in context we can show the coverage
using a `BigWigFileList` as for the Profile Heatmap. If providing a
`BigWigFileList`, separate tracks will drawn for each element of the
object. In more complex scenarios where multiple targets are required a
list of `BigWigFileList`s can be provided and each element will then be
drawn as a track with overlapping coverage within each BigWigFileList.
Colours need to be provided in the identical structure to match the
provided object. Here, we’ve provided a single `BigWigFileList` so we
can provide a simple list matching the names within the
`BigWigFileList.`

``` r
plotHFGC(
  gr, cytobands = grch37.cytobands, 
  coverage = bwfl, linecol = list("Veh" = "grey30", "DHT" = "darkred"),
  zoom = 30
)
```

![*Coverage for the region of interest, shown with a blue highlight.
Note how the counts are lower in the highlighted region for the merged
Veh coverage
track.*](differential_signal_files/figure-gfm/plot-hfgc-coverge-1.png)

### Displaying Genes

Next we might like to add gene models to provide the regulatory context.
These are supplied here in the layout required by the defaults of the
`GeneRegionTrack()` function, with all exons and transcripts annotated.

``` r
gene_models <- read_rds("data/chr10_trans_subset.rds")
plotHFGC(
  gr, cytobands = grch37.cytobands, 
  coverage = bwfl, linecol = list("Veh" = "grey30", "DHT" = "darkred"),
  genes = gene_models, genecol = "wheat",
  zoom = 30
)
```

![*Coverage for our region showing the relationship of signal to
annotated genes*](differential_signal_files/figure-gfm/add-genes-1.png)

### Adding Features

Another useful track to add might be some key features such as promoters
and other annotated regions. Features must **always** be a
`GRangesList`, with each element defining a different type of feature,
as we already have in our `regions_gr` object.

``` r
region_grl <- splitAsList(regions_gr, regions_gr$region)[region_levels]
names(region_col) <- region_levels
plotHFGC(
  gr, cytobands = grch37.cytobands, 
  features = region_grl, featcol = region_col, featstack = "dense",
  coverage = bwfl, linecol = list("Veh" = "grey30", "DHT" = "darkred"),
  genes = gene_models, gene_col = "wheat",
  zoom = 30
)
```

![*The same figure as previously, but with annotated regions added as
features. Any type of feature can be added
here.*](differential_signal_files/figure-gfm/plot-hfgc-1.png)

### Adding HiC Interactions

If long-range interactions are available, these can also be provided as
a GenomicInteractions object, completing all available options for the
HFGC components.

### Adding Annotations To Coverage

An indication of which regions are associated with increased or
decreased ChIP signal can also be a useful annotation to add to plots
such as the above. Although we technically performed no statistical
testing, let’s consider a window with logFC below -1 to be showing
decreased signal.

Similar to the features track, where the names of `GRangesList` elements
denote the different feature types, able to then assigned a colour,
coverage annotation tracks follow these same rules. For each coverage
track being annotated, a `GRangesList` object can denote the ranges
which can be assigned different colours.

``` r
cov_annot <- splitAsList(results_gr, results_gr$status) %>% 
  endoapply(granges)
```

In the above, we have Unchanged and Decreased signal denoted as
annotations. In keeping with the approach of having a matching list
element for every coverage track, we would need to pass this as a list
which matched the coverage track

``` r
plotHFGC(
  gr, 
  features = region_grl, featcol = region_col, featstack = "dense",
  genes = gene_models, genecol = "wheat",
  coverage = bwfl, linecol = list("Veh" = "grey30", "DHT" = "darkred"),
  annotation = cov_annot, 
  annotcol = c(Unchanged = "grey", Decreased = "#3333E6", Increased = "#E6334D"),
  cytobands = grch37.cytobands, zoom = 30
)
```

![*The addition of an annotation track for the coverage tracks shows
which regions were retained during the analysis, as well as those which
were considered as showing changed or unchanged
signal.*](differential_signal_files/figure-gfm/plot-annot-1.png)

Plots are able to be tweaked considerably further via multiple
parameters, however these basic approaches cover the core functionality
of `plotHFCG()` for enabling simple & reproducible plotting across
regions for multiple sites within a larger experiment.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-edgeR2016" class="csl-entry">

Chen, Yunshun, Aaron A T Lun, and Gordon K Smyth. 2016. “From Reads to
Genes to Pathways: Differential Expression Analysis of RNA-Seq
Experiments Using Rsubread and the edgeR Quasi-Likelihood Pipeline.”
*F1000Research* 5: 1438.
<https://doi.org/10.12688/f1000research.8987.2>.

</div>

<div id="ref-Gandolfo2018-oc" class="csl-entry">

Gandolfo, Luke C, and Terence P Speed. 2018. “RLE Plots: Visualizing
Unwanted Variation in High Dimensional Data.” *PLoS One* 13 (2):
e0191629.

</div>

<div id="ref-Law2014-xq" class="csl-entry">

Law, Charity W, Yunshun Chen, Wei Shi, and Gordon K Smyth. 2014. “Voom:
Precision Weights Unlock Linear Model Analysis Tools for <span
class="nocase">RNA-seq</span> Read Counts.” *Genome Biol.* 15 (2): R29.

</div>

<div id="ref-DESeq22014" class="csl-entry">

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15: 550. <https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-csaw2016" class="csl-entry">

Lun, Aaron T L, and Gordon K Smyth. 2016. “Csaw: A Bioconductor Package
for Differential Binding Analysis of ChIP-Seq Data Using Sliding
Windows.” *Nucleic Acids Res.* 44 (5): e45.

</div>

<div id="ref-Lund2012-xo" class="csl-entry">

Lund, Steven P, Dan Nettleton, Davis J McCarthy, and Gordon K Smyth.
2012. “Detecting Differential Expression in <span
class="nocase">RNA-sequence</span> Data Using Quasi-Likelihood with
Shrunken Dispersion Estimates.” *Stat. Appl. Genet. Mol. Biol.* 11 (5).

</div>

<div id="ref-McCarthy2009-qf" class="csl-entry">

McCarthy, Davis J, and Gordon K Smyth. 2009. “Testing Significance
Relative to a Fold-Change Threshold Is a TREAT.” *Bioinformatics* 25
(6): 765–71.

</div>

<div id="ref-DiffBind2012" class="csl-entry">

Ross-Innes, Caryn S., Rory Stark, Andrew E. Teschendorff, Kelly A.
Holmes, H. Raza Ali, Mark J. Dunning, Gordon D. Brown, et al. 2012.
“Differential Oestrogen Receptor Binding Is Associated with Clinical
Outcome in Breast Cancer.” *Nature* 481: –4.
<http://www.nature.com/nature/journal/v481/n7381/full/nature10730.html>.

</div>

<div id="ref-Wilson2019-ln" class="csl-entry">

Wilson, Daniel J. 2019. “The Harmonic Mean *p*-Value for Combining
Dependent Tests.” *Proc. Natl. Acad. Sci. U. S. A.* 116 (4): 1195–1200.

</div>

<div id="ref-Zhang2008-ms" class="csl-entry">

Zhang, Yong, Tao Liu, Clifford A Meyer, Jérôme Eeckhoute, David S
Johnson, Bradley E Bernstein, Chad Nusbaum, et al. 2008. “Model-Based
Analysis of ChIP-Seq (MACS).” *Genome Biol.* 9 (9): R137.

</div>

</div>

<br>

# Session Info

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8    
    ##  [5] LC_MONETARY=en_AU.UTF-8    LC_MESSAGES=en_AU.UTF-8   
    ##  [7] LC_PAPER=en_AU.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggside_0.2.2                scales_1.2.1               
    ##  [3] plyranges_1.18.0            extraChIPs_1.3.10          
    ##  [5] patchwork_1.1.2             edgeR_3.40.2               
    ##  [7] limma_3.54.1                rtracklayer_1.58.0         
    ##  [9] BiocParallel_1.32.5         csaw_1.32.0                
    ## [11] SummarizedExperiment_1.28.0 Biobase_2.58.0             
    ## [13] MatrixGenerics_1.10.0       matrixStats_0.63.0         
    ## [15] Rsamtools_2.14.0            Biostrings_2.66.0          
    ## [17] XVector_0.38.0              GenomicRanges_1.50.2       
    ## [19] GenomeInfoDb_1.34.9         IRanges_2.32.0             
    ## [21] S4Vectors_0.36.1            BiocGenerics_0.44.0        
    ## [23] lubridate_1.9.2             forcats_1.0.0              
    ## [25] stringr_1.5.0               dplyr_1.1.0                
    ## [27] purrr_1.0.1                 readr_2.1.4                
    ## [29] tidyr_1.3.0                 tibble_3.1.8               
    ## [31] ggplot2_3.4.1               tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.4.1            circlize_0.4.15           
    ##   [3] Hmisc_4.8-0                BiocFileCache_2.6.1       
    ##   [5] igraph_1.4.1               lazyeval_0.2.2            
    ##   [7] splines_4.2.3              digest_0.6.31             
    ##   [9] ensembldb_2.22.0           foreach_1.5.2             
    ##  [11] htmltools_0.5.4            fansi_1.0.4               
    ##  [13] magrittr_2.0.3             checkmate_2.1.0           
    ##  [15] EnrichedHeatmap_1.27.2     memoise_2.0.1             
    ##  [17] BSgenome_1.66.3            cluster_2.1.4             
    ##  [19] doParallel_1.0.17          InteractionSet_1.26.1     
    ##  [21] tzdb_0.3.0                 ComplexHeatmap_2.14.0     
    ##  [23] timechange_0.2.0           prettyunits_1.1.1         
    ##  [25] jpeg_0.1-10                colorspace_2.1-0          
    ##  [27] ggrepel_0.9.3              blob_1.2.3                
    ##  [29] rappdirs_0.3.3             xfun_0.37                 
    ##  [31] crayon_1.5.2               RCurl_1.98-1.10           
    ##  [33] VariantAnnotation_1.44.1   survival_3.5-5            
    ##  [35] iterators_1.0.14           glue_1.6.2                
    ##  [37] polyclip_1.10-4            gtable_0.3.1              
    ##  [39] zlibbioc_1.44.0            GetoptLong_1.0.5          
    ##  [41] DelayedArray_0.24.0        shape_1.4.6               
    ##  [43] futile.options_1.0.1       DBI_1.1.3                 
    ##  [45] Rcpp_1.0.10                viridisLite_0.4.1         
    ##  [47] progress_1.2.2             htmlTable_2.4.1           
    ##  [49] clue_0.3-64                foreign_0.8-84            
    ##  [51] bit_4.0.5                  Formula_1.2-5             
    ##  [53] metapod_1.6.0              htmlwidgets_1.6.1         
    ##  [55] httr_1.4.5                 RColorBrewer_1.1-3        
    ##  [57] ellipsis_0.3.2             farver_2.1.1              
    ##  [59] pkgconfig_2.0.3            XML_3.99-0.13             
    ##  [61] Gviz_1.42.1                nnet_7.3-18               
    ##  [63] dbplyr_2.3.1               deldir_1.0-6              
    ##  [65] locfit_1.5-9.7             utf8_1.2.3                
    ##  [67] labeling_0.4.2             tidyselect_1.2.0          
    ##  [69] rlang_1.0.6                AnnotationDbi_1.60.0      
    ##  [71] munsell_0.5.0              tools_4.2.3               
    ##  [73] cachem_1.0.7               cli_3.6.0                 
    ##  [75] generics_0.1.3             RSQLite_2.3.0             
    ##  [77] broom_1.0.3                evaluate_0.20             
    ##  [79] fastmap_1.1.1              yaml_2.3.7                
    ##  [81] knitr_1.42                 bit64_4.0.5               
    ##  [83] AnnotationFilter_1.22.0    KEGGREST_1.38.0           
    ##  [85] formatR_1.14               xml2_1.3.3                
    ##  [87] biomaRt_2.54.0             compiler_4.2.3            
    ##  [89] rstudioapi_0.14            filelock_1.0.2            
    ##  [91] curl_5.0.0                 png_0.1-8                 
    ##  [93] tweenr_2.0.2               stringi_1.7.12            
    ##  [95] highr_0.10                 futile.logger_1.4.3       
    ##  [97] GenomicFeatures_1.50.4     lattice_0.20-45           
    ##  [99] ProtGenerics_1.30.0        Matrix_1.5-1              
    ## [101] vctrs_0.5.2                pillar_1.8.1              
    ## [103] GenomicInteractions_1.32.0 lifecycle_1.0.3           
    ## [105] BiocManager_1.30.20        GlobalOptions_0.1.2       
    ## [107] ComplexUpset_1.3.3         data.table_1.14.8         
    ## [109] bitops_1.0-7               R6_2.5.1                  
    ## [111] BiocIO_1.8.0               latticeExtra_0.6-30       
    ## [113] gridExtra_2.3              codetools_0.2-19          
    ## [115] lambda.r_1.2.4             dichromat_2.0-0.1         
    ## [117] MASS_7.3-58.2              rjson_0.2.21              
    ## [119] withr_2.5.0                GenomicAlignments_1.34.0  
    ## [121] GenomeInfoDbData_1.2.9     parallel_4.2.3            
    ## [123] hms_1.1.2                  VennDiagram_1.7.3         
    ## [125] grid_4.2.3                 rpart_4.1.19              
    ## [127] rmarkdown_2.20             biovizBase_1.46.0         
    ## [129] ggforce_0.4.1              base64enc_0.1-3           
    ## [131] interp_1.1-3               restfulr_0.0.15
