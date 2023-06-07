# Data Preparation for Vignette

The range used for this dataset is chr10:42354900-100000000.
All bam files were subset and indexed using the following example command.

```bash
 samtools view -h -b data/deduplicated/SRR8315180.sorted.bam "chr10:42354900-100000000" > data/ER/SRR8315180.bam
 samtools index data/ER/SRR8315180.bam
```

Peak files were loaded in R using extraChIPs and re-written after subsetting

```r
library(rtracklayer)
library(extraChIPs)
peakFiles <- list.files("data/ER", pattern = "narrowPeak", full.names = TRUE)
peaks <- importPeaks(peakFiles, seqinfo = sq)
subset_peaks <- endoapply(peaks, subsetByOverlaps, GRanges("chr10:42354900-100000000"))
names(subset_peaks) %>% 
  lapply(
    function(x) {
      df <- subset_peaks[[x]] %>% 
        as.data.frame %>% 
        rownames_to_column("name") %>% 
        dplyr::select(
          seqnames, start, end, name, score, strand, ends_with("Value"), peak
        ) %>%
        mutate(start = start - 1, strand = ".")
      fl <- here::here("data/ER", x)
      write.table(
        df, fl, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
      )
    }
  )
```

Coverage BigWig files were also processed in a similar manner

```r
fl <- list.files("output/macs2/ERa", pattern = "pileup.bw", full.names = TRUE)
gr <- GRanges("chr10:42354900-100000000")
bwfl <- BigWigFileList(fl)
names(bwfl) <- c("E2", "E2DHT")
lapply(
  names(bwfl), 
  function(x) {
    cov <- import.bw(bwfl[[x]], which = gr)
    out <- file.path(
      "~/DRMCRL/extraChIPs_vignette/data/ER", paste0(x, "_cov_chr10.bw")
    )
    export.bw(cov, out)
  }
)
fl <- list.files("output/macs2/ERa", pattern = "FE.bw", full.names = TRUE)
bwfl <- BigWigFileList(fl)
names(bwfl) <- c("E2", "E2DHT")
lapply(
  names(bwfl), 
  function(x) {
    cov <- import.bw(bwfl[[x]], which = gr)
    out <- file.path(
      "~/DRMCRL/extraChIPs_vignette/data/ER", paste0(x, "_FE_chr10.bw")
    )
    export.bw(cov, out)
  }
)
```
