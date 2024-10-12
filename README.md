# The LIF-STAT3 Pathway Regulates H3K27me3 Desposition

This repostiory contains the codes and shell commands that used in our manuscript ***The LIF-STAT3 Pathway Regulates H3K27me3 Desposition***. The raw sequencing data have been submitted to GEO under the BioProject [PRJNA1145863](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1145863?reviewer=h6q355m9rgshfs06).

# Procedures

### Step 0: software tools

The following software tools were used in the manuscript:

- [fastp](https://github.com/OpenGene/fastp)
- [hisat2](https://daehwankimlab.github.io/hisat2/)
- [samtools](https://www.htslib.org)
- [sambamba](https://lomereiter.github.io/sambamba/)
- [macs2](https://pypi.org/project/MACS2/)
- [bedtools](https://bedtools.readthedocs.io)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [deepTools](https://deeptools.readthedocs.io)
- [matplotlib-venn](https://pypi.org/project/matplotlib-venn/)

### Step 1: Adapter Trimming and QC

We used `fastp` to trim off the adapter sequences at the 3' end of the read. Some QC information was also returned by the program:

```bash
for i in H3K27me3_Ctrl_rep1 H3K27me3_Ctrl_rep2 H3K27me3_noLIF_rep1 H3K27me3_noLIF_rep2; do
    fastp -l 25 -w 10 -q 20 --detect_adapter_for_pe \
        -i ${i}_R1.fastq.gz -I ${i}_R2.fastq.gz \
        -o 01_fastp/${i}_r1.trim.fastq.gz -O 01_fastp/${i}_r2.trim.fastq.gz \
        -h 01_fastp/${i}.fastp.html -j 01_fastp/${i}.fastp.json
done
```

### Step 2: Read alignment

We used `hisat2` to align the trimmed reads to the mouse reference genome `mm10`. During the alignment, we removed reads that were mapped to the mitochondrial geonme (chrM) and only kept reads that were properly paired and uniquely mapped (MAPQ > 30):

```bash
for i in H3K27me3_Ctrl_rep1 H3K27me3_Ctrl_rep2 H3K27me3_noLIF_rep1 H3K27me3_noLIF_rep2; do
    hisat2 -X 500 -p 10 \
        --no-spliced-alignment --no-temp-splicesite \
        -x /mm10/genome \
        -1 01_fastp/${i}_r1.trim.fastq.gz \
        -2 01_fastp/${i}_r2.trim.fastq.gz \
        --summary-file 02_map/${i}_algn_sum.txt | \
        grep -v chrM | samtools view -h -f 2 -F 3844 -q 30 - | \
        samtools sort - -T ${i}_tmp -o 02_map/${i}.sort.bam
done
```

### Step 3: De-duplication

We used the `sambamba` tool to remove duplicates in the alignment:

```bash
# de-duplicate
for i in H3K27me3_Ctrl_rep1 H3K27me3_Ctrl_rep2 H3K27me3_noLIF_rep1 H3K27me3_noLIF_rep2; do
    sambamba markdup -t 10 -r -p 02_map/${i}.sort.bam \
        03_sambamba/${i}_sambamba.bam
done

# sort
for i in H3K27me3_Ctrl_rep1 H3K27me3_Ctrl_rep2 H3K27me3_noLIF_rep1 H3K27me3_noLIF_rep2; do
    samtools sort -@ 10 03_sambamba/${i}_sambamba.bam \
        -o 03_sambamba/${i%.bam}_sort.bam
done
```

### Step 4: Peak calling

We used `MACS2` for the identification of peak regions. Once the process was done, we also used the UCSC utility [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) to convert the `bedGraph` file to `bigWig` file for the visualisation in the UCSC genome browser:

```bash
# peak calling
for i in H3K27me3_Ctrl_rep1 H3K27me3_Ctrl_rep2 H3K27me3_noLIF_rep1 H3K27me3_noLIF_rep2; do
    macs2 callpeak -t 03_sambamba/${i}_sambamba_sort.bam \
        -c 03_sambamba/${i/H3K27me3/Input}_sambamba_sort.bam \
        -q 0.05 -f BAMPE -B --SPMR \
        -g mm --keep-dup all \
        -n ${i} --outdir 04_callpeak
done

# covnert bdg to bw
for i in 04_callpeak/*_treat_pileup.bdg; do
    sort -k1,1 -k2,2n ${i} > ${i}.sort.bdg
    bedGraphToBigWig ${i}.sort.bdg mm10.chrom.sizes ${i}.sort.bw
done
```

### Step 5: Differential peak identification

We first found "high confidence" peaks for each condition by taking peaks that are present in both biological repilcates. Then we make a union set of the high confidence peaks between control and "no LIF" samples. Then the number of reads from each sample in each peak were counted. Finally, we collected all the counts into an matrix and feed that into `DESeq2`:

```bash
# find high confidence peak
intersectBed -a 04_callpeak/H3K27me3_Ctrl_rep1_peaks.broadPeak \
             -b 04_callpeak/H3K27me3_Ctrl_rep2_peaks.broadPeak \
             -u | cut -f 1-4 > 04_callpeak/H3K27me3_Ctrl_highCon_peaks.bed
intersectBed -a 04_callpeak/H3K27me3_noLIF_rep1_peaks.broadPeak \
             -b 04_callpeak/H3K27me3_noLIF_rep2_peaks.broadPeak \
             -u | cut -f 1-4 > 04_callpeak/H3K27me3_noLIF_highCon_peaks.bed

# generate the union peak set
cat 04_callpeak/H3K27me3_Ctrl_highCon_peaks.bed \
    04_callpeak/H3K27me3_noLIF_highCon_peaks.bed | \
    sort -k1,1 -k2,2n | \
    mergeBed -i - | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' \
    > 04_callpeak/H3K27me3_union_peak.bed

# get counts from each sample and replicate
for i in H3K27me3_Ctrl_rep1 H3K27me3_Ctrl_rep2 H3K27me3_noLIF_rep1 H3K27me3_noLIF_rep2; do
    coverageBed -a 04_callpeak/H3K27me3_union_peak.bed \
        -b 03_sambamba/${i}_sambamba_sort.bam -sorted \
        > 05_coverage/${i}_count.tsv
done

# differential analysis by Deseq2 in R
Rscript run_deseq2.r
```

### Commands to generate figures

To generate the heatmap in **Figures 2C** and **S2C**, the following commands from `deepTools` were used:

```bash
computeMatrix reference-point --referencePoint center \
    -p 10 -b 8000 -a 8000 -R H3K27me3_peaks_up.bed \
    -S LIF-H3K27me3_treat_pileup.bdg.sort.bw Con-H3K27me3_treat_pileup.bdg.sort.bw \
    --skipZeros -o H3k27me3_up.gz
computeMatrix reference-point --referencePoint center \
    -p 10 -b 8000 -a 8000 -R H3K27me3_peaks_down.bed \
    -S LIF-H3K27me3_treat_pileup.bdg.sort.bw Con-H3K27me3_treat_pileup.bdg.sort.bw \
    --skipZeros -o H3k27me3_down.gz

plotProfile -m H3k27me3_up.gz -out H3K27me3_up.pdf \
    --perGroup --colors red blue --plotTitle "H3K27me3_up" \
    --plotHeight 8 --plotWidth 8 --samplesLabel Con noLIF
plotProfile -m H3k27me3_down.gz -out H3K27me3_down.pdf \
    --perGroup --colors red blue --plotTitle "H3K27me3_down" \
    --plotHeight 8 --plotWidth 8 --samplesLabel Con noLIF

plotHeatmap -m H3k27me3_up.gz --colorList 'white,blue' \
    --boxAroundHeatmaps no -out H3K27me3_heatmap_up.pdf \
    --heatmapHeight 9 --heatmapWidth 3 --whatToShow 'heatmap only'
plotHeatmap -m H3k27me3_down.gz --colorList 'white,blue' \
    --boxAroundHeatmaps no -out H3K27me3_heatmap_down.pdf \
    --heatmapHeight 9 --heatmapWidth 3 --whatToShow 'heatmap only'
```

To generate **Figure3G**, the following `python` script was used:

```python
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
my_dpi = 100
plt.figure(figsize=(800/my_dpi, 600/my_dpi), dpi=my_dpi)
venn3(subsets = (6375, 7134, 2207, 19328,2447,979,2357), set_labels = ('Beta-catenin', 'DPF2', 'STAT3'))
plt.savefig('venn_diagram.pdf')
plt.show()
```