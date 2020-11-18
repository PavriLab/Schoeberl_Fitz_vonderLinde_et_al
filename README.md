# tric
repository for Tri-C data analysis related code.

## Overview

### Demultiplexing reads
This task is a very generic one and could in principle also be done by any other dedicated software. Nevertheless, `TriCdemultiplex.py` handles this task at least for the samples generated at the VBCF facilities. The script requires you to pass the BAM file containing the reads to demultiplex, a file containing the barcodes as well as a mapping to the sample name and the number of mismatches allowed in the barcode with the `-m` argument. The `-t` argument specifies the number of threads to use for demultiplexing and `-p` gives the prefix of the demultiplexed fastq files.
```bash
python3 utils/TriCdemultiplex.py -i bams/200109_TriC_1.bam \
                                 -b resource/barcodes.tsv \
                                 -m 1 \
                                 -t 7 \
                                 -p fastqs/TriC_1 \
                                 -s demux_stats_TriC_1.tsv
```
The `barcodes.tsv` is a tab-separated file with the following layout (without the header line)
```
barcode1  barcode2  sampleID  samplename
GACGCATC        GAGGCTGC        125065  priB_d0_3_i8
TGCCGTAG        CCTCGTAG        125066  priB_d0_3_i9
```
|Column |Description|
|---|---|
|barcode1 |first barcode stored in the BC tag of the BAM entry (sample_adaptor_tag) |
|barcode2 |second barcode stored in the B2 tag of the BAM entry (sample_adaptor_secondary_tag)  |
|sampleID |sample id as given by the facility (this is mainly for completeness) |
|samplename |name of the sample used for naming the demultiplexed files (sample_description) |

This information can retrieved via the seqmate site of the VBCF.

### Plotting statistics
The `plotTriCstats.py` is meant to generate a comprehensive view on the statistics the [CCseq](https://github.com/Hughes-Genome-Group/CCseqBasicS) pipeline produces. It requires the report files of the flashed and non-flashed reads of stage 3 (located in the F3 folder) and the the combined stats of stage 6 (located in the F6 folder) as well as the total number of reads entering the pipeline and the sample labels. E.g. for a sequencing run containing 9 samples a valid call of the script would be
```bash
python3 utils/plotTriCstats.py -f6 CCseq/TriC_7_*/F6*/COMBINED_report_CS5.txt \
                               -f3f CCseq/TriC_7_*/F3*/FLASHED_REdig_report_CS5.txt \
                               -f3n CCseq/TriC_7_*/F3*/NONFLASHED_REdig_report_CS5.txt \
                               -rn 4088719 4572065 4859024 5031120 4834634 4942035 6261844 5649615 5372841 \
                               -s mESC_1 mESC_2 mESC_3 priB_d0_1 priB_d0_2 priB_d0_3 priB_d2_1 priB_d2_2 priB_d2_3 \
                               -o plots/TriC_7_Emu_stats.pdf
```
where the wildcard is used to just address all existing folder. **The script is sensitive to the order the folders are given to it, meaning that it requires that all folders and additional information (readcount, sample name) are given in the same order as the folders. E.g. for mESC_1 the folders as well as the readcount have to be at the first place of the list passed to the script. So make sure this is the case, otherwise the results will be not interpretable.**

### Visualizing the results
The final step of the analysis is the visualization of the results of the [CCseq](https://github.com/Hughes-Genome-Group/CCseqBasicS) pipeline. This is done by the `TriCplot.py` script. However, in order to do this we need to generate a few additional files from the pipeline output, namely fetching the COMBINED.bam file of stage 6 and extracting reads that have only 3 or more than 3 fragments, which is necessary to generate the profile plots. The `TriCgetReads.py` script does exactly that. Example commands for this would look as follows:
```bash
# for reads with 3 or more fragments
python3 utils/TriCgetReads.py -b CCseq/${sample}/F6_greenGraphs_combined_sample_CS5/COMBINED_reported_capture_reads_CS5.bam -n 3 -o morethanthree/${sample}_3wayplus.bam --larger
```
The `-n` argument specifies the number of fragments a valid read has. This is an exact process. If one would like to get all reads with exactly n or more fragments, the `--larger` flag has to be set.

The next step is the generation of the read profile over our region of interest for each sample, which is done using the [deeptools](https://deeptools.readthedocs.io/en/develop/) `multiBamSummary`. This requires the BAM files to be position sorted and indexed, which can easily be done with [`samtools`](http://www.htslib.org/doc/samtools.html). The `multiBamSummary` command to use would look like this:
```bash

# 3 or more way reads
for sample in `less samples.txt`;
do
    multiBamSummary bins -p 8 \
                         -b morethanthree/${sample}.sort.bam \
                         --binSize 1000 \
                         --labels ${sample} \
                         --region chr12 \
                         -o profiles/${prefix}_3plus.npz \
                         --outRawCounts profiles/${sample}_3plus.tsv
done
```
The results we need is the `profiles/${sample}_3plus.tsv` file. Note that the used binsize has to be the same as the binsize used for calculating the Tri-C matrices using the `TriC_matrix_simple_MO.py` script of the original publication from the repository of [Marike Oudelaar](https://github.com/oudelaar/TriC). The next step is to collect the Tri-C matrix files, as computed by the CCseq pipeline with the `--tric` flag set from their residing folder (F7*/sample/), which is based on the actual restriction fragment boundaries, and using the above mentioned script to generate a matrix file with a fixed binsize with
```bash
python3 TriC/TriC_matrix_simple_MO.py -f CCseqmatrixfile -c chr12 -l 114435000 -r 114669000 -b 1000 -t 50 -a -o TriCplots
```
where `-t` specifies a cutoff for the value each matrix entry can acquire. The generated files can then be used with the `TriCplot.py` script to generate data visualizations of the pipelines results. An example command is given below:
```bash
python3 TriCplot.py --treatment data/matrices/priB_d2_Emu_*_TriC_interactions_1000_RAW.tab \
                    --treatment_label priB_d2 \
                    --control data/matrices/priB_d0_Emu_*_TriC_interactions_1000_RAW.tab \
                    --control_label priB_d0 \
                    --region "chr12:114435000-114669000" \
                    --capture_bins capture.oligo \
                    --annotation data/annotations/vdj_genes.sort.bed data/annotations/vdj_REs.bed data/annotations/mappability.bw \
                    --annotation_drawstyle Line2D Marker bigwig \
                    --annotation_yMin 0 0 0 \
                    --annotation_yMax 1 1 1 \
                    --annotation_labels Genes REs map \
                    --alternating 0 0 0 \
                    --highlight_annotation 1 \
                    --highlight_features Cg1 4 \
                    --treatment_3plus data/profiles/priB_d2_Emu_3plus.tsv \
                    --control_3plus data/profiles/priB_d0_Emu_3plus.tsv \
                    --profile_yMax 200 \
                    --outputFilePrefix data/priB_Emu
```
A more detailed documentation on the different commandline arguments can be found via `python3 TriCplot.py --help/-h`. Results and test data can be found in the data folder.

### FAQ
#### How do I get the oligos file for my capture probes?
 1. Know the coordinates of the exact mapping positions in your favorite genome build
 2. Go to the UCSC browser and open this genome build
 3. Select the NlaIII restriction enzyme map for this build and add it to the view
 4. Navigate to your probes mapping site
 5. Identify the closest upstream and downstream restriction site
 6. Take the start coordinate of your upstream site and the end coordinate of your downstream site
 7. Make a file with the following format
 ```
 probename  12  123456789 123456900 12  123456789 123456900 1 A
 ```
 The first is the probes name (be careful here because the pipeline is very picky and does only allow \[0-9\], \[a-z\], \[A-Z\] and \[\_\] in the names and will otherwise fail), 12 is the number of the chromosome and 123456789 and 123456900 are start and end positions of the restriction fragment (RF) your probe maps to. The coordinates after this initial block state the exclusion zones which can be used to exclude fragments that map very close to the probes RF (this can usually be the same coordinates a the probes RF if you want to use it a good value would be +- 1kb). The last block 1 and A are SNP positions and SNP base in the RF where 1 and A are a good default for saying we don't have any SNP. A more detailed overview can also be found [here](http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/oligofile.html)
