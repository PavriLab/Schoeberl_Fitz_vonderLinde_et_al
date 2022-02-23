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

### Modified CCseqBasic5.sh script
In order to facilitate filtering of reads in terms of mapping quality to remove reads with a MAPQ score below 30, we had to modify the CCseqBasic script because it does not contain this option by default. The modification is simple and is just done by replacing the SAM to BAM conversion directly after alignment with bowtie with a SAM to BAM conversion which includes a MAPQ filtering. This step is facilitated by the `applyMapqFilter.py` script and simply resets all reads with a MAPQ < 20 or 30 to unaligned. This way the statistics are untouched and can be computed and plotted as usual. Low mapping reads are then part of the unaligned fraction.

Plotting statistics

The plotTriCstats.py is meant to generate a comprehensive view on the statistics the CCseq pipeline produces. It requires the report files of the flashed and non-flashed reads of stage 3 (located in the F3 folder) and the the combined stats of stage 6 (located in the F6 folder) as well as the total number of reads entering the pipeline which are read from the trimming log file (located in the F1 folder) and the sample labels. It suffices to specify the folder(s) where these report folders are located. E.g. for a sequencing run containing 9 samples a valid call of the script would be

python3 plotTriCstats.py -d TriC_7_Emu*
                         -o plots/TriC_7_Emu_stats.pdf

where the wildcard is used to just address all existing folders that are supposed to be analysed.
If these statistics should also be saved in a csv file, the `--dataOut` argument can be specified.
These csv files can be used by the `plotTriCstatsPool.py` script to plot the average statistics across mutiple replicates.

### Summing single capture interaction files
If you are using multiple probes in one sample the pipeline will produce an interaction file for each probe separately. This is usually not what we want and we thus need to recombine them to a single interaction file. This is done with the `sumInteractionFiles.py`, which just sums the counts for the individual restriction fragments and writes the result including the maximum count to a new file. An example command looks like follows:
```bash
python3 sumInteractionFiles.py -i CCseq/${sample}/F6_greenGraphs_combined_sample_CS5/sample_TriC/sample_*_TriC_interactions.txt \
                               -o interactions/${sample}_TriC_interactions.txt
```

### Visualizing the results
The final step of the analysis is the visualization of the results of the [CCseq](https://github.com/Hughes-Genome-Group/CCseqBasicS) pipeline. This is done by the `TriCplot.py` script. However, in order to do this we need to generate a few additional files from the pipeline output, namely fetching the COMBINED.bam file of stage 6 and extracting reads that have only 3 or more than 3 fragments, which is necessary to generate the profile plots. The `TriCgetReads.py` script does exactly that. Example commands for this would look as follows:
```bash
# for reads with 3 or more fragments
python3 TriCgetReads.py -b CCseq/${sample}/F6_greenGraphs_combined_sample_CS5/COMBINED_reported_capture_reads_CS5.bam \
                        -n 3 \
                        --larger \
                        -o morethanthree/${sample}_3wayplus.bam
```
The `-n` argument specifies the number of fragments a valid read has. This is an exact process. If one would like to get all reads with exactly n or more fragments, the `--larger` flag has to be set.

An **optional** next step is the generation of the read profile over our region of interest for each sample, which is done using the [deeptools](https://deeptools.readthedocs.io/en/develop/) `multiBamSummary`. This requires the BAM files to be position sorted and indexed, which can easily be done with [`samtools`](http://www.htslib.org/doc/samtools.html). The `multiBamSummary` command to use would look like this:
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
The results we need is the `profiles/${sample}_3plus.tsv` file. Note that the used binsize has to be the same as the binsize used for calculating the Tri-C matrices using the `TriC_matrix_simple_MO.py` script either of the original publication from the repository of [Marike Oudelaar](https://github.com/oudelaar/TriC) or from this repository (please see section below for more information about the differences).  
If these files are not generated, the profile can be derived directly from the contact matrices during plotting with `TriCplot.py` by adding the argument `--derivedProfile`.

The next step is to collect the Tri-C matrix files, as computed by the CCseq pipeline with the `--tric` flag set from their residing folder (F7*/sample/), which is based on the actual restriction fragment boundaries, and using the above mentioned script to generate a matrix file with a fixed binsize with
```bash
python3 TriC/TriC_matrix_simple_MO.py -f CCseqmatrixfile -c chr12 -l 114435000 -r 114669000 -b 1000 -t 80 -a -o folder
```
where `-t` specifies a cutoff for the value each matrix entry can acquire. The generated files can then be used with the `TriCplot.py` script to generate data visualizations of the pipelines results. An example command is given below:
```bash
python3 TriCplot.py --treatment data/matrices/priB_d2_Emu_*_TriC_interactions_1000_RAW.tab \
                    --treatment_label priB_d2 \
                    --control data/matrices/priB_d0_Emu_*_TriC_interactions_1000_RAW.tab \
                    --control_label priB_d0 \
                    --region "chr12:114435000-114669000" \
                    --binsize 1000 \
            		    --compare_vMax 100 \
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
A more detailed documentation on the different commandline arguments can be found via `python3 TriCplot.py --help/-h`. Results and test data can be found in the data folder. The `data/annotations/vdj_genes.sort.bed` file is a custom bed file with annotations of the region. It is processed by latech, so the display_names can show special characters like greek letters.

### Changes to `TriC_matrix_simple_MO.py`
This script is basically a copy of the original script of [Marieke Oudelaar](https://github.com/oudelaar/TriC), with the addition of the ability to specify a second region which will be used to compute the matrix normalization factor and a new general normalization. Matrices are normalized to 300,000 total interactions. The new parameters are `--normchrom`, `--nstr` and `--nstp` and a typical command would look like this
```bash
python3 TriC/TriC_matrix_simple_MO.py -f CCseqmatrixfile \
                                      -c chr12 -l 114435000 -r 114669000 \
                                      --normchrom chr12 --nstr 114435000 --nstp 114496000 \
                                      -b 1000 \
                                      -t 50 \
                                      -a \
                                      -o TriCplots
```

### Running Environment
All packages that were used can be found in the `TriCenv.yml` conda environment. In addition Latech was installed as techlive 2020. This is necessary for the plotting. 

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
