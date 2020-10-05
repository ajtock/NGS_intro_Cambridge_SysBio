# Introduction to next-generation sequencing data: tools and resources

This [link](https://ajtock.github.io/NGS_intro_Cambridge_SysBio/) will take you to the GitHub Pages rendering of this repository.

Following the success of dye-terminator (Sanger) sequencing in the Human Genome Project, a "next" generation of sequencing technologies became available which were faster, cheaper and produced many more reads in parallel. Despite their higher throughput, the reads were shorter and of lower quality than those produced by Sanger sequencing. One high-throughput sequencing technology became ubiquitous, Illumina sequencing by synthesis (on GAI/GAII/HiSeq/ MiSeq/NextSeq instruments). This technology was developed by Cambridge scientists Shankar Balasubramanian and David Klenerman in the Department of Chemistry, who formed a company called Solexa in 1998. Illumina acquired Solexa in 2007. A third generation of sequencing machines from Pacific Biosciences (PacBio), Oxford Nanopore Technologies (ONT) and Ion Torrent have become available more recently, with PacBio and ONT providing average read lengths of 10--18 kilobases.

This practical aims to familiarise you with Illumina next-generation sequencing (NGS) data and some of the software available for their analysis.

If you have any questions, please email Andy Tock (ajt200@cam.ac.uk).

## The data

The NGS data we are going to analyse are derived from whole-genome sequencing of the Landsberg *erecta* (L*er*) ecotype of the model plant species *Arabidopsis thaliana*, and were published in [Zapata et al. (2016) *PNAS* **113**](https://www.pnas.org/content/113/28/E4052). The data are paired-end reads and so there are two files (`SRR3156163_1.fastq.gz` contains the first read in each pair and `SRR3156163_2.fastq.gz` the second). The reads were downloaded from the the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/SRR3156163), which "provides a comprehensive record of the world's nucleotide sequencing information, covering raw sequencing data, sequence assembly information and functional annotation".

## The pipeline/workflow

Bioinformatics pipelines or workflows consist of sequential data processing and analysis steps that utilize different software tools, with the output file(s) from one step often serving as the input file(s) for the subsequent step(s). These pipelines require input files that conform to standarized data formats.

The goal of our pipeline is to identify DNA sequence differences (variants) in the genome of the L*er* ecotype of *A. thaliana* relative to the reference genome assembly for the Columbia (Col-0) ecotype. To this end, these are the steps in the pipeline that we will work through sequentially:

1. Evaluation of sequencing read quality, including at the level of individual bases
2. Removal of contaminants (e.g., sequencing adapters) and low-quality bases
3. Alignment of reads (L*er*) to a reference genome (Col-0)
4. Filtering of alignments based on the quality of these mappings to the reference genome
5. Detection of DNA sequence differences between the L*er* and Col-0 genomes (variant calling)

## Inspecting the reads in FASTQ format

The sequencing reads are contained in gzip-compressed [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files, a standardized format that NGS data analysis tools have been developed to handle. These files have been downloaded on the computers in the Craik-Marshall Building that you are accessing remotely, so there's no need to download them unless you are working on your own computer. The files are located in the `fastq/` directory (folder). 

Let's first have a look at one of the files to inspect its format. Open a terminal window and use `zcat` and `head` to print out the first eight lines of `SRR3156163_1.fastq.gz` to your screen. We need to use `zcat` here to uncompress the gzip-compressed file. The `|` pipes the output of the `zcat` command to the `head` command.

```
zcat fastq/SRR3156163_1.fastq.gz | head -n 8
```

### Output:
> @SRR3156163.1 1/1
> TTTGCTTGTNNNNNNNNNNNNNTCATCATGAANNNNNNNNNNNNNNNNNNGTCAGATACAANNNNNNNNNNNNNNTTGTGGAAGCAGGAGATGTGGNNGT
> +
> <<<@@????###########################################################################################
> @SRR3156163.2 2/1
> TGATTCGCTTNGNNNNNNNNNGTCGCCACAGCANNNNNNNNNNNNNNNNCGTATAGCATACNNNNNNNNNNNNNNTACGAGCTGCATTAAAGTAGCGCAG
> +
> <<<?@@????#3#########21@=????????###################################################################

The first four lines show data for one read, and the next four lines show data for the subsequent read, each corresponding to the first read in a pair of reads. The second read in each pair is contained in `SRR3156163_2.fastq.gz`.

Data in FASTQ format conform to these standards:

Line | Description
---- | -----------
1 | Begins with '@', followed by information about the read
2 | The DNA sequence
3 | Begins with '+'
4 | A character string of the same length as the sequence, encoding quality scores for each base

