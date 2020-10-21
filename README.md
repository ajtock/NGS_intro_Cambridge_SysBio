# Introduction to next-generation sequencing data: tools and resources

This [link](https://ajtock.github.io/NGS_intro_Cambridge_SysBio/) will take you to the GitHub Pages rendering of this repository.

If you have any questions or comments about the practical, please email Andy Tock at <ajt200@cam.ac.uk>.

## Background

[The cost of DNA sequencing has decreased dramatically since the end of 2007](https://www.genome.gov/about-genomics/fact-sheets/DNA-Sequencing-Costs-Data).
This was the point at which sequencing efforts increasingly moved away from the [Sanger dye-terminator method](https://en.wikipedia.org/wiki/Sanger_sequencing), to a "next" generation of sequencing technologies, which are faster, cheaper and produce many more reads in parallel.
Despite their higher throughput, the reads are shorter and of lower quality than those produced by Sanger sequencing.

One next-generation sequencing (NGS) technology has come to dominate: [Illumina sequencing by synthesis (SBS)](https://www.illumina.com/documents/products/techspotlights/techspotlight_sequencing.pdf) on Genome Analyzer, HiSeq, MiSeq, NextSeq and NovaSeq instruments.
This technology was developed by Cambridge scientists Shankar Balasubramanian and David Klenerman in the Department of Chemistry, who formed a company called Solexa in 1998.
Illumina acquired Solexa in 2007.
The method uses flow cells containing clusters of clonally amplified DNA fragments.
The underlying nucleotide sequence of each fragment is determined based on signal emitted by each cluster under laser excitation of a fluorescently labelled nucleotide that has been added to an extending sequence complementary to the fragment.
Billions of clusters per flow cell are sequenced simultaneously, and an image is captured for each successive base along the DNA fragment, with read length dependent on the method employed. 
This [animation](https://www.youtube.com/watch?v=fCd6B5HRaZ8) shows in detail how Illumina sequencing by synthesis works.

Due in large part to this breakthrough technology, we have seen the cost of sequencing a single person’s genome drop from a prohibitive US$10 million in 2007 to ~US$1,000 by mid-2015!
There have been various refinements to Illumina’s technology, including the use of [patterned flow cells](https://emea.illumina.com/science/technology/next-generation-sequencing/sequencing-technology/patterned-flow-cells.html) that have billions of ordered nanowells.
These patterned flow cells provide optimal "fixed" rather than random spacing between clusters of clonally amplified DNA fragments, for increased cluster density and simplified imaging.

A third generation of sequencing machines from [Pacific Biosciences (PacBio)](https://www.pacb.com/smrt-science/smrt-sequencing/) and [Oxford Nanopore Technologies (ONT)](https://nanoporetech.com/how-it-works) have become available more recently, providing average read lengths of 10–18 kilobases.

A number of considerations come into play when deciding on what kind of sequencing technology to use for a project.
Cost and time requirements differ between technologies.
For example, the costs associated with using newer technology to generate long sequencing reads might not be justifiable or efficient for the intended application.
Some applications require high sequencing depth of coverage of specific parts of the genome (i.e., large numbers of reads overlapping those genomic positions; e.g., identification of protein binding sites using [ChIP-seq](https://en.wikipedia.org/wiki/ChIP_sequencing), or gene expression profiling via [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq)).
Others require greater breadth of coverage of the genome, such as projects aimed at assembling genomes or genomic regions from scratch, termed [*de novo* assembly](https://en.wikipedia.org/wiki/De_novo_sequence_assemblers).
Therefore, the different sequencing read lengths, qualities, abundances and error rates that are associated with different technologies need to be taken into account in order to make the right choice.

This practical aims to familiarise you with Illumina next-generation sequencing data and some of the software available for their analysis.

## The data

The NGS data we are going to analyse are derived from whole-genome sequencing of the Landsberg *erecta* (L*er*) [ecotype](https://en.wikipedia.org/wiki/Ecotype) of the [diploid](https://www.genome.gov/genetics-glossary/Diploid) model plant species [*Arabidopsis thaliana*](https://en.wikipedia.org/wiki/Arabidopsis_thaliana), and were published in [Zapata et al. (2016) *PNAS* **113**](https://www.pnas.org/content/113/28/E4052).
The data are paired-end reads and so there are two files (`SRR3156163_top5M_1.fastq.gz` contains the first read in each pair and `SRR3156163_top5M_2.fastq.gz` the second).
Each read in a pair was sequenced with 100 chemistry cycles on an [Illumina HiSeq 2000](https://www.illumina.com/documents/products/datasheets/datasheet_hiseq2000.pdf), generating 100 consecutive base calls per read (2×100 bp).
The reads were downloaded from the the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/SRR3156163), which "provides a comprehensive record of the world's nucleotide sequencing information, covering raw sequencing data, sequence assembly information and functional annotation".
The top 5 million reads (`top5M`) in each of the two files were extracted in order to reduce time spent on data processing in this practical.

## The pipeline/workflow

Bioinformatics pipelines or workflows consist of sequential data processing and analysis steps that utilise different software tools, with the output file(s) from one step often serving as the input file(s) for the subsequent step(s).
These pipelines require input files that conform to standardised data formats.

The goal of our pipeline is to identify DNA sequence differences (variants) in the genome of the L*er* ecotype of *A. thaliana* relative to the reference genome assembly for the Columbia (Col-0) ecotype.
To this end, these are the steps in the pipeline that we will work through sequentially:

1. Evaluation of sequencing read quality, including at the level of individual bases
2. Removal of technical sequences (e.g., sequencing adapters) and low-quality bases
3. Alignment of reads (L*er*) to a reference genome (Col-0)
4. Filtering of alignments based on the quality of these mappings to the reference genome
5. Detection of DNA sequence differences between the L*er* and Col-0 genomes (variant calling)

## Inspecting the reads in FASTQ format

The sequencing reads are contained in gzip-compressed [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files, a standardised format that NGS data analysis tools have been developed to handle.
These files have been downloaded on the computers in the Craik-Marshall Building that you are accessing remotely, so there's no need to download them unless you are working on your own computer.
The files are located in the `fastq/` directory.

Data in FASTQ format conform to these standards:

| Line | Description |
|:-----|:------------|
| 1    | Begins with '@', followed by information about the read
| 2    | The DNA sequence
| 3    | Begins with '+'
| 4    | A character string of the same length as the sequence, encoding quality scores for each base

Let's first have a look at one of the files to inspect its format.
In a Unix command-line shell, use `zcat` and `head` to print to screen the first eight lines of `SRR3156163_top5M_1.fastq.gz`.
We need to use `zcat` here to uncompress the gzip-compressed file.
The `|` part pipes the output of the `zcat` command to the `head` command.

```
zcat fastq/SRR3156163_top5M_1.fastq.gz | head -n 8
```

### Output:
```
@SRR3156163.1 1/1
TTTGCTTGTNNNNNNNNNNNNNTCATCATGAANNNNNNNNNNNNNNNNNNGTCAGATACAANNNNNNNNNNNNNNTTGTGGAAGCAGGAGATGTGGNNGT
+
<<<@@????###########################################################################################
@SRR3156163.2 2/1
TGATTCGCTTNGNNNNNNNNNGTCGCCACAGCANNNNNNNNNNNNNNNNCGTATAGCATACNNNNNNNNNNNNNNTACGAGCTGCATTAAAGTAGCGCAG
+
<<<?@@????#3#########21@=????????###################################################################
```

The first four lines show data for one read and the next four lines show data for the subsequent read, each corresponding to the first read in a pair of reads.
The second read in each pair is contained in `SRR3156163_top5M_2.fastq.gz`.
The first line for each read contains a unique identifier and, as these are paired-end reads, `/1` indicates that this is the first read in the pair.

The quality score of each base identified in a sequencing read is encoded as a single character on the fourth line.
These represent [Phred quality scores](https://en.wikipedia.org/wiki/Phred_quality_score) that have been [converted into ASCII\_BASE=33 characters](https://drive5.com/usearch/manual/quality_score.html) such that each character encodes a quality score for the corresponding base in the read.
A Phred quality score is logarithmically related to the probability of an incorrect base call *P*, expressed as 1 error in 10<sup>*Q*/10</sup> base calls of *Q* quality, or

> *Q* = -10log<sub>10</sub>*P*  
> *P* = 10<sup>-*Q*/10</sup>  

Accordingly, the ASCII\_BASE 33 character `@` encodes a *Q*-score of 31 and a base-calling error probability of 0.00079.
In the past, Illumina sequencing instruments used the [ASCII\_BASE=64 quality encoding](https://drive5.com/usearch/manual/quality_score.html).

Is the first read composed of mostly high-quality or low-quality base calls?

### Exercise 1

Construct a command that will print to screen the 500th read in `SRR3156163_top5M_1.fastq.gz` in order to inspect its quality.

<details>
  <summary><em><strong>Solution</strong> (click to reveal/hide)</em></summary><p>

  ```
  zcat fastq/SRR3156163_top5M_1.fastq.gz | head -n 2000 | tail -n 4
  ```

  #### Output:
  ```
  @SRR3156163.500 500/1
  GAGGAAGCTTGACGCAGCGGAGGAATCTTTGCTGACCCCATCGGTCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACGGTCGTGGAGCAAGG
  +
  CCCFFFFFHHGHHIJJJJJJGIIFIJIJJJJJJJIJIJJJJJJJCGHFFFFFFDEDD;@CDDCCBCACDDCDEDDBBB8ADEDDD<>?BDDD>BDCBCC?
  ```
</p></details>

Is the 500th read generally better or worse than the first read?

## Step 1. Evaluating read quality using FastQC

Read quality can be evaluated in a more systematic way using dedicated software, such as [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
Let's make sure FastQC is installed by invoking the executable file `fastqc` with an option that will report the version.

```
fastqc --version
```

### Output (abbreviated):
```
FastQC v0.11.9
```

The `--version` and `--help` options are usually, but not always, available for most tools.
The latter generally shows an example of typical usage of the executable file with arguments, along with a list of options for the program.

```
fastqc --help
```

### Output:
```

            FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

	fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of 
    which will help to identify a different potential type of problem in your
    data.
    
    If no files to process are specified on the command line then the program
    will start as an interactive graphical application.  If files are provided
    on the command line then the program will run with no user interaction
    required.  In this mode it is suitable for inclusion into a standardised
    analysis pipeline.
    
    The options for the program as as follows:
    
    -h --help       Print this help file and exit
    
    -v --version    Print the version of the program and exit
    
    -o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the 
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.
                    
    --casava        Files come from raw casava output. Files in the same sample
                    group (differing only by the group number) will be analysed
                    as a set rather than individually. Sequences with the filter
                    flag set in the header will be excluded from the analysis.
                    Files must have the same names given to them by casava
                    (including being gzipped and ending with .gz) otherwise they
                    won't be grouped together correctly.
                    
    --nano          Files come from nanopore sequences and are in fast5 format. In
                    this mode you can pass in directories to process and the program
                    will take in all fast5 files within those directories and produce
                    a single output file from the sequences found in all files.                    
                    
    --nofilter      If running with --casava then don't remove read flagged by
                    casava as poor quality when performing the QC analysis.
                   
    --extract       If set then the zipped output file will be uncompressed in
                    the same directory after it has been created.  By default
                    this option will be set if fastqc is run in non-interactive
                    mode.
                    
    -j --java       Provides the full path to the java binary you want to use to
                    launch fastqc. If not supplied then java is assumed to be in
                    your path.
                   
    --noextract     Do not uncompress the output file after creating it.  You
                    should set this option if you do not wish to uncompress
                    the output when running in non-interactive mode.
                    
    --nogroup       Disable grouping of bases for reads >50bp. All reports will
                    show data for every base in the read.  WARNING: Using this
                    option will cause fastqc to crash and burn if you use it on
                    really long reads, and your plots may end up a ridiculous size.
                    You have been warned!
                    
    --min_length    Sets an artificial lower limit on the length of the sequence
                    to be shown in the report.  As long as you set this to a value
                    greater or equal to your longest read length then this will be
                    the sequence length used to create your read groups.  This can
                    be useful for making directly comaparable statistics from 
                    datasets with somewhat variable read lengths.
                    
    -f --format     Bypasses the normal sequence file format detection and
                    forces the program to use the specified format.  Valid
                    formats are bam,sam,bam_mapped,sam_mapped and fastq
                    
    -t --threads    Specifies the number of files which can be processed
                    simultaneously.  Each thread will be allocated 250MB of
                    memory so you shouldn't run more threads than your
                    available memory will cope with, and not more than
                    6 threads on a 32 bit machine
                  
    -c              Specifies a non-default file which contains the list of
    --contaminants  contaminants to screen overrepresented sequences against.
                    The file must contain sets of named contaminants in the
                    form name[tab]sequence.  Lines prefixed with a hash will
                    be ignored.

    -a              Specifies a non-default file which contains the list of
    --adapters      adapter sequences which will be explicity searched against
                    the library. The file must contain sets of named adapters
                    in the form name[tab]sequence.  Lines prefixed with a hash
                    will be ignored.
                    
    -l              Specifies a non-default file which contains a set of criteria
    --limits        which will be used to determine the warn/error limits for the
                    various modules.  This file can also be used to selectively 
                    remove some modules from the output all together.  The format
                    needs to mirror the default limits.txt file found in the
                    Configuration folder.
                    
   -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                    module. Specified Kmer length must be between 2 and 10. Default
                    length is 7 if not specified.
                    
   -q --quiet       Supress all progress messages on stdout and only report errors.
   
   -d --dir         Selects a directory to be used for temporary files written when
                    generating report images. Defaults to system temp directory if
                    not specified.
                    
BUGS

    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
    or in www.bioinformatics.babraham.ac.uk/bugzilla/
                   
    
```

In the above output, you will see the `-o --outdir` option, which allows you to specify an output directory to which the FastQC results will be written.
However, this directory must exist before running FastQC, so we'd better make it if we want to make use of this option.
The `-p` option in the `mkdir` command below allows us to make a new directory containing subdirectories.

```
mkdir -p results/fastqc/raw_reads
```

Now let's run FastQC on each of our two gzip-compressed FASTQ files located in the `fastq/` directory by using the `*.fastq.gz` wildcard.
Compressed or uncompressed files can be provided as inputs.

```
fastqc --outdir results/fastqc/raw_reads \
       fastq/*.fastq.gz
```

Progress made by FastQC on analysing each file will be printed to screen.

### Output:
```
Started analysis of SRR3156163_top5M_1.fastq.gz
Approx 5% complete for SRR3156163_top5M_1.fastq.gz
Approx 10% complete for SRR3156163_top5M_1.fastq.gz
Approx 15% complete for SRR3156163_top5M_1.fastq.gz
Approx 20% complete for SRR3156163_top5M_1.fastq.gz
Approx 25% complete for SRR3156163_top5M_1.fastq.gz
Approx 30% complete for SRR3156163_top5M_1.fastq.gz
Approx 35% complete for SRR3156163_top5M_1.fastq.gz
Approx 40% complete for SRR3156163_top5M_1.fastq.gz
Approx 45% complete for SRR3156163_top5M_1.fastq.gz
Approx 50% complete for SRR3156163_top5M_1.fastq.gz
Approx 55% complete for SRR3156163_top5M_1.fastq.gz
Approx 60% complete for SRR3156163_top5M_1.fastq.gz
Approx 65% complete for SRR3156163_top5M_1.fastq.gz
Approx 70% complete for SRR3156163_top5M_1.fastq.gz
Approx 75% complete for SRR3156163_top5M_1.fastq.gz
Approx 80% complete for SRR3156163_top5M_1.fastq.gz
Approx 85% complete for SRR3156163_top5M_1.fastq.gz
Approx 90% complete for SRR3156163_top5M_1.fastq.gz
Approx 95% complete for SRR3156163_top5M_1.fastq.gz
Analysis complete for SRR3156163_top5M_1.fastq.gz
Started analysis of SRR3156163_top5M_2.fastq.gz
Approx 5% complete for SRR3156163_top5M_2.fastq.gz
Approx 10% complete for SRR3156163_top5M_2.fastq.gz
Approx 15% complete for SRR3156163_top5M_2.fastq.gz
Approx 20% complete for SRR3156163_top5M_2.fastq.gz
Approx 25% complete for SRR3156163_top5M_2.fastq.gz
Approx 30% complete for SRR3156163_top5M_2.fastq.gz
Approx 35% complete for SRR3156163_top5M_2.fastq.gz
Approx 40% complete for SRR3156163_top5M_2.fastq.gz
Approx 45% complete for SRR3156163_top5M_2.fastq.gz
Approx 50% complete for SRR3156163_top5M_2.fastq.gz
Approx 55% complete for SRR3156163_top5M_2.fastq.gz
Approx 60% complete for SRR3156163_top5M_2.fastq.gz
Approx 65% complete for SRR3156163_top5M_2.fastq.gz
Approx 70% complete for SRR3156163_top5M_2.fastq.gz
Approx 75% complete for SRR3156163_top5M_2.fastq.gz
Approx 80% complete for SRR3156163_top5M_2.fastq.gz
Approx 85% complete for SRR3156163_top5M_2.fastq.gz
Approx 90% complete for SRR3156163_top5M_2.fastq.gz
Approx 95% complete for SRR3156163_top5M_2.fastq.gz
Analysis complete for SRR3156163_top5M_2.fastq.gz
```

Navigate to the output directory and list its contents.

```
cd results/fastqc/raw_reads/
ls -1
```

### Output:
```
SRR3156163_top5M_1_fastqc.html
SRR3156163_top5M_1_fastqc.zip
SRR3156163_top5M_2_fastqc.html
SRR3156163_top5M_2_fastqc.zip
```

### FastQC summary statistics and graphs

Each HTML file contains statistics and graphs summarising the FastQC results:

* [Basic statistics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/1%20Basic%20Statistics.html)
* [Per base sequence quality](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)
* [Per sequence quality scores](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html)
* [Per base sequence content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html)
* [Per sequence GC content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html)
* [Per base N content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html)
* [Sequence length distribution](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)
* [Sequence duplication levels](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)
* [Overrepresented sequences](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html)
* [Adapter content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html)
* [Per tile sequence quality](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/12%20Per%20Tile%20Sequence%20Quality.html)

### Exercise 2

Have a look at the FastQC-generated HTML file for each FASTQ file by opening them in a web browser.

```
firefox SRR3156163_top5M_1_fastqc.html SRR3156163_top5M_2_fastqc.html
```

What quality issues can you see in these reports, and how could they be fixed?

<details>
  <summary><em><strong>Solution</strong> (click to reveal/hide)</em></summary><p>

They are generally high-quality sequencing reads. However, there are some issues that should be addressed:
1. Per-base sequence quality decreases towards the ends of the reads (particularly towards their 3′ ends)
2. There are many duplicated reads, which may have resulted from [PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction) amplification biases
3. Illumina TruSeq adapter sequences are over-represented among reads in `SRR3156163_top5M_1.fastq.gz`

The first and third of these issues can be resolved using software developed to trim off sequencing adapters and low-quality bases.
Duplication can be addressed by discarding either duplicate reads or duplicate alignments to a reference genome.
</p></details>

If you are working with many FASTQ files, [MultiQC](https://multiqc.info/) can be used to aggregate FastQC-generated results and compile one HTML report that's easier to digest than individual reports for each sample. 

## Step 2. Removing technical sequences and low-quality bases using Cutadapt

There are several tools available for filtering and trimming reads to remove technical sequences (e.g., sequencing adapters) and low-quality bases, including [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
Removing these sequences is important because it means that downstream analyses won't be compromised by base calls in which we have low confidence, or by the presence of technical sequences that do not reflect the biology of the sample we have sequenced.
In the case of aligning reads to a reference genome assembly, for example, read cleaning tends to increase the alignment rate.

We're going to use Cutadapt for this step in the pipeline, so let's have a look at a usage example and the available options.

```
cutadapt --help
```

### Output:
```
cutadapt version 2.10

Copyright (C) 2010-2020 Marcel Martin <marcel.martin@scilifelab.se>

cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. All reads from input.fastq will be written to
output.fastq with the adapter sequence removed. Adapter matching is
error-tolerant. Multiple adapter sequences can be given (use further -a
options), but only the best-matching adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.

Citation:

Marcel Martin. Cutadapt removes adapter sequences from high-throughput
sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011.
http://dx.doi.org/10.14806/ej.17.1.200

Run "cutadapt --help" to see all command-line options.
See https://cutadapt.readthedocs.io/ for full documentation.

Options:
  -h, --help            Show this help message and exit
  --version             Show version number and exit
  --debug [{trace}]     Print debug log. 'trace' prints also DP matrices
  -j CORES, --cores CORES
                        Number of CPU cores to use. Use 0 to auto-detect. Default: 1

Finding adapters:
  Parameters -a, -g, -b specify adapters to be removed from each read (or from the first read in a pair if data is paired).
  If specified multiple times, only the best matching adapter is trimmed (but see the --times option). When the special
  notation 'file:FILE' is used, adapter sequences are read from the given FASTA file.

  -a ADAPTER, --adapter ADAPTER
                        Sequence of an adapter ligated to the 3' end (paired data: of the first read). The adapter and
                        subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found
                        if it is a suffix of the read.
  -g ADAPTER, --front ADAPTER
                        Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any
                        preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended
                        ('anchoring'), the adapter is only found if it is a prefix of the read.
  -b ADAPTER, --anywhere ADAPTER
                        Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of the first read). Both
                        types of matches as described under -a und -g are allowed. If the first base of the read is part of the
                        match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed
                        library preparations - do not use if you know which end your adapter was ligated to!
  -e RATE, --error-rate RATE
                        Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching
                        region). Default: 0.1 (=10%)
  --no-indels           Allow only mismatches in alignments. Default: allow both mismatches and indels
  -n COUNT, --times COUNT
                        Remove up to COUNT adapters from each read. Default: 1
  -O MINLENGTH, --overlap MINLENGTH
                        Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
  --match-read-wildcards
                        Interpret IUPAC wildcards in reads. Default: False
  -N, --no-match-adapter-wildcards
                        Do not interpret IUPAC wildcards in adapters.
  --action {trim,mask,lowercase,none}
                        What to do with found adapters. mask: replace with 'N' characters; lowercase: convert to lowercase;
                        none: leave unchanged (useful with --discard-untrimmed). Default: trim
  --rc, --revcomp       Check both the read and its reverse complement for adapter matches. If match is on reverse-complemented
                        version, output that one. Default: check only read

Additional read modifications:
  -u LENGTH, --cut LENGTH
                        Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the
                        beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have
                        different signs. This is applied *before* adapter trimming.
  --nextseq-trim 3'CUTOFF
                        NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G
                        bases.
  -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff [5'CUTOFF,]3'CUTOFF
                        Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both
                        reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated
                        cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
  --quality-base N      Assume that quality values in FASTQ are encoded as ascii(quality + N). This needs to be set to 64 for
                        some old Illumina FASTQ files. Default: 33
  --length LENGTH, -l LENGTH
                        Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at
                        the beginning. This and the following modifications are applied after adapter trimming.
  --trim-n              Trim N's on ends of reads.
  --length-tag TAG      Search for TAG followed by a decimal number in the description field of the read. Replace the decimal
                        number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct
                        fields like 'length=123'.
  --strip-suffix STRIP_SUFFIX
                        Remove this suffix from read names if present. Can be given multiple times.
  -x PREFIX, --prefix PREFIX
                        Add this prefix to read names. Use {name} to insert the name of the matching adapter.
  -y SUFFIX, --suffix SUFFIX
                        Add this suffix to read names; can also include {name}
  --zero-cap, -z        Change negative quality values to zero.

Filtering of processed reads:
  Filters are applied after above read modifications. Paired-end reads are always discarded pairwise (see also --pair-
  filter).

  -m LEN[:LEN2], --minimum-length LEN[:LEN2]
                        Discard reads shorter than LEN. Default: 0
  -M LEN[:LEN2], --maximum-length LEN[:LEN2]
                        Discard reads longer than LEN. Default: no limit
  --max-n COUNT         Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted
                        as a fraction of the read length.
  --max-expected-errors ERRORS, --max-ee ERRORS
                        Discard reads whose expected number of errors (computed from quality values) exceeds ERRORS.
  --discard-trimmed, --discard
                        Discard reads that contain an adapter. Use also -O to avoid discarding too many randomly matching
                        reads.
  --discard-untrimmed, --trimmed-only
                        Discard reads that do not contain an adapter.
  --discard-casava      Discard reads that did not pass CASAVA filtering (header has :Y:).

Output:
  --quiet               Print only error messages.
  --report {full,minimal}
                        Which type of report to print: 'full' or 'minimal'. Default: full
  -o FILE, --output FILE
                        Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. Summary report is sent
                        to standard output. Use '{name}' for demultiplexing (see docs). Default: write to standard output
  --fasta               Output FASTA to standard output even on FASTQ input.
  -Z                    Use compression level 1 for gzipped output files (faster, but uses more space)
  --info-file FILE      Write information about each read and its adapter matches into FILE. See the documentation for the file
                        format.
  -r FILE, --rest-file FILE
                        When the adapter matches in the middle of a read, write the rest (after the adapter) to FILE.
  --wildcard-file FILE  When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE.
                        (Inaccurate with indels.)
  --too-short-output FILE
                        Write reads that are too short (according to length specified by -m) to FILE. Default: discard reads
  --too-long-output FILE
                        Write reads that are too long (according to length specified by -M) to FILE. Default: discard reads
  --untrimmed-output FILE
                        Write reads that do not contain any adapter to FILE. Default: output to same file as trimmed reads

Paired-end options:
  The -A/-G/-B/-U options work like their -a/-b/-g/-u counterparts, but are applied to the second read in each pair.

  -A ADAPTER            3' adapter to be removed from second read in a pair.
  -G ADAPTER            5' adapter to be removed from second read in a pair.
  -B ADAPTER            5'/3 adapter to be removed from second read in a pair.
  -U LENGTH             Remove LENGTH bases from second read in a pair.
  -p FILE, --paired-output FILE
                        Write second read in a pair to FILE.
  --pair-adapters       Treat adapters given with -a/-A etc. as pairs. Either both or none are removed from each read pair.
  --pair-filter (any|both|first)
                        Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to
                        be filtered. Default: any
  --interleaved         Read and/or write interleaved paired-end reads.
  --untrimmed-paired-output FILE
                        Write second read in a pair to this FILE when no adapter was found. Use with --untrimmed-output.
                        Default: output to same file as trimmed reads
  --too-short-paired-output FILE
                        Write second read in a pair to this file if pair is too short. Use also --too-short-output.
  --too-long-paired-output FILE
                        Write second read in a pair to this file if pair is too long. Use also --too-long-output.
```

In the Cutadapt `--help` output above, we can see that the `-a` and `-A` options are used to specify the adapter sequences to be trimmed from the 3’ ends of Read 1 and Read 2 sequences, respectively.
An Illumina TruSeq DNA library preparation kit was used to generate the paired-end sequencing reads we are analysing.
Therefore, we need to consult the [Illumina Adapter Sequences Document](https://emea.support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html?langsel=/gb/) to locate the correct adapter information for inclusion in our Cutadapt command.
For most Illumina read types, including those derived from TruSeq libraries, [adapter trimming is required only at read 3’ ends](https://emea.support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html).

We should first make a subdirectory that will contain the Cutadapt-cleaned reads.

```
mkdir results/cutadapt/
```

### Exercise 3

Based on the options listed above and the [Cutadapt user guide](https://cutadapt.readthedocs.io/en/stable/guide.html#), write a Cutadapt command that will remove:
1. bases with Phred quality scores < 20 ([ASCII_BASE=33](https://drive5.com/usearch/manual/quality_score.html)) at the 3’ end of each read, as we have observed in the FastQC reports that base quality tends to degrade towards the 3’ ends of these reads, which is a general feature of Illumina reads
2. sequences that match a minimum of 4 consecutive bases in Illumina TruSeq adapters (Cutadapt will also remove any bases following [3’ of] a read–adapter match)
3. reads shorted than 30 bases, which will improve alignment performance as the shorter the sequence, the greater the chance that it will align to multiple locations in a reference genome

Cutadapt will print a report of progress and trimming statistics to standard output (stdout), which will be visible on your screen, in the same way FastQC printed its progress to screen.
Try to find a way to redirect the stdout and stderr (error messages, if any) generated by your `cutadapt` command to an appropriately named log file, which will be useful for future reference.

What proportion of read pairs and base calls passed the filters?

<details>
  <summary><em><strong>Solution</strong> (click to reveal/hide)</em></summary><p>

  ```
  (cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --quality-cutoff 20 \
            --overlap 4 \
            --minimum-length 30 \
            --output results/cutadapt/SRR3156163_top5M_1_trimmed.fastq.gz \
            --paired-output results/cutadapt/SRR3156163_top5M_2_trimmed.fastq.gz \
            fastq/SRR3156163_top5M_1.fastq.gz \
            fastq/SRR3156163_top5M_2.fastq.gz) &> results/cutadapt/SRR3156163_top5M_cutadapt_report.txt
  ```

  91% of read pairs and 88.4% of base calls remain after cleaning.
</p></details>

We can now use FastQC to evaluate the quality of the Cutadapt-trimmed reads, so let's make an output directory to contain the FastQC results.

```
mkdir results/fastqc/trimmed_reads
```

### Exercise 4

Write and run a `fastqc` command to evaluate the cleaned reads, specifying the newly created output directory for the FastQC results.
Once this has run to completion, inspect the FastQC-generated HTML reports to see if the cleaning step improved the reads.
Is any further cleaning required?

<details>
  <summary><em><strong>Solution</strong> (click to reveal/hide)</em></summary><p>

```
fastqc --outdir results/fastqc/trimmed_reads \
       results/cutadapt/*.fastq.gz
```

The quality of the reads has improved after cleaning with Cutadapt, but some warnings and so-called "failures" remain.

The warning for the [Per base sequence content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html) module can be safely ignored, as:
> "some types of library will always produce biased sequence composition, normally at the start of the read. Libraries produced by priming using random hexamers (including nearly all RNA-Seq libraries) and those which were fragmented using transposases inherit an intrinsic bias in the positions at which reads start. This bias does not concern an absolute sequence, but instead provides enrichement of a number of different K-mers at the 5’ end of the reads. Whilst this is a true technical bias, it isn't something which can be corrected by trimming and in most cases doesn't seem to adversely affect the downstream analysis."

The [Per sequence GC content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html) may also reflect this technical bias.

The [Sequence length distribution](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html) warning can also be ignored as we expect variable-length sequences after trimming.

"Failure" with regard to [Sequence duplication levels](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html), however, will need attention.
Duplication can be addressed by discarding either duplicate reads or duplicate alignments to a reference genome.
We will use the latter approach.
</p></details>

## Step 3. Aligning the cleaned reads to a reference genome assembly

Now we have removed adapters and low-quality bases from the sequencing reads, we can proceed with the next step in the pipeline: alignment of these cleaned read pairs from the Landsberg *erecta* (L*er*) ecotype to a reference genome assembly for the Columia (Col-0) ecotype.
These two [ecotypes](https://en.wikipedia.org/wiki/Ecotype) are genetically distinct true-breeding lines within the diploid plant species *Arabidopsis thaliana*.
Therefore, while they are genetically different from one another, they are each homozygous at every locus (position) on each chromosome.

Reference genome assemblies are stored in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format), with one of the following file extensions: `.fa`, `.fas` or `.fasta`.
For species with chromosome-level assemblies, each chromosome is represented as a separate nucleotide sequence in the FASTA file.
The first line of each separate sequence in the FASTA file is a description line, beginning with a greater-than symbol followed by a description of the sequence (e.g., `>Chr1`).
The nucleotide sequence is wrapped over the subsequent lines up until the next distinct sequence description line.

The reference genome assembly for *Arabidopsis thaliana* (`TAIR10_chr_all.fa`) is located in the `genome/` directory.
This was previously downloaded from [The Arabidopsis Information Resource (TAIR)](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_chromosome_files), so there's no need to download it for this practical.

### Exercise 5

Run a command that will output the description line of each separate DNA sequence in `genome/TAIR10_chr_all.fa`.
How many distinct sequences corresponding to nuclear chromosomes and non-nuclear DNA-containing organelles are there?

<details>
  <summary><em><strong>Solution</strong> (click to reveal/hide)</em></summary><p>

```
grep '>' genome/TAIR10_chr_all.fa > genome/TAIR10_chr_all_headers.txt
```

#### Output:
```
>1 CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2009-02-02
>2 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02
>3 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02
>4 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02
>5 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02
>mitochondria CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2005-06-03
>chloroplast CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2005-06-03
```

There are 5 nuclear chromosomes and 2 sequences corresponding to mitochondria and chloroplast organelles.
</p></details>


## Alignment using Bowtie 2


Alignment to a reference genome requires index files specific to the alignment software being used.
**You don't need to generate these index files in this case, as they have already been created in the `genome/` directory to save time.**
For future reference, the `bowtie2-build` command that was run to create these index files is:

```
bowtie2-build TAIR10_chr_all.fa TAIR10_chr_all
```

