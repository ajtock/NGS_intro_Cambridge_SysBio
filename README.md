# Introduction to next-generation sequencing data: tools and resources

This [link](https://ajtock.github.io/NGS_intro_Cambridge_SysBio/) will take you to the GitHub Pages rendering of this repository.

Following the success of dye-terminator (Sanger) sequencing in the Human Genome Project, a "next" generation of sequencing technologies became available, which were faster, cheaper and produced many more reads in parallel.
Despite their higher throughput, the reads were shorter and of lower quality than those produced by Sanger sequencing.
One high-throughput sequencing technology became ubiquitous, Illumina sequencing by synthesis (on GAI/GAII/HiSeq/ MiSeq/NextSeq instruments).
This technology was developed by Cambridge scientists Shankar Balasubramanian and David Klenerman in the Department of Chemistry, who formed a company called Solexa in 1998.
Illumina acquired Solexa in 2007.
A third generation of sequencing machines from Pacific Biosciences (PacBio), Oxford Nanopore Technologies (ONT) and Ion Torrent have become available more recently, with PacBio and ONT providing average read lengths of 10–18 kilobases.

This practical aims to familiarise you with Illumina next-generation sequencing (NGS) data and some of the software available for their analysis.

If you have any questions, please email Andy Tock at <ajt200@cam.ac.uk>.

## The data

The NGS data we are going to analyse are derived from whole-genome sequencing of the Landsberg *erecta* (L*er*) ecotype of the model plant species *Arabidopsis thaliana*, and were published in [Zapata et al. (2016) *PNAS* **113**](https://www.pnas.org/content/113/28/E4052).
The data are paired-end reads and so there are two files (`SRR3156163_1.fastq.gz` contains the first read in each pair and `SRR3156163_2.fastq.gz` the second).
The reads were downloaded from the the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/SRR3156163), which "provides a comprehensive record of the world's nucleotide sequencing information, covering raw sequencing data, sequence assembly information and functional annotation".

## The pipeline/workflow

Bioinformatics pipelines or workflows consist of sequential data processing and analysis steps that utilise different software tools, with the output file(s) from one step often serving as the input file(s) for the subsequent step(s).
These pipelines require input files that conform to standardised data formats.

The goal of our pipeline is to identify DNA sequence differences (variants) in the genome of the L*er* ecotype of *A. thaliana* relative to the reference genome assembly for the Columbia (Col-0) ecotype.
To this end, these are the steps in the pipeline that we will work through sequentially:

1. Evaluation of sequencing read quality, including at the level of individual bases
2. Removal of contaminants (e.g., sequencing adapters) and low-quality bases
3. Alignment of reads (L*er*) to a reference genome (Col-0)
4. Filtering of alignments based on the quality of these mappings to the reference genome
5. Detection of DNA sequence differences between the L*er* and Col-0 genomes (variant calling)

## Inspecting the reads in FASTQ format

The sequencing reads are contained in gzip-compressed [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files, a standardised format that NGS data analysis tools have been developed to handle.
These files have been downloaded on the computers in the Craik-Marshall Building that you are accessing remotely, so there's no need to download them unless you are working on your own computer.
The files are located in the `fastq/` directory.

Data in FASTQ format conform to these standards:

Line | Description
---- | -----------
1 | Begins with '@', followed by information about the read
2 | The DNA sequence
3 | Begins with '+'
4 | A character string of the same length as the sequence, encoding quality scores for each base

Let's first have a look at one of the files to inspect its format.
In a Unix command-line shell, use `zcat` and `head` to print to screen the first eight lines of `SRR3156163_1.fastq.gz`.
We need to use `zcat` here to uncompress the gzip-compressed file.
The `|` part pipes the output of the `zcat` command to the `head` command.

```
zcat fastq/SRR3156163_1.fastq.gz | head -n 8
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
The second read in each pair is contained in `SRR3156163_2.fastq.gz`.
The first line for each read contains a unique identifier and, as these are paired-end reads, `/1` indicates that this is the first read in the pair.

The quality score of each base identified in a sequencing read is encoded as a single character on the fourth line.
These represent [Phred quality scores](https://en.wikipedia.org/wiki/Phred_quality_score) that have been [converted into ASCII\_BASE 33 characters](https://drive5.com/usearch/manual/quality_score.html) such that each character encodes a quality score for the corresponding base in the read.
A Phred quality score is logarithmically related to the probability of an incorrect base call *P*, expressed as 1 error in 10<sup>*Q*/10</sup> base calls of *Q* quality, or

> *Q* = -10log<sub>10</sub>*P*  
> *P* = 10<sup>-*Q*/10</sup>  

Accordingly, the ASCII\_BASE 33 character `@` encodes a *Q*-score of 31 and a base-calling error probability of 0.00079.
In the past, Illumina sequencing instruments used the [ASCII\_BASE 64 quality encoding](https://drive5.com/usearch/manual/quality_score.html).

Is the first read composed of mostly high-quality or low-quality base calls?

### Exercise 1

Construct a command that will print to screen the 500th read in `SRR3156163_1.fastq.gz` in order to inspect its quality.

<details>
  <summary><em><strong>Solution</strong> (click to reveal/hide)</em></summary><p>

  ```
  zcat fastq/SRR3156163_1.fastq.gz | head -n 2000 | tail -n 4
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

## Evaluating read quality using FastQC

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
Started analysis of SRR3156163_1.fastq.gz
Approx 5% complete for SRR3156163_1.fastq.gz
Approx 10% complete for SRR3156163_1.fastq.gz
Approx 15% complete for SRR3156163_1.fastq.gz
Approx 20% complete for SRR3156163_1.fastq.gz
Approx 25% complete for SRR3156163_1.fastq.gz
Approx 30% complete for SRR3156163_1.fastq.gz
Approx 35% complete for SRR3156163_1.fastq.gz
Approx 40% complete for SRR3156163_1.fastq.gz
Approx 45% complete for SRR3156163_1.fastq.gz
Approx 50% complete for SRR3156163_1.fastq.gz
Approx 55% complete for SRR3156163_1.fastq.gz
Approx 60% complete for SRR3156163_1.fastq.gz
Approx 65% complete for SRR3156163_1.fastq.gz
Approx 70% complete for SRR3156163_1.fastq.gz
Approx 75% complete for SRR3156163_1.fastq.gz
Approx 80% complete for SRR3156163_1.fastq.gz
Approx 85% complete for SRR3156163_1.fastq.gz
Approx 90% complete for SRR3156163_1.fastq.gz
Approx 95% complete for SRR3156163_1.fastq.gz
Analysis complete for SRR3156163_1.fastq.gz
Started analysis of SRR3156163_2.fastq.gz
Approx 5% complete for SRR3156163_2.fastq.gz
Approx 10% complete for SRR3156163_2.fastq.gz
Approx 15% complete for SRR3156163_2.fastq.gz
Approx 20% complete for SRR3156163_2.fastq.gz
Approx 25% complete for SRR3156163_2.fastq.gz
Approx 30% complete for SRR3156163_2.fastq.gz
Approx 35% complete for SRR3156163_2.fastq.gz
Approx 40% complete for SRR3156163_2.fastq.gz
Approx 45% complete for SRR3156163_2.fastq.gz
Approx 50% complete for SRR3156163_2.fastq.gz
Approx 55% complete for SRR3156163_2.fastq.gz
Approx 60% complete for SRR3156163_2.fastq.gz
Approx 65% complete for SRR3156163_2.fastq.gz
Approx 70% complete for SRR3156163_2.fastq.gz
Approx 75% complete for SRR3156163_2.fastq.gz
Approx 80% complete for SRR3156163_2.fastq.gz
Approx 85% complete for SRR3156163_2.fastq.gz
Approx 90% complete for SRR3156163_2.fastq.gz
Approx 95% complete for SRR3156163_2.fastq.gz
Analysis complete for SRR3156163_2.fastq.gz
```

Navigate to the output directory and list its contents.

```
cd results/fastqc/raw_reads/
ls -1
```

### Output:
```
SRR3156163_1_fastqc.html
SRR3156163_1_fastqc.zip
SRR3156163_2_fastqc.html
SRR3156163_2_fastqc.zip
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

Have a look at the two HTML files by opening them in a web browser.

```
firefox *.html
```

