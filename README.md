# RNA-seq Quality Assessment Assignment - Bi 623 (Summer 2024)

Be sure to upload all relevant materials by the deadline and **double check** to be sure that your offline repository is up-to-date with your online repository. Answers to questions should be included in your final, high-level, report as a `pdf`. This pdf should be generated using Rmarkdown and submitted to Canvas as well as GitHub. Be sure to keep a well-organized, detailed lab notebook!

## Objectives
The objectives of this assignment are to use existing tools for quality assessment and adaptor trimming, compare the quality assessments to those from your own software, and to demonstrate your ability to summarize other important information about this RNA-Seq data set in a high-level report. That is, you should create a cohesive, [well written](https://canvas.uoregon.edu/courses/244873/modules/667769) report for your "PI" about what you've learned about/from your data.

### Data: 
Each of you will be working with 2 of the demultiplexed file pairs. For all steps below, process the two libraries separately. Library assignments are here: ```/projects/bgmp/shared/Bi623/QAA_data_assignments.txt```

The demultiplexed, gzipped `.fastq` files are here: ```/projects/bgmp/shared/2017_sequencing/demultiplexed/```

> [!WARNING]
> Do not move, copy, or unzip these data!
```
______                    _                                                               
|  _  \                  | |                                                              
| | | |___    _ __   ___ | |_   _ __ ___   _____   _____      ___ ___  _ __  _   _        
| | | / _ \  | '_ \ / _ \| __| | '_ ` _ \ / _ \ \ / / _ \    / __/ _ \| '_ \| | | |       
| |/ / (_) | | | | | (_) | |_  | | | | | | (_) \ V /  __/_  | (_| (_) | |_) | |_| |_      
|___/ \___/  |_| |_|\___/ \__| |_| |_| |_|\___/ \_/ \___( )  \___\___/| .__/ \__, ( )     
                                                        |/            | |     __/ |/      
                                                                      |_|    |___/        
                              _         _   _                          _       _        _ 
                             (_)       | | | |                        | |     | |      | |
  ___  _ __   _   _ _ __  _____ _ __   | |_| |__   ___  ___  ___    __| | __ _| |_ __ _| |
 / _ \| '__| | | | | '_ \|_  / | '_ \  | __| '_ \ / _ \/ __|/ _ \  / _` |/ _` | __/ _` | |
| (_) | |    | |_| | | | |/ /| | |_) | | |_| | | |  __/\__ \  __/ | (_| | (_| | || (_| |_|
 \___/|_|     \__,_|_| |_/___|_| .__/   \__|_| |_|\___||___/\___|  \__,_|\__,_|\__\__,_(_)
                               | |                                                        
                               |_|                                                        
```
> [!WARNING]
> Do not move, copy, or unzip these data!

# Part 1 – Read quality score distributions

1. Create a new conda environment called `QAA` and install `FastQC`. Google around if you need a refresher on how to create conda environments. Recommend doing this in an interactive session, not the login node! Record details of how you created this environment in your lab notebook! Make sure you check your installation with:
   - `fastqc --version` (should be 0.12.1)  

2. Using `FastQC` via the command line on Talapas, produce plots of the per-base quality score distributions for R1 and R2 reads. Also, produce plots of the per-base N content, and comment on whether or not they are consistent with the quality score plots.

3. Run your quality score plotting script from your Demultiplexing assignment in Bi622. (Make sure you're using the "running sum" strategy!!) Describe how the `FastQC` quality score distribution plots compare to your own. If different, propose an explanation. Also, does the runtime differ? Mem/CPU usage? If so, why?

4. Comment on the overall data quality of your two libraries. Go beyond per-base qscore distributions. Make and justify a recommendation on whether these data are of high enough quality to use for further analysis.

# Part 2 – Adaptor trimming comparison

5.  In your QAA environment, install `cutadapt` and `Trimmomatic`. Check your installations with:
    - `cutadapt --version` (should be 4.9)
    - `trimmomatic -version` (should be 0.39)

6. Using `cutadapt`, properly trim adapter sequences from your assigned files. Be sure to read how to use `cutadapt`. Use default settings. What proportion of reads (both R1 and R2) were trimmed?

    <details>
    <summary>Try to determine what the adapters are on your own. If you cannot (or if you do, and want to confirm), click here to see the actual adapter sequences used.</summary>
  
    R1: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
    
    R2: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
    </details>

    - *Sanity check*: Use your Unix skills to search for the adapter sequences in your datasets and confirm the expected sequence orientations. Report the commands you used, the reasoning behind them, and how you confirmed the adapter sequences.

7. Use `Trimmomatic` to quality trim your reads. Specify the following, **in this order**:
    - LEADING: quality of 3
    - TRAILING: quality of 3
    - SLIDING WINDOW: window size of 5 and required quality of 15
    - MINLENGTH: 35 bases

    Be sure to output compressed files and clear out any intermediate files.

8. Plot the trimmed read length distributions for both R1 and R2 reads (on the same plot - yes, you will have to use Python or R to plot this. See ICA4 from Bi621). You can produce 2 different plots for your 2 different RNA-seq samples. There are a number of ways you could possibly do this. One useful thing your plot should show, for example, is whether R1s are trimmed more extensively than R2s, or vice versa. Comment on whether you expect R1s and R2s to be adapter-trimmed at different rates and why.

9. CHALLENGE - Run `FastQC` on your trimmed data. Comment on differences you observe between the trimmed and untrimmed data. Include any figures needed to support your conclusions.

# Part 3 – Alignment and strand-specificity
10. Install sofware (record details in lab notebook!!!). In your QAA environment, use conda to install:
    - star
    - numpy
    - matplotlib
    - htseq

11. Find publicly available mouse genome fasta files (Ensemble release 112) and generate an alignment database from them. Align the reads to your mouse genomic database using a splice-aware aligner. Use the settings specified in PS8 from Bi621.

  > [!IMPORTANT]
  > You will need to use gene models to perform splice-aware alignment, see PS8 from Bi621.
    
12. Using your script from PS8 in Bi621, report the number of mapped and unmapped reads from each of your 2 sam files. Make sure that your script is looking at the bitwise flag to determine if reads are primary or secondary mapping (update/fix your script if necessary).

13. Count reads that map to features using `htseq-count`. You should run htseq-count twice: once with `--stranded=yes` and again with `--stranded=reverse`. Use default parameters otherwise.

14. Demonstrate convincingly whether or not the data are from "strand-specific" RNA-Seq libraries. Include any comands/scripts used. Briefly describe your evidence, using quantitative statements (e.g. "I propose that these data are/are not strand-specific, because X% of the reads are y, as opposed to z.").

  > [!TIP]
  > Recall ICA4 from Bi621.

# Challenge (optional!)

Review the [metadata](./metadata) available to you and see if this information leads to any additional insight of your analysis.

# To turn in your work for this assignment

## Upload your:
- [ ] lab notebook,
- [ ] Talapas batch script/code, 
- [ ] FastQC plots, 
- [ ] counts files generated from htseq-count (in a folder would be nice),
- [ ] pdf report (see below), 
- [ ] and any additional plots, code, or code output

to GitHub.
    
## You should create a pdf file (using Rmarkdown) with a high-level report including:
- [ ] all requested plots
- [ ] answers to questions
- [ ] mapped/unmapped read counts from PS8 script (in a nicely formatted table)
    
The three parts of the assignment should be clearly labeled. Be sure to title and write a descriptive figure caption for each image/graph/table you present. 
> [!TIP]
> Think about figure captions you've read and discussed in Journal Club. Find some good examples to model your captions on.

The file should be named `QAA_report.pdf`, and it should be a the top level of your repo AND submitted to Canvas.
