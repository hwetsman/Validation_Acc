# Validation_Acc

## Background

When a lab is doing validation, they must prove their accuracy. To do this they compare their results to the results on the same sample from another lab. In the case of genetics, there are known positive controls from the [1000genomes project](https://www.internationalgenome.org) that are available from the [Coriell Institute](https://www.coriell.org) under license from the US government. The genetic values from these samples were generated by tax dollars and the information is freely available from [NCBI](https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/) or [Ensembl](https://www.ensembl.org).

A lab can buy samples of these positive controls from Coriell (the samples are call Coriells) and download the expected result from public sources. They then report their findings against public sources in their validation document. This validation document is a tedious affair of transcribing numbers and statistics from one spreadsheet to another. This work is usually done by lab validation consultants at a fairly high cost, limiting opening labs to those who are better capitalized. This work of validation computation can be done automatically. This project is the seed that hopes to grow into that goal.

## Inputs
The script will take inputs in the form of .vcf files. For each sample, there are two .vcf files. The first is provided by via bioinformatics pathway from data on the Coriell sample in the 1000genomes project. The second is provided by the lab.

### File Naming
The input files for analysis will all be in the <lab> directory inside a <panel> directory. The files must be named with a strict naming convention. Positive control samples must be named as follows: "<lab>_<panel>_<samplename>_Pos_Control.vcf" where carrots show the inputs for the name. <lab> must match the <lab> directory name. <panel> must match the <panel> directory name. Validation run samples must be named as follows: "<lab>_<panel>_<samplename>_<4digityear>-<2digitmonth>-<2digitday>_<run>.vcf". This will allow the program to identify the positive controls, identify the run files and the order and spacing they are in. It will allow further expansion of the logic to go beyond accuracy to precision and beyond coriells to patient controls.

## Action
The script assigns the 1000genomes vcf the status of gold standard and measures the accuracy of the lab derived .vcf file against it. By the nature of modern methods, much deeper reads are possible with NGS equipment than was available to the 1000genomes project. This means that the lab will read more than what the bioinformatics derived standard shows. The script ignores all findings from the lab that are not in the standard. This is not to say that the script ignores false positives from the lab.

Because of the nature of 1000genome data, the same SNPs were measured in all subjects, and, even if the standard isn't positive for them, the information is available.  This allows us to ignore the reads not present in the standard and still measure both false negative and false positive errors.
