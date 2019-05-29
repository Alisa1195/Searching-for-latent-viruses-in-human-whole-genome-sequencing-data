# Searching for latent viruses in human whole genome sequencing data
__Students:__ Yura Orlov, Nadya Pogodina, Alisa Morshneva (Bioinformatics Institute, Saint-Petersburg)\
__Supervisors:__ Valery Ilinsky, Alexander Rakitko (Genotek, Moscow)


## Project description
  The viral ability to stay in the asymptomatic phase of infection when a virus is not replicating is known as latency. A virus can stay latent for decades without exposing itself, although maintaining a capability to cause acute infections. Hence, viral presence in the human organism can remain unseen. 
    
  Whole-genome human sequencing data also contain sequences of dsDNA-viruses (or integrated RNA-viruses), because they are technically indistinguishable from the host DNA. These sequences stay unaligned and can be identified using databases. There are six viral families that store their genetic material in dsDNA form: Adenoviridae, Herpesviridae, Poxviridae and Polyomaviridae (linear dsDNA), Papovaviridae and Hepadnaviridae (circular dsDNA).
    
  In this project we were searching for viruses in WGS data from private Genotek database and 1000 Genomes project. On 5 population from 1000 Genomes project we performed GWAS (genome-wide association studies) for determine the associations between viral load level (EBV, Epstein-Barr virus) / viral presence (Mastadenoviruses) and SNP in human genome.

## Goal
Characterize viral load and representation in human WGS data and search for possible association between viral load and genetic variations

## Objectives
- [x] Explore literature data and find open WGS databases
- [x] Create a pipeline for WGS data analysis
- [x] Test the pipeline on different samples
- [x] Count viral load in testing data
- [x] Create a table of viral load for samples from 1000 Genomes
- [x] Perform GWAS

## Methods (Юра, проверь, всё ли верно)
#### Searching for viral reads in WGS data
The [bash-script](https://github.com/Alisa1195/Searching-for-latent-viruses-in-human-whole-genome-sequencing-data/blob/master/scripts/processing_script_v6_with_comments.sh) was written in order to extract reads from alignment files (samtools), identify the viral reads using Kraken2 (works with RefSeq database), parse report files to get information about viral representation and count viral load (awk, grep, sed, etc). 

#### GWAS
For performing the GWAS viral data should be represented this way:  

Group | Population | Donor ID | Mastadenovirus load
------------ | ------------- | ------------- | -------------
Case | CEU | NA12342 | 0.268727

It has been done using command line utilities (awk, grep, sed, paste, etc) - the exact commands provided in the [lab report](https://github.com/Alisa1195/Searching-for-latent-viruses-in-human-whole-genome-sequencing-data/blob/master/Lab_report_ILI_Genotek.md)

The final tables for EBV and adenoviruses were used for performing the GWAS (made by Genotek using [Plink](http://zzz.bwh.harvard.edu/plink/))

## Scripts

- __processing_script_v6.sh__ - processes alignment files and provide files with viral load and representation
- __EBV_viral_load_script.R__	- builds plots representing EBV load distribution
- __adeno_viral_load_script.R__	- builds plots representing adenoviral load distribution

all these scripts can be found [here](https://github.com/Alisa1195/Searching-for-latent-viruses-in-human-whole-genome-sequencing-data/tree/master/scripts)

#### Usage examples

Юра, сюда надо кратко примеры команд запуска скриптов

## Results
![The viral load distribution of Epstein-Barr virus in 1000 Genomes data](https://github.com/Alisa1195/Searching-for-latent-viruses-in-human-whole-genome-sequencing-data/tree/master/results/EBV_viral_load.png)

![The viral load distribution of Mastadenoviruses in 1000 Genomes data](https://github.com/Alisa1195/Searching-for-latent-viruses-in-human-whole-genome-sequencing-data/tree/master/results/adenoviruses_viral_load.png)

GWAS reports can be found [here](https://github.com/Alisa1195/Searching-for-latent-viruses-in-human-whole-genome-sequencing-data/tree/master/results/GWAS_report)


## References 
#### Database
- A global reference for human genetic variation, The 1000 Genomes Project Consortium, Nature 526, 68-74 (01 October 2015) doi:10.1038/nature15393

#### Viruses and latency
- Gelderblom, Hans R. 1996. “Structure and Classification of Viruses.” In Medical Microbiology, edited by Samuel Baron, 4th ed. Galveston (TX): University of Texas Medical Branch at Galveston. http://www.ncbi.nlm.nih.gov/books/NBK8174/.
- Grinde, Bjørn. 2013. “Herpesviruses: Latency and Reactivation – Viral Strategies and Host Response.” Journal of Oral Microbiology 5 (October). https://doi.org/10.3402/jom.v5i0.22766.
- Lieberman, Paul M. 2016. “Epigenetics and Genetics of Viral Latency.” Cell Host & Microbe 19 (5): 619–28. https://doi.org/10.1016/j.chom.2016.04.008.
- Mahy, Brian W.J. 2006. “The Diversity of Viruses Infecting Humans.” Biodiversity 7 (1): 34–37. https://doi.org/10.1080/14888386.2006.9712792.
- Moustafa A, Xie C, Kirkness E, Biggs W, Wong E, Turpaz Y, et al. (2017) The blood DNA virome in 8,000 humans. PLoS Pathog 13(3): e1006292. https://doi.org/10.1371/journal.ppat.1006292


#### Searching for viruses, virome studies, GWAS

- Paez-Espino D. et al. Uncovering Earth’s virome. Nature. 2016
- Simmonds P. et al. Virus taxonomy in the age of metagenomics. Nat Rev Microbiol. 2017.
- Edwards R.A., Rohwer F. Viral metagenomics. Nat Rev Microbiol. 2005.
- Nooij S., Schmitz D., Vennema H., Kroneman A., Koopmans M.P.G. Overview of Virus Metagenomic Classification Methods and Their Biological Applications. Front Microbiol. 2018.
- Paez-Espino D. et al. IMG/VR v.2.0: an integrated data management and analysis system for cultivated and environmental viral genomes. Nucleic Acids Res. 2019.
- Paez-Espino D., Pavlopoulos G.A., Ivanova N.N., Kyrpides N.C. Nontargeted virus sequence discovery pipeline and virus clustering for metagenomic data. Nat Protoc. 2017.
- Dehghan, Abbas. 2018. “Genome-Wide Association Studies.” Methods in Molecular Biology (Clifton, N.J.) 1793: 37–49. https://doi.org/10.1007/978-1-4939-7868-7_4.
- A.T.Marees et al, 2018. A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)
- R.Mandage et al, 2017. Genetic factors affecting EBV copy number in lymphoblastoid cell lines derived from the 1000 Genome Project samples (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5487016/)



