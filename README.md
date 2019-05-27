# Searching for latent viruses in human whole genome sequencing data
__Students:__ Yura Orlov, Nadya Pogodina, Alisa Morshneva (Bioinformatics Institute, Saint-Petersburg)\
__Supervisors:__ Valery Ilinsky, Alexander Rakitko (Genotek, Moscow)

1.цели и задачи проекта; - __done__  
2.краткое описание используемых методов;  
3.системные требования для разработанного ПО (требования по памяти/CPU, необходимая версия операционной системы, интерпретатора, библиотек и пр.);  
4.инструкция по запуску разработанного ПО (для консольного приложения - описание ключей запуска, примеры команд с выбранными ключами);  
5.примеры получаемых с помощью ПО результатов (текст, графики, таблицы и пр.);  
6.ссылки на используемую литературу, базы данных и пр.   

## Project description
  The viral ability to stay in the asymptomatic phase of infection when a virus is not replicating is known as latency. A virus can stay latent for decades without exposing itself, although maintaining a capability to cause acute infections. Hence, viral presence in the human organism can remain unseen. 
    
  Whole-genome human sequencing data also contain sequences of dsDNA-viruses (or integrated RNA-viruses), because they are technically indistinguishable from the host DNA. These sequences stay unaligned and can be identified using databases. There are six viral families that store their genetic material in dsDNA form: Adenoviridae, Herpesviridae, Poxviridae and Polyomaviridae (linear dsDNA), Papovaviridae and Hepadnaviridae (circular dsDNA).
    
  In this project we ....

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
про подготовку данных и Plink
