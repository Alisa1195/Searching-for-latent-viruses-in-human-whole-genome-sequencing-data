## Searching for latent viruses in human whole genome sequencing data

**At the start of the project, we faced a challenge of developing a
pipeline for fast viral load and representation assessment in multiple
samples. The test WGS sample (raw reads) of uterine tissue of a woman
with cervical canser provided by Genotek was used for creating and
testing pipeline. We pre-processed the data, aligned reads to the
reference, extracted unmapped reads and then compared different methods
of viruses identification (BLAST, Kraken):**

### Pre-processing

#### Deduplication + adapter trimming <br />

`~/tools/bbmap/clumpify.sh in1=~/data/HPV/raw/mv8970.82B476AA6.1.fastq.gz in2=~/data/HPV/raw/mv8970.82B476AA6.2.fastq.gz out1=F_dedup.fastq.gz out2=R_dedup.fastq.gz dedupe`

##### Summary: <br />

Reads In: 101449352 <br /> Clumps Formed: 26582955 <br /> Duplicates
Found: 5110424 <br /> Reads Out: 96338928 <br /> Bases Out: 9730231728
<br /> Total time: 2892.103 seconds <br />

TruSeq adapter trimming, deleting “N”s on the ends and discarding reads
with total “N” fraction more than 50% <br />

cutadapt <br /> --trim-n <br /> --max-n 0.5 <br /> --minimum-length 20
<br /> -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGCGATAGATCTCGTAT -A <br />
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTCAGAGCCGTGTAGATCT <br /> -o
F\_dedup\_tr1.fastq.gz <br /> -p R\_dedup\_tr1.fastq.gz <br />
F\_dedup.fastq.gz R\_dedup.fastq.gz <br />

Total read pairs processed: 48,169,464 <br /> Read 1 with adapter:
5,087,541 (10.6%) <br /> Read 2 with adapter: 4,083,942 (8.5%) <br />
Pairs that were too short: 1,083,475 (2.2%) <br /> Pairs with too many
N: 1,049 (0.0%) <br /> Pairs written <br /> (passing filters):
47,084,940 (97.7%)<br />

(AT)n sequence trimming <br />

cutadapt -minimum-length 16 -a
ATATATATATATATATATATATATATATATATATATATATATATATATAT -A
ATATATATATATATATATATATATATATATATATATATATATATATATAT <br /> -o
F\_dedup\_tr2.fastq.gz <br /> -p R\_dedup\_tr2.fastq.gz <br />
F\_dedup\_tr1.fastq.gz R\_dedup\_tr1.fastq.gz

Total read pairs processed: 47,084,940 <br /> Read 1 with adapter:
2,564,340 (5.4%) <br /> Read 2 with adapter: 2,182,251 (4.6%) <br />
Pairs that were too short: 1,451,924 (3.1%) <br /> Pairs written
(passing filters): 45,633,016 (96.9%) <br />

One more TruSeq adapter trimming and quality trimming <br />

cutadapt -minimum-length 16 -q 30,30 --pair-filter=any <br /> -a
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGCG <br /> -o F\_dedup\_tr5.fastq.gz
<br /> -p R\_dedup\_tr5.fastq.gz <br /> F\_dedup\_tr2.fastq.gz
R\_dedup\_tr2.fastq.gz <br />

Pairs written (passing filters): 40,198,889 (88.4%) <br />

**Pipeline:** <br /> Adapter trimming and quality trimming
`cutadapt --interleaved --trim-n --max-n 0.5 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGCGATAGATCTCGTAT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTCAGAGCCGTGTAGATCT input.1.fastq.gz input.2.fastq.gz | cutadapt --interleaved -m 20 -a ATATATATATATATATATATATATATATATATATATATATATATATATAT -A ATATATATATATATATATATATATATATATATATATATATATATATATAT - | cutadapt --interleaved -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGCG - | cutadapt --interleaved -m 20 -q 30,30 -o F_output.fastq.gz -p R_output.fastq.gz -`

`/home/orlov239/tools/bbmap/clumpify.sh \ in1=F_output.fastq.gz \ in2=R_output.fastq.gz \ out1=F_dedup_output.fastq.gz \ out2=R_dedup_output.fastq.gz dedupe`

### Alignment of the WGS data to the reference

at first we chose **GRCh38 (2013)** as a reference <br />

Then, there are different strategies: <br /> 1) use GRCh37 - less full
and doesn't include EBV. <br /> One the one hand, we don't need the
precise accuracy of the assembly - we are interested in another part of
the data - unmapped reads. One the other hand, the less reads are left
after alignment, the faster and easier we perform the further steps.
<br /> 2) use one of GRCh38 releases and then extract mapped to EBV
reads <br /> 3) “make” our own reference by exctraction of EBV out of
GRCh38 assembly

#### Indexing of the reference

`bwa index /home/g1195alisa/GRCh38_2013/GCA_000001405.15_GRCh38_genomic.fna`

#### Alignment of paired-end reads to the reference

`bwa mem -t 6 /home/g1195alisa/GRCh38_2013/GCA_000001405.15_GRCh38_genomic.fna /home/orlov239/data/HPV/processed/F_dedup_tr5.fastq.gz /home/orlov239/data/HPV/processed/R_dedup_tr5.fastq.gz`

Alignment took almost 4 hours to be finished (6 threads), resulting in
sam-файл of 25Gb

#### Conversion to bam

`samtools view -S -b alignment_test_3.sam > alignment_test_3.bam`

file has shrinked to 8Gb <br />

**The test alignment statistics** <br />

samtools flagstat alignment\_test\_3.bam <br />

output: <br /> 83078254 + 0 in total (QC-passed reads + QC-failed
reads)<br /> 0 + 0 secondary<br /> 2680476 + 0 supplementary<br /> 0 + 0
duplicates<br /> 80958206 + 0 mapped (97.45% : N/A)<br /> 80397778 + 0
paired in sequencing<br /> 40198889 + 0 read1<br /> 40198889 + 0
read2<br /> 68079812 + 0 properly paired (84.68% : N/A)<br /> 77566460 +
0 with itself and mate mapped<br /> 711270 + 0 singletons (0.88% :
N/A)<br /> 7094860 + 0 with mate mapped to a different chr <br />
4864885 + 0 with mate mapped to a different chr (mapQ&gt;=5) <br />

The percentage of mapped reads is suspiciously high, almost 100%, but if
we look at amount of reads, about 2 millions of reads have been left
unmapped

#### Extraction of unmapped reads

`samtools view -b -f 4 alignment_test_3.bam > unmapped_test.bam (single-end)`

To get paired-end: <br />

**R1 & R2 unmapped** <br />

`samtools view -u -f 12 -F 256 alignment_test_sorted.bam > f1_unmap_unmap.bam`

#### Resulting pipeline:

1.  align + convert to bam + sort <br />

`bwa mem -t 16 /home/g1195alisa/GRCh38_2013/GCA_000001405.15_GRCh38_genomic.fna /home/orlov239/data/HPV/processed/F_dedup_tr5.fastq.gz /home/orlov239/data/HPV/processed/R_dedup_tr5.fastq.gz | samtools view -S -@ 16 -b | samtools sort -@ 16 -o alignment_sorted.bam.gz`

1.  extract unmapped reads <br />

`samtools view -u -f 12 -F 256 alignment_test_sorted.bam > f1_unmap_unmap.bam`

1.  Convert bam to fastq <br /> 

`samtools bam2fq f1_unmap_unmap.bam > unmap_unmap.fastq`

*Summary of the alignment step: the speed needs improvement, some
commands should be combined for convenience.* <br />

### De novo assembly

The study of Moustafa et. al, 2017 compares results of the BLAST and de
novo asembly. We decided to skip the assembly, because it takes a lot of
time, considering the fact we needed to analyse hundreds of genomes

`spades.py --meta --only-assembler -12 /input.fastq -o ./output_folder/`
“scaffolds.fasta”

### Alignment of reads in blast.

We've got unmapped reads (fastq). Next, we should convert them to fasta <br />

`$ paste - - - - < unmapped_test.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > /home/nadya379/unmapped_reads.fasta`
The other convertation method <br />
<https://www.ecseq.com/support/ngs-snippets/convert-fastq-to-fasta-on-the-command-line>
<br />

#### Downloading of the reference, making the database

Link to the reference: NCBI <ftp://ftp.ncbi.nlm.nih.gov/blast/db/>
<br /> For viruses: ref\_viruses\_rep\_genomes. Plus,
ref\_prok\_rep\_genomes and vector should be downloaded. <br /> To get
the result with taxonomical names for each read we need to download the
database and specify the path: <br /> `$ ./update_blastdb.pl taxdb`
`$ mv taxdb* ../blastdb/` `$ cd blastdb`
`$ tar -xzvf taxdb.tar.gz  (распаковка сразу нескольких архивов cat *.tar.gz | tar -xzvf - -i)`
`$ rm taxdb.tar.gz`
`$ export BLASTDB=$BLASTDB:/home/nadya379/ncbi-blast-2.8.1+/blastdb/`

Combine three databases <br />
`$ blastdb_aliastool -dblist "ref_viruses_rep_genomes ref_prok_rep_genomes vector" -dbtype nucl \ -out vir_prok_vec -title "Viruses, prok and vectors"`

Launch BLAST <br /> Info about some flags in blastn
<https://www.ncbi.nlm.nih.gov/books/NBK279684/> Launched the alignment
`$ blastn -query ../unmapped_reads.fasta -db blastdb/vir_prok_vec -out ../align_to_viral_ref -outfmt '6 qseqid sseqid pident nident length mismatch evalue bitscore sscinames sskingdoms’ -evalue 1e-10`

Specify the E-value threshold. Used the one from the article (Moustafa,
2017). If use default 10, get a lot of “0 hits found” <br />

`$ cut -f 10 ../align_to_vir_prok_vec | sort | uniq -c` 2453
Archaea<br /> 3061153 Bacteria <br /> 1983107 N/A <br /> 29028 Viruses
<br />

Extract only viruses, then extract only unique reads and make the
database
`$ grep Viruses ../align_to_vir_prok_vec | sort -u -k1,1 | cut -f 1,11 | awk '{printf(">%s\n%s\n", $1, $2); }' > ../viruses_after_first_align.fasta`

We've got 11620 reads <br /> Get read names we need:
`$paste - - - - < ../viruses_after_first_align.fasta | cut -f 1 | sed 's/^@/>/' | tr "\t" "\n" > ../reads_viruses`

Remove unnecessary characters
`$cat ../reads_viruses.txt | cut -c2- > ../reads_viruses`

Search for these reads in fasta file
`$perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ../reads_viruses ../unmapped_reads.fasta > ../viruses_after_align.fasta`

Make a dataframe for the further filtration
`$ cut -f 1,8,9,10 ../align_to_vir_prok_vec > ../for_R.csv`

Repeat in R, using dplyr.The code: `ССЫЛКА НА ГИТ`

**The point is, we group by reads, then search for the maximum bitscore
in groups and get unique reads, among which we search for the viruses.**
<br /> We've got 10164 reads, weach is less than after a direct search
of all mapped to viruses reads. <br />

We could try to align these reads to rafseq\_genomic, to exclude false
alignments
`$ blastn -query ../viruses_after_align.fasta -db blastdb/refseq_genomic -out ../align_to_all_ref -outfmt '6 qseqid sseqid pident nident length mismatch evalue bitscore sscinames sskingdoms' -evalue 1e-20 -num_threads 8`
Doesn't work... <br />

Try Blat (installation
<http://genomic-identity.wikidot.com/install-blat>) <br /> Blat takes
fasta files. Download: <br />
<https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/> <br /> Merge to
one file: `$ cat *.fna >> all_ref_seq.fasta` Blat doesn't work with this
database either. We decided to use only the first database (prokaryotes,
vectors and viruses)

### Launch of kraken2

<https://ccb.jhu.edu/software/kraken2/>

Used a standart RefSeq database (~40 Gb):

`kraken2-build --standard --threads 7 --db ./databases/standart`

Analysed our data:
`kraken2 --fastq-input --use-names --db ./databases/standart --threads 7 --output ~/data/HPV/reports/kraken2/1st_try/kraken_output --classified-out ~/data/HPV/reports/kraken2/1st_try/kraken_output_clssified_seq.fastq --unclassified-out ~/data/HPV/reports/kraken2/1st_try/kraken_output_unclssified_seq.fastq /home/unmapped_reads/unmapped_test.fastq`

Report: 2120048 sequences (173.79 Mbp) processed in 5.344s (23802.4
Kseq/m, 1951.19 Mbp/m). <br /> 306588 sequences classified (14.46%)
<br /> 1813460 sequences unclassified (85.54%) <br />

Out of 306588 classified sequences about 12000 belong to viral ones
(parsing of kraken\_output file with awk and grep (“vir”), 7000 (awk and
grep either) belong to HPV.

The same for assembled metagenome:
`kraken2 --use-names --db ./databases/standart --threads 7 --output ~/data/HPV/reports/kraken2/1st_try/kraken_output_assembled --classified-out ~/data/HPV/reports/kraken2/1st_try/kraken_output_assembled_classified_seq.fa --unclassified-out ~/data/HPV/reports/kraken2/1st_try/kraken_output_assembled_unclassified_seq.fa /home/orlov239/data/HPV/assembled/1st_attempt/scaffolds.fasta`

9738 sequences (2.97 Mbp) processed in 0.607s (963.2 Kseq/m, 293.92
Mbp/m).<br /> 6230 sequences classified (63.98%)<br /> 3508 sequences
unclassified (36.02%)<br />

10 viral sequences:

4 на HPV type 16<br /> 6 на Human mastadenovirus C<br />

8000 Viromes
------------

tested pipeline at the data from the article Moustafa et.al

**BLAST** `$ cat *.fa > /home/nadya379/all_virome_queries.fa`
`$ blastn -query ../all_virome_queries.fa -db blastdb/vir_prok_vec -out ../Moustafa_output -outfmt '6 qseqid sseqid pident nident length mismatch evalue bitscore sscinames sskingdoms' -evalue 1e-20 -num_threads 8`
`$ cut -f 10 ../Moustafa_output | sort | uniq -c`

1 Archaea <br /> 141 Bacteria <br /> 18 N/A <br /> 6507 Viruses <br />

**Kraken**

<https://genomics.sschmeier.com/ngs-taxonomic-investigation/index.html>

`$ kraken2 --fasta-input --use-names --db ./databases/standart --threads=8 --output /home/g1195alisa/kraken_reports/kraken_output_all --classified-out /home/g1195alisa/kraken_reports/kraken_output_clssified_all.fa --unclassified-out /home/g1195alisa/kraken_reports/kraken_output_unclssified_all.fa --report /home/g1195alisa/kraken_reports/kraken_output_report_all  /home/nadya379/Moustafa/all_virome_queries.fa`

Output files:<br /> kraken\_output\_all<br />
kraken\_output\_clssified\_all<br />
kraken\_output\_unclssified\_all<br /> kraken\_output\_report\_all

**Taxonomic classification in the "report" file
(kraken\_output\_report\_all ):**

1.  Percentage of reads covered by the clade rooted at this taxon<br />
2.  Number of reads covered by the clade rooted at this taxon<br />
3.  Number of reads assigned directly to this taxon<br />
4.  A taxonomy rank code<br />
5.  NCBI taxonomy ID<br />
6.  indented scientific name<br />

We could analyse the representation of different taxones:
`$ cut -f 3 ./kraken_output_all_1 | sort | uniq -c | sort -n -r`

The figures at the end of lines represent what place the virus takes in
original article <br /> 1167 Human gammaherpesvirus 4 (taxid 10376) 2
(1190/8000)\* <br /> 1082 Rhizobium phage RR1-B (taxid 929834)<br /> 62
Human betaherpesvirus 7 (taxid 10372) 1 (1678/8000)<br /> 48
unclassified (taxid 0)<br /> 43 root (taxid 1)<br /> 31 Brucella sp.
09RB8471 (taxid 1149952) <br /> 21 Burkholderia pseudomallei (taxid
28450) ? <br /> 13 Burkholderia virus Bcep22 (taxid 242527) ? <br /> 11
Ralstonia pickettii 12D (taxid 428406) <br /> 10 Torque teno virus 12
(taxid 687351) 3 <br /> 10 Alphatorquevirus (taxid 687331) <br /> 7
Pa6virus (taxid 1982251) <br /> 6 Shinella sp. HZN7 (taxid 879274)
<br /> 6 Ralstonia solanacearum Rs-10-244 (taxid 1457195) <br /> 5
Torque teno virus 2 (taxid 687341) 3 <br /> 5 Torque teno virus 15
(taxid 687354) 3 <br /> 5 Torque teno virus 13 (taxid 687352) 3 <br /> 5
Human herpesvirus 4 type 2 (taxid 12509) 2 <br /> 5 Chromobacterium
rhizoryzae (taxid 1778675) <br /> 5 blood disease bacterium A2-HR MARDI
(taxid 1944648) <br /> 4 Human betaherpesvirus 6B (taxid 32604) 4
(395/8000) <br /> 3 Torque teno virus 3 (taxid 687342) 3 <br /> 3
Ralstonia insidiosa (taxid 190721) <br /> 3 Brucella sp. 141012304
(taxid 1885919) <br /> 2 Sinorhizobium phage PBC5 (taxid 179237) <br />
2 Sindbis virus (taxid 11034) <br /> 2 Pseudomonas (taxid 286) <br /> 2
Propionibacterium phage PHL010M04 (taxid 1235645) <br /> 1 unclassified
Enterobacteriaceae (miscellaneous) (taxid 36866) <br /> 1 Torque teno
virus 6 (taxid 687345) 3 <br /> 1 Torque teno virus 4 (taxid 687343) 3
<br /> 1 Torque teno virus 10 (taxid 687349) 3 <br /> 1 Streptomyces sp.
CMB-StM0423 (taxid 2059884) <br /> 1 Rhizobium sp. IRBG74 (taxid 424182)
<br /> 1 Ralstonia (taxid 48736) <br /> 1 Ralstonia solanacearum (taxid
305) <br /> 1 Ralstonia solanacearum FJAT-1458 (taxid 1130828) <br /> 1
Propionibacterium phage QueenBey (taxid 1654782) <br /> 1
Propionibacterium phage PHL171M01 (taxid 1500827) <br /> 1
Propionibacterium phage PHL025M00 (taxid 1500799) <br /> 1
Propionibacterium phage P104A (taxid 1229787) <br /> 1 Propionibacterium
phage BruceLethal (taxid 1654740)<br /> 1 Moloney murine leukemia virus
(taxid 11801)<br /> 1 Escherichia coli (taxid 562) <br /> 1
Altererythrobacter namhicola (taxid 645517) <br /> 1 Alphaproteobacteria
(taxid 28211)<br />

\*We've got Human gammaherpesvirus 4 at the top, but it's a high-copy
virus, so it's just a smal selection effect <br />

#### Tasks:

1.  Find out how to analyse a lot of genomes in a time: made the bash
    script <br />
2.  Get report zero -use-map-style statistics - done <br />
3.  Ectract only viruses (SRR) - done <br />
4.  Split reads by the launch ID (RG) - done <br />

### Testing of WGS samples from Genotek

extraction of unmapped reads
`samtools view -f 4 ./raw_data/zd3089_clipped.bam > ./zd3089_unmap.bam`
conversion to fastq
`for f in *_unmap.bam; do filename="${f%%.*}"; samtools bam2fq -@ 6 $f > ${fil ename}.fastq; done`
Kraken2
`for f in *_unmap.fastq; do filename="${f%%.*}"; kraken2 --use-names --db /home/orlov239/tools/kraken2/databases/standart --threads=8 --report /home/g1195alisa/genotek_wgs/kraken_reports/${filename}_report  --output /home/g1195alisa/kraken_reports/${filename} $f; done`

**Suspicious results: lots of exotic viruses, while cosmopolitan ones
are not found: need to use some threshold - Kraken has confidence score
option - set 0,65.**

No viruses found in WES samples, which was expected.<br />

### Viral load estimation

counted for HPV<br /> Used fies: `GCA_000001405.15_GRCh38_genomic.fna`
`HPV16_alphaPV9_RefSeq.fna (RefSeq NC_001526, ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Alphapapillomavirus_9/latest_assembly_versions/GCF_000863945.3_ViralProj15505)`
`alignment_sorted.bam.gz` `F_dedup_tr5.fastq.gz` `R_dedup_tr5.fastq.gz`
Формула для подсчета вирусного титра (количество копий вирусного генома
на один клеточный геном): <br />

C=2\*(number of reads mapped to virus genome/virus genom size)/(number
of reads mapped to human genome/human genome size) <br /> или<br />

C=2\*(virus genome coverage/human genome coverage)<br />

**Method 1.** Assess the number of aligned (to human or HPV) reads from
pre-processed fastq <br /> A) Count amount of reads correctly aligned to
human
`samtools view -c -q 10 -f 3 -F 2316 -@ 8 /home/g1195alisa/alignment/alignment_sorted.bam.gz`

-q 10 (MAPQ&gt;10) <br /> -f 3 (read paired (0x1) and read mapped in
proper pair (0x2)) <br /> -F 2316 (read unmapped (0x4), mate unmapped
(0x8), not primary alignment (0x100) and supplementary alignment
(0x800)) <br /> 59,489,753.00 <br />

1.  Align primary fastq files to HPV and count aligned reads
    (MAPQ&gt;10)

`bwa mem -t 8 /home/orlov239/data/reference/viral_RefSeq/HPV16_alphaPV9_RefSeq.fna /home/orlov239/data/HPV/processed/F_dedup_tr5.fastq.gz`
`/home/orlov239/data/HPV/processed/R_dedup_tr5.fastq.gz | samtools view -@ 8 -b | samtools sort -@ 8 -o /home/orlov239/data/HPV/aligned_to_HPV/aligned_to_HPV16.bam`

`samtools view -c -q 10 -f 3 -F 2316 -@ 8 /home/orlov239/data/HPV/aligned_to_HPV/aligned_to_HPV16.bam`
35,629.00

C= 523.72 <br />

**Method 2.** Align unmapped reads to HPV<br />

1.  Extract unmapped single-end reads from the alignment file
    `samtools fastq -@ 8 -f 4 -F 0x900 /home/g1195alisa/alignment/alignment_sorted.bam.gz > /home/orlov239/data/HPV/unmapped/unmapped_to_GRCh38.fastq`

2.  Align them to HPV and count aligned reads
    `bwa mem -t 8 /home/orlov239/data/reference/viral_RefSeq/HPV16_alphaPV9_RefSeq.fna /home/orlov239/data/HPV/unmapped/unmapped_to_GRCh38.fastq | samtools view -@ 10 -b | samtools sort -@ 10 -o`
    `/home/orlov239/data/HPV/aligned_to_HPV/aligned_to_HPV16_from_unmapped_to_GRCh38.bam`

`samtools view -c -q 10 -F 2308 -@ 8 /home/orlov239/data/HPV/aligned_to_HPV/aligned_to_HPV16_from_unmapped_to_GRCh38.bam`
-F 2308 (read unmapped (0x4), not primary alignment (0x100) and
supplementary alignment (0x800))<br /> 7,332.00 <br />

C= 107.77 <br />

**Method 3.** Use number of reads, classified as HPV by Kraken2 <br />

А) Process unmapped reads with Kraken
`/home/orlov239//tools/kraken2/kraken2 --fastq-input --use-names --db ./databases/standart --threads 7 --output ~/data/HPV/reports/kraken2/1st_try/kraken_output --classified-out ~/data/HPV/reports/kraken2/1st_try/kraken_output_clssified_seq.fastq --unclassified-out ~/data/HPV/reports/kraken2/1st_try/kraken_output_unclssified_seq.fastq /home/orlov239/data/HPV/unmapped/unmapped_to_GRCh38.fastq`

В) Count HPV reads
`awk '$1=="C" {print $0}' /home/orlov239/~/data/HPV/reports/kraken2/1st_try/kraken_output | grep -i papillo | cut -f 3 | sort | uniq -c | head`

7,086.00<br />

C= 104.16 <br />

**Methods 4 and 5.** Count viral load using mean and median genome
coverage<br />

Choose random parts of human genome for counting the coverage <br />

Buid index of the .bam file (samtools index)<br /> Create file
GRCh38\_intervals\_4x2000x22.bed, containing intervals' lengths and
their localization in chromosomes:<br /> chr\_name interval\_start
interval\_end<br /> Only the total length maters<br />

Choose random intervals:
`bedtools shuffle -chrom -noOverlapping -i GRCh38_int_rand_4x2000x22_sort.bed -g ~/data/reference/GRCh38/GRCh38_chrName_chrLen > GRCh38_int_rand_4x2000x22.bed`
Sorting:
`sort -k1,1 -k2,2n GRCh38_int_rand_4x2000x22.bed > GRCh38_int_rand_4x2000x22_sort.bed`

Count median coverage of the human genome
`bedtools coverage -hist -sorted -a GRCh38_int_rand_4x2000x22_sort.bed -b /home/g1195alisa/alignment/alignment_sorted.bam -g ~/data/reference/GRCh38/GRCh38_chrName_chrLen | grep "all" | awk 'BEGIN {sum=0}; {sum+=$3}; (sum>=$4/2) {print $2; exit}; END {print sum}`
1<br />

Count mean coverage of the human genome
`samtools depth -a /home/g1195alisa/alignment/alignment_sorted.bam.gz -b GRCh38_intervals.bed | awk 'BEGIN {sum=0; num=0}; {sum+=$3; num+=1}; END {print sum/num}'`
2.10953<br />

Count median coverage of the viral genome
`bedtools genomecov -ibam /home/orlov239/data/HPV/aligned_to_HPV/aligned_to_HPV16_from_unmapped_to_GRCh38.bam | grep "genome" | awk 'BEGIN {sum=0}; {sum+=$3}; (sum>=4953) {print $2; exit}; END {print sum}'`
94<br />

Count mean coverage of the human genome
`bedtools genomecov -d -ibam ~/data/HPV/aligned_to_HPV/aligned_to_HPV16_from_unmapped_to_GRCh38.bam | awk 'BEGIN {sum=0; num=0}; {sum+=$3; num+=1}; END {print sum/num}'`
90.5<br />

Cmed=188 <br />

Cavr=89.12 <br />

**All methods except for the first one have shown similar results**
<br />

We use 3rd method as the simplest one

\_\_So, we've chosen the method to count the viral load, fine-tuned the
Kraken2 parameters, parsed report files to extract the information
needed. Next, we needed to apply our pipeline to several european
populations from 1000 Genomes and provide some statistics for counting
viral lad and representation. Each population involved about a hundred
of samples, so we needed the united cycle, to apply all command for all
sample, choosing from every sample only the reads from the one launch.
We've developed the bash-cycle, the code is also provided in this repo.

1000 Genomes:
-------------

#### 1. Extract reads mapped to decoy sequences

`samtools view -h -@ 8 /home/nadya379/1000_genomes_cram/data_cram/CEU/*.low_coverage.cram chrEBV chrUn_JTFH01000690v1_decoy > /home/nadya379/1000_genomes_cram/data_unmapped_bam/*_decoy.bam`

#### 2. Extract paired unmapped reads

`samtools view -h -f 12 -F 256 -@ 8 /home/nadya379/1000_genomes_cram/data_cram/CEU/*low_coverage.cram > /home/nadya379/1000_genomes_cram/data_unmapped_bam/*_unmapped.bam`

#### 3. Merge two files (1 and 2)

`samtools merge /home/nadya379/1000_genomes_cram/data_unmapped_bam/*.bam /home/nadya379/1000_genomes_cram/data_unmapped_bam/*_decoy.bam /home/nadya379/1000_genomes_cram/data_unmapped_bam/*_unmapped.bam -@ 8`

#### 4. Convert to paired fasta

`samtools fasta -1 /home/nadya379/1000_genomes_cram/data_paired_fasta/*_1.fasta -2 /home/nadya379/1000_genomes_cram/data_paired_fasta/*_2.fasta -@ 8 /home/nadya379/1000_genomes_cram/data_unmapped_bam/*.bam`

#### 5. Launch Kraken2 in paired-end mode

`/home/orlov239/tools/kraken2/kraken2 —use-names —db /home/orlov239/tools/kraken2/databases/standart —threads=8 —report /home/nadya379/1000_genomes_cram/kraken_output_report/*_kraken_report —output /home/nadya379/1000_genomes_cram/kraken_output_report/*_kraken_output —confidence 0.65 —paired /home/nadya379/1000_genomes_cram/data_paired_fasta/*_1.fasta /home/nadya379/1000_genomes_cram/data_paired_fasta/*_2.fasta`

Bash-code for applying the pipeline for all samples:<br/> **should be
launched from the directory with cram-files (!!!)**

**UPD:** current version of the bash-script *ссылка на гитхаб*

### GWAS

#### Preparing the data for GWAS

vcf-files with SNPs for 1000 Genomes were downloaded from .... <br/>

Using the Kraken2 report files we've made the table with the following
columns: <br/> 1) sample ID <br/> 2) population <br/> 3) EBV reads\*
<br/> 4) EBV load\* <br/>

\*adenovirus reads & adenovirus load for adenovirus

We've made separate files fr each population and then merged them:

`cat /home/1000_genomes/kraken_output_report/IBS/viral_load_est_general_2019-05-09-Thu.07\:35\:03 | cut -f 2 | uniq | sed 's/^/ IBS\t /' > ids_IBS > TSI_IDs`

`sudo cat viral_load_est_general_2019-05-07-Tue.21\:48\:22  | grep "Human gammaherpesvirus 4" | cut -f 2,9,10 | sed 's/^/ TSI\t /' > TSI_EBV`

`paste TSI_IDs TSI_EBV > EBV_TSI_report`

`cat EBV_GBR_report EBV_IBS_report EBV_FIN_report EBV_TSI_report_new EBV_TSI_report > EBV_report`

GWAS was perdormed by Genotek using the soft **Plink** <br/>

**Kraken and GWAS reports, plots and scripts can be found in this
repository**
