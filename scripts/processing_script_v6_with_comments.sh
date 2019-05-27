#! /bin/bash
usage="Welcome to the simplest ever script to analyse your data for presense of viral sequences by Kraken2 software.
This script works properly if you sequences are mapped to GRCh_38 containing decoy sequences.
Kraken2 has to be included into PATH variable.

Usage: $0 -id <path to file which contains length of corresponding reference sequences> has to be provided, included in the package
 [-i] <path to the folder with input files>
[-o] <path to the folder with output files>
[-p] <string> string to be added to the folder names (e.g. to specify the samlple more precisely)
[-db] <path to Kraken2 database> to be changed at the beginning of the script for simplicity
[-t] <number of threads to use> default is 8
[-bam] defines .bam format for input files.Ddefault is .cram
[-h] shows this help
Default values for input and output folders are ./

Good luck and have a nice day/night!

"

format_variable=cram
threads=8
input_files_path=./
output_files_path=./
path_to_database=/home/orlov239/tools/kraken2/databases/standart
path_to_length_ID_file=/home/orlov239/scripts/necessary_files/kraken_id_RefSeq_id_lenght.tsv


while [ "$1" != "" ]; do
	case $1 in
		-i | --input_directory)
			shift
			if [ "${1: -1}" != "/" ]; then
				input_files_path=$1/
			else
				input_files_path=$1
			fi
			;;
		-o | --output_directory)
			shift
			output_files_path=$1
			;;
		-p | --folder_prefix)
                        shift
                        prefix=$1
                        ;;
		-id | --id_file)
			shift
			if [ -e $1  ] && [ ! -s $1 ]; then
				echo "File $1 doesn't exist or empty. Cannot extract lengths of virus sequences to count viral load. Interrupting"
				exit 1
			else
				path_to_length_ID_file=$1
			fi
			;;
		-db | --database)
			shift
			path_to_database=$1 
			;;
		-t | --threads)
			shift
			threads=$1;
			;;
		-bam)
			shift
			format_variable=bam
			;;
		-h | --help)
			printf "$usage"
			exit
			;;
		*)
			echo "There isn't such a command $1"
			printf "$usage"
			exit 1
			;;
	esac
	shift
done

##creating output folders if they don't exist

if [ ! -d ${output_files_path}${prefix}/data_bam  ]; then
	mkdir -p ${output_files_path}${prefix}/data_bam
fi

if [ ! -d ${output_files_path}${prefix}/data_fasta  ]; then
        mkdir -p ${output_files_path}${prefix}/data_fasta
fi

if [ ! -d ${output_files_path}${prefix}/kraken_reports  ]; then
        mkdir -p ${output_files_path}${prefix}/kraken_reports
fi

if [[ "${prefix}" != "" ]]; then
	output_files_path=${output_files_path}${prefix}
fi


echo "Input directoty is $input_files_path"

echo "Output directoty is $output_files_path"

viral_load_report_file=${output_files_path}/kraken_reports/viral_load_est_general_$(date +%Y-%m-%d-%a.%H:%M:%S)
echo "Resulting output file for viral load in all samples is ${viral_load_report_file}"

##Iterating through all .cram files in the folder 

for file in $input_files_path*.${format_variable}; do
	
	##extracting filename, sample ID and read group

	filename=$(basename "${file%.*}")
	sample_ID=$(echo ${filename} | cut -d "." -f 1) ###ID of the person
	echo "working on $file"
	sample_RG=$(samtools view -H -@ ${threads} $file | grep RG | cut -f 2 | cut -d : -f 2 | head -1) ###RG in bam/cram file which corresponds to single sequencer run
        echo "Analyzing RG $sample_RG"
	
	##extracting sequences mapped to  decoy EBV and unmapped reads, then merge them and transform to paired-end fasta 

	samtools view -h -@ ${threads} -r ${sample_RG} $file chrEBV chrUn_JTFH01000690v1_decoy > ${output_files_path}/data_bam/${sample_ID}_${sample_RG}_decoy_reads.bam
	samtools view -h -@ ${threads} -r ${sample_RG} -f 12 -F 256 $file > ${output_files_path}/data_bam/${sample_ID}_${sample_RG}_unmapped_reads.bam
	
	samtools merge -@ ${threads} - ${output_files_path}/data_bam/${sample_ID}_${sample_RG}_decoy_reads.bam ${output_files_path}/data_bam/${sample_ID}_${sample_RG}_unmapped_reads.bam | samtools sort \
	-@ ${threads} -n | samtools fasta -@ ${threads} -1 ${output_files_path}/data_fasta/${sample_ID}_${sample_RG}_unmapped_reads_1.fa.gz \
        -2 ${output_files_path}/data_fasta/${sample_ID}_${sample_RG}_unmapped_reads_2.fa.gz \
        -0 /dev/null -s /dev/null -N -F 0x900 -
	
	##classifyng extracteed reads by Kraken2

	kraken2 --use-names --gzip-compressed --db ${path_to_database} --threads=${threads} --report ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_kraken_report \
	--output ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_kraken_output --confidence 0.65 \
	--paired ${output_files_path}/data_fasta/${sample_ID}_${sample_RG}_unmapped_reads_1.fa.gz ${output_files_path}/data_fasta/${sample_ID}_${sample_RG}_unmapped_reads_2.fa.gz
	cat ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_kraken_output | cut -f 3 | sort | uniq -c | sort -n -r  > ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_kraken_output_final
	cat ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_kraken_output | cut -f 3 | sort | grep -i "vir" | uniq -c | sort -n -r  > ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_kraken_output_viral
	echo "analysis for $filename has finished"
	
	##counting reads mapped to human genome

	mapped_to_human=$(samtools view -c -r ${sample_RG} -q 10 -f 3 -F 2316 -@ ${threads} $file)  ##Take into account that it's "rough" estimation as e.g. EBV reads mapped to reference are included

	viral_report_file=${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_kraken_output_viral
	echo "" > ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_virus_load_est 

	##check if Kraken 2 classification file is not empty
	##EM stands for EMpty

if [ -e ${viral_report_file} ] && [ ! -s ${viral_report_file} ]; then
	echo -e "EM\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\tNo reads were classified as viral by Kraken" >> ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_virus_load_est
	echo -e "EM\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\tNo reads were classified as viral by Kraken" >> ${viral_load_report_file} 
	continue
elif [ -s ${viral_report_file} ]; then

##iterate over lines in Kraken2 classification file

while read line; do

	##extracting information about classification
	
        seq_length=0
        info_string=$(echo "$line" | awk 'BEGIN {OFS="\t"} {read_count=$1; \
        tax_id=substr($NF,0,length($NF)-1); \
        $NF=""; $(NF-1)=""; $1=""; \
        tax_name=$0; \
        print read_count, tax_id, tax_name}')
        read_count=$(echo ${info_string} | cut -d ' ' -f 1)
        kraken_tax_id=$(echo ${info_string} | cut -d ' '  -f 2)
        tax_name=$(echo ${info_string} | cut -d ' ' -f 3-)
        seq_length=$( awk -v var=${kraken_tax_id} '($1 == var) {print $0; exit}' \
        ${path_to_length_ID_file} | cut -f 3)
        RefSeq_id=$( awk -v var=${kraken_tax_id} '($1==var) {print $0; exit}' \
        ${path_to_length_ID_file} | cut -f 2)
        
        if [[ ${seq_length} == "" ]]; then
		
		##checking if reads are classified to specie and not a higher rank
		##TE stands for Taxanomic Error

                echo "${tax_name} is not a specie so we cannot extract reference sequence to count viral load"
		echo -e "TE\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\t${kraken_tax_id}\t${tax_name} is not a specie so we cannot extract reference sequence to count viral load" >> ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_virus_load_est
		echo -e "TE\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\t${kraken_tax_id}\t${tax_name} is not a specie so we cannot extract reference sequence to count viral load" >> ${viral_load_report_file} 
                continue
        elif [[ ${seq_length} == "Manual_curation" ]]; then

		##checking if a spicie has a nonsegmented genome
		##MC stands for Multiplicity error

                echo "Ooops for ${tax_name} (Kraken id is ${kraken_tax_id}) you have to do it manually cos it has several segments"
		echo -e "MC\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\t${kraken_tax_id}\t${read_count}\tfor ${tax_name}you have to do it manually cos it has several segments" >> ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_virus_load_est
		echo -e "MC\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\t${kraken_tax_id}\t${read_count}\tfor ${tax_name}you have to do it manually cos it has several segments" >> ${viral_load_report_file}
                continue
        else

		##if everything is OK counting viral load
		##OK stands for OK
		##If you don't want to bother yourself with lots of not very usefull data, output only lines with OK in the 1st field 

                #echo "${tax_name} reference genome length is ${seq_length}"
		virus_load_est=$(awk -v vrc=${read_count} -v vsl=${seq_length} -v mth=${mapped_to_human} -v hsl=3500000000 ' BEGIN {print (2*(vrc/vsl))/(mth/hsl)}')
                echo -e "Finally we got:\n${sample_ID}\t${sample_RG}\t${mapped_to_human}\t${kraken_tax_id}\t${RefSeq_id}\t${tax_name}\t${seq_length}\t${read_count}\t${virus_load_est}"
		echo -e "OK\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\t${kraken_tax_id}\t${RefSeq_id}\t${tax_name}\t${seq_length}\t${read_count}\t${virus_load_est}" >> ${output_files_path}/kraken_reports/${sample_ID}_${sample_RG}_virus_load_est
        	echo -e "OK\t${sample_ID}\t${sample_RG}\t${mapped_to_human}\t${kraken_tax_id}\t${RefSeq_id}\t${tax_name}\t${seq_length}\t${read_count}\t${virus_load_est}" >> ${viral_load_report_file}
	fi
	done < $viral_report_file
else
	echo "Kraken classification file (has to have a name ${viral_report_file}) for ${filename} doesn't exists for For Reasons Unknown"
fi

	
done


