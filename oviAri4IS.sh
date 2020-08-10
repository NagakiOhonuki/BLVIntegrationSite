#!/bin/bash

${oviAri4+FLK-BLV} = path to index
${oviAri4+BLV_LTR+BLV_noLTR} = path to index

# ************ Main Script Begins Here ************ #
#
		echo
		echo "Begin analysis"
		echo "Task started at"
			date
		begin=$(date +%s)

	# Prepare fastq files
		echo "Unzip"
			gunzip *gz
		echo "convert fastq file name"
			mv  *R1* R1.fastq
			mv  *R2* R2.fastq

	# Make analysis directories
		echo
		echo "Create analysis directories"
			mkdir 1_qc
			mkdir 2_map

		echo "Directories 1_qc & 2_map created"
		echo "Begin step 1 in 1_qc directory"
			cd 1_qc

	# Remove adaptor sequence
		echo
		echo "Remove adaptor from read 1"
		echo
			cutadapt -m 1 -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 5 ../R1.fastq -o R1_step1.fastq
		echo
		echo "Remove adaptor from read 2"
		echo
			cutadapt -m 1 -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 ../R2.fastq -o R2_step1.fastq
	# Read QC
		echo "Reads quality check"
		echo
			prinseq-lite.pl -min_len 10 -trim_qual_right 20 -min_qual_mean 20 -fastq R1_step1.fastq -fastq2 R2_step1.fastq -out_good clean	
	# Cleaning intermediated files
		rm R1_step1.fastq
		rm R2_step1.fastq

		echo "Proceed to step 2 in 2_map directory"
		cd ../2_map		

	# Mapping to reference genome
		echo 
		echo "Mapping using bwa"
		echo "Mapping against oviAri4+BLV"
			bwa mem -t 36 -Y -L 0 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" ${oviAri4+FLK-BLV} ../1_qc/clean_1.fastq ../1_qc/clean_2.fastq > oviAri4+virus.sam
		echo 
		echo "Mapping complete - proceed to extract reads"
		echo 
	# Reads extraction and selection
		echo "Begin reads extraction and selection"
		echo 
			samtools view -Sb oviAri4+virus.sam > oviAri4+virus.bam
		echo 
		echo "Filter and accept only first mapped reads"
			samtools view -bh -F 256 -o oviAri4+virus_uniq.bam oviAri4+virus.bam 
			samtools sort oviAri4+virus_uniq.bam -o oviAri4+virus_uniq_sort.bam

		echo "Duplicates removal"
					picard MarkDuplicates INPUT=./oviAri4+virus_uniq_sort.bam OUTPUT=./oviAri4+virus_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt REMOVE_DUPLICATES=TRUE MAX_RECORDS_IN_RAM=null ASSUME_SORTED=true TMP_DIR=./
		echo "Convert bam file to sam file"
			samtools view -h oviAri4+virus_uniq_sort_duprmv.bam > oviAri4+virus_uniq_sort_duprmv.sam
		echo 
		echo "Filter reads which map to BLV"
			awk -F"\t" '($3 ~ /EF600696.1/ || $7 ~/EF600696.1/ || $1~/^@/) {print}' oviAri4+virus_uniq_sort_duprmv.sam > virus_uniq_map_duprmv.sam
		echo 
		echo "Extract reads with soft-clipping"
				awk -F"\t" '($6 ~/S/ || $1~/^@/) {print}' virus_uniq_map_duprmv.sam > total_softclipping.sam
		echo "sam-to-bam"
			samtools view -Sb virus_uniq_map_duprmv.sam > virus_uniq_map_duprmv.bam

	# Cleaning intermediated files
			rm oviAri4+virus.sam
			rm oviAri4+virus.bam
			rm oviAri4+virus_uniq_sort.bam
			rm oviAri4+virus_uniq.bam

	# Unzip fastq file

			gzip ../R*
			
	# Count number of reads 
		echo "Analysis complete - displaying read counts"

		echo "Mapped reads"
			samtools view -F 0x4 ../2_map/oviAri4+virus_uniq_sort_duprmv.bam | cut -f 1 | sort | uniq | wc -l

		echo "Reads mapped to BLV"
			samtools view -F 0x4 ../2_map/virus_uniq_map_duprmv.bam | cut -f 1 | sort | uniq | wc -l

		echo "Soft-clipping Reads mapped to total viral reads before cleaning"
					samtools view -S -F 0x4 ../2_map/total_softclipping.sam | cut -f 1 | sort | uniq | wc -l

		"Congratulations! Data processing complete!"
		echo "Task completed on"
			date
		end=$(date +%s)
		duration=$(($end-$begin))
		echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
		echo

	# make analysis directories
		echo
			cd ..
		echo "Create analysis directories"
			mkdir 3_virus_host_junction
		echo "Directories  3_virus_host_junction created"

		echo "Proceed to mapping step in 3_virus_host_junction directory"
			cd ./3_virus_host_junction

	# Mapping to reference genome (oviAri4+BLV+LTR)
		echo
		echo "Mapping using bwa"
		echo "Mapping against oviAri4+BLV+LTR"
		bwa mem -t 36 -M -R "@RG\tID:sample\tSM:sample\tPL:Illumina" ${oviAri4+BLV_LTR+BLV_noLTR} ../1_qc/clean_1.fastq ../1_qc/clean_2.fastq > oviAri4+virus+LTR.sam
		echo
		echo "Mapping complete - proceed to extract reads"
		echo

	# Reads extraction and selection
		echo "Begin reads extraction and selection"
		echo
			samtools view -Sb oviAri4+virus+LTR.sam > oviAri4+virus+LTR.bam
		echo
		echo "Filter and accept only first mapped reads"
			samtools view -bh -F 256 -o oviAri4+virus+LTR_uniq.bam oviAri4+virus+LTR.bam
		echo
			samtools sort oviAri4+virus+LTR_uniq.bam > oviAri4+virus+LTR_uniq_sort.bam
		echo
			echo "Duplicates removal"
					picard MarkDuplicates INPUT=./oviAri4+virus+LTR_uniq_sort.bam OUTPUT=./oviAri4+virus+LTR_uniq_sort_duprmv.bam METRICS_FILE=marked_dup_metrics.txt TMP_DIR=./tmp REMOVE_DUPLICATES=TRUE MAX_RECORDS_IN_RAM=null ASSUME_SORTED=true

		echo "Convert bam to sam file"
			samtools view -h oviAri4+virus+LTR_uniq_sort_duprmv.bam > oviAri4+virus+LTR_uniq_sort_duprmv.sam
		echo
		echo "Filter reads containing BLV"
			awk -F"\t" '($3 ~ /BLV_/ || $7 ~/BLV_/ || $1~/^@/) {print}' oviAri4+virus+LTR_uniq_sort_duprmv.sam > virus+LTR_duprmv_uniq_sort.sam

	# Remove imtermediate reads 
			echo
			echo
	# making total_chimera_paired_sort.bam file
					awk -F"\t" '($3 ~ /BLV_/ && $7 ~/chr/ ||$3 ~ /chr/ && $7 ~/BLV_/ || $1~/^@/) {print}' virus+LTR_duprmv_uniq_sort.sam > total_chimera_paired.sam
			samtools view -Sb total_chimera_paired.sam > total_chimera_paired.bam

	# IS_analysis by using python
			echo "activate python2.7"
				source activate py2Bio
			echo
		echo "perform IS analysis"
			

			python2.7 ../../shell_scripts/detect_Is.py total_chimera_paired.bam ../../shell_scripts/oviAri4+BLV_LTR+noLTR_chrome_size.txt 1000 20


	# extract Total_IS.sam
		python2.7 ../../shell_scripts/grep_reads.py total_chimera_paired.bam ./total_chimera_paired.bam.out.txt ./

		echo "Congratulations! Data processing complete!"
		echo "Task completed on"
			date
		end=$(date +%s)
		duration=$(($end-$begin))
		echo "Script run time : $(($duration / 3600)) hours $((($duration % 3600) / 60)) minutes $(($duration % 60)) seconds"
		echo
#
# ************ Main Script Ends Here ************ #

