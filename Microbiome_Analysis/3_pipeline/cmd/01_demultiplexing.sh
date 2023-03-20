#!/bin/bash
#SBATCH --qos=1day
#SBATCH --time=12:00:00
#SBATCH --mem=20g
#SBATCH --output=run.out
#SBATCH --error=run.error
#SBATCH --job-name=demultiplexing
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load modules
module load foss/2018b #interpreters
module load FastQC/0.11.8-Java-1.8
module load cutadapt/2.10-foss-2018b-Python-3.6.6
module load xlsx2csv/0.7.4-foss-2018b-Python-3.6.6

 ## --------------------------------------------------------------------
 ## Jan Wälchli | 31.12.2020
 ## --------------------------------------------------------------------

#to be set (sniplet of the sample name to differ bacteria and fungi samples)
bacteria_samples="799"
fungi_samples="ITS"

#running time notification
echo 'Start script'

#convert the design file from xlsx to tab
xlsx2csv ../../1_start/design.xlsx design.csv
cat design.csv | tr ',' '\t' | tail -n +2 > design.tab
rm design.csv

#get the runs
runs=$(awk '{print $3}' design.tab | sort | uniq)

## --------------------------------------------------------------------
## A | Quality Control - FastQC
## --------------------------------------------------------------------

#create a output folder
mkdir ../../4_output 2> /dev/null #suppress error message
rm -r  ../../4_output/qc 2> /dev/null
mkdir ../../4_output/qc

#quality control
fastqc -t 20 -k 0 -q ../../2_data/* -o ../../4_output/qc 2> /dev/null

#remove no longer needed files
rm ../../4_output/qc/*.zip

#running time notification
echo 'A - Quality Control done'


## --------------------------------------------------------------------
## B | Primer Files
## --------------------------------------------------------------------

#create a folder
rm -r  demultiplexed 2> /dev/null
mkdir demultiplexed

#create the primer files
for run in ${runs}; do

	#padding sequence to remove
	padding=$(grep ${run} design.tab | awk '{print $9}' | uniq)
	if [ ${padding} = 'F' ]; then padding=''; fi

	#forward primers
	grep ${run} design.tab | awk '{print $4, $6}' > demultiplexed/${run}_F_primers.txt
	cat demultiplexed/${run}_F_primers.txt | \
	tr ' ' '\n' | tr '_' '-' | \
	sed 's/^/\^/g' | \
	sed 's'/'^\^F'/'>'${run}'-F'/'g' | \
	sed 's'/'^'${padding}/''/'g' > demultiplexed/${run}_F_primers.fasta

	#reverse primers
	grep ${run} design.tab | awk '{print $5, $7}' > demultiplexed/${run}_R_primers.txt
	cat demultiplexed/${run}_R_primers.txt | \
	tr ' ' '\n' | tr '_' '-' | \
	sed 's'/'^R'/'>'${run}'-R'/'g' | \
	sed 's'/'^'${padding}/''/'g' > demultiplexed/${run}_R_primers.fasta


done

#remove no longer needed files
rm demultiplexed/*.txt

#running time notification
echo 'B - Primer Files done'


## ---------------------------------------------------------------------
## C | Demultiplexing
## ---------------------------------------------------------------------

for run in ${runs}; do

	mkdir demultiplexed/${run} 2> /dev/null

	cutadapt \
		-e 0.01 --no-indels \
		-g file:demultiplexed/${run}_F_primers.fasta \
		-G file:demultiplexed/${run}_R_primers.fasta \
		-o demultiplexed/${run}/{name1}-{name2}-r1.fastq \
		-p demultiplexed/${run}/{name1}-{name2}-r2.fastq \
		../../2_data/${run}_r1.fastq.gz ../../2_data/${run}_r2.fastq.gz \
		--discard-untrimmed
done

#running time notification
echo 'C - Demultiplexing done'


## ---------------------------------------------------------------------
## D | Clean up
## ---------------------------------------------------------------------

#move primer files
mkdir demultiplexed/primers
mv demultiplexed/*.fasta demultiplexed/primers

#move samples with primer combinations not in the design file

#forward primer files
awk '{print $3,$4,$3,$5}' design.tab | \
tr ' _' '-' | \
sed 's/$/-r1.fastq/g' | \
uniq > samples_to_keep.txt

#reverse primer files
awk '{print $3,$4,$3,$5}' design.tab | \
tr ' _' '-' | \
sed 's/$/-r2.fastq/g' | \
uniq >> samples_to_keep.txt

#move samples with not used primer combinations
cd demultiplexed
rm -r unsued_samples 2> /dev/null
mkdir unused_samples
for run in ${runs}; do
	cd ${run}
	for sample in *.fastq; do
		if ! grep -qxFe "${sample}" ../../samples_to_keep.txt; then
			mv ${sample} ../unused_samples
		fi
	done
	cd ..
done
cd ..
#not longer required
rm samples_to_keep.txt

#sort files by taxa
rm -r demultiplexed/bacteria demultiplexed/fungi 2> /dev/null

for run in ${runs}; do
	taxa=$(grep ${run} design.tab | awk '{print $2}' | sort | uniq)
	#bacteria
	 if [[ ${taxa} == *'b'* ]]; then
		mkdir demultiplexed/bacteria 2> /dev/null #may already exist
		mkdir demultiplexed/bacteria/${run};
		mv demultiplexed/${run}/*${bacteria_samples}*.fastq demultiplexed/bacteria/${run}
	fi
	#fungi
	 if [[ ${taxa} == *'f'* ]]; then
		mkdir demultiplexed/fungi 2> /dev/null #may already exist
		mkdir demultiplexed/fungi/${run};
		mv demultiplexed/${run}/*${fungi_samples}*.fastq demultiplexed/fungi/${run}
	fi
	rm -r demultiplexed/${run}
done

#running time notification
echo 'D - Clean up done'
echo 'End script'
