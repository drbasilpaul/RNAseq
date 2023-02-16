#!/bin/sh
#$ -S /bin/sh
#!/share/biosoft/java/java

##remark: use perl -p -i -e 's/\r\n$/\n/g'; and chmod 755 to make your bash executable

#definel global variable to save space
hg38_genome=/home/bpaul/ref_files/star/hg38_99
refFasta=/home/bpaul/ref_files/hg38/hg38.fa
#files=/home/bpaul/js-rnaseq859/fastq/fastq-list.txt
datadir=/rsrch3/scratch/anesth/bpaul/js-rnaseq859
outdir=/rsrch3/scratch/anesth/bpaul/js-rnaseq859/aligned
genes=/home/bpaul/ref_files/hg38/hg38.ncbiRefSeq.gtf

###step by step script

##renaming samples to identify uniq sample names
#find . -type f -name 'R859*' | while read FILE ; do
#    newfile="$(echo ${FILE} | sed -e 's/_S/-S/')" ;
#    mv "${FILE}" "${newfile}" ;

##loop
#for sample in $(find ${datadir} -type f -printf "%f\n" | sed 's/_.*//g' | sort -u)

#simple loop through all uniq prefix samples in a directory
#for sample in $(find ${outdir} -type f -maxdepth 1 -exec basename "{}" \; | cut -d'_' -f1 | sort -u)
#ls -1 /rsrch3/scratch/anesth/bpaul/js-rnaseq859 | sed 's/_.*//' | uniq

for sample in $(find ${datadir} -type f -maxdepth 1 -exec basename "{}" \; | cut -d'_' -f1 | sort -u)
#ls -1 /rsrch3/scratch/anesth/bpaul/js-rnaseq859 | uniq
do
#echo "${datadir}/${sample}_L001_R1_001.fastq.gz,${datadir}/${sample}_L002_R1_001.fastq.gz ${datadir}/${sample}_L001_R2_001.fastq.gz,${datadir}/${sample}_L002_R2_001.fastq.gz"

#read from a file using while loop
#while read samples; do
#  echo "$samples"
#done < fastq-list.txt

#copy all the files to a directory
#find "$pwd" -type f -name '*.fastq.gz' -exec sh -c '
#  for f do
#  cp $f /rsrch3/scratch/anesth/bpaul/js-rnaseq859/    

##qc steps
#Step1
#load module to hpc
#  module load fastqc
#module load multiqc

#fastqc $f -f fastq -o ${outdir}
#  echo "$f done"
#  done' sh {} +

#step2
#module load multiqc
#multiqc ${outdir}

##analysis
#step1: alignment 
module load star 
#STAR --runThreadN 124 --runMode genomeGenerate --genomeDir ${hg38_genome} --genomeFastaFiles ${refFasta} --sjdbGTFfile ${genes} --sjdbOverhang 149

#align using STAR
#STAR --runMode alignReads --runThreadN 124 --outSAMtype BAM Unsorted --quantMode GeneCounts --genomeDir ${hg38_genome} --readFilesIn ${datadir}/${sample}_L001_R1_001.fastq,${datadir}/${sample}_L002_R1_001.fastq ${datadir}/${sample}_L001_R2_001.fastq,${datadir}/${sample}_L002_R2_001.fastq --outFileNamePrefix ${outdir}/${sample}_star --sjdbGTFfile ${genes}

#echo "                  ${sample} adukku nannayi, thathwikamaya avalokanam           "

##count junctions using regtools
regtools junctions extract  ${outdir}/${sample}_starAligned.out.sam > ${outdir}/starSJ859/${sample}_junction.bed

#read counts in lincRNAs
#htseq-count -s no -f bam L${sample}_starAligned.sortedByCoord.out.bam ${gtf}gencode.v28lift37.lincRNAs.gtf > L${sample}counts.txt

#convert sam to bam
#samtools view -Sb SRR365860${sample}_starAligned.out.sam -o SRR365860${sample}_starAligned.bam

#echo "                  ${sample} converted to .bam                     "

#Sort: the bamfiles
#samtools sort -@ 24 SRR442342${sample}_starAligned.out.bam -o SRR442342${sample}_sorted.bam

#echo "                  " ${sample} sorted "                   "

#Index: the bamfiles
#samtools index SRR442342${sample}_sorted.bam

#echo "                  " ${sample} chintha saranikal "         "

#calculat coverage using bamcoverage
#bamCoverage -p 24 --ignoreDuplicates -b SRR442342${sample}_sorted.bam -bs 10 -o SRR442342${sample}_sorted.bw

#use Cufflinks: for comparison
#cufflinks -p 8 -G ${genes} -o ${outdir}/cufflinks/rd_S${sample} ${outdir}/tophat/rd_S ${sample}_tophat/accepted_hits_sorted.bam

#echo "((((((((((((((((((((gene kum))))))))))))))))))))"

done
