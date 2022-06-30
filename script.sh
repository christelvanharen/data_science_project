#!/bin/bash

#This code runs some fastq reads through a pipeline to map them to a reference genome and call SNPs

DIR=/home/administrator/resources/
INDEX=/home/administrator/indexes/chr10
REF=/home/administrator/ref_genome/Homo_sapiens.GRCh38.dna.chromosome.10.fa
REFCHR=/home/administrator/new_ref/Homo_sapiens.GRCh38.dna.chromosome.10_chr.fa
ANN=/home/administrator/annotations/gencode.v29.annotation_chr10.gtf
STARDIR=/home/administrator/STAR_output/
GATKDIR=/home/administrator/GATK_output/
CLINVAR=/home/administrator/all_variants/clinvar_chr.vcf
NASGATKDIR=/mnt/NAS/bin3f/GATK_output/

NASGATKDIR2=/mnt/NAS/bin3f/GATK_output2/

for FILE in $DIR*;
do
	if [[ $FILE == *_1.fastq.gz ]]
	then
		FILE1=$FILE
		PREFIX=${FILE%"chr10_1.fastq.gz"}
		FILE2=$PREFIX"chr10_2.fastq.gz"
		echo "$FILE1"
		echo "$FILE2"
		# Create STAR index
		echo "Creatin STAR index..."
		STAR --runMode genomeGenerate --genomeDir $INDEX --genomeFastaFiles $REF --sjdbGTFfile $ANN --sjdbOverhang 50 --outFileNamePrefix chr10
		# MAP the reads with STAR
		PREFIX=${PREFIX#*"/resources/"}
		STAROUTPREFIX=$STARDIR$PREFIX
		GATKOUTPREFIX=$GATKDIR$PREFIX
		echo "$STAROUTPREFIX"
		echo "$GATKOUTPREFIX"
		echo "Creating BAMs..."
		STAR --genomeDir $INDEX --readFilesIn $FILE1 $FILE2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --out>
		# Add a readgroup to BAM file
		echo "Adding readgroup to the BAM file: $PREFIX"
		/home/administrator/gatk/gatk AddOrReplaceReadGroups -I $STAROUTPREFIX"Aligned.sortedByCoord.out.bam" -O $GATKOUTPREFIX"Aligned.sortedByCoor>


		#GATK commands
		BQSR=recalibrator_data.table
		SNCRSUFF=-SNCR.bam
		RECALSUFF=-recalibrated.bam
		/home/administrator/gatk/gatk SplitNCigarReads -I $GATKOUTPREFIX"Aligned.sortedByCoord.out.bam" -R $REFCHR -O $NASGATKDIR2$PREFIX$SNCRSUFF
		/home/administrator/gatk/gatk BaseRecalibrator -I $NASGATKDIR2$PREFIX$SNCRSUFF -R $REFCHR --known-sites $CLINVAR -O $NASGATKDIR2$PREFIX"-"$B>
		/home/administrator/gatk/gatk ApplyBQSR -I $NASGATKDIR2$PREFIX$SNCRSUFF -R $REFCHR -bqsr $NASGATKDIR2$PREFIX"-"$BQSR -O $NASGATKDIR2$PREFIX$>
	fi
done


if [ "$(/bin/ls -A $GATKDIR)" ]; then
	MERGEDBAM=SamplesMergedChr10.bam
	OUTVCF=output_haplotypecaller.vcf.gz
	BQSR=recalibrator_data.table
	MERGEDBAMPRE=${MERGEDBAM%".bam"}
	SNCRSUFF=-SNCR.bam
	RECALSUFF=-recalibrated.bam
	HAPLOTYPEPRE=${OUTVCF%"haplotypecaller.vcf.gz"}
	GENOTYPESUFF=genotypegvcfs.vcf.gz
	echo "$NASGATKDIR$MERGEDBAM"
	echo "$NASGATKDIR$MERGEDBAMPRE$SNCRSUFF"
	echo "$NASGATKDIR$MERGEDBAMPRE$RECALSUFF"
	#merge BAMS
	echo "Merging the BAM files..."

	samtools merge $NASGATKDIR2$MERGEDBAM $NASGATKDIR2*$RECALSUFF
	samtools index $NASGATKDIR2$MERGEDBAM
	echo "Merging complete. Calling genotype..."
	/home/administrator/gatk/gatk HaplotypeCaller -R $REFCHR -I $NASGATKDIR2$MERGEDBAM -O $NASGATKDIR2$OUTVCF -ERC GVCF
	/home/administrator/gatk/gatk GenotypeGVCFs -R $REFCHR -V $NASGATKDIR2$OUTVCF -O $NASGATKDIR2$HAPLOTYPEPRE$GENOTYPESUFF
	echo "Calling genotype complete. You can find the results @ $NASGATKDIR2$HAPLOTYPEPRE$GENOTYPESUFF"
fi
