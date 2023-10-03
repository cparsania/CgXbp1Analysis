#!/bin/bash
#SBATCH --job-name			job
#SBATCH --partition			FHS_NORMAL
#SBATCH --nodes				1
#SBATCH --tasks-per-node		1
#SBATCH --mem				20G
#SBATCH --time				24:00:00
#SBATCH --output			job.%j.out
#SBATCH --error				job.%j.err
#SBATCH --mail-type			FAIL
#SBATCH --mail-user			lakhanp@umac.mo


bowtie2 -p 1 --trim5 8 --local  -x /home/lakhanp/database/C_glabrata/CBS138_s02-m07-r06/bowtie2_index/C_glabrata_CBS138_version_s02-m07-r06_chromosomes.fasta -U /home/hiseq/01/171101_D00691_0098_ACAYWHANXX/Data/Intensities/BaseCalls/CL1019Mix/XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_R1.fastq.gz | samtools view -bS - | samtools sort  -O bam -o XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_bt2.bam


samtools index XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_bt2.bam
samtools flagstat XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_bt2.bam > alignment.stats

mappedReads=`grep -P ' 0 mapped \(' alignment.stats | grep -P -o '^\d+'`
scale=`perl -e "printf('%.3f', 1000000/$mappedReads)"`

##macs2 pileup with 200bp extension
macs2 pileup --extsize 200 -i XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_bt2.bam -o XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_pileup.bdg
error_exit $?

##normalize
printf "Normalizing XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_pileup.bdg with factor %s\n" $scale
macs2 bdgopt -i XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_pileup.bdg -m multiply -p $scale -o temp_normalized.bdg
error_exit $?

##Remove the first line
sed -n '2,$p' temp_normalized.bdg > XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_normalized.bdg
error_exit $?
rm temp_normalized.bdg

##bedSort and bedGraph to bigWig conversion
bedSort XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_normalized.bdg XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_normalized.bdg
bedGraphToBigWig XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_normalized.bdg /home/lakhanp/database/C_glabrata/CBS138_s02-m07-r06/reference/genome.size XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_normalized.bw
error_exit $? 

##Pol-II expression value
perl /home/lakhanp/scripts/miaoScripts/zqWinSGR-v2.pl -feature_file /home/lakhanp/database/C_glabrata/CBS138_s02-m07-r06/annotation/C_glabrata_CBS138_version_s02-m07-r06_CDS_Unique.bed -socre_file XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_normalized.bdg -chrom_column 1 -start_column 2 -end_column 3  -direction_column 6 -bin_count 1 -output_folder /home/lakhanp/Analysis/19_ChIPMix_process/CL2017_CL1019Mix/XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA -outout_name XBPMyc1_THP1_2h_CL1019Mix_GTGGGATA_polii_expr.tab
error_exit $? 

gzip *.bdg

