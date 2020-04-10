#!/bin/bash -l

#$ -pe omp 16

sh STAR.qsub samples/SRR1178008_1.fastq.gz samples/SRR1178008_2.fastq.gz SRR1178008
sh STAR.qsub samples/SRR1178009_1.fastq.gz samples/SRR1178009_2.fastq.gz SRR1178009
sh STAR.qsub samples/SRR1178010_1.fastq.gz samples/SRR1178010_2.fastq.gz SRR1178010
sh STAR.qsub samples/SRR1178014_1.fastq.gz samples/SRR1178014_2.fastq.gz SRR1178014
sh STAR.qsub samples/SRR1178021_1.fastq.gz samples/SRR1178021_2.fastq.gz SRR1178021
sh STAR.qsub samples/SRR1178047_1.fastq.gz samples/SRR1178047_2.fastq.gz SRR1178047
sh STAR.qsub samples/SRR1177981_1.fastq.gz samples/SRR1177981_2.fastq.gz SRR1177981
sh STAR.qsub samples/SRR1177982_1.fastq.gz samples/SRR1177982_2.fastq.gz SRR1177982
sh STAR.qsub samples/SRR1177983_1.fastq.gz samples/SRR1177983_2.fastq.gz SRR1177983
sh STAR.qsub samples/SRR1178050_1.fastq.gz samples/SRR1178050_2.fastq.gz SRR1178050
sh STAR.qsub samples/SRR1178061_1.fastq.gz samples/SRR1178061_2.fastq.gz SRR1178061
sh STAR.qsub samples/SRR1178063_1.fastq.gz samples/SRR1178063_2.fastq.gz SRR1178063
sh STAR.qsub samples/SRR1178004_1.fastq.gz samples/SRR1178004_2.fastq.gz SRR1178004
sh STAR.qsub samples/SRR1178006_1.fastq.gz samples/SRR1178006_2.fastq.gz SRR1178006
sh STAR.qsub samples/SRR1178013_1.fastq.gz samples/SRR1178013_2.fastq.gz SRR1178013
