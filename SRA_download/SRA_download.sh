#!/bin/bash
#$ -S /bin/bash
#$ -q regular.q@cbsubscb01,regular.q@cbsubscb02,regular.q@cbsubscb03,regular.q@cbsubscb04,regular.q@cbsubscb05,regular.q@cbsubscb06,regular.q@cbsubscb07,regular.q@cbsubscb08,regular.q@cbsubscb10,regular.q@cbsubscb11,regular.q@cbsubscb12
#$ -j y
#$ -cwd
#$ -pe bscb 2
#$ -t 1-14:1
#$ -N SRA_download
#$ -l h_rt=06:00:00

d1=$(date +%s)

newdir=$JOB_ID${SGE_TASK_ID}
 
echo $HOSTNAME
echo ${SGE_TASK_ID}
echo $1
echo $newdir

/programs/bin/labutils/mount_server cbsufsrv5 /data1

mkdir -p /workdir/$USER/$newdir
cd /workdir/$USER/$newdir
cp $HOME/GitHub_Repositories/Aegypti.Male.reproductive.tissue.L5/SRA_download/nonEthan.samples.SRA.list .

SRA=$(awk "NR==$SGE_TASK_ID" nonEthan.samples.SRA.list)

fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $SRA

mv fastq/* /fs/cbsufsrv5/data1/ya76/Aedes.Aegypti/Comprehensive.RNAseq/reads/SRA/
cd ..
rm -r ./$newdir


date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

