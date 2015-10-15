#!/bin/bash

if [ "$#" -lt 2 ]
then
    echo -e "Program:\tgrasper.sh (wrapper for GRASPER)"
    echo -e "Contact:\tHeewook Lee  heewlee@indiana.edu"
    echo -e "Version:\t0.1"
    echo -e "\nUsage:\tgrasper.sh <command string : I|G|A|D|S|AD|DS|ADS> <configuration file> [read1 (fastq or fq.gz)] [read2 (fastq or fastq.gz]"
    echo -e "\t\tread parameters are only needed when using \"A\" command";
    echo -e "\nCommand:\tI\tIndexing for BLAST and bwa\n\t\t\t(ONLY needs to be run once for a reference genome)"
    echo -e "\n\t\tG\tRuns pair-wise BLASTN on a given reference\n\t\t\tgenome and construct A-Bruijn graphs\n\t\t\t(ONLY needs to be run once for a reference genome)"
    echo -e "\n\t\tA\t Aligns reads to the reference by invoking bwa-MEM"
    echo -e "\n\t\tD\tDepth serialization, midpoint read sorting,\n\t\t\tobvious discordant pair removal"
    echo -e "\n\t\tS\t SV detection by invoking GRASPER"
    echo -e "\n\nWhile ( I or G ) part can only be run separately, (A, D, and S) parts can be run separtely or combined. For example, if you want to run just \"S\", issue :\n\n> grasper.sh S <config_file>\n\ncommand, if you want to run just \"D\"and \"S\" then you can do so by issueing\n\n> grasper.sh DS <config_file>\n\nor all together by running:\n\n> grasper.sh ADS <config_file> <read1> <read2>\n"
else
    . $2
    CONFIG=$2
    echo ${CONFIG}
#mkdirhier $PROJNAME
#cd $PROJNAME
    
#------------------------------------------
#1. indexing
#------------------------------------------
    if [[ "$1" == I ]]
    then
	echo "[INDEXING] ..."
	formatdb -p F -i ${REFSEQ} 
	bwa index ${REFSEQ} 
	echo -e "[INDEXING] DONE\n"
    fi
    
#------------------------------------------
#2. BLASTN and build A-Bruijn
#------------------------------------------
    if [[ "$1" = G ]]
    then
	echo -e "[Running BLASTN] ..."
	blastall -a ${CPU} -p blastn -d ${REFSEQ} -i ${REFSEQ} -o ${REFSEQ}.blast
	echo -e "\n[BLASTN] DONE\n"
	echo -e "\n[A_l-Bruijn graph construction] ..."
	${REPGRAPH_BIN}/bl2aln -i ${REFSEQ}.blast -o ${REFSEQ}.blast.aln -l ${LEG} -d ${MIN_SEQ_SIM}
	${REPGRAPH_BIN}/repeat_sin -i ${REFSEQ}.blast.aln -s ${REFSEQ} -o ${REFSEQ}.blast.aln.rg -l ${LEG}
	echo -e "\n[A_l-Bruijn graph construction] DONE\n"
    fi
    
#------------------------------------------
#3. Align via BWA
#------------------------------------------
    if [[ "$1" == *A* ]]
    then
	echo -e "[Running BWA-MEM] ...\n"
	bwa mem -M -t 1 ${REFSEQ} $3 $4 > ${SAMPE}
	echo -e "\n[BWA-MEM] DONE\n"
    fi
    
#------------------------------------------
#4. Depth serialization & mid-sorting (some discordant pair removal)
#------------------------------------------
    if [[ "$1" == *D* ]]
    then
	echo -e "\n[Loading depth information from mapped-reads and serializing depth arrays] ...\n"
	java -jar ${GRASPER_BIN}/grasper.jar depth ${REFSEQ}.thread ${SAMPE} ${MEDMADFILE} ${CONFIG} &> ${SAMPE}.depth.log
	echo -e "\n[depth processing] DONE\n"
	echo -e "\n[sorting SAM records based on midpoint and filtering obvious concordant read-pairs] ...\n"
	java -jar ${GRASPER_BIN}/grasper.jar sort ${SAMPE} ${MEDMADFILE} ${CONFIG} &> ${PROJNAME}.sort.log
	echo -e "\n[sorting] DONE\n"
    fi
    
#------------------------------------------
#5. SV Detection
#------------------------------------------
    if [[ "$1" == *S* ]]
    then
	echo -e "\n[GRASPER running to filter out more concordant read-pairs and detect SVs] ...\n"
	java -jar ${GRASPER_BIN}/grasper.jar sv ${REFSEQ}.thread ${SAMPE}.discordant.midsorted ${MEDMADFILE} ${CONFIG} ${SAMPE}.depth > ${SAMPE}.SV 2> ${SAMPE}.SV.log
	echo -e "\n[GRASPER] DONE\n"
    fi
    
#bzip2 ${SAMPE}
    
fi
