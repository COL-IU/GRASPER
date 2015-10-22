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
	if [ $? -eq 0 ]
	then
	    echo -e ""
	else
	    echo -e "\n[INDEXING] error in running formatdb [BLAST].\nSystem Exiting...\n"
	    exit 1
	fi
	bwa index ${REFSEQ} 
	if [ $? -eq 0 ]
	then
	    echo -e "\n[INDEXING] DONE\n"
	else
	    echo -e "\n[INDEXING] error in running bwa index.\nSystem Exiting...\n"
	    exit 1
	fi
    fi
    
#------------------------------------------
#2. BLASTN and build A-Bruijn
#------------------------------------------
    if [[ "$1" = G ]]
    then
	echo -e "[Running BLASTN] ..."
	blastall -a ${CPU} -p blastn -d ${REFSEQ} -i ${REFSEQ} -o ${REFSEQ}.blast &> ${REFSEQ}.blast.log
	if [ $? -eq 0 ]
	then
	    echo -e "\n[BLASTN] DONE\n"
	else
	    echo -e "\n[BLASTN] error in running blastall. See ${REFSEQ}.blast.log.\nSystem Exiting...\n "
	    exit 1
	fi

	echo -e "\n[A_l-Bruijn graph construction] ..."
	${REPGRAPH_BIN}/bl2aln -i ${REFSEQ}.blast -o ${REFSEQ}.blast.aln -l ${LEG} -d ${MIN_SEQ_SIM} &> ${REFSEQ}.bl2aln.log
	if [ $? -eq 0 ]
	then
	    echo -e ""
	else
	    echo -e "\n[A_l-Bruijn graph construction] error in running bl2aln in RepGraph Package. See ${REFSEQ}.bl2aln.log.\nSystem Exiting...\n"
	    exit 1
	fi
	${REPGRAPH_BIN}/repeat_sin -i ${REFSEQ}.blast.aln -s ${REFSEQ} -o ${REFSEQ}.blast.aln.rg -l ${LEG} &> ${REFSEQ}.thread.log
	if [ $? -eq 0 ]
	then
	    echo -e "\n[A_l-Bruijn graph construction] DONE\n"
	else
	    echo -e "\n[A_l-Bruijn graph construction] error in running repeat_sin in RepGraph Package. See ${REFSEQ}.thread.log.\nSystem Exiting...\n"
	    exit 1
	fi
    fi
    
#------------------------------------------
#3. Align via BWA
#------------------------------------------
    if [[ "$1" == *A* ]]
    then
	echo -e "[Running BWA-MEM] ...\n"
	bwa mem -M -t 1 ${REFSEQ} $3 $4 > ${SAMPE}
	if [ $? -eq 0 ]
	then
	    echo -e "\n[BWA-MEM] DONE\n"
	else
	    echo -e "[BWA-MEM] Error in running BWA-MEM.\nSystem Exiting...\n"
	    exit 1
	fi
    fi
    
#------------------------------------------
#4. Depth serialization & mid-sorting (some discordant pair removal)
#------------------------------------------
    if [[ "$1" == *D* ]]
    then
	echo -e "\n[Loading depth information from mapped-reads and serializing depth arrays] ...\n"
	echo -e "java -jar ${GRASPER_BIN}/grasper.jar depth ${REFSEQ}.thread ${SAMPE} ${CONFIG}\n"
	java -jar ${GRASPER_BIN}/grasper.jar depth ${REFSEQ}.thread ${SAMPE} ${CONFIG} &> ${SAMPE}.depth.log
	if [ $? -eq 0 ]
	then
	    echo -e "\n[depth processing] DONE\n"
	else
	    echo -e "[depth processing] Error in running depth command of GRASPER: See ${SAMPE}.depth.log.\nSystem Exiting...\n"
	    exit 1
	fi
	echo -e "\n[sorting SAM records based on midpoint and filtering obvious concordant read-pairs] ...\n"
	java -jar ${GRASPER_BIN}/grasper.jar sort ${SAMPE} ${CONFIG} &> ${PROJNAME}.sort.log
	if [ $? -eq 0 ]
	then
	    echo -e "\n[sorting] DONE\n"
	else
	    echo -e "[sorting] Error in running sort command of GRASPER: See ${PROJNAME}.sort.log.\nSystem Exiting...\n"
	    exit 1
	fi
    fi
    
#------------------------------------------
#5. SV Detection
#------------------------------------------
    if [[ "$1" == *S* ]]
    then
	echo -e "\n[GRASPER running to filter out more concordant read-pairs and detect SVs] ...\n"
	java -jar ${GRASPER_BIN}/grasper.jar sv ${REFSEQ}.thread ${SAMPE}.discordant.midsorted ${CONFIG} ${SAMPE}.depth > ${SAMPE}.SV 2> ${SAMPE}.SV.log
	if [ $? -eq 0 ]
	then
	    echo -e "\n[GRASPER] DONE\n"
	else
	    echo -e "[GRASPER] Error in running sv command of GRASPER: See ${PROJNAME}.SV.log.\nSystem Exiting...\n"
	    exit 1
	fi
	echo -e "\n[GRASPER] DONE\n"
    fi
    
#bzip2 ${SAMPE}
    
fi
