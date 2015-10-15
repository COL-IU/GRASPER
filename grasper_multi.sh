#!/bin/bash
. $1

#-------------------------------
# PREPROCESSION: generate medMAD directory and files
#-------------------------------
if [[ "$2" == P* ]]
then
#cat M275-05/M275-05_K12MG1655_Chr1_stat2.txt | head -n1 | cut -f4,6
    mkdirhier ${MEDMADDIR}
    cat ${LIST} | xargs -I{} sh -c "head -n1 ${BASE}/${JOB}/{}/{}*_stat2.txt | cut -f4,6 > ${MEDMADDIR}/{}.medMAD"
fi

#-------------------------------
# ALIGN TO REFERENCE
#-------------------------------
if [[ "$2" == *1* ]]
then
    mkdirhier ${OUTDIR}/abruijn
    cat ${LIST} | xargs -I{} -P${CPU} sh -c "cd ${OUTDIR}/abruijn;mkdir {};cd {};bwa mem -M -t 1 ${REFSEQ} ${DATADIR}/{}/{}_1.fq.gz ${DATADIR}/{}/{}_2.fq.gz > {}_PE"
fi


#-------------------------------
# UNCOMPRESS _PE.bz2
#-------------------------------
if [[ "$2" == *U* ]]
then
    cat ${LIST} | xargs -I{} -P${CPU} sh -c "cd ${OUTDIR}/abruijn/{};bunzip2 {}_PE.bz2"
fi


#-------------------------------
# DEPTH SERIALIZATION
#-------------------------------
if [[ "$2" == *2* ]]
then
    cat ${LIST} | xargs -I{} -P${CPU} sh -c "cd ${OUTDIR}/abruijn/{};java -jar ${GRASPER_BIN}/grasper.jar depth ${THREAD} {}_PE ${MEDMADDIR}/{}.medMAD ${CONFIG} &> {}.depth.log"
fi


#-------------------------------
# midsort the SAM file and also removes obvious discordant read pairs and unmapped pairs
#-------------------------------
if [[ "$2" == *3* ]]
then
    cat ${LIST} | xargs -I{} -P${CPU} sh -c "cd ${OUTDIR}/abruijn/{};java -jar ${GRASPER_BIN}/grasper.jar sort {}_PE ${MEDMADDIR}/{}.medMAD ${CONFIG} &> {}.sort.log"
fi


#-------------------------------
# SV Detection: takes midsorted discordant readpairs and .depth files, then it filters out concordant read pairs further.
#               Then, it clusters, and assign events to clusters.
#-------------------------------
if [[ "$2" == *4* ]]
then
    cat ${LIST} | xargs -I{} -P${CPU} sh -c "cd ${OUTDIR}/abruijn/{};java -jar ${GRASPER_BIN}/grasper.jar sv ${THREAD} {}_PE.discordant.midsorted ${MEDMADDIR}/{}.medMAD ${CONFIG} {}_PE.depth > {}_PE.SV 2> {}_PE.SV.log"
fi

#-------------------------------
# COMPRESS _PE
#-------------------------------
if [[ "$2" == *C* ]]
then
    cat ${LIST} | xargs -I{} -P${CPU} sh -c "cd ${OUTDIR}/abruijn/{};bzip2 {}_PE"
fi
