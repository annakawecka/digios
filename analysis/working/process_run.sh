#!/bin/sh

if [ $# -eq 0 ] || [ $1 == "-help"  ]; then
  echo "$./process_run_simple.sh #RunNum #download #Merge #EventBld #GeneralSort #Monitor"
  echo "              RunNum = a 3 digit run number, e.g. 001"
  echo "            download = 1/0, is download from DAQ?"
  echo "               Merge = 1/0/-1, is Merge the raw data? -1 is force merge"  
  echo "            EventBld = 1/0/-1, is building event from the meged data? -1 is force build"  
  echo "          GenralSort = 2/1/0/-1, is GeneralSort?  1 = GeneralSort.C, 2 = GeneralSortTrace.C"
  echo "             Monitor = 2/1/0, run ChainMonitor.C?  1 = single run, 2 = all runs"
  exit 1
fi;

source ../../expName.sh #load expName
AnalysisDir=../../analysis

#remote data path
dataloc=/media/DIGIOSDATA3/${expName}/data
daqDir=/home/helios/digios

#===== directory and chat files
GEBDIR=$AnalysisDir/GEBSort
MERGDIR=$AnalysisDir/merged_data
ROOTDIR=$AnalysisDir/root_data
DATADIR=$AnalysisDir/data
MERGECHAT=$AnalysisDir/working/GEBMerge.chat
SORTCHAT=$AnalysisDir/working/GEBSort.chat

RUN=$1

#make RUN to be 3 digit
runLen=${#RUN}
if [ ${runLen} -eq 1 ]; then
   RUN="00"${RUN}
elif [ ${runLen} -eq 2 ]; then
   RUN="0"${RUN}
fi;

isDownload=1
isMerge=1
isSort=1
isGeneralSort=1
isMonitor=1

if [ $# -ge 2 ]; then
 isDownload=$2
fi

if [ $# -ge 3 ]; then
 isMerge=$3
fi

if [ $# -ge 4 ]; then
 isSort=$4
fi

if [ $# -ge 5 ]; then
 isGeneralSort=$5
fi

if [ $# -ge 6 ]; then
 isMonitor=$6
 OverRideMonitor=""
fi

if [ ${isGeneralSort} -eq 2 ]; then
  isMonitor=0;
  OverRideMonitor="overrided by GeneralSortTrace."
fi

echo "============================================="
echo "============ RUN $RUN ========================"
echo "============================================="

echo "Download    : ${isDownload}"
echo "Merge       : ${isMerge}"
echo "Sort        : ${isSort}"
echo "GeneralSort : ${isGeneralSort}"
echo "isMonitor   : ${isMonitor} ${OverRideMonitor}"
echo "============================================="

RED='\033[0;31m'
YELLOW='\033[1;33m'
ORANGE='\033[0;33m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
Cyan='\033[0;36m'
NC='\033[0m'

#=========== Check host name
PCName=$(hostname)
if [ ${PCName} == "digios1" ]; then
    #=== force merge and sort dor digios1, becasue the checking data time is not working
    isMerge=-1
    isSort=-1
fi

#=========== Download raw
if [ ${isDownload} -eq 1 ]; then
  echo -e "${RED}######################### Download raw data${NC}"
  if [ ${PCName} == "digios1" ]; then
      echo "Already in digios1, no need to get data."
  else
      #============ Get the raw data
      IP=192.168.1.2
      echo "RUN $RUN: Get the raw data `date`"
      rsync -avuht --progress "helios@${IP}:${dataloc}/${expName}_run_$RUN.gtd*" ${DATADIR}/.
      rsync -avuht --progress "helios@${IP}:${daqDir}/analysis/working/RunTimeStamp.dat" ${AnalysisDir}/working/.
  fi
fi

echo "============================================="
cat ${AnalysisDir}/working/RunTimeStamp.dat
echo "============================================="

du -hsc $DATADIR/${expName}_run_$RUN*

count=`ls -1 ${DATADIR}/${expName}_run_$RUN* 2>/dev/null | wc -l`
echo "========== Number of Files : ${count}"
if [ ${count} -eq 0 ]; then
    echo "============================================"
    echo "====  RAW Files of RUN-${RUN} not found! "
    echo "============================================"
    exit 
fi


#=========== Merge
if [ $isMerge -eq 1 ]; then
  echo -e "${RED}######################### Merge raw data${NC}"
  
  #==== check Merged_data timeStamp, if newer then raw data, no Merge
  rawDataTime=`stat -f "%Sm" -t "%Y%m%d%H%M%S" ${DATADIR}/${expName}_run_$RUN.gtd* | sort -rn | head -1`
  #==== check if merged data exist
  isMergedDataExist=`ls -1 ${MERGDIR}/GEBMerged_run${RUN}* 2>/dev/null | wc -l`
  if [ ${isMergedDataExist} -gt 0 ]; then
      mergedDataTime=`stat -f "%Sm" -t "%Y%m%d%H%M%S" ${MERGDIR}/GEBMerged_run${RUN}* | sort -rn | head -1`
  else
      mergedDataTime=0
  fi
  
  if [ ${rawDataTime} -ge ${mergedDataTime} ]; then
    echo "RUN $RUN: GEBMerge started at `date`"
    ${GEBDIR}/GEBMerge ${MERGECHAT}  ${MERGDIR}/GEBMerged_run${RUN}.gtd `ls ${DATADIR}/${expName}_run_$RUN.gtd*` > ${MERGDIR}/GEBMerge_run${RUN}.log
    echo "RUN $RUN: GEBMerge DONE at `date`"
  else
    echo -e "${BLUE}Merged data are newer than raw data. No need to merged again.${NC}"
    echo -e "${GREEN}You can Force merging using option -1, ${ORANGE} see ./process_run.sh -help${NC}"
  fi
fi

if [ $isMerge -eq -1 ]; then # force merge
  echo -e "${RED}######################### Merge raw data (force) ${NC}"
  echo "RUN $RUN: GEBMerge started at `date`"
  ${GEBDIR}/GEBMerge ${MERGECHAT}  ${MERGDIR}/GEBMerged_run${RUN}.gtd `ls ${DATADIR}/${expName}_run_$RUN.gtd*` > ${MERGDIR}/GEBMerge_run${RUN}.log
  echo "RUN $RUN: GEBMerge DONE at `date`"
fi

#=========== Sort
if [ $isSort -eq 1 ]; then
  echo -e "${RED}========= GEBSort started sorting run $RUN at `date`${NC}"

  #==== check Merged_data timeStamp, if newer then raw data, no Merge
  mergedDataTime=`stat -f "%Sm" -t "%Y%m%d%H%M%S" ${MERGDIR}/GEBMerged_run${RUN}* | sort -rn | head -1`
  #==== check if root data exist
  isRootDataExist=`ls -1 ${ROOTDIR}/run${RUN}.root 2>/dev/null | wc -l`
  
  if [ ${isRootDataExist} -gt 0 ]; then
    rootDataTime=`stat -f "%Sm" -t "%Y%m%d%H%M%S" ${ROOTDIR}/run${RUN}.root`
  else
    rootDataTime=0
  fi

  if [ ${mergedDataTime} -ge ${rootDataTime} ]; then
    ${GEBDIR}/GEBSort_nogeb -input disk ${MERGDIR}/GEBMerged_run${RUN}.gtd_000 -rootfile ${ROOTDIR}/run${RUN}.root RECREATE -chat ${SORTCHAT} 
    echo "GEBSort DONE at `date`"
    echo -e "========= saved root file --> ${RED} ${ROOTDIR}/run${RUN}.root ${NC} "
    echo "============================================="
    echo "============================================="
  else
    echo -e "${BLUE}Root data are newer than merged data. No need to do again.${NC}"
    echo -e "${GREEN}You can Force event buidling using option -1, ${ORANGE} see ./process_run.sh -help${NC}"
  fi
fi

if [ $isSort -eq -1 ]; then # force Sort
  echo -e "${RED}========= GEBSort started sorting run $RUN at `date` (force)${NC}"
  ${GEBDIR}/GEBSort_nogeb -input disk ${MERGDIR}/GEBMerged_run${RUN}.gtd_000 -rootfile ${ROOTDIR}/run${RUN}.root RECREATE -chat ${SORTCHAT} 
  echo "GEBSort DONE at `date`"
  echo -e "========= saved root file --> {$RED} ${ROOTDIR}/run${RUN}.root ${NC} "
  echo "============================================="
  echo "============================================="
fi

#========== Process_run.C, GeneralSort

# convert to normal number, without zero in front
if [ "${RUN:0:1}" == "0" ] ; then
      runID=${RUN:1:2}
else
      runID=$( printf '%d' $RUN )
fi;


if [ $isGeneralSort -eq 1 ]; then
  echo -e "${RED}######################### GeneralSort.C${NC}"
  root -q -b "process_run.C(${runID})"
fi

if [ $isGeneralSort -eq 2 ]; then
  echo -e "${RED}######################### GeneralSortTrace.C${NC}"
  root -q -b "process_run.C(${runID}, 1)"
fi
root -l ../Armory/runsCheck.C  #check the run Entries, and duration

#=========== Monitor
if [ $isMonitor -eq 1 ] ; then
  root -l "ChainMonitors.C(${runID})"
fi
if [ $isMonitor -eq 2 ] ; then
  root -l ChainMonitors.C
fi;  
  
exit 1
