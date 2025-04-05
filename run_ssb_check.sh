#!/bin/bash

runPeriod="UL2016PreVFP"
StudyName="Testv1"
Channels="MuMu"
#Channels="ElEl"
#Channels="MuEl"
NumEvt=-1
outName="Both-test"

inputlists=("TTbar_Signal_test")
#inputlists=("TTbar_Signal_2")

configdir=""
configch=""

if [ "$Channels" = "MuMu" ]; then
    confch="dimuon.config"
elif [ "$Channels" = "ElEl" ]; then
    confch="dielec.config"
elif [ "$Channels" = "MuEl" ]; then
    confch="muelec.config"
else
    echo "Unknown Channels: $Channels"
    exit 1
fi

config="ULSummer20/${runPeriod}/${confch}"

for i in "${inputlists[@]}"
do
   mkdir -p output/${StudyName}/${runPeriod}/${Channels}
   #./ssb_analysis ${runPeriod}/${i}.list ${StudyName}/${runPeriod}/${Channels}/${i}.root analysis_config.config
   #echo ./ssb_analysis "${runPeriod}/${i}.list" "${StudyName}/${runPeriod}/${Channels}/${i}.root" "$config" "None" ${runPeriod}
   #./ssb_analysis "${i}.list" "${StudyName}/${runPeriod}/${Channels}/${i}.root" "$config" "None" ${runPeriod}
   echo ./ssb_analysis "${i}.list" "${StudyName}/${runPeriod}/${Channels}/${i}_${outName}.root" "$config" "None" ${runPeriod}
   ./ssb_analysis "${runPeriod}/${i}.list" "${StudyName}/${runPeriod}/${Channels}/${i}_${outName}.root" ${config} "None" ${runPeriod} ${NumEvt} > log_${outName}.txt
done

