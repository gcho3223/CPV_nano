#!/bin/bash

runPeriod="UL2016PreVFP"
StudyName="Testv1"
Channels="MuMu"
echo "$runPeriod"

inputlists=("TTbar_Signal_1")

configdir=""
configch=""
echo "Good"

if [ "$Channels" = "MuMu" ]; then
    confch="dimuon.cfg"
elif [ "$Channels" = "ElEl" ]; then
    confch="dielec.cfg"
elif [ "$Channels" = "MuEl" ]; then
    confch="muelec.cfg"
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
   echo ./ssb_analysis "${i}.list" "${StudyName}/${runPeriod}/${Channels}/${i}.root" "ULSummer20/UL2016PreVFP/dimuon.config" "None" ${runPeriod}
   ./ssb_analysis "${runPeriod}/${i}.list" "${StudyName}/${runPeriod}/${Channels}/${i}.root" "ULSummer20/UL2016PreVFP/dimuon.config" "None" ${runPeriod} 10
done

