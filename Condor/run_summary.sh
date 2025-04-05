#!/bin/bash

samplelist=( 
    ### Data & MC Samples ###
    #"Data_DoubleMuon_Run2016B" 
    #"Data_DoubleMuon_Run2016C" 
    #"Data_DoubleMuon_Run2016D" 
    #"Data_DoubleMuon_Run2016E" 
    #"Data_DoubleMuon_Run2016F" 
    #"Data_DoubleMuon_Run2016G" 
    #"Data_DoubleMuon_Run2016HV2" 
    #"Data_DoubleMuon_Run2016HV3" 
    #"Data_SingleMuon_Run2016B" 
    #"Data_SingleMuon_Run2016C" 
    #"Data_SingleMuon_Run2016D" 
    #"Data_SingleMuon_Run2016E" 
    #"Data_SingleMuon_Run2016F" 
    #"Data_SingleMuon_Run2016G" 
    #"Data_SingleMuon_Run2016HV2" 
    #"Data_SingleMuon_Run2016HV3" 
    #"DYJetsToLL_M_10To50" 
    #"DYJetsToLL_M_50" 
    #"ST_tW_antitop" 
    #"ST_tW_top" 
    #"TTbar_WJetToLNu" 
    #"TTbar_WQQ" 
    #"TTbar_ZQQ" 
    #"TTbar_ZToLLNuNu" 
    #"TTJets_others" 
    "TTbar_Signal" 
    #"WJetsToLNu" 
    #"WW" 
    #"WZ" 
    #"ZZ" 
    ### CPV samples ### 
    #"TTJets_Signal_dtG_0" 
    #"TTJets_Signal_dtG_0p5207" 
    #"TTJets_Signal_dtG_1p0415" 
    #"TTJets_Signal_dtG_2p60364" 
    #"TTJets_Signal_dtG_m0p5207" 
    #"TTJets_Signal_dtG_m1p0415" 
    #"TTJets_Signal_dtG_m2p60364" 
)
######################################################################
### before running, should change samples you want to make summary ###
### you can use CPV samples or Dataset samples at the same time    ###
######################################################################
###########################
### version information ###
###########################
version="SelLep"
runPeriod="UL2016PreVFP"
Channels="MuMu"

#########################################################
### set up sample_path                                ###
### calculate total # of log files                    ###
### set up log file path & variables for progress bar ###
#########################################################
### progress bar variables ###
total_sample=0
completed=0
### set up sample_path & calculate total # of log files ###
for i in "${samplelist[@]}"
do
    # CPV 샘플인지 확인
    if [[ $i == TTJets_Signal_dtG_* ]]; then
        sample_path="./Job_Version/${version}/CPV_Sample/${runPeriod}/${Channels}"
    else
        sample_path="./Job_Version/${version}/Dataset/${runPeriod}/${Channels}"
    fi
    # 로그 파일 경로 설정 및 파일 개수 계산
    log_path="${sample_path}/${i}/log_condor/out"
    file_count=$(ls "$log_path"/*.out 2>/dev/null | wc -l)
    total_sample=$((total_sample + file_count))  #sum of total log files
done
# total_sample가 0인지 확인
if (( total_sample == 0 )); then
    echo "No log files found. Exiting."
    echo "Please check the log files in ${sample_path} directory."
    exit 1
fi
#########################################
###        start summary !!           ###
#########################################
for i in "${samplelist[@]}"
do
    ################################################
    ###   check sample name & set up make dir    ###
    ################################################
    # CPV 샘플인지 확인
    if [[ $i == TTJets_Signal_dtG_* ]]; then
        load_path="./Job_Version/${version}/CPV_Sample/${runPeriod}/${Channels}/${i}"
        sample_path="./${version}/CPV_Sample/${runPeriod}/${Channels}/${i}"
    else
        load_path="./Job_Version/${version}/Dataset/${runPeriod}/${Channels}/${i}"
        sample_path="./${version}/Dataset/${runPeriod}/${Channels}/${i}"
    fi
    # 디렉토리 및 파일 생성
    mkdir -p ./Run_Summary/${sample_path}/${i}
    > ./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt
    > ./Run_Summary/${sample_path}/${i}/Memory_Summary_${i}.txt
    ## define log, err, out path ##
    log_path="${load_path}/${i}/log_condor/log"
    out_path="${load_path}/${i}/log_condor/out"
    err_path="${load_path}/${i}/log_condor/err"
    ##########################################  
    ###       calc # of log files         ###
    ##########################################
    file_count=$(ls "$out_path"/*.out 2>/dev/null | wc -l)
    # skip if there is no logfile -> to obtail the number of queue
    if (( file_count == 0 )); then
        echo "No log files found for ${i}" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
        echo "sample path: ${sample_path}" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
        continue
    fi
    failed_evt=()
    ##########################################
    ###      write summary to txt files    ###
    ##########################################
    ## indicate Sample name ##
    echo "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[ ${i} ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
    echo "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[ ${i} ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]" >> "./Run_Summary/${sample_path}/${i}/Memory_Summary_${i}.txt"
    for ((j=0; j<file_count; j++))
    do 
        ##########################################
        ###       memory summary              ###
        ##########################################
        echo ">> $j th-----------------------------------------------------------------------------------" >> "./Run_Summary/${sample_path}/${i}/Memory_Summary_${i}.txt"
        echo "                              	Usage    Requested   Allocated" >> "./Run_Summary/${sample_path}/${i}/Memory_Summary_${i}.txt"
        echo "`grep " Memory (MB)" "$log_path/log_$j.log"`" >> "./Run_Summary/${sample_path}/${i}/Memory_Summary_${i}.txt"""
        ##########################################
        ###       err, out ,log summary        ###
        ##########################################
        if ! grep -q "exit-code 0" "$log_path/log_$j.log" || grep -q "killed" "$err_path/err_$j.err" || ! grep -q "Analysis destructor completed." "$out_path/out_$j.out"; then
            echo ">> $j th-----------------------------------------------------------------------------------" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
            # exit-code 확인 : exit-code 0(정상종료)이 아니면 기록
            if ! grep -q "exit-code 0" "$log_path/log_$j.log"; then
                echo "Log: " >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
                echo "`grep "SlotName" "$log_path/log_$j.log"`" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
                echo "`grep "exit-code" "$log_path/log_$j.log"`" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
                echo "`grep " Memory (MB)" "$log_path/log_$j.log"`" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"""
            fi
            # err 파일 확인 : killed가 있으면 기록
            if grep -q "killed" "$err_path/err_$j.err"; then
                echo "Err: " >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
                echo "`grep "killed" "$err_path/err_$j.err"`" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
            fi
            # out 파일 확인 : End processing(정상종료)이 없으면 기록
            if ! grep -q "End processing" "$out_path/out_$j.out"; then
                    echo "Out: " >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
                    echo "`tail -n 7 "$out_path/out_$j.out"`" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
            fi
            failed_evt[$j]=$j 
        fi
        ################################
        ###       progress bar       ###
        ################################
        ((completed++))
        echo -ne "Progress in $i: ["
        progress=$(( (completed * 100) / total_sample )) 
        for ((k=0; k<=100; k+=2))
        do
            if [ $k -le $progress ]
            then
                echo -ne ">"
            else
                echo -ne " "
            fi
        done
        echo -ne "] $progress% \r"
    done
    echo "------------------------------------ Summary -----------------------------------------------" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
    echo "Number of failed inputlist number: ${#failed_evt[@]}" >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
    echo -n "Failed inputlist number: " >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
    for j in "${failed_evt[@]}"
    do
        echo -n "\"$j\"," >> "./Run_Summary/${sample_path}/${i}/Run_Summary_${i}.txt"
    done
done
echo ""