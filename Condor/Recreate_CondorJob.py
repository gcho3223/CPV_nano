import os
from platform import machine
import sys

samplelist = [
    #["Data_DoubleMuon_Run2016B",[]],
    #["Data_DoubleMuon_Run2016C",[]],
    #["Data_DoubleMuon_Run2016D",[]],
    #["Data_DoubleMuon_Run2016E",[]],
    #["Data_DoubleMuon_Run2016F",[]],
    #["Data_DoubleMuon_Run2016G",[]],
    #["Data_DoubleMuon_Run2016HV2",[]],
    #["Data_DoubleMuon_Run2016HV3",[]],
    #["Data_SingleMuon_Run2016B",[]],
    #["Data_SingleMuon_Run2016C",[]],
    #["Data_SingleMuon_Run2016D",[]],
    #["Data_SingleMuon_Run2016E",[]],
    #["Data_SingleMuon_Run2016F",[]],
    #["Data_SingleMuon_Run2016G",[]],
    #["Data_SingleMuon_Run2016HV2",[]],
    #["Data_SingleMuon_Run2016HV3",[]],
    #["DYJetsToLL_M_10To50",[]],
    #["DYJetsToLL_M_50",[]],
    #["ST_tW_antitop",[]],
    #["ST_tW_top",[]],
    #["TTbar_WJetToLNu",[]],
    #["TTbar_WQQ",[]],
    #["TTbar_ZQQ",[]],
    #["TTbar_ZToLLNuNu",[]],
    #["TTJets_others",[]],
    ["TTbar_Signal",["1","2"]],
    #["WJetsToLNu",[]],
    #["WW",[]],
    #["WZ",[]],
    #["ZZ",[]],
    #["TTJets_Signal_dtG_0",["4","13","16","22"]],
    #["TTJets_Signal_dtG_0p5207",["6","16","28"]],
    #["TTJets_Signal_dtG_1p0415",["9","17","27"]],
    #["TTJets_Signal_dtG_2p60364",["5","13","18","27"]],
    #["TTJets_Signal_dtG_m0p5207",["9","17"]],
    #["TTJets_Signal_dtG_m1p0415",["14"]],
    #["TTJets_Signal_dtG_m2p60364",["3","8","12","19","20","25"]]
]

runPeriod = "UL2016PreVFP"
Channels = "MuMu"
NumEvt = -1
version = "SelLep"

if Channels == "MuMu":
    confch = "ULSummer20/" + runPeriod + "/dimuon.config"
elif Channels == "ElEl":
    confch = "ULSummer20/" + runPeriod + "/dielec.config"
elif Channels == "MuEl":
    confch = "ULSummer20/" + runPeriod + "/muelec.config"

for sample in range(len(samplelist)):
    ### set up jobversion ###
    if "dtG" in samplelist[sample][0]:
        jobversion = "Job_Version/%s/CPV_Sample" %(version)
    else:
        jobversion = "Job_Version/%s/Dataset" %(version)
    ### set up directory ###
    sample_dir = f"{jobversion}/{runPeriod}/{Channels}/{samplelist[sample][0]}"
    ### delete log, err, out files for the failed jobs in log_condor directory ###
    for f in os.listdir(f"{sample_dir}/log_condor/err"):
        os.remove(os.path.join(f"{sample_dir}/log_condor/err", f))
    for f in os.listdir(f"{sample_dir}/log_condor/log"):
        os.remove(os.path.join(f"{sample_dir}/log_condor/log", f))
    for f in os.listdir(f"{sample_dir}/log_condor/out"):
        os.remove(os.path.join(f"{sample_dir}/log_condor/out", f))
    ### created executive condor job: condor_run.sh ###
    with open(f"{sample_dir}/recondor_run.sh", "w") as file:
        file.write("#!/bin/bash \n\n")
        file.write("export SCRAM_ARCH=slc7_amd64_gcc12\n")
        file.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n\n")
        file.write("export LD_PRELOAD=""/usr/lib64/libpdcap.so"" \n")
        file.write("export LD_LIBRARY_PATH=""$LD_LIBRARY_PATH:/usr/lib64/dcap"" \n\n")
        file.write("fileListNum=$((${1}+1)) \n")
        file.write("cd /u/user/gcho/TopPhysics/CPV/NanoAOD/CMSSW_13_3_0/src/SSBNanoAODANCode/ \n")
        file.write("cmsenv \n\n")
        file.write("mkdir -p ./output/%s/\n" %(sample_dir))
        if "dtG" in samplelist[sample][0]:
            file.write("./ssb_analysis CEDM_Sample/%s/%s_${fileListNum}.list /%s/%s_${fileListNum}.root 0 %s \n" % (samplelist[sample][0],samplelist[sample][0],sample_dir,samplelist[sample][0],samplelist[sample][0]))
        else:
            file.write("./ssb_analysis %s/%s_${fileListNum}.list %s/%s_${fileListNum}.root %s None %s %s \n" % (runPeriod,samplelist[sample][0],sample_dir,samplelist[sample][0],confch,runPeriod,NumEvt))
    ### created submit condor job: condor_sub.sub ###
    with open(f"{sample_dir}/recondor_sub.sub", "w") as file:
        file.write("universe = vanilla \n")
        file.write("executable = condor_run.sh \n\n")
        ## set up memory request: data: 6144, MC: 3072, CPV: 2048
        if "Data" in samplelist[sample][0]:
            file.write("request_memory = 6144 MB \n\n") ## data samples
        elif "dtG" in samplelist[sample][0]:
            file.write("request_memory = 2048 MB \n\n") ## CPV samples
        else:
            file.write("request_memory = 3072 MB \n\n") ## MC samples
        file.write("should_transfer_files = YES \n")
        file.write("when_to_transfer_output = ON_EXIT \n\n")
        file.write('Requirements = (machine =!= "cluster291.knu.ac.kr") && (TARGET.Arch == "X86_64") && (TARGET.OpSys == "LINUX") && (TARGET.HasFileTransfer)\n')
        file.write("JobBatchName = %s \n\n" %(samplelist[sample][0]) )
        for job_number in samplelist[sample][1]:
            file.write(f"arguments = {job_number} \n")
            file.write("output = /u/user/gcho/TopPhysics/CPV/NanoAOD/CMSSW_13_3_0/src/SSBNanoAODANCode/Condor/%s/log_condor/out/out_$(Process).out \n" %(sample_dir) )
            file.write("error  = /u/user/gcho/TopPhysics/CPV/NanoAOD/CMSSW_13_3_0/src/SSBNanoAODANCode/Condor/%s/log_condor/err/err_$(Process).err \n" %(sample_dir) )
            file.write("log    = /u/user/gcho/TopPhysics/CPV/NanoAOD/CMSSW_13_3_0/src/SSBNanoAODANCode/Condor/%s/log_condor/log/log_$(Process).log \n\n" %(sample_dir) )
            file.write("queue 1\n\n")
        #file.write("queue %s" %(samplelist[sample][1]) )
    print("%s condor submit file created!!! " %(samplelist[sample][0]))
