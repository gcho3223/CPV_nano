#include <fstream>
#include "TSystem.h"

void inputlistMaker()
{
  // setup dir path
  TString baseDir = "/pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MC/TT_TuneCUETP8M1_13TeV_powheg_80X_V2_forMoriond_backup/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_20200626_145255/200626_060138";
  
  // sub dir
  TString subDirs[4] = {"0000", "0001", "0002", "0003"};
  
  TString outputDir = "TTJets_Signal_backup";
  gSystem->Exec(Form("mkdir -p %s", outputDir.Data()));
  
  int fileCount = 0;
  int listCount = 1;
  ofstream outFile;
  
  for (int i = 0; i < 4; i++)
  {
    TString currentDir = baseDir + "/" + subDirs[i];
    void *dirp = gSystem->OpenDirectory(currentDir.Data());
    if (!dirp) continue;
    
    const char *entry;
    while ((entry = gSystem->GetDirEntry(dirp)))
    {
      TString fileName = entry;
      if (!fileName.EndsWith(".root")) continue;
      
      if (fileCount % 10 == 0) {
        if (outFile.is_open()) outFile.close();
        TString listName = Form("%s/TTJets_Signal_%d.list", outputDir.Data(), listCount++);
        outFile.open(listName.Data());
      }
      
      TString fullPath = "dcap://cluster142.knu.ac.kr/" + currentDir + "/" + fileName;
      outFile << fullPath << endl;
      
      fileCount++;
    }
    
    gSystem->FreeDirectory(dirp);
  }
  
  if (outFile.is_open()) outFile.close();
  
  cout << "Total " << fileCount << " files stored in " << outputDir << "as " << (listCount-1) << " inputfiles" << endl;
}