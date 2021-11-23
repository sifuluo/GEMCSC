#include <iostream>
#include <vector>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TH2.h"

#include "TreeReader.cc"
#include "TreeDataFormat.cc"
#include "Tools.cc"

using namespace std;

void SimplePlot() {
  TChain *tree = new TChain("NtupleMaker/eventTree");
  tree->Add("Privatea_Run4_0.root");

  TreeReader *tr = new TreeReader(tree);
  TFile *out = new TFile("splots.root","RECREATE");
  out->cd();

  vector<TString> LDisks{"D1","D2"};
  vector<TString> LRings{"R1","R2"};

  Long64_t nentries = tree->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    PrintProgress(jentry, nentries,1);
    Long64_t ientry = tree->LoadTree(jentry);
    tr->ReadTree(jentry);
    vector<vector<StationData> > & Stations = tr->Evt.Stations;
    for (unsigned disk = 0; disk < Stations.size(); ++disk) for (unsigned ring = 0; ring > Stations[disk].size(); ++ring) {
      if (disk > 0 || ring >= LRings.size()) continue;
      StationData& ThisStation = Stations[disk][ring];
      for (unsigned ittp = 0; ittp < ThisStation.TPInfos.size(); ++ittp) {
        TPContent& ThisTP = ThisStation.TPInfos[ittp];
        const vector<float> CSCSimHitAve_V{ThisTP.CSCSimHitAve.eta, ThisTP.CSCSimHitAve.phi, ThisTP.CSCSimHitAve.r, ThisTP.CSCSimHitAve.z};
        const vector<float> GEMSimHitAve_V{ThisTP.GEMSimHitAve.eta, ThisTP.GEMSimHitAve.phi, ThisTP.GEMSimHitAve.r, ThisTP.GEMSimHitAve.z};
        bool CanRecoCSC = (ThisTP.NSimHitsCSC > 3);
        bool CanRecoGEM = (ThisTP.NSimHitsGEM > 0);

        //Now begins the analyze for matched information

        
      } // End of TP loop


    } // End of detector loop

  }  // End of event loop
}
