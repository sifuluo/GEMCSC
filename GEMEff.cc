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

void GEMEff() {
  TChain *tree = new TChain("NtupleMaker/eventTree");
  bool BatchRun = true;
  TString outputname;

  if (BatchRun) {
    TString EOSPath = "/eos/user/s/siluo/";
    TString NtuplePath = "MuGun/Ntuple/";
    vector<TString> DataSets{"SMNoPU0"};
    // vector<TString> DataSets{"DMNoPU0"};
    if (DataSets.size() == 1) outputname = DataSets[0];
    else outputname = DataSets[0] + "With" + to_string(DataSets.size());
    for (auto DataSet : DataSets) {
      for (unsigned i = 0; i < 100; ++i) {
        TString DataFile = Form(EOSPath+NtuplePath+DataSet+"/out_%i.root",int(i));
        tree->Add(DataFile);
        cout << "Reading Data File " << DataFile<< endl;
      }
    }
  }
  else {
    TString LocalFile = "../Nov30/out/out.root";
    tree->Add(LocalFile);
    outputname = "plots";
    cout << "Reading Data File " << LocalFile << endl;
  }
  TreeReader *tr = new TreeReader(tree);
  TFile *out = new TFile(outputname + ".root","RECREATE");
  cout << "Saving into " << outputname + ".root"<<endl;
  out->cd();

  // Lists of Disks and Rings
  vector<TString> LDisks{"D0","D1","D2"};
  vector<TString> LRings{"R1"};
  // vector<TString> LRings{"R1","R2","R3","R0"};

  // Branch entry Multiplicity
  // TH2F* Multi = new TH2F("Multi","Branch Multiplicity; Branches; Multiplicity",20,0,20,10,0,10);
  // vector<TString> LBranches{"tp","cscSimHit","gemSimHit"};

  // Efficiency plots
  vector<TString> Eff_T{"CSCReco","CSCMatch","GEMReco","GEMMatch","ClusterReco","ClusterMatch","PadReco", "PadMatch"};
  vector<TString> Eff_V{"Eta","Phi","R","z"};
  vector<TString> Eff_V2{"#eta","#phi","R","z"};
  vector<int> Eff_V_div{60,70,80,220};
  vector<float> Eff_V_low{-3.0,-3.5,0,-1100};
  vector<float> Eff_V_up{3.0,3.5,800,1100};
  vector<vector<vector <TEfficiency*> > > DetEff;
  DetEff.resize(Eff_T.size());
  for (unsigned iET = 0; iET < Eff_T.size(); ++iET) {
    DetEff[iET].resize(Eff_V.size());
    for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
      DetEff[iET][iEV].resize(LDisks.size());
      for (unsigned iES = 0; iES < LDisks.size(); ++iES) {
        TString PlotName = Eff_T[iET] + "EffVs" + Eff_V[iEV] + "_" + LDisks[iES];
        TString PlotTitle = Eff_T[iET] + " Efficiency Vs TP SimHit " + Eff_V2[iEV] + " at " + LDisks[iES] + "; TP SimHit " + Eff_V2[iEV] + "; Efficiency";
        DetEff[iET][iEV][iES] = new TEfficiency(PlotName, PlotTitle, Eff_V_div[iEV], Eff_V_low[iEV], Eff_V_up[iEV]);
        DetEff[iET][iEV][iES]->SetDirectory(out);
      }
    }
  }

  // Closest Digi plots
  bool LookForDigiInMatched = false;
  vector<TString> Close_T{"CSC","GEM","Cluster","Pad"};
  vector<TString> Close_N{"ClosestDigi"};
  vector<int> Close_V_div{   200 , 200 , 60, 80, 70 , 220};//dEta,dPhi,Eta,Phi,R,Z
  vector<float> Close_V_low{0  , 0  , -3, -4, 0  , -1100};
  vector<float> Close_V_up{ 0.05, 0.002, 3 , 4 , 700, 1100};
  vector<vector<vector<vector<TH2F*> > > > TH2Plots;
  TH2Plots.resize(Close_N.size());
  for (auto& vec : TH2Plots) vec.resize(Close_T.size());
  for (unsigned iCT = 0; iCT < Close_T.size(); ++iCT) {
    for (auto& vec : TH2Plots) vec[iCT].resize(LDisks.size());
    for (unsigned iCS = 0; iCS < LDisks.size(); ++iCS) {
      for (auto& vec : TH2Plots) vec[iCT][iCS].resize(LRings.size());
      for (unsigned iCR = 0; iCR < LRings.size(); ++iCR) {
        TH2Plots[0][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "ClosestDigi_" + LDisks[iCS] + LRings[iCR], "Closest Digi in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; #Delta #eta; #Delta #phi", Close_V_div[0], Close_V_low[0], Close_V_up[0], Close_V_div[1], Close_V_low[1], Close_V_up[1] );
        // TH2Plots[1][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "SimHitRVZ_" + LDisks[iCS] + LRings[iCR], "SimHitRVZ in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; z; r", Close_V_div[5], Close_V_low[5], Close_V_up[5], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
        // TH2Plots[2][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "SimHitRVPhi_" + LDisks[iCS] + LRings[iCR], "SimHitRVPhi in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; phi; r", Close_V_div[3], Close_V_low[3], Close_V_up[3], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
        // TH2Plots[3][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "DigiRVZ_" + LDisks[iCS] + LRings[iCR], "DigiRVZ in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; z; r", Close_V_div[5], Close_V_low[5], Close_V_up[5], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
        // TH2Plots[4][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "DigiRVPhi_" + LDisks[iCS] + LRings[iCR], "DigiRVPhi in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; phi; r", Close_V_div[3], Close_V_low[3], Close_V_up[3], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
      }
    }
  }// This initialization is completely rubbish and need to be revised. It will try to resize and repeatedly initialize for nothing if more than 1 Close_N

  // GEM-CSC Matching plots
  // Make the previous plot collections at plots[disk][ring][variable][dettype][name]
  vector<TString> CSCGEMEff_N{"Matched","All"};
  vector<TString> CSCGEMEff_V{"Eta","Phi","R","z","Slope","Chamber"};
  vector<TString> CSCGEMEff_V2{"#eta","#phi","R","z","Slope","Chamber"};
  vector<int> CSCGEMEff_V_div{60,70,400,220,20,2};
  vector<float> CSCGEMEff_V_low{-3.0,-3.5,0,-1100,-10,0};
  vector<float> CSCGEMEff_V_up{3.0,3.5,400,1100,10,2};
  vector<vector<vector<vector<vector<TEfficiency*> > > > >CSCGEMEff;
  vector<vector<vector<vector<vector<TH1F*> > > > >CSCGEMEffH;

  CSCGEMEff.resize(LDisks.size());
  CSCGEMEffH.resize(LDisks.size());
  for (unsigned iED = 0; iED < LDisks.size(); ++iED) {
    CSCGEMEff[iED].resize(LRings.size());
    CSCGEMEffH[iED].resize(LRings.size());
    for (unsigned iER = 0; iER < LRings.size(); ++iER) {
      CSCGEMEff[iED][iER].resize(CSCGEMEff_V.size());
      CSCGEMEffH[iED][iER].resize(CSCGEMEff_V.size());
      for (unsigned iEV = 0; iEV < CSCGEMEff_V.size(); ++iEV) {
        CSCGEMEff[iED][iER][iEV].resize(CSCGEMEff_N.size());
        CSCGEMEffH[iED][iER][iEV].resize(CSCGEMEff_N.size());
        for (unsigned iEN = 0; iEN < CSCGEMEff_N.size(); ++iEN) {
          CSCGEMEff[iED][iER][iEV][iEN].resize(4);
          CSCGEMEff[iED][iER][iEV][iEN][0] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_Either_" + LDisks[iED] + LRings[iER], "CSCGEM(Either) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEff[iED][iER][iEV][iEN][1] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_1_" + LDisks[iED] + LRings[iER], "CSCGEM(Layer1) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEff[iED][iER][iEV][iEN][2] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_2_" + LDisks[iED] + LRings[iER], "CSCGEM(Layer2) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEff[iED][iER][iEV][iEN][3] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_Both_" + LDisks[iED] + LRings[iER], "CSCGEM(Both) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          for (auto eff_ : CSCGEMEff[iED][iER][iEV][iEN]) eff_->SetDirectory(out);
          CSCGEMEffH[iED][iER][iEV][iEN].resize(4);
          CSCGEMEffH[iED][iER][iEV][iEN][0] = new TH1F("CSCGEMVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_Either_" + LDisks[iED] + LRings[iER], "CSCGEM(Either) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEffH[iED][iER][iEV][iEN][1] = new TH1F("CSCGEMVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_1_" + LDisks[iED] + LRings[iER], "CSCGEM(Layer1) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEffH[iED][iER][iEV][iEN][2] = new TH1F("CSCGEMVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_2_" + LDisks[iED] + LRings[iER], "CSCGEM(Layer2) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEffH[iED][iER][iEV][iEN][3] = new TH1F("CSCGEMVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_Both_" + LDisks[iED] + LRings[iER], "CSCGEM(Both) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          for (auto eff_ : CSCGEMEffH[iED][iER][iEV][iEN]) eff_->SetDirectory(out);
        }
      }
    }
  }

  vector<vector<TH2F*> >dPhiVsSlope;
  dPhiVsSlope.resize(LDisks.size());
  for (unsigned iED = 0; iED < LDisks.size(); ++iED) {
    dPhiVsSlope[iED].resize(LRings.size());
    for (unsigned iER = 0; iER < LRings.size(); ++iER) {
      dPhiVsSlope[iED][iER] = new TH2F("dPhiVsSlope_" + LDisks[iED] + LRings[iER],"#Delta #phi Vs Slope at " + LDisks[iED] + " " + LRings[iER] + "; Slope; #Delta #phi", 30,-15,15,200,-1,1);
    }
  }

  int cancscpos(0),cancscneg(0),cangempos(0),cangemneg(0),canbothpos(0),canbothneg(0);
  int tppos(0),tpneg(0),evttppos(0),evttpneg(0);
  int GhostMatchGEM1(0), GhostMatchGEM2;
  vector<int> tpexp, tpreco;
  tpexp.resize(5,0);
  tpreco.resize(5,0);
  // TP reconstruction plots
  //how often is a sim muon matched to ALCT, CLCT, 1GEM cluster, 2GEM clusters? How often do you get any combination? How often do you get an LCT?
  //What is the type of LCT matched to the muon? etc. questions like this need to be answered in order to understand any inefficiencies
  //because there are 5 types of LCTs we can make in ME1/1: ALCT-CLCT, ALCT-CLCT-1GEM, ALCT-CLCT-2GEM, ALCT-2GEM and CLCT-2GEM, we need to know exactly which type was expected to show up in a chamber based on the presence of ALCT/CLCT/GEM, but somehow did not.

  Long64_t nentries = tree->GetEntries();
  int PrintProgressInterval = 1;
  PrintProgress(0, nentries + 1, PrintProgressInterval);
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    // cout << "Finished event init " << jentry <<endl;
    Long64_t ientry = tree->LoadTree(jentry);
    // tree->GetEntry(jentry);
    // cout << "Finished Loading event" <<endl;
    tr->ReadTree(jentry);
    vector<tp*>& tps = tr->Evt.MuonTPs;
    for (tp* thistp : tps) {
      if (thistp->eta > 0) evttppos++;
      else evttpneg++;
    }
    // cout << "READ" <<endl;
    vector<vector<StationData> >& Stations = tr->Evt.Stations;

    for (unsigned disk = 0; disk < Stations.size(); ++disk) for (unsigned ring = 0; ring < Stations[disk].size(); ++ring) {
      if (disk > 2) break; // Considering only first 2 disks;
      if (ring >= LRings.size()) continue;
      StationData& ThisStation = Stations[disk][ring];
      for (unsigned ittp = 0; ittp < ThisStation.TPInfos.size(); ++ittp) {
        TPContent& ThisTP = ThisStation.TPInfos[ittp];
        const vector<float> CSCSimHitAve_V{ThisTP.CSCSimHitAve.eta, ThisTP.CSCSimHitAve.phi, ThisTP.CSCSimHitAve.r, ThisTP.CSCSimHitAve.z};
        const vector<float> GEMSimHitAve_V{ThisTP.GEMSimHitAve.eta, ThisTP.GEMSimHitAve.phi, ThisTP.GEMSimHitAve.r, ThisTP.GEMSimHitAve.z};
        bool CanRecoCSC = (ThisTP.NSimHitsCSC > 3);
        bool CanRecoGEM = (ThisTP.NSimHitsGEM > 0);
        // if (!CanRecoCSC) {
        //   cout << "This TP cannot reconstruct CSCStubs, NSimHitsCSC = " << ThisTP.NSimHitsCSC << ", ThisTP eta = " << ThisTP.TP->eta << ", SimHitAveEta = " << ThisTP.CSCSimHitAve.eta << endl;
        // }
        // if (!CanRecoGEM) {
        //   cout << "This TP cannot reconstruct GEMDigis, NSimHitsGEM = " << ThisTP.NSimHitsGEM << ", ThisTP eta = " << ThisTP.TP->eta << ", SimHitAveEta = " << ThisTP.GEMSimHitAve.eta << endl;
        // }
        if (ThisTP.TP->eta > 0) tppos++;
        else if (ThisTP.TP->eta < 0) tpneg++;
        if (CanRecoCSC) {
          if (ThisTP.CSCSimHitAve.eta > 0) cancscpos++;
          else if (ThisTP.CSCSimHitAve.eta < 0) cancscneg++;
          bool DigiInMatch = ThisTP.MatchCSCStubs.size();
          // TH2Plots[1][0][disk][ring]->Fill(CSCSimHitAve_V[3],CSCSimHitAve_V[2]);
          // TH2Plots[2][0][disk][ring]->Fill(CSCSimHitAve_V[1],CSCSimHitAve_V[2]);
          if (DigiInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[0][iEV][disk]->Fill(1,CSCSimHitAve_V[iEV]); // CSCReco
              DetEff[1][iEV][disk]->Fill(1,CSCSimHitAve_V[iEV]); // CSCMatch
            }
            if (disk == 1 && ThisTP.TP->cross[0]) ThisTP.TP->reco[0] = true;
            if (disk == 2 && ThisTP.TP->cross[1]) ThisTP.TP->reco[1] = true;
          }
          else {
            bool DigiInGeneric = false;
            for (unsigned idigi = 0; idigi < ThisStation.AllCSCStubs.size(); ++idigi) {
              CSCStub* CSCDigi = ThisStation.AllCSCStubs[idigi];
              DigiInGeneric = IsCloseCSC(CSCDigi->eta, CSCSimHitAve_V[0], CSCDigi->phi, CSCSimHitAve_V[1]);
              if (DigiInGeneric) break;
            }

            if (DigiInGeneric) {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[0][iEV][disk]->Fill(1,CSCSimHitAve_V[iEV]); // CSCReco
                DetEff[1][iEV][disk]->Fill(0,CSCSimHitAve_V[iEV]); // CSCMatch
              }
              if (disk == 1 && ThisTP.TP->cross[0]) ThisTP.TP->reco[0] = true;
              if (disk == 2 && ThisTP.TP->cross[1]) ThisTP.TP->reco[1] = true;
            }
            else {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[0][iEV][disk]->Fill(0,CSCSimHitAve_V[iEV]); // CSCReco
              }
            }
          }
          float ClosedR = 999;
          float ClosedEta(999), ClosedPhi(999);
          vector<CSCStub*>& DigisToLook = (LookForDigiInMatched ? ThisTP.MatchCSCStubs : ThisStation.AllCSCStubs);
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            CSCStub* CSCDigi = DigisToLook[idigi];
            if (!(CSCDigi->valid)) continue;
            // TH2Plots[3][0][disk][ring]->Fill(CSCDigi.z,CSCDigi.r);
            // TH2Plots[4][0][disk][ring]->Fill(CSCDigi.phi,CSCDigi.r);
            vector<float> delta = CalcdR(CSCDigi->eta, CSCSimHitAve_V[0], CSCDigi->phi, CSCSimHitAve_V[1]);
            if (delta[0] < ClosedR) {
              ClosedR = delta[0];
              ClosedEta = delta[1];
              ClosedPhi = delta[2];
            }
          }
          if (ClosedR < 999) TH2Plots[0][0][disk][ring]->Fill(ClosedEta,ClosedPhi);
        }

        if (CanRecoGEM) { //GEMDigi
          if (ThisTP.CSCSimHitAve.eta > 0) cangempos++;
          else if (ThisTP.CSCSimHitAve.eta < 0) cangemneg++;
          bool DigiInMatch = ThisTP.MatchGEMDigis.size();
          if (DigiInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[2][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMReco
              DetEff[3][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMMatch
            }
            if (disk == 0 && ThisTP.TP->cross[2]) ThisTP.TP->reco[2]= true;
            if (disk == 1 && ThisTP.TP->cross[3]) ThisTP.TP->reco[3] = true;
            if (disk == 2 && ThisTP.TP->cross[4]) ThisTP.TP->reco[4] = true;
          }
          else {
            bool DigiInGeneric = false;
            for (unsigned idigi = 0; idigi < ThisStation.AllGEMDigis.size(); ++idigi) {
              GEMDigi* GEMDigi = ThisStation.AllGEMDigis[idigi];
              DigiInGeneric = IsCloseGEM(GEMDigi->eta, GEMSimHitAve_V[0], GEMDigi->phi, GEMSimHitAve_V[1]);
              if (DigiInGeneric) break;
            }

            if (DigiInGeneric) {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[2][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMReco
                DetEff[3][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMMatch
              }
              if (disk == 0 && ThisTP.TP->cross[2]) ThisTP.TP->reco[2]= true;
              if (disk == 1 && ThisTP.TP->cross[3]) ThisTP.TP->reco[3] = true;
              if (disk == 2 && ThisTP.TP->cross[4]) ThisTP.TP->reco[4] = true;
            }
            else {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[2][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMReco
              }
            }
          }
          float ClosedR = 999;
          float ClosedEta(999), ClosedPhi(999);
          vector<GEMDigi*>& DigisToLook = (LookForDigiInMatched ? ThisTP.MatchGEMDigis : ThisStation.AllGEMDigis);
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            GEMDigi* GEMDigi = DigisToLook[idigi];
            vector<float> delta = CalcdR(GEMDigi->eta, GEMSimHitAve_V[0], GEMDigi->phi, GEMSimHitAve_V[1]);
            if (delta[0] < ClosedR) {
              ClosedR = delta[0];
              ClosedEta = delta[1];
              ClosedPhi = delta[2];
            }
          }
          if (ClosedR < 999) TH2Plots[0][1][disk][ring]->Fill(ClosedEta,ClosedPhi);
        }

        if (CanRecoGEM) { //GEMPadClusters
          bool ClusterInMatch = ThisTP.MatchGEMPadDigiClusters.size();
          if (ClusterInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[4][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]);
              DetEff[5][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]);
            }
          }
          else {
            bool ClusterInGeneric = false;
            for (unsigned ic = 0; ic < ThisStation.AllGEMPadDigiClusters.size(); ++ic) {
              GEMPadDigiCluster* cl = ThisStation.AllGEMPadDigiClusters[ic];
              // ClusterInGeneric = IsCloseCluster(cl->eta, GEMSimHitAve_V[0], cl->phi, GEMSimHitAve_V[1]);
              ClusterInGeneric = cl->det.SameXY(ThisTP.GEMDetId);
              if (ClusterInGeneric) break;
            }

            if (ClusterInGeneric) {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[4][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMReco
                DetEff[5][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMMatch
              }
            }
            else {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[4][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMReco
              }
            }
          }
          float ClosedR = 999;
          float ClosedEta(999), ClosedPhi(999);
          vector<GEMPadDigiCluster*>& DigisToLook = (LookForDigiInMatched ? ThisTP.MatchGEMPadDigiClusters : ThisStation.AllGEMPadDigiClusters);
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            GEMPadDigiCluster* cl = DigisToLook[idigi];
            vector<float> delta = CalcdR(cl->eta, GEMSimHitAve_V[0], cl->phi, GEMSimHitAve_V[1]);
            if (delta[0] < ClosedR) {
              ClosedR = delta[0];
              ClosedEta = delta[1];
              ClosedPhi = delta[2];
            }
          }
          if (ClosedR < 999) TH2Plots[0][2][disk][ring]->Fill(ClosedEta,ClosedPhi);
          // if (ClosedR < 999) cout << Form("dR = %f, dEta = %f, dPhi = %f", ClosedR, ClosedEta, ClosedPhi);
        }

        if (CanRecoGEM) { //GEMPad
          bool PadInMatch = ThisTP.MatchGEMPadDigis.size();
          if (PadInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[6][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]);
              DetEff[7][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]);
            }
          }
          else {
            bool PadInGeneric = false;
            for (unsigned ip = 0; ip < ThisStation.AllGEMPadDigis.size(); ++ip) {
              GEMPadDigi* pad = ThisStation.AllGEMPadDigis[ip];
              // PadInGeneric = IsClosePad(pad->eta, GEMSimHitAve_V[0], pad->phi, GEMSimHitAve_V[1]);
              PadInGeneric = pad->det.SameXY(ThisTP.GEMDetId);
              if (PadInGeneric) break;
            }

            if (PadInGeneric) {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[6][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMReco
                DetEff[7][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMMatch
              }
            }
            else {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[6][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMReco
              }
            }
          }
          float ClosedR = 999;
          float ClosedEta(999), ClosedPhi(999);
          vector<GEMPadDigi*>& DigisToLook = (LookForDigiInMatched ? ThisTP.MatchGEMPadDigis : ThisStation.AllGEMPadDigis);
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            GEMPadDigi* pad = DigisToLook[idigi];
            vector<float> delta = CalcdR(pad->eta, GEMSimHitAve_V[0], pad->phi, GEMSimHitAve_V[1]);
            if (delta[0] < ClosedR) {
              ClosedR = delta[0];
              ClosedEta = delta[1];
              ClosedPhi = delta[2];
            }
          }
          if (ClosedR < 999) TH2Plots[0][3][disk][ring]->Fill(ClosedEta,ClosedPhi);
          // if (ClosedR < 999) cout << Form("dR = %f, dEta = %f, dPhi = %f", ClosedR, ClosedEta, ClosedPhi);
        }

        if (CanRecoCSC && CanRecoGEM) {
          if (ThisTP.CSCSimHitAve.eta > 0) canbothpos++;
          else if (ThisTP.CSCSimHitAve.eta < 0) canbothneg++;
          vector<float> CSCSimHitAve_VExt = CSCSimHitAve_V;
          CSCSimHitAve_VExt.push_back(0);
          CSCSimHitAve_VExt.push_back(0);
          vector<vector<float> > GEMSimHitV{{},{}};
          bool hasGEMSH1(false), hasGEMSH2(false), overwrittensh(false);
          if (disk > 0) {
            for (auto sh : ThisTP.GEMSimHits) {
              vector<float> GEMSimHitV_{sh->eta, sh->phi, sh->r, sh->z, 0, (float)sh->det.chamber};
              if (sh->det.layer == 1) {
                // if (hasGEMSH1) {
                //   cout << "Overwriting layer1?" <<endl;
                //   overwrittensh = true;
                // }
                GEMSimHitV[0] = GEMSimHitV_;
                hasGEMSH1 = true;
              }
              else if (sh->det.layer == 2) {
                // if (hasGEMSH2) {
                //   cout << "Overwriting layer2?" <<endl;
                //   overwrittensh = true;
                // }
                GEMSimHitV[1] = GEMSimHitV_;
                hasGEMSH2 = true;
              }
            }
          }
          // if (overwrittensh) {
          //   cout << "GEMSimHit Overwritten, GEMSimHits are:" <<endl;
          //   for (auto sh : ThisTP.GEMSimHits) {
          //     cout << *sh <<endl;
          //   }
          // }

          for (unsigned icsc = 0; icsc < ThisTP.MatchCSCStubs.size(); ++icsc) {
            CSCStub* lct = ThisTP.MatchCSCStubs[icsc];
            bool hasGEM1 = (lct->GEM1pad != 255);
            bool hasGEM2 = (lct->GEM2pad != 255);
            // Are the GEM1/GEM2 the same as the GEMPadCluster collection?
            CSCSimHitAve_VExt[4] = lct->slope;
            CSCSimHitAve_VExt[5] = lct->det.chamber;
            for (unsigned iEV = 0; iEV < CSCGEMEff_V.size(); ++iEV) {
              CSCGEMEff[disk][ring][iEV][0][0]->Fill((hasGEM1 || hasGEM2), CSCSimHitAve_VExt[iEV]);
              CSCGEMEff[disk][ring][iEV][0][1]->Fill((hasGEM1), CSCSimHitAve_VExt[iEV]);
              CSCGEMEff[disk][ring][iEV][0][2]->Fill((hasGEM2), CSCSimHitAve_VExt[iEV]);
              CSCGEMEff[disk][ring][iEV][0][3]->Fill((hasGEM1 && hasGEM2), CSCSimHitAve_VExt[iEV]);
              if (disk == 0) continue;
              bool True1 = hasGEMSH1 && hasGEM1;
              bool True2 = hasGEMSH2 && hasGEM2;
              if (hasGEM1 && !hasGEMSH1) {
                cout << "Ghost GEM1 in MatchCSC" << endl;
                GhostMatchGEM1++;
              }
              if (hasGEM2 && !hasGEMSH2) {
                cout << "Ghost GEM2 in MatchCSC" << endl;
                GhostMatchGEM2++;
              }
              if (True1 || True2) CSCGEMEffH[disk][ring][iEV][0][0]->Fill(GEMSimHitV[(True1 ? 0 : 1)][iEV]);
              if (True1) CSCGEMEffH[disk][ring][iEV][0][1]->Fill(GEMSimHitV[0][iEV]);
              if (True2) CSCGEMEffH[disk][ring][iEV][0][2]->Fill(GEMSimHitV[1][iEV]);
              if (True1 && True2) CSCGEMEffH[disk][ring][iEV][0][3]->Fill(GEMSimHitV[0][iEV]);
            }
          }
        }

        if (CanRecoCSC && CanRecoGEM && ThisTP.MatchCSCStubs.size() && ThisTP.MatchGEMDigis.size()) {
          for (unsigned icsc = 0; icsc < ThisTP.MatchCSCStubs.size(); ++icsc) {
            float slope = ThisTP.MatchCSCStubs[icsc]->slope;
            float zsign = (ThisTP.MatchCSCStubs[icsc]->z > 0) - (ThisTP.MatchCSCStubs[icsc]->z < 0);
            for (unsigned igem = 0; igem < ThisTP.MatchGEMDigis.size(); ++igem) {
              float dphi = fabs(ThisTP.MatchGEMDigis[igem]->phi - ThisTP.MatchCSCStubs[icsc]->phi);
              while (dphi > 2 * M_PI) dphi -= 2 * M_PI;
              dPhiVsSlope[disk][ring]->Fill(slope,dphi);
            }
          }
        }
      } // End of tp loop

      if (true) { // Enclose allCSC scope
        for (unsigned icsc = 0; icsc < ThisStation.AllCSCStubs.size(); ++icsc) {
          CSCStub* lct = ThisStation.AllCSCStubs[icsc];
          // if (lct->det.chamber%2 != 0) continue;
          bool hasGEM1 = (lct->GEM1pad != 255);
          bool hasGEM2 = (lct->GEM2pad != 255);
          vector<float> CSCStub_V{lct->eta, lct->phi, lct->r, lct->z, (float)lct->slope, (float)(lct->det.chamber % 2)};
          for (unsigned iEV = 0; iEV < CSCGEMEff_V.size(); ++iEV) {
            CSCGEMEff[disk][ring][iEV][1][0]->Fill((hasGEM1 || hasGEM2), CSCStub_V[iEV]);
            CSCGEMEff[disk][ring][iEV][1][1]->Fill((hasGEM1), CSCStub_V[iEV]);
            CSCGEMEff[disk][ring][iEV][1][2]->Fill((hasGEM2), CSCStub_V[iEV]);
            CSCGEMEff[disk][ring][iEV][1][3]->Fill((hasGEM1 && hasGEM2), CSCStub_V[iEV]);
          }
        }
      }


    } // End of detector loop
    for (tp* thistp : tps) {
      for (unsigned itdet = 0; itdet < 5; ++itdet) {
        if (thistp->cross[itdet]) tpexp[itdet]++;
        if (thistp->reco[itdet]) tpreco[itdet]++;
      }
    }
    PrintProgress(jentry+1, nentries+1,PrintProgressInterval);
  } // End of event loop

  cout << "Can CSC Pos = " << cancscpos << ", Neg = " << cancscneg<< endl;
  cout << "Can GEM Pos = " << cangempos << ", Neg = " << cangemneg<< endl;
  cout << "Can Both Pos = " << canbothpos << ", Neg = " << canbothneg<< endl;
  cout << "TP in +eta : " << tppos << ", TP in -eta : " << tpneg << endl;
  cout << "Event TP in +eta : " << evttppos << ", TP in -eta : " << evttpneg << endl;
  cout << "Ghost GEM1 in MatchCSC: " << GhostMatchGEM1 <<endl;
  cout << "Ghost GEM2 in MatchCSC: " << GhostMatchGEM2 <<endl;
  vector<string> dettags{"ME11","ME21","GE0","GE11","GE21"};
  for (unsigned i = 0; i < 5; ++i) {
    cout << "For " << dettags[i] << " : Expected " << tpexp[i] << " , Reconstructed " << tpreco[i] << " , Efficiency: " << double(tpreco[i])/double(tpexp[i]) << endl;
  }
  out->Write();
  out->Save();
  out->Close();
}
