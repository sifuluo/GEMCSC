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
  // tree->Add("out_Run4.root");
  tree->Add("ntuple/Privatea_Run4_0.root");
  // tree->Add("ntuple/Privatec_Run4.root");
  // tree->Add("pu.root");

  TreeReader *tr = new TreeReader(tree);
  TFile *out = new TFile("plots.root","RECREATE");
  out->cd();

  // Lists of Disks and Rings
  vector<TString> LDisks{"D1","D2"};
  vector<TString> LRings{"R1","R2"};
  // vector<TString> LRings{"R1","R2","R3","R0"};

  // Branch entry Multiplicity
  // TH2F* Multi = new TH2F("Multi","Branch Multiplicity; Branches; Multiplicity",20,0,20,10,0,10);
  // vector<TString> LBranches{"tp","cscSimHit","gemSimHit"};

  // Efficiency plots
  vector<TString> Eff_T{"CSCReco","CSCMatch","GEMReco","GEMMatch","ClusterReco","ClusterMatch"};
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
  bool LookForDigiInMatched = true;
  vector<TString> Close_T{"CSC","GEM","Cluster"};
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
  vector<TString> CSCGEMEff_N{"MatchedSH","AllLCT","MatchedLCT"};
  vector<TString> CSCGEMEff_V{"Eta","Phi","R","z","Slope"};
  vector<TString> CSCGEMEff_V2{"#eta","#phi","R","z","Slope"};
  vector<int> CSCGEMEff_V_div{60,70,80,220,20};
  vector<float> CSCGEMEff_V_low{-3.0,-3.5,0,-1100,-10};
  vector<float> CSCGEMEff_V_up{3.0,3.5,800,1100,10};
  vector<vector<vector<vector<vector<TEfficiency*> > > > >CSCGEMEff;

  CSCGEMEff.resize(LDisks.size());
  for (unsigned iED = 0; iED < LDisks.size(); ++iED) {
    CSCGEMEff[iED].resize(LRings.size());
    for (unsigned iER = 0; iER < LRings.size(); ++iER) {
      CSCGEMEff[iED][iER].resize(CSCGEMEff_V.size());
      for (unsigned iEV = 0; iEV < CSCGEMEff_V.size(); ++iEV) {
        CSCGEMEff[iED][iER][iEV].resize(CSCGEMEff_N.size());
        for (unsigned iEN = 0; iEN < CSCGEMEff_N.size(); ++iEN) {
          CSCGEMEff[iED][iER][iEV][iEN].resize(4);
          CSCGEMEff[iED][iER][iEV][iEN][0] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "GP_Either_" + LDisks[iED] + LRings[iER], "CSCGEM(Either) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEff[iED][iER][iEV][iEN][1] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "GP_1_" + LDisks[iED] + LRings[iER], "CSCGEM(Layer1) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEff[iED][iER][iEV][iEN][2] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "GP_2_" + LDisks[iED] + LRings[iER], "CSCGEM(Layer2) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          CSCGEMEff[iED][iER][iEV][iEN][3] = new TEfficiency("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "GP_Both_" + LDisks[iED] + LRings[iER], "CSCGEM(Both) Efficiency Vs " + CSCGEMEff_V2[iEV] + " for " + CSCGEMEff_N[iEN] + "LCTs at" + LDisks[iED] + " " + LRings[iER] + "; TP SimHit " + CSCGEMEff_V2[iEV] + "; Efficiency", CSCGEMEff_V_div[iEV], CSCGEMEff_V_low[iEV], CSCGEMEff_V_up[iEV]);
          for (auto eff_ : CSCGEMEff[iED][iER][iEV][iEN]) eff_->SetDirectory(out);
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



  // TP reconstruction plots
  //how often is a sim muon matched to ALCT, CLCT, 1GEM cluster, 2GEM clusters? How often do you get any combination? How often do you get an LCT?
  //What is the type of LCT matched to the muon? etc. questions like this need to be answered in order to understand any inefficiencies
  //because there are 5 types of LCTs we can make in ME1/1: ALCT-CLCT, ALCT-CLCT-1GEM, ALCT-CLCT-2GEM, ALCT-2GEM and CLCT-2GEM, we need to know exactly which type was expected to show up in a chamber based on the presence of ALCT/CLCT/GEM, but somehow did not.

  Long64_t nentries = tree->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    PrintProgress(jentry, nentries,1);
    // cout << "Finished event init " << jentry <<endl;
    Long64_t ientry = tree->LoadTree(jentry);
    tree->GetEntry(jentry);
    // cout << "Finished Loading event" <<endl;
    tr->ReadTree();
    vector<tp>& tps = tr->Evt.MuonTPs;
    for (tp thistp : tps) {
      if (thistp.eta > 0) evttppos++;
      else evttpneg++;
    }
    // cout << "READ" <<endl;
    vector<vector<StationData> >& Stations = tr->Evt.Stations;

    for (unsigned disk = 0; disk < Stations.size(); ++disk) for (unsigned ring = 0; ring < Stations[disk].size(); ++ring) {
      if (disk > 1) break; // Considering only first 2 disks;
      if (ring >= LRings.size()) continue;
      StationData& ThisStation = Stations[disk][ring];
      for (unsigned ittp = 0; ittp < ThisStation.TPInfos.size(); ++ittp) {
        TPContent& tp = ThisStation.TPInfos[ittp];
        const vector<float> CSCSimHitAve_V{tp.CSCSimHitAve.eta, tp.CSCSimHitAve.phi, tp.CSCSimHitAve.r, tp.CSCSimHitAve.z};
        const vector<float> GEMSimHitAve_V{tp.GEMSimHitAve.eta, tp.GEMSimHitAve.phi, tp.GEMSimHitAve.r, tp.GEMSimHitAve.z};
        bool CanRecoCSC = (tp.NSimHitsCSC > 3);
        bool CanRecoGEM = (tp.NSimHitsGEM > 0);
        if (!CanRecoCSC) {
          cout << "This TP cannot reconstruct CSCStubs, NSimHitsCSC = " << tp.NSimHitsCSC << ", tp eta = " << tp.TP.eta << ", SimHitAveEta = " << tp.CSCSimHitAve.eta << endl;
        }
        if (!CanRecoGEM) {
          cout << "This TP cannot reconstruct GEMDigis, NSimHitsGEM = " << tp.NSimHitsGEM << ", tp eta = " << tp.TP.eta << ", SimHitAveEta = " << tp.GEMSimHitAve.eta << endl;
        }
        if (tp.TP.eta > 0) tppos++;
        else if (tp.TP.eta < 0) tpneg++;
        if (CanRecoCSC) {
          if (tp.CSCSimHitAve.eta > 0) cancscpos++;
          else if (tp.CSCSimHitAve.eta < 0) cancscneg++;
          bool DigiInMatch = tp.CSCStubs.size();
          // TH2Plots[1][0][disk][ring]->Fill(CSCSimHitAve_V[3],CSCSimHitAve_V[2]);
          // TH2Plots[2][0][disk][ring]->Fill(CSCSimHitAve_V[1],CSCSimHitAve_V[2]);
          if (DigiInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[0][iEV][disk]->Fill(1,CSCSimHitAve_V[iEV]); // CSCReco
              DetEff[1][iEV][disk]->Fill(1,CSCSimHitAve_V[iEV]); // CSCMatch
            }
          }
          else {
            bool DigiInGeneric = false;
            for (unsigned idigi = 0; idigi < ThisStation.CSCStubs.size(); ++idigi) {
              CSCStub CSCDigi = ThisStation.CSCStubs[idigi];
              DigiInGeneric = IsCloseCSC(CSCDigi.eta, CSCSimHitAve_V[0], CSCDigi.phi, CSCSimHitAve_V[1]);
              if (DigiInGeneric) break;
            }

            if (DigiInGeneric) {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[0][iEV][disk]->Fill(1,CSCSimHitAve_V[iEV]); // CSCReco
                DetEff[1][iEV][disk]->Fill(0,CSCSimHitAve_V[iEV]); // CSCMatch
              }
            }
            else {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[0][iEV][disk]->Fill(0,CSCSimHitAve_V[iEV]); // CSCReco
              }
            }
          }
          float ClosedR = 999;
          float ClosedEta(999), ClosedPhi(999);
          vector<CSCStub>& DigisToLook = (LookForDigiInMatched ? tp.CSCStubs : ThisStation.CSCStubs);
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            CSCStub CSCDigi = DigisToLook[idigi];
            if (!(CSCDigi.valid)) continue;
            // TH2Plots[3][0][disk][ring]->Fill(CSCDigi.z,CSCDigi.r);
            // TH2Plots[4][0][disk][ring]->Fill(CSCDigi.phi,CSCDigi.r);
            vector<float> delta = CalcdR(CSCDigi.eta, CSCSimHitAve_V[0], CSCDigi.phi, CSCSimHitAve_V[1]);
            if (delta[0] < ClosedR) {
              ClosedR = delta[0];
              ClosedEta = delta[1];
              ClosedPhi = delta[2];
            }
          }
          if (ClosedR < 999) TH2Plots[0][0][disk][ring]->Fill(ClosedEta,ClosedPhi);
        }

        if (CanRecoGEM) {
          if (tp.CSCSimHitAve.eta > 0) cangempos++;
          else if (tp.CSCSimHitAve.eta < 0) cangemneg++;
          bool DigiInMatch = tp.GEMDigis.size();
          if (DigiInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[2][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMReco
              DetEff[3][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMMatch
            }
          }
          else {
            bool DigiInGeneric = false;
            for (unsigned idigi = 0; idigi < ThisStation.GEMDigis.size(); ++idigi) {
              GEMDigi GEMDigi = ThisStation.GEMDigis[idigi];
              DigiInGeneric = IsCloseGEM(GEMDigi.eta, GEMSimHitAve_V[0], GEMDigi.phi, GEMSimHitAve_V[1]);
              if (DigiInGeneric) break;
            }

            if (DigiInGeneric) {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[2][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // GEMReco
                DetEff[3][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMMatch
              }
            }
            else {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[2][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // GEMReco
              }
            }
          }
          float ClosedR = 999;
          float ClosedEta(999), ClosedPhi(999);
          vector<GEMDigi>& DigisToLook = (LookForDigiInMatched ? tp.GEMDigis : ThisStation.GEMDigis);
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            GEMDigi GEMDigi = DigisToLook[idigi];
            vector<float> delta = CalcdR(GEMDigi.eta, GEMSimHitAve_V[0], GEMDigi.phi, GEMSimHitAve_V[1]);
            if (delta[0] < ClosedR) {
              ClosedR = delta[0];
              ClosedEta = delta[1];
              ClosedPhi = delta[2];
            }
          }
          if (ClosedR < 999) TH2Plots[0][1][disk][ring]->Fill(ClosedEta,ClosedPhi);
        }

        if (CanRecoGEM) {
          bool ClusterInMatch = tp.Clusters.size();
          if (ClusterInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[4][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]);
              DetEff[5][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]);
            }
          }
          else {
            bool ClusterInGeneric = false;
            for (unsigned ic = 0; ic < ThisStation.Clusters.size(); ++ic) {
              GEMPadDigiCluster cl = ThisStation.Clusters[ic];
              ClusterInGeneric = IsCloseCluster(cl.eta, GEMSimHitAve_V[0], cl.phi, GEMSimHitAve_V[1]);
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
          vector<GEMPadDigiCluster>& DigisToLook = (LookForDigiInMatched ? tp.Clusters : ThisStation.Clusters);
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            GEMPadDigiCluster cl = DigisToLook[idigi];
            vector<float> delta = CalcdR(cl.eta, GEMSimHitAve_V[0], cl.phi, GEMSimHitAve_V[1]);
            if (delta[0] < ClosedR) {
              ClosedR = delta[0];
              ClosedEta = delta[1];
              ClosedPhi = delta[2];
            }
          }
          if (ClosedR < 999) TH2Plots[0][2][disk][ring]->Fill(ClosedEta,ClosedPhi);
          // if (ClosedR < 999) cout << Form("dR = %f, dEta = %f, dPhi = %f", ClosedR, ClosedEta, ClosedPhi);
        }

        if (CanRecoCSC && CanRecoGEM) {
          if (tp.CSCSimHitAve.eta > 0) canbothpos++;
          else if (tp.CSCSimHitAve.eta < 0) canbothneg++;
          vector<float> CSCSimHitAve_VExt = CSCSimHitAve_V;
          CSCSimHitAve_VExt.push_back(0);
          for (unsigned icsc = 0; icsc < tp.CSCStubs.size(); ++icsc) {
            CSCStub lct = tp.CSCStubs[icsc];
            bool hasGEM1 = (lct.GEM1pad != 255);
            bool hasGEM2 = (lct.GEM2pad != 255);
            CSCSimHitAve_VExt[4] = lct.slope;
            vector<float> CSCStub_V{lct.eta, lct.phi, lct.r, lct.z, (float) lct.slope};
            for (unsigned iEV = 0; iEV < CSCGEMEff_V.size(); ++iEV) {
              CSCGEMEff[disk][ring][iEV][0][0]->Fill((hasGEM1 || hasGEM2), CSCSimHitAve_VExt[iEV]);
              CSCGEMEff[disk][ring][iEV][0][1]->Fill((hasGEM1), CSCSimHitAve_VExt[iEV]);
              CSCGEMEff[disk][ring][iEV][0][2]->Fill((hasGEM2), CSCSimHitAve_VExt[iEV]);
              CSCGEMEff[disk][ring][iEV][0][3]->Fill((hasGEM1 && hasGEM2), CSCSimHitAve_VExt[iEV]);
              CSCGEMEff[disk][ring][iEV][2][0]->Fill((hasGEM1 || hasGEM2), CSCStub_V[iEV]);
              CSCGEMEff[disk][ring][iEV][2][1]->Fill((hasGEM1), CSCStub_V[iEV]);
              CSCGEMEff[disk][ring][iEV][2][2]->Fill((hasGEM2), CSCStub_V[iEV]);
              CSCGEMEff[disk][ring][iEV][2][3]->Fill((hasGEM1 && hasGEM2), CSCStub_V[iEV]);
            }
          }
        }

        if (CanRecoCSC && CanRecoGEM && tp.CSCStubs.size() && tp.GEMDigis.size()) {
          for (unsigned icsc = 0; icsc < tp.CSCStubs.size(); ++icsc) {
            float slope = tp.CSCStubs[icsc].slope;
            float zsign = (tp.CSCStubs[icsc].z > 0) - (tp.CSCStubs[icsc].z < 0);
            for (unsigned igem = 0; igem < tp.GEMDigis.size(); ++igem) {
              float dphi = fabs(tp.GEMDigis[igem].phi - tp.CSCStubs[icsc].phi);
              while (dphi > 2 * M_PI) dphi -= 2 * M_PI;
              dPhiVsSlope[disk][ring]->Fill(slope,dphi);
            }
          }
        }
      } // End of tp loop


      for (unsigned icsc = 0; icsc < ThisStation.CSCStubs.size(); ++icsc) {
        CSCStub lct = ThisStation.CSCStubs[icsc];
        bool hasGEM1 = (lct.GEM1pad != 255);
        bool hasGEM2 = (lct.GEM2pad != 255);
        vector<float> CSCStub_V{lct.eta, lct.phi, lct.r, lct.z, (float)lct.slope};
        for (unsigned iEV = 0; iEV < CSCGEMEff_V.size(); ++iEV) {
          CSCGEMEff[disk][ring][iEV][1][0]->Fill((hasGEM1 || hasGEM2), CSCStub_V[iEV]);
          CSCGEMEff[disk][ring][iEV][1][1]->Fill((hasGEM1), CSCStub_V[iEV]);
          CSCGEMEff[disk][ring][iEV][1][2]->Fill((hasGEM2), CSCStub_V[iEV]);
          CSCGEMEff[disk][ring][iEV][1][3]->Fill((hasGEM1 && hasGEM2), CSCStub_V[iEV]);
        }
      }
    } // End of detector loop
  } // End of event loop

  cout << "Can CSC Pos = " << cancscpos << ", Neg = " << cancscneg<< endl;
  cout << "Can GEM Pos = " << cangempos << ", Neg = " << cangemneg<< endl;
  cout << "Can Both Pos = " << canbothpos << ", Neg = " << canbothneg<< endl;
  cout << "TP in +eta : " << tppos << ", TP in -eta : " << tpneg << endl;
  cout << "Event TP in +eta : " << evttppos << ", TP in -eta : " << evttpneg << endl;
  out->Write();
  out->Save();
  out->Close();
}
