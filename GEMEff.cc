#include <iostream>
#include <vector>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "TreeReader.cc"
#include "TreeDataFormat.cc"
#include "Tools.cc"


using namespace std;

void GEMEff() {
  TChain *tree = new TChain("NtupleMaker/eventTree");
  tree->Add("out_Run3_2000.root");

  TreeReader tr = TreeReader(tree);
  TFile *out = new TFile("plots.root","RECREATE");

  // Lists of Disks and Rings
  vector<TString> LDisks{"D1","D2"};
  vector<TString> LRings{"R1","R2","R3"};

  // Efficiency plots
  vector<TString> Eff_T{"CSCReco","CSCMatch","GEMReco","GEMMatch"};
  vector<TString> Eff_V{"Eta","Phi","R"};
  vector<int> Eff_V_div{60,70,80};
  vector<double> Eff_V_low{-3.0,-3.5,0};
  vector<double> Eff_V_up{3.0,3.5,800};
  vector<vector<vector <TEfficiency*> > > DetEff;
  DetEff.resize(Eff_T.size());
  for (unsigned iET = 0; iET < Eff_T.size(); ++iET) {
    DetEff[iET].resize(Eff_V.size());
    for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
      DetEff[iET][iEV].resize(LDisks.size());
      for (unsigned iES = 0; iES < LDisks.size(); ++iES) {
        DetEff[iET][iEV][iES] = new TEfficiency(Eff_T[iET] + "EffVs" + Eff_V[iEV] + "_" + LDisks[iES], Eff_T[iET] + " Efficiency Vs TP " + Eff_V[iEV] + " at " + LDisks[iES] + "; TP " + Eff_V[iEV] + "; Efficiency", Eff_V_div[iEV], Eff_V_low[iEV], Eff_V_up[iEV]);
      }
    }
  }

  // Closest Digi plots
  bool LookForDigiInMatched = true;
  vector<TString> Close_T{"CSC","GEM"};
  vector<int> Close_V_div{   100 , 50 , 60, 80, 70 , 110};//dEta,dPhi,Eta,Phi,R,Z
  vector<double> Close_V_low{0  , 0  , -3, -4, 0  , 0};
  vector<double> Close_V_up{ 0.1, 0.05, 3 , 4 , 700, 1100};
  vector<vector<vector<vector<TH2F*> > > > TH2Plots;
  TH2Plots.resize(5);
  for (auto& vec : TH2Plots) vec.resize(Close_T.size());
  // DetClose.resize(Close_T.size());
  for (unsigned iCT = 0; iCT < Close_T.size(); ++iCT) {
    for (auto& vec : TH2Plots) vec[iCT].resize(LDisks.size());
    // DetClose[iCT].resize(LDisks.size());
    for (unsigned iCS = 0; iCS < LDisks.size(); ++iCS) {
      for (auto& vec : TH2Plots) vec[iCT][iCS].resize(LRings.size());
      // DetClose[iCT][iCS].resize(LRings.size());
      for (unsigned iCR = 0; iCR < LRings.size(); ++iCR) {
        TH2Plots[0][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "ClosestDigi_" + LDisks[iCS] + LRings[iCR], "Closest Digi in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; Delta Eta; Delta Phi", Close_V_div[0], Close_V_low[0], Close_V_up[0], Close_V_div[1], Close_V_low[1], Close_V_up[1] );
        TH2Plots[1][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "SimHitRVZ_" + LDisks[iCS] + LRings[iCR], "SimHitRVZ in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; z; r", Close_V_div[5], Close_V_low[5], Close_V_up[5], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
        TH2Plots[2][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "SimHitRVPhi_" + LDisks[iCS] + LRings[iCR], "SimHitRVPhi in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; phi; r", Close_V_div[3], Close_V_low[3], Close_V_up[3], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
        TH2Plots[3][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "DigiRVZ_" + LDisks[iCS] + LRings[iCR], "DigiRVZ in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; z; r", Close_V_div[5], Close_V_low[5], Close_V_up[5], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
        TH2Plots[4][iCT][iCS][iCR] = new TH2F(Close_T[iCT] + "DigiRVPhi_" + LDisks[iCS] + LRings[iCR], "DigiRVPhi in " + Close_T[iCT] + " at " + LDisks[iCS] + " " + LRings[iCR] + "; phi; r", Close_V_div[3], Close_V_low[3], Close_V_up[3], Close_V_div[4], Close_V_low[4], Close_V_up[4] );
      }
    }
  }

  Long64_t nentries = tree->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    // cout << "Finished event init " << jentry <<endl;
    Long64_t ientry = tree->LoadTree(jentry);
    tree->GetEntry(jentry);
    // cout << "Finished Loading event" <<endl;
    tr.ReadTree();
    vector<vector<StationData> >& Stations = tr.Evt.Stations;

    for (unsigned disk = 0; disk < Stations.size(); ++disk) for (unsigned ring = 0; ring < Stations[disk].size(); ++ring) {
      if (disk > 1) break; // Considering only first 2 disks;
      StationData& ThisStation = Stations[disk][ring];
      for (unsigned ittp = 0; ittp < ThisStation.TPInfos.size(); ++ittp) {
        TPContent& tp = ThisStation.TPInfos[ittp];
        const vector<double> CSCSimHitAve_V{tp.CSCSimHitAve.eta, tp.CSCSimHitAve.phi, tp.CSCSimHitAve.r, tp.CSCSimHitAve.z};
        const vector<double> GEMSimHitAve_V{tp.GEMSimHitAve.eta, tp.GEMSimHitAve.phi, tp.GEMSimHitAve.r, tp.GEMSimHitAve.z};
        bool CanRecoCSC = (tp.NSimHitsCSC > 3);
        bool CanRecoGEM = (tp.NSimHitsGEM > 0);

        if (CanRecoCSC) {
          bool DigiInMatch = tp.CSCStubs.size();
          TH2Plots[1][0][disk][ring]->Fill(CSCSimHitAve_V[3],CSCSimHitAve_V[2]);
          TH2Plots[2][0][disk][ring]->Fill(CSCSimHitAve_V[1],CSCSimHitAve_V[2]);
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
              double etadiff = fabs(CSCDigi.eta - CSCSimHitAve_V[0]);
              double phidiff = TVector2::Phi_mpi_pi(CSCDigi.phi - CSCSimHitAve_V[1]);
              double deltar = sqrt(etadiff * etadiff + phidiff * phidiff);
              if (deltar < 0.2) DigiInGeneric = true;
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
          double ClosedR = 999;
          double ClosedEta, ClosedPhi;
          CSCStub CloseDigi;
          vector<CSCStub>& DigisToLook = (LookForDigiInMatched ? tp.CSCStubs : ThisStation.CSCStubs);
          // for (unsigned idigi = 0; idigi < ThisStation.CSCStubs.size(); ++idigi) {
          //   CSCStub CSCDigi = ThisStation.CSCStubs[idigi];
          // for (unsigned idigi = 0; idigi < tp.CSCStubs.size(); ++idigi) {
          //   CSCStub CSCDigi = tp.CSCStubs[idigi];
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            CSCStub CSCDigi = DigisToLook[idigi];
            if (!(CSCDigi.valid)) continue;
            TH2Plots[3][0][disk][ring]->Fill(CSCDigi.z,CSCDigi.r);
            TH2Plots[4][0][disk][ring]->Fill(CSCDigi.phi,CSCDigi.r);
            double etadiff = fabs(CSCDigi.eta - CSCSimHitAve_V[0]);
            double phidiff = TVector2::Phi_mpi_pi(CSCDigi.phi - CSCSimHitAve_V[1]);
            double deltar = sqrt(etadiff * etadiff + phidiff * phidiff);
            if (deltar < ClosedR) {
              CloseDigi = CSCDigi;
              ClosedR = deltar;
              ClosedEta = etadiff;
              ClosedPhi = phidiff;
            }
          }
          if (ClosedR < 999) TH2Plots[0][0][disk][ring]->Fill(ClosedEta,ClosedPhi);
          // if (ClosedPhi < 0.2 && disk == 0) {
          //   bool IsMatched = IsCSCStubMatched(CloseDigi, MatchedCSCstubs, tp);
          //   double IsMatched_d = (IsMatched ? 1. : 0.);
          //   double dZ = CloseDigi.z - CSCSimHitAve_V[3];
          //   if (ClosedPhi < 0.1) { // Close piece
          //     // CSCSta1[0][0]->Fill(IsMatched_d, dZ);
          //   }
          //   else { // Far piece
          //     // CSCSta1[1][0]->Fill(IsMatched_d, dZ);
          //   }
          // }
        }

        if (CanRecoGEM) {
          bool DigiInMatch = tp.GEMDigis.size();
          if (DigiInMatch) {
            for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
              DetEff[2][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // CSCReco
              DetEff[3][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // CSCMatch
            }
          }
          else {
            bool DigiInGeneric = false;
            for (unsigned idigi = 0; idigi < ThisStation.GEMDigis.size(); ++idigi) {
              GEMDigi GEMDigi = ThisStation.GEMDigis[idigi];
              double etadiff = fabs(GEMDigi.eta - GEMSimHitAve_V[0]);
              double phidiff = TVector2::Phi_mpi_pi(GEMDigi.phi - GEMSimHitAve_V[1]);
              double deltar = sqrt(etadiff * etadiff + phidiff * phidiff);
              if (deltar < 0.2) DigiInGeneric = true;
            }

            if (DigiInGeneric) {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[2][iEV][disk]->Fill(1,GEMSimHitAve_V[iEV]); // CSCReco
                DetEff[3][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // CSCMatch
              }
            }
            else {
              for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
                DetEff[2][iEV][disk]->Fill(0,GEMSimHitAve_V[iEV]); // CSCReco
              }
            }
          }
          double ClosedR = 999;
          double ClosedEta, ClosedPhi;
          vector<GEMDigi>& DigisToLook = (LookForDigiInMatched ? tp.GEMDigis : ThisStation.GEMDigis);
          // for (unsigned idigi = 0; idigi < ThisStation.GEMDigis.size(); ++idigi) {
          //   GEMDigi GEMDigi = ThisStation.GEMDigis[idigi];
          // for (unsigned idigi = 0; idigi < tp.GEMDigis.size(); ++idigi) {
          //   GEMDigi GEMDigi = tp.GEMDigis[idigi];
          for (unsigned idigi = 0; idigi < DigisToLook.size(); ++idigi) {
            GEMDigi GEMDigi = DigisToLook[idigi];
            double etadiff = fabs(GEMDigi.eta - GEMSimHitAve_V[0]);
            double phidiff = TVector2::Phi_mpi_pi(GEMDigi.phi - GEMSimHitAve_V[1]);
            double deltar = sqrt(etadiff * etadiff + phidiff * phidiff);
            if (deltar < ClosedR) {
              ClosedR = deltar;
              ClosedEta = etadiff;
              ClosedPhi = phidiff;
            }
          }
          if (ClosedR < 999) TH2Plots[0][1][disk][ring]->Fill(ClosedEta,ClosedPhi);
        }
      }
    }
  }
  out->Write();
  out->Save();
  out->Close();

}
