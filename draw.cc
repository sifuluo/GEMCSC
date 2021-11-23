{
  TString savepath = "../../Report/plots/";
  TFile *f = new TFile("plots.root");
  TCanvas *c6 = new TCanvas("c6","c6",1700,1100);
  c6->Divide(3,2);
  TCanvas *c4 = new TCanvas("c4","c4",1700,1000);
  c4->Divide(2,2);
  TCanvas *c3 = new TCanvas("c3","c3",1700,500);
  c3->Divide(3,1);
  TCanvas *c2 = new TCanvas("c2","c2",1100,500);
  c2->Divide(2,1);
  TCanvas *c1 = new TCanvas("c1","c1",600,500);

  vector<TString> LDisks{"D0","D1","D2"};
  // vector<TString> LRings{"R1","R2","R3","R0"};
  vector<TString> LRings{"R1"};
  vector<TString> Eff_T{"CSCReco","CSCMatch","GEMReco","GEMMatch","ClusterReco","ClusterMatch","PadReco","PadMatch"};
  vector<TString> Eff_V{"Eta","Phi","R","z"};
  for (unsigned iET = 0; iET < Eff_T.size(); ++iET) {
    for (unsigned iES = 0; iES < LDisks.size(); ++iES) {
      for (unsigned iEV = 0; iEV < Eff_V.size(); ++iEV) {
        TEfficiency *eff = (TEfficiency*) f->Get(Eff_T[iET] + "EffVs" + Eff_V[iEV] + "_" + LDisks[iES]);
        c3->cd(iEV+1);
        eff->Draw();
      }
      c3->SaveAs(savepath + Eff_T[iET] + "Eff_" + LDisks[iES] + ".png");
    }
  }

  vector<TString> Close_T{"CSC","GEM","Cluster","Pad"};
  for (unsigned iCT = 0; iCT < Close_T.size(); ++iCT) {
    for (unsigned iCD = 0; iCD < LDisks.size(); ++iCD) {
      for(unsigned iCR = 0; iCR < LRings.size(); ++iCR) {
        TH2F* close = (TH2F*) f->Get(Close_T[iCT] + "ClosestDigi_" + LDisks[iCD] + LRings[iCR]);
        c1->cd();
        close->Rebin2D(5,5); // Original bin width is 0.01
        // close->GetYaxis()->SetRangeUser(0,0.012);
        close->GetYaxis()->SetTitleOffset(1.4);
        close->Draw("colz");
        if (close->GetMean(1) != 0) {
          cout << endl << close->GetName() << " : " <<endl;
          cout << "    eta: " << close->GetMean(1) << " +- " << close->GetRMS(1) << " ; mu + 2 stddev = " << close->GetMean(1) + 2 * close->GetRMS(1) <<endl;
          cout << "    phi: " << close->GetMean(2) << " +- " << close->GetRMS(2) << " ; mu + 2 stddev = " << close->GetMean(2) + 2 * close->GetRMS(2) <<endl;
        }
        c1->SaveAs(savepath + Close_T[iCT] + "ClosestDigi_" + LDisks[iCD] + LRings[iCR] + ".png");
      }
    }
  }
  vector<TString> CSCGEMEff_N{"Matched","All"};
  vector<TString> CSCGEMEff_V{"Eta","Phi","R","z","Slope","Chamber"};
  for (unsigned iED = 0; iED < LDisks.size(); ++iED) {
    for (unsigned iER = 0; iER < LRings.size(); ++iER) {
      for (unsigned iEV = 0; iEV < CSCGEMEff_V.size(); ++iEV) {
        for (unsigned iEN = 0; iEN < CSCGEMEff_N.size(); ++iEN) {
          TEfficiency *eff0 = (TEfficiency*) f->Get("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_Either_" + LDisks[iED] + LRings[iER]);
          TEfficiency *eff1 = (TEfficiency*) f->Get("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_1_" + LDisks[iED] + LRings[iER]);
          TEfficiency *eff2 = (TEfficiency*) f->Get("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_2_" + LDisks[iED] + LRings[iER]);
          TEfficiency *eff3 = (TEfficiency*) f->Get("CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_Both_" + LDisks[iED] + LRings[iER]);
          vector<TEfficiency*> effs{eff0, eff1, eff2, eff3};
          for (unsigned i = 0; i < effs.size(); ++i) {
            c4->cd(i+1);
            effs[i]->Draw();
          }
          c4->SaveAs(savepath + "CSCGEMEffVs" + CSCGEMEff_V[iEV] + CSCGEMEff_N[iEN] + "_" + LDisks[iED] + LRings[iER]+".png");
        }
      }
    }
  }
  for (unsigned iED = 0; iED < LDisks.size(); ++iED) {
    for (unsigned iER = 0; iER < LRings.size(); ++iER) {
      TH2F *h2 = (TH2F*) f->Get("dPhiVsSlope_" + LDisks[iED] + LRings[iER]);
      c1->cd();
      h2->Draw("colz");
      c1->SaveAs(savepath + "dPhiVsSlope" + LDisks[iED] + LRings[iER] + ".png");
    }
  }




  gApplication->Terminate();
}
