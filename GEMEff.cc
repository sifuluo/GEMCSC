#include <iostream>
#include <vector>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

#include "TreeReader.cc"
#include "TreeDataFormat.cc"
#include "Tools.cc"


using namespace std;

void GEMEff() {
  TChain *tree = new TChain("NtupleMaker/eventTree");
  tree->Add("out_Run3.root");

  // TFile *fin = new TFile("out_Run3.root");
  // TTree *tree = (TTree*) fin->Get("NtupleMaker/eventTree");



  TreeReader tr = TreeReader(tree);

  // TFile *out = new TFile("plots.root","RECREATE");

  Long64_t nentries = tree->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    cout << "Finished event init " << jentry <<endl;
    Long64_t ientry = tree->LoadTree(jentry);
    // cout << "Finished Loading Tree" <<endl;
    tree->GetEntry(jentry);
    cout << "Finished Loading event" <<endl;
    // vector<float>* z;
    // tree->SetBranchAddress("allCscStubsLCT_z",&z);
    // if (z->size() == 0) cout << "empty z" <<endl;
    // else cout << "z[0] = " << z->at(0) <<endl;
    tr.ReadTree();
    cout << "Finished ReadTree" <<endl;
    // cout << tr.Evt.Stations[0][0].TPInfos.size() <<endl;
  }


}
