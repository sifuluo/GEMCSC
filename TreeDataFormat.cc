// Data Containers for GEM-CSC
#ifndef TREEDATAFORMAT_CC
#define TREEDATAFORMAT_CC
#include "TROOT.h"
#include "TString.h"
#include "TVector2.h"
#include "TTree.h"
#include "TChain.h"

#include <vector>
#include <utility>
#include <cmath>

#include "Tools.cc"

using namespace std;
//bundles all tp-related information
struct tp{
  float    pt;
  float    eta;
  float    phi;
  float    dxy;
  float    d0;
  float    z0;
  float    d0_prod;
  float    z0_prod;
  int      pdgid;
  int      eventid; //0 for main event, >0 for PU - only available in specific samples
  int      charge;
  unsigned Index; //contains the oriignal unculled tp collection index to facilitate linking
};

//container for SimHit-related information
struct SimHit{
  float    phi;
  float    eta;
  float    r;
  float    z;
  unsigned MatchTp;
};

struct GEMPadDigi{
  float phi;
  float eta;
  float r;
  float z;
};

//container for CSC stub information
struct CSCStub{
  float phi;
  float eta;
  float r;
  float z;
  int   bend;
  int   slope;
  int   strip8;
  int   pattern;
  bool  valid;
  int   MatchTp;
  vector<int> CLCT_hits;
  vector<int> CLCT_positions;
  vector<int> ALCT_hits;
  int layer; //runs from 0 innermost to x outermost, with x=2 for ME2-4 and x=4 for ME1, mind the single layer of ME1/3 being x=3
  GEMPadDigi GEM1;
  GEMPadDigi GEM2;
};

//container for GEM digi information
struct GEMDigi{
  float phi;
  float eta;
  float r;
  float z;
  int   MatchTp;
  int   layer; //0 is innermost, 3 is outermost, given the subdetector in question
};

//container for SimHit averages
struct SimHitAve{
  float phi;
  float eta;
  float r;
  float z;
};

class TPContent{
public:
  // idx = -1 by default
  TPContent(unsigned idx) {
    SortedIndex = RawIndex = idx;
    CSCSimHits.clear();
    GEMSimHits.clear();
    CSCStubs.clear();
    GEMDigis.clear();
  };
  void CalcSimHitAve() {
    float ncsc = NSimHitsCSC = CSCSimHits.size();
    float ngem = NSimHitsGEM = GEMSimHits.size();
    CSCSimHitAve = GEMSimHitAve = SimHitAve{0.,0.,0.,0.};
    for (unsigned i = 0; i < NSimHitsCSC; ++i) {
      CSCSimHitAve.phi += CSCSimHits[i].phi / ncsc;
      CSCSimHitAve.eta += CSCSimHits[i].eta / ncsc;
      CSCSimHitAve.r   += CSCSimHits[i].r / ncsc;
      CSCSimHitAve.z   += CSCSimHits[i].z / ncsc;
    };
    CSCSimHitAve.phi = TVector2::Phi_mpi_pi(CSCSimHitAve.phi);
    for (unsigned i = 0; i < NSimHitsGEM; ++i) {
      GEMSimHitAve.phi += GEMSimHits[i].phi / ngem;
      GEMSimHitAve.eta += GEMSimHits[i].eta / ngem;
      GEMSimHitAve.r   += GEMSimHits[i].r / ngem;
      GEMSimHitAve.z   += GEMSimHits[i].z / ngem;
    };
    GEMSimHitAve.phi = TVector2::Phi_mpi_pi(GEMSimHitAve.phi);
  }

  tp               TP;
  unsigned         RawIndex;
  unsigned         SortedIndex;
  unsigned         NSimHitsCSC;
  unsigned         NSimHitsGEM;
  SimHitAve        CSCSimHitAve;
  SimHitAve        GEMSimHitAve;
  vector<SimHit>   CSCSimHits;
  vector<SimHit>   GEMSimHits;
  vector<CSCStub>  CSCStubs;
  vector<GEMDigi>  GEMDigis;
};

class StationData{
public:
  StationData(TString lable) {
    StationLabel = lable;
    TPInfos.clear();
    CSCStubs.clear();
    GEMDigis.clear();
    GEMPadDigis.clear();
  };
  void AddCSCSimHit(SimHit& d) {
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].CSCSimHits.push_back(d);
      return;
    }
    TPContent tmp((unsigned)d.MatchTp);
    tmp.CSCSimHits.push_back(d);
    TPInfos.push_back(tmp);
  }
  void AddGEMSimHit(SimHit& d) {
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].GEMSimHits.push_back(d);
      return;
    }
    TPContent tmp((unsigned)d.MatchTp);
    tmp.GEMSimHits.push_back(d);
    TPInfos.push_back(tmp);
  }
  void AddCSCStub(CSCStub& d) {
    if (d.MatchTp == -1) {
      CSCStubs.push_back(d);
      return;
    }
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].CSCStubs.push_back(d);
      return;
    }
    TPContent tmp((unsigned)d.MatchTp);
    tmp.CSCStubs.push_back(d);
    TPInfos.push_back(tmp);
  }
  void AddGEMDigi(GEMDigi& d) {
    if (d.MatchTp == -1) {
      GEMDigis.push_back(d);
      return;
    }
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].GEMDigis.push_back(d);
      return;
    }
    TPContent tmp((unsigned)d.MatchTp);
    tmp.GEMDigis.push_back(d);
    TPInfos.push_back(tmp);
  }
  void AddGEMPadDigi(GEMPadDigi& d) {
    GEMPadDigis.push_back(d);
  }
  void FillTP(vector<tp> tps) {
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      bool found = false;
      for (unsigned j = 0; j < tps.size(); ++i) {
        if (TPInfos[i].RawIndex == tps[j].Index){
          TPInfos[i].TP = tps[j];
          TPInfos[i].TP.pt = tps[j].pt;
          TPInfos[i].TP.eta = tps[j].eta;
          TPInfos[i].TP.phi = tps[j].phi;
          TPInfos[i].TP.dxy = tps[j].dxy;
          TPInfos[i].TP.d0 = tps[j].d0;
          TPInfos[i].TP.z0 = tps[j].z0;
          TPInfos[i].TP.d0_prod = tps[j].d0_prod;
          TPInfos[i].TP.z0_prod = tps[j].z0_prod;
          TPInfos[i].TP.pdgid = tps[j].pdgid;
          TPInfos[i].TP.eventid = tps[j].eventid;
          TPInfos[i].TP.charge = tps[j].charge;
          TPInfos[i].TP.Index = tps[j].Index;
          found = true;
          break;
        }
      }
      if (!found) cout << "Cannot find a tp with index = " << TPInfos[i].RawIndex <<endl;
    }
  }
  void CalcSimHitAve() {
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      TPInfos[i].CalcSimHitAve();
    }
  }
  TString            StationLabel;
  vector<TPContent>  TPInfos;
  vector<CSCStub>    CSCStubs;
  vector<GEMDigi>    GEMDigis;
  vector<GEMPadDigi> GEMPadDigis;
};

class EventData{
public:
  EventData() {
    Init();
  }
  void Init() {
    Stations = vector<vector<StationData> >{
      {
        StationData("1-1"),
        StationData("1-2"),
        StationData("1-3")
      },
      {
        StationData("2-1"),
        StationData("2-2")
      },
      {
        StationData("3-1"),
        StationData("3-2")
      },
      {
        StationData("4-1"),
        StationData("4-2")
      }
    };
    MuonTPs.clear();
  }
  StationData& Station(pair<unsigned,unsigned> DaR){
    return Stations[DaR.first][DaR.second];
  }
  void AddTP(tp& muontp) {
    MuonTPs.push_back(muontp);
  }
  void AddCSCSimHit(SimHit& d) {
    pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
    Stations[DaR.first][DaR.second].AddCSCSimHit(d);
  }
  void AddGEMSimHit(SimHit& d) {
    pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
    Stations[DaR.first][DaR.second].AddGEMSimHit(d);
  }
  void AddCSCStub(CSCStub& d) {
    pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
    Stations[DaR.first][DaR.second].AddCSCStub(d);
  }
  void AddGEMDigi(GEMDigi& d) {
    pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
    Stations[DaR.first][DaR.second].AddGEMDigi(d);
  }
  void AddGEMPadDigi(GEMPadDigi& d) {
    pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
    Stations[DaR.first][DaR.second].AddGEMPadDigi(d);
  }
  void FillTP() {
    for (unsigned i = 0; i < Stations.size(); ++i) for(unsigned j = 0; j < Stations[i].size(); ++j) {
      Stations[i][j].FillTP(MuonTPs);
    }
  }
  void CalcSimHitAve() {
    for (unsigned i = 0; i < Stations.size(); ++i) for(unsigned j = 0; j < Stations[i].size(); ++j) {
      Stations[i][j].CalcSimHitAve();
    }
  }

  vector<vector<StationData> > Stations;
  vector<tp> MuonTPs;
};

//container for tp-based information
struct TPmaster{
  unsigned tpIndex;
  vector<vector<bool> > InStation;//default pattern is [disk][ring] for addressing subdetectors; GEM and CSC stations are considered a single entity
};

class Branch_Reader{
public:
  Branch_Reader() {

  };

  void Init(TChain* evttree, TString name_, int data_type_, bool match_ = false) {
    name = name_;
    data_type = data_type_; // 0 for LCT, 1 for ALCT, 2 for CLCT, 3 for GemDigi
    IsMatched = match_;

    // LCT
    if (data_type == 0) {
      InitGP(evttree);
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
      quality = new std::vector<int>;
      detId = new std::vector<int>;
      keywire = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      valid = new std::vector<bool>;
      type = new std::vector<int>;
      evttree->SetBranchAddress(name+"_bend", &bend);
      evttree->SetBranchAddress(name+"_pattern", &pattern);
      evttree->SetBranchAddress(name+"_slope", &slope);
      evttree->SetBranchAddress(name+"_quality", &quality);
      evttree->SetBranchAddress(name+"_detId", &detId);
      evttree->SetBranchAddress(name+"_keywire", &keywire);
      evttree->SetBranchAddress(name+"_strip", &strip);
      evttree->SetBranchAddress(name+"_strip8", &strip8);
      evttree->SetBranchAddress(name+"_valid", &valid);
      evttree->SetBranchAddress(name+"_type", &type);
      if (IsMatched) {
        matchTp = new std::vector<int>;
        evttree->SetBranchAddress(name+"_matchTp", &matchTp);
      }
    }

    // ALCT
    else if (data_type == 1) {
      detId = new std::vector<int>;
      keywire = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      evttree->SetBranchAddress(name+"_detId", &detId);
      evttree->SetBranchAddress(name+"_keywire", &keywire);
      evttree->SetBranchAddress(name+"_hit", &hit);
      evttree->SetBranchAddress(name+"_position", &position);
      evttree->SetBranchAddress(name+"_valid", &valid);
    }

    //CLCT
    else if (data_type == 2) {
      detId = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
      evttree->SetBranchAddress(name+"_detId", &detId);
      evttree->SetBranchAddress(name+"_strip", &strip);
      evttree->SetBranchAddress(name+"_strip8", &strip8);
      evttree->SetBranchAddress(name+"_hit", &hit);
      evttree->SetBranchAddress(name+"_position", &position);
      evttree->SetBranchAddress(name+"_valid", &valid);
      evttree->SetBranchAddress(name+"_bend", &bend);
      evttree->SetBranchAddress(name+"_pattern", &pattern);
      evttree->SetBranchAddress(name+"_slope", &slope);
    }

    // GEM
    else if (data_type == 3) {
      InitGP(evttree);
      detId = new std::vector<int>;
      strip = new std::vector<int>;
      evttree->SetBranchAddress(name+"_detId", &detId);
      evttree->SetBranchAddress(name+"_strip", &strip);
      if (IsMatched) {
        matchTp = new std::vector<int>;
        evttree->SetBranchAddress(name+"_matchTp", &matchTp);
      }
    }

    // GEMPad
    else if (data_type == 4) {
      InitGP(evttree);
      detId = new std::vector<int>;
      pad = new std::vector<int>;
      part = new std::vector<int>;
      evttree->SetBranchAddress(name+"_detId", &detId);
      evttree->SetBranchAddress(name+"_pad", &pad);
      evttree->SetBranchAddress(name+"_part", &part);
      if (IsMatched) {
        matchCSC = new std::vector<int>;
        evttree->SetBranchAddress(name+"_matchCSC", &matchCSC);
      }
    }

    //SimHit
    else if (data_type == 5) {
      if (!IsMatched) cout << name << " not matched? SimHits are always matched!" << endl;
      InitGP(evttree);
      station = new std::vector<int>;
      matchTp = new std::vector<int>;
      evttree->SetBranchAddress(name+"_station", &station);
      evttree->SetBranchAddress(name+"_matchTp", &matchTp);
    }

    // TP
    else if (data_type == 6) {
      if (IsMatched) cout << name << " Matched? TP can't match to anything" << endl;
      pt           = new std::vector<float>;
      eta          = new std::vector<float>;
      phi          = new std::vector<float>;
      dxy          = new std::vector<float>;
      d0           = new std::vector<float>;
      z0           = new std::vector<float>;
      d0_prod      = new std::vector<float>;
      z0_prod      = new std::vector<float>;
      pdgid        = new std::vector<int>;
      nmatch       = new std::vector<int>;
      nloosematch  = new std::vector<int>;
      nstub        = new std::vector<int>;
      eventid      = new std::vector<int>;
      charge       = new std::vector<int>;
      evttree->SetBranchAddress(name+"_pt",          &pt);
      evttree->SetBranchAddress(name+"_eta",         &eta);
      evttree->SetBranchAddress(name+"_phi",         &phi);
      evttree->SetBranchAddress(name+"_dxy",         &dxy);
      evttree->SetBranchAddress(name+"_d0",          &d0);
      evttree->SetBranchAddress(name+"_z0",          &z0);
      evttree->SetBranchAddress(name+"_d0_prod",     &d0_prod);
      evttree->SetBranchAddress(name+"_z0_prod",     &z0_prod);
      evttree->SetBranchAddress(name+"_pdgid",       &pdgid);
      evttree->SetBranchAddress(name+"_nmatch",      &nmatch);
      evttree->SetBranchAddress(name+"_nloosematch", &nloosematch);
      evttree->SetBranchAddress(name+"_nstub",       &nstub);
      evttree->SetBranchAddress(name+"_eventid",     &eventid);
      evttree->SetBranchAddress(name+"_charge",      &charge);
    }
  }

  void Init(TChain* evttree, TString name_, TString data_type_st, bool match_ = false) {
    int data_type_ = -1;
    if      (data_type_st == "LCT") data_type_ = 0;
    else if (data_type_st == "ALCT") data_type_ = 1;
    else if (data_type_st == "CLCT") data_type_ = 2;
    else if (data_type_st == "GEM") data_type_ = 3;
    else if (data_type_st == "GEMPad") data_type_ = 4;
    else if (data_type_st == "SimHit") data_type_ = 5;
    else if (data_type_st == "TP") data_type_ = 6;
    // else if (data_type_st == "RegionalMuon") data_type_ = 7;
    // else if (data_type_st == "MatchMuon") data_type_ = 8;


    if (data_type_ == -1) cout << "Wrong Data Type input for " << name_ << " as " << data_type_st << endl;
    else Init(evttree, name_, data_type_, match_);
  }

  void Reset() {
    if (data_type == 0) {
      ResetGP();
      bend->clear();
      pattern->clear();
      slope->clear();
      quality->clear();
      detId->clear();
      keywire->clear();
      strip->clear();
      strip8->clear();
      valid->clear();
      type->clear();
      if (IsMatched) {
        matchTp->clear();
      }
    }
    else if (data_type == 1) {
      detId->clear();
      keywire->clear();
      hit->clear();
      position->clear();
      valid->clear();
    }
    else if (data_type == 2) {
      detId->clear();
      strip->clear();
      strip8->clear();
      hit->clear();
      position->clear();
      valid->clear();
      bend->clear();
      pattern->clear();
      slope->clear();
    }
    else if (data_type == 3) {
      ResetGP();
      detId->clear();
      strip->clear();
      if (IsMatched) {
        matchTp->clear();
      }
    }
    else if (data_type == 4) {
      ResetGP();
      detId->clear();
      pad->clear();
      part->clear();
      if (IsMatched) {
        matchCSC->clear();
      }
    }
    else if (data_type == 5) {
      ResetGP();
      station->clear();
      matchTp->clear();
    }
    else if (data_type == 6) {
      pt->clear();
      eta->clear();
      phi->clear();
      dxy->clear();
      d0->clear();
      z0->clear();
      d0_prod->clear();
      z0_prod->clear();
      pdgid->clear();
      nmatch->clear();
      nloosematch->clear();
      nstub->clear();
      eventid->clear();
      charge->clear();
    }
  }

  void ResetGP() {
    eta->clear();
    phi->clear();
    z->clear();
    r->clear();
  }

  void InitGP(TChain* evttree) {
    phi = new std::vector<float>;
    eta = new std::vector<float>;
    z = new std::vector<float>;
    r = new std::vector<float>;
    evttree->SetBranchAddress(name+"_phi", &phi);
    evttree->SetBranchAddress(name+"_eta", &eta);
    evttree->SetBranchAddress(name+"_z", &z);
    evttree->SetBranchAddress(name+"_r", &r);
  }

  std::vector<int>*   detId;
  std::vector<int>*   matchTp;
  std::vector<int>*   matchCSC;

  //Matrix extraction
  std::vector<int>*   hit;
  std::vector<int>*   position;

  // GP
  std::vector<float>* phi;
  std::vector<float>* eta;
  std::vector<float>* z;
  std::vector<float>* r;

  //Digi
  std::vector<int>*   bend;
  std::vector<int>*   pattern;
  std::vector<int>*   slope;
  std::vector<int>*   quality;
  std::vector<int>*   keywire;
  std::vector<int>*   strip;
  std::vector<int>*   strip8;
  std::vector<bool>*  valid;
  std::vector<int>*   type;
  std::vector<int>*   part;
  std::vector<int>*   pad;

  //SimHit
  std::vector<int>* station;

  //TP
  std::vector<float>* pt;
  std::vector<float>* dxy;
  std::vector<float>* d0;
  std::vector<float>* z0;
  std::vector<float>* d0_prod;
  std::vector<float>* z0_prod;
  std::vector<int>* pdgid;
  std::vector<int>* nmatch;
  std::vector<int>* nloosematch;
  std::vector<int>* nstub;
  std::vector<int>* eventid;
  std::vector<int>* charge;

private:
  TString name;
  int data_type;
  bool IsMatched;
};

#endif
