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
  int      Index; //contains the oriignal unculled tp collection index to facilitate linking
  vector<TString> DetHit;
};

struct SimHit{
  float    phi;
  float    eta;
  float    r;
  float    z;
  int MatchTp;
};

struct GEMPadDigi{
  int pad;
  int part;
  float phi;
  float eta;
  float r;
  float z;
};

struct GEMPadDigiCluster{
  float phi;
  float eta;
  float r;
  float z;
  int MatchTp;
  vector<int> pads;
};

struct CSCStub{
  float phi;
  float eta;
  float r;
  float z;
  int   bend;
  int   pattern;
  int   slope;
  int   quality;
  int   detId;
  int   keywire;
  int   strip;
  int   strip8;
  bool  valid;
  int   type;
  int   GEM1pad;
  int   GEM1part;
  int   GEM2pad;
  int   GEM2part;
  int   MatchTp;
  vector<int> CLCT_hits;
  vector<int> CLCT_positions;
  vector<int> ALCT_hits;
  int layer; //runs from 0 innermost to x outermost, with x=2 for ME2-4 and x=4 for ME1, mind the single layer of ME1/3 being x=3
  GEMPadDigi GEM1;
  GEMPadDigi GEM2;
};

struct GEMDigi{
  float phi;
  float eta;
  float r;
  float z;
  int   MatchTp;
  int   layer; //0 is innermost, 3 is outermost, given the subdetector in question
};

struct SimHitAve{
  float phi;
  float eta;
  float r;
  float z;
};

// Container for Tracking Particle related information
class TPContent{
public:
  // idx = -1 by default
  TPContent(int idx) {
    SortedIndex = RawIndex = idx;
    CSCSimHits.clear();
    GEMSimHits.clear();
    MatchCSCStubs.clear();
    MatchGEMDigis.clear();
  };
  void CalcSimHitAve() {
    float ncsc = NSimHitsCSC = CSCSimHits.size();
    float ngem = NSimHitsGEM = GEMSimHits.size();
    CSCSimHitAve = GEMSimHitAve = SimHitAve{0.,0.,0.,0.};
    for (unsigned i = 0; i < NSimHitsCSC; ++i) {
      CSCSimHitAve.phi += CSCSimHits[i]->phi / ncsc;
      CSCSimHitAve.eta += CSCSimHits[i]->eta / ncsc;
      CSCSimHitAve.r   += CSCSimHits[i]->r / ncsc;
      CSCSimHitAve.z   += CSCSimHits[i]->z / ncsc;
    };
    CSCSimHitAve.phi = TVector2::Phi_mpi_pi(CSCSimHitAve.phi);
    for (unsigned i = 0; i < NSimHitsGEM; ++i) {
      GEMSimHitAve.phi += GEMSimHits[i]->phi / ngem;
      GEMSimHitAve.eta += GEMSimHits[i]->eta / ngem;
      GEMSimHitAve.r   += GEMSimHits[i]->r / ngem;
      GEMSimHitAve.z   += GEMSimHits[i]->z / ngem;
    };
    GEMSimHitAve.phi = TVector2::Phi_mpi_pi(GEMSimHitAve.phi);
  }
  // void SetTP(tp tp_) {
  //   TP = tp_;
  // }

  tp*               TP;
  int               RawIndex;
  int               SortedIndex;
  int               NSimHitsCSC;
  int               NSimHitsGEM;
  SimHitAve         CSCSimHitAve;
  SimHitAve         GEMSimHitAve;
  vector<SimHit*>   CSCSimHits;
  vector<SimHit*>   GEMSimHits;
  vector<CSCStub*>  MatchCSCStubs;
  vector<GEMDigi*>  MatchGEMDigis;
  vector<GEMPadDigiCluster*> MatchGEMPadDigiClusters;
};

// Container for station information
class StationData{
public:
  StationData(TString lable) {
    StationLabel = lable;
    TPInfos.clear();
    CSCSimHits.clear();
    GEMSimHits.clear();
    MatchCSCStubs.clear();
    AllCSCStubs.clear();
    MatchGEMDigis.clear();
    AllGEMDigis.clear();
    AllGEMPadDigis.clear();
    MatchGEMPadDigiClusters.clear();
    AllGEMPadDigiClusters.clear();
  };
  // void AddCSCSimHit(SimHit& d) {
  //   for (unsigned i = 0; i < TPInfos.size(); ++i) {
  //     if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].CSCSimHits.push_back(d);
  //     return;
  //   }
  //   TPContent tmp((unsigned)d.MatchTp);
  //   tmp.CSCSimHits.push_back(d);
  //   TPInfos.push_back(tmp);
  // }
  // void AddGEMSimHit(SimHit& d) {
  //   for (unsigned i = 0; i < TPInfos.size(); ++i) {
  //     if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].GEMSimHits.push_back(d);
  //     return;
  //   }
  //   TPContent tmp((unsigned)d.MatchTp);
  //   tmp.GEMSimHits.push_back(d);
  //   TPInfos.push_back(tmp);
  // }
  // void AddCSCStub(CSCStub& d) {
  //   if (d.MatchTp == -1) {
  //     CSCStubs.push_back(d);
  //     return;
  //   }
  //   for (unsigned i = 0; i < TPInfos.size(); ++i) {
  //     if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].CSCStubs.push_back(d);
  //     return;
  //   }
  //   TPContent tmp((unsigned)d.MatchTp);
  //   tmp.CSCStubs.push_back(d);
  //   TPInfos.push_back(tmp);
  // }
  // void AddGEMDigi(GEMDigi& d) {
  //   if (d.MatchTp == -1) {
  //     GEMDigis.push_back(d);
  //     return;
  //   }
  //   for (unsigned i = 0; i < TPInfos.size(); ++i) {
  //     if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].GEMDigis.push_back(d);
  //     return;
  //   }
  //   TPContent tmp((unsigned)d.MatchTp);
  //   tmp.GEMDigis.push_back(d);
  //   TPInfos.push_back(tmp);
  // }
  // void AddGEMPadDigi(GEMPadDigi d) {
  //   GEMPadDigis.push_back(d);
  // }
  // void AddGEMPadDigiCluster(GEMPadDigiCluster& d) {
  //   if (d.MatchTp == -1) {
  //     Clusters.push_back(d);
  //     return;
  //   }
  //   for (unsigned i = 0; i < TPInfos.size(); ++i) {
  //     if (TPInfos[i].RawIndex == (unsigned)d.MatchTp) TPInfos[i].Clusters.push_back(d);
  //     return;
  //   }
  //   TPContent tmp((unsigned)d.MatchTp);
  //   tmp.Clusters.push_back(d);
  //   TPInfos.push_back(tmp);
  // }

  int FindTP(int ind) {
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      if (ind == TPInfos[i].RawIndex) return i;
    }
    return -1;
  }

  void SortTP() {
    for (unsigned i = 0; i < CSCSimHits.size(); ++i) {
      int tpindex = FindTP(CSCSimHits[i]->MatchTp);
      if (tpindex == -1) {
        TPContent tmp(CSCSimHits[i]->MatchTp);
        tmp.CSCSimHits.push_back(CSCSimHits[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].CSCSimHits.push_back(CSCSimHits[i]);
    }
    for (unsigned i = 0; i < GEMSimHits.size(); ++i) {
      int tpindex = FindTP(GEMSimHits[i]->MatchTp);
      if (tpindex == -1) {
        TPContent tmp(GEMSimHits[i]->MatchTp);
        tmp.GEMSimHits.push_back(GEMSimHits[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].GEMSimHits.push_back(GEMSimHits[i]);
    }
    for (unsigned i = 0; i < MatchCSCStubs.size(); ++i) {
      int tpindex = FindTP(MatchCSCStubs[i]->MatchTp);
      if (tpindex == -1) {
        TPContent tmp(MatchCSCStubs[i]->MatchTp);
        tmp.MatchCSCStubs.push_back(MatchCSCStubs[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].MatchCSCStubs.push_back(MatchCSCStubs[i]);
    }
    for (unsigned i = 0; i < MatchGEMDigis.size(); ++i) {
      int tpindex = FindTP(MatchGEMDigis[i]->MatchTp);
      if (tpindex == -1) {
        TPContent tmp(MatchGEMDigis[i]->MatchTp);
        tmp.MatchGEMDigis.push_back(MatchGEMDigis[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].MatchGEMDigis.push_back(MatchGEMDigis[i]);
    }
    for (unsigned i = 0; i < MatchGEMPadDigiClusters.size(); ++i) {
      int tpindex = FindTP(MatchGEMPadDigiClusters[i]->MatchTp);
      if (tpindex == -1) {
        TPContent tmp(MatchGEMPadDigiClusters[i]->MatchTp);
        tmp.MatchGEMPadDigiClusters.push_back(MatchGEMPadDigiClusters[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].MatchGEMPadDigiClusters.push_back(MatchGEMPadDigiClusters[i]);
    }
  }

  TString            StationLabel;
  vector<TPContent>  TPInfos;
  vector<SimHit*>     CSCSimHits;
  vector<SimHit*>     GEMSimHits;
  vector<CSCStub*>    MatchCSCStubs;
  vector<CSCStub*>      AllCSCStubs;
  vector<GEMDigi*>    MatchGEMDigis;
  vector<GEMDigi*>      AllGEMDigis;
  vector<GEMPadDigi*> AllGEMPadDigis;
  vector<GEMPadDigiCluster*> MatchGEMPadDigiClusters;
  vector<GEMPadDigiCluster*>   AllGEMPadDigiClusters;
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
        StationData("1-3"),
        StationData("1-0")
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
    CSCSimHits.clear();
    GEMSimHits.clear();
    MatchCSCStubs.clear();
    AllCSCStubs.clear();
    MatchGEMDigis.clear();
    AllGEMDigis.clear();
    AllGEMPadDigis.clear();
    MatchGEMPadDigiClusters.clear();
    AllGEMPadDigiClusters.clear();
    I_MuonTPs.clear();
    I_CSCSimHits.clear();
    I_GEMSimHits.clear();
    I_MatchCSCStubs.clear();
    I_AllCSCStubs.clear();
    I_MatchGEMDigis.clear();
    I_AllGEMDigis.clear();
    I_AllGEMPadDigis.clear();
    I_MatchGEMPadDigiClusters.clear();
    I_AllGEMPadDigiClusters.clear();
  }
  // StationData& Station(pair<unsigned,unsigned> DaR){
  //   return Stations[DaR.first][DaR.second];
  // }
  // void AddTP(tp& muontp) {
  //   MuonTPs.push_back(muontp);
  // }
  // void AddCSCSimHit(SimHit& d) {
  //   pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
  //   Stations[DaR.first][DaR.second].AddCSCSimHit(d);
  // }
  // void AddGEMSimHit(SimHit& d) {
  //   pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
  //   Stations[DaR.first][DaR.second].AddGEMSimHit(d);
  // }
  // void AddCSCStub(CSCStub& d) {
  //   pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
  //   Stations[DaR.first][DaR.second].AddCSCStub(d);
  // }
  // void AddGEMDigi(GEMDigi& d) {
  //   pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
  //   Stations[DaR.first][DaR.second].AddGEMDigi(d);
  // }
  // void AddGEMPadDigi(GEMPadDigi d) {
  //   pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
  //   Stations[DaR.first][DaR.second].AddGEMPadDigi(d);
  // }
  // void AddGEMPadDigiCluster(GEMPadDigiCluster& d) {
  //   pair<unsigned,unsigned> DaR = DiskAndRing(d.r,d.z);
  //   Stations[DaR.first][DaR.second].AddGEMPadDigiCluster(d);
  // }
  void SortByStation() {
    for (unsigned i = 0; i < I_MuonTPs.size(); ++i) {
      tp* p = &(I_MuonTPs[i]);
      MuonTPs.push_back(p);
    }
    for (unsigned i = 0; i < I_CSCSimHits.size(); ++i) {
      SimHit* p = &(I_CSCSimHits[i]);
      CSCSimHits.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].CSCSimHits.push_back(p);
    }
    for (unsigned i = 0; i < I_GEMSimHits.size(); ++i) {
      SimHit* p = &(I_GEMSimHits[i]);
      GEMSimHits.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].GEMSimHits.push_back(p);
    }
    for (unsigned i = 0; i < I_MatchCSCStubs.size(); ++i) {
      CSCStub* p = &(I_MatchCSCStubs[i]);
      MatchCSCStubs.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].MatchCSCStubs.push_back(p);
    }
    for (unsigned i = 0; i < I_AllCSCStubs.size(); ++i) {
      CSCStub* p = &(I_AllCSCStubs[i]);
      AllCSCStubs.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].AllCSCStubs.push_back(p);
    }
    for (unsigned i = 0; i < I_MatchGEMDigis.size(); ++i) {
      GEMDigi* p = &(I_MatchGEMDigis[i]);
      MatchGEMDigis.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].MatchGEMDigis.push_back(p);
    }
    for (unsigned i = 0; i < I_AllGEMDigis.size(); ++i) {
      GEMDigi* p = &(I_AllGEMDigis[i]);
      AllGEMDigis.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].AllGEMDigis.push_back(p);
    }
    for (unsigned i = 0; i < I_AllGEMPadDigis.size(); ++i) {
      GEMPadDigi* p = &(I_AllGEMPadDigis[i]);
      AllGEMPadDigis.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].AllGEMPadDigis.push_back(p);
    }
    for (unsigned i = 0; i < I_MatchGEMPadDigiClusters.size(); ++i) {
      GEMPadDigiCluster* p = &(I_MatchGEMPadDigiClusters[i]);
      MatchGEMPadDigiClusters.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].MatchGEMPadDigiClusters.push_back(p);
    }
    for (unsigned i = 0; i < I_AllGEMPadDigiClusters.size(); ++i) {
      GEMPadDigiCluster* p = &(I_AllGEMPadDigiClusters[i]);
      AllGEMPadDigiClusters.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z);
      Stations[DaR.first][DaR.second].AllGEMPadDigiClusters.push_back(p);
    }
  }
  void SortTP() {
    for (unsigned i = 0; i < Stations.size(); ++i) for(unsigned j = 0; j < Stations[i].size(); ++j) {
      Stations[i][j].SortTP();
      for (unsigned k = 0; k < Stations[i][j].TPInfos.size(); ++k) {
        bool found = false;
        for (unsigned l = 0; l < MuonTPs.size(); ++l) {
          if (Stations[i][j].TPInfos[k].RawIndex == I_MuonTPs[l].Index) {
            I_MuonTPs[l].DetHit.push_back(Stations[i][j].StationLabel);
            Stations[i][j].TPInfos[k].TP = MuonTPs[l];
            found = true;
            break;
          }
        }
        if (!found) cout << "Cannot find a tp with index = " << Stations[i][j].TPInfos[k].RawIndex <<endl;
      }
    }
  }
  void CalcSimHitAve() {
    for (unsigned i = 0; i < Stations.size(); ++i) for(unsigned j = 0; j < Stations[i].size(); ++j) {
      for (unsigned k = 0; k < Stations[i][j].TPInfos.size(); ++k) {
        Stations[i][j].TPInfos[k].CalcSimHitAve();
      }
    }
  }
  void Run() {
    SortByStation();
    SortTP();
    CalcSimHitAve();
  }

  vector<vector<StationData> > Stations;
  vector<tp> I_MuonTPs;
  vector<SimHit> I_CSCSimHits;
  vector<SimHit> I_GEMSimHits;
  vector<CSCStub> I_MatchCSCStubs;
  vector<CSCStub>   I_AllCSCStubs;
  vector<GEMDigi> I_MatchGEMDigis;
  vector<GEMDigi>   I_AllGEMDigis;
  vector<GEMPadDigi>   I_AllGEMPadDigis;
  vector<GEMPadDigiCluster> I_MatchGEMPadDigiClusters;
  vector<GEMPadDigiCluster>   I_AllGEMPadDigiClusters;
  // Make pointers to the instances, so that they are in consistent format with Station and TPContent
  vector<tp*> MuonTPs;
  vector<SimHit*> CSCSimHits;
  vector<SimHit*> GEMSimHits;
  vector<CSCStub*> MatchCSCStubs;
  vector<CSCStub*>   AllCSCStubs;
  vector<GEMDigi*> MatchGEMDigis;
  vector<GEMDigi*>   AllGEMDigis;
  vector<GEMPadDigi*>   AllGEMPadDigis;
  vector<GEMPadDigiCluster*> MatchGEMPadDigiClusters;
  vector<GEMPadDigiCluster*>   AllGEMPadDigiClusters;

  // vector<GEMPadDigi> MatchGEMPadDigis; // This collection is integrated into LCT collections
};



#endif
