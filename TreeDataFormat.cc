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

struct DetId {
  int detId;
  int zendcap;
  int ring;
  int station;
  int layer;
  int chamber;
  int roll; // Exclusive to GEM
  bool SameXY(DetId id) {
    return (zendcap == id.zendcap && (ring == id.ring || abs(ring - id.ring) == 3) && station == id.station && (abs(chamber - id.chamber) < 2 || abs(chamber - id.chamber) == 17 || abs(chamber - id.chamber) == 35) && abs(roll - id.roll) < 2);
  }
  bool operator==(const DetId id) {
    return (detId == id.detId && zendcap == id.zendcap && ring == id.ring && station == id.station && layer == id.layer && chamber == id.chamber && roll == id.roll);
  }
};

pair<int, int> DiskAndRing(double r, double z, DetId det, TString info = "") {
  pair<int, int> DaR = DiskAndRing(r,z);
  if (DaR.first == (det.station - 1 ) && DaR.second == (det.ring - 1)) return DaR; // Consistent case
  else if (det.station == 0 && det.ring == 1 && DaR.first == 0 && DaR.second == 3) return DaR; // ME0 is noted as (0,3) in our code, which is (0,1) in CMSSW
  else if (det.station == 1 && det.ring == 4 && DaR.first == 0 && DaR.second == 0) return DaR; // ME1/1 is divided into ME1a and ME1b, where ME1a is noted as (1,4) in CMSSW
  else {
    cout << info << " ";
    cout << Form("For r = %f, z = %f, DaR gives (%i, %i), while det gives (%i, %i)", r, z, DaR.first + 1, DaR.second + 1, det.station, det.ring)<<endl;
    return DaR;
  }
}

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
  DetId    det;
  float    phi;
  float    eta;
  float    r;
  float    z;
  int MatchIndex;
};

struct GEMPadDigi{
  DetId det;
  int pad;
  int strip;
  int strip8;
  int strip_me1a;
  int strip8_me1a;
  int keywire_min;
  int keywire_max;
  int part;
  float phi;
  float eta;
  float r;
  float z;
  int MatchIndex;
};

struct GEMPadDigiCluster{
  DetId det;
  int pad;
  int strip;
  int strip8;
  int strip_me1a;
  int strip8_me1a;
  int keywire_min;
  int keywire_max;
  int part;
  float phi;
  float eta;
  float r;
  float z;
  int MatchIndex;
  vector<int> pads;
};

struct CSCStub{
  DetId det;
  float phi;
  float eta;
  float r;
  float z;
  int   bend;
  int   pattern;
  int   slope;
  int   quality;
  int   keywire;
  int   strip;
  int   strip8;
  bool  valid;
  int   type;
  int   GEM1pad;
  int   GEM1strip;
  int   GEM1strip8;
  int   GEM1strip_me1a;
  int   GEM1strip8_me1a;
  int   GEM1keywire_min;
  int   GEM1keywire_max;
  int   GEM1roll;
  int   GEM1part;
  int   GEM2pad;
  int   GEM2strip;
  int   GEM2strip8;
  int   GEM2strip_me1a;
  int   GEM2strip8_me1a;
  int   GEM2keywire_min;
  int   GEM2keywire_max;
  int   GEM2roll;
  int   GEM2part;
  int   MatchIndex;
  vector<int> CLCT_hits;
  vector<int> CLCT_positions;
  vector<int> ALCT_hits;
  int layer; //runs from 0 innermost to x outermost, with x=2 for ME2-4 and x=4 for ME1, mind the single layer of ME1/3 being x=3
  GEMPadDigi GEM1;
  GEMPadDigi GEM2;
};

struct GEMDigi{
  DetId det;
  float phi;
  float eta;
  float r;
  float z;
  int   MatchIndex;
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
    MatchGEMPadDigis.clear();
    MatchGEMPadDigiClusters.clear();
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
    // SimHit DetId Verification
    for (unsigned i = 0; i < NSimHitsCSC; ++i) {
      DetId det1 = CSCSimHits[i]->det;
      if (i == 0) CSCDetId = det1;
      else if (!(det1.SameXY(CSCDetId)) ) cout << Form("Inconsistent DetId in CSC: S:%i %i, Z:%d %d, R:%d %d, L:%d %d, C:%d %d",CSCDetId.station, det1.station, CSCDetId.zendcap, det1.zendcap, CSCDetId.ring,det1.ring,CSCDetId.layer,det1.layer,CSCDetId.chamber, det1.chamber) <<endl;
    }
    for (unsigned i = 0; i < NSimHitsGEM; ++i) {
      DetId det1 = GEMSimHits[i]->det;
      if (i == 0) GEMDetId = det1;
      else if (!(det1.SameXY(GEMDetId))) cout << Form("Inconsistent DetId in GEM: S:%i %i, Z:%d %d, R:%d %d, L:%d %d, C:%d %d, Ro:%d %d",GEMDetId.station, det1.station, GEMDetId.zendcap, det1.zendcap, GEMDetId.ring,det1.ring,GEMDetId.layer,det1.layer,GEMDetId.chamber, det1.chamber,GEMDetId.roll, det1.roll) <<endl;
    }
  }

  tp*               TP;
  DetId             GEMDetId;
  DetId             CSCDetId;
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
  vector<GEMPadDigi*>   MatchGEMPadDigis;
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
    MatchGEMPadDigis.clear();
    AllGEMPadDigis.clear();
    MatchGEMPadDigiClusters.clear();
    AllGEMPadDigiClusters.clear();
  };

  int FindTP(int ind) {
    for (unsigned i = 0; i < TPInfos.size(); ++i) {
      if (ind == TPInfos[i].RawIndex) return i;
    }
    return -1;
  }

  void SortTP() {
    for (unsigned i = 0; i < CSCSimHits.size(); ++i) {
      int tpindex = FindTP(CSCSimHits[i]->MatchIndex);
      if (tpindex == -1) {
        TPContent tmp(CSCSimHits[i]->MatchIndex);
        tmp.CSCSimHits.push_back(CSCSimHits[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].CSCSimHits.push_back(CSCSimHits[i]);
    }
    for (unsigned i = 0; i < GEMSimHits.size(); ++i) {
      int tpindex = FindTP(GEMSimHits[i]->MatchIndex);
      if (tpindex == -1) {
        TPContent tmp(GEMSimHits[i]->MatchIndex);
        tmp.GEMSimHits.push_back(GEMSimHits[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].GEMSimHits.push_back(GEMSimHits[i]);
    }
    for (unsigned i = 0; i < MatchCSCStubs.size(); ++i) {
      int tpindex = FindTP(MatchCSCStubs[i]->MatchIndex);
      if (tpindex == -1) {
        TPContent tmp(MatchCSCStubs[i]->MatchIndex);
        tmp.MatchCSCStubs.push_back(MatchCSCStubs[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].MatchCSCStubs.push_back(MatchCSCStubs[i]);
    }
    for (unsigned i = 0; i < MatchGEMDigis.size(); ++i) {
      int tpindex = FindTP(MatchGEMDigis[i]->MatchIndex);
      if (tpindex == -1) {
        TPContent tmp(MatchGEMDigis[i]->MatchIndex);
        tmp.MatchGEMDigis.push_back(MatchGEMDigis[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].MatchGEMDigis.push_back(MatchGEMDigis[i]);
    }
    for (unsigned i = 0; i < MatchGEMPadDigis.size(); ++i) {
      int tpindex = FindTP(MatchGEMPadDigis[i]->MatchIndex);
      if (tpindex == -1) {
        TPContent tmp(MatchGEMPadDigis[i]->MatchIndex);
        tmp.MatchGEMPadDigis.push_back(MatchGEMPadDigis[i]);
        TPInfos.push_back(tmp);
      }
      else TPInfos[tpindex].MatchGEMPadDigis.push_back(MatchGEMPadDigis[i]);
    }
    for (unsigned i = 0; i < MatchGEMPadDigiClusters.size(); ++i) {
      int tpindex = FindTP(MatchGEMPadDigiClusters[i]->MatchIndex);
      if (tpindex == -1) {
        TPContent tmp(MatchGEMPadDigiClusters[i]->MatchIndex);
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
  vector<GEMPadDigi*> MatchGEMPadDigis;
  vector<GEMPadDigi*>   AllGEMPadDigis;
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
        StationData("ME0")
      },
      {
        StationData("1-1"),// Both ME1a and ME1b are filled here. And they have their own copy in the last two StationDatas.
        StationData("1-2"),
        StationData("1-3"),
        StationData("ME1a"),
        StationData("ME1b")
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
    MatchGEMPadDigis.clear();
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
    I_MatchGEMPadDigis.clear();
    I_AllGEMPadDigis.clear();
    I_MatchGEMPadDigiClusters.clear();
    I_AllGEMPadDigiClusters.clear();
  }

  int ME1abCon(DetId det) {
    if (det.station == 1 && det.ring == 4) return 1;
    if (det.station == 1 && det.ring == 1) return 5;
    return 0;
  }

  void SortByStation() {
    for (unsigned i = 0; i < I_MuonTPs.size(); ++i) {
      tp* p = &(I_MuonTPs[i]);
      MuonTPs.push_back(p);
    }
    for (unsigned i = 0; i < I_CSCSimHits.size(); ++i) {
      SimHit* p = &(I_CSCSimHits[i]);
      CSCSimHits.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"CSCSimHits");
      Stations[p->det.station][p->det.ring - 1].CSCSimHits.push_back(p);
      int ME1abRing = ME1abCon(p->det);
      if (ME1abRing) Stations[p->det.station][ME1abRing - 1].CSCSimHits.push_back(p);
    }
    for (unsigned i = 0; i < I_GEMSimHits.size(); ++i) {
      SimHit* p = &(I_GEMSimHits[i]);
      GEMSimHits.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"GEMSimHits");
      Stations[p->det.station][p->det.ring - 1].GEMSimHits.push_back(p);
    }
    for (unsigned i = 0; i < I_MatchCSCStubs.size(); ++i) {
      CSCStub* p = &(I_MatchCSCStubs[i]);
      MatchCSCStubs.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"MatchCSCStubs");
      Stations[p->det.station][p->det.ring - 1].MatchCSCStubs.push_back(p);
      int ME1abRing = ME1abCon(p->det);
      if (ME1abRing) Stations[p->det.station][ME1abRing - 1].MatchCSCStubs.push_back(p);
    }
    for (unsigned i = 0; i < I_AllCSCStubs.size(); ++i) {
      CSCStub* p = &(I_AllCSCStubs[i]);
      AllCSCStubs.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"AllCSCStubs");
      Stations[p->det.station][p->det.ring - 1].AllCSCStubs.push_back(p);
      int ME1abRing = ME1abCon(p->det);
      if (ME1abRing) Stations[p->det.station][ME1abRing - 1].AllCSCStubs.push_back(p);
    }
    for (unsigned i = 0; i < I_MatchGEMDigis.size(); ++i) {
      GEMDigi* p = &(I_MatchGEMDigis[i]);
      MatchGEMDigis.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"MatchGEMDigis");
      Stations[p->det.station][p->det.ring - 1].MatchGEMDigis.push_back(p);
    }
    for (unsigned i = 0; i < I_AllGEMDigis.size(); ++i) {
      GEMDigi* p = &(I_AllGEMDigis[i]);
      AllGEMDigis.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"AllGEMDigis");
      Stations[p->det.station][p->det.ring - 1].AllGEMDigis.push_back(p);
    }
    for (unsigned i = 0; i < I_MatchGEMPadDigis.size(); ++i) {
      GEMPadDigi* p = &(I_MatchGEMPadDigis[i]);
      MatchGEMPadDigis.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"MatchGEMPadDigis");
      Stations[p->det.station][p->det.ring - 1].MatchGEMPadDigis.push_back(p);
    }
    for (unsigned i = 0; i < I_AllGEMPadDigis.size(); ++i) {
      GEMPadDigi* p = &(I_AllGEMPadDigis[i]);
      AllGEMPadDigis.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"AllGEMPadDigis");
      Stations[p->det.station][p->det.ring - 1].AllGEMPadDigis.push_back(p);
    }
    for (unsigned i = 0; i < I_MatchGEMPadDigiClusters.size(); ++i) {
      GEMPadDigiCluster* p = &(I_MatchGEMPadDigiClusters[i]);
      MatchGEMPadDigiClusters.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"MatchGEMPadDigiClusters");
      Stations[p->det.station][p->det.ring - 1].MatchGEMPadDigiClusters.push_back(p);
    }
    for (unsigned i = 0; i < I_AllGEMPadDigiClusters.size(); ++i) {
      GEMPadDigiCluster* p = &(I_AllGEMPadDigiClusters[i]);
      AllGEMPadDigiClusters.push_back(p);
      pair<unsigned,unsigned> DaR = DiskAndRing(p->r,p->z,p->det,"AllGEMPadDigiClusters");
      Stations[p->det.station][p->det.ring - 1].AllGEMPadDigiClusters.push_back(p);
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
  vector<GEMPadDigi>   I_MatchGEMPadDigis;
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
  vector<GEMPadDigi*> MatchGEMPadDigis;
  vector<GEMPadDigi*>   AllGEMPadDigis;
  vector<GEMPadDigiCluster*> MatchGEMPadDigiClusters;
  vector<GEMPadDigiCluster*>   AllGEMPadDigiClusters;
};



#endif
