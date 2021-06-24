#include "TreeDataFormat.cc"

#include "TROOT.h"
#include "TChain.h"
#include "TVector2.h"

#include <utility>
#include <cmath>
#include <iostream>

using namespace std;

class TreeReader{
public:
  TreeReader(TChain* t_) {
    eventTree = t_;
  }
  void Init() {
    br_tp               .Init(eventTree,"tp","TP",false);
    br_cscSimHit        .Init(eventTree,"cscSimHit","SimHit",true);
    br_gemSimHit        .Init(eventTree,"gemSimHit","SimHit",true);
    br_allCscStubsLCT   .Init(eventTree,"allCscStubsLCT","LCT",false);
    br_allCscStubsALCT  .Init(eventTree,"allCscStubsALCT","ALCT",false);
    br_allCscStubsCLCT  .Init(eventTree,"allCscStubsCLCT","CLCT",false);
    br_allALCT          .Init(eventTree,"allALCT","ALCT",false);
    br_allCLCT          .Init(eventTree,"allCLCT","CLCT",false);
    br_allGemDigi       .Init(eventTree,"allGemDigi","GEM",false);
    br_matchCscStubsLCT .Init(eventTree,"matchCscStubsLCT","LCT",true);
    br_matchCscStubsALCT.Init(eventTree,"matchCscStubsALCT","ALCT",false);
    br_matchCscStubsCLCT.Init(eventTree,"matchCscStubsCLCT","CLCT",false);
    br_matchGemDigi     .Init(eventTree,"matchGemDigi","GEM",true);
    br_matchCscGEM1     .Init(eventTree,"matchCscGEM1","GEMPad",true);
    br_matchCscGEM2     .Init(eventTree,"matchCscGEM2","GEMPad",true);
    br_allCscGEM1       .Init(eventTree,"allCscGEM1","GEMPad",true);
    br_allCscGEM2       .Init(eventTree,"allCscGEM2","GEMPad",true);
    br_gemPadDigi       .Init(eventTree,"gemPadDigi","GEMPad",false);
  }

  void ReadTree() {
    Reset();

    //TP
    const unsigned tp_size = br_tp.eta->size();
    for (unsigned i = 0; i < size; ++i) {
      if(abs(br_tp.pdgid->at(i))!=13) continue;
      tp tmp;
      tmp.pt      = br_tp.pt->at(i);
      tmp.eta     = br_tp.eta->at(i);
      tmp.phi     = br_tp.phi->at(i);
      tmp.dxy     = br_tp.dxy->at(i);
      tmp.d0      = br_tp.d0->at(i);
      tmp.z0      = br_tp.z0->at(i);
      tmp.d0_prod = br_tp.d0_prod->at(i);
      tmp.z0_prod = br_tp.z0_prod->at(i);
      tmp.pdgid   = br_tp.pdgid->at(i);
      tmp.eventid = br_tp.eventid->at(i);
      tmp.charge  = br_tp.charge->at(i);
      tmp.index   = i;
      MuonTPs.push_back(tmp);
    };

    //GEM SimHit
    const unsigned GEMSimHitSize = br_gemSimHit.phi->size();
    for (unsigned i = 0; i < GEMSimHitSize; ++i) {
      if (abs(br_tp.pdgid->at(br_gemSimHit.matchTp->at(i))) != 13) continue;
      SimHit tmp;
      tmp.eta     = br_gemSimHit.eta->at(i);
      tmp.phi     = br_gemSimHit.phi->at(i);
      tmp.r       = br_gemSimHit.r->at(i);
      tmp.z       = br_gemSimHit.z->at(i);
      tmp.matchTp = (unsigned)MuonTPindex((unsigned)br_gemSimHit.matchTp->at(i), MuonTPs);
      GEMSimHits.push_back(tmp);
    };

    //CSC SimHit
    const unsigned CSCSimHitSize=br_cscSimHit.phi->size();
    for(unsigned i = 0; i < CSCSimHitSize; ++i) {
      if(abs(br_tp.pdgid->at(br_cscSimHit.matchTp->at(i))) != 13) continue;//remove SimHits not from muons
      SimHit tmp;
      tmp.phi     = br_cscSimHit.phi->at(i);
      tmp.eta     = br_cscSimHit.eta->at(i);
      tmp.r       = br_cscSimHit.r->at(i);
      tmp.z       = br_cscSimHit.z->at(i);
      tmp.matchTp = (unsigned)MuonTPindex((unsigned)br_cscSimHit.matchTp->at(i), MuonTPs);
      CSCSimHits.push_back(tmp);
    };

    //CSCStubs
    const unsigned allCscStubsLCTSize = br_allCscStubsLCT.phi->size();
    if (br_allCscGEM1.phi->size() != allCscStubsLCTSize) cout << "allCscStubsLCT size different from GEM1" <<endl;
    if (br_allCscGEM2.phi->size() != allCscStubsLCTSize) cout << "allCscStubsLCT size different from GEM2" <<endl;
    for(unsigned i = 0; i < allCscStubsLCTSize; ++i) {
      CSCStub tmp;
      tmp.phi     = br_allCscStubsLCT.phi->at(i);
      tmp.eta     = br_allCscStubsLCT.eta->at(i);
      tmp.r       = br_allCscStubsLCT.r->at(i);
      tmp.z       = br_allCscStubsLCT.z->at(i);
      tmp.bend    = br_allCscStubsLCT.bend->at(i);
      tmp.slope   = br_allCscStubsCLCT.slope->at(i);
      if (br_allCscStubsCLCT.slope->at(i) != br_allCscStubsLCT.slope->at(i)) cout << "Inconsistent allCscStubsLCT and allCscStubsCLCT slope" << endl;
      tmp.strip8  = br_allCscStubsLCT.strip8->at(i);
      tmp.pattern = br_allCscStubsLCT.pattern->at(i);
      tmp.valid   = br_allCscStubsLCT.valid->at(i);
      tmp.matchTp = -1;
      tmp.CLCT_hits.clear();
      tmp.CLCT_positions.clear();
      tmp.ALCT_hits.clear();
      for (unsigned j = 0; j < 6; ++j) {
        tmp.CLCT_hits.push_back(br_allCscStubsCLCT.hit->at(6*i+j));
        tmp.CLCT_positions.push_back(br_allCscStubsCLCT.position->at(6*i+j));
        tmp.ALCT_hits,push_back(br_allCscStubsALCT.hit->at(6*i+j));
      }
      tmp.layer = CSClayerIndex(br_allCscStubsLCT.z->at(i));
      if (br_allCscGEM1.matchCSC->at(i) != i) cout << "allCSCGEM1 index inconsistent with allCSCStub" << endl;
      if (br_allCscGEM2.matchCSC->at(i) != i) cout << "allCSCGEM2 index inconsistent with allCSCStub" << endl;
      tmp.GEM1.phi = br_allCscGEM1.phi->at(i);
      tmp.GEM1.eta = br_allCscGEM1.eta->at(i);
      tmp.GEM1.r   = br_allCscGEM1.r->at(i);
      tmp.GEM1.z   = br_allCscGEM1.z->at(i);
      tmp.GEM2.phi = br_allCscGEM2.phi->at(i);
      tmp.GEM2.eta = br_allCscGEM2.eta->at(i);
      tmp.GEM2.r   = br_allCscGEM2.r->at(i);
      tmp.GEM2.z   = br_allCscGEM2.z->at(i);
      CSCStubs.push_back(tmp);
    };

    // MatchedCSCStubs
    const unsigned matchCscStubsLCTSize = br_matchCscStubsLCT.phi->size();
    if (br_matchCscGEM1.phi->size() != matchCscStubsLCTSize) cout << "matchCscStubsLCT size different from GEM1" <<endl;
    if (br_matchCscGEM2.phi->size() != matchCscStubsLCTSize) cout << "matchCscStubsLCT size different from GEM2" <<endl;
    for(unsigned i = 0; i < matchCscStubsLCTSize; ++i) {
      if(abs(br_tp.pdgid->at(br_matchCscStubsLCT.matchTp->at(i)))!=13) continue;
      CSCStub tmp;
      tmp.phi     = br_matchCscStubsLCT.phi->at(i);
      tmp.eta     = br_matchCscStubsLCT.eta->at(i);
      tmp.r       = br_matchCscStubsLCT.r->at(i);
      tmp.z       = br_matchCscStubsLCT.z->at(i);
      tmp.bend    = br_matchCscStubsLCT.bend->at(i);
      tmp.slope   = br_matchCscStubsCLCT.slope->at(i);
      if (br_matchCscStubsCLCT.slope->at(i) != br_matchCscStubsLCT.slope->at(i)) cout << "Inconsistent matchCscStubsLCT and matchCscStubsCLCT slope" << endl;
      tmp.strip8  = br_matchCscStubsLCT.strip8->at(i);
      tmp.pattern = br_matchCscStubsLCT.pattern->at(i);
      tmp.valid   = br_matchCscStubsLCT.valid->at(i);
      tmp.matchTp = MuonTPindex(br_matchCscStubsLCT.matchTp->at(i), MuonTPs);
      tmp.CLCT_hits.clear();
      tmp.CLCT_positions.clear();
      tmp.ALCT_hits.clear();
      for (unsigned j = 0; j < 6; ++j) {
        tmp.CLCT_hits.push_back(br_matchCscStubsCLCT.hit->at(6*i+j));
        tmp.CLCT_positions.push_back(br_matchCscStubsCLCT.position->at(6*i+j));
        tmp.ALCT_hits,push_back(br_matchCscStubsALCT.hit->at(6*i+j));
      }
      tmp.layer = CSClayerIndex(br_matchCscStubsLCT.z->at(i));
      if (br_matchCscGEM1.matchCSC->at(i) != i) cout << "matchCSCGEM1 index inconsistent with matchCSCStub" << endl;
      if (br_matchCscGEM2.matchCSC->at(i) != i) cout << "matchCSCGEM2 index inconsistent with matchCSCStub" << endl;
      tmp.GEM1.phi = br_matchCscGEM1.phi->at(i);
      tmp.GEM1.eta = br_matchCscGEM1.eta->at(i);
      tmp.GEM1.r   = br_matchCscGEM1.r->at(i);
      tmp.GEM1.z   = br_matchCscGEM1.z->at(i);
      tmp.GEM2.phi = br_matchCscGEM2.phi->at(i);
      tmp.GEM2.eta = br_matchCscGEM2.eta->at(i);
      tmp.GEM2.r   = br_matchCscGEM2.r->at(i);
      tmp.GEM2.z   = br_matchCscGEM2.z->at(i);
      MatchedCSCStubs.push_back(tmp);
    };

    //GEMDigis
    const unsigned allGemDigiSize = br_allGemDigi.phi->size();
    for(unsigned i = 0; i < allGemDigiSize; ++i) {
      GEMDigi tmp;
      tmp.phi     = br_allGemDigi.phi->at(i);
      tmp.eta     = br_allGemDigi.eta->at(i);
      tmp.r       = br_allGemDigi.r->at(i);
      tmp.z       = br_allGemDigi.z->at(i);
      tmp.matchTp = -1;
      tmp.layer   = GEMlayerIndex(br_allGemDigi.z->at(i));
      GEMDigis.push_back(tmp);
    };

    //MatchedGEMDigis
    const unsigned matchGemDigiSize = br_matchGemDigi.phi->size();
    for(unsigned i = 0; i < matchGemDigiSize; ++i) {
      if(abs(br_tp.pdgid->at(br_matchGemDigi.matchTp->at(i)))!=13) continue;
      GEMDigi tmp;
      tmp.phi     = br_matchGemDigi.phi->at(i);
      tmp.eta     = br_matchGemDigi.eta->at(i);
      tmp.r       = br_matchGemDigi.r->at(i);
      tmp.z       = br_matchGemDigi.z->at(i);
      tmp.matchTp = MuonTPindex(br_matchGemDigi.matchTp->at(i), MuonTPs);
      tmp.layer   = GEMlayerIndex(br_matchGemDigi.z->at(i));
      MatchedGEMDigis.push_back(tmp);
    };

    //GEMPadDigis
    const unsigned GEMPadDigiSize = br_gemPadDigi.phi->size();
    for(unsigned i = 0; i < GEMPadDigiSize; ++i) {
      GEMPadDigi tmp;
      tmp.phi = br_gemPadDigi.phi->at(i);
      tmp.eta = br_gemPadDigi.eta->at(i);
      tmp.r   = br_gemPadDigi.r->at(i);
      tmp.z   = br_gemPadDigi.z->at(i);
      if (br_gemPadDigi.matchCSC->at(i) != i) cout << "allCSCGEM2 index inconsistent with allCSCStub" << endl;
      GEMPadDigis.push_back(tmp);
    };

  }

  void Reset() {
    MuonTPs.clear();
    GEMSimHits.clear();
    CSCSimHits.clear();
    CSCStubs.clear();
    MatchedCSCStubs.clear();
    GEMDigis.clear();
    MatchedGEMDigis.clear();
    GEMPadDigis.clear();

    vector<TPcontent> TPdummy = {};
    StationInfo = vector<vector<StationContent> >{
      {
        {"1-1", TPdummy, {}, {}},
        {"1-2", TPdummy, {}, {}},
        {"1-3", TPdummy, {}, {}}
      },
      {
        {"2-1", TPdummy, {}, {}},
        {"2-2", TPdummy, {}, {}}
      },
      {
        {"3-1", TPdummy, {}, {}},
        {"3-2", TPdummy, {}, {}}
      },
      {
        {"4-1", TPdummy, {}, {}},
        {"4-2", TPdummy, {}, {}}
      }
    };
  }

  void StationSort() {
    SimHitAve SimHitAve0 = SimHit{0,0,0,0};

    //GEMSimHits
    const unsigned GEMSimHitSize = GEMSimHits.size();
    for (unsigned i = 0; i < GEMSimHitSize; ++i) {
      pair<unsigned, unsigned> DaR = DiskAndRing(GEMSimHits[i].r, GEMSimHits[i].z);
      int SortedIndex = -1;
      vector<TPcontent> &current = StationInfo[DaR.first][DaR.second].TPinfoInStation;
      for (unsigned t = 0; t < current.size(); ++t) {
        if (current[t].tpIndex == GEMSimHits[i].MatchTp) {
          SortedIndex = t;
          break;
        }
      }
      if (SortedIndex == -1) {
        TPcontent tmp;
        tmp.tpIndex = GEMSimHits[i].MatchTp;
        tmp.NSimHitsCSC = 0;
        tmp.NSimHitsGEM = 1;
        tmp.CSCSimHitAve = SimHitAve0;
        tmp.GEMSimHitAve = GEMSimHits[i].ToAve();
        tmp.CSCSimHitIndices.clear();
        tmp.GEMSimHitIndices = vector<unsigned>{i};
        tmp.MatchedCSCStubIndices.clear();
        tmp.MatchedGEMDigiIndices.clear();
        current.push_back(tmp);
      }
      else {
        float N = current[SortedIndex].NSimHitsGEM;
        current[SortedIndex].GEMSimHitAve.phi = TVector2::Phi_mpi_pi((N * current[SortedIndex].GEMSimHitAve.phi + GEMSimHits[i].phi)/(N + 1.));
        current[SortedIndex].GEMSimHitAve.eta = (N * current[SortedIndex].GEMSimHitAve.eta + GEMSimHits[i].eta)/(N + 1.);
        current[SortedIndex].GEMSimHitAve.r   = (N * current[SortedIndex].GEMSimHitAve.r + GEMSimHits[i].r)/(N + 1.);
        current[SortedIndex].GEMSimHitAve.z   = (N * current[SortedIndex].GEMSimHitAve.z + GEMSimHits[i].z)/(N + 1.);
        current[SortedIndex].NSimHitsGEM += 1;
        current[SortedIndex].GEMSimHitIndices.push_back(i);
      }
    }

    //CSCSimHits
    const unsigned CSCSimHitSize = CSCSimHits.size();
    for (unsigned i = 0; i < CSCSimHitSize; ++i) {
      pair<unsigned, unsigned> DaR = DiskAndRing(CSCSimHits[i].r, CSCSimHits[i].z);
      int SortedIndex = -1;
      vector<TPcontent> &current = StationInfo[DaR.first][DaR.second].TPinfoInStation;
      for (unsigned t = 0; t < current.size(); ++t) {
        if (current[t].tpIndex == CSCSimHits[i].MatchTp) {
          SortedIndex = t;
          break;
        }
      }
      if (SortedIndex == -1) {
        TPcontent tmp;
        tmp.tpIndex = CSCSimHits[i].MatchTp;
        tmp.NSimHitsCSC = 1;
        tmp.NSimHitsGEM = 0;
        tmp.CSCSimHitAve = CSCSimHits[i].ToAve();
        tmp.GEMSimHitAve = SimHitAve0;
        tmp.CSCSimHitIndices = vector<unsigned>{i};
        tmp.GEMSimHitIndices.clear();
        tmp.MatchedCSCStubIndices.clear();
        tmp.MatchedGEMDigiIndices.clear();
        current.push_back(tmp);
      }
      else {
        float N = current[SortedIndex].NSimHitsGEM;
        current[SortedIndex].CSCSimHitAve.phi = TVector2::Phi_mpi_pi((N * current[SortedIndex].CSCSimHitAve.phi + CSCSimHits[i].phi)/(N + 1.));
        current[SortedIndex].CSCSimHitAve.eta = (N * current[SortedIndex].CSCSimHitAve.eta + CSCSimHits[i].eta)/(N + 1.);
        current[SortedIndex].CSCSimHitAve.r   = (N * current[SortedIndex].CSCSimHitAve.r + CSCSimHits[i].r)/(N + 1.);
        current[SortedIndex].CSCSimHitAve.z   = (N * current[SortedIndex].CSCSimHitAve.z + CSCSimHits[i].z)/(N + 1.);
        current[SortedIndex].NSimHitsCSC += 1;
        current[SortedIndex].CSCSimHitIndices.push_back(i);
      }
    }

    //MatchedGEMDigis
    const unsigned MatchedGEMDigiSize = MatchedGEMDigis.size();
    for (unsigned i = 0; i < MatchedGEMDigiSize; ++i) {
      pair<unsigned, unsigned> DaR = DiskAndRing(MatchedGEMDigis[i].r, MatchedGEMDigis[i].z);
      int SortedIndex = -1;
      vector<TPcontent> &current = StationInfo[DaR.first][DaR.second].TPinfoInStation;
      for (unsigned t = 0; t < current.size(); ++t) {
        if (current[t].tpIndex == (unsigned)MatchedGEMDigis[i].MatchTp) {
          SortedIndex = t;
          break;
        }
      }
      if (SortedIndex == -1) {
        TPcontent tmp;
        tmp.tpIndex = MatchedGEMDigis[i].MatchTp;
        tmp.NSimHitsCSC = 0;
        tmp.NSimHitsGEM = 0;
        tmp.CSCSimHitAve = SimHitAve0;
        tmp.GEMSimHitAve = SimHitAve0;
        tmp.CSCSimHitIndices.clear();
        tmp.GEMSimHitIndices.clear();
        tmp.MatchedCSCStubIndices.clear();
        tmp.MatchedGEMDigiIndices = vector<unsigned>{i};
        current.push_back(tmp);
      }
      else current[SortedIndex].MatchedGEMDigiIndices.push_back(i);
    }

    //MatchedCSCStubs
    const unsigned MatchedCSCStubSize = MatchedCSCStubs.size();
    for (unsigned i = 0; i < MatchedCSCStubSize; ++i) {
      pair<unsigned, unsigned> DaR = DiskAndRing(MatchedCSCStubs[i].r, MatchedCSCStubs[i].z);
      int SortedIndex = -1;
      vector<TPcontent> &current = StationInfo[DaR.first][DaR.second].TPinfoInStation;
      for (unsigned t = 0; t < current.size(); ++t) {
        if (current[t].tpIndex == (unsigned)MatchedCSCStubs[i].MatchTp) {
          SortedIndex = t;
          break;
        }
      }
      if (SortedIndex == -1) {
        TPcontent tmp;
        tmp.tpIndex = MatchedCSCStubs[i].MatchTp;
        tmp.NSimHitsCSC = 0;
        tmp.NSimHitsGEM = 0;
        tmp.CSCSimHitAve = SimHitAve0;
        tmp.GEMSimHitAve = SimHitAve0;
        tmp.CSCSimHitIndices.clear();
        tmp.GEMSimHitIndices.clear();
        tmp.MatchedCSCStubIndices = vector<unsigned>{i};
        tmp.MatchedGEMDigiIndices.clear();
        current.push_back(tmp);
      }
      else current[SortedIndex].MatchedCSCStubIndices.push_back(i);
    }

    //Generic CSCStubs
    const unsigned CSCStubSize = CSCStubs.size();
    for (unsigned i = 0; i < CSCStubSize; ++i) {
      pair<unsigned, unsigned> DaR = DiskAndRing(CSCStubs[i].r, CSCStubs[i].z);
      StationInfo[DaR.first][DaR.second].CSCStubIndices.push_back(i);
    }

    //Generic GEMDigis
    const unsigned GEMDigiSize = GEMDigis.size();
    for (unsigned i = 0; i < GEMDigiSize; ++i) {
      pair<unsigned, unsigned> DaR = DiskAndRing(GEMDigis[i].r, GEMDigis[i].z);
      StationInfo[DaR.first][DaR.second].GEMDigiIndices.push_back(i);
    }

    const unsigned GEMPadDigiSize = GEMPadDigis.size();
    for (unsigned i = 0; i < GEMPadDigiSize; ++i) {
      pair<unsigned, unsigned> DaR = DiskAndRing(GEMPadDigis[i].r, GEMPadDigis[i].z);
      StationInfo[DaR.first][DaR.second].GEMPadDigiIndices.push_back(i);
    }

  }

  vector<tp> MuonTPs;
  vector<SimHit> GEMSimHits, CSCSimHits;
  vector<CSCStub> CSCStubs, MatchedCSCStubs;
  vector<GEMDigi> GEMDigis, MatchedGEMDigis;
  vector<GEMPadDigi> GEMPadDigis;

  vector<vector<StationContent> > StationInfo;

  Branch_Reader br_tp, br_cscSimHit, br_gemSimHit;
  Branch_Reader br_allCscStubsLCT, br_allCscStubsALCT, br_allCscStubsCLCT;
  Branch_Reader br_allALCT, br_allCLCT, br_allGemDigi;
  Branch_Reader br_matchCscStubsLCT, br_matchCscStubsALCT, br_matchCscStubsCLCT, br_matchGemDigi;
  Branch_Reader br_allCscGEM1, br_allCscGEM2, br_matchCscGEM1, br_matchCscGEM2, br_gemPadDigi;

private:
  TChain *eventTree;
};

int MuonTPindex(unsigned Original_, vector<tp> MuonTPs_){
  const unsigned TPsize=MuonTPs_.size();
  for(unsigned tp=0; tp<TPsize; ++tp){
    if(MuonTPs_[tp].Index==Original_) return tp;
  }
  throw runtime_error("MuonTPindex not found");
}

pair<unsigned,unsigned> DiskAndRing(float r_, float z_){
  pair<unsigned,unsigned> DiskRingPair=make_pair(999,999);
  const float z=fabs(z_);
  const float r=fabs(r_);
  if     (z<640 && z>560 && r<300) DiskRingPair=make_pair(0,0); //1/1
  else if(z<740 && z>660){
    DiskRingPair.first=0;
    if     (r<480 && r>270) DiskRingPair.second=1; //1/2
    else if(r<700 && r>500) DiskRingPair.second=2; //1/3
  }
  else if(z<860 && z>790){
    DiskRingPair.first=1;
    if     (r<360 && r>120) DiskRingPair.second=0; //2/1
    else if(r<700 && r>360) DiskRingPair.second=1; //2/2
  }
  else if(z<960 && z>900){
    DiskRingPair.first=2;
    if     (r<360 && r>120) DiskRingPair.second=0; //3/1
    else if(r<700 && r>360) DiskRingPair.second=1; //3/2
  }
  else if(z<1150 && z>1000){
    DiskRingPair.first=3;
    if     (r<360 && r>120) DiskRingPair.second=0; //4/1
    else if(r<700 && r>360) DiskRingPair.second=1; //4/2
  }
  // else if (z < 560) {
  //   DiskRingPair = make_pair(4,0);
  // }
  if(DiskRingPair.first==999 || DiskRingPair.second==999) {
    cout << "r = " << r << ", z = " << z << endl;
    throw runtime_error("DiskAndRing invalid r-z coordinate");
  }
  return DiskRingPair;
}

//helper to return layer index for GEM detectors
int GEMlayerIndex(float z_){
  const float z=fabs(z_);
  //GE1/1 layers
  if     (z<566 && z>565) return 0;
  else if(z<568.5 && z>567.5) return 1;
  else if(z<570 && z>569) return 2;
  else if(z<572 && z>571) return 3;
  //GE2/1 layers
  else if(z<794 && z>793) return 0;
  else if(z<796.5 && z<795.5) return 1;
  else if(z<797.5 && z>796.5) return 2;
  else if(z<800 && z>799) return 3;
  else return -1;
}

//helper to return layer index for CSC detectors
int CSClayerIndex(float z_){
  const float z=fabs(z_);
  //ME1
  if     (z<588 && z>585) return 0; //ME1/1 innermost
  else if(z<617 && z>614) return 1; //ME1/1 outermost
  else if(z<685 && z>683) return 2; //ME1/2 innermost
  else if(z<695 && z>693) return 3; //ME1/3 single
  else if(z<713 && z>713) return 4; //ME1/2 outermost
  //ME2
  else if(z<816 && z>814) return 0; //ME2/1-2 innermost
  else if(z<841 && z>839) return 1; //ME2/1-2 outermost
  //ME3
  else if(z<925 && z>923) return 0; //ME3/1-2 innermost
  else if(z<950 && z>948) return 1; //ME3/1-2 outermost
  //ME4
  else if(z<1015 && z>1012) return 0; //ME4/1-2 innermost
  else if(z<1040 && z>1037) return 1; //ME4/1-2 outermost
  else return -1;
}

double CalcdR(double eta1, double eta2, double phi1, double phi2) {
  double etadiff = fabs(eta1 - eta2);
  double phidiff = TVector2::Phi_mpi_pi(phi1 - phi2);
  double dr = sqrt(etadiff * etadiff + phidiff * phidiff);
  return dr;
}

double CalcdR(CSCStub stub1, CSCStub stub2) {
  return CalcdR(stub1.eta, stub2.eta, stub1.phi, stub2.phi);
}

SimHitAve ToSimHitAve(SimHit sh) {
  SimHitAve out;
  out.phi = sh.phi;
  out.eta = sh.eta;
  out.r   = sh.r;
  out.z   = sh.z;
  return out;
}
