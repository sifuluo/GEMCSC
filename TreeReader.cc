#ifndef TREEREADER_CC
#define TREEREADER_CC

#include "TreeDataFormat.cc"
#include "BranchReader.cc"
#include "Tools.cc"

#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TVector2.h"
#include "TString.h"

#include <utility>
#include <cmath>
#include <iostream>

using namespace std;

class TreeReader{
public:
  TreeReader(TChain* t_) {
    eventTree = t_;
    DebugMode = false;
    Init();
  }
  void SetDebug(bool dbg) {
    DebugMode = dbg;
  }
  void Init() {
    br_tp                    .Init(eventTree,"tp","TP",false);
    br_cscSimHit             .Init(eventTree,"cscSimHit","SimHit",true);
    br_gemSimHit             .Init(eventTree,"gemSimHit","SimHit",true);
    br_allCscStubsLCT        .Init(eventTree,"allCscStubsLCT","LCT",false);
    br_allCscStubsALCT       .Init(eventTree,"allCscStubsALCT","ALCT",false);
    br_allCscStubsCLCT       .Init(eventTree,"allCscStubsCLCT","CLCT",false);
    br_allALCT               .Init(eventTree,"allALCT","ALCT",false);
    br_allCLCT               .Init(eventTree,"allCLCT","CLCT",false);
    br_allGemDigi            .Init(eventTree,"allGemDigi","GEM",false);
    br_matchCscStubsLCT      .Init(eventTree,"matchCscStubsLCT","LCT",true);
    br_matchCscStubsALCT     .Init(eventTree,"matchCscStubsALCT","ALCT",false);
    br_matchCscStubsCLCT     .Init(eventTree,"matchCscStubsCLCT","CLCT",false);
    br_matchGemDigi          .Init(eventTree,"matchGemDigi","GEM",true);
    br_matchCscGEM1          .Init(eventTree,"matchCscGEM1","GEMPad",true);
    br_matchCscGEM2          .Init(eventTree,"matchCscGEM2","GEMPad",true);
    br_allCscGEM1            .Init(eventTree,"allCscGEM1","GEMPad",true);
    br_allCscGEM2            .Init(eventTree,"allCscGEM2","GEMPad",true);
    br_matchGemPadDigi       .Init(eventTree,"matchGemPadDigi","GEMPad",true);
    br_allGemPadDigi         .Init(eventTree,"allGemPadDigi","GEMPad",false);
    br_matchGemPadDigiCluster.Init(eventTree,"matchGemPadDigiCluster","GEMPadDigiCluster",true);
    br_allGemPadDigiCluster  .Init(eventTree,"allGemPadDigiCluster","GEMPadDigiCluster",false);
  }

  DetId ReadDet(BranchReader br, unsigned i, bool isGEM) {
    DetId tmp;
    tmp.detId   = br.detId->at(i);
    tmp.zendcap = br.zendcap->at(i);
    tmp.ring    = br.ring->at(i);
    tmp.station = br.station->at(i);
    tmp.layer   = br.layer->at(i);
    tmp.chamber = br.chamber->at(i);
    if (isGEM) tmp.roll = br.roll->at(i);
    else tmp.roll = -1;
    return tmp;
  }

  void ReadTree() {
    Evt.Init(); // Reset

    //TP
    const unsigned tp_size = br_tp.eta->size();
    for (unsigned i = 0; i < tp_size; ++i) {
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
      tmp.pdgid   = br_tp.pdgid->at(i);
      tmp.pdgid   = br_tp.pdgid->at(i);
      tmp.eventid = br_tp.eventid->at(i);
      tmp.charge  = br_tp.charge->at(i);
      tmp.Index   = i;
      Evt.I_MuonTPs.push_back(tmp);
    };
    if (DebugMode) {
      cout << "Finished TP, Branch size = " << tp_size << ", Saved size = ";
      if (Evt.MuonTPs.size() != tp_size ) cout << ", Saved size = " << Evt.MuonTPs.size();
      cout <<endl;
    }

    //GEM SimHit
    const unsigned GEMSimHitSize = br_gemSimHit.phi->size();
    for (unsigned i = 0; i < GEMSimHitSize; ++i) {
      if (abs(br_tp.pdgid->at(br_gemSimHit.matchIndex->at(i))) != 13) continue;
      SimHit tmp;
      tmp.det     = br_gemSimHit.ReadDet(i);
      tmp.eta     = br_gemSimHit.eta->at(i);
      tmp.phi     = br_gemSimHit.phi->at(i);
      tmp.r       = br_gemSimHit.r->at(i);
      tmp.z       = br_gemSimHit.z->at(i);
      tmp.MatchIndex = (unsigned)br_gemSimHit.matchIndex->at(i);
      Evt.I_GEMSimHits.push_back(tmp);
    };
    if (DebugMode) cout << "Finished GEMSimHits, Starting CSCSimHits" <<endl;

    //CSC SimHit
    const unsigned CSCSimHitSize = br_cscSimHit.phi->size();
    for(unsigned i = 0; i < CSCSimHitSize; ++i) {
      if(abs(br_tp.pdgid->at(br_cscSimHit.matchIndex->at(i))) != 13) continue;//remove SimHits not from muons
      SimHit tmp;
      tmp.det     = br_cscSimHit.ReadDet(i);
      tmp.phi     = br_cscSimHit.phi->at(i);
      tmp.eta     = br_cscSimHit.eta->at(i);
      tmp.r       = br_cscSimHit.r->at(i);
      tmp.z       = br_cscSimHit.z->at(i);
      tmp.MatchIndex = (unsigned)br_cscSimHit.matchIndex->at(i);
      // Evt.AddCSCSimHit(tmp);
      Evt.I_CSCSimHits.push_back(tmp);
    };
    if (DebugMode) cout << "Finished CSCSimHits, Starting CSCStubs" <<endl;

    //CSCStubs
    const unsigned allCscStubsLCTSize = br_allCscStubsLCT.phi->size();
    for(unsigned i = 0; i < allCscStubsLCTSize; ++i) {
      CSCStub tmp;
      tmp.det     = br_allCscStubsLCT.ReadDet(i);
      tmp.phi     = br_allCscStubsLCT.phi->at(i);
      tmp.eta     = br_allCscStubsLCT.eta->at(i);
      tmp.r       = br_allCscStubsLCT.r->at(i);
      tmp.z       = br_allCscStubsLCT.z->at(i);
      tmp.bend    = br_allCscStubsLCT.bend->at(i);
      tmp.pattern = br_allCscStubsLCT.pattern->at(i);
      tmp.slope   = br_allCscStubsLCT.slope->at(i);
      tmp.quality = br_allCscStubsLCT.quality->at(i);
      tmp.keywire = br_allCscStubsLCT.keywire->at(i);
      tmp.strip   = br_allCscStubsLCT.keywire->at(i);
      tmp.strip8  = br_allCscStubsLCT.strip8->at(i);
      tmp.valid   = br_allCscStubsLCT.valid->at(i);
      tmp.type    = br_allCscStubsLCT.type->at(i);
      tmp.GEM1pad         = br_allCscStubsLCT.GEM1pad->at(i);
      tmp.GEM1strip       = br_allCscStubsLCT.GEM1strip->at(i);
      tmp.GEM1strip8      = br_allCscStubsLCT.GEM1strip8->at(i);
      tmp.GEM1strip_me1a  = br_allCscStubsLCT.GEM1strip_me1a->at(i);
      tmp.GEM1strip8_me1a = br_allCscStubsLCT.GEM1strip8_me1a->at(i);
      tmp.GEM1keywire_min = br_allCscStubsLCT.GEM1keywire_min->at(i);
      tmp.GEM1keywire_max = br_allCscStubsLCT.GEM1keywire_max->at(i);
      tmp.GEM1roll        = br_allCscStubsLCT.GEM1roll->at(i);
      tmp.GEM1part        = br_allCscStubsLCT.GEM1part->at(i);
      tmp.GEM2pad         = br_allCscStubsLCT.GEM2pad->at(i);
      tmp.GEM2strip       = br_allCscStubsLCT.GEM2strip->at(i);
      tmp.GEM2strip8      = br_allCscStubsLCT.GEM2strip8->at(i);
      tmp.GEM2strip_me1a  = br_allCscStubsLCT.GEM2strip_me1a->at(i);
      tmp.GEM2strip8_me1a = br_allCscStubsLCT.GEM2strip8_me1a->at(i);
      tmp.GEM2keywire_min = br_allCscStubsLCT.GEM2keywire_min->at(i);
      tmp.GEM2keywire_max = br_allCscStubsLCT.GEM2keywire_max->at(i);
      tmp.GEM2roll        = br_allCscStubsLCT.GEM2roll->at(i);
      tmp.GEM2part        = br_allCscStubsLCT.GEM2part->at(i);
      tmp.MatchIndex = -1;
      tmp.CLCT_hits.clear();
      tmp.CLCT_positions.clear();
      tmp.ALCT_hits.clear();
      for (unsigned j = 0; j < 6; ++j) {
        tmp.CLCT_hits.push_back(br_allCscStubsCLCT.hit->at(6*i+j));
        tmp.CLCT_positions.push_back(br_allCscStubsCLCT.position->at(6*i+j));
        tmp.ALCT_hits.push_back(br_allCscStubsALCT.hit->at(6*i+j));
      }
      tmp.layer = CSClayerIndex(br_allCscStubsLCT.z->at(i));
      unsigned GEMPadIndex;
      bool found = false;
      for (unsigned j = 0; j < br_allCscGEM1.phi->size(); ++j) {
        if (br_allCscGEM1.matchIndex->at(j) < 0) continue;
        if ((unsigned)br_allCscGEM1.matchIndex->at(j) == i) {
          GEMPadIndex = j;
          found = true;
          break;
        }
      }
      if (found) {
        for (unsigned j = 0; j < br_allCscGEM2.phi->size(); ++j) {
          if (br_allCscGEM2.matchIndex->at(j) < 0) continue;
          if ((unsigned)br_allCscGEM2.matchIndex->at(j) == i) {
            if (GEMPadIndex != j) cout << " allCscStubsLCT GEM1 and GEM2 inconsistent" <<endl;
            break;
          }
        }
        tmp.GEM1.pad = br_allCscGEM1.pad->at(GEMPadIndex);
        tmp.GEM1.part = br_allCscGEM1.part->at(GEMPadIndex);
        // tmp.GEM1.phi = br_allCscGEM1.phi->at(GEMPadIndex);
        // tmp.GEM1.eta = br_allCscGEM1.eta->at(GEMPadIndex);
        // tmp.GEM1.r   = br_allCscGEM1.r->at(GEMPadIndex);
        // tmp.GEM1.z   = br_allCscGEM1.z->at(GEMPadIndex);
        tmp.GEM2.pad = br_allCscGEM2.pad->at(GEMPadIndex);
        tmp.GEM2.part = br_allCscGEM2.part->at(GEMPadIndex);
        // tmp.GEM2.phi = br_allCscGEM2.phi->at(GEMPadIndex);
        // tmp.GEM2.eta = br_allCscGEM2.eta->at(GEMPadIndex);
        // tmp.GEM2.r   = br_allCscGEM2.r->at(GEMPadIndex);
        // tmp.GEM2.z   = br_allCscGEM2.z->at(GEMPadIndex);
      }
      else {
        tmp.GEM1.pad = 255;
        // tmp.GEM1.phi = 0.;
        // tmp.GEM1.eta = 0.;
        // tmp.GEM1.r   = 0.;
        // tmp.GEM1.z   = 0.;
        tmp.GEM2.pad = 255;
        // tmp.GEM2.phi = 0.;
        // tmp.GEM2.eta = 0.;
        // tmp.GEM2.r   = 0.;
        // tmp.GEM2.z   = 0.;
      }
      // Evt.AddCSCStub(tmp);
      Evt.I_AllCSCStubs.push_back(tmp);
    };
    if (DebugMode) cout << "Finished CSCStubs, Starting matchCSCStubs" <<endl;

    // MatchedCSCStubs
    const unsigned matchCscStubsLCTSize = br_matchCscStubsLCT.phi->size();
    for(unsigned i = 0; i < matchCscStubsLCTSize; ++i) {
      if(abs(br_tp.pdgid->at(br_matchCscStubsLCT.matchIndex->at(i)))!=13) continue;
      CSCStub tmp;
      tmp.det     = br_matchCscStubsLCT.ReadDet(i);
      tmp.phi     = br_matchCscStubsLCT.phi->at(i);
      tmp.eta     = br_matchCscStubsLCT.eta->at(i);
      tmp.r       = br_matchCscStubsLCT.r->at(i);
      tmp.z       = br_matchCscStubsLCT.z->at(i);
      tmp.bend    = br_matchCscStubsLCT.bend->at(i);
      tmp.pattern = br_matchCscStubsLCT.pattern->at(i);
      tmp.slope   = br_matchCscStubsLCT.slope->at(i);
      tmp.quality = br_matchCscStubsLCT.quality->at(i);
      tmp.keywire = br_matchCscStubsLCT.keywire->at(i);
      tmp.strip   = br_matchCscStubsLCT.keywire->at(i);
      tmp.strip8  = br_matchCscStubsLCT.strip8->at(i);
      tmp.valid   = br_matchCscStubsLCT.valid->at(i);
      tmp.type    = br_matchCscStubsLCT.type->at(i);
      tmp.GEM1pad = br_matchCscStubsLCT.GEM1pad->at(i);
      tmp.GEM1part= br_matchCscStubsLCT.GEM1part->at(i);
      tmp.GEM2pad = br_matchCscStubsLCT.GEM2pad->at(i);
      tmp.GEM2part= br_matchCscStubsLCT.GEM2part->at(i);
      tmp.MatchIndex = br_matchCscStubsLCT.matchIndex->at(i);
      tmp.CLCT_hits.clear();
      tmp.CLCT_positions.clear();
      tmp.ALCT_hits.clear();
      for (unsigned j = 0; j < 6; ++j) {
        tmp.CLCT_hits.push_back(br_matchCscStubsCLCT.hit->at(6*i+j));
        tmp.CLCT_positions.push_back(br_matchCscStubsCLCT.position->at(6*i+j));
        tmp.ALCT_hits.push_back(br_matchCscStubsALCT.hit->at(6*i+j));
      }
      tmp.layer = CSClayerIndex(br_matchCscStubsLCT.z->at(i));
      unsigned GEMPadIndex;
      bool found = false;
      for (unsigned j = 0; j < br_matchCscGEM1.phi->size(); ++j) {
        if (br_matchCscGEM1.matchIndex->at(j) < 0) continue;
        if ((unsigned)br_matchCscGEM1.matchIndex->at(j) == i) {
          GEMPadIndex = j;
          found = true;
          break;
        }
      }
      if (found) {
        for (unsigned j = 0; j < br_matchCscGEM2.phi->size(); ++j) {
          if (br_matchCscGEM2.matchIndex->at(j) < 0) continue;
          if ((unsigned)br_matchCscGEM2.matchIndex->at(j) == i) {
            if (GEMPadIndex != j) cout << " matchCscStubsLCT GEM1 and GEM2 inconsistent" <<endl;
            break;
          }
        }
        tmp.GEM1.pad = br_matchCscGEM1.pad->at(GEMPadIndex);
        tmp.GEM1.part= br_matchCscGEM1.part->at(GEMPadIndex);
        // tmp.GEM1.phi = br_matchCscGEM1.phi->at(GEMPadIndex);
        // tmp.GEM1.eta = br_matchCscGEM1.eta->at(GEMPadIndex);
        // tmp.GEM1.r   = br_matchCscGEM1.r->at(GEMPadIndex);
        // tmp.GEM1.z   = br_matchCscGEM1.z->at(GEMPadIndex);
        tmp.GEM2.pad = br_matchCscGEM2.pad->at(GEMPadIndex);
        tmp.GEM2.part= br_matchCscGEM2.part->at(GEMPadIndex);
        // tmp.GEM2.phi = br_matchCscGEM2.phi->at(GEMPadIndex);
        // tmp.GEM2.eta = br_matchCscGEM2.eta->at(GEMPadIndex);
        // tmp.GEM2.r   = br_matchCscGEM2.r->at(GEMPadIndex);
        // tmp.GEM2.z   = br_matchCscGEM2.z->at(GEMPadIndex);
      }
      else {
        tmp.GEM1.pad = 255;
        // tmp.GEM1.phi = 0.;
        // tmp.GEM1.eta = 0.;
        // tmp.GEM1.r   = 0.;
        // tmp.GEM1.z   = 0.;
        tmp.GEM2.pad = 255;
        // tmp.GEM2.phi = 0.;
        // tmp.GEM2.eta = 0.;
        // tmp.GEM2.r   = 0.;
        // tmp.GEM2.z   = 0.;
      }
      // Evt.AddCSCStub(tmp);
      Evt.I_MatchCSCStubs.push_back(tmp);
    };
    if (DebugMode) cout << "Finished matchCSCStubs, Starting GEMDigis" <<endl;

    //GEMDigis
    const unsigned allGemDigiSize = br_allGemDigi.phi->size();
    for(unsigned i = 0; i < allGemDigiSize; ++i) {
      GEMDigi tmp;
      tmp.det     = br_allGemDigi.ReadDet(i);
      tmp.phi     = br_allGemDigi.phi->at(i);
      tmp.eta     = br_allGemDigi.eta->at(i);
      tmp.r       = br_allGemDigi.r->at(i);
      tmp.z       = br_allGemDigi.z->at(i);
      tmp.MatchIndex = -1;
      tmp.layer   = GEMlayerIndex(br_allGemDigi.z->at(i));
      // Evt.AddGEMDigi(tmp);
      Evt.I_AllGEMDigis.push_back(tmp);
    };
    if (DebugMode) cout << "Finished GEMDigis, Starting matchGEMDigis" <<endl;

    //MatchedGEMDigis
    const unsigned matchGemDigiSize = br_matchGemDigi.phi->size();
    for(unsigned i = 0; i < matchGemDigiSize; ++i) {
      if(abs(br_tp.pdgid->at(br_matchGemDigi.matchIndex->at(i)))!=13) continue;
      GEMDigi tmp;
      tmp.det     = br_matchGemDigi.ReadDet(i);
      tmp.phi     = br_matchGemDigi.phi->at(i);
      tmp.eta     = br_matchGemDigi.eta->at(i);
      tmp.r       = br_matchGemDigi.r->at(i);
      tmp.z       = br_matchGemDigi.z->at(i);
      tmp.MatchIndex = br_matchGemDigi.matchIndex->at(i);
      tmp.layer   = GEMlayerIndex(br_matchGemDigi.z->at(i));
      // Evt.AddGEMDigi(tmp);
      Evt.I_MatchGEMDigis.push_back(tmp);
    };
    if (DebugMode) cout << "Finished matchGEMDigis, Starting GEMPadDigis" <<endl;

    //GEMPadDigis
    const unsigned allGEMPadDigiSize = br_allGemPadDigi.phi->size();
    for(unsigned i = 0; i < allGEMPadDigiSize; ++i) {
      GEMPadDigi tmp;
      tmp.det = br_allGemPadDigi.ReadDet(i);
      tmp.phi = br_allGemPadDigi.phi->at(i);
      tmp.eta = br_allGemPadDigi.eta->at(i);
      tmp.r   = br_allGemPadDigi.r->at(i);
      tmp.z   = br_allGemPadDigi.z->at(i);
      tmp.part= br_allGemPadDigi.part->at(i);
      tmp.pad = br_allGemPadDigi.pad->at(i);
      // cout << "z = " << tmp.z << ", r = " << tmp.r <<endl;
      // Evt.AddGEMPadDigi(tmp);
      Evt.I_AllGEMPadDigis.push_back(tmp);
    };
    if (DebugMode) cout << "Finished AllGEMPadDigis, Starting MatchGEMPadDigis" <<endl;

    //MatchGEMPadDigis
    const unsigned matchGEMPadDigiSize = br_matchGemPadDigi.phi->size();
    for(unsigned i = 0; i < matchGEMPadDigiSize; ++i) {
      GEMPadDigi tmp;
      tmp.det = br_matchGemPadDigi.ReadDet(i);
      tmp.phi = br_matchGemPadDigi.phi->at(i);
      tmp.eta = br_matchGemPadDigi.eta->at(i);
      tmp.r   = br_matchGemPadDigi.r->at(i);
      tmp.z   = br_matchGemPadDigi.z->at(i);
      tmp.part= br_matchGemPadDigi.part->at(i);
      tmp.pad = br_matchGemPadDigi.pad->at(i);
      tmp.MatchIndex = br_matchGemPadDigi.matchIndex->at(i);
      // cout << "z = " << tmp.z << ", r = " << tmp.r <<endl;
      // Evt.AddGEMPadDigi(tmp);
      Evt.I_MatchGEMPadDigis.push_back(tmp);
    };
    if (DebugMode) cout << "Finished MatchGEMPadDigis, Starting AllGEMPadDigiClusters" <<endl;

    // AllGEMPadDigiClusters
    const unsigned allGemPadDigiClusterSize = br_allGemPadDigiCluster.phi->size();
    unsigned alllensum = 0;
    for (unsigned i = 0; i < allGemPadDigiClusterSize; ++i) alllensum += br_allGemPadDigiCluster.len->at(i);
    if (alllensum != br_allGemPadDigiCluster.pads->size()) cout << "Inconsistent Pads size for allGemPadDigiCluster" <<endl;
    unsigned allGemPadDigiClusterPadIndex = 0;
    for (unsigned i = 0; i < allGemPadDigiClusterSize; ++i) {
      GEMPadDigiCluster tmp;
      tmp.det = br_allGemPadDigiCluster.ReadDet(i);
      tmp.phi = br_allGemPadDigiCluster.phi->at(i);
      tmp.eta = br_allGemPadDigiCluster.eta->at(i);
      tmp.r = br_allGemPadDigiCluster.r->at(i);
      tmp.z = br_allGemPadDigiCluster.z->at(i);
      tmp.MatchIndex = -1;
      tmp.pads.clear();
      for (unsigned k = 0; k < (unsigned)br_allGemPadDigiCluster.len->at(i); ++k) {
        tmp.pads.push_back(br_allGemPadDigiCluster.pads->at(allGemPadDigiClusterPadIndex));
        allGemPadDigiClusterPadIndex++;
      }
      // Evt.AddGEMPadDigiCluster(tmp);
      Evt.I_AllGEMPadDigiClusters.push_back(tmp);
    }
    if (DebugMode) cout << "Finished AllGEMPadDigiClusters, Starting MatchGEMPadDigiClusters" <<endl;

    // MatchedGEMPadDigiClusters
    const unsigned matchGemPadDigiClusterSize = br_matchGemPadDigiCluster.phi->size();
    unsigned matchlensum = 0;
    for (unsigned i = 0; i < matchGemPadDigiClusterSize; ++i) matchlensum += br_matchGemPadDigiCluster.len->at(i);
    if (matchlensum != br_matchGemPadDigiCluster.pads->size()) cout << "Inconsistent Pads size for matchGemPadDigiCluster" <<endl;
    unsigned matchGemPadDigiClusterPadIndex = 0;
    for (unsigned i = 0; i < matchGemPadDigiClusterSize; ++i) {
      if (abs(br_tp.pdgid->at(br_matchGemPadDigiCluster.matchIndex->at(i)))!=13) continue;
      GEMPadDigiCluster tmp;
      tmp.det = br_matchGemPadDigiCluster.ReadDet(i);
      tmp.phi = br_matchGemPadDigiCluster.phi->at(i);
      tmp.eta = br_matchGemPadDigiCluster.eta->at(i);
      tmp.r = br_matchGemPadDigiCluster.r->at(i);
      tmp.z = br_matchGemPadDigiCluster.z->at(i);
      tmp.MatchIndex = br_matchGemPadDigiCluster.matchIndex->at(i);
      if (tmp.MatchIndex == -1)cout << " tmp.MatchIndex = " << tmp.MatchIndex;
      tmp.pads.clear();
      for (unsigned k = 0; k < (unsigned)br_matchGemPadDigiCluster.len->at(i); ++k) {
        tmp.pads.push_back(br_matchGemPadDigiCluster.pads->at(matchGemPadDigiClusterPadIndex));
        matchGemPadDigiClusterPadIndex++;
      }
      // Evt.AddGEMPadDigiCluster(tmp);
      Evt.I_MatchGEMPadDigiClusters.push_back(tmp);
    }
    if (DebugMode) cout << "Finished MatchGEMPadDigiClusters, Starting TPCalc" <<endl;
    Evt.Run();
    if (DebugMode) cout << "Finished Reading" <<endl;
  }

  EventData Evt;

  BranchReader br_tp, br_cscSimHit, br_gemSimHit;
  BranchReader br_allCscStubsLCT, br_allCscStubsALCT, br_allCscStubsCLCT;
  BranchReader br_allALCT, br_allCLCT, br_allGemDigi;
  BranchReader br_matchCscStubsLCT, br_matchCscStubsALCT, br_matchCscStubsCLCT, br_matchGemDigi;
  BranchReader br_allCscGEM1, br_allCscGEM2, br_matchCscGEM1, br_matchCscGEM2, br_allGemPadDigi, br_matchGemPadDigi;
  BranchReader br_matchGemPadDigiCluster, br_allGemPadDigiCluster;

  bool DebugMode;

private:
  TChain *eventTree;
};

bool DaRValidate(double r, double z, DetId det, TString info = "") {
  auto DaR = DiskAndRing(r,z);
  if (DaR.first == (det.station - 1 ) && DaR.second == (det.ring - 1)) return true;
  else {
    cout << info << " ";
    cout << Form("For r = %f, z = %f, DaR gives (%i, %i), while det gives (%i, %i)", r, z, DaR.first, DaR.second, det.station, det.ring)<<endl;
    return false;
  }
}

int MuonTPindex(unsigned Original_, vector<tp> MuonTPs_){
  const unsigned TPsize=MuonTPs_.size();
  for(unsigned tp=0; tp<TPsize; ++tp){
    if(MuonTPs_[tp].Index==Original_) return tp;
  }
  throw runtime_error("MuonTPindex not found");
}

double CalcdR(CSCStub stub1, CSCStub stub2) {
  return CalcdR(stub1.eta, stub2.eta, stub1.phi, stub2.phi)[0];
}

#endif
