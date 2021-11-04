#ifndef BRANCHREADER_CC
#define BRANCHREADER_CC

#include <vector>
#include "TChain.h"
#include "TreeDataFormat.cc"


class BranchReader{
public:
  BranchReader() {

  };

  void Init(TChain* evttree, TString name_, int data_type_, bool match_ = false) {
    name = name_;
    data_type = data_type_; // 0 for LCT, 1 for ALCT, 2 for CLCT, 3 for GemDigi, 4 for GEMPad, 5 for SimHit, 6 for TP, 7 for Cluster
    IsMatched = match_;

    // LCT
    if (data_type == 0) {
      InitGP(evttree);
      InitDet(evttree);
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
      quality = new std::vector<int>;
      keywire = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      valid = new std::vector<bool>;
      type = new std::vector<int>;
      GEM1pad = new std::vector<int>;
      GEM1strip = new std::vector<int>;
      GEM1strip8 = new std::vector<int>;
      GEM1strip_me1a = new std::vector<int>;
      GEM1strip8_me1a = new std::vector<int>;
      GEM1keywire_min = new std::vector<int>;
      GEM1keywire_max = new std::vector<int>;
      GEM1roll = new std::vector<int>;
      GEM1part = new std::vector<int>;
      GEM2pad = new std::vector<int>;
      GEM2strip = new std::vector<int>;
      GEM2strip8 = new std::vector<int>;
      GEM2strip_me1a = new std::vector<int>;
      GEM2strip8_me1a = new std::vector<int>;
      GEM2keywire_min = new std::vector<int>;
      GEM2keywire_max = new std::vector<int>;
      GEM2roll = new std::vector<int>;
      GEM2part = new std::vector<int>;
      evttree->SetBranchAddress(name+"_bend", &bend);
      evttree->SetBranchAddress(name+"_pattern", &pattern);
      evttree->SetBranchAddress(name+"_slope", &slope);
      evttree->SetBranchAddress(name+"_quality", &quality);
      evttree->SetBranchAddress(name+"_keywire", &keywire);
      evttree->SetBranchAddress(name+"_strip", &strip);
      evttree->SetBranchAddress(name+"_strip8", &strip8);
      evttree->SetBranchAddress(name+"_valid", &valid);
      evttree->SetBranchAddress(name+"_type", &type);
      evttree->SetBranchAddress(name+"_GEM1pad", &GEM1pad);
      evttree->SetBranchAddress(name+"_GEM1strip", &GEM1strip);
      evttree->SetBranchAddress(name+"_GEM1strip8", &GEM1strip8);
      evttree->SetBranchAddress(name+"_GEM1strip_me1a", &GEM1strip_me1a);
      evttree->SetBranchAddress(name+"_GEM1strip8_me1a", &GEM1strip8_me1a);
      evttree->SetBranchAddress(name+"_GEM1keywire_min", &GEM1keywire_min);
      evttree->SetBranchAddress(name+"_GEM1keywire_max", &GEM1keywire_max);
      evttree->SetBranchAddress(name+"_GEM1roll", &GEM1roll);
      evttree->SetBranchAddress(name+"_GEM1part", &GEM1part);
      evttree->SetBranchAddress(name+"_GEM2pad", &GEM2pad);
      evttree->SetBranchAddress(name+"_GEM2strip", &GEM2strip);
      evttree->SetBranchAddress(name+"_GEM2strip8", &GEM2strip8);
      evttree->SetBranchAddress(name+"_GEM2strip_me1a", &GEM2strip_me1a);
      evttree->SetBranchAddress(name+"_GEM2strip8_me1a", &GEM2strip8_me1a);
      evttree->SetBranchAddress(name+"_GEM2keywire_min", &GEM2keywire_min);
      evttree->SetBranchAddress(name+"_GEM2keywire_max", &GEM2keywire_max);
      evttree->SetBranchAddress(name+"_GEM2roll", &GEM2roll);
      evttree->SetBranchAddress(name+"_GEM2part", &GEM2part);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->SetBranchAddress(name+"_matchIndex", &matchIndex);
      }
    }

    // ALCT
    else if (data_type == 1) {
      InitDet(evttree);
      keywire = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      evttree->SetBranchAddress(name+"_keywire", &keywire);
      evttree->SetBranchAddress(name+"_hit", &hit);
      evttree->SetBranchAddress(name+"_position", &position);
      evttree->SetBranchAddress(name+"_valid", &valid);
    }

    //CLCT
    else if (data_type == 2) {
      InitDet(evttree);
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
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
      InitDet(evttree);
      strip = new std::vector<int>;
      evttree->SetBranchAddress(name+"_strip", &strip);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->SetBranchAddress(name+"_matchIndex", &matchIndex);
      }
    }

    // GEMPad
    else if (data_type == 4) {
      InitGP(evttree);
      InitDet(evttree);
      pad = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      strip_me1a = new std::vector<int>;
      strip8_me1a = new std::vector<int>;
      keywire_min = new std::vector<int>;
      keywire_max = new std::vector<int>;
      part = new std::vector<int>;
      evttree->SetBranchAddress(name+"_pad", &pad);
      evttree->SetBranchAddress(name+"_strip", &strip);
      evttree->SetBranchAddress(name+"_strip8", &strip8);
      evttree->SetBranchAddress(name+"_strip_me1a", &strip_me1a);
      evttree->SetBranchAddress(name+"_strip8_me1a", &strip8_me1a);
      evttree->SetBranchAddress(name+"_keywire_min", &keywire_min);
      evttree->SetBranchAddress(name+"_keywire_max", &keywire_max);
      evttree->SetBranchAddress(name+"_part", &part);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->SetBranchAddress(name+"_matchIndex", &matchIndex);
      }
    }

    //SimHit
    else if (data_type == 5) {
      if (!IsMatched) cout << name << " not matched? SimHits are always matched!" << endl;
      InitGP(evttree);
      InitDet(evttree);
      matchIndex = new std::vector<int>;
      evttree->SetBranchAddress(name+"_matchIndex", &matchIndex);
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

    // GEMPadDigiCluster
    else if (data_type == 7) {
      InitGP(evttree);
      InitDet(evttree);
      pads = new std::vector<int>;
      pad = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      strip_me1a = new std::vector<int>;
      strip8_me1a = new std::vector<int>;
      keywire_min = new std::vector<int>;
      keywire_max = new std::vector<int>;
      part = new std::vector<int>;
      len = new std::vector<int>;
      evttree->SetBranchAddress(name+"_pads", &pads);
      evttree->SetBranchAddress(name+"_strip", &strip);
      evttree->SetBranchAddress(name+"_strip8", &strip8);
      evttree->SetBranchAddress(name+"_strip_me1a", &strip_me1a);
      evttree->SetBranchAddress(name+"_strip8_me1a", &strip8_me1a);
      evttree->SetBranchAddress(name+"_keywire_min", &keywire_min);
      evttree->SetBranchAddress(name+"_keywire_max", &keywire_max);
      evttree->SetBranchAddress(name+"_part", &part);
      evttree->SetBranchAddress(name+"_len",&len);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->SetBranchAddress(name+"_matchIndex", &matchIndex);
      }
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
    else if (data_type_st == "GEMPadDigiCluster") data_type_ = 7;
    // else if (data_type_st == "MatchMuon") data_type_ = 8;


    if (data_type_ == -1) cout << "Wrong Data Type input for " << name_ << " as " << data_type_st << endl;
    else Init(evttree, name_, data_type_, match_);
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

  void InitDet(TTree* evttree) {
    detId = new std::vector<int>;
    zendcap = new std::vector<int>;
    ring = new std::vector<int>;
    station = new std::vector<int>;
    layer = new std::vector<int>;
    chamber = new std::vector<int>;
    roll = new std::vector<int>;
    evttree->SetBranchAddress(name+"_detId", &detId);
    evttree->SetBranchAddress(name+"_zendcap", &zendcap);
    evttree->SetBranchAddress(name+"_ring", &ring);
    evttree->SetBranchAddress(name+"_station", &station);
    evttree->SetBranchAddress(name+"_layer", &layer);
    evttree->SetBranchAddress(name+"_chamber", &chamber);
    evttree->SetBranchAddress(name+"_roll", &roll);
  }

  DetId ReadDet(unsigned i){
    DetId tmp;
    tmp.detId = detId->at(i);
    tmp.zendcap = zendcap->at(i);
    tmp.ring = ring->at(i);
    tmp.station = station->at(i);
    tmp.layer = layer->at(i);
    tmp.chamber = chamber->at(i);
    if (ring->size() == 0) tmp.ring = -1;
    else tmp.ring = ring->at(i);
    if (roll->size() > 0) tmp.roll = roll->at(i);
    return tmp;
  }

  std::vector<int>*   matchIndex;

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
  std::vector<int>*   GEM1pad;
  std::vector<int>*   GEM1strip;
  std::vector<int>*   GEM1strip8;
  std::vector<int>*   GEM1strip_me1a;
  std::vector<int>*   GEM1strip8_me1a;
  std::vector<int>*   GEM1keywire_min;
  std::vector<int>*   GEM1keywire_max;
  std::vector<int>*   GEM1roll;
  std::vector<int>*   GEM1part;
  std::vector<int>*   GEM2pad;
  std::vector<int>*   GEM2strip;
  std::vector<int>*   GEM2strip8;
  std::vector<int>*   GEM2strip_me1a;
  std::vector<int>*   GEM2strip8_me1a;
  std::vector<int>*   GEM2keywire_min;
  std::vector<int>*   GEM2keywire_max;
  std::vector<int>*   GEM2roll;
  std::vector<int>*   GEM2part;
  std::vector<int>*   part;
  std::vector<int>*   pad;
  std::vector<int>*   strip_me1a;
  std::vector<int>*   strip8_me1a;
  std::vector<int>*   keywire_min;
  std::vector<int>*   keywire_max;
  std::vector<int>*   pads;
  std::vector<int>*   len;

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

  //DetInfo
  // GEMDetId int region(+/-1), int ring, int station, int layer, int chamber, int ieta(roll)
  // CSCDetId int iendcap, int istation, int iring, int ichamber, int ilayer
  std::vector<int>* detId;
  std::vector<int>* zendcap;
  std::vector<int>* ring;
  std::vector<int>* station;
  std::vector<int>* layer;
  std::vector<int>* chamber;
  std::vector<int>* roll;

private:
  TString name;
  int data_type;
  bool IsMatched;
};

#endif
