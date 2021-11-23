#ifndef TOOLS_CC
#define TOOLS_CC

#include "TROOT.h"
#include "TString.h"
#include "TVector2.h"

#include <vector>
#include <utility>
#include <cmath>

using namespace std;

pair<int,int> DiskAndRing(float r_, float z_){
  pair<int,int> DiskRingPair=make_pair(999,999);
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
  else if (z < 560) {
    DiskRingPair = make_pair(0,3);
  }
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

vector<float> CalcdR(float eta1, float eta2, float phi1, float phi2) {
  float etadiff = fabs(eta1 - eta2);
  float phidiff = TVector2::Phi_mpi_pi(phi1 - phi2);
  float dr = sqrt(etadiff * etadiff + phidiff * phidiff);
  vector<float> out{dr,etadiff,phidiff};
  return out;
}

bool IsCloseCSC(float eta1, float eta2, float phi1, float phi2) {
  if (fabs(eta1 - eta2) < 0.0181205 && TVector2::Phi_mpi_pi(phi1 - phi2) < 0.0040559) return true;
  return false;
}

bool IsCloseGEM(float eta1, float eta2, float phi1, float phi2) {
  if (fabs(eta1 - eta2) < 0.0392362 && TVector2::Phi_mpi_pi(phi1 - phi2) < 0.000623863) return true;
  return false;
}

bool IsClosePad(float eta1, float eta2, float phi1, float phi2) {
  if (fabs(eta1 - eta2) < 0.5 && TVector2::Phi_mpi_pi(phi1 - phi2) < 0.623863) return true;
  return false;
}

bool IsCloseCluster(float eta1, float eta2, float phi1, float phi2) {
  // if (fabs(eta1 - eta2) < 0.037425 && TVector2::Phi_mpi_pi(phi1 - phi2) < 0.00141607) return true;
  if (fabs(eta1 - eta2) < 0.37425 && TVector2::Phi_mpi_pi(phi1 - phi2) < 0.141607) return true;
  return false;
}

void PrintProgress (int entry, int entries, int dividen = 100) {
  entries--;
  if (entry % dividen ==0 || entry == entries ) {
    cout << Form("\rProgress: %i / %i ", entry, entries) <<flush;
  }
  if (entry == entries){
    cout <<endl << "Done."<<endl;
  }
}

#endif
