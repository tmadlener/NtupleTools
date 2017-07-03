#ifndef HLTNTUPLEVENT_H__
#define HLTNTUPLEVENT_H__

#include "TROOT.h"
#include "TMath.h"
#include <vector>

class MuonCand {
public:

  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t etaAtVtx;
  Float_t phiAtVtx;
  Int_t   charge;
  Int_t   quality;

  Double_t l1match;
  Double_t l2match;
  Double_t l3match;

  MuonCand() = default;
  virtual ~MuonCand(){};

  ClassDef(MuonCand,2)

};

class hltNtupleEvent_singleMu {
public:
  int runNumber;
  int lumiBlock;
  int eventNumber;

  int nVtx;
  int trueNI;

  double instLumi;

  std::vector<MuonCand> L1muons{};
  std::vector<MuonCand> L2muons{};
  std::vector<MuonCand> L3muons{};
  std::vector<MuonCand> muons{};

  std::vector<int> trigger{};

  hltNtupleEvent_singleMu() = default;
  virtual ~hltNtupleEvent_singleMu() {};

  ClassDef(hltNtupleEvent_singleMu, 1);
};

class MuMuCand {
public:
  double Mu1Pt{};
  double Mu2Pt{};
  double Mu1Eta{};
  double Mu2Eta{};
  double Mu1Phi{};
  double Mu2Phi{};
  double Mu1Ch{};
  double Mu2Ch{};

  double Mu1Px{};
  double Mu1Py{};
  double Mu1Pz{};
  double Mu2Px{};
  double Mu2Py{};
  double Mu2Pz{};

  double MuMuMass{};
  double Pt{};
  double Eta{};
  double Rap{};
  double Phi{};
  double vProb{};

  double Px{};
  double Py{};
  double Pz{};

  double cosThHX{};
  double phiHX{};

  int mu1Match{}; /** index in the L1 muons of hltNtupleEvent. */
  int mu2Match{};

  std::vector<int> triggerDec{};

  MuMuCand() = default;
  virtual ~MuMuCand() {};

  ClassDef(MuMuCand, 1);
};

class L1MuonCand {
public:

  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t etaAtVtx;
  Float_t phiAtVtx;
  Int_t   charge;
  Int_t   quality;

  L1MuonCand(){};
  virtual ~L1MuonCand(){};

  ClassDef(L1MuonCand,1)

};

class hltNtupleEvent {
public:
  int runNumber;
  int lumiBlock;
  int eventNumber;

  int nVtx;
  int trueNI;

  double instLumi;

  std::vector<MuMuCand> muMuCands{};

  std::vector<L1MuonCand> L1muons{};

  hltNtupleEvent() = default;
  virtual ~hltNtupleEvent() {};

  ClassDef(hltNtupleEvent, 1);
};



#endif
