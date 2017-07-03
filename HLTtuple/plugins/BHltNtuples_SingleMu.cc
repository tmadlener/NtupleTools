#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
// #include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

// #include "DataFormats/PatCandidates/interface/Muon.h"

#include "NtupleTools/HLTtuple/src/hltNtupleEvent.h"

#include <map>
#include "TTree.h"
#include <vector>
#include <string>

class BHltNtuples_SingleMu : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  BHltNtuples_SingleMu(const edm::ParameterSet& cfg);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  virtual void endRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  virtual void beginEvent();


private:
  edm::Service<TFileService> m_outfile;

  std::map<std::string, TTree*> m_outTree{};
  hltNtupleEvent_singleMu m_event;

  // edm::InputTag m_puTag;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> m_puToken;
  edm::EDGetTokenT<reco::MuonCollection> m_offlineMuonsToken;
  edm::EDGetTokenT<reco::VertexCollection> m_offlinePVToken;
  // edm::EDGetTokenT<pat::MuonCollection> m_offlineMuonsToken; // relval
  edm::EDGetTokenT<edm::TriggerResults> m_trigResultsToken;
  edm::EDGetTokenT<l1t::MuonBxCollection> m_l1CandToken;
  edm::EDGetTokenT<std::vector<reco::RecoChargedCandidate>> m_l2CollToken;
  edm::EDGetTokenT<std::vector<reco::RecoChargedCandidate>> m_l3CollToken;

  static constexpr double MuMass = 0.106;
  static constexpr double MuMass2 = MuMass*MuMass;

  // Services
  edm::ESHandle<MagneticField>         m_magneticField;
  edm::ESHandle<Propagator>            m_propagator;

  std::vector<std::string> m_trigPathNames;

  double m_matchingdR;

  PropagateToMuon m_muonPropagator;
  // L1MuonMatcherAlgo m_l1MuMatcher;

  // attempts to get more trigger info
  // edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> m_triggerObjectsToken;

  // edm::EDGetTokenT<trigger::TriggerEvent> m_trigEventToken;
  // edm::EDGetTokenT<trigger::TriggerEventWithRefs> m_trigEventWRefsToken;

  void fillL1Muons(const edm::Handle<l1t::MuonBxCollection>& l1cands,
                   const edm::Event& event);

  /*int*/ double closestL1Muon(const TrajectoryStateOnSurface& state, const int otherMatch = -1);

  void fillOfflineMuons(const edm::Handle<reco::MuonCollection>& muons,
                        const edm::Handle<edm::TriggerResults>& trigResults,
                        const edm::Handle<std::vector<reco::RecoChargedCandidate>>& l2Muons,
                        const edm::Handle<std::vector<reco::RecoChargedCandidate>>& l3Muons,
                        const edm::Handle<reco::VertexCollection>& vtxColl,
                        const edm::Event& event);

  std::vector<int> getTriggerResults(const edm::Event& event,
                                     const edm::Handle<edm::TriggerResults>& trigResults);


  void fillMuons(const edm::Handle<std::vector<reco::RecoChargedCandidate>>& muons,
                 std::vector<MuonCand>& candColl);

  /*int*/ double bestMatchingMuon(const reco::Muon& muon, const std::vector<MuonCand>& coll, const double dR);
};

void BHltNtuples_SingleMu::endJob() {}

void BHltNtuples_SingleMu::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}

void BHltNtuples_SingleMu::endRun  (const edm::Run & run, const edm::EventSetup & eventSetup) {}

void BHltNtuples_SingleMu::beginEvent()
{
  m_event.L1muons.clear();
  m_event.L2muons.clear();
  m_event.L3muons.clear();
  m_event.muons.clear();
  m_event.trigger.clear();
  m_event.nVtx = -1;
  m_event.trueNI = -1;
}


BHltNtuples_SingleMu::BHltNtuples_SingleMu(const edm::ParameterSet& cfg) :
  m_puToken(consumes<std::vector<PileupSummaryInfo>>(cfg.getUntrackedParameter<edm::InputTag>("puInfoTag"))),
  m_offlineMuonsToken(consumes<reco::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("OfflineMuonsTag"))),
  m_offlinePVToken(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("offlineVtx"))),
  m_trigResultsToken(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("trigResults"))),
  m_l1CandToken(consumes<l1t::MuonBxCollection>(cfg.getParameter<edm::InputTag>("L1Candidates"))),
  m_l2CollToken(consumes<std::vector<reco::RecoChargedCandidate>>(cfg.getParameter<edm::InputTag>("l2Cands"))),
  m_l3CollToken(consumes<std::vector<reco::RecoChargedCandidate>>(cfg.getParameter<edm::InputTag>("l3Cands"))),
  m_trigPathNames(cfg.getParameter<std::vector<std::string>>("triggerPaths")),
  m_matchingdR(cfg.getParameter<double>("dRcut")),
  m_muonPropagator(cfg.getUntrackedParameter<edm::ParameterSet>("muonPropagator"))// ,
  // m_triggerObjectsToken(consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("trigObjects"))),
  // m_trigEventToken(constumes<trigger::TriggerEvent>(cfg.getParameter<edm::InputTag>("triggerEvent"))),
  // m_trigEventWRefsToken(consumes<trigger::TriggerEventWithRefs>(cfg.getParameter<edm::InputTag>("triggerEventWRefs")))
  // m_l1MuMatcher(cfg)
{
  usesResource("TFileService");
}

void BHltNtuples_SingleMu::beginJob()
{
  m_outTree["ntupleTree"] = m_outfile->make<TTree>("ntupleTree", "ntupleTree");
  m_outTree["ntupleTree"]->Branch("event", &m_event);
}

void BHltNtuples_SingleMu::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  m_muonPropagator.init(eventSetup);
  beginEvent();

  m_event.runNumber = event.id().run();
  m_event.lumiBlock = event.id().luminosityBlock();
  m_event.eventNumber = event.id().event();

  edm::Handle<reco::VertexCollection> vtxColl;
  event.getByToken(m_offlinePVToken, vtxColl);

  edm::Handle<reco::MuonCollection> muons;
  // edm::Handle<pat::MuonCollection> muons; // relval
  event.getByToken(m_offlineMuonsToken, muons);

  edm::Handle<edm::TriggerResults> trigRes;
  event.getByToken(m_trigResultsToken, trigRes);

  edm::Handle<l1t::MuonBxCollection> l1cands;
  event.getByToken(m_l1CandToken, l1cands);

  edm::Handle<std::vector<reco::RecoChargedCandidate>> l3Muons;
  event.getByToken(m_l3CollToken, l3Muons);

  edm::Handle<std::vector<reco::RecoChargedCandidate>> l2Muons;
  event.getByToken(m_l2CollToken, l2Muons);

  fillL1Muons(l1cands, event);
  fillMuons(l2Muons, m_event.L2muons);
  fillMuons(l3Muons, m_event.L3muons);
  fillOfflineMuons(muons, trigRes, l2Muons, l3Muons, vtxColl, event);

  m_outTree["ntupleTree"]->Fill();
}

void BHltNtuples_SingleMu::fillOfflineMuons(const edm::Handle<reco::MuonCollection>& muons,
                                            const edm::Handle<edm::TriggerResults>& trigResults,
                                            const edm::Handle<std::vector<reco::RecoChargedCandidate>>& l2Muons,
                                            const edm::Handle<std::vector<reco::RecoChargedCandidate>>& l3Muons,
                                            const edm::Handle<reco::VertexCollection>& vtxColl,
                                            const edm::Event& event)
{
  const auto trigRes = getTriggerResults(event, trigResults);
  m_event.trigger = trigRes;
  const reco::VertexCollection pvColl  = *(vtxColl.product())  ;

  for (const auto& mu : *muons) {
    if (!(muon::isSoftMuon(mu, pvColl[0]) && mu.pt())) continue;
    auto muStateOnSurf = m_muonPropagator.extrapolate(mu);
    auto bestMu1Match = closestL1Muon(muStateOnSurf, -1);
    // muon HLT uses different dR cuts, but they are matching to gen objects!
    auto bestL2Match = bestMatchingMuon(mu, m_event.L2muons, 0.1);
    auto bestL3Match = bestMatchingMuon(mu, m_event.L3muons, 0.1);

    // store the offline muon as Mu1 in the event
    // store the L3 muon as Mu2 in the event
    // store the L2 information in the DiMu part (neglecting the charge, which should be OK by now)
    // only need offline muons anyway

    MuonCand offlineMuon;

    offlineMuon.pt = mu.pt();
    offlineMuon.eta = mu.eta();
    offlineMuon.phi = mu.phi();
    offlineMuon.charge = mu.charge();
    offlineMuon.l1match = bestMu1Match;
    offlineMuon.l2match = bestL2Match;
    offlineMuon.l3match = bestL3Match;

    std::cout << "offlineMuon.l1match = " << offlineMuon.l1match <<
      ", offlineMuon.l2match = " << offlineMuon.l2match << ", offlineMuon.l3match = " << offlineMuon.l3match << std::endl;

    m_event.muons.push_back(offlineMuon);
  }
}

std::vector<int> BHltNtuples_SingleMu::getTriggerResults(const edm::Event& event,
                                                         const edm::Handle<edm::TriggerResults>& trigResults)
{
  std::vector<int> trigRes{};
  if (trigResults.isValid()) {
    const edm::TriggerNames& triggerNames = event.triggerNames(*trigResults);
    for (const auto& name : m_trigPathNames) {
      auto trigBit = triggerNames.triggerIndex(edm::InputTag(name).label());
      // std::cout << name << " " << trigBit << "/" << trigResults->size() << " -> ";
      if (trigBit < trigResults->size() && trigResults->accept(trigBit) && !trigResults->error(trigBit)) {
        trigRes.push_back(1);
      } else {
        trigRes.push_back(0);
      }
      // std::cout << trigRes.back() << "\n";
    }
  }

  return trigRes;
}


// int
double BHltNtuples_SingleMu::closestL1Muon(const TrajectoryStateOnSurface& state,
                                           const int otherMatch)
{
  if (state.isValid()) {
    const double muEta = state.globalPosition().eta();
    const double muPhi = state.globalPosition().phi();

    // int bestMatch = -1;
    // int secondBestMatch = -1;
    double mindR = 1e20;
    // double secondMindR = 1e20;
    for (size_t i = 0; i < m_event.L1muons.size(); ++i) { // CAUTION: need to be filled at this point!
      const auto& l1Mu = m_event.L1muons[i];
      const double dEta = muEta - l1Mu.eta;
      const double dPhi = muPhi - l1Mu.phi;
      const double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

      if (dR < mindR) {
        // secondBestMatch = bestMatch;
        // bestMatch = i;
        // secondMindR = mindR;
        mindR = dR;
      }
    }

    return mindR;
    // constexpr double maxdR = 0.2;
    //   if (mindR <  maxdR) { // assume that the seoncd best match is also below threshold
    //     // check if the other muon is already matched and return the second best match in this case
    //     // also means that it can only be matched to one, but that is already taken, this one will not be matched at all
    //     if (bestMatch == otherMatch) {
    //       if (secondMindR < 0.2) {
    //         return secondBestMatch;
    //       }
    //     } else {
    //       return bestMatch;
    //     }
    //   }
    // }
    // return -1;
  }
  return -1;
}


void BHltNtuples_SingleMu::fillL1Muons(const edm::Handle<l1t::MuonBxCollection>& l1cands ,
                                       const edm::Event& event)
{

  for (int ibx = l1cands->getFirstBX(); ibx <= l1cands->getLastBX(); ++ibx) {
    if (ibx != 0) continue;
    for (auto it = l1cands->begin(ibx); it != l1cands->end(ibx); it++){

      l1t::MuonRef muon(l1cands, distance(l1cands->begin(l1cands->getFirstBX()),it) );

      MuonCand theL1Mu;

      theL1Mu.pt       = muon -> pt();
      theL1Mu.eta      = muon -> eta();
      theL1Mu.phi      = muon -> phi();
      theL1Mu.etaAtVtx = muon -> etaAtVtx();
      theL1Mu.phiAtVtx = muon -> phiAtVtx();
      theL1Mu.charge   = muon -> charge();
      theL1Mu.quality  = muon -> hwQual();

      m_event.L1muons.push_back(theL1Mu);
    }
  }
}

void BHltNtuples_SingleMu::fillMuons(const edm::Handle<std::vector<reco::RecoChargedCandidate>>& muons,
                                     std::vector<MuonCand>& candColl)
{
  if (muons.isValid()) {
    for (const auto& mu : *muons) {
      MuonCand muCand;
      muCand.pt = mu.pt();
      muCand.eta = mu.eta();
      muCand.phi = mu.phi();
      muCand.charge = mu.charge();

      candColl.push_back(muCand);
    }
  }
}

/*int*/ double BHltNtuples_SingleMu::bestMatchingMuon(const reco::Muon& mu, const std::vector<MuonCand>& coll,
                                                      const double maxdR)
{
  const double muEta = mu.eta();
  const double muPhi = mu.phi();

  double mindR = 1e20;
  // int bestMatch = -1;
  for (size_t i = 0; i < coll.size(); ++i) {
    const double dEta = coll[i].eta - muEta;
    const double dPhi = coll[i].phi - muPhi;
    const double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

    if (dR < mindR) {
      // bestMatch = i;
      mindR = dR;
    }
  }

  return coll.empty() ? -1 : mindR;

  // // check if best matching muon is close enough
  // // constexpr double maxdR = 0.1;
  // if (mindR < maxdR) {
  //   return bestMatch;
  // }

  // return -1;
}

DEFINE_FWK_MODULE(BHltNtuples_SingleMu);
