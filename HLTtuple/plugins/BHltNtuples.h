// #include "FWCore/Framework/interface/EDAnalyzer.h"
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

// #include "DataFormats/PatCandidates/interface/Muon.h"

#include "NtupleTools/HLTtuple/src/hltNtupleEvent.h"

#include <map>
#include "TTree.h"
#include <vector>
#include <string>

class BHltNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  BHltNtuples(const edm::ParameterSet& cfg);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  virtual void endRun(const edm::Run& run, const edm::EventSetup& eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  virtual void beginEvent();


private:
  edm::Service<TFileService> m_outfile;

  std::map<std::string, TTree*> m_outTree{};
  hltNtupleEvent m_event;

  // edm::InputTag m_puTag;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> m_puToken;
  edm::EDGetTokenT<reco::VertexCollection> m_offlinePVToken;
  edm::EDGetTokenT<reco::BeamSpot> m_beamspotToken;
  edm::EDGetTokenT<reco::MuonCollection> m_offlineMuonsToken;
  // edm::EDGetTokenT<pat::MuonCollection> m_offlineMuonsToken; // relval
  edm::EDGetTokenT<edm::TriggerResults> m_trigResultsToken;
  edm::EDGetTokenT<l1t::MuonBxCollection> m_l1CandToken;

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

  void fillMuMu(const edm::Handle<reco::MuonCollection>& muons,
                const edm::Handle<reco::VertexCollection>& vtxColl,
                const edm::Event& event, const edm::EventSetup& eventSetup,
                const edm::Handle<edm::TriggerResults>& trigRes
                // /*const*/ edm::Handle<pat::TriggerObjectStandAloneCollection>& trigObjs
                );

  void fillL1Muons(const edm::Handle<l1t::MuonBxCollection>& l1cands,
                   const edm::Event& event);

  int closestL1Muon(const TrajectoryStateOnSurface& state, const int otherMatch = -1);

};

void BHltNtuples::endJob() {}

void BHltNtuples::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}

void BHltNtuples::endRun  (const edm::Run & run, const edm::EventSetup & eventSetup) {}

void BHltNtuples::beginEvent()
{
  m_event.muMuCands.clear();
  m_event.nVtx = -1;
  m_event.trueNI = -1;
}
