#include "BHltNtuples.h"
#include "BHltNtuples_helper.h"

#include <iostream>
#include <cmath>

BHltNtuples::BHltNtuples(const edm::ParameterSet& cfg) :
  m_puToken(consumes<std::vector<PileupSummaryInfo>>(cfg.getUntrackedParameter<edm::InputTag>("puInfoTag"))),
  m_offlinePVToken(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("offlineVtx"))),
  m_beamspotToken(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamspot"))),
  // m_lumiScalerToken(consumes<LumiScalersCollection>(cfg.getUntrackedParameter<edm::InputTag>("lumiScalerTag"))),
  m_offlineMuonsToken(consumes<reco::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("OfflineMuonsTag"))),
  // m_offlineMuonsToken(consumes<pat::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("OfflineMuonsTag"))), // relval
  m_trigResultsToken(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("trigResults"))),
  m_l1CandToken(consumes<l1t::MuonBxCollection>(cfg.getParameter<edm::InputTag>("L1Candidates"))),
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

void BHltNtuples::beginJob()
{
  m_outTree["ntupleTree"] = m_outfile->make<TTree>("ntupleTree", "ntupleTree");
  m_outTree["ntupleTree"]->Branch("event", &m_event);

  // write the names of the paths into a separate tree to be able to retrieve them
  // in the correct order later easy retrieval
  m_outTree["triggerPaths"]->Branch("pathNames", &m_trigPathNames);
  m_outTree["triggerPaths"]->Fill();
}

void BHltNtuples::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  m_muonPropagator.init(eventSetup);
  beginEvent();

  m_event.runNumber = event.id().run();
  m_event.lumiBlock = event.id().luminosityBlock();
  m_event.eventNumber = event.id().event();

  // std::cout << "EVENT " << m_event.eventNumber << "========================================\n";

  edm::Handle<reco::VertexCollection> vtxColl;
  event.getByToken(m_offlinePVToken, vtxColl);

  size_t nGoodVtx{};
  for (auto it = vtxColl->begin(); it != vtxColl->end(); ++it) {
    if (it->isValid()) nGoodVtx++;
  }
  m_event.nVtx = nGoodVtx;

  edm::Handle<reco::MuonCollection> muons;
  // edm::Handle<pat::MuonCollection> muons; // relval
  event.getByToken(m_offlineMuonsToken, muons);

  edm::Handle<edm::TriggerResults> trigRes;
  event.getByToken(m_trigResultsToken, trigRes);

  edm::Handle<l1t::MuonBxCollection> l1cands;
  event.getByToken(m_l1CandToken, l1cands);

  fillL1Muons(l1cands, event);
  fillMuMu(muons, vtxColl, event, eventSetup, trigRes);// , trigObjects);

  // std::cout << m_event.muMuCands.size() << "\n";
  m_outTree["ntupleTree"]->Fill();
}


void BHltNtuples::fillMuMu (const edm::Handle<reco::MuonCollection> & muons ,
                            const edm::Handle<reco::VertexCollection > & vtxColl,
                            const edm::Event & event, const edm::EventSetup & eventSetup,
                            const edm::Handle<edm::TriggerResults> & trigResults
                            // /*const*/ edm::Handle<pat::TriggerObjectStandAloneCollection>& trigObjs
                            )
{
  //get offline beamspot position
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  event.getByToken(m_beamspotToken,recoBeamSpotHandle);
  const reco::BeamSpot& vertexBeamSpot = *recoBeamSpotHandle;

  //get the transient track builder:
  edm::ESHandle<TransientTrackBuilder> theB;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  const reco::VertexCollection pvColl  = *(vtxColl.product())  ;

  std::vector<int> trigRes{};
  if (trigResults.isValid()) {
    const edm::TriggerNames& triggerNames = event.triggerNames(*trigResults);
    for (const auto& name : m_trigPathNames) {
      auto trigBit = triggerNames.triggerIndex(edm::InputTag(name).label());
      // std::cout << name << " " << trigBit << "/" << trigResults->size() <<"\n";
      if (trigBit < trigResults->size() && trigResults->accept(trigBit) && !trigResults->error(trigBit)) {
        trigRes.push_back(1);
      } else {
        trigRes.push_back(0);
      }
    }
    // std::cout << trigRes.size() << " " << m_trigPathNames.size() << " " << triggerNames.triggerNames().size() << "\n";

  }

  // Loop muon collection
  for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) {
    if( muon::isSoftMuon( (*mu1), pvColl[0] ) && (*mu1).pt()) {

      auto mu1StateOnSurf = m_muonPropagator.extrapolate(*mu1);
      int bestMu1Match = closestL1Muon(mu1StateOnSurf, -1);
      for(std::vector<reco::Muon>::const_iterator mu2=mu1; mu2!=muons->end(); ++mu2) {
        if( mu2 == mu1 ) continue;
        if( muon::isSoftMuon( (*mu2), pvColl[0]) && (*mu2).pt()) {
          if(!( mu1->charge() * mu2->charge() < 0 ))         continue;

          auto mu2StateOnSurf = m_muonPropagator.extrapolate(*mu2);
          int bestMu2Match = closestL1Muon(mu2StateOnSurf, bestMu1Match);

          // std::cout << bestMu1Match << " " << bestMu2Match << std::endl;
          // do jpsi vertex fit
          std::vector<reco::TransientTrack> j_tks;
          j_tks.push_back((*theB).build(mu1->track().get()));
          j_tks.push_back((*theB).build(mu2->track().get()));
          if (j_tks.size()!=2) continue;

          KalmanVertexFitter jkvf;
          TransientVertex jtv = jkvf.vertex(j_tks);
          if (!jtv.isValid()) continue;

          reco::Vertex jpsivertex = jtv;
          float dimuonCL = 0;
          if( (jpsivertex.chi2()>=0.0) && (jpsivertex.ndof()>0) )
            dimuonCL = TMath::Prob(jpsivertex.chi2(), jpsivertex.ndof() );

          // calculate three-track transverse momentum
          math::XYZVector jpperp(mu1->px() + mu2->px() ,
                                 mu1->py() + mu2->py() ,
                                 0.);

          GlobalPoint jVertex = jtv.position();
          GlobalError jerr    = jtv.positionError();

          //calculate decay length  significance w.r.t. the beamspot
          GlobalPoint displacementFromBeamspotJpsi( -1*((vertexBeamSpot.x0() - jVertex.x()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()),
                                                    -1*((vertexBeamSpot.y0() - jVertex.y()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()),
                                                    0);
          reco::Vertex::Point vperpj(displacementFromBeamspotJpsi.x(), displacementFromBeamspotJpsi.y(), 0.);

          MuMuCand theDiMu;

          theDiMu.triggerDec = trigRes;
          theDiMu.mu1Match = bestMu1Match;
          theDiMu.mu2Match = bestMu2Match;

          theDiMu.Mu1Pt    = mu1 -> pt( );
          theDiMu.Mu1Eta   = mu1 -> eta();
          theDiMu.Mu1Phi   = mu1 -> phi();
          theDiMu.Mu1Ch    = mu1 -> charge();
          theDiMu.Mu2Pt    = mu2 -> pt( );
          theDiMu.Mu2Eta   = mu2 -> eta();
          theDiMu.Mu2Phi   = mu2 -> phi();
          theDiMu.Mu2Ch    = mu2 -> charge();

          // 15.06.2017:
          theDiMu.Mu1Px = mu1->px();
          theDiMu.Mu1Py = mu1->py();
          theDiMu.Mu1Pz = mu1->pz();

          theDiMu.Mu2Px = mu2->px();
          theDiMu.Mu2Py = mu2->py();
          theDiMu.Mu2Pz = mu2->pz();

          const auto dimu4Vector = mu1->p4() + mu2->p4();

          theDiMu.MuMuMass = dimu4Vector.mass();
          theDiMu.vProb   = dimuonCL;

          theDiMu.Pt         = jpperp.R();

          theDiMu.Rap        = dimu4Vector.Rapidity();
          theDiMu.Eta = dimu4Vector.Eta();
          theDiMu.Phi = dimu4Vector.Phi();

          theDiMu.Px = dimu4Vector.Px();
          theDiMu.Py = dimu4Vector.Py();
          theDiMu.Pz = dimu4Vector.Pz();

          TLorentzVector jpsi;
          jpsi.SetPtEtaPhiM(theDiMu.Pt, theDiMu.Eta, theDiMu.Phi, theDiMu.MuMuMass);
          TLorentzVector muPlus;
          if (mu1->charge() > 0) {
            muPlus.SetPtEtaPhiM(theDiMu.Mu1Pt, theDiMu.Mu1Eta, theDiMu.Mu1Phi, 0.105658);
          } else {
            muPlus.SetPtEtaPhiM(theDiMu.Mu2Pt, theDiMu.Mu2Eta, theDiMu.Mu2Phi, 0.105658);
          }
          const auto angles = calcAnglesHX(jpsi, muPlus);

          theDiMu.cosThHX = angles.costh;
          theDiMu.phiHX = angles.phi;

          m_event.muMuCands.push_back(theDiMu);
        }
      }
    }
  }
}

int BHltNtuples::closestL1Muon(const TrajectoryStateOnSurface& state,
                               const int otherMatch)
{
  if (state.isValid()) {
    const double muEta = state.globalPosition().eta();
    const double muPhi = state.globalPosition().phi();

    int bestMatch = -1;
    int secondBestMatch = -1;
    double mindR = 1e20;
    for (size_t i = 0; i < m_event.L1muons.size(); ++i) { // CAUTION: need to be filled at this point!
      const auto& l1Mu = m_event.L1muons[i];
      const double dEta = muEta - l1Mu.eta;
      const double dPhi = muPhi - l1Mu.phi;
      const double dR = std::sqrt(dEta*dEta + dPhi*dPhi);

      if (dR < m_matchingdR) {
        if (dR < mindR) {
          secondBestMatch = bestMatch;
          bestMatch = i;
          mindR = dR;
        }
      }
    }

    // check if the other muon is already matched and return the second best match in this case
    // also means that it can only be matched to one, but that is already taken, this one will not be matched at all
    if (bestMatch == otherMatch) bestMatch = secondBestMatch;
    return bestMatch;
  }

  return -1;
}

void BHltNtuples::fillL1Muons(const edm::Handle<l1t::MuonBxCollection>& l1cands ,
                              const edm::Event& event)
{

  for (int ibx = l1cands->getFirstBX(); ibx <= l1cands->getLastBX(); ++ibx) {
    if (ibx != 0) continue;
    for (auto it = l1cands->begin(ibx); it != l1cands->end(ibx); it++){

      l1t::MuonRef muon(l1cands, distance(l1cands->begin(l1cands->getFirstBX()),it) );

      L1MuonCand theL1Mu;

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


DEFINE_FWK_MODULE(BHltNtuples);
