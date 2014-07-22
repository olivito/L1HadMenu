// -*- C++ -*-
//
// Package:    L1TrackTriggerObjectsAnalyzer 
// Class:      L1TrackTriggerObjectsAnalyzer
//
/**\class L1TrackTriggerObjectsAnalyzer L1TrackTriggerObjectsAnalyzer.cc SLHCUpgradeSimulations/L1TrackTriggerObjectsAnalyzer/src/L1TrackTriggerObjectsAnalyzer.cc 
 Description: [one line class summary] 

 Implementation: 
     [Notes on implementation]
*/
//
// Original Author:  Maria
//         Created:  Thu Nov 14 11:22:13 CET 2013 
// $Id$ 
//
//

// system include files 
#include <memory>
//#include <string>

// user include files 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

#include "TFile.h"
#include "TMath.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;

using namespace reco;
using namespace std;

typedef math::XYZTLorentzVector LorentzVector;

//////// 
// class declaration
////

/*
class cmpPtGT {
public:
  bool operator() (const L1TkJetParticle &r1, const L1TkJetParticle &r2) {
    return r1.pt()>r2.pt(); 
  }
};
*/

class cmpPtLorentz {
public:
  bool operator() (const LorentzVector &r1, const LorentzVector &r2) {
    return r1.pt()>r2.pt();
  }
};

/////  

class L1HadMenuDecisionProducer : public edm::EDProducer {
public:

  explicit L1HadMenuDecisionProducer(const edm::ParameterSet&);
  ~L1HadMenuDecisionProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  float dphi(float phi1, float phi2);

  // ----------member data ---------------------------

  const std::string aliasprefix_ = "l1";
                                                                                                                                          
  edm::InputTag L1JetsCentInputTag;
  edm::InputTag L1JetsFwdInputTag;
  edm::InputTag L1EtMissInputTag;
  edm::InputTag L1MHTInputTag;

};

//
// constants, enums and typedefs 
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1HadMenuDecisionProducer::L1HadMenuDecisionProducer(const edm::ParameterSet& iConfig)

{
  L1JetsCentInputTag = iConfig.getParameter<edm::InputTag>("L1JetsCentInputTag");
  L1JetsFwdInputTag = iConfig.getParameter<edm::InputTag>("L1JetsFwdInputTag");
  L1EtMissInputTag = iConfig.getParameter<edm::InputTag>("L1EtMissInputTag");
  L1MHTInputTag= iConfig.getParameter<edm::InputTag>("L1MHTInputTag");

  // output branches
  produces<bool>  (aliasprefix_+"HTT95DoubleJetC30Dphi8bins").setBranchAlias(aliasprefix_+"_HTT95DoubleJetC30Dphi8bins");
  produces<bool>  (aliasprefix_+"HTT125DoubleJetC30Dphi8bins").setBranchAlias(aliasprefix_+"_HTT125DoubleJetC30Dphi8bins");
  produces<bool>  (aliasprefix_+"HTT175").setBranchAlias(aliasprefix_+"_HTT175");
  produces<bool>  (aliasprefix_+"HTT200").setBranchAlias(aliasprefix_+"_HTT200");
  produces<bool>  (aliasprefix_+"DoubleJetC100").setBranchAlias(aliasprefix_+"_DoubleJetC100");
  produces<bool>  (aliasprefix_+"DoubleJetC120").setBranchAlias(aliasprefix_+"_DoubleJetC120");
  produces<bool>  (aliasprefix_+"DoubleJetC60MET60").setBranchAlias(aliasprefix_+"_DoubleJetC60MET60");
  produces<bool>  (aliasprefix_+"QuadJetC60").setBranchAlias(aliasprefix_+"_QuadJetC60");
  produces<bool>  (aliasprefix_+"SingleJet200").setBranchAlias(aliasprefix_+"_SingleJet200");

}


L1HadMenuDecisionProducer::~L1HadMenuDecisionProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1HadMenuDecisionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  produces<bool>  (aliasprefix_+"HTT95DoubleJetC30Dphi8bins").setBranchAlias(aliasprefix_+"_HTT95DoubleJetC30Dphi8bins");
  produces<bool>  (aliasprefix_+"HTT125DoubleJetC30Dphi8bins").setBranchAlias(aliasprefix_+"_HTT125DoubleJetC30Dphi8bins");
  produces<bool>  (aliasprefix_+"HTT175").setBranchAlias(aliasprefix_+"_HTT175");
  produces<bool>  (aliasprefix_+"HTT200").setBranchAlias(aliasprefix_+"_HTT200");
  produces<bool>  (aliasprefix_+"DoubleJetC100").setBranchAlias(aliasprefix_+"_DoubleJetC100");
  produces<bool>  (aliasprefix_+"DoubleJetC120").setBranchAlias(aliasprefix_+"_DoubleJetC120");
  produces<bool>  (aliasprefix_+"DoubleJetC60MET60").setBranchAlias(aliasprefix_+"_DoubleJetC60MET60");
  produces<bool>  (aliasprefix_+"QuadJetC60").setBranchAlias(aliasprefix_+"_QuadJetC60");
  produces<bool>  (aliasprefix_+"SingleJet200").setBranchAlias(aliasprefix_+"_SingleJet200");

  // ----------------------------
  // branch variable defs
  auto_ptr<bool>    HTT95DoubleJetC30Dphi8bins   (new bool);
  auto_ptr<bool>    HTT125DoubleJetC30Dphi8bins   (new bool);
  auto_ptr<bool>    HTT175   (new bool);
  auto_ptr<bool>    HTT200   (new bool);
  auto_ptr<bool>    DoubleJetC100   (new bool);
  auto_ptr<bool>    DoubleJetC120   (new bool);
  auto_ptr<bool>    DoubleJetC60MET60   (new bool);
  auto_ptr<bool>    QuadJetC60   (new bool);
  auto_ptr<bool>    SingleJet200   (new bool);
  // ----------------------------

  //
  // ----------------------------------------------------------------------
  // retrieve the jets objects
  //

  edm::Handle<L1JetParticleCollection> L1JetsCentHandle;
  iEvent.getByLabel(L1JetsCentInputTag, L1JetsCentHandle);
  edm::Handle<L1JetParticleCollection> L1JetsFwdHandle;
  iEvent.getByLabel(L1JetsFwdInputTag, L1JetsFwdHandle);

  std::vector<LorentzVector> jets;
  std::vector<LorentzVector> jetsCent;

  if ( L1JetsCentHandle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsCentHandle -> begin(); jetIter != L1JetsCentHandle->end(); ++jetIter) {
      jets.push_back(jetIter->p4());
      jetsCent.push_back(jetIter->p4());

    } // loop over jets
  }

  // loop also over forward jets
  if ( L1JetsFwdHandle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsFwdHandle -> begin(); jetIter != L1JetsFwdHandle->end(); ++jetIter) {
      // combine with central jets
      jets.push_back(jetIter->p4());
    }
  }

  sort(jets.begin(),jets.end(),cmpPtLorentz());
  sort(jetsCent.begin(),jetsCent.end(),cmpPtLorentz());

  float jet1pt = 0.;

  if (jets.size() > 0) {
    jet1pt = jets.at(0).pt();
  }

  //  float jet1ptCent = 0.;
  float jet2ptCent = 0.;
  //  float jet3ptCent = 0.;
  float jet4ptCent = 0.;
  float dphiJet12Cent = -99.;

  // if (jetsCent.size() > 0) {
  //   jet1ptCent = jetsCent.at(0).pt();
  // }
  if (jetsCent.size() > 1) {
    jet2ptCent = jetsCent.at(1).pt();
    dphiJet12Cent = dphi(jetsCent.at(0).phi(),jetsCent.at(1).phi());
  }
  // if (jetsCent.size() > 2) {
  //   jet3ptCent = jetsCent.at(2).pt();
  // }
  if (jetsCent.size() > 3) {
    jet4ptCent = jetsCent.at(3).pt();
  }

  //
  // ----------------------------------------------------------------------
  // retrieve the caloMHT objects 
  //

  //  float mht = 0.;
  float ht = 0.;

  edm::Handle<L1EtMissParticleCollection> L1MHTHandle;
  iEvent.getByLabel(L1MHTInputTag, L1MHTHandle);
  
  if (L1MHTHandle.isValid() ) {
    for (L1EtMissParticleCollection::const_iterator mhtIter = L1MHTHandle->begin(); mhtIter != L1MHTHandle->end(); ++mhtIter) {
      //      mht = mhtIter -> etMiss();
      ht = mhtIter -> etTotal();
    }
  }
  
  //
  // ----------------------------------------------------------------------
  // retrieve the caloMET objects 2015
  //

  float met=0;

  edm::Handle<L1EtMissParticleCollection> L1EtMissHandle;
  iEvent.getByLabel(L1EtMissInputTag, L1EtMissHandle);

  if (L1EtMissHandle.isValid() ) {
    for (L1EtMissParticleCollection::const_iterator caloEtmIter = L1EtMissHandle->begin(); caloEtmIter != L1EtMissHandle->end(); ++caloEtmIter) {
      met=caloEtmIter -> et();
    }
  }

  *HTT95DoubleJetC30Dphi8bins = bool((ht >= 95.) && (jet2ptCent >= 30.) && (dphiJet12Cent < 2.6));
  *HTT125DoubleJetC30Dphi8bins = bool((ht >= 125.) && (jet2ptCent >= 30.) && (dphiJet12Cent < 2.6));
  *HTT175 = bool(ht >= 175.);
  *HTT200 = bool(ht >= 200.);
  *DoubleJetC100 = bool(jet2ptCent >= 100.);
  *DoubleJetC120 = bool(jet2ptCent >= 120.);
  *DoubleJetC60MET60 = bool((jet2ptCent >= 60.) && (met >= 60.));
  *QuadJetC60 = bool(jet4ptCent >= 60.);
  *SingleJet200 = bool(jet1pt >= 200.);

  // place branches into the event
  iEvent.put(HTT95DoubleJetC30Dphi8bins,aliasprefix_+"HTT95DoubleJetC30Dphi8bins");
  iEvent.put(HTT125DoubleJetC30Dphi8bins,aliasprefix_+"HTT125DoubleJetC30Dphi8bins");
  iEvent.put(HTT175,aliasprefix_+"HTT175");
  iEvent.put(HTT200,aliasprefix_+"HTT200");
  iEvent.put(DoubleJetC100,aliasprefix_+"DoubleJetC100");
  iEvent.put(DoubleJetC120,aliasprefix_+"DoubleJetC120");
  iEvent.put(DoubleJetC60MET60,aliasprefix_+"DoubleJetC60MET60");
  iEvent.put(QuadJetC60,aliasprefix_+"QuadJetC60");
  iEvent.put(SingleJet200,aliasprefix_+"SingleJet200");

}


//____________________________________________________________
float L1HadMenuDecisionProducer::dphi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if (dphi < -3.0) dphi += 2*TMath::Pi();
  else if (dphi > 3.2) dphi -= 2*TMath::Pi();
  //  std::cout << "dphi calc: phi1: " << phi1 << ", phi2: " << phi2 << ", diff: " << phi1-phi2 << ", dphi: " << dphi << std::endl;

  return dphi;
}

// ------------ method called once each job just before starting event loop  ------------
void 
L1HadMenuDecisionProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1HadMenuDecisionProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1HadMenuDecisionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1HadMenuDecisionProducer);
