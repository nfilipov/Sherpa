#ifndef SpinCorrelations_SpinCorrel_SpinCorrel_h
#define SpinCorrelations_SpinCorrel_SpinCorrel_h

//system include files
#include <memory>

// FWCore include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormat includes
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
//#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
//#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//ROOT includes...
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
//#include "TFileService.h"

//#include "Math.h"x
//using namespace reco;

//
// class declaration
//

class SpinCorrel : public edm::EDProducer {
public:
  explicit SpinCorrel(const edm::ParameterSet&);
  ~SpinCorrel();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
 private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void InitTree();

  edm::EDGetTokenT<reco::GenParticleCollection> theGenParticlesToken_;
  //	
  double mommass, mompt, momy,momphi; int momid=-99999, momst=-99999;
  double pmass, ppt, py,pphi;
  double daumass, daupt, daueta,dauphi; int dauid=-99999, daust=-99999;
  double gdaumass, gdaupt, gdaueta,gdauphi; int gdauid=-99999, gdaust=-99999;
  
  TFile * fOut;
  TTree * myTree;

  //  TClonesArr

  int LeptonID[1000];
  int MotherID[1000];
  int DaughterID[1000];
  int GrandDaughterID[1000];

  int LeptonStatus[1000];
  int MotherStatus[1000];
  int DaughterStatus[1000];
  int GrandDaughterStatus[1000];

  double LeptonPlusPt, LeptonPlusEta, LeptonPlusPhi, LeptonMinusMass;
  double LeptonMinusPt, LeptonMinusEta, LeptonMinusPhi, LeptonPlusMass;
  double TauPlusPiPt, TauPlusPiEta,TauPlusPiPhi, TauMinusPiMass;
  double TauMinusPiPt, TauMinusPiEta,TauMinusPiPhi, TauPlusPiMass;

  double PhotonEt, PhotonEta, PhotonPhi, PhotonMass; 
  double PhotonDEt, PhotonDEta, PhotonDPhi, PhotonDMass;

  std::vector<double> PhotonGDEt, PhotonGDEta, PhotonGDPhi, PhotonGDMass;

  const double PiMass=0.13957;
  const double TauMass = 1.777;
  

  
  int nl, nldp, nldm, id, st;

  bool filltree;
  /* int numberofkids; */
  /* int numbergdkids; */

  // edm::InputTag src_;
  // typedef std::vector<int> pdgid;
  // typedef std::vector<reco::GenParticle> recoGenParticles_genParticles__GEN_obj;
  // TBranch        *b_recoGenParticles_genParticles__GEN_obj;   //!
  
};


//
// constants, enums and typedefs
//


//
// static data member definitions
//

#endif
