// -*- C++ -*-
//
// Package:    SpinCorrelations/SpinCorrel
// Class:      SpinCorrel
// 
/**\class SpinCorrel SpinCorrel.cc SpinCorrelations/SpinCorrel/plugins/SpinCorrel.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nicolas Filipovic
//         Created:  Tue, 05 Sep 2017 15:51:32 GMT
//
//

// ...
// most includes are in the header file.
#include "SpinCorrelations/SpinCorrel/interface/SpinCorrel.h"
// class declaration is in the header as well.

// construction
SpinCorrel::SpinCorrel(const edm::ParameterSet& iConfig) :
  theGenParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))){
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  //  src_ = iConfig.getParameter<edm::InputTag>("src");
  //  produces<pat::CompositeCandidateCollection>();
  ///  produces<pdgid>("")
}

// destruction
SpinCorrel::~SpinCorrel()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
SpinCorrel::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  //tiny bit of CMSSW voodoo
  Handle<GenParticleCollection> theGenParticles;
  iEvent.getByToken(theGenParticlesToken_, theGenParticles );

  if (theGenParticles.isValid()){
    for(size_t iGenParticle=0; iGenParticle<theGenParticles->size();++iGenParticle){

	// that's the container for gen-level particles in the event (it is saved at AOD level too, so this will work somewhere else)
	const GenParticle & p = (*theGenParticles)[iGenParticle];
	int id = p.pdgId();
	int st = p.status();

	// another data format for particles in the event, preferably for the decayed candidate. 
	LorentzVector myMom; // that' the Z candidate
	LorentzVector myLepP, myLepM; // e, mu, or tau 
	// TLorentzVector myTau1PCandidates
	
	size_t n =0, nd = 0;
	if ((id == 15 || id == -15) && (st == 11)){
	  const Candidate * mom = p.mother();
	  // std::cout<<"---------"<<std::endl;
	  // std::cout<<" Mother!"<<std::endl;
	  //	  reco::GenParticleRef mom1(theGenParticles,iGenParticle);
	  mommass = mom->mass(); mompt = mom->pt(); momy = mom->rapidity(); momphi = mom->phi();
	  momid = mom->pdgId(); momst = mom->status();

	  n = p.numberOfDaughters(); 
	  pmass = p.mass(); ppt = p.pt(); py = p.eta(); pphi = p.phi();

	  
	    // std::cout<<"---------"<<std::endl;  
	    // std::cout<<"[("<<momid<<","<<momst<<") pt="<< mompt <<"GeV, y="<<momy<<", phi="<<momphi<<", mass="<<mommass<<"];"<<std::endl;
	    // std::cout<<"Tau candidate:";	  
	    // std::cout<<" [("<<id<<","<<st<<") pt="<< ppt <<"GeV, y="<<py<<", phi="<<pphi<<", mass="<<pmass<<"];" << std::endl;
	    // std::cout<<"Here are "<<n<<" tau daughters and "<<nd<<" granddaughters:"<<std::endl;	
	  for (size_t j = 0 ; j<n ; ++j ){
	    const Candidate * d = p.daughter(j);
	    daupt = d->pt(); daueta = d->eta(); dauphi = d->phi(); daumass = d->mass(); 
	    dauid = d->pdgId(); daust = d->status(); 

	    if (daust==2){	
	      //	      std::cout <<"  [("<<dauid<<","<<daust<<") pt="<< daupt <<"GeV, eta="<<daueta<<", phi="<<dauphi<<", mass="<<daumass<<"];"<<std::endl;
	      nd = d->numberOfDaughters();
	      for (size_t jd = 0; jd < nd;++jd){
		const Candidate * gd = d->daughter(jd);
		gdaupt = gd->pt(); gdaueta = gd->eta(); gdauphi = gd->phi(); gdaumass = gd->mass(); 
		gdauid = gd->pdgId(); gdaust = gd->status(); 
		//std::cout <<"    [("<<gdauid<<","<<gdaust<<") pt="<< gdaupt <<"GeV, eta="<<gdaueta<<", phi="<<gdauphi<<", mass="<<gdaumass<<"];"<<std::endl;
		//	  	if(!(gdauid==16 || gdauid==16)
	      }
	    }else{continue;}
	  }
	}
	//if ( id == 15 || id == -15 ) {
	//	std::cout << " tau "<< st << " maman " << mom->pdgId();
	//}
	
      }
  }
   //   typedef Candidate::LorentzVector LorentzVector;
   
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
SpinCorrel::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SpinCorrel::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
SpinCorrel::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SpinCorrel::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SpinCorrel::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SpinCorrel::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SpinCorrel::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SpinCorrel);
