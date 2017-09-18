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

  typedef Candidate::LorentzVector LorentzVector;
  //set things to zero
  LeptonMinusPt = 0; 
  LeptonMinusEta = 0;
  LeptonMinusPhi = 0;
  LeptonMinusMass = 0;
  TauMinusPiPt = 0; 
  TauMinusPiEta = 0;
  TauMinusPiPhi = 0;
  TauMinusPiMass = 0;
  LeptonPlusPt = 0; 
  LeptonPlusEta = 0;
  LeptonPlusPhi = 0;
  LeptonMinusMass = 0;
  TauPlusPiPt = 0;
  TauPlusPiEta = 0;
  TauPlusPiPhi = 0;
  TauPlusPiMass = 0;
  PhotonEt = 0;
  PhotonEta = 0;
  PhotonPhi = 0;
  PhotonMass = 0;
  PhotonGDEt.clear();
  PhotonGDPhi.clear();
  PhotonGDEta.clear();
  PhotonGDMass.clear();

  // const Candidate * Z;
  // const Candidate * Zg;
  //tiny bit of CMSSW voodoo
  Handle<GenParticleCollection> theGenParticles;
  iEvent.getByToken(theGenParticlesToken_, theGenParticles );


  if (theGenParticles.isValid()){
    filltree=  false;
    for(size_t ig=0; ig<theGenParticles->size();++ig){
      
      // that's the container for gen-level particles in the event (it is saved at AOD level too, so this will work somewhere else)
      const GenParticle & p = (*theGenParticles)[ig];
      id = p.pdgId();
      st = p.status();
      
      // tau branch
      // we start with a tau st=3 but off mass shell, radiates whatever to make it on mass shell. st=11, m=1.777
      // then it fragments into a st=2 tau (like hadrons do, afaiu), pt eta phi may change, m=1.777
      // then it decays into whatever (25 pinu, 11 pipipinu, rest lepnu and other hadronic channels)
            
      // tau MINUS
      // if (id == 15){
      //   const Candidate * mom = p.mother();
      // 	if ( mom->status() == 3){
      // 	  filltree=true;
      // 	  LeptonMinusPt = mom->pt(); 
      // 	  LeptonMinusEta = mom->eta(); 
      // 	  LeptonMinusPhi = mom->phi();  /// mass is not taken here because this is off mass shell, wait for status 11
      // 	  // std::cout<<"---------"<<std::endl;  
      // 	  // std::cout<<"[("<<mom->pdgId()<<","<<mom->status()<<") pt="<< mom->pt() <<"GeV, eta="<<mom->rapidity()<<", phi="<<mom->phi()<<", mass="<<mom->mass<<"];"<<std::endl;
      // 	  // std::cout<<"Lepton candidate:";	  
      // 	  // std::cout<<" [("<<id<<","<<st<<") pt="<< LeptonMinusPt <<"GeV, y="<<LeptonMinusEta<<", phi="<<LeptonMinusPhi<<", mass="<<p.mass()<<"];" << std::endl;
      // 	  // std::cout<<"Here are the "<<p.numberOfDaughters()<<" daughters "<<std::endl;	
      // 	}
      // 	if (st == 11){
      // 	  if (id == 15 ) LeptonMinusMass = TauMass;
      // 	  if (id == 13 ) LeptonMinusMass = 0.105;
      // 	  if (id == 11 ) LeptonMinusMass = 0.000511;
      // 	  // mom->setMass(LeptonMinusMass);
      // 	  // lm = mom->p4();
      //  	  if (id == 15 && fabs(LeptonMinusMass - TauMass)>0.01) std::cout << "Lepton mass is not tau mass but "<<LeptonMinusMass<<std::endl;
      // 	}
      // 	  if (st==2){ //that's the decayed one, right ?
      // 	    nldm= p.numberOfDaughters(); // Number of Lepton Daughter Minus
      // 	    if (nldm == 2) {// that's the pinu decay
      // 	      for (int i=0;i<nldm;i++){
      // 		const Candidate * d = p.daughter(i);
      // 		if (d->pdgId()==-211){
      // 		  TauMinusPiPt = d->pt();
      // 		  TauMinusPiEta = d->eta();
      // 		  TauMinusPiPhi = d->phi();
      // 		  TauMinusPiMass = d->mass(); // should be 0.13957 always
      // 		  //		  std::cout << d->pdgId() << ", "<< TauMinusPiPt<<", "<< TauMinusPiEta<<", "<< TauMinusPiPhi<< ", "<< TauMinusPiMass<<std::endl;

      // 		}
      // 	      }	      
      // 	    }
      // 	  }
      // } 
      // //tau PLUS
      //       if (id == -15){
      // 	const Candidate * mom = p.mother();
      // 	if ( mom->status() == 3){
	  // 	  filltree=true;
	 /// mass is not taken here because this is off mass shell, wait for status 11
      // 	  // std::cout<<"---------"<<std::endl;  
      // 	  // std::cout<<"[("<<mom->pdgId()<<","<<mom->status()<<") pt="<< mom->pt() <<"GeV, eta="<<mom->rapidity()<<", phi="<<mom->phi()<<", mass="<<mom->mass<<"];"<<std::endl;
      // 	  // std::cout<<"Lepton candidate:";	  
      // 	  // std::cout<<" [("<<id<<","<<st<<") pt="<< LeptonMinusPt <<"GeV, y="<<LeptonMinusEta<<", phi="<<LeptonMinusPhi<<", mass="<<p.mass()<<"];" << std::endl;
      // 	  // std::cout<<"Here are the "<<p.numberOfDaughters()<<" daughters "<<std::endl;	
      // 	}
	  //	  if (st == 11){
	  //  }
	    // 	    if (id == -13 ) LeptonPlusMass = 0.105;
      // 	    if (id == -11 ) LeptonPlusMass = 0.000511;
      // 	    //	  lp = mom->p4();
      // 	    if (id == -15 && fabs(LeptonPlusMass - TauMass)>0.01) std::cout << "Lepton mass is not tau mass but "<<LeptonPlusMass<<std::endl;
      // 	}
      if (id == -15 && st==2){ //that's the decayed one, right ?7
	LeptonPlusPt = p.pt(); 	    LeptonPlusEta = p.eta();	    LeptonPlusPhi = p.phi(); 	    LeptonPlusMass = p.mass();
	nldp= p.numberOfDaughters(); // Number of Lepton Daughter Plus
	// std::cout <<"tau - status"<< st <<" "<< LeptonPlusPt << ","<<LeptonPlusEta<<","<<LeptonPlusPhi<<","<<LeptonPlusMass<<" # of daughters: "<< p.numberOfDaughters() << std::endl;
	if (nldp == 2) {// that's the pinu decay
	  for (int i=0;i<nldp;i++){
	    const Candidate * d = p.daughter(i);
	    if (d->pdgId()==211){
	      TauPlusPiPt = d->pt();
	      TauPlusPiEta = d->eta();
	      TauPlusPiPhi = d->phi();
	      TauPlusPiMass = d->mass(); // should be 0.13957 always
	      // std::cout <<"    (pi) - status"<< st <<" "<< TauPlusPiPt << ","<<TauPlusPiEta<<","<<TauPlusPiPhi<<","<<TauPlusPiMass<<" # of daughters: "<< d->numberOfDaughters()<< std::endl;
	      //		  std::cout << d->pdgId() << ", "<< TauPlusPiPt<<", "<< TauPlusPiEta<<", "<< TauPlusPiPhi<< ", "<< TauPlusPiMass<<std::endl;
	    }
	  }	      
	}
      }
      
      if (id == 15 && st==2){ //that's the decayed one, right ?7                                                                                                                                                                          
        LeptonMinusPt = p.pt();      LeptonMinusEta = p.eta();        LeptonMinusPhi = p.phi();        LeptonMinusMass = p.mass();
        nldm= p.numberOfDaughters(); // Number of Lepton Daughter Minus                                                                                                                                                                     
	// std::cout <<"tau - status"<< st <<" "<< LeptonMinusPt << ","<<LeptonMinusEta<<","<<LeptonMinusPhi<<","<<LeptonMinusMass<<" # of daughters: "<< p.numberOfDaughters() << std::endl;
        if (nldm == 2) {// that's the pinu decay                                                                                                                                                                                           
          for (int i=0;i<nldm;i++){
            const Candidate * d = p.daughter(i);
            if (d->pdgId()==-211 ){
              TauMinusPiPt = d->pt();
              TauMinusPiEta = d->eta();
              TauMinusPiPhi = d->phi();
              TauMinusPiMass = d->mass(); // should be 0.13957 always                                                                                                                                                                       
	      // std::cout <<"    (pi) - status"<< st <<" "<< TauMinusPiPt << ","<<TauMinusPiEta<<","<<TauMinusPiPhi<<","<<TauMinusPiMass<<" # of daughters: "<< d->numberOfDaughters()<< std::endl;
              //                  std::cout << d->pdgId() << ", "<< TauMinusPiPt<<", "<< TauMinusPiEta<<", "<< TauMinusPiPhi<< ", "<< TauMinusPiMass<<std::endl;                                                                               
            }
          }
        }
      }
      
      //      const Candidate ll;
      
      
      /// photon
      if (id == 22 && st==3 && p.pt()>1 ){ //start from the mother photon (there's always one, that's the sherpa config file here at play)

	PhotonEt = p.pt();
	PhotonEta = p.eta();
	PhotonPhi = p.phi();
	PhotonMass = p.mass();
	int nd =p.numberOfDaughters();
	/*	const Candidate * mom = p.mother();
	if  (!(mom->status()==11 && mom->pdgId()==22)){
	}*/

	// std::cout <<"Photon - status"<< st <<" "<< PhotonEt << ","<<PhotonEta<<","<<PhotonPhi<<","<<PhotonMass<<" # of daughters: "<< nd << std::endl;

	//next, loop over the daughters of this guy, and focus on the status 11 photon(s), if any
      	for (int i = 0 ; i<nd ; i++){
      	  const Candidate * dau = p.daughter(i);
      	  int daughter_id = dau->pdgId();
      	  int daughter_st = dau->status();

      	  PhotonDEt = dau->pt();
      	  PhotonDEta = dau->eta();
      	  PhotonDPhi = dau->phi();
      	  PhotonDMass = dau->mass();	

      	  if (daughter_id == 22 && daughter_st ==11 ) {
      	    // std::cout <<"   daughter - ("<< daughter_id<<","<<daughter_st<<") pt="<< PhotonDEt << ", eta="<<PhotonDEta<<", phi="<<PhotonDPhi<<", m="<<PhotonDMass<<" # of daughters: "<< dau->numberOfDaughters() << std::endl;
      	    int ngd = dau->numberOfDaughters();
      	  
	    // we're looking once at this generation (the granddaughters, and saving all photons with parentage to the hard event.

	    for(int i=0 ; i<ngd ; i++){
      	      const Candidate * gdau = dau->daughter(i);
      	      int gdaughter_id = gdau->pdgId();
	      int gdaughter_st = gdau->status();
	      if(gdaughter_st ==1){
		filltree=true;
		if (gdaughter_id == 22){
		  PhotonGDEt.push_back(gdau->pt());
		  PhotonGDEta.push_back(gdau->eta());
		  PhotonGDPhi.push_back(gdau->phi());
		  PhotonGDMass.push_back(gdau->mass());	
		  // std::cout <<"        grand daughter - ("<< gdaughter_id<<","<<gdaughter_st<<") pt="<< gdau->pt() << ", eta="<< gdau->eta()<<", phi="<< gdau->phi() <<", m="<< gdau->mass() <<" # of daughters: "<< gdau->numberOfDaughters() << std::endl;
		} else if (gdaughter_id == -11 || gdaughter_id == -13){
		  LeptonPlusPt = gdau->pt(); 
		  LeptonPlusEta = gdau->eta(); 
		  LeptonPlusPhi = gdau->phi();
		  LeptonPlusMass = gdau->mass(); 
		  // std::cout <<"        grand daughter - ("<< gdaughter_id<<","<<gdaughter_st<<") pt="<< gdau->pt() << ", eta="<< gdau->eta()<<", phi="<< gdau->phi() <<", m="<< gdau->mass() <<" # of daughters: "<< gdau->numberOfDaughters() << std::endl;
		} else if (gdaughter_id == 11 || gdaughter_id == 13){
                  LeptonMinusPt = gdau->pt();
                  LeptonMinusEta = gdau->eta();
                  LeptonMinusPhi = gdau->phi();
                  LeptonMinusMass = gdau->mass();
		  // std::cout <<"        grand daughter - ("<< gdaughter_id<<","<<gdaughter_st<<") pt="<< gdau->pt() << ", eta="<< gdau->eta()<<", phi="<< gdau->phi() <<", m="<< gdau->mass() <<" # of daughters: "<< gdau->numberOfDaughters() << std::endl;
                }
	      } else if (gdaughter_st == 2){
	      }
	      
      	      //std::cout <<"        grand daughter - ("<< gdaughter_id<<","<<gdaughter_st<<") pt="<< PhotonGDEt << ", eta="<<PhotonGDEta<<", phi="<<PhotonGDPhi<<", m="<<PhotonGDMass<<" # of daughters: "<< gdau->numberOfDaughters() << std::endl;
	      // we got here so this branch can be filled, and we push back the vectors. once this is done, we break from the id == 22 condition.
	    }
	    break; // breaking here should allow to save the good photons only ones.
      	  }
	}	
      }     



      // //mu PLUS
      // if (id == -13){
      // 	const Candidate * mom = p.mother();
      // 	if ( mom->status() == 3){ // take pt before FSR to get real lepton 4momentum (the one that's conserved in the decay..)
      // 	  filltree=true;
      // 	  LeptonPlusPt = mom->pt(); 
      // 	  LeptonPlusEta = mom->eta(); 
      // 	  LeptonPlusPhi = mom->phi();  /// mass is not taken here because this is off mass shell, wait for status 11
      // 	  // std::cout<<"---------"<<std::endl;  
      // 	  // std::cout<<"[("<<mom->pdgId()<<","<<mom->status()<<") pt="<< mom->pt() <<"GeV, eta="<<mom->rapidity()<<", phi="<<mom->phi()<<", mass="<<mom->mass<<"];"<<std::endl;
      // 	  // std::cout<<"Lepton candidate:";	  
      // 	  // std::cout<<" [("<<id<<","<<st<<") pt="<< LeptonMinusPt <<"GeV, y="<<LeptonMinusEta<<", phi="<<LeptonMinusPhi<<", mass="<<p.mass()<<"];" << std::endl;
      // 	}
      // 	if (st == 11){
      // 	  LeptonPlusMass = mom->mass();
      // 	}
      // } 

      // //mu MINUS
      //      if (id == 13){
      // 	const Candidate * mom = p.mother();
      // 	if ( mom->status() == 3){ // take pt before FSR to get real lepton 4momentum (the one that's conserved in the decay..)
      // 	  filltree=true;
      // 	  LeptonMinusPt = mom->pt(); 
      // 	  LeptonMinusEta = mom->eta(); 
      // 	  LeptonMinusPhi = mom->phi();  /// mass is not taken here because this is off mass shell, wait for status 11
      // 	  // std::cout<<"---------"<<std::endl;  
      // 	  // std::cout<<"[("<<mom->pdgId()<<","<<mom->status()<<") pt="<< mom->pt() <<"GeV, eta="<<mom->rapidity()<<", phi="<<mom->phi()<<", mass="<<mom->mass<<"];"<<std::endl;
      // 	  // std::cout<<"Lepton candidate:";	  
      // 	  // std::cout<<" [("<<id<<","<<st<<") pt="<< LeptonMinusPt <<"GeV, y="<<LeptonMinusEta<<", phi="<<LeptonMinusPhi<<", mass="<<p.mass()<<"];" << std::endl;
      // 	}
      // 	if (st == 11){
      // 	  LeptonMinusMass = mom->mass();
      // 	}
      // } 


      //e PLUS
      //     if (id == -11 && st == 1 && p.pt()>3){
	//	const Candidate * mom = p.mother();
      	//	if ( mom->status() == 3){ // take pt before FSR to get real lepton 4momentum (the one that's conserved in the decay..)
	//      	if(p.status()==11 && mom->status()==3){

      /*      if (id == -11 && st == 1 && p.pt()>3){
      	  filltree=true;
      	  LeptonPlusPt = p.pt(); 
      	  LeptonPlusEta = p.eta(); 
      	  LeptonPlusPhi = p.phi();
	  LeptonPlusMass = p.mass();/// mass is not taken here because this is off mass shell, wait for status 11
      	  // std::cout<<"----the mother:"<<std::endl;  
      	  // std::cout<<"[("<<mom->pdgId()<<","<<mom->status()<<") pt="<< mom->pt() <<"GeV, eta="<<mom->eta()<<", phi="<<mom->phi()<<", mass="<<mom->mass()<<"];"<<std::endl;
      	  // std::cout<<"Lepton candidate:";	  
	  //      	  std::cout<<" [("<<id<<","<<st<<"), mom ("<<mom->mother()->pdgId()<<","<<mom->mother()->status()<<"), pt="<< LeptonPlusPt <<"GeV, y="<<LeptonPlusEta<<", phi="<<LeptonPlusPhi<<", mass="<<LeptonPlusMass<<"];" << std::endl;
	  //      	}
	  } 
      */
      
      // if ((p.pdgId()==-11 || p.pdgId()==11) && p.status()<100 && p.pt()>1) 
      // 	{
      // 	  std::cout<<"[("<<p.pdgId()<<","<<p.status()<<") pt="<< p.pt() <<"GeV, eta="<<p.eta()<<", phi="<<p.phi()<<", mass="<<p.mass()<<"];"<<std::endl;
      // 	}
      // if ((p.pdgId()==-13 || p.pdgId()==13) && p.status()<100 && p.pt()>1) 
      // 	{
      // 	  std::cout<<"[("<<p.pdgId()<<","<<p.status()<<") pt="<< p.pt() <<"GeV, eta="<<p.eta()<<", phi="<<p.phi()<<", mass="<<p.mass()<<"];"<<std::endl;
      // 	}
      // if ((p.pdgId()==-15 || p.pdgId()==15) && p.status()<100 && p.pt()>1) 
      // 	{
      // 	  std::cout<<"[("<<p.pdgId()<<","<<p.status()<<") pt="<< p.pt() <<"GeV, eta="<<p.eta()<<", phi="<<p.phi()<<", mass="<<p.mass()<<"];"<<std::endl;
      // 	}

      //e MINUS
	//      if (id == 11 && st == 1 && p.pt()>3){
	//	const Candidate * mom = p.mother();
      // 	//	if ( mom->status() == 3){ // take pt before FSR to get real lepton 4momentum (the one that's conserved in the decay..)
      // 	if(p.status()==11){
	// 	  filltree=true;
	//	  LeptonMinusPt = p.pt(); 
      	//  LeptonMinusEta = p.eta(); 
      	//  LeptonMinusPhi = p.phi();  /// mass is not taken here because this is off mass shell, wait for status 11
	//  LeptonMinusMass = p.mass();      // 	  // std::cout<<"---------"<<std::endl;  
	  //      	  std::cout<<" [("<<id<<","<<st<<"), mom ("<<mom->mother()->pdgId()<<","<<mom->mother()->status()<<"), pt="<< LeptonMinusPt <<"GeV, y="<<LeptonMinusEta<<", phi="<<LeptonMinusPhi<<", mass="<<LeptonMinusMass<<"];" << std::endl;
      // 	  // std::cout<<"[("<<mom->pdgId()<<","<<mom->status()<<") pt="<< mom->pt() <<"GeV, eta="<<mom->rapidity()<<", phi="<<mom->phi()<<", mass="<<mom->mass()<<"];"<<std::endl;
      // 	  // std::cout<<"Lepton candidate:";	  
      // 	  // std::cout<<" [("<<id<<","<<st<<") pt="<< LeptonMinusPt <<"GeV, y="<<LeptonMinusEta<<", phi="<<LeptonMinusPhi<<", mass="<<p.mass()<<"];" << std::endl;
      // 	}
	// } 


     // size_t n =0, nd = 0;

     //  if ((id == 13 || id == -13) && (st == 11)){
     //  	LeptonID[ig] = id;
     //  	const Candidate * mom = p.mother();
     //  	// std::cout<<"---------"<<std::endl;
     //  	// std::cout<<" Mother!"<<std::endl;
     //  	//	  reco::GenParticleRef mom1(theGenParticles,ig);
     //  	mommass = mom->mass(); mompt = mom->pt(); momy = mom->rapidity(); momphi = mom->phi();
     //  	momid = mom->pdgId(); momst = mom->status();

     //  	n = p.numberOfDaughters(); 
     //  	pmass = p.mass(); ppt = p.pt(); py = p.eta(); pphi = p.phi();
     //  	if (filltree){
     //  	  std::cout<<"---------"<<std::endl;  
     //  	  std::cout<<"[("<<momid<<","<<momst<<") pt="<< mompt <<"GeV, y="<<momy<<", phi="<<momphi<<", mass="<<mommass<<"];"<<std::endl;
     //  	  std::cout<<"Tau candidate:";	  
     //  	  std::cout<<" [("<<id<<","<<st<<") pt="<< ppt <<"GeV, y="<<py<<", phi="<<pphi<<", mass="<<pmass<<"];" << std::endl;
     //  	  std::cout<<"Here are "<<n<<" tau daughters and "<<nd<<" granddaughters:"<<std::endl;	
     //  	}
     //  	for (size_t j = 0 ; j<n ; ++j ){
     //  	  const Candidate * d = p.daughter(j);
     //  	  daupt = d->pt(); daueta = d->eta(); dauphi = d->phi(); daumass = d->mass(); 
     //  	  dauid = d->pdgId(); daust = d->status(); 
	    	       
     //  	  if (daust==2){	
     //  	    if(filltree){
     //  	      std::cout <<"  [("<<dauid<<","<<daust<<") pt="<< daupt <<"GeV, eta="<<daueta<<", phi="<<dauphi<<", mass="<<daumass<<"];"<<std::endl;}
     //  	    nd = d->numberOfDaughters();
     //  	    for (size_t jd = 0; jd < nd;++jd){
     //  	      const Candidate * gd = d->daughter(jd);
     //  	      gdaupt = gd->pt(); gdaueta = gd->eta(); gdauphi = gd->phi(); gdaumass = gd->mass(); 
     //  	      gdauid = gd->pdgId(); gdaust = gd->status(); 
     //  	      if (filltree) {std::cout <<"    [("<<gdauid<<","<<gdaust<<") pt="<< gdaupt <<"GeV, eta="<<gdaueta<<", phi="<<gdauphi<<", mass="<<gdaumass<<"];"<<std::endl;}
     //  	    }
     //  	  } else{continue;}
     //  	}
     //  }
	
      }
    }
    if (filltree) myTree->Fill();
  }  

   
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
 


// ------------ method called once each job just before starting event loop  ------------
void 
SpinCorrel::beginJob()
{
  InitTree();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SpinCorrel::endJob() {
  fOut->cd();
  myTree->Write();
}

void SpinCorrel::InitTree(){
  // TFileService
  // edm::Service<TFileService> fs;


  //  myTree=fs->make<TTree>("myTree","My TTree of leptons");
   fOut = new TFile ("/eos/user/n/nfilipov/zg_sherpa/EEG.root","RECREATE");
  //  fOut = new TFile ("./ttg.root","RECREATE");

  myTree=new TTree("myTree","My TTree of leptons");
  myTree->Branch("LeptonID",&LeptonID, "LeptonID[1000]/I");

  myTree->Branch("LeptonMinusPt",&LeptonMinusPt, "LeptonMinusPt/D");
  myTree->Branch("LeptonMinusEta",&LeptonMinusEta, "LeptonMinusEta/D");
  myTree->Branch("LeptonMinusPhi",&LeptonMinusPhi, "LeptonMinusPhi/D");
  myTree->Branch("LeptonMinusMass",&LeptonMinusMass, "LeptonMinusMass/D");

  myTree->Branch("LeptonPlusPt",&LeptonPlusPt, "LeptonPlusPt/D");
  myTree->Branch("LeptonPlusEta",&LeptonPlusEta, "LeptonPlusEta/D");
  myTree->Branch("LeptonPlusPhi",&LeptonPlusPhi, "LeptonPlusPhi/D");
  myTree->Branch("LeptonPlusMass",&LeptonPlusMass, "LeptonPlusMass/D");
 
  myTree->Branch("nldp",&nldp,"nldp/I");
  myTree->Branch("nldm",&nldm,"nldm/I");

  myTree->Branch("PhotonEt",&PhotonGDEt );
  myTree->Branch("PhotonEta",&PhotonGDEta );
  myTree->Branch("PhotonPhi",&PhotonGDPhi );
  myTree->Branch("PhotonMass",&PhotonGDMass);

  myTree->Branch("TauPlusPiPt",&TauPlusPiPt, "TauPlusPiPt/D");
  myTree->Branch("TauPlusPiEta",&TauPlusPiEta, "TauPlusPiEta/D");
  myTree->Branch("TauPlusPiPhi",&TauPlusPiPhi, "TauPlusPiPhi/D");
  myTree->Branch("TauPlusPiMass",&TauPlusPiMass, "TauPlusPiMass/D");

  myTree->Branch("TauMinusPiPt",&TauMinusPiPt, "TauMinusPiPt/D");
  myTree->Branch("TauMinusPiEta",&TauMinusPiEta, "TauMinusPiEta/D");
  myTree->Branch("TauMinusPiPhi",&TauMinusPiPhi, "TauMinusPiPhi/D");
  myTree->Branch("TauMinusPiMass",&TauMinusPiMass, "TauMinusPiMass/D");

  //  myTree->Branch("MotherID",&MotherID, "MotherID[1000]/I");
  /// myTree->Branch("DaughterID",&MotherID, "DaughterID[1000]/I");
  
} 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SpinCorrel::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
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


//define this as a plug-in
DEFINE_FWK_MODULE(SpinCorrel);
