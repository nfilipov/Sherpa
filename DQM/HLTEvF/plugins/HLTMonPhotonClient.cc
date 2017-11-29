#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQM/HLTEvF/interface/HLTMonPhotonClient.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

using namespace edm;

HLTMonPhotonClient::HLTMonPhotonClient(const edm::ParameterSet& iConfig)
{
  
  LogDebug("HLTMonPhotonClient") << "constructor...." ;
  
  logFile_.open("HLTMonPhotonClient.log");
  
  dbe = NULL;
  if (iConfig.getUntrackedParameter < bool > ("DQMStore", false)) {
    dbe = Service < DQMStore > ().operator->();
    dbe->setVerbose(0);
  }
  
  outputFile_ =
    iConfig.getUntrackedParameter <std::string>("outputFile", "");
  if (outputFile_.size() != 0) {
    LogInfo("HLTMonPhotonClient") << "Photon Trigger Monitoring histograms will be saved to " 
			      << outputFile_ ;
  }
  else {
    outputFile_ = "PhotonDQM.root";
  }
  
  bool disable =
    iConfig.getUntrackedParameter < bool > ("disableROOToutput", false);
  if (disable) {
    outputFile_ = "";
  }
  
  sourcetag_=iConfig.getParameter<edm::InputTag>("SourceTag");

  theHLTCollectionLabels = iConfig.getParameter<std::vector<edm::InputTag> >("theHLTCollectionLabels");
  
  dirname_="HLT/HLTMonPhoton/"+iConfig.getParameter<std::string>("@module_label");
  sourcedirname_="HLT/HLTMonPhoton/"+sourcetag_.label();
  
  if (dbe != NULL) {
    dbe->setCurrentFolder(dirname_);
  }
  
}


HLTMonPhotonClient::~HLTMonPhotonClient()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HLTMonPhotonClient::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  TH1F* refhisto = eventCounter->getTH1F();
  if(refhisto->GetBinContent(1) != 0){
    //    relFilterEff->setBinContent(6,refhisto->GetBinContent(5)/refhisto->GetBinContent(2));
    //    cumFilterEff->setBinContent(5,refhisto->GetBinContent(5)/refhisto->GetBinContent(1));
    cumFilterEff->setBinContent(4,refhisto->GetBinContent(4)/refhisto->GetBinContent(1));
    cumFilterEff->setBinContent(3,refhisto->GetBinContent(3)/refhisto->GetBinContent(1));
    cumFilterEff->setBinContent(2,refhisto->GetBinContent(2)/refhisto->GetBinContent(1));
  }else{
    relFilterEff->setBinContent(6,0);
    cumFilterEff->setBinContent(5,0);
    cumFilterEff->setBinContent(4,0);
    cumFilterEff->setBinContent(3,0);
    cumFilterEff->setBinContent(2,0);
  }
  for(int i = 0; i<4; i++){
    if(refhisto->GetBinContent(4-i) != 0){
      relFilterEff->setBinContent(5-i,refhisto->GetBinContent(5-i)/refhisto->GetBinContent(4-i));
    }else{
      relFilterEff->setBinContent(5-i,0);
    }
  }
  if(refhisto->GetBinContent(1) != 0){
    relFilterEff->setBinContent(1,refhisto->GetBinContent(1)/refhisto->GetBinContent(1));
    cumFilterEff->setBinContent(2,refhisto->GetBinContent(2)/refhisto->GetBinContent(1));
    cumFilterEff->setBinContent(1,refhisto->GetBinContent(1)/refhisto->GetBinContent(1));
  }else{
    relFilterEff->setBinContent(1,0);
    cumFilterEff->setBinContent(2,0);
    cumFilterEff->setBinContent(1,0);
  }
  
  TH1F* num;
  TH1F* denom;

  for(int i=0; i<3; i++){   
    denom = pixelhistosEt[0]->getTH1F();
    num = pixelhistosEt[i+1]->getTH1F();
    for(int j=1; j <= pixelhistosEtOut[i]->getNbinsX();j++ ){
      if(denom->GetBinContent(j)!=0){
	pixelhistosEtOut[i]->setBinContent(j,num->GetBinContent(j)/denom->GetBinContent(j));
      }
      else{
	pixelhistosEtOut[i]->setBinContent(j,0.);
      }
    }
    denom = pixelhistosEta[0]->getTH1F();
    num = pixelhistosEta[i+1]->getTH1F();  
    for(int j=1; j <= pixelhistosEtaOut[i]->getNbinsX();j++ ){
      if(denom->GetBinContent(j)!=0)
	pixelhistosEtaOut[i]->setBinContent(j,num->GetBinContent(j)/denom->GetBinContent(j));
      else
	pixelhistosEtaOut[i]->setBinContent(j,0.);
    }
    denom = pixelhistosPhi[0]->getTH1F();
    num = pixelhistosPhi[i+1]->getTH1F();  
    for(int j=1; j <= pixelhistosPhiOut[i]->getNbinsX();j++ ){
      if(denom->GetBinContent(j)!=0)
	pixelhistosPhiOut[i]->setBinContent(j,num->GetBinContent(j)/denom->GetBinContent(j));
      else
	pixelhistosPhiOut[i]->setBinContent(j,0.);
    }
  }

  denom = pixelhistosEt[0]->getTH1F();
  num = pixelhistosEt[3]->getTH1F();
  for(int j=1; j <= pixelhistosEtOut[3]->getNbinsX();j++ ){
    if(denom->GetBinContent(j)!=0){
      pixelhistosEtOut[3]->setBinContent(j,num->GetBinContent(j)/denom->GetBinContent(j));
    }
    else{
      pixelhistosEtOut[3]->setBinContent(j,0.);
    }
  }

  denom = pixelhistosEta[0]->getTH1F();
  num = pixelhistosEta[3]->getTH1F();
  for(int j=1; j <= pixelhistosEtaOut[3]->getNbinsX();j++ ){
    if(denom->GetBinContent(j)!=0){
      pixelhistosEtaOut[3]->setBinContent(j,num->GetBinContent(j)/denom->GetBinContent(j));
    }
    else{
      pixelhistosEtaOut[3]->setBinContent(j,0.);
    }
  }

  denom = pixelhistosPhi[0]->getTH1F();
  num = pixelhistosPhi[3]->getTH1F();
  for(int j=1; j <= pixelhistosPhiOut[3]->getNbinsX();j++ ){
    if(denom->GetBinContent(j)!=0){
      pixelhistosPhiOut[3]->setBinContent(j,num->GetBinContent(j)/denom->GetBinContent(j));
    }
    else{
      pixelhistosPhiOut[3]->setBinContent(j,0.);
    }
  }   
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
HLTMonPhotonClient::beginJob()
{
  
  DQMStore *dbe = 0;
  dbe = Service < DQMStore > ().operator->();
  
  if (dbe) {
    dbe->setCurrentFolder(dirname_);
    dbe->rmdir(dirname_);
  }
    
  if (dbe) {
    dbe->setCurrentFolder(dirname_);

    std::string tmpname = "";
      
    for(int i = 0; i<4; i++){
      LogInfo("HLTMonPhotonClient") << "loop iteration: "<<i ;
      tmpname = sourcedirname_ + "/" + theHLTCollectionLabels[i].label() + " Eta Dist";
      pixelhistosEta[i]=dbe->get(tmpname);
      tmpname = sourcedirname_ + "/" + theHLTCollectionLabels[i].label() + " Et Dist";
      pixelhistosEt[i]=dbe->get(tmpname);
      tmpname = sourcedirname_ + "/" + theHLTCollectionLabels[i].label() + " Phi Dist";
      pixelhistosPhi[i]=dbe->get(tmpname);
    }

    TH1F* refhist;
    int nBins;
    double xMin;
    double xMax;

    refhist = pixelhistosEt[0]->getTH1F();
    nBins = refhist->GetNbinsX();
    xMin = refhist->GetXaxis()->GetXmin();
    xMax = refhist->GetXaxis()->GetXmax();
    pixelhistosEtOut[0]  =dbe->book1D("EcalIsol Eff vs Et","EcalIsol Eff vs Et",nBins,xMin,xMax);
    pixelhistosEtOut[1]  =dbe->book1D("HcalIsol Eff vs Et","HcalIsol Eff vs Et",nBins,xMin,xMax);
    pixelhistosEtOut[2]  =dbe->book1D("Track Isol Eff vs Et","Track Isol Eff vs Et",nBins,xMin,xMax);
    pixelhistosEtOut[3]  =dbe->book1D("Total Eff vs Et","Total Eff vs Et",nBins,xMin,xMax);

    refhist = pixelhistosEta[0]->getTH1F();
    nBins = refhist->GetNbinsX();
    xMin = refhist->GetXaxis()->GetXmin();
    xMax = refhist->GetXaxis()->GetXmax();
    pixelhistosEtaOut[0] =dbe->book1D("EcalIsol Eff vs Eta","EcalIsol Eff vs Eta",nBins,xMin,xMax);
    pixelhistosEtaOut[1] =dbe->book1D("HcalIsol Eff vs Eta","HcalIsol Eff vs Eta",nBins,xMin,xMax);
    pixelhistosEtaOut[2] =dbe->book1D("Track Isol Eff vs Eta","Track Isol Eff vs Eta",nBins,xMin,xMax);
    pixelhistosEtaOut[3] =dbe->book1D("Total Eff vs Eta","Total Eff vs Eta",nBins,xMin,xMax);

    refhist = pixelhistosPhi[0]->getTH1F();
    nBins = refhist->GetNbinsX();
    xMin = refhist->GetXaxis()->GetXmin();
    xMax = refhist->GetXaxis()->GetXmax();
    pixelhistosPhiOut[0] =dbe->book1D("EcalIsol Eff vs Phi","EcalIsol Eff vs Phi",nBins,xMin,xMax);
    pixelhistosPhiOut[1] =dbe->book1D("HcalIsol Eff vs Phi","HcalIsol Eff vs Phi",nBins,xMin,xMax);
    pixelhistosPhiOut[2] =dbe->book1D("Track Isol Eff vs Phi","Track Isol Eff vs Phi",nBins,xMin,xMax);
    pixelhistosPhiOut[3] =dbe->book1D("Total Eff vs Phi","Total Eff vs Phi",nBins,xMin,xMax);

    tmpname = sourcedirname_ + "/Evts Passing Filters";    
    eventCounter = dbe->get(tmpname);

    relFilterEff = dbe->book1D("Relative Filter Effs","Relative Filter Effs",5,0,5);
    relFilterEff->setBinLabel(1,"HLTIsoEtFilterEff");
    relFilterEff->setBinLabel(2,"HLTIsoEcalIsolEff");
    relFilterEff->setBinLabel(3,"HLTIsoHcalIsolEff");
    relFilterEff->setBinLabel(4,"HLTIsoTrackIsolEff");
    //    relFilterEff->setBinLabel(5,"FIXME");
    //    relFilterEff->setBinLabel(6,"FIXME");

    cumFilterEff = dbe->book1D("Cumulative Filter Effs","Cumulative Filter Effs",5,0,5);
    cumFilterEff->setBinLabel(1,"HLTIsoEtFilterEff");
    cumFilterEff->setBinLabel(2,"HLTIsoEcalIsolEff");
    cumFilterEff->setBinLabel(3,"HLTIsoHcalIsolEff");
    cumFilterEff->setBinLabel(4,"HLTIsoTrackIsolEff");
    //    cumFilterEff->setBinLabel(5,"TrackIsolEff");

  } // end "if(dbe)"
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTMonPhotonClient::endJob() {

//     std::cout << "HLTMonPhotonClient: end job...." << std::endl;
 
   if (outputFile_.size() != 0 && dbe)
     dbe->save(outputFile_);
 
   return;
}

DEFINE_FWK_MODULE(HLTMonPhotonClient);
