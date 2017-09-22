// Root stuff
#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"

#include "TMath.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"

#include <TString.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>

#include <fstream>
#include <iostream>

using namespace std;

// miscellaneous  
#include "CMS_lumi.C"
#include "tdrstyle.C"

const TString basedir = "~/Sherpa/";

const TString Trees[4] = {"EEG","MMG","TTG","TTGSPIN"};


// colors
Color_t color1Spp = kCyan+1;
Color_t color2Spp = kCyan+2;
Color_t color3Spp = kCyan+3;
Color_t color1Saa = kGreen+2;
Color_t color2Saa = kGreen+3;
Color_t color3Saa = kGreen+4;
Color_t color1Sraa = kOrange+1;
Color_t color2Sraa = kOrange+4;
Color_t color3Sraa = kOrange+3;
 
const  ofstream ofSyst;

const float gTextSize = 0.04;

void SherpaKinematics()
{
  //setting the style
  setTDRStyle();

  TH1D* hPhotonPt[4]    ;
  TH1D* hPhotonEta[4]   ;
  TH1D* hPhotonPhi[4]   ;
  
  //lepton histograms (p
  TH1D* hLeptonPt[4]    ;
  TH1D* hLeptonEta[4]   ;
  TH1D* hLeptonPhi[4]   ;
			 
  // dilepton histograms
  TH1D* hDileptonPt[4]  ;
  TH1D* hDileptonMass[4];
  TH1D* hDileptonEta[4] ;

  TH1D* hTriobjectPt[4];
  TH1D* hTriobjectMass[4];
  TH1D* hTriobjectEta[4];
  TH1D* hTriobjectPhi[4];
  TH1D* hPiPiMass[2];
  for (int i = 0; i < 4; i++){
    
   TString inputFileName (basedir+Trees[i]+".root");
   TFile *f = TFile::Open(inputFileName,"READ");
   TTree *t = (TTree*) f->Get("myTree");

   TFile plots(basedir+Trees[i]+"_histos.root","RECREATE");
   
   TLorentzVector p,m,g[50],ll,llg[50],pip,pim,pipi;
   double _ppt,_peta,_pphi,_pm,
     _mpt,_meta,_mphi,_mm,
     _gpt,_geta,_gphi,_gm,
     _pippt,_pipeta,_pipphi,_pipm,
     _pimpt,_pimeta,_pimphi,_pimm;

   std::vector<double> *_photonpt=0;
   std::vector<double> *_photoneta=0;
   std::vector<double> *_photonm=0;
   std::vector<double> *_photonphi=0;
   
   //     _llpt,_lleta,_llphi,_llm,
   //  _llgpt,_llgeta,_llgphi,_llgm;

   int nldm, nldp; // number of lepton daughters
   t->SetBranchAddress("nldp",&nldp);
   t->SetBranchAddress("nldm",&nldm);
   
   t->SetBranchAddress("LeptonPlusPt",&_ppt);
   t->SetBranchAddress("LeptonPlusEta",&_peta);
   t->SetBranchAddress("LeptonPlusPhi",&_pphi);
   t->SetBranchAddress("LeptonPlusMass",&_pm);

   t->SetBranchAddress("LeptonMinusPt",&_mpt);
   t->SetBranchAddress("LeptonMinusEta",&_meta);
   t->SetBranchAddress("LeptonMinusPhi",&_mphi);
   t->SetBranchAddress("LeptonMinusMass",&_mm);

   t->SetBranchAddress("PhotonEt",&_photonpt);
   t->SetBranchAddress("PhotonEta",&_photoneta);
   t->SetBranchAddress("PhotonPhi",&_photonphi);
   t->SetBranchAddress("PhotonMass",&_photonm);

   t->SetBranchAddress("TauPlusPiPt",&_pippt);
   t->SetBranchAddress("TauPlusPiEta",&_pipeta);
   t->SetBranchAddress("TauPlusPiPhi",&_pipphi);
   t->SetBranchAddress("TauPlusPiMass",&_pipm);

   t->SetBranchAddress("TauMinusPiPt",&_pimpt);
   t->SetBranchAddress("TauMinusPiEta",&_pimeta);
   t->SetBranchAddress("TauMinusPiPhi",&_pimphi);
   t->SetBranchAddress("TauMinusPiMass",&_pimm);

   hPhotonPt[i]     = new TH1D("hPhotonPt","p_{T}^{#gamma} [GeV/c]",100,0,100);
   hPhotonEta[i]    = new TH1D("hPhotonEta","#eta^{#gamma}",100,-3.,3.);
   hPhotonPhi[i]    = new TH1D("hPhotonPhi","#phi^{#gamma}",100,-3.1415926536,-3.1415926536);
		 
   //lepton histograms (pt, eta, phi, by default done with lepton +)
   hLeptonPt[i]     = new TH1D("hLeptonPt","p_{T}^{lep} [GeV/c]",100,0,200);
   hLeptonEta[i]    = new TH1D("hLeptonEta","#eta^{lep}",100,-3.,3.);
   hLeptonPhi[i]    = new TH1D("hLeptonPhi","#phi^{lep}",100,-3.1415926536,-3.1415926536);
		 
   //dilepton histograms (pt, mass, rapidity)
   hDileptonPt[i]   = new TH1D("hDileptonPt","p_{T}^{ll} [GeV/c]",100,0,200);
   hDileptonMass[i] = new TH1D("hDileptonMass","m^{ll} [GeV/c^{2}]",100,50,150);
   hDileptonEta[i]  = new TH1D("hDileptonEta","y^{ll}",100,-3.,3.);

   // llgamma histograms
   hTriobjectPt[i]   = new TH1D("hTriobjectPt","p_{T}^{ll} [GeV/c]",100,0,200);
   hTriobjectMass[i] = new TH1D("hTriobjectMass","m^{ll} [GeV/c^{2}]",100,0,200);
   hTriobjectEta[i]  = new TH1D("hTriobjectEta","y^{ll#gamma}",100,-3.,3.);
   hTriobjectPhi[i]  = new TH1D("hTriobjectPhi","#phi^{ll#gamma}",100,-3.1415926536,3.1415926536);

   //pipi histogram
   if (i>1){
     hPiPiMass[i-2] = new TH1D("hPiPiMass","m^{#pi#pi} [GeV/c^{2}]",50,0,100);
   }
   //
   //let's-a-go
   int entries = t->GetEntries();
   for (int _i = 0 ; _i < entries ; _i++){
     t->GetEntry(_i);

     //     std::cout << "photon pt = {";
     //looping over photons
     /*     for (vector<double>::iterator it = _photonpt->begin(); it!=_photonpt->end(); ++it){
	    _gpt= *it;
	    _geta = _photoneta->at(*it);
	    hPhotonPt->Fill(_gpt);
	    }*/
     hLeptonPt[i]->Fill(_ppt);
     hLeptonEta[i]->Fill(_peta);
     hLeptonPhi[i]->Fill(_pphi);

     // building the lepton +, lepton -, and leptonlepton objects
     p.SetPtEtaPhiM(_ppt,_peta,_pphi,_pm);
     m.SetPtEtaPhiM(_mpt,_meta,_mphi,_mm);
     ll = p+m;

     hDileptonEta[i]->Fill(ll.Rapidity());
     hDileptonMass[i]->Fill(ll.M());
     hDileptonPt[i]->Fill(ll.Pt());
     // associating the photons to the leptonlepton pair
     std::vector<double>::const_iterator pt;
     std::vector<double>::const_iterator eta;
     std::vector<double>::const_iterator phi;
     std::vector<double>::const_iterator m;

     int order=0;
     
     for( pt = _photonpt->begin(), eta = _photoneta->begin(), phi = _photonphi->begin(), m = _photonm->begin();
	  pt < _photonpt->end() && eta < _photoneta->end() && phi < _photonphi->end() && m < _photonm->end();
	  ++pt, ++eta, ++phi, ++m )
       {
	 g[order].SetPtEtaPhiM(*pt,*eta,*phi,*m);
	 llg[order] = g[order] + ll;

	 hPhotonPt[i]->Fill(*pt);
	 hPhotonEta[i]->Fill(*eta);
	 hPhotonPhi[i]->Fill(*phi);

	 hTriobjectPt[i]->Fill(llg[order].Pt());
	 hTriobjectMass[i]->Fill(llg[order].M());
	 hTriobjectEta[i]->Fill(llg[order].Rapidity());
	 hTriobjectPhi[i]->Fill(llg[order].Phi());
	 
	 order++;
       } // filling objects loop: done

     if(i > 1 && (nldp == 2 && nldm == 2)) // when looking at the tau and tau+spin files, asking for tau(pinu)tau(pinu) decays only
       {
	 pip.SetPtEtaPhiM(_pippt,_pipeta,_pipphi,_pipm);
	 pim.SetPtEtaPhiM(_pimpt,_pimeta,_pimphi,_pimm);
	 pipi = pip+pim;

	 hPiPiMass[i-2]->Fill(pipi.M());
       }
   }
			   
   hPhotonPt[i]->Draw();
   hPhotonEta[i]->Draw();
   hPhotonPhi[i]->Draw();

   hLeptonPt[i]->Draw();
   hLeptonEta[i]->Draw();
   hLeptonPhi[i]->Draw(); 
   
   hDileptonPt[i]->Draw();
   hDileptonMass[i]->Draw();
   hDileptonEta[i]->Draw(); 

   hTriobjectPt[i]->Draw();
   hTriobjectMass[i]->Draw();
   hTriobjectEta[i]->Draw(); 
   hTriobjectPhi[i]->Draw();

   // hPhotonPt[i]->Draw();
   // hPhotonEta[i]->Draw();
   // hPhotonPhi[i]->Draw();

   // hLeptonPt[i]->Draw();
   // hLeptonEta[i]->Draw();
   // hLeptonPhi[i]->Draw(); 
   
   // hDileptonPt[i]->Draw();
   // hDileptonMass[i]->Draw();
   // hDileptonEta[i]->Draw(); 

   // hTriobjectPt[i]->Draw();
   // hTriobjectMass[i]->Draw();
   // hTriobjectEta[i]->Draw(); 
   // hTriobjectPhi[i]->Draw();

   if(i>1){
     hPiPiMass[i-2]->Draw();
   }
   
   //   f->Close();
   //   t->Reset();
   std::cout << "itt vagyok" << std::endl;
   plots.Write();
  }
  std::cout << "Ahhh" << std::endl;
  //  hPhotonPt->Draw();
}

