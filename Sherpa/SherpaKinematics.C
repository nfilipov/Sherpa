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
#include "TLegend.h"

#include <TString.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TTree.h>
#include <TH1F.h>

#include <fstream>
#include <iostream>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;

// miscellaneous
#include "CMS_lumi.C"
#include "tdrstyle.C"

//const TString basedir = "~/postdoc_ELTE/Zgamma/Sherpa/plaots/";  //linux
//const TString basedir = "~/Sherpa/Sherpa/";
const TString basedir = "./";

const TString Trees[4] = {"EEG_100k","MMG_100k","TTG","TTGSPIN"};


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

    //lepton histograms (p
  TH1D* hLeptonPlusPt[4]    ;
  TH1D* hLeptonPlusEta[4]   ;
  TH1D* hLeptonPlusPhi[4]   ;

      //lepton histograms (p
  TH1D* hLeptonMinusPt[4]    ;
  TH1D* hLeptonMinusEta[4]   ;
  TH1D* hLeptonMinusPhi[4]   ;

  // dilepton histograms
  TH1D* hDileptonPt[4]  ;
  TH1D* hDileptonMass[4];
  TH1D* hDileptonEta[4] ;

  TH1D* hTriobjectPt[4];
  TH1D* hTriobjectMass[4];
  TH1D* hTriobjectEta[4];
  TH1D* hTriobjectPhi[4];
  TH1D* hPiPiMass[2];

  TH1D* hDeltaR[4];
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
      _pimpt,_pimeta,_pimphi,_pimm,
      _dRpp, _dRpm; // delta R between photon and lepton(plus/minus)

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

    hPhotonPt[i]     = new TH1D("hPhotonPt",";p_{T}_{#gamma} [GeV/c]",50,10,110);
    hPhotonEta[i]    = new TH1D("hPhotonEta",";#eta_{#gamma}",50,-3.,3.);
    hPhotonPhi[i]    = new TH1D("hPhotonPhi",";#phi_{#gamma}",50,-3.1415926536,-3.1415926536);

    //lepton histograms (pt, eta, phi, by default done with lepton +)
    hLeptonPt[i]     = new TH1D("hLeptonPt",";p_{T}_{lep} [GeV/c]",50,0,200);
    hLeptonEta[i]    = new TH1D("hLeptonEta",";#eta_{lep}",50,-3.,3.);
    hLeptonPhi[i]    = new TH1D("hLeptonPhi",";#phi_{lep}",50,-3.1415926536,-3.1415926536);
    //lepton histograms (pt, eta, phi, by default done with lepton +)
    hLeptonPlusPt[i]     = new TH1D("hLeptonPlusPt",";p_{T}_{lep} [GeV/c]",50,0,200);
    hLeptonPlusEta[i]    = new TH1D("hLeptonPlusEta",";#eta_{lep}",50,-3.,3.);
    hLeptonPlusPhi[i]    = new TH1D("hLeptonPlusPhi",";#phi_{lep}",50,-3.1415926536,-3.1415926536);
    //lepton histograms (pt, eta, phi, by default done with lepton +)
    hLeptonMinusPt[i]     = new TH1D("hLeptonMinusPt",";p_{T}_{lep} [GeV/c]",50,0,200);
    hLeptonMinusEta[i]    = new TH1D("hLeptonMinusEta",";#eta_{lep}",50,-3.,3.);
    hLeptonMinusPhi[i]    = new TH1D("hLeptonMinusPhi",";#phi_{lep}",50,-3.1415926536,-3.1415926536);

    //dilepton histograms (pt, mass, rapidity)
    hDileptonPt[i]   = new TH1D("hDileptonPt",";p_{T}_{ll} [GeV/c]",50,0,200);
    hDileptonMass[i] = new TH1D("hDileptonMass",";m_{ll} [GeV/c^{2}]",50,30,180);
    hDileptonEta[i]  = new TH1D("hDileptonEta",";y_{ll}",50,-3.,3.);

    // llgamma histograms
    hTriobjectPt[i]   = new TH1D("hTriobjectPt",";p_{T}_{ll} [GeV/c]",50,0,200);
    hTriobjectMass[i] = new TH1D("hTriobjectMass",";m_{ll} [GeV/c^{2}]",50,0,200);
    hTriobjectEta[i]  = new TH1D("hTriobjectEta",";y_{ll#gamma}",50,-3.,3.);
    hTriobjectPhi[i]  = new TH1D("hTriobjectPhi",";#phi_{ll#gamma}",50,-3.1415926536,3.1415926536);

    hDeltaR[i] =  new TH1D("hDeltaR",";#DeltaR(l,#gamma)",80,0,8);
    
    //pipi histogram
    if (i>1){
      hPiPiMass[i-2] = new TH1D("hPiPiMass",";m_{#pi#pi} [GeV/c^{2}]",20,0,100);
    }
    //
    //let's-a-go
    int entries = t->GetEntries();
    for (int _i = 0 ; _i < entries ; _i++){
      t->GetEntry(_i);

      // if the lepton + and the lepton - have pt > 10 GeV
      if( true /*_ppt > 10 && _mpt > 10*/){
	
	// building the lepton +, lepton -, and leptonlepton objects

	p.SetPtEtaPhiM(_ppt,_peta,_pphi,_pm);
	m.SetPtEtaPhiM(_mpt,_meta,_mphi,_mm);
	
	ll = p+m;
	
	// associating the photons to the leptonlepton pair

	std::vector<double>::const_iterator pt;
	std::vector<double>::const_iterator eta;
	std::vector<double>::const_iterator phi;
	std::vector<double>::const_iterator m;

	int order=0;

	// looping over the contents of 4 vectors, to build the photon TLorentzVector
	
	for( pt = _photonpt->begin(), eta = _photoneta->begin(), phi = _photonphi->begin(), m = _photonm->begin();
	     pt < _photonpt->end() && eta < _photoneta->end() && phi < _photonphi->end() && m < _photonm->end();
	     ++pt, ++eta, ++phi, ++m )
	  {
	    
	    _dRpp = 0;
	    _dRpm = 0;
	    // Photon TLV
	    g[order].SetPtEtaPhiM(*pt,*eta,*phi,*m);
	    // LLG TLV
	    llg[order] = g[order] + ll;
	    
	    if( _ppt > 15 && _mpt > 15 && *pt > 10){
	    
	      _dRpp = sqrt(pow((_peta - *eta),2) + pow((_pphi - *phi ),2));
	      _dRpm = sqrt(pow((_meta - *eta),2) + pow((_mphi - *phi ),2));

	      if (min(_dRpp, _dRpm) > 0.05){
		hLeptonPt[i]->Fill(_ppt);
		hLeptonEta[i]->Fill(_peta);
		hLeptonPhi[i]->Fill(_pphi);

		hLeptonPlusPt[i]->Fill(_ppt);
		hLeptonPlusEta[i]->Fill(_peta);
		hLeptonPlusPhi[i]->Fill(_pphi);

		hLeptonMinusPt[i]->Fill(_mpt);
		hLeptonMinusEta[i]->Fill(_meta);
		hLeptonMinusPhi[i]->Fill(_mphi);

		hPhotonPt[i]->Fill(*pt);
		hPhotonEta[i]->Fill(*eta);
		hPhotonPhi[i]->Fill(*phi);
		
		hDileptonEta[i]->Fill(ll.Rapidity());
		hDileptonMass[i]->Fill(ll.M());
		hDileptonPt[i]->Fill(ll.Pt());
		
		hTriobjectPt[i]->Fill(llg[order].Pt());
		hTriobjectMass[i]->Fill(llg[order].M());
		hTriobjectEta[i]->Fill(llg[order].Rapidity());
		hTriobjectPhi[i]->Fill(llg[order].Phi());

		hDeltaR[i]->Fill(min(_dRpp,_dRpm));
	      }
	    }
	    order++;
	  } // filling objects loop: done

	if(i > 1 && (nldp == 2 && nldm == 2)) // when looking at the tau and tau+spin files, asking for tau(pinu)tau(pinu) decays only
	  {
	    pip.SetPtEtaPhiM(_pippt,_pipeta,_pipphi,_pipm);
	    pim.SetPtEtaPhiM(_pimpt,_pimeta,_pimphi,_pimm);
	    pipi = pip+pim;
	    TLorentzVector tautau;
	    tautau = ll;
	    TVector3 tautauBoost = tautau.BoostVector();
	    TLorentzVector pipi_comframe = pipi;
	    pipi_comframe.Boost(-tautauBoost);
	    //	    std::cout<<"tautau 4 momentum = ("<<tautau_comframe.E()<<", - "<<tautau_comframe.Px()<<", - "<<tautau_comframe.Py()<<", - "<<tautau_comframe.Pz()<<")"<<std::endl;
	    
	    if (pipi.M()>0) hPiPiMass[i-2]->Fill(pipi_comframe.M());
	  }
      }
    }

   
    hLeptonPlusEta[i]->SetFillStyle(3001);
    hLeptonMinusEta[i]->SetFillStyle(3001);
    hLeptonPlusEta[i]->SetFillColor(kBlue);
    hLeptonMinusEta[i]->SetFillColor(kRed);
    hLeptonPlusEta[i]->SetMarkerSize(0); 
    hLeptonMinusEta[i]->SetMarkerSize(0);

    hLeptonPlusPt[i]->SetFillStyle(3001);
    hLeptonMinusPt[i]->SetFillStyle(3001);
    hLeptonPlusPt[i]->SetFillColor(kBlue);
    hLeptonMinusPt[i]->SetFillColor(kRed);
    hLeptonPlusPt[i]->SetMarkerSize(0);
    hLeptonMinusPt[i]->SetMarkerSize(0);

    hLeptonPlusPhi[i]->SetFillStyle(3001);
    hLeptonMinusPhi[i]->SetFillStyle(3001);
    hLeptonPlusPhi[i]->SetFillColor(kBlue);
    hLeptonMinusPhi[i]->SetFillColor(kRed);
    hLeptonPlusPhi[i]->SetMarkerSize(0);
    hLeptonMinusPhi[i]->SetMarkerSize(0);
    
    hPhotonPt[i]->Draw();
    hPhotonEta[i]->Draw();
    hPhotonPhi[i]->Draw();

    hLeptonPt[i]->Draw();
    hLeptonEta[i]->Draw();
    hLeptonPhi[i]->Draw();

    hLeptonPlusPt[i]->Draw();
    hLeptonPlusEta[i]->Draw();
    hLeptonPlusPhi[i]->Draw();
    hLeptonMinusPt[i]->Draw();
    hLeptonMinusEta[i]->Draw();
    hLeptonMinusPhi[i]->Draw();

    hDileptonPt[i]->Draw();
    hDileptonMass[i]->Draw();
    hDileptonEta[i]->Draw();

    hTriobjectPt[i]->Draw();
    hTriobjectMass[i]->Draw();
    hTriobjectEta[i]->Draw();
    hTriobjectPhi[i]->Draw();

    hDeltaR[i]->Draw();

    if(i>1){
      hPiPiMass[i-2]->Draw();
    }

    //   f->Close();
    //   t->Reset();
    std::cout << "itt vagyok" << std::endl;
    plots.Write();

    TCanvas cLepPMEta; cLepPMEta.cd();
    hLeptonPlusEta[i]->Draw("e2"); hLeptonMinusEta[i]->Draw("e2same");
    //   cLepPMEta.SaveAs("cLepPM_eta_"+Trees[i]+".C");
    cLepPMEta.Close();

    TCanvas cLepPMPt; cLepPMPt.cd();
    hLeptonPlusPt[i]->Draw("e2"); hLeptonMinusPt[i]->Draw("e2same");
    //    cLepPMPt.SaveAs("cLepPM_Pt_"+Trees[i]+".C");
    cLepPMPt.Close();

    TCanvas cLepPMPhi; cLepPMPhi.cd();
    hLeptonPlusPhi[i]->Draw("e2"); hLeptonMinusPhi[i]->Draw("e2same");
    //    cLepPMPhi.SaveAs("cLepPM_Phi_"+Trees[i]+".C");
    cLepPMPhi.Close();
      
  }

  std::cout << "Ahhh" << std::endl;
  //  hPhotonPt->Draw();
}

void drawPlots(){
  setTDRStyle();
  TFile *eppt = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * heppt; heppt = (TH1D*) eppt->Get("hPhotonPt");

  TFile *muppt = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmuppt; hmuppt = (TH1D*) muppt->Get("hPhotonPt");

  TFile *tauppt = TFile::Open("TTG_histos.root","READ");
  TH1D *htauppt; htauppt =(TH1D*)  tauppt->Get("hPhotonPt");

  heppt->SetFillColor(kAzure+1);
  hmuppt->SetFillColor(kRed+1);
  htauppt->SetFillColor(kGreen+1);

  /*  heppt->SetFillStyle(3001);
  hmuppt->SetFillStyle(3001);
  htauppt->SetFillStyle(3001);
  */
  heppt->SetMarkerSize(0);
  hmuppt->SetMarkerSize(0);
  htauppt->SetMarkerSize(0);

  heppt->GetYaxis()->SetRangeUser(10,5000);
  TCanvas cppt;// cppt.Divide(1,2); 
  cppt.cd();

  TPad *ppt = new TPad("ppt","ppt",0,0.3,1.,1.0);
  ppt->SetBottomMargin(0);
  ppt->Draw();
  ppt->SetLogy();
  ppt->cd();


  hmuppt->Draw("e2");
  htauppt->Draw("e2same");
  heppt->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * plpt = new TLegend(0.45,0.4,0.88,0.7,"");
  plpt->SetTextSize(0.032);
  plpt->SetFillStyle(0);
  plpt->SetFillColor(0);
  plpt->SetBorderSize(0);
  plpt->SetTextFont(42);
  plpt->AddEntry(heppt  ,"pp #rightarrow ee#gamma","f");
  plpt->AddEntry(hmuppt ,"pp #rightarrow #mu#mu#gamma","f");
  plpt->AddEntry(htauppt,"pp #rightarrow #tau#tau#gamma","f");
  plpt->Draw();
  */
  cppt.cd();

  TPad *pptr = new TPad("pptr","pptr",0,0.05,1.,0.3);  pptr->SetBottomMargin(0.2);   pptr->Draw();   pptr->cd();

  TH1D * hpptr1; hpptr1 = (TH1D*)hmuppt->Clone("hpptr1");
  hpptr1->Sumw2();
  hpptr1->SetMinimum(0), hpptr1->SetMaximum(2);
  hpptr1->Divide(heppt);

  TH1D * hpptr2; hpptr2 = (TH1D*)htauppt->Clone("hpptr2");
  hpptr2->Sumw2();
  hpptr2->SetMinimum(0), hpptr2->SetMaximum(2);
  hpptr2->Divide(heppt);

  hpptr2->SetFillStyle(3001);
  hpptr1->Draw("e2");
  hpptr2->Draw("e2same");

  //  cppt.SaveAs("cppt.C");
  
  /// lepton pt

  TFile *elpt = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * helpt; helpt = (TH1D*) elpt->Get("hLeptonPt");

  TFile *mulpt = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmulpt; hmulpt = (TH1D*) mulpt->Get("hLeptonPt");

  TFile *taulpt = TFile::Open("TTG_histos.root","READ");
  TH1D * htaulpt; htaulpt =(TH1D*)  taulpt->Get("hLeptonPt");

  helpt->SetFillColor(kAzure+1);
  hmulpt->SetFillColor(kRed+1);
  htaulpt->SetFillColor(kGreen+1);

  /*  helpt->SetFillStyle(3001);
  hmulpt->SetFillStyle(3001);
  htaulpt->SetFillStyle(3001);
  */
  helpt->SetMarkerSize(0);
  hmulpt->SetMarkerSize(0);
  htaulpt->SetMarkerSize(0);

  helpt->GetYaxis()->SetRangeUser(10,5000);
  TCanvas clpt;// clpt.Divide(1,2); 
  clpt.cd();

  TPad *lpt = new TPad("lpt","lpt",0,0.3,1.,1.0);
  lpt->SetBottomMargin(0);
  lpt->Draw();
  lpt->SetLogy();
  lpt->cd();

  hmulpt->Draw("e2");
  htaulpt->Draw("e2same");
  helpt->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * llpt = new TLegend(0.45,0.4,0.88,0.7,"");
  llpt->SetTextSize(0.032);
  llpt->SetFillStyle(0);
  llpt->SetFillColor(0);
  llpt->SetBorderSize(0);
  llpt->SetTextFont(42);
  llpt->AddEntry(helpt  ,"pp #rightarrow ee#gamma","f");
  llpt->AddEntry(hmulpt ,"pp #rightarrow #mu#mu#gamma","f");
  llpt->AddEntry(htaulpt,"pp #rightarrow #tau#tau#gamma","f");
  llpt->Draw();
  */
  clpt.cd();

  TPad *lptr = new TPad("lptr","lptr",0,0.05,1.,0.3);  lptr->SetBottomMargin(0.2);   lptr->Draw();   lptr->cd();

  TH1D * hlptr1; hlptr1 = (TH1D*)hmulpt->Clone("hlptr1");
  hlptr1->Sumw2();
  hlptr1->SetMinimum(0), hlptr1->SetMaximum(2);
  hlptr1->Divide(helpt);

  TH1D * hlptr2; hlptr2 = (TH1D*)htaulpt->Clone("hlptr2");
  hlptr2->Sumw2();
  hlptr2->SetMinimum(0), hlptr2->SetMaximum(2);
  hlptr2->Divide(helpt);

  hlptr2->SetFillStyle(3001);
  hlptr1->Draw("e2");
  hlptr2->Draw("e2same");

  //  clpt.SaveAs("clpt.C");

  ///// eta

  TFile *epeta = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * hepeta; hepeta = (TH1D*) epeta->Get("hPhotonEta");

  TFile *mupeta = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmupeta; hmupeta = (TH1D*) mupeta->Get("hPhotonEta");

  TFile *taupeta = TFile::Open("TTG_histos.root","READ");
  TH1D *htaupeta; htaupeta =(TH1D*)  taupeta->Get("hPhotonEta");

  hepeta->SetFillColor(kAzure+1);
  hmupeta->SetFillColor(kRed+1);
  htaupeta->SetFillColor(kGreen+1);

  /*  hepeta->SetFillStyle(3001);
  hmupeta->SetFillStyle(3001);
  htaupeta->SetFillStyle(3001);
  */
  hepeta->SetMarkerSize(0);
  hmupeta->SetMarkerSize(0);
  htaupeta->SetMarkerSize(0);

  hepeta->GetYaxis()->SetRangeUser(10,5000);
  TCanvas cpeta;// cpeta.Divide(1,2); 
  cpeta.cd();

  TPad *peta = new TPad("peta","peta",0,0.3,1.,1.0);
  peta->SetBottomMargin(0);
  peta->Draw();
  peta->SetLogy();
  peta->cd();


  hmupeta->Draw("e2");
  htaupeta->Draw("e2same");
  hepeta->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * pleta = new TLegend(0.45,0.4,0.88,0.7,"");
  pleta->SetTextSize(0.032);
  pleta->SetFillStyle(0);
  pleta->SetFillColor(0);
  pleta->SetBorderSize(0);
  pleta->SetTextFont(42);
  pleta->AddEntry(hepeta  ,"pp #rightarrow ee#gamma","f");
  pleta->AddEntry(hmupeta ,"pp #rightarrow #mu#mu#gamma","f");
  pleta->AddEntry(htaupeta,"pp #rightarrow #tau#tau#gamma","f");
  pleta->Draw();
  */
  cpeta.cd();

  TPad *petar = new TPad("petar","petar",0,0.05,1.,0.3);  petar->SetBottomMargin(0.2);   petar->Draw();   petar->cd();

  TH1D * hpetar1; hpetar1 = (TH1D*)hmupeta->Clone("hpetar1");
  hpetar1->Sumw2();
  hpetar1->SetMinimum(0), hpetar1->SetMaximum(2);
  hpetar1->Divide(hepeta);

  TH1D * hpetar2; hpetar2 = (TH1D*)htaupeta->Clone("hpetar2");
  hpetar2->Sumw2();
  hpetar2->SetMinimum(0), hpetar2->SetMaximum(2);
  hpetar2->Divide(hepeta);

  hpetar2->SetFillStyle(3001);
  hpetar1->Draw("e2");
  hpetar2->Draw("e2same");

  //  cpeta.SaveAs("cpeta.C");
  
/// lepton Eta

  TFile *eleta = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * heleta; heleta = (TH1D*) eleta->Get("hLeptonEta");

  TFile *muleta = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmuleta; hmuleta = (TH1D*) muleta->Get("hLeptonEta");

  TFile *tauleta = TFile::Open("TTG_histos.root","READ");
  TH1D * htauleta; htauleta =(TH1D*)  tauleta->Get("hLeptonEta");

  heleta->SetFillColor(kAzure+1);
  hmuleta->SetFillColor(kRed+1);
  htauleta->SetFillColor(kGreen+1);

  // heleta->SetFillStyle(3001);
  // hmuleta->SetFillStyle(3001);
  // htauleta->SetFillStyle(3001);

  heleta->SetMarkerSize(0);
  hmuleta->SetMarkerSize(0);
  htauleta->SetMarkerSize(0);

  heleta->GetYaxis()->SetRangeUser(10,5000);
  TCanvas cleta;// cleta.Divide(1,2); 
  cleta.cd();

  TPad *leta = new TPad("leta","leta",0,0.3,1.,1.0);
  leta->SetBottomMargin(0);
  leta->Draw();
  leta->SetLogy();
  leta->cd();

  hmuleta->Draw("e2");
  htauleta->Draw("e2same");
  heleta->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * lleta = new TLegend(0.45,0.4,0.88,0.7,"");
  lleta->SetTextSize(0.032);
  lleta->SetFillStyle(0);
  lleta->SetFillColor(0);
  lleta->SetBorderSize(0);
  lleta->SetTextFont(42);
  lleta->AddEntry(heleta  ,"pp #rightarrow ee#gamma","f");
  lleta->AddEntry(hmuleta ,"pp #rightarrow #mu#mu#gamma","f");
  lleta->AddEntry(htauleta,"pp #rightarrow #tau#tau#gamma","f");
  lleta->Draw();
  */
  cleta.cd();

  TPad *letar = new TPad("letar","letar",0,0.05,1.,0.3);  letar->SetBottomMargin(0.2);   letar->Draw();   letar->cd();

  TH1D * hletar1; hletar1 = (TH1D*)hmuleta->Clone("hletar1");
  hletar1->Sumw2();
  hletar1->SetMinimum(0), hletar1->SetMaximum(2);
  hletar1->Divide(heleta);

  TH1D * hletar2; hletar2 = (TH1D*)htauleta->Clone("hletar2");
  hletar2->Sumw2();
  hletar2->SetMinimum(0), hletar2->SetMaximum(2);
  hletar2->Divide(heleta);

  hletar2->SetFillStyle(3001);
  hletar1->Draw("e2");
  hletar2->Draw("e2same");

  //  cleta.SaveAs("cleta.C");

  
  //// mass dilepton

  TFile *edmass = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * hedmass; hedmass = (TH1D*) edmass->Get("hDileptonMass");

  TFile *mudmass = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmudmass; hmudmass = (TH1D*) mudmass->Get("hDileptonMass");

  TFile *taudmass = TFile::Open("TTG_histos.root","READ");
  TH1D *htaudmass; htaudmass =(TH1D*)  taudmass->Get("hDileptonMass");

  hedmass->SetFillColor(kAzure+1);
  hmudmass->SetFillColor(kRed+1);
  htaudmass->SetFillColor(kGreen+1);

  // hedmass->SetFillStyle(3001);
  // hmudmass->SetFillStyle(3001);
  // htaudmass->SetFillStyle(3001);
  
  hedmass->SetMarkerSize(0);
  hmudmass->SetMarkerSize(0);
  htaudmass->SetMarkerSize(0);

  hedmass->GetYaxis()->SetRangeUser(10,5000);
  TCanvas cdmass;// cdmass.Divide(1,2); 
  cdmass.cd();

  TPad *dmass = new TPad("dmass","dmass",0,0.3,1.,1.0);
  dmass->SetBottomMargin(0);
  dmass->Draw();
  dmass->SetLogy();
  dmass->cd();


  hmudmass->Draw("e2");
  htaudmass->Draw("e2same");
  hedmass->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * ldmass = new TLegend(0.45,0.4,0.88,0.7,"");
  ldmass->SetTextSize(0.032);
  ldmass->SetFillStyle(0);
  ldmass->SetFillColor(0);
  ldmass->SetBorderSize(0);
  ldmass->SetTextFont(42);
  ldmass->AddEntry(hedmass  ,"pp #rightarrow ee#gamma","f");
  ldmass->AddEntry(hmudmass ,"pp #rightarrow #mu#mu#gamma","f");
  ldmass->AddEntry(htaudmass,"pp #rightarrow #tau#tau#gamma","f");
  ldmass->Draw();
  */
  cdmass.cd();

  TPad *dmassr = new TPad("dmassr","dmassr",0,0.05,1.,0.3);  dmassr->SetBottomMargin(0.2);   dmassr->Draw();   dmassr->cd();

  TH1D * hdmassr1; hdmassr1 = (TH1D*)hmudmass->Clone("hdmassr1");
  hdmassr1->Sumw2();
  hdmassr1->SetMinimum(0), hdmassr1->SetMaximum(2);
  hdmassr1->Divide(hedmass);

  TH1D * hdmassr2; hdmassr2 = (TH1D*)htaudmass->Clone("hdmassr2");
  hdmassr2->Sumw2();
  hdmassr2->SetMinimum(0), hdmassr2->SetMaximum(2);
  hdmassr2->Divide(hedmass);

  hdmassr2->SetFillStyle(3001);
  hdmassr1->Draw("e2");
  hdmassr2->Draw("e2same");

  //  cdmass.SaveAs("cdmass.C");
  
  ///// mass triobject

  TFile *etmass = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * hetmass; hetmass = (TH1D*) etmass->Get("hTriobjectMass");

  TFile *mutmass = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmutmass; hmutmass = (TH1D*) mutmass->Get("hTriobjectMass");

  TFile *tautmass = TFile::Open("TTG_histos.root","READ");
  TH1D *htautmass; htautmass =(TH1D*)  tautmass->Get("hTriobjectMass");

  hetmass->SetFillColor(kAzure+1);
  hmutmass->SetFillColor(kRed+1);
  htautmass->SetFillColor(kGreen+1);

  // hetmass->SetFillStyle(3001);
  // hmutmass->SetFillStyle(3001);
  // htautmass->SetFillStyle(3001);
 
  hetmass->SetMarkerSize(0);
  hmutmass->SetMarkerSize(0);
  htautmass->SetMarkerSize(0);

  hetmass->GetYaxis()->SetRangeUser(10,5000);
  TCanvas ctmass;// ctmass.Divide(1,2); 
  ctmass.cd();

  TPad *tmass = new TPad("tmass","tmass",0,0.3,1.,1.0);
  tmass->SetBottomMargin(0);
  tmass->Draw();
  tmass->SetLogy();
  tmass->cd();


  hmutmass->Draw("e2");
  htautmass->Draw("e2same");
  hetmass->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * tlmass = new TLegend(0.45,0.4,0.88,0.7,"");
  tlmass->SetTextSize(0.032);
  tlmass->SetFillStyle(0);
  tlmass->SetFillColor(0);
  tlmass->SetBorderSize(0);
  tlmass->SetTextFont(42);
  tlmass->AddEntry(hetmass  ,"pp #rightarrow ee#gamma","f");
  tlmass->AddEntry(hmutmass ,"pp #rightarrow #mu#mu#gamma","f");
  tlmass->AddEntry(htautmass,"pp #rightarrow #tau#tau#gamma","f");
  tlmass->Draw();
  */
  ctmass.cd();

  TPad *tmassr = new TPad("tmassr","tmassr",0,0.05,1.,0.3);  tmassr->SetBottomMargin(0.2);   tmassr->Draw();   tmassr->cd();

  TH1D * htmassr1; htmassr1 = (TH1D*)hmutmass->Clone("htmassr1");
  htmassr1->Sumw2();
  htmassr1->SetMinimum(0), htmassr1->SetMaximum(2);
  htmassr1->Divide(hetmass);

  TH1D * htmassr2; htmassr2 = (TH1D*)htautmass->Clone("htmassr2");
  htmassr2->Sumw2();
  htmassr2->SetMinimum(0), htmassr2->SetMaximum(2);
  htmassr2->Divide(hetmass);

  htmassr2->SetFillStyle(3001);
  htmassr1->Draw("e2");
  htmassr2->Draw("e2same");

  //  ctmass.SaveAs("ctmass.C");
 
  ///// pt dilepton

    TFile *edpt = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * hedpt; hedpt = (TH1D*) edpt->Get("hDileptonPt");

  TFile *mudpt = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmudpt; hmudpt = (TH1D*) mudpt->Get("hDileptonPt");

  TFile *taudpt = TFile::Open("TTG_histos.root","READ");
  TH1D *htaudpt; htaudpt =(TH1D*)  taudpt->Get("hDileptonPt");

  hedpt->SetFillColor(kAzure+1);
  hmudpt->SetFillColor(kRed+1);
  htaudpt->SetFillColor(kGreen+1);

  // hedpt->SetFillStyle(3001);
  // hmudpt->SetFillStyle(3001);
  // htaudpt->SetFillStyle(3001);
  
  hedpt->SetMarkerSize(0);
  hmudpt->SetMarkerSize(0);
  htaudpt->SetMarkerSize(0);

  hedpt->GetYaxis()->SetRangeUser(10,5000);
  TCanvas cdpt;// cdpt.Divide(1,2); 
  cdpt.cd();

  TPad *dpt = new TPad("dpt","dpt",0,0.3,1.,1.0);
  dpt->SetBottomMargin(0);
  dpt->Draw();
  dpt->SetLogy();
  dpt->cd();


  hmudpt->Draw("e2");
  htaudpt->Draw("e2same");
  hedpt->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * dlpt = new TLegend(0.45,0.4,0.88,0.7,"");
  dlpt->SetTextSize(0.032);
  dlpt->SetFillStyle(0);
  dlpt->SetFillColor(0);
  dlpt->SetBorderSize(0);
  dlpt->SetTextFont(42);
  dlpt->AddEntry(hedpt  ,"pp #rightarrow ee#gamma","f");
  dlpt->AddEntry(hmudpt ,"pp #rightarrow #mu#mu#gamma","f");
  dlpt->AddEntry(htaudpt,"pp #rightarrow #tau#tau#gamma","f");
  dlpt->Draw();
  */
  cdpt.cd();

  TPad *dptr = new TPad("dptr","dptr",0,0.05,1.,0.3);  dptr->SetBottomMargin(0.2);   dptr->Draw();   dptr->cd();

  TH1D * hdptr1; hdptr1 = (TH1D*)hmudpt->Clone("hdptr1");
  hdptr1->Sumw2();
  hdptr1->SetMinimum(0), hdptr1->SetMaximum(2);
  hdptr1->Divide(hedpt);

  TH1D * hdptr2; hdptr2 = (TH1D*)htaudpt->Clone("hdptr2");
  hdptr2->Sumw2();
  hdptr2->SetMinimum(0), hdptr2->SetMaximum(2);
  hdptr2->Divide(hedpt);

  hdptr2->SetFillStyle(3001);
  hdptr1->Draw("e2");
  hdptr2->Draw("e2same");

  //  cdpt.SaveAs("cdpt.C");

  //// pttriobject
  
  TFile *etpt = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * hetpt; hetpt = (TH1D*) etpt->Get("hTriobjectPt");

  TFile *mutpt = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmutpt; hmutpt = (TH1D*) mutpt->Get("hTriobjectPt");

  TFile *tautpt = TFile::Open("TTG_histos.root","READ");
  TH1D *htautpt; htautpt =(TH1D*)  tautpt->Get("hTriobjectPt");

  hetpt->SetFillColor(kAzure+1);
  hmutpt->SetFillColor(kRed+1);
  htautpt->SetFillColor(kGreen+1);

  //  hetpt->SetFillStyle(3001);
  // hmutpt->SetFillStyle(3001);
  // htautpt->SetFillStyle(3001);
  
  hetpt->SetMarkerSize(0);
  hmutpt->SetMarkerSize(0);
  htautpt->SetMarkerSize(0);

  hetpt->GetYaxis()->SetRangeUser(10,5000);
  TCanvas ctpt;// ctpt.Divide(1,2); 
  ctpt.cd();

  TPad *tpt = new TPad("tpt","tpt",0,0.3,1.,1.0);
  tpt->SetBottomMargin(0);
  tpt->Draw();
  tpt->SetLogy();
  tpt->cd();


  hmutpt->Draw("e2");
  htautpt->Draw("e2same");
  hetpt->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * tlpt = new TLegend(0.45,0.4,0.88,0.7,"");
  tlpt->SetTextSize(0.032);
  tlpt->SetFillStyle(0);
  tlpt->SetFillColor(0);
  tlpt->SetBorderSize(0);
  tlpt->SetTextFont(42);
  tlpt->AddEntry(hetpt  ,"pp #rightarrow ee#gamma","f");
  tlpt->AddEntry(hmutpt ,"pp #rightarrow #mu#mu#gamma","f");
  tlpt->AddEntry(htautpt,"pp #rightarrow #tau#tau#gamma","f");
  tlpt->Draw();
  */
  ctpt.cd();

  TPad *tptr = new TPad("tptr","tptr",0,0.05,1.,0.3);  tptr->SetBottomMargin(0.2);   tptr->Draw();   tptr->cd();

  TH1D * htptr1; htptr1 = (TH1D*)hmutpt->Clone("htptr1");
  htptr1->Sumw2();
  htptr1->SetMinimum(0), htptr1->SetMaximum(2);
  htptr1->Divide(hetpt);

  TH1D * htptr2; htptr2 = (TH1D*)htautpt->Clone("htptr2");
  htptr2->Sumw2();
  htptr2->SetMinimum(0), htptr2->SetMaximum(2);
  htptr2->Divide(hetpt);

  htptr2->SetFillStyle(3001);
  htptr1->Draw("e2");
  htptr2->Draw("e2same");

  //  ctpt.SaveAs("ctpt.C");
  //// delta R

  TFile *edr = TFile::Open("EEG_100k_histos.root","READ");
  TH1D * hedr; hedr = (TH1D*) edr->Get("hDeltaR");

  TFile *mudr = TFile::Open("MMG_100k_histos.root","READ");
  TH1D * hmudr; hmudr = (TH1D*) mudr->Get("hDeltaR");

  TFile *taudr = TFile::Open("TTG_histos.root","READ");
  TH1D *htaudr; htaudr =(TH1D*)  taudr->Get("hDeltaR");

  hedr->SetFillColor(kAzure+1);
  hmudr->SetFillColor(kRed+1);
  htaudr->SetFillColor(kGreen+1);

  //   hedr->SetFillStyle(3001);
  // hmudr->SetFillStyle(3001);
  // htaudr->SetFillStyle(3001);
  
  hedr->SetMarkerSize(0);
  hmudr->SetMarkerSize(0);
  htaudr->SetMarkerSize(0);

  hedr->GetYaxis()->SetRangeUser(10,5000);
  TCanvas cdr;// cdr.Divide(1,2); 
  cdr.cd();

  TPad *dr = new TPad("dr","dr",0,0.3,1.,1.0);
  dr->SetBottomMargin(0);
  dr->Draw();
  dr->SetLogy();
  dr->cd();
  


  hmudr->Draw("e2");
  htaudr->Draw("e2same");
  hedr->Draw("e2same");  

  cout << "here" << endl;

  /*  TLegend * pldr = new TLegend(0.45,0.4,0.88,0.7,"");
  pldr->SetX1NDC(0.45),pldr->SetY1NDC(0.4),pldr->SetX2NDC(0.88),pldr->SetY2NDC(0.7);
  pldr->SetTextSize(0.035);
  pldr->SetFillStyle(0);
  //  pldr->SetFillColor(0);
  pldr->SetBorderSize(0);
  pldr->SetTextFont(42);
  pldr->AddEntry(hedr  ,"pp #rightarrow ee#gamma","F");
  pldr->AddEntry(hmudr ,"pp #rightarrow #mu#mu#gamma","F");
  pldr->AddEntry(htaudr,"pp #rightarrow #tau#tau#gamma","F");
  pldr->Draw();
  */
  cdr.cd();

  TPad *drr = new TPad("drr","drr",0,0.05,1.,0.3);  drr->SetBottomMargin(0.2);   drr->Draw();   drr->cd();

  TH1D * hdrr1; hdrr1 = (TH1D*)hmudr->Clone("hdrr1");
  hdrr1->Sumw2();
  hdrr1->SetMinimum(0), hdrr1->SetMaximum(2);
  hdrr1->Divide(hedr);

  TH1D * hdrr2; hdrr2 = (TH1D*)htaudr->Clone("hdrr2");
  hdrr2->Sumw2();
  hdrr2->SetMinimum(0), hdrr2->SetMaximum(2);
  hdrr2->Divide(hedr);

  hdrr1->GetXaxis()->SetTitleSize(0.2);
  hdrr2->SetFillStyle(3001);
  hdrr1->Draw("e2");
  hdrr2->Draw("e2same");

  //  cdr.SaveAs("cdr.C");

  //// end: hPiPiMass

  TFile *tausppm = TFile::Open("TTGSPIN_histos.root","READ");
  TH1D *htausppm; htausppm =(TH1D*)  tausppm->Get("hPiPiMass");

  TFile *tauppm = TFile::Open("TTG_histos.root","READ");
  TH1D *htauppm; htauppm =(TH1D*)  tauppm->Get("hPiPiMass");
  
  
  htausppm->SetFillColorAlpha(kAzure+1,0.7);
  htauppm->SetFillColorAlpha(kRed,0.7);
  // htausppm->SetFillStyle(3001);
  // htauppm->SetFillStyle(3001);

  htauppm->SetMarkerSize(0);
  htausppm->SetMarkerSize(0);

  //  tpt->GetYaxis()->SetRangeUser(10,5000);
  TCanvas ctaus;// ctpt.Divide(1,2); 
  ctaus.cd();
  
  TPad *ttaus = new TPad("ttaus","ttaus",0,0.3,1.,1.0);
  ttaus->SetBottomMargin(0);
  ttaus->Draw();
  ttaus->SetLogy();
  ttaus->cd();

  htauppm->Draw("e2");
  htausppm->Draw("e2same");
  cout << "here" << endl;

  /*  TLegend * tltaus = new TLegend(0.45,0.4,0.88,0.7,"");
  tltaus->SetTextSize(0.032);
  tltaus->SetFillStyle(0);
  tltaus->SetFillColor(0);
  tltaus->SetBorderSize(0);
  tltaus->SetTextFont(42);
  tltaus->AddEntry(htausppm ,"spin correlations disabled","f");
  tltaus->AddEntry(htauppm ,"spin correlations enabled","f");
  */
  ctaus.cd();
 
  TPad *ttausr = new TPad("ttausr","ttausr",0,0.05,1.,0.3);  ttausr->SetBottomMargin(0.2);   ttausr->Draw();   ttausr->cd();

  TH1D * httausr1; httausr1 = (TH1D*)htausppm->Clone("httausr1");
  httausr1->Sumw2();
  httausr1->SetMinimum(0), httausr1->SetMaximum(2);
  httausr1->Divide(htauppm);

  httausr1->Draw("e2");

  //ctaus.SaveAs("ctaus.C");
  
}

void plotMG(){
  setTDRStyle();
  TFile *mg = TFile::Open("ZG_255_1M.root","READ");
  TTreeReader myReader("LHEF",mg);
  
  TFile plots("ZG_madgraph_histos.root","RECREATE");

  TH1D* hPhoton_wePt = new TH1D("hPhoton_wePt",";p_{T}_{#gamma} [GeV/c]",50,10,110);
  TH1D* hPhoton_weEta = new TH1D("hPhoton_weEta",";#eta_{#gamma}",50,-3.,3.);
  TH1D* hPhoton_wePhi = new TH1D("hPhoton_wePhi",";#phi_{#gamma}",50,-3.1415926536,-3.1415926536);
  TH1D* hPhoton_wmPt = new TH1D("hPhoton_wmPt",";p_{T}_{#gamma} [GeV/c]",50,10,110);
  TH1D* hPhoton_wmEta = new TH1D("hPhoton_wmEta",";#eta_{#gamma}",50,-3.,3.);
  TH1D* hPhoton_wmPhi = new TH1D("hPhoton_wmPhi",";#phi_{#gamma}",50,-3.1415926536,-3.1415926536);
  TH1D* hPhoton_wtPt = new TH1D("hPhoton_wtPt",";p_{T}_{#gamma} [GeV/c]",50,10,110);
  TH1D* hPhoton_wtEta = new TH1D("hPhoton_wtEta",";#eta_{#gamma}",50,-3.,3.);
  TH1D* hPhoton_wtPhi = new TH1D("hPhoton_wtPhi",";#phi_{#gamma}",50,-3.1415926536,-3.1415926536);

  TH1D* hTriobject_wePt  = new TH1D("hTriobject_wePt",";p_{T}_{ee#gamma} [GeV/c]",50,0,200);
  TH1D* hTriobject_weMass = new TH1D("hTriobject_weMass",";m_{ee#gamma} [GeV/c^{2}]",50,0,200);
  TH1D* hTriobject_weEta = new TH1D("hTriobject_weEta",";y_{ee#gamma}",50,-3.,3.);
  TH1D* hTriobject_wePhi= new TH1D("hTriobject_wePhi",";#phi_{ee#gamma}",50,-3.1415926536,3.1415926536);
  TH1D* hTriobject_wmPt  = new TH1D("hTriobject_wmPt",";p_{T}_{ee#gamma} [GeV/c]",50,0,200);
  TH1D* hTriobject_wmMass = new TH1D("hTriobject_wmMass",";m_{ee#gamma} [GeV/c^{2}]",50,0,200);
  TH1D* hTriobject_wmEta = new TH1D("hTriobject_wmEta",";y_{ee#gamma}",50,-3.,3.);
  TH1D* hTriobject_wmPhi= new TH1D("hTriobject_wmPhi",";#phi_{ee#gamma}",50,-3.1415926536,3.1415926536);
  TH1D* hTriobject_wtPt  = new TH1D("hTriobject_wtPt",";p_{T}_{ee#gamma} [GeV/c]",50,0,200);
  TH1D* hTriobject_wtMass = new TH1D("hTriobject_wtMass",";m_{ee#gamma} [GeV/c^{2}]",50,0,200);
  TH1D* hTriobject_wtEta = new TH1D("hTriobject_wtEta",";y_{ee#gamma}",50,-3.,3.);
  TH1D* hTriobject_wtPhi= new TH1D("hTriobject_wtPhi",";#phi_{ee#gamma}",50,-3.1415926536,3.1415926536);

  TH1D* hPiPiMass;

  TH1D* hDelta_weR= new TH1D("hDelta_weR",";#DeltaR(e,#gamma)",80,0,8);
  TH1D* hDelta_wmR= new TH1D("hDelta_wmR",";#DeltaR(#mu,#gamma)",80,0,8);
  TH1D* hDelta_wtR= new TH1D("hDelta_wtR",";#DeltaR(#tau,#gamma)",80,0,8);

  TH1D* hElectronPlusPt= new TH1D("hElectronPlusPt",";p_{T}_{e+} [GeV/c]",50,0,200);
  TH1D* hElectronMinusPt= new TH1D("hElectronMinusPt",";p_{T}_{e-} [GeV/c]",50,0,200);
  TH1D* hElectronPlusEta= new TH1D("hElectronPlusEta",";#eta_{e+}",50,-3.,3.);
  TH1D* hElectronMinusEta= new TH1D("hElectronMinusEta",";#eta_{e-}",50,-3.,3.);
  TH1D* hElectronPlusPhi = new TH1D("hElectronPlusPhi",";#phi_{e+}",50,-3.1415926536,-3.1415926536);
  TH1D* hElectronMinusPhi = new TH1D("hElectronMinusPhi",";#phi_{e-}",50,-3.1415926536,-3.1415926536);
  
  TH1D* hMuonPlusPt= new TH1D("hMuonPlusPt",";p_{T}_{#mu+} [GeV/c]",50,0,200);
  TH1D* hMuonMinusPt= new TH1D("hMuonMinusPt",";p_{T}_{#mu-} [GeV/c]",50,0,200);
  TH1D* hMuonPlusEta= new TH1D("hMuonPlusEta",";#eta_{#mu+}",50,-3.,3.);
  TH1D* hMuonMinusEta= new TH1D("hMuonMinusEta",";#eta_{#mu-}",50,-3.,3.);
  TH1D* hMuonPlusPhi = new TH1D("hMuonPlusPhi",";#phi_{#mu+}",50,-3.1415926536,-3.1415926536);
  TH1D* hMuonMinusPhi = new TH1D("hMuonMinusPhi",";#phi_{#mu-}",50,-3.1415926536,-3.1415926536);
  
  TH1D* hTauPlusPt= new TH1D("hTauPlusPt",";p_{T}_{#tau+} [GeV/c]",50,0,200);
  TH1D* hTauMinusPt= new TH1D("hTauMinusPt",";p_{T}_{#tau-} [GeV/c]",50,0,200);
  TH1D* hTauPlusEta= new TH1D("hTauPlusEta",";#eta_{#tau+}",50,-3.,3.);
  TH1D* hTauMinusEta= new TH1D("hTauMinusEta",";#eta_{#tau-}",50,-3.,3.);
  TH1D* hTauPlusPhi = new TH1D("hTauPlusPhi",";#phi_{#tau+}",50,-3.1415926536,-3.1415926536);
  TH1D* hTauMinusPhi = new TH1D("hTauMinusPhi",";#phi_{#tau-}",50,-3.1415926536,-3.1415926536);

  TH1D* hEEPt = new TH1D("hEEPt",";p_{T}_{ee} [GeV/c]",50,0,200);
  TH1D* hEEMass = new TH1D("hEEMass",";m_{ee} [GeV/c^{2}]",50,0,200);
  TH1D* hEERapidity  = new TH1D("hEEEta",";y_{ee}",50,-3.,3.);
  TH1D* hMMPt = new TH1D("hMMPt",";p_{T}_{#mu#mu} [GeV/c]",50,0,200);
  TH1D* hMMMass = new TH1D("hMMMass",";m_{#mu#mu} [GeV/c^{2}]",50,0,200);
  TH1D* hMMRapidity = new TH1D("hMMEta",";y_{#mu#mu}",50,-3.,3.);
  TH1D* hTTPt = new TH1D("hTTPt",";p_{T}_{#tau#tau} [GeV/c]",50,0,200);
  TH1D* hTTMass = new TH1D("hTTMass",";m_{#tau#tau} [GeV/c^{2}]",50,0,200);
  TH1D* hTTRapidity = new TH1D("hTTEta",";y_{#tau#tau}",50,-3.,3.);
  TH1D* hTTPhi = new TH1D("hTTPhi",";#phi_{#tau#tau}",50,-3.1415926536,-3.1415926536);

  //  TTreeReaderValue<Long64_t> rvEventNumber(myReader, "Event.Number");
  TTreeReaderValue<Int_t> rvParticle_size(myReader, "Particle_size");
  TTreeReaderArray<Int_t> raParticlePID(myReader, "Particle.PID");
  TTreeReaderArray<Int_t> raParticleStatus(myReader, "Particle.Status");
  TTreeReaderArray<Int_t> raParticleMom1(myReader, "Particle.Mother1");
  TTreeReaderArray<Int_t> raParticleMom2(myReader, "Particle.Mother2");
  TTreeReaderArray<Double_t> raParticlePx(myReader, "Particle.Px");
  TTreeReaderArray<Double_t> raParticlePy(myReader, "Particle.Py");
  TTreeReaderArray<Double_t> raParticlePz(myReader, "Particle.Pz");
  TTreeReaderArray<Double_t> raParticleE(myReader, "Particle.E");
  TTreeReaderArray<Double_t> raParticleM(myReader, "Particle.M");
  TTreeReaderArray<Double_t> raParticlePT(myReader, "Particle.PT");
  TTreeReaderArray<Double_t> raParticleEta(myReader, "Particle.Eta");
  TTreeReaderArray<Double_t> raParticlePhi(myReader, "Particle.Phi");
  TTreeReaderArray<Double_t> raParticleRapidity(myReader, "Particle.Rapidity");
  TTreeReaderArray<Double_t> raParticleLifeTime(myReader, "Particle.LifeTime");
  TTreeReaderArray<Double_t> raParticleSpin(myReader, "Particle.Spin");
  TLorentzVector muon;	TLorentzVector antimuon;
  TLorentzVector electron;	TLorentzVector antielectron;
  TLorentzVector tau;	TLorentzVector antitau;
  TLorentzVector photon; TLorentzVector Z;
  TLorentzVector dilepton; TLorentzVector triobject;
  TLorentzVector Zgamma;

  double deltaR_p, deltaR_m;

  double photonpt, photoneta, photonphi, electronpt, antielectronpt,electroneta, antielectroneta,electronphi, antielectronphi, muonpt, antimuonpt,muoneta, antimuoneta,muonphi, antimuonphi, taupt, antitaupt,taueta, antitaueta,tauphi, antitauphi ;
    // Loop over all entries of the TTree or TChain.
  while (myReader.Next()) {
    bool event_e = false;
    bool event_m = false;
    bool event_t = false;
    photonpt =0; photoneta=0; photonphi=0; electronpt=0; antielectronpt=0;electroneta=0; antielectroneta=0; electronphi=0; antielectronphi=0; muonpt=0; antimuonpt=0;muoneta=0; antimuoneta=0;muonphi=0; antimuonphi=0; taupt=0; antitaupt =0; taueta=0; antitaueta=0; tauphi=0; antitauphi=0;
    Int_t particlesize = *rvParticle_size;
    for (int ip = 0; ip < particlesize ; ip ++){
      if (raParticlePID[ip] == 11){
	electron.SetPxPyPzE(raParticlePx[ip],raParticlePy[ip],raParticlePz[ip],raParticleE[ip]);
	electronpt = (double) electron.Pt();
	electroneta = (double) electron.Eta();
	electronphi = (double) electron.Phi();
      }
      if (raParticlePID[ip] == -11){
      	antielectron.SetPxPyPzE(raParticlePx[ip],raParticlePy[ip],raParticlePz[ip],raParticleE[ip]);
	antielectronpt = (double) antielectron.Pt();
	antielectroneta = (double) antielectron.Eta();
	antielectronphi = (double) antielectron.Phi();
      }
      if (raParticlePID[ip] == 13){
      	muon.SetPxPyPzE(raParticlePx[ip],raParticlePy[ip],raParticlePz[ip],raParticleE[ip]);
	muonpt = (double) muon.Pt();
	muoneta = (double) muon.Eta();
	muonphi = (double) muon.Phi();
      }
      if (raParticlePID[ip] == -13){
	antimuon.SetPxPyPzE(raParticlePx[ip],raParticlePy[ip],raParticlePz[ip],raParticleE[ip]);
	antimuonpt = (double) antimuon.Pt();
	antimuoneta = (double) antimuon.Eta();
	antimuonphi = (double) antimuon.Phi();
      }
      if (raParticlePID[ip] == 15){
      	tau.SetPxPyPzE(raParticlePx[ip],raParticlePy[ip],raParticlePz[ip],raParticleE[ip]);
	taupt = (double) tau.Pt();
	taueta = (double) tau.Eta();
	tauphi = (double) tau.Phi();
      }
      if (raParticlePID[ip] == -15){
      	antitau.SetPxPyPzE(raParticlePx[ip],raParticlePy[ip],raParticlePz[ip],raParticleE[ip]);
	antitaupt = (double) antitau.Pt();
	antitaueta = (double) antitau.Eta();
	antitauphi = (double) antitau.Phi();
      }
      if (raParticlePID[ip] == 22){
	photon.SetPxPyPzE(raParticlePx[ip],raParticlePy[ip],raParticlePz[ip],raParticleE[ip]);
	photonpt = photon.Pt();
	photoneta = photon.Eta();
	photonphi = photon.Phi();
      }
      //  std::cout<<std::endl;
    }/// end particlesize loop    
    // if ((event_m && event_e) || (event_e && event_t) || (event_t && event_m)){
    //   std::cout <<"at least two flavors in this event!"<< std::endl;
    // }
    if (electronpt > 15 && antielectronpt >15 && photonpt >10) {
      Z = electron + antielectron;
      deltaR_m = sqrt(
      		      pow((photon.Eta() - electron.Eta()),2) +
      		      pow((photon.Phi() - electron.Phi()),2)
      		      );
      deltaR_p = sqrt(
      		      pow((photon.Eta() - antielectron.Eta()),2) +
      		      pow((photon.Phi() - antielectron.Phi()),2)
      		      );
      if(min(deltaR_m,deltaR_p) > 0.05){
	hElectronMinusPt->Fill(electronpt);
	hElectronMinusEta->Fill(electroneta);
	hElectronMinusPhi->Fill(electronphi);
	hElectronPlusPt->Fill(antielectronpt);
	hElectronPlusEta->Fill(antielectroneta);
	hElectronPlusPhi->Fill(antielectronphi);
	Zgamma = Z + photon;
	hPhoton_wePt->Fill((double)photon.Pt());
	hPhoton_weEta->Fill((double)photon.Eta());
	hPhoton_wePhi->Fill((double)photon.Phi());
	hEEPt->Fill((double)Z.Pt());
	hEEMass->Fill((double)Z.M());
	hEERapidity->Fill((double)Z.Rapidity());
	hDelta_weR->Fill(min(deltaR_m,deltaR_p));
	hTriobject_wePt->Fill((double)Zgamma.Pt());
	hTriobject_weEta->Fill((double)Zgamma.Rapidity());
	hTriobject_wePhi->Fill((double)Zgamma.Phi());
	hTriobject_weMass->Fill((double)Zgamma.M());
      }
    }
    if (muonpt > 15 && antimuonpt >15 && photonpt >10) {
      Z = muon + antimuon;
      deltaR_m = sqrt(
		      pow((photon.Eta() - muon.Eta()),2) +
		      pow((photon.Phi() - muon.Phi()),2)
		      );
      deltaR_p = sqrt(
		      pow((photon.Eta() - antimuon.Eta()),2) +
		      pow((photon.Phi() - antimuon.Phi()),2)
		      );
      if(min(deltaR_m,deltaR_p) > 0.05){
	hMuonMinusPt->Fill(muonpt);
	hMuonMinusEta->Fill(muoneta);
	hMuonMinusPhi->Fill(muonphi);
	hMuonPlusPt->Fill(antimuonpt);
	hMuonPlusEta->Fill(antimuoneta);
	hMuonPlusPhi->Fill(antimuonphi);
	Zgamma = Z + photon;
	hPhoton_wmPt->Fill((double)photon.Pt());
	hPhoton_wmEta->Fill((double)photon.Eta());
	hPhoton_wmPhi->Fill((double)photon.Phi());
	hMMPt->Fill((double)Z.Pt());
	hMMMass->Fill((double)Z.M());
	hMMRapidity->Fill((double)Z.Rapidity());
	hDelta_wmR->Fill(min(deltaR_m,deltaR_p));
	hTriobject_wmPt->Fill((double)Zgamma.Pt());
	hTriobject_wmEta->Fill((double)Zgamma.Rapidity());
	hTriobject_wmPhi->Fill((double)Zgamma.Phi());
	hTriobject_wmMass->Fill((double)Zgamma.M());
      }
    }
    if (taupt > 15 && antitaupt >15 && photonpt >10) {
      Z = tau + antitau;
      deltaR_m = sqrt(
		      pow((photon.Eta() - tau.Eta()),2) +
		      pow((photon.Phi() - tau.Phi()),2)
		      );
      deltaR_p = sqrt(
		      pow((photon.Eta() - antitau.Eta()),2) +
		      pow((photon.Phi() - antitau.Phi()),2)
		      );
      if(min(deltaR_m,deltaR_p) > 0.05){
	hTauMinusPt->Fill(taupt);
	hTauMinusEta->Fill(taueta);
	hTauMinusPhi->Fill(tauphi);
	hTauPlusPt->Fill(antitaupt);
	hTauPlusEta->Fill(antitaueta);
	hTauPlusPhi->Fill(antitauphi);
	hTTPhi->Fill((double)Z.Phi());
	Zgamma = Z + photon;
	hPhoton_wtPt->Fill((double)photon.Pt());
	hPhoton_wtEta->Fill((double)photon.Eta());
	hPhoton_wtPhi->Fill((double)photon.Phi());
	hTTPt->Fill((double)Z.Pt());
	hTTMass->Fill((double)Z.M());
	hTTRapidity->Fill((double)Z.Rapidity());
	hDelta_wtR->Fill(min(deltaR_m,deltaR_p));
	hTriobject_wtPt->Fill((double)Zgamma.Pt());
	hTriobject_wtEta->Fill((double)Zgamma.Rapidity());
	hTriobject_wtPhi->Fill((double)Zgamma.Phi());
	hTriobject_wtMass->Fill((double)Zgamma.M());
      }
    }
  }
  plots.Write();
  
  // int entries = t->GetEntries();
  // entries = 3;
  // for (int _i = 0 ; _i < entries ; _i++){
  //   t->GetEntry(_i);
  //   //    for (int _j = 0 ; _j < _particlesize ; _j++){
  //   //      std::cout << _particle[] <<", ";// std::endl;
  //     // }// particle loop
  // }//entry loop
  //}//lepton flavour loop
}
