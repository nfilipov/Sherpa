#ifndef SpinCorrelations_SpinCorrel_SpinCorrel_h
#define SpinCorrelations_SpinCorrel_SpinCorrel_h
// -*- C++ -*-
//
// Package:    SpinCorrelations/SpinCorrel
// Class:      SpinCorrel
// 
/**\class SpinCorrel SpinCorrel.h SpinCorrelations/SpinCorrel/plugins/SpinCorrel.h

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Nicolas Filipovic
//         Created:  Tue, 05 Sep 2017 15:28:26 GMT
//
//
#include <TH1.h>
#include "FWCore/TFWLiteSelector/interface/TFWLiteSelector.h"

//A worker processes the events.  When using PROOF there is one Worker per PROOF CPU Node.
struct SpinCorrelWorker {
  SpinCorrelWorker(const TList*, TList&);
  ~SpinCorrelWorker();
  void process( const edm::Event& iEvent );
  void postProcess(TList&);
  //Place histograms, etc that you want to fill here
  //TH1F* h_a;
};

//Only one Selector is made per job. It gets all the results from each worker.
class SpinCorrel : public TFWLiteSelector<SpinCorrelWorker> {
public :
  SpinCorrel();
  ~SpinCorrel();
  void begin(TList*&);
  void terminate(TList&);
    
private:
    
  SpinCorrel(SpinCorrel const&);
  SpinCorrel operator=(SpinCorrel const&);
  
  ClassDef(SpinCorrel,2)
};
#endif
