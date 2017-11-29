// -*- C++ -*-
//
// Package:    SpinCorrelations/SpinCorrel
// Class:      SpinCorrel
// 
/*

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nicolas Filipovic
//         Created:  Tue, 05 Sep 2017 15:28:26 GMT
//
//


// system include files
#include <memory>
#include <iostream>

#include "TCanvas.h"
// user include files
#include "SpinCorrelations/SpinCorrel/plugins/SpinCorrel.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"

//
// constants shared between the Selector and the Workers
//
//const char* const kA = "a";

//===============================================
//A worker processes the events.
// A new worker is created each time the events are processed
//===============================================

// ------------ constructed for each PROOF Node  ------------
// The arguments are
//   fromSelector: these are copies of values set in the selector and sent to all workers
//            out: these are the items which will be passed back to the selector (e.g. histograms)
SpinCorrelWorker::SpinCorrelWorker(const TList* fromSelector, TList& out ) {
  //h_a  = new TH1F( kA , "a"  , 100,  0, 20 );
  //out.Add(h_a);  
}

SpinCorrelWorker::~SpinCorrelWorker()
{
}
// ------------ method called for each event  ------------
void 
SpinCorrelWorker::process( const edm::Event& iEvent ) {
  using namespace edm;
  
  
//  using namespace edmtest;
//  edm::Handle<ThingCollection> hThings;
//  iEvent.getByLabel("Thing",hThings);
//  for ( ThingCollection::const_iterator it = hThings->begin(); 
//        it != hThings->end(); ++it ) {
//    h_a ->Fill( it->a );
//  }
  
}

// ------------ called after processing the events  ------------
// The argument is the same as for the constructor
void 
SpinCorrelWorker::postProcess(TList& out)
{
}


//===============================================
//Only one Selector is made per job. It gets all the results from each worker.
//===============================================
SpinCorrel::SpinCorrel()
{
}

SpinCorrel::~SpinCorrel()
{
}

// ------------ called just before all workers are constructed  ------------
void SpinCorrel::begin(TList*& toWorkers)
{
}

// ------------ called after all workers have finished  ------------
// The argument 'fromWorkers' contains the accumulated output of all Workers
void SpinCorrel::terminate(TList& fromWorkers) {
  using namespace std;
  std::auto_ptr<TCanvas> canvas( new TCanvas() );
//  {
//    TObject* hist = fromWorkers.FindObject(kA);
//    if(0!=hist) {
//      hist->Draw();
//      canvas->SaveAs( "a.jpg" );
//    } else {
//      cout <<"no '"<<kA<<"' histogram"<<endl;
//    }
 // }

  
}
