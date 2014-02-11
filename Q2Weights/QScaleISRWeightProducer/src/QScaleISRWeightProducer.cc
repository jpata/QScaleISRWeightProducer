// -*- C++ -*-
//
// Package:    QScaleISRWeightProducer
// Class:      QScaleISRWeightProducer
// 
/**\class QScaleISRWeightProducer QScaleISRWeightProducer.cc Q2Weights/QScaleISRWeightProducer/src/QScaleISRWeightProducer.cc

 Description: Module for ttbar parton shower Q-scale reweighting

 Implementation:
     The general idea is that the reweighting recipe build on PDF & alphaS can
     be improved by accounting for ISR effects which have been found to be very sensitive
     to the choosen scale within the parton showering apart from the hard ME.
     This algorithm utilizes the GenPartilce relations to calculate the invariant mass
     of all the particles radiated between the parton from the proton and the parton entering
     the hard scattering. Please note: those particle cascades are never stored by Pythia.
*/
//
// Original Author:  Matthias Komm
//         Created:  Fri Feb  7 15:30:13 CET 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TMath.h"

//
// class declaration
//

class QScaleISRWeightProducer : public edm::EDProducer {
   public:
      explicit QScaleISRWeightProducer(const edm::ParameterSet&);
      ~QScaleISRWeightProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      edm::InputTag _inputTag;
   
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      
      double calculateWeightUp(double isrmass);
      double calculateWeightDown(double isrmass);
      // ----------member data ---------------------------
      
      
      
};

//
// constants, enums and typedefs
//

const double up_p0 = 2.44131;
const double up_p1 = -0.365833;
const double up_p2 = 0.357556;
const double up_p3 = 0.0839954;

const double down_p0 = -0.000786909;
const double down_p1 = -1.71566;
const double down_p2 = 2.07759;

//
// static data member definitions
//

//
// constructors and destructor
//
QScaleISRWeightProducer::QScaleISRWeightProducer(const edm::ParameterSet& iConfig)
{
    _inputTag=iConfig.getParameter<edm::InputTag>("src");


    produces<double>("isrMass");
    produces<double>("weightUp");
    produces<double>("weightDown");
    produces<std::vector<double>>("isrPtList");
}


QScaleISRWeightProducer::~QScaleISRWeightProducer()
{
 
    

}


//
// member functions
//

// ------------ method called to produce the data  ------------

double QScaleISRWeightProducer::calculateWeightUp(double isrmass)
{
    return up_p0+up_p1*TMath::Power(isrmass,up_p2)+up_p3*TMath::Sqrt(isrmass);
}

double QScaleISRWeightProducer::calculateWeightDown(double isrmass)
{
    return TMath::Exp(down_p0*isrmass)*down_p1+down_p2;
}

void
QScaleISRWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle<reco::GenParticleCollection> genParticleCollection;
    iEvent.getByLabel(_inputTag, genParticleCollection);
    
    std::auto_ptr<double> isrMass(new double(0));
    std::auto_ptr<double> weightUp(new double(0));
    std::auto_ptr<double> weightDown(new double(0));
    std::auto_ptr<std::vector<double>> isrPtList(new std::vector<double>(2));
    std::vector<const reco::GenParticle*> afterISRParticles;
    
    for (unsigned iparticle=0; iparticle<genParticleCollection->size();++iparticle)
    {
        const reco::GenParticle& genObject = (*genParticleCollection)[iparticle];
        //find the first ME particles 
        if (genObject.numberOfDaughters()>=2)
        {
            int ntops=0;
            for (unsigned idaughter=0; idaughter<genObject.numberOfDaughters();++idaughter)
            {
                //ME-input-particles have the two tops + add. partons as daughters
                if (fabs(fabs(genObject.daughter(idaughter)->pdgId())-6)<0.5)
                {
                    ++ntops;
                }
            }
            if (ntops==2)
            {
                afterISRParticles.push_back(&genObject);
            }
        }
    }

    if (afterISRParticles.size()==2)
    {
        reco::Candidate::LorentzVector isr1=afterISRParticles[0]->mother(0)->p4()-afterISRParticles[0]->p4();
        reco::Candidate::LorentzVector isr2=afterISRParticles[1]->mother(0)->p4()-afterISRParticles[1]->p4();
        
        *isrMass=(isr1+isr2).mass();
        
        *weightUp=calculateWeightUp(*isrMass);
        *weightDown=calculateWeightDown(*isrMass);
        
        if (isr1.pt()>isr2.pt())
        {
            (*isrPtList)[0]=isr1.pt();
            (*isrPtList)[1]=isr2.pt();
        } 
        else 
        {
            (*isrPtList)[1]=isr1.pt();
            (*isrPtList)[0]=isr2.pt();
        }
    }
    else
    {
        throw cms::Exception("ME-input-particles should only be 2");
    }
    iEvent.put(isrMass,"isrMass");
    iEvent.put(weightUp,"weightUp");
    iEvent.put(weightDown,"weightDown");
    iEvent.put(isrPtList,"isrPtList");
    
}

// ------------ method called once each job just before starting event loop  ------------
void 
QScaleISRWeightProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QScaleISRWeightProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
QScaleISRWeightProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
QScaleISRWeightProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
QScaleISRWeightProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
QScaleISRWeightProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QScaleISRWeightProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(QScaleISRWeightProducer);
