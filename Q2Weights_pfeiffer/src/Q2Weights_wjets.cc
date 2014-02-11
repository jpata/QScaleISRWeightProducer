// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "TMath.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "LHAPDF/LHAPDF.h"

#include "testDown_8var.hpp"
#include "testUp_8var.hpp"

//
// class declaration
//

class Q2WeightsWjets : public edm::EDProducer{
   public:
      explicit Q2WeightsWjets(const edm::ParameterSet&);
      ~Q2WeightsWjets();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:

    struct GenJetHigherPt {
      bool operator() (const reco::GenJet& j1, const reco::GenJet& j2) const {
	return j1.pt() > j2.pt();
      };
    };
    struct GenLeptonHigherPt {
      bool operator() (const math::XYZTLorentzVector l1, const math::XYZTLorentzVector l2) const {
	return l1.pt() > l2.pt();
      };
    };
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  TH1F *hweight_up;
  TH1F *hweight_down;
  
  TH1F* NN_w_up;
  TH1F* NN_w_down;

  testUp_8var tU;
  testDown_8var tD;



 };

//
// constructors and destructor
//
Q2WeightsWjets::Q2WeightsWjets(const edm::ParameterSet& iConfig)
{
  LHAPDF::initPDFSet("cteq6ll.LHpdf");
  LHAPDF::usePDFMember(1,0);

  //TFile *file_weights = new TFile("Wjets_Q2_weights_Jet1_Qscale_pdf.root");
  TFile *file_weights = new TFile(iConfig.getParameter<edm::FileInPath>("fileWeights").fullPath().c_str());
  hweight_up = (TH1F*)file_weights->Get("hweight_up");
  hweight_down = (TH1F*)file_weights->Get("hweight_down"); 

  TFile *file_weight_up = new TFile(iConfig.getParameter<edm::FileInPath>("fileWeightUp").fullPath().c_str());
  TFile *file_weight_down = new TFile(iConfig.getParameter<edm::FileInPath>("fileWeightDown").fullPath().c_str());
  //TFile *file_weight_up = new TFile("simuUp_8var_pdf.root");
  //TFile *file_weight_down = new TFile("simuDown_8var_pdf.root");
  
  NN_w_up = (TH1F*)file_weight_up->Get("bgh");
  NN_w_down = (TH1F*)file_weight_down->Get("bgh");

  produces<std::vector<double> >();
}


Q2WeightsWjets::~Q2WeightsWjets()
{
 
 
}

//
// member functions
//
 double _deltaR(math::XYZTLorentzVector v1,  math::XYZTLorentzVector v2) {
   Double_t DeltaR = 0.0;
   Double_t DeltaPhi = TMath::Abs(v1.Phi()-v2.Phi());
   if (DeltaPhi>TMath::Pi())
     DeltaPhi = 2*TMath::Pi() - DeltaPhi;
   Double_t DeltaEta = v1.Eta() - v2.Eta();
   DeltaR = sqrt(TMath::Power(DeltaPhi, 2) + TMath::Power(DeltaEta, 2));
   return DeltaR;
 }


// ------------ method called for each event  ------------
void
Q2WeightsWjets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  float pt_jet[6];
  float eta_jet[6];
  
  float m_jj[4][4];
  float dR_jj[4][4]; 
  
  float Q_scale;
  float factorization_scale;

  float Nparton;
  
  float x1,x2;
  int id1,id2;

  edm::Handle<reco::GenJetCollection> handle;
  iEvent.getByLabel("kt4GenJets",handle);
  const reco::GenJetCollection& genJets_kt4 = *(handle.product());
  
  reco::GenJetCollection sortedGenJets_kt4 = genJets_kt4;
  
  edm::Handle<GenParticleCollection> genPart;
  iEvent.getByLabel("genParticles", genPart);
  
  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByLabel("generator", genInfo); 
  const GenEventInfoProduct& genEventInfo = *(genInfo.product());
  
  edm::Handle<LHEEventProduct> lheevent;
  iEvent.getByLabel("source", lheevent);
  
  factorization_scale = lheevent->hepeup().SCALUP;
  Q_scale = genEventInfo.qScale();
  
  const gen::PdfInfo* pdf = genEventInfo.pdf();
  id1=pdf->id.first;
  id2=pdf->id.second;
  x1=pdf->x.first;
  x2=pdf->x.second;
  //std::cout << factorization_scale << "  " << Q_scale << "  " << pdf->scalePDF << std::endl;
  
  auto_ptr<vector<math::XYZTLorentzVector> > genLepton ( new vector<math::XYZTLorentzVector> );
  
  Nparton=0;

  for(size_t i = 0; i < genPart->size(); ++ i) {
    //const GenParticle & p = (*genParticles)[i];
    const Candidate & p = (*genPart)[i];
    math::XYZTLorentzVector lvec(p.p4());
    int id = p.pdgId();
    //cout << i << " " <<id << " " << p.status() << " " << p.px() << " " <<p.py() << " " <<p.pz() << " " << p.energy() << " " << p.mass() << std::endl;
    bool is_finalstatus3 = false;
    if ( p.status() == 3) { 
      if(abs(id) == 11 || abs(id)==13 || abs(id)==15) { 
	genLepton->push_back(lvec); 
      }
      
      is_finalstatus3 = true;
      
      for(size_t j=0; j < p.numberOfDaughters(); j++){
	//cout << "  " <<  p.daughter(j)->pdgId() << " " << p.daughter(j)->status() << endl;
	if(p.daughter(j)->status()==3){
	  is_finalstatus3=false;
	}
      }
      if(is_finalstatus3) {
	Nparton++;
      }   
    }
  }  

  for(int i=0; i<6;++i){
    pt_jet[i]=-9999.;
    eta_jet[i]=-9999.;
  }
  
  for(size_t i=0; i<6 && i<sortedGenJets_kt4.size(); ++i){
    bool remove = false;
    //do not use jets matched to leptons
    for(size_t j=0; j<genLepton->size();++j){
      if(_deltaR(genLepton->at(j),sortedGenJets_kt4.at(i).p4())<0.1 ) {
	remove = true;
      }
    }
    if(remove) {
      sortedGenJets_kt4.erase(sortedGenJets_kt4.begin()+i);
      i--;
    }
  }
  
  sort(sortedGenJets_kt4.begin(),sortedGenJets_kt4.end(),GenJetHigherPt());
  
  for(size_t i=0; i<6 && i<sortedGenJets_kt4.size(); ++i){
    
    double pt = sortedGenJets_kt4.at(i).p4().pt();
    if(pt<20) continue;
    double eta = sortedGenJets_kt4.at(i).p4().eta();
    pt_jet[i]=pt;
    eta_jet[i]=eta;
  }
  for(size_t i=0; i<4 && i<sortedGenJets_kt4.size(); ++i){
    for(size_t j=i+1; j<4 && j<sortedGenJets_kt4.size(); ++j){
      dR_jj[i][j] = _deltaR(sortedGenJets_kt4.at(i).p4(),sortedGenJets_kt4.at(j).p4());
      m_jj[i][j] =(sortedGenJets_kt4.at(i).p4()+sortedGenJets_kt4.at(j).p4()).M();
    }
  }

  double alpha = LHAPDF::alphasPDF(Q_scale);
  double alpha_down = LHAPDF::alphasPDF(Q_scale*0.25);
  double alpha_up = LHAPDF::alphasPDF(Q_scale*4.);

  //alhpa_S
  int num=2;
  double weight_up = pow(alpha_up,Nparton-num)/pow(alpha,Nparton-num);
  double weight_down = pow(alpha_down,Nparton-num)/pow(alpha,Nparton-num);
  
  //factorization scale pdf

  double pdf1 = LHAPDF::xfx(1,  x1, factorization_scale,  id1);
  double pdf2 = LHAPDF::xfx(1,  x2, factorization_scale,  id2);
  double pdf1_up = LHAPDF::xfx(1,  x1, factorization_scale*2.,  id1);
  double pdf2_up = LHAPDF::xfx(1,  x2, factorization_scale*2.,  id2);
  double pdf1_down = LHAPDF::xfx(1,  x1, factorization_scale*.5,  id1);
  double pdf2_down = LHAPDF::xfx(1,  x2, factorization_scale*.5,  id2);
  
  weight_up*=pdf1_up/pdf1*pdf2_up/pdf2;
  weight_down*=pdf1_down/pdf1*pdf2_down/pdf2; 

  //first jet pT
  if(pt_jet[0]>0){
    double xpt1 = pt_jet[0];
    if(xpt1>200) xpt1=199.999;
    
    weight_up *= hweight_up->GetBinContent(hweight_up->GetXaxis()->FindBin(xpt1));
    weight_down *= hweight_down->GetBinContent(hweight_down->GetXaxis()->FindBin(xpt1));
  }
 
  if(pt_jet[0]>0 && pt_jet[1]>0){
    double eta1=eta_jet[0];
    double eta2=eta_jet[1];

    double m_12 = m_jj[0][1];
    double dR_12 = dR_jj[0][1];
    
    double m_13=0;
    double m_23=0;
    double dR_13=0;
    double dR_23=0;
    
    if(sortedGenJets_kt4.size()>=3){
      m_13 = m_jj[0][2];
      m_23 = m_jj[1][2];
      dR_13 = dR_jj[0][2];
      dR_23 = dR_jj[1][2];
    }
    
    double NN_up = tU.Value(0,eta1,eta2,m_12,dR_12,m_13,dR_13,m_23,dR_23);
    double NN_down = tD.Value(0,eta1,eta2,m_12,dR_12,m_13,dR_13,m_23,dR_23);
    weight_up *= NN_w_up->GetBinContent(NN_w_up->GetXaxis()->FindBin(NN_up));
    weight_down *= NN_w_down->GetBinContent(NN_w_down->GetXaxis()->FindBin(NN_down));
  }
  
  auto_ptr<std::vector<double> > Q2WeightsWjets (new std::vector<double>);

  Q2WeightsWjets->push_back(weight_up);

  Q2WeightsWjets->push_back(weight_down);

  iEvent.put(Q2WeightsWjets);
}

// ------------ method called once each job just before starting event loop  ------------
void 
Q2WeightsWjets::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Q2WeightsWjets::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
Q2WeightsWjets::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Q2WeightsWjets::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Q2WeightsWjets::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Q2WeightsWjets::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Q2WeightsWjets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Q2WeightsWjets);
