/*class TauValidation
 *  
 *  Class to fill dqm monitor elements from existing EDM file
 *
 */
#include "Validation/EventGenerator/interface/TauValidation.h"
#include "CLHEP/Units/defs.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Validation/EventGenerator/interface/TauDecay_GenParticle.h"
#include "Validation/EventGenerator/interface/PdtPdgMini.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include <iostream>
using namespace edm;

TauValidation::TauValidation(const edm::ParameterSet& iPSet): 
 // wmanager_(iPSet,consumesCollector())
  genparticleCollection_(iPSet.getParameter<edm::InputTag>("genparticleCollection"))
  //,hepmcCollection_(iPSet.getParameter<edm::InputTag>("hepmcCollection"))
  ,NMODEID(TauDecay::NMODEID-1)// fortran to C++ index
  ,zsbins(20)
  ,zsmin(-0.5)
  ,zsmax(0.5)
{    
  //hepmcCollectionToken_=consumes<HepMCProduct>(hepmcCollection_);
  genparticleCollectionToken_=consumes<reco::GenParticleCollection>(genparticleCollection_);


  edm::Service<TFileService> fs;
    histo_pt_ = fs->make<TH1F> ("Pt", "Pt",100,0,100);


}

TauValidation::~TauValidation(){}


void TauValidation::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup){ 
  ///Gathering the reco::GenParticleCollection information
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genparticleCollectionToken_, genParticles );
  
  double weight = 1.;//wmanager_.weight(iEvent);
  //////////////////////////////////////////////
  // find taus#system($cmd);
  //


  for (reco::GenParticleCollection::const_iterator iter = genParticles->begin(); iter != genParticles->end(); ++iter) {
    
    if(abs(iter->pdgId())==PdtPdgMini::Z0 || abs(iter->pdgId())==PdtPdgMini::Higgs0){
     if (!isLastTauinChain(&(*iter))) continue;
      mmc_method(&(*iter),weight);
    }
  }
}//analyze

const reco::GenParticle* TauValidation::GetMother(const reco::GenParticle* tau){
  for (unsigned int i=0;i<tau->numberOfMothers();i++) {
    const reco::GenParticle *mother=static_cast<const reco::GenParticle*>(tau->mother(i));
    if(mother->pdgId() == tau->pdgId()) return GetMother(mother);
    return mother;
  }
  return tau;
}
const reco::GenParticle* TauValidation::GetLastSelf2(const reco::GenParticle* particle){
     for (unsigned int i=0;i<particle->numberOfDaughters();i++) {
          const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(particle->daughter(i));
          if(dau->pdgId() == particle->pdgId()) return GetLastSelf2(dau);
     }
     return particle;
}



const std::vector<const reco::GenParticle*> TauValidation::GetMothers(const reco::GenParticle* boson){
  std::vector<const reco::GenParticle*> mothers;
  for (unsigned int i=0;i<boson->numberOfMothers();i++) {
    const reco::GenParticle *mother=static_cast<const reco::GenParticle*>(boson->mother(i));
    if(mother->pdgId() == boson->pdgId()) return GetMothers(mother);
    mothers.push_back(mother);
  }
  return mothers;
}

int TauValidation::findMother(const reco::GenParticle* tau){
  return TauValidation::GetMother(tau)->pdgId();
}

bool TauValidation::isLastTauinChain(const reco::GenParticle* tau){
  for(unsigned int i = 0; i <tau->numberOfDaughters(); i++){
    if(tau->daughter(i)->pdgId() == tau->pdgId()) return false;
  }
  return true;
}

void TauValidation::findTauList(const reco::GenParticle* tau,std::vector<const reco::GenParticle*> &TauList){
  TauList.insert(TauList.begin(),tau);
  for(unsigned int i=0;i<tau->numberOfMothers();i++) {
    const reco::GenParticle *mother=static_cast<const reco::GenParticle*>(tau->mother(i));
    if(mother->pdgId() == tau->pdgId()){
      findTauList(mother,TauList);
    }
  }
}

void TauValidation::findFSRandBrem(const reco::GenParticle* p, bool doBrem, std::vector<const reco::GenParticle*> &ListofFSR,
				  std::vector<const reco::GenParticle*> &ListofBrem){
  // note this code split the FSR and Brem based one if the tau decays into a tau+photon or not with the Fortran Tauola Interface, this is not 100% correct because photos puts the tau with the regular tau decay products. 
  if(abs(p->pdgId())==15){
    if(isLastTauinChain(p)){ doBrem=true;}
    else{ doBrem=false;}
  }
  int photo_ID=22;
  for(unsigned int i = 0; i <p->numberOfDaughters(); i++){
    const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(p->daughter(i));
    if(abs((dau)->pdgId()) == abs(photo_ID) && !doBrem){ListofFSR.push_back(dau);}
    if(abs((dau)->pdgId()) == abs(photo_ID) && doBrem){ListofBrem.push_back(dau);}
    if(abs((dau)->pdgId()) != 111 && abs((dau)->pdgId()) != 221){ // remove pi0 and eta decays
      findFSRandBrem(dau,doBrem,ListofFSR,ListofBrem);
    }
  }
}

void TauValidation::FindPhotosFSR(const reco::GenParticle* p,std::vector<const reco::GenParticle*> &ListofFSR,double &BosonScale){
  BosonScale=0.0;
  const reco::GenParticle* m=GetMother(p);
  int mother_pid=m->pdgId();
  if(m->pdgId()!=p->pdgId()){
    for(unsigned int i=0; i <m->numberOfDaughters(); i++){
      const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(m->daughter(i));
      if(abs(dau->pdgId()) == 22) {
	ListofFSR.push_back(dau);
      }
    }
  }
  if(abs(mother_pid) == 24) BosonScale=1.0; // W
  if(abs(mother_pid) == 23) BosonScale=2.0; // Z;
  if(abs(mother_pid) == 22) BosonScale=2.0; // gamma;
  if(abs(mother_pid) == 25) BosonScale=2.0; // HSM;
  if(abs(mother_pid) == 35) BosonScale=2.0; // H0;
  if(abs(mother_pid) == 36) BosonScale=2.0; // A0;
  if(abs(mother_pid) == 37) BosonScale=1.0; //Hpm;
}

int TauValidation::tauMother(const reco::GenParticle* tau, double weight){
  if(abs(tau->pdgId()) != 15 ) return -3;
  int mother_pid = findMother(tau);
  if(mother_pid == -2) return -2;
  int label = other;
  if(abs(mother_pid) == 24) label = W;
  if(abs(mother_pid) == 23) label = Z;
  if(abs(mother_pid) == 22) label = gamma;
  if(abs(mother_pid) == 25) label = HSM;
  if(abs(mother_pid) == 35) label = H0;
  if(abs(mother_pid) == 36) label = A0;
  if(abs(mother_pid) == 37) label = Hpm;
  int mother_shortpid=(abs(mother_pid)%10000);
  if(mother_shortpid>500 && mother_shortpid<600 )label = B;
  if(mother_shortpid>400 && mother_shortpid<500)label = D;
  if(label==B || label == D || label == other) return -1;
  return mother_pid;
}

int TauValidation::tauDecayChannel(const reco::GenParticle* tau,int jak_id, unsigned int TauBitMask, double weight){
  int channel = undetermined; 
  if(tau->status() == 1) channel = stable; 
  int allCount   = 0, 
    eCount     = 0, 
    muCount    = 0, 
    pi0Count   = 0, 
    piCount    = 0, 
    rhoCount   = 0, 
    a1Count    = 0, 
    KCount     = 0, 
    KstarCount = 0; 

  countParticles(tau,allCount,eCount,muCount,pi0Count,piCount,rhoCount,a1Count,KCount,KstarCount);

  // resonances   
  if(KCount >= 1)     channel = K; 
  if(KstarCount >= 1) channel = Kstar; 
  if(a1Count >= 1)    channel = a1; 
  if(rhoCount >= 1)   channel = rho; 
 // if(channel!=undetermined && weight!=0.0) TauDecayChannels->Fill(channel,weight); 
  
  // final state products 
  if(piCount == 1 && pi0Count == 0) channel = pi; 
  if(piCount == 1 && pi0Count == 1) channel = pi1pi0; 
  if(piCount == 1 && pi0Count > 1)  channel = pinpi0; 
  if(piCount == 3 && pi0Count == 0) channel = tripi; 
  if(piCount == 3 && pi0Count > 0)  channel = tripinpi0; 
  if(eCount == 1)                   channel = electron; 
  if(muCount == 1)                  channel = muon; 
 // if(weight!=0.0) TauDecayChannels->Fill(channel,weight); 
  return channel; 
} 

void TauValidation::countParticles(const reco::GenParticle* p,int &allCount, int &eCount, int &muCount,
				   int &pi0Count,int &piCount,int &rhoCount,int &a1Count,int &KCount,int &KstarCount){
  for(unsigned int i=0; i<p->numberOfDaughters(); i++){
    const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(p->daughter(i));
    int pid = dau->pdgId();
    allCount++;
    if(abs(pid) == 11)    eCount++;
    if(abs(pid) == 13)    muCount++;
    if(abs(pid) == 111)   pi0Count++;
    if(abs(pid) == 211)   piCount++;
    if(abs(pid) == 213)   rhoCount++;
    if(abs(pid) == 20213) a1Count++;
    if(abs(pid) == 321)   KCount++;
    if(abs(pid) == 323)   KstarCount++;
    countParticles(dau,allCount,eCount,muCount,pi0Count,piCount,rhoCount,a1Count,KCount,KstarCount);
  }
}




void TauValidation::mmc_method(const reco::GenParticle* boson, double weight){
     TLorentzVector taum2(0,0,0,0);
     TLorentzVector taup2(0,0,0,0);



     int ntau(0);
     //  std::cout<<"bosons daughters ";
     for(unsigned int i = 0; i <boson->numberOfDaughters(); i++){
          const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(boson->daughter(i));




          //std::cout<<dau->pdgId()<<" ";
          //  if(ntau==1 && dau->pdgId() == 15)return;
          if(boson->pdgId()!= 15 && abs(dau->pdgId()) == 15) {ntau++;
               //   Particles.clear();
               //  Particles.push_back(static_cast<const reco::GenParticle*>(dau));
               // GetLastSelf(static_cast<const reco::GenParticle*>(dau));
               //  dau=Particles.back();
               dau=GetLastSelf2(dau);

          }
          if ( dau->pdgId() == 15) taup2.SetXYZT(dau->px(),dau->py(),dau->pz(),dau->energy());
          if ( dau->pdgId() == -15) taum2.SetXYZT(dau->px(),dau->py(),dau->pz(),dau->energy());

     }
     //      std::cout<<"bosons daughters end."<<std::endl;

     if(ntau!=2) return; 
     if(abs(boson->pdgId())==PdtPdgMini::Z0 ){
          TLorentzVector tautau(0,0,0,0);
          TLorentzVector pipi(0,0,0,0);
          TLorentzVector taum(0,0,0,0);
          TLorentzVector taup(0,0,0,0);
          TLorentzVector rho_plus,rho_minus,pi_rhominus,pi0_rhominus,pi_rhoplus,pi0_rhoplus,pi_plus,pi_minus;
          //int nSinglePionDecays(0),nSingleMuonDecays(0),nSingleElectronDecays(0);
          //double x1(0),x2(0); 
          TLorentzVector Zboson(boson->px(),boson->py(),boson->pz(),boson->energy());
          for(unsigned int i = 0; i <boson->numberOfDaughters(); i++){
               const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(boson->daughter(i));

               int pid = dau->pdgId();
               //            std::cout<<"hi1................."<<std::endl;

               if(abs(findMother(dau)) != 15 && abs(pid) == 15){

                    dau=GetLastSelf2(dau);
                    /*Particles.clear();
                      Particles.push_back(static_cast<const reco::GenParticle*>(dau));
                      GetLastSelf(static_cast<const reco::GenParticle*>(dau));
                      dau=Particles.back();
                      */
                    TauDecay_GenParticle TD;
                    unsigned int jak_id, TauBitMask;
                    //	std::cout<<"hi1................."<<std::endl;
                    if(TD.AnalyzeTau(dau,jak_id,TauBitMask,false,false)){
                         //	std::cout<<"hi2..................."<<std::endl;
                         std::vector<const reco::GenParticle*> part=TD.Get_TauDecayProducts();




                         if(jak_id==TauDecay::MODE_MUON){


                              std::cout<<"Muonic decay product ids \n";
                              for(unsigned int i=0; i<part.size();i++) std::cout<< part.at(i)->pdgId()<<" ";
                              std::cout<<"End Muonic decay "<<std::endl;


                         }

                         if(jak_id==TauDecay::MODE_ELECTRON){
                         }


                         //TLorentzVector LVtau(dau->px(),dau->py(),dau->pz(),dau->energy());

                         if(jak_id==TauDecay::MODE_PION){
                              std::cout<<"Pionic decay product ids \n";
                              for(unsigned int i=0; i<part.size();i++)  {
                                   std::cout<< part.at(i)->pdgId()<<" ";
                                   histo_pt_->Fill(part.at(i)->pt());
                              }
                              std::cout<<"End Pionic decay "<<std::endl;
                         }


                          


                    }
               }




          }

     }

}


void TauValidation::GetLastSelf(const reco::GenParticle *Particle){
     for (unsigned int i=0; i< Particle->numberOfDaughters(); i++){
          const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(Particle->daughter(i));
          if(Particle->pdgId()==dau->pdgId()){
               Particles.push_back(dau);
               Particle=dau;
               GetLastSelf(Particle);
          }
     }
}


double TauValidation::Zstoa(double zs){
  double a=1-sqrt(fabs(1.0-2*fabs(zs)));
  if(zs<0){
    a*=-1.0;
  }
  return a;
}


double TauValidation::leadingPionMomentum(const reco::GenParticle* tau, double weight){
  return leadingPionP4(tau).P();
}

TLorentzVector TauValidation::leadingPionP4(const reco::GenParticle* tau){
  TLorentzVector p4(0,0,0,0);
  for(unsigned int i = 0; i <tau->numberOfDaughters(); i++){
    const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(tau->daughter(i));
    int pid = dau->pdgId();
    if(abs(pid) == 15) return leadingPionP4(dau);
    if(!(abs(pid)==211 || abs(pid)==13 || abs(pid)==11)) continue;
    if(dau->p() > p4.P()) p4 = TLorentzVector(dau->px(),dau->py(),dau->pz(),dau->energy());
  }
  return p4;
}

TLorentzVector TauValidation::motherP4(const reco::GenParticle* tau){
  const reco::GenParticle* m=GetMother(tau);
  return TLorentzVector(m->px(),m->py(),m->pz(),m->energy());
}

double TauValidation::visibleTauEnergy(const reco::GenParticle* tau){
  TLorentzVector p4(tau->px(),tau->py(),tau->pz(),tau->energy());
  for(unsigned int i = 0; i <tau->numberOfDaughters(); i++){
    const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(tau->daughter(i));
    int pid = dau->pdgId();
    if(abs(pid) == 15) return visibleTauEnergy(dau);
    if(abs(pid) == 12 || abs(pid) == 14 || abs(pid) == 16) {
      p4-=TLorentzVector(dau->px(),dau->py(),dau->pz(),dau->energy());
    }
  }
  return p4.E();
}


