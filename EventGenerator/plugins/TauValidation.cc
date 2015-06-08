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
#include "TMath.h"

#include <iostream>
using namespace edm;

 bool custom_isnan(double var)
{
     volatile double d = var;
     return d != d;

}


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
mass_Z= fs->make<TH1F> ("mass Z", "mass Z", 80,40,150);

mtautau= fs->make<TH1F> ("ditau", "ditau", 80,40,150);
mtautau11= fs->make<TH1F> ("ditau11", "ditau11", 80,40,150);
mtautau12= fs->make<TH1F> ("ditau12", "ditau12", 80,40,150);
mtautau21= fs->make<TH1F> ("ditau21", "ditau21", 80,40,150);
mtautau22= fs->make<TH1F> ("ditau22", "ditau22", 80,40,150);


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

//std::cout<<"begin "<<std::endl;
  for (reco::GenParticleCollection::const_iterator iter = genParticles->begin(); iter != genParticles->end(); ++iter) {
  // std::cout<<iter->pdgId()<<std::endl; 
    if(abs(iter->pdgId())==PdtPdgMini::Z0 || abs(iter->pdgId())==PdtPdgMini::Higgs0){
 //    if (!isLastTauinChain(&(*iter))) continue;
      //std::cout<<" Z ............. "<<std::endl;
      mmc_method_interface(&(*iter),weight);
    }
  }
//std::cout<<"end "<<std::endl;
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




void TauValidation::mmc_method_interface(const reco::GenParticle* boson, double weight){



     int ntau(0);
     bool jak_muon=false;
     bool jak_pion=false;
     //  std::cout<<"bosons daughters ";

TLorentzVector p4vis1_,p4vis2_,p4miss2_,p4miss1_;
TLorentzVector sol_p4miss11_,sol_p4miss12_,solv_p4miss21_,sol_p4miss22_;

for(unsigned int i = 0; i <boson->numberOfDaughters(); i++){
     const reco::GenParticle *dau=static_cast<const reco::GenParticle*>(boson->daughter(i));

     if(ntau==1 && dau->pdgId() == 15)return;

     if(boson->pdgId()!= 15 && abs(dau->pdgId()) == 15) ntau++;


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
//                      std::cout<<"hi1................."<<std::endl;

               if(abs(findMother(dau)) != 15 && abs(pid) == 15){

                    dau=GetLastSelf2(dau);
                    /*Particles.clear();
                      Particles.push_back(static_cast<const reco::GenParticle*>(dau));
                      GetLastSelf(static_cast<const reco::GenParticle*>(dau));
                      dau=Particles.back();
                      */
                    TauDecay_GenParticle TD;
                    unsigned int jak_id, TauBitMask;
  //                  	std::cout<<"hi2................."<<std::endl;
                    if(TD.AnalyzeTau(dau,jak_id,TauBitMask,false,false)){
    //                     	std::cout<<"hi3..................."<<std::endl;
                         std::vector<const reco::GenParticle*> part=TD.Get_TauDecayProducts();




                         if(jak_id==TauDecay::MODE_MUON){
                              jak_muon=true;    
  //                            std::cout<<"Muonic decay product ids \n";
                              for(unsigned int i=0; i<part.size();i++) { //std::cout<< part.at(i)->pdgId()<<" ";
                                   if (fabs(part.at(i)->pdgId())==14 || fabs(part.at(i)->pdgId())==16 ) p4miss2_+=TLorentzVector(part.at(i)->px(),part.at(i)->py(),part.at(i)->pz(),part.at(i)->energy()) ;
                                   else if (fabs(part.at(i)->pdgId())==13) p4vis2_+=TLorentzVector(part.at(i)->px(),part.at(i)->py(),part.at(i)->pz(),part.at(i)->energy());
                              }
    //                          std::cout<<"End Muonic decay "<<std::endl;


                         }

                         if(jak_id==TauDecay::MODE_ELECTRON){
                         }



                         if(jak_id==TauDecay::MODE_PION){
                              jak_pion=true;
                            //  std::cout<<"Pionic decay product ids \n";
                              for(unsigned int i=0; i<part.size();i++)  {
                                   if ( fabs(part.at(i)->pdgId())==16 ) p4miss1_+=TLorentzVector(part.at(i)->px(),part.at(i)->py(),part.at(i)->pz(),part.at(i)->energy()) ;
                                   else if (fabs(part.at(i)->pdgId())==211) p4vis1_+=TLorentzVector(part.at(i)->px(),part.at(i)->py(),part.at(i)->pz(),part.at(i)->energy());
                              
                              
                              
                              //    std::cout<< part.at(i)->pdgId()<<" ";
                                //   histo_pt_->Fill(part.at(i)->pt());
                              }
                            // std::cout<<"End Pionic decay "<<std::endl;
                         }


                          


                    }
               }




          }

     }

     if (jak_muon && jak_pion) {
          mass_Z->Fill((p4vis1_+p4miss1_+p4vis2_+ p4miss2_).M());
          std::cout<<"Boom!..  I found the channel interested \n";
          TH1F mass_dummy_("","",80,50,120);
          if (mmc_method(p4vis1_,p4miss1_,p4vis2_, p4miss2_,mass_dummy_)) {

               double tautau_m=mass_dummy_.GetBinCenter(mass_dummy_.GetMaximumBin());
               //mtautau->Fill(tautau_m);

               if (mass_dummy_.Integral()>0){
//                    mass_dum_count++;
                    double sum_mass=0,sum_weight=0;
                    for (int i=mass_dummy_.GetMaximumBin()-2;i<mass_dummy_.GetMaximumBin()+3;i++){
                         double iweight=mass_dummy_.GetBinContent(i);
                         sum_weight=sum_weight+iweight;
                         double mss=mass_dummy_.GetBinCenter(i);
                         sum_mass=sum_mass+mss*iweight;
                    }


                    mtautau->Fill(sum_mass/sum_weight);

               }






          }else {
               std::cout<<"No mmc solutions \n";
          }

}
}

//bool mmc_method(const TLorentzVector &p4vis1,const TLorentzVector &p4miss1,const TLorentzVector &p4vis2,const TLorentzVector &p4mis2, TLorentzVector &sol_p4miss11, TLorentzVector &sol_p4miss12, TLorentzVector &solv_p4miss21,TLorentzVector &sol_p4miss22){
bool  TauValidation::mmc_method(const TLorentzVector &p4vis1,const TLorentzVector &p4miss1,const TLorentzVector &p4vis2,const TLorentzVector &p4miss2, TH1F &mass_dummy) {
     double p_miss1, phi_miss1, theta_miss1,p_miss2,theta_miss2, phi_miss2;
     double p_vis1,theta_vis1, phi_vis1;
     double m_vis1;
     double m_tau=1.777;
     double p_vis2,theta_vis2,phi_vis2;
     double m_vis2;
     double m_miss2,m_miss1;
     TLorentzVector p4miss11,p4miss12,p4miss21,p4miss22;
     
     
     p_miss1=p4miss1.P();phi_miss1=p4miss1.Phi();theta_miss1=p4miss1.Theta();
     p_miss2=p4miss2.P();phi_miss2=p4miss2.Phi();theta_miss2=p4miss2.Theta();
     p_vis1=p4vis1.P();theta_vis1=p4vis1.Theta(); phi_vis1=p4vis1.Phi();
     p_vis2=p4vis2.P();theta_vis2=p4vis2.Theta(); phi_vis2=p4vis2.Phi();
     m_vis1=p4vis1.M();m_vis2=p4vis2.M();
     m_miss1=0;m_miss2=p4miss2.M();
     
     bool check=false;
     
     for (float kk=0;kk<1.5;kk+= 0.1)
          for (float ii=-TMath::Pi();ii<TMath::Pi();ii+=0.1)
               for (float jj=TMath::Pi();jj>-TMath::Pi();jj-=0.1 ){



                    //m_miss1=0;




                    phi_miss1=ii;
                    phi_miss2=jj;
                    m_miss1=0;
                    m_miss2=kk;

                    double METX= p_miss1*sin(theta_miss1)*cos(phi_miss1)+p_miss2*sin(theta_miss2)*cos(phi_miss2),METY=p_miss1*sin(theta_miss1)*sin(phi_miss1)+p_miss2*sin(theta_miss2)*sin(phi_miss2);
                    // METX=Ran.Gaus(METX,fabs(METX*6/100));
                    //METY=Ran.Gaus(METY,fabs(METY*6/100));

                    double pt_miss1=(-sin(phi_miss2)*METX+cos(phi_miss2)*METY)/sin(phi_miss1-phi_miss2);
                    double pt_miss2=(sin(phi_miss1)*METX-cos(phi_miss1)*METY)/sin(phi_miss1-phi_miss2);




                    //cout<<"En  "<<p_vis1**2+m_vis1**2<<" "<<p_miss1**2+m_miss1**2<<endl;


                    double pz11=(1./(2*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(p_vis1,2))))*(-cos(theta_vis1)*(pow(m_miss1,2)+pow(m_vis1,2)-pow(m_tau,2))*p_vis1+cos(phi_miss1-phi_vis1)*sin(2*theta_vis1)*pow(p_vis1,2)*pt_miss1+sqrt((pow(m_vis1,2)+pow(p_vis1,2))*(pow(m_miss1,4)+pow(pow(m_vis1,2)-pow(m_tau,2),2)+4*     cos(phi_miss1-phi_vis1)*sin(theta_vis1)*(-pow(m_vis1,2)+pow(m_tau,2))*p_vis1*pt_miss1-4*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(sin(phi_miss1-phi_vis1),2)*pow(p_vis1,2))*pow(pt_miss1,2)-2*pow(m_miss1,2)*(pow(m_vis1,2)+pow(m_tau,2)+2*sin(theta_vis1)*p_vis1*(sin(theta_vis1)*p_vis1+cos(phi_miss1-phi_vis1)*pt_miss1)))));

                    double pz12=(-1./(2*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(p_vis1,2))))*(cos(theta_vis1)*(pow(m_miss1,2)+pow(m_vis1,2)-pow(m_tau,2))*p_vis1-cos(phi_miss1-phi_vis1)*sin(2*theta_vis1)*pow(p_vis1,2)*pt_miss1+sqrt((pow(m_vis1,2)+pow(p_vis1,2))*(pow(m_miss1,4)+pow(pow(m_vis1,2)-pow(m_tau,2),2)+4*     cos(phi_miss1-phi_vis1)*sin(theta_vis1)*(-pow(m_vis1,2)+pow(m_tau,2))*p_vis1*pt_miss1-4*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(sin(phi_miss1-phi_vis1),2)*pow(p_vis1,2))*pow(pt_miss1,2)-2*pow(m_miss1,2)*(pow(m_vis1,2)+pow(m_tau,2)+2*sin(theta_vis1)*p_vis1*(sin(theta_vis1)*p_vis1+cos(phi_miss1-phi_vis1)*pt_miss1)))));

                    double pz21=(1./(2*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(p_vis2,2))))*(-cos(theta_vis2)*(pow(m_miss2,2)+pow(m_vis2,2)-pow(m_tau,2))*p_vis2+cos(phi_miss2-phi_vis2)*sin(2*theta_vis2)*pow(p_vis2,2)*pt_miss2+sqrt((pow(m_vis2,2)+pow(p_vis2,2))*(pow(m_miss2,4)+pow(pow(m_vis2,2)-pow(m_tau,2),2)+4*     cos(phi_miss2-phi_vis2)*sin(theta_vis2)*(-pow(m_vis2,2)+pow(m_tau,2))*p_vis2*pt_miss2-4*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(sin(phi_miss2-phi_vis2),2)*pow(p_vis2,2))*pow(pt_miss2,2)-2*pow(m_miss2,2)*(pow(m_vis2,2)+pow(m_tau,2)+2*sin(theta_vis2)*p_vis2*(sin(theta_vis2)*p_vis2+cos(phi_miss2-phi_vis2)*pt_miss2)))));

                    double pz22=(-1./(2*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(p_vis2,2))))*(cos(theta_vis2)*(pow(m_miss2,2)+pow(m_vis2,2)-pow(m_tau,2))*p_vis2-cos(phi_miss2-phi_vis2)*sin(2*theta_vis2)*pow(p_vis2,2)*pt_miss2+sqrt((pow(m_vis2,2)+pow(p_vis2,2))*(pow(m_miss2,4)+pow(pow(m_vis2,2)-pow(m_tau,2),2)+4*     cos(phi_miss2-phi_vis2)*sin(theta_vis2)*(-pow(m_vis2,2)+pow(m_tau,2))*p_vis2*pt_miss2-4*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(sin(phi_miss2-phi_vis2),2)*pow(p_vis2,2))*pow(pt_miss2,2)-2*pow(m_miss2,2)*(pow(m_vis2,2)+pow(m_tau,2)+2*sin(theta_vis2)*p_vis2*(sin(theta_vis2)*p_vis2+cos(phi_miss2-phi_vis2)*pt_miss2)))));

                    //cout<<custom_isnan(pz11)<<" "<<  custom_isnan(pz12)<<" "<<  custom_isnan(pz21) <<" "<<custom_isnan(pz22) <<endl;

                    if (!((custom_isnan(pz11) && custom_isnan(pz12)) || (custom_isnan(pz21) && custom_isnan(pz22)  ))) {

                    //cout<<"Phivalues "<<phi_miss1<<" "<<phi_miss2<<endl;


                    //TLorenzVector p4
                    //cout<<p_vis1<< theta_vis1<<phi_vis1<<m_vis1<<p_miss1<< theta_miss1<<phi_miss1<<m_miss1<<p_vis2<< theta_vis2<<phi_vis2<<m_vis2<<p_miss2<< theta_miss2<<phi_miss2<<m_miss2<<endl;

                    //double zm=sqrt(-2.);
                    // cout<<" isnan "<<custom_isnan(sqrt(-2.))<<" "<<zm<<endl;



                    //double p_miss11_cal=sqrt(z11**2-m_miss1**2);double p_miss12_cal=sqrt(z12**2-m_miss1**2);

                    //cout<<"PZ" <<pz1<<endl;     
                    //   TLorentzVector p4miss11(pt_miss1*cos(phi_miss1),pt_miss1*sin(phi_miss1),pz11,sqrt(pz11**2+pt_miss1**2+m_miss1**2));

 

                    p4miss11.SetXYZT(pt_miss1*cos(phi_miss1),pt_miss1*sin(phi_miss1),pz11,sqrt(pow(pz11,2)+pow(pt_miss1,2)+pow(m_miss1,2)));
                    //TLorentzVector p4miss12(pt_miss1*cos(phi_miss1),pt_miss1*sin(phi_miss1),pz12,sqrt(pz12**2+pt_miss1**2+m_miss1**2));
                    p4miss12.SetXYZT(pt_miss1*cos(phi_miss1),pt_miss1*sin(phi_miss1),pz12,sqrt(pow(pz12,2)+pow(pt_miss1,2)+pow(m_miss1,2)));

                    p4miss21.SetXYZT(pt_miss2*cos(phi_miss2),pt_miss2*sin(phi_miss2),pz21,sqrt(pow(pz21,2)+pow(pt_miss2,2)+pow(m_miss2,2)));
                    //TLorentzVector p4miss12(pt_miss1*cos(phi_miss1),pt_miss1*sin(phi_miss1),pz12,sqrt(pz12**2+pt_miss1**2+m_miss1**2));
                    p4miss22.SetXYZT(pt_miss2*cos(phi_miss2),pt_miss2*sin(phi_miss2),pz22,sqrt(pow(pz22,2)+pow(pt_miss2,2)+pow(m_miss2,2)));

                    //p4vis1.SetXYZT(p_vis1*sin(theta_vis1)*cos(phi_vis1),p_vis1*sin(theta_vis1)*sin(phi_vis1), p_vis1*cos(theta_vis1), sqrt(pow(p_vis1,2)+pow(m_vis1,2)));
                    //p4vis2.SetXYZT(p_vis2*sin(theta_vis2)*cos(phi_vis2),p_vis2*sin(theta_vis2)*sin(phi_vis2), p_vis2*cos(theta_vis2), sqrt(pow(p_vis2,2)+pow(m_vis2,2)));


                    mtautau11->Fill((p4vis1+p4miss11 +p4vis2+ p4miss21).M(),DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11))*DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21)));
                    mtautau12->Fill((p4vis1+p4miss11 +p4vis2+ p4miss22).M(),DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11))*DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22)));
                    mtautau21->Fill((p4vis1+p4miss12 +p4vis2+ p4miss21).M(),DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12))*DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21)));
                    mtautau22->Fill((p4vis1+p4miss12 +p4vis2+ p4miss22).M(),DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12))*DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22)));

                    
                    mass_dummy.Fill((p4vis1+p4miss11 +p4vis2+ p4miss21).M(),DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11))*DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21)));
                    //      if ((p4vis1+p4miss11 +p4vis2+ p4miss22).M()>50 && (p4vis1+p4miss11 +p4vis2+ p4miss22).M()<140)
                     mass_dummy.Fill((p4vis1+p4miss11 +p4vis2+ p4miss22).M(),DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11))*DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22)));
                    //      if ((p4vis1+p4miss12 +p4vis2+ p4miss21).M()>50 && (p4vis1+p4miss12 +p4vis2+ p4miss21).M()<140)
                     mass_dummy.Fill((p4vis1+p4miss12 +p4vis2+ p4miss21).M(),DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12))*DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21)));
                    //      if ((p4vis1+p4miss12 +p4vis2+ p4miss22).M()>50 && (p4vis1+p4miss12 +p4vis2+ p4miss22).M()<140)                
                     mass_dummy.Fill((p4vis1+p4miss12 +p4vis2+ p4miss22).M(),DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12))*DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22)));



                    check=true;
                    }

               }
     return check;
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

TF1 TauValidation::DeltaR_had(double pt){



     TF1 cons2("cons2","pol1",1,66.88);

     // cons2.SetChisquare(7.744005);
     // cons2.SetNDF(9);

     cons2.SetParameter(0,-0.04697871);
     //  cons2.SetParError(0,0.06288041);

     cons2.SetParameter(1,0.03170594);
     //  cons2.SetParError(1,0.001579847);


     TF1 MPV("MPV","expo+pol1(2)",1,92.8);

     //  MPV.SetChisquare(36.90687);
     //  MPV.SetNDF(13);

     MPV.SetParameter(0,-0.5269114);
     //  MPV.SetParError(0,0.3510341);

     MPV.SetParameter(1,-0.09367247);
     // MPV.SetParError(1,0.01833037);

     MPV.SetParameter(2,0.112788);
     //  MPV.SetParError(2,0.009585057);

     MPV.SetParameter(3,-0.0008607203);
     //  MPV.SetParError(3,0.0001189079);




     TF1 sigma2("sigma2","expo+pol1(2)",10.72,77.68);

     //  sigma2.SetChisquare(48.44186);
     //  sigma2.SetNDF(10);

     sigma2.SetParameter(0,-2.376518);
     // sigma2.SetParError(0,0.4077938);

     sigma2.SetParameter(1,-0.1253568);
     //  sigma2.SetParError(1,0.01940982);

     sigma2.SetParameter(2,0.00586322);
     //  sigma2.SetParError(2,0.0005587519);

     sigma2.SetParameter(3,-3.839789e-005);
     //   sigma2.SetParError(3,8.373042e-006);


     TF1 total("DeltaR","gaus(0)+landau(3)",0,0.4);
     Double_t par[6];
     par[0]=0;
     par[1]=0;
     par[2]=0;
     par[3]=cons2.Eval(pt);
     par[4]=MPV.Eval(pt);
     par[5]=sigma2.Eval(pt);

     total.SetParameters(par);



     return total; 


}



TF1 TauValidation::DeltaR_lep(double pt){

     // TCanvas *c1=new TCanvas();c1.Divide(2,3); 
     //  c1.cd(1);
     TF1 cons1 ("cons1","pol1",1,109);

     //  cons1.SetChisquare(18.39828);
     //   cons1.SetNDF(16);

     cons1.SetParameter(0,-0.0004049933);
     //  cons1.SetParError(0,0.002755498);

     cons1.SetParameter(1,0.001609134);
     //  cons1.SetParError(1,6.186914e-005);




     TF1 mean("mean","expo(0)+pol1(2)",10.72,85.24);

     //  mean.SetChisquare(10.05);
     //  mean.SetNDF(11);

     mean.SetParameter(0,-1.319604);
     //  mean.SetParError(0,0.2216975);

     mean.SetParameter(1,-0.0698018);
     //  mean.SetParError(1,0.01949947);

     mean.SetParameter(2,0.05926357);
     //   mean.SetParError(2,0.01851204);

     mean.SetParameter(3,-0.0004089469);
     //  mean.SetParError(3,0.0002275216);




     TF1 sigma1("sigma1","expo+pol1(2)",10.72,109);

     // sigma1.SetChisquare(31.74544);
     // sigma1.SetNDF(14);

     sigma1.SetParameter(0,-2.227225);
     //   sigma1.SetParError(0,0.08536175);

     sigma1.SetParameter(1,-0.04167413);
     //  sigma1.SetParError(1,0.009089547);

     sigma1.SetParameter(2,6.679525e-005);
     //  sigma1.SetParError(2,0.01117666);

     sigma1.SetParameter(3,0.0001051946);
     //  sigma1.SetParError(3,0.0001027565);




     TF1 cons2 ("cons2","pol1",1,109);

     //  cons2.SetChisquare(53.00084);
     //  cons2.SetNDF(16);

     cons2.SetParameter(0,-0.03423635);
     //  cons2.SetParError(0,0.02078572);

     cons2.SetParameter(1,0.008789224);
     //  cons2.SetParError(1,0.000580261);




     TF1 MPV("MPV","expo+pol1(2)",10.72,96.04);


     //  MPV.SetChisquare(4.636014);
     //  MPV.SetNDF(13);

     MPV.SetParameter(0,-0.8407024);
     //  MPV.SetParError(0,0.2137112);

     MPV.SetParameter(1,-0.06564579);
     //   MPV.SetParError(1,0.012763);

     MPV.SetParameter(2,0.07128014);
     //  MPV.SetParError(2,0.01614891);

     MPV.SetParameter(3,-0.0004138105);
     //   MPV.SetParError(3,0.0001800071);




     TF1 sigma2("sigma2","expo+pol1(2)",11.8,92.8);

     //  sigma2.SetChisquare(11.9713);
     //  sigma2.SetNDF(13);

     sigma2.SetParameter(0,-2.364371);
     //  sigma2.SetParError(0,0.369229);

     sigma2.SetParameter(1,-0.09803685);
     //  sigma2.SetParError(1,0.01569372);

     sigma2.SetParameter(2,0.01046975);
     //  sigma2.SetParError(2,0.0007834674);

     sigma2.SetParameter(3,-8.072633e-005);
     // sigma2.SetParError(3,9.291709e-006);


     TF1 total ("DeltaR","gaus(0)+landau(3)",0,1);
     Double_t par[6];
     par[0]=cons1.Eval(pt);
     par[1]=mean.Eval(pt);
     par[2]=sigma1.Eval(pt);
     par[3]=cons2.Eval(pt);
     par[4]=MPV.Eval(pt);
     par[5]=sigma2.Eval(pt);

     total.SetParameters(par);





     return total;



}

