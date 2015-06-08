#ifndef TauValidation_H
#define TauValidation_H

// framework & common header files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

//DQM services
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "TLorentzVector.h"

#include "Validation/EventGenerator/interface/WeightManager.h"

class TauValidation :  public edm::EDAnalyzer 
{
    public:
	// tau decays
        enum  {undetermined,
               electron,
               muon,
               pi,
               rho,
	       a1,
               K,
	       Kstar,
	       pi1pi0,
               pinpi0,
               tripi,
               tripinpi0,
	       stable};
	// tau mother particles 
	enum  {other,
	       B,
	       D,
	       gamma,
	       Z,
	       W,
	       HSM,
	       H0,
	       A0,
	       Hpm};

    public:
	explicit TauValidation(const edm::ParameterSet&);
	virtual ~TauValidation();
	virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
    private:
//	  WeightManager wmanager_;
          void GetLastSelf(const reco::GenParticle *Particle);
           std::vector<const reco::GenParticle *> Particles;
     const reco::GenParticle* GetLastSelf2(const reco::GenParticle* particle);
	int tauMother(const reco::GenParticle*, double weight);
	int tauProngs(const reco::GenParticle*, double weight);
	int tauDecayChannel(const reco::GenParticle* tau,int jak_id,unsigned int TauBitMask, double weight);
	int findMother(const reco::GenParticle*);
	bool isLastTauinChain(const reco::GenParticle* tau);
	void mmc_method_interface(const reco::GenParticle* boson, double weight);
bool mmc_method(const TLorentzVector &p4vis1,const TLorentzVector &p4miss1,const TLorentzVector &p4vis2,const TLorentzVector &p4mis2, TH1F &mass_dummy);
	TF1 DeltaR_had(double pt);
          TF1 DeltaR_lep(double pt);
          
          double leadingPionMomentum(const reco::GenParticle*, double weight);
	double visibleTauEnergy(const reco::GenParticle*);
	TLorentzVector leadingPionP4(const reco::GenParticle*);
	TLorentzVector motherP4(const reco::GenParticle*);
	void photons(const reco::GenParticle*, double weight);
	void findTauList(const reco::GenParticle* tau,std::vector<const reco::GenParticle*> &TauList);
	void findFSRandBrem(const reco::GenParticle* p, bool doBrem, std::vector<const reco::GenParticle*> &ListofFSR,
			   std::vector<const reco::GenParticle*> &ListofBrem);
	void FindPhotosFSR(const reco::GenParticle* p,std::vector<const reco::GenParticle*> &ListofFSR,double &BosonScale);
	const reco::GenParticle* GetMother(const reco::GenParticle* tau);
	const std::vector<const reco::GenParticle*> GetMothers(const reco::GenParticle* boson);
	double Zstoa(double zs);
	void countParticles(const reco::GenParticle* p,int &allCount, int &eCount, int &muCount,
			    int &pi0Count,int &piCount,int &rhoCount,int &a1Count,int &KCount,int &KstarCount);

    	edm::InputTag genparticleCollection_;
	edm::InputTag hepmcCollection_;

  	/// PDT table
  	edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable ;
         TH1F *mass_Z; 
         TH1F *mtautau;
         TH1F *mtautau11;
         TH1F *mtautau12;
         TH1F *mtautau21;
         TH1F *mtautau22;
         
	unsigned int NMODEID;
	int zsbins;
	double zsmin,zsmax;

	edm::EDGetTokenT<reco::GenParticleCollection> genparticleCollectionToken_;
	edm::EDGetTokenT<edm::HepMCProduct> hepmcCollectionToken_;
};

#endif

