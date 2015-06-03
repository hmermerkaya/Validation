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

class TauValidation : public DQMEDAnalyzer
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
	virtual void bookHistograms(DQMStore::IBooker &i, edm::Run const &, edm::EventSetup const &) override;
	virtual void dqmBeginRun(const edm::Run& r, const edm::EventSetup& c) override;
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
	void spinEffectsWHpm(const reco::GenParticle*,int,int,std::vector<const reco::GenParticle*> &part,double weight);
	void spinEffectsZH(const reco::GenParticle* boson, double weight);
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
  
        MonitorElement *nTaus, *nPrimeTaus;
  	MonitorElement *TauPt, *TauEta, *TauPhi, *TauProngs, *TauDecayChannels, *TauMothers, 
	  *TauSpinEffectsW_X, *TauSpinEffectsW_UpsilonRho, *TauSpinEffectsW_UpsilonA1,*TauSpinEffectsW_eX,*TauSpinEffectsW_muX,
	  *TauSpinEffectsHpm_X, *TauSpinEffectsHpm_UpsilonRho, *TauSpinEffectsHpm_UpsilonA1,*TauSpinEffectsHpm_eX,*TauSpinEffectsHpm_muX, 
	  *TauSpinEffectsZ_MVis, *TauSpinEffectsZ_Zs, *TauSpinEffectsZ_Xf, *TauSpinEffectsZ_Xb,
	  *TauSpinEffectsZ_X50to75,*TauSpinEffectsZ_X75to88,*TauSpinEffectsZ_X88to100,*TauSpinEffectsZ_X100to120,*TauSpinEffectsZ_X120UP,
	  *TauSpinEffectsZ_eX, *TauSpinEffectsZ_muX,  *TauSpinEffectsH_X,
	  *TauSpinEffectsH_MVis, *TauSpinEffectsH_Zs, *TauSpinEffectsH_Xf, *TauSpinEffectsH_Xb,
	  *TauSpinEffectsH_eX, *TauSpinEffectsH_muX, *TauSpinEffectsH_rhorhoAcoplanarityplus,  *TauSpinEffectsH_rhorhoAcoplanarityminus,
	  *TauBremPhotonsN,*TauBremPhotonsPt,*TauBremPhotonsPtSum,*TauFSRPhotonsN,*TauFSRPhotonsPt,*TauFSRPhotonsPtSum,
	  *TauSpinEffectsH_pipiAcoplanarity,*TauSpinEffectsH_pipiAcollinearity,*TauSpinEffectsH_pipiAcollinearityzoom, *DecayLength,
	  *LifeTime,
*TauSpinEffectsZhplus_X50to75,*TauSpinEffectsZhplus_nospin_X50to75,*TauSpinEffectsZhplus_Xnew50to75,*TauSpinEffectsZhplus_nospin_Xnew50to75,*TauSpinEffectsZhplus_X75to88,*TauSpinEffectsZhplus_nospin_X75to88,*TauSpinEffectsZhplus_Xnew75to88,*TauSpinEffectsZhplus_nospin_Xnew75to88,*TauSpinEffectsZhplus_X88to100,*TauSpinEffectsZhplus_nospin_X88to100,*TauSpinEffectsZhplus_Xnew88to100,*TauSpinEffectsZhplus_nospin_Xnew88to100,*TauSpinEffectsZhplus_X100to120,*TauSpinEffectsZhplus_nospin_X100to120,*TauSpinEffectsZhplus_Xnew100to120,*TauSpinEffectsZhplus_nospin_Xnew100to120,*TauSpinEffectsZhplus_X120UP,*TauSpinEffectsZhplus_nospin_X120UP,*TauSpinEffectsZhplus_Xnew120UP,*TauSpinEffectsZhplus_nospin_Xnew120UP,
*TauSpinEffectsZhplus_X,
*TauSpinEffectsZhplus_nospin_X,
*TauSpinEffectsZhplus_Xnew,
*TauSpinEffectsZhplus_nospin_Xnew,
*TauSpinEffectsZhminus_X,
*TauSpinEffectsZhminus_nospin_X,
*TauSpinEffectsZhminus_Xnew,
*TauSpinEffectsZhminus_nospin_Xnew,
*TauSpinEffectsZ_unspin_Xnew,
*TauSpinEffectsZ_unspin_X,
*TauSpinEffectsZhminus_X50to75,*TauSpinEffectsZhminus_nospin_X50to75,*TauSpinEffectsZhminus_Xnew50to75,*TauSpinEffectsZhminus_nospin_Xnew50to75,*TauSpinEffectsZhminus_X75to88,*TauSpinEffectsZhminus_nospin_X75to88,*TauSpinEffectsZhminus_Xnew75to88,*TauSpinEffectsZhminus_nospin_Xnew75to88,*TauSpinEffectsZhminus_X88to100,*TauSpinEffectsZhminus_nospin_X88to100,*TauSpinEffectsZhminus_Xnew88to100,*TauSpinEffectsZhminus_nospin_Xnew88to100,*TauSpinEffectsZhminus_X100to120,*TauSpinEffectsZhminus_nospin_X100to120,*TauSpinEffectsZhminus_Xnew100to120,*TauSpinEffectsZhminus_nospin_Xnew100to120,*TauSpinEffectsZhminus_X120UP,*TauSpinEffectsZhminus_nospin_X120UP,*TauSpinEffectsZhminus_Xnew120UP,*TauSpinEffectsZhminus_nospin_Xnew120UP,*TauSpinEffectsZ_X, *TauSpinEffectsZ_nospin_X, *TauSpinEffectsZ_Xnew,  *TauSpinEffectsZ_nospin_Xnew        
           ;

	unsigned int NMODEID;
	MonitorElement *MODEID;
	std::vector<std::vector<MonitorElement *> > MODEInvMass;
double weightplus,weightminus;
bool weightvalid;
	int zsbins;
	double zsmin,zsmax;

	edm::EDGetTokenT<reco::GenParticleCollection> genparticleCollectionToken_;
	edm::EDGetTokenT<edm::HepMCProduct> hepmcCollectionToken_;
};

#endif

