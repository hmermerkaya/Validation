import FWCore.ParameterSet.Config as cms

process = cms.Process("TauTest")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Validation.EventGenerator.TauValidation_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(100)



process.source = cms.Source("PoolSource",
          fileNames = cms.untracked.vstring(

          #     'file:/afs/cern.ch/user/h/hmermerk/cmssw/CMSSW_7_3_1_patch1/src/GeneratorInterface/TauolaInterface/test/DYtoLL_mad_pyth_tauola_13TEV.root'
#	'/store/relval/CMSSW_3_6_0_pre5/RelValHiggs200ChargedTaus/GEN-SIM-DIGI-RAW-HLTDEBUG/START36_V3-v1/0010/E09B3F44-E13D-DF11-995D-00304867BFA8.root'
#	'rfio:/castor/cern.ch/user/s/slehti/HiggsAnalysisData/test_H120_100_1_08t_RAW_RECO.root'
#	'rfio:/castor/cern.ch/user/s/slehti/testData/DYToTauTau_M_20_TuneZ2_7TeV_pythia6_tauola_Fall10_START38_V12_v2_GEN_SIM_RECO_100ev.root'
              # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/020C52E0-1B8B-E411-AE87-002590D0AFE0.root',
              # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/06E9B1AC-1D8B-E411-975F-00259073E4C2.root',
              # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/0A691D86-278B-E411-A8C3-00259073E4FC.root',
              # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/0CCCDE6E-1D8B-E411-A38C-20CF305B059C.root',
              # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/0EA41EC9-2B8B-E411-A855-20CF3019DEED.root',
              # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/102BD642-218B-E411-A575-20CF305616CC.root',
               # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/1048D9B4-258B-E411-B7E2-E0CB4E19F983.root',
               # '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/10D515FD-2C8B-E411-AE01-20CF305B0508.root'
#               '/store/mc/Summer12_DR53X/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-tauola-pythia6/AODSIM/PU_S8_START53_V7A-v1/00000/00744031-9307-E211-B238-0026189438B9.root',
 #              '/store/mc/Summer12_DR53X/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-tauola-pythia6/AODSIM/PU_S8_START53_V7A-v1/00000/00F45E02-A107-E211-A107-003048FFCBB0.root',
 #              '/store/mc/Summer12_DR53X/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-tauola-pythia6/AODSIM/PU_S8_START53_V7A-v1/00000/02A215CA-9F07-E211-8402-001A92971B36.root'
    ) 
)


#process.load("dataset_DYtoJets_8TeV_madgraph_tarball")
process.load("dataset_DYtoTauTau_8TeV_pythia6_tauola")
#process.load("dataset_DYtoTauTau_8TeV_pythia6_v2")
#process.load("dataset_DYtoJets_13TeV_madgraph_pythia8")
#process.load("dataset_DYtoJets_13TeV_madgraph_pythia8_tauola")

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(10000)
     #  input = cms.untracked.int32(-1)    
     )




#Add these 3 lines to put back the summary for timing information at the end of the logfile
#(needed for TimeReport report)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )    
process.TFileService = cms.Service("TFileService",
          fileName = cms.string('dataset_DYtoTauTau_8TeV_pythia6_tauola.root')
          )


process.p = cms.Path(process.tauValidation)




