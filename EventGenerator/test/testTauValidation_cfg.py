import FWCore.ParameterSet.Config as cms

process = cms.Process("TauTest")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Validation.EventGenerator.TauValidation_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
#    input = cms.untracked.int32(-1)    
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

	'file:/afs/cern.ch/user/h/hmermerk/cmssw/CMSSW_7_3_1_patch1/src/GeneratorInterface/TauolaInterface/test/DYtoLL_mad_pyth_tauola_13TEV.root'
#	'/store/relval/CMSSW_3_6_0_pre5/RelValHiggs200ChargedTaus/GEN-SIM-DIGI-RAW-HLTDEBUG/START36_V3-v1/0010/E09B3F44-E13D-DF11-995D-00304867BFA8.root'
#	'rfio:/castor/cern.ch/user/s/slehti/HiggsAnalysisData/test_H120_100_1_08t_RAW_RECO.root'
#	'rfio:/castor/cern.ch/user/s/slehti/testData/DYToTauTau_M_20_TuneZ2_7TeV_pythia6_tauola_Fall10_START38_V12_v2_GEN_SIM_RECO_100ev.root'
    ) 
)



#Add these 3 lines to put back the summary for timing information at the end of the logfile
#(needed for TimeReport report)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )    
process.TFileService = cms.Service("TFileService",
          fileName = cms.string('histo.root')
          )


process.p = cms.Path(process.tauValidation)




