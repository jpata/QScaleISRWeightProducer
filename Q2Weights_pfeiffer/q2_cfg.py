import FWCore.ParameterSet.Config as cms

process = cms.Process("q2weights")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/scratch/hh/lustre/cms/user/peiffer/WjetsTest.root'
    ),
     skipEvents = cms.untracked.uint32(0)
)

process.q2weights = cms.EDProducer('Q2Weights'
)

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('tree.root')
#)

process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('cmssw-out.root'),
#        SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
        outputCommands = cms.untracked.vstring('keep *')
      )


process.p = cms.Path(process.q2weights)
process.outpath = cms.EndPath( process.out )

