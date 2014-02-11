import FWCore.ParameterSet.Config as cms

process = cms.Process("q2weights")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12_DR53X/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/0000/A002DE4D-EE04-E211-A826-003048C68A9E.root'
    ),
     skipEvents = cms.untracked.uint32(0)
)

process.q2weights = cms.EDProducer('Q2WeightsWjets',
    fileWeights=cms.FileInPath("QScaleISRWeightProducer/Q2Weights_pfeiffer/Wjets_Q2_weights_Jet1_Qscale_pdf.root"),
    fileWeightUp=cms.FileInPath("QScaleISRWeightProducer/Q2Weights_pfeiffer/simuUp_8var_pdf.root"),
    fileWeightDown=cms.FileInPath("QScaleISRWeightProducer/Q2Weights_pfeiffer/simuDown_8var_pdf.root"),
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

