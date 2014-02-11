import FWCore.ParameterSet.Config as cms

process = cms.Process("Q2WEIGHTS")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10000/000560C1-FD97-E211-9F33-00304867924E.root"
    )
)

process.prunedGenParticlesForISRWeights = cms.EDProducer("GenParticlePruner",
    src=cms.InputTag("genParticles"),
    select=cms.vstring(
        "drop  *",
        "keep status = 3", #keeps all particles from the hard matrix element
        "+keep abs(pdgId) = 15 & status = 1" #keeps intermediate decaying tau
        )
    )
process.qScaleISRWeights = cms.EDProducer('QScaleISRWeightProducer',
    src=cms.InputTag("prunedGenParticlesForISRWeights")
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.prunedGenParticlesForISRWeights*process.qScaleISRWeights)

process.e = cms.EndPath(process.out)
