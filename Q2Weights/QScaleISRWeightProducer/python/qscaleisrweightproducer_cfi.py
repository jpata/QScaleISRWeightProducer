import FWCore.ParameterSet.Config as cms
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
