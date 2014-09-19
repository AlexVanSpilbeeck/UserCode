import FWCore.ParameterSet.Config as cms

process = cms.Process("SIPIXELDQM")
process.load("Geometry.TrackerSimData.trackerSimGeometryXML_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")


process.load("EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi")
process.siPixelDigis.InputLabel = 'source'
process.siPixelDigis.IncludeErrors = True

process.load("RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi")

process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")

process.load("CalibTracker.SiPixelTools.SiPixelErrorsCalibDigis_cfi")
process.siPixelErrorsDigisToCalibDigis.saveFile=False
process.siPixelErrorsDigisToCalibDigis.SiPixelProducerLabelTag = 'siPixelCalibDigis'
process.load("CalibTracker.SiPixelGainCalibration.SiPixelCalibDigiProducer_cfi")
process.load("CalibTracker.SiPixelSCurveCalibration.SiPixelSCurveCalibrationAnalysis_cfi")
process.siPixelSCurveAnalysis.DetSetVectorSiPixelCalibDigiTag = 'siPixelCalibDigis'
process.siPixelSCurveAnalysis.saveFile = False
process.load("CalibTracker.SiPixelGainCalibration.SiPixelGainCalibrationAnalysis_cfi")
process.siPixelGainCalibrationAnalysis.DetSetVectorSiPixelCalibDigiTag = 'siPixelCalibDigis'
process.siPixelGainCalibrationAnalysis.saveFile = False
process.siPixelGainCalibrationAnalysis.savePixelLevelHists = True 
process.siPixelGainCalibrationAnalysis.vcalHighToLowConversionFac=7
process.siPixelGainCalibrationAnalysis.useVCALHIGH = True 

process.load("DQMServices.Core.DQM_cfg")

process.load("DQMServices.Components.DQMEnvironment_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_31X_GLOBALTAG"
#process.GlobalTag.globaltag = "GR_P_V17::All"
process.GlobalTag.globaltag = "GR_R_71_V7::All"

from CondCore.DBCommon.CondDBCommon_cfi import *
process.siPixelCalibGlobalTag =  cms.ESSource("PoolDBESSource",
                                              CondDBCommon,
                              #               connect = cms.string("frontier://FrontierProd/CMS_COND_31X_PIXEL"),
                                              
                                              toGet = cms.VPSet(
                                                cms.PSet(record = cms.string('SiPixelCalibConfigurationRcd'),
                                                tag = cms.string('GLOBALCALIB_default'))
                                                )
                                              )
process.siPixelCalibGlobalTag.connect = "frontier://FrontierProd/CMS_COND_31X_PIXEL"
process.es_prefer_dbcalib = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.source = cms.Source("SOURCETYPE",
#    debugFlag = cms.untracked.bool(True),
#    debugVebosity = cms.untracked.uint32(1),
    ONEPARAM
    TWOPARAM
    fileNames = cms.untracked.vstring(
    'FILENAME'
       				      )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service("MessageLogger",
    TEXTFILE = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    ),
    destinations = cms.untracked.vstring('TEXTFILE')
)

process.AdaptorConfig = cms.Service("AdaptorConfig")

process.sipixelEDAClient = cms.EDAnalyzer("SiPixelEDAClient",
    EventOffsetForInit = cms.untracked.int32(10),
    ActionOnLumiSection = cms.untracked.bool(False),
    ActionOnRunEnd = cms.untracked.bool(True),
    HighResolutionOccupancy = cms.untracked.bool(True),
    NoiseRateCutValue = cms.untracked.double(-1.)
    
)

process.qTester = cms.EDAnalyzer("QualityTester",
    qtList = cms.untracked.FileInPath('DQM/SiPixelMonitorClient/test/sipixel_qualitytest_config.xml'),
    QualityTestPrescaler = cms.untracked.int32(1),
    getQualityTestsFromFile = cms.untracked.bool(True)
)

process.Reco = cms.Sequence(process.siPixelDigis)
process.Calibration = cms.Sequence(process.siPixelCalibDigis*process.siPixelGainCalibrationAnalysis)
process.DQMmodules = cms.Sequence(process.qTester*process.dqmEnv*process.dqmSaver)
process.p = cms.Path(process.Reco*process.qTester*process.dqmEnv*process.Calibration*process.sipixelEDAClient*process.dqmSaver)

# cms.Path(process.Reco*process.DQMmodules*process.Calibration*process.RAWmonitor*process.DIGImonitprocess.sipixelEDAClient)
#process.DQM.collectorHost = ''
process.dqmSaver.convention = 'Online'
process.dqmSaver.producer = 'DQM'
process.dqmEnv.subSystemFolder = 'Pixel'
process.dqmSaver.dirName = '.'
process.dqmSaver.saveByLumiSection = -1
process.dqmSaver.saveByRun = 1
process.dqmSaver.saveAtJobEnd = True

