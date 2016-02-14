from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'WJets100To200_weight-v'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt','Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt','Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt','Summer15_25nsV2_MC_L1FastJet_AK8PFchs.txt','Summer15_25nsV2_MC_L2Relative_AK8PFchs.txt','Summer15_25nsV2_MC_L3Absolute_AK8PFchs.txt']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =15 
config.Data.totalUnits = -1
config.Data.publication = False

# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'WJets100To200_weight'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
