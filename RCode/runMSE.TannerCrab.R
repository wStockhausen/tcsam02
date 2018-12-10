#Example R script to run the TCSAM02 assessment model in MSE mode. 
#  (NOTE:this requires the rTCSAM02 and related R packages be installed)
#To run:
#  1. Copy this file and the folders "InputFiles.2018AM17","InputFiles.OpModMode", and "InputFiles.EstModMode" 
#      to a folder outside the tcsam02 project (so the runs will NOT be part of the project's git version control).
#  2. Open the copied version of this file in RStudio and change the working directory to be the same as that file.
#  3. Modify the "path2tcsam02" variable below so it agrees with the path to the tcsam02 program
#  4. "Source" this file
#The output will be 5 runs of 10 years each of the operating model and estimation model in folders under "testMSE"
#in the working directory this file was "sourced" in.
require(rTCSAM02);
HCR<-1;
path2tcsam02<-'~/StockAssessments-Crab/AssessmentModelDevelopment/tcsam02/dist/MacOSX-ADMBv12/GNU-MacOSX';
datasets<-list(file="Datasets.inp",
               components=list(
                       Bio="Info",
                       Fishery=c("TCF","SCF","GTF","RKF"),
                       Survey="NMFS",
                       Growth="EBS",
                       ChelaHeights=NULL)
               );
runMSE(os='osx',
       topLevelFolder='./testMSE',
       model='tcsam02',
       path2model=path2tcsam02,
       HCR=HCR,
       numRuns=5,
       firstYr=2018,
       numYrs=10,
       runBaseModel=TRUE,
       baseModelInfo=list(path="./InputFiles.2018AM17",
                          configFile="Configuration.inp",
                          datasets=datasets,
                          pin=TRUE,
                          pinFile="tcsam02.pin",
                          minPhase=5,
                          calcTAC=TRUE,
                          HCR=HCR),
       opModInfo=list(path="./InputFiles.OpModMode",
                       configFile="OpMod.Configuration.inp",
                       optsFile="OpMod.Options.inp",
                       mpiFile="OpMod.ParametersInfo.inp",
                       pinFile="OpModPinFile.pin"),
       estModInfo=list(path="./InputFiles.EstModMode",
                       configFile="EstMod.Configuration.inp",
                       optsFile="EstMod.Options.inp",
                       minPhase=5),
       keepFiles=c("tmp.sh","tcsam02.par"),
       cleanupAll=FALSE,
       test=FALSE,
       verbose=TRUE);
