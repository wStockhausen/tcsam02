//TCSAM02: Bering Sea Tanner crab model
//
// Author: William Stockhausen (william.stockhausen@noaa.gov)
//
// Model units (unless otherwise noted):
//  Individual crab weights are in KG (kilograms).
//  Abundance (numbers) is in MILLIONS (10^6's) of crabs.
//  Biomass (weight) is in 1000's MT (metric tons).
//
//  mxYr denotes the final fishery year (July, mxYr->June, mxYr+1) prior to the assessment.
//  Typically, the final survey is conducted July, mxYr+1.
//  The assessment is conducted for mxYr+1. 
//  Final population abundance is estimated for July 1, mxYr+1.
//
// History:
//  2014-02-11: created
//  2014-05-23: renamed TCSAM2015
//  2014-06-05: 1. moved setDevs() to ) to tcsam::setDevs() in ModelParameterFunctions.hpp/cpp
//              2. moved calcPriors(...) to tcsam::calcPriors(...) in ModelParameterFunctions.hpp/cpp
//  2014-06-09: started review to make sure deep copies are implemented for objects that 
//              should not be modified within a computational unit.
//  2014-06-11: 1. created setFinalVal(...), setFinalVals(...) functionality for ModelParameterInfoTypes
//                 and revised writeToR(...) for parameters to output both initial and final values.
//  2014-06-16: 1. Added sdrLnR_y, sdrSpB_xy, sdrFinlPop_xmsz and MCMC calc.s
//  2014-09-30: 1. Changed to MALE=1, FEMALE=2 ordering of sexes
//              2. Added constant "BioData->recLag" for recruitment lag (set in BioData file)
//  2014-10-30: 1. Corrected sex, maturity, and shell condition dimensioning and loops 
//                  to be independent of ordering of MALE/FEMALE, NEW/OLD_SHELL, and IMMATURE/MATURE
//                  values.
//              2. Added logic to avoid FEMALE-specific logic if model includes only MALEs.
//              3. Updated to use tcsam::ALL_SXs, tcsam::ALL_SCs, tcsam::ALL_MSs
//  2014-11-18: 1. Corrected output of simulated data for retrospective model runs
//              2. Added "-doRetro iRetro" command line option to facilitate retrospective model runs
//              3. Corrected IndexRange behavior for retrospective runs
//              4. MAKE SURE max IndexRanges for TRAWL SURVEYS are set to "-2" (max year index=mxYr+1) 
//                  for correct behavior in retrospective runs
//  2014-11-21: 1. Added ability to simulate data and fit it within a single model run.
//              2. Added "fitSimData" command line option to invoke 1.
//  2014-11-24: 1. Added option to add stochasticity to simulated data fit within the model.
//              2. Added R output for objective function components
//  2014-12-02: 1. Enabled switch to use pin file
//              2. Added 'seed' as input flag to change rng seed (iSeed) for jittering/resampling
//              3. Enabled jitter functionality by modifying setInitVals(...) functions
//  2014-12-08: 1. Added several parameterizations for logistic, double logistic selectivity functions.
//              2. Changed IndexRange to interpret a single input value as a max, not min, index.
//              3. Modified calcNLLs_CatchNatZ to handle more fitting options and not return a value.
//              4. Modified calcMultinomialNLL to write obs, mod values to R list.
//              5. Changed input format for BioData to read in "typical" and "atypical" values
//                  for fishing season midpoint and mating time (the latter by year).
//              6. Changed input format for AggregateCatchData to "long format" from "wide format"
//                  so data is input by sex, maturity, shell_condition factor combinations with rows having
//                      year, value, cv
//  2014-12-22: 1. Finished 6. from above.
//              2. Modified other calcNLLs to write obs, mod, std values to R list.
//              3. Upped version to 01.02 to reflect format changes.
//  2015-02-18: 1. Changed NLL calculations for normal, lognormal error structures
//                  to depend only on the zscores (so they don't include the log(sigma) terms)
//                  so that NLLs are non-negative.
//              2. Added asclogistic50Ln95 selectivity function.
//              3. Implemented ability to read initial value vectors for VectorVector and DevsVector
//                  parameter info classes.
//              4. Revised BioData class and input data file to simplify input.
//              5. Revised aggregateData() and replaceCatchData() functions in AggregateCatchData class.
//              6. Added units_KMT to units_KMT conversion factor (=1) to getConversionFactor(...) functions
//              7. Incremented version to reflect changes.
//  2015-02-19: 1. Added initial size comps calculations based on equilibrium size calcs
//              2. Added gmacs growth function
//  2015-02-26: 1. Added options to use/not use gmacs growth and initial size comps
//              2. Revised model outputs to R to facilitate comparison with rsimTCSAM results
//              3. Fixed problem with fsZ for selectivity functions constrained to be an integer
//  2015-03-01: 1. Finished(?) revisions to R output.
//  2015-03-02: 1. Changed R format for model dimensions (sex, maturity, etc) so all are lower case and "_"s are replaced by spaces
//              2. Changed R format for survey, fishery names so all "_"s are replaced by spaces
//  2015-03-05: 1. Dropped 'm' dimension from prGr_yxmszz (now prGr_yxszz) to match rsimTCSAM (model assumes same growth for all molts)
//              2. Dropped output of summary abundance arrays to decrease size of rep file.
//  2015-03-24: 1. Added "writeParameters" function to write parameter info to csv 
//  2015-04-01: 1. Assigned handling mortality by fishery, year to hmF_fy for output.
//  2015-04-10: 1. Added penalties for final value in devs vectors to achieve consistency with bounds
//              2. Added relative path functionality to data files
//              3. Added debugModelOptions
//  2015-04-20: 1. Added checks on fishery and survey names in DATA_SECTION.
//  2015-04-28: 1. Fixed normalization for model size comps using BY_XMS option: it had been summing
//                 over shell condition as well as size. SHould only have summed over size (now corrected).)
//  2015-05-12: 1. Added fishery effort extrapolation calculations to R output (mp/Eff_list)
//              2. Increased stringency of convergence criteria
//  2015-05-13: 1. Added phase-specific penalties on F-devs that reduce to 0 at some phase
//              2. Added dbgNLLs as debug value. Must consider using std::bitset to implement debug flags
//              3. Incremented internal version to 01.60.
//              4. Can now specify penalty weight and starting phase for last devs to ModelOptions file
//  2015-05-18: 1. Added FIT_BY_XS, FIT_BY_X_ME, FIT_BY_X_SE, FIT_BY_XM_SE options to fit size comps
//              2. Cleaned up application of pS?Devs to pS? so this only happens for required years
//              3. Added mean growth info (mnGrZ_cz, mnGrZ_yxsz) to model rep file
//              4. Reconfigured growth-related model output to mp$T_list$... to match rsimTCSAM
//  2016-01-21: 1. Combined SurveyData, FisheryData classes into a single class: FleetData.
//                  This standardizes data input files, R output lists, and simplifies
//                  adding additional data types (e.g., tagging data, chela height data).
//              2. Updated model version (modVer) to 2016.00. Will create git tag of same to track development.
//  2016-02-29: 1. Model version now YYYY.MM.DD.
//              2. Started to add OFL calculations
//  2016-03-15: 1. Model version incremented.
//              2. OFL calculations now operational
//              3. added debugOFL as command line flag
//  2016-03-30: 1. added csv output file of ALL parameters in last phase Report section
//              2. updated version number. 
//  2016-03-31: 1. Revised NM calc.s so pLnM represents natural mortality rates
//                  on immature male crab (was mature male crab, previously). pLnDMM
//                  and pLnDMXM now reprsent corresponding offsets for mature, not
//                  immature, crab.
//              2. updated version number. 
//              3. added command line flag (-calcOFL) to enable OFL calculations
//              4. now writing initial parameter values to csv file (TCSAM2015.init_params.csv)
//  2016-04-04: 1. writing 2 sets of initial, final parameter values to csv files as
//                  TCSAM2015.params.XXX.init.csv and TCSAM2015.params.XXX.final.csv,
//                  where XXX = 'active' and 'all'.
//              2. Added write to 'jitterInfo.csv' if jittering.
//              3. Changed command line option "-seed" to "-iSeed"(similar to 2013 ver))
//  2016-04-05: 1. Changed exp() to mfexp().
//              2. Incremented version to 2016.04.05.
//              3. Increased convergence criteria strictness in phases 5+.
//  2016-04-06: 1. Added testNaNs() function to identify nan's in calculations.
//              2. Incremented version.
//              3. Added "0" option to skip effort averaging, even if 
//                  effort data is available. Added avg F/Eff to R output.
//  2016-04-11: 1. Added additional cout's to preliminary calcs.
//  2016-04-12: 1. Added additional output in testNaNs().
//              2. Revised penalty function for non-decreasing maturity parameters.
//                  posfun() was generating very large penalties.
//              3. Added option to select penalty function in 2.
//              4. Added VERSION strings to ModelConstants, ModelConfiguration, ModelOptions.
//  2016-04-15: 1. Added calcNoneNLL functions to standardize output when NLL choice is NONE.
//              2. Removed output of T_yxszz in ReportToR_ModelProcesses(). TODO: add option to output.
//              3. Added ctrPrcCallsInPhase to track procedure calls w/in each phase.
//  2016-04-22: 1. Incremented version.
//              2. Recompiled to use logPDF_expnormal() function as priors for
//                  growth parameters.
//              3. revised calcNatMort() to correct error in sequence of adding
//                  parameters to lnM.
//  2016-04-23: 1. Incremented version.
//              2. Added added more diagnostics to deal with exp(prevariable(nan)) problems,
//                  which seem to be linked to asclogistic5095 parameterization.
//  2016-04-25: 1. Added tcsam::FIT_BY_X_MATONLY option to fit aggregated catch data
//              2. Reconfigured fishery and survey capture rate calculations 
//                  somewhat to reflect current arrangement with SEX as an index
//                  in the parameter group (so pLnCDX does not really refer to a
//                  female-specific offset, but is more general).
//  2016-05-18: 1. Changed prMolt2Mat_.. to prM2M_.. in output to R to match R packages.
//  2016-06-20: 1. Mainly cosmetic changes renaming some variables. Revised some output 
//                  to facilitate comparison with TCSAM013 model output.
//  2016-09-30: 1. Renamed project/model tcsam02. 
//              2. renamed input/output files from "TCSAM2015" to "tcsam02".
//--2016-10-13: 1. writing population, survey, fishery numbers/biomass-at-size arrays to rep file
//              2. Incremented version to 20161013.
//--2016-11-06: 1. Refactored OFL calculations.
//              2. Added OFL output to MCMC results.
//--2016-11-08: 1. Revised PopDy and OFL calculations such that z in prM2M refers 
//                  to post-molt size, not pre-molt size.
//--2016-11-12: 1. Revised CatchData types to include weighting factor for
//                  likelihood components (llWgt)
//              2. Revised EffortData to use an IndexBlock to specify ranges
//                  over which to calculate average effort/average F.
//              3. Revised all NLL functions for data to incorporate an input weighting factor.
//--2016-11-14: 1. Added GrowthData and ChelaHeightData classes to handle new
//                  growth and chela height data inputs.
//              2. Modified ModelDatasets to handle new Data classes.
//              3. Removed explicit call to superclass destructor in Tier3_Calculator.
//                  It's called implicitly and resulted in compilation errors on Windows.
//--2016-11-15: 1. Added asclogisticLn50 and dbllogisticLnD50 functions to SelFcns.
//              2. Updated model, model configuration and model options versions
//                  to "2016.11.15".
//--2016-11-16: 1. Revised TCSAM2013 growth calc to use wts::log_gamma_density functions
//                  to try to eliminate NaNs when using optGrowth=0.
//              2. Added output to "GrowthReport.dat" when NaNs detected in growth function.
//              3. Reverted parameterization of growth transition matrix to arithmetic scale.
//              4. Revised growth calc for optGrowth=1 (using cumd_gamma).
//--2016-11-18: 1. Turned off calcNLLs_Recruitment() by setting llWgt = 0 in function
//                  to agree more closely with TCSAM2013 approach. llWgt should be an
//                  input to the ModelControl file.
//              2. Implemented calcNLLs_GrowthData().
//--2016-11-19: 1. Fixed mnZ_n calc.s in calNLLs_GrowthData.
//--2016-11-21: 1. Fixed indexing/dimension problem with grA_xy, grB_xy, and grBeta_xy.
//              2. Fixed problems with writing calcNLLs_GrowthData results to R.
//--2016-11-22: 1. Expanded output to R in calcNLLs_GrowthData.
//--2016-12-01: 1. Fixed problem where sel & ret functions were not being assigned
//                  to output arrays for years where effort extrapolation was used.
//--2016-12-05: 1. Added additional cases and info to exiting from calcGrowth with
//                  problems.
//
// =============================================================================
// =============================================================================
GLOBALS_SECTION
    #include <math.h>
    #include <time.h>
    #include <limits>
    #include <admodel.h>
    #include "TCSAM.hpp"

    adstring model  = tcsam::MODEL;
    adstring modVer = tcsam::VERSION; 
    
    time_t start,finish;
    
    //model objects
    ModelConfiguration*  ptrMC; //ptr to model configuration object
    ModelParametersInfo* ptrMPI;//ptr to model parameters info object
    ModelOptions*        ptrMOs;//ptr to model options object
    ModelDatasets*       ptrMDS;//ptr to model datasets object
    ModelDatasets*       ptrSimMDS;//ptr to simulated model datasets object
    OFLResults           oflResults;//OFL results object for MCMC calculations
        
    //dimensions for R output
    adstring yDms;
    adstring xDms;
    adstring mDms;
    adstring sDms;
    adstring fDms;
    adstring vDms;
    adstring ypDms;
    adstring zbDms;
    adstring zpDms;
    adstring zcDms;
    
    //file streams and filenames
    std::ofstream mcmc;        //stream for mcmc output
    
    //filenames
    adstring fnMCMC = "tcsam02.MCMC.R";
    adstring fnConfigFile;//configuration file
    adstring fnResultFile;//results file name 
    adstring fnPin;       //pin file
    
    //runtime flags (0=false)
    int jitter     = 0;//use jittering for initial parameter values
    int resample   = 0;//use resampling for initial parameter values
    int opModMode  = 0;//run as operating model, no fitting
    int usePin     = 0;//flag to initialize parameter values using a pin file
    int doRetro    = 0;//flag to facilitate a retrospective model run
    int fitSimData = 0;//flag to fit model to simulated data calculated in the PRELIMINARY_CALCs section
    
    int yRetro = 0; //number of years to decrement for retrospective model run
    int iSeed = -1; //default random number generator seed
    random_number_generator rng(iSeed);//random number generator
    int iSimDataSeed = 0;
    random_number_generator rngSimData(-1);//random number generator for data simulation
    
    int doOFL = 0;///<flag (0/1) to do OFL calculations
    
    //debug flags
    int debugModelConfig     = 0;
    int debugModelDatasets   = 0;
    int debugModelParamsInfo = 0;
    int debugModelParams     = 0;
    int debugModelOptions    = 0;
    
    int debugDATA_SECTION    = 0;
    int debugPARAMS_SECTION  = 0;
    int debugPRELIM_CALCS    = 0;
    int debugPROC_SECTION    = 0;
    int debugREPORT_SECTION  = 0;
    
    int showActiveParams = 0;    
    int debugRunModel    = 0;    
    int debugObjFun      = 0;
    int debugOFL         = 0;
    
    int debugMCMC = 0;
    
    //note: consider using std::bitset to implement debug functionality
    int dbgCalcProcs = 10;
    int dbgObjFun = 20;
    int dbgNLLs   = 25;
    int dbgPriors = tcsam::dbgPriors;
    int dbgPopDy  = 70;
    int dbgApply  = 80;
    int dbgDevs   = 90;
    int dbgAll    = tcsam::dbgAll;
    
    int nSXs    = tcsam::nSXs;
    int MALE    = tcsam::MALE;
    int FEMALE  = tcsam::FEMALE;
    int ALL_SXs = tcsam::ALL_SXs;
    
    int nMSs     = tcsam::nMSs;
    int IMMATURE = tcsam::IMMATURE;
    int MATURE   = tcsam::MATURE;
    int ALL_MSs  = tcsam::ALL_MSs;
    
    int nSCs      = tcsam::nSCs;
    int NEW_SHELL = tcsam::NEW_SHELL;
    int OLD_SHELL = tcsam::OLD_SHELL;
    int ALL_SCs   = tcsam::ALL_SCs;
    
    double smlVal = 0.00001;//small value to keep things > 0
    
// =============================================================================
// =============================================================================
DATA_SECTION

 LOCAL_CALCS
    rpt::echo<<"#Starting "<<model<<" (ver "<<modVer<<") Code"<<endl;
    rpt::echo<<"#Starting DATA_SECTION"<<endl;
    cout<<"#Starting "<<model<<" (ver "<<modVer<<") Code"<<endl;
    cout<<"#Starting DATA_SECTION"<<endl;
 END_CALCS

 //Set commandline options
 LOCAL_CALCS
    int on = 0;
    int flg = 0;
    rpt::echo<<"#------Reading command line options---------"<<endl;
    //configFile
    fnConfigFile = "tcsam02.ModelConfig.dat";//default model config filename
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-configFile"))>-1) {
        fnConfigFile = ad_comm::argv[on+1];
        rpt::echo<<"#config file changed to '"<<fnConfigFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //resultsFile
    fnResultFile = "tcsam02.ResultFile";//default results file name
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resultsFile"))>-1) {
        fnResultFile = ad_comm::argv[on+1];
        rpt::echo<<"#results file changed to '"<<fnResultFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-binp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ainp"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //opModMode
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-opModMode"))>-1) {
        opModMode=1;
        rpt::echo<<"#operating model mode turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //calcOFL
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcOFL"))>-1) {
        doOFL=1;
        rpt::echo<<"#OFL calculations turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugOFL
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugOFL"))>-1) {
        debugOFL=1;
        rpt::echo<<"#debugOFL turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelDatasets
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelDatasets"))>-1) {
        debugModelDatasets=1;
        rpt::echo<<"#debugModelDatasets turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParamsInfo
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParamsInfo"))>-1) {
        debugModelParamsInfo=1;
        rpt::echo<<"#debugModelParamsInfo turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelOptions
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelOptions"))>-1) {
        debugModelOptions=1;
        rpt::echo<<"#debugModelOptions turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //doRetro
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-doRetro"))>-1) {
        doRetro=1;
        cout<<"#doRetro turned ON"<<endl;
        rpt::echo<<"#doRetro turned ON"<<endl;
        if (on+1<argc) {
            yRetro=atoi(ad_comm::argv[on+1]);
            cout<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
            rpt::echo<<"#Retrospective model run using yRetro = "<<yRetro<<endl;
        }
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //fitSimData
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-fitSimData"))>-1) {
        fitSimData=1;
        if (on+1<argc) {
            iSimDataSeed=atoi(ad_comm::argv[on+1]);
        } else {
            cout<<"-------------------------------------------"<<endl;
            cout<<"Enter random number seed (0 -> deterministic) for data simulation: ";
            cin>>iSimDataSeed;
        }
        if (iSimDataSeed) rng.reinitialize(iSimDataSeed);
        rpt::echo<<"#Simulating data to fit using "<<iSimDataSeed<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //jitter
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-jitter"))>-1) {
        jitter=1;
        iSeed=(long)start;
        if ((on=option_match(ad_comm::argc,ad_comm::argv,"-iSeed"))>-1) {
            if (on+1<argc) {
                iSeed=atoi(ad_comm::argv[on+1]);
            }
        } 
        rng.reinitialize(iSeed);
        rpt::echo<<"#Jittering for initial parameter values turned ON "<<endl;
        rpt::echo<<iSeed<<"  #iSeed"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //resample
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resample"))>-1) {
        resample=1;
        rpt::echo<<"#Resampling for initial parameter values turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelConfig
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelConfig"))>-1) {
        debugModelConfig=1;
        rpt::echo<<"#debugModelConfig turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugModelParams
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugModelParams"))>-1) {
        debugModelParams=1;
        cout<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#debugModelParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugDATA_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugDATA_SECTION"))>-1) {
        debugDATA_SECTION=1;
        rpt::echo<<"#debugDATA_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPARAMS_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPARAMS_SECTION"))>-1) {
        debugPARAMS_SECTION=1;
        rpt::echo<<"#debugPARAMS_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPRELIM_CALCS
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPRELIM_CALCS"))>-1) {
        debugPRELIM_CALCS=1;
        rpt::echo<<"debugPRELIM_CALCS turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugPROC_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugPROC_SECTION"))>-1) {
        debugPROC_SECTION=1;
        rpt::echo<<"#debugPROC_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugREPORT_SECTION
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugREPORT_SECTION"))>-1) {
        debugREPORT_SECTION=1;
        rpt::echo<<"#debugREPORT_SECTION turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugRunModel
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugRunModel")>-1) {
        debugRunModel=1;
        rpt::echo<<"#debugRunModel turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debugObjFun
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugObjFun")>-1) {
        debugObjFun=1;
        rpt::echo<<"#debugObjFun turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //showActiveParams
    if (option_match(ad_comm::argc,ad_comm::argv,"-showActiveParams")>-1) {
        showActiveParams=1;
        rpt::echo<<"#showActiveParams turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //debuMCMC
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugMCMC")>-1) {
        debugMCMC=1;
        rpt::echo<<"#debugMCMC turned ON"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
 END_CALCS
 
    int nZBs;  //number of model size bins
    int mnYr;  //min model year
    int mxYr;  //max model year
    int mxYrp1;//max model year + 1
    int nSel;  //number of selectivity functions
    int nFsh;  //number of fisheries
    int nSrv;  //number of surveys
 LOCAL_CALCS
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading configuration file '"<<fnConfigFile<<"'"<<endl;
    ad_comm::change_datafile_name(fnConfigFile);
    ptrMC = new ModelConfiguration();
    ptrMC->read(*(ad_comm::global_datafile));
    
    mnYr   = ptrMC->mnYr;
    mxYr   = ptrMC->mxYr;
    if (doRetro){mxYr = mxYr-yRetro; ptrMC->setMaxModelYear(mxYr);}
    if (jitter)   {ptrMC->jitter=1;}
    if (resample) {ptrMC->resample = 1;}
    
    rpt::echo<<"#------------------ModelConfiguration-----------------"<<endl;
    rpt::echo<<(*ptrMC);
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#----finished model configuration---"<<endl;
    rpt::echo<<"#-----------------------------------"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelConfiguration-----------------"<<endl;
        cout<<(*ptrMC);
        cout<<"#-----------------------------------"<<endl;
        cout<<"#----finished model configuration---"<<endl;
        cout<<"#-----------------------------------"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
    
    mxYrp1 = mxYr+1;
    nFsh   = ptrMC->nFsh;
    nSrv   = ptrMC->nSrv;
    nZBs   = ptrMC->nZBs;
 END_CALCS   
    vector zBs(1,nZBs)
    !!zBs  = ptrMC->zMidPts;
    
    //read model parameters info
 LOCAL_CALCS
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading parameters info file '"<<ptrMC->fnMPI<<"'"<<endl;
    if (debugModelParamsInfo) ModelParametersInfo::debug=1;
    ptrMPI = new ModelParametersInfo(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMPI);
    ptrMPI->read(*(ad_comm::global_datafile));
    if (debugModelParamsInfo) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelParamsInfo;
        if (debugModelParamsInfo<0) exit(1);
        ModelParametersInfo::debug=debugModelParamsInfo;
    }
    rpt::echo<<"#----finished reading model parameters info---"<<endl;
    rpt::echo<<"#-----------ModelParametersInfo---------------"<<endl;
    rpt::echo<<(*ptrMPI)<<endl;
    rpt::echo<<"#----finished ModelParametersInfo---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelParametersInfo-----------------"<<endl;
        cout<<(*ptrMPI);
        cout<<"#----finished model parameters info---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
 END_CALCS
        
    //read model data
 LOCAL_CALCS
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading datasets file '"<<ptrMC->fnMDS<<"'"<<endl;
    if (debugModelDatasets) {
        BioData::debug=1;
        FleetData::debug=1;
    }
    ptrMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrMDS->read(*(ad_comm::global_datafile));
    if (debugModelDatasets) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelDatasets;
        if (debugModelDatasets<0) exit(1);
        ModelDatasets::debug=debugModelDatasets;
        BioData::debug=debugModelDatasets;
        FleetData::debug=debugModelDatasets;
    }
    rpt::echo<<"#----finished model datasets---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelDatasets-----------------"<<endl;
        cout<<(*ptrMDS);
        cout<<"#----finished model datasets---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
 END_CALCS
    
    //read model data again to create SimMDS object
 LOCAL_CALCS
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading datasets file again to create SimMDS object '"<<ptrMC->fnMDS<<"'"<<endl;
    ptrSimMDS = new ModelDatasets(ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMDS);
    ptrSimMDS->read(*(ad_comm::global_datafile));
    rpt::echo<<"---SimMDS object after reading datasets---"<<endl;
    rpt::echo<<(*ptrSimMDS);
    rpt::echo<<"#finished SimMDS object"<<endl;
 END_CALCS   
    
    //read model options
 LOCAL_CALCS
    rpt::echo<<"#-----------------------------------"<<endl;
    rpt::echo<<"#Reading model options file '"<<ptrMC->fnMOs<<"'"<<endl;
    if (debugModelOptions) ModelOptions::debug=1;
    ptrMOs = new ModelOptions(*ptrMC);
    ad_comm::change_datafile_name(ptrMC->fnMOs);
    ptrMOs->read(*(ad_comm::global_datafile));
    if (debugModelOptions) {
        cout<<"enter 1 to continue : ";
        cin>>debugModelOptions;
        if (debugModelOptions<0) exit(1);
        ModelOptions::debug=debugModelOptions;
    }
    rpt::echo<<"#------------------ModelOptions-----------------"<<endl;
    rpt::echo<<(*ptrMOs);
    rpt::echo<<"#----finished model options---"<<endl;
    if (debugDATA_SECTION){
        cout<<"#------------------ModelOptions-----------------"<<endl;
        cout<<(*ptrMOs);
        cout<<"#----finished model options---"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debugDATA_SECTION;
        if (debugDATA_SECTION<0) exit(1);
    }
 END_CALCS
        
    //Match up model fisheries with fisheries data
    ivector mapD2MFsh(1,nFsh);
    ivector mapM2DFsh(1,nFsh);
 LOCAL_CALCS
    {
     int idx;
     for (int f=1;f<=nFsh;f++){
         idx = wts::which(ptrMDS->ppFsh[f-1]->name,ptrMC->lblsFsh);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying fishery names and labels in data file and config file."<<endl;
             cout<<"Incorrect fishery name in data file is '"<<ptrMDS->ppFsh[f-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MFsh(f)   = idx;//map from fishery data object f to model fishery idx
         mapM2DFsh(idx) = f;  //map from model fishery idx to fishery data object f
     }
     rpt::echo<<"model fisheries map to fishery data objects: "<<mapM2DFsh<<endl;
     cout<<"model fisheries map to fishery data objects: "<<mapM2DFsh<<endl;
    }
 END_CALCS
        
    //Match up model surveys with surveys data
    ivector mapD2MSrv(1,nSrv);
    ivector mapM2DSrv(1,nSrv);
 LOCAL_CALCS
    {
     int idx;
     for (int v=1;v<=nSrv;v++){
         idx = wts::which(ptrMDS->ppSrv[v-1]->name,ptrMC->lblsSrv);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying survey names and labels in data file and config file."<<endl;
             cout<<"Incorrect survey name in data file is '"<<ptrMDS->ppSrv[v-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MSrv(v)   = idx;//map from survey data object v to model survey idx
         mapM2DSrv(idx) = v;  //map from model survey idx to survey data object v
     }
     rpt::echo<<"model surveys map to survey data objects: "<<mapM2DSrv<<endl;
     cout<<"model surveys map to survey data objects: "<<mapM2DSrv<<endl;
    }
 END_CALCS
 
    //Extract parameter information
    //recruitment parameters
    int npLnR; ivector phsLnR; vector lbLnR; vector ubLnR;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pLnR,npLnR,lbLnR,ubLnR,phsLnR,rpt::echo);
    
    int npLnRCV; ivector phsLnRCV; vector lbLnRCV; vector ubLnRCV;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRCV,npLnRCV,lbLnRCV,ubLnRCV,phsLnRCV,rpt::echo);
    
    int npLgtRX; ivector phsLgtRX; vector lbLgtRX; vector ubLgtRX;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pLgtRX,npLgtRX,lbLgtRX,ubLgtRX,phsLgtRX,rpt::echo);
    
    int npLnRa; ivector phsLnRa; vector lbLnRa; vector ubLnRa;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRa,npLnRa,lbLnRa,ubLnRa,phsLnRa,rpt::echo);
    
    int npLnRb; ivector phsLnRb; vector lbLnRb; vector ubLnRb;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pLnRb,npLnRb,lbLnRb,ubLnRb,phsLnRb,rpt::echo);
    
    int npDevsLnR; ivector mniDevsLnR; ivector mxiDevsLnR; imatrix idxsDevsLnR;
    vector lbDevsLnR; vector ubDevsLnR; ivector phsDevsLnR;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pDevsLnR,npDevsLnR,mniDevsLnR,mxiDevsLnR,idxsDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,rpt::echo);
    
    //natural mortality parameters
    int npLnM; ivector phsLnM; vector lbLnM; vector ubLnM;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pLnM,npLnM,lbLnM,ubLnM,phsLnM,rpt::echo);
    
    int npLnDMT; ivector phsLnDMT; vector lbLnDMT; vector ubLnDMT;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMT,npLnDMT,lbLnDMT,ubLnDMT,phsLnDMT,rpt::echo);
    
    int npLnDMX; ivector phsLnDMX; vector lbLnDMX; vector ubLnDMX;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMX,npLnDMX,lbLnDMX,ubLnDMX,phsLnDMX,rpt::echo);
    
    int npLnDMM; ivector phsLnDMM; vector lbLnDMM; vector ubLnDMM;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMM,npLnDMM,lbLnDMM,ubLnDMM,phsLnDMM,rpt::echo);
    
    int npLnDMXM; ivector phsLnDMXM; vector lbLnDMXM; vector ubLnDMXM;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pLnDMXM,npLnDMXM,lbLnDMXM,ubLnDMXM,phsLnDMXM,rpt::echo);
    
    number zMref;
    !!zMref = ptrMPI->ptrNM->zRef;
    
    //maturity parameters
    int npLgtPrMat; ivector mniLgtPrMat; ivector mxiLgtPrMat; imatrix idxsLgtPrMat;
    vector lbLgtPrMat; vector ubLgtPrMat; ivector phsLgtPrMat;
    !!tcsam::setParameterInfo(ptrMPI->ptrM2M->pLgtPrM2M,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,idxsLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,rpt::echo);
 
    //growth parameters
    int npLnGrA; ivector phsLnGrA; vector lbLnGrA; vector ubLnGrA;
    !!tcsam::setParameterInfo(ptrMPI->ptrGrw->pLnGrA,npLnGrA,lbLnGrA,ubLnGrA,phsLnGrA,rpt::echo);
    
    int npLnGrB; ivector phsLnGrB; vector lbLnGrB; vector ubLnGrB;
    !!tcsam::setParameterInfo(ptrMPI->ptrGrw->pLnGrB,npLnGrB,lbLnGrB,ubLnGrB,phsLnGrB,rpt::echo);
    
    int npLnGrBeta; ivector phsLnGrBeta; vector lbLnGrBeta; vector ubLnGrBeta;
    !!tcsam::setParameterInfo(ptrMPI->ptrGrw->pLnGrBeta,npLnGrBeta,lbLnGrBeta,ubLnGrBeta,phsLnGrBeta,rpt::echo);
    
    //selectivity parameters
    !!nSel = ptrMPI->ptrSel->nPCs;//number of selectivity functions defined
    int npS1; ivector phsS1; vector lbS1; vector ubS1;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pS1,npS1,lbS1,ubS1,phsS1,rpt::echo);
    int npS2; ivector phsS2; vector lbS2; vector ubS2;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pS2,npS2,lbS2,ubS2,phsS2,rpt::echo);
    int npS3; ivector phsS3; vector lbS3; vector ubS3;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pS3,npS3,lbS3,ubS3,phsS3,rpt::echo);
    int npS4; ivector phsS4; vector lbS4; vector ubS4;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pS4,npS4,lbS4,ubS4,phsS4,rpt::echo);
    int npS5; ivector phsS5; vector lbS5; vector ubS5;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pS5,npS5,lbS5,ubS5,phsS5,rpt::echo);
    int npS6; ivector phsS6; vector lbS6; vector ubS6;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pS6,npS6,lbS6,ubS6,phsS6,rpt::echo);
    
    int npDevsS1; ivector mniDevsS1; ivector mxiDevsS1;  imatrix idxsDevsS1;
    vector lbDevsS1; vector ubDevsS1; ivector phsDevsS1;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS1,npDevsS1,mniDevsS1,mxiDevsS1,idxsDevsS1,lbDevsS1,ubDevsS1,phsDevsS1,rpt::echo);
    int npDevsS2; ivector mniDevsS2; ivector mxiDevsS2;  imatrix idxsDevsS2;
    vector lbDevsS2; vector ubDevsS2; ivector phsDevsS2;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS2,npDevsS2,mniDevsS2,mxiDevsS2,idxsDevsS2,lbDevsS2,ubDevsS2,phsDevsS2,rpt::echo);
    int npDevsS3; ivector mniDevsS3; ivector mxiDevsS3;  imatrix idxsDevsS3;
    vector lbDevsS3; vector ubDevsS3; ivector phsDevsS3;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS3,npDevsS3,mniDevsS3,mxiDevsS3,idxsDevsS3,lbDevsS3,ubDevsS3,phsDevsS3,rpt::echo);
    int npDevsS4; ivector mniDevsS4; ivector mxiDevsS4;  imatrix idxsDevsS4;
    vector lbDevsS4; vector ubDevsS4; ivector phsDevsS4;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS4,npDevsS4,mniDevsS4,mxiDevsS4,idxsDevsS4,lbDevsS4,ubDevsS4,phsDevsS4,rpt::echo);
    int npDevsS5; ivector mniDevsS5; ivector mxiDevsS5;  imatrix idxsDevsS5;
    vector lbDevsS5; vector ubDevsS5; ivector phsDevsS5;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS5,npDevsS5,mniDevsS5,mxiDevsS5,idxsDevsS5,lbDevsS5,ubDevsS5,phsDevsS5,rpt::echo);
    int npDevsS6; ivector mniDevsS6; ivector mxiDevsS6;  imatrix idxsDevsS6;
    vector lbDevsS6; vector ubDevsS6; ivector phsDevsS6;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pDevsS6,npDevsS6,mniDevsS6,mxiDevsS6,idxsDevsS6,lbDevsS6,ubDevsS6,phsDevsS6,rpt::echo);
    
        
    //fisheries parameters
    int npHM; ivector phsHM; vector lbHM; vector ubHM;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pHM,npHM,lbHM,ubHM,phsHM,rpt::echo);
    
    int npLnC; ivector phsLnC; vector lbLnC; vector ubLnC;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnC,npLnC,lbLnC,ubLnC,phsLnC,rpt::echo);
    
    int npLnDCT; ivector phsLnDCT; vector lbLnDCT; vector ubLnDCT;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCT,npLnDCT,lbLnDCT,ubLnDCT,phsLnDCT,rpt::echo);
    
    int npLnDCX; ivector phsLnDCX; vector lbLnDCX; vector ubLnDCX;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCX,npLnDCX,lbLnDCX,ubLnDCX,phsLnDCX,rpt::echo);
    
    int npLnDCM; ivector phsLnDCM; vector lbLnDCM; vector ubLnDCM;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCM,npLnDCM,lbLnDCM,ubLnDCM,phsLnDCM,rpt::echo);
    
    int npLnDCXM; ivector phsLnDCXM; vector lbLnDCXM; vector ubLnDCXM;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnDCXM,npLnDCXM,lbLnDCXM,ubLnDCXM,phsLnDCXM,rpt::echo);
    
    int npDevsLnC; ivector mniDevsLnC; ivector mxiDevsLnC; imatrix idxsDevsLnC;
    vector lbDevsLnC; vector ubDevsLnC; ivector phsDevsLnC;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pDevsLnC,npDevsLnC,mniDevsLnC,mxiDevsLnC,idxsDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,rpt::echo);
    
    //surveys parameters
    int npLnQ; ivector phsLnQ; vector lbLnQ; vector ubLnQ;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnQ,npLnQ,lbLnQ,ubLnQ,phsLnQ,rpt::echo);
    
    int npLnDQT; ivector phsLnDQT; vector lbLnDQT; vector ubLnDQT;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQT,npLnDQT,lbLnDQT,ubLnDQT,phsLnDQT,rpt::echo);
    
    int npLnDQX; ivector phsLnDQX; vector lbLnDQX; vector ubLnDQX;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQX,npLnDQX,lbLnDQX,ubLnDQX,phsLnDQX,rpt::echo);
    
    int npLnDQM; ivector phsLnDQM; vector lbLnDQM; vector ubLnDQM;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQM,npLnDQM,lbLnDQM,ubLnDQM,phsLnDQM,rpt::echo);
    
    int npLnDQXM; ivector phsLnDQXM; vector lbLnDQXM; vector ubLnDQXM;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pLnDQXM,npLnDQXM,lbLnDQXM,ubLnDQXM,phsLnDQXM,rpt::echo);

    //other data
    vector dtF_y(mnYr,mxYr);//timing of midpoint of fishing season (by year)
    !!dtF_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
    vector dtM_y(mnYr,mxYr);//timing of mating (by year))
    !!dtM_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);

    //model options
    ivector optsFcAvg(1,nFsh);//option flags for fishery capture rate averaging
    !!optsFcAvg = ptrMOs->optsFcAvg;
    imatrix yrsAvgEff_fy(1,nFsh,mnYr,mxYr);//years corresponding to average effort
    matrix eff_fy(1,nFsh,mnYr,mxYr);       //effort for averaging over fishery-specific time periods
    vector avgEff(1,nFsh);                 //average effort over fishery-specific time periods
    int optGrowth;    //growth function option
    !!optGrowth = ptrMOs->optGrowth;
    int optInitNatZ;  //initial n-at-z option
    !!optInitNatZ = ptrMOs->optInitNatZ;
    int optPenNonDecLgtPrMat;//option for penalty function on non-decreasing prMat parameters
    !!optPenNonDecLgtPrMat = ptrMOs->optPenNonDecLgtPrMat;

    //number of parameter combinations for various processes
    int npcRec;
    !!npcRec = ptrMPI->ptrRec->nPCs;
    int npcNM;
    !!npcNM = ptrMPI->ptrNM->nPCs;
    int npcM2M;
    !!npcM2M = ptrMPI->ptrM2M->nPCs;
    int npcGrw;
    !!npcGrw = ptrMPI->ptrGrw->nPCs;
    int npcSel;
    !!npcSel = ptrMPI->ptrSel->nPCs;
    int npcFsh;
    !!npcFsh = ptrMPI->ptrFsh->nPCs;
    int npcSrv;
    !!npcSrv = ptrMPI->ptrSrv->nPCs;
    
    imatrix idxDevsLnC_fy(1,nFsh,mnYr,mxYr); //matrix to check devs indexing for lnC
    
    //dimensions for output to R
    !!yDms = ptrMC->dimYrsToR;//years (mny:mxy)
    !!xDms = ptrMC->dimSXsToR;//sex
    !!mDms = ptrMC->dimMSsToR;//maturity
    !!sDms = ptrMC->dimSCsToR;//shell condition
    !!fDms = ptrMC->dimFshToR;//fisheries
    !!vDms = ptrMC->dimSrvToR;//surveys
    !!ypDms = ptrMC->dimYrsP1ToR;//years (mny:asy)
    !!zbDms = ptrMC->dimZBsToR;//size bin midpoints
    !!zpDms = ptrMC->dimZPsToR;//size bin midpoints (alternative)
    !!zcDms = ptrMC->dimZCsToR;//size bin cuptoints
    
    //fisheries info
    imatrix hasF_fy(1,nFsh,mnYr,mxYr);//flags indicating fishery activity
    
    //counters for PROCEDURE_SECTION calls
    int ctrProcCalls;       //all calls
    int ctrProcCallsInPhase;//within phase
    !!ctrProcCalls        = 0;
    !!ctrProcCallsInPhase = 0;
    
 LOCAL_CALCS
    rpt::echo<<"#finished DATA_SECTION"<<endl;
    cout<<"#finished DATA_SECTION"<<endl;
//    exit(1);
 END_CALCS
// =============================================================================
// =============================================================================
INITIALIZATION_SECTION

// =============================================================================
// =============================================================================
PARAMETER_SECTION
    !!rpt::echo<<"#Starting PARAMETER_SECTION"<<endl;
    !!cout<<"#Starting PARAMETER_SECTION"<<endl;
 
    //recruitment parameters TODO: implement devs
    init_bounded_number_vector pLnR(1,npLnR,lbLnR,ubLnR,phsLnR);              //mean ln-scale recruitment
    init_bounded_number_vector pLnRCV(1,npLnRCV,lbLnRCV,ubLnRCV,phsLnRCV);    //ln-scale recruitment cv
    init_bounded_number_vector pLgtRX(1,npLgtRX,lbLgtRX,ubLgtRX,phsLgtRX);    //logit-scale male sex ratio
    init_bounded_number_vector pLnRa(1,npLnRa,lbLnRa,ubLnRa,phsLnRa);         //size distribution parameter
    init_bounded_number_vector pLnRb(1,npLnRb,lbLnRb,ubLnRb,phsLnRb);         //size distribution parameter
    init_bounded_vector_vector pDevsLnR(1,npDevsLnR,mniDevsLnR,mxiDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR);//ln-scale rec devs
    matrix devsLnR(1,npDevsLnR,mniDevsLnR,mxiDevsLnR+1);
   
    //natural mortality parameters
    init_bounded_number_vector pLnM(1,npLnM,lbLnM,ubLnM,phsLnM);               //base (immature males)
    init_bounded_number_vector pLnDMT(1,npLnDMT,lbLnDMT,ubLnDMT,phsLnDMT);     //main temporal offsets
    init_bounded_number_vector pLnDMX(1,npLnDMX,lbLnDMX,ubLnDMX,phsLnDMX);     //female offsets
    init_bounded_number_vector pLnDMM(1,npLnDMM,lbLnDMM,ubLnDMM,phsLnDMM);     //mature offsets
    init_bounded_number_vector pLnDMXM(1,npLnDMXM,lbLnDMXM,ubLnDMXM,phsLnDMXM);//female-mature offsets
    
    //growth parameters
    init_bounded_number_vector pLnGrA(1,npLnGrA,lbLnGrA,ubLnGrA,phsLnGrA); //ln-scale mean growth coefficient "a"
    init_bounded_number_vector pLnGrB(1,npLnGrB,lbLnGrB,ubLnGrB,phsLnGrB); //ln-scale mean growth coefficient "b"
    init_bounded_number_vector pLnGrBeta(1,npLnGrBeta,lbLnGrBeta,ubLnGrBeta,phsLnGrBeta);//ln-scale growth scale parameter
    
    //maturity parameters
    init_bounded_vector_vector pLgtPrM2M(1,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat);//logit-scale maturity ogive parameters
    
    //selectivity parameters
    init_bounded_number_vector pS1(1,npS1,lbS1,ubS1,phsS1);
    init_bounded_number_vector pS2(1,npS2,lbS2,ubS2,phsS2);
    init_bounded_number_vector pS3(1,npS3,lbS3,ubS3,phsS3);
    init_bounded_number_vector pS4(1,npS4,lbS4,ubS4,phsS4);
    init_bounded_number_vector pS5(1,npS5,lbS5,ubS5,phsS5);
    init_bounded_number_vector pS6(1,npS6,lbS6,ubS6,phsS6);
    init_bounded_vector_vector pDevsS1(1,npDevsS1,mniDevsS1,mxiDevsS1,lbDevsS1,ubDevsS1,phsDevsS1);
    init_bounded_vector_vector pDevsS2(1,npDevsS2,mniDevsS2,mxiDevsS2,lbDevsS2,ubDevsS2,phsDevsS2);
    init_bounded_vector_vector pDevsS3(1,npDevsS3,mniDevsS3,mxiDevsS3,lbDevsS3,ubDevsS3,phsDevsS3);
    init_bounded_vector_vector pDevsS4(1,npDevsS4,mniDevsS4,mxiDevsS4,lbDevsS4,ubDevsS4,phsDevsS4);
    init_bounded_vector_vector pDevsS5(1,npDevsS5,mniDevsS5,mxiDevsS5,lbDevsS5,ubDevsS5,phsDevsS5);
    init_bounded_vector_vector pDevsS6(1,npDevsS6,mniDevsS6,mxiDevsS6,lbDevsS6,ubDevsS6,phsDevsS6);
    matrix devsS1(1,npDevsS1,mniDevsS1,mxiDevsS1+1);
    matrix devsS2(1,npDevsS2,mniDevsS2,mxiDevsS2+1);
    matrix devsS3(1,npDevsS3,mniDevsS3,mxiDevsS3+1);
    matrix devsS4(1,npDevsS4,mniDevsS4,mxiDevsS4+1);
    matrix devsS5(1,npDevsS5,mniDevsS5,mxiDevsS5+1);
    matrix devsS6(1,npDevsS6,mniDevsS6,mxiDevsS6+1);
    
    //fishing capture rate parameters
    init_bounded_number_vector pHM(1,npHM,lbHM,ubHM,phsHM);                    //handling mortality
    init_bounded_number_vector pLnC(1,npLnC,lbLnC,ubLnC,phsLnC);               //ln-scale base fishing mortality (mature males)
    init_bounded_number_vector pLnDCT(1,npLnDCT,lbLnDCT,ubLnDCT,phsLnDCT);     //ln-scale year-block offsets
    init_bounded_number_vector pLnDCX(1,npLnDCX,lbLnDCX,ubLnDCX,phsLnDCX);     //female offsets
    init_bounded_number_vector pLnDCM(1,npLnDCM,lbLnDCM,ubLnDCM,phsLnDCM);     //immature offsets
    init_bounded_number_vector pLnDCXM(1,npLnDCXM,lbLnDCXM,ubLnDCXM,phsLnDCXM);//female-immature offsets
    init_bounded_vector_vector pDevsLnC(1,npDevsLnC,mniDevsLnC,mxiDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC);//ln-scale deviations
    matrix devsLnC(1,npDevsLnC,mniDevsLnC,mxiDevsLnC+1);
    
    //survey catchbility parameters
    init_bounded_number_vector pLnQ(1,npLnQ,lbLnQ,ubLnQ,phsLnQ);               //base (mature male))
    init_bounded_number_vector pLnDQT(1,npLnDQT,lbLnDQT,ubLnDQT,phsLnDQT);     //main temporal offsets
    init_bounded_number_vector pLnDQX(1,npLnDQX,lbLnDQX,ubLnDQX,phsLnDQX);     //female offsets
    init_bounded_number_vector pLnDQM(1,npLnDQM,lbLnDQM,ubLnDQM,phsLnDQM);     //immature offsets
    init_bounded_number_vector pLnDQXM(1,npLnDQXM,lbLnDQXM,ubLnDQXM,phsLnDQXM);//female-immature offsets
    
    //objective function value
    objective_function_value objFun;
    
    //population-related quantities
    matrix  spB_yx(mnYr,mxYr,1,nSXs);                        //mature (spawning) biomass at mating time
    5darray n_yxmsz(mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//numbers at size, July 1 year y
    5darray nmN_yxmsz(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//natural mortality (numbers) during year)
    5darray tmN_yxmsz(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//total mortality (numbers) during year)
    
    //recruitment-related quantities
    number initMnR;                //mean recruitment for initial size comps
    vector R_y(mnYr,mxYr);         //total number of recruits entering population on July 1, by year
    vector Rx_c(1,npcRec);         //male fraction of recruits by parameter combination
    matrix R_yx(mnYr,mxYr,1,nSXs); //sex-specific fraction of recruits by year
    matrix R_cz(1,npcRec,1,nZBs);  //size distribution of recruits by parameter combination
    matrix R_yz(mnYr,mxYr,1,nZBs); //size distribution of recruits by year
    3darray R_yxz(mnYr,mxYr,1,nSXs,1,nZBs);    //size distribution of recruits by year, sex
    matrix stdvDevsLnR_cy(1,npcRec,mnYr,mxYr); //ln-scale recruitment std. devs by parameter combination and year
    matrix zscrDevsLnR_cy(1,npcRec,mnYr,mxYr); //standardized ln-scale recruitment residuals by parameter combination and year
    
    //natural mortality-related quantities
    3darray M_cxm(1,npcNM,1,nSXs,1,nMSs);                  //natural mortality rate by parameter combination
    5darray M_yxmsz(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//size-specific natural mortality rate
    
    //maturity-related quantities
    matrix  prM2M_cz(1,npcM2M,1,nZBs);         //prob. of immature crab molting to maturity by parameter combination
    3darray prM2M_yxz(mnYr,mxYr,1,nSXs,1,nZBs);//prob. of immature crab molting to maturity given sex x, pre-molt size z
    
    //growth related quantities
    matrix mnGrZ_cz(1,npcGrw,1,nZBs);           //mean post-molt size by parameter combination, pre-molt size
    3darray prGr_czz(1,npcGrw,1,nZBs,1,nZBs);   //prob of growth to z (row=lefthand z index) from zp (col=righthand z index) by parameter combination
    4darray mnGrZ_yxsz(mnYr,mxYr,1,nSXs,1,nSCs,1,nZBs); //mean post-molt size by by year, sex, shell condition, pre-molt size
    5darray prGr_yxszz(mnYr,mxYr,1,nSXs,1,nSCs,1,nZBs,1,nZBs); //prob of growth to z (row=lefthand z index) from zp (col=righthand z index) by year, sex, shell condition
    matrix grA_xy(1,nSXs,mnYr,mxYr+1);    //"a" parameters for growth, by sex and year
    matrix grB_xy(1,nSXs,mnYr,mxYr+1);    //"b" parameters for growth, by sex and year
    matrix grBeta_xy(1,nSXs,mnYr,mxYr+1); //beta parameters for growth, by sex and year
    
    
    //Selectivity (and retention) functions
    matrix sel_cz(1,npcSel,1,nZBs);            //all selectivity functions (fisheries and surveys) by parameter combination (no devs))
    3darray sel_cyz(1,nSel,mnYr,mxYr+1,1,nZBs);//all selectivity functions (fisheries and surveys) by year
    
    //fishery-related quantities
    matrix dvsLnC_fy(1,nFsh,mnYr,mxYr);                      //matrix to capture lnC-devs
    5darray cpF_fxmsy(1,nFsh,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr);//capture rates contributing to averages over a fishery-specific period
    4darray avgFc_fxms(1,nFsh,1,nSXs,1,nMSs,1,nSCs);         //avg capture rate over a fishery-specific period
    4darray avgRatioFc2Eff(1,nFsh,1,nSXs,1,nMSs,1,nSCs);     //ratio of avg capture rate to effort
    matrix  hmF_fy(1,nFsh,mnYr,mxYr);                        //handling mortality
    5darray cpF_fyxms(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs);        //fully-selected fishing capture rates NOT mortality)
    6darray sel_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//selectivity functions
    6darray ret_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//retention functions
    6darray cpF_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//fishing capture rate (NOT mortality) by fishery
    6darray rmF_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//size-specific retained fishing mortality by fishery
    6darray dmF_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//size-specific discard mortality rate  by fishery
    5darray tmF_yxmsz(mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);        //total size-specific fishing mortality rate (over all fisheries)
        
    6darray cpN_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//catch at size in fishery f
    6darray dsN_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discarded catch (numbers, NOT mortality) at size in fishery f
    6darray rmN_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//retained catch (numbers=mortality) at size in fishery f
    6darray dmN_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discard catch mortality at size in fishery f
    
    6darray cpB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//catch at size in fishery f
    6darray dsB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discarded catch (numbers, NOT mortality) at size in fishery f
    6darray rmB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//retained catch (numbers=mortality) at size in fishery f
    6darray dmB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discard catch mortality at size in fishery f
    
    //survey-related quantities
    3darray mb_vyx(1,nSrv,mnYr,mxYr+1,1,nSXs);//mature (spawning) biomass at time of survey
    5darray q_vyxms(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs);        //fully-selected catchability in survey v
    6darray s_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//selectivity functions
    6darray q_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//size-specific catchability in survey v
    6darray n_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//catch abundance at size in survey v
    6darray b_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//catch biomass at size in survey v
    
    //objective function penalties
    vector fPenRecDevs(1,npDevsLnR);//recruitment devs penalties
    
    vector fPenSmoothLgtPrMat(1,npLgtPrMat);//smoothness penalties on pr(mature|z)
    vector fPenNonDecLgtPrMat(1,npLgtPrMat);//non-decreasing penalties on pr(mature|z)
    
    vector fPenDevsS1(1,npDevsS1);//penalties on S1 devs (selectivity parameters)
    vector fPenDevsS2(1,npDevsS2);//penalties on S2 devs (selectivity parameters)
    vector fPenDevsS3(1,npDevsS3);//penalties on S3 devs (selectivity parameters)
    vector fPenDevsS4(1,npDevsS4);//penalties on S4 devs (selectivity parameters)
    vector fPenDevsS5(1,npDevsS5);//penalties on S5 devs (selectivity parameters)
    vector fPenDevsS6(1,npDevsS6);//penalties on S6 devs (selectivity parameters)
    
    vector fPenDevsLnC(1,npDevsLnC);//penalties on LnC devs (fishery capture parameters)
    
    //initial numbers
    vector  R_z(1,nZBs);                         //recruitment size distribution
    3darray S1_msz(1,nMSs,1,nSCs,1,nZBs);        //survival until molting/mating
    matrix  Th_sz(1,nSCs,1,nZBs);                //pr(molt to maturity|pre-molt size, molt)
    3darray T_szz(1,nSCs,1,nZBs,1,nZBs);         //growth matrices (indep. of molt to maturity)
    3darray S2_msz(1,nMSs,1,nSCs,1,nZBs);        //survival after molting/mating
    4darray n_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs); //equilibrium size distribution
    
    //likelihood components
    vector nllRecDevs(1,npcRec);//negative log-likelihoods associated with recruitment
    
    //sdreport variables
    sdreport_vector sdrLnR_y(mnYr,mxYr);
    sdreport_matrix sdrSpB_xy(1,nSXs,mnYr,mxYr);

    
    !!cout<<"#finished PARAMETER_SECTION"<<endl;
    !!rpt::echo<<"#finished PARAMETER_SECTION"<<endl;
    
// =============================================================================
// =============================================================================
PRELIMINARY_CALCS_SECTION
    rpt::echo<<"#Starting PRELIMINARY_CALCS_SECTION"<<endl;
    cout<<"#Starting PRELIMINARY_CALCS_SECTION"<<endl;
    int debug=1;
    
    //set initial values for all parameters
    if (usePin) {
        rpt::echo<<"NOTE: setting initial values for parameters using pin file"<<endl;
    } else {
        rpt::echo<<"NOTE: setting initial values for parameters using setInitVals(...)"<<endl;
        setInitVals();
    }

    cout<<"testing setAllDevs()"<<endl;
    setAllDevs(0,rpt::echo);
        
    {cout<<"writing data to R"<<endl;
     ofstream echo1; echo1.open("ModelData.R", ios::trunc);
     ReportToR_Data(echo1,0,cout);
    }
    
    {cout<<"writing parameters info to R"<<endl;
     ofstream echo1; echo1.open("ModelParametersInfo.R", ios::trunc);
     ptrMPI->writeToR(echo1);
     echo1.close();
     cout<<"finished writing parameters info to R"<<endl;
     
     //write initial parameter values to csv
     cout<<"writing parameters info to csv"<<endl;
     ofstream os1("tcsam02.params.all.init.csv", ios::trunc);
     writeParameters(os1,0,0);//all parameters
     os1.close();
     ofstream os2("tcsam02.params.active.init.csv", ios::trunc);
     writeParameters(os2,0,1);//only parameters that will be active (i.e., phase>0)
     os2.close();
     cout<<"finished writing parameters info to csv"<<endl;
    }
    
    //calculate average effort for fisheries over specified time periods
    cout<<"calculating average effort"<<endl;
    rpt::echo<<"calculating average effort"<<endl;
    rpt::echo<<"mapD2MFsh = "<<mapD2MFsh<<endl;
    avgEff.initialize();;
    for (int fd=1;fd<=nFsh;fd++){//fishery data object
        rpt::echo<<"fd = "<<fd<<endl;
        if (ptrMDS->ppFsh[fd-1]->ptrEff){
            rpt::echo<<"fishery has effort"<<endl;
            IndexBlock* pib = ptrMDS->ppFsh[fd-1]->ptrEff->ptrAvgIB;
            int fm = mapD2MFsh(fd);//index of corresponding model fishery
            int mny = max(mnYr,pib->getMin());
            int mxy = min(mxYr,pib->getMax());
            rpt::echo<<"fm, mny, mxy = "<<fm<<tb<<pib->getMin()<<tb<<pib->getMax()<<endl;
            rpt::echo<<"fm, mny, mxy = "<<fm<<tb<<mny<<tb<<mxy<<endl;
            eff_fy(fm).deallocate();
            eff_fy(fm).allocate(mny,mxy);
            yrsAvgEff_fy(fm).deallocate();
            yrsAvgEff_fy(fm) = pib->getFwdIndexVector();
            for (int rc=0;rc<pib->nRCs;rc++){
                IndexRange* pir = pib->ppIRs[rc];
                mny = max(mnYr,pir->getMin());
                mxy = min(mxYr,pir->getMax());
                eff_fy(fm)(mny,mxy) = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(mny,mxy);
            }
            rpt::echo<<"eff_fy(fm)       = "<<eff_fy(fm)<<endl;
            rpt::echo<<"yrsAvgEff_fy(fm) = "<<yrsAvgEff_fy(fm)<<endl;
            avgEff(fm) = sum(eff_fy(fm))/yrsAvgEff_fy(fm).size();
            rpt::echo<<"avgEff(fm) = "<<avgEff(fm)<<endl;
        }//has effort
    }//fm
    rpt::echo<<"eff_fy       = "<<endl<<eff_fy<<endl;
    rpt::echo<<"yrsAvgEff_fy = "<<endl<<yrsAvgEff_fy<<endl;
    rpt::echo<<"avgEff       = "<<avgEff<<endl;
    rpt::echo<<"finished calculating average effort"<<endl;
    cout<<"finished calculating average effort"<<endl;

    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")<0) {
        cout<<"testing calcRecruitment():"<<endl;
        rpt::echo<<"testing calcRecruitment():"<<endl;
        calcRecruitment(dbgCalcProcs+1,rpt::echo);
        
        cout<<"testing calcNatMort():"<<endl;
        rpt::echo<<"testing calcNatMort():"<<endl;
        calcNatMort(dbgCalcProcs+1,rpt::echo);
        
        cout<<"testing calcGrowth():"<<endl;
        rpt::echo<<"testing calcGrowth():"<<endl;
        calcGrowth(dbgCalcProcs+1,rpt::echo);
        
        cout<<"testing calcMolt2Maturity():"<<endl;
        rpt::echo<<"testing calcMolt2Maturity():"<<endl;
        calcMolt2Maturity(dbgCalcProcs+1,rpt::echo);

        cout<<"testing calcSelectivities():"<<endl;
        rpt::echo<<"testing calcSelectivities():"<<endl;
        calcSelectivities(dbgCalcProcs+1,rpt::echo);

        cout<<"testing calcFisheryFs():"<<endl;
        rpt::echo<<"testing calcFisheryFs():"<<endl;
        calcFisheryFs(dbgCalcProcs+1,rpt::echo);

        cout<<"testing calcSurveyQs():"<<endl;
        rpt::echo<<"testing calcSurveyQs():"<<endl;
        calcSurveyQs(dbgCalcProcs+1,cout);

        cout<<"testing runPopDyMod():"<<endl;
        rpt::echo<<"testing runPopDyMod():"<<endl;
        runPopDyMod(dbgCalcProcs+1,cout);
        rpt::echo<<"n_yxm:"<<endl;
        for (int y=mnYr;y<=(mxYr+1);y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                   rpt::echo<<y<<cc;
                   rpt::echo<<tcsam::getSexType(x)<<cc;
                   rpt::echo<<tcsam::getMaturityType(m)<<cc;
                   rpt::echo<<sum(n_yxmsz(y,x,m))<<endl;
                }
            }
        }
        rpt::echo<<"n_yxmsz:"<<endl;
        for (int y=mnYr;y<=(mxYr+1);y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                       rpt::echo<<y<<cc;
                       rpt::echo<<tcsam::getSexType(x)<<cc;
                       rpt::echo<<tcsam::getMaturityType(m)<<cc;
                       rpt::echo<<tcsam::getShellType(s)<<cc;
                       rpt::echo<<n_yxmsz(y,x,m,s)<<endl;
                    }
                }
            }
        }
        
        if (doOFL&&debugOFL){
            cout<<"Test OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.init.txt", ios::trunc);
            echoOFL<<"----Testing calcOFL()"<<endl;
            calcOFL(mxYr+1,debugOFL,echoOFL);//updates oflResults
            oflResults.writeCSVHeader(echoOFL); echoOFL<<endl;
            oflResults.writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL<<"----Finished testing calcOFL()!"<<endl;
            echoOFL.close();
            cout<<"Finished testing OFL calculations!"<<endl;
        }

        if (fitSimData){
            cout<<"creating sim data to fit in model"<<endl;
            rpt::echo<<"creating sim data to fit in model"<<endl;
            createSimData(1,rpt::echo,iSimDataSeed,ptrMDS);//stochastic if iSimDataSeed<>0
            {cout<<"re-writing data to R"<<endl;
             rpt::echo<<"re-writing data to R"<<endl;
             ofstream echo1; echo1.open("ModelData.R", ios::trunc);
             ReportToR_Data(echo1,0,cout);
            }
        }
        
        cout<<"--Testing calcObjFun()"<<endl;
        rpt::echo<<"Testing calcObjFun()"<<endl;
        calcObjFun(-1,rpt::echo);
        testNaNs(value(objFun),"testing calcObjFun() in PRELIMINARY_CALCS_SECTION");
        rpt::echo<<"--Testing calcObjFun() again"<<endl;
        calcObjFun(dbgAll,rpt::echo);
        
        {cout<<"writing model results to R"<<endl;
            rpt::echo<<"writing initial model report to R"<<endl;
            ofstream echo1; echo1.open("tcsam02.init.rep", ios::trunc);
            ReportToR(echo1,1,cout);
        }
        
        {cout<<"writing model sim data to file"<<endl;
            rpt::echo<<"writing model sim data to file"<<endl;
            createSimData(1,rpt::echo,0,ptrSimMDS);//deterministic
            ofstream echo1; echo1.open("tcsam02.SimData.init.dat", ios::trunc);
            writeSimData(echo1,0,rpt::echo,ptrSimMDS);
        }
        cout<<"#finished PRELIMINARY_CALCS_SECTION"<<endl;
        rpt::echo<<"#finished PRELIMINARY_CALCS_SECTION"<<endl;
        rpt::echo<<"#----------------------------------"<<endl<<endl;
//        int tmp = 1;
//        cout<<"Enter 1 to continue > ";
//        cin>>tmp;
//        if (tmp<0) exit(-1);
    } else {
        writeMCMCHeader();
        cout<<"MCEVAL is on"<<endl;
        rpt::echo<<"MCEVAL is on"<<endl;
    }
    
    
// =============================================================================
// =============================================================================
PROCEDURE_SECTION

    ctrProcCalls++;       //increment procedure section calls counter
    ctrProcCallsInPhase++;//increment in-phase procedure section calls counter
    
    int dbg = 0; //dbgAll;
    if (dbg){
        std::cout<<"---ctrProcCalls = "<<ctrProcCalls<<tb<<ctrProcCallsInPhase<<endl;
        writeParameters(std::cout,0,1);
    }
    
    runPopDyMod(dbg,rpt::echo);

    calcObjFun(dbg,rpt::echo);
    
    if (sd_phase()){
        sdrLnR_y = log(R_y);
        for (int x=1;x<=nSXs;x++){
            for (int y=mnYr; y<=mxYr; y++){
                sdrSpB_xy(x,y) = spB_yx(y,x);
            }
        }
    }
    
    if (mceval_phase()){
        updateMPI(0, cout);
        writeMCMCtoR(mcmc);
    }

//*****************************************
FUNCTION setInitVals
    //recruitment parameters
    setInitVals(ptrMPI->ptrRec->pLnR,    pLnR,    0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRCV,  pLnRCV,  0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLgtRX,  pLgtRX,  0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRa,   pLnRa,   0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pLnRb,   pLnRb,   0,rpt::echo);
    setInitVals(ptrMPI->ptrRec->pDevsLnR,pDevsLnR,0,rpt::echo);

    //natural mortality parameters
    setInitVals(ptrMPI->ptrNM->pLnM,   pLnM,   0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMT, pLnDMT, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMX, pLnDMX, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMM, pLnDMM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrNM->pLnDMXM,pLnDMXM,0,rpt::echo);

    //growth parameters
    setInitVals(ptrMPI->ptrGrw->pLnGrA,   pLnGrA,   0,rpt::echo);
    setInitVals(ptrMPI->ptrGrw->pLnGrB,   pLnGrB,   0,rpt::echo);
    setInitVals(ptrMPI->ptrGrw->pLnGrBeta,pLnGrBeta,0,rpt::echo);

    //maturity parameters
    setInitVals(ptrMPI->ptrM2M->pLgtPrM2M,pLgtPrM2M,0,rpt::echo);

    //selectivity parameters
    setInitVals(ptrMPI->ptrSel->pS1, pS1,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS2, pS2,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS3, pS3,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS4, pS4,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS5, pS5,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pS6, pS6,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS1, pDevsS1,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS2, pDevsS2,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS3, pDevsS3,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS4, pDevsS4,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS5, pDevsS5,0,rpt::echo);
    setInitVals(ptrMPI->ptrSel->pDevsS6, pDevsS6,0,rpt::echo);

    //fully-selected fishing capture rate parameters
    setInitVals(ptrMPI->ptrFsh->pHM,     pHM,     0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnC,    pLnC,    0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCT,  pLnDCT,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCX,  pLnDCX,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCM,  pLnDCM,  0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pLnDCXM, pLnDCXM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrFsh->pDevsLnC,pDevsLnC,0,rpt::echo);

    //survey catchability parameters
    setInitVals(ptrMPI->ptrSrv->pLnQ,   pLnQ,   0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQT, pLnDQT, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQX, pLnDQX, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQM, pLnDQM, 0,rpt::echo);
    setInitVals(ptrMPI->ptrSrv->pLnDQXM,pLnDQXM,0,rpt::echo);

//----------------------------------------------------------------------------------
//write header to MCMC eval file
FUNCTION writeMCMCHeader
    mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::trunc);
    mcmc<<"mcmc=list("<<endl;
    mcmc.close();
    
//******************************************************************************
FUNCTION void writeMCMCtoR(ostream& mcmc,NumberVectorInfo* ptr)
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
//******************************************************************************
FUNCTION void writeMCMCtoR(ostream& mcmc,BoundedNumberVectorInfo* ptr)
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
//******************************************************************************
FUNCTION void writeMCMCtoR(ostream& mcmc,BoundedVectorVectorInfo* ptr)
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
//******************************************************************************
FUNCTION void writeMCMCtoR(ostream& mcmc,DevsVectorVectorInfo* ptr)
    mcmc<<ptr->name<<"="; ptr->writeFinalValsToR(mcmc);
    
//******************************************************************************
FUNCTION void writeMCMCtoR(ofstream& mcmc)
    mcmc.open((char *) fnMCMC, ofstream::out|ofstream::app);
    mcmc<<"list(objFun="<<objFun<<cc<<endl;
    //write parameter values
        //recruitment values
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnR);   mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRCV); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLgtRX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRa);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnRb);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pDevsLnR);  mcmc<<cc<<endl;

        //natural mortality parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMT); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pLnDMXM); mcmc<<cc<<endl;

        //growth parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pLnGrA); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pLnGrB); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pLnGrBeta); mcmc<<cc<<endl;

        //maturity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrM2M->pLgtPrM2M); mcmc<<cc<<endl;

        //selectivity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pS6); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS5); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pDevsS6); mcmc<<cc<<endl;

        //fully-selected fishing capture rate parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pHM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnC); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCT); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCX); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnDCXM); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDevsLnC); mcmc<<cc<<endl;

        //survey catchability parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnQ);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQT);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQX);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQM);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pLnDQXM); mcmc<<cc<<endl;
    
        //write other quantities
        mcmc<<"R_y="; wts::writeToR(mcmc,value(R_y)); mcmc<<cc<<endl;
        ivector bnds = wts::getBounds(spB_yx);
        mcmc<<"MB_xy="; wts::writeToR(mcmc,trans(value(spB_yx)),xDms,yDms); 
        if (doOFL){
            mcmc<<cc<<endl;
            calcOFL(mxYr+1,0,cout);//updates oflresults
            oflResults.writeToR(mcmc,"oflResults",0);//mcm<<cc<<endl;
        }
        
    mcmc<<")"<<cc<<endl;
    mcmc.close();
    
//******************************************************************************
FUNCTION void createSimData(int debug, ostream& cout, int iSimDataSeed, ModelDatasets* ptrSim)
    if (debug)cout<<"simulating model results as data"<<endl;
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vcN_fyxmsz = wts::value(cpN_fyxmsz);
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    for (int f=1;f<=nFsh;f++) {
        if (debug) cout<<"fishery f: "<<f<<endl;
        (ptrSim->ppFsh[f-1])->replaceFisheryCatchData(iSimDataSeed,rngSimData,vcN_fyxmsz(f),vrmN_fyxmsz(f),ptrSim->ptrBio->wAtZ_xmz);
    }
    for (int v=1;v<=nSrv;v++) {
        if (debug) cout<<"survey "<<v<<endl;
        (ptrSim->ppSrv[v-1])->replaceIndexCatchData(iSimDataSeed,rngSimData,vn_vyxmsz(v),ptrSim->ptrBio->wAtZ_xmz);
    }
    if (debug) cout<<"finished simulating model results as data"<<endl;
     
//******************************************************************************
FUNCTION void writeSimData(ostream& os, int debug, ostream& cout, ModelDatasets* ptrSim)
    if (debug)cout<<"writing model results as data"<<endl;
    for (int v=1;v<=nSrv;v++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppSrv[v-1]))<<endl;
    }
    //     cout<<4<<endl;
    for (int f=1;f<=nFsh;f++) {
        os<<"#------------------------------------------------------------"<<endl;
        os<<(*(ptrSim->ppFsh[f-1]))<<endl;
    }
    if (debug) cout<<"finished writing model results as data"<<endl;
     
//******************************************************************************
//* Function: void setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, int debug, ostream& cout)
//* 
//* Description: Sets initial values for a parameter vector.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  pI (NumberVectorInfo*) 
//*     pointer to NumberVectorInfo object
//*  p (param_init_number_vector&)
//*     parameter vector
//* Returns:
//*  void
//* Alters:
//*  p - changes initial values
//******************************************************************************
FUNCTION void setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();//initial values from parameter info
        dvector def = pI->getInitVals();//defaults are initial values
        for (int i=1;i<=np;i++) {
            p(i) = vls(i);  //assign initial value from parameter info
            NumberInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                p(i) = ptrI->drawInitVal(rng,ptrMC->vif);//assign initial value based on resampling prior pdf
            }
        }
        //p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).label()<<": "<<endl;
        rpt::echo<<tb<<"inits  : "<<vls<<endl;
        rpt::echo<<tb<<"default: "<<def<<endl;
        rpt::echo<<tb<<"actual : "<<p<<endl;
        if (debug>=dbgAll) {
            std::cout<<"InitVals for "<<p(1).label()<<": "<<endl;
            std::cout<<tb<<p<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    }
     
//******************************************************************************
//* Function: void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
//* 
//* Description: Sets initial values for a parameter vector.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  pI (BoundedNumberVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_number_vector&)
//*     parameter vector
//* Returns:
//*  void
//* Alters:
//*  p - changes initial values
//******************************************************************************
FUNCTION void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();//initial values from parameter info
        dvector def = 0.5*(pI->getUpperBounds()+pI->getLowerBounds());//defaults are midpoints of ranges
        for (int i=1;i<=np;i++) {
            p(i) = vls(i);  //assign initial value from parameter info
            BoundedNumberInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<"jittering "<<p(i).label()<<endl;
                p(i) = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//only done if parameter phase > 0
            } else 
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                p(i) = ptrI->drawInitVal(rng,ptrMC->vif);
            }
        }
        //p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).label()<<": "<<endl;
        rpt::echo<<tb<<"inits  : "<<vls<<endl;
        rpt::echo<<tb<<"default: "<<def<<endl;
        rpt::echo<<tb<<"actual : "<<p<<endl;
        if (debug>=dbgAll) {
            std::cout<<"InitVals for "<<p(1).label()<<": "<<endl;
            std::cout<<tb<<p<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    }

//******************************************************************************
//* Function: void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//* 
//* Description: Sets initial values for a vector of parameter vectors.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  pI (BoundedVectorVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_vector_vector&)
//*     parameter vector
//* Returns:
//*  void
//* Alters:
//*  p - changes initial values
//******************************************************************************
FUNCTION void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            rpt::echo<<"InitVals "<<p(i).label()<<":"<<endl;
            dvector pns = value(p(i));
            dvector vls = (*pI)[i]->getInitVals();//initial values from parameter info
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=vls(j);
            BoundedVectorInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<tb<<"jittering "<<p(i).label()<<endl;
                dvector rvs = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//get jittered values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                rpt::echo<<tb<<"resampling "<<p(i).label()<<endl;
                dvector rvs = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else {
                rpt::echo<<tb<<"No jittering or resampling "<<p(i).label()<<endl;
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            }
            if (debug>=dbgAll){
                std::cout<<"pns(i) = "<<pns<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
                std::cout<<"p(i)   = "<<p(i)<<endl;
            }
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }

//******************************************************************************
//* Function: void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//* 
//* Description: Sets initial values for a vector of parameter vectors.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  pI (BoundedVectorVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_vector_vector&)
//*     parameter vector
//* Returns:
//*  void
//* Alters:
//*  p - changes initial values
//******************************************************************************
FUNCTION void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            rpt::echo<<"InitVals "<<p(i).label()<<":"<<endl;
            dvector pns = value(p(i));
            dvector vls = (*pI)[i]->getInitVals();//initial values from parameter info
            if (debug>=dbgAll) std::cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<endl;
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=vls(j);
            DevsVectorInfo* ptrI = (*pI)[i];
            if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                rpt::echo<<tb<<"jittering "<<p(i).label()<<endl;
                dvector rvs = wts::jitterParameter(p(i), ptrMC->jitFrac, rng);//get jittered values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else
            if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                rpt::echo<<tb<<"resampling "<<p(i).label()<<endl;
                dvector rvs = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values
                for (int j=p(i).indexmin();j<=p(i).indexmax();j++) p(i,j)=rvs(j);
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"resampled values = "<<rvs<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            } else {
                rpt::echo<<tb<<"No jittering or resampling "<<p(i).label()<<endl;
                rpt::echo<<tb<<"pin values       = "<<pns<<endl;
                rpt::echo<<tb<<"info values      = "<<vls<<endl;
                rpt::echo<<tb<<"final values     = "<<p(i)<<endl;
            }
            if (debug>=dbgAll){
                std::cout<<"pns(i) = "<<pns<<endl;
                std::cout<<"vls(i) = "<<vls<<endl;
                std::cout<<"p(i)   = "<<p(i)<<endl;
            }
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }

//-------------------------------------------------------------------------------------
FUNCTION void setAllDevs(int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"starting setAllDevs()"<<endl;

    if (debug>=dbgAll) cout<<"setDevs() for pLnR"<<endl;
    tcsam::setDevs(devsLnR, pDevsLnR,debug,cout);

    if (debug>=dbgAll) cout<<"setDevs() for pDevsS1"<<endl;
    tcsam::setDevs(devsS1, pDevsS1,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS2"<<endl;
    tcsam::setDevs(devsS2, pDevsS2,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS3"<<endl;
    tcsam::setDevs(devsS3, pDevsS3,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS4"<<endl;
    tcsam::setDevs(devsS4, pDevsS4,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS5"<<endl;
    tcsam::setDevs(devsS5, pDevsS5,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS6"<<endl;
    tcsam::setDevs(devsS6, pDevsS6,debug,cout);
    
    if (debug>=dbgAll) cout<<"setDevs() for pDevsLnC"<<endl;
    tcsam::setDevs(devsLnC, pDevsLnC,debug,cout);
    
    if (debug>=dbgAll) cout<<"finished setAllDevs()"<<endl;


//-------------------------------------------------------------------------------------
//calculate equilibrium size distribution for unmfexploited population
FUNCTION void calcEqNatZF100(dvariable& R, int yr, int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting dvar calcEqNatZF100()"<<endl;

    n_xmsz.initialize();//equilibrium n-at-z
    for (int x=1;x<=nSXs;x++){
        S1_msz.initialize(); //survival until molting/mating
        Th_sz.initialize();  //pr(molt to maturity|pre-molt size, molt)
        T_szz.initialize();  //growth matrices (indep. of molt to maturity)
        S2_msz.initialize(); //survival after molting/mating
        R_z.initialize();    //recruitment size distribution
        R_z = R*R_yx(yr,x)*R_yz(yr);//initial mean recruitment by size
        for (int s=1;s<=nSCs;s++){
            Th_sz(s) = prM2M_yxz(yr,x); //pr(molt to maturity|pre-molt size, molt)
            for (int z=1;z<=nZBs;z++) T_szz(s,z) = prGr_yxszz(yr,x,s,z);//growth matrices
            for (int m=1;m<=nMSs;m++){ 
                S1_msz(m,s) = mfexp(-M_yxmsz(yr,x,m,s)*dtM_y(yr));      //survival until molting/growth/mating
                S2_msz(m,s) = mfexp(-M_yxmsz(yr,x,m,s)*(1.0-dtM_y(yr)));//survival after molting/growth/mating
            }//m
        }//s
        n_xmsz(x) = calcEqNatZ(R_z, S1_msz, Th_sz, T_szz, S2_msz, debug, cout);
    }
    
    if (debug>=dbgPopDy) cout<<"finished dvar calcEqNatZF100()"<<endl;

//-------------------------------------------------------------------------------------
//calculate sex-specific equilibrium size distribution
FUNCTION dvar3_array calcEqNatZ(dvar_vector& R_z,dvar3_array& S1_msz, dvar_matrix& Th_sz, dvar3_array& T_szz, dvar3_array& S2_msz, int debug, ostream& cout)
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgPopDy) cout<<"starting calcEqNatZ()"<<endl;

    //the equilibrium solution
    dvar3_array n_msz(1,nMSs,1,nSCs,1,nZBs); n_msz.initialize();
    
    //create an identity matrix
    dmatrix I = identity_matrix(1,nZBs);
    
    //--calc the state transition matrices
    int i = IMMATURE; 
    int m =   MATURE;
    int n = NEW_SHELL;
    int o = OLD_SHELL;
    //immature new shell crab
    dvar_matrix S2_in = wts::diag(S2_msz(i,n)); //pr(survival|size) for immature new shell crab after molting occurs
    dvar_matrix Tr_in = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab pre-terminal molt
    dvar_matrix Th_in = wts::diag(Th_sz(n));    //pr(molt to maturity|pre-molt size,new shell, molting)
    dvar_matrix Ph_in = identity_matrix(1,nZBs);//pr(molt|pre-molt size, new shell) [assumed that all immature crab molt]
    dvar_matrix S1_in = wts::diag(S1_msz(i,n)); //pr(survival|size) for immature new shell crab before molting occurs
    //immature old shell crab [shouldn't be any of these]
    dvar_matrix S2_io = wts::diag(S2_msz(i,o)); //pr(survival|size) for immature old shell crab after molting occurs (but they didn't molt)
    dvar_matrix Tr_io = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab pre-terminal molt
    dvar_matrix Th_io = wts::diag(Th_sz(o));    //pr(molt to maturity|pre-molt size,old shell, molting)
    dvar_matrix Ph_io = identity_matrix(1,nZBs);//pr(molt|pre-molt size, old shell) [assumed all immature crab molt]
    dvar_matrix S1_io = wts::diag(S1_msz(i,o)); //pr(survival|size) for immature old shell crab before molting occurs
    //mature new shell crab
    dvar_matrix Tr_mn = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab undergoing terminal molt (same as non-terminal molt)
    dvar_matrix Tr_mo = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab undergoing terminal molt (same as non-terminal molt)
    dvar_matrix S2_mn = wts::diag(S2_msz(m,n)); //pr(survival|size) for mature new shell crab after molting occurs
    dvar_matrix S1_mn = wts::diag(S1_msz(m,n)); //pr(survival|size) for mature new shell crab before molting occurs (but they won't molt)
    //mature old shell crab
    dvar_matrix S2_mo = wts::diag(S2_msz(m,o)); //pr(survival|size) for mature old shell crab after molting occurs (but they didn't molt)
    dvar_matrix S1_mo = wts::diag(S1_msz(m,o)); //pr(survival|size) for mature old shell crab before molting occurs (but they won't molt)
    
    //full state transition matrices
    dvar_matrix lA = S2_in * (I-Th_in) * Tr_in * Ph_in * S1_in;//imm, new -> imm, new
    dvar_matrix lB = S2_in * (I-Th_io) * Tr_io * Ph_io * S1_io;//imm, old -> imm, new
    dvar_matrix lC = S2_io * (I-Ph_in) * S1_in;                //imm, new -> imm, old
    dvar_matrix lD = S2_io * (I-Ph_io) * S1_io;                //imm, old -> imm, old
    dvar_matrix lE = S2_mn * Th_in * Tr_mn * Ph_in * S1_in;    //imm, new -> mat, new (terminal molt)
    dvar_matrix lF = S2_mn * Th_io * Tr_mo * Ph_io * S1_io;    //imm, old -> mat, new (terminal molt)
    dvar_matrix lG = S2_mo * S1_mn;                            //mat, new -> mat, old
    dvar_matrix lH = S2_mo * S1_mo;                            //mat, old -> mat, old
    //--done calculating transition matrices
    
    //calculate inverses of matrix quantities
    dvar_matrix iM1 = inv(I - lD);
    dvar_matrix iM2 = inv(I - lA - lB * iM1 * lC);
    dvar_matrix iM3 = inv(I - lH);
    
    //the equilibrium solution is
    n_msz(i,n) = iM2 * R_z;                         //immature, new shell
    n_msz(i,o) = iM1 * lC * n_msz(i,n);             //immature, old shell
    n_msz(m,n) = lE * n_msz(i,n) + lF * n_msz(i,o); //  mature, new shell
    n_msz(m,o) = iM3 * lG * n_msz(m,n);             //  mature, old shell
        
    if (debug>=dbgPopDy) cout<<"finished calcEqNatZ()"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return(n_msz);

//-------------OFL Calculations--------------
FUNCTION void calcOFL(int yr, int debug, ostream& cout)
    if (debug) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcOFL(yr,debug,cout)"<<endl;
        cout<<"year for projection = "<<yr<<endl;
    }

    //1. get initial population for "upcoming" year, yr
    d4_array n_xmsz = wts::value(n_yxmsz(yr));
    if (debug) {cout<<"  males_msz:"<<endl; wts::print(n_xmsz(  MALE),cout,1);}
    if (debug) {cout<<"females_msz:"<<endl; wts::print(n_xmsz(FEMALE),cout,1);}
    
    //2. set yr back one year to get population rates, etc., 
    //   from year prior to projection year
    yr = yr-1;//don't have pop rates, etc. for projection year
    if (debug) cout<<"year for pop rates = "<<yr<<endl;
    
    //3. Determine mean recruitment 
    //   1981 here corresponds to 1982 in TCSAM2013, the year recruitment enters
    //   the model population.
    dvector avgRec_x(1,nSXs);
    if (debug) cout<<"R dims: "<<R_y.indexmin()<<cc<<R_y.indexmax()<<endl;
    for (int x=1;x<=nSXs;x++) 
        avgRec_x(x)= value(mean(elem_prod(R_y(1981,yr),column(R_yx,x)(1981,yr))));
    if (debug) {
        cout<<"R_y(  1981:"<<yr<<")      = "<<R_y(1981,yr)<<endl;
        cout<<"R_yx((1981:"<<yr<<",MALE) = "<<column(R_yx,MALE)(1981,yr)<<endl;
        cout<<"Average recruitment = "<<avgRec_x<<endl;
    }

    //4. Determine population rates for next year, using yr
    double dtF = dtF_y(yr);//time at which fisheries occur
    double dtM = dtM_y(yr);//time at which mating occurs
    
    PopDyInfo* pPIM = new PopDyInfo(nZBs,dtM);//  males info
    pPIM->R_z   = value(R_yz(yr));
    pPIM->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(MALE);
    pPIM->M_msz = value(M_yxmsz(yr,MALE));
    pPIM->T_szz = value(prGr_yxszz(yr,MALE));
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = value(prM2M_yxz(yr,MALE));
    
    PopDyInfo* pPIF = new PopDyInfo(nZBs,dtM);//females info
    pPIF->R_z   = value(R_yz(yr));
    pPIF->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(FEMALE);
    pPIF->M_msz = value(M_yxmsz(yr,FEMALE));
    pPIF->T_szz = value(prGr_yxszz(yr,FEMALE));
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = value(prM2M_yxz(yr,FEMALE));
    if (debug) cout<<"calculated pPIM, pPIF."<<endl;
    
    //5. Determine fishery conditions for next year based on averages for recent years
        int oflAvgPeriodYrs = 1;
        //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
        //a. Calculate average handling mortality, retention curves and capture rates
        int ny;   //number of years fishery is active
        dvector avgHM_f(1,nFsh);
        avgHM_f.initialize();
        for (int f=1;f<=nFsh;f++){
            ny = 0;
            for (int y=yr-oflAvgPeriodYrs+1;y<=yr;y++){
                ny         += hasF_fy(f,y);
                avgHM_f(f) += value(hmF_fy(f,y));
            }
            avgHM_f(f) /= 1.0*ny;
        }
        if (debug) cout<<"avgHm_f = "<<avgHM_f<<endl;

        d5_array avgCapF_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged capture mortality
        d5_array avgRFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged retention function
        avgCapF_xfmsz.initialize();
        avgRFcn_xfmsz.initialize();
        for (int x=1;x<=nSXs;x++){
            for (int f=1;f<=nFsh;f++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) {
                        for (int z=1;z<=nZBs;z++){
                            ny = 0;
                            for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                                ny += hasF_fy(f,y);
                                avgCapF_xfmsz(x,f,m,s,z) += value(cpF_fyxmsz(f,y,x,m,s,z));
                                avgRFcn_xfmsz(x,f,m,s,z) += value(ret_fyxmsz(f,y,x,m,s,z));
                            }
                            avgCapF_xfmsz(x,f,m,s,z) /= 1.0*ny;
                            avgRFcn_xfmsz(x,f,m,s,z) /= 1.0*ny;
                        }
                    }
                }
            }
        }
        if (debug){
            for (int f=1;f<=nFsh;f++){
                cout<<"avgCapF_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgCapF_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
            }
        }
        
        CatchInfo* pCIM = new CatchInfo(nZBs,nFsh,dtF);//male catch info
        pCIM->setCaptureRates(avgCapF_xfmsz(MALE));
        pCIM->setRetentionFcns(avgRFcn_xfmsz(MALE));
        pCIM->setHandlingMortality(avgHM_f);
        double maxCapF = pCIM->findMaxTargetCaptureRate(cout);
        if (debug) cout<<"maxCapF = "<<maxCapF<<endl;
        
        CatchInfo* pCIF = new CatchInfo(nZBs,nFsh,dtF);//female catch info
        pCIF->setCaptureRates(avgCapF_xfmsz(FEMALE));
        pCIF->setRetentionFcns(avgRFcn_xfmsz(FEMALE));
        pCIF->setHandlingMortality(avgHM_f);
        pCIF->maxF = maxCapF;//need to set this for females
        
    //6. Create PopProjectors
        PopProjector* pPPM = new PopProjector(pPIM,pCIM);
        PopProjector* pPPF = new PopProjector(pPIF,pCIF);
        if (debug) cout<<"created pPPs."<<endl;
        
    //7. Create Equilibrium_Calculators
        Equilibrium_Calculator* pECM = new Equilibrium_Calculator(pPPM);
        Equilibrium_Calculator* pECF = new Equilibrium_Calculator(pPPF);
        if (debug) cout<<"created pECs."<<endl;
        
    //8. Define OFL_Calculator   
        OFL_Calculator*  pOC;
        if (debug) cout<<"declared pOC."<<endl;
        
    //9. Determine TIER LEVEL, define Tier_Calculators, calculate OFL
        int tier = 3;
        if (tier==3){
            //5. Determine Fmsy and Bmsy
            Tier3_Calculator* pT3CM = new Tier3_Calculator(0.35,pECM);
            Tier3_Calculator* pT3CF = new Tier3_Calculator(0.35,pECF);
            if (debug) cout<<"created pT3Cs."<<endl;
            pOC = new OFL_Calculator(pT3CM,pT3CF);
            if (debug) {
                cout<<"created pOC."<<endl;
//                OFL_Calculator::debug=1;
//                Tier3_Calculator::debug=1;
//                Equilibrium_Calculator::debug=0;
            }
            oflResults = pOC->calcOFLResults(avgRec_x,n_xmsz,cout);
            if (debug) {
                cout<<"calculated oflResults."<<endl;
                oflResults.writeCSVHeader(cout); cout<<endl;
                oflResults.writeToCSV(cout); cout<<endl;
//                OFL_Calculator::debug=0;
//                Tier3_Calculator::debug=0;
//                Equilibrium_Calculator::debug=0;
            }
            
            if (debug){
                oflResults.writeCSVHeader(cout); cout<<endl;
                oflResults.writeToCSV(cout); cout<<endl;
            }
        }//Tier 3 calculation
    
    if (debug) {
        int n = 100;
        MultiYearPopProjector* pMYPPM = new MultiYearPopProjector(pPPM);
        MultiYearPopProjector* pMYPPF = new MultiYearPopProjector(pPPF);
        pMYPPM->projectUnFished(n,avgRec_x(  MALE),n_xmsz(  MALE),cout);
        double myPP_B0 = pMYPPM->matBio_y(n);
        pMYPPM->project(n,avgRec_x(  MALE),oflResults.Fmsy,n_xmsz(  MALE),cout);
        pMYPPF->project(n,avgRec_x(FEMALE),oflResults.Fmsy,n_xmsz(FEMALE),cout);
        double myPP_Bmsy = pMYPPM->matBio_y(n);
        double myPP_prjB = pMYPPM->matBio_y(1);
        double myPP_MSY  = pMYPPM->totCM_y(n)+pMYPPF->totCM_y(n);
        double myPP_OFL  = pMYPPM->totCM_y(1)+pMYPPF->totCM_y(1);
        cout<<"#------------------------"<<endl;
        cout<<"MYPP B0     = "<<myPP_B0  <<". B0     = "<<oflResults.B0  <<endl;
        cout<<"MYPP Bmsy   = "<<myPP_Bmsy<<". Bmsy   = "<<oflResults.Bmsy<<endl;
        cout<<"MYPP prjMMB = "<<myPP_prjB<<". prjMMB = "<<oflResults.prjB<<endl;
        cout<<"MYPP MSY    = "<<myPP_MSY <<". MSY    = "<<oflResults.MSY <<endl;
        cout<<"MYPP OFL    = "<<myPP_OFL <<". OFL    = "<<oflResults.OFL <<endl;
        cout<<"#------------------------"<<endl;
        cout<<"finished calcOFL(yr,debug,cout)"<<endl<<endl<<endl;
    }
        
//-------------END OFL Calculations----------   
//-------------------------------------------------------------------------------------
FUNCTION void initPopDyMod(int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting initPopDyMod()"<<endl;
    
    spB_yx.initialize();
    n_yxmsz.initialize();
    nmN_yxmsz.initialize();
    tmN_yxmsz.initialize();
       
    setAllDevs(debug,cout);//set devs vectors
    
    calcRecruitment(debug,cout);//calculate recruitment
    calcNatMort(debug,cout);    //calculate natural mortality rates
    calcGrowth(debug,cout);     //calculate growth transition matrices
    calcMolt2Maturity(debug,cout);   //calculate maturity ogives
    
    calcSelectivities(debug,cout); //calculate selectivity functions
    calcFisheryFs(debug,cout);     //calculate fishery F's
    calcSurveyQs(debug,cout);      //calculate survey Q's
    
    if (optInitNatZ==0){
        //will build up population from recruitment (like TCSAM2013)
        //do nothing, because n_yxmsz has already been initialized to 0
    } else if (optInitNatZ==1){
        //use equilibrium calculation to set initial n-at-z (like gmacs)
        //assumes no fishing occurs before model start
        calcEqNatZF100(initMnR,mnYr,debug,cout);//calculate n_xmsz
        n_yxmsz(mnYr) = n_xmsz;
    } else {
        cout<<"Unrecognized option for initial n-at-z: "<<optInitNatZ<<endl;
        cout<<"Terminating!"<<endl;
        exit(-1);
    }

    
    if (debug>=dbgPopDy) cout<<"finished initPopDyMod()"<<endl;

//-------------------------------------------------------------------------------------
FUNCTION void runPopDyMod(int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting runPopDyMod()"<<endl;
    //initialize population model
    initPopDyMod(debug, cout);
    
    //run population model
    for (int y=mnYr;y<=mxYr;y++){
        doSurveys(y,debug,cout);
        runPopDyModOneYear(y,debug,cout);        
    }
    doSurveys(mxYr+1,debug,cout);//do final surveys
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyMod()"<<endl;
    
//-------------------------------------------------------------------------------------
//calculate surveys.
FUNCTION void doSurveys(int y,int debug,ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting doSurveys("<<y<<")"<<endl;

    for (int v=1;v<=nSrv;v++){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    n_vyxmsz(v,y,x,m,s) = elem_prod(q_vyxmsz(v,y,x,m,s),n_yxmsz(y,x,m,s));
                }
            }
        }
    }
    for (int v=1;v<=nSrv;v++){
        mb_vyx(v,y) = calcSpB(n_vyxmsz(v,y),y,debug,cout);
    }
    if (debug>=dbgPopDy) cout<<"finished doSurveys("<<y<<")"<<endl;

//-------------------------------------------------------------------------------------
FUNCTION void runPopDyModOneYear(int yr, int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"Starting runPopDyModOneYear("<<yr<<")"<<endl;

    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n2_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n3_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n4_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n5_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    
    if (dtF_y(yr)<=dtM_y(yr)){//fishery occurs BEFORE molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs BEFORE molting/growth/maturity"<<endl;
        //apply natural mortality before fisheries
        n1_xmsz = applyNatMort(n_yxmsz(yr),yr,dtF_y(yr),debug,cout);
        //conduct fisheries
        n2_xmsz = applyFshMort(n1_xmsz,yr,debug,cout);
        //apply natural mortality from fisheries to molting/growth/maturity
        if (dtF_y(yr)==dtM_y(yr)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,yr,dtM_y(yr)-dtF_y(yr),debug,cout);
        }
        //calc mature (spawning) biomass at time of mating, but BEFORE growth/maturity (TODO: does this make sense??)
        spB_yx(yr) = calcSpB(n3_xmsz,yr,debug,cout);
        //apply molting, growth and maturation
        n4_xmsz = applyMGM(n3_xmsz,yr,debug,cout);
        //apply natural mortality to end of year
        if (dtM_y(yr)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,yr,1.0-dtM_y(yr),debug,cout);
        }
    } else {              //fishery occurs AFTER molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs AFTER molting/growth/maturity"<<endl;
        //apply natural mortality before molting/growth/maturity
        n1_xmsz = applyNatMort(n_yxmsz(yr),yr,dtM_y(yr),debug,cout);
        //calc mature (spawning) biomass at time of mating (TODO: does this make sense??)
        spB_yx(yr) = calcSpB(n1_xmsz,yr,debug,cout);
        //apply molting, growth and maturation
        n2_xmsz = applyMGM(n1_xmsz,yr,debug,cout);
        //apply natural mortality from molting/growth/maturity to fisheries
        if (dtM_y(yr)==dtF_y(yr)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,yr,dtF_y(yr)-dtM_y(yr),debug,cout);
        }
        //conduct fisheries
        n4_xmsz = applyFshMort(n3_xmsz,yr,debug,cout);
        //apply natural mortality to end of year
        if (dtF_y(yr)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,yr,1.0-dtF_y(yr),debug,cout);
        }
    }
    
    //advance surviving individuals to next year
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n_yxmsz(yr+1,x,m,s) = n5_xmsz(x,m,s);
            }
        }
    }
    //add in recruits (NOTE: R_y(y) here corresponds to R_y(y+1) in TCSAM2013)
    for (int x=1;x<=nSXs;x++) n_yxmsz(yr+1,x,IMMATURE,NEW_SHELL) += R_yxz(yr,x);
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyModOneYear("<<yr<<")"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION dvar_vector calcSpB(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
    if (debug>dbgApply) cout<<"starting calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector spb(1,nSXs); spb.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int s=1;s<=nSCs;s++) spb(x) += n0_xmsz(x,MATURE,s)*ptrMDS->ptrBio->wAtZ_xmz(x,MATURE);//dot product here
    }
    if (debug>dbgApply) cout<<"finished calcSpB("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return spb;
    
//-------------------------------------------------------------------------------------
FUNCTION dvar4_array applyNatMort(dvar4_array& n0_xmsz, int y, double dt, int debug, ostream& cout)
    if (debug>dbgApply) cout<<"starting applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n1_xmsz(x,m,s) = elem_prod(mfexp(-M_yxmsz(y,x,m,s)*dt),n0_xmsz(x,m,s));//survivors
                nmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
                tmN_yxmsz(y,x,m,s) += n0_xmsz(x,m,s)-n1_xmsz(x,m,s); //natural mortality
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyNatMort("<<y<<cc<<dt<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
//-------------------------------------------------------------------------------------
FUNCTION dvar4_array applyFshMort(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
    if (debug>dbgApply) cout<<"starting applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar_vector tm_z(1,nZBs);//total mortality (numbers) by size
    dvar_vector tvF_z(1,nZBs);//total fishing mortality rate by size, for use in calculating fishing rate components
    dvector     tdF_z(1,nZBs);//total fishing mortality rate by size, for use in calculating fishing rate components
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);//numbers surviving fisheries
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                tmF_yxmsz(y,x,m,s) = 0.0;//total fishing mortality rate
                for (int f=1;f<=nFsh;f++) tmF_yxmsz(y,x,m,s) += rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s);
                n1_xmsz(x,m,s) = elem_prod(mfexp(-tmF_yxmsz(y,x,m,s)),n0_xmsz(x,m,s));//numbers surviving all fisheries
                tm_z = n0_xmsz(x,m,s)-n1_xmsz(x,m,s);  //numbers killed by all fisheries
                tmN_yxmsz(y,x,m,s) += tm_z;            //add in numbers killed by all fisheries to total killed
                
                //calculate fishing rate components (need to ensure NOT dividing by 0)
                tdF_z = value(tmF_yxmsz(y,x,m,s));
                tvF_z = elem_prod(1-wts::isEQ(tdF_z,0.0),tmF_yxmsz(y,x,m,s)) + 
                                  wts::isEQ(tdF_z,0.0);
//                cout<<"y,x,m,s = "<<y<<tb<<x<<tb<<m<<tb<<s<<endl;
//                cout<<"tdF_z      = "<<tdF_z<<endl;
//                cout<<"tdF_z==0.0 = "<<wts::isEQ(tdF_z,0.0)<<endl;
//                cout<<"tvF_z      = "<<tvF_z<<endl;
//                int tmp; cout<<"Enter 1 to continue > "; cin>>tmp; if (!tmp) exit(-1);
                for (int f=1;f<=nFsh;f++){                   
                    cpN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div( cpF_fyxmsz(f,y,x,m,s),tvF_z),tm_z);//numbers captured in fishery f
                    rmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(rmF_fyxmsz(f,y,x,m,s),tvF_z),tm_z); //retained mortality in fishery f (numbers)
                    dmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_div(dmF_fyxmsz(f,y,x,m,s),tvF_z),tm_z); //discards mortality in fishery f (numbers)
                    dsN_fyxmsz(f,y,x,m,s) = cpN_fyxmsz(f,y,x,m,s)-rmN_fyxmsz(f,y,x,m,s);//discarded catch (NOT mortality) in fishery f (numbers)                    
                }
            }
        }
    }
    if (debug>dbgApply) cout<<"finished applyFshMort("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
//------------------------------------------------------------------------------
//Apply molting/growth/maturity to population numbers    
FUNCTION dvar4_array applyMGM(dvar4_array& n0_xmsz, int y, int debug, ostream& cout)
    if (debug>dbgApply) cout<<"starting applyMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        n1_xmsz(x,IMMATURE,NEW_SHELL) = elem_prod(1.0-prM2M_yxz(y,x),prGr_yxszz(y,x,NEW_SHELL)*n0_xmsz(x,IMMATURE,NEW_SHELL));
        n1_xmsz(x,IMMATURE,OLD_SHELL) = 0.0;
        n1_xmsz(x,MATURE,NEW_SHELL)   = elem_prod(    prM2M_yxz(y,x),prGr_yxszz(y,x,NEW_SHELL)*n0_xmsz(x,IMMATURE,NEW_SHELL));
        n1_xmsz(x,MATURE,OLD_SHELL)   = n0_xmsz(x,MATURE,NEW_SHELL)+n0_xmsz(x,MATURE,OLD_SHELL);
    }
    if (debug>dbgApply) cout<<"finished applyNatMGM("<<y<<")"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return n1_xmsz;
    
//-------------------------------------------------------------------------------------
//calculate recruitment.
FUNCTION void calcRecruitment(int debug, ostream& cout)
    if (debug>dbgCalcProcs) cout<<"starting calcRecruitment()"<<endl;

    RecruitmentInfo* ptrRI = ptrMPI->ptrRec;
    
    R_y.initialize();
    Rx_c.initialize();
    R_yx.initialize();
    R_cz.initialize();
    R_yz.initialize();
    R_yxz.initialize();
    stdvDevsLnR_cy.initialize();
    zscrDevsLnR_cy.initialize();
    
    int k; int y;
    dvector dzs = zBs+(zBs[2]-zBs[1])/2.0-zBs[1];
    for (int pc=1;pc<=ptrRI->nPCs;pc++){
        ivector pids = ptrRI->getPCIDs(pc);
        k=ptrRI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvariable mnLnR    = pLnR(pids[k++]);
        dvariable lnRCV    = pLnRCV(pids[k++]);
        dvariable lgtRX    = pLgtRX(pids[k++]);
        dvariable lnRa     = pLnRa(pids[k++]);
        dvariable lnRb     = pLnRb(pids[k++]);
        if (debug>dbgCalcProcs){
            cout<<"pids  = "<<pids<<endl;
            cout<<"mnLnR = "<<mnLnR<<endl;
            cout<<"lnRCV = "<<lnRCV<<endl;
            cout<<"lgtRX = "<<lgtRX<<endl;
            cout<<"lnRa  = "<<lnRa<<endl;
            cout<<"lnRb  = "<<lnRb<<endl;
        }

        int useDevs = pids[k]; k++;
        dvariable mnR;   //mean recruitment
        dvariable varLnR;//ln-scale variance in recruitment
        dvar_vector dvsLnR;
        ivector idxDevsLnR;
        varLnR = log(1.0+mfexp(2.0*lnRCV));//ln-scale variance
        mnR    = mfexp(mnLnR+varLnR/2.0);  //mean recruitment
        if (useDevs) {
            dvsLnR     = devsLnR(useDevs);
            idxDevsLnR = idxsDevsLnR(useDevs);
            if (debug>dbgCalcProcs) {
                cout<<"lims(dvsLnR) = "<<dvsLnR.indexmin()<<cc<<dvsLnR.indexmax()<<endl;
                cout<<"idx(dvsLnR) = "<<idxDevsLnR<<endl;
                cout<<"dvsLnR = "<<dvsLnR<<endl;
            }
        }
        
        Rx_c(pc) = 1.0/(1.0+mfexp(-lgtRX));
        R_cz(pc) = elem_prod(pow(dzs,mfexp(lnRa-lnRb)-1.0),mfexp(-dzs/mfexp(lnRb)));
        R_cz(pc) /= sum(R_cz(pc));//normalize to sum to 1

        imatrix idxs = ptrRI->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if (y==mnYr) initMnR = mnR;
            if ((mnYr<=y)&&(y<=mxYr)){
                if (debug>dbgCalcProcs+10) cout<<"y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
                if (useDevs){
                    R_y(y) = mfexp(mnLnR+dvsLnR[idxDevsLnR[y]]);
                } else {
                    R_y(y) = mnR;
                }
                if (debug>dbgCalcProcs+10) cout<<"R_y(y)="<<R_y(y)<<tb;
                if (MALE==nSXs){
                    R_yx(y,MALE) = 1.0;//only tracking males
                } else {
                    R_yx(y,MALE)   = Rx_c(pc);
                    R_yx(y,FEMALE) = 1.0-R_yx(y,MALE);
                    if (debug>dbgCalcProcs+10) cout<<R_yx(y,MALE)<<endl;
                }

                R_yz(y) = R_cz(pc);
                if (debug>dbgCalcProcs+10) cout<<"R_yz(y)="<<R_yz(y)<<endl;
                
                for (int x=1;x<=nSXs;x++) R_yxz(y,x) = R_y(y)*R_yx(y,x)*R_yz(y);

                stdvDevsLnR_cy(pc,y) = sqrt(varLnR); //ln-scale std dev
                zscrDevsLnR_cy(pc,y) = dvsLnR(idxDevsLnR(y))/stdvDevsLnR_cy(pc,y);//standardized ln-scale rec devs
            } else {
                if (debug>dbgCalcProcs) cout<<"skipping y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
            }
        }//idx
    }//pc
    
    if (debug>dbgCalcProcs) {
        cout<<"R_y = "<<R_y<<endl;
        cout<<"R_yx(MALE) = "<<column(R_yx,MALE)<<endl;
        cout<<"R_yz = "<<endl<<R_yz<<endl;
        cout<<"zscr = "<<zscrDevsLnR_cy<<endl;
        cout<<"finished calcRecruitment()"<<endl;
    }

//******************************************************************************
//* Function: void calcNatMort(void)
//* 
//* Description: Calculates natural mortality rates for all years.
//* 
//* Inputs:
//*  none
//* Returns:
//*  void
//* Alters:
//*  M_yxmsz - year/sex/maturity state/shell condition/size-specific natural mortality rate
//******************************************************************************
FUNCTION void calcNatMort(int debug, ostream& cout)  
    if(debug>dbgCalcProcs) cout<<"Starting calcNatMort()"<<endl;
    
    NaturalMortalityInfo* ptrNM = ptrMPI->ptrNM;
    
    dvar_matrix lnM(1,nSXs,1,nMSs);
    dvar3_array M_xmz(1,nSXs,1,nMSs,1,nZBs);
    
    M_cxm.initialize();
    M_yxmsz.initialize();

    int y; 
    for (int pc=1;pc<=ptrNM->nPCs;pc++){
        if (debug>dbgCalcProcs) cout<<"pc = "<<pc<<endl;
        lnM.initialize();
        ivector pids = ptrNM->getPCIDs(pc);
        int k=ptrNM->nIVs+1;//1st parameter variable column
        if (debug>dbgCalcProcs) cout<<"pids = "<<pids(k,pids.indexmax())<<endl;
        //add in base (ln-scale) natural mortality (immature males)
        if (pids[k]) {
            if (debug>dbgCalcProcs) cout<<"Adding pLnM["<<pids[k]<<"]: "<<pLnM(pids[k])<<endl;
            for (int x=1;x<=nSXs;x++) lnM(x) += pLnM(pids[k]);
        }   k++;
        //add in main offset from base
        if (pids[k]) {
            if (debug>dbgCalcProcs) cout<<"Adding pLnDMT["<<pids[k]<<"]: "<<pLnDMT(pids[k])<<endl;
            for (int x=1;x<=nSXs;x++) lnM(x) += pLnDMT(pids[k]);
        } k++;
        if (FEMALE<=nSXs){
            //add in offset for females
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding pLnDMX["<<pids[k]<<"]: "<<pLnDMX(pids[k])<<endl;
                lnM(FEMALE) += pLnDMX(pids[k]);
            }          k++;
        }
        //add in offset for mature crab
        if (pids[k]) {
            if (debug>dbgCalcProcs) cout<<"Adding pLnMM["<<pids[k]<<"]: "<<pLnDMM(pids[k])<<endl;
            for (int x=1;x<=nSXs;x++) lnM(x,MATURE) += pLnDMM(pids[k]);
        } k++;
        if (FEMALE<=nSXs){
            //add in offset for mature females
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding pLnDMXM["<<pids[k]<<"]: "<<pLnDMXM(pids[k])<<endl;
                lnM(FEMALE,MATURE) += pLnDMXM(pids[k]);
            }  k++; //advance k to zScaling in pids
        }
        
        //convert from ln-scale to arithmetic scale
        M_cxm(pc) = mfexp(lnM);
        if (debug>dbgCalcProcs){
            cout<<"lnM:"<<endl;
            for (int x=1;x<=nSXs;x++) cout<<"lnM("<<x<<")= "<<lnM(x)<<endl;
            cout<<"M_xm:"<<endl;
            for (int x=1;x<=nSXs;x++) cout<<"M_xm("<<x<<")= "<<M_cxm(pc,x)<<endl;
        }
        
        //add in size-scaling, if requested
        M_xmz.initialize();
        if (pids[k]&&(current_phase()>=pids[k])) {
            if (debug>dbgCalcProcs) cout<<"adding size scaling"<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++) M_xmz(x,m) = M_cxm(pc,x,m)*(zMref/zBs);//factor in size dependence
            }
        } else {
            if (debug>dbgCalcProcs) cout<<"not adding size scaling"<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++) M_xmz(x,m) = M_cxm(pc,x,m);//no size dependence
            }
        }
        if (debug>dbgCalcProcs) cout<<"finished scaling"<<endl;
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrNM->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//only model index for natural mortality is year
            if ((mnYr<=y)&&(y<=mxYr)){
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) M_yxmsz(y,x,m,s) = M_xmz(x,m);
                    }
                }
            }
        }
    }//loop over pcs
    if (debug>dbgCalcProcs) {
        for (int y=mnYr;y<=mxYr;y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) cout<<"M_yxmsz("<<y<<tb<<x<<tb<<m<<tb<<s<<")="<<M_yxmsz(y,x,m,s)<<endl;
                }
            }
        }
        cout<<"Finished calcNatMort()"<<endl;
    }
    
//-------------------------------------------------------------------------------------
//calculate Pr(maturity-at-size)
FUNCTION void calcMolt2Maturity(int debug, ostream& cout)
    if (debug>dbgCalcProcs) cout<<"starting calcMolt2Maturity()"<<endl;

    Molt2MaturityInfo* ptrM2MI = ptrMPI->ptrM2M;
    
    prM2M_cz.initialize();
    prM2M_yxz.initialize();
    
    int k; int y; int x;
    for (int pc=1;pc<=ptrM2MI->nPCs;pc++){
        ivector pids = ptrM2MI->getPCIDs(pc);
        k=ptrM2MI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvar_vector lgtPrM2M = pLgtPrM2M(pids[k++]);
        int vmn = lgtPrM2M.indexmin();
        int vmx = lgtPrM2M.indexmax();
        if (debug>dbgCalcProcs){
            cout<<"pc = "<<pc<<". mn = "<<vmn<<", mx = "<<vmx<<endl;
            cout<<"lgtPrM2M = "<<lgtPrM2M<<endl;
        }

        prM2M_cz(pc) = 1.0;//default is 1
        prM2M_cz(pc)(vmn,vmx) = 1.0/(1.0+mfexp(-lgtPrM2M));
        if (debug>dbgCalcProcs){
            cout<<"pc = "<<pc<<". mn = "<<vmn<<", mx = "<<vmx<<endl;
            cout<<"prM2M = "<<prM2M_cz(pc)<<endl;
        }
        
        imatrix idxs = ptrM2MI->getModelIndices(pc);
        if (debug>dbgCalcProcs) cout<<"molt-to-maturity indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if ((mnYr<=y)&&(y<=mxYr)){
                x = idxs(idx,2);
                if (debug>dbgCalcProcs) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(x)<<endl;
                prM2M_yxz(y,x) = prM2M_cz(pc);//note: this change made a difference, but not sure why!
            }
        }
    }
    
    if (debug>dbgCalcProcs) {
        for (int y=mnYr;y<=mxYr;y++){
            for (int x=1;x<=nSXs;x++){
                cout<<"prM2M_yxz("<<y<<tb<<x<<")="<<prM2M_yxz(y,x)<<endl;
            }
        }
        cout<<"finished calcMolt2Maturity()"<<endl;
    }

//******************************************************************************
//* Function: void calcGrowth(void)
//* 
//* Description: Calculates growth transition matrices for all years.
//* 
//* Inputs:
//*  none
//* Returns:
//*  void
//* Alters:
//*  prGr_yxszz - year/sex/maturity state/size-specific growth transition matrices
//******************************************************************************
FUNCTION void calcGrowth(int debug, ostream& cout)  
    if(debug>dbgCalcProcs) cout<<"Starting calcGrowth()"<<endl;
    
    GrowthInfo* ptrGrI = ptrMPI->ptrGrw;
    
    dvariable grA;
    dvariable grB;
    dvariable grBeta;
    
    mnGrZ_cz.initialize();
    prGr_czz.initialize();
    mnGrZ_yxsz.initialize();
    prGr_yxszz.initialize();
    grA_xy.initialize();
    grB_xy.initialize();
    grBeta_xy.initialize();
    
    dvar_matrix prGr_zz(1,nZBs,1,nZBs);

    int y; int x;
    for (int pc=1;pc<=ptrGrI->nPCs;pc++){
        ivector pids = ptrGrI->getPCIDs(pc);
        int k=ptrGrI->nIVs+1;//1st parameter column
//        grA = mfexp(pLnGrA(pids[k])); k++; //"a" coefficient for mean growth
//        grB = mfexp(pLnGrB(pids[k])); k++; //"b" coefficient for mean growth
//        grBeta = mfexp(pLnGrBeta(pids[k])); k++; //shape factor for gamma function growth transition
        grA = pLnGrA(pids[k]); k++; //"a" coefficient for mean growth
        grB = pLnGrB(pids[k]); k++; //"b" coefficient for mean growth
        grBeta = pLnGrBeta(pids[k]); k++; //scale factor for gamma function growth transition
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<". grB:"<<tb<<grB<<". grBeta:"<<grBeta<<endl;
        }
        
        //compute growth transition matrix for this pc
        prGr_zz.initialize();
        dvariable invBeta = 1.0/grBeta;           //inverse scale for gamma density function
        dvar_vector mnZs = mfexp(grA)*pow(zBs,grB);//mean post-molt sizes with zBs as pre-molt sizes
        dvar_vector mnIs = mnZs - zBs;              //mean molt increments
        dvar_vector alIs = mnIs/grBeta;             //gamma density alpha (location) parameters
        //check all mean molt increments are > 0
        if (isnan(value(sum(sqrt(mnIs))))){
            ofstream os("GrowthReport."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
            std::cout<<"##Found negative growth increments!"<<endl;
            std::cout<<"jitter seed = "<<iSeed<<endl;
            std::cout<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"##Found negative growth increments!"<<endl;
            os<<"jitter seed = "<<iSeed<<endl;
            os<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<" grB:"<<tb<<grB<<" grBeta:"<<grBeta<<endl;
            os<<"zBs   = "<<zBs<<endl;
            os<<"mnZs  = "<<mnZs<<endl;
            os<<"mnIs  = "<<mnIs<<endl;
            os<<"sdIs  = "<<sqrt(mnIs*grBeta)<<endl;
            os<<"alIs  = "<<alIs<<endl;
            os.close();
            exit(-1);
        }
        if (optGrowth==0) {
            //old style (TCSAM2013)
            for (int z=1;z<nZBs;z++){//pre-molt growth bin
                dvector dZs =  zBs(z,nZBs) - zBs(z);//realized growth increments (note non-neg. growth only)
                if (debug) cout<<"dZs: "<<dZs.indexmin()<<":"<<dZs.indexmax()<<endl;
                //dvar_vector prs = elem_prod(pow(dZs,alZ(z)-1.0),mfexp(-dZs/grBeta)); //pr(dZ|z)
                dvar_vector prs = wts::log_gamma_density(dZs,alIs(z),invBeta);
                prs = mfexp(prs);//gamma pdf
                if (debug) cout<<"prs: "<<prs.indexmin()<<":"<<prs.indexmax()<<endl;
                if (prs.size()>10) prs(z+10,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                if (debug) cout<<prs<<endl;
                prs = prs/sum(prs);//normalize to sum to 1
                if (debug) cout<<prs<<endl;
                prGr_zz(z)(z,nZBs) = prs;
            }//zs
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else if (optGrowth==1){
            //using cumd_gamma function like gmacs
            for (int z=1;z<nZBs;z++){
                dvar_vector sclIs = (ptrMC->zCutPts(z+1,nZBs+1)-zBs(z))/grBeta;//scaled increments at size bin cut points
                dvar_vector cprs(z,nZBs); cprs.initialize();
                dvar_vector prs(z,nZBs); prs.initialize();
                cprs(z) = cumd_gamma(sclIs(z+1),alIs(z));
                prs(z)  = cprs(z);
                for (int zp=z+1;zp<=nZBs;zp++){
                    cprs(zp) = cumd_gamma(sclIs(zp+1),alIs(z));
                    prs(zp)  = cprs(zp)-cprs(zp-1);//cumulative pr from zCs(zp) to zCs(zp+1)
                }
                prs(nZBs)   += 1.0 - cprs(nZBs);//treat final size bin as accumulator
                //cout<<"prs indices: "<<prs.indexmin()<<"  "<<prs.indexmax()<<endl;
                //check sum = 1
                if (sfabs(1.0-sum(prs))>1.0e-10){
                    std::cout<<"##Errors in calculating growth transition matrix: sum NE 1"<<endl;
                    std::cout<<"jitter seed = "<<iSeed<<endl;
                    std::cout<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
                    ofstream os("GrowthReportSum1."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
                    os<<"Errors in calculating growth transition matrix: sum NE 1"<<endl;
                    os<<"jitter seed = "<<iSeed<<endl;
                    os<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
                    os<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<" grB:"<<tb<<grB<<" grBeta:"<<grBeta<<endl;
                    os<<"zBs   = "<<zBs<<endl;
                    os<<"mnZs  = "<<mnZs<<endl;
                    os<<"mnIs  = "<<mnIs<<endl;
                    os<<"sdIs  = "<<sqrt(mnIs*grBeta)<<endl;
                    os<<"alIs  = "<<alIs<<endl;
                    os<<"z = "<<z<<tb<<"zB = "<<zBs(z)<<endl;
                    os<<"cutpts = "<<ptrMC->zCutPts(z+1,nZBs+1)<<endl;
                    os<<"sclIs = "<<sclIs<<endl;
                    os<<"sum(prs) = "<<sum(prs)<<endl;
                    os<<"prs  = "<<prs<<endl;
                    os<<"cprs = "<<cprs<<endl;
                    os.close();
                    exit(-1);
                }
                if (prs.size()>10) prs(z+10,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                prs = prs/sum(prs);//normalize to sum to 1
                if (debug) cout<<prs<<endl;
                prGr_zz(z)(z,nZBs) = prs;
            }//zs
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else {
            cout<<"Unrecognized growth option: "<<optGrowth<<endl;
            cout<<"Terminating!"<<endl;
            exit(-1);
        }
        
        mnGrZ_cz(pc) = mnZs;        
        prGr_czz(pc) = trans(prGr_zz);//transpose so rows are post-molt (i.e., lefthand z index is post-molt, or "to") z's so n+ = prGr_zz*n
        
        if (isnan(value(sum(prGr_zz)))){
            ofstream os("GrowthReportPrZZ."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
            std::cout<<"##Found NaN in prGz_zz in calcGrowth!"<<endl;
            std::cout<<"jitter seed = "<<iSeed<<endl;
            std::cout<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"##Found NaN in prGz_zz in calcGrowth!"<<endl;
            os<<"jitter seed = "<<iSeed<<endl;
            os<<"---Phase = "<<current_phase()<<". num PROC calls = "<<ctrProcCalls<<". num PROC calls in phase = "<<ctrProcCallsInPhase<<endl;
            os<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<" grB:"<<tb<<grB<<" grBeta:"<<grBeta<<endl;
            os<<"zBs   = "<<zBs<<endl;
            os<<"mnZs  = "<<mnZs<<endl;
            os<<"mnIs  = "<<mnIs<<endl;
            os<<"sdIs  = "<<sqrt(mnIs*grBeta)<<endl;
            os<<"alIs  = "<<alIs<<endl;
            if (optGrowth==0) {
                //old style (TCSAM2013)
                for (int z=1;z<nZBs;z++){//pre-molt growth bin
                    dvector dZs =  zBs(z,nZBs) - zBs(z);//realized growth increments (note non-neg. growth only)
                    os<<"Zpre = "<<zBs(z)<<endl;
                    os<<"dZs: "<<dZs<<endl;
                    os<<"Zbs: "<<zBs(z,nZBs)<<endl;
                    //dvar_vector prs = elem_prod(pow(dZs,alZ(z)-1.0),mfexp(-dZs/grBeta)); //pr(dZ|z)
                    dvar_vector prs = wts::log_gamma_density(dZs,alIs(z),invBeta);
                    os<<"log(prs): "<<prs<<endl;
                    prs = mfexp(prs);//gamma pdf
                    os<<"prs     : "<<prs<<endl;
                    if (prs.size()>10) prs(z+10,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                    os<<"normalization factor = "<<sum(prs)<<endl;
                    os<<"prs     : "<<prs<<endl;
                }
            } else if (optGrowth==1){
                //using cumd_gamma function like gmacs
                for (int z=1;z<nZBs;z++){
                    dvar_vector sclIs = (ptrMC->zCutPts(z+1,nZBs+1)-zBs(z))/grBeta;//scaled increments at size bin cut points
                    dvar_vector prs(z,nZBs);
                    prs(z) = cumd_gamma(sclIs(z+1),alIs(z));
                    for (int zp=z+1;zp<=nZBs;zp++){
                        prs(zp) = cumd_gamma(sclIs(zp+1),alIs(z))-prs(zp-1);//cumulative pr from zCs(zp) to zCs(zp+1)
                    }
                    prs(nZBs) += 1.0 - cumd_gamma(sclIs(nZBs+1),alIs(z));//treat final size bin as accumulator
                    cout<<"zB = "<<zBs(z)<<tb<<"sum(prs) = "<<sum(prs)<<endl;
                    cout<<"prs = "<<prs<<endl;
                    if (prs.size()>10) prs(z+10,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                    os<<"normalization factor = "<<sum(prs)<<endl;
                    os<<"prs     : "<<prs<<endl;
                }//zs
            }//optGrowth
            
            d3_array val = value(prGr_czz);
            os<<"prGr_czz = "<<endl; wts::print(val, os, 1); os<<endl;
            os.close();
            exit(-1);
            //testNaNs(value(sum(prGr_zz)),"Calculating growth");
        }
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrGrI->getModelIndices(pc);
        if (debug) cout<<"growth indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1); //year index
            if ((mnYr<=y)&&(y<=mxYr)){
                x = idxs(idx,2); //sex index
                grA_xy(x,y) = grA;
                grB_xy(x,y) = grB;
                grBeta_xy(x,y)  = grBeta;
                for (int s=1;s<=nSCs;s++){
                    mnGrZ_yxsz(y,x,s) = mnGrZ_cz(pc);
                    prGr_yxszz(y,x,s) = prGr_czz(pc);
                    //for (int z=1;z<=nZBs;z++) prGr_yxszz(y,x,s,z) = prGr_czz(pc,z);
                }//s
            }
        }//idx
    }//pc
    
    //set values for mxYr+1 to those for mxYr
    for (int x=1;x<=nSXs;x++){
        grA_xy(x,mxYr+1) = grA_xy(x,mxYr);
        grB_xy(x,mxYr+1) = grB_xy(x,mxYr);
        grBeta_xy(x,mxYr+1) = grBeta_xy(x,mxYr);
    }
    if (debug>dbgCalcProcs) cout<<"finished calcGrowth()"<<endl;

//******************************************************************************
//* Function: void calcSelectivities(int debug=0, ostream& cout=std::cout)
//* 
//* Description: Calculates all selectivity functions.
//* 
//* Required inputs:
//*  none
//* Returns:
//*  void
//* Alters:
//*  sel_cyz - selectivity array
//******************************************************************************
FUNCTION void calcSelectivities(int debug, ostream& cout)  
    if(debug>dbgCalcProcs) cout<<"Starting calcSelectivities()"<<endl;
    
    SelectivityInfo* ptrSel = ptrMPI->ptrSel;

    double fsZ;             //fully selected size
    int idSel;              //selectivity function id
//    int idxFSZ = ptrSel->nIVs+ptrSel->nPVs+1;//index for fsZ in pids vector below
    int idxFSZ = 1;//index for fsZ in pXDs vector below
    
    ivector mniSelDevs(1,6);//min indices of devs vectors
    ivector mxiSelDevs(1,6);//max indices of devs vectors
    dvar_vector params(1,6);//vector for number_vector params
    dvar_vector paramsp(1,6);//vector for number_vector params + dev offsets
        
    sel_cz.initialize();//selectivities w/out deviations
    sel_cyz.initialize();//selectivity array

    int y;
    for (int pc=1;pc<=ptrSel->nPCs;pc++){
        params.initialize();
        ivector pids = ptrSel->getPCIDs(pc);
        dvector pXDs = ptrSel->getPCXDs(pc);
        //extract the number parameters
        int k=ptrSel->nIVs+1;//1st parameter variable column
        if (pids[k]) {params[1] = pS1(pids[k]);}   k++;
        if (pids[k]) {params[2] = pS2(pids[k]);}   k++;
        if (pids[k]) {params[3] = pS3(pids[k]);}   k++;
        if (pids[k]) {params[4] = pS4(pids[k]);}   k++;
        if (pids[k]) {params[5] = pS5(pids[k]);}   k++;
        if (pids[k]) {params[6] = pS6(pids[k]);}   k++;
        if (debug>dbgCalcProcs) {
            cout<<"pc: "<<pc<<tb<<"pids = "<<pids<<endl;
            cout<<tb<<"params:"<<tb<<params<<endl;
        }
        
        int useDevsS1=pids[k++];
        dvar_vector dvsS1; ivector idxDevsS1;
        if (useDevsS1){
            dvsS1 = devsS1(useDevsS1);
            idxDevsS1 = idxsDevsS1(useDevsS1);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS1) = "<<idxDevsS1<<endl;
                cout<<"dvsS1      = "<<dvsS1<<endl;
            }
        }
        int useDevsS2=pids[k++];
        dvar_vector dvsS2; ivector idxDevsS2;
        if (useDevsS2){
            dvsS2 = devsS2(useDevsS2);
            idxDevsS2 = idxsDevsS2(useDevsS2);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS2) = "<<idxDevsS2<<endl;
                cout<<"dvsS2      = "<<dvsS2<<endl;
            }
        }
        int useDevsS3=pids[k++];
        dvar_vector dvsS3; ivector idxDevsS3;
        if (useDevsS3){
            dvsS3 = devsS3(useDevsS3);
            idxDevsS3 = idxsDevsS3(useDevsS3);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS3) = "<<idxDevsS3<<endl;
                cout<<"dvsS3      = "<<dvsS3<<endl;
            }
        }
        int useDevsS4=pids[k++];
        dvar_vector dvsS4; ivector idxDevsS4;
        if (useDevsS4){
            dvsS4 = devsS4(useDevsS4);
            idxDevsS4 = idxsDevsS4(useDevsS4);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS4) = "<<idxDevsS4<<endl;
                cout<<"dvsS4      = "<<dvsS4<<endl;
            }
        }
        int useDevsS5=pids[k++];
        dvar_vector dvsS5; ivector idxDevsS5;
        if (useDevsS5){
            dvsS5 = devsS5(useDevsS5);
            idxDevsS5 = idxsDevsS5(useDevsS5);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS5) = "<<idxDevsS5<<endl;
                cout<<"dvsS5      = "<<dvsS5<<endl;
            }
        }
        int useDevsS6=pids[k++];
        dvar_vector dvsS6; ivector idxDevsS6;
        if (useDevsS6){
            dvsS6 = devsS6(useDevsS6);
            idxDevsS6 = idxsDevsS6(useDevsS6);
            if (debug>dbgCalcProcs){
                cout<<"idx(dvsS6) = "<<idxDevsS6<<endl;
                cout<<"dvsS6      = "<<dvsS6<<endl;
            }
        }

//        fsZ   = pids[idxFSZ];
        fsZ   = pXDs[idxFSZ];
        idSel = pids[ptrSel->nIVs+ptrSel->nPVs+idxFSZ+1];
        if (debug>dbgCalcProcs) cout<<tb<<"fsZ: "<<fsZ<<tb<<"idSel"<<tb<<idSel<<tb<<SelFcns::getSelFcnID(idSel)<<tb<<params<<endl;

        sel_cz(pc) = SelFcns::calcSelFcn(idSel, zBs, params, fsZ);
        if (debug>dbgCalcProcs) cout<<tb<<"pc = "<<pc<<tb<<"sel_cz: "<<sel_cz(pc)<<endl;
            
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSel->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//year
            if ((mnYr<=y)&&(y<=mxYr+1)){
                paramsp = params;//set paramsp equal to base params
                k=ptrSel->nIVs+1+6;//1st devs vector variable column
                if (useDevsS1){
                    if (idxDevsS1[y]){
                        if (debug>dbgCalcProcs) cout<<tb<<idx<<tb<<y<<tb<<useDevsS1<<tb<<idxDevsS1[y]<<tb<<paramsp[1]<<tb<<devsS1(useDevsS1,idxDevsS1[y])<<endl;
                        paramsp[1] += devsS1(useDevsS1,idxDevsS1[y]);
                    }
                }
                if (useDevsS2) if (idxDevsS2[y]){paramsp[2] += devsS2(useDevsS2,idxDevsS2[y]);}
                if (useDevsS3) if (idxDevsS3[y]){paramsp[3] += devsS3(useDevsS3,idxDevsS3[y]);}
                if (useDevsS4) if (idxDevsS4[y]){paramsp[4] += devsS4(useDevsS4,idxDevsS4[y]);}
                if (useDevsS5) if (idxDevsS5[y]){paramsp[5] += devsS5(useDevsS5,idxDevsS5[y]);}
                if (useDevsS6) if (idxDevsS6[y]){paramsp[6] += devsS6(useDevsS6,idxDevsS6[y]);}
                sel_cyz(pc,y) = SelFcns::calcSelFcn(idSel, zBs, paramsp, fsZ);
                if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"paramsp: "<<paramsp<<tb<<"sel: "<<sel_cyz(pc,y)<<endl;
            } else {
                if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"y outside model range--skipping year!"<<endl;
            }
        }//idx
    }//pc
    if (debug>dbgCalcProcs) cout<<"finished calcSelectivities()"<<endl;

//******************************************************************************
//* Function: void calcFisheryFs(int debug, ostream& cout)
//* 
//* Description: Calculates fishery F's for all years.
//* 
//* Inputs:
//*  none
//* Returns:
//*  void
//* Alters:
//*  cpF_fyxms  - fully-selected fishery/year/sex/maturity state/shell condition-specific capture rate
//*  cpF_fyxmsz - fishery/year/sex/maturity state/shell condition/size-specific capture rate
//*  rmF_fyxmsz - fishery/year/sex/maturity state/shell condition/size-specific retained mortality rate
//*  dmF_fyxmsz - fishery/year/sex/maturity state/shell condition/size-specific discard mortality rate
//******************************************************************************
FUNCTION void calcFisheryFs(int debug, ostream& cout)  
    if(debug>dbgCalcProcs) cout<<"Starting calcFisheryFs()"<<endl;
    
    FisheriesInfo* ptrFsh = ptrMPI->ptrFsh;
    
    dvariable hm;           //handling mortality
    dvar_vector lnC_m(1,nMSs);//ln-scale capture rate
    dvar_vector arC_m(1,nMSs);//arithmetic-scale capture rate
    
    dvsLnC_fy.initialize();
    for (int f=1;f<=nFsh;f++) idxDevsLnC_fy(f) = -1;
    
    hasF_fy.initialize();   //flags indicating whether or not fishery occurs
    hmF_fy.initialize();    //handling mortality
    cpF_fyxms.initialize(); //fully-selected capture rate
    sel_fyxmsz.initialize();//selectivity functions
    ret_fyxmsz.initialize();//retention functions
    cpF_fyxmsz.initialize();//size-specific capture rate
    rmF_fyxmsz.initialize();//retention rate
    dmF_fyxmsz.initialize();//discard mortality rate
    tmF_yxmsz.initialize(); //total mortality rate
    
    /******************************************************\n
     * Fully-selected annual capture rates are calculated  \n
     * using 2 approaches:                                 \n
     *  1. from parameter values (if useER=0 below)        \n
     *  2. based on effort and the average ratio between   \n
     *     effort and capture rates over some time period  \n
     *     in the fishery        (if useER=1 below)        \n
     * Consequently, calculating all the above quantities  \n
     * requires 2 passes through parameter combinations.   \n
    ******************************************************/

    int idxER = ptrFsh->idxUseER;//index into pids below for flag to use effort ratio
    int y; int f; int x; int idSel; int idRet; int useER; int useDevs;
    //Pass 1: calculations based on parameter values
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids: "<<pids<<endl;
        useER = pids[idxER];//flag to use effort ratio
        if (!useER){//calculate capture rates from parameters
            lnC_m.initialize();
            arC_m.initialize();
            int k=ptrFsh->nIVs+1;//1st parameter variable column
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[k]) {hm = pHM(pids[k]);}        k++;
            //set base (ln-scale) capture rate (mature males)
            if (pids[k]) {lnC_m += pLnC(pids[k]);}   k++;
            //add in main offset 1 (temporal, perhaps)
            if (pids[k]) {lnC_m += pLnDCT(pids[k]);} k++;
            //add in main offset 2 (for females, perhaps)
            if (pids[k]) {lnC_m += pLnDCX(pids[k]);} k++;
            //add in immature offset 1
            if (pids[k]) {lnC_m(IMMATURE) += pLnDCM(pids[k]);} k++;
            //add in immature offset 2 (for immature females, perhaps)
            if (pids[k]) {lnC_m(IMMATURE) += pLnDCXM(pids[k]);} k++;

            //extract devs vector
            useDevs = pids[k]; k++;
            dvar_vector dvsLnC;             
            ivector idxDevsLnC;
            if (useDevs) {
                dvsLnC     = devsLnC(useDevs);
                idxDevsLnC = idxsDevsLnC(useDevs);
                if (debug>dbgCalcProcs){
                    cout<<"y   idx    devsLnC"<<endl;
                    for (int i=idxDevsLnC.indexmin();i<=idxDevsLnC.indexmax();i++) {
                        cout<<i<<tb<<idxDevsLnC(i)<<tb;
                        if (idxDevsLnC(i)) cout<<dvsLnC[idxDevsLnC(i)];
                        cout<<endl;
                    }//i
                }
            } else {
                arC_m = mfexp(lnC_m);
            }
            
            k = ptrFsh->nIVs+ptrFsh->nPVs+1;//1st extra variable column
            idSel = pids[k++];//selectivity function id
            idRet = pids[k++];//retention function id
        
            //convert from ln-scale to arithmetic scale
            if (debug>dbgCalcProcs){
                cout<<"pc: "<<pc<<". idSel = "<<idSel<<". idRet = "<<idRet<<tb<<"lnC:"<<lnC_m<<endl;
                if (useDevs) {
                    cout<<tb<<tb<<"dvsLnC["<<dvsLnC.indexmin()<<cc<<dvsLnC.indexmax()<<"] = "<<dvsLnC<<endl;
                } else {
                    cout<<tb<<tb<<"arC_m:"<<arC_m<<endl;
                }
            }

            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                x = idxs(idx,3);//sex
                if ((mnYr<=y)&&(y<=mxYr)){
                    hasF_fy(f,y) = 1;//flag indicating occurrence of fishery in year y
                    hmF_fy(f,y) = hm;//save discard mortality rate
                    if (debug>dbgCalcProcs) cout<<"f,y,x,useDevs = "<<f<<cc<<y<<cc<<x<<cc<<useDevs<<endl;
                    if (useDevs) {
                        idxDevsLnC_fy(f,y) = idxDevsLnC[y];
                        dvsLnC_fy(f,y)     = dvsLnC[idxDevsLnC[y]];
                        arC_m = mfexp(lnC_m+dvsLnC[idxDevsLnC[y]]);//recalculate C_xm w/ devs
                        if (debug>dbgCalcProcs) {
                            cout<<"idxDevsLnC[y],dvsLnC[idxDevsLnC[y]], arC_m: "<<idxDevsLnC[y]<<tb<<dvsLnC[idxDevsLnC[y]]<<tb<<arC_m<<endl;
                        }
                    }
                    for (int m=1;m<=nMSs;m++){
                        cpF_fyxms(f,y,x,m)  = arC_m(m); //fully-selected capture rate (independent of shell condition)
                        for (int s=1;s<=nSCs;s++){
                            sel_fyxmsz(f,y,x,m,s) = sel_cyz(idSel,y);          //selectivity
                            cpF_fyxmsz(f,y,x,m,s) = cpF_fyxms(f,y,x,m,s)*sel_cyz(idSel,y);//size-specific capture rate
                            if (idRet){//fishery has retention
                                ret_fyxmsz(f,y,x,m,s) = sel_cyz(idRet,y);      //retention curves
                                rmF_fyxmsz(f,y,x,m,s) = elem_prod(sel_cyz(idRet,y),         cpF_fyxmsz(f,y,x,m,s));//retention mortality
                                dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_cyz(idRet,y)),cpF_fyxmsz(f,y,x,m,s));//discard mortality
                            } else {//discard only
                                dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality
                            }
                        }//s
                    }//m
                }//(mnYr<=y)&&(y<=mxYr)
            }//idx
        }//useER=FALSE
    }//pc
    if (debug) cout<<"finished pass 1"<<endl;
    
    //calculate ratio of average capture rate to effort
    if (debug>dbgCalcProcs) cout<<"calculating avgRatioFc2Eff"<<endl;
    dvariable tot;
    avgRatioFc2Eff.initialize();
    for (int f=1;f<=nFsh;f++){//model fishery objects
        int fd = mapM2DFsh(f);//index of corresponding fishery data object
        if (ptrMDS->ppFsh[fd-1]->ptrEff){
            IndexBlock* pib = ptrMDS->ppFsh[fd-1]->ptrEff->ptrAvgIB;
            int mny = max(mnYr,pib->getMin());//adjust for model min
            int mxy = min(mxYr,pib->getMax());//adjust for model max
            if (debug>dbgCalcProcs) cout<<"f,mny,mxy = "<<f<<tb<<mny<<tb<<mxy<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++){
                        tot.initialize();
                        cpF_fxmsy(f,x,m,s).deallocate();
                        cpF_fxmsy(f,x,m,s).allocate(mny,mxy);//TODO: do this in preliminary calcs so it happens once!
                        cpF_fxmsy(f,x,m,s).initialize();
                        switch (optsFcAvg(f)){
                            case 0:
                                //won't use effort to fill in
                                avgFc_fxms(f,x,m,s) = -1;
                            case 1:
                                for (int y=yrsAvgEff_fy(f).indexmin();y<=yrsAvgEff_fy(f).indexmax();y++) 
                                    cpF_fxmsy(f,x,m,s,yrsAvgEff_fy(f,y)) = cpF_fyxms(f,yrsAvgEff_fy(f,y),x,m,s);
                                avgFc_fxms(f,x,m,s) = sum(cpF_fxmsy(f,x,m,s))/yrsAvgEff_fy(f).size(); break;
                            case 2:
                                for (int y=yrsAvgEff_fy(f).indexmin();y<=yrsAvgEff_fy(f).indexmax();y++) 
                                    cpF_fxmsy(f,x,m,s,yrsAvgEff_fy(f,y)) = 1.0-mfexp(-cpF_fyxms(f,yrsAvgEff_fy(f,y),x,m,s));
                                avgFc_fxms(f,x,m,s) = sum(cpF_fxmsy(f,x,m,s))/yrsAvgEff_fy(f).size(); break;
                            case 3:
                                for (int y=yrsAvgEff_fy(f).indexmin();y<=yrsAvgEff_fy(f).indexmax();y++)
                                    cpF_fxmsy(f,x,m,s,yrsAvgEff_fy(f,y)) = mean(cpF_fyxmsz(f,yrsAvgEff_fy(f,y),x,m,s));
                                avgFc_fxms(f,x,m,s) = sum(cpF_fxmsy(f,x,m,s))/yrsAvgEff_fy(f).size(); break;
                            default:
                                cout<<"optsFcAvg("<<f<<") = "<<optsFcAvg(f)<<" to calculate average Fc is invalid."<<endl;
                                cout<<"Aborting..."<<endl;
                                exit(-1);
                        }
                        if (optsFcAvg(f)) avgRatioFc2Eff(f,x,m,s) = avgFc_fxms(f,x,m,s)/avgEff(f);
                    }//s
                }//m
            }//x
        }//has ptrEff
    }//f
    if (debug>dbgCalcProcs) cout<<"calculated avgRatioFc2Eff"<<endl;
    
    //Pass 2: calculations based on effort and average effort:capture rate ratios
    int fd; double eff;
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        useER = pids[idxER];//flag to use effort ratio
        if (useER){//calculate capture rates from parameters
            int k=ptrFsh->nIVs+1;//1st parameter variable column
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[k]) {hm = pHM(pids[k]);} k++;
            
            k = ptrFsh->nIVs+ptrFsh->nPVs+1;//1st extra variable column
            idSel = pids[k++];   //selectivity function id
            idRet = pids[k++];   //retention function id
            
            if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<". hm = "<<hm<<". idSel = "<<idSel<<". idRet = "<<idRet<<". Using ER"<<endl;

            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                x = idxs(idx,3);//sex
                if ((mnYr<=y)&&(y<=mxYr)){
                    fd = mapM2DFsh(f);//index of corresponding fishery data object
                    eff = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(y);
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            //fully-selected capture rate
                            switch(optsFcAvg(f)) {
                                case 0:
                                    cpF_fyxms(f,y,x,m,s) = -1;
                                    break; //do nothing
                                case 1:
                                    cpF_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                                case 2:
                                    cpF_fyxms(f,y,x,m,s) = -log(1.0-avgRatioFc2Eff(f,x,m,s)*eff); break;
                                case 3:
                                    cpF_fyxms(f,y,x,m,s) = avgRatioFc2Eff(f,x,m,s)*eff; break;
                            }
                            if (debug>dbgCalcProcs) cout<<"f, y, x, m, s, eff, avgRatFcp2E, cpF = "<<f<<tb<<y<<tb<<x<<tb<<eff<<tb<<avgRatioFc2Eff(f,x,m,s)<<tb<<cpF_fyxms(f,y,x,m,s)<<endl;
                            if (!debug) testNaNs(value(cpF_fyxms(f,y,x,m,s)),"calcFisheryFs: 2nd pass");
                            sel_fyxmsz(f,y,x,m,s) = sel_cyz(idSel,y);
                            cpF_fyxmsz(f,y,x,m,s) = cpF_fyxms(f,y,x,m,s)*sel_cyz(idSel,y);//size-specific capture rate
                            if (idRet){//fishery has retention
                                ret_fyxmsz(f,y,x,m,s) = sel_cyz(idRet,y);
                                rmF_fyxmsz(f,y,x,m,s) = elem_prod(sel_cyz(idRet,y),         cpF_fyxmsz(f,y,x,m,s));//retention mortality rate
                                dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-sel_cyz(idRet,y)),cpF_fyxmsz(f,y,x,m,s));//discard mortality rate
                            } else {//discard only
                                dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality rate
                            }
                        }//s
                    }//m
                }
            }
        }//useER=TRUE
    }
    if (debug>dbgCalcProcs) {
        for (int f=1;f<=nFsh;f++){
            for (int y=mnYr;y<=mxYr;y++){
                cout<<"cpF_fyxmsz("<<f<<tb<<y<<tb<<"  MALE,m,s) = "<<cpF_fyxmsz(f,y,  MALE,MATURE,NEW_SHELL)<<endl;
                cout<<"cpF_fyxmsz("<<f<<tb<<y<<tb<<"FEMALE,m,s) = "<<cpF_fyxmsz(f,y,FEMALE,MATURE,NEW_SHELL)<<endl;
            }
        }
        cout<<"finished calcFisheryFs()"<<endl;
    }

//******************************************************************************
//* Function: void calcSurveyQs(int debug, ostream& cout)
//* 
//* Description: Calculates survey catchabilities for all years.
//* 
//* Inputs:
//*  none
//* Returns:
//*  void
//* Alters:
//*  q_vyxms  - fully-selected survey/year/sex/maturity state/shell condition-specific catchability
//*  q_vyxmsz - survey/year/sex/maturity state/shell condition/size-specific catchability
//******************************************************************************
FUNCTION void calcSurveyQs(int debug, ostream& cout)  
    if(debug>dbgCalcProcs) cout<<"Starting calcSurveyQs()"<<endl;
    
    SurveysInfo* ptrSrv = ptrMPI->ptrSrv;
    
    dvar_vector lnQ_m(1,nMSs);//ln-scale Q's
    dvar_vector arQ_m(1,nMSs);//arithmetic-scale Q's
    
    q_vyxms.initialize();
    q_vyxmsz.initialize();

    int y; int v; int x; int idSel;
    for (int pc=1;pc<=ptrSrv->nPCs;pc++){
        lnQ_m.initialize();
        arQ_m.initialize();
        ivector pids = ptrSrv->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids = "<<pids<<endl;
        
        int k=ptrSrv->nIVs+1;//1st parameter variable column
        //add in base (ln-scale) catchability (mature males)
        if (pids[k]) {lnQ_m += pLnQ(pids[k]);}   k++;
        //add in main offset 1 (for time period, for example)
        if (pids[k]) {lnQ_m += pLnDQT(pids[k]);} k++;
        //add in main offset 2 (for females, for example)
        if (pids[k]) {lnQ_m += pLnDQX(pids[k]);} k++;
        //add in immature offset 1
        if (pids[k]) {lnQ_m(IMMATURE) += pLnDQM(pids[k]);}  k++;
        //add in immature offset 2 (for females, for example)
        if (pids[k]) {lnQ_m(IMMATURE) += pLnDQXM(pids[k]);} k++; 
        
        idSel = pids[k];//selectivity function id
        
        //convert from ln-scale to arithmetic scale
        arQ_m = mfexp(lnQ_m);
        if (debug>dbgCalcProcs){
            cout<<"lnQ_m:"<<endl<<lnQ_m<<endl;
            cout<<"arQ_m:"<<endl<<arQ_m<<endl;
        }
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSrv->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            v = idxs(idx,1);//survey
            y = idxs(idx,2);//year
            x = idxs(idx,3);//sex
            if ((mnYr<=y)&&(y<=mxYrp1)){
                for (int m=1;m<=nMSs;m++){
                    q_vyxms(v,y,x,m) = arQ_m(m);//indep of s currently
                    for (int s=1;s<=nSCs;s++){
                        s_vyxmsz(v,y,x,m,s) = sel_cyz(idSel,y);
                        q_vyxmsz(v,y,x,m,s) = arQ_m(m)*sel_cyz(idSel,y);
                    }//s
                }//m
            }//(mnYr<=y)&&(y<=mxYrp1)
        }//idx
    }//pc
    if (debug>dbgCalcProcs) cout<<"finished calcSurveyQs()"<<endl;
    
//-------------------------------------------------------------------------------------
//Calculate penalties for objective function. TODO: finish
FUNCTION void calcPenalties(int debug, ostream& cout)
    if (debug>=dbgObjFun) cout<<"Started calcPenalties()"<<endl;
    if (debug<0) cout<<"list("<<endl;//start list of penalties by category

    if (debug<0) cout<<tb<<"maturity=list("<<endl;//start of maturity penalties list
    //smoothness penalties on maturity parameters (NOT maturity ogives)
    double penWgtSmthLgtPrMat = ptrMOs->wgtSmthLgtPrMat;
    fPenSmoothLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
    for (int i=1;i<npLgtPrMat;i++){
        dvar_vector v; v = 1.0*pLgtPrM2M(i);
        fPenSmoothLgtPrMat(i) = norm2(calc2ndDiffs(v));
        objFun += penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i)<<"),"<<endl;
    }
    {
        int i = npLgtPrMat;
        dvar_vector v; v = 1.0*pLgtPrM2M(i);
        fPenSmoothLgtPrMat(i) = norm2(calc2ndDiffs(v));
        objFun += penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat*fPenSmoothLgtPrMat(i)<<")"<<endl;
    }
    if (debug<0) cout<<tb<<tb<<")"<<cc<<endl;//end of smoothness penalties list

    //non-decreasing penalties on maturity parameters (NOT maturity ogives)
    double penWgtNonDecLgtPrMat = ptrMOs->wgtNonDecLgtPrMat;
    fPenNonDecLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"nondecreasing=list(";//start of non-decreasing penalties list
    for (int i=1;i<npLgtPrMat;i++){
        dvar_vector v; v = calc1stDiffs(pLgtPrM2M(i));
        if (optPenNonDecLgtPrMat==0){
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            } 
        } else if (optPenNonDecLgtPrMat==1){
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        }
        objFun += penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i)<<"),";
    }
    {
        int i = npLgtPrMat;
        dvar_vector v; v = calc1stDiffs(pLgtPrM2M(i));
        if (optPenNonDecLgtPrMat==0){
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            }
        } else if (optPenNonDecLgtPrMat==1){
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        }
        objFun += penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat*fPenNonDecLgtPrMat(i)<<")";
    }
    if (debug<0) cout<<tb<<tb<<")"<<endl;//end of non-decreasing penalties list    
    if (debug<0) cout<<tb<<"),";//end of maturity penalties list
    
    //penalties on final value of dev vectors
    double penWgt = 0.0;
    if (current_phase()>=ptrMOs->phsLastDevsPen) penWgt = ptrMOs->wgtLastDevsPen;
    if (debug<0) cout<<tb<<"final.devs=list("<<endl;//start of devs penalties list
    //recruitment devs
    if (ptrMPI->ptrRec->pDevsLnR->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsLnR=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnR,devsLnR);
        if (debug<0) cout<<cc<<endl;
    }
    //S1 devs
    if (ptrMPI->ptrSel->pDevsS1->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS1=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS1,devsS1);
        if (debug<0) cout<<cc<<endl;
    }
    //S2 devs
    if (ptrMPI->ptrSel->pDevsS2->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS2=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS2,devsS2);
        if (debug<0) cout<<cc<<endl;
    }
    //S3 devs
    if (ptrMPI->ptrSel->pDevsS3->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS3=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS3,devsS3);
        if (debug<0) cout<<cc<<endl;
    }
    //S4 devs
    if (ptrMPI->ptrSel->pDevsS4->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS4=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS4,devsS4);
        if (debug<0) cout<<cc<<endl;
    }
    //S5 devs
    if (ptrMPI->ptrSel->pDevsS5->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS5=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS5,devsS5);
        if (debug<0) cout<<cc<<endl;
    }
    //S6 devs
    if (ptrMPI->ptrSel->pDevsS6->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsS6=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS6,devsS6);
        if (debug<0) cout<<cc<<endl;
    }
    //capture rate devs
    if (ptrMPI->ptrFsh->pDevsLnC->getSize()){
        if (debug<0) cout<<tb<<tb<<"pDevsLnC=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnC,devsLnC);
        if (debug<0) cout<<cc<<endl;
    }
    
    if (debug<0) cout<<tb<<"NULL)"<<endl;//end of devs penalties list
    if (debug<0) cout<<")"<<cc<<endl;//end of penalties list
    
    //Apply penalties to F-devs
    {
        if (debug<0) cout<<"penFDevs=list("<<endl;
        double penWgt = 1.0/log(1.0+square(ptrMOs->cvFDevsPen));
        if (ptrMOs->phsDecrFDevsPen<=current_phase()){
            double scl = max(1.0,(double)(ptrMOs->phsZeroFDevsPen-ptrMOs->phsDecrFDevsPen));
            penWgt *= max(0.0,(double)(ptrMOs->phsZeroFDevsPen-current_phase()))/scl;
        }
        double effCV = std::numeric_limits<double>::infinity();
        if (penWgt>0) {
            effCV = sqrt(mfexp(1.0/penWgt)-1.0);
            if (debug<0) rpt::echo<<"phase: "<<current_phase()<<"; penWgt = "<<penWgt<<"; effCV = "<<effCV<<endl;
        } else {
            if (debug<0) rpt::echo<<"phase: "<<current_phase()<<"; penWgt = "<<penWgt<<"; effCV = Inf"<<endl;
        }
        for (int i=pDevsLnC.indexmin();i<=pDevsLnC.indexmax();i++){
            dvariable fpen = 0.5*norm2(pDevsLnC(i));
            objFun += penWgt*fpen;
            if (debug<0) {
                rpt::echo<<tb<<i<<": pen="<<fpen<<cc<<" objfun="<<penWgt*fpen<<endl;
                if (penWgt>0) {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV="<<effCV<<cc<<"pen="<<fpen<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                } else {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV=Inf"    <<cc<<"pen="<<fpen<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                }
            }
        }
        if (debug<0) cout<<tb<<tb<<"NULL)"<<endl;//end of penFDevs lists
    }
    
    if (!debug) testNaNs(value(objFun),"in calcPenalties()");
    
    if (debug>=dbgObjFun) cout<<"Finished calcPenalties()"<<endl;

//-------------------------------------------------------------------------------------
//Calculate penalties on final values in devs vectors
FUNCTION void calcDevsPenalties(int debug, ostream& cout, double penWgt, param_init_bounded_vector_vector& pDevs, dvar_matrix devs)    
    double scale = 1.0e-2;//scale for posfun calculation
    int idx;
    double lower; double upper;
    dvariable fPenLower;
    dvariable fPenUpper;
    if (debug<0) cout<<"list(";//start of list
    for (int i=pDevs.indexmin();i<=pDevs.indexmax();i++){
        if (pDevs(i).get_phase_start()){
            idx = devs(i).indexmax();//index of final dev
            lower = pDevs(i).get_minb();
            upper = pDevs(i).get_maxb();
            fPenLower.initialize();
            fPenUpper.initialize();
            posfun2(devs(i,idx)-lower,scale,fPenLower);
            posfun2(upper-devs(i,idx),scale,fPenUpper);
            objFun += penWgt*(fPenLower+fPenUpper);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"pen="<<fPenLower+fPenUpper<<cc<<"objfun="<<penWgt*(fPenLower+fPenUpper)<<cc<<
                                                          "val="<<devs(i,idx)<<cc<<"minb="<<lower<<cc<<"maxb="<<upper<<cc<<"fPenL="<<fPenLower<<cc<<"fPenU="<<fPenUpper<<"),"<<endl;
        }
    }
    if (!debug) testNaNs(value(objFun),"in calcDevsPenalties()");
    if (debug<0) cout<<tb<<tb<<"NULL)";//end of penalties list
    
//-------------------------------------------------------------------------------------
//Calculate 1st differences of vector
FUNCTION dvar_vector calc1stDiffs(const dvar_vector& d)
//    cout<<"Starting calc1stDiffs"<<endl;
    RETURN_ARRAYS_INCREMENT();
    int mn = d.indexmin();
    int mx = d.indexmax();
    dvar_vector cp; cp = d;
    dvar_vector r(mn,mx-1);
    r = cp(mn+1,mx).shift(mn)-cp(mn,mx-1);
    RETURN_ARRAYS_DECREMENT();
//    cout<<"Finished calc1stDiffs"<<endl;
    return r;

//-------------------------------------------------------------------------------------
//Calculate 2nd differences of vector
FUNCTION dvar_vector calc2ndDiffs(const dvar_vector& d)
//    cout<<"Starting calc2ndDiffs"<<endl;
    RETURN_ARRAYS_INCREMENT();
    int mn = d.indexmin();
    int mx = d.indexmax();
    dvar_vector r = calc1stDiffs(calc1stDiffs(d));
    RETURN_ARRAYS_DECREMENT();
//    cout<<"Finished calc2ndDiffs"<<endl;
    return r;

//-------------------------------------------------------------------------------------
//Calculate recruitment components in the likelihood.
FUNCTION void calcNLLs_Recruitment(int debug, ostream& cout)
    if (debug>=dbgObjFun) cout<<"Starting calcNLLs_Recruitment"<<endl;
    double nllWgtRecDevs = 0.0;//TODO: read in from input file (as vector?))
    nllRecDevs.initialize();
    if (debug<0) cout<<"list("<<endl;
    if (debug<0) cout<<tb<<"recDevs=list("<<endl;
    for (int pc=1;pc<npcRec;pc++){
        nllRecDevs(pc) = 0.5*norm2(zscrDevsLnR_cy(pc));
        for (int y=mnYr;y<=mxYr;y++) if (value(stdvDevsLnR_cy(pc,y))>0) {nllRecDevs(pc) += log(stdvDevsLnR_cy(pc,y));}
        objFun += nllWgtRecDevs*nllRecDevs(pc);
        if (debug<0){
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs(pc)<<cc<<"objfun="<<nllWgtRecDevs*nllRecDevs(pc)<<cc;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<cc;
            cout<<"stdvs="; wts::writeToR(cout,value(stdvDevsLnR_cy(pc))); cout<<")"<<cc<<endl;
        }
    }//pc
    {
        int pc = npcRec;
        nllRecDevs(pc) = 0.5*norm2(zscrDevsLnR_cy(pc));
        for (int y=mnYr;y<=mxYr;y++) if (value(stdvDevsLnR_cy(pc,y))>0) {nllRecDevs(pc) += log(stdvDevsLnR_cy(pc,y));}
        objFun += nllWgtRecDevs*nllRecDevs(pc);
        if (debug<0){
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs(pc)<<cc<<"objfun="<<nllWgtRecDevs*nllRecDevs(pc)<<cc;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<cc;
            cout<<"stdvs="; wts::writeToR(cout,value(stdvDevsLnR_cy(pc))); cout<<")"<<endl;
        }
    }//pc
    if (debug<0) cout<<tb<<")";//recDevs
    if (debug<0) cout<<")";
    
    if (!debug) testNaNs(value(objFun),"in calcNLLs_Recruitment()");
    
    if (debug>=dbgObjFun) cout<<"Finished calcNLLs_Recruitment"<<endl;

//-------------------------------------------------------------------------------------
//Calculate objective function TODO: finish
FUNCTION void calcObjFun(int debug, ostream& cout)
    if ((debug>=dbgObjFun)||(debug<0)) cout<<"Starting calcObjFun"<<endl;

    objFun.initialize();
    
    //objective function penalties
    calcPenalties(debug,cout);

    //prior likelihoods
    calcAllPriors(debug,cout);

    //recruitment component
    calcNLLs_Recruitment(debug,cout);
    
    //data components
    calcNLLs_Fisheries(debug,cout);
    calcNLLs_Surveys(debug,cout);
    calcNLLs_GrowthData(debug,cout);
    //calcNLLs_ChelaHeightData(debug,cout);  //TODO: implement this!
    
    if ((debug>=dbgObjFun)||(debug<0)){
        cout<<"proc call          = "<<ctrProcCalls<<endl;
        cout<<"proc call in phase = "<<ctrProcCallsInPhase<<endl;
        cout<<"total objFun       = "<<objFun<<endl;
        cout<<"Finished calcObjFun"<<endl<<endl;
    }
    
//******************************************************************************
//* Function: void calcNLLs_GrowthData
//* 
//* Description: Calculates NLLs for growth data
//* 
//* Inputs:
//*  none
//* Returns:
//*  void
//* Alters:
//*  objFun
//******************************************************************************
FUNCTION void calcNLLs_GrowthData(int debug, ostream& cout)  
    if(debug>dbgObjFun) cout<<"Starting calcNLLs_GrowthData()"<<endl;
    
    if (debug<0) cout<<"list("<<endl;
    for (int i=0;i<ptrMDS->nGrw;i++){
        GrowthData* pGD = ptrMDS->ppGrw[i];
        d3_array gd_xcn = ptrMDS->ppGrw[i]->inpData_xcn;
        if (debug<0) cout<<"GrowthData."<<i+1<<"=list("<<endl;
        for (int x=1;x<=nSXs;x++){
            int nObs = ptrMDS->ppGrw[i]->nObs_x(x);
            if (nObs>0) {
                double wgt = ptrMDS->ppGrw[i]->llWgt;
                /* observation year */
                ivector year_n = ptrMDS->ppGrw[i]->obsYears_xn(x);
                /* pre-molt size, by observation */
                dvar_vector zpre_n = ptrMDS->ppGrw[i]->inpData_xcn(x,2);
                /* post-molt size, by observation */
                dvar_vector zpst_n = ptrMDS->ppGrw[i]->inpData_xcn(x,3);
                /* molt increment, by observation */
                dvar_vector incZ_n = zpst_n - zpre_n;
                /* mean post-molt size, by observation */
                dvar_vector mnZ_n = elem_prod(mfexp(grA_xy(x)(year_n)),pow(zpre_n,grB_xy(x)(year_n)));
                /* multiplicative scale factor, by observation */
                dvar_vector ibeta_n = 1.0/grBeta_xy(x)(year_n);
                /* location factor, by observation */
                dvar_vector alpha_n = elem_prod(mnZ_n-zpre_n,ibeta_n);
                dvar_vector nlls_n(1,nObs); nlls_n.initialize();
                nlls_n = -wts::log_gamma_density(incZ_n,alpha_n,ibeta_n);
                dvariable nll = sum(nlls_n);
                if (isnan(value(nll))){
                    ofstream os("GrowthData.NLLs.NanReport.dat");
                    os<<"phase = "<<current_phase()<<endl;
                    os<<"sex   = "<<tcsam::getSexType(x)<<endl;
                    os<<"nll   = "<<nll<<endl;
                    os<<"nObs  = "<<nObs<<endl;
                    os<<"years   = "<<year_n<<endl;
                    os<<"grA     = "<<grA_xy(x)(year_n)<<endl;
                    os<<"grB     = "<<grB_xy(x)(year_n)<<endl;
                    os<<"zpre_n  = "<<zpre_n<<endl;
                    os<<"zpst_n  = "<<zpst_n<<endl;
                    os<<"mnZ_n   = "<<mnZ_n<<endl;
                    os<<"incZ    = "<<zpst_n-zpre_n<<endl;
                    os<<"mnInc   = "<<mnZ_n-zpre_n<<endl;
                    os<<"ibeta_n = "<<ibeta_n<<endl;
                    os<<"alpha_n = "<<alpha_n<<endl;
                    os<<"nlls_n  = "<<nlls_n<<endl;
                    dvar_vector zscrs = elem_div((zpst_n-mnZ_n),sqrt(elem_prod(mnZ_n,grB_xy(x)(year_n))));
                    os<<"zscrs   = "<<zscrs<<endl;
                    os.close();
                    exit(-1);
                }
                objFun += wgt*nll;
                if (debug<0) {
                    dvar_vector zscrs = elem_div((zpst_n-mnZ_n),sqrt(elem_prod(mnZ_n,grB_xy(x)(year_n))));
                    cout<<tcsam::getSexType(x)<<"=list(type='normal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl;
                    cout<<"years="; wts::writeToR(cout,year_n);         cout<<cc<<endl;
                    cout<<"zPre=";  wts::writeToR(cout,value(zpre_n));  cout<<cc<<endl;
                    cout<<"zPst=";  wts::writeToR(cout,value(zpst_n));  cout<<cc<<endl;
                    cout<<"grA=";   wts::writeToR(cout,value(grA_xy(x)(year_n))); cout<<cc<<endl;
                    cout<<"grB=";   wts::writeToR(cout,value(grB_xy(x)(year_n))); cout<<cc<<endl;
                    cout<<"mnZ =";  wts::writeToR(cout,value(mnZ_n));   cout<<cc<<endl;
                    cout<<"ibeta="; wts::writeToR(cout,value(ibeta_n)); cout<<cc<<endl;
                    cout<<"alpha="; wts::writeToR(cout,value(alpha_n)); cout<<cc<<endl;
                    cout<<"nlls=";  wts::writeToR(cout,value(nlls_n));  cout<<cc<<endl;
                    cout<<"zscrs="; wts::writeToR(cout,value(zscrs));   cout<<"),"<<endl;
                }
            }//nObs>0
        }//x
        if (debug<0) cout<<"NULL),";
    }//datasets
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>dbgObjFun) cout<<"finished calcNLLs_GrowthData()"<<endl;

//-------------------------------------------------------------------------------------
FUNCTION void testNaNs(double v, adstring str) 
    if (isnan(v)){
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"----NaN detected: "<<str<<"---"<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        ofstream os("NaNReport.rep");
        os<<"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        os<<"#----NaN detected: "<<str<<"---"<<endl;
        os<<"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        updateMPI(0,cout);
        os<<"nanRep=list(phase="<<current_phase()<<",pcCtr="<<ctrProcCalls<<cc<<",pcCtrInPhs="<<ctrProcCallsInPhase<<cc<<endl;
        ReportToR_Params(os,0,cout);         os<<","<<endl;
        ReportToR_ModelProcesses(os,0,cout); os<<","<<endl;
        ReportToR_ModelResults(os,0,cout);   os<<","<<endl;
        ReportToR_ModelFits(os,-1,cout);     os<<endl;
        os<<")"<<endl;
        os.close();
        exit(-1);
    }

//-------------------------------------------------------------------------------------
//Calculate norm2 NLL contribution to objective function
FUNCTION void calcNorm2NLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcNorm2NLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    for (int i=1;i<=yrs.size();i++){
        y = yrs(i);
        if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
            zscr(y) = (obs[i]-mod[y]);
        }
    }
    nll += 0.5*norm2(zscr);
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='norm2',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNorm2NLL()"<<endl;
    
//-------------------------------------------------------------------------------------
//Calculate normal NLL contribution to objective function
FUNCTION void calcNormalNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
   if (debug>=dbgAll) cout<<"Starting calcNormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (obs[i]-mod[y])/stdv[i];
                //  nll += log(stdv[i]);
            }
        }
    }
    nll += 0.5*norm2(zscr);
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='normal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcNormalNLL()"<<endl;
    
//-------------------------------------------------------------------------------------
//Calculate lognormal NLL contribution to objective function
FUNCTION void calcLognormalNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcLognormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (log(obs[i]+smlVal)-log(mod[y]+smlVal))/stdv[i];
                //nll += log(stdv[i]);
            }
        }
    }
    nll += 0.5*norm2(zscr);
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='lognormal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcLognormalNLL()"<<endl;
    
//-------------------------------------------------------------------------------------
//Pass-through for no NLL contribution to objective function for aggregated catch
FUNCTION void calcNoneNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcNoneNLL(agg catch)"<<endl;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='none',wgt="<<0<<cc<<"nll="<<0<<cc<<"objfun="<<0<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNoneNLL(agg catch)"<<endl;
    
//-------------------------------------------------------------------------------------
//Calculate multinomial NLL contribution to objective function
FUNCTION void calcMultinomialNLL(double wgt, dvar_vector& mod, dvector& obs, double& ss, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcMultinomialNLL()"<<endl;
    dvariable nll = -ss*(obs*(log(mod+smlVal)-log(obs+smlVal)));//note dot-product sums
    objFun += wgt*nll;
    if (debug<0){
        dvector vmod = value(mod);
        dvector nlls = -ss*(elem_prod(obs,log(vmod+smlVal)-log(obs+smlVal)));
        dvector zscrs = elem_div(obs-vmod,sqrt(elem_prod((vmod+smlVal),1.0-(vmod+smlVal))/ss));//pearson residuals
        double effN = 0.0;
        if (ss>0) effN = (vmod*(1.0-vmod))/norm2(obs-vmod);
        cout<<"list(nll.type='multinomial',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<"ss="<<ss<<cc<<"effN="<<effN<<cc<<endl; 
        adstring dzbs = "size=c("+ptrMC->csvZBs+")";
        cout<<"nlls=";  wts::writeToR(cout,nlls, dzbs); cout<<cc<<endl;
        cout<<"obs=";   wts::writeToR(cout,obs,  dzbs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,vmod, dzbs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,zscrs,dzbs); cout<<endl;
        cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcMultinomialNLL()"<<endl;
 
//-------------------------------------------------------------------------------------
//Pass through for no NLL contribution to objective function for size comps
FUNCTION void calcNoneNLL(double wgt, dvar_vector& mod, dvector& obs, double& ss, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcNoneNLL(size comps)"<<endl;
    if (debug<0){
        cout<<"list(nll.type='none',wgt="<<wgt<<cc<<"nll="<<0.0<<cc<<"objfun="<<0.0<<cc<<"ss="<<0<<cc<<"effN="<<0<<cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNoneNLL(size comps)"<<endl;
 
//-------------------------------------------------------------------------------------
//Calculate time series contribution to objective function
FUNCTION void calcNLL(int llType, double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
    switch (llType){
        case tcsam::LL_NONE:
            calcNoneNLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_LOGNORMAL:
            calcLognormalNLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORMAL:
            calcNormalNLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        case tcsam::LL_NORM2:
            calcNorm2NLL(wgt,mod,obs,stdv,yrs,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(1)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }    

//-------------------------------------------------------------------------------------
//Calculate size frequency contribution to objective function
FUNCTION void calcNLL(int llType, double wgt, dvar_vector& mod, dvector& obs, double& ss, int debug, ostream& cout)
    switch (llType){
        case tcsam::LL_NONE:
            calcNoneNLL(wgt,mod,obs,ss,debug,cout);
            break;
        case tcsam::LL_MULTINOMIAL:
            calcMultinomialNLL(wgt,mod,obs,ss,debug,cout);
            break;
        default:
            cout<<"Unrecognized likelihood type in calcNLL(2)"<<endl;
            cout<<"Input type was "<<llType<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
    }
   
//-------------------------------------------------------------------------------------
//Calculate aggregate catch (abundance or biomass) components to objective function
FUNCTION void calcNLLs_AggregateCatch(AggregateCatchData* ptrAB, dvar5_array& mA_yxmsz, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcNLLs_AggregateCatch()"<<endl;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvar_vector tAB_y(mny,mxy);
    int isBio = ptrAB->type==AggregateCatchData::KW_BIOMASS_DATA;
    if (debug>=dbgAll) cout<<"isBio="<<isBio<<tb<<"type="<<ptrAB->type<<endl;
    if (debug<0) cout<<"list(fit.type='"<<tcsam::getFitType(ptrAB->optFit)<<"',fits=list("<<endl;
    if (ptrAB->optFit==tcsam::FIT_BY_TOT){
        tAB_y.initialize();
        if (isBio){
            for (int x=1;x<=nSXs;x++) {
                for (int m=1;m<=nMSs;m++) {
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                }//m
            }//x
        } else {
            for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y));//sum over x,m,s,z
        }
        if (debug>=dbgAll) cout<<"FIT_BY_TOT: "<<tAB_y<<endl;
        if (debug<0) {
            cout<<"list(";
            cout<<"x="<<qt<<tcsam::getSexType(ALL_SXs)     <<qt<<cc;
            cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
            cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)   <<qt<<cc;
            cout<<"nll=";
        }
        calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(ALL_SXs,ALL_MSs,ALL_SCs), ptrAB->sd_xmsy(ALL_SXs,ALL_MSs,ALL_SCs), ptrAB->yrs, debug, cout);                
        if (debug<0) cout<<")";
        if (debug<0) cout<<")";
    } else if (ptrAB->optFit==tcsam::FIT_BY_X){
        for (int x=1;x<=nSXs;x++){
            tAB_y.initialize();
            if (isBio){
                for (int m=1;m<=nMSs;m++) {
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                }//m
            } else {
                for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x));//sum over m,s,z
            }
            if (debug>=dbgAll) cout<<"FIT_BY_X("<<x<<"): "<<tAB_y<<endl;
            if (debug<0) {
                cout<<"list(";
                cout<<"x="<<qt<<tcsam::getSexType(x)           <<qt<<cc;
                cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
                cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)   <<qt<<cc;
                cout<<"nll=";
            }
            calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,ALL_MSs,ALL_SCs), ptrAB->sd_xmsy(x,ALL_MSs,ALL_SCs), ptrAB->yrs, debug, cout); 
            if (debug<0) cout<<"),"<<endl;
        }//x
        if (debug<0) cout<<"NULL)";
    } else if ((ptrAB->optFit==tcsam::FIT_BY_XM)||(ptrAB->optFit==tcsam::FIT_BY_X_MATONLY)){
        int mnM = IMMATURE;
        if (ptrAB->optFit==tcsam::FIT_BY_X_MATONLY) mnM = MATURE;
        for (int x=1;x<=nSXs;x++){
            for (int m=mnM;m<=nMSs;m++){
                tAB_y.initialize();
                if (isBio){
                    if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                    for (int s=1;s<=nSCs;s++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//s
                } else {
                    for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x,m));//sum over s,z
                }
                if (debug<0) {
                    cout<<"list(";
                    cout<<"x="<<qt<<tcsam::getSexType(x)        <<qt<<cc;
                    cout<<"m="<<qt<<tcsam::getMaturityType(m)   <<qt<<cc;
                    cout<<"s="<<qt<<tcsam::getShellType(ALL_SCs)<<qt<<cc;
                    cout<<"nll=";
                }
                calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,m,ALL_SCs), ptrAB->sd_xmsy(x,m,ALL_SCs), ptrAB->yrs, debug, cout); 
                if (debug<0) cout<<"),"<<endl;
            }//m
        }//x
        if (debug<0) cout<<"NULL)";
    } else if (ptrAB->optFit==tcsam::FIT_BY_XS){
        for (int x=1;x<=nSXs;x++){
            for (int s=1;s<=nSCs;s++){
                tAB_y.initialize();
                if (isBio){
                    for (int m=1;m<=nMSs;m++) {
                        if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                        }//y
                    }//m
                } else {
                    for (int m=1;m<=nMSs;m++) {
                        for (int y=mny;y<=mxy;y++) {
                            tAB_y(y) += sum(mA_yxmsz(y,x,m,s));//sum over m,z
                        }//y
                    }//m
                }
                if (debug<0) {
                    cout<<"list(";
                    cout<<"x="<<qt<<tcsam::getSexType(x)           <<qt<<cc;
                    cout<<"m="<<qt<<tcsam::getMaturityType(ALL_MSs)<<qt<<cc;
                    cout<<"s="<<qt<<tcsam::getShellType(s)         <<qt<<cc;
                    cout<<"nll=";
                }
                calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,ALL_MSs,s), ptrAB->sd_xmsy(x,ALL_MSs,s), ptrAB->yrs, debug, cout); 
                if (debug<0) cout<<"),"<<endl;
            }//s
        }//x
        if (debug<0) cout<<"NULL)";
    } else if (ptrAB->optFit==tcsam::FIT_BY_XMS){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    tAB_y.initialize();
                    if (isBio){
                        if (debug>=dbgAll) cout<<"w("<<x<<cc<<m<<") = "<<ptrMDS->ptrBio->wAtZ_xmz(x,m)<<endl;
                        for (int y=mny;y<=mxy;y++) tAB_y(y) += mA_yxmsz(y,x,m,s)*ptrMDS->ptrBio->wAtZ_xmz(x,m);
                    } else {
                        for (int y=mny;y<=mxy;y++) tAB_y(y) += sum(mA_yxmsz(y,x,m,s));//sum over z
                    }
                    if (debug<0) {
                        cout<<"list(";
                        cout<<"x="<<qt<<tcsam::getSexType(x)     <<qt<<cc;
                        cout<<"m="<<qt<<tcsam::getMaturityType(m)<<qt<<cc;
                        cout<<"s="<<qt<<tcsam::getShellType(s)   <<qt<<cc;
                        cout<<"nll=";
                    }
                    calcNLL(ptrAB->llType, ptrAB->llWgt, tAB_y, ptrAB->C_xmsy(x,m,s), ptrAB->sd_xmsy(x,m,s), ptrAB->yrs, debug, cout); 
                    if (debug<0) cout<<"),"<<endl;
                }//s
            }//m
        }//x
        if (debug<0) cout<<"NULL)";
    } else {
        std::cout<<"Calling calcNLLs_AggregateCatch with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrAB->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<")";
    if (debug>=dbgAll){
        cout<<"Finished calcNLLs_AggregateCatch()"<<endl;
    }

//-------------------------------------------------------------------------------------
//Calculate catch size frequencies components to objective function
FUNCTION void calcNLLs_CatchNatZ(SizeFrequencyData* ptrZFD, dvar5_array& mA_yxmsz, int debug, ostream& cout)
    if (debug>=dbgNLLs) cout<<"Starting calcNLLs_CatchNatZ()"<<endl;
    if (ptrZFD->optFit==tcsam::FIT_NONE) return;
    ivector yrs = ptrZFD->yrs;
    int y;
    double ss;
    dvariable nT;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvector     oP_z;//observed size comp.
    dvar_vector mP_z;//model size comp.
    if (ptrZFD->optFit==tcsam::FIT_BY_XE){
        oP_z.allocate(1,nSXs*nZBs);
        mP_z.allocate(1,nSXs*nZBs);
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_X_ME){
        oP_z.allocate(1,nMSs*nZBs);
        mP_z.allocate(1,nMSs*nZBs);
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X_SE){
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_XME){
        oP_z.allocate(1,nSXs*nMSs*nZBs);
        mP_z.allocate(1,nSXs*nMSs*nZBs);
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XM_SE){
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
    } else 
    {
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
    }
    if (debug<0) cout<<"list("<<endl;
    for (int iy=1;iy<=yrs.size();iy++) {
        y = yrs[iy];
        if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
        if ((mny<=y)&&(y<=mxy)) {
            if (ptrZFD->optFit==tcsam::FIT_BY_TOT){
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=ALL_SXs;x++){
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                    }
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    for (int x=1;x<=nSXs;x++){
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++)mP_z += mA_yxmsz(y,x,m,s);
                        }
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    if (debug<0) {
                        cout<<"'"<<y<<"'=list(";
                        cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                        cout<<"y="<<y<<cc;
                        cout<<"x='"<<tcsam::getSexType(ALL_SXs)<<"'"<<cc;
                        cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                        cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                        cout<<"fit=";
                    }
                    calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,debug,cout);
                    if (debug<0) cout<<")"<<cc<<endl;
                }
                //FIT_BY_TOT
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_X){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"y="<<y<<cc;
                            cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                            cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                            cout<<"fit=";
                        }
                        calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//nT>0
                }//x
                //FIT_BY_X
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XM){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        ss = 0;
                        nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int s=1;s<=nSCs;s++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XM
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XS){
                for (int x=1;x<=nSXs;x++) {
                    for (int s=1;s<=nSCs;s++){
                        ss = 0;
                        nT.initialize();
                        for (int m=1;m<=nMSs;m++) nT += sum(mA_yxmsz(y,x,m,s));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int m=1;m<=ALL_MSs;m++) {
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            for (int m=1;m<=nMSs;m++) mP_z += mA_yxmsz(y,x,m,s);
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XS
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XMS){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) {
                            ss = 0;
                            nT = sum(mA_yxmsz(y,x,m,s));//=0 if not calculated
                            if (value(nT)>0){
                                oP_z.initialize();//observed size comp.
                                mP_z.initialize();//model size comp.                            
                                ss   += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                if (sum(oP_z)>0) oP_z /= sum(oP_z);
                                if (debug>=dbgNLLs){
                                    cout<<"ss = "<<ss<<endl;
                                    cout<<"oP_Z = "<<oP_z<<endl;
                                }
                                mP_z += mA_yxmsz(y,x,m,s);
                                mP_z /= nT;//normalize model size comp
                                if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"y="<<y<<cc;
                                    cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//nT>0
                        }//s
                    }//m
                }//x
                //FIT_BY_XMS
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XE){
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        for (int m=1;m<=ALL_MSs;m++) {
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                        }
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                        }
                    }//x
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        int mnz = 1+(x-1)*nZBs;
                        int mxz = x*nZBs;
                        dvar_vector mPt = mP_z(mnz,mxz);
                        dvector oPt = oP_z(mnz,mxz);
                        if (debug<0) {
                            cout<<"'"<<y<<"'=list(";
                            cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                            cout<<"y="<<y<<cc;
                            cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                            cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                            cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                            cout<<"fit=";
                        }
                        calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//x
                }//nT>0
                //FIT_BY_XE
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_X_ME){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs;
                            int mxz = m*nZBs;
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }//s
                            for (int s=1;s<=nSCs;s++) {
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//s
                        }//m
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs;
                            int mxz = m*nZBs;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//nT>0
                }//x
                //FIT_BY_X_ME
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_X_SE){
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int s=1;s<=nSCs;s++) {
                            int mnz = 1+(s-1)*nZBs;
                            int mxz = s*nZBs;
                            for (int m=1;m<=ALL_MSs;m++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }//m
                            for (int m=1;m<=nMSs;m++) {
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//m
                        }//s
                        if (sum(oP_z)>0) oP_z /= sum(oP_z);
                        if (debug>=dbgNLLs){
                            cout<<"ss = "<<ss<<endl;
                            cout<<"oP_Z = "<<oP_z<<endl;
                        }
                        mP_z /= nT;//normalize model size comp
                        if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                        for (int s=1;s<=nSCs;s++) {
                            int mnz = 1+(s-1)*nZBs;
                            int mxz = s*nZBs;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(ALL_MSs)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//s
                    }//nT>0
                }//x
                //FIT_BY_X_SE
            } else
            if (ptrZFD->optFit==tcsam::FIT_BY_XME){
                ss = 0;
                nT = sum(mA_yxmsz(y));//=0 if not calculated
                if (value(nT)>0){
                    oP_z.initialize();//observed size comp.
                    mP_z.initialize();//model size comp.
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            for (int s=1;s<=ALL_SCs;s++) {
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                            }
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            if (m<=nMSs) {for (int s=1;s<=nSCs;s++) mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);}
                        }//m
                    }//x
                    if (sum(oP_z)>0) oP_z /= sum(oP_z);
                    if (debug>=dbgNLLs){
                        cout<<"ss = "<<ss<<endl;
                        cout<<"oP_Z = "<<oP_z<<endl;
                    }
                    mP_z /= nT;//normalize model size comp
                    if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                    for (int x=1;x<=nSXs;x++) {
                        for (int m=1;m<=nMSs;m++) {
                            int mnz = 1+(m-1)*nZBs+(x-1)*nMSs*nZBs;
                            int mxz = mnz+nZBs-1;
                            dvar_vector mPt = mP_z(mnz,mxz);
                            dvector oPt = oP_z(mnz,mxz);
                            if (debug<0) {
                                cout<<"'"<<y<<"'=list(";
                                cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                cout<<"y="<<y<<cc;
                                cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                cout<<"s='"<<tcsam::getShellType(ALL_SCs)<<"'"<<cc;
                                cout<<"fit=";
                            }
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//x
                }//nT>0
                //FIT_BY_XME
            } else 
            if (ptrZFD->optFit==tcsam::FIT_BY_XM_SE){
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++) {
                        ss = 0;
                        nT = sum(mA_yxmsz(y,x,m));//=0 if not calculated
                        if (value(nT)>0){
                            oP_z.initialize();//observed size comp.
                            mP_z.initialize();//model size comp.
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs;
                                int mxz = s*nZBs;
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
                                mP_z(mnz,mxz).shift(1) += mA_yxmsz(y,x,m,s);
                            }//s
                            if (sum(oP_z)>0) oP_z /= sum(oP_z);
                            if (debug>=dbgNLLs){
                                cout<<"ss = "<<ss<<endl;
                                cout<<"oP_Z = "<<oP_z<<endl;
                            }
                            mP_z /= nT;//normalize model size comp
                            if (debug>=dbgNLLs) cout<<"mP_z = "<<mP_z<<endl;
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs;
                                int mxz = s*nZBs;
                                dvar_vector mPt = mP_z(mnz,mxz);
                                dvector oPt = oP_z(mnz,mxz);
                                if (debug<0) {
                                    cout<<"'"<<y<<"'=list(";
                                    cout<<"fit.type='"<<tcsam::getFitType(ptrZFD->optFit)<<"'"<<cc;
                                    cout<<"y="<<y<<cc;
                                    cout<<"x='"<<tcsam::getSexType(x)<<"'"<<cc;
                                    cout<<"m='"<<tcsam::getMaturityType(m)<<"'"<<cc;
                                    cout<<"s='"<<tcsam::getShellType(s)<<"'"<<cc;
                                    cout<<"fit=";
                                }
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//s
                        }//nT>0
                    }//m
                }//x
                //FIT_BY_XM_SE
            } else 
            {
                std::cout<<"Calling calcNLLs_CatchNatZ with invalid fit option."<<endl;
                std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrZFD->optFit)<<qt<<endl;
                std::cout<<"Aborting..."<<endl;
                exit(-1);
            }
        } //if ((mny<=y)&&(y<=mxy))
    } //loop over iy
    if (debug<0) cout<<"NULL)";
    if (debug>=dbgNLLs) cout<<"Finished calcNLLs_CatchNatZ()"<<endl;

//-------------------------------------------------------------------------------------
//Calculate fishery components to objective function
FUNCTION void calcNLLs_Fisheries(int debug, ostream& cout)
//    if (debug>0) debug = dbgAll+10;
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Fisheries()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int f=1;f<=nFsh;f++){
        if (debug>=dbgAll) cout<<"calculating NLLs for fishery "<<ptrMC->lblsFsh[f]<<endl;
        if (debug<0) cout<<ptrMC->lblsFsh[f]<<"=list("<<endl;
        FleetData* ptrObs = ptrMDS->ppFsh[f-1];
        if (ptrObs->hasRCD){//retained catch data
            if (debug<0) cout<<"retained.catch=list("<<endl;
            if (ptrObs->ptrRCD->hasN && ptrObs->ptrRCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---retained catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrN,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasB && ptrObs->ptrRCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---retained catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrB,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---retained catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrRCD->ptrZFD,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasTCD){//observed total catch data
            if (debug<0) cout<<"total.catch=list("<<endl;
            if (ptrObs->ptrTCD->hasN && ptrObs->ptrTCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---total catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrN,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasB && ptrObs->ptrTCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---total catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrB,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---total catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrTCD->ptrZFD,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (ptrObs->hasDCD){//observed discard catch data
            if (debug<0) cout<<"discard.catch=list("<<endl;
            if (ptrObs->ptrDCD->hasN && ptrObs->ptrDCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---discard catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrN,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasB && ptrObs->ptrDCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---discard catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrB,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---discard catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrDCD->ptrZFD,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//fisheries
    if (debug<0) cout<<"NULL)"<<endl;
    if (!debug) testNaNs(value(objFun),"in calcNLLs_Fisheries()");
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Fisheries()"<<endl;

//-------------------------------------------------------------------------------------
//Calculate survey components to objective function
FUNCTION void calcNLLs_Surveys(int debug, ostream& cout)
//    if (debug>0) debug = dbgAll+10;
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Surveys()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int v=1;v<=nSrv;v++){
        if (debug>=dbgAll) cout<<"calculating NLLs for survey "<<ptrMC->lblsSrv[v]<<endl;
        if (debug<0) cout<<ptrMC->lblsSrv[v]<<"=list("<<endl;
        FleetData* ptrObs = ptrMDS->ppSrv[v-1];
        if (ptrObs->hasICD){//index catch data
            if (debug<0) cout<<"index.catch=list("<<endl;
            if (ptrObs->ptrICD->hasN && ptrObs->ptrICD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---index catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrN,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasB && ptrObs->ptrICD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---index catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrB,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---index catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                calcNLLs_CatchNatZ(ptrObs->ptrICD->ptrZFD,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (debug<0) cout<<"NULL),"<<endl;
        }
        if (debug<0) cout<<"NULL),"<<endl;
    }//surveys loop
    if (debug<0) cout<<"NULL)"<<endl;
    if (!debug) testNaNs(value(objFun),"in calcNLLs_Surveys()");
    if (debug>=dbgAll) cout<<"Finished calcNLLs_Surveys()"<<endl;

//-------------------------------------------------------------------------------------
//Calculate contributions to objective function from all priors                                         
FUNCTION void calcAllPriors(int debug, ostream& cout)
    if (debug>=dbgPriors) cout<<"Starting calcAllPriors()"<<endl;
    if (debug<0) cout<<"list("<<endl;

    //recruitment parameters
    if (debug<0) cout<<tb<<"recruitment=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnR,  pLnR,  debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRCV,pLnRCV,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLgtRX,pLgtRX,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRa, pLnRa, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnRb, pLnRb, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pDevsLnR,devsLnR,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
   
    //natural mortality parameters
    if (debug<0) cout<<tb<<"'natural mortality'=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnM,   pLnM,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMT, pLnDMT, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMX, pLnDMX, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMM, pLnDMM, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pLnDMXM,pLnDMXM,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //growth parameters
    if (debug<0) cout<<tb<<"growth=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pLnGrA,   pLnGrA,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pLnGrB,   pLnGrB,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pLnGrBeta,pLnGrBeta,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //maturity parameters
    if (debug<0) cout<<tb<<"maturity=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrM2M->pLgtPrM2M,pLgtPrM2M,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //selectivity parameters
    if (debug<0) cout<<tb<<"'selectivity functions'=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS1,pS1,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS2,pS2,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS3,pS3,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS4,pS4,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS5,pS5,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pS6,pS6,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS1,devsS1,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS2,devsS2,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS3,devsS3,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS4,devsS4,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS5,devsS5,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSel->pDevsS6,devsS6,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //fishing mortality parameters
    if (debug<0) cout<<tb<<"fisheries=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnC,    pLnC,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCT,  pLnDCT, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCX,  pLnDCX, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCM,  pLnDCM, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnDCXM, pLnDCXM,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDevsLnC,devsLnC,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
   
    //survey catchability parameters
    if (debug<0) cout<<tb<<"surveys=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnQ,    pLnQ,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQT,  pLnDQT, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQX,  pLnDQX, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQM,  pLnDQM, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pLnDQXM, pLnDQXM,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<endl;
    
    if (debug<0) cout<<")"<<endl;
    
    if (!debug) testNaNs(value(objFun),"in calcAllPriors()");
    if (debug>=dbgPriors) cout<<"Finished calcAllPriors()"<<endl;

//-------------------------------------------------------------------------------------
//Write data to file as R list
FUNCTION void ReportToR_Data(ostream& os, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_Data(...)"<<endl;
    ptrMDS->writeToR(os,"data",0);
    if (debug) cout<<"Finished ReportToR_Data(...)"<<endl;

//-------------------------------------------------------------------------------------
//Write parameter values to file as R list
FUNCTION void ReportToR_Params(ostream& os, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_Params(...)"<<endl;
    ptrMPI->writeToR(os);
    if (debug) cout<<"Finished ReportToR_Params(...)"<<endl;

//-------------------------------------------------------------------------------------
//Write model processes to file as R list
FUNCTION void ReportToR_ModelProcesses(ostream& os, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_ModelProcesses(...)"<<endl;
//    wts::adstring_matrix aM2M = tcsam::convertPCs(ptrMPI->ptrM2M);
//    wts::adstring_matrix aGr  = tcsam::convertPCs(ptrMPI->ptrGrw);
    os<<"mp=list("<<endl;
        os<<"M_cxm    ="; wts::writeToR(os,value(M_cxm),   adstring("pc=1:"+str(npcNM )),xDms,mDms);   os<<cc<<endl;
        os<<"prM2M_cz ="; wts::writeToR(os,value(prM2M_cz),adstring("pc=1:"+str(npcM2M)),zbDms);       os<<cc<<endl;
        os<<"sel_cz   ="; wts::writeToR(os,value(sel_cz),  adstring("pc=1:"+str(npcSel)),zbDms);       os<<cc<<endl;
        os<<"M_yxmsz  ="; wts::writeToR(os,wts::value(M_yxmsz),   yDms,xDms,mDms,sDms,zbDms);  os<<cc<<endl;
        os<<"prM2M_yxz =";wts::writeToR(os,     value(prM2M_yxz), yDms,xDms,zbDms);            os<<cc<<endl;
        os<<"T_list=list("<<endl;
            os<<"mnZAM_cz   ="; wts::writeToR(os,value(mnGrZ_cz),adstring("pc=1:"+str(npcGrw )),zbDms);       os<<cc<<endl;
            os<<"T_czz      ="; wts::writeToR(os,value(prGr_czz),adstring("pc=1:"+str(npcGrw )),zbDms,zpDms); os<<cc<<endl;
            os<<"mnZAM_yxsz =";  wts::writeToR(os,wts::value(mnGrZ_yxsz),yDms,xDms,sDms,zbDms);              os<<cc<<endl;
            if (0){//TDODO: develop option to output
                os<<"T_yxszz    =";  wts::writeToR(os,wts::value(prGr_yxszz),yDms,xDms,sDms,zbDms,zpDms);    os<<endl;
            } else {
                os<<"T_yxszz    =NULL"<<endl;
            }
        os<<")"<<cc<<endl;
        os<<"R_list=list("<<endl;
            os<<"R_y  ="; wts::writeToR(os,value(R_y), yDms);                                os<<cc<<endl;
            os<<"R_yx ="; wts::writeToR(os,value(R_yx),yDms,xDms);                           os<<cc<<endl;
            os<<"R_yz ="; wts::writeToR(os,value(R_yz),yDms,zbDms);                          os<<cc<<endl;
            os<<"Rx_c ="; wts::writeToR(os,value(Rx_c),adstring("pc=1:"+str(npcRec)));       os<<cc<<endl;
            os<<"R_cz ="; wts::writeToR(os,value(R_cz),adstring("pc=1:"+str(npcRec)),zbDms); os<<endl;
        os<<")"<<cc<<endl;
        d6_array tmF_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
        for (int f=1;f<=nFsh;f++){
            for (int y=mnYr;y<=mxYr;y++){
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++){
                            tmF_fyxmsz(f,y,x,m,s) = value(rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s));
                        }
                    }
                }
            }
        }
        os<<"Eff_list=list("<<endl;
        for (int f=1;f<=nFsh;f++){
            int fd = mapM2DFsh(f);//index of corresponding fishery data object
            os<<"`"<<ptrMDS->ppFsh[fd-1]->name<<"`=list("<<endl;
            if (ptrMDS->ppFsh[fd-1]->ptrEff) {
                adstring yDmsp = "y="+str(eff_fy(f).indexmin())+":"+str(eff_fy(f).indexmax());
                os<<"optFcAvg="<<optsFcAvg(f)<<cc;
                os<<"avgEff="<<avgEff(f)<<cc<<endl;
                os<<"eff_y ="; wts::writeToR(os,eff_fy(f),yDmsp); os<<cc<<endl;
                os<<"avgRatioFc2Eff ="; wts::writeToR(os,     value(avgRatioFc2Eff(f)),xDms,mDms,sDms); os<<cc<<endl;
                os<<"cpF_xmsy =";       wts::writeToR(os,wts::value(cpF_fxmsy(f)),     xDms,mDms,sDms,yDmsp); os<<cc<<endl;
                os<<"avgFc_xms =";      wts::writeToR(os,     value(avgFc_fxms(f)),    xDms,mDms,sDms); os<<endl;
            }
            os<<")"<<cc<<endl;
        }
        os<<"NULL)"<<cc<<endl;
        os<<"F_list=list("<<endl;
            os<<"hm_fy     ="; wts::writeToR(os,     value(hmF_fy),    fDms,yDms);                     os<<cc<<endl;
            os<<"cpF_fyxms ="; wts::writeToR(os,wts::value(cpF_fyxms ),fDms,yDms,xDms,mDms,sDms);      os<<cc<<endl;
            os<<"sel_fyxmsz="; wts::writeToR(os,wts::value(sel_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"ret_fyxmsz="; wts::writeToR(os,wts::value(ret_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmF_yxmsz ="; wts::writeToR(os,wts::value(tmF_yxmsz),      yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"cpF_fyxmsz="; wts::writeToR(os,wts::value(cpF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmF_fyxmsz="; wts::writeToR(os,wts::value(rmF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmF_fyxmsz="; wts::writeToR(os,wts::value(dmF_fyxmsz),fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmF_fyxmsz="; wts::writeToR(os,            tmF_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;
        os<<"S_list=list("<<endl;
            os<<"sel_vyxmsz="; wts::writeToR(os,wts::value(s_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"Q_vyxms   ="; wts::writeToR(os,wts::value(q_vyxms), vDms,ypDms,xDms,mDms,sDms);       os<<cc<<endl;
            os<<"Q_vyxmsz  ="; wts::writeToR(os,wts::value(q_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")";
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelProcesses(...)"<<endl;

//-------------------------------------------------------------------------------------
//Write model results to file as R list
FUNCTION void ReportToR_ModelResults(ostream& os, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_ModelResults(...)"<<endl;
    
    d3_array ones(1,nSXs,1,nMSs,1,nZBs);//weighting for simple sums
    for (int x=1;x<=nSXs;x++) ones(x) = 1.0;
    
    //population numbers, biomass
    d5_array vn_yxmsz = wts::value(n_yxmsz);
    d5_array vb_yxmsz = tcsam::calcBiomass(vn_yxmsz,ptrMDS->ptrBio->wAtZ_xmz);
//    d4_array n_yxms   = tcsam::calcYXMSfromYXMSZ(vn_yxmsz,ones);
//    d4_array b_yxms   = tcsam::calcYXMSfromYXMSZ(vn_yxmsz,ptrMDS->ptrBio->wAtZ_xmz);
        
    //numbers, biomass captured (NOT mortality)
    d6_array vcpN_fyxmsz = wts::value(cpN_fyxmsz);
    d6_array vcpB_fyxmsz = tcsam::calcBiomass(vcpN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
//    d5_array cpN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vcpN_fyxmsz,ones);
//    d5_array cpB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vcpN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass discards (NOT mortality)
    d6_array vdsN_fyxmsz = wts::value(dsN_fyxmsz);
    d6_array vdsB_fyxmsz = tcsam::calcBiomass(vdsN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
//    d5_array dsN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdsN_fyxmsz,ones);
//    d5_array dsB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdsN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass retained (mortality)
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    d6_array vrmB_fyxmsz = tcsam::calcBiomass(vrmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
//    d5_array rmN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vrmN_fyxmsz,ones);
//    d5_array rmB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vrmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass discard mortality
    d6_array vdmN_fyxmsz = wts::value(dmN_fyxmsz);
    d6_array vdmB_fyxmsz = tcsam::calcBiomass(vdmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
//    d5_array dmN_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdmN_fyxmsz,ones);
//    d5_array dmB_fyxms   = tcsam::calcIYXMSfromIYXMSZ(vdmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //survey numbers, biomass
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vb_vyxmsz = tcsam::calcBiomass(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
//    d5_array n_vyxms   = tcsam::calcIYXMSfromIYXMSZ(vn_vyxmsz,ones);
//    d5_array b_vyxms   = tcsam::calcIYXMSfromIYXMSZ(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    os<<"mr=list("<<endl;
        os<<"iN_xmsz ="; wts::writeToR(os,vn_yxmsz(mnYr),xDms,mDms,sDms,zbDms); os<<cc<<endl;
        os<<"P_list=list("<<endl;
            os<<"MB_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms);                       os<<cc<<endl;
//            os<<"B_yxms   ="; wts::writeToR(os,       b_yxms,ypDms,xDms,mDms,sDms);             os<<cc<<endl;
            os<<"N_yxmsz  ="; wts::writeToR(os,     vn_yxmsz,ypDms,xDms,mDms,sDms,zbDms);        os<<cc<<endl;
            os<<"B_yxmsz  ="; wts::writeToR(os,     vb_yxmsz,ypDms,xDms,mDms,sDms,zbDms);        os<<cc<<endl;
            os<<"nmN_yxmsz="; wts::writeToR(os,wts::value(nmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmN_yxmsz="; wts::writeToR(os,wts::value(tmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;    
        os<<"F_list=list("<<endl;
//            os<<"cpB_fyxms ="; wts::writeToR(os,cpB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
//            os<<"dsB_fyxms ="; wts::writeToR(os,dsB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
//            os<<"rmB_fyxms ="; wts::writeToR(os,rmB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
//            os<<"dmB_fyxms ="; wts::writeToR(os,dmB_fyxms,fDms,yDms,xDms,mDms,sDms); os<<cc<<endl;
            os<<"cpN_fyxmsz="; wts::writeToR(os,vcpN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dsN_fyxmsz="; wts::writeToR(os,vdsN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmN_fyxmsz="; wts::writeToR(os,vrmN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmN_fyxmsz="; wts::writeToR(os,vdmN_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"cpB_fyxmsz="; wts::writeToR(os,vcpB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dsB_fyxmsz="; wts::writeToR(os,vdsB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"rmB_fyxmsz="; wts::writeToR(os,vrmB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"dmB_fyxmsz="; wts::writeToR(os,vdmB_fyxmsz,fDms,yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;
        os<<"S_list=list("<<endl;
           os<<"MB_vyx  ="; wts::writeToR(os,value(mb_vyx),vDms,ypDms,xDms);            os<<cc<<endl;
//           os<<"B_vyxms ="; wts::writeToR(os,      b_vyxms,vDms,ypDms,xDms,mDms,sDms);  os<<cc<<endl;
           os<<"N_vyxmsz="; wts::writeToR(os,    vn_vyxmsz,vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
           os<<"B_vyxmsz="; wts::writeToR(os,    vb_vyxmsz,vDms,ypDms,xDms,mDms,sDms,zbDms); os<<endl;
       os<<")";
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelResults(...)"<<endl;
    
//-------------------------------------------------------------------------------------
//Write quantities related to model fits to file as R list
FUNCTION void ReportToR_ModelFits(ostream& os, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_ModelFits(...)"<<endl;
    //recalc objective function components and and write results to os
    os<<"model.fits=list("<<endl;
        os<<tb<<"penalties="; calcPenalties(-1,os);      os<<cc<<endl;
        os<<tb<<"priors=";    calcAllPriors(-1,os);      os<<cc<<endl;
        os<<tb<<"components=list("<<endl;
            os<<tb<<tb<<"recruitment="; calcNLLs_Recruitment(-1,os); os<<endl;
        os<<tb<<")"<<cc<<endl;
        os<<tb<<"fisheries=";  calcNLLs_Fisheries(-1,os);  os<<cc<<endl; 
        os<<tb<<"surveys=";    calcNLLs_Surveys(-1,os);    os<<cc<<endl;  
        os<<tb<<"growthdata="; calcNLLs_GrowthData(-1,os); os<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelFits(...)"<<endl;

//-------------------------------------------------------------------------------------
//Update MPI for current parameter values (mainly for export))
FUNCTION void updateMPI(int debug, ostream& cout)
    if (debug) cout<<"Starting updateMPI(...)"<<endl;
//    NumberVectorInfo::debug=1;
//    BoundedVectorInfo::debug=1;
//    DevsVectorInfo::debug=1;
//    DevsVectorVectorInfo::debug=1;
    //recruitment parameters
    ptrMPI->ptrRec->pLnR->setFinalVals(pLnR);
    ptrMPI->ptrRec->pLnRCV->setFinalVals(pLnRCV);
    ptrMPI->ptrRec->pLgtRX->setFinalVals(pLgtRX);
    ptrMPI->ptrRec->pLnRa->setFinalVals(pLnRa);
    ptrMPI->ptrRec->pLnRb->setFinalVals(pLnRb);
    //cout<<"setting final vals for pDevsLnR"<<endl;
    for (int p=1;p<=ptrMPI->ptrRec->pDevsLnR->getSize();p++) (*ptrMPI->ptrRec->pDevsLnR)[p]->setFinalVals(pDevsLnR(p));
     
    //natural mortality parameters
    ptrMPI->ptrNM->pLnM->setFinalVals(pLnM);
    ptrMPI->ptrNM->pLnDMT->setFinalVals(pLnDMT);
    ptrMPI->ptrNM->pLnDMX->setFinalVals(pLnDMX);
    ptrMPI->ptrNM->pLnDMM->setFinalVals(pLnDMM);
    ptrMPI->ptrNM->pLnDMXM->setFinalVals(pLnDMXM);
    
    //growth parameters
    ptrMPI->ptrGrw->pLnGrA->setFinalVals(pLnGrA);
    ptrMPI->ptrGrw->pLnGrB->setFinalVals(pLnGrB);
    ptrMPI->ptrGrw->pLnGrBeta->setFinalVals(pLnGrBeta);
    
    //maturity parameters
    //cout<<"setting final vals for pLgtPrM2M"<<endl;
    for (int p=1;p<=npLgtPrMat;p++) (*ptrMPI->ptrM2M->pLgtPrM2M)[p]->setFinalVals(pLgtPrM2M(p));
    
    //selectivity parameters
    ptrMPI->ptrSel->pS1->setFinalVals(pS1);
    ptrMPI->ptrSel->pS2->setFinalVals(pS2);
    ptrMPI->ptrSel->pS3->setFinalVals(pS3);
    ptrMPI->ptrSel->pS4->setFinalVals(pS4);
    ptrMPI->ptrSel->pS5->setFinalVals(pS5);
    ptrMPI->ptrSel->pS6->setFinalVals(pS6);
    //cout<<"setting final vals for pDevsS1"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS1->getSize();p++) (*ptrMPI->ptrSel->pDevsS1)[p]->setFinalVals(pDevsS1(p));
    //cout<<"setting final vals for pDevsS2"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS2->getSize();p++) (*ptrMPI->ptrSel->pDevsS2)[p]->setFinalVals(pDevsS2(p));
    //cout<<"setting final vals for pDevsS3"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS3->getSize();p++) (*ptrMPI->ptrSel->pDevsS3)[p]->setFinalVals(pDevsS3(p));
    //cout<<"setting final vals for pDevsS4"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS4->getSize();p++) (*ptrMPI->ptrSel->pDevsS4)[p]->setFinalVals(pDevsS4(p));
    //cout<<"setting final vals for pDevsS5"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS5->getSize();p++) (*ptrMPI->ptrSel->pDevsS5)[p]->setFinalVals(pDevsS5(p));
    //cout<<"setting final vals for pDevsS6"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS6->getSize();p++) (*ptrMPI->ptrSel->pDevsS6)[p]->setFinalVals(pDevsS6(p));
     
    //fully-selected fishing capture rate parameters
    ptrMPI->ptrFsh->pHM->setFinalVals(pHM);
    ptrMPI->ptrFsh->pLnC->setFinalVals(pLnC);
    ptrMPI->ptrFsh->pLnDCT->setFinalVals(pLnDCT);
    ptrMPI->ptrFsh->pLnDCX->setFinalVals(pLnDCX);
    ptrMPI->ptrFsh->pLnDCM->setFinalVals(pLnDCM);
    ptrMPI->ptrFsh->pLnDCXM->setFinalVals(pLnDCXM);
    //cout<<"setting final vals for pDevsLnC"<<endl;
    for (int p=1;p<=ptrMPI->ptrFsh->pDevsLnC->getSize();p++) (*ptrMPI->ptrFsh->pDevsLnC)[p]->setFinalVals(pDevsLnC(p));
    
    //survey catchability parameters
    ptrMPI->ptrSrv->pLnQ->setFinalVals(pLnQ);
    ptrMPI->ptrSrv->pLnDQT->setFinalVals(pLnDQT);
    ptrMPI->ptrSrv->pLnDQX->setFinalVals(pLnDQX);
    ptrMPI->ptrSrv->pLnDQM->setFinalVals(pLnDQM);
    ptrMPI->ptrSrv->pLnDQXM->setFinalVals(pLnDQXM);
    
    if (debug) cout<<"Finished updateMPI(...)"<<endl;

//-------------------------------------------------------------------------------------
//Write results to file as R list
FUNCTION void ReportToR(ostream& os, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR(...)"<<endl;

    updateMPI(debug,cout);
        
    os<<"res=list("<<endl;
        //model configuration
        ptrMC->writeToR(os,"mc",0); os<<","<<endl;
        
        //model data
        ptrMDS->writeToR(os,"data",0); os<<","<<endl;
        
        //parameter values
        ReportToR_Params(os,debug,cout); os<<","<<endl;
        
        //model processes
        ReportToR_ModelProcesses(os,debug,cout); os<<","<<endl;
        
        //model results
        ReportToR_ModelResults(os,debug,cout); os<<","<<endl;

        //model fit quantities
        ReportToR_ModelFits(os,debug,cout); os<<","<<endl;
        
        //simulated model data
        createSimData(debug, cout, 0, ptrSimMDS);//deterministic
        ptrSimMDS->writeToR(os,"sim.data",0); 
        
        //do OFL calculations
        if (doOFL){
            cout<<"ReportToR: starting OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
            calcOFL(mxYr+1,1,echoOFL);//updates oflResults
            oflResults.writeCSVHeader(echoOFL); echoOFL<<endl;
            oflResults.writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL.close();
            os<<","<<endl;
            oflResults.writeToR(os,"oflResults",0);
            cout<<"ReportToR: finished OFL calculations"<<endl;
        }

    os<<")"<<endl;
    if (debug) cout<<"Finished ReportToR(...)"<<endl;

//----------------------------------------------------------------------
//Write parameter information to file
FUNCTION void writeParameters(ostream& os,int toR, int willBeActive)           
    os<<"index, phase, idx.mn, idx.mx, min, max, value, name, type"<<endl;
    //recruitment parameters
    wts::writeParameter(os,pLnR,toR,willBeActive);      
    wts::writeParameter(os,pLnRCV,toR,willBeActive);      
    wts::writeParameter(os,pLgtRX,toR,willBeActive);      
    wts::writeParameter(os,pLnRa,toR,willBeActive);      
    wts::writeParameter(os,pLnRb,toR,willBeActive);      
    wts::writeParameter(os,pDevsLnR,toR,willBeActive);      
    
    //natural mortality parameters
    wts::writeParameter(os,pLnM,toR,willBeActive);      
    wts::writeParameter(os,pLnDMT,toR,willBeActive);      
    wts::writeParameter(os,pLnDMX,toR,willBeActive);      
    wts::writeParameter(os,pLnDMM,toR,willBeActive);      
    wts::writeParameter(os,pLnDMXM,toR,willBeActive);      
    
    //growth parameters
    wts::writeParameter(os,pLnGrA,toR,willBeActive);      
    wts::writeParameter(os,pLnGrB,toR,willBeActive);      
    wts::writeParameter(os,pLnGrBeta,toR,willBeActive);      
    
    //maturity parameters
    wts::writeParameter(os,pLgtPrM2M,toR,willBeActive);      
    
    //selectivity parameters
    wts::writeParameter(os,pS1,toR,willBeActive);      
    wts::writeParameter(os,pS2,toR,willBeActive);      
    wts::writeParameter(os,pS3,toR,willBeActive);      
    wts::writeParameter(os,pS4,toR,willBeActive);      
    wts::writeParameter(os,pS5,toR,willBeActive);      
    wts::writeParameter(os,pS6,toR,willBeActive);      
    wts::writeParameter(os,pDevsS1,toR,willBeActive);      
    wts::writeParameter(os,pDevsS2,toR,willBeActive);      
    wts::writeParameter(os,pDevsS3,toR,willBeActive);      
    wts::writeParameter(os,pDevsS4,toR,willBeActive);      
    wts::writeParameter(os,pDevsS5,toR,willBeActive);      
    wts::writeParameter(os,pDevsS6,toR,willBeActive);      
    
    //fishery parameters
    wts::writeParameter(os,pHM,toR,willBeActive);      
    wts::writeParameter(os,pLnC,toR,willBeActive);      
    wts::writeParameter(os,pLnDCT,toR,willBeActive);      
    wts::writeParameter(os,pLnDCX,toR,willBeActive);      
    wts::writeParameter(os,pLnDCM,toR,willBeActive);      
    wts::writeParameter(os,pLnDCXM,toR,willBeActive);      
    wts::writeParameter(os,pDevsLnC,toR,willBeActive);      
    
    //survey parameters
    wts::writeParameter(os,pLnQ,toR,willBeActive);      
    wts::writeParameter(os,pLnDQT,toR,willBeActive);      
    wts::writeParameter(os,pLnDQX,toR,willBeActive);      
    wts::writeParameter(os,pLnDQM,toR,willBeActive);      
    wts::writeParameter(os,pLnDQXM,toR,willBeActive);      
    
// =============================================================================
// =============================================================================
REPORT_SECTION
        
    //write active parameters to rpt::echo
    rpt::echo<<"Finished phase "<<current_phase()<<endl;
    {
        //write objective function components only
        ofstream os0("tcsam02.ModelFits."+itoa(current_phase(),10)+".R", ios::trunc);
        ReportToR_ModelFits(os0,0,cout);
        os0.close();
    }
    if (last_phase()) {
        //write report as R file
        ReportToR(report,1,rpt::echo);
        //write parameter values to csv
        ofstream os1("tcsam02.params.all.final.csv", ios::trunc);
        writeParameters(os1,0,0);
        os1.close();
        //write parameter values to csv
        ofstream os2("tcsam02.params.active.final.csv", ios::trunc);
        writeParameters(os2,0,1);
        os2.close();

        if (option_match(ad_comm::argc,ad_comm::argv,"-jitter")>-1) {
            ofstream fs("jitterInfo.csv");
            fs<<"seed"<<cc<<"objfun"<<endl;
            fs<<iSeed<<cc<<objFun<<endl;
        }
    }
    
// =============================================================================
// =============================================================================
BETWEEN_PHASES_SECTION
    rpt::echo<<endl<<endl<<"#---------------------"<<endl;
    rpt::echo<<"Starting phase "<<current_phase()<<" of "<<initial_params::max_number_phases<<endl;
    ctrProcCallsInPhase=0;//reset in-phase counter

// =============================================================================
// =============================================================================
FINAL_SECTION
    {cout<<"writing model sim data to file"<<endl;
        ofstream echo1; echo1.open("ModelSimData.dat", ios::trunc);
        writeSimData(echo1,0,cout,ptrSimMDS);
    }

    if (option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1) {
        mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::app);
        mcmc<<"NULL)"<<endl;
        mcmc.close();
    }
    
    long hour,minute,second;
    double elapsed_time;
    
    time(&finish); 
    elapsed_time = difftime(finish,start);
    
    hour = long(elapsed_time)/3600;
    minute = long(elapsed_time)%3600/60;
    second = (long(elapsed_time)%3600)%60;
    cout << endl << endl << "Starting time: " << ctime(&start);
    cout << "Finishing time: " << ctime(&finish);
    cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;
    
// =============================================================================
// =============================================================================
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 5000,5000,5000,5000,5000,5000,10000
//  convergence_criteria 0.1,0.1,.01,.001,.001,.001,1e-3,1e-4
  convergence_criteria 0.5,0.1,.01,.001,1e-4,1e-5,1e-6

// =============================================================================
// =============================================================================
TOP_OF_MAIN_SECTION
  arrmblsize = 1000000000; //must be smaller than 2,147,483,647
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(40000000); // this may be incorrect in the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1500000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(4000);
  time(&start);

