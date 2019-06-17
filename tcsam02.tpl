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
//  2016-03-31: 1. Revised NM calc.s so pM represents natural mortality rates
//                  on immature male crab (was mature male crab, previously). pDM3
//                  and pDM4 now reprsent corresponding offsets for mature, not
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
//--2016-12-13: 1. optGrowth=0 now calculates growth matrices identically to TCSAM2013
//                  by shifting calculation of density to lower size bin cutpoint.
//              2. Changed width of growth matrices from 10 bins to 11 bins 
//                  to match TCSAM2013.
//--2016-12-16: 1. Modified memory allocation for output biomass arrays in
//                  calcBiomass functions in SummaryFunctions.cpp. Previous
//                  declarations used default constructors, which only copied
//                  the reference to the input abundance array (i.e., used a 
//                  shallow copy); thus, the output biomass array referenced the 
//                  same memory addresses as the input abundance array.
//--2017-01-03: 1. Added call to ReportToR_ModelFits() in PRELIM_CALCS section
//                  to get initial likelihood values.
//              2. Added "yr" to likelihood info for size comps.
//              3. Added options for penalizing smoothness and non-decreasing
//--2017-01-04: 1. revised calc1stDiffs & calc2ndDiffs to use ADMB first_difference().
//              2. revised prM2M smoothness penalties to use 0.5*norm2(...).
//              3. calcNLL functions now change value of smlVal to reflect
//                  parallel values in TCSAM2013.
//--2017-01-08: 1. Added jitter seed, if used, to model report file as "jitter".
//              2. Added objFun to ReportToR_ModelFits as "objfun".
//--2017-01-09: 1. Added checkParams functions.
//              2. added ctrDebugParams commandline option.
//--2017-02-01: 1. Increased precision on all output to 12.
//--2017-02-06: 1. changed mfexp to exp in calcGrowth for prsp to match TCSAM2013.
//              2. Increment tcsam::VERSION to "2017.02.06".
//--2017-02-13: 1. Revised calcGrowth for TCSAM2013 growth to match more closely.
//                  Apparently fixed nan issue when estimating TCSAM2013-type growth.
//              2. Incremented tcsam::VERSION to "2017.02.13".
//              3. Changed penalty weights for prM2M so that they are vectors,
//                  with weights assigned for each parameter combination.
//--2017-02-20: 1. Adjusted setInitVals for jittered devs vectors to ensure all values
//                  were within bounds.
//              2. Set jittering fraction for devs to 0.1*jitFrac
//              3. Added maxGrad to ReportToR() and ReportToR_ModelFits()
//--2017-02-22: 1. Added equilibrium unfished size distribution to OFL output.
//              2. Incremented model version to "2012.02.22".
//--2017-02-23: 1. Added equilibrium fished size distribution to OFL output.
//              2. Added option to parameterize natural mortality on arithmetic
//                  scale as in TCSAM2013 in ModelOptions. 
//              3. Incremented model options version to "2012.02.23".
//              4. Incremented model version to "2012.02.23".
//--2017-02-27: 1. Added model options to calculate OFL using averaged selectivities
//                  and max capture rates, rather than averaged capture rates
//              2. Added options to use external values for max capture rates
//                  in OFL calculations (ala SCF in TCSAM2013)
//              3. Added options for extrapolating effort to capture rates
//              4. Incremented ModelOptions version to "2012.02.27".
//--2017-03-27: 1. Fixed implementation of new Model Options.
//              2. Working on implementation of effort extrapolation
//              3. Incremented model version to 2017-03-27.
//--2017-04-03: 1. Implemented effort extrapolation inputs for ModelOptions. 
//              2. Implemented effort extrapolation logic in tcsam02.tpl.
//              3. Revised FisheriesInfo parameter combinations format to specify 
//                  effort extrapolation scenario.
//              4. Added "version" to ModelParametersInfo file as consistency check.
//              5. Working on output for effort extrapolation.
//--2017-04-05: 1. Revised debug output for OFL calculations.
//--2017-04-07: 1. Revised effort extrapolation to use observed, predicted effort,
//                  not "observed", "predicted" capture rates in the likelihood.
//--2017-04-13: 1. Fixed error in WriteToR for an IndexBlock.
//              2. Added output from calcNLLs_ExtrapolatedEffort to rep file 
//                  via ReportToR_ModelFits.
//--2017-04-14: 1. Added zscrEffX_nxmsy and nllEffX_nxms arrays, updated 
//                  R output in calcNLLs_ExtrapolatedEffort().
//--2017-05-10: 1. Added labels to I/O for GrowthData datasets
//              2. Added labels to parameter combinations and parameters via
//                  ModelParametersInfo file.
//--2017-05-15: 1. Added pLgtRet as parameter for retained catch fraction (for old shell crab)
//                  in directed fisheries.
//              2. Refactored survey LnDQ parameter names from 
//                  pLnDQT,pLnDQX,pLnDQM,pLnDQXM
//                 to
//                  pDQ1,pDQ2,pDQ3,pDQ4
//                 to reflect more generic usage
//              3. Refactored natural mortality pLnDM parameter names from 
//                  pDMT,pDMX,pDMM,pDMXM
//                 to
//                  pDM1,pDM2,pDM3,pDM4
//                 to reflect more generic usage
//              4. Added logic to apply pLgtRet as retFrac.
//--2017-06-05: 1. Refactored survey pLnDC parameter names from 
//                  pDC1,pDC2,pDC3,pDC4 and associated idx's
//                 to
//                  pDC1,pDC2,pDC3,pDC4
//                 to reflect more generic usage
//              2. Refactored calcFisheryFs: replaced lnC_m, arC_m with lnC, arC
//                  so function uses input sex/maturity/shell condition categories
//--2017-06-06: 1. revised wts::writeParameter functions to tcsam::writeParameters
//                  to include two categories and parameter labels in output
//              2. Incremented model version to 2017.06.06.
//--2017-06-12: 1. revised calcOFL to avoid division by 0 when fisheries
//                  are not conducted within the averaging periods for
//                  handling mortality, fishery capture rates, selectivity and retention curves
//--2017-06-19: 1. added commandline variable "phsItsRewgt" to set min phase to calculate 
//                  effective weights for size compositions
//              2. added Mc-I harmonic mean and Francis weights calculations based
//                  on Punt, 2017
//--2017-07-08: 1. added cout's in parameter section to help with deviant pin files 
//--2017-07-17: 1. changed single quotes to double quotes in csv-type output in
//                  all tcsam::writeParameter(...) functions so Excel does the correct thing.
//--2017-07-19: 1. Added calcCohortProgression(...) and call to it in report section.
//--2017-08-01: 1. Added final year MMB, OFL info to jitterInfo.csv output
//              2. calcOFL ONLY called in last phase now, unless debugOFL has been set
//--2017-08-02: 1. Added likeprof_numbers lkMMB, lkQM, lkQF, lkNMMM, lkNMMF, lkNMIM, lkNMIF
//--2017-08-12: 1. Added selectivity function "asclogistic5099", with parameters z50 and z99.
//--2017-08-23: 1. Corrected problem with cohort progression when directed fishery was closed 
//                  (i.e., when maxF=0).
//              2. Added selectivity function "asclogistic95Ln50" with parameters z95 and ln(z95-z50)
//--2017-08-24: 1. Added objective function penalties to keep growth parameters in range
//              2. Added option for re-parameterized mean growth functions in terms of min,max post-molt sizes
//              3. Dropped "Ln" notation from all growth parameters
//              4. Incremented ModelParametersInfo version to "2017.08.24"
//              5. Incremented ModelOptions version to "2017.08.24"
//              6. refactored getPCXDs(pc) from inheriting classes to base class ParameterGroupInfo
//              7. updated calcNLLs_GrowthData to reflect to growth options
//              8. Outputting name assigned to growth dataset as name of R list in model fits
//--2017-08-25: 1. Corrected z-score calculation for growth data to use grBeta_xy, not grB_xy,
//                   in variance calculation
//--2017-08-30: 1. Changed fishery averaging period for OFL calculations from 1 year to 5.
//                  TODO: this should be an input in the Model Options file.
//--2017-09-18: 1. Changes to OFL calculations to try to output quantities for old projection model.
//--2017-10-03: 1. Modifying code to incorporate iterative re-weighting of size comps
//                  by McAllister-Ianelli or Francis methods
//--2017-10-15: 1. Finished modifying code to incorporate iterative re-weighting.
//--2017-10-20: 1. Added ascending normal and 4- and 6-parameter versions of the 
//                  double normal selectivity functions
//--2017-10-23: 1. Implemented new version of DevsVectors such that a vector has
//                  all elements bounded as a BoundedVector but the sum-to-1 constraint
//                  must be enforced in the likelihood. Previously, the sum-to-1
//                  constraint was identically satisfied by setting the final element
//                  to the negative sum of the previous elements, while the bounds on
//                  the final element were enforced in the likelihood.
//--2017-10-31: 1. Revised ParameterInfoTypes to incorporate inputs to parameters 
//                  on non-arithmetic scales. Changed growth parameters to implement this.
//              3. Stopped printing debugging info for new stuff by default.
//--2017-10-31: 1. Implemented non-arithmetic input scales for all devs and for 
//                  all selectivity parameters.
//              2. Implemented non-arithmetic input scales for recruitment parameters
//                  except pLnR, and for pLnM (now pM) and pLnQ (now pQ). 
//              3. Added rmse's to output for likelihood calculations.
//--2017-11-08: 1. Refactored ModelParameterInfoTypes to better handle parameter scaling.
//                  BounderVectorInfo now inherits from VectorInfo and method names 
//                  were changed or added to reflect whether inputs/outputs are
//                  on the parameter scale. 
//--2017-11-13: 1. Continued refactoring ModelParameterInfoTypes to better handle parameter scaling.
//              2. Changed MPI version to 20171113.
//--2017-11-14: 1. More debugging info being printed, documentation improved.
//              2. Added scripts under "docs/scripts" to produce html documentation
//                  using doxygen.
//--2017-11-15: 1. Revised csv output from writeParameters(...) to include min/max/value
//                  on both arithmetic and parameter scales.
//              2. TODO: change jittering to occur on ARITHMETIC scales, not parameter
//                  scales (logit-scale parameters tend to end up very near bounds when
//                  jittering occurs on logit-scale).
//--2017-11-12: 1. Finished 2. above. Fixed problems with jittering devs.
//              2. Updated VERSION to "2017.11.16".
//--2017-11-19: 1. Fixed problem with BoundedVectorInfo::writeToR(...)
//--2017-11-25: 1. Fixed problem with BoundedVectorVectorInfo::write(...)
//              2. Fixed problem with BoundedVectorInfo::calcParamScaleVals(...)
//              3. Fixed problem with calcNatMort(...) not converting params to
//                  arithmetic scale.
//--2017-11-27: 1. Increased penalty on growth increments approaching 0 and
//                  added debugging output to rpt::echo when penalty is invoked.
//--2017-11-29: 1. Revised ModelOptions so file specifies the likelihood weight and "eps"
//                  value in posfun(...) for the penalty on approaching negative 
//                  growth increments.
//              2. Revised ModelOptions to reflect change in devs vectors info 
//                  implementation from old approach (last dev = -sum(other devs)
//                  but penalized if not within bounds)
//                  to new approach (square(sum(devs)) penalized for > 0)
//              3. ModelOptions version updated to 2017.11.29.
//--2017-12-01: 1. Revised output of "GrowthData.NLLs.NanReport.dat"
//              2. ModelOptions revised to specify minimum/max pre-molt CWs at which
//                  the growth increment must be positive. Needed to do this because
//                  it is possible to have growth data from crab sizes outside the
//                  model range that should be included in growth estimation. These
//                  are enforced by applying posfun-type likelihood penalties when 
//                  growth increments start to approach negative values.
//              3. Updated ModelOptions version to "2017.12.01".
//              4. Changed 'final.devs' to 'devsSumSq' in R list for likelihood penalties.
//--2017-12-05: 1. Revised ALL selectivity functions to use mfexp() rather than exp().
//                  Using exp() in the ...normal() functions led to model instability!
//              2. Added save_params() function, which will save parameter values
//                  to a file "tcsam02.PP.XX.par" when called, with PP replaced by the 
//                  current estimation phase and XX by the number of PROCEDURE_SECTION
//                  function calls within the phase. Should only be used for debugging
//                  convergence problems.
//              3. A csv version of 2 is available (but commented out) at th end of the
//                  PROCEDURE_SECTION code section.
//              4. Updated tcsam::VERSION to "2017.12.05".
//--2017-12-06: 1. Implemented first cut at incorporating chela height data (as male
//                  maturity ogives) into likelihood calculations.
//              2. Updated tcsam::VERSION to "2017.12.06".
//--2018-01-30: 1. Implemented a Dirichlet NLL for size comps. TODO: Need to figure out
//                  how to implement associated effective N's as parameters.
//--2018-02-12: 1. TODO: Need to figure out how to skip growth data from years > assmt year.
//--2018-02-15: 1. Added tb to end of VectorInfo::writePart1() so output could be read back in correctly.
//--2018-02-21: 1. writing gradients to csv file in report section.
//              2. set nllWgtRecDevs = 1 (was 0) in calcNLLs_Recruitment()
//--2018-02-22: 1. Revised rmse calculation in calcNLLs_Recruitment() to use actual number
//                  of z-scores, not size of z-scores vector (because vector runs mnYr:mxYr).
//              2. added nllWgt by parameter combination to recruitment parameters info.
//              3. updated MPI version to "2018.02.22".
//              4. revised calculations in calcNLLs_Recruitment.
//--2018-02-26: 1. Revised calcNLLs_ChelaHeightData to turn "on" nll calculation in objective function
//                  only when estimation is active.
//              2. Added output to report file reflecting predicted male maturity ogives for
//                  population and surveys
//-2018-03-05:  1. Revised calcNLLs_Recruitment because zscores were being recalculated within
//                  the function, which might be called multiple times at the REPORT stage
//                  leading to changes in the zscores, associated NLLs, and objective function 
//                  even though no parameter changes had occurred.
//-2018-03-13:  1. Revised approach to jittering because old version resulted in too high
//                  a proportion of values near a limit if the default value was near that limit
//-2018-03-15:  1. Eliminated writes to terminal screen in PRELIMINARY_CALCS
//-2018-03-19:  1. Now sets usePin to 2 if "-mcpin" is a commandline argument (i.e., running NUTS MCMC)
//-2018-03-21:  1. Adding commandline option "calcDynB0" to run dynamic B0 calculations after final phase
//              2. Updated model version to 2018.03.21.
//-2018-03-28:  1. Revised calcPrM2M to correctly handle small size bins being fixed at 0.
//-2018-04-05:  1. Added new growth parameterization option based on post-molt size at specified 
//                  pre-molt size and ln-scale slope.
//              2. Incremented ModelOptions.VERSION to 2018.04.05.
//              3. Added print statements in PROCEDURE_SECTION to diagnose nuts behavior.
//              4. Changed modVer to be specified in tpl (easier to update).
//              5. modVer incremented to 2018.04.05.
//              6. Added macros PRINT2B1 and PRINT2B2 to facilitate printing to
//                  both std::cout and rpt::echo. Converted lots of writes.
//-2018-04-08:  1. Modified calcNLLs_ChelaHeightData so it always adds wgt*nll to
//                  the objective function (previously, this only occurred when
//                  the prM2M parameters were active--thus never in the
//                  PRELIMINARY_CALCS section).
//-2018-04-10:  1. Fix to output to R in calcNoneNLLs for size comps.
//              2. changed command line option from -dynB0 to -calcDynB0
//              3. changed lists that use survey or fishery names to surround
//                  them with `s to be valid names in R
//-2018-04-11:  1. Revised calculations in calcNLLs_GrowthData and calcNLLs_MaturityData
//                  so that likelihoods are calculated for datasets with likelihood weights
//                  set to 0 (so the fit is not included in the parameter optimization)
//                  ONLY when output to R is requested (i.e., when debug<0 in call to function).
//              2. changed lists that use growth dataset or maturity dataset names to surround
//                  them with `s to be valid names in R
//              3. Added 'GAMMA' (LL_GAMMA and STR_LL_GAMMA) as likelihood type
//              4. Incremented tcsam::VERSION and modVer to 2018.04.11.
//-2018-04-24:  1. Revised some print statements (need to surround PRINT2B macros
//                  with brackets to include in a 1-line if statement.
//-2018-04-29:  1. Revised output for iterative reweighting to write to "effectiveWeights.R".
//-2018-08-24:  1. Incremented version to "2018.08.24".
//              2. Added FIT_BY_X_MSE option to fit size compositions
//-2018-08-27:  1. Corrected error in normalization factor for asclogistic5095 function.
//              2. Incremented version to "2018.08.27".
//              3. Fixed effN=nan problem with output to R when size comp component 
//                  is identically 0 but ss>0 (possible when fitting an 'extended' comp).
//-2018-08-29:  1. Incremented version to "2018.08.29".
//              2. Corrected missing calculation of penalty on approaching negative growth
//                  increments when growth parameterization option 2 was selected 
//                  (did not affect optimization because it simply added a very
//                  large constant to the objective function).
//-2018-08-30:  1. Corrected (again) missing calculation of penalty on approaching negative growth
//                  increments when growth parameterization option 2 was selected 
//                  (did not affect optimization because it simply added a very
//                  large constant to the objective function).
//-2018-09-02:  1. Revised format for mcmc output to R.
//              2. Revised threshold on discrepancy between msy and totCM that
//                  triggers diagnostic printout in OFL_Calculator::calcMSY(...).
//              3. Model version now taken from tcsam::VERSION in ModelConstants.hpp.
//-2018-09-21:  1. Revising code to accommodate running in an operating mode and an 
//                  estimation mode to facilitate MSEs which will require feedback
//                  between the two modes and extension to subsequent years.
//-2018-10-29:  1. Added pvNPSel parameter vector implementation to be able to estimate
//                  BSFRF availability. Modified tpl and ModelParametersInfo.
//              2. Added pMSE_LnC
//-2018-11-06:  1. Added "runAlt" commandline flag to run the alternative pop dy model
//              2. Implemented runAltPopDyMod based on PopDyInfo, etc. classes converted
//                  to dvar-type variables.
//-2018-11-07:  1. Added HarvestStrategies for the MSE; revised ModelOptions
//                  to incorporate MSE-related options.
//-2018-11-13:  1.Changes to allow pMSE_LnC to be estimated correctly. Now need to
//                  figure out how to write output files to step between MSE OpMod 
//                  and EstMod iterations.
//-2018-11-18: 1. Lots of changes to tcsam02, ModelParametersInfo, ModelParameterInfoTypes,
//                  ModelData, CatchData, and SummaryFunctions to incorporate
//                  necessary MSE OpMod gyrations (not completed yet, but getting there).
//-2018-11-20: 1. Revised OpModMode so it reads in "state" from a file and only
//                  calculates the projected year. Also created MSE_OpModInfo class.
//-2018-11-28: 1. Revised a number of classes to be able to update an input MPI file
//                  for "next" year and write it and an associated pin file for use
//                  in an MSE.
//-2018-11-30: 1. Added check in the PRELIMIINARY_CALCS for input TAC=0 (directed fishery closed) 
//                  when running the operating model (mseOpModMode==1) so that pMSE_LnC is not estimated.
//             2. Added function projectPopForZeroTAC to project population when TAC=0
//                  (directed fishery is closed).
//             3. Moved opModMode calculations in FINAL_SECTION to "finishOpModMode" function
//                  so that these could be called from PRELIMINARY_CALCS when TAC=0. 
//             4. Modified opModMode output of data files to handle closed
//                  directed fishery correctly.
//-2018-12-03: 1. Corrected problems with OFL calculations associated with code changes to
//                  incorporate MSE calculations. OFL results now agree (again) with results
//                  from the 2018 assessment.
//-2018-12-26:  1. Corrected some problems with writing the report file (missing commas, etc).
//              2. Added ReportToR_OpModMode().
//-2018-12-31:  1. Corrected some problems with ReportToR_OpModMode(). Needs to be fleshed out.
//              2. Corrected how R_y is allocated in MSEClasses.cpp from mnYr,mxYr to mnYr,mxYr-1
//                  so that input file is read correctly.
//-2019-01-03:  1. Resolved issues with OpMod files and input. In process, revised a lot of associated
//                  variables and functions. NOTE: don't use "\t" in writing output that will be
//                  read back in to ADMB--ADMB reads it as a bunch of zeros, not a tab character.
//              2. Expanded ReportToR_OpModMode(); Removed sim data from ReportToR().
//              3. Expanded ModelOptions to include min/max years for OpMod recruitment statistics.
//                  ModelOptions version now 2019.01.03.
//-2019-01-28:  1. Expanded MSE_OpModInfo to include time series of recruitment, TACs, and OFLs
//                  and write them to R in the "results" list written by ReportToR_OpModMode()
//                  to facilitate further analysis of model results.
//-2019-01-30:  1. Added boiler plate to calculate inputs for Tanner crab MSE HCR 6.
//-2019-02-13:  1. Starting work on revisions to incorporate SBS data. Need to:
//                  1. Add a survey max availability parameter pA to SurveysInfo
//                  2. Add a constant selectivity function (BSFRF sel assumed =1 for all sizes)
//                  3. Incorporate survey max availability parameter into calcSurveyQs
//              2. Added a constant selectivity function (=1 for all sizes) to ModelSelectivities.
//              3. Added survey max availability parameter pA to SurveysInfo and tpl.
//              4. Incorporated max availability parameter into calcSurveyQs.
//-2019-02-20:  1. Cohort progression now output as separate R files in PRELIMINARY_CALCS
//                  and FINAL sections if run in "normal" mode.
//              2. Added inputs to calcCohortProgression and ReportToR_CohortProgression
//                  to specify number of size bins at recruitment to include and whether
//                  or not to include natural mortality and/or fishing mortality on 
//                  cohort progression.
//              3. Now outputting cohort progressions to file using:
//                  a. 3 size bins, M+F in rep files
//                  b. 1 size bin, no M, no F in "CohortProgression.noMF.R"
//                  c. 1 size bin, M, no F    in "CohortProgression.MnoF.R"
//                  d. 1 size bin, no M, F    in "CohortProgression.FnoM.R"
//                  e. 1 size bin, M+F        in "CohortProgression.MF.R"
//                  f. 1 size bin, M+F        in "CohortProgression.init.R"
//-2019-02-21:  1. Added model option maxGrowthZBEx to be able to change the extent
//                  that immature crab can grow (had been hard-wired at 11).
//              2. Incremented VERSION in ModelOptions to "2019.02.21".
//-2019-04-02:  1. Corrected some problems associated with indexing in the 
//                  nonparametric selectivity functions.
//              2. Corrected some problems with reporting penalties on NPSel functions
//                  to R.
//              3. Fixed "debug" setting in calcCohortProgression call in 
//                  reportToR_CohortProgression.
//-2019-06-17:  1. Added MaturityOgiveData (as opposed to ChelaHeightData) datasets 
//
// =============================================================================
// =============================================================================
//--Commandline Options
//  -configFile  fnConfigFile      default model config filename="tcsam02.ModelConfig.dat"
//  -pin fnPin                     use ascii  pin file fnPin to set parameter values. if fnPin not given, fnPin = "tcsam02.pin"
//  -binp fnPin                    use binary pin file fnPin to set parameter values
//  -ainp fnPin                    use ascii  pin file fnPin to set parameter values
//  -mceval                        flag that mceval is ON (i.e., mcevalOn=1)
//  -mcpin fnPin                   pin file to run NUTS MCMC (i.e., usePin=2)
//  -runAlt                        flag to run alternative pop dy equations
//  -calcOFL                       flag to turn on OFL calculations
//  -calcTAC hcr                   flag to calculate TAC using the indicated harvest control rule (hcr)
//  -calcDynB0                     flag to turn on dynamic B0 calculations
//  -doRetro  yRetro               flag to turn on retrospective calculations, dropping yRetro years
//  -jitter                        flag to turn on parameter jittering
//  -resample                      flag to turn on parameter resampling
//  -iSeed  val                    use val as random number generator (RNG) seed
//  -fitSimData iSimDataSeed       flag to fit simulated data using RNG seed iSimDataSeed
//  -mseMode type                  flag to use mseMode (type = "mseOpModMode" or "mseEstModMode")
///////----flags to print debugging info-----
//  -debugModelConfig
//  -debugModelParams
//  -ctrDebugParams ctr            flag to turn on debugParams when PROCEDURE counter reaches ctr
//  -debugDATA_SECTION
//  -debugPARAMS_SECTION
//  -debugPRELIM_CALCS
//  -debugPROC_SECTION
//  -debugREPORT_SECTION
//  -debugOFL
//  -debugModelDatasets
//  -debugModelParamsInfo
//  -debugModelOptions
//  -debugRunModel
//  -debugObjFun
//  -debugMCMC
//  -showActiveParams
// =============================================================================
// =============================================================================
GLOBALS_SECTION
    #include <math.h>
    #include <time.h>
    #include <limits>
    #include <admodel.h>
    #include "TCSAM.hpp"

    #ifdef PRINT2B1
        #undef PRINT2B1
    #endif
    #define PRINT2B1(o) std::cout<<(o)<<std::endl; rpt::echo<<(o)<<std::endl;
    #ifdef PRINT2B2
        #undef PRINT2B2
    #endif
    #define PRINT2B2(t,o) std::cout<<(t)<<(o)<<std::endl; rpt::echo<<(t)<<(o)<<std::endl;

    adstring model  = tcsam::MODEL;
    adstring modVer = tcsam::VERSION; 
    
    time_t start,finish;
    
    //model objects
    ModelConfiguration*  ptrMC; //ptr to model configuration object
    ModelParametersInfo* ptrMPI;//ptr to model parameters info object
    ModelOptions*        ptrMOs;//ptr to model options object
    ModelDatasets*       ptrMDS;//ptr to model datasets object
    ModelDatasets*       ptrSimMDS;//ptr to simulated model datasets object
    OFLResults*          ptrOFLResults;//ptr to OFL results object for MCMC calculations
    
    //MSE objects
    MSE_OpModInfo* ptrOMI=0; //ptr to MSE OpModInfo object
    
    //population dynamics objects for alternative calculations
    int runAlt = 0;
    PopDyInfo*    pPDI=0;//  pointer to population dynamics info
    CatchInfo*    pCDI=0;//  pointer to catch info
    PopProjector* pPPr=0;//  pointer to population projector
    
    //MSE OpMod objects
    PopProjector** ppOpModPPr_x=0;//pointer to sex-specific array of population projectors 
          
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
    long ctrMCMC = 0;    //counter for mcmc output
    std::ofstream mcmc;  //stream for mcmc output
    
    //filenames
    adstring fnMCMC = "tcsam02.MCMC.R";
    adstring fnConfigFile;//configuration file
    adstring fnPin;       //pin file
    
    //runtime flags (0=false)
    int jitter     = 0;//use jittering for initial parameter values
    int resample   = 0;//use resampling for initial parameter values
    int mcevalOn   = 0;//flag indicating model is being run in mceval phase
    int mseMode    = 0;//flag indicating model is being run in an MSE
    int mseOpModMode  = 0;//flag indicating model is being run in an MSE in operating model mode
    int mseEstModMode = 0;//flag indicating model is being run in an MSE in estimation model mode 
    int usePin     = 0;//flag to initialize parameter values using a pin file
    int doRetro    = 0;//flag to facilitate a retrospective model run
    int fitSimData = 0;//flag to fit model to simulated data calculated in the PRELIMINARY_CALCs section
    int doOFL      = 0;///<flag (0/1) to do OFL calculations
    int doTAC      = 0;///<calculate TAC using harvest control rule indicated by value of doTAC
    int doDynB0    = 0;//flag to run dynamic B0 calculations after final phase
        
    int yRetro = 0; //number of years to decrement for retrospective model run
    int iSeed = -1; //default random number generator seed
    random_number_generator rng(iSeed);//random number generator
    int iSimDataSeed = 0;
    random_number_generator rngSimData(-1);//random number generator for data simulation
    
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
    
    int ctrDebugParams = 0;//PROCEDURE_SECTION call counter value to start debugging output
    
    int showActiveParams = 0;    
    int debugRunModel    = 0;    
    int debugObjFun      = 0;
    int debugOFL         = 0;
    int debugMCMC        = 0;
    
    //note: consider using std::bitset to implement debug functionality
    int dbgCalcProcs = 10;
    int dbgObjFun = 20;
    int dbgNLLs   = 25;
    int dbgPriors = tcsam::dbgPriors;
    int dbgPopDy  = 70;
    int dbgApply  = 80;
    int dbgDevs   = 90;
    int dbgAll    = tcsam::dbgAll;
    
    int phsItsRewgt  = 1000;//min phase to calculate effective weights for size comps
    int maxItsRewgt = 0;    //maximum number of terations for re-weighting
    int numItsRewgt = 0;    //number of re-weighting iterations completed
        
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
    {
        adstring msg = "#Starting "+model+" (ver "+modVer+") Code";
        PRINT2B1(msg)
        PRINT2B1("#Starting DATA_SECTION")
    }
 END_CALCS

 //Set commandline options
 LOCAL_CALCS
    int on = 0;
    int flg = 0;
    PRINT2B1("#------Reading command line options---------")
    //configFile
    fnConfigFile = "tcsam02.ModelConfig.dat";//default model config filename
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-configFile"))>-1) {
        fnConfigFile = ad_comm::argv[on+1];
        rpt::echo<<"#config file changed to '"<<fnConfigFile<<"'"<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg=1;
    }
    //parameter input file
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-pin"))>-1) {
        usePin = 1;
        fnPin = ad_comm::argv[on+1];
        if (std::fstream(fnPin)){
            ad_comm::change_pinfile_name(fnPin);
            rpt::echo<<"#Initial parameter values from pin file: "<<fnPin<<endl;
        } else {
            rpt::echo<<"#Initial parameter values from pin file: tcsam02.pin"<<endl;
        }
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
    //mceval phase is on
    if ((option_match(ad_comm::argc,ad_comm::argv,"-mceval")>-1)){
        mcevalOn = 1;
        rpt::echo<<"#mceval is ON."<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //parameter input file for running NUTS MCMC 
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-mcpin"))>-1) {
        usePin = 2;
        fnPin = ad_comm::argv[on+1];
        rpt::echo<<"#Initial parameter values for running NUTS MCMC from pin file: "<<fnPin<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
    }
    //runAlt
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-runAlt"))>-1) {
        runAlt=1;
        rpt::echo<<"#run alternative population dynamics functions"<<endl;
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
    //calcTAC
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcTAC"))>-1) {
        doTAC = atoi(ad_comm::argv[on+1]);
        doOFL = 1;
        rpt::echo<<"#Calculating TAC using harvest control rule "<<doTAC<<endl;
        rpt::echo<<"#OFL calculations also turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //run dynamic B0 calculations after final phase
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-calcDynB0"))>-1) {
        doDynB0 = 1;
        rpt::echo<<"#Running dynamic B0 calculations after final phase."<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
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
    //ctrDebugParams
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-ctrDebugParams"))>-1) {
        ctrDebugParams=atoi(ad_comm::argv[on+1]);
        rpt::echo<<"Starting debugParams at counter "<<ctrDebugParams<<endl;
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
        ofstream fs("jitterInfo.dat");
        fs<<"seed = "<<iSeed<<endl;
        fs.close();
        flg = 1;
    }
    //resample
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-resample"))>-1) {
        resample=1;
        rpt::echo<<"#Resampling for initial parameter values turned ON "<<endl;
        rpt::echo<<"#-------------------------------------------"<<endl;
        flg = 1;
    }
    //mseMode
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-mseMode"))>-1) {
        mseMode=1;
        adstring type = ad_comm::argv[on+1];
        if (type=="mseOpModMode"){
            mseOpModMode = 1;
            rpt::echo<<"#MSE operating model mode turned ON"<<endl;
            iSeed=(long)start;
            if ((on=option_match(ad_comm::argc,ad_comm::argv,"-iSeed"))>-1) {
                if (on+1<argc) {
                    iSeed=atoi(ad_comm::argv[on+1]);
                }
            } 
            rng.reinitialize(iSeed);
            rpt::echo<<tb<<iSeed<<"  #iSeed used for random recruitment"<<endl;
        } else if (type=="mseEstModMode"){
            mseEstModMode = 1;
            rpt::echo<<"#MSE estimation model mode turned ON"<<endl;
        }
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
    //debugOFL
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-debugOFL"))>-1) {
        debugOFL=1;
        rpt::echo<<"#debugOFL turned ON"<<endl;
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
    //debuMCMC
    if (option_match(ad_comm::argc,ad_comm::argv,"-debugMCMC")>-1) {
        debugMCMC=1;
        rpt::echo<<"#debugMCMC turned ON"<<endl;
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
    if (mseOpModMode&&mseEstModMode){
        PRINT2B1("");
        PRINT2B1("--------ERROR!-----------")
        PRINT2B1("mseOpModMode and mseEstModMode cannot both be 'on'.")
        PRINT2B1("Terminating model run!!")
        PRINT2B1("--------ERROR!-----------")
        PRINT2B1("")
        exit(1);
    }
 END_CALCS
 
    int nZBs;  //number of model size bins
    int mnYr;  //min model year
    int mxYr;  //max model year
    int mxYrp1;//max model year + 1
    int nFsh;  //number of fisheries
    int nSrv;  //number of surveys
 LOCAL_CALCS
    PRINT2B1("#-----------------------------------")
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#--Reading configuration file ",fnConfigFile)
    ad_comm::change_datafile_name(fnConfigFile);
    ptrMC = new ModelConfiguration();
    ptrMC->read(*(ad_comm::global_datafile));
    PRINT2B1("#--Finished reading configuration file")
    
    mnYr   = ptrMC->mnYr;
    mxYr   = ptrMC->mxYr;
    if (doRetro){mxYr = mxYr-yRetro; ptrMC->setMaxModelYear(mxYr);}
    if (jitter)   {ptrMC->jitter   = 1;}
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
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#Reading parameters info file ",ptrMC->fnMPI)
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
    PRINT2B1("#----finished reading model parameters info---")
//    rpt::echo<<"#-----------ModelParametersInfo---------------"<<endl;
//    rpt::echo<<(*ptrMPI)<<endl;
//    rpt::echo<<"#----finished ModelParametersInfo---"<<endl;
//    if (debugDATA_SECTION){
//        cout<<"#------------------ModelParametersInfo-----------------"<<endl;
//        cout<<(*ptrMPI);
//        cout<<"#----finished model parameters info---"<<endl;
//        cout<<"enter 1 to continue : ";
//        cin>>debugDATA_SECTION;
//        if (debugDATA_SECTION<0) exit(1);
//    }
//    ptrMPI->project(0);//TODO:testing only!
//    ofstream ofs; ofs.open("projectedMPI.txt");
//    ptrMPI->write(ofs);
//    ofs.close();
//    cout<<"wrote out ptojectedMPI!"<<endl;
//    exit(1);
 END_CALCS
        
    //read model data
 LOCAL_CALCS
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#--Reading datasets file ",ptrMC->fnMDS);
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
    PRINT2B1("#----finished model datasets---")
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
//    rpt::echo<<"---SimMDS object after reading datasets---"<<endl;
//    rpt::echo<<(*ptrSimMDS);
//    rpt::echo<<"#finished SimMDS object"<<endl;
 END_CALCS   
    
    //read model options
 LOCAL_CALCS
    PRINT2B1("#-----------------------------------")
    PRINT2B2("#--Reading model options file ",ptrMC->fnMOs)
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
    PRINT2B1("#--finished reading model options file---")
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
    //check commandline options for iterative reweighting overrides
    //min phase in which to calculate effective weights for size compositions
    phsItsRewgt = ptrMOs->phsIterativeReweighting;
    maxItsRewgt  = ptrMOs->maxIterations;
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-phsIterativeReweighing"))>-1) {
        if (on+1<argc) {
            phsItsRewgt=atoi(ad_comm::argv[on+1]);
        }
        flg = 1;
    }
    if ((on=option_match(ad_comm::argc,ad_comm::argv,"-maxIterationsReweighting"))>-1) {
        if (on+1<argc) {
            phsItsRewgt=atoi(ad_comm::argv[on+1]);
        }
        flg = 1;
    }
    rpt::echo<<phsItsRewgt<<tb<<"#phsIterativeReweighing"<<endl;
    rpt::echo<<maxItsRewgt<<tb<<"#maxIterationsReweighting"<<endl;
    rpt::echo<<"#-------------------------------------------"<<endl;
 END_CALCS
        
    //Match up model fisheries with fisheries datasets
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
     PRINT2B2("model fisheries map to fishery data objects: ",mapM2DFsh)
    }
 END_CALCS
        
    //Match up model surveys with surveys datasets
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
     PRINT2B2("model surveys map to survey data objects: ",mapM2DSrv)
    }
 END_CALCS
 
    //Match up model surveys with chela height datasets
    //Assign size bin indices to observations
    ivector mapD2MChd(1,ptrMDS->nCHD);
 LOCAL_CALCS
    {
     int idx;
     for (int v=1;v<=ptrMDS->nCHD;v++){
         idx = wts::which(ptrMDS->ppCHD[v-1]->survey,ptrMC->lblsSrv);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying survey names and labels in CH data file and config file."<<endl;
             cout<<"Incorrect survey name in CH data file is '"<<ptrMDS->ppCHD[v-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MChd(v)   = idx;//map from chela height dataset object v to model survey index
         ptrMDS->ppCHD[v-1]->calcSizeBinIndices(ptrMC->zCutPts);
     }
     PRINT2B2("chela height datasets map to model surveys: ",mapD2MChd)
    }
 END_CALCS
 
    //Match up model surveys with maturity ogive datasets
    //Assign size bin indices to observations
    ivector mapD2MMOd(1,ptrMDS->nMOD);
 LOCAL_CALCS
    {
     int idx;
     for (int v=1;v<=ptrMDS->nMOD;v++){
         idx = wts::which(ptrMDS->ppMOD[v-1]->survey,ptrMC->lblsSrv);
         if (idx<1){
             cout<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
             cout<<"Error specifying survey names and labels in MO data file and config file."<<endl;
             cout<<"Incorrect survey name in CH data file is '"<<ptrMDS->ppMOD[v-1]<<"'"<<endl;
             cout<<"Please fix names in files!!"<<endl;
             exit(-1);
         }
         mapD2MMOd(v)   = idx;//map from maturity ogive dataset object v to model survey index
         ptrMDS->ppMOD[v-1]->calcSizeBinRemapper(ptrMC->zCutPts);
     }
     PRINT2B2("maturity ogive datasets map to model surveys: ",mapD2MMOd)
    }
 END_CALCS
 
    //Extract parameter information
    //recruitment parameters
    int npLnR; ivector phsLnR; vector lbLnR; vector ubLnR;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pLnR,npLnR,lbLnR,ubLnR,phsLnR,rpt::echo);
    
    int npRCV; ivector phsRCV; vector lbRCV; vector ubRCV;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pRCV,npRCV,lbRCV,ubRCV,phsRCV,rpt::echo);
    
    int npRX; ivector phsRX; vector lbRX; vector ubRX;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pRX,npRX,lbRX,ubRX,phsRX,rpt::echo);
    
    int npRa; ivector phsRa; vector lbRa; vector ubRa;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pRa,npRa,lbRa,ubRa,phsRa,rpt::echo);
    
    int npRb; ivector phsRb; vector lbRb; vector ubRb;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pRb,npRb,lbRb,ubRb,phsRb,rpt::echo);
    
    int npDevsLnR; ivector mniDevsLnR; ivector mxiDevsLnR; imatrix idxsDevsLnR;
    vector lbDevsLnR; vector ubDevsLnR; ivector phsDevsLnR;
    !!tcsam::setParameterInfo(ptrMPI->ptrRec->pDevsLnR,npDevsLnR,mniDevsLnR,mxiDevsLnR,idxsDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR,rpt::echo);
 LOCAL_CALCS    
    if (mseOpModMode) {
        phsLnR = -1;
        phsRCV = -1;
        phsRX  = -1;
        phsRa  = -1;
        phsRb  = -1;
        phsDevsLnR = -1;
    }
 END_CALCS    
    
    //natural mortality parameters
    int npM; ivector phsM; vector lbM; vector ubM;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pM,npM,lbM,ubM,phsM,rpt::echo);
    
    int npDM1; ivector phsDM1; vector lbDM1; vector ubDM1;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pDM1,npDM1,lbDM1,ubDM1,phsDM1,rpt::echo);
    
    int npDM2; ivector phsDM2; vector lbDM2; vector ubDM2;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pDM2,npDM2,lbDM2,ubDM2,phsDM2,rpt::echo);
    
    int npDM3; ivector phsDM3; vector lbDM3; vector ubDM3;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pDM3,npDM3,lbDM3,ubDM3,phsDM3,rpt::echo);
    
    int npDM4; ivector phsDM4; vector lbDM4; vector ubDM4;
    !!tcsam::setParameterInfo(ptrMPI->ptrNM->pDM4,npDM4,lbDM4,ubDM4,phsDM4,rpt::echo);
    
    number zMref;
    !!zMref = ptrMPI->ptrNM->zRef;
 LOCAL_CALCS    
    if (mseOpModMode) {
        phsM = -1;
        phsDM1 = -1;
        phsDM2 = -1;
        phsDM3 = -1;
        phsDM4 = -1;
    }
 END_CALCS    
        
    
    //maturity parameters
    int npLgtPrMat; ivector mniLgtPrMat; ivector mxiLgtPrMat; imatrix idxsLgtPrMat;
    vector lbLgtPrMat; vector ubLgtPrMat; ivector phsLgtPrMat;
    !!tcsam::setParameterInfo(ptrMPI->ptrM2M->pvLgtPrM2M,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,idxsLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat,rpt::echo);
 LOCAL_CALCS    
    if (mseOpModMode) {
        phsLgtPrMat = -1;
    }
 END_CALCS    
        
 
    //growth parameters
    int npGrA; ivector phsGrA; vector lbGrA; vector ubGrA;
    !!tcsam::setParameterInfo(ptrMPI->ptrGrw->pGrA,npGrA,lbGrA,ubGrA,phsGrA,rpt::echo);
    
    int npGrB; ivector phsGrB; vector lbGrB; vector ubGrB;
    !!tcsam::setParameterInfo(ptrMPI->ptrGrw->pGrB,npGrB,lbGrB,ubGrB,phsGrB,rpt::echo);
    
    int npGrBeta; ivector phsGrBeta; vector lbGrBeta; vector ubGrBeta;
    !!tcsam::setParameterInfo(ptrMPI->ptrGrw->pGrBeta,npGrBeta,lbGrBeta,ubGrBeta,phsGrBeta,rpt::echo);
 LOCAL_CALCS    
    if (mseOpModMode) {
        phsGrA = -1;
        phsGrB = -1;
        phsGrBeta = -1;
    }
 END_CALCS    
        
    
    //selectivity parameters
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
    
    int npNPSel; ivector mniNPSel; ivector mxiNPSel; imatrix idxsNPSel;
    vector lbNPSel; vector ubNPSel; ivector phsNPSel;
    !!tcsam::setParameterInfo(ptrMPI->ptrSel->pvNPSel,npNPSel,mniNPSel,mxiNPSel,idxsNPSel,lbNPSel,ubNPSel,phsNPSel,rpt::echo);
 LOCAL_CALCS    
    if (mseOpModMode) {
        phsS1 = -1;
        phsS2 = -1;
        phsS3 = -1;
        phsS4 = -1;
        phsS5 = -1;
        phsS6 = -1;
        phsDevsS1 = -1;
        phsDevsS2 = -1;
        phsDevsS3 = -1;
        phsDevsS4 = -1;
        phsDevsS5 = -1;
        phsDevsS6 = -1;
        phsNPSel  = -1;
    }
 END_CALCS    
        
    //fisheries parameters
    int npHM; ivector phsHM; vector lbHM; vector ubHM;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pHM,npHM,lbHM,ubHM,phsHM,rpt::echo);
    
    int npLnC; ivector phsLnC; vector lbLnC; vector ubLnC;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnC,npLnC,lbLnC,ubLnC,phsLnC,rpt::echo);
    
    int npDC1; ivector phsDC1; vector lbDC1; vector ubDC1;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC1,npDC1,lbDC1,ubDC1,phsDC1,rpt::echo);
    
    int npDC2; ivector phsDC2; vector lbDC2; vector ubDC2;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC2,npDC2,lbDC2,ubDC2,phsDC2,rpt::echo);
    
    int npDC3; ivector phsDC3; vector lbDC3; vector ubDC3;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC3,npDC3,lbDC3,ubDC3,phsDC3,rpt::echo);
    
    int npDC4; ivector phsDC4; vector lbDC4; vector ubDC4;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pDC4,npDC4,lbDC4,ubDC4,phsDC4,rpt::echo);
    
    int npDevsLnC; ivector mniDevsLnC; ivector mxiDevsLnC; imatrix idxsDevsLnC;
    vector lbDevsLnC; vector ubDevsLnC; ivector phsDevsLnC;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pDevsLnC,npDevsLnC,mniDevsLnC,mxiDevsLnC,idxsDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC,rpt::echo);
    
    int npLnEffX; ivector phsLnEffX; vector lbLnEffX; vector ubLnEffX;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLnEffX,npLnEffX,lbLnEffX,ubLnEffX,phsLnEffX,rpt::echo);
    
    int npLgtRet; ivector phsLgtRet; vector lbLgtRet; vector ubLgtRet;
    !!tcsam::setParameterInfo(ptrMPI->ptrFsh->pLgtRet,npLgtRet,lbLgtRet,ubLgtRet,phsLgtRet,rpt::echo);
    
 LOCAL_CALCS    
    if (mseOpModMode) {
        phsHM = -1;
        phsLnC = -1;
        phsDC1 = -1;
        phsDC2 = -1;
        phsDC3 = -1;
        phsDC4 = -1;
        phsDevsLnC = -1;
        phsLnEffX = -1;
        phsLgtRet = -1;
    }
 END_CALCS    
         
    //surveys parameters
    int npQ; ivector phsQ; vector lbQ; vector ubQ;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pQ,npQ,lbQ,ubQ,phsQ,rpt::echo);
    
    int npDQ1; ivector phsDQ1; vector lbDQ1; vector ubDQ1;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ1,npDQ1,lbDQ1,ubDQ1,phsDQ1,rpt::echo);
    
    int npDQ2; ivector phsDQ2; vector lbDQ2; vector ubDQ2;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ2,npDQ2,lbDQ2,ubDQ2,phsDQ2,rpt::echo);
    
    int npDQ3; ivector phsDQ3; vector lbDQ3; vector ubDQ3;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ3,npDQ3,lbDQ3,ubDQ3,phsDQ3,rpt::echo);
    
    int npDQ4; ivector phsDQ4; vector lbDQ4; vector ubDQ4;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pDQ4,npDQ4,lbDQ4,ubDQ4,phsDQ4,rpt::echo);
    
    int npA; ivector phsA; vector lbA; vector ubA;
    !!tcsam::setParameterInfo(ptrMPI->ptrSrv->pA,npA,lbA,ubA,phsA,rpt::echo);
    
 LOCAL_CALCS    
    if (mseOpModMode) {
        phsQ = -1;
        phsDQ1 = -1;
        phsDQ2 = -1;
        phsDQ3 = -1;
        phsDQ4 = -1;
        phsA   = -1;
    }
 END_CALCS    
         
    //MSE-related parameters
    int npMSE_LnC; ivector phsMSE_LnC; vector lbMSE_LnC; vector ubMSE_LnC;
    !!tcsam::setParameterInfo(ptrMPI->ptrMSE->pMSE_LnC,npMSE_LnC,lbMSE_LnC,ubMSE_LnC,phsMSE_LnC,rpt::echo);
    !!if (!mseOpModMode) {phsMSE_LnC = -1;}
    !!if ( mseOpModMode) {phsMSE_LnC =  1;}
     
    //other data
    vector dtF_y(mnYr,mxYr);//timing of midpoint of fishing season (by year)
    !!dtF_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);
    vector dtM_y(mnYr,mxYr);//timing of mating (by year))
    !!dtM_y = ptrMDS->ptrBio->fshTiming_y(mnYr,mxYr);

    //Model Options
    //--effort averaging scenarios
    int nEASs;//number of effort averaging scenarios
    !!nEASs = ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios->nAvgs;
    imatrix yrsAvgEff_ny(1,nEASs,mnYr,mxYr);//years for effort averaging scenarios
    matrix eff_ny(1,nEASs,mnYr,mxYr);       //effort for averaging over fishery-specific time periods
    vector avgEff_n(1,nEASs);               //average effort over fishery-specific time periods
    //--capture rate averaging scenarios
    int nCRASs;//number of effort averaging scenarios
    !!nCRASs = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios->nAvgs;
    5darray obsEff_nxmsy(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr);         //observed effort for averaging by capture rate averaging scenario

    //number of parameter combinations for various processes
    int npcRec;   //number of recruitment parameter combinations
    !!npcRec = ptrMPI->ptrRec->nPCs;
    int npcNM;    //number of natural mortality parameter combinations
    !!npcNM = ptrMPI->ptrNM->nPCs;
    int npcM2M;   //number of molt-to-maturity parameter combinations
    !!npcM2M = ptrMPI->ptrM2M->nPCs;
    int npcGrw;   //number of growth parameter combinations
    !!npcGrw = ptrMPI->ptrGrw->nPCs;
    int npcSel;   //number of selectivity functions (by parameter combination)
    !!npcSel = ptrMPI->ptrSel->nPCs;
    int npcFsh;   //number of fishery parameter combinations
    !!npcFsh = ptrMPI->ptrFsh->nPCs;
    int npcSrv;   //number of survey parameter combinations
    !!npcSrv = ptrMPI->ptrSrv->nPCs;
    
    ivector nDevsLnR_c(1,npcRec);              //number of recruit devs by parameter combination
    
    //growth arrays
    matrix zGrA_xy(1,nSXs,mnYr,mxYr);//pre-molt size corresponding to pGrA in alt growth parameterization
    matrix zGrB_xy(1,nSXs,mnYr,mxYr);//pre-molt size corresponding to pGrB in alt growth parameterization
    
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
    
    //MSE-related variables
    number opModMnLnR;
    number opModSdLnR;
    number opModDevLnR;
    number prjR;
    number inpOFL; //OFL for upcoming year
    number inpTAC; //TAC for upcoming year
 LOCAL_CALCS
    if (mseOpModMode){
        PRINT2B1("#--Creating ptrOMI")
        ptrOMI = new MSE_OpModInfo(ptrMC);
        ad_comm::change_datafile_name("OpModStateFile.txt");
        PRINT2B1("#--Reading OpModStateFile.txt")
        (*ad_comm::global_datafile)>>(*ptrOMI);
        ptrOMI->write(rpt::echo);        
        PRINT2B1("#--Finished reading OpModStateFile.txt")
        prjR = 0.0;
        PRINT2B1("#--Reading TAC.txt")
        ad_comm::change_datafile_name("TAC.txt");
        (*ad_comm::global_datafile)>>inpTAC;
        (*ad_comm::global_datafile)>>inpOFL;
        PRINT2B1("#--Finished reading TAC.txt")
        rpt::echo<<"TAC, OFL is: "<<inpTAC<<tb<<inpOFL<<endl;
        cout<<"TAC, OFL is: "<<inpTAC<<tb<<inpOFL<<endl;
    }
 END_CALCS
 
    !!PRINT2B1("#finished DATA_SECTION")
            
// =============================================================================
// =============================================================================
INITIALIZATION_SECTION

// =============================================================================
// =============================================================================
PARAMETER_SECTION
    !!PRINT2B1("#Starting PARAMETER_SECTION")
 
    //recruitment parameters
    init_bounded_number_vector pLnR(1,npLnR,lbLnR,ubLnR,phsLnR);    //mean ln-scale recruitment
    !!cout<<"pLnR = "<<pLnR<<endl;
    init_bounded_number_vector pRCV(1,npRCV,lbRCV,ubRCV,phsRCV);    //ln-scale recruitment cv
    !!cout<<"pRCV = "<<pRCV<<endl;
    init_bounded_number_vector pRX(1,npRX,lbRX,ubRX,phsRX);         //logit-scale male sex ratio
    !!cout<<"pRX = "<<pRX<<endl;
    init_bounded_number_vector pRa(1,npRa,lbRa,ubRa,phsRa);         //size distribution parameter
    !!cout<<"pRa = "<<pRa<<endl;
    init_bounded_number_vector pRb(1,npRb,lbRb,ubRb,phsRb);         //size distribution parameter
    !!cout<<"pRb = "<<pRb<<endl;
    init_bounded_vector_vector pDevsLnR(1,npDevsLnR,mniDevsLnR,mxiDevsLnR,lbDevsLnR,ubDevsLnR,phsDevsLnR);//ln-scale rec devs
    !!for (int p=1;p<=npDevsLnR;p++) cout<<"pDevsLnR["<<p<<"] = "<<pDevsLnR[p]<<endl;
    matrix devsLnR(1,npDevsLnR,mniDevsLnR,mxiDevsLnR+1);
    !!cout<<"got past recruitment parameters"<<endl;
   
    //natural mortality parameters
    init_bounded_number_vector pM(1,npM,lbM,ubM,phsM);//base ln-scale
    !!cout<<"pM = "<<pM<<endl;
    init_bounded_number_vector pDM1(1,npDM1,lbDM1,ubDM1,phsDM1);//offset 1s
    !!cout<<"pDM1 = "<<pDM1<<endl;
    init_bounded_number_vector pDM2(1,npDM2,lbDM2,ubDM2,phsDM2);//offset 2s
    !!cout<<"pDM2 = "<<pDM2<<endl;
    init_bounded_number_vector pDM3(1,npDM3,lbDM3,ubDM3,phsDM3);//offset 3s
    !!cout<<"pDM3 = "<<pDM3<<endl;
    init_bounded_number_vector pDM4(1,npDM4,lbDM4,ubDM4,phsDM4);//offset 4s
    !!cout<<"pDM4 = "<<pDM4<<endl;
    !!cout<<"got past natural mortality parameters"<<endl;
    
    //growth parameters
    init_bounded_number_vector pGrA(1,npGrA,lbGrA,ubGrA,phsGrA); //ln-scale mean growth coefficient "a"
    !!cout<<"pGrA = "<<pGrA<<endl;
    init_bounded_number_vector pGrB(1,npGrB,lbGrB,ubGrB,phsGrB); //ln-scale mean growth coefficient "b"
    !!cout<<"pGrB = "<<pGrB<<endl;
    init_bounded_number_vector pGrBeta(1,npGrBeta,lbGrBeta,ubGrBeta,phsGrBeta);//ln-scale growth scale parameter
    !!cout<<"pGrBeta = "<<pGrBeta<<endl;
    !!cout<<"got past growth parameters"<<endl;
    
    //maturity parameters
    init_bounded_vector_vector pvLgtPrM2M(1,npLgtPrMat,mniLgtPrMat,mxiLgtPrMat,lbLgtPrMat,ubLgtPrMat,phsLgtPrMat);//logit-scale maturity ogive parameters
    !!for (int p=1;p<=npLgtPrMat;p++) cout<<"pvLgtPrM2M["<<p<<"] = "<<pvLgtPrM2M[p]<<endl;
    !!cout<<"got past maturity parameters"<<endl;
    
    //selectivity parameters
    init_bounded_number_vector pS1(1,npS1,lbS1,ubS1,phsS1);
    !!cout<<"pS1 = "<<pS1<<endl;
    init_bounded_number_vector pS2(1,npS2,lbS2,ubS2,phsS2);
    !!cout<<"pS2 = "<<pS2<<endl;
    init_bounded_number_vector pS3(1,npS3,lbS3,ubS3,phsS3);
    !!cout<<"pS3 = "<<pS3<<endl;
    init_bounded_number_vector pS4(1,npS4,lbS4,ubS4,phsS4);
    !!cout<<"pS4 = "<<pS4<<endl;
    init_bounded_number_vector pS5(1,npS5,lbS5,ubS5,phsS5);
    !!cout<<"pS5 = "<<pS5<<endl;
    init_bounded_number_vector pS6(1,npS6,lbS6,ubS6,phsS6);
    !!cout<<"pS6 = "<<pS6<<endl;
    init_bounded_vector_vector pDevsS1(1,npDevsS1,mniDevsS1,mxiDevsS1,lbDevsS1,ubDevsS1,phsDevsS1);
    !!for (int p=1;p<=npDevsS1;p++) cout<<"pDevsS1["<<p<<"] = "<<pDevsS1[p]<<endl;
    init_bounded_vector_vector pDevsS2(1,npDevsS2,mniDevsS2,mxiDevsS2,lbDevsS2,ubDevsS2,phsDevsS2);
    !!for (int p=1;p<=npDevsS2;p++) cout<<"pDevsS2["<<p<<"] = "<<pDevsS2[p]<<endl;
    init_bounded_vector_vector pDevsS3(1,npDevsS3,mniDevsS3,mxiDevsS3,lbDevsS3,ubDevsS3,phsDevsS3);
    !!for (int p=1;p<=npDevsS3;p++) cout<<"pDevsS3["<<p<<"] = "<<pDevsS3[p]<<endl;
    init_bounded_vector_vector pDevsS4(1,npDevsS4,mniDevsS4,mxiDevsS4,lbDevsS4,ubDevsS4,phsDevsS4);
    !!for (int p=1;p<=npDevsS4;p++) cout<<"pDevsS4["<<p<<"] = "<<pDevsS4[p]<<endl;
    init_bounded_vector_vector pDevsS5(1,npDevsS5,mniDevsS5,mxiDevsS5,lbDevsS5,ubDevsS5,phsDevsS5);
    !!for (int p=1;p<=npDevsS5;p++) cout<<"pDevsS5["<<p<<"] = "<<pDevsS5[p]<<endl;
    init_bounded_vector_vector pDevsS6(1,npDevsS6,mniDevsS6,mxiDevsS6,lbDevsS6,ubDevsS6,phsDevsS6);
    !!for (int p=1;p<=npDevsS6;p++) cout<<"pDevsS6["<<p<<"] = "<<pDevsS6[p]<<endl;
    matrix devsS1(1,npDevsS1,mniDevsS1,mxiDevsS1+1);
    matrix devsS2(1,npDevsS2,mniDevsS2,mxiDevsS2+1);
    matrix devsS3(1,npDevsS3,mniDevsS3,mxiDevsS3+1);
    matrix devsS4(1,npDevsS4,mniDevsS4,mxiDevsS4+1);
    matrix devsS5(1,npDevsS5,mniDevsS5,mxiDevsS5+1);
    matrix devsS6(1,npDevsS6,mniDevsS6,mxiDevsS6+1);
    
    init_bounded_vector_vector pvNPSel(1,npNPSel,mniNPSel,mxiNPSel,lbNPSel,ubNPSel,phsNPSel);
    !!for (int p=1;p<=npNPSel;p++) cout<<"pvNPSel["<<p<<"] = "<<pvNPSel[p]<<endl;
    !!cout<<"got past selectivity parameters"<<endl;
    
    //fishing capture rate parameters
    init_bounded_number_vector pHM(1,npHM,lbHM,ubHM,phsHM);     //handling mortality
    !!cout<<"pHM = "<<pHM<<endl;
    init_bounded_number_vector pLnC(1,npLnC,lbLnC,ubLnC,phsLnC);//ln-scale base fishing mortality (e.g., mature males)
    !!cout<<"pLnC = "<<pLnC<<endl;
    init_bounded_number_vector pDC1(1,npDC1,lbDC1,ubDC1,phsDC1);//ln-offset 1 (e.g., year-block offsets)
    !!cout<<"pDC1 = "<<pDC1<<endl;
    init_bounded_number_vector pDC2(1,npDC2,lbDC2,ubDC2,phsDC2);//ln-offset 2 (e.g., female offsets)
    !!cout<<"pDC2 = "<<pDC2<<endl;
    init_bounded_number_vector pDC3(1,npDC3,lbDC3,ubDC3,phsDC3);//ln-offset 3 (e.g., immature offsets)
    !!cout<<"pDC3 = "<<pDC3<<endl;
    init_bounded_number_vector pDC4(1,npDC4,lbDC4,ubDC4,phsDC4);//ln-offset 4 (e.g., female-immature offsets)
    !!cout<<"pDC4 = "<<pDC4<<endl;
    init_bounded_vector_vector pDevsLnC(1,npDevsLnC,mniDevsLnC,mxiDevsLnC,lbDevsLnC,ubDevsLnC,phsDevsLnC);//ln-scale deviations
    !!for (int p=1;p<=npDevsLnC;p++) cout<<"pDevsLnC["<<p<<"] = "<<pDevsLnC[p]<<endl;
    !!cout<<npLnEffX<<tb<<lbLnEffX<<tb<<ubLnEffX<<tb<<phsLnEffX<<endl;
    init_bounded_number_vector pLnEffX(1,npLnEffX,lbLnEffX,ubLnEffX,phsLnEffX);//ln-scale effort extrapolation parameters
    !!cout<<"pLnEffX = "<<pLnEffX<<endl;
    !!cout<<npLgtRet<<tb<<lbLgtRet<<tb<<ubLgtRet<<tb<<phsLgtRet<<endl;
    init_bounded_number_vector pLgtRet(1,npLgtRet,lbLgtRet,ubLgtRet,phsLgtRet);//lgt-scale retained fraction parameters
    !!cout<<"pLgtRet = "<<pLgtRet<<endl;
    matrix devsLnC(1,npDevsLnC,mniDevsLnC,mxiDevsLnC+1);
    !!cout<<"got past capture rate parameters"<<endl;
   
    //survey catchability parameters
    init_bounded_number_vector pQ(1,npQ,lbQ,ubQ,phsQ);//base (e.g., mature male)
    !!cout<<"pQ = "<<pQ<<endl;
    init_bounded_number_vector pDQ1(1,npDQ1,lbDQ1,ubDQ1,phsDQ1);//ln-offset 1 (e.g., main temporal offsets)
    !!cout<<"pDQ1 = "<<pDQ1<<endl;
    init_bounded_number_vector pDQ2(1,npDQ2,lbDQ2,ubDQ2,phsDQ2);//ln-offset 2 (e.g., female offsets)
    !!cout<<"pDQ2 = "<<pDQ2<<endl;
    init_bounded_number_vector pDQ3(1,npDQ3,lbDQ3,ubDQ3,phsDQ3);//ln-offset 3 (e.g., immature offsets)
    !!cout<<"pDQ3 = "<<pDQ3<<endl;
    init_bounded_number_vector pDQ4(1,npDQ4,lbDQ4,ubDQ4,phsDQ4);//ln-offset 4 (e.g., female-immature offsets)
    !!cout<<"pDQ4 = "<<pDQ4<<endl;
    init_bounded_number_vector pA(1,npA,lbA,ubA,phsA);//max availability
    !!cout<<"pA = "<<pA<<endl;
    !!cout<<"got past survey catchability parameters"<<endl;
   
    //MSE-related parameters
    init_bounded_number_vector pMSE_LnC(1,npMSE_LnC,lbMSE_LnC,ubMSE_LnC,phsMSE_LnC);//
    !!cout<<"pMSE_LnC = "<<pMSE_LnC<<endl;
    !!cout<<"got past MSE parameters"<<endl;
   
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
    vector stdvDevsLnR_c(1,npcRec);            //ln-scale recruitment std. devs by parameter combination
    matrix devsLnR_cy(1,npcRec,mnYr,mxYr);     //unstandardized ln-scale recruitment residuals, by parameter combination and year
    matrix zscrDevsLnR_cy(1,npcRec,mnYr,mxYr); //standardized ln-scale recruitment residuals, by parameter combination and year
    
    //natural mortality-related quantities
    vector M_c(1,npcNM);                                   //natural mortality rate by parameter combination
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
    matrix npSel_cz(1,npNPSel,1,nZBs);           //nonparametric  selectivity functions (fisheries and surveys) by parameter combination)
    matrix sel_cz(1,npcSel,1,nZBs);              //all selectivity functions (fisheries and surveys) by parameter combination (no devs))
    3darray sel_cyz(1,npcSel,mnYr,mxYr+1,1,nZBs);//all selectivity functions (fisheries and surveys) by year
    
    //fishery-related quantities
    4darray avgFc_nxms(1,nCRASs,1,nSXs,1,nMSs,1,nSCs);         //avg capture rates for capture rate averaging scenarios
    4darray avgFc2Eff_nxms(1,nCRASs,1,nSXs,1,nMSs,1,nSCs);     //ratios of avg capture rate to effort for capture rate averaging scenarios
    5darray obsFc_nxmsy(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr);    //"observed" capture rates for capture rate averaging scenarios
    5darray prdFc_nxmsy(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr);    //"predicted" capture rates for capture rate averaging scenarios
    5darray prdEff_nxmsy(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr);   //"predicted" effort for capture rate averaging scenarios
    5darray zscrEffX_nxmsy(1,nCRASs,1,nSXs,1,nMSs,1,nSCs,mnYr,mxYr); //effort z-scores for effort extrapolation/capture rate averaging scenarios
    4darray nllEffX_nxms(1,nCRASs,1,nSXs,1,nMSs,1,nSCs);             //nlls for effort extrapolation/capture rate averaging scenarios
    
    matrix  hmF_fy(1,nFsh,mnYr,mxYr);                        //handling mortality
    matrix dvsLnC_fy(1,nFsh,mnYr,mxYr);                      //matrix to capture lnC-devs
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
    !!cout<<"got here 12"<<endl;
    
    6darray cpB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//catch at size in fishery f
    6darray dsB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discarded catch (numbers, NOT mortality) at size in fishery f
    6darray rmB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//retained catch (numbers=mortality) at size in fishery f
    6darray dmB_fyxmsz(1,nFsh,mnYr,mxYr,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//discard catch mortality at size in fishery f
    !!cout<<"got here 13"<<endl;
    
    //survey-related quantities
    3darray mb_vyx(1,nSrv,mnYr,mxYr+1,1,nSXs);//mature (spawning) biomass at time of survey
    5darray a_vyxms(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs);        //max availability in survey v
    5darray q_vyxms(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs);        //fully-selected catchability in survey v
    6darray a_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//availability functions
    6darray s_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//selectivity functions
    6darray q_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//size-specific catchability in survey v (includes availability)
    6darray n_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//catch abundance at size in survey v
    6darray b_vyxmsz(1,nSrv,mnYr,mxYr+1,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//catch biomass at size in survey v

    //MSE-related quantities
    vector prj_hmF_f(1,nFsh);                         //projected handling mortality
    4darray prj_capF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//projected fishery capture rates
    4darray prj_retF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//projected fishery retention functions
    4darray prj_selF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//projected fishery selectivity functions
    vector prj_spB_x(1,nSXs);                                 //projected spawning biomass
    4darray prj_n_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);          //projected population size at beginning of new year
    5darray prj_cpN_fxmsz(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//projected capture abundance
    5darray prj_rmN_fxmsz(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//projected retained mortality (abundance)
    5darray prj_dmN_fxmsz(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//projected discards mortality (abundance)
    5darray prj_dsN_fxmsz(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);//projected discards (abundance)
    matrix prjRetCatchMortBio_fx(1,nFsh,1,nSXs);              //projected retained catch mortality (biomass))
    matrix prjDscCatchMortBio_fx(1,nFsh,1,nSXs);              //projected discard catch mortality (biomass)
    matrix prjTotCatchMortBio_fx(1,nFsh,1,nSXs);              //projected total catch mortality (biomass))
    5darray prj_n_vxmsz(1,nSrv,1,nSXs,1,nMSs,1,nSCs,1,nZBs);  //projected surveys abundance

    //objective function penalties
    vector fPenRecDevs(1,npDevsLnR);//recruitment devs penalties
    
    vector fPenSmoothLgtPrMat(1,npLgtPrMat);//smoothness penalties on pr(mature|z)
    vector fPenNonDecLgtPrMat(1,npLgtPrMat);//non-decreasing penalties on pr(mature|z)
    
    vector fPenSmoothNPSel(1,npNPSel);//smoothness penalties on nonparametric selectivities
    
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
    
    //sdreport variables
    sdreport_vector sdrLnR_y(mnYr,mxYr);
    sdreport_matrix sdrSpB_xy(1,nSXs,mnYr,mxYr);
    
    //likelihood profile numbers
    likeprof_number lkMMB; 
    likeprof_number lkQM; 
    likeprof_number lkQF; 
    likeprof_number lkNMIM; 
    likeprof_number lkNMIF;
    likeprof_number lkNMMM; 
    likeprof_number lkNMMF;

    
    !!PRINT2B1("#finished PARAMETER_SECTION")
    
// =============================================================================
// =============================================================================
PRELIMINARY_CALCS_SECTION
    PRINT2B1("#Starting PRELIMINARY_CALCS_SECTION")
    int debug=1;
    
//---------PRELIMINARY_CALCS: ALT POP DY CALCULATIONS---------------------------
    if (runAlt){
        //create population projection objects
        PRINT2B1("#--Creating population projection objects")
        pPDI = new PopDyInfo(nZBs);         //  generic population dynamics info
        pCDI = new CatchInfo(nZBs,nFsh);    //  generic catch info
        pPPr = new PopProjector(pPDI,pCDI); //  generic population projector
        PRINT2B1("#--Created population projection objects")
    }

//---------PRELIMINARY_CALCS: ALL MODEL RUNS------------------------------------    
    PRINT2B1("PRELIMINARY_CALCS: ALL MODEL RUNS")
    //set initial values for all parameters
    if (usePin) {
        PRINT2B1("NOTE: setting initial values for parameters using pin file")
    } else {
        PRINT2B1("NOTE: setting initial values for parameters using MPI")
    }
    setInitVals(0,rpt::echo);
    
    //calculate average effort for fisheries over specified time periods and 
    //allocate associated arrays
    PRINT2B1(" ")
    PRINT2B1("--calculating average effort")
    rpt::echo<<"mapD2MFsh = "<<mapD2MFsh<<endl;    
    rpt::echo<<"nEASs     = "<<nEASs<<endl;    
    if (nEASs){
        EffAvgScenarios* ptrEASs = ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios;
        for (int n=1;n<=nEASs;n++){//effort averaging scenarios
            EffAvgScenario* ptrEAS = ptrEASs->ppEASs[n-1];
            int fm = ptrEAS->f;    //index for fishery associated with this averaging scenario
            int fd = mapM2DFsh(fm);//index for corresponding fishery data object
            rpt::echo<<"n = "<<n<<". fm = "<<fm<<". fd = "<<fd<<endl;
            if (!ptrMDS->ppFsh[fd-1]->ptrEff){
                cout<<"---------------------------------------------------"<<endl;
                cout<<"No effort data given for "<<ptrMC->lblsFsh[fm]<<","<<endl;
                cout<<"but effort averaging requested! Aborting..."<<endl;
                exit(-1);
            }
            rpt::echo<<"fishery has effort"<<endl;
            //extract years for effort averaging, keeping model limits in mind
            yrsAvgEff_ny(n).deallocate();
            yrsAvgEff_ny(n) = wts::extractVector(mnYr,mxYr,ptrEAS->ptrIB->getFwdIndexVector());
            //assign values for effort averaging
            eff_ny(n).deallocate();
            eff_ny(n) = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(yrsAvgEff_ny(n));
            rpt::echo<<"eff_ny(n)       = "<<eff_ny(n)<<endl;
            rpt::echo<<"yrsAvgEff_ny(n) = "<<yrsAvgEff_ny(n)<<endl;
            avgEff_n(n) = sum(eff_ny(n))/yrsAvgEff_ny(n).size();
            rpt::echo<<"avgEff_n(n) = "<<avgEff_n(n)<<endl;
        }//n
        rpt::echo<<"eff_ny       = "<<endl<<eff_ny<<endl;
        rpt::echo<<"yrsAvgEff_ny = "<<endl<<yrsAvgEff_ny<<endl;
        rpt::echo<<"avgEff_n     = "<<avgEff_n<<endl;
        rpt::echo<<"finished calculating average effort scenarios"<<endl;
        cout<<"finished calculating average effort scenarios"<<endl;

        PRINT2B1(" ")
        PRINT2B1("--starting allocation of arrays for average capture rates and effort extrapolation")
        obsEff_nxmsy.initialize();
        for (int n=1;n<=nCRASs;n++){
            CapRateAvgScenario* ptrCRAS = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios->ppCRASs[n-1];
            int idEAS = ptrCRAS->idEffAvgInfo;
            int mnx, mxx, mnm, mxm, mns, mxs;
            mnx = mxx = ptrCRAS->x;//sex index
            if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
            mnm = mxm = ptrCRAS->m;//maturity index
            if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
            mns = mxs = ptrCRAS->s;//shell index
            if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
            for (int x=mnx;x<=mxx;x++){
                for (int m=mnm;m<=mxm;m++){
                    for (int s=mns;s<=mxs;s++) {
                        for (int iy=yrsAvgEff_ny(idEAS).indexmin();iy<=yrsAvgEff_ny(idEAS).indexmax();iy++)
                            obsEff_nxmsy(n,x,m,s)(yrsAvgEff_ny(idEAS,iy)) = eff_ny(idEAS,iy);
                        rpt::echo<<"n = "<<n<<". idEAS = "<<idEAS<<endl;
                        rpt::echo<<"yrsAvgEff_ny(idEAS) = "<<yrsAvgEff_ny(idEAS)<<endl;
                        rpt::echo<<"eff_ny(idEAS)       ="<<eff_ny(idEAS)<<endl;
                        rpt::echo<<"obsEff_nxmsy("<<n<<cc<<x<<cc<<m<<cc<<s<<") = "<<obsEff_nxmsy(n,x,m,s)<<endl;
                    }//s
                }//m
            }//x
        }//n
        PRINT2B1("--finished allocation of arrays for average capture rates and effort extrapolation")
    } else {
        PRINT2B1("--NO effort averaging scenarios defined!")
    }

    if (!mseMode){
//---------PRELIMINARY_CALCS: NON-MSE MODEL RUNS ONLY---------------------------    
        PRINT2B1("--PRELIMINARY_CALCS: NON-MSE MODEL RUNS ONLY")
        {
            PRINT2B1("writing effective MPI after setInitVals")
            ofstream os; os.open("effectiveMPI.dat", ios::trunc);
            os.precision(12);
            ptrMPI->setToWriteVectorInitialValues(true);
            os<<(*ptrMPI)<<endl;
            os.close();
            PRINT2B1("finished writing effective MPI after setInitVals")
        }

        PRINT2B1("testing setAllDevs()")
        setAllDevs(tcsam::dbgAll,rpt::echo);
        PRINT2B1("finished testing setAllDevs()")

        {
            PRINT2B1("writing data to R")
            ofstream os; os.open("ModelData.R", ios::trunc);
            os.precision(12);
            ReportToR_Data(os,0,cout);
            PRINT2B1("finished writing data to R")
        }

        {
            PRINT2B1("writing parameters info to R")
            ofstream os; os.open("ModelParametersInfo.R", ios::trunc);
            os.precision(12);
            ptrMPI->writeToR(os);
            os.close();
            PRINT2B1("finished writing parameters info to R")

            //write initial parameter values to csv
            PRINT2B1("writing parameters info to csv")
            ofstream os1("tcsam02.params.all.init.csv", ios::trunc);
            os1.precision(12);
            writeParameters(os1,0,0);//all parameters
            os1.close();
            ofstream os2("tcsam02.params.active.init.csv", ios::trunc);
            os2.precision(12);
            writeParameters(os2,0,1);//only parameters that will be active (i.e., phase>0)
            os2.close();
            PRINT2B1("finished writing parameters info to csv")
        }

        if (!mcevalOn) {
//-----------PRELIMINARY_CALCS: "ORDINARY" MODEL RUNS ONLY---------------------- 
            //this section runs for an "ordinary" model run
            int dbgLevel = 0; //set at dbgCalcProcs+1 to print to terminal 
            PRINT2B1("testing calcRecruitment():")
            calcRecruitment(dbgLevel,rpt::echo);

            PRINT2B1("testing calcNatMort():")
            calcNatMort(dbgLevel,rpt::echo);

            PRINT2B1("testing calcGrowth():")
            calcGrowth(dbgLevel,rpt::echo);

            PRINT2B1("testing calcPrM2M():")
            calcPrM2M(dbgLevel,rpt::echo);

            PRINT2B1("testing calcSelectivities():")
            calcSelectivities(dbgLevel,rpt::echo);

            PRINT2B1("testing calcFisheryFs():")
            calcFisheryFs(dbgLevel,rpt::echo);

            PRINT2B1("testing calcSurveyQs():")
            calcSurveyQs(dbgLevel,cout);

            if (!runAlt){
                PRINT2B1("testing runPopDyMod():")
                runPopDyMod(dbgLevel,cout);
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
            } else {
                PRINT2B1("testing runAltPopDyMod():")
                runAltPopDyMod(dbgLevel,cout);
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
            }

            if (doOFL&&debugOFL){
                PRINT2B1("Testing OFL calculations")
                ofstream echoOFL; echoOFL.open("calcOFL.init.txt", ios::trunc);
                echoOFL.precision(12);
                echoOFL<<"----Testing calcOFL()"<<endl;
                calcOFL(mxYr+1,debugOFL,echoOFL);//updates ptrOFLResults
                ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
                ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
                echoOFL<<"----Finished testing calcOFL()!"<<endl;
                echoOFL.close();
                PRINT2B1("Finished testing OFL calculations!")
            }

            if (fitSimData){
                int dbgLevel = 0;
                PRINT2B1("creating sim data to fit in model")
                createSimData(dbgLevel,rpt::echo,iSimDataSeed,ptrMDS);//stochastic if iSimDataSeed<>0
                {
                    PRINT2B1("re-writing data to R")
                    ofstream echo1; echo1.open("ModelData.R", ios::trunc);
                    echo1.precision(12);
                    ReportToR_Data(echo1,0,cout);
                }
            }

            {
                PRINT2B1("--Testing calcObjFun()")
                if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
                calcObjFun(dbgAll,rpt::echo);
                PRINT2B2("--Finished testing calcObjFun(): ",objFun)
            }

            {
                //write objective function components only
                PRINT2B1("Writing model fits to R")
                ofstream os0("tcsam02.ModelFits.init.R", ios::trunc);
                os0.precision(12);
                ReportToR_ModelFits(os0,-1.0,0,cout);
                os0.close();
                PRINT2B1("Finished writing model fits to R")
            }

            {
                PRINT2B1("writing model sim data to file")
                int dbgLevel = 0;
                createSimData(dbgLevel,rpt::echo,0,ptrSimMDS);//deterministic
                ofstream echo1; echo1.open("tcsam02.SimData.init.dat", ios::trunc);
                echo1.precision(12);
                writeSimData(echo1,0,rpt::echo,ptrSimMDS);
                echo1.close();
                PRINT2B1("finished writing model sim data to file")
            }
            
            {
                PRINT2B1("writing cohort progression to R")
                ofstream echo1; echo1.open("CohortProgression.init.R", ios::trunc);
                echo1.precision(12);
                ReportToR_CohortProgression(echo1,1,1,1,1,0,cout);
                echo1.close();
                PRINT2B1("finished writing cohort progression to file")
            }
            
            {
                //must do this last because call to calcDynB0 in ReportToR
                //sets fishing mortality to 0 across all fleets
                //and re-runs population model
                PRINT2B1("writing initial model report to R")
                ofstream echo1; echo1.open("tcsam02.init.rep", ios::trunc);
                echo1.precision(12);
                ReportToR(echo1,-1.0,1,cout);
                echo1.close();
                PRINT2B1("finished writing model report to R")
            }        
        } else {
//-----------PRELIMINARY_CALCS: mceval MODEL RUNS ONLY--------------------------   
            writeMCMCHeader();
            PRINT2B1("MCEVAL is on")
        }//if mcevalOn
        PRINT2B1("")
        PRINT2B2("obj fun = ",objFun)
    } else {
//-------PRELIMINARY_CALCS: MSE MODEL RUNS ONLY---------------------------------   
         if (mseOpModMode){
//-----------PRELIMINARY_CALCS: MSE OpMod RUNS ONLY-----------------------------
            PRINT2B1("PRELIMINARY_CALCS: MSE OpModMode")
            dvector vLnR_y = log(ptrOMI->R_y(ptrMOs->opModRecStatsMinYr,ptrMOs->opModRecStatsMaxYr));
            cout<<"vLnR_y = "<<vLnR_y<<endl;
            opModMnLnR = mean(vLnR_y);
            opModSdLnR = sqrt(wts::variance(vLnR_y));
            opModDevLnR = wts::drawSampleNormal(rng, 0.0, opModSdLnR);
            prjR = mfexp(opModMnLnR+opModDevLnR);
            PRINT2B2("mean recruitment: ",opModMnLnR);
            PRINT2B2("stdv recruitment: ",opModSdLnR);
            PRINT2B2("recruitment dev : ",opModDevLnR);
            PRINT2B2("expected total recruitment: ",mfexp(opModMnLnR+square(opModSdLnR)/2.0));
            PRINT2B2("projected total recruitment: ",prjR);
            //create population projection objects
            PRINT2B1("#--MSE OpMod: creating population projection objects")
            dvariable maxCapF = 0.0;
            ppOpModPPr_x = new PopProjector*[nSXs];
            for (int x=1;x<=nSXs;x++){
                //set population info
                PopDyInfo* ptrPDI = new PopDyInfo(nZBs); //generic population dynamics info
                ptrPDI->w_mz  = 1.0*ptrOMI->wAtZ_xmz(x);  //assign weight-at-length
                ptrPDI->R_z   = 1.0*ptrOMI->R_z;          //relative recruitment-at-size
                ptrPDI->M_msz = 1.0*ptrOMI->M_xmsz(x);    //rates of natural mortality
                ptrPDI->T_szz = 1.0*ptrOMI->prGr_xszz(x); //growth transition matrices
                //pr(terminal molt|size)
                for (int s=1;s<=nSCs;s++) ptrPDI->Th_sz(s) = 1.0*ptrOMI->prM2M_xz(x);
                
                //set catch info
                CatchInfo* ptrCDI = new CatchInfo(nZBs,nFsh); //generic catch info
                prj_hmF_f = 1.0*ptrOMI->hmF_f;
                prj_capF_fmsz.initialize();
                prj_retF_fmsz.initialize();
                prj_selF_fmsz.initialize();
                for (int f=1;f<=nFsh;f++){
                    for (int m=1;m<=nMSs;m++){
                        prj_capF_fmsz(f,m) = 1.0*ptrOMI->cpF_fxmsz(f,x,m);
                        prj_retF_fmsz(f,m) = 1.0*ptrOMI->ret_fxmsz(f,x,m);
                        prj_selF_fmsz(f,m) = 1.0*ptrOMI->sel_fxmsz(f,x,m);
                    }//m
                }//f
                ptrCDI->setHandlingMortality(prj_hmF_f);
                ptrCDI->setCaptureRates(prj_capF_fmsz);
                ptrCDI->setRetentionFcns(prj_retF_fmsz);
                ptrCDI->setSelectivityFcns(prj_selF_fmsz);
                if (x==  MALE) maxCapF = ptrCDI->findMaxTargetCaptureRate(cout);//TODO: how to handle FEMALEs
                if (x==FEMALE) ptrCDI->maxF = maxCapF;//need to set for females based on males
                PopProjector* ptrPPr = new PopProjector(ptrPDI,ptrCDI); //generic population projector
                ptrPPr->dtF = 1.0*ptrOMI->dtF; //time at which fisheries occur
                ptrPPr->dtM = 1.0*ptrOMI->dtM; //time at which mating occurs

                ppOpModPPr_x[x-1] = ptrPPr;
            }
            PRINT2B1("#--MSE OpMod: created population projection objects")
            if (inpTAC>0.0){
                dvariable mseCapF = mfexp(pMSE_LnC[1]);
                projectPopForTAC(mseCapF,dbgAll,rpt::echo);
                calcObjFunForTAC(dbgAll,rpt::echo);
            } else {
                PRINT2B1("#--TAC is 0: directed fishery is closed!!");
                PRINT2B1("#--RUNNING OpMod without fitting for TAC.")
                finishOpModMode();
                PRINT2B1("#--Exiting model without fitting for TAC.")
                exit(0);
            }
        } else if (mseEstModMode){
//-----------PRELIMINARY_CALCS: MSE EstMod RUNS ONLY-----------------------------   
             PRINT2B1("PRELIMINARY_CALCS: MSE EstModMode")
        }
    }
    PRINT2B1("#finished PRELIMINARY_CALCS_SECTION")
    PRINT2B1("#----------------------------------")
    
// =============================================================================
// =============================================================================
PROCEDURE_SECTION
    int dbg = 0; //dbgAll;
    
    ctrProcCalls++;       //increment procedure section calls counter
    ctrProcCallsInPhase++;//increment in-phase procedure section calls counter
    if (dbg>=dbgObjFun){
        PRINT2B1(" ")
        PRINT2B1("--PROCEDURE_SECTION----------------")
    }

//    {
//        adstring fn = "tcsam02.active."+str(current_phase())+"."+str(ctrProcCalls)+".csv";
//        ofstream os; os.open((char*) fn, ios::trunc);
//        writeParameters(os,0,1);
//    }
    
    if (mc_phase()&&(option_match(ad_comm::argc,ad_comm::argv,"-nuts")>0)){
        adstring msg;
        PRINT2B1("--Running NUTS MCMC-------")
        msg = "---phase="+str(current_phase())+". ctrProcCalls = "+str(ctrProcCalls)+tb+str(ctrProcCallsInPhase);
        PRINT2B1(msg)
        msg = "number of active params = "+str(initial_params::nvarcalc())+tb+str(initial_params::num_active_calc());
        PRINT2B1(msg);
        rpt::echo<<"number of active params = "<<initial_params::nvarcalc()<<cc<<initial_params::num_active_calc()<<endl;
        //writeParameters(std::cout,0,1);
        writeParameters(rpt::echo,0,1);
    }
    if (ctrDebugParams&&(ctrProcCalls>=ctrDebugParams)){
        adstring msg;
        PRINT2B1("--writing parameters to file for debugging")
        msg = "---phase="+str(current_phase())+". ctrProcCalls = "+str(ctrProcCalls)+tb+str(ctrProcCallsInPhase);
        PRINT2B1(msg)
        writeParameters(std::cout,0,1);
        writeParameters(rpt::echo,0,1);
        adstring fn = "tcsam02.Debug."+str(current_phase())+"."+str(ctrProcCalls)+".rep";
        ofstream os; os.open((char*) fn, ios::trunc);
        os.precision(12);
        PRINT2B1("--writing report to R file for debugging")
        ReportToR(os,-1.0,1,cout);
        PRINT2B1("--finished writing report to R file for debugging")
    }
    
    if (checkParams(0,std::cout)){
        adstring msg;
        PRINT2B1("--checking parameters for debugging")
        msg = "---phase="+str(current_phase())+". ctrProcCalls = "+str(ctrProcCalls)+tb+str(ctrProcCallsInPhase);
        PRINT2B1(msg)
        checkParams(1,std::cout);
        checkParams(1,rpt::echo);
        writeParameters(std::cout,0,1);
        writeParameters(rpt::echo,0,1);
        PRINT2B1("--exiting...")
        ad_exit(-1);
    }
    
    if (mseOpModMode){
        dvariable mseCapF = mfexp(pMSE_LnC[1]);//multiplier on capture rate in directed fishery
        projectPopForTAC(mseCapF,0,cout);
        calcObjFunForTAC(dbg,cout);
    } else {
        if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
        calcObjFun(dbg,rpt::echo);
    }
    
    if ((!mseOpModMode)&&(ctrProcCallsInPhase==1)){
        //write objective function components only
        adstring fn = "tcsam02.ModelFits."+itoa(current_phase(),10)+"-"+str(ctrProcCallsInPhase)+"."+"R";
        ofstream os0(fn, ios::trunc);
        os0.precision(12);
        ReportToR_ModelFits(os0,-1.0,0,cout);
        os0.close();
        PRINT2B2("just after ReportToR_ModelFits. obj fun = ",objFun)
        PRINT2B1(" ")            
    }
    
    //assign values to likelihood profile numbers
    lkMMB = spB_yx(mxYr,MALE); 
    lkQM  = q_vyxms(1,mxYr+1,  MALE,MATURE,NEW_SHELL); 
    lkQF  = q_vyxms(1,mxYr+1,FEMALE,MATURE,NEW_SHELL); 
    lkNMIM = M_yxmsz(mxYr,  MALE,IMMATURE,NEW_SHELL)(1); 
    lkNMIF = M_yxmsz(mxYr,FEMALE,IMMATURE,NEW_SHELL)(1);
    lkNMMM = M_yxmsz(mxYr,  MALE,MATURE,NEW_SHELL)(1); 
    lkNMMF = M_yxmsz(mxYr,FEMALE,MATURE,NEW_SHELL)(1);
    
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
    
//    if ((current_phase()==1)&&(ctrProcCallsInPhase>450)){
//     //cout<<"writing parameters info to csv"<<endl;
//     adstring fn = "tcsam02.params.all."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".csv";
//     ofstream os1(fn, ios::trunc);
//     os1.precision(12);
//     os1<<"objFun = "<<cc<<objFun<<cc<<"phase ="<<cc<<current_phase()<<cc<<"proc calls in phase = "<<cc<<ctrProcCallsInPhase<<endl;
//     writeParameters(os1,0,0);//all parameters
//     os1.close();
//    }
//
    if (dbg>=dbgObjFun) {PRINT2B1("--END PROCEDURE_SECTION----------------")}
            
//-------------------------------------------------------------------------------------
FUNCTION void runAltPopDyMod(int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting runAltPopDyMod()"<<endl;
    //initialize population model
    initAltPopDyMod(debug, cout);
    
    //run population model
    for (int y=mnYr;y<=mxYr;y++){
        doSurveys(y,debug,cout);
       if (debug>=dbgPopDy) cout<<"--year = "<<y<<endl;        
         for (int x=1;x<=nSXs;x++){
            if (debug>=dbgPopDy) cout<<"----sex = "<<x<<endl;
            runAltPopDyModOneYear(y,x,debug,cout);  
            if (debug>=dbgPopDy) {cout<<"------R_z    = "<<endl; cout<<R_yxz(y,x)<<endl;}
        }
        if (debug>=dbgPopDy) cout<<endl;
    }
    doSurveys(mxYr+1,debug,cout);//do final surveys
    
    if (debug>=dbgPopDy) cout<<"finished runAltPopDyMod()"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION void initAltPopDyMod(int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting initAltPopDyMod()"<<endl;
    
    spB_yx.initialize();   //annual spawning biomass
    n_yxmsz.initialize();  //annual population abundance
    nmN_yxmsz.initialize();//number of crab killed by natural mortality
    tmN_yxmsz.initialize();//total number of crab killed 
    tmF_yxmsz.initialize();//total fishing mortality rate 
       
    setAllDevs(debug,cout);//set devs vectors
    
    calcRecruitment(debug,cout);//calculate recruitment
    calcNatMort(debug,cout);    //calculate natural mortality rates
    calcGrowth(debug,cout);     //calculate growth transition matrices
    calcPrM2M(debug,cout);      //calculate maturity ogives
    
    calcSelectivities(debug,cout); //calculate selectivity functions
    calcFisheryFs(debug,cout);     //calculate fishery F's
    calcSurveyQs(debug,cout);      //calculate survey Q's
    
    if (ptrMOs->optInitNatZ==0){
        //will build up population from recruitment (like TCSAM2013)
        //do nothing, because n_yxmsz has already been initialized to 0
    } else if (ptrMOs->optInitNatZ==1){
        //use equilibrium calculation to set initial n-at-z (like gmacs)
        //assumes no fishing occurs before model start
        calcEqNatZF100(initMnR,mnYr,debug,cout);//calculate n_xmsz
        n_yxmsz(mnYr) = n_xmsz;
    } else {
        cout<<"Unrecognized option for initial n-at-z: "<<ptrMOs->optInitNatZ<<endl;
        cout<<"Terminating!"<<endl;
        exit(-1);
    }

    if (debug>=dbgPopDy) cout<<"finished initAltPopDyMod()"<<endl;

//-------------------------------------------------------------------------------------
FUNCTION void runAltPopDyModOneYear(int y, int x, int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"Starting runAltPopDyModOneYear("<<y<<cc<<x<<")"<<endl;

   //1. Update population rates, based on year and sex
    pPDI->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(x);//assign weight-at-length
    pPDI->R_z   = R_yz(y);                   //relative recruitment-at-size
    pPDI->M_msz = M_yxmsz(y,x);              //rates of natural mortality
    pPDI->T_szz = prGr_yxszz(y,x);           //growth transition matrices
    for (int s=1;s<=nSCs;s++) pPDI->Th_sz(s) = prM2M_yxz(y,x);//pr(terminal molt|size)
    //if (debug) cout<<"updated pPDI for sex "<<x<<" in "<<y<<endl;
    
    //2. Update fishery conditions based on year and sex
    //TODO: move declarations to PARAMETER_SECTION
    dvar_vector hmF_f(1,nFsh);//handling mortality
    dvar4_array capF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//fishery capture rates
    dvar4_array retF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//fishery retention functions
    dvar4_array selF_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);//fishery selectivity functions
    hmF_f.initialize();
    capF_fmsz.initialize();
    retF_fmsz.initialize();
    selF_fmsz.initialize();
    for (int f=1;f<=nFsh;f++){
        hmF_f(f) = hmF_fy(f,y);
            for (int m=1;m<=nMSs;m++){
//                for (int s=1;s<=nSCs;s++) {
                            capF_fmsz(f,m) = cpF_fyxmsz(f,y,x,m);
                            retF_fmsz(f,m) = ret_fyxmsz(f,y,x,m);
                            selF_fmsz(f,m) = sel_fyxmsz(f,y,x,m);
//                }//s
            }//m
    }//f

    //update CatchInfo objects
    pCDI->setCaptureRates(capF_fmsz);
    pCDI->setRetentionFcns(retF_fmsz);
    pCDI->setSelectivityFcns(selF_fmsz);
    pCDI->setHandlingMortality(hmF_f);
    if (debug) cout<<"updated pCDI for sex "<<x<<" in "<<y<<endl;

    //3. run PopProjector
    pPPr->dtF = dtF_y(y); //time at which fishery occurs
    pPPr->dtM = dtM_y(y); //time at which mating occurs
    dvariable dirF = -1.0;//don't change scale on directed fishery F's
//    if (debug>=dbgPopDy) PopProjector::debug=1000;   
    n_yxmsz(y+1,x)       = pPPr->project(dirF,n_yxmsz(y,x),cout);
//    if (debug>=dbgPopDy) PopProjector::debug=0;
    n_yxmsz(y+1,x,IMMATURE,NEW_SHELL) += R_yxz(y,x);
    spB_yx(y,x)          = pPPr->getMatureBiomassAtMating();
    dvar4_array cpN_fmsz = pPPr->getFisheriesCaptureAbundance();
    dvar4_array rmN_fmsz = pPPr->getRetainedCatchMortality();
    dvar4_array dmN_fmsz = pPPr->getDiscardCatchMortality();
    for (int f=1;f<=nFsh;f++){
        for (int m=1;m<=nMSs;m++){
            cpN_fyxmsz(f,y,x,m) = cpN_fmsz(f,m);                          //capture abundance
            rmN_fyxmsz(f,y,x,m) = rmN_fmsz(f,m);                          //retained mortality abundance
            dmN_fyxmsz(f,y,x,m) = dmN_fmsz(f,m);                          //discard mortality abundance
            dsN_fyxmsz(f,y,x,m) = cpN_fyxmsz(f,y,x,m)-rmN_fyxmsz(f,y,x,m);//total discard abundance
        }//m
    }//f
    if (debug) cout<<"ran pPPr for sex "<<x<<" in "<<y<<endl;
    
    if (debug>=dbgPopDy) cout<<"finished runAltPopDyModOneYear("<<y<<cc<<x<<")"<<endl;
    
///**
// * MSE OpModMode ONLY:
// * Project the population from July 1, year to July 1, year+1 based on
// * directed fishing for males at a rate that will yield the TAC when the
// * model has converged.
// * 
// * This function is directly based on runAltPopDyModForOneYear.
// * 
// */
FUNCTION void projectPopForTAC(dvariable& mseCapF, int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"projectPopForTAC()"<<endl;

    dvar4_array cpN_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array rmN_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array dmN_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
        
    prjRetCatchMortBio_fx.initialize();
    prjDscCatchMortBio_fx.initialize();
    prjTotCatchMortBio_fx.initialize();
    prj_cpN_fxmsz.initialize(); //capture abundance
    prj_rmN_fxmsz.initialize(); //retained mortality abundance
    prj_dmN_fxmsz.initialize(); //discard mortality abundance
    prj_dsN_fxmsz.initialize(); //total discard abundance
    
    for (int x=1;x<=nSXs;x++){   
        //project population under directed fishing male capture rate mseCapF
        dvariable dirF = 1.0*mseCapF;
//        if (debug>=dbgPopDy) PopProjector::debug=1000;
        PopProjector* pPPr = ppOpModPPr_x[x-1];
        dvar3_array n_msz  = 1.0*ptrOMI->n_xmsz(x);//sex "x" population abundance at start of year
        prj_n_xmsz(x)      = 1.0*pPPr->project(dirF,n_msz,cout);
//        if (debug>=dbgPopDy) PopProjector::debug=0;
        prj_spB_x(x)       = 1.0*pPPr->getMatureBiomassAtMating();
        prj_n_xmsz(x,IMMATURE,NEW_SHELL) += prjR*ptrOMI->R_x(x)*ptrOMI->R_z;
        
        //extract catch info
        cpN_fmsz.initialize();
        rmN_fmsz.initialize();
        dmN_fmsz.initialize();
        cpN_fmsz = pPPr->getFisheriesCaptureAbundance();
        rmN_fmsz = pPPr->getRetainedCatchMortality();
        dmN_fmsz = pPPr->getDiscardCatchMortality();
//        if (debug){
//            cout<<"#--check on fishery rates used"<<endl;
//            for (int f=1;f<=nFsh;f++){
//               for (int m=1;m<=nMSs;m++){
//                    for (int s=1;s<=nSCs;s++){
//                        cout<<"#----x,f,m,s = "<<x<<cc<<f<<cc<<m<<cc<<s<<endl;
//                        cout<<"diff(cpF_fxmsz) = "<<pPPr->pCI->cpF_fmsz(f,m,s)-ptrOMI->cpF_fxmsz(f,x,m,s)<<endl;
//                        cout<<"diff(rmF_fxmsz) = "<<pPPr->pCI->retF_fmsz(f,m,s)-ptrOMI->ret_fxmsz(f,x,m,s)<<endl;
//                    }
//                }
//            }
//        }
        for (int f=1;f<=nFsh;f++){
            for (int m=1;m<=nMSs;m++){
                prj_cpN_fxmsz(f,x,m) = cpN_fmsz(f,m); //capture abundance
                prj_rmN_fxmsz(f,x,m) = rmN_fmsz(f,m); //retained mortality abundance
                prj_dmN_fxmsz(f,x,m) = dmN_fmsz(f,m); //discard mortality abundance
                //total discard abundance
                prj_dsN_fxmsz(f,x,m) = prj_cpN_fxmsz(f,x,m)-prj_rmN_fxmsz(f,x,m);
                for (int s=1;s<=nSCs;s++){
                    //retained catch mortality (biomass)
                    prjRetCatchMortBio_fx(f,x) += prj_rmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                    //discarded catch mortality (biomass)
                    prjDscCatchMortBio_fx(f,x) += prj_dmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                }//s
            }//m
            prjTotCatchMortBio_fx(f,x) = prjRetCatchMortBio_fx(f,x) + prjDscCatchMortBio_fx(f,x);
        }//f
        if (debug) cout<<"ran pPPr for sex "<<x<<endl;
    }//x
    
    if (debug>=dbgPopDy){
        cout<<"#--projected population with directed fishery F = "<<mseCapF<<endl;
        cout<<"#----total catch mortality (biomass) = "<<sum(prjTotCatchMortBio_fx)<<endl;
        cout<<"#----retained catch mortality (biomass): "<<endl<<prjRetCatchMortBio_fx<<endl;
        cout<<"#----discards catch mortality (biomass): "<<endl<<prjDscCatchMortBio_fx<<endl;
    }
    if (debug>=dbgPopDy) cout<<"finished projectPopForTAC()"<<endl;
    
//-------------------------------------------------------------------------------------
//Calculate objective function TODO: finish
FUNCTION void calcObjFunForTAC(int debug, ostream& cout)
    if ((debug>=dbgObjFun)||(debug<0)) {PRINT2B1("----Starting calcObjFunForTAC()")}

    int k = 0;
    dvector objFunV(0,20);

    //reset objective function
    objFun.initialize(); objFunV[k++] = value(objFun);
    
    //Fit OpMod catch to TAC
    dvariable prdTAC = 0.0;
    dvariable prdTot = 0.0;
    for (int f=1;f<=nFsh;f++){
        prdTAC += sum(prjRetCatchMortBio_fx(f));
        prdTot += sum(prjTotCatchMortBio_fx(f));
    }
    
    dvariable nllForTAC = square(inpTAC-prdTAC);
    
    dvariable nllForOFL = 0.0; 
    posfun(inpOFL-prdTot,1.0e-2,nllForOFL);
    nllForOFL *= 1000.0;
    
    objFun += (nllForTAC + nllForOFL);   objFunV[k++] = value(objFun);
    
    if ((debug>=dbgObjFun)||(debug<0)){
        int k=0;
        PRINT2B2("proc call          = ",ctrProcCalls)
        PRINT2B2("proc call in phase = ",ctrProcCallsInPhase)
        PRINT2B2("pMSE_LnC =",pMSE_LnC[1]);
        PRINT2B2("mseCapF  =",mfexp(pMSE_LnC[1]));
        PRINT2B2("after initialization. objFun =",objFunV[k++])
        PRINT2B2("inpTAC    = ",inpTAC)    
        PRINT2B2("prdTAC    = ",prdTAC)    
        PRINT2B2("nllForTAC = ",nllForTAC)    
        PRINT2B2("inpOFL    = ",inpOFL)    
        PRINT2B2("prdTot    = ",prdTot)    
        for (int f=1;f<=nFsh;f++) {
            adstring s = str(f)+" = ";
            for (int x=1;x<=nSXs;x++) s = s + str(value(prjTotCatchMortBio_fx(f,x)))+"\t";
            PRINT2B1(s);
        }
        PRINT2B2("nllForOFL = ",nllForOFL)    
        PRINT2B2("after MSE OpMod call objFun =",value(objFun))
        PRINT2B1("----Finished calcObjFunForTAC()")
    }
    
///**
// * MSE OpModMode ONLY:
// * Project the population from July 1, year to July 1, year+1 given that
// * the TAC is 0 and the directed fishery is closed.
// * 
// * This function is directly based on projectPopForTAC.
// * 
// */
FUNCTION void projectPopForZeroTAC(int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"projectPopForZeroTAC()"<<endl;

    dvar4_array cpN_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array rmN_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array dmN_fmsz(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
        
    prjRetCatchMortBio_fx.initialize();
    prjDscCatchMortBio_fx.initialize();
    prjTotCatchMortBio_fx.initialize();
    prj_cpN_fxmsz.initialize(); //capture abundance
    prj_rmN_fxmsz.initialize(); //retained mortality abundance
    prj_dmN_fxmsz.initialize(); //discard mortality abundance
    prj_dsN_fxmsz.initialize(); //total discard abundance
    
    for (int x=1;x<=nSXs;x++){   
        //project population under no directed fishing
        dvariable dirF = 0.0;
//        if (debug>=dbgPopDy) PopProjector::debug=1000;
        PopProjector* pPPr = ppOpModPPr_x[x-1];
        dvar3_array n_msz  = 1.0*ptrOMI->n_xmsz(x);//sex "x" population abundance at start of year
        prj_n_xmsz(x)      = 1.0*pPPr->project(dirF,n_msz,cout);
//        if (debug>=dbgPopDy) PopProjector::debug=0;
        prj_spB_x(x)       = 1.0*pPPr->getMatureBiomassAtMating();
        prj_n_xmsz(x,IMMATURE,NEW_SHELL) += prjR*ptrOMI->R_x(x)*ptrOMI->R_z;
        
        //extract catch info
        cpN_fmsz.initialize();
        rmN_fmsz.initialize();
        dmN_fmsz.initialize();
        cpN_fmsz = pPPr->getFisheriesCaptureAbundance();
        rmN_fmsz = pPPr->getRetainedCatchMortality();
        dmN_fmsz = pPPr->getDiscardCatchMortality();
//        if (debug){
//            cout<<"#--check on fishery rates used"<<endl;
//            for (int f=1;f<=nFsh;f++){
//               for (int m=1;m<=nMSs;m++){
//                    for (int s=1;s<=nSCs;s++){
//                        cout<<"#----x,f,m,s = "<<x<<cc<<f<<cc<<m<<cc<<s<<endl;
//                        cout<<"diff(cpF_fxmsz) = "<<pPPr->pCI->cpF_fmsz(f,m,s)-ptrOMI->cpF_fxmsz(f,x,m,s)<<endl;
//                        cout<<"diff(rmF_fxmsz) = "<<pPPr->pCI->retF_fmsz(f,m,s)-ptrOMI->ret_fxmsz(f,x,m,s)<<endl;
//                    }
//                }
//            }
//        }
        for (int f=1;f<=nFsh;f++){
            for (int m=1;m<=nMSs;m++){
                prj_cpN_fxmsz(f,x,m) = cpN_fmsz(f,m); //capture abundance
                prj_rmN_fxmsz(f,x,m) = rmN_fmsz(f,m); //retained mortality abundance
                prj_dmN_fxmsz(f,x,m) = dmN_fmsz(f,m); //discard mortality abundance
                //total discard abundance
                prj_dsN_fxmsz(f,x,m) = prj_cpN_fxmsz(f,x,m)-prj_rmN_fxmsz(f,x,m);
                for (int s=1;s<=nSCs;s++){
                    //retained catch mortality (biomass)
                    prjRetCatchMortBio_fx(f,x) += prj_rmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                    //discarded catch mortality (biomass)
                    prjDscCatchMortBio_fx(f,x) += prj_dmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
                }//s
            }//m
            prjTotCatchMortBio_fx(f,x) = prjRetCatchMortBio_fx(f,x) + prjDscCatchMortBio_fx(f,x);
        }//f
        if (debug) cout<<"ran pPPr for sex "<<x<<endl;
    }//x
    
    if (debug>=dbgPopDy){
        cout<<"#--projected population directed fishery closed (TAC=0)"<<endl;
        cout<<"#----total catch mortality (biomass) = "<<sum(prjDscCatchMortBio_fx)+sum(prjDscCatchMortBio_fx)<<endl;
        cout<<"#----retained catch mortality (biomass): "<<endl<<prjRetCatchMortBio_fx<<endl;
        cout<<"#----discards catch mortality (biomass): "<<endl<<prjDscCatchMortBio_fx<<endl;
    }
    if (debug>=dbgPopDy) cout<<"finished projectPopForZeroTAC()"<<endl;
    
/////**
//// * MSE OpModMode ONLY:
//// * Project the population from July 1, year to July 1, year+1 given that
//// * the TAC is 0 and the directed fishery is closed.
//// * 
//// * This function is directly based on projectPopForTAC.
//// * 
//// */
//FUNCTION void projectPopForZeroTAC(int debug, ostream& cout)
//    if (debug>=dbgPopDy) cout<<"projectPopForZeroTAC()"<<endl;
//
//    prjRetCatchMortBio_fx.initialize();
//    prjDscCatchMortBio_fx.initialize();
//    prjTotCatchMortBio_fx.initialize();
//    prj_cpN_fxmsz.initialize(); //capture abundance
//    prj_rmN_fxmsz.initialize(); //retained mortality abundance
//    prj_dmN_fxmsz.initialize(); //discard mortality abundance
//    prj_dsN_fxmsz.initialize(); //total discard abundance
//    
//    dvariable maxCapF = 0.0;
//    
//    for (int x=1;x<=nSXs;x++){   
//        //1. Update population rates, based on year and sex
//        pPDI->w_mz  = ptrOMI->wAtZ_xmz(x);   //assign weight-at-length
//        pPDI->R_z   = ptrOMI->R_z;           //relative recruitment-at-size
//        pPDI->M_msz = ptrOMI->M_xmsz(x);     //rates of natural mortality
//        pPDI->T_szz = ptrOMI->prGr_xszz(x);  //growth transition matrices
//        //pr(terminal molt|size)
//        for (int s=1;s<=nSCs;s++) pPDI->Th_sz(s) = ptrOMI->prM2M_xz(x);
//        //if (debug) cout<<"updated pPDI for sex "<<x<<" in "<<y<<endl;
//
//        //2. Update fishery conditions based on year and sex
//        prj_hmF_f = ptrOMI->hmF_f;
//        prj_capF_fmsz.initialize();
//        prj_retF_fmsz.initialize();
//        prj_selF_fmsz.initialize();
//        for (int f=1;f<=nFsh;f++){
//            for (int m=1;m<=nMSs;m++){
//                prj_capF_fmsz(f,m) = ptrOMI->cpF_fxmsz(f,x,m);
//                prj_retF_fmsz(f,m) = ptrOMI->ret_fxmsz(f,x,m);
//                prj_selF_fmsz(f,m) = ptrOMI->sel_fxmsz(f,x,m);
//            }//m
//        }//f
//
//        
//        //directed fishery is closed,so set capture rates to 0
//        prj_capF_fmsz(1) = 0.0*prj_selF_fmsz(1);
//
//        //update CatchInfo objects
//        pCDI->setHandlingMortality(prj_hmF_f);
//        pCDI->setCaptureRates(prj_capF_fmsz);
//        pCDI->setRetentionFcns(prj_retF_fmsz);
//        pCDI->setSelectivityFcns(prj_selF_fmsz);
//        if (debug) cout<<"updated pCDI for sex "<<x<<endl;
//
//        //3. run PopProjector
//        pPPr->dtF = ptrOMI->dtF; //time at which fishery occurs
//        pPPr->dtM = ptrOMI->dtM; //time at which mating occurs
//        dvariable dirF = -1.0;//don't change scale on directed fishery F's
////        if (debug>=dbgPopDy) PopProjector::debug=1000;
//        dvar3_array n_msz   = ptrOMI->n_xmsz(x);//sex "x" population abundance at start of year
//        prj_n_xmsz(x)       = pPPr->project(dirF,n_msz,cout);
////        if (debug>=dbgPopDy) PopProjector::debug=0;
//        prj_n_xmsz(x,IMMATURE,NEW_SHELL) += prjR*ptrOMI->R_x(x)*ptrOMI->R_z;
//        prj_spB_x(x)         = pPPr->getMatureBiomassAtMating();
//        dvar4_array cpN_fmsz = pPPr->getFisheriesCaptureAbundance();
//        dvar4_array rmN_fmsz = pPPr->getRetainedCatchMortality();
//        dvar4_array dmN_fmsz = pPPr->getDiscardCatchMortality();
//        for (int f=1;f<=nFsh;f++){
//            for (int m=1;m<=nMSs;m++){
//                prj_cpN_fxmsz(f,x,m) = cpN_fmsz(f,m); //capture abundance
//                prj_rmN_fxmsz(f,x,m) = rmN_fmsz(f,m); //retained mortality abundance
//                prj_dmN_fxmsz(f,x,m) = dmN_fmsz(f,m); //discard mortality abundance
//                //total discard abundance
//                prj_dsN_fxmsz(f,x,m) = prj_cpN_fxmsz(f,x,m)-prj_rmN_fxmsz(f,x,m);
//                for (int s=1;s<=nSCs;s++){
//                    //retained catch mortality (biomass)
//                    prjRetCatchMortBio_fx(f,x) += prj_rmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
//                    //discarded catch mortality (biomass)
//                    prjDscCatchMortBio_fx(f,x) += prj_dmN_fxmsz(f,x,m,s)*ptrOMI->wAtZ_xmz(x,m);
//                }//s
//            }//m
//            prjTotCatchMortBio_fx(f,x) = prjRetCatchMortBio_fx(f,x) + prjDscCatchMortBio_fx(f,x);
//        }//f
//        if (debug) cout<<"ran pPPr for sex "<<x<<endl;
//    }//x
//    
//    if (debug>=dbgPopDy) cout<<"finished projectPopForZeroTAC()"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION void runPopDyMod(int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting runPopDyMod()"<<endl;
    //initialize population model
    initPopDyMod(0, cout);
    
    //run population model
    for (int y=mnYr;y<=mxYr;y++){
        doSurveys(y,0,cout);
        runPopDyModOneYear(y,debug,cout);        
    }
    doSurveys(mxYr+1,0,cout);//do final surveys
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyMod()"<<endl;
    
//-------------------------------------------------------------------------------------
FUNCTION void initPopDyMod(int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting initPopDyMod()"<<endl;
    
    spB_yx.initialize();
    n_yxmsz.initialize();
    nmN_yxmsz.initialize();
    tmN_yxmsz.initialize();
    tmF_yxmsz.initialize();
       
    setAllDevs(debug,cout);//set devs vectors
    
    calcRecruitment(debug,cout);//calculate recruitment
    calcNatMort(debug,cout);    //calculate natural mortality rates
    calcGrowth(debug,cout);     //calculate growth transition matrices
    calcPrM2M(debug,cout);      //calculate maturity ogives
    
    calcSelectivities(debug,cout); //calculate selectivity functions
    calcFisheryFs(debug,cout);     //calculate fishery F's
    calcSurveyQs(debug,cout);      //calculate survey Q's
    
    if (ptrMOs->optInitNatZ==0){
        //will build up population from recruitment (like TCSAM2013)
        //do nothing, because n_yxmsz has already been initialized to 0
    } else if (ptrMOs->optInitNatZ==1){
        //use equilibrium calculation to set initial n-at-z (like gmacs)
        //assumes no fishing occurs before model start
        calcEqNatZF100(initMnR,mnYr,debug,cout);//calculate n_xmsz
        n_yxmsz(mnYr) = n_xmsz;
    } else {
        cout<<"Unrecognized option for initial n-at-z: "<<ptrMOs->optInitNatZ<<endl;
        cout<<"Terminating!"<<endl;
        exit(-1);
    }

    if (debug>=dbgPopDy) cout<<"finished initPopDyMod()"<<endl;

//-------------------------------------------------------------------------------------
FUNCTION void runPopDyModOneYear(int y, int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"Starting runPopDyModOneYear("<<y<<")"<<endl;

    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n2_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n3_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n4_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    dvar4_array n5_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    
    n1_xmsz.initialize();
    n2_xmsz.initialize();
    n3_xmsz.initialize();
    n4_xmsz.initialize();
    n5_xmsz.initialize();
    
    if (dtF_y(y)<=dtM_y(y)){//fishery occurs BEFORE molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs BEFORE molting/growth/maturity"<<endl;
        //apply natural mortality before fisheries
        n1_xmsz = applyNatMort(n_yxmsz(y),y,dtF_y(y),debug,cout);
        //conduct fisheries
        n2_xmsz = applyFshMort(n1_xmsz,y,debug,cout);
        //apply natural mortality from fisheries to molting/growth/maturity
        if (dtF_y(y)==dtM_y(y)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,y,dtM_y(y)-dtF_y(y),debug,cout);
        }
        //calc mature (spawning) biomass at time of mating, but BEFORE growth/maturity (TODO: does this make sense??)
        spB_yx(y) = calcSpB(n3_xmsz,y,debug,cout);
        //apply molting, growth and maturation
        n4_xmsz = applyMGM(n3_xmsz,y,debug,cout);
        //apply natural mortality to end of year
        if (dtM_y(y)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,y,1.0-dtM_y(y),debug,cout);
        }
    } else {              //fishery occurs AFTER molting/growth/maturity
        if (debug>=dbgPopDy) cout<<"Fishery occurs AFTER molting/growth/maturity"<<endl;
        //apply natural mortality before molting/growth/maturity
        n1_xmsz = applyNatMort(n_yxmsz(y),y,dtM_y(y),debug,cout);
        //calc mature (spawning) biomass at time of mating (TODO: does this make sense??)
        spB_yx(y) = calcSpB(n1_xmsz,y,debug,cout);
        //apply molting, growth and maturation
        n2_xmsz = applyMGM(n1_xmsz,y,debug,cout);
        //apply natural mortality from molting/growth/maturity to fisheries
        if (dtM_y(y)==dtF_y(y)) {
            n3_xmsz = n2_xmsz;
        } else {
            n3_xmsz = applyNatMort(n2_xmsz,y,dtF_y(y)-dtM_y(y),debug,cout);
        }
        //conduct fisheries
        n4_xmsz = applyFshMort(n3_xmsz,y,debug,cout);
        //apply natural mortality to end of year
        if (dtF_y(y)==1.0) {
            n5_xmsz = n4_xmsz;
        } else {
            n5_xmsz = applyNatMort(n4_xmsz,y,1.0-dtF_y(y),debug,cout);
        }
    }
    
    //advance surviving individuals to next year
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                n_yxmsz(y+1,x,m,s) = n5_xmsz(x,m,s);
            }
        }
    }
    
    //add in recruits (NOTE: R_y(y) here corresponds to R_y(y+1) in TCSAM2013)
    for (int x=1;x<=nSXs;x++) n_yxmsz(y+1,x,IMMATURE,NEW_SHELL) += R_yxz(y,x);
    
    if (debug>=dbgPopDy){
        cout<<"----year = "<<y<<endl;
        for (int x=1;x<=nSXs;x++){
            cout<<"----sex = "<<x<<endl;
            cout<<"------n0_msz = "<<endl; wts::print(n_yxmsz(y,x),cout,1);
            cout<<"------n1_msz = "<<endl; wts::print(n1_xmsz(x),cout,1);
            for (int f=1;f<=nFsh;f++){
                cout<<"------fishery = "<<f<<endl;
                cout<<"--------hmF_f = "<<hmF_fy(f,y)<<endl;
                cout<<"--------cpF_msz = "<<endl; wts::print(cpF_fyxmsz(f,y,x),cout,1);
                cout<<"--------rmF_msz = "<<endl; wts::print(rmF_fyxmsz(f,y,x),cout,1);
                cout<<"--------dmF_msz = "<<endl; wts::print(dmF_fyxmsz(f,y,x),cout,1);
//                cout<<"--------cpN_msz = "<<endl; wts::print(cpN_fyxmsz(f,y,x),cout,1);
//                cout<<"--------rmN_msz = "<<endl; wts::print(rmN_fyxmsz(f,y,x),cout,1);
//                cout<<"--------dmN_msz = "<<endl; wts::print(dmN_fyxmsz(f,y,x),cout,1);
            }//f
            cout<<"------n2_msz = "<<endl; wts::print(n2_xmsz(x),cout,1);
            cout<<"------n3_msz = "<<endl; wts::print(n3_xmsz(x),cout,1);
            cout<<"------n4_msz = "<<endl; wts::print(n4_xmsz(x),cout,1);
            cout<<"------n5_msz = "<<endl; wts::print(n5_xmsz(x),cout,1);
            cout<<"------R_z    = "<<endl; cout<<R_yxz(y,x)<<endl;
        }
        cout<<endl;
    }
    
    if (debug>=dbgPopDy) cout<<"finished runPopDyModOneYear("<<y<<")"<<endl;
    
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
    dvector     tdF_z(1,nZBs);//for use in calculating fishing rate components
    dvar_vector tm_z(1,nZBs); //total mortality (numbers) by size
    dvar_vector tvF_z(1,nZBs);//for use in calculating fishing rate components
    dvar_vector tfF_z(1,nZBs);//for use in calculating fishing rate components
    dvar4_array n1_xmsz(1,nSXs,1,nMSs,1,nSCs,1,nZBs);//numbers surviving fisheries
    n1_xmsz.initialize();
    for (int x=1;x<=nSXs;x++){
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nSCs;s++){
                //tmF_yxmsz(y,x,m,s).initialize();//total fishing mortality rate
                for (int f=1;f<=nFsh;f++) tmF_yxmsz(y,x,m,s) += rmF_fyxmsz(f,y,x,m,s)+dmF_fyxmsz(f,y,x,m,s);
                n1_xmsz(x,m,s) = elem_prod(mfexp(-tmF_yxmsz(y,x,m,s)),n0_xmsz(x,m,s));//numbers surviving all fisheries
                tm_z = n0_xmsz(x,m,s)-n1_xmsz(x,m,s);                                 //numbers killed by all fisheries
                tmN_yxmsz(y,x,m,s) += tm_z;            //add in numbers killed by all fisheries to total killed
                
                //calculate fishing rate components (need to ensure NOT dividing by 0)
                tdF_z = value(tmF_yxmsz(y,x,m,s));
                tvF_z = elem_prod(1-wts::isEQ(tdF_z,0.0),tmF_yxmsz(y,x,m,s)) + wts::isEQ(tdF_z,0.0);//= tmF_yxmsz(y,x,m,s,z) if tmF_yxmsz(y,x,m,s,z) > 0, else = 1
                tfF_z = elem_div(1.0-mfexp(-tmF_yxmsz(y,x,m,s)),tvF_z);//= (1-exp(-tmF_yxmsz(y,x,m,s,z)))/tmF_yxmsz(y,x,m,s,z) if tmF_yxmsz(y,x,m,s,z) > 0, else = 1
//                cout<<"y,x,m,s = "<<y<<tb<<x<<tb<<m<<tb<<s<<endl;
//                cout<<"tdF_z      = "<<tdF_z<<endl;
//                cout<<"tdF_z==0.0 = "<<wts::isEQ(tdF_z,0.0)<<endl;
//                cout<<"tvF_z      = "<<tvF_z<<endl;
//                int tmp; cout<<"Enter 1 to continue > "; cin>>tmp; if (!tmp) exit(-1);
                for (int f=1;f<=nFsh;f++){                   
                    cpN_fyxmsz(f,y,x,m,s) = elem_prod(elem_prod(cpF_fyxmsz(f,y,x,m,s),tfF_z),n0_xmsz(x,m,s)); //numbers captured in fishery f
                    rmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_prod(rmF_fyxmsz(f,y,x,m,s),tfF_z),n0_xmsz(x,m,s)); //retained mortality in fishery f (numbers)
                    dmN_fyxmsz(f,y,x,m,s) = elem_prod(elem_prod(dmF_fyxmsz(f,y,x,m,s),tfF_z),n0_xmsz(x,m,s)); //discards mortality in fishery f (numbers)
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
        dvar_vector np_z = prGr_yxszz(y,x,NEW_SHELL)*n0_xmsz(x,IMMATURE,NEW_SHELL);
        n1_xmsz(x,IMMATURE,NEW_SHELL) = elem_prod(1.0-prM2M_yxz(y,x),np_z);
        n1_xmsz(x,IMMATURE,OLD_SHELL) = 0.0;
        n1_xmsz(x,MATURE,NEW_SHELL)   = elem_prod(    prM2M_yxz(y,x),np_z);
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
    nDevsLnR_c.initialize();
    stdvDevsLnR_c.initialize();
    devsLnR_cy.initialize();
    zscrDevsLnR_cy.initialize();
    
    int k; int y;
    dvector dzs = zBs+(zBs[2]-zBs[1])/2.0-zBs[1];
    dvar_vector ptLnR  = ptrRI->pLnR->calcArithScaleVals(pLnR);
    dvar_vector ptRCV  = ptrRI->pRCV->calcArithScaleVals(pRCV);
    dvar_vector ptRX   = ptrRI->pRX->calcArithScaleVals(pRX);
    dvar_vector ptRa   = ptrRI->pRa->calcArithScaleVals(pRa);
    dvar_vector ptRb   = ptrRI->pRb->calcArithScaleVals(pRb);
    for (int pc=1;pc<=ptrRI->nPCs;pc++){
        ivector pids = ptrRI->getPCIDs(pc);
        k=ptrRI->nIVs+1;//first parameter variable column in ParameterComnbinations
        dvariable mnLnR = ptLnR(pids[k++]);
        dvariable cvR   = ptRCV(pids[k++]);
        dvariable xR    = ptRX(pids[k++]);
        dvariable aR    = ptRa(pids[k++]);
        dvariable bR    = ptRb(pids[k++]);
        if (debug>dbgCalcProcs){
            cout<<"pids  = "<<pids<<endl;
            cout<<"mnLnR = "<<mnLnR<<endl;
            cout<<"cvR   = "<<cvR<<endl;
            cout<<"xR    = "<<xR<<endl;
            cout<<"aR    = "<<aR<<endl;
            cout<<"bR    = "<<bR<<endl;
        }

        int useDevs = pids[k]; k++;
        dvariable mnR;   //mean recruitment
        dvariable varLnR;//ln-scale variance in recruitment
        dvar_vector dvsLnR;
        ivector idxDevsLnR;
        varLnR = log(1.0+square(cvR));    //ln-scale variance
        stdvDevsLnR_c(pc) = sqrt(varLnR); //ln-scale std dev
        mnR    = mfexp(mnLnR+varLnR/2.0); //mean recruitment
        if (useDevs) {
            dvsLnR     = devsLnR(useDevs);
            idxDevsLnR = idxsDevsLnR(useDevs);
            if (debug>dbgCalcProcs) {
                cout<<"lims(dvsLnR) = "<<dvsLnR.indexmin()<<cc<<dvsLnR.indexmax()<<endl;
                cout<<"idx(dvsLnR) = "<<idxDevsLnR<<endl;
                cout<<"dvsLnR = "<<dvsLnR<<endl;
            }
        }
        
        Rx_c(pc) = xR;
        R_cz(pc) = elem_prod(pow(dzs,(aR/bR)-1.0),mfexp(-dzs/bR));
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

                nDevsLnR_c(pc)++;//increment count of devs for this pc
                devsLnR_cy(pc,y) = dvsLnR(idxDevsLnR(y));//rec devs by pc and year
            } else {
                if (debug>dbgCalcProcs) cout<<"skipping y,i = "<<y<<tb<<idxDevsLnR(y)<<endl;
            }
        }//idx
        zscrDevsLnR_cy(pc) = devsLnR_cy(pc)/stdvDevsLnR_c(pc);//standardized zscores (assuming mean=0)
    }//pc
    
    if (debug>dbgCalcProcs) {
        cout<<"R_y = "<<R_y<<endl;
        cout<<"R_yx(MALE) = "<<column(R_yx,MALE)<<endl;
        cout<<"R_yz  = "<<endl<<R_yz<<endl;
        cout<<"nDevs = "<<nDevsLnR_c<<endl;
        cout<<"sigR  = "<<stdvDevsLnR_c<<endl;
        cout<<"devs  = "<<devsLnR_cy<<endl;
        cout<<"zscr  = "<<zscrDevsLnR_cy<<endl;
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
    
    dvariable lnM;
    dvar_vector M_z(1,nZBs);
    
    M_c.initialize();
    M_yxmsz.initialize();

    int y; int mnx; int mxx; int mnm; int mxm; int mns; int mxs;
    dvar_vector ptM  = ptrNM->pM->calcArithScaleVals(pM);
    dvar_vector ptDM1  = ptrNM->pDM1->calcArithScaleVals(pDM1);
    dvar_vector ptDM2  = ptrNM->pDM2->calcArithScaleVals(pDM2);
    dvar_vector ptDM3  = ptrNM->pDM3->calcArithScaleVals(pDM3);
    dvar_vector ptDM4  = ptrNM->pDM4->calcArithScaleVals(pDM4);
    for (int pc=1;pc<=ptrNM->nPCs;pc++){
        if (debug>dbgCalcProcs) cout<<"pc = "<<pc<<endl;
        lnM.initialize();
        ivector pids = ptrNM->getPCIDs(pc);
        int k=ptrNM->nIVs+1;//1st parameter variable column
        if (debug>dbgCalcProcs) cout<<"pids = "<<pids(k,pids.indexmax())<<endl;
        if (ptrMOs->optParamNM==0){
            //add in base (arithmetic-scale) natural mortality
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptM["<<pids[k]<<"]: "<<ptM(pids[k])<<endl;
                lnM += log(ptM(pids[k]));
            } k++;
            //add in ln-scale offset 1
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM1["<<pids[k]<<"]: "<<ptDM1(pids[k])<<endl;
                lnM += ptDM1(pids[k]);
            } k++;
            //add in ln-scale offset 2
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM2["<<pids[k]<<"]: "<<ptDM2(pids[k])<<endl;
                lnM += ptDM2(pids[k]);
            } k++;
            //add in ln-scale offset 3
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM3["<<pids[k]<<"]: "<<ptDM3(pids[k])<<endl;
                lnM += ptDM3(pids[k]);
            } k++;
            //add in ln-scale offset 4
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptDM4["<<pids[k]<<"]: "<<ptDM4(pids[k])<<endl;
                lnM += ptDM4(pids[k]);
            }  k++; //advance k to zScaling in pids
        } else if (ptrMOs->optParamNM==1){
            //TCSAM2013 parameterization: arithmetic scale multipliers to base
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Adding ptM["<<pids[k]<<"]: "<<ptM(pids[k])<<endl;
                lnM += log(ptM(pids[k]));
            } k++;
            //multiply offset 1 (for immature crab)
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM1["<<pids[k]<<"]: "<<ptDM1(pids[k])<<endl;
                lnM += log(ptDM1(pids[k]));//add on ln-scale
            } k++;
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM2["<<pids[k]<<"]: "<<ptDM2(pids[k])<<endl;
                lnM += log(ptDM2(pids[k]));
            } k++;
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM3["<<pids[k]<<"]: "<<ptDM3(pids[k])<<endl;
                lnM += log(ptDM3(pids[k]));
            } k++;
            //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
            if (pids[k]) {
                if (debug>dbgCalcProcs) cout<<"Multiplying by offset ptDM4["<<pids[k]<<"]: "<<ptDM4(pids[k])<<endl;
                lnM += log(ptDM4(pids[k]));
            } k++;
        }//optParamNM
        
        //convert from ln-scale to arithmetic scale
        M_c(pc) = mfexp(lnM);
        if (debug>dbgCalcProcs) cout<<"lnM= "<<lnM<<tb<<"M_c = "<<M_c(pc)<<endl;
        
        //add in size-scaling, if requested
        M_z.initialize();
        //cout<<"k = "<<k<<tb<<"pids[k] = "<<pids[k]<<endl;
        if (pids[k]&&(current_phase()>=pids[k])) {
            if (debug>dbgCalcProcs) cout<<"adding size scaling"<<endl;
            M_z = M_c(pc)*(zMref/zBs);//factor in size dependence
        } else {
            if (debug>dbgCalcProcs) cout<<"not adding size scaling"<<endl;
            M_z = M_c(pc);//no size dependence
        }
        if (debug>dbgCalcProcs) cout<<"finished scaling"<<endl;
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrNM->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);//year
            if ((mnYr<=y)&&(y<=mxYr)){
                mnx = mxx = idxs(idx,2);//sex index
                if (mnx==tcsam::ALL_SXs){mnx = 1; mxx = tcsam::nSXs;}
                mnm = mxm = idxs(idx,3);//maturity state
                if (mnm==tcsam::ALL_MSs){mnm = 1; mxm = tcsam::nMSs;}
                mns = mxs = idxs(idx,4);//shell condition
                if (mns==tcsam::ALL_SCs){mns = 1; mxs = tcsam::nSCs;}
                for (int x=mnx;x<=mxx;x++){
                    for (int m=mnm;m<=mxm;m++){
                        for (int s=mns;s<=mxs;s++) M_yxmsz(y,x,m,s) = 1.0*M_z;
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
FUNCTION void calcPrM2M(int debug, ostream& cout)
    if (debug>dbgCalcProcs) cout<<"starting calcPrM2M()"<<endl;

    Molt2MaturityInfo* ptrM2MI = ptrMPI->ptrM2M;
    
    prM2M_cz.initialize();
    prM2M_yxz.initialize();
    
    int k; int y; int mnx; int mxx;
    for (int pc=1;pc<=ptrM2MI->nPCs;pc++){
        ivector pids = ptrM2MI->getPCIDs(pc);
        k=ptrM2MI->nIVs+1;//first parameter variable column in ParameterComnbinations
        BoundedVectorInfo* pBVI = (*ptrM2MI->pvLgtPrM2M)[pids[k]];
        dvar_vector lgtPrM2M = pBVI->calcArithScaleVals(pvLgtPrM2M(pids[k++]));
        ivector idxf = pBVI->getFwdIndices();//map indices to model size bins
        int imn = idxf.indexmin();
        int imx = idxf.indexmax();
        if (debug>dbgCalcProcs){
            cout<<"pc = "<<pc<<". mn = "<<imn<<", mx = "<<imx<<endl;
            cout<<"lgtPrM2M = "<<lgtPrM2M<<endl;
            cout<<"fwd indices = "<<idxf<<endl;
            ivector idxr = pBVI->getRevIndices();
            cout<<"rev indices min, max = "<<idxr.indexmin()<<cc<<idxr.indexmax()<<endl;
            cout<<"rev indices = "<<idxr<<endl;
        }

        prM2M_cz(pc) = 1.0;                                          //1 for size classes larger than max estimated
        if (idxf[1]>1) prM2M_cz(pc)(1,idxf[1]-1) = 0.0;              //0 for size classes smaller than min estimated
        dvar_vector prM2Mp = 1.0/(1.0+mfexp(-lgtPrM2M));             //logistic otherwise
        for (int i=imn;i<=imx;i++) prM2M_cz(pc)(idxf[i]) = prM2Mp[i];//assign to model size bins         
        if (debug>dbgCalcProcs){
            cout<<"prM2M = "<<prM2M_cz(pc)<<endl;
        }
        
        imatrix idxs = ptrM2MI->getModelIndices(pc);
        if (debug>dbgCalcProcs) cout<<"molt-to-maturity indices"<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1);
            if ((mnYr<=y)&&(y<=mxYr)){
                mnx = mxx = idxs(idx,2);//sex index
                if (mnx==tcsam::ALL_SXs){mnx = 1; mxx = tcsam::nSXs;}
                if (debug>dbgCalcProcs) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(mnx)<<": "<<tcsam::getSexType(mnx)<<endl;
                for (int x=mnx;x<=mxx;x++) prM2M_yxz(y,x) = prM2M_cz(pc);//note: this change made a difference, but not sure why!
            }
        }
    }
    
    if (debug>dbgCalcProcs) {
        for (int y=mnYr;y<=mxYr;y++){
            for (int x=1;x<=nSXs;x++){
                cout<<"prM2M_yxz("<<y<<tb<<x<<")="<<prM2M_yxz(y,x)<<endl;
            }
        }
        cout<<"finished calcPrM2M()"<<endl;
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
    
    GrowthInfo* ptrGrw = ptrMPI->ptrGrw;
    
    int maxZBEx = ptrMOs->maxGrowthZBEx;//maximum size bin extent for growth (was hard-wired as 11)
    
    dvariable grA;
    dvariable grB;
    dvariable grBeta;
    
    double zGrA = 0.0;
    double zGrB = 0.0;
    
    mnGrZ_cz.initialize();
    prGr_czz.initialize();
    mnGrZ_yxsz.initialize();
    prGr_yxszz.initialize();
    grA_xy.initialize();
    grB_xy.initialize();
    grBeta_xy.initialize();
    
    zGrA_xy.initialize();
    zGrB_xy.initialize();

    dvar_matrix prGr_zz(1,nZBs,1,nZBs);

    int y; int mnx; int mxx;
    dvar_vector ptGrA    = ptrGrw->pGrA->calcArithScaleVals(pGrA);
    dvar_vector ptGrB    = ptrGrw->pGrB->calcArithScaleVals(pGrB);
    dvar_vector ptGrBeta = ptrGrw->pGrBeta->calcArithScaleVals(pGrBeta);
    if (debug>dbgCalcProcs){
        cout<<"pGrA:  "<<pGrA<<endl;
        cout<<"ptGrA: "<<ptGrA<<endl;
        cout<<"pGrB:  "<<pGrB<<endl;
        cout<<"ptGrB: "<<ptGrB<<endl;
        cout<<"pGrBeta:  "<<pGrBeta<<endl;
        cout<<"ptGrBeta: "<<ptGrBeta<<endl;
    }
    for (int pc=1;pc<=ptrGrw->nPCs;pc++){
        ivector pids = ptrGrw->getPCIDs(pc);
        int k=ptrGrw->nIVs+1;//1st parameter column
        grA = ptGrA(pids[k]); k++; //"a" coefficient for mean growth
        grB = ptGrB(pids[k]); k++; //"b" coefficient for mean growth
        grBeta = ptGrBeta(pids[k]); k++; //scale factor for gamma function growth transition
        
        if (debug>dbgCalcProcs){
            cout<<"pc: "<<pc<<tb<<"grA:"<<tb<<grA<<". grB:"<<tb<<grB<<". grBeta:"<<grBeta<<endl;
        }
        
        //compute growth transition matrix for this pc
        dvar_vector mnZs(1,nZBs); mnZs.initialize();  //mean post-molt sizes with zBs as pre-molt sizes
        if (ptrMOs->optGrowthParam==0){
            //TCSAM2013 parameterization with ln-scale intercept, slope
            mnZs = mfexp(grA+grB*log(zBs));
        } else if (ptrMOs->optGrowthParam==1){
            //parameterization at min, max pre-molt sizes
            dvector pXDs = ptrGrw->getPCXDs(pc);
            zGrA = pXDs[1]; //pre-molt size corresponding to pGrA as mean post-molt size
            zGrB = pXDs[2]; //pre-molt size corresponding to pGrB as mean post-molt size
            if (debug>dbgCalcProcs) cout<<"growth parameterization 1. "<<"zGrA:"<<tb<<zGrA<<". zGrB:"<<tb<<zGrB<<endl;
            mnZs = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(zBs/zGrA));
            if (debug>dbgCalcProcs) cout<<"mnZs:"<<tb<<mnZs<<endl;
        } else if (ptrMOs->optGrowthParam==2){
            //parameterization at min pre-molt size and ln-scale slope
            dvector pXDs = ptrGrw->getPCXDs(pc);
            zGrA = pXDs[1]; //pre-molt size corresponding to pGrA as mean post-molt size
            if (debug>dbgCalcProcs) cout<<"growth parameterization 2. "<<"zGrA:"<<tb<<zGrA<<endl;
            mnZs = grA*mfexp(grB*log(zBs/zGrA));
            if (debug>dbgCalcProcs) cout<<"mnZs:"<<tb<<mnZs<<endl;
        } else {
            //throw error
            PRINT2B1(" ")
            PRINT2B1("#---------------------")
            PRINT2B2("Invalid growth parameterization option ",ptrMOs->optGrowthParam)
            PRINT2B1("Please select a valid option and re-run.")
            ad_exit(-1);
        }

        dvar_vector mnIs = mnZs - zBs;              //mean molt increments
        dvariable invBeta = 1.0/grBeta;             //inverse scale for gamma density function
        dvar_vector alIs = mnIs*invBeta;            //gamma density alpha (location) parameters
        dvar_vector mnpIs = mnZs - (zBs - 2.5);     //mean molt increment (adjusted to start of size bin)
        dvar_vector alpIs = mnpIs*invBeta;          //gamma density alpha (location) parameters
        //check all mean molt increments are > 0
        if (isnan(value(sum(sqrt(mnIs))))){
            ofstream os("GrowthReport."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
            os.precision(12);
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
        prGr_zz.initialize();
        if (ptrMOs->optGrowthPDF==0) {
            //old style (TCSAM2013)
            for (int z=1;z<=nZBs;z++){//pre-molt growth bin
                dvector dpZs =  zBs(z,nZBs) - (zBs(z)-2.5);//realized growth increments (note non-neg. growth only)
                dvar_vector prs  = elem_prod(pow(dpZs,alpIs(z)-1.0),exp(-dpZs/grBeta)); //pr(dZ|z): use exp like TCSAM2013
                int mxIZ = min(nZBs,z+maxZBEx-1);//limit growth range
                prGr_zz(z)(z,mxIZ) = prs(z,mxIZ);
                prGr_zz(z) /= sum(prGr_zz(z));
            }//zs
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else if (ptrMOs->optGrowthPDF==1){
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
                    os.precision(12);
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
                if (prs.size()>maxZBEx) prs(z+maxZBEx,nZBs) = 0.0;//limit growth range
                prs = prs/sum(prs);//normalize to sum to 1
                if (debug) cout<<prs<<endl;
                prGr_zz(z)(z,nZBs) = prs;
            }//zs
            prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
        } else {
            cout<<"Unrecognized growth option: "<<ptrMOs->optGrowthPDF<<endl;
            cout<<"Terminating!"<<endl;
            exit(-1);
        }
        
        mnGrZ_cz(pc) = mnZs;        
        prGr_czz(pc) = trans(prGr_zz);//transpose so rows are post-molt (i.e., lefthand z index is post-molt, or "to") z's so n+ = prGr_zz*n
        
        if (isnan(value(sum(prGr_zz)))){
            ofstream os("GrowthReportPrZZ."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".dat");
            os.precision(12);
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
            if (ptrMOs->optGrowthPDF==0) {
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
                    if (prs.size()>11) prs(z+11,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
                    os<<"normalization factor = "<<sum(prs)<<endl;
                    os<<"prs     : "<<prs<<endl;
                }
            } else if (ptrMOs->optGrowthPDF==1){
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
                    if (prs.size()>11) prs(z+11,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
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
        imatrix idxs = ptrGrw->getModelIndices(pc);
        if (debug) cout<<"growth indices for pc "<<pc<<endl<<idxs<<endl;
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            y = idxs(idx,1); //year index
            if ((mnYr<=y)&&(y<=mxYr)){
                mnx = mxx = idxs(idx,2);//sex index
                if (mnx==tcsam::ALL_SXs){mnx = 1; mxx = tcsam::nSXs;}
                if (debug>dbgCalcProcs) cout<<"y = "<<y<<tb<<"sex = "<<tcsam::getSexType(mnx)<<": "<<tcsam::getSexType(mnx)<<endl;
                for (int x=mnx;x<=mxx;x++){
                    grA_xy(x,y) = grA;
                    grB_xy(x,y) = grB;
                    zGrA_xy(x,y) = zGrA;
                    zGrB_xy(x,y) = zGrB;
                    grBeta_xy(x,y)  = grBeta;
                    for (int s=1;s<=nSCs;s++){
                        mnGrZ_yxsz(y,x,s) = mnGrZ_cz(pc);
                        prGr_yxszz(y,x,s) = prGr_czz(pc);
                        //for (int z=1;z<=nZBs;z++) prGr_yxszz(y,x,s,z) = prGr_czz(pc,z);
                    }//s
                }//x
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
//*  npSel_cz - nonparametric selectivities by size
//*  sel_cz   - by parameter combination and size
//*  sel_cyz - selectivity array by year and size
//******************************************************************************
FUNCTION void calcSelectivities(int debug, ostream& cout)  
    if(debug>dbgCalcProcs) cout<<"Starting calcSelectivities()"<<endl;
    
    SelectivityInfo* ptrSel = ptrMPI->ptrSel;

    double fsZ;    //fully selected size
    int idSel;     //selectivity function id
    int idxFSZ = 1;//index for fsZ in pXDs vector below
    
    ivector mniSelDevs(1,6);//min indices of devs vectors
    ivector mxiSelDevs(1,6);//max indices of devs vectors
    dvar_vector params(1,6);//vector for number_vector params
    dvar_vector paramsp(1,6);//vector for number_vector params + dev offsets
        
    npSel_cz.initialize();//nonparameteric selectivities
    sel_cz.initialize();//selectivities w/out deviations
    sel_cyz.initialize();//selectivity array

    int y;
    dvar_vector ptS1 = ptrSel->pS1->calcArithScaleVals(pS1);
    dvar_vector ptS2 = ptrSel->pS2->calcArithScaleVals(pS2);
    dvar_vector ptS3 = ptrSel->pS3->calcArithScaleVals(pS3);
    dvar_vector ptS4 = ptrSel->pS4->calcArithScaleVals(pS4);
    dvar_vector ptS5 = ptrSel->pS5->calcArithScaleVals(pS5);
    dvar_vector ptS6 = ptrSel->pS6->calcArithScaleVals(pS6);
    for (int pc=1;pc<=ptrSel->nPCs;pc++){
        ivector pids = ptrSel->getPCIDs(pc);
        dvector pXDs = ptrSel->getPCXDs(pc);
        fsZ   = pXDs[idxFSZ];
        idSel = pids[ptrSel->nIVs+ptrSel->nPVs+idxFSZ+1];
        if (debug) cout<<"pc = "<<pc<<tb<<"idSel = "<<idSel<<tb<<SelFcns::getSelFcnID(idSel)<<endl;
        if (idSel==SelFcns::ID_CONSTANT){
            //calculate constant selectivity function
            sel_cz(pc)     = SelFcns::constant(zBs);
            if (debug>dbgCalcProcs) cout<<"pc = "<<pc<<tb<<"sel_cz: "<<sel_cz(pc)<<endl;
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrSel->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                y = idxs(idx,1);//year
                if ((mnYr<=y)&&(y<=mxYr+1)) {
                    sel_cyz(pc,y) = sel_cz(pc);
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<endl<<"constant: "<<paramsp<<endl<<"sel: "<<sel_cyz(pc,y)<<endl;
                } else {
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"y outside model range--skipping year!"<<endl;
                }
            }
        } else if (idSel==SelFcns::ID_NONPARAMETRIC){
            //calculate nonparametric selectivity function
            int pcNP = pids[ptrSel->nIVs+ptrSel->nPVs];
            if (debug>dbgCalcProcs) cout<<tb<<"pcNP = "<<pcNP<<endl;
            dvar_vector nonParParams = pvNPSel(pcNP);
            if (debug>dbgCalcProcs) cout<<tb<<"nonParParams = "<<pvNPSel(pcNP)<<endl;
            int idZ = (int) fsZ;
            if (debug>dbgCalcProcs) cout<<tb<<"idZ = "<<idZ<<endl;
            npSel_cz(pcNP) = SelFcns::calcSelFcn(idSel, zBs, nonParParams, idZ);
            if (debug>dbgCalcProcs) cout<<tb<<"npSel_cz(pcNP) = "<<npSel_cz(pcNP)<<endl;
            sel_cz(pc)     = npSel_cz(pcNP);
            if (debug>dbgCalcProcs) cout<<tb<<"pcNP = "<<pcNP<<tb<<"pc = "<<pc<<tb<<"sel_cz: "<<sel_cz(pc)<<endl;
            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrSel->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                y = idxs(idx,1);//year
                if ((mnYr<=y)&&(y<=mxYr+1)) {
                    sel_cyz(pc,y) = npSel_cz(pcNP);
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<endl<<"nonParParams: "<<paramsp<<endl<<"sel: "<<sel_cyz(pc,y)<<endl;
                } else {
                    if (debug>dbgCalcProcs) cout<<tb<<"y = "<<y<<tb<<"y outside model range--skipping year!"<<endl;
                }
            }
        } else {
            //extract the number parameters
            params.initialize();
            int k=ptrSel->nIVs+1;//1st parameter variable column
            if (pids[k]) {params[1] = ptS1(pids[k]);}   k++;
            if (pids[k]) {params[2] = ptS2(pids[k]);}   k++;
            if (pids[k]) {params[3] = ptS3(pids[k]);}   k++;
            if (pids[k]) {params[4] = ptS4(pids[k]);}   k++;
            if (pids[k]) {params[5] = ptS5(pids[k]);}   k++;
            if (pids[k]) {params[6] = ptS6(pids[k]);}   k++;
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

            if (debug>dbgCalcProcs) cout<<tb<<"fsZ: "<<fsZ<<tb<<"idSel"<<tb<<idSel<<tb<<SelFcns::getSelFcnID(idSel)<<tb<<params<<endl;
            //calc selectivity function WITHOUT any annual devs
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
        }    
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
    
    dvariable hm; //handling mortality
    dvariable lnC;//ln-scale capture rate
    dvariable arC;//arithmetic-scale capture rate
    
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
    //tmF_yxmsz.initialize(); //total mortality rate
    
    /*********************************************************\n
     * Fully-selected annual capture rates are calculated     \n
     * using 2 approaches:                                    \n
     *  1. directly from parameter values (if useEX=0 below)  \n
     *  2. based on effort extrapolated to                    \n
     *      capture rate        (if useEX>1 below)            \n
     * Consequently, calculating all the above quantities     \n
     * requires 2 passes through parameter combinations.      \n
    ***********************************************************/

    int idxEX = ptrFsh->idxUseEX;//index into pids below for flag to use effort extrapolation
    int y, f, k, mnx, mxx, mnm, mxm, mns, mxs; 
    int idSel, idRet, useDevs;
    //Pass 1: calculations based on parameter values
    if (debug>dbgCalcProcs) cout<<"starting pass 1"<<endl;
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids: "<<pids<<endl;
        int useEX = pids[ptrFsh->idxUseEX];//flag to use effort extrapolation
        if (!useEX){//calculate capture rates from parameters
            lnC.initialize();
            arC.initialize();
             //get handling mortality (default to 1)
            hm = 1.0;
            k = FisheriesInfo::idxHM;
            if (pids[k]) {hm = pHM(pids[k]);}
            //set base (ln-scale) capture rate
            k = FisheriesInfo::idxLnC;
            if (pids[k]) {lnC += pLnC(pids[k]);}
            //add in offset 1 (temporal, perhaps)
            k = FisheriesInfo::idxDC1;
            if (pids[k]) {lnC += pDC1(pids[k]);}
            //add in offset 2 (for females, perhaps)
            k = FisheriesInfo::idxDC2;
            if (pids[k]) {lnC += pDC2(pids[k]);}
            //add in offset 3 (for immature crab, perhaps)
            k = FisheriesInfo::idxDC3;
            if (pids[k]) {lnC += pDC3(pids[k]);}
            //add in offset 4 (for immature females, perhaps)
            k = FisheriesInfo::idxDC4;
            if (pids[k]) {lnC += pDC4(pids[k]);}

            //extract devs vector
            k = FisheriesInfo::idxDevsLnC;
            useDevs = pids[k];
            dvar_vector dvsLnC;             
            ivector idxDevsLnC;
            if (useDevs) {
                dvsLnC     = devsLnC(useDevs);
                idxDevsLnC = idxsDevsLnC(useDevs);
            } else {
                arC = mfexp(lnC);
            }
            
            //extract logistic scale retention fraction for old shell crab
            k = FisheriesInfo::idxLgtRet;
            dvariable retFrac = 1.0;
            if (pids[k]) {
                retFrac = 1.0/(1.0+mfexp(-pLgtRet[pids[k]]));
                if (debug>dbgCalcProcs) {
                    cout<<"pc: "<<pc<<". retFrac = "<<retFrac<<endl;
                }
            } k++;

            idSel = pids[FisheriesInfo::idxSelFcn];//selectivity function id
            idRet = pids[FisheriesInfo::idxRetFcn];//retention function id

            //convert from ln-scale to arithmetic scale
            if (debug>dbgCalcProcs){
                cout<<"pc: "<<pc<<". idSel = "<<idSel<<". idRet = "<<idRet<<". lnC = "<<lnC<<". retFrac = "<<retFrac<<endl;
                if (useDevs) {
                    cout<<tb<<tb<<"dvsLnC["<<dvsLnC.indexmin()<<cc<<dvsLnC.indexmax()<<"] = "<<dvsLnC<<endl;
                } else {
                    cout<<tb<<tb<<"arC:"<<arC<<endl;
                }
            }

            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                if ((mnYr<=y)&&(y<=mxYr)){
                    hasF_fy(f,y) = 1;//flag indicating occurrence of fishery in year y
                    hmF_fy(f,y) = hm;//save discard mortality rate
                    mnx = mxx = idxs(idx,3);//sex index
                    if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
                    mnm = mxm = idxs(idx,4);//maturity index
                    if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
                    mns = mxs = idxs(idx,5);//shell condition index
                    if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
                    for (int x=mnx;x<=mxx;x++){
                        if (debug>dbgCalcProcs) cout<<"f,y,x,useDevs = "<<f<<cc<<y<<cc<<x<<cc<<useDevs<<endl;
                        if (useDevs) {
                            idxDevsLnC_fy(f,y) = idxDevsLnC[y];
                            dvsLnC_fy(f,y)     = dvsLnC[idxDevsLnC[y]];
                            arC = mfexp(lnC+dvsLnC[idxDevsLnC[y]]);//recalculate arC w/ devs
                            if (debug>dbgCalcProcs) {
                                cout<<"idxDevsLnC[y],dvsLnC[idxDevsLnC[y]], arC: "<<idxDevsLnC[y]<<tb<<dvsLnC[idxDevsLnC[y]]<<tb<<arC<<endl;
                            }
                        }
                        for (int m=mnm;m<=mxm;m++){
                            cpF_fyxms(f,y,x,m)  = arC; //fully-selected capture rate (independent of shell condition)
                            for (int s=mns;s<=mxs;s++){
                                sel_fyxmsz(f,y,x,m,s) = sel_cyz(idSel,y);          //selectivity
                                cpF_fyxmsz(f,y,x,m,s) = cpF_fyxms(f,y,x,m,s)*sel_cyz(idSel,y);//size-specific capture rate
                                if (idRet){//fishery has retention
                                    ret_fyxmsz(f,y,x,m,s) = retFrac*sel_cyz(idRet,y);      //retention curves
                                    rmF_fyxmsz(f,y,x,m,s) = elem_prod(ret_fyxmsz(f,y,x,m,s),         cpF_fyxmsz(f,y,x,m,s));//retention mortality rate
                                    dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-ret_fyxmsz(f,y,x,m,s)),cpF_fyxmsz(f,y,x,m,s));//discard mortality rate
                                } else {//discard only
                                    dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality rate
                                }
                            }//s
                        }//m
                    }//x
                }//(mnYr<=y)&&(y<=mxYr)
            }//idx
        }//useEX=FALSE
    }//pc
    if (debug>dbgCalcProcs) cout<<"finished pass 1"<<endl;
    
    //calculate ratio of average capture rate to effort
    if (debug>dbgCalcProcs) cout<<"--calculating average capture rates"<<endl;
    avgFc_nxms.initialize();
    avgFc2Eff_nxms.initialize();
    obsFc_nxmsy.initialize();
    prdFc_nxmsy.initialize();
    prdEff_nxmsy.initialize();
    CapRateAvgScenarios* pCRASs = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios;
    int nCRASs = pCRASs->nAvgs;
    for (int n=1;n<=nCRASs;n++){//capture rate averaging scenarios
        CapRateAvgScenario* ptrCRAS = pCRASs->ppCRASs[n-1];
        int idEAS = ptrCRAS->idEffAvgInfo;//index to associated average effort
        int idPar = ptrCRAS->idParam;     //index to associated extrapolation parameter
        int f     = ptrCRAS->f;  //fishery
        int fd    = mapM2DFsh(f);//index of corresponding fishery data object
        dvector eff_y = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y;//corresponding effort time series
        mnx = mxx = ptrCRAS->x;//sex index
        if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
        mnm = mxm = ptrCRAS->m;//maturity index
        if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
        mns = mxs = ptrCRAS->s;//shell index
        if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
        for (int x=mnx;x<=mxx;x++){
            for (int m=mnm;m<=mxm;m++){
                for (int s=mns;s<=mxs;s++){
                    if (debug>dbgCalcProcs)cout<<"capture rate averaging for n="<<n<<"idEAS="<<idEAS<<cc<<" avgEff_n(idEAS)="<<avgEff_n(idEAS)<<cc
                                               <<" f="<<f<<cc<<" x="<<x<<cc<<" m="<<m<<cc<<" s="<<s<<endl;
                    //loop over years to extract "observed" F's and calculate averages
                    for (int iy=yrsAvgEff_ny(idEAS).indexmin();iy<=yrsAvgEff_ny(idEAS).indexmax();iy++){
                        //cout<<iy<<cc<<yrsAvgEff_ny(idEAS,iy)<<endl;
                        obsFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = cpF_fyxms(f,yrsAvgEff_ny(idEAS,iy),x,m,s);
                    }
                    avgFc_nxms(n,x,m,s) = sum(obsFc_nxmsy(n,x,m,s))/yrsAvgEff_ny(idEAS).size();
                    avgFc2Eff_nxms(n,x,m,s) = avgFc_nxms(n,x,m,s)/avgEff_n(idEAS);
                    
                    //loop over year again to calculate "predicted" F's based on effort
                    for (int iy=yrsAvgEff_ny(idEAS).indexmin();iy<=yrsAvgEff_ny(idEAS).indexmax();iy++){
                        //fully-selected capture rate
                        if (idPar==0){
                            //extrapolation based on effort extrapolation ratio only
                            prdFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) =                      
                                                          avgFc2Eff_nxms(n,x,m,s)*eff_y(yrsAvgEff_ny(idEAS,iy)); 
                            //predicted effort from "observed" capture rates
                            prdEff_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = 
                                    obsFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy))/avgFc2Eff_nxms(n,x,m,s);
                        } else {
                            //extrapolation based on effort extrapolation parameters, as well as ratio
                            prdFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = 
                                    mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(n,x,m,s)*eff_y(yrsAvgEff_ny(idEAS,iy));
                            //predicted effort from "observed" capture rates
                            prdEff_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy)) = 
                                    obsFc_nxmsy(n,x,m,s,yrsAvgEff_ny(idEAS,iy))/(mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(n,x,m,s));
                        }
                    }
                    if (debug>dbgCalcProcs){
                        cout<<"yrsAvgEff_ny(idEAS)     = "<<yrsAvgEff_ny(idEAS)<<endl;
                        cout<<"avgFc_nxms(n,x,m,s)     = "<<avgFc_nxms(n,x,m,s)<<endl;
                        cout<<"avgFc2Eff_nxms(n,x,m,s) = "<<avgFc2Eff_nxms(n,x,m,s)<<endl;
                        cout<<"obsFc_nxmsy(n,x,m,s)    = "<<obsFc_nxmsy(n,x,m,s)<<endl;
                        cout<<"prdFc_nxmsy(n,x,m,s)    = "<<prdFc_nxmsy(n,x,m,s)<<endl;
                        cout<<"obsEff_nxmsy(n,x,m,s)   = "<<obsEff_nxmsy(n,x,m,s)<<endl;
                        cout<<"prdEff_nxmsy(n,x,m,s)   = "<<prdEff_nxmsy(n,x,m,s)<<endl;
                    }
                }//s
            }//m
        }//x
    }//n - capture rate averaging
    if (debug>dbgCalcProcs) cout<<"calculated avgFc2Eff_nxms"<<endl;
    
    //Pass 2: calculations based on effort and effort extrapolation parameters or average effort:capture rate ratios
    if (debug>dbgCalcProcs) cout<<"Starting pass 2"<<endl;
    int fd; double eff;
    for (int pc=1;pc<=ptrFsh->nPCs;pc++){
        ivector pids = ptrFsh->getPCIDs(pc);
        int useEX   = pids[FisheriesInfo::idxUseEX];  //flag to use direct effort extrapolation (+ index to EX scenario [i.e., "n"])
        int idPar   = pids[FisheriesInfo::idxLnEffX]; //index to parameters for effort extrapolation
        if (useEX){//calculate capture rates from effort extrapolation
            //get handling mortality (default to 1)
            hm = 1.0;
            if (pids[ptrFsh->idxHM]) {hm = pHM(pids[ptrFsh->idxHM]);}
            
            //extract logistic scale retention fraction for old shell crab
            dvariable retFrac = 1.0;
            k = FisheriesInfo::idxLgtRet;
            if (pids[k]) {
                retFrac = 1.0/(1.0+mfexp(-pLgtRet[pids[k]]));
                if (debug>dbgCalcProcs) {
                    cout<<"pc: "<<pc<<". retFrac = "<<retFrac<<endl;
                }
            }

            idSel = pids[FisheriesInfo::idxSelFcn];//selectivity function id
            idRet = pids[FisheriesInfo::idxRetFcn];//retention function id

            //loop over model indices as defined in the index blocks
            imatrix idxs = ptrFsh->getModelIndices(pc);
            for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
                f = idxs(idx,1);//fishery
                y = idxs(idx,2);//year
                if ((mnYr<=y)&&(y<=mxYr)){
                    hasF_fy(f,y) = 1;//flag indicating occurrence of fishery in year y
                    hmF_fy(f,y) = hm;//save discard mortality rate
                    fd = mapM2DFsh(f);//index of corresponding fishery data object
                    eff = ptrMDS->ppFsh[fd-1]->ptrEff->eff_y(y);
                    mnx = mxx = idxs(idx,3);//sex index
                    if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
                    mnm = mxm = idxs(idx,4);//maturity index
                    if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
                    mns = mxs = idxs(idx,5);//shell index
                    if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
                    if (debug>dbgCalcProcs) cout<<"f, y n, x, m, s, pLnEffX avgFc2Eff, eff cpF  = "<<endl;
                    for (int x=mnx;x<=mxx;x++){
                        for (int m=mnm;m<=mxm;m++){
                            for (int s=mns;s<=mxs;s++){
                                //fully-selected capture rate
                                if (idPar==0){
                                        //extrapolation based on effort extrapolation ratio only
                                    //prdFc_nxmsy(useEX,x,m,s,y) = avgFc2Eff_nxms(useEX,x,m,s)*eff; 
                                    //cpF_fyxms(f,y,x,m,s)   = prdFc_nxmsy(useEX,x,m,s,y); 
                                    cpF_fyxms(f,y,x,m,s) = avgFc2Eff_nxms(useEX,x,m,s)*eff; 
                                } else {
                                    //extrapolation based on effort extrapolation parameters, as well as ratio
                                    //prdFc_nxmsy(useEX,x,m,s,y) = mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(useEX,x,m,s)*eff;
                                    //cpF_fyxms(f,y,x,m,s) = prdFc_nxmsy(useEX,x,m,s,y);
                                    cpF_fyxms(f,y,x,m,s) = mfexp(pLnEffX(idPar))*avgFc2Eff_nxms(useEX,x,m,s)*eff;
                                }
                                if (debug>dbgCalcProcs) {
                                    if (idPar==0) cout<<f<<tb<<y<<useEX<<tb<<x<<tb<<0             <<avgFc2Eff_nxms(useEX,x,m,s)<<tb<<eff<<tb<<cpF_fyxms(f,y,x,m,s)<<endl;
                                    if (idPar)    cout<<f<<tb<<y<<useEX<<tb<<x<<tb<<pLnEffX(idPar)<<avgFc2Eff_nxms(useEX,x,m,s)<<tb<<eff<<tb<<cpF_fyxms(f,y,x,m,s)<<endl;
                                }
                                if (!debug) testNaNs(value(cpF_fyxms(f,y,x,m,s)),"calcFisheryFs: 2nd pass");
                                sel_fyxmsz(f,y,x,m,s) = sel_cyz(idSel,y);
                                cpF_fyxmsz(f,y,x,m,s) = cpF_fyxms(f,y,x,m,s)*sel_cyz(idSel,y);//size-specific capture rate
                                if (idRet){//fishery has retention
                                    ret_fyxmsz(f,y,x,m,s) = retFrac*sel_cyz(idRet,y);
                                    rmF_fyxmsz(f,y,x,m,s) = elem_prod(ret_fyxmsz(f,y,x,m,s),         cpF_fyxmsz(f,y,x,m,s));//retention mortality rate
                                    dmF_fyxmsz(f,y,x,m,s) = elem_prod(hm*(1.0-ret_fyxmsz(f,y,x,m,s)),cpF_fyxmsz(f,y,x,m,s));//discard mortality rate
                                } else {//discard only
                                    dmF_fyxmsz(f,y,x,m,s) = hm*cpF_fyxmsz(f,y,x,m,s);//discard mortality rate
                                }
                            }//s
                        }//m
                    }//x
                }//(mnYr<=y)&&(y<=mxYr)
            }//idx
        }//useEX=TRUE
    }//pc
    if (debug>dbgCalcProcs) cout<<"Finished pass 2"<<endl;
    
    if (debug>dbgCalcProcs) {
        for (int f=1;f<=nFsh;f++){
            cout<<"cpF_fyxmsz("<<f<<",mnYr:mxYr,  MALE,MATURE,NEW SHELL) = ";
            for (int y=mnYr;y<=mxYr;y++) {cout<<cpF_fyxmsz(f,y,  MALE,MATURE,NEW_SHELL)<<tb;} cout<<endl;
            cout<<"cpF_fyxmsz("<<f<<",mnYr:mxYr,FEMALE,MATURE,NEW SHELL) = ";
            for (int y=mnYr;y<=mxYr;y++) {cout<<cpF_fyxmsz(f,y,FEMALE,MATURE,NEW_SHELL)<<tb;} cout<<endl;
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
//*  q_vyxms  - fully-selected catchability by survey/year/sex/maturity state/shell condition
//*  q_vyxmsz - size-specific catchability by survey/year/sex/maturity state/shell condition
//*  a_vyxms  - max availability by survey/year/sex/maturity state/shell condition
//*  a_vyxmsz - size-specific availability by survey/year/sex/maturity state/shell condition
//******************************************************************************
FUNCTION void calcSurveyQs(int debug, ostream& cout)  
    if(debug>dbgCalcProcs) cout<<"Starting calcSurveyQs()"<<endl;
    
    SurveysInfo* ptrSrv = ptrMPI->ptrSrv;
    
    dvariable lnA;//ln-scale A's
    dvariable arA;//arithmetic-scale A's
    dvariable lnQ;//ln-scale Q's
    dvariable arQ;//arithmetic-scale Q's
    
    a_vyxms.initialize();
    a_vyxmsz.initialize();
    q_vyxms.initialize();
    q_vyxmsz.initialize();

    int y; int v; int mnx; int mxx; int mnm; int mxm; int mns; int mxs;
    int idAvl; int idSel;
    dvar_vector ptA  = ptrSrv->pA->calcArithScaleVals(pA);
    dvar_vector ptQ  = ptrSrv->pQ->calcArithScaleVals(pQ);
    for (int pc=1;pc<=ptrSrv->nPCs;pc++){
        lnA.initialize();
        arA.initialize();
        lnQ.initialize();
        arQ.initialize();
        ivector pids = ptrSrv->getPCIDs(pc);
        if (debug>dbgCalcProcs) cout<<"pc: "<<pc<<tb<<"pids = "<<pids<<endl;
        
        int k=ptrSrv->nIVs+1;//1st parameter variable column
        //add in base catchability (e.g., for mature male)
        if (pids[k]) {lnQ += log(ptQ(pids[k]));} k++;
        //add in ln-scale offset 1            (e.g., for time period)
        if (pids[k]) {lnQ += pDQ1(pids[k]);} k++;
        //add in ln-scale offset 2            (e.g., for females)
        if (pids[k]) {lnQ += pDQ2(pids[k]);} k++;
        //add in ln-scale offset 3            (e.g., for immature crab)
        if (pids[k]) {lnQ += pDQ3(pids[k]);} k++;
        //add in ln-scale offset 4            (e.g., for immature females)
        if (pids[k]) {lnQ += pDQ4(pids[k]);} k++; 
        
        //get log-scale max availability
        if (pids[k]) {lnA += log(ptA(pids[k]));} k++;
        
        idAvl = pids[k++];//availability function id
        idSel = pids[k++];//selectivity function id
        
        //convert from ln-scale to arithmetic scale
        arA = mfexp(lnA);
        arQ = mfexp(lnQ);
        if (debug>dbgCalcProcs){
            cout<<"lnA:"<<lnA<<tb<<"arA:"<<endl<<arA<<endl;
            cout<<"lnQ:"<<lnQ<<tb<<"arQ:"<<endl<<arQ<<endl;
        }
        
        //loop over model indices as defined in the index blocks
        imatrix idxs = ptrSrv->getModelIndices(pc);
        for (int idx=idxs.indexmin();idx<=idxs.indexmax();idx++){
            v = idxs(idx,1);//survey
            y = idxs(idx,2);//year
            if ((mnYr<=y)&&(y<=mxYrp1)){
                mnx = mxx = idxs(idx,3);//sex index
                if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
                mnm = mxm = idxs(idx,4);//maturity state index
                if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
                mns = mxs = idxs(idx,5);//shell condition index
                if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
                for (int x=mnx;x<=mxx;x++){
                    for (int m=mnm;m<=mxm;m++){
                        a_vyxms(v,y,x,m) = arA;//max availability
                        q_vyxms(v,y,x,m) = arQ;//fully-selected catchability
                        for (int s=mns;s<=mxs;s++){
                            a_vyxmsz(v,y,x,m,s) = arA;//default: availability = 1
                            if (idAvl>0) a_vyxmsz(v,y,x,m,s) = sel_cyz(idAvl,y);
                            if (idSel>0) s_vyxmsz(v,y,x,m,s) = sel_cyz(idSel,y);
                            q_vyxmsz(v,y,x,m,s) = arA*arQ*elem_prod(a_vyxmsz(v,y,x,m,s),s_vyxmsz(v,y,x,m,s));
                        }//s
                    }//m
                }//x
            }//(mnYr<=y)&&(y<=mxYrp1)
        }//idx
    }//pc
    
    if (debug>dbgCalcProcs) cout<<"finished calcSurveyQs()"<<endl;
    
//-------------------------------------------------------------------------------------
//calculate equilibrium size distribution for unexploited population
FUNCTION void calcEqNatZF100(dvariable& R, int yr, int debug, ostream& cout)
    if (debug>=dbgPopDy) cout<<"starting void calcEqNatZF100(R,yr)"<<endl;

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
    
    if (debug>=dbgPopDy) cout<<"finished void calcEqNatZF100(R,yr)"<<endl;

//-------------------------------------------------------------------------------------
//calculate sex-specific equilibrium size distribution
FUNCTION dvar3_array calcEqNatZ(dvar_vector& R_z,dvar3_array& S1_msz, dvar_matrix& Th_sz, dvar3_array& T_szz, dvar3_array& S2_msz, int debug, ostream& cout)
    RETURN_ARRAYS_INCREMENT();
    if (debug>=dbgPopDy) cout<<"starting dvar3_array calcEqNatZ(...)"<<endl;

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
        
    if (debug>=dbgPopDy) cout<<"finished dvar3_array calcEqNatZ(...)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return(n_msz);

//-------------cohort progression calculations--------------
///**
// * Calculate cohort progression.
// * 
// * @param yr  - year determining population process rates
// * @param nzp - number of size bins to include for initial recruitment
// * @param includeM - flag (0,1) to include natural mortality
// * @param includeF - flag (0,1) to include fishing mortality
// * @param debug - flag to print debugging info
// * @param cout - output stream to write to
// * 
// * @return 5d array of cohort abundance by yxmsz (n_yxmsz)
// */
FUNCTION d5_array calcCohortProgression(int yr, int nzp, int includeM, int includeF, int debug, ostream& cout)
    if (debug) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcCohortProgression(...)"<<endl;
        cout<<"year for progression = "<<yr<<endl;
    }

    int nyp = 20;//number of years to track cohorts
    //1. set initial cohort abundance
    dvar5_array n_yxmsz(0,nyp,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n_yxmsz.initialize();
    n_yxmsz(0,  MALE,IMMATURE,NEW_SHELL)(1,nzp) = R_yz(yr)(1,nzp)/sum(R_yz(yr)(1,nzp));//set initial abundance-at-size
    n_yxmsz(0,FEMALE,IMMATURE,NEW_SHELL)(1,nzp) = R_yz(yr)(1,nzp)/sum(R_yz(yr)(1,nzp));//set initial abundance-at-size
    
    //2. Determine population rates, based on yr
    double dtF = dtF_y(yr);//time at which fisheries occur
    double dtM = dtM_y(yr);//time at which mating occurs
    
    PopDyInfo* pPIM = new PopDyInfo(nZBs);//  males info
    pPIM->R_z   = R_yz(yr);//NOTE: this is not used in calculation
    pPIM->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(MALE);
    pPIM->M_msz = includeM*M_yxmsz(yr,MALE);
    pPIM->T_szz = prGr_yxszz(yr,MALE);
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = prM2M_yxz(yr,MALE);
    
    PopDyInfo* pPIF = new PopDyInfo(nZBs);//females info
    pPIF->R_z   = R_yz(yr);
    pPIF->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(FEMALE);
    pPIF->M_msz = includeM*M_yxmsz(yr,FEMALE);
    pPIF->T_szz = prGr_yxszz(yr,FEMALE);
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = prM2M_yxz(yr,FEMALE);
    if (debug) cout<<"calculated pPIM, pPIF."<<endl;
    
    //3. Determine fishery conditions based on yr
    dvar5_array avgCapF_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery capture rates
    dvar5_array avgRFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery retention functions
    dvar5_array avgSFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery selectivity functions
    avgCapF_xfmsz.initialize();
    avgRFcn_xfmsz.initialize();
    avgSFcn_xfmsz.initialize();
    dvar_vector avgHM_f(1,nFsh);
    avgHM_f.initialize();
    if (includeF){
        int avgPeriodYrs = 1;  
        //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
        //a. Calculate average handling mortality, retention curves and capture rates
        int ny;   //number of years fishery is active
        for (int f=1;f<=nFsh;f++){
            ny = 0;
            for (int y=yr-avgPeriodYrs+1;y<=yr;y++){
                ny         += hasF_fy(f,y);
                avgHM_f(f) += hmF_fy(f,y);
            }
            avgHM_f(f) /= wts::max(1.0,1.0*ny);
        }
        if (debug) cout<<"avgHm_f = "<<avgHM_f<<endl;

        for (int f=1;f<=nFsh;f++){
            avgPeriodYrs = 1;
            if (debug) cout<<"avgPeriodYrs("<<f<<") = "<<avgPeriodYrs<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) {
                        for (int z=1;z<=nZBs;z++){
                            ny = 0;
                            for (int y=(yr-avgPeriodYrs+1);y<=yr;y++) {
                                //if (debug) cout<<"y = "<<y<<endl;
                                ny += hasF_fy(f,y);
                                avgCapF_xfmsz(x,f,m,s,z) += cpF_fyxmsz(f,y,x,m,s,z);
                                avgRFcn_xfmsz(x,f,m,s,z) += ret_fyxmsz(f,y,x,m,s,z);
                                avgSFcn_xfmsz(x,f,m,s,z) += sel_fyxmsz(f,y,x,m,s,z);
                            }//y
                            double fac = 1.0/max(1.0,1.0*ny);
                            avgCapF_xfmsz(x,f,m,s,z) *= fac;
                            avgRFcn_xfmsz(x,f,m,s,z) *= fac;
                            avgSFcn_xfmsz(x,f,m,s,z) *= fac;
                        }//z
                    }//s
                }//m
            }//x
        }//f
        if (debug){
            for (int f=1;f<=nFsh;f++){
                cout<<"avgCapF_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgCapF_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgSFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgSFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
            }
        }
    }

    CatchInfo* pCIM = new CatchInfo(nZBs,nFsh);//male catch info
    pCIM->setCaptureRates(avgCapF_xfmsz(MALE));
    pCIM->setRetentionFcns(avgRFcn_xfmsz(MALE));
    pCIM->setHandlingMortality(avgHM_f);
    dvariable maxCapF = pCIM->findMaxTargetCaptureRate(cout);
    if (debug) cout<<"maxCapF = "<<maxCapF<<endl;

    CatchInfo* pCIF = new CatchInfo(nZBs,nFsh);//female catch info
    pCIF->setCaptureRates(avgCapF_xfmsz(FEMALE));
    pCIF->setRetentionFcns(avgRFcn_xfmsz(FEMALE));
    pCIF->setHandlingMortality(avgHM_f);
    pCIF->maxF = maxCapF;//need to set this for females

    //4. Create PopProjectors
    PopProjector* pPPM = new PopProjector(pPIM,pCIM);
    pPPM->dtF = dtF;
    pPPM->dtM = dtM;
    PopProjector* pPPF = new PopProjector(pPIF,pCIF);
    if (debug) cout<<"created pPPs."<<endl;
    pPPF->dtF = dtF;
    pPPF->dtM = dtM;

    //5. Create multi-year population projectors
    MultiYearPopProjector* pMYPPM = new MultiYearPopProjector(pPPM);
    MultiYearPopProjector* pMYPPF = new MultiYearPopProjector(pPPF);
    if (debug) cout<<"created pMPPs."<<endl;
    
    //6.project with no recruitment
    pMYPPM->project(nyp,0.0,maxCapF,n_yxmsz(0,  MALE),cout);
    pMYPPF->project(nyp,0.0,maxCapF,n_yxmsz(0,FEMALE),cout);
    if (debug) cout<<"projected cohort."<<endl;
    for (int y=0;y<=nyp;y++){
        n_yxmsz(y,  MALE) = pMYPPM->n_ymsz(y);
        n_yxmsz(y,FEMALE) = pMYPPF->n_ymsz(y);
    }
    d5_array vn_yxmsz = wts::value(n_yxmsz);
    
    cout<<"finished calcCohortProgression(...)"<<endl<<endl<<endl;
    return vn_yxmsz;
//-------------END cohort progression calculations----------   
        
//-------------OFL Calculations--------------
FUNCTION void calcOFL(int yr, int debug, ostream& cout)
    if (debug) {
        cout<<endl<<endl<<"#------------------------"<<endl;
        cout<<"starting calcOFL(yr,debug,cout)"<<endl;
        cout<<"year for projection = "<<yr<<endl;
    }

    //1. get initial population for "upcoming" year, yr
    dvar4_array n_xmsz = n_yxmsz(yr);
    if (debug) {cout<<"  males_msz:"<<endl; wts::print(n_xmsz(  MALE),cout,1);}
    if (debug) {cout<<"females_msz:"<<endl; wts::print(n_xmsz(FEMALE),cout,1);}
    
    //2. set yr back one year to get population rates, etc., 
    //   from year prior to projection year
    yr = yr-1;//don't have pop rates, etc. for projection year
    if (debug) cout<<"year for pop rates = "<<yr<<endl;
    
    //3. Determine mean recruitment 
    //   1981 here corresponds to 1982 in TCSAM2013, the year recruitment enters
    //   the model population.
    dvar_vector avgRec_x(1,nSXs);
    if (debug) cout<<"R dims: "<<R_y.indexmin()<<cc<<R_y.indexmax()<<endl;
    for (int x=1;x<=nSXs;x++) 
        avgRec_x(x)= mean(elem_prod(R_y(1981,yr),column(R_yx,x)(1981,yr)));
    if (debug) {
        cout<<"R_y(  1981:"<<yr<<")      = "<<R_y(1981,yr)<<endl;
        cout<<"R_yx((1981:"<<yr<<",MALE) = "<<column(R_yx,MALE)(1981,yr)<<endl;
        cout<<"Average recruitment = "<<avgRec_x<<endl;
    }

    //4. Determine population rates for next year, using yr
    double dtF = dtF_y(yr);//time at which fisheries occur
    double dtM = dtM_y(yr);//time at which mating occurs
    
    PopDyInfo* pPIM = new PopDyInfo(nZBs);//  males info
    pPIM->R_z   = R_yz(yr);
    pPIM->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(MALE);
    pPIM->M_msz = M_yxmsz(yr,MALE);
    pPIM->T_szz = prGr_yxszz(yr,MALE);
    for (int s=1;s<=nSCs;s++) pPIM->Th_sz(s) = prM2M_yxz(yr,MALE);
    
    PopDyInfo* pPIF = new PopDyInfo(nZBs);//females info
    pPIF->R_z   = R_yz(yr);
    pPIF->w_mz  = ptrMDS->ptrBio->wAtZ_xmz(FEMALE);
    pPIF->M_msz = M_yxmsz(yr,FEMALE);
    pPIF->T_szz = prGr_yxszz(yr,FEMALE);
    for (int s=1;s<=nSCs;s++) pPIF->Th_sz(s) = prM2M_yxz(yr,FEMALE);
    if (debug) cout<<"calculated pPIM, pPIF."<<endl;
    
    //5. Determine fishery conditions for next year based on averages for recent years
        int oflAvgPeriodYrs = 5;  //TODO: this should be an input
        //assumption here is that ALL fisheries EXCEPT the first are bycatch fisheries
        //a. Calculate average handling mortality, retention curves and capture rates
        int ny;   //number of years fishery is active
        dvar_vector avgHM_f(1,nFsh);
        avgHM_f.initialize();
        for (int f=1;f<=nFsh;f++){
            ny = 0;
            for (int y=yr-oflAvgPeriodYrs+1;y<=yr;y++){
                ny         += hasF_fy(f,y);
                avgHM_f(f) += hmF_fy(f,y);
            }
            avgHM_f(f) /= wts::max(1.0,1.0*ny);
        }
        if (debug) cout<<"avgHm_f = "<<avgHM_f<<endl;

        dvar4_array avgCapF_xfms(1,nSXs,1,nFsh,1,nMSs,1,nSCs);//averaged max fishery capture rates
        dvar5_array avgCapF_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery capture rates
        dvar5_array avgRFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery retention functions
        dvar5_array avgSFcn_xfmsz(1,nSXs,1,nFsh,1,nMSs,1,nSCs,1,nZBs);//averaged fishery selectivity functions
        avgCapF_xfms.initialize();
        avgCapF_xfmsz.initialize();
        avgRFcn_xfmsz.initialize();
        avgSFcn_xfmsz.initialize();
        for (int f=1;f<=nFsh;f++){
            oflAvgPeriodYrs = ptrMOs->oflNumYrsForAvgCapRate(f);
            if (debug) cout<<"oflAvgPeriodYrs("<<f<<") = "<<oflAvgPeriodYrs<<endl;
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) {
                        for (int z=1;z<=nZBs;z++){
                            ny = 0;
                            for (int y=(yr-oflAvgPeriodYrs+1);y<=yr;y++) {
                                //if (debug) cout<<"y = "<<y<<endl;
                                ny += hasF_fy(f,y);
                                avgCapF_xfms(x,f,m,s) += cpF_fyxms(f,y,x,m,s);
                                avgCapF_xfmsz(x,f,m,s,z) += cpF_fyxmsz(f,y,x,m,s,z);
                                avgRFcn_xfmsz(x,f,m,s,z) += ret_fyxmsz(f,y,x,m,s,z);
                                avgSFcn_xfmsz(x,f,m,s,z) += sel_fyxmsz(f,y,x,m,s,z);
                            }//y
                            double fac = 1.0/max(1.0,1.0*ny);
                            avgCapF_xfms(x,f,m,s) *= fac;
                            avgCapF_xfmsz(x,f,m,s,z) *= fac;
                            avgRFcn_xfmsz(x,f,m,s,z) *= fac;
                            avgSFcn_xfmsz(x,f,m,s,z) *= fac;
                        }//z
                    }//s
                }//m
            }//x
        }//f
        if (debug){
            for (int f=1;f<=nFsh;f++){
                cout<<"avgCapF_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgCapF_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgCapF_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgRFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgRFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgSFcn_xfmsz(  MALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(  MALE,f,MATURE,NEW_SHELL)<<endl;
                cout<<"avgSFcn_xfmsz(FEMALE,"<<f<<",MATURE,NEW_SHELL) = "<<endl<<tb<<avgSFcn_xfmsz(FEMALE,f,MATURE,NEW_SHELL)<<endl;
            }
        }
        
        //determine average capture rates
        for (int f=1;f<=nFsh;f++){
            if (ptrMOs->optOFLAvgCapRate(f)==0){
                //use averaged selectivity functions, max capture rates
                if (debug) cout<<"Using average selectivity functions for fishery "<<f<<" for OFL calculations"<<endl;
                for (int x=1;x<=nSXs;x++){
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) avgCapF_xfmsz(x,f,m,s) = avgCapF_xfms(x,f,m,s)*avgSFcn_xfmsz(x,f,m,s);
                    }//m
                }//x
            } else {
                if (debug) cout<<"Using size-specific average capture rates for fishery "<<f<<" for OFL calculations"<<endl;
                //use size-specific averaged capture rates
                //do nothing
            }
        }//f
        
        CatchInfo* pCIM = new CatchInfo(nZBs,nFsh);//male catch info
        pCIM->setCaptureRates(avgCapF_xfmsz(MALE));
        pCIM->setCaptureRates(avgCapF_xfms(MALE));
        pCIM->setSelectivityFcns(avgSFcn_xfmsz(MALE));
        pCIM->setRetentionFcns(avgRFcn_xfmsz(MALE));
        pCIM->setHandlingMortality(avgHM_f);
        dvariable maxCapF = pCIM->findMaxTargetCaptureRate(cout);
        if (debug) cout<<"maxCapF = "<<maxCapF<<endl;
        
        CatchInfo* pCIF = new CatchInfo(nZBs,nFsh);//female catch info
        pCIF->setCaptureRates(avgCapF_xfmsz(FEMALE));
        pCIF->setCaptureRates(avgCapF_xfms(FEMALE));
        pCIF->setSelectivityFcns(avgSFcn_xfmsz(FEMALE));
        pCIF->setRetentionFcns(avgRFcn_xfmsz(FEMALE));
        pCIF->setHandlingMortality(avgHM_f);
        pCIF->maxF = maxCapF;//need to set this for females
        
    //6. Create PopProjectors
        PopProjector* pPPM = new PopProjector(pPIM,pCIM);
        pPPM->dtF = dtF;
        pPPM->dtM = dtM;
        PopProjector* pPPF = new PopProjector(pPIF,pCIF);
        pPPF->dtF = dtF;
        pPPF->dtM = dtM;
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
                OFL_Calculator::debug=1;
                Tier3_Calculator::debug=1;
                Equilibrium_Calculator::debug=0;
                cout<<"Calculating ptrOFLResults"<<endl;
            }
            ptrOFLResults = pOC->calcOFLResults(avgRec_x,n_xmsz,cout);
            if (debug) {
                cout<<"calculated ptrOFLResults->"<<endl;
                ptrOFLResults->writeCSVHeader(cout); cout<<endl;
                ptrOFLResults->writeToCSV(cout); cout<<endl;
                ptrOFLResults->writeToR(cout,ptrMC,"oflResults",0); cout<<endl;
                OFL_Calculator::debug=0;
                Tier3_Calculator::debug=0;
                Equilibrium_Calculator::debug=0;
            }
        }//Tier 3 calculation
    
    if (debug) {
        int n = 100;
        MultiYearPopProjector* pMYPPM = new MultiYearPopProjector(pPPM);
        MultiYearPopProjector* pMYPPF = new MultiYearPopProjector(pPPF);
        pMYPPM->projectUnFished(n,avgRec_x(  MALE),n_xmsz(  MALE),cout);
        dvariable myPP_B0 = pMYPPM->matBio_y(n);
        pMYPPM->project(n,avgRec_x(  MALE),ptrOFLResults->Fmsy,n_xmsz(  MALE),cout);
        pMYPPF->project(n,avgRec_x(FEMALE),ptrOFLResults->Fmsy,n_xmsz(FEMALE),cout);
        dvariable myPP_Bmsy = pMYPPM->matBio_y(n);
        dvariable myPP_prjB = pMYPPM->matBio_y(1);
        dvariable myPP_MSY  = pMYPPM->totCM_y(n)+pMYPPF->totCM_y(n);
        dvariable myPP_OFL  = pMYPPM->totCM_y(1)+pMYPPF->totCM_y(1);
        cout<<"#------------------------"<<endl;
        cout<<"MYPP B0     = "<<myPP_B0  <<". B0     = "<<ptrOFLResults->B0  <<endl;
        cout<<"MYPP Bmsy   = "<<myPP_Bmsy<<". Bmsy   = "<<ptrOFLResults->Bmsy<<endl;
        cout<<"MYPP prjMMB = "<<myPP_prjB<<". prjMMB = "<<ptrOFLResults->prjB<<endl;
        cout<<"MYPP MSY    = "<<myPP_MSY <<". MSY    = "<<ptrOFLResults->MSY <<endl;
        cout<<"MYPP OFL    = "<<myPP_OFL <<". OFL    = "<<ptrOFLResults->OFL <<endl;
        cout<<"#------------------------"<<endl;
        cout<<"finished calcOFL(yr,debug,cout)"<<endl<<endl<<endl;
    }
//-------------END OFL Calculations----------   
        
//-------------------------------------------------------------------------------------
//Calculate penalties for objective function. TODO: finish
FUNCTION void calcPenalties(int debug, ostream& cout)
    if (debug>dbgObjFun) cout<<"Started calcPenalties()"<<endl;
    if (debug<0) cout<<"list("<<endl;//start list of penalties by category

    //growth-related penalties
    if (debug>dbgObjFun) cout<<"calculating growth-related penalties"<<endl;
    if (debug<0) cout<<tb<<"growth=list(negativeGrowth=list("<<endl;//start of growth penalties list
    GrowthInfo* ptrGrw = ptrMPI->ptrGrw;
    //calculate growth parameters on arithmetic scale
    dvar_vector ptGrA = ptrGrw->pGrA->calcArithScaleVals(pGrA);
    dvar_vector ptGrB = ptrGrw->pGrB->calcArithScaleVals(pGrB);
    //loop over parameter combinations
    dvariable grA, grB;
    dvector zBsp(1,nZBs+2);
    dvar_vector dZ(1,nZBs+2);
    zBsp(1,nZBs) = zBs;
    zBsp(nZBs+1) = ptrMOs->minGrowthCW;
    zBsp(nZBs+2) = ptrMOs->maxGrowthCW;
    for (int pc=1;pc<=(ptrGrw->nPCs-1);pc++){
        ivector pids = ptrGrw->getPCIDs(pc);
        int k=ptrGrw->nIVs+1;//1st parameter column
        grA = ptGrA(pids[k]); k++; //"a" coefficient for mean growth
        grB = ptGrB(pids[k]); k++; //"b" coefficient for mean growth
                
        //add objective function penalties to keep mean size increments positive
        dvariable pen; pen.initialize();
        dZ.initialize();
        if (ptrMOs->optGrowthParam==0){
            dZ = mfexp(grA+grB*log(zBsp)) - zBsp;
        } else if (ptrMOs->optGrowthParam==1){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(zBsp/zGrA)) - zBsp;
        } else if (ptrMOs->optGrowthParam==2){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(grB*log(zBsp/zGrA)) - zBsp;
        }
        posfun(dZ,ptrMOs->epsNegGrowth,pen);
        if (pen>0.0){
            rpt::echo<<"--Growth Increments Approaching 0 for pc = "<<pc<<". pen = "<<pen<<endl;
            rpt::echo<<"params = "<<grA<<cc<<grB<<endl;
            rpt::echo<<"zBsp = "<<zBsp<<endl;
            rpt::echo<<"dZ   = "<<dZ<<endl;
            rpt::echo<<"--------"<<endl;
        }
        objFun += ptrMOs->wgtNegGrowth*pen;
        if (debug<0) {
            cout<<tb<<tb<<tb<<"'"<<pc<<"'=list(wgt="<<ptrMOs->wgtNegGrowth<<cc<<"pen="<<pen<<cc<<"objfun="<<ptrMOs->wgtNegGrowth*pen<<"),"<<endl;
        }
    }
    {
        int pc = ptrGrw->nPCs;
        ivector pids = ptrGrw->getPCIDs(pc);
        int k=ptrGrw->nIVs+1;//1st parameter column
        grA = (*ptrMPI->ptrGrw->pGrA)[pids[k]]->calcArithScaleVal(pGrA(pids[k])); k++; //"a" coefficient for mean growth
        grB = (*ptrMPI->ptrGrw->pGrB)[pids[k]]->calcArithScaleVal(pGrB(pids[k])); k++; //"b" coefficient for mean growth
        
        //add objective function penalty to keep mean size increments positive
        dvariable pen; pen.initialize();
        dZ.initialize();
        if (ptrMOs->optGrowthParam==0){
            dZ = mfexp(grA+grB*log(zBsp)) - zBsp;
        } else if (ptrMOs->optGrowthParam==1){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(zBsp/zGrA)) - zBsp;
        } else if (ptrMOs->optGrowthParam==2){
            dvector pXDs = ptrGrw->getPCXDs(pc);
            double zGrA = pXDs[1];
            double zGrB = pXDs[2];
            dZ = grA*mfexp(grB*log(zBsp/zGrA)) - zBsp;
        }
        posfun(dZ,ptrMOs->epsNegGrowth,pen);
        if (pen>0.0){
            rpt::echo<<"--Growth Increments Approaching 0 for pc = "<<pc<<". pen = "<<pen<<endl;
            rpt::echo<<"params = "<<grA<<cc<<grB<<endl;
            rpt::echo<<"zBsp = "<<zBsp<<endl;
            rpt::echo<<"dZ   = "<<dZ<<endl;
            rpt::echo<<"--------"<<endl;
        }
        objFun += ptrMOs->wgtNegGrowth*pen;
        if (debug<0) {
            cout<<tb<<tb<<tb<<"'"<<pc<<"'=list(wgt="<<ptrMOs->wgtNegGrowth<<cc<<"pen="<<pen<<cc<<"objfun="<<ptrMOs->wgtNegGrowth*pen<<")"<<endl;
        }
    }
    if (debug<0) cout<<tb<<tb<<"))"<<cc<<"#end of growth-related penalties"<<endl;//end of growth penalties list
    
    //maturity-related penalties
    if (debug>dbgObjFun) cout<<"calculating maturity-related penalties"<<endl;
    if (debug<0) cout<<tb<<"maturity=list("<<endl;//start of maturity penalties list
    //smoothness penalties
    dvector penWgtSmthLgtPrMat = ptrMOs->wgtPenSmthPrM2M;
    if (ptrMOs->optPenSmthPrM2M==0){
        if (debug>dbgObjFun) cout<<"--calculating smoothness penalties on parameters"<<endl;
        //smoothness penalties on maturity PARAMETERS (NOT maturity ogives)
        fPenSmoothLgtPrMat.initialize();
        if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
        for (int i=1;i<npLgtPrMat;i++){
            dvar_vector v = 1.0*pvLgtPrM2M(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<"),"<<endl;
        }
        {
            int i = npLgtPrMat;
            dvar_vector v = 1.0*pvLgtPrM2M(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<")"<<endl;
        }
        if (debug<0) cout<<tb<<tb<<")"<<cc<<"#end of smoothness penalties"<<endl;//end of smoothness penalties list
    } else if (ptrMOs->optPenSmthPrM2M==1){
        //smoothness penalties on maturity OGIVES (NOT maturity parameters)
        if (debug>dbgObjFun) cout<<"calculating smoothness penalties on ogives"<<endl;
        fPenSmoothLgtPrMat.initialize();
        if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
        for (int i=1;i<ptrMPI->ptrM2M->nPCs;i++){
            dvar_vector v = 1.0*prM2M_cz(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<"),"<<endl;
        }
        {
            int i = ptrMPI->ptrM2M->nPCs;
            dvar_vector v = 1.0*prM2M_cz(i);
            fPenSmoothLgtPrMat(i) = 0.5*norm2(calc2ndDiffs(v));
            objFun += penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i);
            if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthLgtPrMat(i)<<cc<<"pen="<<fPenSmoothLgtPrMat(i)<<cc<<"objfun="<<penWgtSmthLgtPrMat(i)*fPenSmoothLgtPrMat(i)<<")"<<endl;
        }
        if (debug<0) cout<<tb<<tb<<")"<<cc<<"#end of smoothness penalties"<<endl;//end of smoothness penalties list
    }
    //penalties for decreasing maturity parameters/ogives
    dvector penWgtNonDecLgtPrMat = ptrMOs->wgtPenNonDecPrM2M;
    fPenNonDecLgtPrMat.initialize();
    if (debug<0) cout<<tb<<tb<<"nondecreasing=list(";//start of non-decreasing penalties list
    int np;
    if (ptrMOs->optPenNonDecPrM2M==0||ptrMOs->optPenNonDecPrM2M==1) np = npLgtPrMat;
    else if (ptrMOs->optPenNonDecPrM2M==2||ptrMOs->optPenNonDecPrM2M==3) np = ptrMPI->ptrM2M->nPCs;
    for (int i=1;i<np;i++){
        dvar_vector v; 
        if (ptrMOs->optPenNonDecPrM2M==0){
            v = calc1stDiffs(pvLgtPrM2M(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            } 
        } else if (ptrMOs->optPenNonDecPrM2M==1){
            v = calc1stDiffs(pvLgtPrM2M(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        } else if (ptrMOs->optPenNonDecPrM2M==2){
            v = calc1stDiffs(prM2M_cz(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            } 
        } else if (ptrMOs->optPenNonDecPrM2M==3){
            v = calc1stDiffs(prM2M_cz(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        }
        objFun += penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat(i)<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i)<<"),";
    }
    {
        int i = np;
        dvar_vector v; 
        if (ptrMOs->optPenNonDecPrM2M==0){
            v = calc1stDiffs(pvLgtPrM2M(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            }
        } else if (ptrMOs->optPenNonDecPrM2M==1){
            v = calc1stDiffs(pvLgtPrM2M(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        } else if (ptrMOs->optPenNonDecPrM2M==2){
            v = calc1stDiffs(prM2M_cz(i));
            for (int iv=v.indexmin();iv<=v.indexmax();iv++){
                posfun2(v(iv),1.0E-2,fPenNonDecLgtPrMat(i));
            } 
        } else if (ptrMOs->optPenNonDecPrM2M==3){
            v = calc1stDiffs(prM2M_cz(i));
            fPenNonDecLgtPrMat(i) = sum(mfexp(-10.0*v));
        }
        objFun += penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i);
        if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtNonDecLgtPrMat(i)<<cc<<"pen="<<fPenNonDecLgtPrMat(i)<<cc<<"objfun="<<penWgtNonDecLgtPrMat(i)*fPenNonDecLgtPrMat(i)<<")";
    }
    if (debug<0) cout<<tb<<tb<<") "<<"#end of non-decreasing penalties"<<endl;//end of non-decreasing penalties list    
    if (debug<0) cout<<tb<<")"<<cc<<"#end of maturity penalties"<<endl;//end of maturity penalties list
    
    //penalties on nonparametric selectivity functions
    if (debug<0) cout<<tb<<"nonParSelFcns=list("<<endl;//start of nonparametric selectivity function penalties list
    //smoothness penalties
    fPenSmoothNPSel.initialize();
    if (ptrMPI->ptrSel->pvNPSel->getSize()>0){
        if (debug<0) cout<<tb<<tb<<"smoothness=list(";//start of smoothness penalties list
        dvector penWgtSmthNPSel = ptrMOs->wgtPenSmthNPSel;
        if (ptrMOs->optPenSmthNPSel==0){
            //smoothness penalties on selectivity PARAMETERS (NOT selectivity functions)
            for (int i=1;i<npNPSel;i++){
                dvar_vector v = 1.0*pvNPSel(i);
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(v));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<"),"<<endl;
            }
            {
                int i = npNPSel;
                dvar_vector v = 1.0*pvNPSel(i);
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(v));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<")"<<endl;
            }
        } else if (ptrMOs->optPenSmthNPSel==1){
            //smoothness penalties on nonparametric selectivity curves (NOT parameter vectors)
            fPenSmoothNPSel.initialize();
            for (int i=ptrMOs->wgtPenSmthNPSel.indexmin();i<ptrMOs->wgtPenSmthNPSel.indexmax();i++){
                dvar_vector v = 1.0*pvNPSel(i);
                dvar_vector s = 1.0*npSel_cz(i)(v.indexmin(),v.indexmax());//only want to smooth over estimated part of sel function
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(s));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<"),"<<endl;
            }
            {
                int i = ptrMOs->wgtPenSmthNPSel.indexmax();
                dvar_vector v = 1.0*pvNPSel(i);
                dvar_vector s = 1.0*npSel_cz(i)(v.indexmin(),v.indexmax());
                fPenSmoothNPSel(i) = 0.5*norm2(calc2ndDiffs(s));
                objFun += penWgtSmthNPSel(i)*fPenSmoothNPSel(i);
                if (debug<0) cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgtSmthNPSel(i)<<cc<<"pen="<<fPenSmoothNPSel(i)<<cc<<"objfun="<<penWgtSmthNPSel(i)*fPenSmoothNPSel(i)<<")"<<endl;
            }
        }
        if (debug<0) cout<<tb<<tb<<") "<<"#end of smoothness penalties"<<endl;//end of smoothness penalties list
    }
    if (debug<0) cout<<tb<<")"<<cc<<"#end of nonparameteric selectivities penalties"<<endl;//end of nonparametric selectivity penalties list
    
    //penalties on sums of dev vectors to enforce sum-to-zero
    double penWgt = 0.0;
    if (current_phase()>=ptrMOs->phsSqSumDevsPen) penWgt = ptrMOs->wgtSqSumDevsPen;
    if (debug<0) rpt::echo<<"#----Check on sum(devs) in calcPenalties "<<current_phase()<<tb<<ctrProcCallsInPhase<<endl;
    if (debug<0) cout<<tb<<"devsSumSq=list("<<endl;//start of devs penalties list
    //recruitment devs
    if (ptrMPI->ptrRec->pDevsLnR->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsLnR = ";
            for (int i=1;i<=pDevsLnR.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsLnR[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsLnR=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnR,devsLnR);        
        if (debug<0) cout<<cc<<endl;
    }
    //S1 devs
    if (ptrMPI->ptrSel->pDevsS1->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS1 = ";
            for (int i=1;i<=pDevsS1.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS1[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS1=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS1,devsS1);
        if (debug<0) cout<<cc<<endl;
    }
    //S2 devs
    if (ptrMPI->ptrSel->pDevsS2->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS2 = ";
            for (int i=1;i<=pDevsS2.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS2[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS2=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS2,devsS2);
        if (debug<0) cout<<cc<<endl;
    }
    //S3 devs
    if (ptrMPI->ptrSel->pDevsS3->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS3 = ";
            for (int i=1;i<=pDevsS3.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS3[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS3=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS3,devsS3);
        if (debug<0) cout<<cc<<endl;
    }
    //S4 devs
    if (ptrMPI->ptrSel->pDevsS4->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS4 = ";
            for (int i=1;i<=pDevsS4.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS4[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS4=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS4,devsS4);
        if (debug<0) cout<<cc<<endl;
    }
    //S5 devs
    if (ptrMPI->ptrSel->pDevsS5->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS5 = ";
            for (int i=1;i<=pDevsS5.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS5[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS5=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS5,devsS5);
        if (debug<0) cout<<cc<<endl;
    }
    //S6 devs
    if (ptrMPI->ptrSel->pDevsS6->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsS6 = ";
            for (int i=1;i<=pDevsS6.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsS6[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsS6=";
        calcDevsPenalties(debug,cout,penWgt,pDevsS6,devsS6);
        if (debug<0) cout<<cc<<endl;
    }
    //capture rate devs
    if (ptrMPI->ptrFsh->pDevsLnC->getSize()){
        if (debug<0) {
            rpt::echo<<"pDevsLnC = ";
            for (int i=1;i<=pDevsLnC.indexmax();i++) rpt::echo<<"["<<i<<"]="<<sum(devsLnC[i])<<tb;
            rpt::echo<<endl;
        }
        if (debug<0) cout<<tb<<tb<<"pDevsLnC=";
        calcDevsPenalties(debug,cout,penWgt,pDevsLnC,devsLnC);
        if (debug<0) cout<<cc<<endl;
    }
    if (debug<0) rpt::echo<<"#------"<<endl<<endl;
    
    if (debug<0) cout<<tb<<"NULL)"<<endl;//end of devs penalties list
    if (debug<0) cout<<")"<<cc<<endl;//end of penalties list
    
    //Apply diminishing penalties to variances of F-devs
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
                double rmse = sqrt(value(norm2(pDevsLnC(i)))/pDevsLnC(i).size());
                rpt::echo<<tb<<i<<": pen="<<fpen<<cc<<" objfun="<<penWgt*fpen<<endl;
                if (penWgt>0) {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV="<<effCV<<cc
                            <<"pen="<<fpen<<cc<<"rmse="<<rmse<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                } else {
                    cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"effCV=Inf"    <<cc
                            <<"pen="<<fpen<<cc<<"rmse="<<rmse<<cc<<"objfun="<<penWgt*fpen<<"),"<<endl;
                }
            }
        }
        if (debug<0) cout<<tb<<tb<<"NULL)"<<endl;//end of penFDevs lists
    }
    
    if (!debug) testNaNs(value(objFun),"in calcPenalties()");
    
    if (debug>dbgObjFun) cout<<"Finished calcPenalties()"<<endl;

//-------------------------------------------------------------------------------------
//Calculate penalties on sums of devs vectors
FUNCTION void calcDevsPenalties(int debug, ostream& cout, double penWgt, param_init_bounded_vector_vector& pDevs, dvar_matrix devs)    
    dvariable fPen;
    if (debug<0) cout<<"list(";//start of list
    for (int i=pDevs.indexmin();i<=pDevs.indexmax();i++){
        if (pDevs(i).get_phase_start()){
            fPen.initialize();
            fPen = square(sum(devs(i)));
            objFun += penWgt*fPen;
            if (debug<0) {
                double rmse = sqrt(value(norm2(devs(i)))/devs(i).size());
                cout<<tb<<tb<<tb<<"'"<<i<<"'=list(wgt="<<penWgt<<cc<<"pen="<<fPen<<cc<<"objfun="<<penWgt*fPen<<cc<<"val="<<sum(devs(i))<<cc<<
                                                  "rmse="<<rmse<<"),"<<endl;
            }
        }
    }
    if (!debug) testNaNs(value(objFun),"in calcDevsPenalties()");
    if (debug<0) cout<<tb<<tb<<"NULL)";//end of penalties list
    
//-------------------------------------------------------------------------------------
//Calculate 1st differences of vector
FUNCTION dvar_vector calc1stDiffs(const dvar_vector& d)
//    cout<<"Starting calc1stDiffs"<<endl;
    RETURN_ARRAYS_INCREMENT();
//    int mn = d.indexmin();
//    int mx = d.indexmax();
//    dvar_vector cp = 1.0*d;
//    dvar_vector r(mn,mx-1);
//    r = cp(mn+1,mx).shift(mn)-cp(mn,mx-1);
    dvar_vector r = first_difference(d);
    RETURN_ARRAYS_DECREMENT();
//    cout<<"Finished calc1stDiffs"<<endl;
    return r;

//-------------------------------------------------------------------------------------
//Calculate 2nd differences of vector
FUNCTION dvar_vector calc2ndDiffs(const dvar_vector& d)
//    cout<<"Starting calc2ndDiffs"<<endl;
    RETURN_ARRAYS_INCREMENT();
//    int mn = d.indexmin();
//    int mx = d.indexmax();
//    dvar_vector r = calc1stDiffs(calc1stDiffs(d));
    dvar_vector r = first_difference(first_difference(d));
    RETURN_ARRAYS_DECREMENT();
//    cout<<"Finished calc2ndDiffs"<<endl;
    return r;

//-------------------------------------------------------------------------------------
//Calculate recruitment components in the likelihood.
FUNCTION void calcNLLs_Recruitment(int debug, ostream& cout)
    if (debug>dbgObjFun) cout<<"Starting calcNLLs_Recruitment"<<endl;
    dvariable nllRecDevs;
    if (debug<0) cout<<"list("<<endl;
    if (debug<0) cout<<tb<<"recDevs=list("<<endl;
    for (int pc=1;pc<npcRec;pc++){
        double nllWgtRecDevs = ptrMPI->ptrRec->xd(pc,1);
        nllRecDevs.initialize();
        nllRecDevs = 0.5*norm2(zscrDevsLnR_cy(pc));
        nllRecDevs += nDevsLnR_c(pc)*log(stdvDevsLnR_c(pc));
        objFun += nllWgtRecDevs*nllRecDevs;
        if (debug<0){
            double rmse = sqrt(value(norm2(devsLnR_cy(pc)))/nDevsLnR_c(pc));
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs<<cc
                    <<"objfun="<<nllWgtRecDevs*nllRecDevs<<cc
                    <<"n="<<nDevsLnR_c(pc)<<cc<<"rmse="<<rmse<<cc<<"sigmaR="<<stdvDevsLnR_c(pc)<<cc<<endl;
            cout<<"devs="; wts::writeToR(cout,value(devsLnR_cy(pc))); cout<<cc<<endl;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<")"<<cc<<endl;
        }
    }//pc
    {
        int pc = npcRec;
        double nllWgtRecDevs = ptrMPI->ptrRec->xd(pc,1);
        nllRecDevs.initialize();
        nllRecDevs = 0.5*norm2(zscrDevsLnR_cy(pc));
        nllRecDevs += nDevsLnR_c(pc)*log(stdvDevsLnR_c(pc));
        objFun += nllWgtRecDevs*nllRecDevs;
        if (debug<0){
            double rmse = sqrt(value(norm2(devsLnR_cy(pc)))/nDevsLnR_c(pc));
            cout<<tb<<tb<<"'"<<pc<<"'=list(type='normal',wgt="<<nllWgtRecDevs<<cc<<"nll="<<nllRecDevs<<cc
                    <<"objfun="<<nllWgtRecDevs*nllRecDevs<<cc
                    <<"n="<<nDevsLnR_c(pc)<<cc<<"rmse="<<rmse<<cc<<"sigmaR="<<stdvDevsLnR_c(pc)<<cc<<endl;
            cout<<"devs="; wts::writeToR(cout,value(devsLnR_cy(pc))); cout<<cc<<endl;
            cout<<"zscrs="; wts::writeToR(cout,value(zscrDevsLnR_cy(pc))); cout<<")"<<endl;
        }
    }//pc
    if (debug<0) cout<<tb<<")";//recDevs
    if (debug<0) cout<<")";
    
    if (!debug) testNaNs(value(objFun),"in calcNLLs_Recruitment()");
    
    if (debug>dbgObjFun) cout<<"Finished calcNLLs_Recruitment"<<endl;

//-------------------------------------------------------------------------------------
//Calculate objective function TODO: finish
FUNCTION void calcObjFun(int debug, ostream& cout)
    if ((debug>=dbgObjFun)||(debug<0)) {PRINT2B1("----Starting calcObjFun()")}

    int k = 0;
    dvector objFunV(0,20);

    //reset objective function
    objFun.initialize();                     objFunV[k++] = value(objFun);
    //objective function penalties
    calcPenalties(debug,cout);               objFunV[k++] = value(objFun);
    //prior likelihoods
    calcAllPriors(debug,cout);               objFunV[k++] = value(objFun);
    //recruitment component
    calcNLLs_Recruitment(debug,cout);        objFunV[k++] = value(objFun);    
    //effort extrapolation
    calcNLLs_ExtrapolatedEffort(debug,cout); objFunV[k++] = value(objFun);
    
    //data components
    calcNLLs_Fisheries(debug,cout);         objFunV[k++] = value(objFun);
    calcNLLs_Surveys(debug,cout);           objFunV[k++] = value(objFun);
    calcNLLs_GrowthData(debug,cout);        objFunV[k++] = value(objFun);
    calcNLLs_ChelaHeightData(debug,cout);   objFunV[k++] = value(objFun);
    calcNLLs_MaturityOgiveData(debug,cout); objFunV[k++] = value(objFun);
    if ((debug>=dbgObjFun)||(debug<0)){
        PRINT2B2("proc call          = ",ctrProcCalls)
        PRINT2B2("proc call in phase = ",ctrProcCallsInPhase)
        PRINT2B2("after initialization. objFun =",objFunV[0])
        PRINT2B2("after calcPenalties.  objFun =",objFunV[1])
        PRINT2B2("after calcAllPriors.  objFun =",objFunV[2])
        PRINT2B2("after calcNLLs_Recruitment.        objFun =",objFunV[3])
        PRINT2B2("after calcNLLs_ExtrapolatedEffort. objFun =",objFunV[4])
        PRINT2B2("after calcNLLs_Fisheries.          objFun =",objFunV[5])
        PRINT2B2("after calcNLLs_Surveys.            objFun =",objFunV[6])
        PRINT2B2("after calcNLLs_GrowthData.         objFun =",objFunV[7])
        PRINT2B2("after calcNLLs_ChelaHeightData.    objFun =",objFunV[8])    
        PRINT2B2("after calcNLLs_MaturityOgiveData.  objFun =",objFunV[9])    
        PRINT2B1("----Finished calcObjFun()")
    }
    
//******************************************************************************
//* Function: void calcNLLs_MaturityOgiveData
//* 
//* Description: Calculates NLLs for maturity ogive data
//* 
//* Inputs:
//*  debug - flag to print debugging info 
//* Returns:
//*  void
//* Alters:
//*  objFun
//******************************************************************************
FUNCTION void calcNLLs_MaturityOgiveData(int debug, ostream& cout)  
    //debug = dbgObjFun+1;
    if(debug>dbgObjFun) cout<<"Starting calcNLLs_MaturityOgiveData()"<<endl;
    
    if (debug<0) cout<<"list("<<endl;
    for (int i=0;i<ptrMDS->nMOD;i++){
        MaturityOgiveData* pMOD = ptrMDS->ppMOD[i];
        if ((pMOD->llWgt>0.0)||(debug<0)){
            if (debug<0) cout<<"`"<<pMOD->name<<"`=list("<<endl;
            const int nObs = pMOD->nObs;
            if (debug>dbgObjFun) cout<<"nObs= "<<nObs<<tb;
            if (nObs>0) {
                /* index of model survey corresponding to dataset */
                int v = mapD2MMOd(i+1);//
                if (debug>dbgObjFun) cout<<"v= "<<v<<tb<<"survey: "<<ptrMC->lblsSrv[v]<<endl;
                /* sex corresponding to this dataset */
                const int SEX = pMOD->sex;
                /* likelihood multiplier for this dataset */
                const double wgt = pMOD->llWgt;
                /* year corresponding to observed fractions */
                const ivector y_n = pMOD->obsYear_n;
               /* observation size bin corresponding to observation */
                const ivector obsZ_n = pMOD->obsSize_n;
               /* observation size bin index for size corresponding to observation */
                const ivector obsIZ_n = pMOD->obsZBI_n;
                /* sample sizes for observed fractions */
                const dvector ss_n = pMOD->obsSS_n;
                /* observed fractions of new shell mature crab at size */
                const dvector obsPM_n = pMOD->obsPrMat_n;
                
                /* calculate model-predicted ratios by year for 
                 * all observed size bins */
                int mny = n_vyxmsz(v).indexmin();
                int mxy = n_vyxmsz(v).indexmax();
                dvar_matrix modPrMat_yzp(mny,mxy,1,pMOD->nZBs); modPrMat_yzp.initialize();
                dvar_vector vmMat_z(1,nZBs);
                dvar_vector vmTot_z(1,nZBs);
                for (int y=mny;y<=mxy;y++){
                    vmMat_z  = n_vyxmsz(v,y,SEX,MATURE,NEW_SHELL);
                    vmTot_z  = vmMat_z + n_vyxmsz(v,y,SEX,IMMATURE,NEW_SHELL);
                    modPrMat_yzp(y) = elem_div(pMOD->zbRemapper*vmMat_z,
                                               pMOD->zbRemapper*vmTot_z+1.0e-10);//small value added
                    if (debug>dbgObjFun) cout<<y<<tb<<modPrMat_yzp(y)<<endl;
                }
                    
                //compare model predictions with observations
                dvar_vector modPM_n(1,nObs);  modPM_n.initialize();
                dvar_vector nlls_n(1,nObs);   nlls_n.initialize();
                dvector zscrs_n(1,nObs);      zscrs_n.initialize();
                if (debug>dbgObjFun) cout<<"n"<<tb<<"y"<<tb<<"z"<<tb<<"iz"<<tb<<"ss"<<tb<<"obsPM"<<tb<<"modPM"<<tb<<"nll"<<tb<<"zscr"<<endl;
                for (int n=1;n<=nObs;n++){
                    if (debug>dbgObjFun) cout<<n<<tb<<y_n(n)<<tb<<obsZ_n(n)<<tb<<obsIZ_n(n)<<tb<<ss_n(n)<<tb<<obsPM_n(n)<<tb;
//                    if ((obsIZ(n)>0)&(obsIZ(n)<=nZBs)){
//                        if (debug>dbgObjFun) cout<<y_n(n)<<tb;
                        modPM_n(n) = modPrMat_yzp(y_n(n),obsIZ_n(n));
                        if (debug>dbgObjFun) cout<<modPM_n(n)<<tb;
                        if ((modPM_n(n)>0.0)&(modPM_n(n)<1.0)){
                            if (obsPM_n(n)>0.0) nlls_n(n) -= ss_n(n)*obsPM_n(n)*(log(modPM_n(n))-log(obsPM_n(n)));
                            if (obsPM_n(n)<1.0) nlls_n(n) -= ss_n(n)*(1.0-obsPM_n(n))*(log(1.0-modPM_n(n))-log(1.0-obsPM_n(n)));
                            if (debug>dbgObjFun) cout<<nlls_n(n)<<tb;
                            double modPMv = value(modPM_n(n));
                            zscrs_n(n) = (obsPM_n(n)-modPMv)/sqrt(modPMv*(1.0-modPMv)/ss_n(n));
                            if (debug>dbgObjFun) cout<<zscrs_n(n)<<endl;
                        }
//                    }
                }//loop over n
                dvariable nll = sum(nlls_n);
                objFun += wgt*nll;
                if (debug>dbgObjFun) cout<<tb<<"sum(nlls) = "<<nll<<endl;
                if (debug<0) {
                    cout<<"fleet="<<qt<<ptrMC->lblsSrv[v]<<qt<<cc;
                    cout<<"sex="<<qt<<tcsam::getSexType(SEX)<<qt<<cc<<"type='binomial'"<<cc
                            <<"wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl;
                    cout<<"y=";     wts::writeToR(cout,y_n);             cout<<cc<<endl;
                    cout<<"n=";     wts::writeToR(cout,ss_n);            cout<<cc<<endl;
                    cout<<"z=";     wts::writeToR(cout,obsZ_n);          cout<<cc<<endl;
                    cout<<"i=";     wts::writeToR(cout,obsIZ_n);         cout<<cc<<endl;
                    cout<<"obsPM="; wts::writeToR(cout,obsPM_n);         cout<<cc<<endl;
                    cout<<"modPM="; wts::writeToR(cout,value(modPM_n));  cout<<cc<<endl;
                    cout<<"nlls=";  wts::writeToR(cout,value(nlls_n));   cout<<cc<<endl;
                    cout<<"zscrs="; wts::writeToR(cout,zscrs_n);         cout<<cc<<endl;
                    cout<<"rmse="<<sqrt(norm2(zscrs_n)/zscrs_n.size())<<"),"<<endl;
                }
            }//nObs>0
        }//((pMOD->llWgt>0)||(debug<0))
    }//datasets (i)
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>dbgObjFun) cout<<"finished calcNLLs_MaturityOgiveData()"<<endl;

///******************************************************************************
//* Function: void calcNLLs_ChelaHeightData
//* 
//* Description: Calculates NLLs for chela height data
//* 
//* Inputs:
//*  debug - flag to print debugging info 
//* Returns:
//*  void
//* Alters:
//*  objFun
//******************************************************************************
FUNCTION void calcNLLs_ChelaHeightData(int debug, ostream& cout)  
    if(debug>dbgObjFun) cout<<"Starting calcNLLs_ChelaHeightData()"<<endl;
    
    if (debug<0) cout<<"list("<<endl;
    for (int i=0;i<ptrMDS->nCHD;i++){
        ChelaHeightData* pCHD = ptrMDS->ppCHD[i];
        if ((pCHD->llWgt>0.0)||(debug<0)){
            if (debug<0) cout<<"`"<<pCHD->name<<"`=list("<<endl;
            int nObs = pCHD->nObs;
            //cout<<"nObs= "<<nObs<<tb;
            if (nObs>0) {
                /* index of model survey corresponding to dataset */
                int v = mapD2MChd(i+1);//
                //cout<<"v= "<<v<<endl;
                /* likelihood multiplier for this dataset */
                double wgt = pCHD->llWgt;
                /* year corresponding to observed fractions */
                ivector y_n = pCHD->obsYear_n;
                /* sample sizes from observations */
                dvector ss_n = pCHD->obsSS_n;
               /* model size bin index for size corresponding to observation */
                ivector obsIZ = pCHD->obsSizeBinIndex_n;
                /* observed fractions of new shell mature crab at size */
                dvector obsPM = 1.0*pCHD->obsPrMat_n;

                dvar_vector modPM(1,nObs);  modPM.initialize();
                dvar_vector nlls_n(1,nObs); nlls_n.initialize();
                dvector zscrs_n(1,nObs);    zscrs_n.initialize();
                for (int n=1;n<=nObs;n++){
                    //cout<<"n="<<n<<tb<<"obsIZ="<<tb<<obsIZ(n)<<tb;
                    if ((obsIZ(n)>0)&(obsIZ(n)<=nZBs)){
                        //cout<<y_n(n)<<tb;
                        dvariable nMat = n_vyxmsz(v,y_n(n),MALE,MATURE,NEW_SHELL,obsIZ(n));
                        dvariable nTot = nMat + n_vyxmsz(v,y_n(n),MALE,IMMATURE,NEW_SHELL,obsIZ(n));
                        modPM(n) = nMat/nTot;
                        //cout<<nMat<<tb<<nTot<<tb<<modPM(n)<<tb;
                        if ((modPM(n)>0.0)&(modPM(n)<1.0)){
                            if (obsPM(n)>0.0) nlls_n(n) -= ss_n(n)*obsPM(n)*(log(modPM(n))-log(obsPM(n)));
                            if (obsPM(n)<1.0) nlls_n(n) -= ss_n(n)*(1.0-obsPM(n))*(log(1.0-modPM(n))-log(1.0-obsPM(n)));
                            //cout<<nlls_n(n)<<tb;
                            double modPMv = value(modPM(n));
                            zscrs_n(n) = (obsPM(n)-modPMv)/sqrt(modPMv*(1.0-modPMv)/ss_n(n));
                            //cout<<"zscrs_n(n)="<<zscrs_n(n)<<endl;
                        }
                    }
                }//loop over n
                dvariable nll = sum(nlls_n);
                objFun += wgt*nll;
                if (debug<0) {
                    cout<<"fleet="<<qt<<ptrMC->lblsSrv[v]<<qt<<cc;
                    cout<<"sex="<<qt<<"MALE"<<qt<<cc<<"type='binomial'"<<cc
                            <<"wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl;
                    cout<<"y=";     wts::writeToR(cout,y_n);             cout<<cc<<endl;
                    cout<<"n=";     wts::writeToR(cout,ss_n);            cout<<cc<<endl;
                    cout<<"z=";     wts::writeToR(cout,pCHD->obsSize_n); cout<<cc<<endl;
                    cout<<"i=";     wts::writeToR(cout,obsIZ);           cout<<cc<<endl;
                    cout<<"obsPM="; wts::writeToR(cout,obsPM);           cout<<cc<<endl;
                    cout<<"modPM="; wts::writeToR(cout,value(modPM));    cout<<cc<<endl;
                    cout<<"nlls=";  wts::writeToR(cout,value(nlls_n));   cout<<cc<<endl;
                    cout<<"zscrs="; wts::writeToR(cout,zscrs_n);         cout<<cc<<endl;
                    cout<<"rmse="<<sqrt(norm2(zscrs_n)/zscrs_n.size())<<"),"<<endl;
                }
            }//nObs>0
        }//((pCHD->llWgt>0)||(debug<0))
    }//datasets (i)
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>dbgObjFun) cout<<"finished calcNLLs_ChelaHeightData()"<<endl;

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
        if ((pGD->llWgt>0.0)||(debug<0)){
            if (debug<0) cout<<"`"<<ptrMDS->ppGrw[i]->name<<"`=list("<<endl;
            for (int x=1;x<=nSXs;x++){
                int nObs = pGD->nObs_x(x);
                if (nObs>0) {
                    double wgt = pGD->llWgt;
                    /* observation year */
                    ivector year_n = pGD->obsYears_xn(x);
                    /* pre-molt size, by observation */
                    dvar_vector zpre_n = pGD->inpData_xcn(x,2);
                    /* post-molt size, by observation */
                    dvar_vector zpst_n = ptrMDS->ppGrw[i]->inpData_xcn(x,3);
                    /* molt increment, by observation */
                    dvar_vector incZ_n = zpst_n - zpre_n;
                    /* mean post-molt size, by observation */
                    //dvar_vector mnZ_n = elem_prod(mfexp(grA_xy(x)(year_n)),pow(zpre_n,grB_xy(x)(year_n)));
                    dvar_vector mnZ_n(zpre_n.indexmin(),zpre_n.indexmax());
                    if (ptrMOs->optGrowthParam==0){
                        mnZ_n = mfexp(grA_xy(x)(year_n)+elem_prod(grB_xy(x)(year_n),log(zpre_n)));
                    } else if (ptrMOs->optGrowthParam==1){
                        mnZ_n = elem_prod(
                                    grA_xy(x)(year_n),
                                    mfexp(
                                        elem_prod(
                                            elem_div(log(elem_div( grB_xy(x)(year_n), grA_xy(x)(year_n))),
                                                     log(elem_div(zGrB_xy(x)(year_n),zGrA_xy(x)(year_n)))),
                                            log(elem_div(zpre_n,zGrA_xy(x)(year_n)))
                                        )
                                    )
                                );
                    } else if (ptrMOs->optGrowthParam==2){
                        mnZ_n = elem_prod(
                                    grA_xy(x)(year_n),
                                    mfexp(
                                        elem_prod(
                                            grB_xy(x)(year_n),
                                            log(elem_div(zpre_n,zGrA_xy(x)(year_n)))
                                        )
                                    )
                                );
                    } else {
                        //throw error
                        PRINT2B1(" ")
                        PRINT2B1("#---------------------")
                        PRINT2B2("Unknown growth parameterization option",ptrMOs->optGrowthParam)
                        PRINT2B1("Terminating model run. Please correct.")
                        ad_exit(-1);
                    }
                    /* multiplicative scale factor, by observation */
                    dvar_vector ibeta_n = 1.0/grBeta_xy(x)(year_n);
                    /* location factor, by observation */
                    dvar_vector alpha_n = elem_prod(mnZ_n-zpre_n,ibeta_n);
                    dvar_vector nlls_n(1,nObs); nlls_n.initialize();
                    nlls_n = -wts::log_gamma_density(incZ_n,alpha_n,ibeta_n);
                    dvariable nll = sum(nlls_n);
                    if (isnan(value(nll))){
                        dvar_vector zscrs = elem_div((zpst_n-mnZ_n),sqrt(elem_prod(mnZ_n,grBeta_xy(x)(year_n))));
                        ofstream os("GrowthData.NLLs.NanReport.dat");
                        os.precision(12);
                        os<<"phase = "<<current_phase()<<endl;
                        os<<"sex   = "<<tcsam::getSexType(x)<<endl;
                        os<<"nll   = "<<nll<<endl;
                        os<<"nObs  = "<<nObs<<endl;
                        os<<"year  grA   zGrA    grB    zGrB   zpre_n  zpst_n  mnZ_n   incZ   mnInc  ibeta_n alpha_n nll_n  zscr"<<endl;
                        for (int n=1;n<=nObs;n++){
                            os<<year_n(n)<<tb<<
                                    grA_xy(x)(year_n(n))<<tb<<grB_xy(x)(year_n(n))<<tb<<
                                    zGrA_xy(x)(year_n(n))<<tb<<zGrB_xy(x)(year_n(n))<<tb<<
                                    zpre_n(n)<<tb<<zpst_n(n)<<tb<<mnZ_n(n)<<tb<<
                                    zpst_n(n)-zpre_n(n)<<tb<<mnZ_n(n)-zpre_n(n)<<tb<<
                                    ibeta_n(n)<<tb<<alpha_n(n)<<tb<<nlls_n(n)<<tb<<zscrs(n)<<endl;
                        }
                        os.close();
                        exit(-1);
                    }
                    objFun += wgt*nll;
                    if (debug<0) {
                        dvar_vector zscrs = elem_div((zpst_n-mnZ_n),sqrt(elem_prod(mnZ_n,grBeta_xy(x)(year_n))));
                        cout<<tcsam::getSexType(x)<<"=list(type='gamma',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl;
                        cout<<"years="; wts::writeToR(cout,year_n);         cout<<cc<<endl;
                        cout<<"zPre=";  wts::writeToR(cout,value(zpre_n));  cout<<cc<<endl;
                        cout<<"zPst=";  wts::writeToR(cout,value(zpst_n));  cout<<cc<<endl;
                        cout<<"grA=";   wts::writeToR(cout,value(grA_xy(x)(year_n))); cout<<cc<<endl;
                        cout<<"grB=";   wts::writeToR(cout,value(grB_xy(x)(year_n))); cout<<cc<<endl;
                        cout<<"mnZ =";  wts::writeToR(cout,value(mnZ_n));   cout<<cc<<endl;
                        cout<<"ibeta="; wts::writeToR(cout,value(ibeta_n)); cout<<cc<<endl;
                        cout<<"alpha="; wts::writeToR(cout,value(alpha_n)); cout<<cc<<endl;
                        cout<<"nlls=";  wts::writeToR(cout,value(nlls_n));  cout<<cc<<endl;
                        cout<<"zscrs="; wts::writeToR(cout,value(zscrs));   cout<<cc<<endl;
                        cout<<"rmse="<<sqrt(value(norm2(zscrs))/zscrs.size())<<"),"<<endl;
                    }
                }//nObs>0
            }//x
            if (debug<0) cout<<"NULL),";
        }//(pGD->llWgt>0.0)||(debug<0)
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
        os.precision(12);
        os<<"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        os<<"#----NaN detected: "<<str<<"---"<<endl;
        os<<"#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        updateMPI(0,cout);
        os<<"nanRep=list(phase="<<current_phase()<<",pcCtr="<<ctrProcCalls<<cc<<",pcCtrInPhs="<<ctrProcCallsInPhase<<cc<<endl;
        ReportToR_Params(os,0,cout);         os<<","<<endl;
        ReportToR_ModelProcesses(os,0,cout); os<<","<<endl;
        ReportToR_ModelResults(os,0,cout);   os<<","<<endl;
        ReportToR_ModelFits(os,-1.0,-1,cout);     os<<endl;
        os<<")"<<endl;
        os.close();
        exit(-1);
    }

//-------------------------------------------------------------------------------------
//Calculate norm2 NLL contribution to objective function
FUNCTION void calcNorm2NLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcNorm2NLL()"<<endl;
    int y;
    double rmse = 0.0; int cnt = 0;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    for (int i=1;i<=yrs.size();i++){
        y = yrs(i);
        if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
            zscr(y) = (obs[i]-mod[y]); cnt++;
        }
    }
    if (cnt>0) rmse = sqrt(value(norm2(zscr))/cnt);
    nll += 0.5*norm2(zscr);
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='norm2',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<cc<<endl;
        cout<<"rmse="<<rmse<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcNorm2NLL()"<<endl;
    
//-------------------------------------------------------------------------------------
//Calculate normal NLL contribution to objective function
FUNCTION void calcNormalNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
   if (debug>=dbgAll) cout<<"Starting calcNormalNLL()"<<endl;
    int y;
    double rmse = 0.0; int cnt = 0;
    dvariable nll = 0.0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (obs[i]-mod[y])/stdv[i]; cnt++;
            }
        }
        if (cnt>0) rmse = sqrt(value(norm2(zscr))/cnt);
        nll += 0.5*norm2(zscr);
    }
    objFun += wgt*nll;
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='normal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<cc<<endl;
        cout<<"rmse="<<rmse<<")";
    }
   if (debug>=dbgAll) cout<<"Finished calcNormalNLL()"<<endl;
    
//-------------------------------------------------------------------------------------
//Calculate lognormal NLL contribution to objective function
FUNCTION void calcLognormalNLL(double wgt, dvar_vector& mod, dvector& obs, dvector& stdv, ivector& yrs, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcLognormalNLL()"<<endl;
    int y;
    dvariable nll = 0.0;
    double rmse = 0.0; int cnt = 0;
    dvar_vector zscr(mod.indexmin(),mod.indexmax());
    zscr.initialize();
    if (sum(stdv)>0){
        for (int i=1;i<=yrs.size();i++){
            y = yrs(i);
            if ((zscr.indexmin()<=y)&&(y<=zscr.indexmax())) {
                zscr(y) = (log(obs[i]+smlVal)-log(mod[y]+smlVal))/stdv[i]; cnt++;
            }
        }
        if (cnt>0) {
            rmse = sqrt(value(norm2(zscr))/cnt);
            nll = 0.5*norm2(zscr);
            objFun += wgt*nll;
        }
    }
    if (debug<0){
        adstring obsyrs = wts::to_qcsv(yrs);
        adstring modyrs = str(mod.indexmin())+":"+str(mod.indexmax());
        cout<<"list(nll.type='lognormal',wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<endl; 
        cout<<"obs=";   wts::writeToR(cout,obs,        obsyrs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,value(mod), modyrs); cout<<cc<<endl;
        cout<<"stdv=";  wts::writeToR(cout,stdv,       obsyrs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,value(zscr),modyrs); cout<<cc<<endl;
        cout<<"rmse="<<rmse<<")";
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
FUNCTION void calcMultinomialNLL(double wgt, dvar_vector& mod, dvector& obs, double& ss, int& yr, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcMultinomialNLL()"<<endl;
    dvariable nll = -ss*(obs*(log(mod+smlVal)-log(obs+smlVal)));//note dot-product sums
    objFun += wgt*nll;
    if (debug<0){
        dvector vmod = value(mod);
        dvector nlls = -ss*(elem_prod(obs,log(vmod+smlVal)-log(obs+smlVal)));
        dvector zscrs = elem_div(obs-vmod,sqrt(elem_prod((vmod+smlVal),1.0-(vmod+smlVal))/ss));//pearson residuals
        double effN = 0.0;
        if ((ss>0)&&(norm2(obs-vmod)>0)) effN = (vmod*(1.0-vmod))/norm2(obs-vmod);
        cout<<"list(nll.type='multinomial',yr="<<yr<<cc<<"wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<"ss="<<ss<<cc<<"effN="<<effN<<cc<<endl; 
        adstring dzbs = "size=c("+ptrMC->csvZBs+")";
        cout<<"nlls=";  wts::writeToR(cout,nlls, dzbs); cout<<cc<<endl;
        cout<<"obs=";   wts::writeToR(cout,obs,  dzbs); cout<<cc<<endl;
        cout<<"mod=";   wts::writeToR(cout,vmod, dzbs); cout<<cc<<endl;
        cout<<"zscrs="; wts::writeToR(cout,zscrs,dzbs); cout<<endl;
        cout<<")";
    }
    if (debug>=dbgAll) cout<<"Finished calcMultinomialNLL()"<<endl;
 
//-------------------------------------------------------------------------------------
////Calculate Dirichlet NLL contribution to objective function
//FUNCTION void calcDirichletNLL(double wgt, dvar_vector& mod, dvector& obs, double& ss, int& yr, int debug, ostream& cout)
//    if (debug>=dbgAll) cout<<"Starting calcDirichletNLL()"<<endl;
//    dvariable nll = -ss*(obs*(log(mod+smlVal)-log(obs+smlVal)));//note dot-product sums
//    objFun += wgt*nll;
//    
//    dvector o = obs/sum(obs);
//    int c1 = o.indexmin();
//    int c2 = o.indexmax();
//    dvar_vector p = mod;
//    dvariable vn = ???;
//    
//    dvariable lmnB = 0.0;
//    dvariable sj = 0.0;
//    dvariable alpha0 = 0.0;
//    dvar_vector alpha = vn * p/sum(p);
//    for ( int j = c1; j <= c2; j++ ) {
//        aj = alpha(j);
//        alpha0 += aj;
//        lmnB += gammln(aj);
//        sj += (aj - 1.0) * log(1e-10 + obs(j));
//    }
//    lmnB -= gammln(alpha0);
//    nll -= sj - lmnB;
//    
//    if (debug<0){
//        dvector vmod = value(mod);
//        dvector nlls = -ss*(elem_prod(obs,log(vmod+smlVal)-log(obs+smlVal)));
//        dvector zscrs = elem_div(obs-vmod,sqrt(elem_prod((vmod+smlVal),1.0-(vmod+smlVal))/ss));//pearson residuals
//        double effN = 0.0;
//        if (ss>0) effN = (vmod*(1.0-vmod))/norm2(obs-vmod);
//        cout<<"list(nll.type='Dirichlet',yr="<<yr<<cc<<"wgt="<<wgt<<cc<<"nll="<<nll<<cc<<"objfun="<<wgt*nll<<cc<<"ss="<<ss<<cc<<"effN="<<effN<<cc<<endl; 
//        adstring dzbs = "size=c("+ptrMC->csvZBs+")";
//        cout<<"nlls=";  wts::writeToR(cout,nlls, dzbs); cout<<cc<<endl;
//        cout<<"obs=";   wts::writeToR(cout,obs,  dzbs); cout<<cc<<endl;
//        cout<<"mod=";   wts::writeToR(cout,vmod, dzbs); cout<<cc<<endl;
//        cout<<"zscrs="; wts::writeToR(cout,zscrs,dzbs); cout<<endl;
//        cout<<")";
//    }
//    if (debug>=dbgAll) cout<<"Finished calcDirchletNLL()"<<endl;
 
//-------------------------------------------------------------------------------------
//Pass through for no NLL contribution to objective function for size comps
FUNCTION void calcNoneNLL(double wgt, dvar_vector& mod, dvector& obs, double& ss, int& yr, int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcNoneNLL(size comps)"<<endl;
    if (debug<0){
        cout<<"list(nll.type='none',yr="<<yr<<cc<<"wgt="<<wgt<<cc<<"nll="<<0.0<<cc<<"objfun="<<0.0<<cc<<"ss="<<0<<cc<<"effN="<<0<<")";
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
FUNCTION void calcNLL(int llType, double wgt, dvar_vector& mod, dvector& obs, double& ss, int& yr, int debug, ostream& cout)
    switch (llType){
        case tcsam::LL_NONE:
            calcNoneNLL(wgt,mod,obs,ss,yr,debug,cout);
            break;
        case tcsam::LL_MULTINOMIAL:
            calcMultinomialNLL(wgt,mod,obs,ss,yr,debug,cout);
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
FUNCTION d5_array calcNLLs_CatchNatZ(SizeFrequencyData* ptrZFD, dvar5_array& mA_yxmsz, int debug, ostream& cout)
    if (debug>=dbgNLLs) cout<<"Starting calcNLLs_CatchNatZ()"<<endl;
    d5_array effWgtComps_xmsyn;
    if (ptrZFD->optFit==tcsam::FIT_NONE) return effWgtComps_xmsyn;
    ivector yrs = ptrZFD->yrs;
    int y;
    double ss;
    dvariable nT;
    int mny = mA_yxmsz.indexmin();
    int mxy = mA_yxmsz.indexmax();//may NOT be mxYr
    dvector     oP_z;//observed size comp.
    dvar_vector mP_z;//model size comp.
    if (debug<0) cout<<"list("<<endl;
    if (ptrZFD->optFit==tcsam::FIT_BY_TOT){
        effWgtComps_xmsyn.allocate(tcsam::ALL_SXs,tcsam::ALL_SXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                    calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                    effWgtComps_xmsyn(tcsam::ALL_SXs,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    if (debug<0) cout<<")"<<cc<<endl;
                }//value(nT)>0
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_TOT
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                        calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XM){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               1,nMSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                            effWgtComps_xmsyn(x,m,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//m
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XM
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XS){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               1,nSCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                            effWgtComps_xmsyn(x,tcsam::ALL_MSs,s,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//nT>0
                    }//s
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XS
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XMS){
        effWgtComps_xmsyn.allocate(1,nSXs,1,nMSs,1,nSCs,mny,mxy,1,3);
        oP_z.allocate(1,nZBs);
        mP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mP_z,oP_z,ss,y,debug,cout);
                                effWgtComps_xmsyn(x,m,s,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//nT>0
                        }//s
                    }//m
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XMS
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XE){
        effWgtComps_xmsyn.allocate(tcsam::ALL_SXs,tcsam::ALL_SXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSXs*nZBs);
        mP_z.allocate(1,nSXs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                        calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                        if (debug<0) cout<<")"<<cc<<endl;
                    }//x
                    effWgtComps_xmsyn(tcsam::ALL_SXs,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                }//nT>0
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XE
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X_ME){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nMSs*nZBs);
        mP_z.allocate(1,nMSs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X_ME
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_X_SE){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//s
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X_SE
    } else
    if (ptrZFD->optFit==tcsam::FIT_BY_XME){
        effWgtComps_xmsyn.allocate(tcsam::ALL_SXs,tcsam::ALL_SXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSXs*nMSs*nZBs);
        mP_z.allocate(1,nSXs*nMSs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                            calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                            if (debug<0) cout<<")"<<cc<<endl;
                        }//m
                    }//x
                    effWgtComps_xmsyn(tcsam::ALL_SXs,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                }//nT>0
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XME
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_XM_SE){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               1,nMSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nSCs*nZBs);
        mP_z.allocate(1,nSCs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
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
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//s
                            effWgtComps_xmsyn(x,m,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                        }//nT>0
                    }//m
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_XM_SE
    } else 
    if (ptrZFD->optFit==tcsam::FIT_BY_X_MSE){
        effWgtComps_xmsyn.allocate(1,nSXs,
                               tcsam::ALL_MSs,tcsam::ALL_MSs,
                               tcsam::ALL_SCs,tcsam::ALL_SCs,
                               mny,mxy,1,3);
        oP_z.allocate(1,nMSs*nSCs*nZBs);
        mP_z.allocate(1,nMSs*nSCs*nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug>=dbgNLLs) cout<<"y = "<<y<<endl;
            if ((mny<=y)&&(y<=mxy)) {
                for (int x=1;x<=nSXs;x++) {
                    ss = 0;
                    nT = sum(mA_yxmsz(y,x));//=0 if not calculated
                    if (value(nT)>0){
                        oP_z.initialize();//observed size comp.
                        mP_z.initialize();//model size comp.
                        for (int m=1;m<=nMSs;m++) {
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs+(m-1)*nSCs*nZBs;
                                int mxz =       s*nZBs+(m-1)*nSCs*nZBs;
                                ss += ptrZFD->ss_xmsy(x,m,s,iy);
                                oP_z(mnz,mxz).shift(1) += ptrZFD->PatZ_xmsyz(x,m,s,iy);
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
                            for (int s=1;s<=nSCs;s++) {
                                int mnz = 1+(s-1)*nZBs+(m-1)*nSCs*nZBs;
                                int mxz =       s*nZBs+(m-1)*nSCs*nZBs;
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
                                calcNLL(ptrZFD->llType,ptrZFD->llWgt,mPt,oPt,ss,y,debug,cout);
                                if (debug<0) cout<<")"<<cc<<endl;
                            }//s
                        }//m
                        effWgtComps_xmsyn(x,tcsam::ALL_MSs,tcsam::ALL_SCs,y) = calcEffWgtComponents(ss, oP_z, mP_z, debug, cout);
                    }//nT>0
                }//x
            } //if ((mny<=y)&&(y<=mxy))
        } //loop over iy
        //FIT_BY_X_MSE
    } else 
    {
        std::cout<<"Calling calcNLLs_CatchNatZ with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(ptrZFD->optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug<0) cout<<"NULL)";
    if (debug>=dbgNLLs) cout<<"Finished calcNLLs_CatchNatZ()"<<endl;
    return effWgtComps_xmsyn;
    
//-------------------------------------------------------------------------------------
//Calculate effective weight components for a size composition 
FUNCTION dvector calcEffWgtComponents(double ss, dvector& obs, dvar_vector& mod, int debug, ostream& cout)
    //if (debug) cout<<"Starting calcEffWgtComponents(...)"<<endl;
    dvector effWgts(1,3); effWgts = 0.0;
    dvector vmd = value(mod);
    //set counter for valid comps
    effWgts(1) = 1.0;
    //McAllister-Ianelli E_y/N_y ala Punt, 2017
    effWgts(2) = ((vmd*(1.0-vmd))/norm2(obs-vmd))/ss;
    //Francis weights ala Punt, 2017
    int N = obs.size();
    dvector L(1,N); L.fill_seqadd(1.0,1.0);
    double obsMnL = L*obs;
    double prdMnL = L*vmd;
    double stdMnL = sqrt(vmd*elem_prod(L-prdMnL,L-prdMnL)/N);
    effWgts(3) = (obsMnL-prdMnL)/stdMnL;
    //if (debug) cout<<"Starting calcEffWgtComponents(...)"<<endl;
    return effWgts;
    
//-------------------------------------------------------------------------------------
//Calculate fishery components to objective function
FUNCTION void calcNLLs_Fisheries(int debug, ostream& cout)
//    if (debug>0) debug = dbgAll+10;
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Fisheries()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int f=1;f<=nFsh;f++){
        if (debug>=dbgAll) cout<<"calculating NLLs for fishery "<<ptrMC->lblsFsh[f]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsFsh[f]<<"`=list("<<endl;
        int fd = mapM2DFsh(f);//get index for fishery data corresponding to model fishery f
        FleetData* ptrObs = ptrMDS->ppFsh[fd-1];
        if (ptrObs->hasRCD){//retained catch data
            if (debug<0) cout<<"retained.catch=list("<<endl;
            if (ptrObs->ptrRCD->hasN && ptrObs->ptrRCD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---retained catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrN,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasB && ptrObs->ptrRCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---retained catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrRCD->ptrB,rmN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---retained catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
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
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrN,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasB && ptrObs->ptrTCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---total catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrTCD->ptrB,cpN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---total catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
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
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrN,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasB && ptrObs->ptrDCD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---discard catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.0001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrDCD->ptrB,dsN_fyxmsz(f),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---discard catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
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
//Calculate extrapolated effort components for the objective function
FUNCTION void calcNLLs_ExtrapolatedEffort(int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"Starting calcNLLs_ExtrapolatedEffort()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    CapRateAvgScenarios* ptrCRASs = ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios;
    double stdv = sqrt(log(1.0+square(0.1)));//TODO: check this!!
    int mnx, mxx, mnm, mxm, mns, mxs;
    nllEffX_nxms.initialize();
    zscrEffX_nxmsy.initialize();
    for (int n=1;n<=ptrCRASs->nAvgs;n++){
        CapRateAvgScenario* ptrCRAS = ptrCRASs->ppCRASs[n-1];
        int idEAS = ptrCRAS->idEffAvgInfo;//index to associated average effort
        int f     = ptrCRAS->f;//fishery
        mnx = mxx = ptrCRAS->x;//sex index
        if (mnx==tcsam::ALL_SXs) {mnx=1; mxx=tcsam::nSXs;}
        mnm = mxm = ptrCRAS->m;//maturity index
        if (mnm==tcsam::ALL_MSs) {mnm=1; mxm=tcsam::nMSs;}
        mns = mxs = ptrCRAS->s;//shell index
        if (mns==tcsam::ALL_SCs) {mns=1; mxs=tcsam::nSCs;}
        for (int x=mnx;x<=mxx;x++){
            for (int m=mnm;m<=mxm;m++){
                for (int s=mns;s<=mxs;s++){
                    zscrEffX_nxmsy(n,x,m,s) = (log(obsEff_nxmsy(n,x,m,s)+smlVal)-log(prdEff_nxmsy(n,x,m,s)+smlVal))/stdv;
                    nllEffX_nxms(n,x,m,s)   = norm2(zscrEffX_nxmsy(n,x,m,s));
                    if (debug>dbgAll){
                        cout<<"n    x   m    s      stdv   wgt     nll"<<endl;
                        cout<<n<<tb<<x<<tb<<m<<tb<<s<<tb<<stdv<<tb<<ptrCRAS->llWgt<<tb<<nllEffX_nxms(n,x,m,s)<<endl;
                        cout<<"obsEff_nxmsy(n,x,m,s) = "<<obsEff_nxmsy(n,x,m,s)<<endl;
                        cout<<"prdEff_nxmsy(n,x,m,s) = "<<prdEff_nxmsy(n,x,m,s)<<endl;
                        cout<<"zscores               = "<<zscrEffX_nxmsy(n,x,m,s)<<endl;
                    }
                    if ((ptrCRAS->idParam)&&(active(pLnEffX(ptrCRAS->idParam)))){
                        objFun += (ptrCRAS->llWgt)*nllEffX_nxms(n,x,m,s);
                    } else {
                        objFun += 0.0;//add nothing, for now
                    }
                }//s
            }//m
        }//x
        if (debug<0){
            cout<<"`"<<n<<"`=list(n="<<n<<",f='"<<ptrMC->lblsFsh(f)<<"',stdv="<<stdv<<",wgt="<<ptrCRAS->llWgt<<cc<<endl;
            cout<<"obsFc=";   wts::writeToR(cout,obsFc_nxmsy(n),   xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"prdFc=";   wts::writeToR(cout,prdFc_nxmsy(n),   xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"obsEff=";  wts::writeToR(cout,obsEff_nxmsy(n),  xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"prdEff=";  wts::writeToR(cout,prdEff_nxmsy(n),  xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"zscores="; wts::writeToR(cout,zscrEffX_nxmsy(n),xDms,mDms,sDms,yDms);cout<<cc<<endl;
            cout<<"nlls=";    wts::writeToR(cout,nllEffX_nxms(n),  xDms,mDms,sDms);cout<<endl;
            cout<<")"<<cc<<endl;
        }
    }//n
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug>=dbgAll) cout<<"Finished calcNLLs_ExtrapolatedEffort()"<<endl;
        
//-------------------------------------------------------------------------------------
//Calculate survey components to objective function
FUNCTION void calcNLLs_Surveys(int debug, ostream& cout)
//    if (debug>0) debug = dbgAll+10;
    if (debug>=dbgAll) cout<<"Starting calcNLLs_Surveys()"<<endl;
    if (debug<0) cout<<"list("<<endl;
    for (int v=1;v<=nSrv;v++){
        if (debug>=dbgAll) cout<<"calculating NLLs for survey "<<ptrMC->lblsSrv[v]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsSrv[v]<<"`=list("<<endl;
        FleetData* ptrObs = ptrMDS->ppSrv[v-1];
        if (ptrObs->hasICD){//index catch data
            if (debug<0) cout<<"index.catch=list("<<endl;
            if (ptrObs->ptrICD->hasN && ptrObs->ptrICD->ptrN->optFit){
                if (debug>=dbgAll) cout<<"---index catch abundance"<<endl;
                if (debug<0) cout<<"abundance="<<endl;
                smlVal = 0.000001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrN,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasB && ptrObs->ptrICD->ptrB->optFit){
                if (debug>=dbgAll) cout<<"---index catch biomass"<<endl;
                if (debug<0) cout<<"biomass="<<endl;
                smlVal = 0.000001;//match to TCSAM2013
                calcNLLs_AggregateCatch(ptrObs->ptrICD->ptrB,n_vyxmsz(v),debug,cout);
                if (debug<0) cout<<","<<endl;
            }
            if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                if (debug>=dbgAll) cout<<"---index catch size frequencies"<<endl;
                if (debug<0) cout<<"n.at.z="<<endl;
                smlVal = 0.001;//match to TCSAM2013
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
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pLnR,pLnR,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRCV,pRCV,debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRX,pRX,  debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRa, pRa, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pRb, pRb, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrRec->pDevsLnR,devsLnR,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
   
    //natural mortality parameters
    if (debug<0) cout<<tb<<"'natural mortality'=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pM,   pM,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM1, pDM1, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM2, pDM2, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM3, pDM3, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrNM->pDM4, pDM4, debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //growth parameters
    if (debug<0) cout<<tb<<"growth=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pGrA,   pGrA,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pGrB,   pGrB,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrGrw->pGrBeta,pGrBeta,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
    
    //maturity parameters
    if (debug<0) cout<<tb<<"maturity=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrM2M->pvLgtPrM2M,pvLgtPrM2M,debug,cout); if (debug<0){cout<<endl;}
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
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pLnC, pLnC, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC1, pDC1, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC2, pDC2, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC3, pDC3, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDC4, pDC4, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrFsh->pDevsLnC,devsLnC,debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<cc<<endl;
   
    //survey catchability parameters
    if (debug<0) cout<<tb<<"surveys=list("<<endl;
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pQ,   pQ,   debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ1, pDQ1, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ2, pDQ2, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ3, pDQ3, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pDQ4, pDQ4, debug,cout); if (debug<0){cout<<cc<<endl;}
    if (debug<0) {cout<<tb;} tcsam::calcPriors(objFun,ptrMPI->ptrSrv->pA,   pA,   debug,cout); if (debug<0){cout<<endl;}
    if (debug<0) cout<<tb<<")"<<endl;
    
    if (debug<0) cout<<")"<<endl;
    
    if (!debug) testNaNs(value(objFun),"in calcAllPriors()");
    if (debug>=dbgPriors) cout<<"Finished calcAllPriors()"<<endl;

//*****************************************
FUNCTION void setInitVals(int debug, ostream& os)
    //recruitment parameters
    setInitVals(ptrMPI->ptrRec->pLnR, pLnR, usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRCV, pRCV, usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRX,  pRX,  usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRa,  pRa,  usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pRb,  pRb,  usePin, debug, os);
    setInitVals(ptrMPI->ptrRec->pDevsLnR,pDevsLnR,usePin, debug, os);

    //natural mortality parameters
    setInitVals(ptrMPI->ptrNM->pM,   pM,   usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM1, pDM1, usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM2, pDM2, usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM3, pDM3, usePin, debug, os);
    setInitVals(ptrMPI->ptrNM->pDM4, pDM4, usePin, debug, os);

    //growth parameters
    setInitVals(ptrMPI->ptrGrw->pGrA,   pGrA,   usePin, debug, os);
    setInitVals(ptrMPI->ptrGrw->pGrB,   pGrB,   usePin, debug, os);
    setInitVals(ptrMPI->ptrGrw->pGrBeta,pGrBeta,usePin, debug, os);

    //maturity parameters
    setInitVals(ptrMPI->ptrM2M->pvLgtPrM2M,pvLgtPrM2M,usePin, debug, os);

    //selectivity parameters
    setInitVals(ptrMPI->ptrSel->pS1, pS1,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS2, pS2,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS3, pS3,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS4, pS4,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS5, pS5,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pS6, pS6,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS1, pDevsS1,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS2, pDevsS2,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS3, pDevsS3,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS4, pDevsS4,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS5, pDevsS5,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pDevsS6, pDevsS6,usePin, debug, os);
    setInitVals(ptrMPI->ptrSel->pvNPSel, pvNPSel,usePin, debug, os);

    //fully-selected fishing capture rate parameters
    setInitVals(ptrMPI->ptrFsh->pHM,  pHM,  usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pLnC, pLnC, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC1, pDC1, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC2, pDC2, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC3, pDC3, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDC4, pDC4, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pDevsLnC,pDevsLnC,usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pLnEffX, pLnEffX, usePin, debug, os);
    setInitVals(ptrMPI->ptrFsh->pLgtRet, pLgtRet, usePin, debug, os);

    //survey catchability parameters
    setInitVals(ptrMPI->ptrSrv->pQ,   pQ,   usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ1, pDQ1, usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ2, pDQ2, usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ3, pDQ3, usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pDQ4, pDQ4, usePin, debug, os);
    setInitVals(ptrMPI->ptrSrv->pA,   pA,   usePin, debug, os);

    //MSE parameters    
    setInitVals(ptrMPI->ptrMSE->pMSE_LnC, pMSE_LnC, usePin, debug, os);

//*****************************************
FUNCTION int checkParams(int debug, ostream& os)
    int res = 0;
    //recruitment parameters
    res += checkParams(ptrMPI->ptrRec->pLnR, pLnR, debug,os);
    res += checkParams(ptrMPI->ptrRec->pRCV, pRCV, debug,os);
    res += checkParams(ptrMPI->ptrRec->pRX,  pRX,  debug,os);
    res += checkParams(ptrMPI->ptrRec->pRa,  pRa,  debug,os);
    res += checkParams(ptrMPI->ptrRec->pRb,  pRb,  debug,os);
    res += checkParams(ptrMPI->ptrRec->pDevsLnR,pDevsLnR,debug,os);

    //natural mortality parameters
    res += checkParams(ptrMPI->ptrNM->pM, pM, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM1, pDM1, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM2, pDM2, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM3, pDM3, debug,os);
    res += checkParams(ptrMPI->ptrNM->pDM4, pDM4, debug,os);

    //growth parameters
    res += checkParams(ptrMPI->ptrGrw->pGrA,   pGrA,   debug,os);
    res += checkParams(ptrMPI->ptrGrw->pGrB,   pGrB,   debug,os);
    res += checkParams(ptrMPI->ptrGrw->pGrBeta,pGrBeta,debug,os);

    //maturity parameters
    res += checkParams(ptrMPI->ptrM2M->pvLgtPrM2M,pvLgtPrM2M,debug,os);

    //selectivity parameters
    res += checkParams(ptrMPI->ptrSel->pS1, pS1,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS2, pS2,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS3, pS3,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS4, pS4,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS5, pS5,debug,os);
    res += checkParams(ptrMPI->ptrSel->pS6, pS6,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS1, pDevsS1,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS2, pDevsS2,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS3, pDevsS3,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS4, pDevsS4,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS5, pDevsS5,debug,os);
    res += checkParams(ptrMPI->ptrSel->pDevsS6, pDevsS6,debug,os);
    res += checkParams(ptrMPI->ptrSel->pvNPSel, pvNPSel,debug,os);

    //fully-selected fishing capture rate parameters
    res += checkParams(ptrMPI->ptrFsh->pHM,     pHM,     debug,os);
    res += checkParams(ptrMPI->ptrFsh->pLnC,    pLnC,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC1,    pDC1,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC2,    pDC2,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC3,    pDC3,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDC4,    pDC4,    debug,os);
    res += checkParams(ptrMPI->ptrFsh->pDevsLnC,pDevsLnC,debug,os);
    res += checkParams(ptrMPI->ptrFsh->pLnEffX, pLnEffX, debug,os);
    res += checkParams(ptrMPI->ptrFsh->pLgtRet, pLgtRet, debug,os);

    //survey catchability parameters
    res += checkParams(ptrMPI->ptrSrv->pQ,   pQ,   debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ1, pDQ1, debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ2, pDQ2, debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ3, pDQ3, debug,os);
    res += checkParams(ptrMPI->ptrSrv->pDQ4, pDQ4, debug,os);
    res += checkParams(ptrMPI->ptrSrv->pA,   pA, debug,os);

    //MSE parameters
    res += checkParams(ptrMPI->ptrMSE->pMSE_LnC,pMSE_LnC,debug,os);
    
    return res;

//******************************************************************************
//* Function: void checkParams(NumberVectorInfo* pI, param_init_number_vector& p, int debug, ostream& cout)
//* 
//* Description: Check values for a parameter vector.
//*
//* Inputs:
//*  pI (NumberVectorInfo*) 
//*     pointer to NumberVectorInfo object
//*  p (param_init_number_vector&)
//*     parameter vector
//* Returns:
//*  int - 0 = all good, 1 = nan detected
//* Alters:
//*  none
//******************************************************************************
FUNCTION int checkParams(NumberVectorInfo* pI, param_init_number_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting checkParams(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            if (isnan(value(p(i)))){
                r++;
                if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<"]"<<endl;
            }
        }
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;
     
//******************************************************************************
//* Function: void checkParams(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
//* 
//* Description: Check values for a parameter vector.
//*
//* Inputs:
//*  pI (BoundedNumberVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_number_vector&)
//*     parameter vector
//* Returns:
//*  int - 0 = all good, 1 = nan detected
//* Alters:
//*  none
//******************************************************************************
FUNCTION int checkParams(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting checkParams(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            if (isnan(value(p(i)))){
                r++;
                if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<"]"<<endl;
            }
        }
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;

//******************************************************************************
//* Function: void checkParams(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//* 
//* Description: Checks values for a vector of parameter vectors.
//*
//* Inputs:
//*  pI (BoundedVectorVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_vector_vector&)
//*     parameter vector
//* Returns:
//*  int - 0 = all good, 1 = nan detected
//* Alters:
//*  none
//******************************************************************************
FUNCTION int checkParams(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting checkParams(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) {
                if (isnan(value(p(i,j)))){
                    r++;
                    if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<cc<<j<<"]"<<endl;
                }
            }
        }
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;

//******************************************************************************
//* Function: void checkParams(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//* 
//* Description: Checks values for a vector of parameter vectors.
//*
//* Inputs:
//*  pI (DevsVectorVectorInfo*) 
//*     pointer to BoundedNumberVectorInfo object
//*  p (param_init_bounded_vector_vector&)
//*     parameter vector
//* Returns:
//*  int - 0 = all good, 1 = nan detected
//* Alters:
//*  none
//******************************************************************************
FUNCTION int checkParams(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, ostream& cout)
//    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting checkParams(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    int r = 0;
    if (np){
        for (int i=1;i<=np;i++) {
            for (int j=p(i).indexmin();j<=p(i).indexmax();j++) {
                if (isnan(value(p(i,j)))){
                    r++;
                    if (debug) cout<<"NaN detected in checkParams() for "<<p(1).label()<<"["<<i<<cc<<j<<"]"<<endl;
                }
            }
        }
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Finished checkParams(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }
    return r;

//----------------------------------------------------------------------------------
//write header to MCMC eval file
FUNCTION writeMCMCHeader
    mcmc.open((char*)(fnMCMC),ofstream::out|ofstream::trunc);
    ctrMCMC = 0;
    mcmc<<"mcmc<-list();"<<endl;
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
    ctrMCMC+=1;
    std::cout<<"writing mcmc iteration "<<ctrMCMC<<endl;
    mcmc<<"mcmc[["<<ctrMCMC<<"]]<-list(objFun="<<objFun<<cc<<endl;
    //write parameter values
        //recruitment values
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pLnR);     mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRCV);     mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRX);      mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRa);      mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pRb);      mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrRec->pDevsLnR); mcmc<<cc<<endl;

        //natural mortality parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pM);   mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrNM->pDM4); mcmc<<cc<<endl;

        //growth parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pGrA);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pGrB);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrGrw->pGrBeta); mcmc<<cc<<endl;

        //maturity parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrM2M->pvLgtPrM2M); mcmc<<cc<<endl;

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
        writeMCMCtoR(mcmc,ptrMPI->ptrSel->pvNPSel); mcmc<<cc<<endl;

        //fully-selected fishing capture rate parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pHM);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnC); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC1); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC2); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC3); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDC4); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pDevsLnC); mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLnEffX);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrFsh->pLgtRet);  mcmc<<cc<<endl;

        //survey catchability parameters
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pQ);    mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ1);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ2);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ3);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pDQ4);  mcmc<<cc<<endl;
        writeMCMCtoR(mcmc,ptrMPI->ptrSrv->pA);    mcmc<<cc<<endl;
    
        //write other quantities
        mcmc<<"R_y="; wts::writeToR(mcmc,value(R_y)); mcmc<<cc<<endl;
        ivector bnds = wts::getBounds(spB_yx);
        mcmc<<"MB_xy="; wts::writeToR(mcmc,trans(value(spB_yx)),xDms,yDms); 
        if (doOFL){
            mcmc<<cc<<endl;
            calcOFL(mxYr+1,0,cout);//updates oflresults
            ptrOFLResults->writeToR(mcmc,ptrMC,"ptrOFLResults",0);//mcm<<cc<<endl;
        }
        
    mcmc<<");"<<endl;
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
//* Function: void setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, int usePin, int debug, ostream& os)
//* 
//* Description: Sets initial values for a vector of parameters from the associated NumberInfo object.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//*  @param pI - pointer to a NumberVectorInfo object
//*  @param p - reference to a param_init_number_vector
//*  @param usePin - flag to use values 
//*  @param debug - flag to print debugging info
//*  @param os - stream for debugging output
//*
//*  @alters pI - if usePin=1 or jittering or resampling occurs, then the initial values will be updated 
//*  @alters p - if usePin=0, the initial values will be updated 
//******************************************************************************
FUNCTION void setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, int usePin, int debug, ostream& os)
//    debug=dbgAll;
    if (debug>=dbgAll) os<<"Starting setInitVals(NumberVectorInfo* pI, param_init_number_vector& p, usePin, debug, os) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        //parameters have been defined
        dvector aovls = 1.0*pI->getInitVals();             //original initial values on arithmetic scales from parameter info
        dvector povls = 1.0*pI->getInitValsOnParamScales();//original initial values on parameter scales from parameter info
        dvector afvls = 1.0*pI->getInitVals();             //final initial values on arithmetic scales from parameter info
        if (usePin) {
            //use values from pinfile assigned to p as initial values
            os<<"Using pin file to set initial values for "<<p(1).label()<<endl;
            pI->setInitValsFromParamVals(p);//update initial values in pI to those from pinfile
            afvls = 1.0*pI->getInitVals();
            os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
            for (int i=1;i<=np;i++) os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
        } else {
            //use values from pI as initial values
            for (int i=1;i<=np;i++) {
                p(i) = povls(i);  //assign original initial value from parameter info
                NumberInfo* ptrI = (*pI)[i];
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    //assign final initial value based on resampling prior pdf
                    os<<"Using resampling to set initial values for "<<p(i).label()<<endl;
                    afvls(i) = ptrI->drawInitVal(rng,ptrMC->vif); //get resampled initial value on arithmetic scale
                    p(i) = ptrI->calcParamScaleVal(afvls(i));     //calc initial parameter value
                    ptrI->setInitVal(afvls(i));                  //update ptrI
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                }
            }
        }
    } else {
        //parameters have not been defined
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        os<<"Finished setInitVals(NumberVectorInfo* pI, param_init_number_vector& p) for "<<p(1).label()<<endl; 
    }
     
//******************************************************************************
//* Function: void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int usePin, int debug, ostream& os)
//* 
//* Sets initial values for a vecotr of bounded parameters fro the associated BoundedNumberVectorInfo object.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//*  @param pI - pointer to a BoundedNumberVectorInfo object
//*  @param p - reference to a param_init_number_vector
//*  @param usePin - flag to use values 
//*  @param debug - flag to print debugging info
//*  @param os - stream for debugging output
//*
//*  @alters pI - if usePin=1 or jittering or resampling occurs, then the initial values will be updated 
//*  @alters p - if usePin=0, the initial values will be updated 
//******************************************************************************
FUNCTION void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int usePin, int debug, ostream& os)
//    debug=dbgAll;
    if (debug>=dbgAll) os<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        //parameters have been defined
        dvector aovls = 1.0*pI->getInitVals();             //original initial values on arithmetic scales from parameter info
        dvector povls = 1.0*pI->getInitValsOnParamScales();//original initial values on parameter scales from parameter info
        dvector afvls = 1.0*pI->getInitVals();             //original initial values on arithmetic scales from parameter info
        if (usePin) {
            //use values from pinfile assigned to p as initial values
            os<<"Using pin file to set initial values for "<<p(1).label()<<endl;
            pI->setInitValsFromParamVals(p);//update initial values in pI to those from pinfile
            afvls = 1.0*pI->getInitVals();
            os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
            for (int i=1;i<=np;i++) os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
        } else {
            //use values from pI as initial values
            for (int i=1;i<=np;i++) {
                p(i) = povls(i);  //assign original initial value from parameter info
                BoundedNumberInfo* ptrI = (*pI)[i];
                if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                    //assign final initial value based on jittering
                    os<<"Using jittering to set initial values for "<<p(i).label()<<endl;
                    afvls(i) = ptrI->jitterInitVal(rng,ptrMC->jitFrac);//get jittered initial value on arithmetic scale
                    p(i) = ptrI->calcParamScaleVal(afvls(i));          //calc initial parameter value
                    ptrI->setInitVal(afvls(i));                       //update ptrI
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                } else 
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    //assign final initial value based on resampling prior pdf
                    os<<"Using resampling to set initial values for "<<p(i).label()<<endl;
                    afvls(i) = ptrI->drawInitVal(rng,ptrMC->vif); //get resampled initial value on arithmetic scale
                    p(i) = ptrI->calcParamScaleVal(afvls(i));     //calc initial parameter value
                    ptrI->setInitVal(afvls(i));                  //update ptrI
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                    os<<"param  aovalue   povalue    afvalue   pfvalue"<<endl;
                    os<<tb<<p(i).label()<<": "<<aovls(i)<<tb<<povls(i)<<tb<<afvls(i)<<tb<<value(p(i))<<endl;
                }
            }
        }
    } else {
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).label()<<endl; 
    }

//******************************************************************************
//* Function: void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int usePin, int debug, ostream& os)
//* 
//* Sets initial values for a vector of parameter vectors.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  @param pI - pointer to BoundedNumberVectorInfo object
//*  @param p - reference to a param_init_bounded_vector_vector
//*  @param usePin - flag to use values 
//*  @param debug - flag to print debugging info
//*  @param os - stream for debugging output
//* 
//*  @alters pI - if usePin=1 or jittering or resampling occurs, then the initial values will be updated 
//*  @alters p - if usePin=0, the initial values will be updated 
//******************************************************************************
FUNCTION void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int usePin, int debug, ostream& os)
    //debug=dbgAll;
    if (debug>=dbgAll) os<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        if (usePin){
            //use values from pinfile assigned to p as initial values
            for (int i=1;i<=np;i++) {
                os<<"Using pin file to set initial values for "<<p(i).label()<<endl;
                BoundedVectorInfo* ptrI = (*pI)[i];
                dvector pnvls = value(p(i));                       //original initial values on parameter scale from pin file
                dvector aovls = 1.0*ptrI->getInitVals();            //original initial values on arithmetic scale from parameter info
                dvector povls = 1.0*ptrI->getInitValsOnParamScale();//original initial values on parameter scale from parameter info
                ptrI->setInitValsFromParamVals(p(i));               //set final initial values on arithmetic scale for parameter info
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                os<<tb<<"orig param inits : "<<povls<<endl;
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        } else {
            //use values based on pI as initial values for p
            for (int i=1;i<=np;i++) {
                os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                dvector pnvls = value(p(i));                            //original initial values on parameter scale from pin file
                if (debug>=dbgAll) os<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<pnvls.indexmin()<<tb<<pnvls.indexmax()<<endl;
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                BoundedVectorInfo* ptrI = static_cast<DevsVectorInfo*>((*pI)[i]);
                dvector aovls = 1.0*ptrI->getInitVals();            //initial values on arithmetic scale from parameter info
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                dvector povls = 1.0*ptrI->getInitValsOnParamScale();//initial values on parameter scales from parameter info
                os<<tb<<"orig param inits : "<<povls<<endl;
                p(i)=povls;//set initial values on parameter scales from parameter info
                if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                    os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->jitterInitVals(rng,ptrMC->jitFrac);//get jittered values on arithmetic scale
                    ptrI->setInitVals(afvls);                               //set jittered values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();                   //get jittered values on param scale as initial parameter values
                } else
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values on arithmetic scale
                    ptrI->setInitVals(afvls);                          //set resampled values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();              //get resampled values on param scale as initial parameter values
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                }
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        }
    } else {
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        os<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }

//******************************************************************************
//* Function: void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int usePin, int debug, ostream& os)
//* 
//* Sets initial values for a vector of devs parameter vectors based on the associated DevsVectorVectorInfo object.
//*
//* Note: this function MUST be declared/defined as a FUNCTION in the tpl code
//*     because the parameter assignment is a private method but the model_parameters 
//*     class has friend access.
//* 
//* Inputs:
//*  @param pI - pointer to BoundedNumberVectorInfo object
//*  @param p - reference to parameter_init_bounded_vector_vector acting as devs_vector_vector
//*  @param usePin - flag to use values 
//*  @param debug - flag to print debugging info
//*  @param os - stream for debugging output
//* 
//*  @alters pI - if usePin=1 or jittering or resampling occurs, then the initial values will be updated 
//*  @alters p - if usePin=0, the initial values will be updated 
//******************************************************************************
FUNCTION void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int usePin, int debug, ostream& os)
    //debug=dbgAll;
    if (debug>=dbgAll) os<<"Starting setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    int np = pI->getSize();
    if (np){
        if (usePin){
            //use values from pinfile assigned to p as initial values
            for (int i=1;i<=np;i++) {
                os<<"Using pin file to set initial values for "<<p(i).label()<<endl;
                DevsVectorInfo* ptrI = static_cast<DevsVectorInfo*>((*pI)[i]);
                dvector pnvls = value(p(i));                       //original initial values on parameter scale from pin file
                dvector aovls = 1.0*ptrI->getInitVals();            //original initial values on arithmetic scale from parameter info
                dvector povls = 1.0*ptrI->getInitValsOnParamScale();//original initial values on parameter scale from parameter info
                ptrI->setInitValsFromParamVals(p(i));               //set final initial values on arithmetic scale for parameter info
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                os<<tb<<"orig param inits : "<<povls<<endl;
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        } else {
            //use values based on pI as initial values for p
            for (int i=1;i<=np;i++) {
                dvector pnvls = value(p(i));                            //original initial values on parameter scale from pin file
                dvector aovls = 1.0*(*pI)[i]->getInitVals();            //initial values on arithmetic scale from parameter info
                dvector povls = 1.0*(*pI)[i]->getInitValsOnParamScale();//initial values on parameter scales from parameter info
                if (debug>=dbgAll) os<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<pnvls.indexmin()<<tb<<pnvls.indexmax()<<endl;
                p(i)=povls;//set initial values on parameter scales from parameter info
                DevsVectorInfo* ptrI = (*pI)[i];
                if ((p(i).get_phase_start()>0)&&(ptrMC->jitter)&&(ptrI->jitter)){
                    os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using jittering to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->jitterInitVals(rng,ptrMC->jitFrac);//get jittered values on arithmetic scale
                    ptrI->setInitVals(afvls);                               //set jittered values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();                   //get jittered values on param scale as initial parameter values
                } else
                if ((p(i).get_phase_start()>0)&&(ptrMC->resample)&&(ptrI->resample)){
                    os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    if (debug>=dbgAll) os<<"Using resampling to set initial values for "<<p(i).label()<<" : "<<endl;
                    dvector afvls = ptrI->drawInitVals(rng,ptrMC->vif);//get resampled values on arithmetic scale
                    ptrI->setInitVals(afvls);                          //set resampled values as initial values on arithmetic scale
                    p(i)=ptrI->getInitValsOnParamScale();              //get resampled values on param scale as initial parameter values
                } else {
                    os<<"Using MPI to set initial values for "<<p(i).label()<<endl;
                }
                os<<tb<<"pinfile    inits : "<<pnvls<<endl;
                os<<tb<<"orig arith inits : "<<aovls<<endl;
                os<<tb<<"orig param inits : "<<povls<<endl;
                os<<tb<<"final arith inits: "<<ptrI->getInitVals()<<endl;//final initial values on arithmetic scale from parameter info
                os<<tb<<"final param inits: "<<value(p(i))<<endl;        //final initial values on parameter scale from parameter info
            }
        }
    } else {
        os<<"InitVals for "<<p(1).label()<<" not defined because np = "<<np<<endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        os<<"Finished setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).label()<<endl; 
    }

//-------------------------------------------------------------------------------------
FUNCTION void setAllDevs(int debug, ostream& cout)
    if (debug>=dbgAll) cout<<"starting setAllDevs()"<<endl;

    if (debug>=dbgAll) cout<<"setDevs() for pLnR"<<endl;
    tcsam::setDevs(devsLnR, pDevsLnR, ptrMPI->ptrRec->pDevsLnR,debug,cout);

    if (debug>=dbgAll) cout<<"setDevs() for pDevsS1"<<endl;
    tcsam::setDevs(devsS1, pDevsS1, ptrMPI->ptrSel->pDevsS1, debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS2"<<endl;
    tcsam::setDevs(devsS2, pDevsS2, ptrMPI->ptrSel->pDevsS2,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS3"<<endl;
    tcsam::setDevs(devsS3, pDevsS3, ptrMPI->ptrSel->pDevsS3,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS4"<<endl;
    tcsam::setDevs(devsS4, pDevsS4, ptrMPI->ptrSel->pDevsS4,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS5"<<endl;
    tcsam::setDevs(devsS5, pDevsS5, ptrMPI->ptrSel->pDevsS5,debug,cout);
    if (debug>=dbgAll) cout<<"setDevs() for pDevsS6"<<endl;
    tcsam::setDevs(devsS6, pDevsS6, ptrMPI->ptrSel->pDevsS6,debug,cout);
    
    if (debug>=dbgAll) cout<<"setDevs() for pDevsLnC"<<endl;
    tcsam::setDevs(devsLnC, pDevsLnC, ptrMPI->ptrFsh->pDevsLnC,debug,cout);
    
    if (debug>=dbgAll) cout<<"finished setAllDevs()"<<endl;

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
        os<<"M_c      ="; wts::writeToR(os,value(M_c),     adstring("pc=1:"+str(npcNM )));      os<<cc<<endl;
        os<<"prM2M_cz ="; wts::writeToR(os,value(prM2M_cz),adstring("pc=1:"+str(npcM2M)),zbDms);os<<cc<<endl;
        os<<"sel_cz   ="; wts::writeToR(os,value(sel_cz),  adstring("pc=1:"+str(npcSel)),zbDms);os<<cc<<endl;
        os<<"M_yxmsz  ="; wts::writeToR(os,wts::value(M_yxmsz),   yDms,xDms,mDms,sDms,zbDms);   os<<cc<<endl;
        os<<"prM2M_yxz =";wts::writeToR(os,     value(prM2M_yxz), yDms,xDms,zbDms);             os<<cc<<endl;
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
        
        if (nEASs){
        os<<"EffX_list=list("<<endl;
            os<<"effAvgScenarios="; ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios->writeToR(os); os<<cc<<endl;
            os<<"capRateScenarios="; ptrMOs->ptrEffXtrapScenarios->ptrCapRateAvgScenarios->writeToR(os); os<<cc<<endl;
            os<<"avgEff=list("<<endl;
            for (int n=eff_ny.indexmin();n<=eff_ny.indexmax();n++){
                os<<"`"<<n<<"`=list(avgEff="<<avgEff_n(n)<<cc;
                os<<"eff="; wts::writeToR(os,eff_ny(n),ptrMOs->ptrEffXtrapScenarios->ptrEffAvgScenarios->ppEASs[n-1]->getYDimsForR()); os<<"), "<<endl;
            }
            os<<"NULL)"<<endl;
        os<<")"<<cc<<endl;
        } else {
            os<<"EffX_list=NULL"<<cc<<endl;
        }
        
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
            os<<"avl_vyxmsz="; wts::writeToR(os,wts::value(a_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"sel_vyxmsz="; wts::writeToR(os,wts::value(s_vyxmsz),vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"A_vyxms   ="; wts::writeToR(os,wts::value(a_vyxms), vDms,ypDms,xDms,mDms,sDms);       os<<cc<<endl;
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
    
    //population male maturity ogives
    dmatrix vPM_yz(mnYr,mxYrp1,1,nZBs);
    for (int y=mnYr;y<=mxYrp1;y++){
        vPM_yz(y) = elem_div(vn_yxmsz(y,MALE,MATURE,NEW_SHELL),
                            (1.0e-5)+vn_yxmsz(y,MALE,IMMATURE,NEW_SHELL)+vn_yxmsz(y,MALE,MATURE,NEW_SHELL));
    }
        
    //numbers, biomass captured (NOT mortality)
    d6_array vcpN_fyxmsz = wts::value(cpN_fyxmsz);
    d6_array vcpB_fyxmsz = tcsam::calcBiomass(vcpN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass discards (NOT mortality)
    d6_array vdsN_fyxmsz = wts::value(dsN_fyxmsz);
    d6_array vdsB_fyxmsz = tcsam::calcBiomass(vdsN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass retained (mortality)
    d6_array vrmN_fyxmsz = wts::value(rmN_fyxmsz);
    d6_array vrmB_fyxmsz = tcsam::calcBiomass(vrmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //numbers, biomass discard mortality
    d6_array vdmN_fyxmsz = wts::value(dmN_fyxmsz);
    d6_array vdmB_fyxmsz = tcsam::calcBiomass(vdmN_fyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    
    //survey numbers, biomass
    d6_array vn_vyxmsz = wts::value(n_vyxmsz);
    d6_array vb_vyxmsz = tcsam::calcBiomass(vn_vyxmsz,ptrMDS->ptrBio->wAtZ_xmz);
    d3_array vPM_vyz(1,nSrv,mnYr,mxYrp1,1,nZBs);
    for (int v=1;v<=nSrv;v++){
        for (int y=mnYr;y<=mxYrp1;y++){
            vPM_vyz(v,y) = elem_div(vn_vyxmsz(v,y,MALE,MATURE,NEW_SHELL),
                                    (1.0e-5)+vn_vyxmsz(v,y,MALE,IMMATURE,NEW_SHELL)+vn_vyxmsz(v,y,MALE,MATURE,NEW_SHELL));
        }
    }
    
    os<<"mr=list("<<endl;
        os<<"iN_xmsz ="; wts::writeToR(os,vn_yxmsz(mnYr),xDms,mDms,sDms,zbDms); os<<cc<<endl;
        os<<"P_list=list("<<endl;
            os<<"MB_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms);                        os<<cc<<endl;
            os<<"N_yxmsz  ="; wts::writeToR(os,vn_yxmsz,ypDms,xDms,mDms,sDms,zbDms);             os<<cc<<endl;
            os<<"B_yxmsz  ="; wts::writeToR(os,vb_yxmsz,ypDms,xDms,mDms,sDms,zbDms);             os<<cc<<endl;
            os<<"prM_yz   ="; wts::writeToR(os,vPM_yz,ypDms,zbDms);                              os<<cc<<endl;
            os<<"nmN_yxmsz="; wts::writeToR(os,wts::value(nmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"tmN_yxmsz="; wts::writeToR(os,wts::value(tmN_yxmsz),yDms,xDms,mDms,sDms,zbDms); os<<endl;
        os<<")"<<cc<<endl;    
        os<<"F_list=list("<<endl;
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
           os<<"MB_vyx  ="; wts::writeToR(os,value(mb_vyx),vDms,ypDms,xDms);                 os<<cc<<endl;
           os<<"N_vyxmsz="; wts::writeToR(os,    vn_vyxmsz,vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
           os<<"B_vyxmsz="; wts::writeToR(os,    vb_vyxmsz,vDms,ypDms,xDms,mDms,sDms,zbDms); os<<cc<<endl;
           os<<"prM_vyz ="; wts::writeToR(os,    vPM_vyz,vDms,ypDms,zbDms);                  os<<endl;
       os<<")";
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelResults(...)"<<endl;
    
//-------------------------------------------------------------------------------------
//Write cohort progression quantities to file as R list
///**
// * Write cohort progression quantities to file as R list.
// * 
// * @param os  - output stream to  write to
// * @param nzp - number of size bins to include for initial recruitment
// * @param includeM - flag (0,1) to include natural mortality
// * @param includeF - flag (0,1) to include fishing mortality
// * @param xtra - flag (0,1) to include "extra" info (prM2M and growth arrays)
// * @param debug - flag to print debugging info
// * @param cout - output stream to write to
// * 
// * @return 5d array of cohort abundance by yxmsz (n_yxmsz)
// */
FUNCTION void ReportToR_CohortProgression(ostream& os, int nzp, int includeM, int includeF, int xtra, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_CohortProgression(...)"<<endl;
    d5_array n_yxmsz = calcCohortProgression(mxYr,nzp,includeM,includeF,debug,cout);
    int ny = n_yxmsz.indexmax();
    if (xtra){
        dmatrix  prM2M_xz  = value(prM2M_yxz(mxYr));
        d3_array mnZAM_xsz = value(mnGrZ_yxsz(mxYr));
        ivector pd3(1,3); pd3(1)=2; pd3(2)=1; pd3(3)=3;
        dmatrix mnZAM_xz   = wts::permuteDims(pd3,mnZAM_xsz)(NEW_SHELL);//new shelll same as old shell
        d4_array  T_xszz   = wts::value(prGr_yxszz(mxYr));
        ivector pd4(1,4); pd4(1)=2; pd4(2)=1; pd4(3)=3; pd4(4)=4;
        d3_array T_xzz   = wts::permuteDims(pd4,T_xszz)(NEW_SHELL);//new shell same as old shell
        os<<"cohortprogression=list("<<endl;
            os<<"n_yxmsz  ="; wts::writeToR(os,n_yxmsz,adstring("y=0:"+str(ny)),xDms,mDms,sDms,zbDms); os<<cc<<endl;
            os<<"prM2M_xz ="; wts::writeToR(os,prM2M_xz,xDms,zbDms);    os<<cc<<endl;
            os<<"mnZAM_xz ="; wts::writeToR(os,mnZAM_xz,xDms,zbDms);    os<<cc<<endl;
            os<<"T_xzz    ="; wts::writeToR(os,T_xzz,xDms,zbDms,zpDms); os<<endl;
        os<<")";
    } else {
        os<<"cohortprogression=list("<<endl;
            os<<"n_yxmsz="; wts::writeToR(os,n_yxmsz,adstring("y=0:"+str(ny)),xDms,mDms,sDms,zbDms); os<<endl;
        os<<")";
    }
    if (debug) cout<<"Finished ReportToR_CohortProgression(...)"<<endl;
    
//-------------------------------------------------------------------------------------
//Write quantities related to model fits to file as R list
FUNCTION void ReportToR_ModelFits(ostream& os, double maxGrad, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_ModelFits(...)"<<endl;
    os<<"model.fits=list("<<endl;
        os<<tb<<"objfun="<<value(objFun)<<cc<<"maxGrad="<<maxGrad<<cc<<endl;
        //recalc objective function components and and write results to os
        objFun.initialize();
        os<<tb<<"penalties="; calcPenalties(-1,os);      os<<cc<<endl;
        os<<tb<<"#end of penalties"<<endl;
        os<<tb<<"priors=";    calcAllPriors(-1,os);      os<<cc<<endl;
        os<<tb<<"#end of priors"<<endl;
        os<<tb<<"components=list("<<endl;
            os<<tb<<tb<<"recruitment="; calcNLLs_Recruitment(-1,os); os<<endl;
        os<<tb<<")"<<cc<<endl;
        os<<tb<<"#end of components"<<endl;
        os<<tb<<"fisheries=";  calcNLLs_Fisheries(-1,os);          os<<cc<<endl; 
        os<<tb<<"#end of fisheries"<<endl;
        os<<tb<<"surveys=";    calcNLLs_Surveys(-1,os);            os<<cc<<endl;  
        os<<tb<<"#end of surveys"<<endl;
        os<<tb<<"growthdata="; calcNLLs_GrowthData(-1,os);         os<<cc<<endl;
        os<<tb<<"#end of growthdata"<<endl;
        os<<tb<<"chelaheightdata="; calcNLLs_ChelaHeightData(-1,os);  os<<cc<<endl;
        os<<tb<<"#end of chelaheightdata"<<endl;
        os<<tb<<"maturityogivedata="; calcNLLs_MaturityOgiveData(-1,os);  os<<cc<<endl;
        os<<tb<<"#end of maturityogivedata"<<endl;
        os<<tb<<"effortdata="; calcNLLs_ExtrapolatedEffort(-1,os); os<<endl;
        os<<tb<<"#end of effortdata"<<endl;
    os<<")";
    if (debug) cout<<"Finished ReportToR_ModelFits(...)"<<endl;

//-------------------------------------------------------------------------------------
//Write quantities related to model fits to file as R list
FUNCTION void ReportToR_OFLResults(ostream& os, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_OFLResults(...)"<<endl;
    os<<tb<<"#--OFLResults:"<<endl;
    ptrOFLResults->writeToR(os,ptrMC,"oflResults",debug);
    os<<tb<<"#--end of OFLResults"<<endl;
    if (debug) cout<<"Finished ReportToR_OFLResults(...)"<<endl;

//-------------------------------------------------------------------------------------
//Update MPI for current parameter values (mainly for export))
FUNCTION void updateMPI(int debug, ostream& cout)
    if (debug) cout<<"Starting updateMPI(...)"<<endl;

    //recruitment parameters
    if (debug) cout<<"starting recruitment parameters"<<endl;
    ptrMPI->ptrRec->pLnR->setFinalValsFromParamVals(pLnR);
    ptrMPI->ptrRec->pRCV->setFinalValsFromParamVals(pRCV);
    ptrMPI->ptrRec->pRX->setFinalValsFromParamVals(pRX);
    ptrMPI->ptrRec->pRa->setFinalValsFromParamVals(pRa);
    ptrMPI->ptrRec->pRb->setFinalValsFromParamVals(pRb);
    if (debug) cout<<"setting final vals for pDevsLnR"<<endl;
    for (int p=1;p<=ptrMPI->ptrRec->pDevsLnR->getSize();p++) 
        (*ptrMPI->ptrRec->pDevsLnR)[p]->setFinalValsFromParamVals(pDevsLnR(p));
    if (debug) cout<<"finished recruitment parameters"<<endl;
     
    //natural mortality parameters
    if (debug) cout<<"starting natural mortality parameters"<<endl;
    ptrMPI->ptrNM->pM->setFinalValsFromParamVals(pM);
    ptrMPI->ptrNM->pDM1->setFinalValsFromParamVals(pDM1);
    ptrMPI->ptrNM->pDM2->setFinalValsFromParamVals(pDM2);
    ptrMPI->ptrNM->pDM3->setFinalValsFromParamVals(pDM3);
    ptrMPI->ptrNM->pDM4->setFinalValsFromParamVals(pDM4);
    if (debug) cout<<"finished natural mortality parameters"<<endl;
    
    //growth parameters
    if (debug) cout<<"starting growth parameters"<<endl;
    ptrMPI->ptrGrw->pGrA->setFinalValsFromParamVals(pGrA);
    ptrMPI->ptrGrw->pGrB->setFinalValsFromParamVals(pGrB);
    ptrMPI->ptrGrw->pGrBeta->setFinalValsFromParamVals(pGrBeta);
    if (debug) cout<<"finished growth parameters"<<endl;
    
    //maturity parameters
    //cout<<"setting final vals for pvLgtPrM2M"<<endl;
    if (debug) cout<<"starting prM2M parameters"<<endl;
    for (int p=1;p<=ptrMPI->ptrRec->pDevsLnR->getSize();p++) 
        (*ptrMPI->ptrM2M->pvLgtPrM2M)[p]->setFinalValsFromParamVals(pvLgtPrM2M(p));
    if (debug) cout<<"finished prM2M parameters"<<endl;
    
    //selectivity parameters
    if (debug) cout<<"starting selectivities"<<endl;
    ptrMPI->ptrSel->pS1->setFinalValsFromParamVals(pS1);
    ptrMPI->ptrSel->pS2->setFinalValsFromParamVals(pS2);
    ptrMPI->ptrSel->pS3->setFinalValsFromParamVals(pS3);
    ptrMPI->ptrSel->pS4->setFinalValsFromParamVals(pS4);
    ptrMPI->ptrSel->pS5->setFinalValsFromParamVals(pS5);
    ptrMPI->ptrSel->pS6->setFinalValsFromParamVals(pS6);
    if (debug) cout<<"setting final vals for pDevsS1"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS1->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS1)[p]->setFinalValsFromParamVals(pDevsS1(p));
    if (debug) cout<<"setting final vals for pDevsS2"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS2->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS2)[p]->setFinalValsFromParamVals(pDevsS2(p));
    if (debug) cout<<"setting final vals for pDevsS3"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS3->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS3)[p]->setFinalValsFromParamVals(pDevsS3(p));
    if (debug) cout<<"setting final vals for pDevsS4"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS4->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS4)[p]->setFinalValsFromParamVals(pDevsS4(p));
    if (debug) cout<<"setting final vals for pDevsS5"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS5->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS5)[p]->setFinalValsFromParamVals(pDevsS5(p));
    if (debug) cout<<"setting final vals for pDevsS6"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pDevsS6->getSize();p++) 
        (*ptrMPI->ptrSel->pDevsS6)[p]->setFinalValsFromParamVals(pDevsS6(p));
    if (debug) cout<<"setting final vals for pvNPSel"<<endl;
    for (int p=1;p<=ptrMPI->ptrSel->pvNPSel->getSize();p++) 
        (*ptrMPI->ptrSel->pvNPSel)[p]->setFinalValsFromParamVals(pvNPSel(p));
    if (debug) cout<<"finished selectivities"<<endl;
     
    //fully-selected fishing capture rate parameters
    if (debug) cout<<"starting fisheries parameters"<<endl;
    ptrMPI->ptrFsh->pHM->setFinalValsFromParamVals(pHM);
    ptrMPI->ptrFsh->pLnC->setFinalValsFromParamVals(pLnC);
    ptrMPI->ptrFsh->pDC1->setFinalValsFromParamVals(pDC1);
    ptrMPI->ptrFsh->pDC2->setFinalValsFromParamVals(pDC2);
    ptrMPI->ptrFsh->pDC3->setFinalValsFromParamVals(pDC3);
    ptrMPI->ptrFsh->pDC4->setFinalValsFromParamVals(pDC4);
    if (debug) cout<<"setting final vals for pDevsLnC"<<endl;
    for (int p=1;p<=ptrMPI->ptrFsh->pDevsLnC->getSize();p++) 
        (*ptrMPI->ptrFsh->pDevsLnC)[p]->setFinalValsFromParamVals(pDevsLnC(p));
    ptrMPI->ptrFsh->pLnEffX->setFinalValsFromParamVals(pLnEffX);
    ptrMPI->ptrFsh->pLgtRet->setFinalValsFromParamVals(pLgtRet);
    if (debug) cout<<"finished fisheries parameters"<<endl;
    
    //survey catchability parameters
    if (debug) cout<<"starting surveys parameters"<<endl;
    ptrMPI->ptrSrv->pQ->setFinalValsFromParamVals(pQ);
    ptrMPI->ptrSrv->pDQ1->setFinalValsFromParamVals(pDQ1);
    ptrMPI->ptrSrv->pDQ2->setFinalValsFromParamVals(pDQ2);
    ptrMPI->ptrSrv->pDQ3->setFinalValsFromParamVals(pDQ3);
    ptrMPI->ptrSrv->pDQ4->setFinalValsFromParamVals(pDQ4);
    ptrMPI->ptrSrv->pA->setFinalValsFromParamVals(pA);
    if (debug) cout<<"finished surveys parameters"<<endl;
    
    //MSE-related parameters
    ptrMPI->ptrMSE->pMSE_LnC->setFinalValsFromParamVals(pMSE_LnC);
    if (debug) cout<<"finished MSE parameters"<<endl;
    
    if (debug) cout<<"Finished updateMPI(...)"<<endl;

//-------------------------------------------------------------------------------------
//Write results to file as R list
FUNCTION void ReportToR(ostream& os, double maxGrad, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR(...)"<<endl;

    updateMPI(debug,cout);
        
    os<<"res=list("<<endl;
        os<<"mode='estModMode'"<<cc<<endl;
        if (jitter) os<<"jitter="<<iSeed<<cc<<endl;
        os<<"objFun="<<value(objFun)<<cc<<"maxGrad="<<maxGrad<<cc<<endl;
        //model configuration
        ptrMC->writeToR(os,"mc",0); os<<","<<endl;
        os<<tb<<"#end of mc"<<endl;
        
        //model data
        ptrMDS->writeToR(os,"data",0); os<<","<<endl;
        os<<tb<<"#end of data"<<endl;
        
        //parameter values
        ReportToR_Params(os,debug,cout); os<<","<<endl;
        os<<tb<<"#end of params"<<endl;
        
        //model processes
        ReportToR_ModelProcesses(os,debug,cout); os<<","<<endl;
        os<<tb<<"#end of modelprocesses"<<endl;
        
        //model results
        ReportToR_ModelResults(os,debug,cout); os<<","<<endl;
        os<<tb<<"#end of modelresults"<<endl;

        //model fit quantities
        ReportToR_ModelFits(os,maxGrad,debug,cout); os<<","<<endl;
        os<<tb<<"#end of modelfits"<<endl;
        
        //simulated model data
//        createSimData(debug, cout, 0, ptrSimMDS);//deterministic
//        ptrSimMDS->writeToR(os,"sim.data",0); os<<","<<endl;
//        os<<tb<<"#end of sim.data"<<endl;
        
        //cohort projections (recruits in 3 size bins, M+F included, no extra info)
        ReportToR_CohortProgression(os,3,1,1,0,debug,cout);
        os<<tb<<"#end of cohortprogression"<<endl;
        
        //do OFL calculations
        if (doOFL&&last_phase()){
            cout<<"ReportToR: starting OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
            echoOFL.precision(12);
            calcOFL(mxYr+1,1,echoOFL);//updates ptrOFLResults
            ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
            ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL.close();
            //pick up writing ReportTotR
            os<<","<<endl;
            ptrOFLResults->writeToR(os,ptrMC,"ptrOFLResults",0);
            os<<tb<<"#end of ptrOFLResults"<<endl;
            cout<<"ReportToR: finished OFL calculations"<<endl;
        }

        //do dynamic B0 calculations
        if (doDynB0&&last_phase()){
            cout<<"ReportToR: starting dynamic B0 calculations"<<endl;
            os<<","<<endl;
            calcDynB0(-1,os);
            os<<tb<<"#end of dynamic B0 results"<<endl;
            cout<<"ReportToR: finished dynamic B0 calculations"<<endl;
        }

    os<<")"<<endl;
    if (debug) cout<<"Finished ReportToR(...)"<<endl;

//-------------------------------------------------------------------------------------
//Write results to file as R list
FUNCTION void ReportToR_OpModMode(ostream& os, double maxGrad, int debug, ostream& cout)
    if (debug) cout<<"Starting ReportToR_OpModMode(...)"<<endl;

    updateMPI(debug,cout);
        
    os<<"res=list("<<endl;
        os<<"mode='opModMode'"<<cc<<"objFun="<<value(objFun)<<cc<<"maxGrad="<<maxGrad<<cc<<endl;
        //model configuration
        ptrMC->writeToR(os,"mc",0); os<<","<<endl;
        os<<tb<<"#--end of mc"<<endl;
        
        //parameter values
        ReportToR_Params(os,debug,cout); os<<cc<<endl;
        os<<tb<<"#--end of params"<<endl;
        
        //OpMod info
        os<<"info="; ptrOMI->writeToR(os); os<<cc<<endl;
        os<<tb<<"#--end of OpMod info"<<endl;
        
        //OpMod results
        os<<"results=list("<<endl;
            os<<"catchStats=list("<<endl;
                os<<"inpTAC="<<inpTAC<<cc<<"inpOFL="<<inpOFL<<cc<<endl;
        adstring tacDims = "y=";
        if (ptrOMI->nTACs==0){
            tacDims = tacDims+str(ptrOMI->mxYr)+":"+str(ptrOMI->mxYr);
            dvector TAC_y(ptrOMI->mxYr,ptrOMI->mxYr);
            dvector OFL_y(ptrOMI->mxYr,ptrOMI->mxYr);
            TAC_y(ptrOMI->mxYr) = inpTAC;
            OFL_y(ptrOMI->mxYr) = inpOFL;
        } else {
            tacDims = tacDims+str(ptrOMI->yrsTAC(1))+":"+str(ptrOMI->yrsTAC(ptrOMI->nTACs)+1);
            dvector TAC_y(ptrOMI->TAC_y.indexmin(),ptrOMI->TAC_y.indexmax()+1);
            dvector OFL_y(ptrOMI->OFL_y.indexmin(),ptrOMI->OFL_y.indexmax()+1);
            TAC_y(ptrOMI->TAC_y.indexmin(),ptrOMI->TAC_y.indexmax()) = ptrOMI->TAC_y;
            OFL_y(ptrOMI->OFL_y.indexmin(),ptrOMI->OFL_y.indexmax()) = ptrOMI->OFL_y;
            TAC_y(ptrOMI->TAC_y.indexmax()+1) = inpTAC;
            OFL_y(ptrOMI->OFL_y.indexmax()+1) = inpOFL;
                os<<"TAC_y=";     wts::writeToR(os,TAC_y,tacDims); os<<cc<<endl;
                os<<"OFL_y=";     wts::writeToR(os,OFL_y,tacDims); os<<cc<<endl;
        }
                os<<"capF="<<mfexp(pMSE_LnC[1])<<cc;
                os<<"retCatchMort="<<sum(prjRetCatchMortBio_fx)<<cc;
                os<<"totCatchMort="<<sum(prjTotCatchMortBio_fx)<<endl;
            os<<")"<<cc<<endl;
        adstring dimRecYrsToR = "y="+str(ptrOMI->R_y.indexmin())+":"+str(ptrOMI->R_y.indexmax()+1);
        dvector dR_y(ptrOMI->R_y.indexmin(),ptrOMI->R_y.indexmax()+1);
        dR_y(ptrOMI->R_y.indexmin(),ptrOMI->R_y.indexmax()) = ptrOMI->R_y;
        dR_y(ptrOMI->R_y.indexmax()+1)                      = prjR;
            os<<"recStats=list("<<endl;
                os<<"minYr="<<ptrMOs->opModRecStatsMinYr<<cc<<"maxYr="<<ptrMOs->opModRecStatsMaxYr<<cc<<endl;
                os<<"mnLnR="<<opModMnLnR<<cc<<"sdLnR="<<opModSdLnR<<cc<<"devLnR="<<opModDevLnR<<cc<<"prjR="<<prjR<<cc<<endl;
                os<<"R_y="; wts::writeToR(os,dR_y,dimRecYrsToR); 
            os<<")"<<cc<<endl;
            os<<"n_xmsz="; wts::writeToR(os,prj_n_xmsz,
                                         ptrMC->dimSXsToR,
                                         ptrMC->dimMSsToR,
                                         ptrMC->dimSCsToR,
                                         ptrMC->dimZBsToR); os<<endl;
            os<<")"<<endl<<tb<<"#--end of om results"<<endl;
        os<<")"<<endl<<"#end of OpMod report"<<endl;
        
    if (debug) cout<<"Finished ReportToR_OpModMode(...)"<<endl;

//----------------------------------------------------------------------
//Write parameter information to file
FUNCTION void writeParameters(ostream& os,int toR, int willBeActive)      
    adstring ctg1, ctg2;
    if (!toR) tcsam::writeCSVHeaderForParametersFile(os);
    //recruitment parameters
    ctg1="population processes";
    ctg2="recruitment";
    tcsam::writeParameters(os,pLnR,ctg1,ctg2,ptrMPI->ptrRec->pLnR,toR,willBeActive);      
    tcsam::writeParameters(os,pRCV,ctg1,ctg2,ptrMPI->ptrRec->pRCV,toR,willBeActive);      
    tcsam::writeParameters(os,pRX, ctg1,ctg2,ptrMPI->ptrRec->pRX, toR,willBeActive);      
    tcsam::writeParameters(os,pRa, ctg1,ctg2,ptrMPI->ptrRec->pRa, toR,willBeActive);      
    tcsam::writeParameters(os,pRb, ctg1,ctg2,ptrMPI->ptrRec->pRb, toR,willBeActive);      
    tcsam::writeParameters(os,pDevsLnR,ctg1,ctg2,ptrMPI->ptrRec->pDevsLnR,toR,willBeActive);      
    
    //natural mortality parameters
    ctg1="population processes";
    ctg2="natural mortality";
    tcsam::writeParameters(os,pM,ctg1,ctg2,ptrMPI->ptrNM->pM,toR,willBeActive);      
    tcsam::writeParameters(os,pDM1,ctg1,ctg2,ptrMPI->ptrNM->pDM1,toR,willBeActive);      
    tcsam::writeParameters(os,pDM2,ctg1,ctg2,ptrMPI->ptrNM->pDM2,toR,willBeActive);      
    tcsam::writeParameters(os,pDM3,ctg1,ctg2,ptrMPI->ptrNM->pDM3,toR,willBeActive);      
    tcsam::writeParameters(os,pDM4,ctg1,ctg2,ptrMPI->ptrNM->pDM4,toR,willBeActive);      
    
    //growth parameters
    ctg1="population processes";
    ctg2="growth";
    tcsam::writeParameters(os,pGrA,ctg1,ctg2,ptrMPI->ptrGrw->pGrA,toR,willBeActive);      
    tcsam::writeParameters(os,pGrB,ctg1,ctg2,ptrMPI->ptrGrw->pGrB,toR,willBeActive);      
    tcsam::writeParameters(os,pGrBeta,ctg1,ctg2,ptrMPI->ptrGrw->pGrBeta,toR,willBeActive);      
    
    //maturity parameters
    ctg1="population processes";
    ctg2="maturity";
    tcsam::writeParameters(os,pvLgtPrM2M,ctg1,ctg2,ptrMPI->ptrM2M->pvLgtPrM2M,toR,willBeActive);      
    
    //selectivity parameters
    ctg1="selectivity";
    ctg2="selectivity";
    tcsam::writeParameters(os,pS1,ctg1,ctg2,ptrMPI->ptrSel->pS1,toR,willBeActive);      
    tcsam::writeParameters(os,pS2,ctg1,ctg2,ptrMPI->ptrSel->pS2,toR,willBeActive);      
    tcsam::writeParameters(os,pS3,ctg1,ctg2,ptrMPI->ptrSel->pS3,toR,willBeActive);      
    tcsam::writeParameters(os,pS4,ctg1,ctg2,ptrMPI->ptrSel->pS4,toR,willBeActive);      
    tcsam::writeParameters(os,pS5,ctg1,ctg2,ptrMPI->ptrSel->pS5,toR,willBeActive);      
    tcsam::writeParameters(os,pS6,ctg1,ctg2,ptrMPI->ptrSel->pS6,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS1,ctg1,ctg2,ptrMPI->ptrSel->pDevsS1,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS2,ctg1,ctg2,ptrMPI->ptrSel->pDevsS2,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS3,ctg1,ctg2,ptrMPI->ptrSel->pDevsS3,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS4,ctg1,ctg2,ptrMPI->ptrSel->pDevsS4,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS5,ctg1,ctg2,ptrMPI->ptrSel->pDevsS5,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsS6,ctg1,ctg2,ptrMPI->ptrSel->pDevsS6,toR,willBeActive);      
    tcsam::writeParameters(os,pvNPSel,ctg1,ctg2,ptrMPI->ptrSel->pvNPSel,toR,willBeActive);      
    
    //fishery parameters
    ctg1="fisheries";
    ctg2="fisheries";
    tcsam::writeParameters(os,pHM, ctg1,ctg2,ptrMPI->ptrFsh->pHM, toR,willBeActive);      
    tcsam::writeParameters(os,pLnC,ctg1,ctg2,ptrMPI->ptrFsh->pLnC,toR,willBeActive);      
    tcsam::writeParameters(os,pDC1,ctg1,ctg2,ptrMPI->ptrFsh->pDC1,toR,willBeActive);      
    tcsam::writeParameters(os,pDC2,ctg1,ctg2,ptrMPI->ptrFsh->pDC2,toR,willBeActive);      
    tcsam::writeParameters(os,pDC3,ctg1,ctg2,ptrMPI->ptrFsh->pDC3,toR,willBeActive);      
    tcsam::writeParameters(os,pDC4,ctg1,ctg2,ptrMPI->ptrFsh->pDC4,toR,willBeActive);      
    tcsam::writeParameters(os,pDevsLnC,ctg1,ctg2,ptrMPI->ptrFsh->pDevsLnC,toR,willBeActive);      
    tcsam::writeParameters(os,pLnEffX, ctg1,ctg2,ptrMPI->ptrFsh->pLnEffX, toR,willBeActive);      
    tcsam::writeParameters(os,pLgtRet, ctg1,ctg2,ptrMPI->ptrFsh->pLgtRet, toR,willBeActive);      
    
    //survey parameters
    ctg1="surveys";
    ctg2="surveys";
    tcsam::writeParameters(os,pQ,  ctg1,ctg2,ptrMPI->ptrSrv->pQ,  toR,willBeActive);      
    tcsam::writeParameters(os,pDQ1,ctg1,ctg2,ptrMPI->ptrSrv->pDQ1,toR,willBeActive);      
    tcsam::writeParameters(os,pDQ2,ctg1,ctg2,ptrMPI->ptrSrv->pDQ2,toR,willBeActive);      
    tcsam::writeParameters(os,pDQ3,ctg1,ctg2,ptrMPI->ptrSrv->pDQ3,toR,willBeActive);      
    tcsam::writeParameters(os,pDQ4,ctg1,ctg2,ptrMPI->ptrSrv->pDQ4,toR,willBeActive);
    tcsam::writeParameters(os,pA,  ctg1,ctg2,ptrMPI->ptrSrv->pA,  toR,willBeActive);      
    
    //MSE directed fishery
    if (mseOpModMode){
        ctg1="MSE-related";
        ctg2="ln-scale directed fishery capture rate";
        tcsam::writeParameters(os,pMSE_LnC,ctg1,ctg2,ptrMPI->ptrMSE->pMSE_LnC,toR,willBeActive);
    }
    
// =============================================================================
// =============================================================================
REPORT_SECTION
    PRINT2B1(" ")
    PRINT2B1("--REPORT_SECTION---------")
        
    //max gradient
    double maxGrad = max(fabs(gradients));

    //write active parameters to rpt::echo
    PRINT2B2("Finished phase ",current_phase())
    {
        //write final gradients to file
        ofstream os0("tcsam02.Gradients."+itoa(current_phase(),10)+".csv", ios::trunc);
        for (int i=gradients.indexmin();i<=gradients.indexmax();i++){
            os0<<i<<cc<<gradients[i]<<endl;
        }
        os0.close();
    }
    if (!mseOpModMode) {
        //write objective function components only
        ofstream os0("tcsam02.ModelFits."+itoa(current_phase(),10)+".R", ios::trunc);
        os0.precision(12);
        ReportToR_ModelFits(os0,maxGrad,0,cout);
        os0.close();
    }
    if (last_phase()) {
        PRINT2B1("--ReportToR: last phase ")
        PRINT2B2("----last phase objFun =",value(objFun))
        if (!mseOpModMode){
            //write report as R file
            report.precision(12);
            ReportToR(report,maxGrad,1,rpt::echo);
            //write parameter values to csv
            ofstream os1("tcsam02.params.all.final.csv", ios::trunc);
            os1.precision(12);
            writeParameters(os1,0,0);
            os1.close();
            //write parameter values to csv
            ofstream os2("tcsam02.params.active.final.csv", ios::trunc);
            os2.precision(12);
            writeParameters(os2,0,1);
            os2.close();

            if (option_match(ad_comm::argc,ad_comm::argv,"-jitter")>-1) {
                ofstream fs("jitterInfo.csv");
                fs.precision(20);
                fs<<"seed"<<cc<<"objfun"<<cc<<"maxGrad"<<cc<<"MMB";
                if (doOFL) fs<<cc<<"B0"<<cc<<"Bmsy"<<cc<<"Fmsy"<<cc<<"OFL"<<cc<<"curB";
                fs<<endl;
                fs<<iSeed<<cc<<value(objFun)<<cc<<maxGrad<<cc<<spB_yx(mxYr,MALE);
                if (doOFL) fs<<cc<<ptrOFLResults->B0<<cc<<ptrOFLResults->Bmsy<<cc<<ptrOFLResults->Fmsy<<cc<<ptrOFLResults->OFL<<cc<<ptrOFLResults->curB;
                fs<<endl;
                fs.close();
            }
        } else if (mseOpModMode){
            PRINT2B1("#--Write REPORT FILE for OpModMode--")
            //write report as R file
            report.precision(12);
            ReportToR_OpModMode(report,maxGrad,1,rpt::echo);
            //write parameter values to csv
            ofstream os1("tcsam02.params.all.final.csv", ios::trunc);
            os1.precision(12);
            writeParameters(os1,0,0);
            os1.close();
            //write parameter values to csv
            ofstream os2("tcsam02.params.active.final.csv", ios::trunc);
            os2.precision(12);
            writeParameters(os2,0,1);
            os2.close();
        } 
    }
    PRINT2B1("--FINISHED REPORT_SECTION---------")

// =============================================================================
// =============================================================================
BETWEEN_PHASES_SECTION
    PRINT2B1(" ")
    PRINT2B1("#--BETWEEN_PHASES_SECTION---------------------")
    adstring msg = "#----Starting phase "+str(current_phase())+" of "+str(initial_params::max_number_phases);
    PRINT2B1(msg)
    if (mseOpModMode){
        dvariable mseCapF = mfexp(pMSE_LnC[1]);
        projectPopForTAC(mseCapF,0,cout);
    } else {
        if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
        calcObjFun(dbgObjFun,cout);
        for (int n=1;n<=npLnEffX;n++){
            rpt::echo<<"#----is pLnEffX["<<n<<"] active? "<<active(pLnEffX[n])<<endl;
        }
        if ((current_phase()>=phsItsRewgt)){
            PRINT2B1("#--Calculating effective weights for size compositions")
            //note that the following modifies the objective function value
            ofstream os; 
            if ((current_phase()==phsItsRewgt)){
                os.open("effectiveWeights.R", ios::trunc);
                os<<"effWgts=list("<<endl;
            } else {
                os.open("effectiveWeights.R", ios::app);
            }
            os<<"#--BETWEEN_PHASES_SECTION: Calculating effective weights for size compositions at start of "<<str(current_phase())<<endl;
            os<<"effWgts."<<current_phase()-phsItsRewgt+1<<"=list("<<endl;
            os<<"surveys="  <<endl; calcWeightsForSurveySizeComps( -1,os); os<<","<<endl;
            os<<"fisheries="<<endl; calcWeightsForFisherySizeComps(-1,os); os<<endl;
            os<<"),"<<endl;
            os<<"#--BETWEEN_PHASES_SECTION: Finished calculating effective weights for size compositions"<<endl<<endl;
            if ((ptrMOs->optIterativeReweighting)&&(numItsRewgt<maxItsRewgt)){
                PRINT2B1("#--Re-weighting size compositions using reWeightXXXSizeComps")
                os<<"#--BETWEEN_PHASES_SECTION: Re-weighting size compositions using reWeightXXXSizeComps"<<endl;
                //reWeightSurveySizeComps(-1,os);
                reWeightFisherySizeComps(-1,os);
                os<<"#--BETWEEN_PHASES_SECTION: Finished r-weighting size compositions using reWeightXXXSizeComps"<<endl<<endl;
                PRINT2B1("#--Finished re-weighting size compositions using reWeightXXXSizeComps")
                numItsRewgt++;
                calcObjFun(dbgObjFun,cout);
            }
            os.close();
        }
    }
    ctrProcCallsInPhase=0;//reset in-phase counter
        
    PRINT2B1("#---END BETWEEN_PHASES_SECTION")

//----------------------------------------------------------------------
//re-weight survey size comps
FUNCTION void reWeightSurveySizeComps(int debug,ostream& cout)
    if (debug) cout<<"#--Starting reWeightSurveySizeComps()"<<endl;
    if (ptrMOs->optIterativeReweighting){
        for (int v=1;v<=nSrv;v++){
            if (debug) cout<<"#--reweighting size comps for survey "<<ptrMC->lblsSrv[v]<<endl;
            int vd = mapM2DSrv(v);//get index for survey data corresponding to model survey v
            FleetData* ptrObs = ptrMDS->ppSrv[vd-1];
            if (ptrObs->hasICD){//index catch data
                 if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                    if (debug) cout<<"#---index catch size frequencies"<<endl;
                    ptrObs->ptrICD->ptrZFD->applyReWeightingFactors();
                }
            }
        }//surveys loop
    }
    if (debug) cout<<"#--Finished reWeightSurveySizeComps()"<<endl<<endl;

//----------------------------------------------------------------------
//re-weight fishery size comps
FUNCTION void reWeightFisherySizeComps(int debug,ostream& cout)
    if (debug) cout<<"#--Starting reWeightsForFisherySizeComps()"<<endl;
    if (ptrMOs->optIterativeReweighting){
        for (int f=1;f<=nFsh;f++){
            if (debug) cout<<"#---reweighting size comps for fishery "<<ptrMC->lblsFsh[f]<<endl;
            int fd = mapM2DFsh(f);//get index for fishery data corresponding to model fishery f
            FleetData* ptrObs = ptrMDS->ppFsh[fd-1];
            if (ptrObs->hasRCD){//retained catch data
                if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                    if (debug) cout<<"#----retained catch size frequencies"<<endl;
                    ptrObs->ptrRCD->ptrZFD->applyReWeightingFactors();
                }
            }
            if (ptrObs->hasTCD){//observed total catch data
                 if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                    if (debug) cout<<"#----total catch size frequencies"<<endl;
                    ptrObs->ptrTCD->ptrZFD->applyReWeightingFactors();
                }
            }
            if (ptrObs->hasDCD){//observed discard catch data
                if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                    if (debug) cout<<"#----discard catch size frequencies"<<endl;
                    ptrObs->ptrDCD->ptrZFD->applyReWeightingFactors();
                }
            }
        }//fisheries
    }
    if (debug) cout<<"#--Finished reWeightFisherySizeComps()"<<endl;

//-------------------------------------------------------------------------------------
//**
// * Calculate effective weights for all survey size comps
// * 
// * NOTE: As a side effect, this function modifies the value of the objective function
// * because it calls calcNLLs_CatchNatZ() to calculate the effective weights. THE USER
// * MUST RECALCULATE objFun using calcObjFun() after calling this function to obtain a valid
// * value for the objective function.
// * 
// * @param debug - integer indicating debug level (<0 writes R list to cout)
// * @param cout - output stream to write debugging info to
// * 
// * @details If ptrMOs->optIterativeReweighting > 0, then re-weighting factors for the 
// * size comps are calculated (but not applied) using the method corresponding to its value 
// * (1=McAllister-Ianelli, 2=Francis). Re-weighting is accomplished using reWeightSurveySizeComps().
// * 
// * @return none.
// */
FUNCTION void calcWeightsForSurveySizeComps(int debug, ostream& cout)
    if (debug) cout<<"#--Starting calcWeightsForSurveySizeComps()"<<endl;
    int opt = ptrMOs->optIterativeReweighting;
    adstring nDims = "c('N','McAllister-Ianelli','Francis')";
    if (debug<0) cout<<"list("<<endl;
    for (int v=1;v<=nSrv;v++){
        if (debug) cout<<"#--calculating size comps weights for survey "<<ptrMC->lblsSrv[v]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsSrv[v]<<"`=list("<<endl;
        int vd = mapM2DSrv(v);//get index for survey data corresponding to model survey v
        FleetData* ptrObs = ptrMDS->ppSrv[vd-1];
        if (ptrObs->hasICD){//index catch data
            if (debug<0) cout<<"index.catch=list("<<endl;
            if (ptrObs->ptrICD->hasZFD && ptrObs->ptrICD->ptrZFD->optFit){
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrICD->ptrZFD,n_vyxmsz(v),0,cout);//don't debug
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Index catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrICD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForSurveySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrICD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<")"<<endl;
        }//ptrObs->hasICD
        if (debug<0) cout<<"),"<<endl;
    }//surveys loop
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug) cout<<"#--Finished calcWeightsForSurveySizeComps()"<<endl;

//-------------------------------------------------------------------------------------
//**
// * Calculate effective weights for fishery size comps
// * 
// * NOTE: As a side effect, this function modifies the value of the objective function
// * because it calls calcNLLs_CatchNatZ() to calculate the effective weights. THE USER
// * MUST RECALCULATE objFun using calcObjFun() after calling this function to obtain a valid
// * value for the objective function.
// * 
// * @param debug - integer indicating debug level (<0 writes R list to cout)
// * @param cout - output stream to write debugging info to
// * 
// * @details If ptrMOs->optIterativeReweighting > 0, then re-weighting factors for the 
// * size comps are calculated (but not applied) using the method corresponding to its value 
// * (1=McAllister-Ianelli, 2=Francis). Re-weighting is accomplished using reWeightFisherySizeComps().
// * 
// * @return none.
// */
FUNCTION void calcWeightsForFisherySizeComps(int debug, ostream& cout)
    if (debug) cout<<"#--Starting calcWeightsForFisherySizeComps()"<<endl;
    int opt = ptrMOs->optIterativeReweighting;
    adstring nDims = "c('N','McAllister-Ianelli','Francis')";
    if (debug<0) cout<<"list("<<endl;
    for (int f=1;f<=nFsh;f++){
        if (debug) cout<<"#--calculating effective weights for fishery "<<ptrMC->lblsFsh[f]<<endl;
        if (debug<0) cout<<"`"<<ptrMC->lblsFsh[f]<<"`=list("<<endl;
        int fd = mapM2DFsh(f);//get index for fishery data corresponding to model fishery f
        FleetData* ptrObs = ptrMDS->ppFsh[fd-1];
        if (ptrObs->hasRCD){//retained catch data
            if (debug<0) cout<<"retained.catch=list("<<endl;
            if (ptrObs->ptrRCD->hasZFD && ptrObs->ptrRCD->ptrZFD->optFit){
                if (debug) cout<<"#---retained catch size frequencies"<<endl;
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrRCD->ptrZFD,rmN_fyxmsz(f),0,cout);
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Retained catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrRCD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForFisherySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrRCD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<"),"<<endl;
        }//ptrObs->hasRCD
        if (ptrObs->hasTCD){//observed total catch data
            if (debug<0) cout<<"total.catch=list("<<endl;
            if (ptrObs->ptrTCD->hasZFD && ptrObs->ptrTCD->ptrZFD->optFit){
                if (debug) cout<<"#---total catch size frequencies"<<endl;
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrTCD->ptrZFD,cpN_fyxmsz(f),0,cout);
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Total catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrTCD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForFisherySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrTCD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<"),"<<endl;
        }//ptrObs->hasTCD
        if (ptrObs->hasDCD){//observed discard catch data
            if (debug<0) cout<<"discard.catch=list("<<endl;
            if (ptrObs->ptrDCD->hasZFD && ptrObs->ptrDCD->ptrZFD->optFit){
                if (debug) cout<<"#---discard catch size frequencies"<<endl;
                smlVal = 0.001;//match to TCSAM2013
                d5_array effWgtComps_xmsyn = calcNLLs_CatchNatZ(ptrObs->ptrDCD->ptrZFD,dsN_fyxmsz(f),0,cout);
                if (effWgtComps_xmsyn.size()>0){
                    if (debug>0) cout<<"#--Discard catch size comps weights using "<<tcsam::getFitType(ptrObs->ptrDCD->ptrZFD->optFit)<<endl;
                    d4_array effWgts_nxms = calcEffWgts(effWgtComps_xmsyn,debug,cout);
                    if (debug>0) {cout<<"effWgts_nxms ="<<endl; wts::print(effWgts_nxms,cout,0); cout<<endl;}
                    if (debug<0) {
                        ivector bnds = wts::getBounds(effWgts_nxms);
                        adstring x = tcsamDims::getSXsForR(bnds(3),bnds(4));
                        adstring m = tcsamDims::getMSsForR(bnds(5),bnds(6));
                        adstring s = tcsamDims::getSCsForR(bnds(7),bnds(8));
                        cout<<"effWgts="<<endl; wts::writeToR(cout, effWgts_nxms,nDims,x,m,s); cout<<cc<<endl;
                    }
                    if (opt) {
                        cout<<"#----Calculating reweighting size comps in call to calcWeightsForFisherySizeComps!!"<<endl;
                        if (debug<0) {
                            cout<<"opt="<<opt<<cc<<endl;
                            cout<<"applied="<<endl;
                        }
                        ptrObs->ptrDCD->ptrZFD->calcReWeightingFactors(effWgts_nxms(opt),debug,cout);
                    }
                }//if effWgtComps_xmsyn.size()>0
            }
            if (debug<0) cout<<"),"<<endl;
        }//ptrObs->hasDCD
        if (debug<0) cout<<"NULL),"<<endl;
    }//fisheries
    if (debug<0) cout<<"NULL)"<<endl;
    if (debug) cout<<"#--Finished calcWeightsForFisherySizeComps()"<<endl;

//-------------------------------------------------------------------------------------
//**
// * Calculate effective weights for size comps
// * 
// * This function calculates effective weights for size compositions 
// * using the the McAllister-Ianelli (Punt 2017, 1.B) and 
// * Francis (Punt 2017, 1.C) methods.
// * 
// * @param effWgtComps - 5d-array (xmsyn) from call to calcNLLs_CatchNatZ()
// * @debug 
// * @cout
// * 
// * @return 4d-array (nxms) with effective weights by xms
// * 
// * @details 
// *  n=0 is the number of actual size comps processed
// *  n=1 is the harmonic mean of the individual (annual) McAllister-Ianelli tuning weights
// *  n=2 is the Francis tuning weight
// * 
// */
FUNCTION d4_array calcEffWgts(d5_array& effWgtComps,int debug, ostream& cout)
    if (debug>0) cout<<"Starting calcEffWgts()"<<endl;
    ivector d = wts::getBounds(effWgtComps);
    d4_array effWgts_nxms(0,2,d[1],d[2],d[3],d[4],d[5],d[6]);
    
    //calculate McAllister-Ianelli weights using the harmonic mean (1.B in Punt, 2017)
    for (int x=d[1];x<=d[2];x++){
        for (int m=d[3];m<=d[4];m++){
            for (int s=d[5];s<=d[6];s++){
                double N = sum(column(effWgtComps(x,m,s),1));
                double harmn = 0.0;
                for (int y=effWgtComps(x,m,s).indexmin();y<=effWgtComps(x,m,s).indexmax();y++){
                    if (effWgtComps(x,m,s,y,1)>0.0) harmn += 1.0/effWgtComps(x,m,s,y,2);
                }
                effWgts_nxms(0,x,m,s) = N;      //number of actual size comps
                if (N>0){//to avoid nan's
                    effWgts_nxms(1,x,m,s) = N/harmn;//harmonic mean  of annual Mc-I tuning weights 
                } else effWgts_nxms(1,x,m,s) = 1.0;
            }
        }
    }
    
    //calculate Francis weights  (1.C in Punt, 2017))
    for (int x=d[1];x<=d[2];x++){
        for (int m=d[3];m<=d[4];m++){
            for (int s=d[5];s<=d[6];s++){
                double N = sum(column(effWgtComps(x,m,s),1));   //number of actual size comps
                int Np   = column(effWgtComps(x,m,s),1).size(); //total number of years
                dvector effWgt_y = column(effWgtComps(x,m,s),3);//z-scores for Francis tuning weight
                effWgts_nxms(0,x,m,s) = N;
                if (N>0){//to avoid nan's
                    effWgts_nxms(2,x,m,s) = 1.0/(wts::variance(effWgt_y)*Np/N);//Francis tuning weight
                } else effWgts_nxms(2,x,m,s) = 1.0;
            }
        }
    }
    
    if (debug>0) wts::print(effWgts_nxms,cout,1);
    
    if (debug>0) cout<<"Finished calcEffWgts()"<<endl;
    return effWgts_nxms;

//-------------------------------------------------------------------------------------
FUNCTION save_params
    adstring fn = "tcsam02."+str(current_phase())+"."+str(ctrProcCallsInPhase)+".par";
    ofstream os1(fn, ios::trunc);
    os1.precision(12);
    initial_params::save(os1, 12);
    os1.close();

//-------------------------------------------------------------------------------------
//--NOTE: As a side effect, this function sets all fishery-related F's  to 0.
FUNCTION void calcDynB0(int debug, ostream& cout)
    if (debug>0) cout<<"starting calcDynB0()"<<endl;
    //open file for output
    adstring fn = "DynamicB0.R";
    ofstream os(fn, ios::trunc);
    os.precision(12);
    
    //write MMB from final phase to R
    os<<"res=list("<<endl;
    os<<"MB_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms); os<<cc<<endl;
    if (debug<0){
        cout<<"res=list("<<endl;
        cout<<"MB_yx    ="; wts::writeToR(cout,value(spB_yx), yDms,xDms); cout<<cc<<endl;
    }
    
    //initialize population model
    if (!runAlt) initPopDyMod(0, cout); else initAltPopDyMod(0,cout);
    
    //reset fishery F's to 0
    hasF_fy.initialize();   //flags indicating whether or not fishery occurs
    hmF_fy.initialize();    //handling mortality
    cpF_fyxms.initialize(); //fully-selected capture rate
    cpF_fyxmsz.initialize();//size-specific capture rate
    rmF_fyxmsz.initialize();//retained mortality rate
    dmF_fyxmsz.initialize();//discard mortality rate
    tmF_yxmsz.initialize(); //total mortality rate
    
    //run population model
    if (!runAlt){ 
        for (int y=mnYr;y<=mxYr;y++) runPopDyModOneYear(y,0,cout);
    } else {
        for (int y=mnYr;y<=mxYr;y++) 
            for (int x=1;x<=nSXs;x++) runAltPopDyModOneYear(y,x,0,cout);
    }
    
    //write MMB from dynamic B0 calculations to R
    os<<"dB0_yx    ="; wts::writeToR(os,value(spB_yx), yDms,xDms); os<<endl;
    os<<")"<<endl;
    if (debug<0){
        cout<<"dB0_yx    ="; wts::writeToR(cout,value(spB_yx), yDms,xDms); cout<<endl;
        cout<<")"<<endl;
    }
    os.close();
    
    if (debug>0) cout<<"finished calcDynB0()"<<endl;
    
//-----------------------------------    
FUNCTION void writeStateForOpMod(int y,ostream& os)
    if (!mseOpModMode){
        os<<"#--OpMod state for projecting year "<<y+1<<"--"<<endl;
        os<<mnYr    <<tb<<"#start year for recruitment"  <<endl;
        os<<y+1     <<tb<<"#year for projection"<<endl;
        os<<dtF_y(y)<<tb<<"#dtF"<<endl;
        os<<dtM_y(y)<<tb<<"#dtM"<<endl;
        os<<endl;
        os<<0<<tb<<"#number of years for TAC/OFLs"<<endl;
        os<<endl;
        os<<"#wAtZ_xmz:"<<endl; wts::print(ptrMDS->ptrBio->wAtZ_xmz,os,1); os<<endl;
        os<<"#R_y:"<<endl<<R_y(mnYr,y)<<endl;
        os<<"#R_x:"<<endl<<R_yx(y)<<endl;
        os<<"#R_z:"<<endl<<R_yz(y)<<endl;
        os<<endl;
        os<<"#M_xmsz:"   <<endl; wts::print(M_yxmsz(y),os,1);    os<<endl;
        os<<"#prGr_xszz:"<<endl; wts::print(prGr_yxszz(y),os,1); os<<endl;
        os<<"#prM2M_xz:" <<endl; wts::print(prM2M_yxz(y),os,1);  os<<endl;
        os<<"#hmF_f:"<<endl<<column(hmF_fy,y)<<endl<<endl;
        os<<"#cpF_fxmsz:"<<endl;
        for (int f=1;f<=nFsh;f++) {os<<tb<<"#"<<f<<endl; wts::print(cpF_fyxmsz(f,y),os,2); os<<endl;}
        os<<"#ret_fxmsz:"<<endl;
        for (int f=1;f<=nFsh;f++) {os<<tb<<"#"<<f<<endl; wts::print(ret_fyxmsz(f,y),os,2); os<<endl;}
        os<<"#sel_fxmsz:"<<endl;
        for (int f=1;f<=nFsh;f++) {os<<tb<<"#"<<f<<endl; wts::print(sel_fyxmsz(f,y),os,2); os<<endl;}
        os<<"#q_vxmsz:"<<endl;
        for (int v=1;v<=nSrv;v++) {os<<tb<<"#"<<v<<endl; wts::print(q_vyxmsz(v,y),os,2); os<<endl;}
        os<<"#n_xmsz (pop state at start of year to be projected):"<<endl; wts::print(n_yxmsz(y+1),os,1); os<<endl;
    } else {
        os<<"#--OpMod state for projecting year "<<ptrOMI->mxYr+1<<"--"<<endl;
        os<<ptrOMI->mnYr  <<tb<<"#start year for recruitment"<<endl;
        os<<ptrOMI->mxYr+1<<tb<<"#year for projection"       <<endl;
        os<<ptrOMI->dtF<<tb<<"#dtF"<<endl;
        os<<ptrOMI->dtM<<tb<<"#dtM"<<endl;
        os<<endl;
        os<<ptrOMI->nTACs+1<<tb<<"#number of years for TAC/OFLs"<<endl;
        os<<"#year TAC  OFL"<<endl;
        for (int y=1;y<=ptrOMI->nTACs;y++)
            os<<ptrOMI->yrsTAC(y)<<tb<<ptrOMI->TAC_y(y)<<tb<<ptrOMI->OFL_y(y)<<endl;
        os<<ptrOMI->mxYr<<tb<<inpTAC<<tb<<inpOFL<<endl;
        os<<endl;
        os<<"#wAtZ_xmz:"<<endl; wts::print(ptrOMI->wAtZ_xmz,os,1); os<<endl; 
        os<<"#R_y:"<<endl<<ptrOMI->R_y<<tb<<prjR<<endl;
        os<<"#R_x:"<<endl<<ptrOMI->R_x<<endl;
        os<<"#R_z:"<<endl<<ptrOMI->R_z<<endl;
        os<<endl;
        os<<"#M_xmsz:"   <<endl; wts::print(ptrOMI->M_xmsz,   os,1); os<<endl;
        os<<"#prGr_xszz:"<<endl; wts::print(ptrOMI->prGr_xszz,os,1); os<<endl;
        os<<"#prM2M_xz:" <<endl; wts::print(ptrOMI->prM2M_xz, os,1); os<<endl;
        os<<"#hmF_f:"    <<endl<<ptrOMI->hmF_f<<endl;
        os<<endl;
        os<<"#cpF_fxmsz:"<<endl; wts::print(ptrOMI->cpF_fxmsz,os,1); os<<endl;
        os<<"#ret_fxmsz:"<<endl; wts::print(ptrOMI->ret_fxmsz,os,1); os<<endl;
        os<<"#sel_fxmsz:"<<endl; wts::print(ptrOMI->sel_fxmsz,os,1); os<<endl;
        os<<"#q_vxmsz:"  <<endl; wts::print(ptrOMI->q_vxmsz,  os,1); os<<endl;
        os<<"#n_xmsz (pop state at start of year to be projected):"<<endl; 
          wts::print(prj_n_xmsz,os,1); os<<endl;
    }
    
//-----------------------------------    
FUNCTION void writeEstModPinFile(int closed, ostream& os)
    updateMPI(0,cout);
    ptrMPI->addNextYearToInfo(closed);
    os<<"#####---EstMod Pin File"<<endl;
    ptrMPI->writePin(os);
    
//-------------------------------------
FUNCTION int calcTAC(int hcr, double OFL)  
    int closed = 1;
    double TAC = 0.0;
    adstring info;
    if (hcr==1){
        dmatrix vspB_yx = value(spB_yx);
        double MMB = vspB_yx(mxYr,  MALE);
        double MFB = vspB_yx(mxYr,FEMALE);
        ivector perm(1,2); perm[1]=2;perm[2]=1;
        dmatrix vspB_xy = wts::permuteDims(perm,vspB_yx);
        double aveMFB = mean(vspB_xy(FEMALE)(ptrMOs->HCR1_avgMinYr,ptrMOs->HCR1_avgMaxYr));
        TAC = HarvestStrategies::HCR1_FemaleRamp(MFB, aveMFB, MMB);
        info = "#--HCR1: MFB = "+str(MFB)+cc+"aveMFB = "+str(aveMFB)+cc+"ratio = "+str(MFB/aveMFB)+cc+"TAC = "+str(TAC);
    }
    if (hcr==6){
        double Fmsy        = value(ptrOFLResults->Fmsy);
        d3_array selF_msz  = value(ptrOFLResults->pCIM->selF_fmsz(1));
        d3_array M_msz     = value(ptrOFLResults->pPDIM->M_msz);
        d3_array n_msz     = value(this->n_vyxmsz(1,mxYr,MALE));
        dmatrix w_mz       = value(ptrMDS->ptrBio->wAtZ_xmz(MALE));
        dvector cpB_z(20,32); cpB_z.initialize();
        for (int m=1;m<=nMSs;m++){
            for (int s=1;s<=nMSs;s++){
                cpB_z += elem_prod(selF_msz(m,s)(20,32),
                                   elem_prod(exp(-M_msz(m,s)(20,32)),
                                             elem_prod(n_msz(m,s)(20,32),w_mz(m)(20,32))
                                            )
                                   )*(1-exp(-Fmsy));       
            }
        }
        double CWmsy = sum(cpB_z);
    }
    if (TAC>0.0) closed=0;

    //--save TAC and OFL to file for OpMod to read
    adstring fn = "TAC_"+str(mxYr+1)+".txt";
    ofstream os; os.open(fn, ios::trunc);
    os<<"#---TAC, OFL for MSE OpMod"<<endl;
    os<<"# TAC      OFL"<<endl;
    os<<TAC<<tb<<OFL<<endl;
    os<<endl;
    os<<info<<endl;
    os.close();
    
    return(closed);
    
//----------------------
FUNCTION finishOpModMode
    PRINT2B1("#----MSE OpModMode: Recalculating population dynamics")
    if (inpTAC>0.0){
        dvariable mseCapF = mfexp(pMSE_LnC[1]);
        projectPopForTAC(mseCapF,0,cout);
        calcObjFunForTAC(dbgObjFun,cout);
        PRINT2B2("#--Final obj fun = ",objFun)
    } else {
        projectPopForZeroTAC(0,cout);
    }

    //--define parent folder for subsequent files
    ptrMC->mxYr   = mxYr+1;
    ptrMC->asYr   = ptrMC->asYr+1;
    adstring nxtOpModFolder  = "../"+str(ptrMC->asYr)+".OpMod";
    adstring nxtEstModFolder = "../"+str(ptrMC->asYr)+".EstMod";

    //--write model configuration file for OpMod in mxYr+1
    adstring nwMC = wts::concatenateFilePaths(nxtOpModFolder,
                                              "OpMod.Configuration.inp");
    PRINT2B2("writing MC to ",nwMC)
    ptrMC->write(nwMC);

    //--write model datasets file and dataset files for estimation model in mxYr+1
    //--bio data
    {
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnBioData);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ptrBio->write(ofs);
        ofs.close();
    }
    //--fishery data
    for (int i=1;i<=ptrMDS->nFsh;i++){
        if ((i!=1)||(inpTAC>0.0)){
            //update fishery i (skipping directed fishery if TAC=0)
            adstring name  = ptrMDS->ppFsh[i-1]->name;
            //update fishery data for year mxYr+1
            FleetData::debug=0;
            ptrMDS->ppFsh[i-1]->addFisheryCatchData(ptrMC->mxYr,
                                                    prj_cpN_fxmsz(i),
                                                    prj_rmN_fxmsz(i),
                                                    ptrMDS->ptrBio->wAtZ_xmz,
                                                    0.0,
                                                    0.0,
                                                    rng);
        }
        //write new fishery data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsFisheryData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppFsh[i-1]->write(ofs);
        ofs.close();
        FleetData::debug=0;
    }
    //--survey data
    //----calculate survey data for year mxYr+1
    cout<<"Calculating survey data for for year "<<ptrMC->mxYr+1<<endl;
    prj_n_vxmsz.initialize();
    for (int v=1;v<=nSrv;v++){
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++){
                    //NOTE: using survey characteristics from final year of "real" data
                    prj_n_vxmsz(v,x,m,s) = elem_prod(ptrOMI->q_vxmsz(v,x,m,s),prj_n_xmsz(x,m,s));
                }
            }
        }
    }
    cout<<"Calculated survey data for for year "<<ptrMC->mxYr+1<<endl;
    //----write survey data files
    for (int i=1;i<=ptrMDS->nSrv;i++){
        adstring name  = ptrMDS->ppSrv[i-1]->name;
        //update survey data for year mxYr+1
        FleetData::debug=0;
        ptrMDS->ppSrv[i-1]->addIndexCatchData(ptrMC->mxYr+1,
                                              prj_n_vxmsz(i),
                                              ptrMDS->ptrBio->wAtZ_xmz,
                                              0.0,
                                              0.0,
                                              rng);

        //write new survey data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsSurveyData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppSrv[i-1]->write(ofs);
        ofs.close();
        FleetData::debug=0;
    }
    //--growth data
    for (int i=1;i<=ptrMDS->nGrw;i++){
        adstring name  = ptrMDS->ppGrw[i-1]->name;
        //update growth data for year mxYr+1
        //NOTHING TO UPDATE

        //write new growth data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsGrowthData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppGrw[i-1]->write(ofs);
        ofs.close();
    }
    //--chela height data
    for (int i=1;i<=ptrMDS->nCHD;i++){
        adstring name  = ptrMDS->ppCHD[i-1]->name;
        //update chela height data for year mxYr+1
        //NOTHING TO UPDATE

        //write new chela height data
        adstring nwF = wts::concatenateFilePaths(nxtEstModFolder,ptrMDS->fnsChelaHeightData[i]);
        ofstream ofs; ofs.open(nwF,ios::trunc);
        ptrMDS->ppCHD[i-1]->write(ofs);
        ofs.close();
    }
    adstring nwMDS = wts::concatenateFilePaths(nxtEstModFolder,ptrMC->fnMDS);
    rpt::echo<<"writing MDS to '"<<nwMDS<<"'"<<endl;
    ptrMDS->write(nwMDS);

    //write OpMod "state" to file
    {
    adstring nwF = "OpModStateFile_"+str(mxYr+2)+".txt";//writing out for July 1, mxYr+2 (yes, this is correct!)
    ofstream ofs; ofs.open(nwF,ios::trunc);
    writeStateForOpMod(0,ofs);
    ofs.close();
    }

// =============================================================================
// =============================================================================
FINAL_SECTION
    PRINT2B1(" ")
    PRINT2B1("#--Starting FINAL_SECTION")
        
    if (!mseMode){
        PRINT2B1("#--mseMode is OFF")
        if (mcevalOn) {
            PRINT2B1("#--mceval is ON")
            PRINT2B1("#----Closing mcmc file")
            mcmc.open((char*)(fnMCMC),ios::app);
            mcmc.precision(12);
            //mcmc<<"NULL)"<<endl;
            mcmc.close();
            PRINT2B1(" ")
        }

        if (!mcevalOn){
            PRINT2B1("#--mceval is OFF")
            PRINT2B2("obj fun = ",objFun)
            PRINT2B1(" ")
            {
                PRINT2B1("#----Writing sim data to file")
                ofstream echo1; echo1.open("ModelSimData.dat", ios::trunc);
                echo1.precision(12);
                writeSimData(echo1,0,cout,ptrSimMDS);
                PRINT2B2("obj fun = ",objFun)
                PRINT2B1(" ")
            }

            {
                PRINT2B1("#----Calculating final effective weights for size compositions")
                PRINT2B1("#----Note that the value of objFun is not valid now!!")
                //note that this modifies the value of the objective function!
                ofstream os; os.open("effectiveWeights.R", ios::app);
                os<<"#--Calculating effective weights in FINAL_PHASE--"<<endl;
                os<<"effWgts.final=list("<<endl;
                os<<"surveys=";   calcWeightsForSurveySizeComps( -1,os); os<<","<<endl;
                os<<"fisheries="; calcWeightsForFisherySizeComps(-1,os); os<<endl;
                os<<")"<<endl;
                os<<"#--Finished calculating effective weights in FINAL_PHASE--"<<endl;
                os<<")"<<endl;
                os.close();
                PRINT2B2("#--obj fun = ",objFun)
                PRINT2B1(" ")
            }

            if (doDynB0>0){
                PRINT2B1("#----Calculating dynamic B0")
                PRINT2B1("#----Note that the value of objFun is not valid now!!")
                calcDynB0(1,rpt::echo);
                PRINT2B2("#--obj fun = ",objFun)
                PRINT2B1(" ")
            }

            PRINT2B1("#----Recalculating final objective function value")
            if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
            calcObjFun(dbgObjFun,cout);
            PRINT2B2("#--Final obj fun = ",objFun)
                    
            {
                PRINT2B1("#----writing cohort progression to R")
                ofstream echo1; echo1.precision(12);
                echo1.open("CohortProgression.noMF.R", ios::trunc);                
                ReportToR_CohortProgression(echo1,1,0,0,1,0,cout);
                echo1.close();
                echo1.open("CohortProgression.FnoM.R", ios::trunc);                
                ReportToR_CohortProgression(echo1,1,0,1,1,0,cout);
                echo1.close();
                echo1.open("CohortProgression.MnoF.R", ios::trunc);                
                ReportToR_CohortProgression(echo1,1,1,0,1,0,cout);
                echo1.close();
                echo1.open("CohortProgression.MF.R", ios::trunc);                
                ReportToR_CohortProgression(echo1,1,1,1,1,0,cout);
                echo1.close();
                PRINT2B1("#----finished writing cohort progression to file")
            }
            
           //do OFL calculations
            if (doOFL){
                PRINT2B1("#----Starting OFL calculations")
                ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
                echoOFL.precision(12);
                calcOFL(mxYr+1,1,echoOFL);//updates ptrOFLResults
                ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
                ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
                echoOFL.close();
                PRINT2B1("#----Finished OFL calculations")
            }
            
            //do TAC calculations (doTAC = HCR number)
            PRINT2B2("#----doTAC = ",doTAC)
            if (doTAC>0) {
                PRINT2B1("#----Starting TAC calculations")
                int closed = calcTAC(doTAC, value(ptrOFLResults->OFL));
                PRINT2B1("#----Finished TAC calculations")
                
                //write state for operating model
                {
                    PRINT2B1("#----Writing OpMod state file")
                    adstring nwF = "OpModStateFile_"+str(mxYr+1)+".txt";
                    //nwF = wts::concatenateFilePaths(parent,nwF);
                    ofstream ofs; ofs.open(nwF,ios::trunc);
                    writeStateForOpMod(mxYr,ofs);
                    ofs.close();
                    PRINT2B1("#----Finished writing OpMod state file")
                }
            
                //write MPI for EstMod for upcoming year
                //----devs for upcoming year are zero
                {
                    ptrMPI->addNextYearToInfo(closed);
                    ptrMPI->setToWriteVectorInitialValues(false);
                    adstring nwF = "EstMod.ParametersInfo."+str(mxYr+1)+".inp";
                    ofstream ofs; ofs.open(nwF, ios::trunc);
                    ptrMPI->write(ofs);
                    ofs.close();
                }

                //--write pin for EstMod for upcoming year
                {
                    adstring nwF = "EstModPinFile_"+str(mxYr+1)+".txt";
                    ofstream ofs; ofs.open(nwF, ios::trunc);
                    ptrMPI->writePin(ofs);
                    ofs.close();
                }
            }//doTAC>0
        }//!mcevalOn
    } else if (mseEstModMode){
        //running in estModMode
        PRINT2B1("#----MSE EstModMode: Recalculating population dynamics")
        if (!runAlt) runPopDyMod(0,cout); else runAltPopDyMod(0,cout);
        calcObjFun(dbgObjFun,cout);
        PRINT2B2("#--Final obj fun = ",objFun)
        
        //--calculate OFL
        {
            cout<<"#----Starting OFL calculations"<<endl;
            ofstream echoOFL; echoOFL.open("calcOFL.final.txt", ios::trunc);
            echoOFL.precision(12);
            calcOFL(mxYr+1,1,echoOFL);//updates ptrOFLResults
            ptrOFLResults->writeCSVHeader(echoOFL); echoOFL<<endl;
            ptrOFLResults->writeToCSV(echoOFL);     echoOFL<<endl;
            echoOFL.close();
            cout<<"#----Finished OFL calculations"<<endl;
        }
        
        //--calculate TAC for upcoming year using harvest control rule
        int closed = calcTAC(doTAC, value(ptrOFLResults->OFL));
                
        //write MPI for EstMod for upcoming year
        //----devs for upcoming year are zero
        {
            ptrMPI->addNextYearToInfo(closed);
            ptrMPI->setToWriteVectorInitialValues(false);
            adstring nwF = "EstMod.ParametersInfo."+str(mxYr+1)+".inp";
            ofstream ofs; ofs.open(nwF, ios::trunc);
            ptrMPI->write(ofs);
            ofs.close();
        }

        //--write pin for EstMod for upcoming year
        {
            adstring nwF = "EstModPinFile_"+str(mxYr+1)+".txt";
            ofstream ofs; ofs.open(nwF, ios::trunc);
            ptrMPI->writePin(ofs);
            ofs.close();
        }
                
    } else if (mseOpModMode){
        finishOpModMode();
    }//mseOpModMode
    
    int hour,minute,second;
    double elapsed_time;
    
    time(&finish); 
    elapsed_time = difftime(finish,start);
    
    hour   = (int) long(elapsed_time)/3600;
    minute = (int) long(elapsed_time)%3600/60;
    second = (int) (long(elapsed_time)%3600)%60;
    PRINT2B1("")
    PRINT2B1("#------------------")
    PRINT2B2("Starting time : ",ctime(&start))
    PRINT2B2("Finishing time: ",ctime(&finish))
    adstring hms = "This run took: "+str(hour)+" hours, "+str(minute)+" minutes, "+str(second)+" seconds.";
    PRINT2B1(hms)
    PRINT2B1("#--Finished FINAL_SECTION")
    PRINT2B1("#------------------------")
    PRINT2B1("#------------------------")
    PRINT2B1("")
        
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
  arrmblsize = 2147000000; //must be smaller than 2,147,483,647
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(40000000); // this may be incorrect in the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1500000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(7000);
  time(&start);

