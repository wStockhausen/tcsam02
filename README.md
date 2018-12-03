#TCSAM02
##Introduction
TCSAM02 is the current (Septemtber 2018) modeling framework for the 
Bering Sea and Aleutian Islands Tanner crab assessment. It supersedes the
TCSAM2013 modeling framework, which was used for the Tanner crab assessment from 
September 2014 to May, 2017. The two modeling frameworks were demonstrated to
produce identical ("exactly equivalent") results at the May 2017 Crab Plan Team (CPT)
meeting in Juneau. The TCSAM02 framework was subsequently recommended by the CPT, and 
approved by the North Pacific Fishery Management Council's Science and Statistical Committee,
for use in the September 2017 Tanner crab assessment and has been used in subsequent assessments.


##Required libraries
The admb and wtsADMB C++ libraries are required to compile TCSAM02. The relevant ADMB library can be found at http://www.admb-project.org. The
current version of TCSAM02 uses ADMB version 12. wtsADMB is a library of ADMB-RELATED C++ functions available as soure code on GitHub
at https://github.com/wStockhausen/wtsADMB. 


##Setting up TCSAM02 as a Netbeans Project
Netbeans should be configured with the C++ modules installed. 
1. Use Team/git/clone to clone the repository into TopFolder/tcsam02, where "TopFolder" is an arbitrary directory.
2. Create a new Netbeans C++ application called "tcsam02" in the directory TopFolder (so TopFolder/tcsam02 will be
the top-level folder in the Netbeans project). Using the Project Creator:
    a. create a "new" C++ application project in TopFolder
    b. uncheck the "Create Main File" option
    c. use "tcsam02" as the project name
    d. save the project
3. Open the "tcsam02" project in Netbeans.
    a. add the files under the "include"subfolder to the "Header Files" (right click the icon, select "Add existing file")
    b. add the files under the "src" subfolder to the "Source Files"
    c. add tcsam02.tpl to the "Important Files"
    d. add the cloned Makefile to the "Important Files"
4. Right-click the "tcsam02" project icon, select "Properties", and set up  the compiler and linker options per your system.


##TCSAM02 commandline options
* -configFile filename : specify full path to model configuration file
* -resultsFile filename : specify full path to model results (.rep?) file
* -pin filename : specify pin file for initial parameter values 
* -mcpin filename : specify pin file for running NUTS mcmc
* -calcOFL : flag to calculate OFL-related values after final phase
* -calcDynB0 : flag to calculate dynamic B0 after final phase
* -opModMode : flag to run in operating model mode (run model without parameter estimation)
* -mseMode type : flag to run model in MSE 'type' mode (type = 'mseOpModMode' or 'mseEstModMode')
* -doRetro yRetro : do retrospective model run for year assYr-yRetro
* -fitSimData iSimDataSeed : fit simulated data with random error based on a random number seed 
* -ctrDebugParams ctr : flag to print diagnostic information related to parameters starting at counter ctr
* -jitter iSeed : flag to jitter initial parameter values using a random number seed
* -resample : flag to turn on resampling for initial parameter values
* -debugModelConfig : flag to print diagnostic information related to OFL calculations
* -debugModelDatasets : flag to print diagnostic information related to OFL calculations
* -debugModelParamsInfo : flag to print diagnostic information related to OFL calculations
* -debugModelOptions : flag to print diagnostic information related to OFL calculations
* -debugModelParams : flag to print diagnostic information related to model parameters
* -debugDATA_SECTION : flag to print diagnostic information related to the DATA_SECTION
* -debugPARAMS_SECTION : flag to print diagnostic information related to the PARAMETERS_SECTION
* -debugPRELIM_CALCS : flag to print diagnostic information related to the PRELIMINARY_CALCS_SECTION
* -debugPROC_SECTION : flag to print diagnostic information related to the PROCEDURE_SECTION
* -debugREPORT_SECTION : flag to print diagnostic information related to the REPORT_SECTION
* -debugRunModel : flag to print diagnostic information related to running the model
* -debugObjFun : flag to print diagnostic information related to calculating the objective function
* -debugMCMC : flag to print diagnostic information related to an MCMC run
* -debugOFL : flag to print diagnostic information related to OFL calculations
* -showActiveParams : flag to print information related to the active parameters

