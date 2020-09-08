# TCSAM02

## Introduction
TCSAM02 (Tanner Crab Stock Assessment Model, version 2) is the current (September 2020) 
modeling framework for the Bering Sea and Aleutian Islands Tanner crab assessment. It supersedes the
TCSAM2013 modeling framework, which was used for the Tanner crab assessment from 
September 2014 to May, 2017. The two modeling frameworks were demonstrated to
produce identical ("exactly equivalent") results at the May 2017 Crab Plan Team (CPT)
meeting. The TCSAM02 framework was subsequently recommended by the CPT, and 
approved by the North Pacific Fishery Management Council's Science and Statistical Committee,
for use in the September 2017 Tanner crab assessment and has been used in subsequent assessments.


## Required libraries
The admb and wtsADMB C++ libraries are required to compile TCSAM02. The relevant ADMB library can be found at http://www.admb-project.org. The
current version of TCSAM02 uses ADMB version 12. wtsADMB is a library of ADMB-RELATED C++ functions available as soure code on GitHub
at https://github.com/wStockhausen/wtsADMB. 


## Setting up TCSAM02 as a Netbeans Project
Netbeans should be configured with the C++ modules installed. 
1. Use Team/git/clone to clone the repository into TopFolder/tcsam02, where "TopFolder" is an arbitrary directory.
2. Create a new Netbeans C++ application called "tcsam02" in the directory TopFolder (so TopFolder/tcsam02 will be
the top-level folder in the Netbeans project). Using the Project Creator:
    * create a "new" C++ application project in TopFolder
    * uncheck the "Create Main File" option
    * use "tcsam02" as the project name
    * save the project
3. Open the "tcsam02" project in Netbeans.
    * edit nbproject/configurations.xml to replace "tcsam02-Makefile.mk" with "Makefile" in both instances
    * delete "tcsam02-Makefile.mk"
    * Note that "Makefile" is now under the "Important Files" icon in the Project view
    * edit the Makefile variables PLATFORM and ADMB_HOME_WIN or ADMB_HOME_MAC under "#Environment" to conform to your ADMB installation
    * add tcsam02.tpl to the "Important Files"
    * right-click the Makefile under "Important Files", Select "Make Target" and run the .tpl2cpp target (may need to "Add" it first)
    * assuming the tpl2cpp ran successfully (if not, check your path to tpl2cpp.exe), the files tcsam02.htp and tcsam02.cpp will have been created
        * add tcsam02.htp to the "Header Files" (right click the icon, select "Add existing file")
        * add tcsam02.cpp to the "Source Files" (right click the icon, select "Add existing file")
    * add the files under the "include" subfolder to the "Header Files" (right click the icon, select "Add existing file")
    * add the files under the "src" subfolder to the "Source Files" (right click the icon, select "Add existing file")
4. Right-click the "tcsam02" project icon, select "Properties", and set up  the compiler and 
linker options per your system.
    * Windows (64-bit ADMB, MinGW toolset, for non-optimized compilation)
        * Under "Build/C++ Compiler", 
            * "Include Directories": add ".", "include", and the paths to the wtsADMB include directory, the ADMB include directory, and the ADMB include/contrib directory
            * "Additional Options": add "-std=c++11 -O3 -D_FILE_OFFSET_BITS=64" (no quotes)
        * Under Linker", file paths under "Libraries" should point to the libwtsadmb.a and libadmb-contrib.a libraries
    * OSX (64-bit ADMB, GNU toolset, for non-optimized compilation)
        * Under "Build/C++ Compiler", 
            * "Include Directories": add ".", "include", and the paths to the wtsADMB include directory, the ADMB include directory, and the ADMB include/contrib directory
            * "Additional Options": none
        * Under Linker", 
            * file paths under "Libraries" should point to the libwtsadmb.a and libadmb-contrib.a libraries
            * "Additional Options": add "-g" (no quotes) to include debugging symbols

## TCSAM02 commandline options
### file options
* -configFile filename : flag to specify full path to model configuration file
* -pin filename : flag to specify pin file for initial parameter values 
* -binp fnPin : flag to use binary pin file fnPin to set parameter values
* -ainp fnPin : flag to use ascii  pin file fnPin to set parameter values
* -mcpin filename : flag to specify pin file for running NUTS mcmc

### run configuration options
* -mceval : flag to run in mceval mode (i.e., mcevalOn=1)
* -runAlt : flag to run alternative pop dy equations
* -calcOFL : flag to calculate OFL-related values after final phase
* -calcTAC hcr : flag to calculate TAC using the indicated harvest control rule (hcr)
* -calcDynB0 : flag to calculate dynamic B0 after final phase
* -doRetro yRetro : flag to do retrospective model run for year assYr-yRetro
* -fitSimData iSimDataSeed : flag to fit simulated data with random error based on a random number seed 
* -jitter iSeed : flag to jitter initial parameter values using a random number seed
* -resample : flag to turn on resampling for initial parameter values

### commandline flags to print debugging info
* -ctrDebugParams ctr : flag to print diagnostic information related to parameters starting at counter ctr
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

