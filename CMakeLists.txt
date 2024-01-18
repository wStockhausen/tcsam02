cmake_minimum_required(VERSION 3.26)

project(
    tcsam02
    LANGUAGES CXX
    VERSION 2017.06.06
    DESCRIPTION "CMake script to build tcsam02, the assessment model for Tanner crab"
)

#--set some local variables
set(ADMB_HOME "/Users/williamstockhausen/Work/Programming/admb/build/admb")
set(WTSADMB_HOME "/Users/williamstockhausen/Work/StockAssessments-Crab/AssessmentModelDevelopment/VS_Code/wtsADMB")
set(TPL "tcsam02")
set(TPL_BASE "../${TPL}.tpl")
set(${TPL}_INCLUDE "include")
set(${TPL}_SRC "src")
set(myCXX_FLAGS "-DUSE_ADMB_CONTRIBS")

#--add external executables; set associated locations
add_executable(tpl2cpp IMPORTED)
set_property(
    TARGET tpl2cpp
    PROPERTY IMPORTED_LOCATION 
        "${ADMB_HOME}/bin/tpl2cpp"
)
add_executable(tpl2rem IMPORTED)
set_property(
    TARGET tpl2rem 
    PROPERTY IMPORTED_LOCATION 
        "${ADMB_HOME}/bin/tpl2rem"
)

#--add external libraries
add_library(libadmb STATIC IMPORTED)
set_property(
    TARGET libadmb 
    PROPERTY IMPORTED_LOCATION 
        "${ADMB_HOME}/lib/libadmb-contrib.a"
)
add_library(libwtsADMB STATIC IMPORTED)
set_property(
    TARGET libwtsADMB 
    PROPERTY IMPORTED_LOCATION 
        "${WTSADMB_HOME}/_build/libwtsADMB.a"
)

#--add "global" include folders
#----relative paths are interpreted as relative to the current source folder
include_directories(
    "_build"
    "${WTSADMB_HOME}/include"
    "${ADMB_HOME}/include"
    "${ADMB_HOME}/include/contrib"
    "${${TPL}_INCLUDE}"
)

#--add command to create .ccp and .htp files from tpl target (TODO?: add logic for tpl2rem)
#----copy tpl to _build folder
add_custom_command(
    OUTPUT ${TPL}.tpl
    COMMAND cp  ${TPL_BASE} ${TPL}.tpl        
    DEPENDS "${TPL_BASE}"
            "${TPL_PERSONAL}"
    COMMENT "Copying ${TPL}.tpl from ${TPL_BASE}."
)
#----run tpl2cpp in _build folder to create ${TPL}.cpp and ${TPL}.htp files there
add_custom_command(
    OUTPUT ${TPL}.cpp ${TPL}.htp
    COMMAND tpl2cpp "${TPL}"
    DEPENDS "${TPL}.tpl"
    COMMENT "Creating ${TPL}.cpp and ${TPL}.htp from ${TPL}.tpl using tpl2cpp."
)

#--add command to compile and create executable
add_executable(${TPL} ${TPL}.cpp ${TPL}.htp)
#--add _build folder to target-specific include folder
# target_include_directories(${TPL}
#     "."
#     "../${${TPL}_INCLUDE}"
# )
#--add preprocessor flags (no -D) for target
target_compile_definitions(${TPL} 
    PUBLIC USE_ADMB_CONTRIBS
    PUBLIC _USE_MATH_DEFINES
)
#--add compiler flags for target
target_compile_options(${TPL} 
    PUBLIC -Wall
    PUBLIC -O3
)
#--add compiler features
target_compile_features(${TPL} PUBLIC cxx_std_14)
#--add target include directories 
#--NOTE: (${GMACS_INCLUDE} already in global include_directories)
# target_include_directories(${TPL} PUBLIC xxx)
#--add additional sources to target 
#--TODO?: could move this to CMakeLists file in ${GMACS_SRC}??
target_sources(${TPL} 
    PUBLIC "${${TPL}_SRC}/CatchData.cpp"
    PUBLIC "${${TPL}_SRC}/HarvestStrategies.cpp"
    PUBLIC "${${TPL}_SRC}/IndivData.cpp"
    PUBLIC "${${TPL}_SRC}/MSEClasses.cpp"
    PUBLIC "${${TPL}_SRC}/ModelConfiguration.cpp"
    PUBLIC "${${TPL}_SRC}/ModelConstants.cpp"
    PUBLIC "${${TPL}_SRC}/ModelData.cpp"
    PUBLIC "${${TPL}_SRC}/ModelFunctions.cpp"
    PUBLIC "${${TPL}_SRC}/ModelIndexBlocks.cpp"
    PUBLIC "${${TPL}_SRC}/ModelIndexFunctions.cpp"
    PUBLIC "${${TPL}_SRC}/ModelOptions.cpp"
    PUBLIC "${${TPL}_SRC}/ModelParameterFunctions.cpp"
    PUBLIC "${${TPL}_SRC}/ModelParameterInfoTypes.cpp"
    PUBLIC "${${TPL}_SRC}/ModelParameterInfoVectorTypes.cpp"
    PUBLIC "${${TPL}_SRC}/ModelParametersInfo.cpp"
    PUBLIC "${${TPL}_SRC}/ModelPopDyClasses.cpp"
    PUBLIC "${${TPL}_SRC}/ModelSelectivities.cpp"
    PUBLIC "${${TPL}_SRC}/OFLCalcs.cpp"
    PUBLIC "${${TPL}_SRC}/SummaryFunctions.cpp"
)
#--add libraries to link for target 
#----(not sure this is necessary since these are also "global" libraries defined above)
target_link_libraries(${TPL} 
    PUBLIC libwtsADMB
    PUBLIC libadmb
)