#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=MacOSX-ADMBv12
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/CatchData.o \
	${OBJECTDIR}/src/HarvestStrategies.o \
	${OBJECTDIR}/src/IndivData.o \
	${OBJECTDIR}/src/MSEClasses.o \
	${OBJECTDIR}/src/ModelConfiguration.o \
	${OBJECTDIR}/src/ModelConstants.o \
	${OBJECTDIR}/src/ModelData.o \
	${OBJECTDIR}/src/ModelFunctions.o \
	${OBJECTDIR}/src/ModelIndexBlocks.o \
	${OBJECTDIR}/src/ModelIndexFunctions.o \
	${OBJECTDIR}/src/ModelOptions.o \
	${OBJECTDIR}/src/ModelParameterFunctions.o \
	${OBJECTDIR}/src/ModelParameterInfoTypes.o \
	${OBJECTDIR}/src/ModelParametersInfo.o \
	${OBJECTDIR}/src/ModelPopDyClasses.o \
	${OBJECTDIR}/src/ModelSelectivities.o \
	${OBJECTDIR}/src/OFLCalcs.o \
	${OBJECTDIR}/src/SummaryFunctions.o \
	${OBJECTDIR}/tcsam02.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/dist/MacOSX-ADMBv12/GNU-MacOSX/libwtsadmb.a /Users/WilliamStockhausen/Programs/admb_v12/build/dist/lib/libadmb-contrib.a /Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/dist/MacOSX-ADMBv12/GNU-MacOSX/libwtsadmb.a /Users/WilliamStockhausen/Programs/admb_v12/build/dist/lib/libadmb-contrib.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02: /Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/dist/MacOSX-ADMBv12/GNU-MacOSX/libwtsadmb.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02: /Users/WilliamStockhausen/Programs/admb_v12/build/dist/lib/libadmb-contrib.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02: /Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/dist/MacOSX-ADMBv12/GNU-MacOSX/libwtsadmb.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02: /Users/WilliamStockhausen/Programs/admb_v12/build/dist/lib/libadmb-contrib.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02 ${OBJECTFILES} ${LDLIBSOPTIONS} -g

${OBJECTDIR}/src/CatchData.o: src/CatchData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/CatchData.o src/CatchData.cpp

${OBJECTDIR}/src/HarvestStrategies.o: src/HarvestStrategies.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/HarvestStrategies.o src/HarvestStrategies.cpp

${OBJECTDIR}/src/IndivData.o: src/IndivData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/IndivData.o src/IndivData.cpp

${OBJECTDIR}/src/MSEClasses.o: src/MSEClasses.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/MSEClasses.o src/MSEClasses.cpp

${OBJECTDIR}/src/ModelConfiguration.o: src/ModelConfiguration.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelConfiguration.o src/ModelConfiguration.cpp

${OBJECTDIR}/src/ModelConstants.o: src/ModelConstants.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelConstants.o src/ModelConstants.cpp

${OBJECTDIR}/src/ModelData.o: src/ModelData.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelData.o src/ModelData.cpp

${OBJECTDIR}/src/ModelFunctions.o: src/ModelFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelFunctions.o src/ModelFunctions.cpp

${OBJECTDIR}/src/ModelIndexBlocks.o: src/ModelIndexBlocks.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelIndexBlocks.o src/ModelIndexBlocks.cpp

${OBJECTDIR}/src/ModelIndexFunctions.o: src/ModelIndexFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelIndexFunctions.o src/ModelIndexFunctions.cpp

${OBJECTDIR}/src/ModelOptions.o: src/ModelOptions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelOptions.o src/ModelOptions.cpp

${OBJECTDIR}/src/ModelParameterFunctions.o: src/ModelParameterFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelParameterFunctions.o src/ModelParameterFunctions.cpp

${OBJECTDIR}/src/ModelParameterInfoTypes.o: src/ModelParameterInfoTypes.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelParameterInfoTypes.o src/ModelParameterInfoTypes.cpp

${OBJECTDIR}/src/ModelParametersInfo.o: src/ModelParametersInfo.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelParametersInfo.o src/ModelParametersInfo.cpp

${OBJECTDIR}/src/ModelPopDyClasses.o: src/ModelPopDyClasses.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelPopDyClasses.o src/ModelPopDyClasses.cpp

${OBJECTDIR}/src/ModelSelectivities.o: src/ModelSelectivities.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ModelSelectivities.o src/ModelSelectivities.cpp

${OBJECTDIR}/src/OFLCalcs.o: src/OFLCalcs.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/OFLCalcs.o src/OFLCalcs.cpp

${OBJECTDIR}/src/SummaryFunctions.o: src/SummaryFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/SummaryFunctions.o src/SummaryFunctions.cpp

${OBJECTDIR}/tcsam02.o: tcsam02.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I. -Iinclude -I/Users/WilliamStockhausen/StockAssessments-Crab/AssessmentModelDevelopment/wtsADMB/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/include -I/Users/WilliamStockhausen/Programs/admb_v12/build/dist/contrib/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/tcsam02.o tcsam02.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tcsam02

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
