#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelFunctions.hpp"
#include "ModelIndexFunctions.hpp"
#include "ModelParameterInfoTypes.hpp"
#include "ModelParametersInfo.hpp"
#include "ModelSelectivities.hpp"

/*----------------------------------------------------------------------------\n
**  Includes:
**      ParametersGroupInfo\n
**      RecruitmentInfo
**      NaturalMortalityInfo\n
**      GrowthInfo\n
**      MaturityInfo\n
**      SelectivityInfo\n
**      FisheriesInfo\n
**      SurveysInfo\n
**      ModelParametersInfo\n
**----------------------------------------------------------------------------*/
int ParameterGroupInfo::debug   = 0;
int RecruitmentInfo::debug      = 0;
int NaturalMortalityInfo::debug = 0;
int GrowthInfo::debug           = 0;
int Molt2MaturityInfo::debug    = 0;
int SelectivityInfo::debug      = 0;
int FisheriesInfo::debug        = 0;
int SurveysInfo::debug          = 0;
int MSE_Info::debug             = 0;
int ModelParametersInfo::debug  = 0;
const adstring ModelParametersInfo::version = "2020.02.02";
    
/*----------------------------------------------------------------------------*/
/**
 * Function to convert parameter combinations to R dimensions 
 * @param pgi
 * @return wts::adstring_matrix
 */
wts::adstring_matrix tcsam::convertPCs(ParameterGroupInfo * pgi){
    //cout<<"starting tcsam::convertPCs for "<<pgi->name<<endl;
    int nIVs = pgi->nIVs;
    int nPCs = pgi->nPCs;
    //cout<<"nIVs = "<<nIVs<<"; nPCs = "<<nPCs<<endl;
    //wts::adstring_matrix::debug = 1;
    wts::adstring_matrix a(1,nIVs,0,nPCs);
    //wts::adstring_matrix::debug = 0;
    //cout<<"allocated a"<<endl;
    int ibsIdx = 1;
    for (int i=1;i<=nIVs;i++){
        if (pgi->lblIVs(i)==tcsam::STR_SEX)             {a(i,0)=pgi->lblIVs(i);} else
        if (pgi->lblIVs(i)==tcsam::STR_MATURITY_STATE)  {a(i,0)=pgi->lblIVs(i);} else
        if (pgi->lblIVs(i)==tcsam::STR_SHELL_CONDITION) {a(i,0)=pgi->lblIVs(i);} else 
        if ((pgi->ibsIdxs.allocated())&&(i==pgi->ibsIdxs(ibsIdx))){
            a(i,0)=pgi->ppIBSs[ibsIdx-1]->getType();
            if (ibsIdx<pgi->nIBSs) ibsIdx++;//increment to next
        }
    }
    for (int r=1;r<=nPCs;r++){//loop over rows
        //cout<<"r = "<<r<<endl;
        int ibsIdx = 1;
        for (int i=1;i<=nIVs;i++){//loop over index variables
            //cout<<tb<<"i = "<<i<<" "<<pgi->lblIVs(i)<<tb<<pgi->in(r,i)<<endl;
            if (pgi->lblIVs(i)==tcsam::STR_SEX)             {a(i,r)=tcsam::getSexType(pgi->in(r,i));}      else
            if (pgi->lblIVs(i)==tcsam::STR_MATURITY_STATE)  {a(i,r)=tcsam::getMaturityType(pgi->in(r,i));} else
            if (pgi->lblIVs(i)==tcsam::STR_SHELL_CONDITION) {a(i,r)=tcsam::getShellType(pgi->in(r,i));}    else 
            if ((pgi->ibsIdxs.allocated())&&(i==pgi->ibsIdxs(ibsIdx))){
                a(i,r)=pgi->ppIBSs[ibsIdx-1]->getIndexBlock(r)->asString();
                if (ibsIdx<pgi->nIBSs) ibsIdx++;//increment to next
            } else {a(i,r)=str(pgi->in(r,i));}
        }
    }
    //cout<<a<<endl;
    //cout<<"finished tcsam::convertPCs for "<<pgi->name<<endl;
    return a;
}

/**
 * Adjusts parameter phases to be consistent with maxYr.
 * 
 * The phase for any parameter that is originally set to be estimated ONLY in a 
 * time block later than maxYr is set to negative so it will not be estimated.
 * 
 * @param pI - pointer to NumberVectorInfo for vector of parameters of interest
 * @param pgi - pointer to ParameterGroupInfo for the parameter group of interest
 * @param pid - index to parameter in ParameterGroupInfo
 * @param maxYr - max year for estimated parameters
 * @param debug - flag to print debugging info
 */
void tcsam::adjustParamPhase(NumberVectorInfo* pI, ParameterGroupInfo* pgi, int pid, int maxYr, int debug){
    if (debug) rpt::echo<<"Starting adjustParamPhase(NumberVectorInfo*...) for "<<pI->name<<" for max year "<<maxYr<<" using pid "<<pid<<endl;
    int np = pI->getSize();//number of parameters
    ivector phases = pI->getPhases();
    if (debug) rpt::echo<<"--initial phases: "<<phases<<endl;
    ivector done(1,np); done=0;
    for (int p=1;p<=np;p++){ //loop over parameters
        if (debug) rpt::echo<<"--checking parameter "<<p<<endl;
        if ((done[p]==0)&&(phases[p]>0)){ //check only parameters that haven't been determined and are set to be estimated
            int nPCs = pgi->nPCs;//number of parameter combinations in the associated parameter group
            bool estParam = false; 
            for (int pc=1;pc<=nPCs;pc++){ //loop over parameter combinations
                ivector pcids = pgi->getPCIDs(pc);
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter "<<p<<". and pcid = "<<pcids<<endl;
                int pcid = pcids[pid];//get id (index) of parameter used in this pc
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter "<<p<<" and got pcid "<<pcid<<endl;
                if (pcid==p){//make check for this parameter, since it's included in this pc
                    ivector years = pgi->getYearsForPC(pc);
                    if (debug) rpt::echo<<"----pc = "<<pc<<tb<<"years = "<<years<<endl;
                    for (int iy=years.indexmin();iy<=years.indexmax();iy++){
                        if (years[iy]<=maxYr){
                            //parameter is in pc that pertains to model year interval
                            estParam=true;//so parameter should be estimated
                            done[p]=1;//don't need to check further pc's
                            if (debug) rpt::echo<<"----estParam=TRUE and done=1"<<endl;
                            break;
                        }
                    }//--iy loop
                    if (estParam) break;
                }//--pcid==p
            }//--loop over pc's
            if (debug && !estParam) rpt::echo<<"----estParam=FALSE"<<endl;
            phases[p] = estParam ? phases[p] : -phases[p];//set phase to negative if parameter should NOT be estimated
        }//--(!done[p])&&(phases[p]>0)
    }//--loop over parameters
    if (debug) rpt::echo<<"--final phases: "<<phases<<endl;
    pI->setPhases(phases);
    if (debug) rpt::echo<<"Finished adjustParamPhase(NumberVectorInfo*...) for "<<pI->name<<" for max year "<<maxYr<<endl<<endl;
}

/**
 * Adjusts parameter vector phases to be consistent with maxYr.
 * 
 * Estimation phases are defined for each associated parameter vector (not 
 * each element). The phase for any parameter vector that is originally set to be 
 * estimated ONLY in a time block later than maxYr is set to negative so it will 
 * not be estimated.
 * 
 * @param pVVI  - pointer to VectorVectorInfo for vector of parameters of interest
 * @param pgi   - pointer to ParameterGroupInfo for the parameter group of interest
 * @param pid   - index to parameter vector in ParameterGroupInfo
 * @param maxYr - max year for estimated parameter vectors
 * @param debug - flag to print debugging info
 */
void tcsam::adjustParamPhase(VectorVectorInfo* pVVI, ParameterGroupInfo* pgi, int pid, int maxYr, int debug){
    if (debug) rpt::echo<<"Starting adjustParamPhase(VectorVectorInfo*...) for "<<pVVI->name<<" for max year "<<maxYr<<" using pid "<<pid<<endl;
    int np = pVVI->getSize();//number of parameter vectors
    ivector done(1,np); done=0;
    ivector phases(1,np);
    for (int p=1;p<=np;p++){ //loop over parameter vectors
        int phase = (*pVVI)[p]->getPhase();
        if (debug) rpt::echo<<"--checking parameter vector "<<p<<" with phase "<<phase<<endl;
        if ((done[p]==0)&&(phase>0)){ //check only parameter vectors that haven't been determined and are set to be estimated
            int nPCs = pgi->nPCs;//number of parameter combinations in the associated parameter group
            bool estParam = false; 
            for (int pc=1;pc<=nPCs;pc++){ //loop over parameter combinations
                ivector pcids = pgi->getPCIDs(pc);
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter vector "<<p<<". and pcid = "<<pcids<<endl;
                int pcid = pcids[pid];//get id (index) of parameter vector used in this pc
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter vector "<<p<<" and got pcid "<<pcid<<endl;
                if (pcid==p){//make check for this parameter vector, since it's included in this pc
                    ivector years = pgi->getYearsForPC(pc);
                    if (debug) rpt::echo<<"----pc = "<<pc<<tb<<"years = "<<years<<endl;
                    for (int iy=years.indexmin();iy<=years.indexmax();iy++){
                        if (years[iy]<=maxYr){
                            //parameter is in pc that pertains to model year interval
                            estParam=true;//so parameter should be estimated
                            done[p]=1;//don't need to check further pc's
                            if (debug) rpt::echo<<"----estParam=TRUE and done=1"<<endl;
                            break;
                        }
                    }//--iy loop
                    if (estParam) break;
                }//--pcid==p
            }//--loop over pc's
            if (debug && !estParam) rpt::echo<<"----estParam=FALSE"<<endl;
            phase = estParam ? phase : -phase;//set phase to negative if parameter should NOT be estimated
        }//--(!done[p])&&(phase>0)
        (*pVVI)[p]->setPhase(phase);
        phases[p] = phase;
    }//--loop over parameter vectors
    if (debug) rpt::echo<<"--final phases: "<<phases<<endl;
    if (debug) rpt::echo<<"Finished adjustParamPhase(VectorVectorInfo*...) for "<<pVVI->name<<" for max year "<<maxYr<<endl<<endl;
}

/**
 * Adjusts devs vector phases to be consistent with maxYr.
 * 
 * For each devs vector represented by the DevsVectorVectorInfo object,
 * the phase for each element associated with a year later than maxYr 
 * is set to negative so it will not be estimated.
 * 
 * @param pDVVI - pointer to DevsVectorVectorInfo for the devs vectors of interest
 * @param maxYr - max year for estimated parameter vectors
 * @param debug - flag to print debugging info
 */
void tcsam::adjustParamPhase(DevsVectorVectorInfo* pDVVI, int maxYr, int debug){
    if (debug) rpt::echo<<"Starting adjustParamPhase(DevsVectorVectorInfo*...) for "<<pDVVI->name<<" for max year "<<maxYr<<endl;
    int np = pDVVI->getSize();//number of devs vectors
    for (int p=1;p<=np;p++){ //loop over devs vectors
        DevsVectorInfo* pDVI = (*pDVVI)[p];
        ivector phases = pDVI->getPhases();
        ivector years  = pDVI->getFwdIndices();
        if (debug){
            rpt::echo<<"--p = "<<p<<endl;
            rpt::echo<<"----original phases = "<<phases<<endl;
            rpt::echo<<"----years           ="<<years<<endl;
        }
        for (int i=years.indexmin();i<=years.indexmax();i++){
            if (years[i]>maxYr) phases[i] = -abs(phases[i]);//set estimation phase to negative
        }
        pDVI->setPhases(phases);
        if (debug) rpt::echo<<"----final phases: "<<pDVI->getPhases()<<endl;
    }//--loop over devs vectors
    if (debug) {
        rpt::echo<<"--final phases for all parameters: "<<pDVVI->getParameterPhases()<<endl;
        rpt::echo<<"Finished adjustParamPhase(DevsVectorVectorInfo*...) for "<<pDVVI->name<<" for max year "<<maxYr<<endl<<endl;
    }
}
/**
 * Adjusts parameter phases to be consistent with maxYr.
 * 
 * The phase for any parameter that is originally set to be estimated ONLY in a 
 * time block later than maxYr is set to negative so it will not be estimated.
 * 
 * @param pI - pointer to NumberVectorInfo for vector of parameters of interest
 * @param pgi - pointer to ParameterGroupInfo for the parameter group of interest
 * @param pid - index to parameter in ParameterGroupInfo
 * @param maxYr - max year for estimated parameters
 * @param sfs - ivector indicating fishery or survey adjustment
 * @param debug - flag to print debugging info
 */
void tcsam::adjustParamPhase(NumberVectorInfo* pI, ParameterGroupInfo* pgi, int pid, int maxYr, const ivector& sfs, int debug){
    if (debug) rpt::echo<<"Starting adjustParamPhase(NumberVectorInfo*...) for "<<pI->name<<" for max year "<<maxYr<<" using pid "<<pid<<endl;
    int np = pI->getSize();//number of parameters
    ivector phases = pI->getPhases();
    if (debug) rpt::echo<<"--initial phases: "<<phases<<endl;
    ivector done(1,np); done=0;
    for (int p=1;p<=np;p++){ //loop over parameters
        if (debug) rpt::echo<<"--checking parameter "<<p<<endl;
        if ((done[p]==0)&&(phases[p]>0)){ //check only parameters that haven't been determined and are set to be estimated
            int nPCs = pgi->nPCs;//number of parameter combinations in the associated parameter group
            bool estParam = false; 
            for (int pc=1;pc<=nPCs;pc++){ //loop over parameter combinations
                ivector pcids = pgi->getPCIDs(pc);
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter "<<p<<". and pcid = "<<pcids<<endl;
                int pcid = pcids[pid];//get id (index) of parameter used in this pc
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter "<<p<<" and got pcid "<<pcid<<endl;
                if (pcid==p){//make check for this parameter, since it's included in this pc
                    ivector years = pgi->getYearsForPC(pc);
                    if (debug) rpt::echo<<"----pc = "<<pc<<tb<<"years = "<<years<<endl;
                    for (int iy=years.indexmin();iy<=years.indexmax();iy++){
                        if (years[iy]<=maxYr+sfs[pc]){//adjust for survey or fishery
                            //parameter is in pc that pertains to model year interval
                            estParam=true;//so parameter should be estimated
                            done[p]=1;//don't need to check further pc's
                            if (debug) rpt::echo<<"----estParam=TRUE and done=1"<<endl;
                            break;
                        }
                    }//--iy loop
                    if (estParam) break;
                }//--pcid==p
            }//--loop over pc's
            if (debug && !estParam) rpt::echo<<"----estParam=FALSE"<<endl;
            phases[p] = estParam ? phases[p] : -phases[p];//set phase to negative if parameter should NOT be estimated
        }//--(!done[p])&&(phases[p]>0)
    }//--loop over parameters
    if (debug) rpt::echo<<"--final phases: "<<phases<<endl;
    pI->setPhases(phases);
    if (debug) rpt::echo<<"Finished adjustParamPhase(NumberVectorInfo*...) for "<<pI->name<<" for max year "<<maxYr<<endl<<endl;
}

/**
 * Adjusts parameter vector phases to be consistent with maxYr.
 * 
 * Estimation phases are defined for each associated parameter vector (not 
 * each element). The phase for any parameter vector that is originally set to be 
 * estimated ONLY in a time block later than maxYr is set to negative so it will 
 * not be estimated.
 * 
 * @param pVVI  - pointer to VectorVectorInfo for vector of parameters of interest
 * @param pgi   - pointer to ParameterGroupInfo for the parameter group of interest
 * @param pid   - index to parameter vector in ParameterGroupInfo
 * @param maxYr - max year for estimated parameter vectors
 * @param sfs - ivector indicating fishery or survey adjustment
 * @param debug - flag to print debugging info
 */
void tcsam::adjustParamPhase(VectorVectorInfo* pVVI, ParameterGroupInfo* pgi, int pid, int maxYr, const ivector& sfs, int debug){
    if (debug) rpt::echo<<"Starting adjustParamPhase(VectorVectorInfo*...) for "<<pVVI->name<<" for max year "<<maxYr<<" using pid "<<pid<<endl;
    int np = pVVI->getSize();//number of parameter vectors
    ivector done(1,np); done=0;
    ivector phases(1,np);
    for (int p=1;p<=np;p++){ //loop over parameter vectors
        int phase = (*pVVI)[p]->getPhase();
        if (debug) rpt::echo<<"--checking parameter vector "<<p<<" with phase "<<phase<<endl;
        if ((done[p]==0)&&(phase>0)){ //check only parameter vectors that haven't been determined and are set to be estimated
            int nPCs = pgi->nPCs;//number of parameter combinations in the associated parameter group
            bool estParam = false; 
            for (int pc=1;pc<=nPCs;pc++){ //loop over parameter combinations
                ivector pcids = pgi->getPCIDs(pc);
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter vector "<<p<<". and pcid = "<<pcids<<endl;
                int pcid = pcids[pid];//get id (index) of parameter vector used in this pc
                if (debug) rpt::echo<<"checking pc "<<pc<<" for parameter vector "<<p<<" and got pcid "<<pcid<<endl;
                if (pcid==p){//make check for this parameter vector, since it's included in this pc
                    ivector years = pgi->getYearsForPC(pc);
                    if (debug) rpt::echo<<"----pc = "<<pc<<tb<<"years = "<<years<<endl;
                    for (int iy=years.indexmin();iy<=years.indexmax();iy++){
                        if (years[iy]<=maxYr+sfs[pc]){//adjust for survey or fishery
                            //parameter is in pc that pertains to model year interval
                            estParam=true;//so parameter should be estimated
                            done[p]=1;//don't need to check further pc's
                            if (debug) rpt::echo<<"----estParam=TRUE and done=1"<<endl;
                            break;
                        }
                    }//--iy loop
                    if (estParam) break;
                }//--pcid==p
            }//--loop over pc's
            if (debug && !estParam) rpt::echo<<"----estParam=FALSE"<<endl;
            phase = estParam ? phase : -phase;//set phase to negative if parameter should NOT be estimated
        }//--(!done[p])&&(phase>0)
        (*pVVI)[p]->setPhase(phase);
        phases[p] = phase;
    }//--loop over parameter vectors
    if (debug) rpt::echo<<"--final phases: "<<phases<<endl;
    if (debug) rpt::echo<<"Finished adjustParamPhase(VectorVectorInfo*...) for "<<pVVI->name<<" for max year "<<maxYr<<endl<<endl;
}

/**
 * Adjusts devs vector phases to be consistent with maxYr.
 * 
 * For each devs vector represented by the DevsVectorVectorInfo object,
 * the phase for each element associated with a year later than maxYr 
 * is set to negative so it will not be estimated.
 * 
 * @param pDVVI - pointer to DevsVectorVectorInfo for the devs vectors of interest
 * @param maxYr - max year for estimated parameter vectors
 * @param sfs - ivector indicating fishery or survey adjustment
 * @param debug - flag to print debugging info
 */
void tcsam::adjustParamPhase(DevsVectorVectorInfo* pDVVI, ParameterGroupInfo* pgi, int pid, int maxYr, const ivector& sfs, int debug){
    if (debug) rpt::echo<<"Starting adjustParamPhase(DevsVectorVectorInfo*...) for "<<pDVVI->name<<" for max year "<<maxYr<<endl;
    int np = pDVVI->getSize();//number of devs vectors
    for (int p=1;p<=np;p++){ //loop over devs vectors
        DevsVectorInfo* pDVI = (*pDVVI)[p];
        ivector phases = pDVI->getPhases();
        ivector years  = pDVI->getFwdIndices();
        if (debug){
            rpt::echo<<"--p = "<<p<<endl;
            rpt::echo<<"----original phases = "<<phases<<endl;
            rpt::echo<<"----years           ="<<years<<endl;
        }
        int adj = 0;
        for (int pc=1;pc<=pgi->nPCs;pc++){
            if (pgi->getPCIDs(pc)[pid]==p){
                if (debug) rpt::echo<<"devs vector "<<p<<" is in pc "<<pc<<endl;
                if (sfs[pc]>0) {
                    adj=1;//adjustment for surveys
                    if (debug) rpt::echo<<"--adjusting max year for survey"<<endl;
                }
            }
        }
        for (int i=years.indexmin();i<=years.indexmax();i++){
            if (years[i]>maxYr+adj) phases[i] = -abs(phases[i]);//set estimation phase to negative
        }
        pDVI->setPhases(phases);
        if (debug) rpt::echo<<"----final phases: "<<pDVI->getPhases()<<endl;
    }//--loop over devs vectors
    if (debug) {
        rpt::echo<<"--final phases for all parameters: "<<pDVVI->getParameterPhases()<<endl;
        rpt::echo<<"Finished adjustParamPhase(DevsVectorVectorInfo*...) for "<<pDVVI->name<<" for max year "<<maxYr<<endl<<endl;
    }
}

/*----------------------------------------------------------------------------\n
 * ParameterGroupInfo\n
 -----------------------------------------------------------------------------*/
/**
 * Constructor.
 */
ParameterGroupInfo::ParameterGroupInfo(){
    nIVs=0; nPVs=0; nXIs=0; nIBSs=0; nPCs=0;
    ppIBSs=0; ppIdxs=0;
}
/**
 * Destructor.
 */
ParameterGroupInfo::~ParameterGroupInfo(){
    if (ppIBSs) {
        for (int i=0;i<nIBSs;i++) delete ppIBSs[i]; 
        delete ppIBSs;
        nIBSs=0;
    }
    if (ppIdxs) {
        for (int p=0;p<nPCs;p++) delete ppIdxs[p]; 
        delete ppIdxs;
        ppIdxs=0;
    }
}

/**
 * Tests whether a given year is associated with a parameter combination.
 * 
 * @param y - the year in question
 * @param pc - the pc in question
 * 
 * @return - int 1/0 if y is/not associated with the parameter combination
 */
int ParameterGroupInfo::isYearInPC(int y, int pc){
    if (debug) cout<<"starting ParameterGroupInfo::isYearInPC("<<y<<cc<<pc<<") for "<<name<<endl;
    IndexBlockSet* pIBS = getIndexBlockSet("YEAR_BLOCK");
    ivector iv = pIBS->getFwdIndexVector(pc);
    if (debug) cout<<pc<<": "<<iv<<endl;
    for (int i=iv.indexmin();i<=iv.indexmax();i++) 
        if (iv(i)==y) {
            if (debug) cout<<"year is in pc "<<endl;
            return(1);
        }
    if (debug) cout<<"year is NOT in pc "<<endl;
    return(0);
}

/**
 * Tests whether a given pc year is associated a year.
 * 
 * @param pc - the pc in question
 * @param y - the year in question
 * 
 * @return - int 1/0 if the parameter combination is/not associated with the year
 */
int ParameterGroupInfo::isPCinYearBlock(int pc, const ivector& yb){
    if (debug) cout<<"starting ParameterGroupInfo::isPCInYearBlock("<<pc<<cc<<endl<<yb<<") for "<<name<<endl;
    IndexBlockSet* pIBS = getIndexBlockSet("YEAR_BLOCK");
    ivector iv = pIBS->getFwdIndexVector(pc);
    if (debug) cout<<pc<<": "<<iv<<endl;
    for (int i=iv.indexmin();i<=iv.indexmax();i++) 
        for (int y=yb.indexmin();y<=yb.indexmax();y++) 
            if (iv(i)==y) {
                if (debug) cout<<"pc is in year block "<<endl;
                return(1);
            }
    if (debug) cout<<"pc is NOT in year block "<<endl;
    return(0);
}

/**
 * Finds the number of parameter combinations that correspond to a given year.
 * 
 * @param y - the year to find the PCs for
 * @return - number of corresponding PCs
 */
int ParameterGroupInfo::getNumPCsForYear(int y){
    if (debug) cout<<"starting ParameterGroupInfo::getNumPCsForYear("<<y<<") for "<<name<<endl;
    IndexBlockSet* pIBS = getIndexBlockSet("YEAR_BLOCK");
    if (debug) cout<<"nPCs = "<<nPCs<<endl;
    //count pc's that include year
    int k = 0;
    for (int pc=1;pc<=nPCs;pc++){
        ivector iv = pIBS->getFwdIndexVector(pc);
        if (debug) cout<<pc<<": "<<iv<<endl;
        for (int i=iv.indexmin();i<=iv.indexmax();i++) 
            if (iv(i)==y) k++;
    }
    if (debug) cout<<"Found "<<k<<" PCs that include "<<y<<endl;
    return(k);
}

/**
 * Set flags to include/excude (1/0) parameter combinations in likelihood
 * given a set of model years.
 * 
 * @param model_years - the model years in question
 * 
 * @return - nothing
 * 
 * @details Determines results of calls to isPCinLikelihood(pc).
 */
void ParameterGroupInfo::setLikelihoodFlagsForPCs(const ivector& model_years){
    for (int pc=1;pc<=nPCs;pc++){
        includePCinLikelihood(pc) = isPCinYearBlock(pc,model_years);
    }
}

/**
 * Finds the pc indices that corresponds to a given year.
 * 
 * @param y - the year to find the PCs for
 * @return - ivector of the corresponding PCs (or unallocated, if no corresponding PC found)
 */
ivector ParameterGroupInfo::getPCsForYear(int y){
    if (debug) cout<<"starting ParameterGroupInfo::getPCforYear("<<y<<") for "<<name<<endl;
    ivector pcs; //output vector of pc indices
    int k = getNumPCsForYear(y);//count pc's that include year
    if (debug) cout<<"nPCs = "<<nPCs<<endl;
    if (debug) cout<<"Found "<<k<<" PCs that include "<<y<<endl;
    if (k) {
        IndexBlockSet* pIBS = getIndexBlockSet("YEAR_BLOCK");
        //now collect the pc indices
        pcs.allocate(1,k);
        k = 1;//reset counter
        for (int pc=1;pc<=nPCs;pc++){
            ivector iv = pIBS->getFwdIndexVector(pc);
            if (debug) cout<<pc<<": "<<iv<<endl;
            for (int i=iv.indexmin();i<=iv.indexmax();i++) 
                if (iv(i)==y) pcs(k++) = pc;
        }
    }
    if (debug) {
        if (pcs.allocated()) 
            cout<<"pcs corresponding to "<<y<<" are "<<pcs<<" for "<<name<<endl; 
        else
            cout<<"No pcs corresponding to "<<y<<" were found for "<<name<<endl; 
        cout<<"finished ParameterGroupInfo::getPCforYear("<<y<<") for "<<name<<endl;
    }
    return(pcs);
}

/**
 * Finds the years that a parameter combination (PC) is defined for
 * 
 * @param pc - the PC of interest
 * 
 * @return - ivector the corresponding years (or unallocated, if no corresponding PC found)
 */
ivector ParameterGroupInfo::getYearsForPC(int pc){
    if (debug) cout<<"starting ParameterGroupInfo::getYearsForPC("<<pc<<") for "<<name<<endl;
    IndexBlockSet* pIBS = ParameterGroupInfo::getIndexBlockSet("YEAR_BLOCK");
    ivector iv = pIBS->getFwdIndexVector(pc);
    if (debug) {
        cout<<"years for pc: "<<iv<<endl;
        cout<<"Finished ParameterGroupInfo::getYearsForPC("<<pc<<") for "<<name<<endl;
    }
    return iv;
}

/**
 * Add a year to a parameter combination.
 * 
 * @param pc - parameter combination index to add year to
 * @param y - year to add
 */
void ParameterGroupInfo::addYearToPC(int pc, int y){
    if (debug) cout<<"starting ParameterGroupInfo::addYearToPC("<<pc<<cc<<y<<") for "<<name<<endl;
    IndexBlockSet* pIBS = getIndexBlockSet("YEAR_BLOCK");
    IndexBlock* pIB  = pIBS->getIndexBlock(pc);
    if (debug) IndexBlock::debug=1;
    pIB->addElement(y);
    if (debug) IndexBlock::debug=0;
    createIndices();
    if (debug) cout<<"finished ParameterGroupInfo::addYearToPC("<<pc<<cc<<y<<") for "<<name<<endl;
}

/**
 * Gets indices for parameter combination pc.
 * @param pc - id for desired parameter combination
 * 
 * @return - ivector with corresponding parameter indices
 */
ivector ParameterGroupInfo::getPCIDs(int pc){
    if (debug) cout<<"starting ParameterGroupInfo::getPCIDs(int pc)"<<endl;
    if (debug) {
        cout<<"pcids for "<<pc<<"th parameter combination: "<<in(pc)<<endl;
        cout<<"finished ParameterGroupInfo::getPCIDs(int pc)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    return in(pc);
}

/**
 * Gets "extra" values for parameter combination pc.
 * 
 * @param pc - id for desired parameter combination
 * 
 * @return - dvector of corresponding values
 */
dvector ParameterGroupInfo::getPCXDs(int pc){
    if (debug) cout<<"starting ParameterGroupInfo::getPCXDs(int pc)"<<endl;
    if (debug) {
        cout<<"pcids for "<<pc<<"th parameter combination: "<<xd(pc)<<endl;
        cout<<"finished ParameterGroupInfo::getPCXDs(int pc)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    return xd(pc);
}

/**
 * Gets a matrix of model indices corresponding 
 * to the given parameter combination pc.
 * 
 * @param pc
 * 
 * @return - imatrix of model indices
 */
imatrix ParameterGroupInfo::getModelIndices(int pc){
    if (debug) cout<<"starting ParameterGroupInfo::getModelIndices(int pc)"<<endl;
    if (pc>nPCs) {
        cout<<"Error in ParameterGroupInfo::getModelIndices(int pc) for "<<name<<endl;
        cout<<"pc was "<<pc<<" but max defined is "<<nPCs<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    imatrix idxs = (*ppIdxs[pc-1]);
    if (debug) {
        cout<<"imatrix idxs has rows "<<idxs.indexmin()<<":"<<idxs.indexmax()<<endl;
        cout<<"finished ParameterGroupInfo::getModelIndices(int pc)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    return idxs;
}

/**
 * Creates the pointer to the IndexBlockSets object 
 * and the individual IndexBlockSet objects corresponding
 * to the blocks defined for model indices for this 
 * ParameterGroup object.
 */
void ParameterGroupInfo::createIndexBlockSets(){
    if (debug) cout<<"starting ParameterGroupInfo::createIndexBlockSets() "<<nIBSs<<endl;
    if(ppIBSs) delete(ppIBSs);
    if (nIBSs){
        ppIBSs = new IndexBlockSet*[nIBSs];
        for (int i=1;i<=nIBSs;i++){
            ppIBSs[i-1] = new IndexBlockSet();
        }
    }
    if (debug) cout<<"finished ParameterGroupInfo::createIndexBlockSets()"<<endl;
}

/* 
 * Returns a pointer to the index block set identified by "type".
 * Inputs:
 *  adstring type:  "type" identifying index block set to return
 * Returns:
 *  pointer to the identified IndexBlockSet
 */
IndexBlockSet* ParameterGroupInfo::getIndexBlockSet(adstring type){
    if (debug) cout<<"starting  IndexBlockSets::getIndexBlockSet("<<type<<")"<<endl;
    IndexBlockSet* p = 0;
    int s=0;
    while(s<nIBSs){
        p = ppIBSs[s];
        if (p->getType()==type) break;
        s++;
    }
    if (debug) cout<<"finished  IndexBlockSets::getIndexBlockSet("<<type<<")"<<endl;
    return p;
}

/**
 * Creates the set of imatrix's representing the model indices 
 * corresponding to each parameter combination.
 */
void ParameterGroupInfo::createIndices(void){
    if (debug) cout<<"starting void ParameterGroupInfo::createIndices(void)"<<endl;
    if (ppIdxs) delete(ppIdxs);
    if (nPCs>0){
        //create pointer array
        ppIdxs = new imatrix*[nPCs];
        //loop over parameter combinations and create an indices imatrix for each
        for (int p=1;p<=nPCs;p++){
            int nc=1;//number of rows for indices imatrix
            imatrix tmp(1,nIVs);//matrix of values of indices
            for (int i=1;i<=nIVs;i++){//loop over index variables
                if (lblIVs(i)==tcsam::STR_SEX){
                    tmp(i).allocate(1,1);
                    tmp(i,1) = in(p,i);
                    nc *= 1;//just for clarity
                } else
                if (lblIVs(i)==tcsam::STR_MATURITY_STATE){
                    tmp(i).allocate(1,1);
                    tmp(i,1) = in(p,i);
                    nc *= 1;//just for clarity
                } else
                if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {
                    tmp(i).allocate(1,1);
                    tmp(i,1) = in(p,i);
                    nc *= 1;//just for clarity
                } else
                if (lblIVs(i)==tcsam::STR_FISHERY) {
                    tmp(i).allocate(1,1);
                    tmp(i,1) = in(p,i);
                    nc *= 1;//just for clarity
                } else
                if (lblIVs(i)==tcsam::STR_SURVEY) {
                    tmp(i).allocate(1,1);
                    tmp(i,1) = in(p,i);
                    nc *= 1;//just for clarity
                } else {
                    if (debug) cout<<"Parsing label '"<<lblIVs(i)<<"' for IndexBlockSet dim type."<<endl;
                    //parse label (should be of form 'type_BLOCK')
                    int n = lblIVs(i).pos("_");
                    if (n) { //can parse label correctly
                        adstring type = lblIVs(i)(1,n-1);
                        int id = in(p,i);//block index in IndexBlockSet
                        IndexBlockSet* pIBS = 0;
                        if (debug) cout<<"Checking local IndexBlockSets using type '"<<type<<"'."<<endl;
    //                    pIBS = ptrIBSs->getIndexBlockSet(type);
                        pIBS = ppIBSs[in(p,i)-1];
                        if (pIBS){
                            ivector idxs = pIBS->getFwdIndexVector(p);
                            tmp(i).allocate(idxs.indexmin(),idxs.indexmax());
                            tmp(i) = idxs;
                            nc *= idxs.size();
                            if (debug) cout<<"indices for block = "<<tmp(i)<<endl;
                        } else {
                            cout<<"Error in ParameterGroupInfo::createIndices(void) for "<<name<<endl;
                            cout<<"Could not find IndexBlockSet for dim type '"<<type<<endl;
                            cout<<"Dim label was '"<<lblIVs(i)<<"'"<<endl;
                            cout<<"Aborting..."<<endl;
                            exit(-1);
                        }
                    } else {
                        cout<<"Error in ParameterGroupInfo::createIndices(void) for "<<name<<endl;
                        cout<<"Could not parse index variable label '"<<lblIVs(i)<<"' correctly"<<endl;
                        cout<<"Should have form 'type_BLOCK', where 'type' is dimension type (e.g., 'YEAR')."<<endl;
                        cout<<"Aborting..."<<endl;
                        exit(-1);
                    }
                }            
            }
            ppIdxs[p-1] = new imatrix();
            (*ppIdxs[p-1]).allocate(1,nc,1,nIVs);
            imatrix m(1,nIVs,1,nc);
            for (int i=1;i<=nIVs;i++){//loop over index variables
                ivector itmp = tmp(i);
                int sz = itmp.size();
                int cmn = 1;
                while ((cmn+sz-1)<=nc){
                    m(i)(cmn,cmn+sz-1) = itmp.shift(cmn);
                    cmn += sz;
                }
            }
            (*ppIdxs[p-1]) = trans(m);
            if (debug) {
                cout<<"Index matrix for pc "<<p<<endl;
                cout<<"#";
                for (int i=1;i<=nIVs;i++) cout<<lblIVs(i)<<tb; cout<<endl;
                cout<<(*ppIdxs[p-1])<<endl;  
            }
        }//p
    }//nPCs>0
    if (debug) {
        cout<<"finished void ParameterGroupInfo::createIndices(void)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Reads info for a BoundedNumberVector object from an input filestream.
 * @param is
 * @param lbl
 * @param pBNVI
 * @return 
 */
BoundedNumberVectorInfo* ParameterGroupInfo::read(cifstream& is, adstring& lbl, BoundedNumberVectorInfo* pBNVI){
    adstring param;
    is>>param;
    if (param==lbl){
        if (debug) cout<<"Reading "<<lbl<<endl;
        pBNVI = new BoundedNumberVectorInfo(lbl);
        is>>(*pBNVI);
        if (debug) cout<<"ptr to "<<lbl<<": "<<pBNVI<<endl<<(*pBNVI)<<endl;
    } else {
        tcsam::readError(is,lbl,param);
        exit(-1);
    }
    return pBNVI;
}

/**
 * 
 * @param is
 * @param lbl
 * @param pBVVI
 * @return 
 */
BoundedVectorVectorInfo* ParameterGroupInfo::read(cifstream& is, adstring& lbl, BoundedVectorVectorInfo* pBVVI){
    if (debug) cout<<"Starting BVVI* PGI::read(is,lbl,BVVI*) for "<<lbl<<endl;
    adstring param;
    is>>param;
    if (param==lbl){
        if (debug) cout<<"Reading "<<lbl<<endl;
        pBVVI = new BoundedVectorVectorInfo(lbl);
        is>>(*pBVVI);
        if (debug) cout<<"ptr to "<<lbl<<": "<<pBVVI<<endl<<(*pBVVI)<<endl;
    } else {
        tcsam::readError(is,lbl,param);
        exit(-1);
    }
    if (debug) cout<<"Finished BVVI* PGI::read(is,lbl,BVVI*) for "<<lbl<<endl;
    return pBVVI;
}

/**
 * 
 * @param is
 * @param lbl
 * @param pDVVI
 * @return 
 */
DevsVectorVectorInfo* ParameterGroupInfo::read(cifstream& is, adstring& lbl, DevsVectorVectorInfo* pDVVI){
    adstring param;
    is>>param;
    if (param==lbl){
        if (debug) cout<<"Reading "<<lbl<<endl;
        pDVVI = new DevsVectorVectorInfo(lbl);
        is>>(*pDVVI);
        if (debug) cout<<"ptr to "<<lbl<<": "<<pDVVI<<endl<<(*pDVVI)<<endl;
    } else {
        tcsam::readError(is,lbl,param);
        exit(-1);
    }
    return pDVVI;
}

/**
 * Reads info for the parameter group from an input filestream.
 * @param is
 */
void ParameterGroupInfo::read(cifstream& is){
    if (debug) cout<<"starting void ParameterGroupInfo::read(is)"<<endl;
    adstring str; 
    //read parameter combinations
    is>>str;
    if (debug) cout<<str<<tb<<"#Required keyword (PARAMETER_COMBINATIONS)"<<endl;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETER_COMBINATIONS)"<<endl;
    if (str=="PARAMETER_COMBINATIONS"){
        is>>nPCs;
        if (debug) cout<<nPCs<<tb<<"#number of parameter combinations"<<endl;
        rpt::echo<<nPCs<<tb<<"#number of parameter combinations"<<endl;
        if (nPCs>0){
            includePCinLikelihood.allocate(1,nPCs);
            includePCinLikelihood = 1;//default  to include all
            if (debug){
                cout<<"#id  ";
                for (int i=1;i<=nIVs;i++) cout<<lblIVs(i)<<tb; 
                for (int i=1;i<=nPVs;i++) cout<<lblPVs(i)<<tb; 
                for (int i=1;i<=nXIs;i++) cout<<lblXIs(i)<<tb; 
            }
            rpt::echo<<"#id  "; 
            for (int i=1;i<=nIVs;i++) rpt::echo<<lblIVs(i)<<tb; 
            for (int i=1;i<=nPVs;i++) rpt::echo<<lblPVs(i)<<tb; 
            for (int i=1;i<=nXIs;i++) rpt::echo<<lblXIs(i)<<tb; 
            if (debug) cout<<tb<<"label"<<tb<<endl;
            rpt::echo<<tb<<"label"<<tb<<endl;
            if (debug) cout<<"About to allocate ppIBSs"<<endl;
            for (int i=0;i<nIBSs;i++) ppIBSs[i]->allocate(nPCs);
            if (debug) cout<<"Finished allocating ppIBSs"<<endl;
            //read parameters combinations definition matrix
            int ibsIdx=1;
            in.allocate(1,nPCs,1,nIVs+nPVs+nXIs);
            if (nXIs) xd.allocate(1,nPCs,1,nXIs);
            pcLabels.allocate(1,nPCs);
            for (int r=1;r<=nPCs;r++){//loop over rows
                is>>str; //read id
                rpt::echo<<str<<tb;
                if (debug) cout<<"pc = "<<r<<tb<<"looping over index variables"<<endl;
                for (int i=1;i<=nIVs;i++){//loop over index variables
                    is>>str;
                    if (debug) cout<<tb<<"index variable"<<tb<<i<<tb<<str<<endl;
                    rpt::echo<<str<<tb;
                    if (lblIVs(i)==tcsam::STR_SEX)             {in(r,i) = tcsam::getSexType(str);}      else
                    if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {in(r,i) = tcsam::getMaturityType(str);} else
                    if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {in(r,i) = tcsam::getShellType(str);}    else
                    if ((ibsIdxs.allocated())&&(i==ibsIdxs(ibsIdx))){//variable is a block
                        ppIBSs[ibsIdx-1]->getIndexBlock(r)->parse(str);
                        in(r,i)=ibsIdx;//index to associated IndexBlockSet
                        if (ibsIdx<nIBSs) ibsIdx++;//increment to next IBS
                    } else {in(r,i)=::atoi(str);}
                }
                if (debug) cout<<"looping over parameter indices"<<endl;
                for (int p=1;p<=nPVs;p++) {is>>in(r,nIVs+p); rpt::echo<<in(r,nIVs+p)<<tb;}  
                if (debug) cout<<"looping over extra variables"<<endl;
                for (int x=1;x<=nXIs;x++) {//loop over "extra" variables
                    is>>str;
                    if (debug) cout<<str<<endl;
                    rpt::echo<<str<<tb;
                    if (lblXIs(x)==tcsam::STR_SELFCN){
                        //identify selectivity function and return function index
                        in(r,nIVs+nPVs+x) = SelFcns::getSelFcnID(str);
                    } else {
                        in(r,nIVs+nPVs+x)=::atoi(str);
                        xd(r,x)=::atof(str);
                    }
                    //rpt::echo<<"x = "<<x<<tb<<"str = "<<str<<tb<<"in() = "<<in(r,nIVs+nPVs+x)<<tb<<"xd() = "<<xd(r,x)<<endl;
                }
                is>>pcLabels(r);
                if (debug) cout<<pcLabels(r)<<endl;
                rpt::echo<<pcLabels(r)<<tb;
                rpt::echo<<endl;
                if (debug) {
                    cout<<"pc row "<<r<<": "<<in(r);
                    if (nXIs) cout<<tb<<"xd: "<<xd(r);
                    cout<<endl;
                }
            }//--r
            //revert reading back to sub-class to read values for parameters
        }//nPCs>0
    } else {
        cout<<"Reading ParameterGroupInfo for '"<<name<<"' from file "<<is.get_file_name()<<endl;
        cout<<"Expected 'PARAMETER_COMBINATIONS' keyword but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug) ParameterGroupInfo::write(cout);        
    
    if (nIVs) createIndices();
    
    //now back to sub-class to read values for parameters
    if (debug) cout<<"finished void ParameterGroupInfo::read(is)"<<endl;
}

//void ParameterGroupInfo::setToWriteVectorInitialValues(bool flag){
//    //do nothing [not sure why have to provide a concrete
//}

void ParameterGroupInfo::write(std::ostream& os){
    if (debug) cout<<"starting ParameterGroupInfo::write(std::ostream& os)"<<endl;
//    os<<(*ptrIBSs)<<endl;
    os<<"PARAMETER_COMBINATIONS"<<endl;
    os<<nPCs<<tb<<"#number of parameter combinations"<<endl;
    if (nPCs>0){
        os<<"#id  "; 
        for (int i=1;i<=nIVs;i++) os<<lblIVs(i)<<tb; 
        for (int i=1;i<=nPVs;i++) os<<lblPVs(i)<<tb; 
        for (int i=1;i<=nXIs;i++) os<<lblXIs(i)<<tb; 
        os<<"label"<<tb;
        os<<endl;
        for (int r=1;r<=nPCs;r++){//loop over rows
            os<<r<<tb;
            int ibsIdx = 1;
            for (int i=1;i<=nIVs;i++){//loop over index variables           
                if (lblIVs(i)==tcsam::STR_SEX)             {os<<tcsam::getSexType(in(r,i))<<tb;} else
                if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {os<<tcsam::getMaturityType(in(r,i))<<tb;} else
                if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {os<<tcsam::getShellType(in(r,i))<<tb;} else 
                if ((ibsIdxs.allocated())&&(i==ibsIdxs(ibsIdx))){
                    os<<(*ppIBSs[ibsIdx-1]->getIndexBlock(r))<<tb;
                    if (ibsIdx<nIBSs) ibsIdx++;//increment to next
                } else {os<<in(r,i)<<tb;}
            }
            for (int p=1;p<=nPVs;p++) os<<in(r,nIVs+p)<<tb;      //loop over parameter variables
            for (int x=1;x<=nXIs;x++) { //loop over extra indices
                if (lblXIs(x)==tcsam::STR_SELFCN) os<<SelFcns::getSelFcnID(in(r,nIVs+nPVs+x))<<tb; else
                os<<xd(r,x)<<tb;
        //        os<<in(r,nIVs+nPVs+x)<<tb;

            }
            os<<pcLabels(r)<<tb;
            os<<endl;
        }//--r
    }//nPCs>0
    if (debug) cout<<"finished ParameterGroupInfo::write(std::ostream& os)"<<endl;
}
/**
 * This should be called from sub-classes to write generic information to R.
 * @param os
 */
void ParameterGroupInfo::writeToR(std::ostream& os){
    adstring lbls = "";
//    if (nIVs) lbls += wts::to_qcsv(lblIVs);
//    if (nPVs) {if (lbls=="") lbls += wts::to_qcsv(lblPVs); else lbls += cc+wts::to_qcsv(lblPVs);}
    if (nPVs) lbls += wts::to_qcsv(lblPVs);
    if (nXIs) {if (lbls=="") lbls += wts::to_qcsv(lblXIs); else lbls += cc+wts::to_qcsv(lblXIs);}
    os<<"pgi=list(name="<<qt<<name<<qt;
    if (nPCs>0) {
        os<<cc<<endl;
        os<<"pcs=list("<<endl;
            for (int p=1;p<=nPCs;p++){
                os<<qt<<p<<qt<<"=list(";
                os<<"label='"<<pcLabels(p)<<"',"<<endl;
                int ibsIdx = 1;
                for (int i=1;i<=nIVs;i++){//loop over index variables
                    os<<lblIVs(i)<<"='";
                    if (lblIVs(i)==tcsam::STR_SEX)             {os<<tcsam::getSexType(in(p,i))     <<"',";} else
                    if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {os<<tcsam::getMaturityType(in(p,i))<<"',";} else
                    if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {os<<tcsam::getShellType(in(p,i))   <<"',";} else 
                    if ((ibsIdxs.allocated())&&(i==ibsIdxs(ibsIdx))){
                        os<<(*ppIBSs[ibsIdx-1]->getIndexBlock(p))<<"',";
                        if (ibsIdx<nIBSs) ibsIdx++;//increment to next
                    } else {os<<in(p,i)<<"',";}
                }
                ivector iv = getPCIDs(p);
                imatrix im = getModelIndices(p);
                os<<"ids.PC="; wts::writeToR(os,iv(nIVs+1,nIVs+nPVs+nXIs),lbls); os<<cc<<endl;
                adstring ids = "index=c(1:"+str(im.indexmax())+")";
                adstring tps = "type=c("+wts::to_qcsv(lblIVs)+")";
                os<<"ids.Mod="; wts::writeToR(os,im,ids,tps); os<<"),"<<endl;
            }
        os<<"NULL)"<<endl;
    }//nPCs>0
    os<<")";
}
/*------------------------------------------------------------------------------
 * RecruitmentInfo
 -----------------------------------------------------------------------------*/
adstring RecruitmentInfo::NAME = "recruitment";
RecruitmentInfo::RecruitmentInfo(){
    if (debug) cout<<"starting RecruitmentInfo::RecruitmentInfo()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 1;
    lblIVs.allocate(1,nIVs);
    lblIVs(k=1) = "YEAR_BLOCK";
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs = 6;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pLnR";    dscPVs(k++) = "ln-scale mean recruitment";
    lblPVs(k) = "pRCV";    dscPVs(k++) = "recruitment cv's";
    lblPVs(k) = "pRX";     dscPVs(k++) = "fraction males at recruitment";
    lblPVs(k) = "pRa";     dscPVs(k++) = "size-at-recruitment parameter a";
    lblPVs(k) = "pRb";     dscPVs(k++) = "size-at-recruitment parameter b";    
    lblPVs(k) = "pDevsLnR";dscPVs(k++) = "ln-scale recruitment devs";    
    k=1;
    pLnR = new BoundedNumberVectorInfo(lblPVs(k++));
    pRCV = new BoundedNumberVectorInfo(lblPVs(k++));
    pRX  = new BoundedNumberVectorInfo(lblPVs(k++));
    pRa  = new BoundedNumberVectorInfo(lblPVs(k++));
    pRb  = new BoundedNumberVectorInfo(lblPVs(k++));
    pDevsLnR = new DevsVectorVectorInfo(lblPVs(k++));
    
    nXIs=1;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "nllWgt";
}

RecruitmentInfo::~RecruitmentInfo(){
    if (pLnR) delete pLnR;  pLnR = 0;
    if (pRCV) delete pRCV;  pRCV = 0;
    if (pRX)  delete pRX;   pRX  = 0;
    if (pRa)  delete pRa;   pRa  = 0;
    if (pRb)  delete pRb;   pRb  = 0;
    if (pDevsLnR) delete pDevsLnR; pDevsLnR = 0;
}

/**
 * Set the maximum year for estimating parameters.
 * 
 * @param mxYr - the max year to allow estimated parameters
 */
void RecruitmentInfo::setMaxYear(int mxYr){
    int k = nIVs+1;
    if (pLnR)     tcsam::adjustParamPhase(pLnR,this,k++,mxYr,1);
    if (pRCV)     tcsam::adjustParamPhase(pRCV,this,k++,mxYr,1);
    if (pRX)      tcsam::adjustParamPhase(pRX, this,k++,mxYr,1);
    if (pRa)      tcsam::adjustParamPhase(pRa, this,k++,mxYr,1);
    if (pRb)      tcsam::adjustParamPhase(pRb, this,k++,mxYr,1);
    if (pDevsLnR) tcsam::adjustParamPhase(pDevsLnR,mxYr,1);
}

/**
 * update RecruitmentInfo for a 1-year projected scenario.
 * 
 * @param closed - flag indicating whether directed fishery is closed
 */
void RecruitmentInfo::addNextYearToInfo(int closed){
    if (debug) cout<<"starting void RecruitmentInfo::project("<<closed<<")"<<endl;
    int idxDevsLnR = nIVs+6;
    //find pcs corresponding to mxYr
    int mxYr = ModelConfiguration::mxYr;
    ivector pcs = ParameterGroupInfo::getPCsForYear(mxYr);
    if (debug) {
        if (pcs.size()) cout<<"PCs "<<pcs<<" apply to "<<mxYr<<endl;
        else cout<<"No PCs apply to <<"<<mxYr<<endl;
    }
    //define ivector to track which vectors in pDevsLnR have been modified
    ivector idsDevsLnR(pcs.indexmin(),pcs.indexmax());
    idsDevsLnR = 0;
    for (int i=pcs.indexmin();i<=pcs.indexmax();i++){
        ParameterGroupInfo::addYearToPC(pcs[i],mxYr+1);

        ivector ids = ParameterGroupInfo::getPCIDs(pcs[i]);
        if (debug) cout<<"PC IDs for pc "<<pcs[i]<<" = "<<ids<<endl;
        if (pDevsLnR->getSize()){
            if (debug) cout<<"extending pDevsLnR["<<ids[idxDevsLnR]<<"]"<<endl<<(*pDevsLnR)[ids[idxDevsLnR]]->getInitVals()<<endl;
            //Need to make sure not to addValue multiple times to same devs vector.
            int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
            for (int j=pcs.indexmin();j<=pcs.indexmax();j++) sumdsp += (idsDevsLnR[j]==ids[idxDevsLnR]);
            if (!sumdsp){
                if (debug) cout<<"Adding value to pDevsLnR["<<ids[idxDevsLnR]<<"]"<<endl;
                idsDevsLnR[i] = ids[idxDevsLnR];
                DevsVectorInfo* pDVI = (*pDevsLnR)[ids[idxDevsLnR]];
                if (debug) DevsVectorInfo::debug=1;
                pDVI->addValueOnParameterScale(0.0,mxYr+1);
                if (debug) DevsVectorInfo::debug=0;
            } else {
                if (debug) cout<<"Already added value to pDevsLnR["<<ids[idxDevsLnR]<<"]"<<endl;
            }
        }
    }
    if (debug) cout<<"finished void RecruitmentInfo::project("<<closed<<")"<<endl; 
}

void RecruitmentInfo::read(cifstream & is){
    if (debug) cout<<"starting void RecruitmentInfo::read(cifstream & is)"<<endl;    
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in RecruitmentInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    if (nPCs>0){
        is>>str;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            pLnR   = ParameterGroupInfo::read(is,lblPVs(k),pLnR);   
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnR)<<endl;   k++;
            pRCV = ParameterGroupInfo::read(is,lblPVs(k),pRCV); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRCV)<<endl; k++;
            pRX = ParameterGroupInfo::read(is,lblPVs(k),pRX); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRX)<<endl; k++;
            pRa  = ParameterGroupInfo::read(is,lblPVs(k),pRa);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRa)<<endl;  k++;
            pRb  = ParameterGroupInfo::read(is,lblPVs(k),pRb);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRb)<<endl;  k++;
            pDevsLnR = ParameterGroupInfo::read(is,lblPVs(k),pDevsLnR); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsLnR)<<endl; 
        } else {
            cout<<"Error reading RecruitmentInfo from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }//nPCs>0
    
    if (debug){
        cout<<"RecruitmentInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished RecruitmentInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write estimation phases for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write estimation phases to file
 */
void RecruitmentInfo::setToWriteVectorEstimationPhases(bool flag){
    if (pDevsLnR){
        for (int i=1;i<=pDevsLnR->getSize();i++){
            DevsVectorInfo* vi = (*pDevsLnR)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void RecruitmentInfo::setToWriteVectorInitialValues(bool flag){
    if (pDevsLnR){
        for (int i=1;i<=pDevsLnR->getSize();i++){
            DevsVectorInfo* vi = (*pDevsLnR)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
}

/**
 * Write RecruitmentInfo to output stream in ADMB format
 * 
 * @param os - output stream
 */
void RecruitmentInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    if (nPCs>0){
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pLnR)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pRCV)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pRX)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pRa)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pRb)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsLnR)<<endl;
    }
}

void RecruitmentInfo::writeToPin(std::ostream & os){
    os<<"#--recruitment parameters"<<endl;
    pLnR->writeToPin(os);
    pRCV->writeToPin(os);
    pRX->writeToPin(os);
    pRa->writeToPin(os);
    pRb->writeToPin(os);
    pDevsLnR->writeToPin(os);
}

/**
 * Write RecruitmentInfo to output stream in R format
 * 
 * @param os - output stream
 */
void RecruitmentInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"rec=list("<<endl;
        ParameterGroupInfo::writeToR(os);   
        if (nPCs) {
            os<<cc<<endl;
            pLnR->writeToR(os,"pLnR",indent+1); os<<cc<<endl;
            pRCV->writeToR(os,"pRCV",indent+1); os<<cc<<endl;
            pRX->writeToR(os, "pRX", indent+1); os<<cc<<endl;
            pRa->writeToR(os, "pRa", indent+1); os<<cc<<endl;
            pRb->writeToR(os, "pRb", indent+1); os<<cc<<endl;
            pDevsLnR->writeToR(os,"pDevsLnR",indent+1); os<<endl;
        }
    os<<")";
}

/*------------------------------------------------------------------------------
 * NaturalMortalityInfo
 -----------------------------------------------------------------------------*/
adstring NaturalMortalityInfo::NAME = "natural_mortality";
NaturalMortalityInfo::NaturalMortalityInfo(){
    if (debug) cout<<"starting NaturalMortalityInfo::NaturalMortalityInfo()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    if (debug) cout<<1<<endl;
    
    int k;
    //define index variables for parameters
    nIVs = 4;
    lblIVs.allocate(1,nIVs);
    k = 1;
    lblIVs(k++) = "YEAR_BLOCK";
    lblIVs(k++) = tcsam::STR_SEX;
    lblIVs(k++) = tcsam::STR_MATURITY_STATE;
    lblIVs(k++) = tcsam::STR_SHELL_CONDITION;
    if (debug) cout<<2<<endl;
    //define index block sets for index variables
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    if (debug) cout<<3<<endl;
    
    //define parameters
    nPVs=5;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pM";   dscPVs(k++) = "base natural mortality rate";
    lblPVs(k) = "pDM1"; dscPVs(k++) = "offset 1";
    lblPVs(k) = "pDM2"; dscPVs(k++) = "offset 2";
    lblPVs(k) = "pDM3"; dscPVs(k++) = "offset 3";
    lblPVs(k) = "pDM4"; dscPVs(k++) = "offset 4";
    if (debug) cout<<3<<endl;
    k=1;
    pM   = new BoundedNumberVectorInfo(lblPVs(k++));
    pDM1 = new BoundedNumberVectorInfo(lblPVs(k++));
    pDM2 = new BoundedNumberVectorInfo(lblPVs(k++));
    pDM3 = new BoundedNumberVectorInfo(lblPVs(k++));
    pDM4 = new BoundedNumberVectorInfo(lblPVs(k++));
    if (debug) cout<<4<<endl;
    
    //define "extra" indices
    nXIs=1;
    lblXIs.allocate(1,nXIs);
    lblXIs(k=1) = "zScaling";    
    if (debug) cout<<5<<endl;
    if (debug) cout<<"finished NaturalMortalityInfo::NaturalMortalityInfo()"<<endl;
}

NaturalMortalityInfo::~NaturalMortalityInfo(){
    if (pM)   delete pM;    pM =0;
    if (pDM1) delete pDM1;  pDM1 =0;
    if (pDM2) delete pDM2;  pDM2 =0;
    if (pDM3) delete pDM3;  pDM3 =0;
    if (pDM4) delete pDM4;  pDM4 =0;
}

/**
 * Set the maximum year for estimating parameters.
 * 
 * @param mxYr - the max year to allow estimated parameters
 */
void NaturalMortalityInfo::setMaxYear(int mxYr){
    int k = nIVs+1;
    if (pM)   tcsam::adjustParamPhase(pM,  this,k++,mxYr,1);
    if (pDM1) tcsam::adjustParamPhase(pDM1,this,k++,mxYr,1);
    if (pDM2) tcsam::adjustParamPhase(pDM2,this,k++,mxYr,1);
    if (pDM3) tcsam::adjustParamPhase(pDM3,this,k++,mxYr,1);
    if (pDM4) tcsam::adjustParamPhase(pDM4,this,k++,mxYr,1);
}

void NaturalMortalityInfo::addNextYearToInfo(int closed){
    if (debug) cout<<"starting void NaturalMortalityInfo::project("<<closed<<")"<<endl;
    //find pcs corresponding to mxYr
    int mxYr = ModelConfiguration::mxYr;
    ivector pcs = ParameterGroupInfo::getPCsForYear(mxYr);
    if (debug) {
        if (pcs.size()) cout<<"PCs "<<pcs<<" apply to "<<mxYr<<endl;
        else cout<<"No PCs apply to <<"<<mxYr<<endl;
    }
    for (int i=1;i<=pcs.size();i++){
        ParameterGroupInfo::addYearToPC(pcs[i],mxYr+1);
    }
    if (debug) cout<<"finished void NaturalMortalityInfo::project("<<closed<<")"<<endl;
}

void NaturalMortalityInfo::read(cifstream & is){
    if (debug) cout<<"starting NaturalMortalityInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in NaturalMortalityInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    is>>zRef;    
    rpt::echo<<"zRef = "<<zRef<<endl;
    if (debug) ParameterGroupInfo::debug=1;
    ParameterGroupInfo::read(is);
    if (debug) ParameterGroupInfo::debug=0;
    
    if (nPCs>0){
        is>>str;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            pM    = ParameterGroupInfo::read(is,lblPVs(k),pM);    
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pM)<<endl;    k++;
            pDM1  = ParameterGroupInfo::read(is,lblPVs(k),pDM1);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM1)<<endl;  k++;
            pDM2  = ParameterGroupInfo::read(is,lblPVs(k),pDM2);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM2)<<endl;  k++;
            pDM3  = ParameterGroupInfo::read(is,lblPVs(k),pDM3);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM3)<<endl;  k++;
            pDM4 = ParameterGroupInfo::read(is,lblPVs(k),pDM4); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM4)<<endl; k++;
        } else {
            cout<<"Error reading NaturalMortalityInfo from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }//nPCs>0
    
    if (debug){
        cout<<"NMInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished NaturalMortalityInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

void NaturalMortalityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    os<<zRef<<tb<<"#reference size for scaling"<<endl;
    ParameterGroupInfo::write(os);
    
    if (nPCs>0){
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pM)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDM1)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDM2)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDM3)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDM4)<<endl;
    }
}

void NaturalMortalityInfo::writeToPin(std::ostream & os){
    os<<"#--natural mortality parameters"<<endl;
    pM->writeToPin(os);
    pDM1->writeToPin(os);
    pDM2->writeToPin(os);
    pDM3->writeToPin(os);
    pDM4->writeToPin(os);
}

void NaturalMortalityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"nm=list("<<endl;
        ParameterGroupInfo::writeToR(os);    
        if (nPCs>0){
            os<<cc<<endl;
            pM  ->writeToR(os,"pM", indent+1);   os<<cc<<endl;
            pDM1->writeToR(os,"pDM1", indent+1); os<<cc<<endl;
            pDM2->writeToR(os,"pDM2", indent+1); os<<cc<<endl;
            pDM3->writeToR(os,"pDM3", indent+1); os<<cc<<endl;
            pDM4->writeToR(os,"pDM4", indent+1); os<<endl;
        }
    os<<")";
}
        
/*------------------------------------------------------------------------------
 * GrowthInfo
 -----------------------------------------------------------------------------*/
adstring GrowthInfo::NAME = "growth";
GrowthInfo::GrowthInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 2;
    lblIVs.allocate(1,nIVs);
    lblIVs(1) = "YEAR_BLOCK";
    lblIVs(2) = tcsam::STR_SEX;
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs=3;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pGrA";    dscPVs(k++) = "mean growth coefficient 'a'";
    lblPVs(k) = "pGrB";    dscPVs(k++) = "mean growth coefficient 'b'";
    lblPVs(k) = "pGrBeta"; dscPVs(k++) = "growth transition matrix scale factor";
    k=1;
    pGrA    = new BoundedNumberVectorInfo(lblPVs(k++));
    pGrB    = new BoundedNumberVectorInfo(lblPVs(k++));
    pGrBeta = new BoundedNumberVectorInfo(lblPVs(k++));
    
    //define "extra" indices
    nXIs=2;
    lblXIs.allocate(1,nXIs);
    lblXIs(k=1) = "zScaleGrA"; //pre-molt size corresponding to pGrA in alt parameterization
    lblXIs(k=2) = "zScaleGrB"; //pre-molt size corresponding to pGrA in alt parameterization
    
}

GrowthInfo::~GrowthInfo(){
    if (pGrA)    delete pGrA;    pGrA    =0;
    if (pGrB)    delete pGrB;    pGrB    =0;
    if (pGrBeta) delete pGrBeta; pGrBeta =0;
}

/**
 * Set the maximum year for estimating parameters.
 * 
 * @param mxYr - the max year to allow estimated parameters
 */
void GrowthInfo::setMaxYear(int mxYr){
    int k = nIVs+1;
    if (pGrA)    tcsam::adjustParamPhase(pGrA,   this,k++,mxYr,1);
    if (pGrB)    tcsam::adjustParamPhase(pGrB,   this,k++,mxYr,1);
    if (pGrBeta) tcsam::adjustParamPhase(pGrBeta,this,k++,mxYr,1);
}

void GrowthInfo::addNextYearToInfo(int closed){
    if (debug) cout<<"starting void GrowthInfo::project("<<closed<<")"<<endl;
    //find pcs corresponding to mxYr
    int mxYr = ModelConfiguration::mxYr;
    ivector pcs = ParameterGroupInfo::getPCsForYear(mxYr);
    if (debug) {
        if (pcs.size()) cout<<"PCs "<<pcs<<" apply to "<<mxYr<<endl;
        else cout<<"No PCs apply to <<"<<mxYr<<endl;
    }
    for (int i=1;i<=pcs.size();i++){
        ParameterGroupInfo::addYearToPC(pcs[i],mxYr+1);
    }
    if (debug) cout<<"finished void GrowthInfo::project("<<closed<<")"<<endl;
}

void GrowthInfo::read(cifstream & is){
    if (debug) cout<<"starting GrowthInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in GrowthInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    if (nPCs>0){
        is>>str;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            pGrA     = ParameterGroupInfo::read(is,lblPVs(k),pGrA);    
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pGrA)<<endl;    k++;
            pGrB     = ParameterGroupInfo::read(is,lblPVs(k),pGrB);    
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pGrB)<<endl;    k++;
            pGrBeta  = ParameterGroupInfo::read(is,lblPVs(k),pGrBeta); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pGrBeta)<<endl; k++;
         } else {
            cout<<"Error reading GrowthInfo from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }//nPCs>0
    
    if (debug){
        cout<<"GrowthInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished GrowthInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

void GrowthInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
     ParameterGroupInfo::write(os);
    
    if (nPCs>0){
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pGrA)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pGrB)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pGrBeta)<<endl;
    }//nPCs
}

void GrowthInfo::writeToPin(std::ostream & os){
    os<<"#--growth parameters"<<endl;
    pGrA->writeToPin(os);
    pGrB->writeToPin(os);
    pGrBeta->writeToPin(os);
 }

void GrowthInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"grw=list("<<endl;
        ParameterGroupInfo::writeToR(os);         os<<cc<<endl;
        pGrA->writeToR(os,   "pGrA",   indent+1); os<<cc<<endl;
        pGrB->writeToR(os,   "pGrB",   indent+1); os<<cc<<endl;
        pGrBeta->writeToR(os,"pGrBeta",indent+1); os<<endl;
    os<<")";
}
        
/*------------------------------------------------------------------------------
 * MaturityInfo
 -----------------------------------------------------------------------------*/
adstring Molt2MaturityInfo::NAME = "molt_to_maturity";
Molt2MaturityInfo::Molt2MaturityInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 2;
    lblIVs.allocate(1,nIVs);
    lblIVs(1) = "YEAR_BLOCK";
    lblIVs(2) = tcsam::STR_SEX;
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs = 1;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pvLgtPrM2M";dscPVs(k++) = "logit-scale parameter vectors for Pr(molt-to-maturity|size)"; 
    k=1;
    pvLgtPrM2M = new BoundedVectorVectorInfo(lblPVs(k));
    
    nXIs = 0;    
}

Molt2MaturityInfo::~Molt2MaturityInfo(){
    if (pvLgtPrM2M) delete pvLgtPrM2M; pvLgtPrM2M = 0;
}

/**
 * Set the maximum year for estimating parameters.
 * 
 * @param mxYr - the max year to allow estimated parameters
 */
void Molt2MaturityInfo::setMaxYear(int mxYr){
    int k = nIVs+1;
    if (pvLgtPrM2M) tcsam::adjustParamPhase(pvLgtPrM2M,this,k++,mxYr,1);
}

void Molt2MaturityInfo::addNextYearToInfo(int closed){
    if (debug) cout<<"starting void Molt2MaturityInfo::project("<<closed<<")"<<endl;
    //find pcs corresponding to mxYr
    int mxYr = ModelConfiguration::mxYr;
    ivector pcs = ParameterGroupInfo::getPCsForYear(mxYr);
    if (debug) {
        if (pcs.size()) cout<<"PCs "<<pcs<<" apply to "<<mxYr<<endl;
        else cout<<"No PCs apply to <<"<<mxYr<<endl;
    }
    for (int i=1;i<=pcs.size();i++){
        ParameterGroupInfo::addYearToPC(pcs[i],mxYr+1);
    }
    if (debug) cout<<"finished void Molt2MaturityInfo::project("<<closed<<")"<<endl;
}

void Molt2MaturityInfo::read(cifstream & is){
    if (debug) cout<<"starting void MaturityInfo::read(cifstream & is)"<<endl;    
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in MaturityInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    if (nPCs>0){
        is>>str;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            pvLgtPrM2M = ParameterGroupInfo::read(is,lblPVs(k),pvLgtPrM2M);
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pvLgtPrM2M)<<endl;  
        } else {
            cout<<"Error reading MaturityInfo from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }
    
    if (debug){
        cout<<"MaturityInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished MaturityInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void Molt2MaturityInfo::setToWriteVectorInitialValues(bool flag){
    if (pvLgtPrM2M){
        for (int i=1;i<=pvLgtPrM2M->getSize();i++){
            BoundedVectorInfo* vi = (*pvLgtPrM2M)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
}

void Molt2MaturityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    if (nPCs>0){
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pvLgtPrM2M)<<endl;
    }
}

void Molt2MaturityInfo::writeToPin(std::ostream & os){
    os<<"#--molt-to-maturity parameters"<<endl;
    pvLgtPrM2M->writeToPin(os);
}

void Molt2MaturityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"mat=list("<<endl;
        ParameterGroupInfo::writeToR(os);
        if (nPCs>0){
            os<<cc<<endl;
            pvLgtPrM2M->writeToR(os,"pvLgtPrM2M",indent+1); os<<endl;
        }
    os<<")";
}

/*------------------------------------------------------------------------------
 * SelectivityInfo
 -----------------------------------------------------------------------------*/
adstring SelectivityInfo::NAME = "selectivities";
const int SelectivityInfo::nIVs =  1;
const int SelectivityInfo::nPVs = 13;
const int SelectivityInfo::nXIs =  2;
const int SelectivityInfo::idxS1=SelectivityInfo::nIVs+ 1;//parameter combinations index for pS1
const int SelectivityInfo::idxS2=SelectivityInfo::nIVs+ 2;//parameter combinations index for pS2
const int SelectivityInfo::idxS3=SelectivityInfo::nIVs+ 3;//parameter combinations index for pS3
const int SelectivityInfo::idxS4=SelectivityInfo::nIVs+ 4;//parameter combinations index for pS4
const int SelectivityInfo::idxS5=SelectivityInfo::nIVs+ 5;//parameter combinations index for pS5
const int SelectivityInfo::idxS6=SelectivityInfo::nIVs+ 6;//parameter combinations index for pS6
const int SelectivityInfo::idxDevsS1=SelectivityInfo::nIVs+ 7;//parameter combinations index for pDevsS1
const int SelectivityInfo::idxDevsS2=SelectivityInfo::nIVs+ 8;//parameter combinations index for pDevsS2
const int SelectivityInfo::idxDevsS3=SelectivityInfo::nIVs+ 9;//parameter combinations index for pDevsS3
const int SelectivityInfo::idxDevsS4=SelectivityInfo::nIVs+10;//parameter combinations index for pDevsS4
const int SelectivityInfo::idxDevsS5=SelectivityInfo::nIVs+11;//parameter combinations index for pDevsS5
const int SelectivityInfo::idxDevsS6=SelectivityInfo::nIVs+12;//parameter combinations index for pDevsS6
const int SelectivityInfo::idxNPSel =SelectivityInfo::nIVs+13;//parameter combinations index for pvNPSel
SelectivityInfo::SelectivityInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    ParameterGroupInfo::nIVs = SelectivityInfo::nIVs;
    lblIVs.allocate(1,nIVs);
    lblIVs(1) = "YEAR_BLOCK";
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    ParameterGroupInfo::nPVs=SelectivityInfo::nPVs;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pS1"; dscPVs(k++) = "1st input to selectivity function";
    lblPVs(k) = "pS2"; dscPVs(k++) = "2nd input to selectivity function";
    lblPVs(k) = "pS3"; dscPVs(k++) = "3rd input to selectivity function";
    lblPVs(k) = "pS4"; dscPVs(k++) = "4th input to selectivity function";
    lblPVs(k) = "pS5"; dscPVs(k++) = "5th input to selectivity function";
    lblPVs(k) = "pS6"; dscPVs(k++) = "6th input to selectivity function";
    lblPVs(k) = "pDevsS1"; dscPVs(k++) = "devs to 1st input to selectivity function";
    lblPVs(k) = "pDevsS2"; dscPVs(k++) = "devs to 2nd input to selectivity function";
    lblPVs(k) = "pDevsS3"; dscPVs(k++) = "devs to 3rd input to selectivity function";
    lblPVs(k) = "pDevsS4"; dscPVs(k++) = "devs to 4th input to selectivity function";
    lblPVs(k) = "pDevsS5"; dscPVs(k++) = "devs to 5th input to selectivity function";
    lblPVs(k) = "pDevsS6"; dscPVs(k++) = "devs to 6th input to selectivity function";
    lblPVs(k) = "pvNPSel"; dscPVs(k++) = "non-parametric selectivity functions";
    
    k=1;
    pS1 = new BoundedNumberVectorInfo(lblPVs(k++));
    pS2 = new BoundedNumberVectorInfo(lblPVs(k++));
    pS3 = new BoundedNumberVectorInfo(lblPVs(k++));
    pS4 = new BoundedNumberVectorInfo(lblPVs(k++));
    pS5 = new BoundedNumberVectorInfo(lblPVs(k++));
    pS6 = new BoundedNumberVectorInfo(lblPVs(k++));
    pDevsS1 = new DevsVectorVectorInfo(lblPVs(k++));
    pDevsS2 = new DevsVectorVectorInfo(lblPVs(k++));
    pDevsS3 = new DevsVectorVectorInfo(lblPVs(k++));
    pDevsS4 = new DevsVectorVectorInfo(lblPVs(k++));
    pDevsS5 = new DevsVectorVectorInfo(lblPVs(k++));
    pDevsS6 = new DevsVectorVectorInfo(lblPVs(k++));
    pvNPSel = new BoundedVectorVectorInfo(lblPVs(k++));
    
    ParameterGroupInfo::nXIs=SelectivityInfo::nXIs;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "fsZ";
    lblXIs(k++) = "selFcn";
    
}

SelectivityInfo::~SelectivityInfo(){
    if (pS1) delete pS1; pS1=0;
    if (pS2) delete pS2; pS2=0;
    if (pS3) delete pS3; pS3=0;
    if (pS4) delete pS4; pS4=0;
    if (pS5) delete pS5; pS5=0;
    if (pS6) delete pS6; pS6=0;
    if (pDevsS1) delete pDevsS1; pDevsS1=0;
    if (pDevsS2) delete pDevsS2; pDevsS2=0;
    if (pDevsS3) delete pDevsS3; pDevsS3=0;
    if (pDevsS4) delete pDevsS4; pDevsS4=0;
    if (pDevsS5) delete pDevsS5; pDevsS5=0;
    if (pDevsS6) delete pDevsS6; pDevsS6=0;
    if (pvNPSel) delete pvNPSel; pvNPSel=0;
}

/**
 * Set max year for selectivity function identified by pc
 * 
 * @param mxYr - max year to allow
 * @param sfs - ivector indicating fishery, survey, or both
 */
void SelectivityInfo::setMaxYear(int mxYr, const ivector& sfs){
    if (pS1) tcsam::adjustParamPhase(pS1,this,idxS1,mxYr,sfs,1);
    if (pS2) tcsam::adjustParamPhase(pS2,this,idxS2,mxYr,sfs,1);
    if (pS3) tcsam::adjustParamPhase(pS3,this,idxS3,mxYr,sfs,1);
    if (pS4) tcsam::adjustParamPhase(pS4,this,idxS4,mxYr,sfs,1);
    if (pS5) tcsam::adjustParamPhase(pS5,this,idxS5,mxYr,sfs,1);
    if (pS6) tcsam::adjustParamPhase(pS6,this,idxS6,mxYr,sfs,1);
    if (pDevsS1) tcsam::adjustParamPhase(pDevsS1,this,idxDevsS1,mxYr,sfs,1);
    if (pDevsS2) tcsam::adjustParamPhase(pDevsS2,this,idxDevsS2,mxYr,sfs,1);
    if (pDevsS3) tcsam::adjustParamPhase(pDevsS3,this,idxDevsS3,mxYr,sfs,1);
    if (pDevsS4) tcsam::adjustParamPhase(pDevsS4,this,idxDevsS4,mxYr,sfs,1);
    if (pDevsS5) tcsam::adjustParamPhase(pDevsS5,this,idxDevsS5,mxYr,sfs,1);
    if (pDevsS6) tcsam::adjustParamPhase(pDevsS6,this,idxDevsS6,mxYr,sfs,1);
    if (pvNPSel) tcsam::adjustParamPhase(pvNPSel,this,idxNPSel,mxYr,sfs,1);
}

void SelectivityInfo::addNextYear1(adstring flt_type, int y, imatrix& fltpcs){
    //debug=1;
    if (debug) cout<<endl<<"starting SelectivityInfo::addNextYear1('"<<flt_type<<"',"<<y<<")"<<endl;
    
    //define vector for selectivity function pcs that have been updated already
    ivector updated(1,nPCs); int selpc; int fltpc;
    //define ivector to track which vectors in pDevsS1-pDevsS6 have been modified
    ivector idsDevsS1(1,nPCs);
    ivector idsDevsS2(1,nPCs);
    ivector idsDevsS3(1,nPCs);
    ivector idsDevsS4(1,nPCs);
    ivector idsDevsS5(1,nPCs);
    ivector idsDevsS6(1,nPCs);
    idsDevsS1 = 0;
    idsDevsS2 = 0;
    idsDevsS3 = 0;
    idsDevsS4 = 0;
    idsDevsS5 = 0;
    idsDevsS6 = 0;
    int k = 0;//index into "updated" vector
    updated = -1;//set all elements to -1 (i.e., not updated)
    for (int idxS=3;idxS<=4;idxS++){//loop over indices for selectivity/retention or availability/selectivity
        for (int i=fltpcs.indexmin();i<=fltpcs.indexmax();i++){//loop over rows of input imatrix
            fltpc = fltpcs(i,1);   //fleet pc
            selpc = fltpcs(i,idxS);//pc for idxS-type function for ith fleet pc
            if (selpc){
                //check 
                int sumu = 0;
                for (int j=1;j<=nPCs;j++) sumu += (updated(j)==fltpcs(i,idxS));
                if (!sumu){
                    //selectivity function identified by selpc has not been updated yet
                    ivector ids = getPCIDs(selpc);
                    if (debug) {
                        cout<<"will update function pc "<<selpc<<" for function type "<<idxS<<" for fleet pc "<<fltpc<<endl;
                        cout<<"ids for function parameters: "<<ids<<endl;
                    }
                    //add year to YEAR_BLOCK for this function pc
                    addYearToPC(selpc,y);
                    //deal with devs vectors
                    if (pDevsS1->getSize()&&ids[idxDevsS1]){
                        if (debug) cout<<"extending pDevsS1["<<ids[idxDevsS1]<<"]"<<endl<<(*pDevsS1)[ids[idxDevsS1]]->getInitVals()<<endl;
                        //Need to make sure not to addValue multiple times to same devs vector.
                        int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
                        for (int j=1;j<=nPCs;j++) sumdsp += (idsDevsS1[j]==ids[idxDevsS1]);
                        if (!sumdsp){
                            if (debug) cout<<"Adding value to pDevsS1["<<ids[idxDevsS1]<<"]"<<endl;
                            idsDevsS1[i] = ids[idxDevsS1];
                            DevsVectorInfo* pDVI = (*pDevsS1)[ids[idxDevsS1]];
                            if (debug) DevsVectorInfo::debug=1;
                            pDVI->addValueOnParameterScale(0.0,y);
                            if (debug) DevsVectorInfo::debug=0;
                        } else {
                            if (debug) cout<<"Already added value to pDevsS1["<<ids[idxDevsS1]<<"]"<<endl;
                        }
                    }
                    if (pDevsS2->getSize()&&ids[idxDevsS2]){
                        if (debug) cout<<"extending pDevsS2["<<ids[idxDevsS2]<<"]"<<endl<<(*pDevsS2)[ids[idxDevsS2]]->getInitVals()<<endl;
                        //Need to make sure not to addValue multiple times to same devs vector.
                        int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
                        for (int j=1;j<=nPCs;j++) sumdsp += (idsDevsS2[j]==ids[idxDevsS2]);
                        if (!sumdsp){
                            if (debug) cout<<"Adding value to pDevsS2["<<ids[idxDevsS2]<<"]"<<endl;
                            idsDevsS2[i] = ids[idxDevsS2];
                            DevsVectorInfo* pDVI = (*pDevsS2)[ids[idxDevsS2]];
                            if (debug) DevsVectorInfo::debug=1;
                            pDVI->addValueOnParameterScale(0.0,y);
                            if (debug) DevsVectorInfo::debug=0;
                        } else {
                            if (debug) cout<<"Already added value to pDevsS2["<<ids[idxDevsS2]<<"]"<<endl;
                        }
                    }
                    if (pDevsS3->getSize()&&ids[idxDevsS3]){
                        if (debug) cout<<"extending pDevsS3["<<ids[idxDevsS3]<<"]"<<endl<<(*pDevsS3)[ids[idxDevsS3]]->getInitVals()<<endl;
                        //Need to make sure not to addValue multiple times to same devs vector.
                        int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
                        for (int j=1;j<=nPCs;j++) sumdsp += (idsDevsS3[j]==ids[idxDevsS3]);
                        if (!sumdsp){
                            if (debug) cout<<"Adding value to pDevsS3["<<ids[idxDevsS3]<<"]"<<endl;
                            idsDevsS3[i] = ids[idxDevsS3];
                            DevsVectorInfo* pDVI = (*pDevsS3)[ids[idxDevsS3]];
                            if (debug) DevsVectorInfo::debug=1;
                            pDVI->addValueOnParameterScale(0.0,y);
                            if (debug) DevsVectorInfo::debug=0;
                        } else {
                            if (debug) cout<<"Already added value to pDevsS3["<<ids[idxDevsS3]<<"]"<<endl;
                        }
                    }
                    if (pDevsS4->getSize()&&ids[idxDevsS4]){
                        if (debug) cout<<"extending pDevsS4["<<ids[idxDevsS4]<<"]"<<endl<<(*pDevsS4)[ids[idxDevsS4]]->getInitVals()<<endl;
                        //Need to make sure not to addValue multiple times to same devs vector.
                        int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
                        for (int j=1;j<=nPCs;j++) sumdsp += (idsDevsS4[j]==ids[idxDevsS4]);
                        if (!sumdsp){
                            if (debug) cout<<"Adding value to pDevsS4["<<ids[idxDevsS4]<<"]"<<endl;
                            idsDevsS4[i] = ids[idxDevsS4];
                            DevsVectorInfo* pDVI = (*pDevsS4)[ids[idxDevsS4]];
                            if (debug) DevsVectorInfo::debug=1;
                            pDVI->addValueOnParameterScale(0.0,y);
                            if (debug) DevsVectorInfo::debug=0;
                        } else {
                            if (debug) cout<<"Already added value to pDevsS4["<<ids[idxDevsS4]<<"]"<<endl;
                        }
                    }
                    if (pDevsS5->getSize()&&ids[idxDevsS5]){
                        if (debug) cout<<"extending pDevsS5["<<ids[idxDevsS5]<<"]"<<endl<<(*pDevsS5)[ids[idxDevsS5]]->getInitVals()<<endl;
                        //Need to make sure not to addValue multiple times to same devs vector.
                        int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
                        for (int j=1;j<=nPCs;j++) sumdsp += (idsDevsS5[j]==ids[idxDevsS5]);
                        if (!sumdsp){
                            if (debug) cout<<"Adding value to pDevsS5["<<ids[idxDevsS5]<<"]"<<endl;
                            idsDevsS5[i] = ids[idxDevsS5];
                            DevsVectorInfo* pDVI = (*pDevsS5)[ids[idxDevsS5]];
                            if (debug) DevsVectorInfo::debug=1;
                            pDVI->addValueOnParameterScale(0.0,y);
                            if (debug) DevsVectorInfo::debug=0;
                        } else {
                            if (debug) cout<<"Already added value to pDevsS5["<<ids[idxDevsS5]<<"]"<<endl;
                        }
                    }
                    if (pDevsS6->getSize()&&ids[idxDevsS6]){
                        if (debug) cout<<"extending pDevsS6["<<ids[idxDevsS6]<<"]"<<endl<<(*pDevsS6)[ids[idxDevsS6]]->getInitVals()<<endl;
                        //Need to make sure not to addValue multiple times to same devs vector.
                        int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
                        for (int j=1;j<=nPCs;j++) sumdsp += (idsDevsS6[j]==ids[idxDevsS6]);
                        if (!sumdsp){
                            if (debug) cout<<"Adding value to pDevsS6["<<ids[idxDevsS6]<<"]"<<endl;
                            idsDevsS6[i] = ids[idxDevsS6];
                            DevsVectorInfo* pDVI = (*pDevsS6)[ids[idxDevsS6]];
                            if (debug) DevsVectorInfo::debug=1;
                            pDVI->addValueOnParameterScale(0.0,y);
                            if (debug) DevsVectorInfo::debug=0;
                        } else {
                            if (debug) cout<<"Already added value to pDevsS6["<<ids[idxDevsS6]<<"]"<<endl;
                        }
                    }
                    //add selpc to list of updated functions
                    updated(++k) = selpc;
                } else {
                    //already been updated
                    if (debug) cout<<"Function pc "<<selpc<<"for fleet pc "<<fltpc<<" has already been updated "<<endl;
                }
            } else {
                //function of idxS-type not defined for ith fleet pc
                if (debug) cout<<"Function type "<<idxS<<" was not defined for fleet pc "<<fltpc<<endl;
            }
        }//--i
        if (debug) cout<<endl;
    }//--idxS
    if (debug) cout<<"updated SelectivityInfo pcs:"<<endl<<tb<<updated<<endl;
    if (debug) cout<<"finished SelectivityInfo::addNextYear1('"<<flt_type<<"',"<<y<<")"<<endl<<endl;
}

void SelectivityInfo::read(cifstream & is){
    if (debug) cout<<"starting SelectivityInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in SelectivityInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    
    int debugPGI = ParameterGroupInfo::debug;
    if (debug) ParameterGroupInfo::debug=1;
    ParameterGroupInfo::read(is);
    if (debug) ParameterGroupInfo::debug=debugPGI;
    
    if (nPCs>0){
        is>>str;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            pS1 = ParameterGroupInfo::read(is,lblPVs(k),pS1); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS1)<<endl;  k++;
            pS2 = ParameterGroupInfo::read(is,lblPVs(k),pS2); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS2)<<endl;  k++;
            pS3 = ParameterGroupInfo::read(is,lblPVs(k),pS3); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS3)<<endl;  k++;
            pS4 = ParameterGroupInfo::read(is,lblPVs(k),pS4); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS4)<<endl;  k++;
            pS5 = ParameterGroupInfo::read(is,lblPVs(k),pS5); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS5)<<endl;  k++;
            pS6 = ParameterGroupInfo::read(is,lblPVs(k),pS6); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS6)<<endl;  k++;
            pDevsS1 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS1); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS1)<<endl;  k++;
            pDevsS2 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS2); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS2)<<endl;  k++;
            pDevsS3 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS3); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS3)<<endl;  k++;
            pDevsS4 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS4); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS4)<<endl;  k++;
            pDevsS5 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS5); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS5)<<endl;  k++;
            pDevsS6 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS6); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS6)<<endl;  k++;
            pvNPSel = ParameterGroupInfo::read(is,lblPVs(k),pvNPSel); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pvNPSel)<<endl;  k++;
        } else {
            cout<<"Error reading SelectivityInfo from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }//nPCs>0
    
    if (debug){
        cout<<"SelectivityInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished SelectivityInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void SelectivityInfo::setToWriteVectorInitialValues(bool flag){
    if (pDevsS1){
        for (int i=1;i<=pDevsS1->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS1)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS2){
        for (int i=1;i<=pDevsS2->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS2)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS3){
        for (int i=1;i<=pDevsS3->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS3)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS4){
        for (int i=1;i<=pDevsS4->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS4)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS5){
        for (int i=1;i<=pDevsS5->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS5)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS6){
        for (int i=1;i<=pDevsS6->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS6)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pvNPSel){
        for (int i=1;i<=pvNPSel->getSize();i++){
            BoundedVectorInfo* vi = (*pvNPSel)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
}

/**
 * Sets the flags to write estimations phases for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write estimation phases to file
 */
void SelectivityInfo::setToWriteVectorEstimationPhases(bool flag){
    if (pDevsS1){
        for (int i=1;i<=pDevsS1->getSize();i++){
            DevsVectorInfo* vi = (*pDevsS1)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS2){
        for (int i=1;i<=pDevsS2->getSize();i++){
            DevsVectorInfo* vi = (*pDevsS2)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS3){
        for (int i=1;i<=pDevsS3->getSize();i++){
            DevsVectorInfo* vi = (*pDevsS3)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS4){
        for (int i=1;i<=pDevsS4->getSize();i++){
            DevsVectorInfo* vi = (*pDevsS4)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS5){
        for (int i=1;i<=pDevsS5->getSize();i++){
            DevsVectorInfo* vi = (*pDevsS5)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
    if (pDevsS6){
        for (int i=1;i<=pDevsS6->getSize();i++){
            DevsVectorInfo* vi = (*pDevsS6)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
}

void SelectivityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    if (nPCs>0) {
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pS1)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pS2)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pS3)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pS4)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pS5)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pS6)<<endl;

        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsS1)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsS2)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsS3)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsS4)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsS5)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsS6)<<endl;

        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pvNPSel)<<endl;
    }
 }

void SelectivityInfo::writeToPin(std::ostream & os){
    os<<"#--selectivity parameters"<<endl;
    pS1->writeToPin(os);
    pS2->writeToPin(os);
    pS3->writeToPin(os);
    pS4->writeToPin(os);
    pS5->writeToPin(os);
    pS6->writeToPin(os);
    
    pDevsS1->writeToPin(os);
    pDevsS2->writeToPin(os);
    pDevsS3->writeToPin(os);
    pDevsS4->writeToPin(os);
    pDevsS5->writeToPin(os);
    pDevsS6->writeToPin(os);
    
    pvNPSel->writeToPin(os);
 }

void SelectivityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"sel=list("<<endl;
        ParameterGroupInfo::writeToR(os); 
        if (nPCs>0) {
            os<<cc<<endl;
            pS1->writeToR(os,"pS1",indent+1); os<<cc<<endl;
            pS2->writeToR(os,"pS2",indent+1); os<<cc<<endl;
            pS3->writeToR(os,"pS3",indent+1); os<<cc<<endl;
            pS4->writeToR(os,"pS4",indent+1); os<<cc<<endl;
            pS5->writeToR(os,"pS5",indent+1); os<<cc<<endl;
            pS6->writeToR(os,"pS6",indent+1); os<<cc<<endl;
            pDevsS1->writeToR(os,"pDevsS1",indent+1); os<<cc<<endl;
            pDevsS2->writeToR(os,"pDevsS2",indent+1); os<<cc<<endl;
            pDevsS3->writeToR(os,"pDevsS3",indent+1); os<<cc<<endl;
            pDevsS4->writeToR(os,"pDevsS4",indent+1); os<<cc<<endl;
            pDevsS5->writeToR(os,"pDevsS5",indent+1); os<<cc<<endl;
            pDevsS6->writeToR(os,"pDevsS6",indent+1); os<<cc<<endl;
            pvNPSel->writeToR(os,"pvNPSel",indent+1); os<<endl;
        }//nPCs>0
    os<<")";
}

/*------------------------------------------------------------------------------
 * FisheriesInfo
 -----------------------------------------------------------------------------*/
const adstring FisheriesInfo::NAME = "fisheries";
const int FisheriesInfo::nIVs=5;//number of index variables
const int FisheriesInfo::nPVs=9;//number of parameter variables
const int FisheriesInfo::nXIs=3;//number of "extra" variables
const int FisheriesInfo::idxHM     = FisheriesInfo::nIVs+1;//column in parameter combinations matrix with parameter index for column in parameter combinations matrix indicating handling mortality parameters
const int FisheriesInfo::idxLnC    = FisheriesInfo::nIVs+2;//column in parameter combinations matrix with parameter index for ln-scale base mean capture rate (mature males)
const int FisheriesInfo::idxDC1    = FisheriesInfo::nIVs+3;//column in parameter combinations matrix with parameter index for main year_block ln-scale offsets
const int FisheriesInfo::idxDC2    = FisheriesInfo::nIVs+4;//column in parameter combinations matrix with parameter index for ln-scale female offsets
const int FisheriesInfo::idxDC3    = FisheriesInfo::nIVs+5;//column in parameter combinations matrix with parameter index for ln-scale immature offsets
const int FisheriesInfo::idxDC4    = FisheriesInfo::nIVs+6;//column in parameter combinations matrix with parameter index for ln-scale female-immature offsets 
const int FisheriesInfo::idxDevsLnC = FisheriesInfo::nIVs+7;//column in parameter combinations matrix with parameter index for annual ln-scale devs w/in year_blocks
const int FisheriesInfo::idxLnEffX = FisheriesInfo::nIVs+8;//column in parameter combinations matrix with parameter index for ln-scale effort extrapolation 
const int FisheriesInfo::idxLgtRet = FisheriesInfo::nIVs+9;//column in parameter combinations matrix with parameter index for logit-scale retained fraction (for old shell crab)
const int FisheriesInfo::idxSelFcn = FisheriesInfo::nIVs+FisheriesInfo::nPVs+1;//column in parameter combinations matrix indicating selectivity function index
const int FisheriesInfo::idxRetFcn = FisheriesInfo::nIVs+FisheriesInfo::nPVs+2;//column in parameter combinations matrix indicating retention function index
const int FisheriesInfo::idxUseEX  = FisheriesInfo::nIVs+FisheriesInfo::nPVs+3;//column in parameter combinations matrix indicating effort extrapolation use
/*------------------------------------------------------------------------------
 * FisheriesInfo\n
 * Encapsulates the following fishery-related parameters:\n
 *   pHM    : handling mortality (0-1)
 *   pLnC   : ln-scale base mean capture rate (mature males)
 *   pDC1   : ln-scale offset 1 (e.g., main year_block ln-scale offsets)
 *   pDC2   : ln-scale offset 2 (e.g., female offsets)
 *   pDC3   : ln-scale offset 3 (e.g., immature offsets)
 *   pDC4   : ln-scale offset 4 (e.g., female-immature offsets)   
 *
 *   pDevsLnC : annual ln-scale devs w/in year_blocks
 * 
 *   pLnEffX: ln-scale effort extrapolation parameters
 *   pLgtRet: logit-scale max retention parameter
 *
 * Notes:
 *  1. FISHERY, YEAR_BLOCK, SEX, MATURITY_STATE, SHELL_CONDITION are the index variables for the parameters
*----------------------------------------------------------------------------*/
FisheriesInfo::FisheriesInfo(){
    if (debug) cout<<"starting FisheriesInfo::FisheriesInfo()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    //create "independent variables" for parameter group assignment
    ParameterGroupInfo::nIVs = FisheriesInfo::nIVs;//number of independent variables
    lblIVs.allocate(1,nIVs);
    k=1;
    lblIVs(k++) = tcsam::STR_FISHERY;
    lblIVs(k++) = "YEAR_BLOCK";
    lblIVs(k++) = tcsam::STR_SEX;
    lblIVs(k++) = tcsam::STR_MATURITY_STATE;
    lblIVs(k++) = tcsam::STR_SHELL_CONDITION;
    //create index block sets for "BLOCKS" (e.g. year blocks)
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 2; //YEAR_BLOCK
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    ParameterGroupInfo::nPVs=FisheriesInfo::nPVs;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pHM";      dscPVs(k++) = "handling mortality (0-1)";
    lblPVs(k) = "pLnC";     dscPVs(k++) = "ln-scale base mean capture rate (mature male crab)";
    lblPVs(k) = "pDC1";     dscPVs(k++) = "ln-scale capture rate offset 1";
    lblPVs(k) = "pDC2";     dscPVs(k++) = "ln-scale capture rate offset 2";
    lblPVs(k) = "pDC3";     dscPVs(k++) = "ln-scale capture rate offset 3";
    lblPVs(k) = "pDC4";     dscPVs(k++) = "ln-scale capture rate offset 4";
    lblPVs(k) = "pDevsLnC"; dscPVs(k++) = "ln-scale annual capture rate devs";
    lblPVs(k) = "pLnEffX";  dscPVs(k++) = "ln-scale effort extrapolation parameters";
    lblPVs(k) = "pLgtRet";  dscPVs(k++) = "logit-scale retained fraction parameters";
    k=1;
    pHM   = new BoundedNumberVectorInfo(lblPVs(k++));
    pLnC  = new BoundedNumberVectorInfo(lblPVs(k++));    
    pDC1  = new BoundedNumberVectorInfo(lblPVs(k++));    
    pDC2  = new BoundedNumberVectorInfo(lblPVs(k++));    
    pDC3  = new BoundedNumberVectorInfo(lblPVs(k++));    
    pDC4  = new BoundedNumberVectorInfo(lblPVs(k++));    
    pDevsLnC = new DevsVectorVectorInfo(lblPVs(k++));    
    pLnEffX  = new BoundedNumberVectorInfo(lblPVs(k++)); 
    pLgtRet  = new BoundedNumberVectorInfo(lblPVs(k++)); 
    
    //create "extra indices"
    ParameterGroupInfo::nXIs=FisheriesInfo::nXIs;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "idx.SelFcn";
    lblXIs(k++) = "idx.RetFcn";
    lblXIs(k++) = "useEffX";   
    
    if (debug) cout<<"finished FisheriesInfo::FisheriesInfo()"<<endl;
}

FisheriesInfo::~FisheriesInfo(){
    if (pHM)  delete pHM;   pHM=0;
    if (pLnC) delete pLnC;  pLnC=0;
    if (pDC1) delete pDC1;  pDC1=0;
    if (pDC2) delete pDC2;  pDC2=0;
    if (pDC3) delete pDC3;  pDC3=0;
    if (pDC4) delete pDC4;  pDC4=0;
    if (pDevsLnC) delete pDevsLnC; pDevsLnC=0;
    if (pLnEffX)  delete pLnEffX;  pLnEffX=0;
    if (pLgtRet)  delete pLgtRet;  pLgtRet=0;
}

/**
 * Set the maximum year for estimating parameters.
 * 
 * @param mxYr - the max year to allow estimated parameters
 */
void FisheriesInfo::setMaxYear(int mxYr){
    int k = nIVs+1;
    if (pHM)  tcsam::adjustParamPhase(pHM,this,k++,mxYr,1);
    if (pLnC) tcsam::adjustParamPhase(pLnC,this,k++,mxYr,1);
    if (pDC1) tcsam::adjustParamPhase(pDC1,this,k++,mxYr,1);
    if (pDC2) tcsam::adjustParamPhase(pDC2,this,k++,mxYr,1);
    if (pDC3) tcsam::adjustParamPhase(pDC3,this,k++,mxYr,1);
    if (pDC4) tcsam::adjustParamPhase(pDC4,this,k++,mxYr,1);
    if (pDevsLnC) tcsam::adjustParamPhase(pDevsLnC,mxYr,1); k++;
    if (pLnEffX)  tcsam::adjustParamPhase(pLnEffX,this,k++,mxYr,1);
    if (pLgtRet)  tcsam::adjustParamPhase(pLgtRet,this,k++,mxYr,1);
}

imatrix FisheriesInfo::addNextYear1(int closed){
    //debug=1;
    if (debug) cout<<"starting void FisheriesInfo::addNextYear1(closed="<<closed<<")"<<endl;
    //find pcs corresponding to mxYr
    int mxYr = ModelConfiguration::mxYr;
    ivector pcs = ParameterGroupInfo::getPCsForYear(mxYr);
    if (debug) {
        if (pcs.size()) cout<<pcs.size()<<" PCs ("<<pcs<<") apply to "<<mxYr<<endl;
        else cout<<"No PCs apply to <<"<<mxYr<<endl;
    }
    //define ivector to track which vectors in pDevsLnC have been modified
    ivector idsDevsLnC(pcs.indexmin(),pcs.indexmax());
    idsDevsLnC = 0;
    //define imatrix to track associated selectivity and retention functions
    imatrix selpcs(pcs.indexmin(),pcs.indexmax(),1,4); selpcs.initialize();
    for (int i=pcs.indexmin();i<=pcs.indexmax();i++){
        ivector ids = ParameterGroupInfo::getPCIDs(pcs[i]);
        if ((!closed)||(ids[1]!=1)){
            //directed fishery is not closed, or this pc does not pertain to directed fishery
            ParameterGroupInfo::addYearToPC(pcs[i],mxYr+1);
            if (debug) {
                cout<<endl;
                cout<<"Extending PC "<<pcs[i]<<", which pertains to fishery = "<<ids[1]<<endl;
                cout<<"PC IDs for pc "<<pcs[i]<<" = "<<ids<<endl;
                cout<<"will extend pDevsLnC["<<ids[idxDevsLnC]<<"]"<<endl<<(*pDevsLnC)[ids[idxDevsLnC]]->getInitVals()<<endl;
            }
            if (pDevsLnC->getSize()){
                //Need to make sure not to addValue multiple times to same devs vector.
                int sumdsp = 0;//if sumdsp>0, then already added value to indexed devs vector
                for (int j=1;j<=pcs.size();j++) sumdsp += (idsDevsLnC[j]==ids[idxDevsLnC]);
                if (!sumdsp){
                    if (debug) cout<<"Adding value to pDevsLnC["<<ids[idxDevsLnC]<<"]"<<endl;
                    idsDevsLnC[i] = ids[idxDevsLnC];
                    DevsVectorInfo* pDVI = (*pDevsLnC)[ids[idxDevsLnC]];
                    if (debug) DevsVectorInfo::debug=1;
                    pDVI->addValueOnParameterScale(0.0,mxYr+1);
                    if (debug) DevsVectorInfo::debug=0;
                } else {
                    if (debug) cout<<"Already added value to pDevsLnC["<<ids[idxDevsLnC]<<"]"<<endl;
                }
            }
            //add selectivity/retention info to selpcs
            selpcs(i,1) = pcs[i];//fishery pc
            selpcs(i,2) = ids[1];//fishery id
            selpcs(i,3) = ids[idxSelFcn];//pc identifying selectivity function
            selpcs(i,4) = ids[idxRetFcn];//pc identifying retention function
        } else {
            if (debug) {
                cout<<endl;
                cout<<"Not extending PC "<<pcs[i]<<", which pertains to fishery = "<<ids[1]<<endl;
                cout<<"PC IDs for pc "<<pcs[i]<<" = "<<ids<<endl;
                cout<<"This is the directed fishery"<<endl;
                cout<<"Not adding additional year for this pc or pDevsLnC["<<ids[idxDevsLnC]<<"]."<<endl;
            }
        }
    }
    if (1) cout<<endl<<"selpcs: "<<endl<<selpcs<<endl<<endl;
    if (debug) cout<<"finished void FisheriesInfo::addNextYear1(closed="<<closed<<")"<<endl;
    return(selpcs);
}

void FisheriesInfo::read(cifstream & is){
    if (debug) cout<<"starting FisheriesInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in FisheriesInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    if (nPCs>0){
        is>>str;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            pHM     = ParameterGroupInfo::read(is,lblPVs(k),pHM);     
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pHM)<<endl;  k++;
            pLnC    = ParameterGroupInfo::read(is,lblPVs(k),pLnC);    
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnC)<<endl;  k++;
            pDC1  = ParameterGroupInfo::read(is,lblPVs(k),pDC1);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC1)<<endl;  k++;
            pDC2  = ParameterGroupInfo::read(is,lblPVs(k),pDC2);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC2)<<endl;  k++;
            pDC3  = ParameterGroupInfo::read(is,lblPVs(k),pDC3);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC3)<<endl;  k++;
            pDC4 = ParameterGroupInfo::read(is,lblPVs(k),pDC4); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC4)<<endl;  k++;
            pDevsLnC = ParameterGroupInfo::read(is,lblPVs(k),pDevsLnC); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsLnC)<<endl;  k++;
            pLnEffX = ParameterGroupInfo::read(is,lblPVs(k),pLnEffX); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnEffX)<<endl;  k++;
            pLgtRet = ParameterGroupInfo::read(is,lblPVs(k),pLgtRet); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLgtRet)<<endl;  k++;
        } else {
            cout<<"Error reading FisheriesInfo from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }//nPCs>0
    
    if (debug){
        cout<<"FisheriesInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished FisheriesInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write estimation phases for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write estimation phases to file
 */
void FisheriesInfo::setToWriteVectorEstimationPhases(bool flag){
    if (pDevsLnC){
        for (int i=1;i<=pDevsLnC->getSize();i++){
            DevsVectorInfo* vi = (*pDevsLnC)[i];
            vi->readPhases = flag ? INT_TRUE : INT_FALSE;        
        }
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void FisheriesInfo::setToWriteVectorInitialValues(bool flag){
    if (pDevsLnC){
        for (int i=1;i<=pDevsLnC->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsLnC)[i];
            vi->readVals = flag ? INT_TRUE : INT_FALSE;        
        }
    }
}

void FisheriesInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    if (nPCs>0){
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pHM)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pLnC)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDC1)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDC2)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDC3)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDC4)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDevsLnC)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pLnEffX)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pLgtRet)<<endl;
    }//nPCs>0
 }

void FisheriesInfo::writeToPin(std::ostream & os){
    os<<"#--fisheries parameters"<<endl;
    pHM->writeToPin(os);
    pLnC->writeToPin(os);
    pDC1->writeToPin(os);
    pDC2->writeToPin(os);
    pDC3->writeToPin(os);
    pDC4->writeToPin(os);
    
    pDevsLnC->writeToPin(os);
    
    pLnEffX->writeToPin(os);
    pLgtRet->writeToPin(os);
}

void FisheriesInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"fsh=list("<<endl;
        ParameterGroupInfo::writeToR(os);
        if (nPCs>0) {
            os<<cc<<endl;
            pLnC->writeToR(os, "pLnC", indent++); os<<cc<<endl;
            pDC1->writeToR(os, "pDC1", indent++); os<<cc<<endl;
            pDC2->writeToR(os, "pDC2", indent++); os<<cc<<endl;
            pDC3->writeToR(os, "pDC3", indent++); os<<cc<<endl;
            pDC4->writeToR(os, "pDC4", indent++); os<<cc<<endl;
            pLnEffX->writeToR(os, "pLnEffX", indent++); os<<cc<<endl;
            pLgtRet->writeToR(os, "pLgtRet", indent++); os<<cc<<endl;
            pDevsLnC->writeToR(os,"pDevsLnC",indent++); os<<endl;
        }
    os<<")";
}

/*------------------------------------------------------------------------------
 * SurveysInfo\n
 * Encapsulates the following recruitment-related parameters:\n
 *   pQ   : base q (mature males)
 *   pDQ1 : ln-scale offset 1 (e.g., main temporal offset)
 *   pDQ2 : ln-scale offset 1 (e.g., female offsets)
 *   pDQ3 : ln-scale offset 1 (e.g., immature offsets)
 *   pDQ4 : ln-scale offset 1 (e.g., female-immature offsets)
 * Notes:
 *  1. index variables for parameters:
 *      a. SURVEY
 *      b. YEAR_BLOCK
 *      c. SEX
 *      d. MATURITY
 *      e. SHELL
*----------------------------------------------------------------------------*/
const adstring SurveysInfo::NAME = "surveys";
const int SurveysInfo::nIVs=5;//number of index variables
const int SurveysInfo::nPVs=6;//number of parameter variables
const int SurveysInfo::nXIs=2;//number of "extra" variables
const int SurveysInfo::idxAvlFcn = SurveysInfo::nIVs+SurveysInfo::nPVs+1;
const int SurveysInfo::idxSelFcn = SurveysInfo::nIVs+SurveysInfo::nPVs+2;
SurveysInfo::SurveysInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    ParameterGroupInfo::nIVs = SurveysInfo::nIVs;
    lblIVs.allocate(1,nIVs);
    k=1;
    lblIVs(k++) = tcsam::STR_SURVEY;
    lblIVs(k++) = "YEAR_BLOCK";
    lblIVs(k++) = tcsam::STR_SEX;
    lblIVs(k++) = tcsam::STR_MATURITY_STATE;
    lblIVs(k++) = tcsam::STR_SHELL_CONDITION;
    
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 2;//index for YEAR_BLOCK
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    ParameterGroupInfo::nPVs=SurveysInfo::nPVs;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pQ";   dscPVs(k++) = "base catchability (e.g. mature male crab)";
    lblPVs(k) = "pDQ1"; dscPVs(k++) = "offset 1 for ln-scale catchability";
    lblPVs(k) = "pDQ2"; dscPVs(k++) = "offset 2 for ln-scale catchability";
    lblPVs(k) = "pDQ3"; dscPVs(k++) = "offset 3 for ln-scale catchability";
    lblPVs(k) = "pDQ4"; dscPVs(k++) = "offset 4 for ln-scale catchability";
    lblPVs(k) = "pA";   dscPVs(k++) = "max availability";
    k=1;
    pQ   = new BoundedNumberVectorInfo(lblPVs(k++));
    pDQ1 = new BoundedNumberVectorInfo(lblPVs(k++));
    pDQ2 = new BoundedNumberVectorInfo(lblPVs(k++));
    pDQ3 = new BoundedNumberVectorInfo(lblPVs(k++));
    pDQ4 = new BoundedNumberVectorInfo(lblPVs(k++));
    pA   = new BoundedNumberVectorInfo(lblPVs(k++));
    
    ParameterGroupInfo::nXIs=SurveysInfo::nXIs;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "idx.AvlFcn";
    lblXIs(k++) = "idx.SelFcn";
    
}

SurveysInfo::~SurveysInfo(){
    if (pQ)   delete pQ;   pQ   =0;
    if (pDQ1) delete pDQ1; pDQ1 =0;
    if (pDQ2) delete pDQ2; pDQ2 =0;
    if (pDQ3) delete pDQ3; pDQ3 =0;
    if (pDQ4) delete pDQ4; pDQ4 =0;
    if (pA)   delete pA;   pA =0;
}

/**
 * Set the maximum year for estimating parameters.
 * 
 * @param mxYr - the max year to allow estimated parameters
 */
void SurveysInfo::setMaxYear(int mxYr){
    int k = nIVs+1;
    if (pQ)   tcsam::adjustParamPhase(pQ,this,k++,mxYr,1);
    if (pDQ1) tcsam::adjustParamPhase(pDQ1,this,k++,mxYr,1);
    if (pDQ2) tcsam::adjustParamPhase(pDQ2,this,k++,mxYr,1);
    if (pDQ3) tcsam::adjustParamPhase(pDQ3,this,k++,mxYr,1);
    if (pDQ4) tcsam::adjustParamPhase(pDQ4,this,k++,mxYr,1);
    if (pA)   tcsam::adjustParamPhase(pA,this,k++,mxYr,1);
    
    //need to do something here for selectivity functions
}

imatrix SurveysInfo::addNextYear1(int closed){
    //debug=1;
    if (debug) cout<<"starting void SurveysInfo::addNextYear1("<<closed<<")"<<endl;
    //find pcs corresponding to mxYr
    int mxYr = ModelConfiguration::mxYr+1;
    ivector pcs = ParameterGroupInfo::getPCsForYear(mxYr);
    if (debug) {
        if (pcs.size()) cout<<"PCs "<<pcs<<" apply to "<<mxYr<<endl;
        else cout<<"No PCs apply to <<"<<mxYr<<endl;
    }
    //define imatrix to track associated selectivity and retention functions
    imatrix selpcs(pcs.indexmin(),pcs.indexmax(),1,4); selpcs.initialize();
    for (int i=pcs.indexmin();i<=pcs.indexmax();i++){
        ivector ids = ParameterGroupInfo::getPCIDs(pcs[i]);
        ParameterGroupInfo::addYearToPC(pcs[i],mxYr+1);
        if (debug) {
            cout<<endl;
            cout<<"Extending PC "<<pcs[i]<<", which pertains to survey = "<<ids[1]<<endl;
            cout<<"PC IDs for pc "<<pcs[i]<<" = "<<ids<<endl;
        }
        //add selectivity/retention info to selpcs
        selpcs(i,1) = pcs[i];//survey pc
        selpcs(i,2) = ids[1];//survey id
        selpcs(i,3) = ids[idxAvlFcn];//pc identifying associated availability function
        selpcs(i,4) = ids[idxSelFcn];//pc identifying associated selectivity function
    }
    if (debug) cout<<endl<<"selpcs: "<<endl<<selpcs<<endl<<endl;
    if (debug) cout<<"finished void SurveysInfo::addNextYear1("<<closed<<")"<<endl;
    return(selpcs);
}

void SurveysInfo::read(cifstream & is){
    if (debug) cout<<"starting SurveysInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in SurveysInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    if (nPCs>0){
        is>>str;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            pQ = ParameterGroupInfo::read(is,lblPVs(k),pQ);    
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pQ)<<endl;  k++;
            pDQ1 = ParameterGroupInfo::read(is,lblPVs(k),pDQ1);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ1)<<endl;  k++;
            pDQ2 = ParameterGroupInfo::read(is,lblPVs(k),pDQ2);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ2)<<endl;  k++;
            pDQ3 = ParameterGroupInfo::read(is,lblPVs(k),pDQ3);  
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ3)<<endl;  k++;
            pDQ4 = ParameterGroupInfo::read(is,lblPVs(k),pDQ4); 
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ4)<<endl;  k++;
            pA = ParameterGroupInfo::read(is,lblPVs(k),pA);    
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pA)<<endl;  k++;
        } else {
            cout<<"Error reading SurveysInfo from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }
    
    if (debug){
        cout<<"SurveysInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished SurveysInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

void SurveysInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    if (nPCs>0){
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pQ)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDQ1)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDQ2)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDQ3)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pDQ4)<<endl;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pA)<<endl;
    }//nPCs>0
 }

void SurveysInfo::writeToPin(std::ostream & os){
    os<<"#--surveys parameters"<<endl;
    pQ->writeToPin(os);
    pDQ1->writeToPin(os);
    pDQ2->writeToPin(os);
    pDQ3->writeToPin(os);
    pDQ4->writeToPin(os);
    pA->writeToPin(os);
 }

void SurveysInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"srv=list("<<endl;
        ParameterGroupInfo::writeToR(os);
        if (nPCs>0) {
            os<<cc<<endl;
            pQ  ->writeToR(os,"pQ",indent++);   os<<cc<<endl;
            pDQ1->writeToR(os,"pDQ1",indent++); os<<cc<<endl;
            pDQ2->writeToR(os,"pDQ2",indent++); os<<cc<<endl;
            pDQ3->writeToR(os,"pDQ3",indent++); os<<cc<<endl;
            pDQ4->writeToR(os,"pDQ4",indent++); os<<cc<<endl;
            pA  ->writeToR(os,"pA",indent++);   os<<endl;
        }
    os<<")";
}

/*------------------------------------------------------------------------------
 * MSE_Info
 -----------------------------------------------------------------------------*/
adstring MSE_Info::NAME = "MSE";
int MSE_Info::idxMSE_LnC  = 1+1;//column in parameter combinations matrix with parameter index for capture rate for MSE op mod 
/*------------------------------------------------------------------------------
 * MSE_Info\n
 * Encapsulates the following MSE-related parameters:\n
 *   pMSE_LnC : ln-scale capture rate for MSE op mod   
 *
 * Notes:
 *  1. FISHERY is the only index variable for the parameters
*----------------------------------------------------------------------------*/
MSE_Info::MSE_Info(){
    if (debug) cout<<"starting MSE_Info::MSE_Info()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    //create "independent variables" for parameter group assignment
    nIVs = 1;//number of independent variables
    lblIVs.allocate(1,nIVs);
    k=1;
    lblIVs(k++) = tcsam::STR_FISHERY;
    //create index block sets for "BLOCKS" (e.g. year blocks)
    nIBSs = 0;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs=1;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pMSE_LnC"; dscPVs(k++) = "ln-scale capture rate for directed fishery in MSE op mod";
    k=1;
    pMSE_LnC = new BoundedNumberVectorInfo(lblPVs(k));
    
    //create "extra indices"
    nXIs=0;
    
    if (debug) cout<<"finished MSE_Info::MSE_Info()"<<endl;
}

MSE_Info::~MSE_Info(){
    if (pMSE_LnC) delete pMSE_LnC;  pMSE_LnC=0;
}

void MSE_Info::addNextYearToInfo(int closed){
    //no year block to project from
}

void MSE_Info::read(cifstream & is){
    if (debug) cout<<"starting MSE_Info::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in MSE_Info::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    if (debug) {
        cout<<"MSE_Info starting ParameterGroupInfo::read(is)"<<endl;
        ParameterGroupInfo::debug=1;
    }
    ParameterGroupInfo::read(is);
    if (debug) {
        ParameterGroupInfo::debug=0;
        cout<<"MSE_Info finished ParameterGroupInfo::read(is)"<<endl;
    }
    
    if (nPCs>0) {
        is>>str;
        if (debug) cout<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
        if (str=="PARAMETERS"){
            int k=1;
            if (debug) cout<<"MSE_Info::read--starting PGI::read(is,lblPVs(k),pMSE_LnC"<<endl;
            pMSE_LnC = ParameterGroupInfo::read(is,lblPVs(k),pMSE_LnC);    
            if (debug) {
                cout<<"MSE_Info::read--finished PGI::read(is,lblPVs(k),pMSE_LnC"<<endl;
                cout<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; cout<<(*pMSE_LnC)<<endl;
            }
            rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pMSE_LnC)<<endl;  k++;
        } else {
            cout<<"Error reading MSE_Info from "<<is.get_file_name()<<endl;
            cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }//nPCs
    if (debug) cout<<"finished MSE_Info::read(is)"<<endl;
}

void MSE_Info::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    if (nPCs>0){
        os<<"PARAMETERS"<<endl;

        int k=1;
        os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
        os<<(*pMSE_LnC)<<endl;
    }
 }

void MSE_Info::writeToPin(std::ostream & os){
    os<<"#--MSE parameters"<<endl;
    pMSE_LnC->writeToPin(os);
 }

void MSE_Info::writeToR(std::ostream & os){
    int indent=0;
    os<<"mse=list("<<endl;
        ParameterGroupInfo::writeToR(os);
        if (nPCs>0) {
            os<<cc<<endl;
            pMSE_LnC->writeToR(os, "pMSE_LnC", indent++); os<<endl;
        }
    os<<")";
}

/*------------------------------------------------------------------------------
 * ModelParametersInfo
 -----------------------------------------------------------------------------*/
ModelParametersInfo::ModelParametersInfo(ModelConfiguration& mc){
    ptrMC=&mc;
    ptrRec=0;
    ptrNM =0;
    ptrGrw=0;
    ptrM2M=0;
    ptrSel=0;
    ptrFsh=0;
    ptrSrv=0;
    ptrMSE=0;
}

ModelParametersInfo::~ModelParametersInfo(){
    if (ptrRec) {delete ptrRec;     ptrRec  = 0;}
    if (ptrNM)  {delete ptrNM;      ptrNM   = 0;}
    if (ptrGrw) {delete ptrGrw;     ptrGrw  = 0;}
    if (ptrM2M) {delete ptrM2M;     ptrM2M  = 0;}
    if (ptrSel) {delete ptrSel;     ptrSel  = 0;}
    if (ptrFsh) {delete ptrFsh;     ptrFsh  = 0;}
    if (ptrSrv) {delete ptrSrv;     ptrSrv  = 0;}
    if (ptrMSE) {delete ptrMSE;     ptrMSE  = 0;}
}

/**
 * Set the maximum year for estimating parameters.
 * 
 * @param mxYr - the max fishery year to allow estimated parameters
 */
void ModelParametersInfo::setMaxYear(int mxYr){
    if (ptrRec) ptrRec->setMaxYear(mxYr);
    if (ptrNM)  ptrNM ->setMaxYear(mxYr);
    if (ptrGrw) ptrGrw->setMaxYear(mxYr+1); //need +1 in case one has growth info from current surveys
    if (ptrM2M) ptrM2M->setMaxYear(mxYr+1); //need +1 in case one has maturity info from current surveys
    if (ptrFsh) ptrFsh->setMaxYear(mxYr);
    if (ptrSrv) ptrSrv->setMaxYear(mxYr+1); //need +1 for current surveys
    if (ptrSel){
        int nSelPCs = ptrSel->nPCs;
        ivector sfs(1,nSelPCs); sfs.initialize();
        //only need to loop over survey-related sel/avl functions
        for (int i=1;i<=ptrSrv->nPCs;i++){
            int selpc = ptrSrv->getPCIDs(i)[SurveysInfo::idxSelFcn];
            int avlpc = ptrSrv->getPCIDs(i)[SurveysInfo::idxAvlFcn];
            if (selpc) sfs[selpc] = 1;
            if (avlpc) sfs[avlpc] = 1;
        }
        rpt::echo<<"sfs = "<<sfs<<endl;
        ptrSel->setMaxYear(mxYr, sfs);
    }
    //if (ptrSel) ptrSel->setMaxYear(mxYr); //handled by FisheriesInfo, SurveysInfo
    //if (ptrMSE) ptrMSE->setMaxYear(mxYr); //not needed for MSE_Info
}

/**
 * Sets the flags to write estimation phases for all vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void ModelParametersInfo::setToWriteVectorEstimationPhases(bool flag){
    if (ptrRec) ptrRec->setToWriteVectorEstimationPhases(flag);
    if (ptrNM)  ptrNM ->setToWriteVectorEstimationPhases(flag); //not needed
    if (ptrGrw) ptrGrw->setToWriteVectorEstimationPhases(flag); //not needed
    if (ptrM2M) ptrM2M->setToWriteVectorEstimationPhases(flag); //not needed
    if (ptrSel) ptrSel->setToWriteVectorEstimationPhases(flag);
    if (ptrFsh) ptrFsh->setToWriteVectorEstimationPhases(flag);
    if (ptrSrv) ptrSrv->setToWriteVectorEstimationPhases(flag); //not needed for SurveysInfo
    if (ptrMSE) ptrMSE->setToWriteVectorEstimationPhases(flag); //not needed for MSE_Info
}

/**
 * Sets the flags to write initial values for all vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void ModelParametersInfo::setToWriteVectorInitialValues(bool flag){
    if (ptrRec) ptrRec->setToWriteVectorInitialValues(flag);
    if (ptrNM)  ptrNM ->setToWriteVectorInitialValues(flag);
    if (ptrGrw) ptrGrw->setToWriteVectorInitialValues(flag);
    if (ptrM2M) ptrM2M->setToWriteVectorInitialValues(flag);
    if (ptrSel) ptrSel->setToWriteVectorInitialValues(flag);
    if (ptrFsh) ptrFsh->setToWriteVectorInitialValues(flag);
    if (ptrSrv) ptrSrv->setToWriteVectorInitialValues(flag); //not needed for SurveysInfo
    if (ptrMSE) ptrMSE->setToWriteVectorInitialValues(flag); //not needed for MSE_Info
}

void ModelParametersInfo::read(cifstream & is){
    if (debug) cout<<"starting void ModelParametersInfo::read(cifstream & is)"<<endl;
    
    rpt::echo<<"#------Reading Model Parameters Info file-----"<<endl;
    rpt::echo<<"file is '"<<is.get_file_name()<<"'."<<endl;
    
    adstring str;
    is>>str;
    if (str!=ModelParametersInfo::version){
        cout<<"-----------------------------------------------------------"<<endl;
        cout<<"Model Parameters Info file '"<<is.get_file_name()<<"'"<<endl;
        cout<<"has incorrect version number!!!"<<endl;
        cout<<"Expected '"<<ModelParametersInfo::version<<"'. Got '"<<str<<"'."<<endl;
        cout<<"Terminating run..."<<endl;
        exit(-1);
    }
    
    //read recruitment parameters
    rpt::echo<<"#---reading  Recruitment Info"<<endl;
    ptrRec = new RecruitmentInfo();
    is>>(*ptrRec);
    rpt::echo<<"#---created  RecruitmentInfo object"<<endl;
    
    //read natural mortality parameters
    rpt::echo<<"#---reading Natural Mortality Info"<<endl;
    ptrNM = new NaturalMortalityInfo();
    is>>(*ptrNM);
    if (debug) cout<<"created NaturalMortalityInfo object"<<endl;
    
    //read growth parameters
    rpt::echo<<"#---reading Growth Info"<<endl;
    ptrGrw = new GrowthInfo();
    is>>(*ptrGrw);
    if (debug) cout<<"created GrowthInfo object"<<endl;
    
    //read maturity parameters
    rpt::echo<<"#---reading Molt2Maturity Info"<<endl;
    ptrM2M = new Molt2MaturityInfo();
    is>>(*ptrM2M);
    if (debug) cout<<"created Molt2MaturityInfo object"<<endl;
    
    //read selectivity function parameters
    rpt::echo<<"#---reading Selectivity Info"<<endl;
    ptrSel = new SelectivityInfo();
    is>>(*ptrSel);
    if (debug) cout<<"created SelectivityInfo object"<<endl;
    
    //read fisheries parameters
    rpt::echo<<"#---reading Fisheries Info"<<endl;
    ptrFsh = new FisheriesInfo();
    is>>(*ptrFsh);
    if (debug) cout<<"created FisheriesInfo object"<<endl;
    
    //read surveys parameters
    rpt::echo<<"#---reading Surveys Info"<<endl;
    ptrSrv = new SurveysInfo();
    is>>(*ptrSrv);
    rpt::echo<<"#---read Surveys Info"<<endl;
    
    //read MSE-related parameters
    rpt::echo<<"#---reading MSE Info"<<endl;
    ptrMSE = new MSE_Info();
    is>>(*ptrMSE);
    rpt::echo<<"#---read MSE Info"<<endl;
    
    if (debug) cout<<"finished void ModelParametersInfo::read(cifstream & is)"<<endl;
}
    
void ModelParametersInfo::write(std::ostream & os){
    os<<version<<tb<<"#model parameters info version"<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Recruitment parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrRec)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Natural mortality parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrNM)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Growth parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrGrw)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Molt-to-maturity parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrM2M)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Selectivity parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrSel)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Fisheries parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrFsh)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Surveys parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrSrv)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# MSE-related parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrMSE)<<endl;
}

void ModelParametersInfo::writePin(std::ostream & os){
    os<<"#-------------------------------"<<endl;
    os<<"# Recruitment parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrRec->writeToPin(os); os<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Natural mortality parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrNM->writeToPin(os); os<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Growth parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrGrw->writeToPin(os); os<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Molt-to-maturity parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrM2M->writeToPin(os); os<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Selectivity parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrSel->writeToPin(os); os<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Fisheries parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrFsh->writeToPin(os); os<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Surveys parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrSrv->writeToPin(os); os<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# MSE-related parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    ptrMSE->writeToPin(os); os<<endl;
}

void ModelParametersInfo::addNextYearToInfo(int closed){
    //recruitment info
    ptrRec->addNextYearToInfo(closed);
    
    //natural mortality info
    ptrNM->addNextYearToInfo(closed);
    
    //growth info
    ptrGrw->addNextYearToInfo(closed);
    
    //molt-to-maturity info
    ptrM2M->addNextYearToInfo(closed);
    
    //fisheries info
    imatrix fshpcs = ptrFsh->addNextYear1(closed);
    //selectivity info for fisheries
    ptrSel->addNextYear1("fisheries",ModelConfiguration::mxYr+1,fshpcs);
    
    //surveys info
    imatrix srvpcs = ptrSrv->addNextYear1(closed);    
    //selectivity info for surveys
    ptrSel->addNextYear1("surveys",ModelConfiguration::mxYr+2,srvpcs);
    
    //MSE info
    ptrMSE->addNextYearToInfo(closed);
}

void ModelParametersInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"mpi=list(version='"<<version<<"'"<<cc<<endl;
    ptrRec->writeToR(os); os<<cc<<endl;
    ptrNM->writeToR(os);  os<<cc<<endl;
    ptrGrw->writeToR(os); os<<cc<<endl;
    ptrM2M->writeToR(os); os<<cc<<endl;
    ptrSel->writeToR(os); os<<cc<<endl;
    ptrFsh->writeToR(os); os<<cc<<endl;
    ptrSrv->writeToR(os); os<<cc<<endl;
    ptrMSE->writeToR(os); os<<endl;
    os<<")";
}
        
