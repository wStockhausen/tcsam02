#include <limits>
#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelPopDyClasses.hpp"

using namespace tcsam;

////////////////////////////////////////////////////////////////////////////////
//PopDyInfo
////////////////////////////////////////////////////////////////////////////////
/**flag to print debug info*/
int PopDyInfo::debug = 0;
/**
 * Constructor.
 * 
 * @param npZBs - number of size bins
 */
PopDyInfo::PopDyInfo(int npZBs){
    nZBs = npZBs;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    R_z.allocate(1,nZBs);
    w_mz.allocate(1,nMSs,1,nZBs);
    Th_sz.allocate(1,nSCs,1,nZBs);
    M_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    T_szz.allocate(1,nSCs,1,nZBs,1,nZBs);
    
    np_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    S_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
}

PopDyInfo& PopDyInfo::operator=(const PopDyInfo& o) {
    nZBs = o.nZBs;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    R_z = o.R_z;
    w_mz = o.w_mz;
    Th_sz = o.Th_sz;
    M_msz = o.M_msz;
    T_szz = o.T_szz;
    
    np_msz = o.np_msz;
    S_msz  = o.S_msz;
    
    return *this;
}

/**
 * Calculate single-sex mature biomass, given population abundance.
 * 
 * @param n_msz - single-sex population abundance
 * 
 * @return mature biomass (1000's t)
 */
dvariable PopDyInfo::calcMatureBiomass(dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting double PopDyInfo::calcMMB(n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvariable mmb = 0.0; 
    for (int s=1;s<=nSCs;s++) mmb += n_msz(MATURE,s)*w_mz(MATURE);//dot product here
    if (debug) cout<<"finished dvariable PopDyInfo::calcMMB(n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return mmb;
}

/**
 * Calculate total single-sex biomass, given population abundance.
 * 
 * @param n_msz - single-sex population abundance
 * 
 * @return total biomass (1000's t)
 */
dvariable PopDyInfo::calcTotalBiomass(dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::calcBiomass()"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvariable bio = 0.0;
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            bio += w_mz(m)*n_msz(m,s);
        }//m
    }//s
    if (debug) cout<<"finished PopDyInfo::calcBiomass()"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return bio;
}

/**
 * Calculate single-sex survival probabilities, given natural mortality rates,
 * over a given period of time (dt)
 * 
 * @param dt - period of time (in years)
 * 
 * @return d3_array S_msz
 */
dvar3_array PopDyInfo::calcSurvival(double dt, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::calcSurvival(dt)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    S_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            S_msz(m,s) = exp(-M_msz(m,s)*dt); //survival over dt
        }//m
    }//s  
    if (debug) cout<<"finished PopDyInfo::calcSurvival(dt)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return S_msz;
}

/**
 * Apply natural mortality rates over a given period of time (dt) to
 * a single sex component of a population.
 * 
 * @param dt - period of time (in years)
 * 
 * @return d3_array - final population abundance
 */
dvar3_array PopDyInfo::applyNM(double dt, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::applyNM(dt,n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    np_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            np_msz(m,s) = elem_prod(exp(-M_msz(m,s)*dt),n_msz(m,s)); //survival over dt
        }//m
    }//s  
    if (debug) cout<<"finished PopDyInfo::applyNM(dt,n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return np_msz;
}

/** Apply molting/growth to population */
/**
 * 
 * @param n_msz
 * @return 
 */
dvar3_array PopDyInfo::applyMG(dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::applyMG(n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    np_msz.initialize();
    np_msz(IMMATURE,NEW_SHELL) = elem_prod(1.0-Th_sz(NEW_SHELL),T_szz(NEW_SHELL)*n_msz(IMMATURE,NEW_SHELL));
//    np_msz(IMMATURE,OLD_SHELL) = 0.0;
    np_msz(MATURE,NEW_SHELL)   = elem_prod(    Th_sz(NEW_SHELL),T_szz(NEW_SHELL)*n_msz(IMMATURE,NEW_SHELL));
    np_msz(MATURE,OLD_SHELL)   = n_msz(MATURE,NEW_SHELL)+n_msz(MATURE,OLD_SHELL);
    if (debug) cout<<"finished PopDyInfo::applyMG(n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return np_msz;
}

/**
 * Add single-sex recruitment to population abundance for one sex.
 * 
 * @param R - total single-sex recruitment
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 * @return single-sex abundance with recruitment
 */
dvar3_array PopDyInfo::addRecruitment(dvariable R, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::addRecruitment(R, n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    np_msz = 1.0*n_msz;
    np_msz(IMMATURE,NEW_SHELL) += R*R_z;
    if (debug) cout<<"finished PopDyInfo::addRecruitment(R, n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return np_msz;
}

/**
 * Write important quantities to output stream as an R list.
 * 
 * @param os - the output stream
 * @param ptrMC - pointer to the model configuration
 * @param name  - name for list in R
 * @param debug - flag to print debugging info
 * 
 * @return void
 */
void PopDyInfo::writeToR(ostream& os, ModelConfiguration* ptrMC, adstring name, int dbg){
    if (dbg) std::cout<<"Starting PopDyInfo::writeToR()"<<endl;
    adstring mDms  = ptrMC->dimMSsToR;//maturity
    adstring sDms  = ptrMC->dimSCsToR;//shell condition
    adstring zDms  = ptrMC->dimZBsToR;//size bin midpoints
    adstring zpDms = ptrMC->dimZPsToR;//size bin midpoints
    os<<name<<"=list(";
    os<<"R_z="; wts::writeToR(os,value(R_z),zDms);                os<<cc<<endl;
    os<<"w_mz="; wts::writeToR(os,w_mz,mDms,zDms);                os<<cc<<endl;
    os<<"M_msz="; wts::writeToR(os,value(M_msz),mDms,sDms,zDms);  os<<cc<<endl;
    os<<"Th_sz="; wts::writeToR(os,value(Th_sz),sDms,zDms);       os<<cc<<endl;
    os<<"T_szz="; wts::writeToR(os,value(T_szz),sDms,zDms,zpDms); os<<endl;
    os<<")";
    if (dbg) std::cout<<"Finished PopDyInfo::writeToR()"<<endl;
}
        
////////////////////////////////////////////////////////////////////////////////
//CatchInfo
////////////////////////////////////////////////////////////////////////////////
/** flag to print debug info */
int CatchInfo::debug = 0;
/**
 * Constructor.
 * 
 * @param npZBs - number of size bins
 * @param npFsh - number of fisheries
 */
CatchInfo::CatchInfo(int npZBs, int npFsh){
    nZBs = npZBs;
    nFsh = npFsh;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    maxF = 1.0;//default scale    
    
    hm_f.allocate(1,nFsh);
    cpF_fms.allocate(1,nFsh,1,nMSs,1,nSCs);
    cpF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    selF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    retF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    
    rmF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dmF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    
    cmN_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    cpN_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    rmN_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dmN_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    
    totFM.allocate(1,nZBs);
    S_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
}

CatchInfo& CatchInfo::operator =(const CatchInfo& o){
    nZBs = o.nZBs;
    nFsh = o.nFsh;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    maxF = o.maxF;//default scale
    
    hm_f      = o.hm_f;
    cpF_fms   = o.cpF_fms;
    cpF_fmsz  = o.cpF_fmsz;
    selF_fmsz = o.selF_fmsz;
    retF_fmsz = o.retF_fmsz;
    
    rmF_fmsz = o.rmN_fmsz;
    dmF_fmsz = o.dmN_fmsz;
    
    cmN_msz  = o.cmN_msz;
    cpN_fmsz = o.cpN_fmsz;
    rmN_fmsz = o.rmN_fmsz;
    dmN_fmsz = o.dmN_fmsz;
    
    totFM = o.totFM;
    S_msz = o.S_msz;
    
    return *this;
}

/**
 * Find max fishing mortality rate for target fishery (f=1)
 * 
 * Modifies: maxF
 * 
 * @return max fishing mortality rate
 */
dvariable CatchInfo::findMaxTargetCaptureRate(ostream& cout){
    maxF = wts::max(cpF_fmsz(1));
    return maxF;
}

/**
 * Calculate single-sex catch abundance (cm_msz, cp_fmsz, rm_fmsz, dm_fmsz) and 
 * post-fisheries abundance based on target fishing mortality rate 'dirF' and 
 * pre-fisheries population abundance n_msz.
 * 
 * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
 * 
 * Modifies: 
 *      rmF_fmsz - retained mortality rate, by fishery
 *      dmF_fmsz - discard mortality rate, by fishery
 *      cmN_msz  - total fishing mortality (abundance)
 *      cpN_fmsz - fishery captures, by fishery (abundance)
 *      rmN_fmsz - retained mortality, by fishery (abundance)
 *      dmN_fmsz - discard mortality, by fishery (abundance)
 * 
 * @param dirF - directed fishery fishing mortality rate
 * @param n_msz - pre-fisheries population size
 * 
 * @return np_msz - post-fisheries population size
 * 
 */
dvar3_array CatchInfo::applyFM(dvariable dirF, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting Catch_Calculator::calcCatch(double dirF, d3_array n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvariable ratF = 1.0;        //default target fishery (f=1) scaling ratio
    if ((dirF>=0.0)&&(maxF>0.0)) 
        ratF = dirF/maxF;        //target fishery (f=1) scaling ratio, if target fishery is open   (maxF>0)
    dvar3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();//number surviving fisheries
    rmF_fmsz.initialize();//retained catch mortality rates
    dmF_fmsz.initialize();//discards catch mortality rates
    cmN_msz.initialize();//total catch mortality (abundance)
    cpN_fmsz.initialize();//capture abundance, by fishery
    rmN_fmsz.initialize();//retained catch mortality (abundance)
    dmN_fmsz.initialize();//discards catch mortality (abundance)
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM.initialize();
            //directed fishery
            rmF_fmsz(1,m,s) = elem_prod(retF_fmsz(1,m,s),              ratF*cpF_fmsz(1,m,s));
            dmF_fmsz(1,m,s) = elem_prod(hm_f(1)*(1.0-retF_fmsz(1,m,s)),ratF*cpF_fmsz(1,m,s));
            totFM += rmF_fmsz(1,m,s)+dmF_fmsz(1,m,s);
            //bycatch fisheries
            for (int f=2;f<=nFsh;f++){
                rmF_fmsz(f,m,s) = elem_prod(retF_fmsz(f,m,s),              cpF_fmsz(f,m,s));
                dmF_fmsz(f,m,s) = elem_prod(hm_f(f)*(1.0-retF_fmsz(f,m,s)),cpF_fmsz(f,m,s));
                totFM += rmF_fmsz(f,m,s)+dmF_fmsz(f,m,s);
            }
            
            np_msz(m,s) = elem_prod(mfexp(-totFM),n_msz(m,s));//survival after all fisheries
            cmN_msz(m,s) = n_msz(m,s)-np_msz(m,s);             //total catch mortality, all fisheries
            
            //total capture abundance, directed fishery
            cpN_fmsz(1,m,s) = elem_prod(elem_div(ratF*cpF_fmsz(1,m,s),totFM),cmN_msz(m,s));
            //retained catch mortality (abundance), directed fishery
            rmN_fmsz(1,m,s) = elem_prod(elem_div(rmF_fmsz(1,m,s),totFM),cmN_msz(m,s));
            //discard catch mortality (abundance), directed fishery
            dmN_fmsz(1,m,s) = elem_prod(elem_div(dmF_fmsz(1,m,s),totFM),cmN_msz(m,s));
            //bycatch fisheries
            for (int f=2;f<=nFsh;f++){
                //total capture abundance
                cpN_fmsz(f,m,s) = elem_prod(elem_div(cpF_fmsz(f,m,s),totFM),cmN_msz(m,s));
                //retained catch mortality (abundance)
                rmN_fmsz(f,m,s) = elem_prod(elem_div(rmF_fmsz(f,m,s),totFM),cmN_msz(m,s));
                //discard catch mortality (abundance)
                dmN_fmsz(f,m,s) = elem_prod(elem_div(dmF_fmsz(f,m,s),totFM),cmN_msz(m,s));
            }//f
        }//m
    }//s
    if (debug) cout<<"finished Catch_Calculator::calcCatch(double dirF, d3_array n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return np_msz;
}

/**
 * Calculates probabilities of surviving fisheries, given directed 
 * fishing capture rate 'dirF'.
 * 
 * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
 * 
 * @param dirF - capture rate for target fishery (f=1)
 * 
 * @return d3_array of survival probabilities S_msz.
 */
dvar3_array CatchInfo::calcSurvival(dvariable dirF, ostream& cout){
    RETURN_ARRAYS_INCREMENT();
    dvariable ratF = 1.0;     //default target fishery scaling ratio
    if ((dirF>=0)&&(maxF>0)) 
        ratF = dirF/maxF;     //target fishery (f=1) scaling ratio
    S_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM.initialize();
            totFM += elem_prod(ratF*cpF_fmsz(1,m,s),
                              retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
            for (int f=2;f<=nFsh;f++)
                totFM += elem_prod(cpF_fmsz(f,m,s),
                                  retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
            S_msz(m,s) = exp(-totFM);                //survival of fisheries
        }//m
    }//s
    RETURN_ARRAYS_DECREMENT();
    return S_msz;
}

/**
 * Set single-sex size-specific fishery capture rates.
 * 
 * @param capF_fmsz - input capture rates
 */
void CatchInfo::setCaptureRates(dvar4_array& capFp_fmsz){
    cpF_fmsz = capFp_fmsz;
}

/**
 * Set single-sex fishery capture rates.
 * 
 * @param capF_fmsz - input capture rates
 */
void CatchInfo::setCaptureRates(dvar3_array& capFp_fms){
    cpF_fms = capFp_fms;
}

/**
 * Set single-sex selectivity functions.
 * 
 * @param retF_fmsz - input selectivity functions
 */
void CatchInfo::setSelectivityFcns(dvar4_array& selFp_fmsz){
    selF_fmsz = selFp_fmsz;
}

/**
 * Set single-sex retention functions.
 * 
 * @param retF_fmsz - input retention functions
 */
void CatchInfo::setRetentionFcns(dvar4_array& retFp_fmsz){
    retF_fmsz = retFp_fmsz;
}

/**
 * Set handling mortality rates, by fishery.
 * 
 * @param pHM_f - dvector of handling mortality rates
 */
void CatchInfo::setHandlingMortality(dvar_vector& pHM_f){
    hm_f = pHM_f;
}

/**
 * Write important quantities to output stream as an R list.
 * 
 * @param os - the output stream
 * @param ptrMC - pointer to the model configuration
 * @param name  - name for list in R
 * @param debug - flag to print debugging info
 * 
 * @return void
 */
void CatchInfo::writeToR(ostream& os, ModelConfiguration* ptrMC, adstring name, int debug){
    if (debug) std::cout<<"Starting CatchInfo::writeToR()"<<endl;
    adstring xDms = ptrMC->dimSXsToR;//sex
    adstring mDms = ptrMC->dimMSsToR;//maturity
    adstring sDms = ptrMC->dimSCsToR;//shell condition
    adstring zDms = ptrMC->dimZBsToR;//size bin midpoints
    adstring fDms = ptrMC->dimFshToR;//fisheries
    os<<name<<"=list(";
    os<<"hm_f="; wts::writeToR(os,value(hm_f),fDms); os<<cc<<endl;
    os<<"capF_fms="; wts::writeToR(os,value(cpF_fms),fDms,mDms,sDms); os<<cc<<endl;
    os<<"capF_fmsz="; wts::writeToR(os,wts::value(cpF_fmsz),fDms,mDms,sDms,zDms); os<<cc<<endl;
    os<<"selF_fmsz="; wts::writeToR(os,wts::value(selF_fmsz),fDms,mDms,sDms,zDms); os<<cc<<endl;
    os<<"retF_fmsz="; wts::writeToR(os,wts::value(retF_fmsz),fDms,mDms,sDms,zDms); os<<endl;
    os<<")";
    if (debug) std::cout<<"Finished CatchInfo::writeToR()"<<endl;
}
        
////////////////////////////////////////////////////////////////////////////////
//PopProjector
////////////////////////////////////////////////////////////////////////////////
/** flag to print debug info */
int PopProjector::debug = 0;
/**
 * Constructor
 * 
 * @param npZBs - number of size bins
 */
PopProjector::PopProjector(PopDyInfo* pPIp, CatchInfo* pCIp){
    pPI = pPIp;
    pCI = pCIp;
    
    nMSs = pPI->nMSs;
    nSCs = pPI->nSCs;
    nZBs = pPI->nZBs;
    nFsh = pCI->nFsh;  
    
    n1_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n2_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n3_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n4_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n5_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
}

/**
 * Project sex-specific component of population ahead one year, WITHOUT recruitment.
 * 
 * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
 * 
 * Also calculates:
 *      matBio - spawning biomass at mating time
 *      pCI elements:
 *          rmF_fmsz - retained mortality rate, by fishery
 *          dmF_fmsz - discard mortality rate, by fishery 
 *          cmN_msz - total fishing mortality          (abundance)
 *          cpN_fmsz - fishery captures, by fishery    (abundance)
 *          rmN_fmsz - retained mortality, by fishery  (abundance)
 *          dmN_fmsz - discard mortality, by fishery   (abundance)
 * 
 * @param n_msz - initial numbers-at-maturity state/shell condition/size
 * @param dirF - multiplier on fishing mortality rate in directed fishery
 * 
 * @return final numbers-at-maturity state/shell condition/size, WITHOUT recruitment
 */
dvar3_array PopProjector::project(dvariable dirF, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::project(dirF, n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    n1_msz.initialize();
    n2_msz.initialize();
    n3_msz.initialize();
    n4_msz.initialize();
    n5_msz.initialize();
    //if (debug){PopDyInfo::debug=1; CatchInfo::debug=1;}
    if (dtF<=dtM){ //fisheries occur BEFORE molting/growth/maturity 
        if (debug) cout<<"dtF(<=dtM) = "<<dtF<<endl;
        //apply natural mortality BEFORE fisheries
        n1_msz = pPI->applyNM(dtF, n_msz,cout);
        //if (debug) cout<<1<<endl;
        //apply fisheries
        n2_msz = pCI->applyFM(dirF, n1_msz,cout);
        //if (debug) cout<<2<<endl;
        //apply natural mortality after fisheries but before molting/growth
        if (dtF==dtM){
            n3_msz = n2_msz;
            //if (debug) cout<<3<<endl;
        } else {
            n3_msz = pPI->applyNM(dtM-dtF, n2_msz,cout);
            //if (debug) cout<<4<<endl;
        }        
        //apply molting/growth
        n4_msz = pPI->applyMG(n3_msz,cout);
        //if (debug) cout<<5<<endl;
        //apply natural mortality AFTER molting/growth
        n5_msz = pPI->applyNM(1.0-dtM, n4_msz,cout);
        //if (debug) cout<<6<<endl;
        
        //calculate mature biomass-at-mating from pre-molting/growth abundance
        matBio = pPI->calcMatureBiomass(n3_msz,cout);
        //if (debug) cout<<7<<endl;
    } else { //fisheries occur AFTER molting/growth/maturity 
        if (debug) cout<<"dtM<dtF"<<endl;
        //apply natural mortality BEFORE molting/growth
        n1_msz = pPI->applyNM(dtM, n_msz,cout);        
        //apply molting/growth
        n2_msz = pPI->applyMG(n1_msz,cout);
        //apply natural mortality after molting/growth but before fisheries
        if (dtF==dtM){
            n3_msz = n2_msz;
        } else {
            n3_msz = pPI->applyNM(dtF-dtM, n2_msz,cout);
        }
        //apply fisheries
        n4_msz = pCI->applyFM(dirF, n3_msz,cout);
        //apply natural mortality AFTER fisheries
        n5_msz = pPI->applyNM(1.0-dtF, n4_msz,cout);
        
        //calculate mature biomass-at-mating from pre-molting/growth abundance
        matBio = pPI->calcMatureBiomass(n1_msz,cout);
    }
    
    if (debug){
        cout<<"------n0_msz = "<<endl; wts::print(n_msz, cout,1);
        cout<<"------n1_msz = "<<endl; wts::print(n1_msz,cout,1);
        for (int f=1;f<=nFsh;f++){
            cout<<"----fishery = "<<f<<endl;
            cout<<"----hmF = "<<pCI->hm_f(f)<<endl;
            cout<<"------cpF_msz = "<<endl; wts::print(pCI->cpF_fmsz(f),cout,1);
            cout<<"------rmF_msz = "<<endl; wts::print(pCI->rmF_fmsz(f),cout,1);
            cout<<"------dmF_msz = "<<endl; wts::print(pCI->dmF_fmsz(f),cout,1);
//            cout<<"------cpN_msz = "<<endl; wts::print(pCI->cpN_fmsz(f),cout,1);
//            cout<<"------rmN_msz = "<<endl; wts::print(pCI->rmN_fmsz(f),cout,1);
//            cout<<"------dmN_msz = "<<endl; wts::print(pCI->dmN_fmsz(f),cout,1);
        }//f
        cout<<"------n2_msz = "<<endl; wts::print(n2_msz,cout,1);
        cout<<"------n3_msz = "<<endl; wts::print(n3_msz,cout,1);
        cout<<"------n4_msz = "<<endl; wts::print(n4_msz,cout,1);
        cout<<"------n5_msz = "<<endl; wts::print(n5_msz,cout,1);
    }
    if (debug){PopDyInfo::debug=0; CatchInfo::debug=0;}
    if (debug) cout<<"finished PopProjector::project(dirF, n_msz)"<<endl;   
    dvar3_array np_msz = 1.0*n5_msz;
    RETURN_ARRAYS_DECREMENT();
    return np_msz;
}

/**
 * Project unfished single-sex population abundance forward one year,
 * based on single-sex population abundance on July 1. Recruitment
 * is NOT added in.
 * 
 * Also calculates:
 *      matBio - mature biomass at mating time
 *      pCI elements:
 *          cm_msz - total fishing mortality         (abundance) [=0]
 *          cp_fmsz - fishery captures, by fishery   (abundance) [=0]
 *          rm_fmsz - retained mortality, by fishery (abundance) [=0]
 *          dm_fmsz - discard mortality, by fishery  (abundance) [=0]
 * 
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 * @return final sex-specific abundance, without recruitment
 */
dvar3_array PopProjector::projectUnFished(dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::projectUnFished(dirF, n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    n1_msz.initialize();
    n2_msz.initialize();
    n3_msz.initialize();
    pCI->cmN_msz.initialize(); //set to 0's
    pCI->cpN_fmsz.initialize();//set to 0's
    pCI->rmN_fmsz.initialize();//set to 0's
    pCI->dmN_fmsz.initialize();//set to 0's
    //d3_array n1_msz(1,nMSs,1,nSCs,1,nZBs);
    //apply natural mortality BEFORE molting/growth/maturity
    n1_msz = pPI->applyNM(dtM, n_msz,cout);
    //apply molting/growth
    n2_msz = pPI->applyMG(n1_msz,cout);
    //apply natural mortality AFTER molting/growth
    n3_msz = pPI->applyNM(1.0-dtM, n2_msz,cout);

    //calculate mature biomass-at-mating from pre-molting/growth abundance
    matBio = pPI->calcMatureBiomass(n1_msz,cout);

    if (debug) cout<<"finished PopProjector::projectUnFished(n_msz)"<<endl;   
    dvar3_array np_msz = 1.0*n3_msz;
    RETURN_ARRAYS_DECREMENT();
    return np_msz;
}

/**
 * Add single-sex recruitment to population abundance for one sex.
 * 
 * @param R - total single-sex recruitment
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 * @return final single-sex abundance with recruitment
 */
dvar3_array PopProjector::addRecruitment(dvariable R, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::addRecruitment(R, n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvar3_array np_msz = pPI->addRecruitment(R,n_msz,cout);
    if (debug) cout<<"finished PopProjector::addRecruitment(R, n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return np_msz;
}

/**
 * Calculate mature biomass-at-mating by projecting population forward to mating time.
 * 
 * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
 * 
 * @param n_msz - initial (July 1) numbers-at-maturity state/shell condition/size
 * 
 * @return MMB-at-mating
 */
dvariable PopProjector::projectMatureBiomassAtMating(dvariable dirF, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::projectMatureBiomassAtMating(dirF, n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"dtF = "<<dtF<<"; dtM = "<<dtM<<endl;
    n1_msz.initialize();
    n2_msz.initialize();
    n3_msz.initialize();
    if (dtF<=dtM){ //fisheries occur BEFORE molting/growth/maturity 
        if (debug) cout<<"dtF<=dtM"<<endl;
        //apply natural mortality BEFORE fisheries
        n1_msz = pPI->applyNM(dtF,n_msz, cout);
        if (debug) {cout<<"n1_msz ="<<endl; wts::print(n1_msz,cout,1);}
        //apply fisheries
        n2_msz = pCI->applyFM(dirF, n1_msz, cout);
        if (debug) {cout<<"n2_msz ="<<endl; wts::print(n2_msz,cout,1);}
        //apply natural mortality after fisheries but before molting/growth
        if (dtF==dtM){
            if (debug) cout<<"dtF=dtM"<<endl;
            n3_msz = n2_msz;
        } else {
            if (debug) cout<<"dtF<dtM"<<endl;
            n3_msz = pPI->applyNM(dtM-dtF,n2_msz,cout);
        }
        if (debug) {cout<<"n3_msz ="<<endl; wts::print(n3_msz,cout,1);}
        
        //calculate mature biomass at mating
        matBio = pPI->calcMatureBiomass(n3_msz,cout);
    } else { //fisheries occur AFTER molting/growth/maturity 
        if (debug) cout<<"dtF>dtM"<<endl;
        //apply natural mortality BEFORE molting/growth
        n1_msz = pPI->applyNM(dtM,n_msz,cout);
        if (debug) {cout<<"n1_msz ="<<endl; wts::print(n1_msz,cout,1);}
        
        //calculate mature biomass at mating
        matBio = pPI->calcMatureBiomass(n1_msz,cout);
    }
    
    if (debug) {
        cout<<"matBio = "<<matBio<<endl;
        cout<<"finished PopProjector::projectMatureBiomassAtMating(dirF, n_msz)"<<endl;   
    }
    RETURN_ARRAYS_DECREMENT();
    return matBio;
}
////////////////////////////////////////////////////////////////////////////////
//MultiYearPopProjector
////////////////////////////////////////////////////////////////////////////////
/** flag to print debug info */
int MultiYearPopProjector::debug = 0;
/**
 * Project multiple years at constant recruitment and directed F.
 * 
 * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
 * 
 * @param n - number of years to project
 * @param R - (constant) single-sex recruitment
 * @param dirF - (constant) directed F
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 */
void MultiYearPopProjector::project(int n, dvariable R, dvariable dirF, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"Starting MultiYearPopProjector::project(n,R,n_msz)"<<endl;
    if (debug) cout<<"nFsh = "<<pPP->nFsh<<endl;
    n_ymsz.allocate(0,n,1,pPP->nMSs,1,pPP->nSCs,1,pPP->nZBs);
    matBio_y.allocate(1,n);
    totCM_y.allocate(1,n);
    cp_yf.allocate(0,n,1,pPP->nFsh);
    rm_yf.allocate(0,n,1,pPP->nFsh);
    dm_yf.allocate(0,n,1,pPP->nFsh);
    n_ymsz.initialize();
    matBio_y.initialize();
    totCM_y.initialize();
    cp_yf.initialize();
    rm_yf.initialize();
    dm_yf.initialize();
    n_ymsz(0) = n_msz;
    if (debug){PopProjector::debug=1;PopDyInfo::debug=1;}
    for (int y=1;y<=n;y++){
        dvar3_array np_msz = pPP->project(dirF,n_ymsz(y-1),cout);
        matBio_y(y) = pPP->matBio;
        for (int f=1;f<=pPP->nFsh;f++) {
            cp_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->cpN_fmsz(f),cout);
            rm_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->rmN_fmsz(f),cout);
            dm_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->dmN_fmsz(f),cout);
        }
        totCM_y(y) = pPP->pPI->calcTotalBiomass(pPP->pCI->cmN_msz,cout);
        n_ymsz(y)  = pPP->addRecruitment(R,np_msz,cout);
    }
    if (debug){PopProjector::debug=0;PopDyInfo::debug=0;}
    if (debug){
        cout<<"matBio = "<<endl<<matBio_y<<endl;
        cout<<"totCM_y = "<<endl<<totCM_y<<endl;
        cout<<"cp_yf = "<<endl<<cp_yf<<endl;
        cout<<"rm_yf = "<<endl<<rm_yf<<endl;
        cout<<"dm_yf = "<<endl<<dm_yf<<endl;
    }
    if (debug) cout<<"Finished MultiYearPopProjector::project(n,R,n_msz)"<<endl;
}
/**
 * Project multiple years at constant recruitment and directed F.
 * 
 * @param n - number of years to project
 * @param R - (constant) single-sex recruimtent
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 */
void MultiYearPopProjector::projectUnFished(int n, dvariable R, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"Starting MultiYearPopProjector::projectUnFished(n,R,n_msz)"<<endl;
    n_ymsz.allocate(0,n,1,pPP->nMSs,1,pPP->nSCs,1,pPP->nZBs);
    matBio_y.allocate(1,n);
    rm_yf.allocate(0,n,1,pPP->nFsh);
    dm_yf.allocate(0,n,1,pPP->nFsh);
    totCM_y.allocate(1,n);
    n_ymsz.initialize();
    matBio_y.initialize();
    rm_yf.initialize();
    dm_yf.initialize();
    totCM_y.initialize();
    n_ymsz(0) = n_msz;
    for (int y=1;y<=n;y++){
        dvar3_array np_msz = pPP->projectUnFished(n_ymsz(y-1),cout);
        matBio_y(y) = pPP->matBio;
        n_ymsz(y) = pPP->addRecruitment(R,np_msz,cout);
    }
    if (debug){
        cout<<"matBio = "<<endl<<matBio_y<<endl;
        cout<<"Finished MultiYearPopProjector::projectUnFished(n,R,n_msz)"<<endl;
    }
}

