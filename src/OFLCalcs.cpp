#include <limits>
#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "OFLCalcs.hpp"

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
 * @params dtMp - time at which mating occurs (as fraction of year)
 */
PopDyInfo::PopDyInfo(int npZBs, double dtMp){
    dtM = dtMp;
    nZBs = npZBs;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    R_z.allocate(1,nZBs);
    w_mz.allocate(1,nMSs,1,nZBs);
    Th_sz.allocate(1,nSCs,1,nZBs);
    M_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    T_szz.allocate(1,nSCs,1,nZBs,1,nZBs);
}

/**
 * Calculate single-sex mature biomass, given population abundance.
 * 
 * @param n_msz - single-sex population abundance
 * 
 * @return mature biomass (1000's t)
 */
double PopDyInfo::calcMatureBiomass(d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting double PopDyInfo::calcMMB(n_msz)"<<endl;
    double mmb = 0.0; 
    for (int s=1;s<=nSCs;s++) mmb += n_msz(MATURE,s)*w_mz(MATURE);//dot product here
    if (debug) cout<<"finished double PopDyInfo::calcMMB(n_msz)"<<endl;
    return mmb;
}

/**
 * Calculate total single-sex biomass, given population abundance.
 * 
 * @param n_msz - single-sex population abundance
 * 
 * @return total biomass (1000's t)
 */
double PopDyInfo::calcTotalBiomass(d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::calcBiomass()"<<endl;
    double bio = 0.0;
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            bio += w_mz(m)*n_msz(m,s);
        }//m
    }//s
    if (debug) cout<<"finished PopDyInfo::calcBiomass()"<<endl;
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
d3_array PopDyInfo::calcSurvival(double dt, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::calcSurvival(dt)"<<endl;
    d3_array S_msz(1,nMSs,1,nSCs,1,nZBs);
    S_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            S_msz(m,s) = exp(-M_msz(m,s)*dt); //survival over dt
        }//m
    }//s  
    if (debug) cout<<"finished PopDyInfo::calcSurvival(dt)"<<endl;
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
d3_array PopDyInfo::applyNM(double dt, d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::applyNM(dt,n_msz)"<<endl;
    d3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            np_msz(m,s) = elem_prod(exp(-M_msz(m,s)*dt),n_msz(m,s)); //survival over dt
        }//m
    }//s  
    if (debug) cout<<"finished PopDyInfo::applyNM(dt,n_msz)"<<endl;
    return np_msz;
}

/** Apply molting/growth to population */
/**
 * 
 * @param n_msz
 * @return 
 */
d3_array PopDyInfo::applyMG(d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::applyMG(n_msz)"<<endl;
    d3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();
    np_msz(IMMATURE,NEW_SHELL) = elem_prod(1.0-Th_sz(NEW_SHELL),T_szz(NEW_SHELL)*n_msz(IMMATURE,NEW_SHELL));
//    np_msz(IMMATURE,OLD_SHELL) = 0.0;
    np_msz(MATURE,NEW_SHELL)   = elem_prod(    Th_sz(NEW_SHELL),T_szz(NEW_SHELL)*n_msz(IMMATURE,NEW_SHELL));
    np_msz(MATURE,OLD_SHELL)   = n_msz(MATURE,NEW_SHELL)+n_msz(MATURE,OLD_SHELL);
    if (debug) cout<<"finished PopDyInfo::applyMG(n_msz)"<<endl;
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
d3_array PopDyInfo::addRecruitment(double R, d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopDyInfo::addRecruitment(R, n_msz)"<<endl;
    d3_array np_msz = 1.0*n_msz;
    np_msz(IMMATURE,NEW_SHELL) += R*R_z;
    if (debug) cout<<"finished PopDyInfo::addRecruitment(R, n_msz)"<<endl;
    return np_msz;
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
 * @param dtFp - time at which fisheries occur (as fraction of year)
 */
CatchInfo::CatchInfo(int npZBs, int npFsh, double dtFp){
    dtF  = dtFp;
    nZBs = npZBs;
    nFsh = npFsh;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    maxF = 1.0;//default scale
    cm_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    cp_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    rm_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    dm_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
}

/**
 * Find max fishing mortality rate for target fishery (f=1)
 * 
 * Modifies: maxF
 * 
 * @return max fishing mortality rate
 */
double CatchInfo::findMaxTargetCaptureRate(ostream& cout){
    maxF = wts::max(capF_fmsz(1));
    return maxF;
}

/**
 * Calculate single-sex catch abundance (cm_msz, cp_fmsz, rm_fmsz, dm_fmsz) and 
 * post-fisheries abundance based on target fishing mortality rate 'dirF' and 
 * pre-fisheries population abundance n_msz.
 * 
 * Modifies: 
 *      cm_msz  - total fishing mortality (abundance)
 *      cp_fmsz - fishery captures, by fishery (abundance)
 *      rm_fmsz - retained mortality, by fishery (abundance)
 *      dm_fmsz - discard mortality, by fishery (abundance)
 * 
 * @param dirF - directed fishery fishing mortality rate
 * @param n_msz - pre-fisheries population size
 * 
 * @return np_msz - post-fisheries population size
 * 
 */
d3_array CatchInfo::applyFM(double dirF, d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting Catch_Calculator::calcCatch(double dirF, d3_array n_msz)"<<endl;
    double ratF = dirF/maxF;//target fishery (f=1) scaling ratio
    dvector totFM(1,nZBs);
    d3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();
    cm_msz.initialize();
    cp_fmsz.initialize();
    rm_fmsz.initialize();
    dm_fmsz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM.initialize();
            //directed fishery
            totFM += elem_prod(ratF*capF_fmsz(1,m,s),
                               retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
            //bycatch fisheries
            for (int f=2;f<=nFsh;f++)
                totFM += elem_prod(capF_fmsz(f,m,s),
                                   retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
            np_msz(m,s) = elem_prod(exp(-totFM),n_msz(m,s));//survival after all fisheries
            cm_msz(m,s) = n_msz(m,s)-np_msz(m,s);           //total catch mortality, all fisheries
            
            //total capture abundance, directed fishery
            cp_fmsz(1,m,s) = elem_prod(elem_div(ratF*capF_fmsz(1,m,s),totFM),cm_msz(m,s));
            //retained catch mortality (abundance), directed fishery
            rm_fmsz(1,m,s) =        elem_prod(     retF_fmsz(1,m,s),cp_fmsz(1,m,s));
            //discard catch mortality (abundance), directed fishery
            dm_fmsz(1,m,s) = hm_f(1)*elem_prod(1.0-retF_fmsz(1,m,s),cp_fmsz(1,m,s));
            //bycatch fisheries
            for (int f=2;f<=nFsh;f++){
                //total capture abundance
                cp_fmsz(f,m,s) = elem_prod(elem_div(capF_fmsz(f,m,s),totFM),cm_msz(m,s));
                //retained catch mortality (abundance)
                rm_fmsz(f,m,s) =        elem_prod(     retF_fmsz(f,m,s),cp_fmsz(f,m,s));
                //discard catch mortality (abundance)
                dm_fmsz(f,m,s) = hm_f(f)*elem_prod(1.0-retF_fmsz(f,m,s),cp_fmsz(f,m,s));
            }//f
        }//m
    }//s
    if (debug) cout<<"finished Catch_Calculator::calcCatch(double dirF, d3_array n_msz)"<<endl;
    return np_msz;
}

/**
 * Calculates probabilities of surviving fisheries, given directed 
 * fishing capture rate 'dirF'.
 * 
 * @param dirF - capture rate for target fishery (f=1)
 * 
 * @return d3_array of survival probabilities S_msz.
 */
d3_array CatchInfo::calcSurvival(double dirF, ostream& cout){
    double ratF = dirF/maxF;//target fishery (f=1) scaling ratio
    dvector totFM(1,nZBs);
    d3_array S_msz(1,nMSs,1,nSCs,1,nZBs);
    S_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM.initialize();
            totFM += elem_prod(ratF*capF_fmsz(1,m,s),
                              retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
            for (int f=2;f<=nFsh;f++)
                totFM += elem_prod(capF_fmsz(f,m,s),
                                  retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
            S_msz(m,s) = exp(-totFM);                //survival of fisheries
        }//m
    }//s
    return S_msz;
}
/**
 * Set single-sex fishery capture rates.
 * 
 * @param capF_fmsz - input capture rates
 */
void CatchInfo::setCaptureRates(d4_array& capFp_fmsz){
    capF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    capF_fmsz = capFp_fmsz;
}

/**
 * Set single-sex retention functions.
 * 
 * @param retF_fmsz - input retention functions
 */
void CatchInfo::setRetentionFcns(d4_array& retFp_fmsz){
    retF_fmsz.allocate(1,nFsh,1,nMSs,1,nSCs,1,nZBs);
    retF_fmsz = retFp_fmsz;
}

/**
 * Set handling mortality rates, by fishery.
 * 
 * @param pHM_f - dvector of handling mortality rates
 */
void CatchInfo::setHandlingMortality(dvector& pHM_f){
    hm_f.allocate(1,nFsh);
    hm_f = pHM_f;
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
    
    dtM = pPI->dtM;
    dtF = pCI->dtF;
}

/**
 * Project sex-specific component of population ahead one year, WITHOUT recruitment.
 * Also calculates:
 *      matBio - spawning biomass at mating time
 *      pCI elements:
 *          cm_msz - total fishing mortality          (abundance)
 *          cp_fmsz - fishery captures, by fishery    (abundance)
 *          rm_fmsz - retained mortality, by fishery  (abundance)
 *          dm_fmsz - discard mortality, by fishery   (abundance)
 * 
 * @param n_msz - initial numbers-at-maturity state/shell condition/size
 * @param dirF - multiplier on fishing mortality rate in directed fishery
 * 
 * @return final numbers-at-maturity state/shell condition/size, WITHOUT recruitment
 */
d3_array PopProjector::project(double dirF, d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::project(dirF, n_msz)"<<endl;
    d3_array n1_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n2_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n3_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n4_msz(1,nMSs,1,nSCs,1,nZBs);
    d3_array n5_msz(1,nMSs,1,nSCs,1,nZBs);
    if (debug){PopDyInfo::debug=1; CatchInfo::debug=1;}
    if (dtF<=dtM){ //fisheries occur BEFORE molting/growth/maturity 
        if (debug) cout<<"dtF<=dtM"<<endl;
        //apply natural mortality BEFORE fisheries
        n1_msz = pPI->applyNM(dtF, n_msz,cout);
        if (debug) cout<<1<<endl;
        //apply fisheries
        n2_msz = pCI->applyFM(dirF, n1_msz,cout);
        if (debug) cout<<2<<endl;
        //apply natural mortality after fisheries but before molting/growth
        if (dtF==dtM){
            n3_msz = n2_msz;
            if (debug) cout<<3<<endl;
        } else {
            n3_msz = pPI->applyNM(dtM-dtF, n2_msz,cout);
            if (debug) cout<<4<<endl;
        }        
        //apply molting/growth
        n4_msz = pPI->applyMG(n3_msz,cout);
        if (debug) cout<<5<<endl;
        //apply natural mortality AFTER molting/growth
        n5_msz = pPI->applyNM(1.0-dtM, n4_msz,cout);
        if (debug) cout<<6<<endl;
        
        //calculate mature biomass-at-mating from pre-molting/growth abundance
        matBio = pPI->calcMatureBiomass(n3_msz,cout);
        if (debug) cout<<7<<endl;
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
    
    if (debug){PopDyInfo::debug=0; CatchInfo::debug=0;}
    if (debug) cout<<"finished PopProjector::project(dirF, n_msz)"<<endl;   
    return n5_msz;
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
d3_array PopProjector::projectUnFished(d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::projectUnFished(dirF, n_msz)"<<endl;
    pCI->cm_msz.initialize(); //set to 0's
    pCI->cp_fmsz.initialize();//set to 0's
    pCI->rm_fmsz.initialize();//set to 0's
    pCI->dm_fmsz.initialize();//set to 0's
    //d3_array n1_msz(1,nMSs,1,nSCs,1,nZBs);
    //apply natural mortality BEFORE molting/growth/maturity
    d3_array n1_msz = pPI->applyNM(dtM, n_msz,cout);
    //apply molting/growth
    d3_array n2_msz = pPI->applyMG(n1_msz,cout);
    //apply natural mortality AFTER molting/growth
    d3_array n3_msz = pPI->applyNM(1.0-dtM, n2_msz,cout);

    //calculate mature biomass-at-mating from pre-molting/growth abundance
    matBio = pPI->calcMatureBiomass(n1_msz,cout);

    if (debug) cout<<"finished PopProjector::projectUnFished(n_msz)"<<endl;   
    return n3_msz;
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
d3_array PopProjector::addRecruitment(double R, d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::addRecruitment(R, n_msz)"<<endl;
    d3_array n1_msz = pPI->addRecruitment(R,n_msz,cout);
    if (debug) cout<<"finished PopProjector::addRecruitment(R, n_msz)"<<endl;
    return n1_msz;
}

/**
 * Calculate mature biomass-at-mating by projecting population forward to mating time.
 * 
 * @param n_msz - initial (July 1) numbers-at-maturity state/shell condition/size
 * 
 * @return MMB-at-mating
 */
double PopProjector::projectMatureBiomassAtMating(double dirF, d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting PopProjector::projectMatureBiomassAtMating(dirF, n_msz)"<<endl;
    if (debug) cout<<"dtF = "<<dtF<<"; dtM = "<<dtM<<endl;
    if (dtF<=dtM){ //fisheries occur BEFORE molting/growth/maturity 
        if (debug) cout<<"dtF<=dtM"<<endl;
        d3_array n1_msz(1,nMSs,1,nSCs,1,nZBs);
        d3_array n2_msz(1,nMSs,1,nSCs,1,nZBs);
        d3_array n3_msz(1,nMSs,1,nSCs,1,nZBs);
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
        d3_array n1_msz(1,nMSs,1,nSCs,1,nZBs);
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
 * @param n - number of years to project
 * @param R - (constant) single-sex recruitment
 * @param dirF - (constant) directed F
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 */
void MultiYearPopProjector::project(int n, double R, double dirF, d3_array& n_msz, ostream& cout){
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
        d3_array np_msz = pPP->project(dirF,n_ymsz(y-1),cout);
        matBio_y(y) = pPP->matBio;
        for (int f=1;f<=pPP->nFsh;f++) {
            cp_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->cp_fmsz(f),cout);
            rm_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->rm_fmsz(f),cout);
            dm_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->dm_fmsz(f),cout);
        }
        totCM_y(y) = pPP->pPI->calcTotalBiomass(pPP->pCI->cm_msz,cout);
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
void MultiYearPopProjector::projectUnFished(int n, double R, d3_array& n_msz, ostream& cout){
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
        d3_array np_msz = pPP->projectUnFished(n_ymsz(y-1),cout);
        matBio_y(y) = pPP->matBio;
        n_ymsz(y) = pPP->addRecruitment(R,np_msz,cout);
    }
    if (debug){
        cout<<"matBio = "<<endl<<matBio_y<<endl;
        cout<<"Finished MultiYearPopProjector::projectUnFished(n,R,n_msz)"<<endl;
    }
}
////////////////////////////////////////////////////////////////////////////////
//Equilibrium_Calculator
////////////////////////////////////////////////////////////////////////////////
/** flag to print debug info */
int Equilibrium_Calculator::debug = 0;
/**
 * Calculate equilibrium (longterm) population abundance on July 1.
 * 
 * NOTE: these equations assume prM2M(z) is based on post-molt size, not pre-molt size.
 * 
 * @param R_z - longterm recruitment at size
 * @param S1_msz - survival prior to molting/growth
 * @param Th_sz - pr(molt-to-maturity|post-size) for previously-immature crab after growth
 * @param T_szz - growth transition matrix for molting crab
 * @param S2_msz - survival following molting/growth
 * @param cout - output stream for debug info
 * 
 * @return equlibrium (longterm) population abundance on July 1 (as d3_array)
 */
d3_array Equilibrium_Calculator::calcEqNatZ(dvector& R_z, d3_array& S1_msz, 
                                            dmatrix& Th_sz, d3_array& T_szz, 
                                            d3_array& S2_msz, ostream& cout){
    if (debug) cout<<"starting d3_array Equilibrium_Calculator::calcEqNatZ(...)"<<endl;

    int nMSs = pPP->pPI->nMSs;
    int nSCs = pPP->pPI->nSCs;
    int nZBs = pPP->pPI->nZBs;
    
    //the equilibrium solution
    d3_array n_msz(1,nMSs,1,nSCs,1,nZBs);
    
    //create an identity matrix
    dmatrix I = identity_matrix(1,nZBs);
    
    //--calc the state transition matrices
    int i = IMMATURE; 
    int m = MATURE;
    int n = NEW_SHELL;
    int o = OLD_SHELL;
    //immature new shell crab
    dmatrix S2_in = wts::diag(S2_msz(i,n)); //pr(survival|size) for immature new shell crab after molting occurs
    dmatrix Tr_in = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab pre-terminal molt
    dmatrix Th_in = wts::diag(Th_sz(n));    //pr(molt to maturity|post-molt size,new shell, molting)
    dmatrix Ph_in = identity_matrix(1,nZBs);//pr(molt|pre-molt size, new shell) [assumed that all immature crab molt]
    dmatrix S1_in = wts::diag(S1_msz(i,n)); //pr(survival|size) for immature new shell crab before molting occurs
    //immature old shell crab [shouldn't be any of these]
    dmatrix S2_io = wts::diag(S2_msz(i,o)); //pr(survival|size) for immature old shell crab after molting occurs (but they didn't molt)
    dmatrix Tr_io = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab pre-terminal molt
    dmatrix Th_io = wts::diag(Th_sz(o));    //pr(molt to maturity|post-molt size,old shell, molting)
    dmatrix Ph_io = identity_matrix(1,nZBs);//pr(molt|pre-molt size, old shell) [assumed all immature crab molt]
    dmatrix S1_io = wts::diag(S1_msz(i,o)); //pr(survival|size) for immature old shell crab before molting occurs
    //mature new shell crab
    dmatrix Tr_mn = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab undergoing terminal molt (same as non-terminal molt)
    dmatrix Tr_mo = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab undergoing terminal molt (same as non-terminal molt)
    dmatrix S2_mn = wts::diag(S2_msz(m,n)); //pr(survival|size) for mature new shell crab after molting occurs
    dmatrix S1_mn = wts::diag(S1_msz(m,n)); //pr(survival|size) for mature new shell crab before molting occurs (but they won't molt)
    //mature old shell crab
    dmatrix S2_mo = wts::diag(S2_msz(m,o)); //pr(survival|size) for mature old shell crab after molting occurs (but they didn't molt)
    dmatrix S1_mo = wts::diag(S1_msz(m,o)); //pr(survival|size) for mature old shell crab before molting occurs (but they won't molt)
    
    //full state transition matrices
    dmatrix lA = S2_in * (I-Th_in) * Tr_in * Ph_in * S1_in;//imm, new -> imm, new
    dmatrix lB = S2_in * (I-Th_io) * Tr_io * Ph_io * S1_io;//imm, old -> imm, new
    dmatrix lC = S2_io * (I-Ph_in) * S1_in;                //imm, new -> imm, old
    dmatrix lD = S2_io * (I-Ph_io) * S1_io;                //imm, old -> imm, old
    dmatrix lE = S2_mn * Th_in * Tr_mn * Ph_in * S1_in;    //imm, new -> mat, new (terminal molt)
    dmatrix lF = S2_mn * Th_io * Tr_mo * Ph_io * S1_io;    //imm, old -> mat, new (terminal molt)
    dmatrix lG = S2_mo * S1_mn;                            //mat, new -> mat, old
    dmatrix lH = S2_mo * S1_mo;                            //mat, old -> mat, old
    //--done calculating transition matrices
    
    //calculate inverses of matrix quantities
    dmatrix iM1 = inv(I - lD);
    dmatrix iM2 = inv(I - lA - lB * iM1 * lC);
    dmatrix iM3 = inv(I - lH);
    
    //the equilibrium solution is
    n_msz.initialize();
    n_msz(i,n) = iM2 * R_z;                         //immature, new shell
    n_msz(i,o) = iM1 * lC * n_msz(i,n);             //immature, old shell
    n_msz(m,n) = lE * n_msz(i,n) + lF * n_msz(i,o); //  mature, new shell
    n_msz(m,o) = iM3 * lG * n_msz(m,n);             //  mature, old shell
        
    if (debug) cout<<"finished d3_array Equilibrium_Calculator::calcEqNatZ(...)"<<endl;
    return n_msz;
}

/**
 * Calculate equilibrium unfished abundance on July 1 when 
 * longterm (average) male recruitment is R.
 * 
 * @param R - longterm (average) male recruitment
 * @param cout - output stream for debug info
 * 
 * @return equilibrium (male) population abundance on July 1
 */
d3_array Equilibrium_Calculator::calcEqNatZF0(double R, ostream& cout){
    if (debug) {
        cout<<"starting d3_array Equilibrium_Calculator::calcEqNatZF0(double R)"<<endl;
        cout<<"R = "<<R<<endl;
    }

    d3_array S1_msz = pPP->pPI->calcSurvival(pPP->dtM,cout);
    d3_array S2_msz = pPP->pPI->calcSurvival(1.0-pPP->dtM,cout);
    
    dvector R_zp = R*pPP->pPI->R_z;
    
    d3_array n_msz = calcEqNatZ(R_zp, S1_msz, pPP->pPI->Th_sz, pPP->pPI->T_szz, S2_msz,cout);
    
    if (debug) {
        cout<<"n_msz = "<<endl; wts::print(n_msz,cout,1);
        cout<<"finished d3_array Equilibrium_Calculator::calcEqNatZF0(double R)"<<endl;
    }
    return n_msz;
}

/**
 * Calculate equilibrium abundance on July 1 when dirF is the 
 * directed fishery capture rate and the longterm (average) 
 * male recruitment is R.
 * 
 * @param R - longterm (average) male recruitment
 * @param dirF - directed fishery capture rate
 * @param cout - output stream for debug info
 * 
 * @return equilibrium (male) population abundance on July 1
 */
d3_array Equilibrium_Calculator::calcEqNatZFM(double R, double dirF, ostream& cout){
    if (debug) {
        cout<<"starting d3_array Equilibrium_Calculator::calcEqNatZFM(double R)"<<endl;
        cout<<"R = "<<R<<". dirF = "<<dirF<<endl;
    }
    int nMSs = pPP->pPI->nMSs;
    int nSCs = pPP->pPI->nSCs;
    int nZBs = pPP->pPI->nZBs;
    
    double dtM = pPP->dtM;
    double dtF = pPP->dtF;
    
    d3_array S1_msz(1,nMSs,1,nSCs,1,nZBs);//survival until molting/mating
    d3_array S2_msz(1,nMSs,1,nSCs,1,nZBs);//survival after molting/mating
    d3_array n_msz(1,nMSs,1,nSCs,1,nZBs); //equilibrium size distribution on July 1
    
    if (dtF<=dtM){
        //fisheries occur BEFORE molting/growth/maturity 
        d3_array S1a_msz = pPP->pPI->calcSurvival(dtF,cout);     //survival prior to fisheries
        d3_array S1F_msz = pPP->pCI->calcSurvival(dirF,cout);    //survival of fisheries
        d3_array S1b_msz = pPP->pPI->calcSurvival(dtM-dtF,cout); //survival after fisheries, before mating
        S1_msz = elem_prod(S1b_msz,elem_prod(S1F_msz,S1a_msz));//total survival before mating
        S2_msz = pPP->pPI->calcSurvival(1.0-dtM,cout);//survival after mating/molting/growth
    } else {
        //fisheries occur AFTER molting/growth/maturity 
        S1_msz = pPP->pPI->calcSurvival(dtM,cout);                   //survival before mating/molting/growth
        d3_array S2a_msz = pPP->pPI->calcSurvival(dtF-dtM,cout);     //survival afterMGM, before fisheries
        d3_array S2F_msz = pPP->pCI->calcSurvival(dirF,cout);        //survival of fisheries
        d3_array S2b_msz = pPP->pPI->calcSurvival(1.0-dtF,cout);     //survival after fisheries, to year end
        S2_msz = elem_prod(S2b_msz,elem_prod(S2F_msz,S2a_msz)); //total survival after MGM to year end
    }
    dvector R_zp = R*pPP->pPI->R_z;
    n_msz = calcEqNatZ(R_zp, S1_msz, pPP->pPI->Th_sz, pPP->pPI->T_szz, S2_msz, cout);
    
    if (debug) {
        cout<<"n_msz = "<<endl; wts::print(n_msz,cout,1);
        cout<<"finished d3_array Equilibrium_Calculator::calcEqNatZFM(double R)"<<endl;
    }
    return n_msz; 
}
/**
 * Calculate single-sex equilibrium mature biomass-at-mating
 * for an unfished stock.
 * 
 * @param R - longterm (average) recruitment (male-only, millions)
 * @param cout - output stream for debug info
 * 
 * @return mature biomass-at-mating for unfished population (1000's t)
 */
double Equilibrium_Calculator::calcEqMatureBiomassAtMatingF0(double R, ostream& cout){
    if (debug) {
        cout<<"starting double Equilibrium_Calculator::calcEqMatureBiomassAtMatingF0(double R)"<<endl;
        cout<<"R = "<<R<<endl;
    }
    //calculate July 1 unfished size distribution
    d3_array n_msz = calcEqNatZF0(R,cout);
    
    //advance to mating in unfished population
    d3_array nm_msz = pPP->pPI->applyNM(pPP->dtM,n_msz,cout);
    
    //calculate mature biomass at time of mating
    double eqMB = pPP->pPI->calcMatureBiomass(nm_msz,cout);
    if (debug) {
        cout<<"eqMB = "<<eqMB<<endl;
        cout<<"finished double Equilibrium_Calculator::calcEqMatureBiomassAtMatingF0(double R)"<<endl;
    }
    return eqMB;
}

/**
 * Calculate equilibrium single-sex mature biomass-at-mating, 
 * when dirF is the directed fishery capture rate
 * and longterm (average) male recruitment is R.
 * 
 * @param R - longterm (average) male recruitment
 * @param dirF - directed fishery capture rate
 * @param cout - output stream for debug info
 * 
 * @return equilibrium spawning biomass-at-mating 
 */
double Equilibrium_Calculator::calcEqMatureBiomassAtMatingFM(double R, double dirF, ostream& cout){
    if (debug) {
        cout<<"starting Equilibrium_Calculator::calcEqMatureBiomassAtMatingFM(double R, double dirF)"<<endl;
        cout<<"R = "<<R<<". dirF = "<<dirF<<endl;
    }
    //calculate equilibrium size distribution on July 1
    d3_array n_msz = calcEqNatZFM(R,dirF,cout);
    
    //equilibrium MMB
    double eqMB = 0;//dummy value    
    //advance to population to time of mating
    if (pPP->dtF<=pPP->dtM){ //fisheries occur BEFORE molting/growth/maturity 
        d3_array n3_msz(1,pPP->nMSs,1,pPP->nSCs,1,pPP->nZBs);
        if (debug) {cout<<"dtF<=dtM"<<endl;}
        //apply natural mortality BEFORE fisheries
        d3_array n1_msz = pPP->pPI->applyNM(pPP->dtF,n_msz,cout);
        if (debug) {cout<<"n1_msz ="<<endl; wts::print(n1_msz,cout,1);}
        //apply fisheries
        d3_array n2_msz = pPP->pCI->applyFM(dirF, n1_msz, cout);
        if (debug) {cout<<"n2_msz ="<<endl; wts::print(n2_msz,cout,1);}
        //apply natural mortality after fisheries but before molting/growth
        if (pPP->dtF==pPP->dtM){
            if (debug) cout<<"dtF=dtM"<<endl;
            n3_msz = n2_msz;
        } else {
            if (debug) cout<<"dtF<dtM"<<endl;
            n3_msz = pPP->pPI->applyNM(pPP->dtM-pPP->dtF,n2_msz,cout);
        }
        if (debug) {cout<<"n3_msz ="<<endl; wts::print(n3_msz,cout,1);}
        
        //calculate mature biomass at mating
        eqMB = pPP->pPI->calcMatureBiomass(n3_msz,cout);
    } else { //fisheries occur AFTER molting/growth/maturity 
        if (debug) cout<<"dtF>dtM"<<endl;
        //apply natural mortality BEFORE molting/growth
        d3_array n1_msz = pPP->pPI->applyNM(pPP->dtM,n_msz,cout);
        if (debug) {cout<<"n1_msz ="<<endl; wts::print(n1_msz,cout,1);}
        
        //calculate mature biomass at mating
        eqMB = pPP->pPI->calcMatureBiomass(n1_msz,cout);
    }
    if (debug) {
        cout<<"eqMB = "<<eqMB<<endl;
        cout<<"finished Equilibrium_Calculator::calcEqMatureBiomassAtMatingFM(double R, double dirF)"<<endl;
    }
    return eqMB;
}

////////////////////////////////////////////////////////////////////////////////
//Tier3_Calculator
////////////////////////////////////////////////////////////////////////////////
/** flag to print debug info */
int Tier3_Calculator::debug = 0;
/**
 * Constructor.
 * 
 * @params XX - target SPR rate for MSY calculations
 * @param npZBs - number of size bins
 */
Tier3_Calculator::Tier3_Calculator(double XXp, Equilibrium_Calculator* pECp) : Tier_Calculator(pECp) {
    XX   = XXp;
}

/**
 * Calculate B100, which is equilibrium mature male biomass (MMB)
 * for an unfished stock.
 * 
 * Also calculates class members B100 and B0 (=B100).
 * 
 * @param R - longterm (average) recruitment (male-only, millions)
 * @param cout - output stream for debug info
 * 
 * @return MMB for unfished population (1000's t)
 */
double Tier3_Calculator::calcB100(double R, ostream& cout){
    if (debug) cout<<"starting d3_array Tier3_Calculator::calcB100(double R)"<<endl;
    //calculate mature biomass at time of mating
    B100 = pEC->calcEqMatureBiomassAtMatingF0(R,cout);
    B0 = B100;
    if (debug) cout<<"finished double Tier3_Calculator::calcB100(double R)"<<endl;
    return B100;
}

/**
 * Calculate equilibrium mature male biomass-at-mating (MMB) when directed fishery
 * is fished at longterm capture rate Fmsy. For a Tier 3 stock, 
 * Bmsy = B35% = 0.35*B100.
 * 
 * Also calculates class members B100 and B0.
* 
 * @param R - longterm (average) recruitment (male-only, millions)
 * @param cout - output stream for debug info
 * 
 * @return equilibrium MMB for population fished at Fmsy (1000's t)
 */
double Tier3_Calculator::calcBmsy(double R, ostream& cout){
    if (debug) cout<<"starting double Tier3_Calculator::calcBmsy(double R)"<<endl;
    Bmsy = XX*calcB100(R,cout);
    if (debug) cout<<"finished double Tier3_Calculator::calcBmsy(double R)"<<endl;
    return Bmsy;
}

/**
 * Calculate directed fishery capture rate, Fmsy, that results in 
 * equilibrium (longterm) MMB = Bmsy. For a Tier 3 stock, Fmsy and 
 * Bmsy are given by SPR-based proxies F35% nd B35%, where F35% is
 * the directed fishery capture rate that results in equilibrium
 * MMB = B35% = 0.35*MMB.
 * 
 * Also calculates class members: B100, B0 (=B100), and Bmsy (= XX*B100)
 * 
 * @param R - longterm (average) recruitment (male-only)
 * @param cout - output stream for debug info
 * 
 * @return Fmsy
 */
double Tier3_Calculator::calcFmsy(double R, ostream& cout){
    if (debug) cout<<"starting Tier3_Calculator::calcFmsy(double R)"<<endl;
    //calculate unfished mmb (B0)
    B0   = calcB100(R,cout);
    Bmsy = XX*B0;
    if (debug){cout<<"B100 ="<<B0<<". Bmsy = "<<Bmsy<<endl;}
    
    //From initial guess for FXX, iterate to improve FXX    
    double dF   = 0.0001;//"delta" F for derivative calculation
    double FXX  = 0.5;   //initial guess for FXX
    double dFXX = 1.0;   //"delta" FXX to update (initial value is a dummy)
    double mmbp, mmb, mmbm, dMMBdF, XXp;
    int i=0;
    while ((i<20)&&(sfabs(dFXX)>0.00001)){
        mmbp = pEC->calcEqMatureBiomassAtMatingFM(R,FXX+dF,cout);
        mmb  = pEC->calcEqMatureBiomassAtMatingFM(R,FXX,cout);
        mmbm = pEC->calcEqMatureBiomassAtMatingFM(R,FXX-dF,cout);
        dMMBdF = 0.5*(mmbp-mmbm)/dF;//derivative of mmb wrto F
        XXp   = mmb/B0;           //ratio of mmb for current FXX relative to unfished 
        dFXX  = (Bmsy - mmb)/dMMBdF;
        FXX  += dFXX;
        if (debug) {
            cout<<"--iteration = "<<++i<<endl;
            cout<<"----mmb  = "<<mmb<<". Bmsy = "<<Bmsy<<endl;
            cout<<"----XXp  = "<<XXp<<". dFXX = "<<dFXX<<endl;
            cout<<"----FXX  = "<<FXX<<endl;
        }
    }//i loop
    
    if (debug) cout<<"finished Tier3_Calculator::calcFmsy(double R)"<<endl;
    return FXX;//Fmsy for Tier 3 stock
}

////////////////////////////////////////////////////////////////////////////////
//OFL_Calculator
////////////////////////////////////////////////////////////////////////////////
/** flag to print debug info */
int OFL_Calculator::debug = 0;
/**
 * Constructor.
 * 
 * @param pTCM - pointer to a Tier_Calculator object for males
 * @param pTCF - pointer to a Tier_Calculator object for females
 */
OFL_Calculator::OFL_Calculator(Tier_Calculator* pTCMp, Tier_Calculator* pTCFp){
    //inputs
    pTCM = pTCMp;
    pTCF = pTCFp;
    
    //other constants
    alpha = 0.1; 
    beta  = 0.25;
}

/**
 * Calculate Fofl using Harvest Control Rule (HCR) for "current" MMB.
 * 
 * @param currMMB - "current" MMB used to determine B/Bmsy
 * @param Bmsy - Bmsy
 * @param Fmsy - Fmsy
 * @param cout - output stream for debug info
 * 
 * @return Fofl derived from HCR
 */
double OFL_Calculator::calcHCR(double currMMB, double Bmsy, double Fmsy, ostream& cout){
    if (debug) cout<<"starting OFL_Calculator::calcHCR(double currMMB, double Bmsy, double Fmsy)"<<endl;
    double Fofl = 0.0;
    double ratio  = currMMB/Bmsy;
    if (ratio < beta){Fofl = 0.0;} else
    if (ratio < 1.0){
        Fofl = Fmsy*(ratio-alpha)/(1-alpha);
    } else Fofl = Fmsy;
    
    if (debug) cout<<"finished OFL_Calculator::calcHCR(double currMMB, double Bmsy, double Fmsy)"<<endl;
    return Fofl;    
}

/**
 * Calculate Fofl using the Harvest Control Rule (HCR), based on Bmsy,
 * Fmsy, and the initial (July 1) male abundance n_msz. This is an 
 * iterative process because Fofl is given by the HCR based on the
 * projected MMB, which itself depends on what the target fishery capture
 * rate (F) was.
 * 
 * Also calculates prjMMB;
 * 
 * @param Bmsy - equilibrium MMB when population fished at Fmsy
 * @param Fmsy - directed fishery capture rate yielding MSY
 * @param n_msz - initial (July 1) male abundance
 * @param cout - output stream for debug info
 * 
 * @return Fofl
 */
double OFL_Calculator::calcFofl(double Bmsy, double Fmsy, d3_array& n_msz, ostream& cout){
    if (debug) {
        cout<<"starting double OFL_Calculator::calcFofl(Bmsy, Fmsy,n_msz)"<<endl;
        cout<<"Bmsy = "<<Bmsy<<"; Fmsy = "<<Fmsy<<endl;
        cout<<"n_msz = "<<endl;
        for (int m=1;m<=tcsam::nMSs;m++) {cout<<"m = "<<m<<endl; cout<<n_msz(m)<<endl;}
    }
    
    PopProjector* pPrjM = pTCM->pEC->pPP;
    
    //start with guess for Fofl based on currMMB
    double currMMB = pPrjM->projectMatureBiomassAtMating(Fmsy,n_msz,cout);
    if (debug) cout<<"init currMMB = "<<currMMB<<"; B/Bmsy = "<<currMMB/Bmsy<<endl;
    double Fofl = calcHCR(currMMB,Bmsy,Fmsy,cout);
    if (debug) cout<<"init Fofl = "<<Fofl<<endl;
    //now iterate until Fofl yields currMMB
    double Foflp = 0.0;
    double criF = 0.001;
    double delF = std::numeric_limits<double>::infinity();
    while (sfabs(delF)>criF){
        //calculate currMMB based on Fofl
        currMMB = pPrjM->projectMatureBiomassAtMating(Fofl,n_msz,cout);
        if (debug) cout<<"--updated currMMB = "<<currMMB<<"; B/Bmsy = "<<currMMB/Bmsy<<endl;
        //update Fofl based on currMMB
        Foflp = Fofl;
        Fofl  = calcHCR(currMMB,Bmsy,Fmsy,cout);
        //check convergence
        delF = Fofl - Foflp;
        if (debug) cout<<"--updated Fofl = "<<Fofl<<"; delF = "<<delF<<endl;
    }
    prjMMB = currMMB;
    if (debug) cout<<"finished double OFL_Calculator::calcFofl(Bmsy, Fmsy,n_msz)"<<endl;
    return Fofl;
}
/**
 * Calculate MSY, the maximum sustainable yield, based on the 
 * longterm sex-specific recruitment level (R_y) and Fmsy. 
 * 
 * @param R_x - longterm (average) sex-specific recruitment
 * @param Fmsy - directed fishery F rate at MSY
 * @param cout - output stream for debug info
 * 
 * @return the MSY
 */
double OFL_Calculator::calcMSY(dvector R_x, double Fmsy, ostream& cout){
    if (debug) cout<<"starting double OFL_Calculator::calcMSY(R_x,Fmsy)"<<endl;
    int nFsh = pTCM->pEC->pPP->nFsh;
    if (debug) cout<<"nFsh = "<<nFsh<<endl;
    
    msy_fx.allocate(0,nFsh,1,tcsam::nSXs);
    msy_fx.initialize();
    
    double totCM;
        
    //do males
    Equilibrium_Calculator* pECM = pTCM->pEC;
    PopProjector*           pPPM = pTCM->pEC->pPP;
    //calculate equilibrium numbers-at-size for males
    d3_array n_msz = pECM->calcEqNatZFM(R_x(MALE),Fmsy,cout);
    //calculate associated catch
    pPPM->project(Fmsy,n_msz,cout);
    //calc retained catch biomass (directed fishery only) from catch abundance for males
    msy_fx(0,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->rm_fmsz(1),cout);
    //calc discard mortality biomass from catch abundance for males
    for (int f=1;f<=nFsh;f++){
        msy_fx(f,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->dm_fmsz(f),cout);
    }
    totCM = pPPM->pPI->calcTotalBiomass(pPPM->pCI->cm_msz,cout);    
    
    if (tcsam::nSXs>1){
        //do females
        Equilibrium_Calculator* pECF = pTCF->pEC;
        PopProjector*           pPPF = pTCF->pEC->pPP;
        //calculate equilibrium numbers-at-size for females
        d3_array n_msz = pECF->calcEqNatZFM(R_x(FEMALE),Fmsy,cout);
        //calculate associated catch
        pPPF->project(Fmsy,n_msz,cout);
        //calc retained catch biomass (directed fishery only) from catch abundance for females
        msy_fx(0,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->rm_fmsz(1),cout);
        //calc discard mortality biomass from catch abundance for females
        for (int f=1;f<=nFsh;f++){
            msy_fx(f,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->dm_fmsz(f),cout);
        }
        totCM += pPPF->pPI->calcTotalBiomass(pPPF->pCI->cm_msz,cout);    
    }
    
    double msy = sum(msy_fx);
    if (debug||(msy!=totCM)) {
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"in double OFL_Calculator::calcMSY(R_x,Fmsy,n_xmsz)"<<endl;
        cout<<"MSY   = "<<msy<<endl;
        cout<<"totCM = "<<totCM<<endl;
        cout<<"msy_fx: "<<endl<<msy_fx<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"finished double OFL_Calculator::calcMSY(R_x,Fmsy)"<<endl;
    }
    return(msy);
}

/**
 * Calculate the total OFL (retained+discard mortality) taken
 * when fishing on population starting with initial abundance (July 1)
 * n_xmsz at directed fishery capture rate of Fofl.
 * 
 * Note: catch (biomass) by fishery and sex is available in ofl_fx after function exits:
 *      ofl_fx(0,x) - retained catch (biomass) in directed fishery, by sex x
 *      ofl_fx(f,x) - discard mortality (biomass) for fishery f, sex x
 * 
 * @param Fofl - directed fishery capture rate
 * @param n_xmsz - initial population abundance
 * @param cout - output stream for debug info
 * 
 * @return total OFL (biomass)
 */
double OFL_Calculator::calcOFL(double Fofl, d4_array& n_xmsz, ostream& cout){
    if (debug) cout<<"starting double OFL_Calculator::calcOFL(Fofl,n_xmsz)"<<endl;
    
    int nFsh = pTCM->pEC->pPP->nFsh;
    if (debug) cout<<"nFsh = "<<nFsh<<endl;
    
    ofl_fx.allocate(0,nFsh,1,tcsam::nSXs);
    ofl_fx.initialize();

    double totCM;
    
    //do males
    PopProjector* pPPM = pTCM->pEC->pPP;
    //calc catch abundance of males at Fofl
    if (debug) PopProjector::debug=1;
    pPPM->project(Fofl,n_xmsz(  MALE),cout);
    //calc retained catch biomass (directed fishery only) from catch abundance for males
    ofl_fx(0,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->rm_fmsz(1),cout);
    //calc discard mortality biomass from catch abundance for males
    for (int f=1;f<=nFsh;f++){
        ofl_fx(f,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->dm_fmsz(f),cout);
    }
    totCM = pPPM->pPI->calcTotalBiomass(pPPM->pCI->cm_msz,cout);    
    
    if (tcsam::nSXs>1){
        //do females
        PopProjector* pPPF = pTCF->pEC->pPP;
        //calc catch abundance of females at Fofl
        pPPF->project(Fofl,n_xmsz(FEMALE),cout);
        //calc retained catch biomass (directed fishery only) from catch abundance for females
        ofl_fx(0,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->rm_fmsz(1),cout);
        //calc discard mortality biomass from catch abundance for females
        for (int f=1;f<=nFsh;f++){
            ofl_fx(f,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->dm_fmsz(f),cout);
        }
        totCM += pPPF->pPI->calcTotalBiomass(pPPF->pCI->cm_msz,cout);    
    }
    if (debug) PopProjector::debug=0;
    
    double ofl = sum(ofl_fx);
    if (debug||(ofl!=totCM)) {
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"in double OFL_Calculator::calcOFL(Fofl,n_xmsz)"<<endl;
        cout<<"OFL   = "<<ofl<<endl;
        cout<<"totCM = "<<totCM<<endl;
        cout<<"ofl_fx: "<<endl<<ofl_fx<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"finished double OFL_Calculator::calcOFL(Fofl,n_xmsz)"<<endl;
    }
    return ofl;
}
/**
 * Calculate the (projected) MMB when fished at Fofl.
 * 
 * @param Fofl - target fishery capture rate
 * @param n_msz - initial (July 1) population
 * @param cout - output stream for debug info
 * 
 * @return - the (projected) MMB
 */
double OFL_Calculator::calcPrjMMB(double Fofl, d3_array& n_msz, ostream& cout){
    if (debug) cout<<"starting double OFL_Calculator::calcPrjMMB(Fofl,n_msz)"<<endl;
    PopProjector* pPrjM = pTCM->pEC->pPP;
    double mmb = pPrjM->projectMatureBiomassAtMating(Fofl,n_msz,cout);
    if (debug) cout<<"finished double OFL_Calculator::calcPrjMMB(Fofl,n_msz)"<<endl;
    return mmb;
}

/**
 * Calculate all results associated with the OFL.
 * 
 * @param R - assumed sex-specific equilibrium recruitment
 * @param n_xmsz - d4_array with initial population abundance
 * @param cout - output stream for debug info
 * 
 * @return OFLResults object.
 */
OFLResults OFL_Calculator::calcOFLResults(dvector R, d4_array& n_xmsz, ostream& cout){
    if (debug) cout<<"starting OFLResults OFL_Calculator::calcOFLResults(R,n_xmsz,cout)"<<endl;
    int nFsh = pTCM->pEC->pPP->pCI->nFsh;
    OFLResults res;
    
    res.avgRec_x = R;
    res.curB     = pTCM->pEC->pPP->pPI->calcMatureBiomass(n_xmsz(MALE),cout);
    if (debug) cout<<"calcOFLResults: calculated curB"<<endl;
    
    res.eqNatZF0_xmsz.allocate(1,tcsam::nSXs,
                               1,tcsam::nMSs,
                               1,tcsam::nSCs,
                               1,pTCM->pEC->pPP->nZBs);
    if (debug) {
        cout<<"calcOFLResults: allocated eq NatZ for F=0"<<endl;
        Equilibrium_Calculator::debug=1;
    }
    res.eqNatZF0_xmsz(MALE) = pTCM->pEC->calcEqNatZF0(R(MALE),cout);
    if (debug) cout<<"calcOFLResults: calculated eq NatZ(MALE) for F=0"<<endl;
    if (tcsam::nSXs>1) {
        res.eqNatZF0_xmsz(FEMALE) = pTCF->pEC->calcEqNatZF0(R(FEMALE),cout);
        if (debug) cout<<"calcOFLResults: calculated eq NatZ(FEMALE) for F=0"<<endl;
    }
    if (debug) Equilibrium_Calculator::debug=0;
            
    res.Fmsy = pTCM->calcFmsy(R(MALE),cout);//also calculates B0 and Bmsy
    res.B0   = pTCM->B0;
    res.Bmsy = pTCM->Bmsy;
    res.MSY  = calcMSY(R,res.Fmsy,cout);
    if (debug) cout<<"calcOFLResults: calculated MSY"<<endl;
    
    res.eqNatZFM_xmsz.allocate(1,tcsam::nSXs,
                               1,tcsam::nMSs,
                               1,tcsam::nSCs,
                               1,pTCM->pEC->pPP->nZBs);
    if (debug) {
        cout<<"calcOFLResults: allocated eq NatZ for F=Fmsy"<<endl;
        Equilibrium_Calculator::debug=1;
    }
    res.eqNatZFM_xmsz(MALE) = pTCM->pEC->calcEqNatZFM(R(MALE),res.Fmsy,cout);
    if (debug) cout<<"calcOFLResults: calculated eq NatZ(MALE) for F=Fmsy"<<endl;
    if (tcsam::nSXs>1) {
        res.eqNatZFM_xmsz(FEMALE) = pTCF->pEC->calcEqNatZFM(R(FEMALE),res.Fmsy,cout);
        if (debug) cout<<"calcOFLResults: calculated eq NatZ(FEMALE) for F=Fmsy"<<endl;
    }
    if (debug) Equilibrium_Calculator::debug=0;
            
    res.Fofl = calcFofl(res.Bmsy,res.Fmsy,n_xmsz(MALE),cout);
    res.prjB = prjMMB;
    res.OFL  = calcOFL(res.Fofl,n_xmsz,cout);
    if (debug) cout<<"calcOFLResults: calculated OFL"<<endl;
    
    res.ofl_fx.allocate(0,nFsh,1,tcsam::nSXs);
    res.ofl_fx = ofl_fx;
    if (debug) cout<<"finished OFLResults OFL_Calculator::calcOFLResults(R,n_xmsz,cout)"<<endl;
    return res;
}
////////////////////////////////////////////////////////////////////////////////
//OFLResults
////////////////////////////////////////////////////////////////////////////////
int OFLResults::debug = 0;
/**
 * Assignment operator for OFLResults class.
 * 
 * @param OFLResults object to copy.
 * @return reference to the copied object
 */
OFLResults& OFLResults::operator=(const OFLResults& o){
    if (debug) std::cout<<"starting OFLResults::operator=(const OFLResults&)"<<endl;
    avgRec_x = o.avgRec_x;//average recruitment, by sex
    B0       = o.B0;      //equilibrium MMB for unfished population
    Fmsy     = o.Fmsy;    //equilibrium F on directed fishery for males resulting in MSY
    Bmsy     = o.Bmsy;    //equilibrium MMB when fished at Fmsy 
    MSY      = o.MSY;     //equilibrium yield (1000's t) when fished at Fmsy
    Fofl     = o.Fofl;    //F on directed fishery for males resulting in the OFL
    OFL      = o.OFL;     //total OFL (1000's t)
    ofl_fx   = o.ofl_fx;  //fishery/sex-specific mortality components to OFL (f=0 is retained catch, f>0 is total catch mortality)
    prjB     = o.prjB;    //projected MMB for projection year when current population is fished at Fofl.
    curB     = o.curB;    //"current" MMB at beginning of projection year
    eqNatZF0_xmsz.deallocate(); //unfished equilibrium size distribution
    eqNatZF0_xmsz.allocate(o.eqNatZF0_xmsz);
    if (debug) std::cout<<"got here 1"<<endl;
    eqNatZF0_xmsz = o.eqNatZF0_xmsz;
    eqNatZFM_xmsz.deallocate(); //unfished equilibrium size distribution
    eqNatZFM_xmsz.allocate(o.eqNatZF0_xmsz);
    if (debug) std::cout<<"got here 2"<<endl;
    eqNatZFM_xmsz = o.eqNatZF0_xmsz;
    if (debug) std::cout<<"finished OFLResults::operator=(const OFLResults&)"<<endl;
    return *this;
}
/**
 * Write csv header for OFL results to output stream.
 *  
 * @param os - output stream to write to
 */
void OFLResults::writeCSVHeader(ostream& os){
    os<<"avgRecM";
    if (tcsam::nSXs>1) os<<", AvgRecF";
    os<<", B0, Bmsy, Fmsy, MSY, Fofl, OFL, prjB, prjB/Bmsy, curB, prjB/curB";
}

/**
 * Write values to output stream in csv format
 * 
 * @param cout - output stream to write to
 */
void OFLResults::writeToCSV(ostream& os){
    os<<avgRec_x(  MALE);
    if (tcsam::nSXs>1) os<<cc<<avgRec_x(FEMALE);
    os<<cc<<B0<<cc<<Bmsy<<cc<<Fmsy<<cc<<MSY<<cc<<Fofl<<cc<<OFL<<cc;
    os<<prjB<<cc<<prjB/Bmsy<<cc<<curB<<cc<<prjB/curB;
}
/**
 * Write values as R list to output stream
 * 
 * @param os - output stream to write to
 * @param ptrMC - pointer to ModelConfiguration object
 * @param name - name for R list
 * @param debug - flag to print debugging info
 */
void OFLResults::writeToR(ostream& os, ModelConfiguration* ptrMC, adstring name, int debug){
    adstring xDms = ptrMC->dimSXsToR;//sex
    adstring mDms = ptrMC->dimMSsToR;//maturity
    adstring sDms = ptrMC->dimSCsToR;//shell condition
    adstring zDms = ptrMC->dimZBsToR;//size bin midpoints
    os<<name<<"=list(OFL="<<OFL<<cc<<"Fofl="<<Fofl<<cc<<"prjB="<<prjB<<cc<<"curB="<<curB;
    os<<cc<<"Fmsy="<<Fmsy<<cc<<"Bmsy="<<Bmsy<<cc<<"MSY="<<MSY;
    os<<cc<<"B100="<<B0<<cc<<"avgRecM="<<avgRec_x(MALE);
    if (tcsam::nSXs>1) os<<cc<<"avgRecF="<<avgRec_x(FEMALE);
    os<<cc<<"eqNatZF0_xmsz="; wts::writeToR(os,eqNatZF0_xmsz,xDms,mDms,sDms,zDms); 
    os<<cc<<"eqNatZFM_xmsz="; wts::writeToR(os,eqNatZFM_xmsz,xDms,mDms,sDms,zDms); 
    os<<")";
}
