#include <limits>
#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelPopDyClasses.hpp"

using namespace tcsam;

/**flags to print debug info*/
int PopDyInfo::debug = 0;
int CatchInfo::debug = 0;
int PopProjector::debug = 0;
int MultiYearPopProjector::debug = 0;
////////////////////////////////////////////////////////////////////////////////
//PopDyInfo
////////////////////////////////////////////////////////////////////////////////
void PopDyInfo::allocate(){
    if (debug) cout<<"starting PopDyInfo::allocate()"<<endl;
    R_z.deallocate();
    w_mz.deallocate();
    Th_sz.deallocate();
    M_msz.deallocate();
    T_szz.deallocate();
    
    np_msz.deallocate();
    S_msz.deallocate();
    
    R_z.allocate(1,nZBs);
    w_mz.allocate(1,nMSs,1,nZBs);
    Th_sz.allocate(1,nSCs,1,nZBs);
    M_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    T_szz.allocate(1,nSCs,1,nZBs,1,nZBs);
    
    np_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    S_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    if (debug) cout<<"finished PopDyInfo::allocate()"<<endl;
}
/**
 * Constructor.
 * 
 * @param npZBs - number of size bins
 */
PopDyInfo::PopDyInfo(int npZBs){
    nZBs = npZBs;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    allocate();
}

PopDyInfo::PopDyInfo(const PopDyInfo& o): PopDyInfo(o.nZBs) {
//    nZBs = o.nZBs;
//    
//    nMSs = tcsam::nMSs;
//    nSCs = tcsam::nSCs;
    if (debug) cout<<"starting PopDyInfo::PopDyInfo(const PopDyInfo& o): PopDyInfo(o.nZBs)"<<endl;
    allocate();
    
    R_z   = 1.0*o.R_z;
    w_mz  = 1.0*o.w_mz;
    Th_sz = 1.0*o.Th_sz;
    M_msz = 1.0*o.M_msz;
    T_szz = 1.0*o.T_szz;
    
    np_msz = 1.0*o.np_msz;
    S_msz  = 1.0*o.S_msz;
    if (debug) cout<<"this = "<<this<<"; o"<<&o<<endl;
    if (debug) cout<<"finished PopDyInfo::PopDyInfo(const PopDyInfo& o): PopDyInfo(o.nZBs)"<<endl;
}

PopDyInfo& PopDyInfo::operator=(const PopDyInfo& o) {
    nZBs = o.nZBs;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    allocate();
    
    R_z   = 1.0*o.R_z;
    w_mz  = 1.0*o.w_mz;
    Th_sz = 1.0*o.Th_sz;
    M_msz = 1.0*o.M_msz;
    T_szz = 1.0*o.T_szz;
    
    np_msz = 1.0*o.np_msz;
    S_msz  = 1.0*o.S_msz;
    
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
    if (debug) cout<<"starting PopDyInfo::calcSurvival(dt); dt="<<dt<<endl;
    RETURN_ARRAYS_INCREMENT();
    S_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            S_msz(m,s) = exp(-M_msz(m,s)*dt); //survival over dt
        }//m
    }//s  
    if (debug) cout<<"finished PopDyInfo::calcSurvival(dt)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return 1.0*S_msz;
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
    if (debug) cout<<"starting PopDyInfo::applyNM(dt,n_msz); dt="<<dt<<endl;
    RETURN_ARRAYS_INCREMENT();
    np_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            np_msz(m,s) = elem_prod(exp(-M_msz(m,s)*dt),n_msz(m,s)); //survival over dt
        }//m
    }//s  
    if (debug) {
        cout<<"np_msz = "<<endl; wts::print(np_msz,cout,1); cout<<endl;
    }
    if (debug) cout<<"finished PopDyInfo::applyNM(dt,n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return 1.0*np_msz;
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
    if (debug) {
        cout<<"np_msz = "<<endl; wts::print(np_msz,cout,1); cout<<endl;
    }
    if (debug) cout<<"finished PopDyInfo::applyMG(n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return 1.0*np_msz;
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
    if (debug) {
        cout<<"#--n_msz  = "<<&n_msz<<endl;
        cout<<"#--np_msz = "<<&np_msz<<endl;
        cout<<"np_msz = "<<endl; wts::print(np_msz,cout,1); cout<<endl;
    }

    if (debug) cout<<"finished PopDyInfo::addRecruitment(R, n_msz)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return 1.0*np_msz;
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
    os<<"S_msz="; wts::writeToR(os,value(S_msz),mDms,sDms,zDms);  os<<cc<<endl;
    os<<"Th_sz="; wts::writeToR(os,value(Th_sz),sDms,zDms);       os<<cc<<endl;
    os<<"T_szz="; wts::writeToR(os,value(T_szz),sDms,zDms,zpDms); os<<endl;
    os<<")";
    if (dbg) std::cout<<"Finished PopDyInfo::writeToR()"<<endl;
}
        
////////////////////////////////////////////////////////////////////////////////
//CatchInfo
////////////////////////////////////////////////////////////////////////////////
void CatchInfo::allocate(){
    if (debug) cout<<"starting CatchInfo::allocate()"<<endl;
    hm_f.deallocate();
    cpF_fms.deallocate();
    cpF_fmsz.deallocate();
    selF_fmsz.deallocate();
    retF_fmsz.deallocate();
    
    rmF_fmsz.deallocate();
    dmF_fmsz.deallocate();
    
    cmN_msz.deallocate();
    cpN_fmsz.deallocate();
    rmN_fmsz.deallocate();
    dmN_fmsz.deallocate();
    
    totFM_z.deallocate();
    S_msz.deallocate();
    
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
    
    totFM_z.allocate(1,nZBs);
    S_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    if (debug) cout<<"finished CatchInfo::allocate()"<<endl;
}
/**
 * Constructor.
 * 
 * @param npZBs - number of size bins
 * @param npFsh - number of fisheries
 */
CatchInfo::CatchInfo(int npZBs, int npFsh){
    if (debug) cout<<"starting CatchInfo::CatchInfo(int npZBs, int npFsh)"<<endl;
    nZBs = npZBs;
    nFsh = npFsh;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    allocate();
    
    maxF = 1.0;//default scale    
    
    if (debug) cout<<"finished CatchInfo::CatchInfo(int npZBs, int npFsh)"<<endl;
}

    /**
     * Class copy constructor.
     * 
     * @param o - CatchInfo object
     */
CatchInfo::CatchInfo(const CatchInfo& o): CatchInfo(o.nZBs,o.nFsh){
//    nZBs = o.nZBs;
//    nFsh = o.nFsh;
//    
//    nMSs = tcsam::nMSs;
//    nSCs = tcsam::nSCs;
    if (debug) cout<<"starting CatchInfo::CatchInfo(const CatchInfo& o): CatchInfo(o.nZBs,o.nFsh)"<<endl;
    allocate();
    
    maxF = o.maxF;//default scale
    
    hm_f      = 1.0*o.hm_f;
    cpF_fms   = 1.0*o.cpF_fms;
    for (int f=1;f<=nFsh;f++){
        cpF_fmsz(f)  = 1.0*o.cpF_fmsz(f);
        selF_fmsz(f) = 1.0*o.selF_fmsz(f);
        retF_fmsz(f) = 1.0*o.retF_fmsz(f);

        rmF_fmsz(f) = 1.0*o.rmN_fmsz(f);
        dmF_fmsz(f) = 1.0*o.dmN_fmsz(f);

        cpN_fmsz(f) = 1.0*o.cpN_fmsz(f);
        rmN_fmsz(f) = 1.0*o.rmN_fmsz(f);
        dmN_fmsz(f) = 1.0*o.dmN_fmsz(f);
    }
    
    totFM_z = 1.0*o.totFM_z;
    S_msz   = 1.0*o.S_msz;
    cmN_msz = 1.0*o.cmN_msz;
    if (debug) cout<<"this = "<<this<<"; o"<<&o<<endl;
    if (debug) cout<<"starting CatchInfo::CatchInfo(const CatchInfo& o): CatchInfo(o.nZBs,o.nFsh)"<<endl;
}

CatchInfo& CatchInfo::operator =(const CatchInfo& o){
    nZBs = o.nZBs;
    nFsh = o.nFsh;
    
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    
    allocate();
    
    maxF = o.maxF;//default scale
    
    hm_f      = 1.0*o.hm_f;
    cpF_fms   = 1.0*o.cpF_fms;
    for (int f=1;f<=nFsh;f++){
        cpF_fmsz(f)  = 1.0*o.cpF_fmsz(f);
        selF_fmsz(f) = 1.0*o.selF_fmsz(f);
        retF_fmsz(f) = 1.0*o.retF_fmsz(f);

        rmF_fmsz(f) = 1.0*o.rmN_fmsz(f);
        dmF_fmsz(f) = 1.0*o.dmN_fmsz(f);

        cpN_fmsz(f) = 1.0*o.cpN_fmsz(f);
        rmN_fmsz(f) = 1.0*o.rmN_fmsz(f);
        dmN_fmsz(f) = 1.0*o.dmN_fmsz(f);
    }
    
    totFM_z = 1.0*o.totFM_z;
    S_msz   = 1.0*o.S_msz;
    cmN_msz = 1.0*o.cmN_msz;
    
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
    if (debug) cout<<"starting CatchInfo::applyFM(dvariable dirF, dvar3_array n_msz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvariable ratF = 1.0;        //default target fishery (f=1) scaling ratio
    if ((dirF>=0.0)&&(maxF>0.0)) 
        ratF = dirF/maxF;        //target fishery (f=1) scaling ratio, if target fishery is open   (maxF>0)
    dvar3_array np_msz(1,nMSs,1,nSCs,1,nZBs);
    np_msz.initialize();//number surviving fisheries
    cmN_msz.initialize();//total catch mortality (abundance)
    if (nFsh){
        rmF_fmsz.initialize();//retained catch mortality rates
        dmF_fmsz.initialize();//discards catch mortality rates
        cpN_fmsz.initialize();//capture abundance, by fishery
        rmN_fmsz.initialize();//retained catch mortality (abundance)
        dmN_fmsz.initialize();//discards catch mortality (abundance)

        dvector     tdF_z(1,nZBs);//for use in calculating fishing rate components
        dvar_vector tvF_z(1,nZBs);//for use in calculating fishing rate components
        dvar_vector tfF_z(1,nZBs);//for use in calculating fishing rate components
        for (int s=1;s<=nSCs;s++){
            for (int m=1;m<=nMSs;m++){ 
                //calculate fishery rates
                totFM_z.initialize();
                //--directed fishery
                rmF_fmsz(1,m,s) = elem_prod(retF_fmsz(1,m,s),              ratF*cpF_fmsz(1,m,s));
                dmF_fmsz(1,m,s) = elem_prod(hm_f(1)*(1.0-retF_fmsz(1,m,s)),ratF*cpF_fmsz(1,m,s));
                totFM_z += rmF_fmsz(1,m,s)+dmF_fmsz(1,m,s);
                //--bycatch fisheries
                for (int f=2;f<=nFsh;f++){
                    rmF_fmsz(f,m,s) = elem_prod(retF_fmsz(f,m,s),              cpF_fmsz(f,m,s));
                    dmF_fmsz(f,m,s) = elem_prod(hm_f(f)*(1.0-retF_fmsz(f,m,s)),cpF_fmsz(f,m,s));
                    totFM_z        += rmF_fmsz(f,m,s)+dmF_fmsz(f,m,s);
                }

                //calculate numbers surviving and numbers killed
                np_msz(m,s) = elem_prod(mfexp(-totFM_z),n_msz(m,s));//survival after all fisheries
                cmN_msz(m,s) = n_msz(m,s)-np_msz(m,s);              //total catch mortality, all fisheries

                //calculate capture abundance, retained abundance, and discard abundance mortality
                tdF_z = value(totFM_z);
                tvF_z = elem_prod(1-wts::isEQ(tdF_z,0.0),totFM_z) + wts::isEQ(tdF_z,0.0);//= totFM_z(z)                     if totFM_z(z) > 0, else = 1
                tfF_z = elem_div(1.0-mfexp(-totFM_z),tvF_z);                             //= (1-exp(-totFM_z(z))/totFM_z(z) if totFM_z(z) > 0, else = 1
                //--directed fishery
                cpN_fmsz(1,m,s) = elem_prod(elem_prod(ratF*cpF_fmsz(1,m,s),tfF_z),n_msz(m,s));
                //retained catch mortality (abundance), directed fishery
                rmN_fmsz(1,m,s) = elem_prod(elem_prod(     rmF_fmsz(1,m,s),tfF_z),n_msz(m,s));
                //discard catch mortality (abundance), directed fishery
                dmN_fmsz(1,m,s) = elem_prod(elem_prod(     dmF_fmsz(1,m,s),tfF_z),n_msz(m,s));
                //--bycatch fisheries
                for (int f=2;f<=nFsh;f++){
                    //total capture abundance
                    cpN_fmsz(f,m,s) = elem_prod(elem_prod(cpF_fmsz(f,m,s),tfF_z),n_msz(m,s));
                    //retained catch mortality (abundance)
                    rmN_fmsz(f,m,s) = elem_prod(elem_prod(rmF_fmsz(f,m,s),tfF_z),n_msz(m,s));
                    //discard catch mortality (abundance)
                    dmN_fmsz(f,m,s) = elem_prod(elem_prod(dmF_fmsz(f,m,s),tfF_z),n_msz(m,s));
                }//f
            }//m
        }//s
    } else {
        if (debug) cout<<"no fisheries defined"<<endl;
        for (int s=1;s<=nSCs;s++){
            for (int m=1;m<=nMSs;m++) np_msz(m,s) = n_msz(m,s);
        }
    }
    if (debug) {
        cout<<"#--In CatchInfo::calcCatch"<<endl;
        cout<<"dirF, maxCapF,ratF    = "<<dirF<<cc<<maxF<<cc<<ratF<<endl;
        for (int s=1;s<=nSCs;s++){
            for (int m=1;m<=nMSs;m++){ 
                cout<<"m,s     = "<<m<<cc<<s<<endl;
                cout<<"n_msz   = "<<n_msz(m,s)<<endl;
                cout<<"cmN_msz = "<<cmN_msz(m,s)<<endl;
                for (int f=1;f<=nFsh;f++){
                    cout<<"f = "<<f<<endl;
                    cout<<"cpF_fmsz = "<<cpF_fmsz(f,m,s)<<endl;
                    cout<<"cpN_fmsz = "<<cpN_fmsz(f,m,s)<<endl;
                    cout<<"ret_fmsz = "<<retF_fmsz(f,m,s)<<endl;
                    cout<<"rmF_fmsz = "<<rmF_fmsz(f,m,s)<<endl;
                    cout<<"rmN_fmsz = "<<rmN_fmsz(f,m,s)<<endl;
                    cout<<"dmF_fmsz = "<<dmF_fmsz(f,m,s)<<endl;
                    cout<<"dmN_fmsz = "<<dmN_fmsz(f,m,s)<<endl;
                }
            }
        }
    }
    if (debug) cout<<"finished CatchInfo::calcCatch(dvariable dirF, dvar3_array n_msz)"<<endl;
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
    if (debug) cout<<"starting CatchInfo::calcSurvival(double dirF)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvariable ratF = 1.0;     //default target fishery scaling ratio
    if ((dirF>=0)&&(maxF>0)) 
        ratF = dirF/maxF;     //target fishery (f=1) scaling ratio
    S_msz.initialize();
    for (int s=1;s<=nSCs;s++){
        for (int m=1;m<=nMSs;m++){ 
            totFM_z.initialize();
            totFM_z += elem_prod(ratF*cpF_fmsz(1,m,s),
                                 retF_fmsz(1,m,s) + hm_f(1)*(1.0-retF_fmsz(1,m,s)));
            for (int f=2;f<=nFsh;f++)
                totFM_z += elem_prod(cpF_fmsz(f,m,s),
                                     retF_fmsz(f,m,s) + hm_f(f)*(1.0-retF_fmsz(f,m,s)));
            S_msz(m,s) = exp(-totFM_z);                //survival of fisheries
        }//m
    }//s
    if (debug) cout<<"#--S_msz = "<<&S_msz<<endl;
    if (debug) cout<<"finished CatchInfo::calcSurvival(double dirF)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return 1.0*S_msz;
}

/**
 * Set single-sex size-specific fishery capture rates.
 * 
 * @param capF_fmsz - input capture rates
 */
void CatchInfo::setCaptureRates(dvar4_array& capFp_fmsz){
    if (debug) cout<<"starting CatchInfo::setCaptureRates(dvar4_array& capFp_fmsz)"<<endl;
    for (int f=1;f<=nFsh;f++) cpF_fmsz(f) = 1.0*capFp_fmsz(f);
    if (debug) cout<<"#--cpF_fmsz = "<<&cpF_fmsz<<endl<<"capFp_fmsz = "<<&capFp_fmsz<<endl;
    if (debug) cout<<"finished CatchInfo::setCaptureRates(dvar4_array& capFp_fmsz)"<<endl;
}

/**
 * Set single-sex fishery capture rates.
 * 
 * @param capFp_fms - input capture rates
 */
void CatchInfo::setCaptureRates(dvar3_array& capFp_fms){
    if (debug) cout<<"starting CatchInfo::setCaptureRates(dvar3_array& capFp_fms)"<<endl;
    cpF_fms = 1.0*capFp_fms;
    if (debug) cout<<"#--cpF_fms = "<<&cpF_fms<<endl<<"capFp_fms = "<<&capFp_fms<<endl;
    if (debug) cout<<"finished CatchInfo::setCaptureRates(dvar3_array& capFp_fms)"<<endl;
}

/**
 * Set single-sex selectivity functions.
 * 
 * @param selFp_fmsz - input selectivity functions
 */
void CatchInfo::setSelectivityFcns(dvar4_array& selFp_fmsz){
    if (debug) cout<<"starting CatchInfo::setSelectivityFcns(dvar3_array& selFp_fmsz)"<<endl;
    for (int f=1;f<=nFsh;f++) selF_fmsz(f) = 1.0*selFp_fmsz(f);
    if (debug) cout<<"#--selF_fmsz = "<<&selF_fmsz<<endl<<"selFp_fmsz = "<<&selFp_fmsz<<endl;
    if (debug) cout<<"finished CatchInfo::setSelectivityFcns(dvar3_array& selFp_fmsz)"<<endl;
}

/**
 * Set single-sex retention functions.
 * 
 * @param retFp_fmsz - input retention functions
 */
void CatchInfo::setRetentionFcns(dvar4_array& retFp_fmsz){
    if (debug) cout<<"finished CatchInfo::setRetentionFcns(dvar3_array& retFp_fmsz)"<<endl;
    for (int f=1;f<=nFsh;f++) retF_fmsz(f) = retFp_fmsz(f);
    if (debug) cout<<"#--retF_fmsz = "<<&retF_fmsz<<endl<<"retFp_fmsz = "<<&retFp_fmsz<<endl;
    if (debug) cout<<"finished CatchInfo::setRetentionFcns(dvar3_array& retFp_fmsz)"<<endl;
}

/**
 * Set handling mortality rates, by fishery.
 * 
 * @param pHM_f - dvector of handling mortality rates
 */
void CatchInfo::setHandlingMortality(dvar_vector& pHM_f){
    if (debug) cout<<"starting CatchInfo::setHandlingMortality(dvar_vector& pHM_f)"<<endl;
    hm_f = 1.0*pHM_f;
    if (debug) cout<<"#--hm_f = "<<&hm_f<<endl<<"pHM_f = "<<&pHM_f<<endl;
    if (debug) cout<<"finished CatchInfo::setHandlingMortality(dvar_vector& pHM_f)"<<endl;
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
/**
 * Constructor
 * 
 * @param pPIp - pointer to PopDyInfo object
 * @param pCIp - pointer to CatchInfo object
 * 
 * NOTE: this does NOT set dtF or dtM.
 */
PopProjector::PopProjector(PopDyInfo* pPIp, CatchInfo* pCIp){
    if (debug) cout<<"Starting PopProjector::PopProjector(PopDyInfo* pPIp, CatchInfo* pCIp)"<<endl;
    pPI = new PopDyInfo(*pPIp);
    pCI = new CatchInfo(*pCIp);
    if (debug) {
        cout<<"got 1"<<endl;
        cout<<"pPI="<<pPI<<cc<<"pPIp="<<pPIp<<endl;
        cout<<"pCI="<<pCI<<cc<<"pCIp="<<pCIp<<endl;
    }
    
    nMSs = pPI->nMSs;
    nSCs = pPI->nSCs;
    nZBs = pPI->nZBs;
    nFsh = pCI->nFsh;  
    if (debug) cout<<"got 2"<<endl;
    
    n1_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n2_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n3_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n4_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n5_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    if (debug) cout<<"Finished PopProjector::PopProjector(PopDyInfo* pPIp, CatchInfo* pCIp)"<<endl;
}

/**
 * Copy constructor
 * 
 * @param o - PopProjector object
 */
PopProjector::PopProjector(const PopProjector& o){
    if (debug) cout<<"Starting PopProjector::PopProjector(const PopProjector& o)"<<endl;
    pPI = new PopDyInfo(*(o.pPI));
    pCI = new CatchInfo(*(o.pCI));
    if (debug) {
        cout<<"got 1"<<endl;
        cout<<"pPI="<<pPI<<cc<<"o.pPI="<<o.pPI<<endl;
        cout<<"pCI="<<pCI<<cc<<"o.pCI="<<o.pCI<<endl;
    }
    dtF = o.dtF;//time of fishery
    dtM = o.dtM;//time of mating/growth
    
    nMSs = pPI->nMSs;
    nSCs = pPI->nSCs;
    nZBs = pPI->nZBs;
    nFsh = pCI->nFsh;  
    
    n1_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n2_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n3_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n4_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    n5_msz.allocate(1,nMSs,1,nSCs,1,nZBs);
    if (debug) cout<<"Finished PopProjector::PopProjector(const PopProjector& o)"<<endl;
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
        if (debug) cout<<1<<endl;
        //apply fisheries
        n2_msz = pCI->applyFM(dirF, n1_msz,cout);
        if (debug) cout<<2<<endl;
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
        cout<<"n_msz  = "<<&n_msz<<endl;
        cout<<"n1_msz = "<<&n1_msz<<endl;
        cout<<"n2_msz = "<<&n2_msz<<endl;
        cout<<"n3_msz = "<<&n3_msz<<endl;
        cout<<"n4_msz = "<<&n4_msz<<endl;
        cout<<"n5_msz = "<<&n5_msz<<endl;
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
//    if (debug){PopDyInfo::debug=0; CatchInfo::debug=0;}
    if (debug) cout<<"finished PopProjector::project(dirF, n_msz)"<<endl;   
    RETURN_ARRAYS_DECREMENT();
    return 1.0*n5_msz;
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

    if (debug) {
        cout<<"matBio = "<<matBio<<endl;
        cout<<"finished PopProjector::projectUnFished(n_msz)"<<endl;   
    }
    RETURN_ARRAYS_DECREMENT();
    return 1.0*n3_msz;
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
    if (debug) {
        cout<<"n_msz  = "<<&n_msz<<endl;
        cout<<"np_msz = "<<&np_msz<<endl;
    }
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
    if (dtF<=dtM){ //fisheries occur BEFORE molting/growth/maturity 
        n2_msz.initialize();
        n3_msz.initialize();
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
/**
 * Class constructor.
 * 
 * @param pPPp - pointer to a PopProjector object
 */
MultiYearPopProjector::MultiYearPopProjector(PopProjector* pPPp){
    pPP = new PopProjector(*pPPp);
}
/**
 * Copy constructor.
 * 
 * @param o - MultiYearPopProjector object
 */
MultiYearPopProjector::MultiYearPopProjector(const MultiYearPopProjector& o){
    pPP = new PopProjector(*(o.pPP));
}
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
    n_ymsz.deallocate(); matBio_y.deallocate(); rm_yf.deallocate(); dm_yf.deallocate(); totCM_y.deallocate();
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
        if (debug) cout<<"y = "<<y<<endl;
        dvar3_array np_msz = pPP->project(dirF,n_ymsz(y-1),cout);
        if (debug) cout<<"n_ymsz(y-1) = "<<&(n_ymsz(y-1))<<endl<<"np_msz = "<<&np_msz<<endl;
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
        cout<<"matBio_y = "<<endl<<matBio_y<<endl;
        cout<<"totCM_y  = "<<endl<<totCM_y<<endl;
        cout<<"cp_yf    = "<<endl<<cp_yf<<endl;
        cout<<"rm_yf    = "<<endl<<rm_yf<<endl;
        cout<<"dm_yf    = "<<endl<<dm_yf<<endl;
    }
    if (debug) cout<<"Finished MultiYearPopProjector::project(n,R,n_msz)"<<endl;
}
/**
 * Project multiple years at constant recruitment and directed F.
 * 
 * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
 * 
 * @param n - number of years to project
 * @param R - single-sex recruitment vector
 * @param dirF - (constant) directed F
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 */
void MultiYearPopProjector::project(dvar_vector R, dvariable dirF, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"Starting MultiYearPopProjector::project(R,n_msz)"<<endl;
    if (debug) cout<<"nFsh = "<<pPP->nFsh<<endl;
    n_ymsz.deallocate(); matBio_y.deallocate(); rm_yf.deallocate(); dm_yf.deallocate(); totCM_y.deallocate();
    int n = R.size();
    n_ymsz.allocate(0,n,1,pPP->nMSs,1,pPP->nSCs,1,pPP->nZBs);
    matBio_y.allocate(1,n);
    totCM_y.allocate(1,n);
    int npFsh = pPP->nFsh;
    if (!npFsh) npFsh = 1;
    cp_yf.allocate(0,n,1,npFsh);
    rm_yf.allocate(0,n,1,npFsh);
    dm_yf.allocate(0,n,1,npFsh);
    n_ymsz.initialize();
    matBio_y.initialize();
    totCM_y.initialize();
    cp_yf.initialize();
    rm_yf.initialize();
    dm_yf.initialize();
    n_ymsz(0) = n_msz;
    if (debug){PopProjector::debug=1;PopDyInfo::debug=1;}
    for (int y=1;y<=n;y++){
        if (debug) cout<<"y = "<<y<<endl;
        dvar3_array np_msz = pPP->project(dirF,n_ymsz(y-1),cout);
        if (debug) cout<<"n_ymsz(y-1) = "<<&(n_ymsz(y-1))<<endl<<"np_msz = "<<&np_msz<<endl;
        matBio_y(y) = pPP->matBio;
        for (int f=1;f<=pPP->nFsh;f++) {
            cp_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->cpN_fmsz(f),cout);
            rm_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->rmN_fmsz(f),cout);
            dm_yf(y,f) = pPP->pPI->calcTotalBiomass(pPP->pCI->dmN_fmsz(f),cout);
        }
        totCM_y(y) = pPP->pPI->calcTotalBiomass(pPP->pCI->cmN_msz,cout);
        n_ymsz(y)  = pPP->addRecruitment(R[y],np_msz,cout);
    }
    if (debug){PopProjector::debug=0;PopDyInfo::debug=0;}
    if (debug){
        cout<<"matBio_y = "<<endl<<matBio_y<<endl;
        cout<<"totCM_y  = "<<endl<<totCM_y<<endl;
        cout<<"cp_yf    = "<<endl<<cp_yf<<endl;
        cout<<"rm_yf    = "<<endl<<rm_yf<<endl;
        cout<<"dm_yf    = "<<endl<<dm_yf<<endl;
    }
    if (debug) cout<<"Finished MultiYearPopProjector::project(R,n_msz)"<<endl;
}
/**
 * Project multiple years at constant recruitment and directed F.
 * 
 * @param n - number of years to project
 * @param R - (constant) single-sex recruitment
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 */
void MultiYearPopProjector::projectUnFished(int n, dvariable R, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"Starting MultiYearPopProjector::projectUnFished(n,R,n_msz)"<<endl;
    n_ymsz.deallocate(); matBio_y.deallocate(); rm_yf.deallocate(); dm_yf.deallocate(); totCM_y.deallocate();
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
        if (debug) cout<<"n_ymsz(y-1) = "<<&(n_ymsz(y-1))<<endl<<"np_msz = "<<&np_msz<<endl;
        matBio_y(y) = pPP->matBio;
        n_ymsz(y) = pPP->addRecruitment(R,np_msz,cout);
        if (debug) cout<<"n_ymsz(y) = "<<&(n_ymsz(y))<<endl;
    }
    if (debug){
        cout<<"matBio = "<<endl<<matBio_y<<endl;
        cout<<"Finished MultiYearPopProjector::projectUnFished(n,R,n_msz)"<<endl;
    }
}
/**
 * Project multiple years at constant recruitment and directed F.
 * 
 * @param R - single-sex recruitment vector (size is the number of years to project)
 * @param n_msz - initial abundance
 * @param cout - output stream for debug info
 * 
 */
void MultiYearPopProjector::projectUnFished(dvar_vector R, dvar3_array& n_msz, ostream& cout){
    if (debug) cout<<"Starting MultiYearPopProjector::projectUnFished(R,n_msz)"<<endl;
    n_ymsz.deallocate(); matBio_y.deallocate(); rm_yf.deallocate(); dm_yf.deallocate(); totCM_y.deallocate();
    int n = R.size();
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
        if (debug) cout<<"n_ymsz(y-1) = "<<&(n_ymsz(y-1))<<endl<<"np_msz = "<<&np_msz<<endl;
        matBio_y(y) = pPP->matBio;
        n_ymsz(y) = pPP->addRecruitment(R[y],np_msz,cout);
        if (debug) cout<<"n_ymsz(y) = "<<&(n_ymsz(y))<<endl;
    }
    if (debug){
        cout<<"matBio = "<<endl<<matBio_y<<endl;
        cout<<"Finished MultiYearPopProjector::projectUnFished(R,n_msz)"<<endl;
    }
}

