#include <limits>
#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelPopDyClasses.hpp"
#include "OFLCalcs.hpp"

using namespace tcsam;

/** flags to print debug info */
int Equilibrium_Calculator::debug = 0;
int Tier3_Calculator::debug = 0;
int OFL_Calculator::debug = 0;
int OFLResults::debug = 0;
////////////////////////////////////////////////////////////////////////////////
//Equilibrium_Calculator
////////////////////////////////////////////////////////////////////////////////
/**
 * Class constructor.
 * 
 * @param pPPp - pointer to a PopProjector object to base equilibrium calculations on
 */
Equilibrium_Calculator::Equilibrium_Calculator(PopProjector* pPPp){
    pPP = pPPp;
    nMSs = pPP->pPI->nMSs;
    nSCs = pPP->pPI->nSCs;
    nZBs = pPP->pPI->nZBs;
    I_z = identity_matrix(1,nZBs);
}
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
dvar3_array Equilibrium_Calculator::calcEqNatZ(dvar_vector& R_z, 
                                               dvar3_array& S1_msz, 
                                               dvar_matrix& Th_sz, 
                                               dvar3_array& T_szz, 
                                               dvar3_array& S2_msz, 
                                               ostream& cout){
    if (debug) cout<<"starting dvar3_array Equilibrium_Calculator::calcEqNatZ(...)"<<endl;
    if (debug) {
        cout<<"R_z = "<<R_z<<endl;
        cout<<"S1_msz = "<<endl; wts::print(S1_msz,cout,1); cout<<endl;
        cout<<"Th_sz = "<<endl;  wts::print(Th_sz,cout,1);  cout<<endl;
        cout<<"T_szz = "<<endl;  wts::print(T_szz,cout,1);  cout<<endl;
        cout<<"S2_msz = "<<endl; wts::print(S2_msz,cout,1); cout<<endl;
    }
    RETURN_ARRAYS_INCREMENT();

    //the equilibrium solution
    dvar3_array n_msz(1,nMSs,1,nSCs,1,nZBs);
        
    //--calc the state transition matrices
    int i = IMMATURE; 
    int m = MATURE;
    int n = NEW_SHELL;
    int o = OLD_SHELL;
    //immature new shell crab
    dvar_matrix S2_in = wts::diag(S2_msz(i,n)); //pr(survival|size) for immature new shell crab after molting occurs
    dvar_matrix Tr_in = T_szz(n);               //pr(Z_post|Z_pre, molt) for immature new shell crab pre-terminal molt
    dvar_matrix Th_in = wts::diag(Th_sz(n));    //pr(molt to maturity|post-molt size,new shell, molting)
    dvar_matrix Ph_in = identity_matrix(1,nZBs);//pr(molt|pre-molt size, new shell) [assumed that all immature crab molt]
    dvar_matrix S1_in = wts::diag(S1_msz(i,n)); //pr(survival|size) for immature new shell crab before molting occurs
    //immature old shell crab [shouldn't be any of these]
    dvar_matrix S2_io = wts::diag(S2_msz(i,o)); //pr(survival|size) for immature old shell crab after molting occurs (but they didn't molt)
    dvar_matrix Tr_io = T_szz(o);               //pr(Z_post|Z_pre, molt) for immature old shell crab pre-terminal molt
    dvar_matrix Th_io = wts::diag(Th_sz(o));    //pr(molt to maturity|post-molt size,old shell, molting)
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
    dvar_matrix lA = S2_in * (I_z-Th_in) * Tr_in * Ph_in * S1_in;//imm, new -> imm, new
    dvar_matrix lB = S2_in * (I_z-Th_io) * Tr_io * Ph_io * S1_io;//imm, old -> imm, new
    dvar_matrix lC = S2_io * (I_z-Ph_in) * S1_in;                //imm, new -> imm, old
    dvar_matrix lD = S2_io * (I_z-Ph_io) * S1_io;                //imm, old -> imm, old
    dvar_matrix lE = S2_mn * Th_in * Tr_mn * Ph_in * S1_in;      //imm, new -> mat, new (terminal molt)
    dvar_matrix lF = S2_mn * Th_io * Tr_mo * Ph_io * S1_io;      //imm, old -> mat, new (terminal molt)
    dvar_matrix lG = S2_mo * S1_mn;                              //mat, new -> mat, old
    dvar_matrix lH = S2_mo * S1_mo;                              //mat, old -> mat, old
    //--done calculating transition matrices
    
    //calculate inverses of matrix quantities
    dvar_matrix iM1 = inv(I_z - lD);
    dvar_matrix iM2 = inv(I_z - lA - lB * iM1 * lC);
    dvar_matrix iM3 = inv(I_z - lH);
    
    //the equilibrium solution is
    n_msz.initialize();
    n_msz(i,n) = iM2 * R_z;                         //immature, new shell
    n_msz(i,o) = iM1 * lC * n_msz(i,n);             //immature, old shell
    n_msz(m,n) = lE * n_msz(i,n) + lF * n_msz(i,o); //  mature, new shell
    n_msz(m,o) = iM3 * lG * n_msz(m,n);             //  mature, old shell
        
    if (debug) cout<<"finished dvar3_array Equilibrium_Calculator::calcEqNatZ(...)"<<endl;
    RETURN_ARRAYS_DECREMENT();
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
dvar3_array Equilibrium_Calculator::calcEqNatZF0(dvariable R, ostream& cout){
    if (debug) {
        cout<<"starting dvar3_array Equilibrium_Calculator::calcEqNatZF0(double R)"<<endl;
        cout<<"R = "<<R<<endl;
    }
    RETURN_ARRAYS_INCREMENT();

    dvar3_array S1_msz = pPP->pPI->calcSurvival(pPP->dtM,cout);
    dvar3_array S2_msz = pPP->pPI->calcSurvival(1.0-pPP->dtM,cout);
    if (debug){
        cout<<"#--S1_msz = "<<&S1_msz<<endl;
        cout<<"#--S2_msz = "<<&S2_msz<<endl;
    }
    
    dvar_vector R_zp = R*pPP->pPI->R_z;
    
    dvar3_array n_msz = calcEqNatZ(R_zp, S1_msz, pPP->pPI->Th_sz, pPP->pPI->T_szz, S2_msz,cout);
    
    if (debug) {
        cout<<"Equilibrium_Calculator::calcEqNatZF0: n_msz = "<<endl; wts::print(n_msz,cout,1); cout<<endl;
        cout<<"finished dvar3_array Equilibrium_Calculator::calcEqNatZF0(double R)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
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
dvar3_array Equilibrium_Calculator::calcEqNatZFM(dvariable R, dvariable dirF, ostream& cout){
    if (debug) {
        cout<<"starting dvar3_array Equilibrium_Calculator::calcEqNatZFM(double R)"<<endl;
        cout<<"R = "<<R<<". dirF = "<<dirF<<endl;
    }
    RETURN_ARRAYS_INCREMENT();
    
    double dtM = pPP->dtM;//time to mating
    double dtF = pPP->dtF;//time to fisheries
    
    dvar3_array S1_msz(1,nMSs,1,nSCs,1,nZBs);//survival until molting/mating
    dvar3_array S2_msz(1,nMSs,1,nSCs,1,nZBs);//survival after molting/mating
    dvar3_array n_msz(1,nMSs,1,nSCs,1,nZBs); //equilibrium size distribution on July 1
    
    if (dtF<=dtM){
        //fisheries occur BEFORE molting/growth/maturity 
        dvar3_array S1a_msz = pPP->pPI->calcSurvival(dtF,cout);     //survival prior to fisheries
        dvar3_array S1F_msz = pPP->pCI->calcSurvival(dirF,cout);    //survival of fisheries
        dvar3_array S1b_msz = pPP->pPI->calcSurvival(dtM-dtF,cout); //survival after fisheries, before mating
        if (debug){
            cout<<"#--S1a_msz = "<<&S1a_msz<<endl;
            cout<<"#--S1F_msz = "<<&S1F_msz<<endl;
            cout<<"#--S1b_msz = "<<&S1b_msz<<endl;
        }
        S1_msz = elem_prod(S1b_msz,elem_prod(S1F_msz,S1a_msz));//total survival before mating
        S2_msz = pPP->pPI->calcSurvival(1.0-dtM,cout);//survival after mating/molting/growth
    } else {
        //fisheries occur AFTER molting/growth/maturity 
        S1_msz = pPP->pPI->calcSurvival(dtM,cout);                   //survival before mating/molting/growth
        dvar3_array S2a_msz = pPP->pPI->calcSurvival(dtF-dtM,cout);     //survival afterMGM, before fisheries
        dvar3_array S2F_msz = pPP->pCI->calcSurvival(dirF,cout);        //survival of fisheries
        dvar3_array S2b_msz = pPP->pPI->calcSurvival(1.0-dtF,cout);     //survival after fisheries, to year end
        if (debug){
            cout<<"#--S2a_msz = "<<&S2a_msz<<endl;
            cout<<"#--S2F_msz = "<<&S2F_msz<<endl;
            cout<<"#--S2b_msz = "<<&S2b_msz<<endl;
        }
        S2_msz = elem_prod(S2b_msz,elem_prod(S2F_msz,S2a_msz)); //total survival after MGM to year end
    }
    dvar_vector R_zp = R*pPP->pPI->R_z;
    n_msz = calcEqNatZ(R_zp, S1_msz, pPP->pPI->Th_sz, pPP->pPI->T_szz, S2_msz, cout);
    
    if (debug) {
        cout<<"n_msz = "<<endl; wts::print(n_msz,cout,1);
        cout<<"finished dvar3_array Equilibrium_Calculator::calcEqNatZFM(double R)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
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
dvariable Equilibrium_Calculator::calcEqMatureBiomassAtMatingF0(dvariable R, ostream& cout){
    if (debug) {
        cout<<"starting double Equilibrium_Calculator::calcEqMatureBiomassAtMatingF0(double R)"<<endl;
        cout<<"R = "<<R<<endl;
    }
    RETURN_ARRAYS_INCREMENT();
    //calculate July 1 unfished size distribution
    dvar3_array n_msz = calcEqNatZF0(R,cout);
    
    //advance to mating in unfished population
    dvar3_array nm_msz = pPP->pPI->applyNM(pPP->dtM,n_msz,cout);
    
    //calculate mature biomass at time of mating
    dvariable eqMB = pPP->pPI->calcMatureBiomass(nm_msz,cout);
    if (debug) {
        cout<<"eqMB = "<<eqMB<<endl;
        cout<<"finished double Equilibrium_Calculator::calcEqMatureBiomassAtMatingF0(double R)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
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
dvariable Equilibrium_Calculator::calcEqMatureBiomassAtMatingFM(dvariable R, dvariable dirF, ostream& cout){
    if (debug) {
        cout<<"starting Equilibrium_Calculator::calcEqMatureBiomassAtMatingFM(double R, double dirF)"<<endl;
        cout<<"R = "<<R<<". dirF = "<<dirF<<endl;
    }
    RETURN_ARRAYS_INCREMENT();
    //calculate equilibrium size distribution on July 1
    dvar3_array n_msz = calcEqNatZFM(R,dirF,cout);
    
    //equilibrium MMB
    dvariable eqMB = 0;//dummy value    
    //advance to population to time of mating
    if (pPP->dtF<=pPP->dtM){ //fisheries occur BEFORE molting/growth/maturity 
        if (debug) {cout<<"dtF<=dtM"<<endl;}
        //apply natural mortality BEFORE fisheries
        dvar3_array n1_msz = pPP->pPI->applyNM(pPP->dtF,n_msz,cout);
        if (debug) {cout<<"n1_msz ="<<endl; wts::print(n1_msz,cout,1);}
        //apply fisheries
        dvar3_array n2_msz = pPP->pCI->applyFM(dirF, n1_msz, cout);
        if (debug) {cout<<"n2_msz ="<<endl; wts::print(n2_msz,cout,1);}
        //apply natural mortality after fisheries but before molting/growth
        dvar3_array n3_msz(1,nMSs,1,nSCs,1,nZBs);
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
        dvar3_array n1_msz = pPP->pPI->applyNM(pPP->dtM,n_msz,cout);
        if (debug) {cout<<"n1_msz ="<<endl; wts::print(n1_msz,cout,1);}
        
        //calculate mature biomass at mating
        eqMB = pPP->pPI->calcMatureBiomass(n1_msz,cout);
    }
    if (debug) {
        cout<<"eqMB = "<<eqMB<<endl;
        cout<<"finished Equilibrium_Calculator::calcEqMatureBiomassAtMatingFM(dvariable R, dvariable dirF)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return eqMB;
}

////////////////////////////////////////////////////////////////////////////////
//Tier3_Calculator
////////////////////////////////////////////////////////////////////////////////
/** flag to print debug info */
/**
 * Constructor.
 * 
 * @params XX - target SPR rate for MSY calculations
 * @param pECp - pointer to an Equilibrium_Calculator object
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
dvariable Tier3_Calculator::calcB100(dvariable R, ostream& cout){
    if (debug) cout<<"starting dvariable Tier3_Calculator::calcB100(dvariable R)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    //calculate mature biomass at time of mating
    B100 = pEC->calcEqMatureBiomassAtMatingF0(R,cout);
    B0 = B100;
    if (debug) {
        cout<<"B100 = "<<B100<<endl;
        cout<<"finished dvariable Tier3_Calculator::calcB100(dvariable R)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return 1.0*B100;
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
dvariable Tier3_Calculator::calcBmsy(dvariable R, ostream& cout){
    if (debug) cout<<"starting double Tier3_Calculator::calcBmsy(dvariable R)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    Bmsy = XX*calcB100(R,cout);
    if (debug) {
        cout<<"Bmsy = "<<Bmsy<<endl;
        cout<<"finished double Tier3_Calculator::calcBmsy(dvariable R)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return 1.0*Bmsy;
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
dvariable Tier3_Calculator::calcFmsy(dvariable R, ostream& cout){
    if (debug) cout<<"starting Tier3_Calculator::calcFmsy(R)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    //calculate unfished mmb (B0)
    B0   = calcB100(R,cout);
    Bmsy = XX*B0;
    if (debug){cout<<"in Tier3_Calculator::calcFmsy(R): B100 ="<<B0<<". Bmsy = "<<Bmsy<<endl;}
    
    //From initial guess for FXX, iterate to improve FXX    
    dvariable dF   = 0.0001;//"delta" F for derivative calculation
    dvariable FXX  = 0.5;   //initial guess for FXX
    dvariable dFXX = 1.0;   //"delta" FXX to update (initial value is a dummy)
    dvariable mmbp, mmb, mmbm, dMMBdF, XXp;
    int i=0;
    while ((i++ < maxIts)){
        mmbp = pEC->calcEqMatureBiomassAtMatingFM(R,FXX+dF,cout);
        mmb  = pEC->calcEqMatureBiomassAtMatingFM(R,FXX,cout);
        mmbm = pEC->calcEqMatureBiomassAtMatingFM(R,FXX-dF,cout);
        dMMBdF = 0.5*(mmbp-mmbm)/dF;//derivative of mmb wrto F
        XXp   = mmb/B0;           //ratio of mmb for current FXX relative to unfished 
        dFXX  = (Bmsy - mmb)/dMMBdF;
        FXX  += dFXX;
        if (debug&&(i==maxIts)) {
            cout<<"--max iteration = "<<i<<endl;
            cout<<"----mmb  = "<<mmb<<". Bmsy = "<<Bmsy<<endl;
            cout<<"----XXp  = "<<XXp<<". dFXX = "<<dFXX<<endl;
            cout<<"----FXX  = "<<FXX<<endl;
        }
    }//i loop
    
    if (sfabs(value(dFXX))>0.00001){
        cout<<endl<<"-------ERROR!!-------"<<endl;
        cout<<"Convergence in calFmsy failed: |dfXX| = "<<sfabs(value(dFXX))<<endl;
        cout<<"-------ERROR!!-------"<<endl<<endl;
    }
    
    if (debug) {
        cout<<"Fmsy = "<<FXX<<endl;
        cout<<"finished Tier3_Calculator::calcFmsy(R)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return FXX;//Fmsy for Tier 3 stock
}

////////////////////////////////////////////////////////////////////////////////
//OFL_Calculator
////////////////////////////////////////////////////////////////////////////////
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
dvariable OFL_Calculator::calcHCR(dvariable currMMB, dvariable Bmsy, dvariable Fmsy, ostream& cout){
    if (debug) cout<<"starting OFL_Calculator::calcHCR(currMMB, Bmsy, Fmsy)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    dvariable Fofl = 0.0;
    double ratio  = value(currMMB/Bmsy);
    if (ratio < beta){Fofl = 0.0;} else
    if (ratio < 1.0){
        Fofl = Fmsy*(ratio-alpha)/(1-alpha);
    } else Fofl = Fmsy;
    
    if (debug) {
        cout<<"Fofl = "<<Fofl<<endl;
        cout<<"finished OFL_Calculator::calcHCR(currMMB, Bmsy, Fmsy)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
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
dvariable OFL_Calculator::calcFofl(dvariable Bmsy, dvariable Fmsy, dvar3_array& n_msz, ostream& cout){
    if (debug) {
        cout<<"starting double OFL_Calculator::calcFofl(Bmsy, Fmsy,n_msz)"<<endl;
        cout<<"Bmsy = "<<Bmsy<<"; Fmsy = "<<Fmsy<<endl;
        cout<<"n_msz = "<<endl;wts::print(n_msz,cout,1); cout<<endl;        
    }
    RETURN_ARRAYS_INCREMENT();
    
    PopProjector* pPrjM = pTCM->pEC->pPP;
    
    //start with guess for Fofl based on currMMB
    dvariable currMMB = pPrjM->projectMatureBiomassAtMating(Fmsy,n_msz,cout);
    if (debug) cout<<"init currMMB = "<<currMMB<<"; B/Bmsy = "<<currMMB/Bmsy<<endl;
    dvariable Fofl = calcHCR(currMMB,Bmsy,Fmsy,cout);
    if (debug) cout<<"init Fofl = "<<Fofl<<endl;
    //now iterate until Fofl yields currMMB
    dvariable Foflp = 0.0;
    int i = 1;
    while (i++ < maxIts){
        if (debug&&(i==maxIts)) cout<<"HCR iteration for Fofl = "<<i<<endl;
        //calculate currMMB based on Fofl
        currMMB = pPrjM->projectMatureBiomassAtMating(Fofl,n_msz,cout);
        if (debug) cout<<"--updated prjMMB = "<<currMMB<<"; B/Bmsy = "<<currMMB/Bmsy<<endl;
        //update Fofl based on currMMB
        Foflp = Fofl;
        Fofl  = calcHCR(currMMB,Bmsy,Fmsy,cout);
        if (debug) cout<<"--updated Fofl = "<<Fofl<<"; delF = "<<Fofl - Foflp<<endl;
    }
    double criF = 0.001;
    if (sfabs(value(Fofl - Foflp))>criF){
        cout<<endl<<"-------ERROR!!-------"<<endl;
        cout<<"Convergence in calcFofl failed: |delF| = "<<sfabs(value(Fofl - Foflp))<<" > "<<criF<<endl;
        cout<<"-------ERROR!!-------"<<endl<<endl;
    }
    prjMMB = currMMB;
    if (debug) {
        cout<<"Fofl   = "<<Fofl<<endl;
        cout<<"prjMMB = "<<prjMMB<<endl;
        cout<<"finished double OFL_Calculator::calcFofl(Bmsy, Fmsy,n_msz)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
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
dvariable OFL_Calculator::calcMSY(dvar_vector R_x, dvariable Fmsy, ostream& cout){
    if (debug) {
        cout<<"starting double OFL_Calculator::calcMSY(R_x,Fmsy)"<<endl;
        cout<<"R_x = "<<R_x<<tb<<"Fmsy = "<<Fmsy<<endl;
    }
    RETURN_ARRAYS_INCREMENT();
    int nFsh = pTCM->pEC->pPP->nFsh;
    if (debug) cout<<"nFsh = "<<nFsh<<endl;
    
    msy_fx.allocate(0,nFsh,1,tcsam::nSXs);
    msy_fx.initialize();
    
    dvariable totCM;
        
    //do males
    if (debug) cout<<"Calculating MALE portion of MSY"<<endl;
    Equilibrium_Calculator* pECM = pTCM->pEC;
    PopProjector*           pPPM = pTCM->pEC->pPP;
    //calculate equilibrium numbers-at-size for males
    dvar3_array n_msz = pECM->calcEqNatZFM(R_x(MALE),Fmsy,cout);
    //calculate associated catch
    pPPM->project(Fmsy,n_msz,cout);
    //calc retained catch biomass (directed fishery only) from catch abundance for males
    msy_fx(0,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->rmN_fmsz(1),cout);
    //calc discard mortality biomass from catch abundance for males
    for (int f=1;f<=nFsh;f++){
        msy_fx(f,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->dmN_fmsz(f),cout);
    }
    totCM = pPPM->pPI->calcTotalBiomass(pPPM->pCI->cmN_msz,cout);    
    
    if (tcsam::nSXs>1){
        //do females
        if (debug) cout<<"Calculating FEMALE portion of MSY"<<endl;
        Equilibrium_Calculator* pECF = pTCF->pEC;
        PopProjector*           pPPF = pTCF->pEC->pPP;
        //calculate equilibrium numbers-at-size for females
        dvar3_array n_msz = pECF->calcEqNatZFM(R_x(FEMALE),Fmsy,cout);
        //calculate associated catch
        pPPF->project(Fmsy,n_msz,cout);
        //calc retained catch biomass (directed fishery only) from catch abundance for females
        msy_fx(0,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->rmN_fmsz(1),cout);
        //calc discard mortality biomass from catch abundance for females
        for (int f=1;f<=nFsh;f++){
            msy_fx(f,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->dmN_fmsz(f),cout);
        }
        totCM += pPPF->pPI->calcTotalBiomass(pPPF->pCI->cmN_msz,cout);    
    }
    
    dvariable msy = sum(msy_fx);
    if (debug||abs(value(0.5*(msy-totCM)/(msy+totCM)))>0.01) {
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"in double OFL_Calculator::calcMSY(R_x,Fmsy,n_xmsz)"<<endl;
        cout<<"MSY   = "<<msy<<endl;
        cout<<"totCM = "<<totCM<<endl;
        cout<<"msy_fx: "<<endl<<msy_fx<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"finished double OFL_Calculator::calcMSY(R_x,Fmsy)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
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
dvariable OFL_Calculator::calcOFL(dvariable Fofl, dvar4_array& n_xmsz, ostream& cout){
    if (debug) cout<<"starting double OFL_Calculator::calcOFL(Fofl,n_xmsz)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    
    int nFsh = pTCM->pEC->pPP->nFsh;
    if (debug) cout<<"nFsh = "<<nFsh<<endl;
    
    ofl_fx.allocate(0,nFsh,1,tcsam::nSXs);
    ofl_fx.initialize();

    dvariable totCM;
    
    //do males
    if (debug) cout<<"calculating OFL portion for MALEs"<<endl;
    PopProjector* pPPM = pTCM->pEC->pPP;
    //calc catch abundance of males at Fofl
    pPPM->project(Fofl,n_xmsz(  MALE),cout);
    //calc retained catch biomass (directed fishery only) from catch abundance for males
    ofl_fx(0,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->rmN_fmsz(1),cout);
    //calc discard mortality biomass from catch abundance for males
    for (int f=1;f<=nFsh;f++){
        ofl_fx(f,MALE) = pPPM->pPI->calcTotalBiomass(pPPM->pCI->dmN_fmsz(f),cout);
    }
    totCM = pPPM->pPI->calcTotalBiomass(pPPM->pCI->cmN_msz,cout);    
    
    if (tcsam::nSXs>1){
        //do females
        if (debug) cout<<"calculating OFL portion for FEMALEs"<<endl;
        PopProjector* pPPF = pTCF->pEC->pPP;
        //calc catch abundance of females at Fofl
        pPPF->project(Fofl,n_xmsz(FEMALE),cout);
        //calc retained catch biomass (directed fishery only) from catch abundance for females
        ofl_fx(0,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->rmN_fmsz(1),cout);
        //calc discard mortality biomass from catch abundance for females
        for (int f=1;f<=nFsh;f++){
            ofl_fx(f,FEMALE) = pPPF->pPI->calcTotalBiomass(pPPF->pCI->dmN_fmsz(f),cout);
        }
        totCM += pPPF->pPI->calcTotalBiomass(pPPF->pCI->cmN_msz,cout);    
    }
    
    dvariable ofl = sum(ofl_fx);
    if (debug||(ofl!=totCM)) {
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"in double OFL_Calculator::calcOFL(Fofl,n_xmsz)"<<endl;
        cout<<"OFL   = "<<ofl<<endl;
        cout<<"totCM = "<<totCM<<endl;
        cout<<"ofl_fx: "<<endl<<ofl_fx<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"finished double OFL_Calculator::calcOFL(Fofl,n_xmsz)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
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
dvariable OFL_Calculator::calcPrjMMB(dvariable Fofl, dvar3_array& n_msz, ostream& cout){
    if (debug) {
        cout<<"starting dvariable OFL_Calculator::calcPrjMMB(Fofl,n_msz)"<<endl;
        cout<<"Fofl = "<<Fofl<<endl;
    }
    RETURN_ARRAYS_INCREMENT();
    PopProjector* pPrjM = pTCM->pEC->pPP;
    dvariable mmb = pPrjM->projectMatureBiomassAtMating(Fofl,n_msz,cout);
    if (debug) {
        cout<<"mmb = "<<mmb<<endl;
        cout<<"finished dvariable OFL_Calculator::calcPrjMMB(Fofl,n_msz)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return mmb;
}

/**
 * Calculate all results associated with the OFL.
 * 
 * @param R - assumed sex-specific equilibrium recruitment
 * @param n_xmsz - d4_array with initial population abundance
 * @param cout - output stream for debug info
 * 
 * @return pointer to OFLResults object.
 */
OFLResults* OFL_Calculator::calcOFLResults(dvar_vector R, dvar4_array& n_xmsz, ostream& cout){
    if (debug) cout<<"starting OFLResults OFL_Calculator::calcOFLResults(R,n_xmsz,cout)"<<endl;
    RETURN_ARRAYS_INCREMENT();
    int nFsh = pTCM->pEC->pPP->pCI->nFsh;
    OFLResults* res = new OFLResults();
    
    res->avgRec_x = R;
    res->finlNatZ_xmsz.allocate(1,tcsam::nSXs,
                                1,tcsam::nMSs,
                                1,tcsam::nSCs,
                                1,pTCM->pEC->pPP->nZBs);
    res->finlNatZ_xmsz = n_xmsz;
    res->pPDIM = pTCM->pEC->pPP->pPI;
    res->pCIM  = pTCM->pEC->pPP->pCI;
    if (tcsam::nSXs>1){
        res->pPDIF = pTCF->pEC->pPP->pPI;
        res->pCIF  = pTCF->pEC->pPP->pCI;
    }
    res->curB     = pTCM->pEC->pPP->pPI->calcMatureBiomass(n_xmsz(MALE),cout);
    if (debug) cout<<"calcOFLResults: calculated curB"<<endl;
        
    res->eqNatZF0_xmsz.allocate(1,tcsam::nSXs,
                                1,tcsam::nMSs,
                                1,tcsam::nSCs,
                                1,pTCM->pEC->pPP->nZBs);
    res->eqNatZF0_xmsz(MALE) = pTCM->pEC->calcEqNatZF0(R(MALE),cout);
    if (debug) cout<<"calcOFLResults: calculated eq NatZ(MALE) for F=0"<<endl;
    if (tcsam::nSXs>1) {
        res->eqNatZF0_xmsz(FEMALE) = pTCF->pEC->calcEqNatZF0(R(FEMALE),cout);
        if (debug) cout<<"calcOFLResults: calculated eq NatZ(FEMALE) for F=0"<<endl;
    }
            
    res->Fmsy = pTCM->calcFmsy(R(MALE),cout);//also calculates B0 and Bmsy
    res->B0   = pTCM->B0;
    res->Bmsy = pTCM->Bmsy;
    res->MSY  = calcMSY(R,res->Fmsy,cout);
    if (debug) cout<<"calcOFLResults: calculated MSY"<<endl;
    
    res->eqNatZFM_xmsz.allocate(1,tcsam::nSXs,
                                1,tcsam::nMSs,
                                1,tcsam::nSCs,
                                1,pTCM->pEC->pPP->nZBs);
    if (debug) {
        cout<<"calcOFLResults: allocated eq NatZ for F=Fmsy"<<endl;
    }
    res->eqNatZFM_xmsz(MALE) = pTCM->pEC->calcEqNatZFM(R(MALE),res->Fmsy,cout);
    if (debug) cout<<"calcOFLResults: calculated eq NatZ(MALE) for F=Fmsy"<<endl;
    if (tcsam::nSXs>1) {
        res->eqNatZFM_xmsz(FEMALE) = pTCF->pEC->calcEqNatZFM(R(FEMALE),res->Fmsy,cout);
        if (debug) cout<<"calcOFLResults: calculated eq NatZ(FEMALE) for F=Fmsy"<<endl;
    }
            
    res->Fofl = calcFofl(res->Bmsy,res->Fmsy,n_xmsz(MALE),cout);
    if (debug) cout<<"calcOFLResults: calculated Fofl"<<endl;
    res->prjB = prjMMB;
    res->OFL  = calcOFL(res->Fofl,n_xmsz,cout);
    if (debug) cout<<"calcOFLResults: calculated OFL"<<endl;
    
    res->ofl_fx.allocate(0,nFsh,1,tcsam::nSXs);
    res->ofl_fx = ofl_fx;
    RETURN_ARRAYS_DECREMENT();
    if (debug) cout<<"finished OFLResults* OFL_Calculator::calcOFLResults(R,n_xmsz,cout)"<<endl;
    return res;
}
////////////////////////////////////////////////////////////////////////////////
//OFLResults
////////////////////////////////////////////////////////////////////////////////
OFLResults::OFLResults(){
    pPDIM=0;
    pPDIF=0;
    pCIM=0;
    pCIF=0;
}

OFLResults::~OFLResults(){
    std::cout<<"Deleting OFLResults object."<<endl;
    if (pPDIM) delete(pPDIM); pPDIM=0;
    if (pPDIF) delete(pPDIF); pPDIF=0;
    if (pCIM)  delete(pCIM);  pCIM=0;
    if (pCIF)  delete(pCIF);  pCIF=0;
    std::cout<<"Deleted OFLResults object."<<endl;
}

/**
 * Assignment operator for OFLResults class.
 * 
 * @param OFLResults object to copy.
 * @return reference to the copied object
 */
OFLResults& OFLResults::operator=(const OFLResults& o){
    if (debug) std::cout<<"starting OFLResults::operator=(const OFLResults&)"<<endl;
    avgRec_x = 1.0*o.avgRec_x;//average recruitment, by sex
    B0       = 1.0*o.B0;      //equilibrium MMB for unfished population
    Fmsy     = 1.0*o.Fmsy;    //equilibrium F on directed fishery for males resulting in MSY
    Bmsy     = 1.0*o.Bmsy;    //equilibrium MMB when fished at Fmsy 
    MSY      = 1.0*o.MSY;     //equilibrium yield (1000's t) when fished at Fmsy
    Fofl     = 1.0*o.Fofl;    //F on directed fishery for males resulting in the OFL
    OFL      = 1.0*o.OFL;     //total OFL (1000's t)
    ofl_fx   = 1.0*o.ofl_fx;  //fishery/sex-specific mortality components to OFL (f=0 is retained catch, f>0 is total catch mortality)
    prjB     = 1.0*o.prjB;    //projected MMB for projection year when current population is fished at Fofl.
    curB     = 1.0*o.curB;    //"current" MMB at beginning of projection year
    finlNatZ_xmsz.deallocate(); //final pop state in assessment model
    finlNatZ_xmsz.allocate(o.finlNatZ_xmsz);
    if (debug) std::cout<<"got here 0"<<endl;
    eqNatZF0_xmsz.deallocate(); //unfished equilibrium size distribution
    eqNatZF0_xmsz.allocate(o.eqNatZF0_xmsz);
    if (debug) std::cout<<"got here 1"<<endl;
    eqNatZF0_xmsz = o.eqNatZF0_xmsz;
    eqNatZFM_xmsz.deallocate(); //unfished equilibrium size distribution
    eqNatZFM_xmsz.allocate(o.eqNatZF0_xmsz);
    if (debug) std::cout<<"got here 2"<<endl;
    for (int x=1;x<=tcsam::nSXs;x++) eqNatZFM_xmsz(x) = 1.0*o.eqNatZF0_xmsz(x);
    
    pPDIM = o.pPDIM;
    pPDIF = o.pPDIF;
    pCIM  = o.pCIM;
    pCIF  = o.pCIF;
    
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
    adstring fDms = "f=c('directed fishery',"+wts::to_qcsv(ptrMC->lblsFsh)+")";//fisheries
    os<<name<<"=list(OFL="<<OFL<<cc<<"Fofl="<<Fofl<<cc<<"prjB="<<prjB<<cc<<"curB="<<curB;
    os<<cc<<"Fmsy="<<Fmsy<<cc<<"Bmsy="<<Bmsy<<cc<<"MSY="<<MSY;
    os<<cc<<"B100="<<B0<<cc<<"avgRec="<<sum(avgRec_x)<<cc<<"avgRecM="<<avgRec_x(MALE);
    if (tcsam::nSXs>1) os<<cc<<"avgRecF="<<avgRec_x(FEMALE);
    os<<cc<<endl;
    os<<"ofl_fx="; wts::writeToR(os,value(ofl_fx),fDms,xDms); os<<cc<<endl;
    os<<"finlNatZ_xmsz="; wts::writeToR(os,finlNatZ_xmsz,xDms,mDms,sDms,zDms); os<<cc<<endl;
    pPDIM->writeToR(os,ptrMC,"popDyInfoM",0); os<<cc<<endl;
    if (tcsam::nSXs>1) {pPDIF->writeToR(os,ptrMC,"popDyInfoF",0); os<<cc<<endl;}
    pCIM->writeToR(os,ptrMC,"catchInfoM",0); os<<cc<<endl;
    if (tcsam::nSXs>1) {pCIF->writeToR(os,ptrMC,"catchInfoF",0); os<<cc<<endl;}
    os<<"eqNatZF0_xmsz="; wts::writeToR(os,eqNatZF0_xmsz,xDms,mDms,sDms,zDms); os<<cc<<endl; 
    os<<"eqNatZFM_xmsz="; wts::writeToR(os,eqNatZFM_xmsz,xDms,mDms,sDms,zDms); os<<endl;
    os<<")";
}
