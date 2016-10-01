/* 
 * File:   OFLCalcs.hpp
 * Author: WilliamStockhausen
 *
 * Created on March 4, 2016, 10:16 AM
 */

#ifndef OFLCALCS_HPP
#define	OFLCALCS_HPP

/**
 * Class to facilitate single-sex population dynamics calculations.
 */
class PopDyInfo {
    public:
        static int debug;
        int nMSs;//number of maturity states
        int nSCs;//number of shell conditions
        int nZBs;//number of size bins
        
    public:
        dvector  R_z;      //recruitment size distribution
        dmatrix  w_mz;     //weight-at-size
        d3_array M_msz;    //natural mortality
        dmatrix  Th_sz;    //pr(molt to maturity|pre-molt size, molt)
        d3_array T_szz;    //growth matrices (indep. of molt to maturity)
        
    public:
        /**
         * Class constructor.
         * 
         * @param npZBs - number of size bins
         */
        PopDyInfo(int npZBs);
        /**
         * Class destructor.
         */
        ~PopDyInfo(){}
        /**
         * Calculate single-sex mature biomass, given population abundance.
         * 
         * @param n_msz - single-sex population abundance
         * @param cout - stream for writing debug info
         * 
         * @return mature biomass (units??)
         */
        double calcMatureBiomass(d3_array& n_msz, ostream& cout);
        /**
         * Calculate total single-sex biomass, given population abundance.
         * 
         * @param n_msz - single-sex population abundance
         * @param cout - stream for writing debug info
         * 
         * @return total biomass (units??)
         */
        double calcTotalBiomass(d3_array& n_msz, ostream& cout);
        /**
         * Calculate single-sex survival probabilities, given natural mortality rates,
         * over a given period of time (dt)
         * 
         * @param dt - period of time (in years)
         * @param cout - stream for writing debug info
         * 
         * @return d3_array S_msz
         */
        d3_array calcSurvival(double dt, ostream& cout);
        /**
         * Calculate effects of natural mortality on population abundance.
         * 
         * @param dt - time interval (fraction of year)
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final population abundance (after applying natural mortality)
         */
        d3_array applyNM(double dt, d3_array& n_msz, ostream& cout);
        /**
         * Calculate effects of molting/growth on population abundance.
         * 
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final population abundance (after applying molting/growth)
         */
        d3_array applyMG(d3_array& n_msz, ostream& cout);
};//PopDyInfo

/**
 * Class to facilitate fishery catch calculations for one sex.
 */
class CatchInfo {
    public:
        static int debug;
        int nMSs;//number of maturity states
        int nSCs;//number of shell conditions
        int nZBs;//number of size bins
        int nFsh;//number of fisheries
        
    public:
        double   maxF;     //maximum capture rate in target fishery
        d4_array capF_fmsz;//fishery capture rates
        d4_array retF_fmsz;//retention function
        dvector  hm_f;     //discard mortality rates
        
    public:
        d3_array ct_msz; //catch mortality (abundance)
        d4_array rm_fmsz;//retained mortality (abundance)
        d4_array dm_fmsz;//discard mortality (abundance)
        
    public:
        /**
         * Class constructor.
         * 
         * @param npZBs - number of size bins
         * @param npFsh - number of fisheries
         */
        CatchInfo(int npZBs, int npFsh);
        /**
         * Class destructor.
         */
        ~CatchInfo(){}
        /**
         * Calculate maximum fishing capture rate for target fishery (f=1)
         * 
         * Modifies: maxF
         * 
         * @param cout - stream for writing debug info
         * 
         * @return maxF
         */
        double findMaxTargetCaptureRate(ostream& cout);
        /**
         * Calculate catch abundance (ct_msz, rm_fmsz, dm_fmsz) and post-fisheries
         * abundance based on target fishery capture rate 'dirF' and initial population 
         * abundance n_msz.
         * 
         * Modifies: ct_msz, rm_fmsz, dm_fmsz
         * 
         * @param dirF - fully-selected directed fishery capture rate
         * @param n_msz - pre-fisheries population size
         * @param cout - stream for writing debug info
         * 
         * @return np_msz - post-fisheries population size (d3_array)
         * 
         */
        d3_array applyFM(double dirF, d3_array& n_msz, ostream& cout);
        /**
         * Calculates probabilities of surviving fisheries, given directed 
         * fishing capture rate 'dirF'.
         * 
         * @param dirF - fully-selected capture rate for target fishery (f=1)
         * @param cout - stream for writing debug info
         * 
         * @return survival probabilities S_msz as d3_array.
         */
        d3_array calcSurvival(double dirF, ostream& cout);
        /**
         * Set sex-specific fishery capture rates.
         * 
         * @param x - sex to select
         * @param capF_fxmsz - input capture rates
         */
        void setCaptureRates(int x, d5_array& capF_fxmsz);
        /**
         * Set sex-specific retention functions.
         * 
         * @param x - sex to select
         * @param retF_fxmsz - input retention functions
         */
        void setRetentionFcns(int x, d5_array& retF_fxmsz);
        /**
         * Set handling mortality rates, by fishery.
         * 
         * @param pHM_f - dvector of handling mortality rates
         */
        void setHandlingMortality(dvector& pHM_f);
};//CatchInfo

class Tier3_Calculator {
    public:
        static int debug;
        int nMSs;//number of maturity states
        int nSCs;//number of shell conditions
        int nZBs;//number of size bins
        int nFsh;//number of fisheries
        double dtF; //time to fisheries
        double dtM; //time to molting/growth
        
    public:
        PopDyInfo* pPI;//pointer to PopDyInfo object for males
        CatchInfo* pCI;//pointer to CatchInfo object for males
    public:
        double XX;  //SPR level for Fmsy, Bmsy
        double B100;//unfished MMB
        double Fmsy;//=FXX; directed fishery capture F resulting in MSY
        double Bmsy;//=XX*B100; longterm B resulting from directed fishing at F
        
    public:
        /**
         * Class constructor.
         * 
         * @param npZBs - number of size bins
         * @param npFsh - number of fisheries
         */
        Tier3_Calculator(int npZBs, int npFsh);
        /**
         * Class destructor.
         */
        ~Tier3_Calculator(){delete pPI; delete pCI;}
        /**
         * Calculate B100, which is equilibrium mature male biomass (MMB)
         * for an unfished stock.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return MMB for unfished population (1000's t)
         */
        double calcB100(double R, ostream& cout);
        /**
         * Calculate equilibrium spawning biomass (as MMB) when directed fishery
         * is fished at longterm capture rate Fmsy. For a Tier 3 stock, 
         * Bmsy = B35% = 0.35*B100.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return equilibrium MMB for population fished at Fmsy (1000's t)
         */
        double calcBmsy(double R, ostream& cout);
        /**
         * Calculate directed fishery capture rate, Fmsy, that results in 
         * equilibrium (longterm) MMB = Bmsy. For a Tier 3 stock, Fmsy and 
         * Bmsy are given by SPR-based proxies F35% nd B35%, where F35% is
         * the directed fishery capture rate that results in equilibrium
         * MMB = B35% = 0.35*MMB.
         * 
         * @param R - longterm (average) recruitment (male-only)
         * @param cout - output stream for debug info
         * 
         * @return Fmsy
         */
        double calcFmsy(double R, ostream& cout);
        
    public:
        /**
         * Calculate equilibrium (longterm) population abundance on July 1.
         * 
         * @param R_z - longterm recruitment at size
         * @param S1_msz - survival prior to molting/growth
         * @param Th_sz - pr(molt-to-maturity|size) for immature crab
         * @param T_szz - growth transition matrix for molting crab
         * @param S2_msz - survival following molting/growth
         * @param cout - output stream for debug info
         * 
         * @return equlibrium (longterm) population abundance on July 1 (as d3_array)
         */
        d3_array calcEqNatZ(dvector& R_z,d3_array& S1_msz, 
                            dmatrix& Th_sz, d3_array& T_szz, 
                            d3_array& S2_msz, ostream& cout);
        /**
         * Calculate equilibrium unfished abundance on July 1 when 
         * longterm (average) male recruitment is R.
         * 
         * @param R - longterm (average) male recruitment
         * @param cout - output stream for debug info
         * 
         * @return equilibrium (male) population abundance on July 1
         */
        d3_array calcEqNatZF0(double R, ostream& cout);
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
        d3_array calcEqNatZFM(double R, double dirF, ostream& cout);
        /**
         * Calculate equilibrium MMB when dirF is the directed fishery capture
         * rate and longterm (average) male recruitment is R.
         * 
         * @param R - longterm (average) male recruitment
         * @param dirF - directed fishery capture rate
         * @param cout - output stream for debug info
         * 
         * @return equilibrium MMB 
         */
        double   calcEqMMBatF(double R, double dirF, ostream& cout);
};

/**
 * Class to enable sex-specific population projections.
 */
class PopProjector{
    public:
        static int debug;
        int nMSs;//number of maturity states
        int nSCs;//number of shell conditions
        int nZBs;//number of size bins
        int nFsh;//number of fisheries
        double dtF; //time to fisheries
        double dtM; //time to molting/growth
        
    public:
        double matBio;       //mature biomass at time of mating
    
    public:
        PopDyInfo* pPI;//pointer to single sex PopDyInfo object
        CatchInfo* pCI;//pointer to single sex CatchInfo object
    
    public:
        /**
         * Class constructor.
         * 
         * @param npZBs - number of size bins
         * @param npFsh - number of fisheries
         */
        PopProjector(int npZBs, int npFsh);
        ~PopProjector(){delete pPI; delete pCI;}
        /**
         * Project sex-specific population abundance forward one year,
         * based on sex-specific population abundance on July 1.
         * Also calculates:
         *      spB - spawning (mature) biomass at mating time
         *      ct_msz - total fishing mortality (abundance)
         *      rm_fmsz - retained mortality, by fishery
         *      dm_fmsz - discard mortality, by fishery
         * 
         * @param n_msz - initial abundance
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param cout - output stream for debug info
         * 
         * @return final sex-specific abundance
         */
        d3_array project(double dirF, d3_array& n_msz, ostream& cout);
        /**
         * Calculate mature biomass based on single-sex population
         * abundance on July 1.
         * 
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return mature biomass (units??)
         */
        double projectMatureBiomass(double dirF, d3_array& n_msz, ostream& cout);
        /**
         * Project sex-specific population abundance forward by time interval "dt", 
         * applying only natural mortality.
         * 
         * @param dt - time interval
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return - final abundance (after natural mortality)
         */
        d3_array applyNM(double dt, d3_array& n_msz, ostream& cout);
        /**
         * Project sex-specific population abundance forward through 
         * pulse fisheries.
         * 
         * @param n_msz - initial abundance
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param cout - output stream for debug info
         * 
         * @return - final abundance (after fisheries)
         */
        d3_array applyFM(double dirF, d3_array& n_msz, ostream& cout);
        /**
         * Project sex-specific population abundance forward through 
         * molting/growth.
         * 
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return - final abundance (after molting/growth)
         */
        d3_array applyMG(d3_array& n_msz, ostream& cout);
};

class OFL_Calculator{
    public:
        static int debug;
        int nSXs;//number of sexes
        int nFsh;//number of fisheries

    public:
        double alpha; //HCR intercept on (B/Bmsy) axis
        double beta;  //(B/Bmsy) limit below which the directed fishery is closed

    public:
        double prjMMB;  //projected ("current") MMB when OFL is taken
        dmatrix ofl_fx; //catch by fishery/sex when OFL is taken
        
    public:
        Tier3_Calculator* pT3C;//pointer to a Tier3_Calculator
        PopProjector*     pPrjM;//pointer to PopProjector for males
        PopProjector*     pPrjF;//pointer to PopProjector for females
    public:
        /**
         * Class constructor.
         * 
         * @param npFsh - number of fisheries
         */
        OFL_Calculator(int npFsh);
        /**
         * Class destructor.
         */
        ~OFL_Calculator(){delete pT3C;delete pPrjM;delete pPrjF;}
        /**
         * Calculate the Bmsy. For crab stocks, Bmsy is the longterm average
         * MMB (mature male biomass) obtained when the capture rate for males 
         * in the directed fishery is Fmsy.
         *  
         * @param tier - stock Tier level for status determination and OFL setting
         * @param R - longterm (average) MALE-only recruitment
         * @param cout - output stream for debug info
         * 
         * @return Bmsy
         */
        double calcBmsy(int tier, double R, ostream& cout);
        /**
         * Calculate Fmsy, the directed fishery male capture rate, based on the 
         * stock Tier level and the longterm male-only recruitment level (R). 
         * 
         * @param tier - stock Tier level for status determination and OFL setting
         * @param R - longterm (average) MALE-only recruitment
         * @param cout - output stream for debug info
         * 
         * @return the F that yields MSY
         */
        double calcFmsy(int tier, double R, ostream& cout);
        /**
         * Calculate Fofl using the Harvest Control Rule (HCR), based on Bmsy,
         * Fmsy, and the initial (July 1) male abundance n_msz. This is an 
         * iterative process because Fofl is given by the HCR based on the
         * projected MMB, which itself depends on what the target fishery capture
         * rate (F) was.
         * 
         * @param Bmsy - equilibrium MMB when population fished at Fmsy
         * @param Fmsy - directed fishery capture rate yielding MSY
         * @param n_msz - initial (July 1) male abundance
         * @param cout - output stream for debug info
         * 
         * @return Fofl
         */
        double calcFofl(double Bmsy, double Fmsy, d3_array& n_msz, ostream& cout);
        /**
         * Calculate the total OFL (retained+discard mortality) taken
         * when fishing on population starting with initial abundance (July 1)
         * n_xmsz at directed fishery capture rate of Fofl.
         * 
         * Note: catch (biomass) by fishery and sex is available in ofl_fx:
         *      ofl_fx(0,MALE) - retained catch (biomass) in directed fishery
         *      ofl_fx(f,x) - discard mortality (biomass) for fishery f, sex x
         * 
         * @param Fofl - directed fishery capture rate
         * @param n_xmsz - initial population abundance
         * @param cout - output stream for debug info
         * 
         * @return total OFL (biomass)
         */
        double calcOFL(double Fofl, d4_array& n_xmsz, ostream& cout);
        /**
         * Calculate the (projected) MMB when fished at Fofl.
         * 
         * @param Fofl - target fishery capture rate
         * @param n_msz - initial (July 1) population
         * @param cout - output stream for debug info
         * 
         * @return - the (projected) MMB
         */
        double calcPrjMMB(double Fofl, d3_array& n_msz, ostream& cout);
        /**
         * Calculate Fofl using Harvest Control Rule (HCR).
         * 
         * @param currMMB - "current" MMB used to determine B/Bmsy
         * @param Bmsy - Bmsy
         * @param Fmsy - Fmsy
         * @param cout - output stream for debug info
         * 
         * @return Fofl derived from HCR
         */
        double calcHCR(double prjMMB, double Bmsy, double Fmsy, ostream& cout);
};

/**
 * Convenience class encapsulating OFL results for one model (MCMC) instance
 */
class OFLResults {
    public:
        double OFL;     //total OFL (1000's t)
        //dmatrix OFL_fx; //fishery mortality components to OFL (f=0 is retained catch, f>0 is total catch mortality)
        double Fofl;    //F on directed fishery for males resulting in the OFL
        double prjB;    //projected MMB for coming year when current population fished at Fofl.
        double Fmsy;    //equilibrium F on directed fishery for males resulting in MSY
        double Bmsy;    //equilibrium MMB when fished at Fmsy 
        double B100;    //equilibrium MMB for unfished population
        
    public:
        OFLResults(){}
        ~OFLResults(){}
        
    public:
        /**
         * Write csv header for OFL results to output stream
         *  
         * @param os - output stream to write to
         */
        void writeHeaderToCSV(ostream& os);
        /**
         * Write values in csv format to output stream
         * 
         * @param os - output stream to write to
         */
        void writeToCSV(ostream& os);
        /**
         * Write values as R list to output stream
         * 
         * @param os - output stream to write to
         * @param name - name for R list
         * @param debug - flag to print debugging info
         */
        void writeToR(ostream& os, adstring name, int debug);
};

#endif	/* OFLCALCS_HPP */

