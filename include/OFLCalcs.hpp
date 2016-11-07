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
        double dtM; //time at which mating occurs (as fraction of year)
        
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
         * @params dtMp - time at which mating occurs (as fraction of year)
         */
        PopDyInfo(int npZBs, double dtMp);
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
         * @return mature single-sex biomass (units??)
         */
        double calcMatureBiomass(d3_array& n_msz, ostream& cout);
        /**
         * Calculate total single-sex biomass, given population abundance.
         * 
         * @param n_msz - single-sex population abundance
         * @param cout - stream for writing debug info
         * 
         * @return total single-sex biomass (units??)
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
         * Calculate effects of natural mortality on population abundance for single sex.
         * 
         * @param dt - time interval (fraction of year)
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final single-sex population abundance (after applying natural mortality)
         */
        d3_array applyNM(double dt, d3_array& n_msz, ostream& cout);
        /**
         * Calculate effects of molting/growth on population abundance for single sex.
         * 
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final single-sex population abundance (after applying molting/growth)
         */
        d3_array applyMG(d3_array& n_msz, ostream& cout);
        
        /**
         * Calculate population abundance after recruitment for single sex.
         * 
         * @param R - single-sex recruitment
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final single-sex population abundance after applying recruitment
         */
        d3_array addRecruitment(double R, d3_array& n_msz, ostream& cout);
};//PopDyInfo

/**
 * Class to facilitate fishery single-sex catch calculations.
 */
class CatchInfo {
    public:
        static int debug;
        int nMSs;//number of maturity states
        int nSCs;//number of shell conditions
        int nZBs;//number of size bins
        int nFsh;//number of fisheries
        double dtF; //time at which fisheries occur (as fraction of year)
        
    public:
        double   maxF;     //maximum capture rate in target fishery
        d4_array capF_fmsz;//fishery capture rates
        d4_array retF_fmsz;//retention function
        dvector  hm_f;     //discard mortality rates
        
    public:
        d3_array cp_msz; //total captured (abundance)
        d3_array cm_msz; //catch mortality (abundance)
        d4_array rm_fmsz;//retained mortality (abundance)
        d4_array dm_fmsz;//discard mortality (abundance)
        
    public:
        /**
         * Class constructor.
         * 
         * @param npZBs - number of size bins
         * @param npFsh - number of fisheries
         * @param dtFp - time at which fisheries occur (as fraction of year)
         */
        CatchInfo(int npZBs, int npFsh, double dtFp);
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
         * Set single-sex fishery capture rates.
         * 
         * @param capFp_fmsz - input capture rates
         */
        void setCaptureRates(d4_array& capFp_fmsz);
        /**
         * Set single-sex retention functions.
         * 
         * @param retFp_fxmsz - input retention functions
         */
        void setRetentionFcns(d4_array& retFp_fmsz);
        /**
         * Set handling mortality rates, by fishery.
         * 
         * @param pHM_f - dvector of handling mortality rates
         */
        void setHandlingMortality(dvector& pHM_f);
};//CatchInfo

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
         * @param pPIp - pointer to single sex PopDyInfo object
         * @param pCIp - pointer to single sex CatchInfo object
         */
        PopProjector(PopDyInfo* pPIp, CatchInfo* pCIp);
        /**
         * Class destructor
         */
        ~PopProjector(){delete pPI; delete pCI;}
        /**
         * Project single-sex population abundance forward one year,
         * based on single-sex population abundance on July 1. Recruitment
         * is NOT added in.
         * 
         * Also calculates:
         *      spB - spawning (mature) biomass at mating time
         *      cp_msz - total fishery captures (abundance)
         *      cm_msz - total fishing mortality (abundance)
         *      rm_fmsz - retained mortality, by fishery
         *      dm_fmsz - discard mortality, by fishery
         * 
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return final sex-specific abundance, without recruitment
         */
        d3_array project(double dirF, d3_array& n_msz, ostream& cout);
        /**
         * Add single-sex recruitment to population abundance for one sex.
         * 
         * @param R - total single-sex recruitment
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return final single-sex abundance with recruitment
         */
        d3_array addRecruitment(double R, d3_array& n_msz, ostream& cout);
        /**
         * Calculate mature biomass-at-mating based on single-sex population
         * abundance on July 1.
         * 
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return mature biomass-at-mating (1000's t)
         */
        double projectMatureBiomassAtMating(double dirF, d3_array& n_msz, ostream& cout);
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

class Equilibrium_Calculator {
    public:
        static int debug;
        
    public:
        PopProjector* pPP;//pointer to single-sex population projector
        
    public:
        /**
         * Class constructor.
         * 
         * @param pPPp - pointer to a PopProjector object to base equilibrium calculations on
         */
        Equilibrium_Calculator(PopProjector* pPPp){pPP = pPPp;}
        /**
         * Class destructor.
         */
        ~Equilibrium_Calculator(){delete pPP;}
        /**
         * Calculate equilibrium (longterm) population single-sex abundance on July 1.
         * 
         * @param R_z - longterm single-sex recruitment at size
         * @param S1_msz - survival prior to molting/growth
         * @param Th_sz - pr(molt-to-maturity|size) for immature crab
         * @param T_szz - growth transition matrix for molting crab
         * @param S2_msz - survival following molting/growth
         * @param cout - output stream for debug info
         * 
         * @return equilibrium (longterm) population single-sex abundance on July 1 (as d3_array)
         */
        d3_array calcEqNatZ(dvector& R_z,d3_array& S1_msz, 
                            dmatrix& Th_sz, d3_array& T_szz, 
                            d3_array& S2_msz, ostream& cout);
        /**
         * Calculate equilibrium unfished single-sex abundance on July 1 when 
         * longterm (average) single-sex recruitment is R.
         * 
         * @param R - longterm (average) male recruitment
         * @param cout - output stream for debug info
         * 
         * @return equilibrium population single-sex abundance on July 1
         */
        d3_array calcEqNatZF0(double R, ostream& cout);
        /**
         * Calculate equilibrium single-sex abundance on July 1 when dirF is the 
         * directed fishery capture rate and the longterm (average) 
         * single-sex recruitment is R.
         * 
         * @param R - longterm (average) single-sex recruitment
         * @param dirF - directed fishery capture rate
         * @param cout - output stream for debug info
         * 
         * @return equilibrium population single-sex abundance on July 1
         */
        d3_array calcEqNatZFM(double R, double dirF, ostream& cout);
        /**
         * Calculate single-sex equilibrium mature biomass-at-mating
         * for an unfished stock.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return mature biomass-at-mating for unfished population (1000's t)
         */
        double calcEqMatureBiomassAtMatingF0(double R, ostream& cout);
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
        double calcEqMatureBiomassAtMatingFM(double R, double dirF, ostream& cout);
};

/**
 * Base class for Tier calculations to determine B0, Bmsy, and Fmsy prior to
 * OFL calculations.
 */
class Tier_Calculator {
    public:
        Equilibrium_Calculator* pEC;//pointer to Equilibrium_Calculator object for males
        
    public:
        double B0;//unfished MMB
        double Fmsy;//=FXX; directed fishery capture F resulting in MSY
        double Bmsy;//=XX*B100; longterm B resulting from directed fishing at F
        
    public:
        /**
         * Class constructor.
         * 
         * @param pEC - pointer to Equilibrium_Calculator object
         */
        Tier_Calculator(Equilibrium_Calculator* pECp){pEC = pECp;}
        /**
         * Class destructor.
         */
        virtual ~Tier_Calculator(){delete pEC;}
        /**
         * Calculate B0, which is equilibrium mature male biomass (MMB)
         * for an unfished stock.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return MMB for unfished population (1000's t)
         */
        virtual double calcB0(double R, ostream& cout)=0;
        /**
         * Calculate equilibrium spawning biomass (as MMB) when directed fishery
         * is fished at longterm capture rate Fmsy.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return Bmsy for population fished at Fmsy (1000's t)
         */
        virtual double calcBmsy(double R, ostream& cout)=0;
        /**
         * Calculate directed fishery capture rate, Fmsy, that results in 
         * equilibrium (longterm) MMB = Bmsy.
         * 
         * @param R - longterm (average) recruitment (male-only)
         * @param cout - output stream for debug info
         * 
         * @return Fmsy
         */
        virtual double calcFmsy(double R, ostream& cout)=0;
};

/**
 * Tier_Calculator class to do Tier 3 calculations for B0, Bmsy, and Fmsy
 * based on SPR considerations (B0=B100, Bmsy=B35, Fmsy=F35).
 */
class Tier3_Calculator : public Tier_Calculator {
    public:
        static int debug;
    public:
        double XX;  //SPR level for Fmsy, Bmsy
        double B100;//unfished MMB
        
    public:
        /**
         * Class constructor.
         * 
         * @params XX - target SPR rate for MSY calculations
         * @param pEC - pointer to Equilibrium_Calculator object
         */
        Tier3_Calculator(double XX, Equilibrium_Calculator* pECp);
        /**
         * Class destructor
         */
        ~Tier3_Calculator(){Tier_Calculator::~Tier_Calculator();}
        /**
         * Calculate B0, which is equilibrium mature male biomass (MMB)
         * for an unfished stock.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return MMB for unfished population (1000's t)
         */
        double calcB0(double R, ostream& cout){return calcB100(R,cout);}
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
};


/**
 * Convenience class encapsulating OFL results for one model (MCMC) instance
 */
class OFLResults {
    public:
        dvector avgRec_x; //average recruitment, by sex
        double B0;       //equilibrium MMB for unfished population
        double Fmsy;     //equilibrium F on directed fishery for males resulting in MSY
        double Bmsy;     //equilibrium MMB when fished at Fmsy 
        double Fofl;     //F on directed fishery for males resulting in the OFL
        double OFL;      //total OFL (1000's t)
        //dmatrix OFL_fx; //fishery mortality components to OFL (f=0 is retained catch, f>0 is total catch mortality)
        double prjB;     //projected MMB for coming year when current population fished at Fofl.
        
    public:
        OFLResults(){}
        ~OFLResults(){}
        
    public:
        /**
         * Write csv header for OFL results to output stream
         *  
         * @param os - output stream to write to
         */
        void writeCSVHeader(ostream& os);
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
};//OFLResults
/**
 * Class to calculate OFL results.
 */
class OFL_Calculator{
    public:
        static int debug;

    public:
        int tier;     //Tier for status determination and OFL calculation
        double alpha; //HCR intercept on (B/Bmsy) axis
        double beta;  //(B/Bmsy) limit below which the directed fishery is closed

    public:
        double prjMMB;  //projected ("current") MMB when OFL is taken
        dmatrix ofl_xf; //catch by sex/fishery when OFL is taken
        
    public:
        Tier_Calculator*  pTC;  //pointer to a Tier_Calculator object
        PopProjector*     pPrjF;//pointer to PopProjector for females
    public:
        /**
         * Constructor.
         * 
         * @param pTC   - pointer to a Tier_Calculator object
         * @param pPrjF - pointer to PopProjector for females
         */
        OFL_Calculator(Tier_Calculator* pTC, PopProjector* pPrjF);
        /**
         * Class destructor.
         */
        ~OFL_Calculator(){delete pTC; delete pPrjF;}
        /**
         * Calculate Fofl using the Harvest Control Rule (HCR).
         * 
         * @param currMMB - "current" MMB used to determine B/Bmsy
         * @param Bmsy - Bmsy
         * @param Fmsy - Fmsy
         * @param cout - output stream for debug info
         * 
         * @return Fofl derived from HCR
         */
        double calcHCR(double currMMB, double Bmsy, double Fmsy, ostream& cout);
        /**
         * Calculate the Bmsy. For crab stocks, Bmsy is the longterm average
         * MMB (mature male biomass) obtained when the capture rate for males 
         * in the directed fishery is Fmsy.
         *  
         * @param R - longterm (average) MALE-only recruitment
         * @param cout - output stream for debug info
         * 
         * @return Bmsy
         */
        double calcBmsy(double R, ostream& cout){return pTC->calcBmsy(R,cout);}
        /**
         * Calculate Fmsy, the directed fishery male capture rate, based on the 
         * stock Tier level and the longterm male-only recruitment level (R). 
         * 
         * @param R - longterm (average) MALE-only recruitment
         * @param cout - output stream for debug info
         * 
         * @return the F that yields MSY
         */
        double calcFmsy(double R, ostream& cout){return pTC->calcFmsy(R,cout);}
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
         * Calculate all results associated with the OFL.
         * 
         * @param R - assumed equilibrium recruitment, by sex
         * @param n_xmsz - d4_array with initial population abundance
         * @param cout - output stream for debug info
         * 
         * @return OFLResults object.
         */
        OFLResults calcOFLResults(dvector R, d4_array& n_xmsz, ostream& cout);
};//OFL_Calculator

#endif	/* OFLCALCS_HPP */

