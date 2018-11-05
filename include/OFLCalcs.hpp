/* 
 * File:   OFLCalcs.hpp
 * Author: WilliamStockhausen
 *
 * Created on March 4, 2016, 10:16 AM
 */

#ifndef OFLCALCS_HPP
#define	OFLCALCS_HPP

/**
 * Class to calculate single-sex equilibrium conditions.
 */
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
         * @return equilibrium (longterm) population single-sex abundance on July 1 (as dvar3_array)
         */
        dvar3_array calcEqNatZ(dvar_vector& R_z,dvar3_array& S1_msz, 
                            dvar_matrix& Th_sz, dvar3_array& T_szz, 
                            dvar3_array& S2_msz, ostream& cout);
        /**
         * Calculate equilibrium unfished single-sex abundance on July 1 when 
         * longterm (average) single-sex recruitment is R.
         * 
         * @param R - longterm (average) male recruitment
         * @param cout - output stream for debug info
         * 
         * @return equilibrium population single-sex abundance on July 1
         */
        dvar3_array calcEqNatZF0(dvariable R, ostream& cout);
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
        dvar3_array calcEqNatZFM(dvariable R, dvariable dirF, ostream& cout);
        /**
         * Calculate single-sex equilibrium mature biomass-at-mating
         * for an unfished stock.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return mature biomass-at-mating for unfished population (1000's t)
         */
        dvariable calcEqMatureBiomassAtMatingF0(dvariable R, ostream& cout);
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
        dvariable calcEqMatureBiomassAtMatingFM(dvariable R, dvariable dirF, ostream& cout);
};

/**
 * Base class for Tier calculations to determine B0, Bmsy, and Fmsy prior to
 * OFL calculations.
 */
class Tier_Calculator {
    public:
        Equilibrium_Calculator* pEC;//pointer to Equilibrium_Calculator object for males
        
    public:
        dvariable B0;//unfished MMB
        dvariable Fmsy;//=FXX; directed fishery capture F resulting in MSY
        dvariable Bmsy;//=XX*B100; longterm B resulting from directed fishing at F
        
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
        virtual dvariable calcB0(dvariable R, ostream& cout)=0;
        /**
         * Calculate equilibrium spawning biomass (as MMB) when directed fishery
         * is fished at longterm capture rate Fmsy.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return Bmsy for population fished at Fmsy (1000's t)
         */
        virtual dvariable calcBmsy(dvariable R, ostream& cout)=0;
        /**
         * Calculate directed fishery capture rate, Fmsy, that results in 
         * equilibrium (longterm) MMB = Bmsy.
         * 
         * @param R - longterm (average) recruitment (male-only)
         * @param cout - output stream for debug info
         * 
         * @return Fmsy
         */
        virtual dvariable calcFmsy(dvariable R, ostream& cout)=0;
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
        dvariable B100;//unfished MMB
        
    public:
        /**
         * Class constructor.
         * 
         * @params XX - target SPR rate for MSY calculations
         * @param pEC - pointer to Equilibrium_Calculator object
         */
        Tier3_Calculator(double XX, Equilibrium_Calculator* pECp);
        /**
         * Class destructor (calls superclass destructor).
         */
        ~Tier3_Calculator(){}
        /**
         * Calculate B0, which is equilibrium mature male biomass (MMB)
         * for an unfished stock.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return MMB for unfished population (1000's t)
         */
        dvariable calcB0(dvariable R, ostream& cout){return calcB100(R,cout);}
        /**
         * Calculate B100, which is equilibrium mature male biomass (MMB)
         * for an unfished stock.
         * 
         * @param R - longterm (average) recruitment (male-only, millions)
         * @param cout - output stream for debug info
         * 
         * @return MMB for unfished population (1000's t)
         */
        dvariable calcB100(dvariable R, ostream& cout);
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
        dvariable calcBmsy(dvariable R, ostream& cout);
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
        dvariable calcFmsy(dvariable R, ostream& cout);
};


/**
 * Convenience class encapsulating OFL results for one model (MCMC) instance
 */
class OFLResults {
    public:
        static int debug;//flag to print debugging info
    public:
        /** average recruitment, by sex */
        dvar_vector avgRec_x;
        /** equilibrium MMB for unfished population */
        dvariable B0;
        /** equilibrium F on directed fishery for males resulting in MSY */
        dvariable Fmsy;
        /** equilibrium MMB when fished at Fmsy */
        dvariable Bmsy;
        /** equilibrium yield (1000's t) when fished at Fmsy */
        dvariable MSY;
        /** F on directed fishery for males resulting in the OFL */
        dvariable Fofl;
        /** total OFL (1000's t) */
        dvariable OFL;
        /** fishery/sex-specific mortality components to OFL (f=0 is retained catch, f>0 is total catch mortality) */
        dvar_matrix ofl_fx;
        /** projected MMB for projection year when current population is fished at Fofl. */
        dvariable prjB;
        /** "current" MMB at beginning of projection year */
        dvariable curB;
        /** "final" size distribution from assessment model */
        dvar4_array finlNatZ_xmsz; 
        /** unfished equilibrium size distribution */
        dvar4_array eqNatZF0_xmsz;
        /** fished equilibrium size distribution */
        dvar4_array eqNatZFM_xmsz;
        
        /** pointer to PopDyInfo for males */
        PopDyInfo* pPDIM;
        /** pointer to PopDyInfo for females */
        PopDyInfo* pPDIF;
        /** pointer to CatchInfo for males */
        CatchInfo* pCIM;
        /** pointer to CatchInfo for females */
        CatchInfo* pCIF;
        
        
    public:
        /** constructor */
        OFLResults();
        /** destructor */
        ~OFLResults();
        
    public:
        /**
         * Assignment operator for OFLResults class.
         * 
         * @param o - OFLResults object to copy.
         * @return - reference to the copied object
         */
        OFLResults& operator=(const OFLResults& o);
        
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
         * @param ptrMC - pointer to ModelConfiguration object
         * @param name - name for R list
         * @param debug - flag to print debugging info
         */
        void writeToR(ostream& os, ModelConfiguration* ptrMC, adstring name, int debug);
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
        dvariable prjMMB;  //projected ("current") MMB when OFL is taken
        dvar_matrix ofl_fx; //catch by sex/fishery when OFL is taken
        dvar_matrix msy_fx; //MSY catch by sex/fishery
        
    public:
        Tier_Calculator*  pTCM;  //pointer to a Tier_Calculator object for males
        Tier_Calculator*  pTCF;  //pointer to a Tier_Calculator object for females
    public:
        /**
         * Constructor.
         * 
         * @param pTCM - pointer to a Tier_Calculator object for males
         * @param pTCF - pointer to a Tier_Calculator object for females
         */
        OFL_Calculator(Tier_Calculator* pTCMp, Tier_Calculator* pTCFp);
        /**
         * Class destructor.
         */
        ~OFL_Calculator(){delete pTCM; delete pTCF;}
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
        dvariable calcHCR(dvariable currMMB, dvariable Bmsy, dvariable Fmsy, ostream& cout);
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
        dvariable calcBmsy(dvariable R, ostream& cout){return pTCM->calcBmsy(R,cout);}
        /**
         * Calculate Fmsy, the directed fishery male capture rate, based on the 
         * stock Tier level and the longterm male-only recruitment level (R). 
         * 
         * @param R - longterm (average) MALE-only recruitment
         * @param cout - output stream for debug info
         * 
         * @return the F that yields MSY
         */
        dvariable calcFmsy(dvariable R, ostream& cout){return pTCM->calcFmsy(R,cout);}
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
        dvariable calcMSY(dvar_vector R_x, dvariable Fmsy, ostream& cout);
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
        dvariable calcFofl(dvariable Bmsy, dvariable Fmsy, dvar3_array& n_msz, ostream& cout);
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
        dvariable calcOFL(dvariable Fofl, dvar4_array& n_xmsz, ostream& cout);
        /**
         * Calculate the (projected) MMB when fished at Fofl.
         * 
         * @param Fofl - target fishery capture rate
         * @param n_msz - initial (July 1) population
         * @param cout - output stream for debug info
         * 
         * @return - the (projected) MMB
         */
        dvariable calcPrjMMB(dvariable Fofl, dvar3_array& n_msz, ostream& cout);
        /**
         * Calculate all results associated with the OFL.
         * 
         * @param R - assumed equilibrium recruitment, by sex
         * @param n_xmsz - dvar4_array with initial population abundance
         * @param cout - output stream for debug info
         * 
         * @return pointer to OFLResults object.
         */
        OFLResults* calcOFLResults(dvar_vector R, dvar4_array& n_xmsz, ostream& cout);
};//OFL_Calculator

#endif	/* OFLCALCS_HPP */

