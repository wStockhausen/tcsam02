/* 
 * File:   ModelPopDyClasses.hpp
 * Author: WilliamStockhausen
 *
 * Created on November 3, 2018, 2:45 PM
 */

#ifndef MODELPOPDYCLASSES_HPP
#define MODELPOPDYCLASSES_HPP

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
        dmatrix  w_mz;     //weight-at-size
        dvar_vector  R_z;      //recruitment size distribution
        dvar3_array  M_msz;    //natural mortality
        dvar_matrix  Th_sz;    //pr(molt to maturity|pre-molt size, molt)
        dvar3_array  T_szz;    //growth matrices (indep. of molt to maturity)
        
    private:
        dvar3_array np_msz;
        dvar3_array S_msz;
    
        
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
         * Assignment operator for PopDyInfo class.
         * 
         * @param o - object to copy.
         * @return - reference to the copied object
         */
        PopDyInfo& operator=(const PopDyInfo& o);
        /**
         * Calculate single-sex mature biomass, given population abundance.
         * 
         * @param n_msz - single-sex population abundance
         * @param cout - stream for writing debug info
         * 
         * @return mature single-sex biomass (units??)
         */
        dvariable calcMatureBiomass(dvar3_array& n_msz, ostream& cout);
        /**
         * Calculate total single-sex biomass, given population abundance.
         * 
         * @param n_msz - single-sex population abundance
         * @param cout - stream for writing debug info
         * 
         * @return total single-sex biomass (units??)
         */
        dvariable calcTotalBiomass(dvar3_array& n_msz, ostream& cout);
        /**
         * Calculate single-sex survival probabilities, given natural mortality rates,
         * over a given period of time (dt)
         * 
         * @param dt - period of time (in years)
         * @param cout - stream for writing debug info
         * 
         * @return dvar3_array S_msz
         */
        dvar3_array calcSurvival(double dt, ostream& cout);
        /**
         * Calculate effects of natural mortality on population abundance for single sex.
         * 
         * @param dt - time interval (fraction of year)
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final single-sex population abundance (after applying natural mortality)
         */
        dvar3_array applyNM(double dt, dvar3_array& n_msz, ostream& cout);
        /**
         * Calculate effects of molting/growth on population abundance for single sex.
         * 
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final single-sex population abundance (after applying molting/growth)
         */
        dvar3_array applyMG(dvar3_array& n_msz, ostream& cout);
        
        /**
         * Calculate population abundance after recruitment for single sex.
         * 
         * @param R - single-sex recruitment
         * @param n_msz - initial population abundance
         * @param cout - stream for writing debug info
         * 
         * @return final single-sex population abundance after applying recruitment
         */
        dvar3_array addRecruitment(dvariable R, dvar3_array& n_msz, ostream& cout);
        
        /**
         * Write important quantities to output stream as an R list.
         * 
         * @param os - the output stream
         * @param ptrMC - pointer to the model configuration
         * @param name  - name for list in R
         * @param dbg - flag to print debugging info
         * 
         * @return void
         */
        void writeToR(ostream& os, ModelConfiguration* ptrMC, adstring name, int dbg);
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
        
    public:
        /** maximum capture rate in target fishery */
        dvariable   maxF;
        /** fishery capture rates */
        dvar3_array cpF_fms;
        /** size-specific fishery capture rates */
        dvar4_array cpF_fmsz;
        /** size-specific fishery retained mortality rates */
        dvar4_array rmF_fmsz;
        /** size-specific fishery discard mortality rates */
        dvar4_array dmF_fmsz;
        /** selectivity functions */
        dvar4_array selF_fmsz;
        /** retention functions */
        dvar4_array retF_fmsz;
        /** discard mortality rates */
        dvar_vector  hm_f; 
        
    public:
        /** total catch mortality (abundance) */
        dvar3_array cmN_msz; 
        /** total captured, by fishery (abundance) */
        dvar4_array cpN_fmsz;//
        /** retained mortality, by fishery (abundance) */
        dvar4_array rmN_fmsz;//
        /** discard mortality, by fishery (abundance) */
        dvar4_array dmN_fmsz;//
        
    private:
        dvar_vector totFM;//total fishing mortality
        dvar3_array S_msz;//survival following fisheries
          
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
         * Assignment operator for CatchInfo class.
         * 
         * @param o - object to copy.
         * @return - reference to the copied object
         */
        CatchInfo& operator=(const CatchInfo& o);
        /**
         * Calculate maximum fishing capture rate for target fishery (f=1)
         * 
         * Modifies: maxF
         * 
         * @param cout - stream for writing debug info
         * 
         * @return maxF
         */
        dvariable findMaxTargetCaptureRate(ostream& cout);
        /**
         * Calculate single-sex catch abundance (cm_msz, cp_fmsz, rm_fmsz, dm_fmsz) and 
         * post-fisheries abundance based on target fishery capture rate 'dirF' 
         * and pre-fisheries population abundance n_msz.
         * 
         * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
         * 
         * Modifies: 
         *      cm_msz  - total fishing mortality (abundance)
         *      cp_fmsz - fishery captures, by fishery (abundance)
         *      rm_fmsz - retained mortality, by fishery (abundance)
         *      dm_fmsz - discard mortality, by fishery
         * 
         * @param dirF - fully-selected directed fishery capture rate
         * @param n_msz - pre-fisheries population size
         * @param cout - stream for writing debug info
         * 
         * @return np_msz - post-fisheries population size (dvar3_array)
         * 
         */
        dvar3_array applyFM(dvariable dirF, dvar3_array& n_msz, ostream& cout);
        /**
         * Calculates probabilities of surviving fisheries, given directed 
         * fishing capture rate 'dirF'.
         * 
         * @param dirF - fully-selected capture rate for target fishery (f=1)
         * @param cout - stream for writing debug info
         * 
         * @return survival probabilities S_msz as dvar3_array.
         */
        dvar3_array calcSurvival(dvariable dirF, ostream& cout);
        /**
         * Set single-sex size-specific fishery capture rates.
         * 
         * @param capFp_fmsz - input capture rates
         */
        void setCaptureRates(dvar4_array& capFp_fmsz);
        /**
         * Set single-sex fully-selected fishery capture rates.
         * 
         * @param capFp_fms - input capture rates
         */
        void setCaptureRates(dvar3_array& capFp_fms);
        /**
         * Set single-sex selectivity functions.
         * 
         * @param selFp_fmsz - input selectivity functions
         */
        void setSelectivityFcns(dvar4_array& selFp_fmsz);
        /**
         * Set single-sex retention functions.
         * 
         * @param retFp_fmsz - input retention functions
         */
        void setRetentionFcns(dvar4_array& retFp_fmsz);
        /**
         * Set handling mortality rates, by fishery.
         * 
         * @param pHM_f - dvar_vector of handling mortality rates
         */
        void setHandlingMortality(dvar_vector& pHM_f);
        /**
         * Write important quantities to output stream as an R list.
         * 
         * @param os - the output stream
         * @param ptrMC - pointer to the model configuration
         * @param name  - name for list in R
         * @param dbg - flag to print debugging info
         * 
         * @return void
         */
        void writeToR(ostream& os, ModelConfiguration* ptrMC, adstring name, int dbg);
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
        dvariable matBio;       //mature biomass at time of mating
    
    public:
        PopDyInfo* pPI;//pointer to single sex PopDyInfo object
        CatchInfo* pCI;//pointer to single sex CatchInfo object
        
    private:
        dvar3_array n1_msz;
        dvar3_array n2_msz;
        dvar3_array n3_msz;
        dvar3_array n4_msz;
        dvar3_array n5_msz;
    
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
         * based on single-sex population abundance on July 1. 
         * 
         * NOTE: Must set dtF and dtM prior to calling this method.
         * 
         * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
         * NOTE: Recruitment is NOT added in.
         * 
         * Also calculates:
         *      matBio - mature biomass at mating time
         *      pCI elements:
         *          cm_msz  - total fishing mortality        (abundance)
         *          cp_fmsz - fishery captures, by fishery   (abundance)
         *          rm_fmsz - retained mortality, by fishery (abundance)
         *          dm_fmsz - discard mortality, by fishery  (abundance)
         * 
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return final sex-specific abundance, WITHOUT RECRUITMENT
         */
        dvar3_array project(dvariable dirF, dvar3_array& n_msz, ostream& cout);
        
        /**
         * Project unfished single-sex population abundance forward one year,
         * based on single-sex population abundance on July 1. Recruitment
         * is NOT added in.
         * 
         * Also calculates:
         *      matBio - mature biomass at mating time
         *      pCI elements:
         *          cm_msz  - total fishing mortality (abundance) [=0]
         *          cp_fmsz - fishery captures (abundance)        [=0]
         *          rm_fmsz - retained mortality, by fishery      [=0]
         *          dm_fmsz - discard mortality, by fishery       [=0]
         * 
         * 
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return final sex-specific abundance, without recruitment
         */
        dvar3_array projectUnFished(dvar3_array& n_msz, ostream& cout);
        
        /**
         * Add single-sex recruitment to population abundance for one sex.
         * 
         * @param R - total single-sex recruitment
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return final single-sex abundance with recruitment
         */
        dvar3_array addRecruitment(dvariable R, dvar3_array& n_msz, ostream& cout);
        
        /**
         * Calculate mature biomass-at-mating based on single-sex population
         * abundance on July 1.
         * 
         * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
         * 
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return mature biomass-at-mating (1000's t)
         */
        dvariable projectMatureBiomassAtMating(dvariable dirF, dvar3_array& n_msz, ostream& cout);
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
        dvar3_array applyNM(dvariable dt, dvar3_array& n_msz, ostream& cout);
        /**
         * Project sex-specific population abundance forward through 
         * pulse fisheries.
         * 
         * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
         * 
         * @param n_msz - initial abundance
         * @param dirF - multiplier on fishing mortality rate in directed fishery
         * @param cout - output stream for debug info
         * 
         * @return - final abundance (after fisheries)
         */
        dvar3_array applyFM(dvariable dirF, dvar3_array& n_msz, ostream& cout);
        /**
         * Project sex-specific population abundance forward through 
         * molting/growth.
         * 
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         * @return - final abundance (after molting/growth)
         */
        dvar3_array applyMG(dvar3_array& n_msz, ostream& cout);
        /**
         * Get the mature biomass-at-mating. Must have called \code{project}
         * first.
         * 
         * @return - (dvariable) mature biomass-at-mating for the projected year
         */
        dvariable getMatureBiomassAtMating(){return matBio;}
        /**
         * Get the single-sex fishery capture abundance, by fishery,
         * maturity state, shell condition and size.
         *  
         * NOTE: Must have called \code{project} first.
         * 
         * @return - (dvar4_array) capture abundance by fishery (fmsz)
         */
        dvar4_array getFisheriesCaptureAbundance(){return pCI->cpN_fmsz;}
        /**
         * Get the single-sex fishery retained catch mortality, by fishery,
         * maturity state, shell condition and size.
         *  
         * NOTE: Must have called \code{project} first.
         * 
         * @return - (dvar4_array) retained catch mortality by fishery (fmsz)
         */
        dvar4_array getRetainedCatchMortality(){return pCI->rmN_fmsz;}
        /**
         * Get the single-sex fishery discard catch mortality, by fishery,
         * maturity state, shell condition and size.
         *  
         * NOTE: Must have called \code{project} first.
         * 
         * @return - (dvar4_array) discard catch mortality by fishery (fmsz)
         */
        dvar4_array getDiscardCatchMortality(){return pCI->dmN_fmsz;}
        /**
         * Get the single-sex total catch mortality, by 
         * maturity state, shell condition and size.
         *  
         * NOTE: Must have called \code{project} first.
         * 
         * @return - (dvar3_array) total catch mortality (msz)
         */
        dvar3_array getTotalCatchMortality(){return pCI->cmN_msz;}
};

/**
 * Class to make multi-year single-sex population projections.
 */
class MultiYearPopProjector {
    public:
        static int debug;
        
    public:
        PopProjector* pPP;//pointer to single-sex population projector
        
    public:
        dvar4_array n_ymsz;//numbers-at-size
        dvar_vector matBio_y;//mature biomass
        dvar_matrix cp_yf;   //captured biomass, by fishery
        dvar_matrix rm_yf;   //retained catch biomass, by fishery
        dvar_matrix dm_yf;   //discard catch mortality, by fishery
        dvar_vector totCM_y; //total catch biomass
        
    public:
        /**
         * Class constructor.
         * 
         * @param pPPp - pointer to a PopProjector object to base equilibrium calculations on
         */
        MultiYearPopProjector(PopProjector* pPPp){pPP = pPPp;}
        /**
         * Class destructor.
         */
        ~MultiYearPopProjector(){delete pPP;}
        /**
         * Project multiple years at constant recruitment and directed F.
         * 
         * NOTE: If dirF &lt 0, then directed fishing mortality is not rescaled.
         * 
         * Calculates:
         *      n_ymsz   - numbers-at-size
         *      matBio_y - mature biomass
         *      cp_yf    - captured biomass, by fishery
         *      rm_yf    - retained catch biomass, by fishery
         *      dm_yf    - discard catch mortality, by fishery
         *      totCM_y  - total catch biomass
         * 
         * @param n - number of years to project
         * @param R - (constant) single-sex recruimtent
         * @param dirF - (constant) directed F
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         */
        void project(int n, dvariable R, dvariable dirF, dvar3_array& n_msz, ostream& cout);
        /**
         * Project multiple years at constant recruitment, with no fishing.
         * 
         * Calculates:
         *      n_ymsz   - numbers-at-size
         *      matBio_y - mature biomass
         *      cp_yf    - captured biomass, by fishery        [=0]
         *      rm_yf    - retained catch biomass, by fishery  [=0]
         *      dm_yf    - discard catch mortality, by fishery [=0]
         *      totCM_y  - total catch biomass                 [=0]
         * 
         * @param n - number of years to project
         * @param R - (constant) single-sex recruitment
         * @param n_msz - initial abundance
         * @param cout - output stream for debug info
         * 
         */
        void projectUnFished(int n, dvariable R, dvar3_array& n_msz, ostream& cout);
};



#endif /* MODELPOPDYCLASSES_HPP */

