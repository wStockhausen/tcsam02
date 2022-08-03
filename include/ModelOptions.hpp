/* 
 * ModelOptions.hpp
 */

#ifndef MODELOPTIONS_HPP
#define MODELOPTIONS_HPP

#include <admodel.h>
#include "ModelConfiguration.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelParameterInfoTypes.hpp"

/**
 * Class encapsulating information on a single empirical selectivity function
 * 
 * This class encapsulates information on a single empirical selectivity function.
 */
class EmpiricalSelFcn{
    public:
        /** flag to print debugging info */
        static int debug;
    private:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
    public:
        /** id associated with the empirical selectivity function */
        int id;
        /* sizes at which the empirical sel function was evaluated */
        dvector zBs;
        /** values for empirical sel function */
        dvector esf;
    public:
        /**
         * Class instantiator
         * 
         * @param mc - pointer to the model configuration object
         */
        EmpiricalSelFcn(ModelConfiguration& mc);
        /**
         * Class destructor
         */
        ~EmpiricalSelFcn();
        
        /**
         * Read from portion of ModelOptions input stream
         * 
         * @param is - input stream
         */
        void read(cifstream & is);
        /**
         * Write to stream
         * 
         * @param os - output stream
         */
        void write(std::ostream & os);
        /**
         * Write to stream in R format
         * 
         * @param os - output stream
         */
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EmpiricalSelFcn & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EmpiricalSelFcn & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on ALL empirical selectivity functions
 * 
 * This class encapsulates information on ALL empirical selectivity functions.
 */
class EmpiricalSelFcns{
    public:
        /** flag to print debugging info */
        static int debug;
    private:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
    public:
        /** number of empirical selectivity functions */
        int nESFs; 
        /** pointer to array of pointers to individual empirical selectivity functions */
        EmpiricalSelFcn** ppESFs;
    public:
        EmpiricalSelFcns(ModelConfiguration& mc);
        ~EmpiricalSelFcns();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EmpiricalSelFcns & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EmpiricalSelFcns & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on a single empirical selectivity function prior
 * 
 * This class encapsulates information on a single empirical selectivity function prior.
 */
class EmpiricalSelFcnPrior{
    public:
        /** flag to print debugging info */
        static int debug;
    private:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
    public:
        /** id associated with the empirical selectivity function prior */
        int id;
        /** id of selectivity function to apply prior to */
        int sel_id;
        /* weight to assign to prior probability */
        double priorWgt;
        /* pdf type for prior */
        adstring priorType;
        /* vector of size bins for likelihood function */
        dvector zBs;
        /* vector of first parameter values for likelihood function */
        dvector p1;
        /* vector of second parameter values for likelihood function */
        dvector p2;
        /* */
        ModelPDFInfo* pMPI;
    public:
        /**
         * Class instantiator
         * 
         * @param mc - pointer to the model configuration object
         */
        EmpiricalSelFcnPrior(ModelConfiguration& mc);
        /**
         * Class destructor
         */
        ~EmpiricalSelFcnPrior();
        /**
         * Gets the multiplicative weight set on the prior probability in the likelihood.
         * 
         * @return - the weight
         */
        double getPriorWgt(){return priorWgt;}
        /**
         * Gets the prior type, as an adstring.
         * 
         * @return - the prior type, as an adstring object 
         */
        adstring getPriorType(){return priorType;}
        /**
         * Sets the prior type, based on an adstring value.
         * 
         * @param prior - the prior type, as an adstring
         */
        void setPriorType(adstring & prior);
        
        /**
         * Read from portion of ModelOptions input stream
         * 
         * @param is - input stream
         */
        void read(cifstream & is);
        /**
         * Write to stream
         * 
         * @param os - output stream
         */
        void write(std::ostream & os);
        /**
         * Write to stream in R format
         * 
         * @param os - output stream
         */
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EmpiricalSelFcnPrior & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EmpiricalSelFcnPrior & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on ALL empirical selectivity function priors
 * 
 * This class encapsulates information on ALL empirical selectivity function priors.
 */
class EmpiricalSelFcnPriors{
    public:
        /** flag to print debugging info */
        static int debug;
    private:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
    public:
        /** number of empirical selectivity function priors */
        int nESPs; 
        /** pointer to array of pointers to individual empirical selectivity function priors */
        EmpiricalSelFcnPrior** ppESPs;
    public:
        EmpiricalSelFcnPriors(ModelConfiguration& mc);
        ~EmpiricalSelFcnPriors();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EmpiricalSelFcnPriors & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EmpiricalSelFcnPriors & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on a single effort averaging scenario
 * 
 * This class encapsulates information on a single fishery effort averaging scenario 
 * used to extrapolate fully-selected fishing mortality rates to time periods
 * when effort data is available but catch data is not.
 */
class EffAvgScenario{
    public:
        /** flag to print debugging info */
        static int debug;
    private:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
    public:
        /** id associated with the effort averaging scenario */
        int id;
        /** index of associated fishery */
        int f;//
        /** pointer to an IndexBlock object */
        IndexBlock* ptrIB;//
    private:
        /** names for index variables */
        adstring_array strVals;//
    public:
        EffAvgScenario(ModelConfiguration& mc);
        ~EffAvgScenario();
        
        /**
         * Sets the max year for the averaging period (for retrospective analysis).
         * 
         * @param mxYr - max year for averaging 
         */
        void setMaxYearForAveraging(int mxYr);
        
        /**
         * Gets the time block for the effort averaging scenario as an R array dimension
         * 
         * @return - an adstring object
         */
        adstring getYDimsForR(){return ptrIB->getAsRDim();}
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EffAvgScenario & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EffAvgScenario & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on ALL effort averaging scenarios
 * 
 * This class encapsulates information on ALL fishery effort averaging scenarios 
 * used to extrapolate fully-selected fishing mortality rates to time periods
 * when effort data is available but catch data is not.
 */
class EffAvgScenarios{
    public:
        /** flag to print debugging info */
        static int debug;
    private:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
    public:
        /** number of effort averaging scenarios */
        int nAvgs; 
        /** pointer to array of pointers to individual effort averaging scenarios */
        EffAvgScenario** ppEASs;
    public:
        EffAvgScenarios(ModelConfiguration& mc);
        ~EffAvgScenarios();
        
        /**
         * Sets the max year for all averaging periods (for retrospective analysis).
         * 
         * @param mxYr - max year for averaging 
         */
        void setMaxYearForAveraging(int mxYr);
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EffAvgScenarios & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EffAvgScenarios & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on a single capture rate averaging scenario
 * 
 *  This class encapsulates information on a single fishery capture rate averaging scenario 
 *  used in the OFL calculations.
 */
class CapRateAvgScenario{
    public:
        /** flag to print debugging info */
        static int debug;
    public:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
        int id;
        int f;//fishery index
        int x;//sex index
        int m;//maturity state index
        int s;//shell condition index
        int idParam;//id of associated pLnEffX parameter in ModelParametersInfo file
        int idEffAvgInfo;//id of associated effort averaging info object
        int optAvg;      //averaging option
        double llWgt;    //likelihood weight
    private:
        adstring_array strVals;//names for index variables
    public:
        CapRateAvgScenario(ModelConfiguration& mc);
        ~CapRateAvgScenario();
        
//        /**
//         * Sets the max year for the averaging period (for retrospective analysis).
//         * 
//         * @param mxYr - max year for averaging 
//         */
//        void setMaxYearForAveraging(int mxYr);
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, CapRateAvgScenario & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, CapRateAvgScenario & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on all capture rate averaging scenarios for OFL calculations.
 */
class CapRateAvgScenarios{
    public:
        /** flag to print debugging info */
        static int debug;
    public:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
        int nAvgs; //number of capture rate averaging info objects
        CapRateAvgScenario** ppCRASs;//pointer to array of CapRateAvgInfo pointers 
    public:
        CapRateAvgScenarios(ModelConfiguration& mc);
        ~CapRateAvgScenarios();
        
//        /**
//         * Sets the max year for all averaging periods (for retrospective analysis).
//         * 
//         * @param mxYr - max year for averaging 
//         */
//        void setMaxYearForAveraging(int mxYr);
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, CapRateAvgScenarios & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, CapRateAvgScenarios & obj){obj.write(os); return os;}
};

/**
 * Class encapsulating information on effort extrapolation scenarios.
 */
class EffXtrapScenarios{
    public:
        /** flag to print debugging info */
        static int debug;
    private:
        /** pointer to the global ModelConfiguration object */
        ModelConfiguration* ptrMC;
    public:
        /** pointer to effort averaging scenarios object */
        EffAvgScenarios* ptrEffAvgScenarios;//pointer to effort averaging scenarios object 
        /** pointer to capture rate averaging scenarios object */
        CapRateAvgScenarios* ptrCapRateAvgScenarios;//pointer to capture rate averaging scenarios object 
    public:
        /**
         * Class constructor.
         * 
         * @param mc - reference to the ModelConfigurtion object
         */
        EffXtrapScenarios(ModelConfiguration& mc);
        /**
         * Class destructor for effort extrapolation scenarios.
         */
        ~EffXtrapScenarios();
        
        /**
         * Sets the max year for all averaging periods (for retrospective analysis).
         * 
         * @param mxYr - max year for averaging 
         */
        void setMaxYearForAveraging(int mxYr);
        
        /**
         * Gets the years defining the id'th averaging time period
         * 
         * @param id - the id of the desired time period
         * @return - an ivector with the years corresponding to the time period
         */
        ivector getTimePeriodForCapRateAveraging(int id);
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EffXtrapScenarios & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EffXtrapScenarios & obj){obj.write(os); return os;}
};

/**
 * ProjectionOptions class definition.
 */
class ProjectionOptions {
public:
    /** flag to print debugging info */
    static int debug;
    
public:
    /** number of repetitions to project*/
    int nReps;
    /** number of years to project */
    int nYrs;
    /** number of F's to project */
    int nFs;
    /** F's to project */
    dvector Fs;
    /** number of F multipliers to project */
    int nFMs;
    /** F multipliers to project */
    dvector FMs;
        /**
         * Class constructor.
         */
        ProjectionOptions();
        /**
         * Class destructor for projection scenarios.
         */
        ~ProjectionOptions();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, ProjectionOptions & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, ProjectionOptions & obj){obj.write(os); return os;}
};

/**
 * CatchDataSimOptions class
 */
class CatchDataSimOptions{
public:
    static int debug; //debugging flag
    adstring name;    //fleet name
    long rngSeed;     //seed for random number generator 
    double expFacAbd; //expansion factor for catch abundance data
    double expFacBio; //expansion factor for catch biomass data
    double expFacZCs; //expansion factor for size composition data
    
    CatchDataSimOptions();
    ~CatchDataSimOptions();
    void read(cifstream & is);
    void write(std::ostream & os);
    void writeToR(std::ostream & os);    
    
    friend cifstream& operator >>(cifstream & is, CatchDataSimOptions & obj){obj.read(is); return is;}
    friend std::ostream& operator <<(std::ostream & os, CatchDataSimOptions & obj){obj.write(os); return os;}
};

/**
 * SimOptions class definition.
 */
class SimOptions {
public:
    /** flag to print debugging info */
    static int debug;
    
public:
    /** pointer to ModelConfiguration object */
    int nRetCatch;
    int nTotCatch;
    int nDscCatch;
    int nIdxCatch;
    CatchDataSimOptions** ppRetCatch;
    CatchDataSimOptions** ppTotCatch;
    CatchDataSimOptions** ppDscCatch;
    CatchDataSimOptions** ppIdxCatch;
    
    /** random number generator flag for growth data */
    int grwRngSeed;//--if not 0, reset rng seed to this when generating random numbers
    /* CV multiplier for growth data */
    double grwMultFac;
    /** random number generator flag for maturity ogive data */
    int modRngSeed;//--if not 0, reset rng seed to this when generating random numbers
    /* sample size divisor for maturity ogive data */
    double modDivFac;
    /**
     * Class constructor.
     */
    SimOptions();
    /**
     * Class destructor for projection scenarios.
     */
    ~SimOptions();

    void read(cifstream & is);
    void write(std::ostream & os);
    void writeToR(std::ostream & os);

    friend cifstream& operator >>(cifstream & is, SimOptions & obj){obj.read(is); return is;}
    friend std::ostream& operator <<(std::ostream & os, SimOptions & obj){obj.write(os); return os;}
};

/**
 * ModelOptions class definition.
 */
    class ModelOptions {
    public:
        /** version string for class */
        const static adstring VERSION;
        /** flag to print debug info */
        static int debug;        
    public:
        /** pointer to model configuration object */
        ModelConfiguration* ptrMC; 
        
        /** labels for initial n-at-z options */
        adstring_array optsInitNatZ;
        /** selected option for initial n-at-z calculations */
        int optInitNatZ;               
        
        /** labels for options for natural mortality parameterization */
        adstring_array optsParamNM;
        /** integer indicating option for natural mortality parameterization */
        int optParamNM;
        
        /** labels for growth parameterization options */
        adstring_array optsGrowthParam;  
        /** selected option for growth parameterization */
        int optGrowthParam;                 
        
        /** labels for growth pdf options */
        adstring_array optsGrowthPDF;  
        /** selected option for growth pdf */
        int optGrowthPDF;   
        
        /** max extent of size bins for growth probabilities */
        int maxGrowthZBEx;
        
        /** min pre-molt size to apply penalty to prevent negative growth increments */
        double minGrowthCW;
        /** max pre-molt size to apply penalty to prevent negative growth increments */
        double maxGrowthCW;
        /** likelihood weight on penalties to prevent negative growth increments */
        double wgtNegGrowth;
        /** eps in penalty function using posfun to prevent negative growth increments */
        double epsNegGrowth;
        
        /** labels for options for penalties on M2M parameters or ogives smoothness */
        adstring_array optsPenSmthPrM2M;
        /** integer indicating option for penalty on M2M parameters or ogives smoothness */
        int optPenSmthPrM2M;
        /** weight for penalties on M2M parameters or ogives smoothness */
        dvector wgtPenSmthPrM2M;      
        /** labels for options for penalties on non-decreasing M2M parameters or ogives */
        adstring_array optsPenNonDecPrM2M;
        /** integer indicating option for penalty on non-decreasing M2M parameters or ogives */
        int optPenNonDecPrM2M;
        /** weight for penalties on maturity non-decreasing M2M parameters or ogives */
        dvector wgtPenNonDecPrM2M;    
        
        /** labels for options for smoothing penalties on nonparametric sel functions */
        adstring_array optsPenSmthNPSel;
        /** integer indicating option for smoothing penalties on nonparametric sel functionss */
        int optPenSmthNPSel;
        /** weight for smoothing penalties on nonparametric sel functions */
        dvector wgtPenSmthNPSel;      
        
        /** pointer to empirical selectivity functions object */
        EmpiricalSelFcns* ptrEmpiricalSelFcns;
        
        /** pointer to empirical selectivity function priors object */
        EmpiricalSelFcnPriors* ptrEmpiricalSelFcnPriors;
        
        /** pointer to effort extrapolation scenarios object */
        EffXtrapScenarios* ptrEffXtrapScenarios;
        
        /** likelihood penalty weight (as CVs) for F-devs */
        double cvFDevsPen;              
        /** phase to start decreasing penalty weights on F-devs */
        int phsDecrFDevsPen;            
        /** phase at which to turn off penalties on F-devs */
        int phsZeroFDevsPen;     
        
        /** likelihood penalty weight for squared sum of devs */
        double wgtSqSumDevsPen;          
        /** phase to start the penalties on squared sum of devs */
        int phsSqSumDevsPen;             
        
        /** options for OFL calculations regarding capture rate averaging */
        adstring_array optsOFLAvgCapRate;
        /** integer vector indicating option for OFL capture rate averaging, by fishery */
        ivector optOFLAvgCapRate;
        /** integer vector indicating number of years for OFL capture rate averaging, by fishery */
        ivector oflNumYrsForAvgCapRate;
        /** average capture rate info for OFL calculations*/
        dvector oflAvgCapRateInfo;
        
        /** options for iterative re-weighting of size compositions */
        adstring_array optsIterativeReweighting;
        /** type of iterative re-weighting for size compositions */
        int optIterativeReweighting;
        /** phase to start iterative re-weighting for size compositions */
        int phsIterativeReweighting;
        /** max number of iterations */
        int maxIterations;
        
        /** options for MSEs*/
        adstring_array optsMSE;
        /** min year for recruitment statistics for OpMod */
        int opModRecStatsMinYr;
        /** max year for recruitment statistics for OpMod */
        int opModRecStatsMaxYr;
        /** harvest control rule to use */
        int HCR; 
        /** min year for averaging in HCR 1*/
        int HCR1_avgMinYr;
        /** max year for averaging in HCR 1*/
        int HCR1_avgMaxYr;
         /** id of ramp to use in HCR 2 */
        int HCR2_rampID;
        
        /** pointer to projection options for non-MCMC runs */
        ProjectionOptions* ptrProjOpts;
        /** pointer to projection options for MCMC runs */
        ProjectionOptions* ptrProjOptsMCMC;
        
        /** pointer to simulation options */
        SimOptions* ptrSimOpts;
 
    public:
        /**
         * Class constructor
         * 
         * @param mc - ModelConfiguration object
         */
        ModelOptions(ModelConfiguration& mc);
        /**
         * Class destructor.
         */
        ~ModelOptions(){}
        
        /**
         * Read from input stream in ADMB format.
         * 
         * @param is - input stream
         */
        void read(cifstream & is);        //read file in ADMB format
        /**
         * Write to output stream in ADMB format
         * 
         * @param os - output stream
         */
        void write(std::ostream & os);         //write object to file in ADMB format
        /**
         * Write object to R file as a list.
         * 
         * @param os - output stream to write to
         * @param nm - name for R object
         * @param indent - number of tabs to indent
         */
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read from input filestream in ADMB format
         */
        friend cifstream&    operator >>(cifstream & is, ModelOptions & obj){obj.read(is);return is;}
        /**
         * Operator to write to output stream in ADMB format
         */
        friend std::ostream& operator <<(std::ostream & os,   ModelOptions & obj){obj.write(os);;return os;}
    };


#endif /* MODELOPTIONS_HPP */

