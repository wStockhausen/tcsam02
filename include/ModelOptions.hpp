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

