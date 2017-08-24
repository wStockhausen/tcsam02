/* 
 * File:   ModelOptions.hpp
 * Author: WilliamStockhausen
 *
 * Created on February 27, 2017, 4:24 AM
 */

#ifndef MODELOPTIONS_HPP
#define MODELOPTIONS_HPP

#include <admodel.h>
#include "ModelConfiguration.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelParameterInfoTypes.hpp"

/*------------------------------------------------------------------------------\n
 * EffAvgScenario: class encapsulating effort averaging information \n
 *----------------------------------------------------------------------------*/
class EffAvgScenario{
    public:
        static int debug;
    private:
        ModelConfiguration* ptrMC;
    public:
        int id;
        int f;//fishery index
        IndexBlock* ptrIB;//pointer to an IndexBlock object
    private:
        adstring_array strVals;//names for index variables
    public:
        EffAvgScenario(ModelConfiguration& mc);
        ~EffAvgScenario();
        
        /**
         * Get the time block as an R array dimension
         * @return 
         */
        adstring getYDimsForR(){return ptrIB->getAsRDim();}
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EffAvgScenario & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EffAvgScenario & obj){obj.write(os); return os;}
};

/*------------------------------------------------------------------------------\n
 * EffAvgOptions: scenarios for effort averaging \n
 *----------------------------------------------------------------------------*/
class EffAvgScenarios{
    public:
        static int debug;
    private:
        ModelConfiguration* ptrMC;
    public:
        int nAvgs; //number of effort averaging info objects
        EffAvgScenario** ppEASs;//pointer to array of EffAvgInfo pointers 
    public:
        EffAvgScenarios(ModelConfiguration& mc);
        ~EffAvgScenarios();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, EffAvgScenarios & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, EffAvgScenarios & obj){obj.write(os); return os;}
};

/*------------------------------------------------------------------------------\n
 * CapRateAvgScenario: fishery capture rate averaging info \n
 *----------------------------------------------------------------------------*/
class CapRateAvgScenario{
    public:
        static int debug;
    public:
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

/*------------------------------------------------------------------------------\n
 * CapRateAvgScenarios: scenarios for capture rate averaging \n
 *----------------------------------------------------------------------------*/
class CapRateAvgScenarios{
    public:
        static int debug;
    public:
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

/*------------------------------------------------------------------------------\n
 * EffXtrapScenarios: effort extrapolation scenarios\n
 *----------------------------------------------------------------------------*/
class EffXtrapScenarios{
    public:
        static int debug;
    private:
        ModelConfiguration* ptrMC;
    public:
        EffAvgScenarios* ptrEffAvgScenarios;//pointer to effort averaging scenarios object 
        CapRateAvgScenarios* ptrCapRateAvgScenarios;//pointer to capture rate averaging scenarios object 
    public:
        EffXtrapScenarios(ModelConfiguration& mc);
        ~EffXtrapScenarios();
        
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
        /* version string for class */
        const static adstring VERSION;
        /* flag to print debug info */
        static int debug;  //flag to print debug info        
    public:
        /* pointer to model configuration object */
        ModelConfiguration* ptrMC;      //pointer to model configuration object
        
        /* labels for initial n-at-z options */
        adstring_array optsInitNatZ;
        /* selected option for initial n-at-z calculations */
        int optInitNatZ;               
        
        /* labels for options for penalties on non-decreasing M2M parameters/ogives */
        adstring_array optsParamNM;
        /* integer indicating option for natural mortality parameterization */
        int optParamNM;
        
        /* labels for growth parameterization options */
        adstring_array optsGrowthParam;  
        /* selected option for growth parameterization */
        int optGrowthParam;                 
        
        /* labels for growth pdf options */
        adstring_array optsGrowthPDF;  
        /* selected option for growth pdf */
        int optGrowthPDF;                 
        
        /* labels for options for penalties on M2M parameters/ogives smoothness */
        adstring_array optsPenSmthPrM2M;
        /* integer indicating option for penalty on M2M parameters/ogives smoothness */
        int optPenSmthPrM2M;
        /* weight for penalties on M2M parameters/ogives smoothness */
        dvector wgtPenSmthPrM2M;      
        /* labels for options for penalties on non-decreasing M2M parameters/ogives */
        adstring_array optsPenNonDecPrM2M;
        /* integer indicating option for penalty on non-decreasing M2M parameters/ogives */
        int optPenNonDecPrM2M;
        /* weight for penalties on maturity non-decreasing M2M parameters/ogives */
        dvector wgtPenNonDecPrM2M;    
        
        /* pointer to effort extrapolation scenarios object */
        EffXtrapScenarios* ptrEffXtrapScenarios;
        
        /* likelihood penalties for F-devs */
        double cvFDevsPen;              
        /* phase to start decreasing fpenCV */
        int phsDecrFDevsPen;            
        /* phase at which to turn off penalties on F-devs */
        int phsZeroFDevsPen;     
        
        /* likelihood penalties for last dev in each devs vector */
        double wgtLastDevsPen;          
        /* phase to start the penalty on the last devs */
        int phsLastDevsPen;             
        
        /* options for OFL calculations regarding capture rate averaging */
        adstring_array optsOFLAvgCapRate;
        /* integer vector indicating option for OFL capture rate averaging, by fishery */
        ivector optOFLAvgCapRate;
        /* integer vector indicating number of years for OFL capture rate averaging, by fishery */
        ivector oflNumYrsForAvgCapRate;
        /* average capture rate info for OFL calculations*/
        dvector oflAvgCapRateInfo;

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

