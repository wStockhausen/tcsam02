/* 
 * File:   ModelOptions.hpp
 * Author: WilliamStockhausen
 *
 * Created on February 27, 2017, 4:24 AM
 */

#ifndef MODELOPTIONS_HPP
#define MODELOPTIONS_HPP

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
        
        /* labels for growth options */
        adstring_array optsGrowth;  
        /* selected option for growth calculations */
        int optGrowth;                 
        
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
        
        /* labels for effort extrapolation estimation */
        adstring_array optsEffXtrEst;       
        /* selected options (by fishery) for average capture rate estimation */
        ivector optEffXtrEst;               
        /* labels for effort extrapolation estimation */
        adstring_array optsEffXtrAvgFc;       
        /* selected options (by fishery) for average capture rate estimation */
        ivector optEffXtrAvgFc;               
        
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

