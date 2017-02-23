/* 
 * File:   ModelConfiguration.hpp
 * Author: william.stockhausen
 *
 * Created on March 12, 2013, 8:59 AM
 * 
 * Includes:
 *  class ModelConfiguration
 *  class ModelOptions
 * 
 * History:
 * 2014-10-30: 1. Removed nSXs, nSCs, nMSs as static variables that were set through
 *                the model configuration file. These should now be changed in ModelConstants.hpp
 *                and the model recompiled for different configurations.
 * 2014-12-03: 1. Added asYr (assessment year) as input, replacing mxYr. Now, mxYr = asYr-1.
 * 2015-03-01: 1. Added "dims" adstring variables to facilitate output to R.
 *             2. Added optsGrowth and optsInitNatZ options to ModelOptions.
 * 2015-05-12: 1. Added cvFDevsPen, phsDecrFDevsPen,phsZeroFDevsPen for F-devs penalties to ModelOptions
 * 2015-05-13: 1. Added phsLastDevPen, wgtLastDevPen to ModelOptions
 * 2015-05-26: 1. Added penWgtSmthLgtPrMat, penWgtNonDecLgtPrMat to ModelOptions
 * 2016-04-13: 1. Added optPenNonDecLgtPrMat to ModelOptions
 *             2. Added version strings to ModelConfiguration, ModelOptions
 *             3. Updated documentation
 */

#ifndef MODELCONFIGURATION_HPP
    #define	MODELCONFIGURATION_HPP

//--------------------------------------------------------------------------------
//          ModelConfiguration
//--------------------------------------------------------------------------------
    class ModelConfiguration {
    public:
        /* version string for class */
        const static adstring VERSION;
        
        static int debug;  //flag to print debug info
        
        /* numbr of size bins */
        static int nZBs;//number of size bins 
        /* min model year */
        static int mnYr;//min model year
        /* assessment year (final pop numbers given for July 1, asYr) */
        static int asYr;//assessment year (final pop numbers given for July 1, asYr)
        /* max model year (mxYr = asYr-1) */
        static int mxYr;//max model year (mxYr = asYr-1)
        /* number of fisheries */
        static int nFsh;//number of fisheries 
        /* number of surveys */
        static int nSrv;//number of surveys 
        
        /* flag to jitter initial parameter values */
        static int    jitter;  //flag to jitter initial parameter values
        /* fraction to jitter bounded parameter values */
        static double jitFrac; //fraction to jitter bounded parameter values
        /* flag to resample initial parameter values */
        static int    resample;//flag to resample initial parameter values
        /* variance inflation factor for resampling parameter values */
        static double vif;     //variance inflation factor for resampling parameter values
    public:
        /* model configuration name */
        adstring cfgName;//model configuration name

        /* flag to run operating model only */
        int runOpMod;    //flag to run operating model only
        /* flag to fit to priors */
        int fitToPriors; //flag to fit to priors
        
        /* size bin midpoints (CW in mm) */
        dvector zMidPts;     //size bin midpoints (CW in mm)
        /* size bin cutpoints (CW in mm) */
        dvector zCutPts;     //size bin cutpoints (CW in mm)
        /* vector of 1's same size as zMidPts */
        dvector onesZMidPts; //vector of 1's same size as zMidPts

        /* model datasets file name */
        adstring fnMDS; //model datasets file name
        /* model parameters info file name */
        adstring fnMPI; //model parameters info file name
        /* model options file name */
        adstring fnMOs; //model options file name
        
        /* labels for fisheries */
        adstring_array lblsFsh;//labels for fisheries
        /* labels for surveys */
        adstring_array lblsSrv;//labels for surveys
        
        adstring csvYrs;  //csv string of model years (mnYr:mxYr)
        adstring csvYrsP1;//csv string of model years+1 (mnYr:mxYr+1)
        adstring csvSXs;//csv string of sexes
        adstring csvMSs;//csv string of maturity states
        adstring csvSCs;//csv string of shell conditions
        adstring csvZCs;//csv string of size bin cutpoints
        adstring csvZBs;//csv string of size bin midpoints
        adstring csvFsh;//csv string of fishery labels
        adstring csvSrv;//csv string of survey labels

        adstring dimYrsToR;  //R dim string of model years (mnYr:mxYr)
        adstring dimYrsP1ToR;//R dim string of model years+1 (mnYr:mxYr+1)
        adstring dimSXsToR;//R dim string of sexes
        adstring dimMSsToR;//R dim string of maturity states
        adstring dimSCsToR;//R dim string of shell conditions
        adstring dimZCsToR;//R dim string of size bin cutpoints
        adstring dimZBsToR;//R dim string of size bin midpoints
        adstring dimZPsToR;//R dim string of size bin midpoints (alternative)
        adstring dimFshToR;//R dim string of fishery labels
        adstring dimSrvToR;//R dim string of survey labels

    public:
        ModelConfiguration();
        ~ModelConfiguration();
        ModelConfiguration& operator =(const ModelConfiguration & rhs);
        
        /**
         * Set max model year (for retrospective model runs).
         * 
         * @param yr - new max model year
         */
        void setMaxModelYear(int yr);
        /**
         * Tests if mnYr <= yr <= mxYr. 
         * 
         * @param yr
         * @return 1 if true, 0 if false
         */
        int isModelYear(int yr){if ((mnYr<=yr)&&(yr<=mxYr)) return 1; return 0;}
        /**
         * Read input file in ADMB format.
         * 
         * @param fn - name of file to read
         */
        void read(const adstring & fn);   //read file in ADMB format
        /**
         * Read from input file stream in ADMB format.
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);        //read file in ADMB format
        /**
         * Write object to file in ADMB format.
         * 
         * @param fn - name of file to write
         */
        void write(const adstring & fn);  //write object to file in ADMB format
        /**
         * Write object to output stream in ADMB format.
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
        void writeToR(std::ostream& os, std::string nm, int indent=0);
        /**
         * Operator to read from input filestream in ADMB format
         */
        friend cifstream&    operator >>(cifstream & is, ModelConfiguration & obj){obj.read(is);return is;}
        /**
         * Operator to write to output stream in ADMB format
         */
        friend std::ostream& operator <<(std::ostream & os,   ModelConfiguration & obj){obj.write(os);;return os;}
    };

//--------------------------------------------------------------------------------
//          ModelOptions
//--------------------------------------------------------------------------------
    class ModelOptions {
    public:
        /* version string for class */
        const static adstring VERSION;
        /* flag to print debug info */
        static int debug;  //flag to print debug info        
    public:
        /* pointer to model configuration object */
        ModelConfiguration* ptrMC;      //pointer to model configuration object
        /* labels for capture rate averaging options */
        adstring_array optsFcAvg;       //labels for capture rate averaging options
        /* selected options for averaging capture rate */
        ivector optFcAvg;               //selected options for averaging capture rate
        /* labels for growth options */
        adstring_array optsGrowth;  
        /* selected option for growth calculations */
        int optGrowth;                 
        /* labels for initial n-at-z options */
        adstring_array optsInitNatZ;
        /* selected option for initial n-at-z calculations */
        int optInitNatZ;               
        /* penalty for F-devs */
        double cvFDevsPen;              
        /* phase to start decreasing fpenCV */
        int phsDecrFDevsPen;            
        /* phase at which to turn off penalties on F-devs */
        int phsZeroFDevsPen;            
        /* penalty for last dev in each devs vector */
        double wgtLastDevsPen;          
        /* phase to start the penalty on the last devs */
        int phsLastDevsPen;             
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
        /* labels for options for penalties on non-decreasing M2M parameters/ogives */
        adstring_array optsParamNM;
        /* integer indicating option for natural mortality parameterization */
        int optParamNM;

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
    
#endif	/* MODELCONFIGURATION_HPP */

