/* 
 * File:   ModelConfiguration.hpp
 * Author: william.stockhausen
 *
 * Created on March 12, 2013, 8:59 AM
 * 
 * Includes:
 *  class ModelConfiguration
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
 * 2017-02-27: 1. Extracted ModelOptions class to ModelOptions.hpp
 * 2020-01-29: 1. Added yRetro 
 * 2020-07-08: 1. Added mnYrAvgRec, mxYrOffsetAvgRec.
 *             2. Incremented VERSION to '2020.07.08'.
 * 2021-04-09: 1. Added maxZs dvector for sex-specific max sizes
 *             2. Incremented VERSION to '2021.04.09'.
 * 2021-04-10: 1. Added maxZs dvector for sex-specific max sizes
 *             2. Incremented VERSION to '2021.04.10'.
 * 2023-03-07: 1. Replaced vector of model size bin cutpoints with min cutpoint and bin size
 *             2. Incremented VERSION to '2023.03.13'.
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
        
        /* vector of max sizes, by sex */
        static dvector maxZs;
        /* vector of indices to max size bin, by sex */
        static ivector maxZBs;
        /* max size possible at recruitment */
        static double maxZRec;
        /* index to max possible size bin for recruitment */
        static int maxZBRec;
        /* number of size bins */
        static int nZBs;
        /* maximum size bin cutpoint */
        static double maxZC;
        /* minimum size bin cutpoint */
        static double minZC;
        /* bin size */
        static double delZ;
        /* size bin midpoints (CW in mm) */
        static dvector zMidPts;
        /* size bin cutpoints (CW in mm) */
        static dvector zCutPts; 
        /* min model year */
        static int mnYr;//min model year
        /* assessment year (final pop numbers given for July 1, asYr) */
        static int asYr;//assessment year (final pop numbers given for July 1, asYr)
        /* max model year (mxYr = asYr-1) */
        static int mxYr;//max model year (mxYr = asYr-1)
        /* min year for OFL average recruitment calculation */
        static int mnYrAvgRec;
        /* offset to max year for OFL average recruitment calculation */
        static int mxYrOffsetAvgRec;
        /* number of retrospective years to incorporate */
        static int yRetro;
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
         * Set number of retrospective years to peel off.
         * 
         * @param _yRetro - number of retrospective years to peel off
         */
        void setNumRetroYears(int _yRetro){yRetro = _yRetro;}
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
    
#endif	/* MODELCONFIGURATION_HPP */

