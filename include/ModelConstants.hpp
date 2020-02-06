/* 
 * File:   ModelConstants.hpp
 * Author: william.stockhausen
 * 
 * History:
 * 20140605: added tcsam::dgbAll, dbgPriors
 * 20140918: renamed string constants of form "name_STR" to follow format "STR_name"
 * 20140930: switched MALE, FEMALE so MALE=1, FEMALE=2
 * 20141030: switched ANY to ALL
 * 20150113: added FIT_BY_XSE, FIT_BY_XMSE
 * 20150302: added RSTR_... constants and converter
 * 20150417: revised STR_FIT_BY_ string values
 * 20150518: revised FIT_BY_ and STR_FIT_BY_ values
 * 20160413: added tcsam::VERSION string to indicate model version
 * 20160424: added tcsam:FIT_BY_X_MATONLY to fit by sex for mature crab only
 * 20170206: incremented tcsam::VERSION to "2017.02.06"
 * 20170606: incremented tcsam::VERSION to "2017.06.06"
 * 20171025: incremented tcsam::VERSION to "2017.10.25"
 * 20171116: incremented tcsam::VERSION to "2017.11.16"
 * 20171205: incremented tcsam::VERSION to "2017.12.05"
 * 20171206: incremented tcsam::VERSION to "2017.12.06"
 * 20180411: incremented tcsam::VERSION to "2018.04.11"
 *           added LL_GAMMA and STR_LL_GAMMA
 * 20180824: incremented tcsam::VERSION to "2018.08.24"
 *           added FIT_BY_X_MSE string and int constants
 * 20180827: incremented tcsam::VERSION to "2018.08.27"
 *           reflecting changes to tpl file.
 * 20180829: incremented tcsam::VERSION to "2018.08.29"
 *           reflecting changes to tpl file.
 * 20180902: incremented tcsam::VERSION to "2018.09.02"
 *           reflecting changes to tpl file.
 * 20200202: incremented tcsam::VERSION to "2020.02.02"
 *           reflecting changes to tpl file.
 */

#pragma once
#ifndef MODELCONSTANTS_HPP
    #define	MODELCONSTANTS_HPP
   
class rpt{
    public:
        /* Global output filestream */
        static std::ofstream echo;
};

class tcsamDims{
    public:
        static adstring formatForR(const adstring& s);
        static adstring getSXsForR(int mn,int mx);
        static adstring getMSsForR(int mn,int mx);
        static adstring getSCsForR(int mn,int mx);
};

namespace tcsam{
    /* adstring indicating model name */
    const adstring MODEL = "tcsam02";
    /* adstring indicating model version */
    const adstring VERSION = "2020.02.02";
    
    /* minimum debugging level that will print ALL debug info */
    const int dbgAll = 100;
    /* minimum debugging level that will print debug info for prior calcs */
    const int dbgPriors = 30;
    
    /* Model constant for generic "ALL" */
    const adstring STR_ALL = "ALL";
    /* Model constant for generic "NONE" */
    const adstring STR_NONE = "NONE";
    
    /* Model dimension name for sex */
    const adstring STR_SEX = "SEX";
    /* Model dimension name for maturity state */
    const adstring STR_MATURITY_STATE = "MATURITY_STATE";
    /* Model dimension name for shell condition */
    const adstring STR_SHELL_CONDITION = "SHELL_CONDITION";
    /* Model dimension name for size (bins) */
    const adstring STR_SIZE = "SIZE";
    /* Model dimension name for year */
    const adstring STR_YEAR = "YEAR";
    /* Model dimension name for fisheries */
    const adstring STR_FISHERY = "FISHERY";
    /* Model dimension name for surveys */
    const adstring STR_SURVEY = "SURVEY";
    /* Model flag name for selectivity functions */
    const adstring STR_SELFCN = "selFcn";
    
    //sexes
    /* number of sexes in model */
    const int nSXs      = 2;//number of model sexes
    /* integer indicating sex is male */
    const int MALE      = 1;//integer indicating sex=male
    /* integer indicating sex is female */
    const int FEMALE    = 2;//integer indicating sex=female
    /* integer indicating all sexes combined */
    const int ALL_SXs   = nSXs+1*(nSXs>1);//integer indicating all sexes combined
    /* adstring indicating sex is male */
    const adstring STR_MALE    = "MALE";   //string indicating male
    /* adstring indicating sex is female */
    const adstring STR_FEMALE  = "FEMALE"; //string indicating female
    /* adstring indicating all sexes combined */
    const adstring STR_ALL_SXs = "ALL_SEX";//string indicating all sexes combined
    
    //maturity states
    const int nMSs     = 2; //number of model maturity states
    const int IMMATURE = 1; //integer indicating immature
    const int MATURE   = 2; //integer indicating mature
    const int ALL_MSs = nMSs+1*(nMSs>1); //integer indicating all maturity states combined
    const adstring STR_IMMATURE = "IMMATURE";//string indicating immature
    const adstring STR_MATURE   = "MATURE";  //string indicating mature
    const adstring STR_ALL_MSs  = "ALL_MATURITY"; //string indicating all maturity states combined
        
    //shell conditions
    const int nSCs = 2;     //number of model shell conditions
    const int NEW_SHELL = 1;//integer indicating new shell condition
    const int OLD_SHELL = 2;//integer indicating old shell condition
    const int ALL_SCs = nSCs+1*(nSCs>1);//integer indicating all shell conditions combined
    const adstring STR_NEW_SHELL = "NEW_SHELL"; //string indicating new shell condition
    const adstring STR_OLD_SHELL = "OLD_SHELL"; //string indicating old shell condition
    const adstring STR_ALL_SCs   = "ALL_SHELL"; //string indicating all shell conditions combined
    
    //objective function fitting option types
    const adstring STR_FIT_NONE    = "NONE";
    const adstring STR_FIT_BY_TOT  = "BY_TOTAL";
    const adstring STR_FIT_BY_X    = "BY_X";
    const adstring STR_FIT_BY_XM   = "BY_XM";
    const adstring STR_FIT_BY_XS   = "BY_XS";
    const adstring STR_FIT_BY_XMS  = "BY_XMS";
    const adstring STR_FIT_BY_XE   = "BY_XE";
    const adstring STR_FIT_BY_X_ME = "BY_X_ME";
    const adstring STR_FIT_BY_X_SE = "BY_X_SE";
    const adstring STR_FIT_BY_XME  = "BY_XME";
    const adstring STR_FIT_BY_XM_SE = "BY_XM_SE";
    const adstring STR_FIT_BY_X_MSE = "BY_X_MSE";
    const adstring STR_FIT_BY_X_MATONLY = "BY_X_MATONLY";
    const int FIT_NONE    = 0;
    const int FIT_BY_TOT  = 1;
    const int FIT_BY_X    = 2;
    const int FIT_BY_XM   = 3;
    const int FIT_BY_XS   = 4;
    const int FIT_BY_XMS  = 5;
    const int FIT_BY_XE   = 6;
    const int FIT_BY_X_ME = 7;
    const int FIT_BY_X_SE = 8;
    const int FIT_BY_XME  = 9;
    const int FIT_BY_XM_SE = 10;
    const int FIT_BY_X_MATONLY = 11;
    const int FIT_BY_X_MSE = 12;
    
    /** adstring constant indicating likelihood type 'NONE' */
    const adstring STR_LL_NONE        = "NONE";
    /** adstring constant indicating likelihood type 'NORM2' */
    const adstring STR_LL_NORM2       = "NORM2";
    /** adstring constant indicating likelihood type 'NORMAL' */
    const adstring STR_LL_NORMAL      = "NORMAL";
    /** adstring constant indicating likelihood type 'LOGNORMAL' */
    const adstring STR_LL_LOGNORMAL   = "LOGNORMAL";
    /** adstring constant indicating likelihood type 'MULTINOMIAL' */
    const adstring STR_LL_MULTINOMIAL = "MULTINOMIAL";
    /** adstring constant indicating likelihood type 'BINOMIAL' */
    const adstring STR_LL_BINOMIAL    = "BINOMIAL";
    /** adstring constant indicating likelihood type 'GAMMA' */
    const adstring STR_LL_GAMMA    = "GAMMA";
    /** integer constant indicating likelihood type 'NONE' */
    const int LL_NONE        = 0;
    /** integer constant indicating likelihood type 'NORM2' */
    const int LL_NORM2       = 1;
    /** integer constant indicating likelihood type 'NORMAL' */
    const int LL_NORMAL      = 2;
    /** integer constant indicating likelihood type 'LOGNORMAL' */
    const int LL_LOGNORMAL   = 3;
    /** integer constant indicating likelihood type 'MULTINOMIAL' */
    const int LL_MULTINOMIAL = 4;
    /** integer constant indicating likelihood type 'BINOMIAL' */
    const int LL_BINOMIAL    = 5;
    /** integer constant indicating likelihood type 'GAMMA' */
    const int LL_GAMMA    = 6;
    
    //Stock-recruit function types
    const adstring STR_CONSTANT = "CONSTANT";
    const adstring STR_BEVHOLT  = "BEVHOLT";
    const adstring STR_RICKER   = "RICKER";
    const int SRTYPE_CONSTANT = 1;//integer indicating a constant recruitment SRR
    const int SRTYPE_BEVHOLT  = 2;//integer indicating a Beverton-Holt SRR
    const int SRTYPE_RICKER   = 3;//integer indicating a Ricker SRR
    
    const adstring STR_VAR = "VARIANCE";
    const adstring STR_STD = "STD_DEV";
    const adstring STR_CV  = "CV";
    const int SCLTYPE_VAR = 0;//integer indicating variances are given
    const int SCLTYPE_STD = 1;//integer indicating std devs are given
    const int SCLTYPE_CV  = 2;//integer indicating cv's are given
    
    //units
    const adstring UNITS_ONES      = "ONES";//"ONES"
    const adstring UNITS_THOUSANDS = "THOUSANDS";//"THOUSANDS"
    const adstring UNITS_MILLIONS  = "MILLIONS";//"MILLIONS"
    const adstring UNITS_BILLIONS  = "BILLIONS";//"BILLIONS"
    const adstring UNITS_GM        = "GM";//"GM"
    const adstring UNITS_KG        = "KG";//"KG"
    const adstring UNITS_MT        = "MT";//"MT"
    const adstring UNITS_KMT       = "THOUSANDS_MT";//"THOUSANDS_MT"
    const adstring UNITS_LBS       = "LBS";//"LBS"
    const adstring UNITS_MLBS      = "MILLIONS_LBS";//"MILLIONS_LBS"
    //unit conversions
    const double CONV_KGtoLBS = 2.20462262; //multiplier conversion from kg to lbs
    
    //transform types from arithmetic space to parameter space
    /* alternative adstring constant indicating parameter scale is arithmetic */
    const adstring STR_SCALE_NONE      = "NONE";
    /* adstring constant indicating parameter values are on the arithmetic (untransformed) scale */
    const adstring STR_SCALE_ARITHM    = "ARITHMETIC";
    /* adstring constant indicating parameter values are on the logit scale */
    const adstring STR_SCALE_LOGIT     = "LOGIT";
    /* adstring constant indicating parameter values are on the log scale */
    const adstring STR_SCALE_LOG       = "LOG";
    /* adstring constant indicating parameter values are on the probit scale */
    const adstring STR_SCALE_PROBIT    = "PROBIT";
    /* int constant indicating parameter values are on the arithmetic (untransformed) scale */
    const int SCALE_ARITHM    = 0;
    /* int constant indicating parameter values are on the logit scale */
    const int SCALE_LOGIT     = 1;
    /* int constant indicating parameter values are on the natural log scale */
    const int SCALE_LOG       = 2;
    /* int constant indicating parameter values are on the probit scale */
    const int SCALE_PROBIT    = 3;
    
    int getMaturityType(adstring s);
    adstring getMaturityType(int i);

    int getSexType(adstring s);
    adstring getSexType(int i);

    int getShellType(adstring s);
    adstring getShellType(int i);

    //Stock-recruit function types
    int getSRType(adstring s);
    adstring getSRType(int i);

    /**
     * Translate error scale type from adstring  to int version.
     * 
     * @param fitType
     * @return 
     */
    int getErrorScaleType(adstring sclType);
    /**
     * Translate error scale type from int to adstring version.
     * 
     * @param fitType
     * @return 
     */
    adstring getErrorScaleType(int sclFlg);
    
    /**
     * Translate fit type from adstring to int version.
     * 
     * @param fitType
     * @return 
     */
    int getFitType(adstring fitType);
    
    /**
     * Translate fit type from integer to adstring version.
     * 
     * @param i
     * @return 
     */
    adstring getFitType(int i);

    double convertToStdDev(double sclVal, double mnVal, int sclFlg);
    dvector convertToStdDev(dvector sclVal, dvector mnVal, int sclFlg);
    
    /**
     * Gets multiplicative conversion factor from "from" units to "to" units. 
     * 
     * @param from - adstring UNITS_ keyword
     * @param to   - adstring UNITS_ keyword
     * @return - multiplicative factor: to_units  = factor*from_units
     */
    double getConversionMultiplier(adstring from,adstring to);
        
    /**
     * Translate from adstring likelihood type to int version.
     * 
     * @param llType - adstring likelihood type
     * @return - corresponding integer version
     */
    int getLikelihoodType(adstring llType);
    /**
     * Translate from integer likelihood type to adstring version.
     * 
     * @param llType - integer likelihood type
     * @return - corresponding adstring version
     */
    adstring getLikelihoodType(int llType);
        
    /**
     * Translate from adstring scale type to int version.
     * 
     * @param sclType - adstring scale type
     * @return -integer scale type
     */
    int getScaleType(adstring sclType);
    /**
     * Translate from integer scale type to adstring version.
     * 
     * @param sclType - integer scale type
     * @return  - corresponding adstring version
     */
    adstring getScaleType(int sclType);
}


#endif	/* MODELCONSTANTS_HPP */

