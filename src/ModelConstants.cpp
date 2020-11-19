#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
/**
 * Changes:
 * 20141203: 1. Changed std::exit(-1) to exit(-1) for MinGW compatibility.
 * 20150113: 1. Added FIT_BY_XSE, FIT_BY_XMSE
 * 20150210: 1. Added conversion from UNITS_KMT to other units
 * 20150302: 1. Added tcsamDims::formatForR(adstring).
 *           2. Revised getDDsForR(...) functions (DD=SX,SC,MS) to use formatForR(...)
 * 20171024: 1. Added transform types
 * 20201117: 1. Added LL_MULINOMIALP, LL_DIRICHLET, LL_DIRICHLETP
*/

using namespace tcsam;

std::ofstream rpt::echo("EchoData.dat",std::ios::trunc);

/**
 * Format for R output: lower case, replace "_"'s with spaces
 * @param s
 * @return 
 */
adstring tcsamDims::formatForR(const adstring& s){
//    rpt::echo<<"In formatForR() with '"<<s<<"'"<<endl;
    adstring tmp; tmp = s; tmp.to_lower();
//    rpt::echo<<"tmp = "<<tmp<<endl;
    int p = tmp.pos('_');
//    rpt::echo<<"p = "<<p<<endl;
    while (p>0){
        tmp(p)=' ';
        p = tmp.pos('_');
//        rpt::echo<<"p = "<<p<<endl;
    }
//    rpt::echo<<"original: '"<<s<<"'. Re-formatted: '"<<tmp<<"'"<<endl;
    return tmp;
}

/**
 * Convert a vector of model sexes (defined by the input indices) to an
 * R character vector formatted as part of an R list, using 'x'
 * as the name.
 * 
 * @param mn
 * @param mx
 * @return 
 */
adstring tcsamDims::getSXsForR(int mn,int mx){
    adstring dms = "x=c("+qt+formatForR(getSexType(mn))+qt;
    for (int i=(mn+1);i<=mx;i++) {
        dms += cc+qt+formatForR(getSexType(i))+qt;
    }
    dms += ")";
    return dms;
}

/**
 * Convert a vector of model maturity states (defined by the input indices) to an
 * R character vector formatted as part of an R list, using 'm'
 * as the name.
 * 
 * @param mn
 * @param mx
 * @return 
 */
adstring tcsamDims::getMSsForR(int mn,int mx){
    adstring dms = "m=c("+qt+formatForR(getMaturityType(mn))+qt;
    for (int i=(mn+1);i<=mx;i++) {
        dms += cc+qt+formatForR(getMaturityType(i))+qt;
    }
    dms += ")";
    return dms;
}

/**
 * Convert a vector of model shell conditions (defined by the input indices) to an
 * R character vector formatted as part of an R list, using 's'
 * as the name.
 * 
 * @param mn
 * @param mx
 * @return 
 */
adstring tcsamDims::getSCsForR(int mn,int mx){
    adstring dms = "s=c("+qt+formatForR(getShellType(mn))+qt;
    for (int i=(mn+1);i<=mx;i++) {
        dms += cc+qt+formatForR(getShellType(i))+qt;
    }
    dms += ")";
    return dms;
}

int tcsam::getMaturityType(adstring s){
    s.to_upper();
    if (s==STR_IMMATURE) return IMMATURE; else
    if (s==STR_MATURE)   return MATURE;   else
    if (s==STR_ALL_MSs)  return ALL_MSs;  else
    if (s==STR_ALL)      return ALL_MSs;  else
    std::cout<<"Unrecognized MaturityType '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getMaturityType(int i){
    if (i<=ALL_SXs){
        if (i==IMMATURE) return STR_IMMATURE; else
        if (i==MATURE)   return STR_MATURE;   else
        if (i==ALL_MSs)  return STR_ALL_MSs; 
    }
    std::cout<<"Unrecognized or inappropriate MaturityType '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getSexType(adstring s){
    s.to_upper();
    if (s==STR_MALE)    return MALE;    else
    if (s==STR_FEMALE)  return FEMALE;  else
    if (s==STR_ALL_SXs) return ALL_SXs; else
    if (s==STR_ALL)     return ALL_SXs; else
    std::cout<<"Unrecognized SexType '"<<s<<"'"<<std::endl;
    std::cout<<"Potential values are:"<<endl;
    std::cout<<"'"<<STR_MALE<<"'"<<std::endl;
    std::cout<<"'"<<STR_FEMALE<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getSexType(int i){
    if (i<=ALL_SXs){
        if (i==MALE)    return STR_MALE;    else
        if (i==FEMALE)  return STR_FEMALE;  else
        if (i==ALL_SXs) return STR_ALL_SXs; 
    } 
    std::cout<<"Unrecognized or inappropriate SexType '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getShellType(adstring s){
    s.to_upper();
    if (s==STR_NEW_SHELL) return NEW_SHELL; else
    if (s==STR_OLD_SHELL) return OLD_SHELL; else
    if (s==STR_ALL_SCs)   return ALL_SCs;   else
    if (s==STR_ALL)       return ALL_SCs;   else
    std::cout<<"Unrecognized ShellType '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getShellType(int i){
    if (i<=ALL_SCs){
        if (i==NEW_SHELL) return STR_NEW_SHELL; else
        if (i==OLD_SHELL) return STR_OLD_SHELL; else
        if (i==ALL_SCs)   return STR_ALL_SCs;
    }
    std::cout<<"Unrecognized or inappropriate ShellType '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

int tcsam::getSRType(adstring s){
    s.to_lower();
    if (s==STR_CONSTANT) return SRTYPE_CONSTANT;
    if (s==STR_BEVHOLT)  return SRTYPE_BEVHOLT;
    if (s==STR_RICKER)   return SRTYPE_RICKER;
    std::cout<<"Unrecognized Stock-Recruit type (SRType) '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
adstring tcsam::getSRType(int i){
    if (i==SRTYPE_CONSTANT) return STR_CONSTANT;
    if (i==SRTYPE_BEVHOLT)  return STR_BEVHOLT;
    if (i==SRTYPE_RICKER)   return STR_RICKER;
    std::cout<<"Unrecognized Stock-Recruit type (SRType) '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}
/**
 * Converts the error scale type from adstring to integer version.
 * @param s - the adstring version
 * @return  - the integer version
 */
int tcsam::getErrorScaleType(adstring s){
    if (s==STR_VAR) return SCLTYPE_VAR;
    if (s==STR_STD) return SCLTYPE_STD;
    if (s==STR_CV)  return SCLTYPE_CV;
    std::cout<<"Unrecognized scale type '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return 0;
}
/**
 * Converts the error scale type from integer to adstring
 * @param i - the integer type
 * @return  - the adstring version
 */
adstring tcsam::getErrorScaleType(int i){
    if (i==SCLTYPE_VAR) return STR_VAR;
    if (i==SCLTYPE_STD) return STR_STD;
    if (i==SCLTYPE_CV)  return STR_CV;
    std::cout<<"Unrecognized scale type '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

/**
 * Translate from adstring fit type to int version.
 * @param s - likelihood fit type
 * @return 
 */
int tcsam::getFitType(adstring s){
    if (s==STR_FIT_NONE)    return FIT_NONE;
    if (s==STR_FIT_BY_TOT)  return FIT_BY_TOT;
    if (s==STR_FIT_BY_X)    return FIT_BY_X;
    if (s==STR_FIT_BY_XM)   return FIT_BY_XM;
    if (s==STR_FIT_BY_XS)   return FIT_BY_XS;
    if (s==STR_FIT_BY_XMS)  return FIT_BY_XMS;
    if (s==STR_FIT_BY_XE)   return FIT_BY_XE;
    if (s==STR_FIT_BY_X_ME) return FIT_BY_X_ME;
    if (s==STR_FIT_BY_X_SE) return FIT_BY_X_SE;
    if (s==STR_FIT_BY_XME)  return FIT_BY_XME;
    if (s==STR_FIT_BY_XM_SE) return FIT_BY_XM_SE;
    if (s==STR_FIT_BY_X_MSE) return FIT_BY_X_MSE;
    if (s==STR_FIT_BY_X_MATONLY) return FIT_BY_X_MATONLY;
    std::cout<<"Unrecognized fit type '"<<s<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return -1;
}

/**
 * Translate from int fit type to adstring version.
 * @param fitType
 * @return 
 */
adstring tcsam::getFitType(int i){
    if (i==FIT_NONE)    return STR_FIT_NONE;
    if (i==FIT_BY_TOT)  return STR_FIT_BY_TOT;
    if (i==FIT_BY_X)    return STR_FIT_BY_X;
    if (i==FIT_BY_XM)   return STR_FIT_BY_XM;
    if (i==FIT_BY_XS)   return STR_FIT_BY_XS;
    if (i==FIT_BY_XMS)  return STR_FIT_BY_XMS;
    if (i==FIT_BY_XE)   return STR_FIT_BY_XE;
    if (i==FIT_BY_X_ME) return STR_FIT_BY_X_ME;
    if (i==FIT_BY_X_SE) return STR_FIT_BY_X_SE;
    if (i==FIT_BY_XME)  return STR_FIT_BY_XME;
    if (i==FIT_BY_XM_SE) return STR_FIT_BY_X_MSE;
    if (i==FIT_BY_X_MSE) return STR_FIT_BY_XM_SE;
    if (i==FIT_BY_X_MATONLY) return STR_FIT_BY_X_MATONLY;
    std::cout<<"Unrecognized fit type '"<<i<<"'"<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return adstring("");
}

/***************************************************************
*   Converts from variance type indicated by sclFlg to         *
*   standard deviation.                                        *
***************************************************************/
double tcsam::convertToStdDev(double sclVal, double mnVal, int sclFlg){
    double sdv = 0.0;
    if (sclFlg==SCLTYPE_CV) {
        sdv = sclVal*mnVal;
    } else if (sclFlg==SCLTYPE_STD) {
        sdv = sclVal;
    } else if (sclFlg==SCLTYPE_VAR) {
        sdv = sqrt(sclVal);
    }
    return sdv;
}

/***************************************************************
*   Converts from variance type indicated by sclFlg to         *
*   standard deviation.                                        *
***************************************************************/
dvector tcsam::convertToStdDev(dvector sclVal, dvector mnVal, int sclFlg){
    dvector sdv = 0.0*sclVal;
    if (sclFlg==SCLTYPE_CV) {
        sdv = elem_prod(sclVal,mnVal);
    } else if (sclFlg==SCLTYPE_STD) {
        sdv = sclVal;
    } else if (sclFlg==SCLTYPE_VAR) {
        sdv = sqrt(sclVal);
    }
    return sdv;
}
    
/**
 * Gets multiplicative conversion factor from "from" units to "to" units. 
 * @param from - adstring UNITS_ keyword
 * @param to   - adstring UNITS_ keyword
 * @return - multiplicative factor: to_units  = factor*from_units
 */
double tcsam::getConversionMultiplier(adstring from,adstring to){
    double cnv = 1.0;
    if (from==UNITS_ONES){
        if (to==UNITS_ONES)      return 1.0E0;
        if (to==UNITS_THOUSANDS) return 1.0E-3;
        if (to==UNITS_MILLIONS)  return 1.0E-6;
        if (to==UNITS_BILLIONS)  return 1.0E-9;
    }
    if (from==UNITS_THOUSANDS){
        if (to==UNITS_ONES)      return 1.0E+3;
        if (to==UNITS_THOUSANDS) return 1.0E+0;
        if (to==UNITS_MILLIONS)  return 1.0E-3;
        if (to==UNITS_BILLIONS)  return 1.0E-6;
    }
    if (from==UNITS_MILLIONS){
        if (to==UNITS_ONES)      return 1.0E6;
        if (to==UNITS_THOUSANDS) return 1.0E3;
        if (to==UNITS_MILLIONS)  return 1.0E0;
        if (to==UNITS_BILLIONS)  return 1.0E-3;
    }
    if (from==UNITS_BILLIONS){
        if (to==UNITS_ONES)      return 1.0E-9;
        if (to==UNITS_THOUSANDS) return 1.0E-6;
        if (to==UNITS_MILLIONS)  return 1.0E-3;
        if (to==UNITS_BILLIONS)  return 1.0E0;
    }
    if (from==UNITS_GM){
        if (to==UNITS_GM)   return 1.0E+0;
        if (to==UNITS_KG)   return 1.0E-3;
        if (to==UNITS_MT)   return 1.0E-6;
        if (to==UNITS_KMT)  return 1.0E-9;
        if (to==UNITS_LBS)  return 1.0E-3 * CONV_KGtoLBS;
        if (to==UNITS_MLBS) return 1.0E-9 * CONV_KGtoLBS;
    }
    if (from==UNITS_KG){
        if (to==UNITS_GM)   return 1.0E+3;
        if (to==UNITS_KG)   return 1.0E+0;
        if (to==UNITS_MT)   return 1.0E-3;
        if (to==UNITS_KMT)  return 1.0E-6;
        if (to==UNITS_LBS)  return 1.0E+0 * CONV_KGtoLBS;
        if (to==UNITS_MLBS) return 1.0E-6 * CONV_KGtoLBS;
    }
    if (from==UNITS_MT){
        if (to==UNITS_GM)   return 1.0E+6;
        if (to==UNITS_KG)   return 1.0E+3;
        if (to==UNITS_MT)   return 1.0E+0;
        if (to==UNITS_KMT)  return 1.0E-3;
        if (to==UNITS_LBS)  return 1.0E+3 * CONV_KGtoLBS;
        if (to==UNITS_MLBS) return 1.0E-3 * CONV_KGtoLBS;
    }
    if (from==UNITS_KMT){
        if (to==UNITS_GM)   return 1.0E+9;
        if (to==UNITS_KG)   return 1.0E+6;
        if (to==UNITS_MT)   return 1.0E+3;
        if (to==UNITS_KMT)  return 1.0E+0;
        if (to==UNITS_LBS)  return 1.0E+6 * CONV_KGtoLBS;
        if (to==UNITS_MLBS) return 1.0E+0 * CONV_KGtoLBS;
    }
    if (from==UNITS_LBS){
        if (to==UNITS_GM)   return 1.0E+3/CONV_KGtoLBS;
        if (to==UNITS_KG)   return 1.0E+0/CONV_KGtoLBS;
        if (to==UNITS_MT)   return 1.0E-3/CONV_KGtoLBS;
        if (to==UNITS_KMT)  return 1.0E-6/CONV_KGtoLBS;
        if (to==UNITS_LBS)  return 1.0E+0;
        if (to==UNITS_MLBS) return 1.0E-6;
    }
    if (from==UNITS_MLBS){
        if (to==UNITS_GM)   return 1.0E+9/CONV_KGtoLBS;
        if (to==UNITS_KG)   return 1.0E+6/CONV_KGtoLBS;
        if (to==UNITS_MT)   return 1.0E+3/CONV_KGtoLBS;
        if (to==UNITS_KMT)  return 1.0E+0/CONV_KGtoLBS;
        if (to==UNITS_LBS)  return 1.0E+6;
        if (to==UNITS_MLBS) return 1.0E+0;
    }
    std::cout<<"Error converting units!"<<std::endl;
    std::cout<<"No conversion defined from '"<<from<<"' to '"<<to<<"'."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return cnv;
}

/**
 * Translate from integer likelihood type to adstring version.
 * 
 * @param llType - integer likelihood type
 * @return - corresponding adstring version
 */
adstring tcsam::getLikelihoodType(int llType){
    adstring type = STR_LL_NONE;
    if (llType==LL_BINOMIAL)     return STR_LL_BINOMIAL;
    if (llType==LL_GAMMA)        return STR_LL_GAMMA;
    if (llType==LL_LOGNORMAL)    return STR_LL_LOGNORMAL;
    if (llType==LL_MULTINOMIAL)  return STR_LL_MULTINOMIAL;
    if (llType==LL_MULTINOMIALP) return STR_LL_MULTINOMIALP;
    if (llType==LL_DIRICHLET)    return STR_LL_DIRICHLET;
    if (llType==LL_DIRICHLETP)   return STR_LL_DIRICHLETP;
    if (llType==LL_NONE)         return STR_LL_NONE;
    if (llType==LL_NORM2)        return STR_LL_NORM2;
    if (llType==LL_NORMAL)       return STR_LL_NORMAL;
    std::cout<<"Likelihood type integer '"<<llType<<"' not recognized."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return type;
}

/**
 * Translate from adstring likelihood type to int version.
 * 
 * @param llType - adstring likelihood type
 * @return - corresponding integer version
 */
int tcsam::getLikelihoodType(adstring llType){
    int type = 0;
    if (llType==STR_LL_BINOMIAL)     return LL_BINOMIAL;
    if (llType==STR_LL_GAMMA)        return LL_GAMMA;
    if (llType==STR_LL_LOGNORMAL)    return LL_LOGNORMAL;
    if (llType==STR_LL_MULTINOMIAL)  return LL_MULTINOMIAL;
    if (llType==STR_LL_MULTINOMIALP) return LL_MULTINOMIALP;
    if (llType==STR_LL_DIRICHLET)    return LL_DIRICHLET;
    if (llType==STR_LL_DIRICHLETP)   return LL_DIRICHLETP;
    if (llType==STR_LL_NONE)         return LL_NONE;
    if (llType==STR_LL_NORM2)        return LL_NORM2;
    if (llType==STR_LL_NORMAL)       return LL_NORMAL;
    std::cout<<"Likelihood type keyword '"<<llType<<"' not recognized."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    exit(-1);
    return type;
}

/**
 * Translate from adstring scale type to int version.
 * 
 * @param sclType - adstring scale type
 * @return -integer transform type
 */
int tcsam::getScaleType(adstring sclType){
    int type = -1;
    sclType.to_upper();
    if (sclType==STR_SCALE_NONE)   return SCALE_ARITHM;
    if (sclType==STR_SCALE_ARITHM) return SCALE_ARITHM;
    if (sclType==STR_SCALE_LOGIT)  return SCALE_LOGIT;
    if (sclType==STR_SCALE_LOG)    return SCALE_LOG;
    if (sclType==STR_SCALE_PROBIT) return SCALE_PROBIT;
    std::cout<<"transform type keyword '"<<sclType<<"' not recognized."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    ad_exit(-1);
    return type;
}

/**
 * Translate from integer scale type to adstring version.
 * 
 * @param sclType - integer scale type
 * @return  - corresponding adstring version
 */
adstring tcsam::getScaleType(int sclType){
    adstring type = "<UNDEFINED>";
    if (sclType==SCALE_ARITHM) return STR_SCALE_ARITHM;
    if (sclType==SCALE_LOGIT)  return STR_SCALE_LOGIT;
    if (sclType==SCALE_LOG)    return STR_SCALE_LOG;
    if (sclType==SCALE_PROBIT) return STR_SCALE_PROBIT;
    std::cout<<"transform type integer '"<<sclType<<"' not recognized."<<std::endl;
    std::cout<<"Aborting..."<<std::endl;
    ad_exit(-1);
    return type;
}
