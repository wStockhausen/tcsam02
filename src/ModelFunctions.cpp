#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelFunctions.hpp"

/**
 * Calculates a jittered value.
 * 
 * @param i - default initial value
 * @param l - lower bound
 * @param u - upper bound
 * @param jitFrac - jittering fraction (0-1)
 * @param r - uniform random deviate
 * 
 * @return jittered value
 */
double tcsam::jitterIt(double i, double l, double u, double jitFrac, double r){
    double d = u - l;
    l = l+0.001*d;//shrink lower bound
    u = u-0.001*d;//shrink upper bound
    d = u - l;    //shrink interval
    double lp = i - 0.5*d*jitFrac;
    double up = i + 0.5*d*jitFrac;
    double rp = i + (r-0.5)*d*jitFrac;
    if (rp>u)      {rp = lp - (rp-u);}
    else if (rp<l) {rp = up + (l-rp);}
    return rp;
}
/**
 * Prints a file read error.
 * 
 * @param is - input filestream
 * @param expP - the expected read value
 * @param gotP - the obtained read value
 */
void tcsam::readError(cifstream & is, const char * expP, adstring gotP){
    cout<<"Reading "<<is.get_file_name()<<endl;
    cout<<"Expected parameter name '"<<expP<<"' but got '"<<gotP<<"'."<<endl;
    cout<<"Aborting..."<<endl;
}
//----------------------------------------------------------------------
//          Model Functions
//----------------------------------------------------------------------
/***********************************************************************
modfcn_constant
    f(x) = c
    params(1) = c (constant value)
***********************************************************************/
dvar_vector tcsam::modfcn_constant(dvector& x, dvar_vector params, dvector& consts){
    RETURN_ARRAYS_INCREMENT();
    dvar_vector res = params(1)+0.0*x;//has same size as x
    RETURN_ARRAYS_DECREMENT();
    return res;
}

/***********************************************************************
modfcn_exp
    f(x) = exp(alpha+beta*x)
    params(1) = alpha
    params(2) = beta
***********************************************************************/
dvar_vector tcsam::modfcn_exp(dvector& x, dvar_vector params, dvector& consts){
    RETURN_ARRAYS_INCREMENT();
    dvar_vector res = mfexp(params(1)+params(2)*x);//has same size as x
    RETURN_ARRAYS_DECREMENT();
    return res;
}

/***********************************************************************
modfcn_logistic
    f(x) = 1/{1+exp(-beta*[x-x50])}
    params(1) = x50
    params(2) = beta (slope)
***********************************************************************/
dvar_vector tcsam::modfcn_logistic(dvector& x, dvar_vector params, dvector& consts){
    RETURN_ARRAYS_INCREMENT();
    dvar_vector res = 1.0/(1.0+mfexp(-params(2)*(x-params(1))));//has same size as x
    RETURN_ARRAYS_DECREMENT();
    return res;
}

/***********************************************************************
modfcn_logistic95
    f(x) = 1/{1+exp[-2*ln(19)*(x-x05)/(x95-x05)+ln(19)]}
    params(1) = x05
    params(2) = offset from x05 to x95 (x95-x05)
***********************************************************************/
dvar_vector tcsam::modfcn_logistic95(dvector& x, dvar_vector params, dvector& consts){
    RETURN_ARRAYS_INCREMENT();
    double s1 = log(1.0/0.05 - 1.0);
    double s2 = log(1.0/0.95 - 1.0);
    dvar_vector res = 1.0/(1.0+mfexp((s2-s1)*(x-params(1))/params(2)+s1));//normalized selfcn
    RETURN_ARRAYS_DECREMENT();
    return res;
}

/***********************************************************************
modfcn_none
    note that this function has no parameters
    and returns a 0-value vector
***********************************************************************/
dvar_vector tcsam::modfcn_none(dvector& x, dvar_vector params, dvector& consts){
    RETURN_ARRAYS_INCREMENT();
    dvar_vector res = 0.0*x;//has same size as x
    RETURN_ARRAYS_DECREMENT();
    return res;
}


//----------------------------------------------------------------------
//          ModelFunctionInfo
//----------------------------------------------------------------------
const adstring ModelFunctionInfo::FCNNAME_NONE         = "none";
const adstring ModelFunctionInfo::FCNNAME_BETA         = "beta";
const adstring ModelFunctionInfo::FCNNAME_BH           = "bh";
const adstring ModelFunctionInfo::FCNNAME_CAUCHY       = "cauchy";
const adstring ModelFunctionInfo::FCNNAME_CHISQ        = "chisq";
const adstring ModelFunctionInfo::FCNNAME_CONSTANT     = "constant";
const adstring ModelFunctionInfo::FCNNAME_DELTAFCN     = "deltafcn";
const adstring ModelFunctionInfo::FCNNAME_DOUBLENORMAL = "doublenormal";
const adstring ModelFunctionInfo::FCNNAME_EXP          = "exp";
const adstring ModelFunctionInfo::FCNNAME_GAMMA        = "gamma";
const adstring ModelFunctionInfo::FCNNAME_LOG          = "log";
const adstring ModelFunctionInfo::FCNNAME_LOGISTIC     = "logistic";
const adstring ModelFunctionInfo::FCNNAME_LOGISTIC95   = "logistic95";
const adstring ModelFunctionInfo::FCNNAME_LOGISTICEVP  = "logisticEVP";
const adstring ModelFunctionInfo::FCNNAME_LOGNORMAL    = "lognormal";
const adstring ModelFunctionInfo::FCNNAME_LORENZEN     = "lorenzen";
const adstring ModelFunctionInfo::FCNNAME_MULTINOM     = "multinom";
const adstring ModelFunctionInfo::FCNNAME_NORMAL       = "normal";
const adstring ModelFunctionInfo::FCNNAME_POISSON      = "poisson";
const adstring ModelFunctionInfo::FCNNAME_RICKER       = "ricker";
const adstring ModelFunctionInfo::FCNNAME_SLOT         = "slot";
const adstring ModelFunctionInfo::FCNNAME_T            = "t";
const adstring ModelFunctionInfo::FCNNAME_VBCEVP       = "VBcevp";
const adstring ModelFunctionInfo::FCNNAME_WEIBULL      = "weibull";

/***************************************************************
*   instance creation                                          *
***************************************************************/
ModelFunctionInfo* ModelFunctionInfo::getInfo(adstring name){
    ModelFunctionInfo* pMFI = 0;
    if (name==FCNNAME_NONE) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 0;
        pMFI->nmsVar.allocate(1,1);
        pMFI->nmsVar(1) = "none";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = tcsam::modfcn_none;
    } else
    if (name==FCNNAME_BH) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 2;
        pMFI->nmsVar.allocate(1,2);
        pMFI->nmsVar(1) = "alpha";
        pMFI->nmsVar(2) = "beta";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = 0;
    } else
    if (name==FCNNAME_CONSTANT) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nmsVar.allocate(1,1);
        pMFI->nmsVar(1) = "constant";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = tcsam::modfcn_constant;
    } else
    if (name==FCNNAME_EXP) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 2;
        pMFI->nmsVar.allocate(1,2);
        pMFI->nmsVar(1) = "alpha";
        pMFI->nmsVar(2) = "beta";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = tcsam::modfcn_exp;
    } else
    if (name==FCNNAME_LOGISTIC) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 2;
        pMFI->nmsVar.allocate(1,2);
        pMFI->nmsVar(1) = "X50";
        pMFI->nmsVar(2) = "slope";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = tcsam::modfcn_logistic;
    } else
    if (name==FCNNAME_LOGISTIC95) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 2;
        pMFI->nmsVar.allocate(1,2);
        pMFI->nmsVar(1) = "X05";
        pMFI->nmsVar(2) = "DX95";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = tcsam::modfcn_logistic95;
    } else
    if (name==FCNNAME_LOGNORMAL) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 2;
        pMFI->nmsVar.allocate(1,2);
        pMFI->nmsVar(1) = "med";
        pMFI->nmsVar(2) = "cv";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = 0;
    } else
    if (name==FCNNAME_LORENZEN) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 1;
        pMFI->nmsVar.allocate(1,1);
        pMFI->nmsVar(1) = "ref. rate";
        pMFI->nFxd = 1;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "ref. size";
        pMFI->pFcn = 0;
    } else
    if (name==FCNNAME_NORMAL) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 2;
        pMFI->nmsVar.allocate(1,2);
        pMFI->nmsVar(1) = "mean";
        pMFI->nmsVar(2) = "stdv";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = 0;
    } else
    if (name==FCNNAME_RICKER) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 2;
        pMFI->nmsVar.allocate(1,2);
        pMFI->nmsVar(1) = "alpha";
        pMFI->nmsVar(2) = "beta";
        pMFI->nFxd = 0;
        pMFI->nmsFxd.allocate(1,1);
        pMFI->nmsFxd(1) = "none";
        pMFI->pFcn = 0;
    } else
    if (name==FCNNAME_VBCEVP) {
        pMFI = new ModelFunctionInfo();
        pMFI->fcnName = name;
        pMFI->nVar = 6;
        pMFI->nmsVar.allocate(1,pMFI->nVar);
        pMFI->nmsVar(1) = "L1";
        pMFI->nmsVar(2) = "dL1M";
        pMFI->nmsVar(3) = "dLM2";
        pMFI->nmsVar(4) = "stdDev1";
        pMFI->nmsVar(5) = "stdDev2";
        pMFI->nmsVar(6) = "stdDev3";
        pMFI->nFxd = 4;
        pMFI->nmsFxd.allocate(1,pMFI->nFxd);
        pMFI->nmsFxd(1) = "age1";
        pMFI->nmsFxd(2) = "age2";
        pMFI->nmsFxd(3) = "sdmType";
        pMFI->nmsFxd(4) = "errType";
        pMFI->pFcn = 0;
    } else
    {
        cout<<"Model function name '"<<name<<"' was not recognized."<<endl;
        int test;
        cout<<"Enter 1 to continue: ";
        cin>>test;
        if (test<0) exit(-1);
    }
    return pMFI;
}

/***************************************************************
*   Get names for function constants.                          *
***************************************************************/
adstring ModelFunctionInfo::getStringForConstsNames(){
    adstring nms = "";
    for (int p=1;p<=nFxd;p++) nms += nmsFxd(p)+"  ";
    return nms;
}

/***************************************************************
*   Get names for function parameters.                         *
***************************************************************/
adstring ModelFunctionInfo::getStringForParamsNames(){
    adstring nms = "";
    for (int p=1;p<=nVar;p++) nms += nmsVar(p)+"  ";
    return nms;
}

/***************************************************************
*   Calculate the function.                                    *
***************************************************************/
dvar_vector ModelFunctionInfo::calc(dvector& x, dvar_vector params, dvector& consts){
    if (pFcn) return (*pFcn)(x,params,consts);
    RETURN_ARRAYS_INCREMENT();
    dvar_vector res;
    RETURN_ARRAYS_DECREMENT();
    return res;
}

//----------------------------------------------------------------------
//          ModelPDFInfo
//----------------------------------------------------------------------
const adstring ModelPDFInfo::PDFTYPE_NONE         = "none";
const adstring ModelPDFInfo::PDFTYPE_AR1_NORMAL   ="ar1_normal";
//const adstring ModelPDFInfo::PDFTYPE_BETA         = "beta";
const adstring ModelPDFInfo::PDFTYPE_CAUCHY       = "cauchy";
const adstring ModelPDFInfo::PDFTYPE_CHISQ        = "chisquare";
const adstring ModelPDFInfo::PDFTYPE_CONSTANT     = "constant";
//const adstring ModelPDFInfo::PDFTYPE_DELTAFCN     = "deltafcn";
const adstring ModelPDFInfo::PDFTYPE_EXPNORMAL    = "expnormal";
const adstring ModelPDFInfo::PDFTYPE_EXPONENTIAL  = "exponential";
const adstring ModelPDFInfo::PDFTYPE_GAMMA        = "gamma";
const adstring ModelPDFInfo::PDFTYPE_INVCHISQ     = "invchisquare";
const adstring ModelPDFInfo::PDFTYPE_INVGAMMA     = "invgamma";
const adstring ModelPDFInfo::PDFTYPE_INVGAUSSIAN  = "invgaussian";
//const adstring ModelPDFInfo::PDFTYPE_LOGISTIC     = "logistic";
const adstring ModelPDFInfo::PDFTYPE_LOGNORMAL    = "lognormal";
const adstring ModelPDFInfo::PDFTYPE_LOGSCALE_NORMAL   = "logscale_normal";
//const adstring ModelPDFInfo::PDFTYPE_MULTINOMIAL       = "multinomial";
const adstring ModelPDFInfo::PDFTYPE_NORMAL            = "normal";
//const adstring ModelPDFInfo::PDFTYPE_POISSON           = "poisson";
const adstring ModelPDFInfo::PDFTYPE_SCALED_INVCHISQ   = "scaled_invchisquare";
const adstring ModelPDFInfo::PDFTYPE_SCALEDCV_INVCHISQ = "scaledCV_invchisquare";
const adstring ModelPDFInfo::PDFTYPE_T                 = "t";
const adstring ModelPDFInfo::PDFTYPE_TRUNCATED_NORMAL  ="truncated_normal";
//const adstring ModelPDFInfo::PDFTYPE_WEIBULL      = "weibull";

/***************************************************************
*   instance creation                                          *
***************************************************************/
ModelPDFInfo* ModelPDFInfo::getInfo(adstring type){
    ModelPDFInfo* pMPI = 0;
    if (type==PDFTYPE_NONE) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_NONE;
        pMPI->nVar = 0;
        pMPI->nmsVar.allocate(1,1);
        pMPI->nmsVar(1) = "none";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_constant;
        pMPI->pSmplr = 0;
    } else
    if (type==PDFTYPE_CONSTANT) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_CONSTANT;
        pMPI->nVar = 1;
        pMPI->nmsVar.allocate(1,1);
        pMPI->nmsVar(1) = "constant";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_constant;
        pMPI->pSmplr = 0;
    } else
    if (type==PDFTYPE_AR1_NORMAL) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_AR1_NORMAL;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "mean";
        pMPI->nmsVar(2) = "stdv";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->vpPDF = wts::logPDF_AR1_normal;
        pMPI->vpSmplr = wts::samplePDF_AR1_normal;
    } else
//     if (type==PDFTYPE_BETA) {
//         pMPI = new ModelPDFInfo();
//         pMPI->pdfType = PDFTYPE_BETA;
//         pMPI->nVar = 2;
//         pMPI->nmsVar.allocate(1,2);
//         pMPI->nmsVar(1) = "alpha";
//         pMPI->nmsVar(2) = "beta";
//         pMPI->nFxd = 0;
//         pMPI->nmsFxd.allocate(1,1);
//         pMPI->nmsFxd(1) = "none";
//         pMPI->pPDF = wts::logPDF_beta;
//         pMPI->pSmplr = wts::samplePDF_beta;
//     } else
    if (type==PDFTYPE_CAUCHY) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_CAUCHY;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "x0";
        pMPI->nmsVar(2) = "gamma";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_cauchy;
        pMPI->pSmplr = wts::samplePDF_cauchy;
    } else
    if (type==PDFTYPE_CHISQ) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_CHISQ;
        pMPI->nVar = 1;
        pMPI->nmsVar.allocate(1,1);
        pMPI->nmsVar(1) = "nu";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_chisquare;
        pMPI->pSmplr = wts::samplePDF_chisquare;
    } else
    if (type==PDFTYPE_EXPNORMAL) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_EXPNORMAL;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "expmean";
        pMPI->nmsVar(2) = "stdv";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_expnormal;
        pMPI->pSmplr = wts::samplePDF_expnormal;
    } else
    if (type==PDFTYPE_EXPONENTIAL) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_EXPONENTIAL;
        pMPI->nVar = 1;
        pMPI->nmsVar.allocate(1,1);
        pMPI->nmsVar(1) = "lambda";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_exponential;
        pMPI->pSmplr = wts::samplePDF_exponential;
    } else
    if (type==PDFTYPE_GAMMA) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_GAMMA;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "r";
        pMPI->nmsVar(2) = "mu";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_gamma;
        pMPI->pSmplr = wts::samplePDF_gamma;
    } else
    if (type==PDFTYPE_INVCHISQ) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_INVCHISQ;
        pMPI->nVar = 1;
        pMPI->nmsVar.allocate(1,1);
        pMPI->nmsVar(1) = "nu";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_invchisquare;
        pMPI->pSmplr = wts::samplePDF_invchisquare;
    } else
    if (type==PDFTYPE_INVGAMMA) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_INVGAMMA;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "r";
        pMPI->nmsVar(2) = "mu";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_invgamma;
        pMPI->pSmplr = wts::samplePDF_invgamma;
    } else
    if (type==PDFTYPE_INVGAUSSIAN) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_INVGAUSSIAN;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "mu";
        pMPI->nmsVar(2) = "lambda";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_invgaussian;
        pMPI->pSmplr = wts::samplePDF_invgaussian;
    } else
    if (type==PDFTYPE_LOGNORMAL) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_LOGNORMAL;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "med";
        pMPI->nmsVar(2) = "cv";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_lognormal;
        pMPI->pSmplr = wts::samplePDF_lognormal;
    } else
    if (type==PDFTYPE_LOGSCALE_NORMAL) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_LOGSCALE_NORMAL;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "med";
        pMPI->nmsVar(2) = "cv";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_logscale_normal;
        pMPI->pSmplr = wts::samplePDF_lognormal;//just use lognormal sampler
    } else
    if (type==PDFTYPE_NORMAL) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_NORMAL;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "mean";
        pMPI->nmsVar(2) = "stdv";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_normal;
        pMPI->pSmplr = wts::samplePDF_normal;
    } else
    if (type==PDFTYPE_SCALED_INVCHISQ) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_SCALED_INVCHISQ;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "nu";
        pMPI->nmsVar(2) = "s";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_scaled_invchisquare;
        pMPI->pSmplr = wts::samplePDF_scaled_invchisquare;
    } else
    if (type==PDFTYPE_SCALEDCV_INVCHISQ) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_SCALEDCV_INVCHISQ;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "nu";
        pMPI->nmsVar(2) = "cv";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_scaledCV_invchisquare;
        pMPI->pSmplr = wts::samplePDF_scaledCV_invchisquare;
    } else
    if (type==PDFTYPE_T) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_T;
        pMPI->nVar = 1;
        pMPI->nmsVar.allocate(1,1);
        pMPI->nmsVar(1) = "nu";
        pMPI->nFxd = 0;
        pMPI->nmsFxd.allocate(1,1);
        pMPI->nmsFxd(1) = "none";
        pMPI->pPDF = wts::logPDF_t;
        pMPI->pSmplr = wts::samplePDF_t;
    } else
    if (type==PDFTYPE_TRUNCATED_NORMAL) {
        pMPI = new ModelPDFInfo();
        pMPI->pdfType = PDFTYPE_TRUNCATED_NORMAL;
        pMPI->nVar = 2;
        pMPI->nmsVar.allocate(1,2);
        pMPI->nmsVar(1) = "mean";
        pMPI->nmsVar(2) = "stdv";
        pMPI->nFxd = 2;
        pMPI->nmsFxd.allocate(1,2);
        pMPI->nmsFxd(1) = "min";
        pMPI->nmsFxd(2) = "max";
        pMPI->pPDF = wts::logPDF_truncated_normal;
        pMPI->pSmplr = wts::samplePDF_truncated_normal;
    } else
    {
        cout<<"Model pdf type '"<<type<<"' was not recognized."<<endl;
        int test;
        cout<<"Enter 1 to continue: ";
        cin>>test;
        if (test<0) exit(-1);
    }
    return pMPI;
}

/***************************************************************
*   get function constants                                     *
***************************************************************/
adstring ModelPDFInfo::getStringForConstsNames(){
    adstring nms = "";
    for (int p=1;p<=nFxd;p++) nms += nmsFxd(p)+"  ";
    return nms;
}

/***************************************************************
*   get function parameters                                    *
***************************************************************/
adstring ModelPDFInfo::getStringForParamsNames(){
    adstring nms = "";
    for (int p=1;p<=nVar;p++) nms += nmsVar(p)+"  ";
    return nms;

}

/***************************************************************
*   calculate the log(pdf(x))                                  *
***************************************************************/
dvariable ModelPDFInfo::calcLogPDF(prevariable& val,dvar_vector params,dvector& consts){
    RETURN_ARRAYS_INCREMENT();
    dvariable res = (*pPDF)(val,params,consts);
    RETURN_ARRAYS_DECREMENT();
    return res;
}

dvar_vector ModelPDFInfo::calcLogPDF(dvar_vector& vals,dvar_vector params,dvector& consts){
    RETURN_ARRAYS_INCREMENT();
    dvar_vector res = (*vpPDF)(vals,params,consts);
    RETURN_ARRAYS_DECREMENT();
    return res;
}

/***************************************************************
*   get a sample from the pdf(x)                               *
***************************************************************/
double ModelPDFInfo::drawSample(random_number_generator& rng, dvector& params, dvector& consts){
    double val = 0.0;
    if (pSmplr) val = (*pSmplr)(rng,params,consts);
    return val;
}

//----------------------------------------------------------------------
//          ModelTransformInfo
//----------------------------------------------------------------------
const adstring ModelTransformInfo::TYPE_NONE         = "none";
const adstring ModelTransformInfo::TYPE_COS          = "cos";
const adstring ModelTransformInfo::TYPE_EXP          = "exp";
const adstring ModelTransformInfo::TYPE_EXPNEG       = "expneg";
const adstring ModelTransformInfo::TYPE_LOG          = "log";
const adstring ModelTransformInfo::TYPE_LOGNEG       = "logneg";
const adstring ModelTransformInfo::TYPE_LOGISTIC     = "logistic";
const adstring ModelTransformInfo::TYPE_LOGIT        = "logit";
const adstring ModelTransformInfo::TYPE_SIN          = "sin";
const adstring ModelTransformInfo::TYPE_SQRT         = "sqrt";
const adstring ModelTransformInfo::TYPE_SQR          = "sqr";
const adstring ModelTransformInfo::TYPE_TAN          = "tan";

/***************************************************************
*   instance creation                                          *
***************************************************************/
ModelTransformInfo* ModelTransformInfo::getInfo(adstring type){
    ModelTransformInfo* pMTI = 0;
    if (type==TYPE_NONE) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::none;
        pMTI->pInv     = wts::none;
        pMTI->pFcnDvar = wts::none;
        pMTI->pInvDvar = wts::none;
    } else
    if (type==TYPE_COS) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::cos;
        pMTI->pInv     = wts::acos;
        pMTI->pFcnDvar = wts::cos;
        pMTI->pInvDvar = wts::acos;
    } else
    if (type==TYPE_EXP) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::exp;
        pMTI->pInv     = wts::log;
        pMTI->pFcnDvar = wts::exp;
        pMTI->pInvDvar = wts::log;
    } else
    if (type==TYPE_EXPNEG) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::expneg;
        pMTI->pInv     = wts::logneg;
        pMTI->pFcnDvar = wts::expneg;
        pMTI->pInvDvar = wts::logneg;
    } else
    if (type==TYPE_LOG) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::log;
        pMTI->pInv     = wts::exp;
        pMTI->pFcnDvar = wts::log;
        pMTI->pInvDvar = wts::exp;
    } else
    if (type==TYPE_LOGNEG) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::logneg;
        pMTI->pInv     = wts::expneg;
        pMTI->pFcnDvar = wts::logneg;
        pMTI->pInvDvar = wts::expneg;
    } else
    if (type==TYPE_LOGISTIC) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 2;
        pMTI->consts.allocate(1,2);
        pMTI->nmsFxd.allocate(1,2);
        pMTI->nmsFxd(1)="min";
        pMTI->nmsFxd(2)="max";
        pMTI->pFcn     = wts::logistic;
        pMTI->pInv     = wts::logit;
        pMTI->pFcnDvar = wts::logistic;
        pMTI->pInvDvar = wts::logit;
    } else
    if (type==TYPE_LOGIT) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 2;
        pMTI->consts.allocate(1,2);
        pMTI->nmsFxd.allocate(1,2);
        pMTI->nmsFxd(1)="min";
        pMTI->nmsFxd(2)="max";
        pMTI->pFcn     = wts::logit;
        pMTI->pInv     = wts::logistic;
        pMTI->pFcnDvar = wts::logit;
        pMTI->pInvDvar = wts::logistic;
    } else
    if (type==TYPE_SIN) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::sin;
        pMTI->pInv     = wts::asin;
        pMTI->pFcnDvar = wts::sin;
        pMTI->pInvDvar = wts::asin;
    } else
    if (type==TYPE_SQR) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::square;
        pMTI->pInv     = wts::sqrt;
        pMTI->pFcnDvar = wts::square;
        pMTI->pInvDvar = wts::sqrt;
    } else
    if (type==TYPE_SQRT) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::sqrt;
        pMTI->pInv     = wts::square;
        pMTI->pFcnDvar = wts::sqrt;
        pMTI->pInvDvar = wts::square;
    } else
    if (type==TYPE_TAN) {
        pMTI = new ModelTransformInfo();
        pMTI->fcnType = type;
        pMTI->nC = 0;
        pMTI->consts.allocate(1,1);
        pMTI->pFcn     = wts::tan;
        pMTI->pInv     = wts::atan;
        pMTI->pFcnDvar = wts::tan;
        pMTI->pInvDvar = wts::atan;
    } else
    {
        cout<<"Model transform type '"<<type<<"' was not recognized."<<endl;
        int test;
        cout<<"Enter 1 to continue: ";
        cin>>test;
        if (test<0) exit(-1);
    }
    return pMTI;
}

/***************************************************************
*   get function constants                                     *
***************************************************************/
adstring ModelTransformInfo::getStringForConstsNames(){
    adstring nms = "";
    for (int p=1;p<=nC;p++) nms += nmsFxd(p)+"  ";
    return nms;
}

