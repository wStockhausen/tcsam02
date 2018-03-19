//--------------------------------------------------------------------------------
//  Includes:
//      ModelConsts
//      ModelFunctionInfo
//      ModelPDFInfo
//      ModelTransformInfo
//--------------------------------------------------------------------------------
#pragma once
#ifndef MODELFUNCTIONS_HPP
    #define MODELFUNCTIONS_HPP

#include <admodel.h>

//----------------------------------------------------------------------
//          Model Functions
//----------------------------------------------------------------------
typedef dvar_vector (*ModFcnPtr)          (dvector&,                dvar_vector,       dvector&);
typedef double      (*TransformFcnPtr)    (double,                  const dvector&);
typedef dvariable   (*dvarTransformFcnPtr)(const prevariable&,      const dvector&);
typedef dvariable   (*PdfFcnPtr)          (const prevariable&,      const dvar_vector&,const dvector&);
typedef double      (*PdfSamplerPtr)      (random_number_generator&,const dvector&,    const dvector&);
typedef dvar_vector (*vPdfFcnPtr)         (const dvar_vector&,      const dvar_vector&,const dvector&);
typedef dvector     (*vPdfSamplerPtr)     (int n, random_number_generator&,const dvector&,    const dvector&);

namespace tcsam{
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
    double jitterIt(double i, double l, double u, double jitFrac, double r);
    /**
     * Prints a file read error.
     * 
     * @param is - input filestream
     * @param expP - the expected read value
     * @param gotP - the obtained read value
     */
    void readError(cifstream & is, const char * expP, adstring gotP);
    
    /***********************************************************************
    modfcn_constant
        f(x) = c
        params(1) = c (constant value)
    ***********************************************************************/
    dvar_vector modfcn_constant(dvector& x, dvar_vector params, dvector& consts);
    /***********************************************************************
    modfcn_exp
        f(x) = exp(alpha+beta*x)
        params(1) = alpha
        params(2) = beta
    ***********************************************************************/
    dvar_vector modfcn_exp(dvector& x, dvar_vector params, dvector& consts);
    /***********************************************************************
    modfcn_logistic
        f(x) = 1/{1+exp(-beta*[x-x50])}
        params(1) = x50
        params(2) = beta (slope)
    ***********************************************************************/
    dvar_vector modfcn_logistic(dvector& x, dvar_vector params, dvector& consts);
    /***********************************************************************
    modfcn_logistic95
        f(x) = 1/{1+exp[-2*ln(19)*(x-x05)/(x95-x05)+ln(19)]}
        params(1) = x05
        params(2) = offset from x05 to x95 (x95-x05)
    ***********************************************************************/
    dvar_vector modfcn_logistic95(dvector& x, dvar_vector params, dvector& consts);
    /***********************************************************************
    modfcn_none
        note that this function has no parameters
        and returns a 0-value vector
    ***********************************************************************/
    dvar_vector modfcn_none(dvector& x, dvar_vector params, dvector& consts);
}  //namespace tcsam
//----------------------------------------------------------------------
//          ModelFunctionInfo
//----------------------------------------------------------------------
    class ModelFunctionInfo {
    public:
        const static adstring FCNNAME_NONE;
        const static adstring FCNNAME_BETA;
        const static adstring FCNNAME_BH;
        const static adstring FCNNAME_CAUCHY;
        const static adstring FCNNAME_CHISQ;
        const static adstring FCNNAME_CONSTANT;
        const static adstring FCNNAME_DELTAFCN;
        const static adstring FCNNAME_DOUBLENORMAL;
        const static adstring FCNNAME_EXP;
        const static adstring FCNNAME_GAMMA;
        const static adstring FCNNAME_LOG;
        const static adstring FCNNAME_LOGISTIC;
        const static adstring FCNNAME_LOGISTIC95;
        const static adstring FCNNAME_LOGISTICEVP;
        const static adstring FCNNAME_LOGNORMAL;
        const static adstring FCNNAME_LORENZEN;
        const static adstring FCNNAME_MULTINOM;
        const static adstring FCNNAME_NORMAL;
        const static adstring FCNNAME_POISSON;
        const static adstring FCNNAME_RICKER;
        const static adstring FCNNAME_SLOT;
        const static adstring FCNNAME_T;
        const static adstring FCNNAME_VBCEVP;
        const static adstring FCNNAME_WEIBULL;

        public:
            static ModelFunctionInfo* getInfo(adstring name);
        protected:
            adstring fcnName;
            int nVar;//number of variable parameters
            adstring_array nmsVar;//array of names for variable parameters
            int nFxd;//number of fixed parameters
            adstring_array nmsFxd;//array of names for fixed parameters
            ModFcnPtr pFcn;//model function
        protected:
            ModelFunctionInfo(){pFcn=0;}
        public:
            ~ModelFunctionInfo(){pFcn=0;}
            adstring getFunctionName(void){return fcnName;}
            int getNumConsts(void){return nFxd;}
            int getNumParams(void){return nVar;}
            adstring_array getNamesForConsts(){return nmsFxd;}
            adstring_array getNamesForParams(){return nmsVar;}
            adstring getStringForConstsNames();
            adstring getStringForParamsNames();
            dvar_vector calc(dvector& x, dvar_vector params, dvector& consts);//calc model function
    };

//----------------------------------------------------------------------
//          ModelPDFInfo
//----------------------------------------------------------------------
    class ModelPDFInfo {
    public:
        const static adstring PDFTYPE_NONE;
        const static adstring PDFTYPE_AR1_NORMAL;
//        const static adstring PDFTYPE_BETA;
        const static adstring PDFTYPE_CAUCHY;
        const static adstring PDFTYPE_CHISQ;
        const static adstring PDFTYPE_CONSTANT;
        const static adstring PDFTYPE_EXPNORMAL;
        const static adstring PDFTYPE_EXPONENTIAL;
        const static adstring PDFTYPE_GAMMA;
        const static adstring PDFTYPE_INVCHISQ;
        const static adstring PDFTYPE_INVGAMMA;
        const static adstring PDFTYPE_INVGAUSSIAN;
        const static adstring PDFTYPE_LOGNORMAL;
        const static adstring PDFTYPE_LOGSCALE_NORMAL;
        const static adstring PDFTYPE_NORMAL;
        const static adstring PDFTYPE_SCALED_INVCHISQ;
        const static adstring PDFTYPE_SCALEDCV_INVCHISQ;
        const static adstring PDFTYPE_T;
        const static adstring PDFTYPE_TRUNCATED_NORMAL;
        

        public:
            static ModelPDFInfo* getInfo(adstring pdfType);
        protected:
            adstring pdfType;     //pdf type as string
            int nVar;             //number of variable parameters
            adstring_array nmsVar;//array of names for variable parameters
            int nFxd;             //number of fixed parameters
            adstring_array nmsFxd;//array of names for fixed parameters
            PdfFcnPtr pPDF;       //pointer to log pdf function
            PdfSamplerPtr pSmplr; //ptr to function to sample pdf
            vPdfFcnPtr vpPDF;       //pointer to log vector pdf function
            vPdfSamplerPtr vpSmplr; //ptr to function to sample vector pdf
        protected:
            ModelPDFInfo(){pPDF=0;pSmplr=0;vpPDF=0;vpSmplr=0;}
        public:
            ~ModelPDFInfo(){pPDF=0;pSmplr=0;vpPDF=0;vpSmplr=0;}
            int canCalcLogPDF(prevariable& val){if (pPDF)  return 1; else return 0;}
            int canCalcLogPDF(dvar_vector& val){if (vpPDF) return 1; else return 0;}
            dvariable calcLogPDF(prevariable& val,dvar_vector params,dvector& consts);
            dvar_vector calcLogPDF(dvar_vector& val,dvar_vector params,dvector& consts);
            int canSample(void){if ((pSmplr)||(vpSmplr)) return 1; else return 0;}
            double drawSample(random_number_generator& rng,dvector& params, dvector& consts);
            dvector drawSample(int n,random_number_generator& rng,dvector& params, dvector& consts);
            adstring getPDFType(void){return pdfType;}
            int getNumConsts(void){return nFxd;}
            int getNumParams(void){return nVar;}
            adstring_array getNamesForConsts(){return nmsFxd;}
            adstring_array getNamesForParams(){return nmsVar;}
            adstring getStringForConstsNames();
            adstring getStringForParamsNames();
    };

//----------------------------------------------------------------------
//          ModelTransformInfo
//----------------------------------------------------------------------
    class ModelTransformInfo {
    public:
        const static adstring TYPE_NONE;
        const static adstring TYPE_EXP;
        const static adstring TYPE_EXPNEG;
        const static adstring TYPE_COS;
        const static adstring TYPE_LOG;
        const static adstring TYPE_LOGNEG;
        const static adstring TYPE_LOGIT;
        const static adstring TYPE_LOGISTIC;
        const static adstring TYPE_SQRT;
        const static adstring TYPE_SQR;
        const static adstring TYPE_SIN;
        const static adstring TYPE_TAN;

        public:
            static ModelTransformInfo* getInfo(adstring type);
            dvector consts;        //values of transform constants
        protected:
            adstring fcnType;      //transform type
            int nC;                //number of transform constants
            adstring_array nmsFxd; //names for consts
            TransformFcnPtr pFcn;
            TransformFcnPtr pInv;
            dvarTransformFcnPtr pFcnDvar;
            dvarTransformFcnPtr pInvDvar;
        protected:
            ModelTransformInfo(){pFcn=0;pInv=0;pFcnDvar=0;pInvDvar=0;}
        public:
            ~ModelTransformInfo(){pFcn=0;pInv=0;pFcnDvar=0;pInvDvar=0;}
            adstring getTransformType(void){return fcnType;}
            int getNumConsts(void){return nC;}
            adstring_array getNamesForConsts(){return nmsFxd;}
            adstring getStringForConstsNames();
            double calcTransform(double val){return (*pFcn)(val,consts);}
            double calcInvTransform(double val){return (*pInv)(val,consts);}
            dvariable calcTransform(const prevariable& val){RETURN_ARRAYS_INCREMENT();dvariable res = (*pFcnDvar)(val,consts);RETURN_ARRAYS_DECREMENT();return res;}
            dvariable calcInvTransform(const prevariable& val){RETURN_ARRAYS_INCREMENT();dvariable res = (*pInvDvar)(val,consts);RETURN_ARRAYS_DECREMENT();return res;}
    };

#endif