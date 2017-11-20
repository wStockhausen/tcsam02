/*------------------------------------------------------------------------------
 *  Includes:
 *      NumberInfo
 *      BoundedNumberInfo
 *      VectorInfo
 *      BoundedVectorInfo
 *      DevsVectorInfo
 *      NumberVectorInfo
 *      BoundedNumberVectorInfo
 *      VectorVectorInfo
 *      BoundedVectorVectorInfo
 *      DevsVectorVectorInfo
 *----------------------------------------------------------------------------*/
#include <admodel.h>
#include <typeinfo>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelFunctions.hpp"
#include "ModelConfiguration.hpp"
#include "ModelParameterInfoTypes.hpp"

int NumberInfo::debug              = 0;
int BoundedNumberInfo::debug       = 0;
int VectorInfo::debug              = 0;
int BoundedVectorInfo::debug       = 0;
int DevsVectorInfo::debug          = 0;
int NumberVectorInfo::debug        = 0;
int BoundedNumberVectorInfo::debug = 0;
int VectorVectorInfo::debug        = 0;
int BoundedVectorVectorInfo::debug = 0;
int DevsVectorVectorInfo::debug    = 0;
/*------------------------------------------------------------------------------
 *  NumberInfo
 *----------------------------------------------------------------------------*/
/**
 * Calculate arithmetic-scale value corresponding to parameter-scale value.
 * 
 * @param x - parameter-scale value
 * @return - arithmetic-scale value
 */
double NumberInfo::calcArithScaleVal(double x){
    if (debug) rpt::echo<<"starting NumberInfo::calcArithScaleVal(double&) "<<this<<endl;
    double z;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            z = 1.0/(1.0+mfexp(-x));//scaled to [0,1]
            break;
        default:
            cout<<"NumberInfo::calcArithScaleVal(double) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"starting NumberInfo::calcArithScaleVal(double&) "<<this<<endl;
    return z;
}

/**
 * Calculates the arithmetic-scale value corresponding to given parameter-scale value.
 * 
 * @param x - parameter-scale value
 * @return - arithmetic-scale value
 */
dvariable NumberInfo::calcArithScaleVal(const prevariable& x){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting NumberInfo::calcArithScaleVal(prevariable&) "<<this<<endl;
    dvariable z;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);//scaled to [0,1]
            break;
        case tcsam::SCALE_LOGIT:
            z = 1.0/(1.0+mfexp(-x));//scaled to [0,1]
            break;
        default:
            cout<<"NumberInfo::calcArithScaleVal(dvariable&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished NumberInfo::calcArithScaleVal(prevariable&) "<<this<<endl;
    RETURN_ARRAYS_DECREMENT();
    return z;
}

/**
 * Calculate parameter-scale value corresponding to given arithmetic-scale value.
 * 
 * @param x - arithmetic-scale value
 * @return - parameter-scale value
 */
double NumberInfo::calcParamScaleVal(double x){
    if (debug) rpt::echo<<"starting NumberInfo::calcParamScaleVal(double&) "<<this<<endl;
    double y,z;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = log(x);//on log scale
            break;
        case tcsam::SCALE_LOGIT:
            z = log(x/(1.0-x));//on logit scale
            break;
        default:
            cout<<"NumberInfo::calcParamScaleVal(double) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized for NumberInfo!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedNumberInfo::calcParamScaleVal(double&) "<<this<<endl;
    return z;
}
            
/**
 * Draws a random value on the arithmetic scale based on the specified prior probability distribution.
 * 
 * Note that this DOES NOT update initVal.
 * If the parameter estimation phase is \< 0, the value for initVal is returned.
 * 
 * @param rng - a reference to a random number generator
 * @param vif - the variance inflation factor
 * 
 * @return - the random value, or value of initVal
 */
double NumberInfo::drawInitVal(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting NumberInfo::drawSample(random_number_generator& rng, vif)"<<this<<endl;
    double smpl  = initVal;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        smpl = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        rpt::echo<<phase<<tb<<pMPI->canSample()<<tb<<initVal<<tb<<smpl<<endl;
        rpt::echo<<"finished NumberInfo::drawSample(random_number_generator& rng,vif)"<<this<<endl;
    }
    return smpl;
}

/**
 * Calculates the log-scale prior probability based on an input value.
 * 
 * The input value is on the arithmetic scale, but is presumably from the associated parameter value.
 * 
 * @param x - the arithmetic-scale input value, presumably based on the associated parameter
 * 
 * @return - the log-scale prior, as a dvariable
 */
dvariable NumberInfo::calcLogPrior(prevariable& x){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting NumberInfo::calcLogPrior(prevariable& x)"<<this<<endl;
    dvariable val = pMPI->calcLogPDF(x,priorParams,priorConsts);
    if (debug) {
        if (pMPI->getNumParams()) rpt::echo<<"priorParams = "<<priorParams<<tb;
        if (pMPI->getNumConsts()) rpt::echo<<"priorConsts = "<<priorConsts<<tb;
        rpt::echo<<endl;
        rpt::echo<<"x = "<<x<<"; ln(pdf(x)) = "<<val<<endl;
        rpt::echo<<"finished NumberInfo::calcLogPrior(prevariable& x)"<<this<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return val;
}

/**
 * Sets the prior type, based on an adstring value.
 * 
 * @param prior - the prior type, as an adstring
 */
void NumberInfo::setPriorType(adstring & prior){
    if (debug) rpt::echo<<"starting NumberInfo::setPriorType(adstring prior)"<<this<<endl;
    if (pMPI) delete pMPI;
    pMPI = ModelPDFInfo::getInfo(prior);
    if (pMPI) {
        if (pMPI->getNumParams()>0) priorParams.allocate(1,pMPI->getNumParams());
        if (pMPI->getNumConsts()>0) priorConsts.allocate(1,pMPI->getNumConsts());
    }
    if (debug) rpt::echo<<"finished NumberInfo::setPriorType(adstring prior)"<<this<<endl;
}

/**
 * Reads the parameter info in ADMB format from an input stream.
 * 
 * The read order for the parameter info is:
 * <ul>
 *  <li> initVal
 *  <li> scaleType
 *  <li> phase
 *  <li> resample
 *  <li> priorWgt
 *  <li> priorType
 *  <li> priorParams
 *  <li> priorConsts
 *  <li> label
 * </ul>
 * 
 * @param is - the input stream
 * 
 */
void NumberInfo::read(cifstream & is){
    if (debug) rpt::echo<<"Starting NumberInfo::read(cifstream & is) "<<this<<endl;
    adstring str;
    is>>initVal;
    is>>str; scaleType=tcsam::getScaleType(str);
    std::cout<<name<<tb<<str<<tb<<scaleType<<tb<<tcsam::getScaleType(scaleType)<<endl;
    is>>phase;
    is>>str; resample=wts::getOnOffType(str);
    is>>priorWgt;
    is>>priorType;//prior type
    if (debug){
        rpt::echo<<initVal<<tb<<"#initVal"<<endl;
        rpt::echo<<phase  <<tb<<"#phase"<<endl;
        rpt::echo<<wts::getOnOffType(resample)<<tb<<"#resample"<<endl;
        rpt::echo<<priorWgt<<tb<<"#priorWgt"<<endl;
        rpt::echo<<priorType<<tb<<"#prior type"<<endl;
    }
    setPriorType(priorType);
    if (pMPI->getNumParams()) {
        is>>priorParams;
        if (debug) rpt::echo<<priorParams<<tb<<"#prior params"<<endl;
    } else {
        is>>str; str.to_upper();
        if (str!=tcsam::STR_NONE){
            rpt::echo<<"no prior params, so input should be 'none', but got '"<<str<<"'"<<endl;
            rpt::echo<<"Please fix!!"<<endl;
            exit(-1);
        }
        if (debug) rpt::echo<<"#no prior params"<<endl;
    }
    if (pMPI->getNumConsts()) {
        is>>priorConsts;
        if (debug) rpt::echo<<priorConsts<<tb<<"#prior consts"<<endl;
    } else {
        is>>str; str.to_upper();
        if (str!=tcsam::STR_NONE){
            rpt::echo<<"no prior consts, so input should be 'none', but got '"<<str<<"'"<<endl;
            rpt::echo<<"Please fix!!"<<endl;
            exit(-1);
        }
        if (debug) rpt::echo<<"#no prior consts"<<endl;
    }
    is>>label;
    if (debug) rpt::echo<<label<<tb<<"#label"<<endl;
    if (debug) rpt::echo<<"Done NumberInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
 * Write to ofstream object.\n
 * Write order is:\n
 *   initVal, phase, resample, priorWgt, \n
 *   priorType, priorParams, priorConsts\n
 ***************************************************************/
/**
 * Writes info to output stream.
 * 
 * Write order is in same sense as the read order:
 * <ul>
 *  <li> initVal
 *  <li> scaleType
 *  <li> phase
 *  <li> resample
 *  <li> priorWgt
 *  <li> priorType
 *  <li> priorParams
 *  <li> priorConsts
 *  <li> label
 * </ul>
 * 
 * @param os - the output stream to write to
 */
void NumberInfo::write(ostream & os){
    os<<initVal<<tb;
    os<<tcsam::getScaleType(scaleType)<<tb;
    os<<phase<<tb;
    os<<wts::getOnOffType(resample)<<tb;
    os<<priorWgt<<tb;
    os<<priorType<<tb;
    if (pMPI->getNumParams()) os<<priorParams<<tb<<tb; else os<<"none"<<tb<<tb;
    if (pMPI->getNumConsts()) os<<priorConsts<<tb<<tb; else os<<"none"<<tb<<tb;
    os<<label<<tb<<tb;
    os<<"#";
    if (pMPI->getNumParams()){
        os<<"prior params are ";
        adstring_array pNames = pMPI->getNamesForParams();
        for (int p=pNames.indexmin();p<=pNames.indexmax()-1;p++) os<<pNames(p)<<cc;
        os<<pNames(pNames.indexmax())<<". ";
    } else {os<<"no prior params. ";}
    if (pMPI->getNumConsts()){
        os<<"prior consts are ";
        adstring_array cNames = pMPI->getNamesForConsts();
        for (int c=cNames.indexmin();c<=cNames.indexmax();c++) os<<cNames(c)<<tb;
        os<<cNames(cNames.indexmax())<<". ";
    } else {os<<"no prior consts. ";}
}

/**
 * Writes the parameter info to an output stream as an R list.
 * 
 * @param os - the output stream object, by reference
 */
void NumberInfo::writeToR(ostream& os){
    os<<"list(";
    writeToR1(os);
    os<<")";
}

/**
 * Writes the parameter info to an output stream as members of an R list.
 * 
 * @param os - the output stream object, by reference
 */
void NumberInfo::writeToR1(ostream& os){
    os<<"label='"<<label<<qt<<cc;
    os<<"scaleType='"<<tcsam::getScaleType(scaleType)<<qt<<cc;
    os<<"initVal="<<initVal<<cc;
    os<<"finalVal="<<finlVal<<cc;
    os<<"phase="<<phase<<cc;
    os<<"resample="<<qt<<wts::getOnOffType(resample)<<qt<<cc;
    os<<"priorWgt="<<priorWgt<<cc;
    os<<"pdfType=list(type='"<<priorType<<"',";
    if (pMPI->getNumParams()) {
        int N=pMPI->getNumParams();
        adstring_array names=pMPI->getNamesForParams();
        os<<"params=list("; for(int n=1;n<N;n++) os<<names(n)<<"="<<priorParams(n)<<",";os<<names(N)<<"="<<priorParams(N)<<"),";
    } else {os<<"params=NULL,";}
    if (pMPI->getNumConsts()) {
        int N=pMPI->getNumConsts();
        adstring_array names=pMPI->getNamesForConsts();
        os<<"consts=list("; for(int n=1;n<N;n++) os<<names(n)<<"="<<priorConsts(n)<<",";os<<names(N)<<"="<<priorConsts(N)<<")";
    } else {os<<"consts=NULL)";}
}
////////////////////////////NumberInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedNumberInfo
 *----------------------------------------------------------------------------*/
/**
 * Sets the initial value for the parameter.
 * 
 * @param x - the initial value to set
 */
void BoundedNumberInfo::setInitVal(double x){
    if (debug) rpt::echo<<"starting BoundedNumberInfo::setInitVal(double x)"<<this<<endl;
    if (x<lower) {initVal = lower+(upper-lower)/1000000.0;} else
    if (x>upper) {initVal = upper-(upper-lower)/1000000.0;} else
    {initVal=x;}
    if (debug) rpt::echo<<"finished BoundedNumberInfo::setInitVal(double x)"<<this<<endl;
}

/**
 * Calculate bounded arithmetic-scale value corresponding to parameter-scale value.
 * 
 * @param x - parameter-scale value
 * @return - arithmetic-scale value
 */
double BoundedNumberInfo::calcArithScaleVal(double x){
    if (debug) rpt::echo<<"starting BoundedNumberInfo::calcArithScaleVal(double&) "<<this<<endl;
    double y, z;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            y = 1.0/(1.0+exp(-x));//scaled to [0,1]
            z = lower+(upper-lower)*y;
            break;
        default:
            cout<<"BoundedNumberInfo::calcArithScaleVal(double&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"starting BoundedNumberInfo::calcArithScaleVal(double&) "<<this<<endl;
    return z;
}

/**
 * Calculates the arithmetic-scale value corresponding to given parameter-scale value.
 * 
 * @param x - the parameter-scale value
 * @return - the corresponding arithmetic-scale value, as a dvariable
 * 
 * @overrides NumberInfo::calcArithScaleVal(const prevariable& x)
 */
dvariable BoundedNumberInfo::calcArithScaleVal(const prevariable& x){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting BoundedNumberInfo::calcArithScaleVal(prevariable&) "<<this<<endl;
    dvariable y, z;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            y = 1.0/(1.0+exp(-x));//scaled to [0,1]
            z = lower+(upper-lower)*y;
            break;
        default:
            cout<<"BoundedNumberInfo::calcArithScaleVal(dvariable&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedNumberInfo::calcArithScaleVal(prevariable&) "<<this<<endl;
    RETURN_ARRAYS_DECREMENT();
    return z;
}

/**
 * Calculate logit-scale value corresponding to given bounded arithmetic-scale value.
 * 
 * @param x - arithmetic-scale value
 * @return - logit-scale value
 */
double BoundedNumberInfo::calcParamScaleVal(double x){
    if (debug) rpt::echo<<"starting BoundedNumberInfo::calcParamScaleVal(double&) "<<this<<endl;
    double y,z;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = log(x);
            break;
        case tcsam::SCALE_LOGIT:
            if (x<=lower) return -25.0;
            if (x>=upper) return 25.0;
            y = (x-lower)/(upper-lower);//scaled from [lower,upper] to [0,1]
            z = log(y/(1.0-y));//on logit scale
            break;
        default:
            cout<<"BoundedNumberInfo::calcParamScaleVal(double) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedNumberInfo::calcParamScaleVal(double&) "<<this<<endl;
    return z;
}
            
/**
 * Draws a jittered random value on the arithmetic scale based on the 
 * bounds set and the fraction of the range to jitter (jitFrac).
 * 
 * Note that this DOES NOT update initVal.
 * If the estimation phase is \< 0, the value of initVal is returned
 * 
 * @param rng - the random number generator
 * @param jitFrac - the fraction of the range across which to jitter
 * 
 * @return - the random number on the arithmetic scale, or the value of initVal
 */
double BoundedNumberInfo::jitterInitVal(random_number_generator& rng, double jitFrac){
    double vp = initVal;
    if (phase>0){
        double d = upper-lower;
        double r = rng.better_rand();
        vp = min(max(lower+0.0001*d,initVal+wts::min(1.0,jitFrac)*(r-0.5)*d),upper-0.0001*d);
    }
    return vp;
}

/**
 * Reads the parameter info in ADMB format from an input stream.
 * 
 * The read order for the parameter info is:
 * <ul>
 *  <li> jitter
 *  <li> lower bound, on arithmetic scale
 *  <li> upper bound, on arithmetic scale
 *  <li> initVal
 *  <li> scaleType (as an adstring)
 *  <li> phase
 *  <li> resample
 *  <li> priorWgt
 *  <li> priorType
 *  <li> priorParams
 *  <li> priorConsts
 *  <li> label
 * </ul>
 * 
 * @param is - the input stream
 * 
 */
void BoundedNumberInfo::read(cifstream & is){
    if (debug) {
        rpt::echo<<"Starting BoundedNumberInfo::read(cifstream & is) "<<this<<endl;
        NumberInfo::debug=1;
    }
    adstring str;
    is>>str;   if (debug) rpt::echo<<str<<tb<<"#jitter"<<endl; jitter=wts::getOnOffType(str);
    is>>lower; if (debug) rpt::echo<<lower<<tb<<"#lower"<<endl;
    is>>upper; if (debug) rpt::echo<<upper<<tb<<"#upper"<<endl;
    NumberInfo::read(is);//finish reading
    if (debug) {
        NumberInfo::debug=0;
        rpt::echo<<"Done BoundedNumberInfo::read(cifstream & is) "<<this<<endl;
    }
}

/**
 * Write the parameter information to an output stream in ADMB format.
 * 
 * @param os - the output stream object
 */
void BoundedNumberInfo::write(ostream & os){
    os<<wts::getOnOffType(jitter)<<tb;
    os<<lower<<tb;
    os<<upper<<tb;
    NumberInfo::write(os);//finish writing
}

/**
 * Write the parameter information to an output stream as an R list.
 * 
 * @param os - the output stream object
 */
void BoundedNumberInfo::writeToR(ostream& os){
    os<<"list(";
    writeToR1(os);
    os<<")";
}

/**
 * Write the parameter information to an output stream as members of an R list.
 * 
 * @param os - the output stream object
 */
void BoundedNumberInfo::writeToR1(ostream& os){
    os<<"scaleType='"<<tcsam::getScaleType(scaleType)<<qt<<cc;
    os<<"lower="<<lower<<cc;
    os<<"upper="<<upper<<cc;
    os<<"jitter="<<qt<<wts::getOnOffType(jitter)<<qt<<cc;
    NumberInfo::writeToR1(os);
}
////////////////////////////BoundedNumberInfo/////////////////////////////////
/*------------------------------------------------------------------------------
 *  VectorInfo
 *-----------------------------------------------------------------------------*/
/**
 * Sets the prior type, based on an adstring value.
 * 
 * @param prior - the prior type, as an adstring
 */
void VectorInfo::setPriorType(adstring & prior){
    if (debug) rpt::echo<<"starting NumberInfo::setPriorType(adstring prior)"<<this<<endl;
    if (pMPI) delete pMPI;
    pMPI = ModelPDFInfo::getInfo(prior);
    if (pMPI) {
        if (pMPI->getNumParams()>0) priorParams.allocate(1,pMPI->getNumParams());
        if (pMPI->getNumConsts()>0) priorConsts.allocate(1,pMPI->getNumConsts());
    }
    if (debug) rpt::echo<<"finished NumberInfo::setPriorType(adstring prior)"<<this<<endl;
}

/**
 * Calculates arithmetic-scale value corresponding to the input parameter-scale value.
 * 
 * @param x - parameter-scale values as a double
 * @return - arithmetic-scale values as a double
 */
double VectorInfo::calcArithScaleVals(double x){
    if (debug) rpt::echo<<"starting VectorInfo::calcArithScaleVal(double) "<<this<<endl;
    double z = 0.0;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            z = 1.0/(1.0+mfexp(-1.0*x));
            break;
        default:
            cout<<"VectorInfo::calcArithScaleVal(double) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished VectorInfo::calcArithScaleVal(double) "<<this<<endl;
    return z;
}

/**
 * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
 * 
 * @param x - parameter-scale values as a dvector
 * @return - arithmetic-scale values as a dvector
 */
dvector VectorInfo::calcArithScaleVals(const dvector& x){
    if (debug) rpt::echo<<"starting VectorInfo::calcArithScaleVal(dvector&) "<<this<<endl;
    dvector z(x.indexmin(),x.indexmax());
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            z = 1.0/(1.0+mfexp(-1.0*x));
            break;
        default:
            cout<<"VectorInfo::calcArithScaleVal(dvector&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished VectorInfo::calcArithScaleVal(dvector&) "<<this<<endl;
    return z;
}

/**
 * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
 * 
 * @param x - parameter-scale values as dvar_vector
 * @return - arithmetic-scale values as dvar_vector
 */
dvar_vector VectorInfo::calcArithScaleVals(const dvar_vector& x){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting VectorInfo::calcArithScaleVal(dvar_vector&) "<<this<<endl;
    dvar_vector z(x.indexmin(),x.indexmax());
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            z = 1.0/(1.0+mfexp(-1.0*x));
            break;
        default:
            cout<<"VectorInfo::calcArithScaleVal(dvar_vector&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedVectorInfo::calcArithScaleVal(dvar_vector&) "<<this<<endl;
    RETURN_ARRAYS_DECREMENT();
    return z;
}
/**
 * Calculate parameter-scale value corresponding to given bounded arithmetic-scale value.
 * 
 * @param x - arithmetic-scale double
 * @return - parameter-scale double
 */
double VectorInfo::calcParamScaleVals(double x){
    if (debug) rpt::echo<<"starting VectorInfo::calcParamScaleVal(double) "<<this<<endl;
    double z = 0;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = log(x);//on log scale
            break;
        case tcsam::SCALE_LOGIT:
            z = log(x/(1.0-x));//on logit scale
            break;
        default:
            cout<<"VectorInfo::calcParamScaleVal(double) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished VectorInfo::calcParamScaleVal(double) "<<this<<endl;
    return z;
}

/**
 * Calculate parameter-scale values corresponding to given bounded arithmetic-scale values.
 * 
 * @param x - arithmetic-scale dvector
 * @return - parameter-scale dvector
 */
dvector VectorInfo::calcParamScaleVals(dvector& x){
    if (debug) rpt::echo<<"starting VectorInfo::calcParamScaleVal(dvector&) "<<this<<endl;
    dvector z(x.indexmin(),x.indexmax());
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = log(x);//on log scale
            break;
        case tcsam::SCALE_LOGIT:
            z = log(elem_div(x,1.0-x));//on logit scale
            break;
        default:
            cout<<"VectorInfo::calcParamScaleVal(dvector&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished VectorInfo::calcParamScaleVal(dvector&) "<<this<<endl;
    return z;
}
/**
 * Sets the vector of initial values to the input dvector element-by-element.
 * 
 * @param [in] x - a dvector of initial values to set. Indices should run 1:N
 */
void VectorInfo::setInitVals(dvector& x){initVals=1.0*x;}     

/**
 * Draws initial values based on resampling the prior.
 * 
 * If the parameter estimation phase is \< 0, a copy of the vector of initial 
 * values (initVals) is returned.
 * 
 * @param [in] rng - random number generator object
 * @param [in] vif - variance inflation factor
 * 
 * @return - a dvector of random (or initial) values. Indices run 1:N
 */
dvector VectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting VectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl  = 1.0*initVals;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        for (int i=1;i<=N;i++) smpl(i) = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        rpt::echo<<phase<<tb<<pMPI->canSample()<<endl;
        rpt::echo<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        rpt::echo<<"finished VectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    }
    return smpl;
}

/**
 * Calculates a vector of the log-scale prior probability based on an input vector of values.
 * 
 * The input values are presumably from the associated parameter values, 
 * but should be on the "arithmetic" scale.
 * 
 * @param x - the input dvar_vector on the arithmetic scale, presumably based on the associated parameters
 * 
 * @return - the log-scale prior, as a dvar_vector
 */
dvar_vector VectorInfo::calcLogPrior(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting VectorInfo::calcLogPrior(pv) for "<<name<<endl;
    dvar_vector lps;
    if (pMPI->canCalcLogPDF(pv)) {
        lps = pMPI->calcLogPDF(pv,priorParams,priorConsts);
    } else {
        lps.allocate(pv.indexmin(),pv.indexmax());
        lps.initialize();
        dvariable x;
        for (int i=pv.indexmin();i<=pv.indexmax();i++) {
            x = pv(i);//have to copy pv(i) )
            lps(i) = pMPI->calcLogPDF(x,priorParams,priorConsts);
        }
    }
    if (debug) {
        rpt::echo<<"logPrior: "<<lps<<endl;
        rpt::echo<<"finished VectorInfo::calcLogPriors(pv) for "<<name<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/**
 * Reads the parameter info from an input filestream.
 * The read order is:
 * <ul>
 * <li> idxType - the index type (as an adstring)
 * <li> the index block defining the vector indices (which determines N), as an IndexBlock
 * <li> readVals - an adstring flag to subsequently read a vector of initial values
 * <li> the values for a NumberInfo object (i.e., appropriate to NumberInfo::read(is))
 * </ul>
 * 
 * @param [in] is - the filestream to read from
 */
void VectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"Starting VectorInfo::read(cifstream & is) for "<<name<<endl;
    readPart1(is);
    readPart2(is);
}

/**
 * Reads the index type and index block for the VectorInfo.
 * 
 * @param is - the stream to read from
 */
void VectorInfo::readPart1(cifstream & is){
    if (debug) rpt::echo<<"Starting VectorInfo::readPart1(cifstream & is) for "<<name<<endl;
    adstring str;
    is>>idxType;
    if (debug) rpt::echo<<idxType<<tb<<"#index type"<<endl;
    int imn; int imx; 
    tcsam::getIndexLimits(idxType,imn,imx);
    ptrIB = new IndexBlock(imn,imx);
    is>>(*ptrIB);
    if (debug) rpt::echo<<(*ptrIB)<<tb<<"#index block"<<endl;
    N = ptrIB->getSize();
    initVals.allocate(1,N);
    is>>str; readVals = wts::getBooleanType(str);
    if (debug) {
        rpt::echo<<"idxType    = "<<idxType<<endl;
        rpt::echo<<"IndexBlock = "<<(*ptrIB)<<endl;
        rpt::echo<<"N          = "<<N<<endl;
        rpt::echo<<"readVals   = "<<wts::getBooleanType(readVals)<<endl;
        rpt::echo<<"Finished VectorInfo::readPart1(cifstream & is) for "<<name<<endl;
    }
}

void VectorInfo::readPart2(cifstream & is){
    if (debug) rpt::echo<<"Starting VectorInfo::readPart2(cifstream & is) for "<<name<<endl;
    adstring str;
    is>>initVal; if (debug) rpt::echo<<initVal<<tb<<"#initVal"<<endl;
    is>>str; scaleType=tcsam::getScaleType(str);
    is>>phase;
    is>>str; resample=wts::getOnOffType(str);
    is>>priorWgt;
    is>>priorType;//prior type
    if (debug){
        rpt::echo<<initVal<<tb<<"#initVal"<<endl;
        rpt::echo<<phase  <<tb<<"#phase"<<endl;
        rpt::echo<<wts::getOnOffType(resample)<<tb<<"#resample"<<endl;
        rpt::echo<<priorWgt<<tb<<"#priorWgt"<<endl;
        rpt::echo<<priorType<<tb<<"#prior type"<<endl;
    }
    setPriorType(priorType);
    if (pMPI->getNumParams()) {
        is>>priorParams; if (debug) rpt::echo<<priorParams<<tb<<"#prior params"<<endl;
    } else {
        is>>str; str.to_upper();
        if (str!=tcsam::STR_NONE){
            rpt::echo<<"no prior params, so input should be 'none', but got '"<<str<<"'"<<endl;
            rpt::echo<<"Please fix!!"<<endl;
            exit(-1);
        }
        if (debug) rpt::echo<<"#no prior params"<<endl;
    }
    if (pMPI->getNumConsts()) {
        is>>priorConsts; if (debug) rpt::echo<<priorConsts<<tb<<"#prior consts"<<endl;
    } else {
        is>>str; str.to_upper();
        if (str!=tcsam::STR_NONE){
            rpt::echo<<"no prior consts, so input should be 'none', but got '"<<str<<"'"<<endl;
            rpt::echo<<"Please fix!!"<<endl;
            exit(-1);
        }
        if (debug) rpt::echo<<"#no prior consts"<<endl;
    }
    is>>label; if (debug) rpt::echo<<label<<tb<<"#label"<<endl;
    initVals = initVal;
    if (debug) {
        rpt::echo<<"idxType = "<<idxType<<endl;
        rpt::echo<<"IndexBlock = "<<(*ptrIB)<<endl;
        rpt::echo<<"N     = "<<N<<endl;
        rpt::echo<<"init values = "<<initVals<<endl;
        rpt::echo<<"Finished VectorInfo::readPart2(cifstream & is) for "<<name<<endl;
    }
}

/**
 * Writes the parameter info to an output stream, in ADMB format.
 * 
 * @param [in] os - the output stream to write to
 */
void VectorInfo::write(std::ostream & os){
    writePart1(os);
    writePart2(os);
}

void VectorInfo::writePart1(std::ostream& os){
    os<<idxType<<tb;
    os<<(*ptrIB)<<tb;
    os<<wts::getBooleanType(readVals);
}

void VectorInfo::writePart2(std::ostream& os){
    os<<tcsam::getScaleType(scaleType)<<tb;
    os<<initVal<<tb;
    os<<phase<<tb;
    os<<wts::getOnOffType(resample)<<tb;
    os<<priorWgt<<tb;
    os<<priorType<<tb;
    if (pMPI->getNumParams()) os<<priorParams<<tb<<tb; else os<<"none"<<tb<<tb;
    if (pMPI->getNumConsts()) os<<priorConsts<<tb<<tb; else os<<"none"<<tb<<tb;
    os<<label<<tb<<tb;
    os<<"#";
    if (pMPI->getNumParams()){
        os<<"prior params are ";
        adstring_array pNames = pMPI->getNamesForParams();
        for (int p=pNames.indexmin();p<=pNames.indexmax()-1;p++) os<<pNames(p)<<cc;
        os<<pNames(pNames.indexmax())<<". ";
    } else {os<<"no prior params. ";}
    if (pMPI->getNumConsts()){
        os<<"prior consts are ";
        adstring_array cNames = pMPI->getNamesForConsts();
        for (int c=cNames.indexmin();c<=cNames.indexmax();c++) os<<cNames(c)<<tb;
        os<<cNames(cNames.indexmax())<<". ";
    } else {os<<"no prior consts. ";}
}
/**
 * Writes the parameter info to an output stream as an R list.
 * 
 * @param [in] os - the output stream to write to
 */
void VectorInfo::writeToR(ostream& os){
    if (debug) rpt::echo<<"VectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
    writeToR1(os);
    os<<")";
}

void VectorInfo::writeToR1(ostream& os){
    os<<"label='"<<label<<qt<<cc;
    os<<"scaleType='"<<tcsam::getScaleType(scaleType)<<qt<<cc;
    os<<"initVal="<<initVal<<cc;
    os<<"finalVal="<<finlVal<<cc;
    os<<"phase="<<phase<<cc;
    os<<"resample="<<qt<<wts::getOnOffType(resample)<<qt<<cc;
    os<<"priorWgt="<<priorWgt<<cc;
    os<<"pdfType=list(type='"<<priorType<<"',";
    if (pMPI->getNumParams()) {
        int N=pMPI->getNumParams();
        adstring_array names=pMPI->getNamesForParams();
        os<<"params=list("; for(int n=1;n<N;n++) os<<names(n)<<"="<<priorParams(n)<<",";os<<names(N)<<"="<<priorParams(N)<<"),";
    } else {os<<"params=NULL,";}
    if (pMPI->getNumConsts()) {
        int N=pMPI->getNumConsts();
        adstring_array names=pMPI->getNamesForConsts();
        os<<"consts=list("; for(int n=1;n<N;n++) os<<names(n)<<"="<<priorConsts(n)<<",";os<<names(N)<<"="<<priorConsts(N)<<")"<<cc<<endl;
    } else {os<<"consts=NULL)"<<cc<<endl;}
    os<<"initVals=";  wts::writeToR(os,initVals,wts::to_qcsv(ptrIB->getFwdIndexVector())); os<<cc<<endl;
    os<<"finalVals="; wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
}

/**
 * Writes the final values to an output stream as an R structure.
 * 
 * @param [in] os - the output stream to write to
 */
void VectorInfo::writeFinalValsToR(ostream& os){
    if (debug) rpt::echo<<"VectorInfo::writeFinalValsToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
}
////////////////////////////VectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedVectorInfo
 *----------------------------------------------------------------------------*/
/**
 * Calculates arithmetic-scale values corresponding to the input parameter-scale value.
 * 
 * @param x - parameter-scale value as a double
 * @return - arithmetic-scale value as a double
 * 
 * @overrides VectorInfo::calcArithScaleVals(double)
 */
double BoundedVectorInfo::calcArithScaleVals(double x){
   if (debug) rpt::echo<<"starting BoundedVectorInfo::calcArithScaleVals(double) "<<this<<endl;
    double z=0.0;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            z = lower+(upper-lower)/(1.0+mfexp(-1.0*x));
            break;
        default:
            cout<<"BoundedVectorInfo::calcArithScaleVal(double) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedVectorInfo::calcArithScaleVals(double) "<<this<<endl;
    return z;
}
/**
 * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
 * 
 * @param x - parameter-scale values as a dvector
 * @return - arithmetic-scale values as a dvector
 */
dvector BoundedVectorInfo::calcArithScaleVals(const dvector& x){
    if (debug) rpt::echo<<"starting BoundedVectorInfo::calcArithScaleVal(dvector&) "<<this<<endl;
    dvector z(x.indexmin(),x.indexmax());
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            z = lower+(upper-lower)/(1.0+mfexp(-1.0*x));
            break;
        default:
            cout<<"BoundedVectorInfo::calcArithScaleVal(dvector&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedVectorInfo::calcArithScaleVal(dvector&) "<<this<<endl;
    return z;
}

/**
 * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
 * 
 * @param x - parameter-scale values as dvar_vector
 * @return - arithmetic-scale values as dvar_vector
 */
dvar_vector BoundedVectorInfo::calcArithScaleVals(const dvar_vector& x){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting BoundedVectorInfo::calcArithScaleVal(dvar_vector&) "<<this<<endl;
    dvar_vector z(x.indexmin(),x.indexmax());
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = mfexp(x);
            break;
        case tcsam::SCALE_LOGIT:
            z = lower+(upper-lower)/(1.0+mfexp(-1.0*x));
            break;
        default:
            cout<<"VectorInfo::calcArithScaleVal(dvar_vector&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedVectorInfo::calcArithScaleVal(dvar_vector&) "<<this<<endl;
    RETURN_ARRAYS_DECREMENT();
    return z;
}
/**
 * Calculate parameter-scale value corresponding to given bounded arithmetic-scale value.
 * 
 * @param x - arithmetic-scale double
 * @return - parameter-scale double
 */
double BoundedVectorInfo::calcParamScaleVals(double x){
    if (debug) rpt::echo<<"starting BoundedVectorInfo::calcParamScaleVal(double) "<<this<<endl;
    double z = 0;
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = log(x);
            break;
        case tcsam::SCALE_LOGIT:
            {
                //copy x
                double xp = x;
                //replace x "near" bounds with values equivalent to z = -25 or +25
                if (xp<=lower) xp = lower+(upper-lower)/(1+mfexp( 25.0));
                if (xp>=upper) xp = lower+(upper-lower)/(1+mfexp(-25.0));
                double y = (xp-lower)/(upper-lower);//scaled from [lower,upper] to [0,1]
                z = log(y/(1.0-y));//on logit scale
                break;
            }
        default:
            cout<<"BoiundedVectorInfo::calcParamScaleVal(double) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedVectorInfo::calcParamScaleVal(double) "<<this<<endl;
    return z;
}
/**
 * Calculate parameter-scale values corresponding to given bounded arithmetic-scale values.
 * 
 * @param x - arithmetic-scale dvector
 * @return - parameter-scale dvector
 */
dvector BoundedVectorInfo::calcParamScaleVals(dvector& x){
    if (debug) rpt::echo<<"starting BoundedVectorInfo::calcParamScaleVal(dvector&) "<<this<<endl;
    dvector z(x.indexmin(),x.indexmax());
    switch (scaleType){
        case tcsam::SCALE_ARITHM:
            z = x;
            break;
        case tcsam::SCALE_LOG:
            z = log(x);
            break;
        case tcsam::SCALE_LOGIT:
            {
                //copy x
                dvector xp(x.indexmin(),x.indexmax()); xp = x;
                //replace x "near" bounds with values equivalent to z = -25 or +25
                for (int i=xp.indexmin();xp.indexmax();i++){
                    if (xp(i)<=lower) xp(i) = lower+(upper-lower)/(1+mfexp( 25.0));
                    if (xp(i)>=upper) xp(i) = lower+(upper-lower)/(1+mfexp(-25.0));
                }
                dvector y = (xp-lower)/(upper-lower);//scaled from [lower,upper] to [0,1]
                z = log(elem_div(y,1.0-y));//on logit scale
                break;
            }
        default:
            cout<<"BoundedVectorInfo::calcParamScaleVal(dvector&) "<<this<<endl;
            cout<<"Parameter scale type "<<scaleType<<" not recognized!!"<<endl;
            ad_exit(1);
            break;
    }
    if (debug) rpt::echo<<"finished BoundedVectorInfo::calcParamScaleVal(dvector&) "<<this<<endl;
    return z;
}
/**
 * Set initial values on the arithmetic scale.
 * 
 * @param x - double for initial values
 */
void BoundedVectorInfo::setInitVals(double x){
    if (debug) rpt::echo<<"starting BoundedVectorInfo::setInitVals(double x) for "<<name<<endl;
    initVals = x;
    for (int i=1;i<=N;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
    if (debug) {
        rpt::echo<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        rpt::echo<<"finished BoundedVectorInfo::setInitVals(double x) for "<<name<<endl;
    }
}
/**
 * Set initial values on the arithmetic scale.
 * 
 * @param x - dvector of initial values
 */
void BoundedVectorInfo::setInitVals(dvector& x){
    if (debug) {
        rpt::echo<<"starting BoundedVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
        rpt::echo<<"input  x index limits: "<<x.indexmin()<<cc<<x.indexmax()<<endl;
        rpt::echo<<"initVals index limits: "<<initVals.indexmin()<<cc<<initVals.indexmax()<<endl;
    }
    initVals = 1.0*x;
    for (int i=1;i<=N;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
    if (debug) {
        rpt::echo<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        rpt::echo<<"finished BoundedVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    }
}

/**
 * Read initial values on the arithmetic scale from a stream.
 * 
 * @param is - the stream to read from
 */
void BoundedVectorInfo::readInitVals(cifstream & is){
    is>>initVals;
    for (int i=1;i<=N;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
}

/**
 * Draw a random sample from a prior defined on the arithmetic scale. \n
 * If phase<0, return initVals rather than resampling.
 * 
 * @param rng - random number generator
 * @param vif - variance inflation factor
 * 
 * @return - dvector of resampled values on the arithmetic scale
*/
dvector BoundedVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting BoundedVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl  = 1.0*initVals;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        for (int i=1;i<=N;i++) smpl(i) = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        rpt::echo<<phase<<tb<<pMPI->canSample()<<endl;
        rpt::echo<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        rpt::echo<<"finished BoundedVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    }
    return smpl;
}

/**
 * Calculates a vector of the log-scale prior probability based on an input vector of values.
 * 
 * The input values are presumably from the associated parameter values, 
 * but should be on the "arithmetic" scale.
 * 
 * @param x - the input dvar_vector on the arithmetic scale, presumably based on the associated parameters
 * 
 * @return - the log-scale prior, as a dvar_vector
 */
dvar_vector BoundedVectorInfo::calcLogPrior(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting BoundedVectorInfo::calcLogPrior(pv) for "<<name<<endl;
    dvar_vector lps;
    if (pMPI->canCalcLogPDF(pv)) {
        lps = pMPI->calcLogPDF(pv,priorParams,priorConsts);
    } else {
        lps.allocate(pv.indexmin(),pv.indexmax());
        lps.initialize();
        dvariable x;
        for (int i=pv.indexmin();i<=pv.indexmax();i++) {
            x = pv(i);//have to copy pv(i) )
            lps(i) = pMPI->calcLogPDF(x,priorParams,priorConsts);
        }
    }
    if (debug) {
        rpt::echo<<"logPrior: "<<lps<<endl;
        rpt::echo<<"finished BoundedVectorInfo::calcLogPriors(pv) for "<<name<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/**
 * Draws a vector of jittered random values on the arithmetic scale based on the 
 * bounds set and the fraction of the range to jitter (jitFrac).
 * 
 * Note that this DOES NOT update initVals.
 * If the estimation phase is \< 0, the value of initVals is returned
 * 
 * @param rng - the random number generator
 * @param jitFrac - the fraction of the range across which to jitter
 * 
 * @return - dvector of random numbers on the arithmetic scale, or initVals
 */
dvector BoundedVectorInfo::jitterInitVals(random_number_generator& rng, double jitFrac){
    dvector vp = 1.0*initVals;
    if (phase>0){
        double d = upper-lower;
        for (int i=1;i<=N;i++){
            double r = rng.better_rand();
            vp(i) = min(max(lower+0.0001*d,initVals(i)+wts::min(1.0,jitFrac)*(r-0.5)*d),upper-0.0001*d);
        }
    }
    return vp;
}

/**
 * Reads the parameter info from an input filestream.
 * The read order is:
 * <ul>
 *  <li> idxType - the index type (as an adstring)
 *  <li> the index block defining the vector indices (which determines N), as an IndexBlock
 *  <li> readVals - an adstring flag to subsequently read a vector of initial values
 *  <li> lower bound - the lower bound on the arithmetic scale
 *  <li> upper bound - the upper bound on the arithmetic scale
 *  <li> jitter - adstring flag to "jitter" initial values
 *  <li> scaleType
 *  <li> initVal
 *  <li> phase
 *  <li> resample
 *  <li> priorWgt
 *  <li> priorType
 *  <li> priorParams
 *  <li> priorConsts
 *  <li> label
 * </ul>
 * 
 * @param [in] is - the filestream to read from
 * 
 * @overrides VectorInfo::read(cifstream& is)
 */
void BoundedVectorInfo::read(cifstream & is){
    if (debug) {
        rpt::echo<<"Starting BoundedVectorInfo::read(cifstream & is) for "<<name<<endl;
        VectorInfo::debug = 1;
    }
    VectorInfo::readPart1(is);
    adstring str;
    is>>str; jitter=wts::getOnOffType(str);
    is>>lower;
    is>>upper;
    if (debug){
        rpt::echo<<wts::getOnOffType(jitter)<<tb<<"#jitter"<<endl;
        rpt::echo<<lower<<tb<<"#lower"<<endl;
        rpt::echo<<upper<<tb<<"#upper"<<endl;
    }
    VectorInfo::readPart2(is);
    initVals = initVal;
    if (debug) {
        VectorInfo::debug = 0;
        rpt::echo<<"Done BoundedVectorInfo::read(cifstream & is) for "<<name<<endl;
    }
}

/**
 * Write parameter information to output stream in ADMB format.
 * 
 * @param os - output stream to write to
 */
void BoundedVectorInfo::write(ostream & os){
    VectorInfo::writePart1(os);
    os<<wts::getOnOffType(jitter)<<tb;
    os<<lower<<tb;
    os<<upper<<tb;
    VectorInfo::writePart2(os);
}

/**
 * Write parameter information to output stream in R format.
 * 
 * @param os - output stream
 */
void BoundedVectorInfo::writeToR(ostream& os){
    if (debug) rpt::echo<<"BoundedVectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
        os<<"jitter='"<<wts::getOnOffType(jitter)<<"',lower="<<lower<<",upper="<<upper<<cc;
        VectorInfo::writeToR1(os);
    os<<")";
}

/**
 * Writes final values to an output stream as an R vector.
 * 
 * @param os - the output stream.
 */
void BoundedVectorInfo::writeFinalValsToR(ostream& os){
    if (debug) rpt::echo<<"starting BoundedVectorInfo::writeFinalValsToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
    if (debug) rpt::echo<<"finished BoundedVectorInfo::writeFinalValsToR for "<<this->name<<endl;
}

////////////////////////////BoundedVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  DevsVectorInfo
 *----------------------------------------------------------------------------*/
/**
 * Set initial values on the arithmetic scale.
 * 
 * @param x - dvector of initial values
 */
void DevsVectorInfo::setInitVals(dvector& x){
    if (debug) {
        rpt::echo<<"starting DevsVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
        BoundedVectorInfo::debug=1;
    }
    BoundedVectorInfo::setInitVals(x);
    initVals -= mean(initVals);
    if (max(fabs(initVals))>max(fabs(lower),fabs(upper))) initVals *= max(fabs(lower),fabs(upper))/max(fabs(initVals));
    if (debug) {
        BoundedVectorInfo::debug=0;
        rpt::echo<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        rpt::echo<<"finished BoundedVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    }
}

/**
 *   Draw a random sample on the arithmetic scale from the prior.\n 
 *   If phase<0, return initVals rather than resampling. \n
 * 
 * @param rng - random number generator
 * @param vif - variance inflation factor
*/
dvector DevsVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting DevsVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl = initVals;
    if (phase>0){
        smpl = BoundedVectorInfo::drawInitVals(rng,vif);
        smpl -= mean(smpl);//ensure sum of devs = 0
        if (max(fabs(smpl))>max(fabs(lower),fabs(upper))) smpl *= max(fabs(lower),fabs(upper))/max(fabs(smpl));
    }
    if (debug) {
        rpt::echo<<phase<<tb<<pMPI->canSample()<<endl;
        rpt::echo<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        rpt::echo<<"finished DevsVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    }
    return smpl;
}

/**
 * Draws a vector of jittered random values on the arithmetic scale based on the 
 * bounds set and the fraction of the range to jitter (jitFrac).
 * 
 * Note that this DOES NOT update initVals.
 * If the estimation phase is \< 0, the value of initVals is returned
 * 
 * @param rng - the random number generator
 * @param jitFrac - the fraction of the range across which to jitter
 * 
 * @return - the vector of jittered values on the arithmetic scale, or the value of initVals
 */
dvector DevsVectorInfo::jitterInitVals(random_number_generator& rng, double jitFrac){
    if (debug){
        rpt::echo<<"--in DevsVectorInfo::jitterInitVals"<<endl;
        rpt::echo<<"initVals   = "<<initVals<<endl;
    }
    dvector vp = 1.0*initVals;
    if (phase>0){
        if (debug) rpt::echo<<"lower = "<<lower<<tb<<"upper = "<<upper<<endl;
        vp = BoundedVectorInfo::jitterInitVals(rng,jitFrac);
        if (debug) rpt::echo<<"vp         = "<<vp<<endl;
        vp -= mean(vp);//ensure sum of devs = 0
        if (debug) rpt::echo<<"vp-mean(vp) = "<<vp<<endl<<"sum(vpp) = "<<sum(vp)<<endl;
        if (max(fabs(vp))>max(fabs(lower),fabs(upper))) {
            double rescl = max(fabs(lower),fabs(upper))/max(fabs(vp));
            if (debug) rpt::echo<<"rescaling using "<<rescl<<endl;
            vp *= rescl;//rescale to within bounds
        }        
    }
    if (debug) {
        rpt::echo<<"final vp = "<<vp<<endl;
        rpt::echo<<"------------"<<endl;
    }
    return vp;
}

/*------------------------------------------------------------------------------
 *  NumberVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   deallocation                                               *
***************************************************************/
void NumberVectorInfo::deallocate(void){
    if (debug) rpt::echo<<"starting NumberVectorInfo::deallocate(void) "<<this<<endl;
    if (ppNIs) {
        for (int p=0;p<nNIs;p++) if (ppNIs[p]!=0) delete ppNIs[p];
        delete[] ppNIs;
        ppNIs = 0;
    }
    if (debug) rpt::echo<<"finished NumberVectorInfo::deallocate(void) "<<this<<endl;
}

/**
 * Gets the parameter scale types for the associated parameters.
 * 
 * @return - an adstring_array of scale types
 */
adstring_array NumberVectorInfo::getScaleTypes(void){
    if (debug) rpt::echo<<"starting NumberVectorInfo::getScaleTypes(void) "<<this<<endl;
    adstring_array types(1,nNIs);
    if (ppNIs) for (int i=1;i<=nNIs;i++) types(i) = ppNIs[i-1]->getScaleType();
    if (debug) rpt::echo<<"finished NumberVectorInfo::getScaleTypes(void) "<<this<<endl;
    return types;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector NumberVectorInfo::getPhases(void){
    if (debug) rpt::echo<<"starting NumberVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nNIs);
    phases.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) phases(i) = ppNIs[i-1]->getPhase();
    if (debug) rpt::echo<<"finished NumberVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector NumberVectorInfo::getPriorWgts(){
    if (debug) rpt::echo<<"starting NumberVectorInfo::getPriorWgts()"<<this<<endl;
    dvector wts(1,nNIs);
    for (int i=1;i<=nNIs;i++) wts(i) = ppNIs[i-1]->getPriorWgt();
    if (debug) rpt::echo<<"finished NumberVectorInfo::getPriorWgts()"<<this<<endl;
    return wts;
}

/**
* Calculates arithmetic-scale values corresponding to the input parameter-scale values.
* 
* @param x - parameter-scale values as a dvector
* @return - arithmetic-scale values as a dvector
*/
dvector NumberVectorInfo::calcArithScaleVals(const dvector& x){
    if (debug) rpt::echo<<"starting NumberVectorInfo::calcArithScaleVals(dvector&) "<<this<<endl;
    dvector asv(1,nNIs);
    asv.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) asv(i) = (ppNIs[i-1])->calcArithScaleVal(x(i));
    if (debug) rpt::echo<<"finished NumberVectorInfo::calcArithScaleVals(dvector&) "<<this<<endl;
    return asv;
}

/**
* Calculates arithmetic-scale values corresponding to the input parameter-scale values.
* 
* @param x - parameter-scale values as dvar_vector
* @return - arithmetic-scale values as dvar_vector
*/
dvar_vector NumberVectorInfo::calcArithScaleVals(const dvar_vector& x){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting NumberVectorInfo::calcArithScaleVals(dvar_vector&) "<<this<<endl;
    dvar_vector asv(1,nNIs);
    asv.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) asv(i) = (ppNIs[i-1])->calcArithScaleVal(x(i));
    if (debug) rpt::echo<<"finished NumberVectorInfo::calcArithScaleVals(dvar_vector&) "<<this<<endl;
    RETURN_ARRAYS_DECREMENT();
    return asv;
}

/**
* Calculates parameter-scale values corresponding to the input arithmetic-scale values.
* 
* @param x - the arithmetic-scale values as dvector
* @return - the parameter-scale values as dvector
*/
dvector NumberVectorInfo::calcParamScaleVals(dvector& x){
    if (debug) rpt::echo<<"starting NumberVectorInfo::calcParamScaleVals(dvector&) "<<this<<endl;
    dvector asv(1,nNIs);
    asv.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) asv(i) = (ppNIs[i-1])->calcParamScaleVal(x(i));
    if (debug) rpt::echo<<"finished NumberVectorInfo::calcParamScaleVals(dvector&) "<<this<<endl;
    return asv;
}

/**
 * Gets a dvector of the initial values for each parameter in the 
 * associated param_init_number_vector.
 * 
 * @return - a dvector
 */
dvector NumberVectorInfo::getInitVals(){
    if (debug) rpt::echo<<"starting NumberVectorInfo::getInitVals()"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = ppNIs[i-1]->getInitVal();
    if (debug) rpt::echo<<"finished NumberVectorInfo::getInitVals()"<<this<<endl;
    return initVals;
}

/**
 * Gets a dvector of the initial values for each parameter in the 
 * associated param_init_number_vector.
 * 
 * @return - a dvector
 */
dvector NumberVectorInfo::getInitValsOnParamScales(){
    if (debug) rpt::echo<<"starting NumberVectorInfo::getInitValsOnParamScales()"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = ppNIs[i-1]->getInitValOnParamScale();
    if (debug) rpt::echo<<"finished NumberVectorInfo::getInitValsOnParamScales()"<<this<<endl;
    return initVals;
}

/**
 * Sets initial values on arithmetic scale from corresponding parameter values, which may
 * not be one the same scale. Primarily used when parameter values are set
 * using a pin file.
 * 
 * @param x - a dvar_vector with the parameter-scale values to set
 */
void NumberVectorInfo::setInitValsFromParamVals(const dvar_vector& x){
    if (debug) rpt::echo<<"starting NumberVectorInfo::setInitValsFromParamVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) ppNIs[i-1]->setInitValFromParamVal(x(i));
    if (debug) rpt::echo<<"finished NumberVectorInfo::setInitValsFromParamVals(x)"<<this<<endl;
}

/**
 * Sets final values on arithmetic scale from corresponding parameter values, which may
 * not be one the same scale. Primarily used for output to R.
 * 
 * @param x - a dvar_vector with the parameter-scale values to set
 */
void NumberVectorInfo::setFinalValsFromParamVals(const dvar_vector& x){
    if (debug) rpt::echo<<"starting NumberVectorInfo::setFinalVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) ppNIs[i-1]->setFinalValFromParamVal(x(i));
    if (debug) rpt::echo<<"finished NumberVectorInfo::setFinalVals(x)"<<this<<endl;
}

/***************************************************************
*   Draw initial parameter values.                             *
***************************************************************/
dvector NumberVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting NumberVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = ppNIs[i-1]->drawInitVal(rng,vif);
    if (debug) rpt::echo<<"finished NumberVectorInfo::drawInitVals(random_number_generator& rng)"<<this<<endl;
    return initVals;
}

/**
 * Calculates a vector of the log-scale prior probability based on an input vector of values.
 * 
 * The input values are presumably from the associated parameter values, 
 * but should be on the "arithmetic" scale.
 * 
 * @param x - the input dvar_vector on the arithmetic scale, presumably based on the associated parameters
 * 
 * @return - the log-scale prior, as a dvar_vector
 */
dvar_vector NumberVectorInfo::calcLogPriors(dvar_vector & pv){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"starting NumberVectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
    dvar_vector lps(pv.indexmin(),pv.indexmax());
    lps.initialize();
    if (ppNIs) {
        dvariable x;
        for (int i=pv.indexmin();i<=pv.indexmax();i++) {
            x = pv(i);//have to copy pv(i) )
            lps(i) = ppNIs[i-1]->calcLogPrior(x);
        }
    }
    if (debug) {
        rpt::echo<<"logPriors: "<<lps<<endl;
        rpt::echo<<"finished NumberVectorInfo::calcLogPriors(pv)"<<this<<tb<<name<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void NumberVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"starting NumberVectorInfo::read(cifstream & is) "<<this<<endl;
    is>>nNIs;
    if (debug) rpt::echo<<"reading info for parameter vector"<<tb<<name<<tb<<nNIs<<endl;
    if (ppNIs) deallocate();
    if (nNIs>0) {
        ppNIs = new NumberInfo*[nNIs];
        int idx;
        for (int p=0;p<nNIs;p++) {
            is>>idx;
            ppNIs[idx] = new NumberInfo();
            is>>(*ppNIs[idx]);
        }
        if (debug) {
            for (int p=0;p<nNIs;p++) rpt::echo<<(*ppNIs[p])<<endl;
        }
    }
    if (debug) rpt::echo<<"finished NumberVectorInfo::read(cifstream & is) "<<this<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void NumberVectorInfo::write(ostream & os){
    os<<tb<<nNIs<<"  #number of parameters"<<endl;
    os<<"#     param initial       resample  prior   prior  prior   prior"<<endl;
    os<<"#id   scale  value  phase  values?  weight  type   params  consts  label"<<endl;
    if (nNIs){
        for (int p=0;p<(nNIs-1);p++) os<<(p+1)<<tb<<(*ppNIs[p])<<endl;
        os<<nNIs<<tb<<(*ppNIs[nNIs-1]);
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void NumberVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nNIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nNIs;p++) {os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nNIs;               os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}

/**
 * Writes final values to stream as an R list.
 * @param os - the stream to write to
 */
void NumberVectorInfo::writeFinalValsToR(ostream& os){
    if (nNIs){
        os<<"list("<<endl;
        for (int p=1;p<nNIs;p++) {os<<tb<<"'"<<p<<"'="<<ppNIs[p-1]->getFinalVal()<<","<<endl;}
        int p=nNIs;               os<<tb<<"'"<<p<<"'="<<ppNIs[p-1]->getFinalVal()<<")";
    } else {
        os<<"NULL";
    }
}

/**
 * Function to return labels for individual parameters as an
 * adstring_array.
 * 
 * @return adstring_array
 */
adstring_array NumberVectorInfo::getLabels(void){
    adstring_array arr(1,1);
    if (nNIs){
        adstring_array arr(1,nNIs);
        for (int i=1;i<=nNIs;i++) arr(i) = ppNIs[i-1]->label;
        return(arr);
    }
    return(arr);
}
////////////////////////////NumberVectorInfo/////////////////////////////////

/**----------------------------------------------------------------------------/n
 *  BoundedNumberVectorInfo
 *----------------------------------------------------------------------------*/
/**
 * Gets a pointer to the ith BoundedNumberInfo element in the vector. 
 * i must be in the interval 1<=i<=nNIs.
 * 
 * @param i - index (1-based) to BoundedNumberInfo
 * @return pointer to ith BoundedNumberInfo object
 */
BoundedNumberInfo* BoundedNumberVectorInfo::operator[](int i){
    if (ppNIs&&(i<=nNIs)) return (static_cast<BoundedNumberInfo*>(ppNIs[i-1])); 
    return 0;
}
            
/**
 * Get the lower bounds on the arithmetic scale for each parameter in the associated 
 * param_init_bounded_number_vector.
 * 
 * @return - a dvector with the lower bound on the arithmetic scale for each parameter
*/
dvector BoundedNumberVectorInfo::getLowerBounds(void){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector lbs(1,nNIs);
    lbs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) lbs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getLowerBound();
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::getLowerBounds(void) "<<this<<endl;
    return lbs;
}

/**
 * Gets the lower bound on the parameter scale for each parameter in the associated 
 * param_init_bounded_number_vector.
 * 
 * @return - a dvector with the lower bounds on the parameter scale for each parameter
*/
dvector BoundedNumberVectorInfo::getLowerBoundsOnParamScales(void){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::getLowerBoundsOnParamScales(void) "<<this<<endl;
    dvector lbs(1,nNIs);
    lbs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) lbs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getLowerBoundOnParamScale();
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::getLowerBoundsOnParamScales(void) "<<this<<endl;
    return lbs;
}

/**
 * Get the upper bounds on the arithmetic scale for each parameter in the associated 
 * param_init_bounded_number_vector.
 * 
 * @return - a dvector with the upper bound on the arithmetic scale for each parameter
*/
dvector BoundedNumberVectorInfo::getUpperBounds(void){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nNIs);
    ubs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) ubs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getUpperBound();
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

/**
 * Get the upper bounds on the parameter scale for each parameter in the associated 
 * param_init_bounded_number_vector.
 * 
 * @return - a dvector with the upper bound on the parameter scale for each parameter
*/
dvector BoundedNumberVectorInfo::getUpperBoundsOnParamScales(void){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nNIs);
    ubs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) ubs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getUpperBoundOnParamScale();
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

/**
 * Read parameters info from an input stream in ADMB format.
 * 
 * @param is - input stream
 */
void BoundedNumberVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::read(cifstream & is) "<<name<<endl;
    is>>nNIs; if (debug) rpt::echo<<"nNIs ="<<tb<<nNIs<<endl;
    if (ppNIs) deallocate();
    if (nNIs>0) {
        int idx;
        ppNIs = new NumberInfo*[nNIs];
        for (int p=0;p<nNIs;p++) {
            is>>idx; if (debug) rpt::echo<<"reading idx = "<<idx<<endl;
            if (idx<=nNIs){
                BoundedNumberInfo* pBNI = new BoundedNumberInfo();
                if (debug) BoundedNumberInfo::debug = 1;
                is>>(*pBNI);
                if (debug) BoundedNumberInfo::debug = 0;
                ppNIs[idx-1] = pBNI;
            } else {
                rpt::echo<<"Error reading file "<<is.get_file_name()<<endl;
                rpt::echo<<"for bounded parameter '"<<name<<"' defined for "<<nNIs<<" values."<<endl;
                rpt::echo<<"Tried to read "<<idx<<"th bounded number";
                cout<<"Aborting. See EchoOut.dat for details."<<endl;
                exit(-1.0);
            }
        }
        if (debug) {
            rpt::echo<<"TESTING read for BoundedNumberVectorInfo "<<name<<endl;
            for (int p=0;p<nNIs;p++) rpt::echo<<typeid(*ppNIs[p]).name()<<tb<<p+1<<tb<<(*ppNIs[p])<<endl;
        }
    }
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::read(cifstream & is) "<<name<<endl;
}

/**
 * Write parameters info to an output stream in ADMB format.
 * 
 * @param os - the output stream
 */
void BoundedNumberVectorInfo::write(ostream & os){
    if (debug) rpt::echo<<"Starting BoundedNumberVectorInfo::write(ostream & os) for "<<name<<endl;
    os<<tb<<nNIs<<"  #number of bounded parameters"<<endl;
    os<<"#           lower  upper  initial param        resample  prior   prior  prior   prior"<<endl;
    os<<"#id jitter? bounds bounds  value  scale  phase  values?  weight  type   params  consts  label"<<endl;
    if (nNIs){
        for (int p=0;p<(nNIs-1);p++) os<<(p+1)<<tb<<(*ppNIs[p])<<endl;
        os<<nNIs<<tb<<(*ppNIs[nNIs-1]);
    }
    if (debug) rpt::echo<<endl<<"Finished BoundedNumberVectorInfo::write(ostream & os)"<<endl;
}

/**
 * Write parameters info to stream as an R list.
 * 
 * @param os - the output stream
 * @param nm - the name of the R list
 * @param indent - number of tabs to indent
 */
void BoundedNumberVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nNIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nNIs;p++) {os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nNIs;               os<<tb<<"'"<<p<<"'="; ppNIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}
////////////////////////////BoundedNumberVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  VectorVectorInfo
 *----------------------------------------------------------------------------*/
/**
 * Deallocate memory by deleting all pointers, etc.
 */
void VectorVectorInfo::deallocate(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::deallocate(void) for "<<name<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) rpt::echo<<"finished VectorVectorInfo::deallocate(void) for "<<name<<endl;
}

/**
 * Get the minimum index for each parameter vector as an ivector.
 * 
 * @return - ivector of the minimum index for each vector
 */
ivector VectorVectorInfo::getMinIndices(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getMinIndices(void) for "<<name<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) idxs=1;
    if (debug) rpt::echo<<"finished VectorVectorInfo::getMinIndices(void) for "<<name<<endl;
    return idxs;
}

/**
 * Get the maximum index for each parameter vector as an ivector.
 * 
 * @return - ivector of the maximum index for each vector
 */
ivector VectorVectorInfo::getMaxIndices(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getMaxIndices(void) for "<<name<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getSize();
    if (debug) rpt::echo<<"finished VectorVectorInfo::getMaxIndices(void) for "<<name<<endl;
    return idxs;
}

/**
 * Gets the parameter-scale types for each associated parameter vector.
 * 
 * @return - an adstring_array
 */
adstring_array VectorVectorInfo::getScaleTypes(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getScaleTypes(void) for "<<name<<endl;
    adstring_array v(1,nVIs);
    if (ppVIs) for (int i=1;i<=nVIs;i++) v(i) = ppVIs[i-1]->getScaleType();
    if (debug) rpt::echo<<"finished VectorVectorInfo::getScaleTypes(void) for "<<name<<endl;
    return v;
}

/**
 * Get the starting estimation phase for each associated parameter vector, as an ivector.
 * 
 * @return - ivector of the start phase for each associated parameter vector
 */
ivector VectorVectorInfo::getPhases(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getPhases(void) for "<<name<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) rpt::echo<<"finished VectorVectorInfo::getPhases(void) for "<<name<<endl;
    return phases;
}

/**
 * Get the likelihood weights for the log prior probability associated with each parameter vector
 * 
 * @return - the vector of weights
 */
dvector VectorVectorInfo::getPriorWgts(){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getPriorWgts() for "<<name<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) rpt::echo<<"finished VectorVectorInfo::getPriorWgts() for "<<name<<endl;
    return wts;
}

/**
 * Function to return labels for individual parameters as an
 * adstring_array.
 * 
 * @return adstring_array
 */
adstring_array VectorVectorInfo::getLabels(void){
    adstring_array arr(1,1);
    if (nVIs){
        adstring_array arr(1,nVIs);
        for (int i=1;i<=nVIs;i++) arr(i) = (*this)[i]->label;
        return(arr);
    }
    return(arr);
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void VectorVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"starting VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
    is>>nVIs;
    if (debug) rpt::echo<<"reading info for parameter vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        ppVIs = new VectorInfo*[nVIs];
        int idx;
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            if ((idx>0)&&(idx<=nVIs)){
                ppVIs[idx-1] = new VectorInfo();
                is>>(*ppVIs[idx-1]);
            } else {
                rpt::echo<<"Error in VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
                rpt::echo<<"Error reading '"<<name<<"' from "<<is.get_file_name()<<endl;
                rpt::echo<<"Invalid index "<<idx<<". Must be >0 and <="<<nVIs<<endl;
                rpt::echo<<"Aborting..."<<endl;
                cout<<"Aborting. See EchoOut.dat for details."<<endl;
                exit(-1);
            }
        }
        //read vector values, if flagged
        for (int p=0;p<nVIs;p++) {
            VectorInfo* pVI = ppVIs[p];
            if (pVI->readVals) {
                if (debug) rpt::echo<<"Reading vector values for "<<p+1<<"th vector"<<endl;
                pVI->readInitVals(is);
                if (debug) rpt::echo<<pVI->getInitVals()<<endl;
            }
        }
        if (debug) {for (int p=0;p<nVIs;p++) rpt::echo<<(*ppVIs[p])<<endl;}
    }
    if (debug) rpt::echo<<"finished VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void VectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of parameters"<<endl;
    os<<"#     index index   read   initial param        resample  prior   prior  prior   prior"<<endl;
    os<<"#id   type  block  values? value   scale  phase  values?  weight  type   params  consts  label"<<endl;
    if (nVIs){
        for (int p=0;p<=(nVIs-1);p++) os<<(p+1)<<tb<<(*ppVIs[p])<<endl;
        os<<"#--initial values read in (index  values):";
        for (int p=0;p<=(nVIs-1);p++) {
            os<<endl<<(p+1)<<tb<<(*ppVIs[p]);
        }
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void VectorVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nVIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}

/**
 * Writes final values to an output stream as an R list of vector.
 * 
 * @param os - the output stream
 */
void VectorVectorInfo::writeFinalValsToR(ostream& os){
    if (nVIs){
        os<<"list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<")";
    } else {
        os<<"NULL";
    }
}
////////////////////////////VectorVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedVectorVectorInfo
 *----------------------------------------------------------------------------*/
/**
 * Gets a vector of the lower bounds (on the arithmetic scale) for the associated BoundedVectorInfo instances.
 * 
 * @return - a dvector of the lower bounds, on the arithmetic scale
 */
dvector BoundedVectorVectorInfo::getLowerBounds(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector bnds(1,nVIs);
    bnds.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) bnds(i) = (static_cast<BoundedVectorInfo*>(ppVIs[i-1]))->getLowerBound();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    return bnds;
}

/**
 * Gets a vector of the lower bounds (on the parameter vector scales) for the associated BoundedVectorInfo instances.
 * 
 * @return - a dvector of the lower bounds, on the parameter vector scales
 */
dvector BoundedVectorVectorInfo::getLowerBoundsOnParamScales(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getLowerBoundsOnParamScales(void) "<<this<<endl;
    dvector bnds(1,nVIs);
    bnds.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) bnds(i) = (static_cast<BoundedVectorInfo*>(ppVIs[i-1]))->getLowerBoundOnParamScale();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getLowerBoundsOnParamScales(void) "<<this<<endl;
    return bnds;
}

/**
 * Gets a vector of the upper bounds (on the arithmetic scale) for the associated BoundedVectorInfo instances.
 * 
 * @return - a dvector of the upper bounds
 */
dvector BoundedVectorVectorInfo::getUpperBounds(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector bnds(1,nVIs);
    bnds.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) bnds(i) = (static_cast<BoundedVectorInfo*>(ppVIs[i-1]))->getUpperBound();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    return bnds;
}

/**
 * Gets a vector of the upper bounds (on the parameter vector scales) for the associated BoundedVectorInfo instances.
 * 
 * @return - a dvector of the upper bounds, on the parameter vector scales
 */
dvector BoundedVectorVectorInfo::getUpperBoundsOnParamScales(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getUpperBoundsOnParamScales(void) "<<this<<endl;
    dvector bnds(1,nVIs);
    bnds.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) bnds(i) = (static_cast<BoundedVectorInfo*>(ppVIs[i-1]))->getUpperBoundOnParamScale();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getUpperBoundsOnParamScales(void) "<<this<<endl;
    return bnds;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void BoundedVectorVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::read(cifstream & is) "<<name<<endl;
    is>>nVIs;
    if (debug) rpt::echo<<"reading info for parameter vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        int idx;
        ppVIs = new VectorInfo*[nVIs];
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            if ((idx>0)&&(idx<=nVIs)){
                BoundedVectorInfo* pBVI = new BoundedVectorInfo();
                is>>(*pBVI);
            ppVIs[idx-1] = pBVI;
            } else {
                rpt::echo<<"Error in BoundedVectorVectorInfo::read(cifstream & is)"<<endl;
                rpt::echo<<"Error reading '"<<name<<"' from "<<is.get_file_name()<<endl;
                rpt::echo<<"Invalid index "<<idx<<". Must be >0 and <="<<nVIs<<endl;
                rpt::echo<<"Aborting..."<<endl;
                cout<<"Aborting. See EchoOut.dat for details"<<endl;
                exit(-1);
            }
        }
        //read vector values, if flagged
        for (int p=0;p<nVIs;p++) {
            BoundedVectorInfo* pVI = static_cast<BoundedVectorInfo*>(ppVIs[p]);
            if (pVI->readVals) {
                if (debug) rpt::echo<<"Reading vector values for "<<p+1<<"th bounded vector"<<endl;
                pVI->readInitVals(is);
                if (debug) rpt::echo<<pVI->getInitVals()<<endl;
            }
        }
        if (debug) {for (int p=0;p<nVIs;p++) rpt::echo<<p+1<<tb<<(*ppVIs[p])<<endl;}
    }
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::read(cifstream & is) "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedVectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of bounded parameters"<<endl;
    os<<"#     index index   read           lower  upper  initial param       resample  prior   prior  prior   prior"<<endl;
    os<<"#id   type  block  values? jitter? bounds bounds  value  scale phase  values?  weight  type   params  consts  label"<<endl;
    if (nVIs){
        for (int p=0;p<=(nVIs-1);p++) os<<(p+1)<<tb<<(*ppVIs[p])<<endl;
        os<<"#--initial values read in (index  values):";
        for (int p=0;p<=(nVIs-1);p++) {
            os<<endl<<(p+1)<<tb<<(*ppVIs[p]);
        }
    }
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
void BoundedVectorVectorInfo::writeToR(ostream& os, adstring nm, int indent){
    if (nVIs){
        os<<nm<<"=list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeToR(os); os<<")"<<endl;
    } else {
        os<<nm<<"=NULL";
    }
}

/**
 * Writes final values to an output stream as an R list of vectors.
 * 
 * @param os - the output stream
 */
void BoundedVectorVectorInfo::writeFinalValsToR(ostream& os){
    if (nVIs){
        os<<"list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<")";
    } else {
        os<<"NULL";
    }
}

////////////////////////////BoundedVectorVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  DevsVectorVectorInfo
 *----------------------------------------------------------------------------*/

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void DevsVectorVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::read(cifstream & is) "<<name<<endl;
    is>>nVIs;
    if (debug) rpt::echo<<"reading info for devs vector"<<tb<<name<<tb<<nVIs<<endl;
    if (ppVIs) deallocate();
    if (nVIs>0) {
        int idx;
        ppVIs = new VectorInfo*[nVIs];
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            if (debug) rpt::echo<<"reading index "<<idx<<endl;
            if ((idx>0)&&(idx<=nVIs)){
                DevsVectorInfo* pDVI = new DevsVectorInfo();
                if (debug) {DevsVectorInfo::debug=1; BoundedVectorInfo::debug=1;}
                is>>(*pDVI);
                if (debug) {DevsVectorInfo::debug=0; BoundedVectorInfo::debug=0;}
                ppVIs[idx-1] = pDVI;
            } else {
                rpt::echo<<"Error in DevsVectorVectorInfo::read(cifstream & is)"<<endl;
                rpt::echo<<"Error reading '"<<name<<"' from "<<is.get_file_name()<<endl;
                rpt::echo<<"Invalid index "<<idx<<". Must be >0 and <="<<nVIs<<endl;
                rpt::echo<<"Aborting..."<<endl;
                cout<<"Aborting. See EchoOut.dat for details"<<endl;
                exit(-1);
            }
        }
        //read vector values, if flagged
        for (int p=0;p<nVIs;p++) {
            DevsVectorInfo* pVI = static_cast<DevsVectorInfo*>(ppVIs[p]);
            if (pVI->readVals) {
                if (debug) rpt::echo<<"Reading vector values for "<<p+1<<"th devs vector"<<endl;
                pVI->readInitVals(is);
                if (debug) rpt::echo<<pVI->getInitVals()<<endl;
            }
        }
        if (debug) {
            for (int p=0;p<nVIs;p++) rpt::echo<<p+1<<tb<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::read(cifstream & is) "<<name<<endl;
}

////////////////////////////DevsVectorVectorInfo/////////////////////////////////
