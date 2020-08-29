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
int NumberVectorInfo::debug        = 0;
int BoundedNumberVectorInfo::debug = 0;
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
    is>>initVal; finlVal=initVal;
    is>>str; scaleType=tcsam::getScaleType(str);
    if (debug) rpt::echo<<name<<tb<<str<<tb<<scaleType<<tb<<tcsam::getScaleType(scaleType)<<endl;
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

/**
 * Writes info to output stream.
 * 
 * Write order is in same sense as the read order:
 * <ul>
 *  <li> finlVal
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
    os<<finlVal<<tb;
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
    finlVal=initVal;
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
 * Calculate parameter-scale value corresponding to given bounded arithmetic-scale value.
 * 
 * @param x - arithmetic-scale value
 * @return - parameter-scale value
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
        //double d = upper-lower;
        double r = rng.better_rand();
        //vp = min(max(lower+0.0001*d,initVal+wts::min(1.0,jitFrac)*(r-0.5)*d),upper-0.0001*d);
        vp = tcsam::jitterIt(initVal, lower, upper, jitFrac, r);
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

/**
 * Gets an ivector of the estimation phase for each parameter in the
 * associated param_init_number_vector.
 * 
 * @return - an ivector
 */
ivector NumberVectorInfo::getPhases(void){
    if (debug) rpt::echo<<"starting NumberVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nNIs);
    phases.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) phases(i) = ppNIs[i-1]->getPhase();
    if (debug) rpt::echo<<"finished NumberVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/**
 * Sets an ivector for the estimation phase for each parameter in the
 * associated param_init_number_vector.
 * 
 * @param - an ivector of estimation phases
 * 
 * @return - nothing
 */
void NumberVectorInfo::setPhases(const ivector& _phases){
    if (debug) rpt::echo<<"starting NumberVectorInfo::setPhases(...) "<<this<<endl;
    if (ppNIs) for (int i=1;i<=nNIs;i++) ppNIs[i-1]->setPhase(_phases[i]);
    if (debug) rpt::echo<<"finished NumberVectorInfo::setPhases(...) "<<this<<endl;
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
 * not be on the same scale. Primarily used when parameter values are set
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

/**
 * Read from input stream in TCSAM02 format.
 * 
 * @param is
 */
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

/**
 * Write to output stream in TCSAM02 format.
 * 
 * @param os
 */
void NumberVectorInfo::write(ostream & os){
    os<<tb<<nNIs<<"  #number of parameters"<<endl;
    os<<"#     param initial       resample  prior   prior  prior   prior"<<endl;
    os<<"#id   scale  value  phase  values?  weight  type   params  consts  label"<<endl;
    if (nNIs){
        for (int p=0;p<(nNIs-1);p++) os<<(p+1)<<tb<<(*ppNIs[p])<<endl;
        os<<nNIs<<tb<<(*ppNIs[nNIs-1]);
    }
}

/**
 * Write to output stream in ADMB pin-file format.
 * 
 * @param os
 */
void NumberVectorInfo::writeToPin(ostream & os){
    if (nNIs){
        for (int p=0;p<nNIs;p++) {
            os<<"#"<<name<<"["<<p+1<<"]:"<<endl;
            double tmp = ppNIs[p]->getFinalVal();
            os<<ppNIs[p]->calcParamScaleVal(tmp)<<endl;
        }
    } else {
        os<<"#"<<name<<"[0]:"<<endl<<0.00<<endl;
    }
}

/**
 *   Write parameters info to stream in R format. 
 * 
 * @param os
 * @param nm
 * @param indent
 */
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

