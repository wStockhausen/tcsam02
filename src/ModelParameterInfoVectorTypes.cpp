//ModelParameterInfoVectorTypes.cpp
/*------------------------------------------------------------------------------
 *  Includes:
 *      VectorInfo
 *      BoundedVectorInfo
 *      DevsVectorInfo
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
#include "ModelParameterVectorInfoTypes.hpp"

int VectorInfo::debug              = 0;
int BoundedVectorInfo::debug       = 0;
int DevsVectorInfo::debug          = 0;
int VectorVectorInfo::debug        = 0;
int BoundedVectorVectorInfo::debug = 0;
int DevsVectorVectorInfo::debug    = 0;

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
 * Add a value to the end of the init and final vectors.
 * 
 * @param val - value to add
 * @param ibVal - corresponding IndexBlock value to add
 */
void VectorInfo::addValueOnArithmeticScale(double val, int ibVal){
    if (debug) cout<<"starting VectorInfo::addValueOnArithmeticScale("<<val<<cc<<ibVal<<")"<<endl;
    ptrIB->addElement(ibVal);
    dvector tmp = 1.0*initVals;
    N++;
    initVals.deallocate(); initVals.allocate(1,N);
    initVals(1,N-1) = tmp;
    initVals(N) = val;
    if (debug){
        cout<<"initVals orig: "<<tmp<<endl;
        cout<<"initVals updt: "<<initVals<<endl;
    }
    if (finlVals.allocated()){
        dvector tmp = 1.0*finlVals;
        finlVals.deallocate(); finlVals.allocate(1,N);
        finlVals(1,N-1) = tmp;
        finlVals(N) = val;
        if (debug){
            cout<<"finlVals orig: "<<tmp<<endl;
            cout<<"finlVals updt: "<<finlVals<<endl;
        }
    }
    if (debug) cout<<"finished VectorInfo::addValueOnArithmeticScale("<<val<<cc<<ibVal<<")"<<endl;
}

/**
 * Add a value to the end of the init and final vectors.
 * 
 * @param val - value (on parameter scale) to add
 * @param ibVal - corresponding IndexBlock value to add
 */
void VectorInfo::addValueOnParameterScale(double val, int ibVal){
    if (debug) cout<<"starting VectorInfo::addValueOnParameterScale("<<val<<cc<<ibVal<<")"<<endl;
    double aval = calcArithScaleVals(val);
    addValueOnArithmeticScale(aval,ibVal);
    if (debug) cout<<"finished VectorInfo::addValueOnParameterScale("<<val<<cc<<ibVal<<")"<<endl;
}

/**
 * Reads the parameter info from an input filestream.
 * The read order is:
 * <ul>
 * <li> idxType - the index type (as an adstring)
 * <li> the index block defining the vector indices (which determines N), as an IndexBlock
 * <li> readPhases - an adstring flag to subsequently read a vector of estimation phases
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
    if (debug) rpt::echo<<"Finished VectorInfo::read(cifstream & is) for "<<name<<endl;
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
    phases.allocate(1,N);
    initVals.allocate(1,N);
    is>>str; readPhases = wts::getBooleanType(str);
    is>>str; readVals   = wts::getBooleanType(str);
    if (debug) {
        rpt::echo<<"idxType    = "<<idxType<<endl;
        rpt::echo<<"IndexBlock = "<<(*ptrIB)<<endl;
        rpt::echo<<"N          = "<<N<<endl;
        rpt::echo<<"readPhases = "<<wts::getBooleanType(readPhases)<<endl;
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
    phases   = phase;
    initVals = initVal;//set default
    if (debug) {
        rpt::echo<<"idxType = "<<idxType<<endl;
        rpt::echo<<"IndexBlock = "<<(*ptrIB)<<endl;
        rpt::echo<<"N     = "<<N<<endl;
        rpt::echo<<"init phases = "<<phases<<endl;
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
    if (debug) cout<<"Starting VectorInfo::write()"<<endl;
    writePart1(os);
    writePart2(os);
    if (debug) cout<<"Finished VectorInfo::write()"<<endl;
}

void VectorInfo::writePart1(std::ostream& os){
    os<<idxType<<tb;
    os<<(*ptrIB)<<tb;
    os<<wts::getBooleanType(readPhases)<<tb;
    os<<wts::getBooleanType(readVals)  <<tb;
}

void VectorInfo::writePart2(std::ostream& os){
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
            cout<<"BoundedVectorInfo::calcArithScaleVal(dvar_vector&) "<<this<<endl;
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
            cout<<"BoundedVectorInfo::calcParamScaleVal(double) "<<this<<endl;
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
                for (int i=xp.indexmin();i<=xp.indexmax();i++){
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
            //vp(i) = min(max(lower+0.0001*d,initVals(i)+wts::min(1.0,jitFrac)*(r-0.5)*d),upper-0.0001*d);
            vp(i) = tcsam::jitterIt(initVals(i), lower, upper, jitFrac, r);
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
 *  <li> readPhases - an adstring flag to subsequently read a vector of phases
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
 * This calls VectorInfo::readPart1(...) and VectorInfo::readPart2(...).
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
    //initVals = initVal;  <-already done in readPart2
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
////////////////////////////DevsVectorInfo/////////////////////////////////

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

/**
 * Gets the ivector of phases for the associated param_init_number_vector (a
 * concatenation of the phases ivectors for all associated parameter vectors).
 * 
 * @return ivector(1,npT) of phases for each element of the associated param_init_number_vector
 */
ivector VectorVectorInfo::getParameterPhases(void){
    ivector phases(1,npT);
    int ctr = 1;
    for (int v=1;v<=nVIs;v++) {
        ivector phsv = (*this)[v]->getPhases();
        for (int j=phsv.indexmin();j<=phsv.indexmax();j++) phases(ctr++) = phsv(j);
    }
    return phases;
}

/**
 * Convenience method to set estimation phases for elements of all associated 
 * parameter vectors (and thus the associated param_init_number_vector).
 * 
 * @param phases - ivector(1,npT) to use to set phases
 */
void VectorVectorInfo::setParameterPhases(const ivector& phases){
    int ctr = 1;
    for (int v=1;v<=nVIs;v++) {
        ivector phsv = (*this)[v]->getPhases();
        for (int j=phsv.indexmin();j<=phsv.indexmax();j++) phsv(j) = phases(ctr++);
        (*this)[v]->setPhases(phsv);
    }
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
                ppVIs[idx-1]->read(is);
            } else {
                rpt::echo<<"Error in VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
                rpt::echo<<"Error reading '"<<name<<"' from "<<is.get_file_name()<<endl;
                rpt::echo<<"Invalid index "<<idx<<". Must be >0 and <="<<nVIs<<endl;
                rpt::echo<<"Aborting..."<<endl;
                cout<<"Aborting. See EchoOut.dat for details."<<endl;
                exit(-1);
            }
        }

        //read vectors for phases and/or initial values, if flagged
        readOptional(is);

        if (debug) {for (int p=0;p<nVIs;p++) rpt::echo<<(*ppVIs[p])<<endl;}
        
        //allocate mnIdxs, mxIdxs and determine npT
        allocateIndices();
    }
    if (debug) rpt::echo<<"finished VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
}

void VectorVectorInfo::readOptional(cifstream & is){
    //read vectors for phases and/or initial values, if flagged
    int idx;
    for (int p=0;p<nVIs;p++) {
        VectorInfo* pVI = ppVIs[p];
        if (pVI->readPhases) {
            if (debug) rpt::echo<<"Reading phases vector for "<<p+1<<"th vector"<<endl;
            is>>idx; 
            if (idx!=p+1){
               rpt::echo<<"Error reading initial vector phases for "<<p+1<<"th vector in "<<name<<endl;
               rpt::echo<<"Expected idx="<<p+1<<" but got "<<idx<<"."<<endl;
               rpt::echo<<"Check MPI file!"<<endl;
               exit(0);
            }
            pVI->readPhaseVector(is);
            if (debug) rpt::echo<<pVI->getPhases()<<endl;
        }
        if (pVI->readVals) {
            if (debug) rpt::echo<<"Reading vector values for "<<p+1<<"th vector"<<endl;
            is>>idx; 
            if (idx!=p+1){
                rpt::echo<<"Error reading initial vector values for "<<p+1<<"th vector in "<<name<<endl;
                rpt::echo<<"Expected idx="<<p+1<<" but got "<<idx<<"."<<endl;
                rpt::echo<<"Check MPI file!"<<endl;
                exit(0);
            }
            pVI->readInitVals(is);
            if (debug) rpt::echo<<pVI->getInitVals()<<endl;
        }
    }
}

void VectorVectorInfo::allocateIndices(){
    //allocate and set mnIdxs, mxIdxs
    mnIdxs.allocate(1,nVIs);
    mxIdxs.allocate(1,nVIs);
    for (int v=1;v<=nVIs;v++){
        mnIdxs[v] = 1;
        mxIdxs[v] = ppVIs[v-1]->getSize();
    }
    
    //determine npT, the total number of parameters
    npT = 0;//reset to zero, in case otherwise
    for (int v=1;v<=nVIs;v++) npT += mxIdxs(v)-mnIdxs(v)+1;        
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void VectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of parameters"<<endl;
    os<<"#     index index   read    read  initial  param        resample  prior   prior  prior   prior"<<endl;
    os<<"#id   type  block  phases? values? value   scale  phase  values?  weight  type   params  consts  label"<<endl;
    if (nVIs){
        for (int p=0;p<nVIs;p++) os<<(p+1)<<tb<<(*ppVIs[p])<<endl;
        os<<"#--id + phases and/or id + initial values:";
        for (int p=0;p<nVIs;p++) {
            if ((ppVIs[p])->readPhases) os<<endl<<p+1<<tb<<(ppVIs[p])->getPhases();
            if ((ppVIs[p])->readVals)   os<<endl<<p+1<<tb<<(ppVIs[p])->getInitVals();
        }
    }
}

/**
 * Write to output stream in ADMB pin-file format.
 * 
 * @param os
 */
void VectorVectorInfo::writeToPin(ostream & os){
    if (nVIs){
        for (int p=0;p<nVIs;p++) {
            os<<"#"<<name<<"["<<p+1<<"]:"<<endl;
            dvector tmp = ppVIs[p]->getFinalVals();
            os<<ppVIs[p]->calcParamScaleVals(tmp)<<endl;
        }
    } else {
        os<<"#"<<name<<"[0]:"<<endl<<0.00<<endl;
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
        //read vectors for phases and/or initial values, if flagged
        readOptional(is);
        
        if (debug) {for (int p=0;p<nVIs;p++) rpt::echo<<p+1<<tb<<(*ppVIs[p])<<endl;}
        
        //allocate mnIdxs, mxIdxs and determine npT
        allocateIndices();
    }
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::read(cifstream & is) "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedVectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of bounded parameters"<<endl;
    os<<"#     index index   read     read            lower upper initial  param        resample  prior   prior  prior   prior"<<endl;
    os<<"#id   type  block  phases?  values? jitter?  bound bound value   scale  phase  values?  weight  type   params  consts  label"<<endl;
    if (nVIs){
        for (int p=0;p<nVIs;p++) os<<(p+1)<<tb<<(*ppVIs[p])<<endl;
        os<<"#--id + phases and/or id + initial values:";
        for (int p=0;p<nVIs;p++) {
            if ((ppVIs[p])->readPhases) os<<endl<<p+1<<tb<<(ppVIs[p])->getPhases();
            if ((ppVIs[p])->readVals)   os<<endl<<(p+1)<<tb<<(ppVIs[p])->getInitVals();
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

        //read vectors for phases and/or initial values, if flagged
        readOptional(is);

        if (debug) for (int p=0;p<nVIs;p++) rpt::echo<<p+1<<tb<<(*ppVIs[p])<<endl;
        
        //allocate mnIdxs, mxIdxs and determine npT
        allocateIndices();
    }
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::read(cifstream & is) "<<name<<endl;
}

////////////////////////////DevsVectorVectorInfo/////////////////////////////////
