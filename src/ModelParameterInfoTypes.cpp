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
/***************************************************************
*   Draw a random sample from the prior.                       * 
*   If phase<0, return initVal rather than resampling.         *
***************************************************************/
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

/******************************************************************
*   Calculate log prior probability.                              *
******************************************************************/
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

/***************************************************************
*   Set prior type.                                            *
***************************************************************/
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

/*********************************************\n
 * Read from the cifstream object.
 * Read order is:
 *  initVal, phase, resample, priorWgt, 
 *  priorType, priorParams, priorConsts
**********************************************/
void NumberInfo::read(cifstream & is){
    if (debug) rpt::echo<<"Starting NumberInfo::read(cifstream & is) "<<this<<endl;
    adstring str;
    is>>initVal;
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
void NumberInfo::write(ostream & os){
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
    os<<"label='"<<label<<"'"<<cc;
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
 * Draw a random sample based on the parameter's prior.
 * If the parameter's phase is < 0, the input initial value is returned.
 * 
 * @param rng - reference to random number generator
 * @param vif - variance inflation factor
 * 
 * @return random value, or the initial value
 */
double BoundedNumberInfo::drawInitVal(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting BoundedNumberInfo::drawSample(random_number_generator& rng)"<<this<<endl;
    double smpl  = initVal;
    if (resample&&(phase>0)&&(pMPI->canSample())) {
        smpl = pMPI->drawSample(rng,priorParams,priorConsts);
    }
    if (debug) {
        rpt::echo<<phase<<tb<<pMPI->canSample()<<tb<<initVal<<tb<<smpl<<endl;
        rpt::echo<<"finished BoundedNumberInfo::drawSample(random_number_generator& rng)"<<this<<endl;
    }
    return smpl;
}

/**
 * Read the parameter information from an input filestream in ADMB format.
 * 
 * @param is - the input filestream object
 */
void BoundedNumberInfo::read(cifstream & is){
    if (debug) rpt::echo<<"Starting BoundedNumberInfo::read(cifstream & is) "<<this<<endl;
    adstring str;
    is>>lower;
    is>>upper;
    is>>str; jitter=wts::getOnOffType(str);
    if (debug){
        rpt::echo<<lower<<tb<<"#lower"<<endl;
        rpt::echo<<upper<<tb<<"#upper"<<endl;
        rpt::echo<<wts::getOnOffType(jitter)<<tb<<"#jitter"<<endl;
    }
    NumberInfo::read(is);//finish reading
    if (debug) rpt::echo<<"Done BoundedNumberInfo::read(cifstream & is) "<<this<<endl;
}

/**
 * Write the parameter information to an output stream in ADMB format.
 * 
 * @param os - the output stream object
 */
void BoundedNumberInfo::write(ostream & os){
    os<<lower<<tb;
    os<<upper<<tb;
    os<<wts::getOnOffType(jitter)<<tb;
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
    os<<"lower="<<lower<<cc;
    os<<"upper="<<upper<<cc;
    os<<"jitter="<<qt<<wts::getOnOffType(jitter)<<qt<<cc;
    NumberInfo::writeToR1(os);
}
////////////////////////////BoundedNumberInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  VectorInfo\n
 *----------------------------------------------------------------------------*/
/***************************************************************
*   Draw a random sample from the prior.                \n 
*   If phase<0, return initVals rather than resampling. \n
***************************************************************/
dvector VectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting VectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl  = initVals;
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

/***************************************************************
*   Calculate log prior probability for each element of the\n
 *  parameter vector.                                      \n
***************************************************************/
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
            lps(i) = NumberInfo::calcLogPrior(x);
        }
    }
    if (debug) {
        rpt::echo<<"logPrior: "<<lps<<endl;
        rpt::echo<<"finished VectorInfo::calcLogPriors(pv) for "<<name<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*  Read from cifstream object.\n
 * Read order is:
 * mni, mxi, readVals + NumberInfo::read(is)
***************************************************************/
void VectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"Starting VectorInfo::read(cifstream & is) for "<<name<<endl;
    adstring str;
    is>>idxType;
    int imn; int imx; 
    tcsam::getIndexLimits(idxType,imn,imx);
    ptrIB = new IndexBlock(imn,imx);
    is>>(*ptrIB);
    N = ptrIB->getSize();
    initVals.allocate(1,N);
    is>>str; readVals = wts::getBooleanType(str);
    NumberInfo::read(is);
    initVals = initVal;
    if (debug) {
        rpt::echo<<"idxType = "<<idxType<<endl;
        rpt::echo<<"IndexBlock = "<<(*ptrIB)<<endl;
        rpt::echo<<"N     = "<<N<<endl;
        rpt::echo<<"init values = "<<initVals<<endl;
        rpt::echo<<"Finished VectorInfo::read(cifstream & is) for "<<name<<endl;
    }
}

/***************************************************************
*   write                                                      *
***************************************************************/
void VectorInfo::write(ostream & os){
    os<<idxType<<tb;
    os<<(*ptrIB)<<tb;
    os<<wts::getBooleanType(readVals);
    NumberInfo::write(os);
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void VectorInfo::writeToR(ostream& os){
    if (debug) rpt::echo<<"VectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
    NumberInfo::writeToR1(os); os<<cc<<endl;
    os<<"initVals=";  wts::writeToR(os,initVals,wts::to_qcsv(ptrIB->getFwdIndexVector())); os<<cc<<endl;
    os<<"finalVals="; wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
    os<<")";
}

/**
 * Writes final values to an output stream as an R vector
 * @param os - the output stream.
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
/***************************************************************\n
*   Sets initial values.       \n
***************************************************************/
void BoundedVectorInfo::setInitVals(dvector& x){
    if (debug) {
        rpt::echo<<"starting BoundedVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
        rpt::echo<<"input  x index limits: "<<x.indexmin()<<cc<<x.indexmax()<<endl;
        rpt::echo<<"initVals index limits: "<<initVals.indexmin()<<cc<<initVals.indexmax()<<endl;
    }
    initVals = x;
    for (int i=1;i<=N;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
    if (debug) {
        rpt::echo<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        rpt::echo<<"finished BoundedVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    }
}

/***************************************************************
*   read initial values. \n
***************************************************************/
void BoundedVectorInfo::readInitVals(cifstream & is){
    is>>initVals;
    for (int i=1;i<=N;i++) {
        if (initVals(i)<=lower) initVals(i) = lower+(upper-lower)/1000000.0; else
        if (initVals(i)>=upper) initVals(i) = upper-(upper-lower)/1000000.0; 
    }
}

/***************************************************************
*   Draw a random sample from the prior.                \n 
*   If phase<0, return initVals rather than resampling. \n
***************************************************************/
dvector BoundedVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting BoundedVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl  = initVals;
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

/***************************************************************
*   Calculate log prior probability for each element of the\n
 *  parameter vector.                                      \n
***************************************************************/
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
            lps(i) = BoundedNumberInfo::calcLogPrior(x);
        }
    }
    if (debug) {
        rpt::echo<<"logPrior: "<<lps<<endl;
        rpt::echo<<"finished BoundedVectorInfo::calcLogPriors(pv) for "<<name<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return lps;
}

/***************************************************************
*  Read from cifstream object.\n
 * Read order is:
 * mni, mxi, readVals + BoundedNumberInfo::read(is)
***************************************************************/
void BoundedVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"Starting BoundedVectorInfo::read(cifstream & is) for "<<name<<endl;
    adstring str;
    is>>idxType;
    int imn; int imx; 
    tcsam::getIndexLimits(idxType,imn,imx);
    ptrIB = new IndexBlock(imn,imx);
    is>>(*ptrIB);
    N = ptrIB->getSize();
    initVals.allocate(1,N);
    is>>str; 
    readVals = wts::getBooleanType(str);
    BoundedNumberInfo::read(is);
    initVals = initVal;
    if (debug) {
        rpt::echo<<"idxType    = "<<idxType<<endl;
        rpt::echo<<"IndexBlock = "<<(*ptrIB)<<endl;
        rpt::echo<<"N          = "<<N<<endl;
        rpt::echo<<"Done BoundedVectorInfo::read(cifstream & is) for "<<name<<endl;
    }
}

/***************************************************************
*   write                                                      *
***************************************************************/
void BoundedVectorInfo::write(ostream & os){
    os<<idxType<<tb;
    os<<(*ptrIB)<<tb;
    os<<wts::getBooleanType(readVals);
    BoundedNumberInfo::write(os);
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void BoundedVectorInfo::writeToR(ostream& os){
    if (debug) rpt::echo<<"BoundedVectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
    BoundedNumberInfo::writeToR1(os); os<<cc<<endl;
    os<<"initVals=";  wts::writeToR(os,initVals,wts::to_qcsv(ptrIB->getFwdIndexVector())); os<<cc<<endl;
    os<<"finalVals="; wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
    os<<")";
}

/**
 * Writes final values to an output stream as an R vector.
 * @param os - the output stream.
 */
void BoundedVectorInfo::writeFinalValsToR(ostream& os){
    if (debug) rpt::echo<<"BoundedVectorInfo::writeFinalValsToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
}

////////////////////////////BoundedVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  DevsVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************\n
*   Calculates devs.       \n
***************************************************************/
void DevsVectorInfo::calcDevs(void){
    initVals(N) = -sum(initVals(1,(N-1)));
}
/***************************************************************\n
*   Sets initial values.       \n
***************************************************************/
void DevsVectorInfo::setInitVals(dvector& x){
    if (debug) rpt::echo<<"starting DevsVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    BoundedVectorInfo::setInitVals(x);//use parent class
    calcDevs();//ensure devs
    if (debug) {
        rpt::echo<<"initVals: "<<initVals<<endl<<"vector x: "<<x<<endl;
        rpt::echo<<"finished DevsVectorInfo::setInitVals(dvector& x) for "<<name<<endl;
    }
}
/**
 * Sets initial values 1:(N-1) to those of the vector x, but
 * sets the value for element N to -sum(initVals(1,N-1)) so
 * the sum over all elements is 0. x may have size N-1.
 * 
 * @param x - param_init_bounded_vector of initial values
 */
void DevsVectorInfo::setInitVals(param_init_bounded_vector & x){
    if (debug) {
        rpt::echo<<"starting DevsVectorInfo::setInitVals(param_init_bounded_vector & x) for "<<name<<endl;
        rpt::echo<<"input x  index limits: "<<x.indexmin()<<cc<<x.indexmax()<<endl;
        rpt::echo<<"initVals index limits: "<<1<<cc<<N<<endl;
    }
    initVals(1,N-1) = value(x);
    calcDevs();
    if (debug) rpt::echo<<"finished DevsVectorInfo::setInitVals(param_init_bounded_vector & x) for "<<name<<endl;
}     

/**
 * Sets final values 1:(N-1) to those of the vector x, but
 * sets the value for element N to -sum(initVals(1,N-1)) so
 * the sum over all elements is 0. x may have size N-1.
 * 
 * @param x - param_init_bounded_vector of final values
 */
void DevsVectorInfo::setFinalVals(param_init_bounded_vector & x){
    if (debug) {
        rpt::echo<<"starting DevsVectorInfo::setFinalVals(param_init_bounded_vector & x) for "<<name<<endl;
        rpt::echo<<"input x  index limits: "<<x.indexmin()<<cc<<x.indexmax()<<endl;
        rpt::echo<<"initVals index limits: "<<1<<cc<<N<<endl;
    }
    if (!finlVals.allocated()) finlVals.allocate(initVals.indexmin(),initVals.indexmax());
    finlVals(1,N-1) = value(x);
    finlVals(N) = -sum(finlVals(1,N-1));
    if (debug) rpt::echo<<"finished DevsVectorInfo::setFinalVals(param_init_bounded_vector & x) for "<<name<<endl;
}     

/***************************************************************
*   read initial values. \n
***************************************************************/
void DevsVectorInfo::readInitVals(cifstream & is){
    BoundedVectorInfo::readInitVals(is);
    calcDevs();
}

/***************************************************************
*   Draw a random sample from the prior.                \n 
*   If phase<0, return initVals rather than resampling. \n
***************************************************************/
dvector DevsVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting DevsVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    dvector smpl = BoundedVectorInfo::drawInitVals(rng,vif);
    smpl(N) = -sum(smpl(1,(N-1)));
    if (debug) {
        rpt::echo<<phase<<tb<<pMPI->canSample()<<endl;
        rpt::echo<<"initVals: "<<initVals<<endl<<"samples: "<<smpl<<endl;
        rpt::echo<<"finished DevsVectorInfo::drawInitVals(random_number_generator& rng) for "<<name<<endl;
    }
    return smpl;
}

/***************************************************************
*  Read from cifstream object.\n
 * Read order is:
 * mni, mxi, readVals + BoundedNumberInfo::read(is)
***************************************************************/
void DevsVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"Starting DevsVectorInfo::read(cifstream & is) for "<<name<<endl;
    BoundedVectorInfo::read(is);
    initVals = 0.0;
    if (debug) {
        rpt::echo<<"idxType    = "<<idxType<<endl;
        rpt::echo<<"IndexBlock = "<<(*ptrIB)<<endl;
        rpt::echo<<"N          = "<<N<<endl;
        rpt::echo<<"Done DevsVectorInfo::read(cifstream & is) for "<<name<<endl;
    }
}

/***************************************************************
*   writeToR                                                    *
***************************************************************/
void DevsVectorInfo::writeToR(ostream& os){
    if (debug) rpt::echo<<"DevsVectorInfo::writeToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    os<<"list(";
    BoundedNumberInfo::writeToR1(os); os<<cc<<endl;
    os<<"initVals=";  wts::writeToR(os,initVals,wts::to_qcsv(ptrIB->getFwdIndexVector())); os<<cc<<endl;
    os<<"finalVals="; wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
    os<<")";
}

/**
 * Writes final values to an output stream as an R vector.
 * @param os - the output stream.
 */
void DevsVectorInfo::writeFinalValsToR(ostream& os){
    if (debug) rpt::echo<<"DevsVectorInfo::writeFinalValsToR for "<<this->name<<endl;
    if (!finlVals.allocated()) {
        finlVals.allocate(initVals.indexmin(),initVals.indexmax()); 
        finlVals = initVals;
    }
    wts::writeToR(os,finlVals,wts::to_qcsv(ptrIB->getFwdIndexVector()));
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

/***************************************************************
*   Get initial parameter values.                              *
***************************************************************/
dvector NumberVectorInfo::getInitVals(){
    if (debug) rpt::echo<<"starting NumberVectorInfo::getInitVals()"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = ppNIs[i-1]->getInitVal();
    if (debug) rpt::echo<<"finished NumberVectorInfo::getInitVals()"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Set initial parameter values.                              *
***************************************************************/
void NumberVectorInfo::setInitVals(param_init_number_vector& x){
    if (debug) rpt::echo<<"starting NumberVectorInfo::setInitVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) ppNIs[i-1]->setInitVal(x(i));
    if (debug) rpt::echo<<"finished NumberVectorInfo::setInitVals(x)"<<this<<endl;
}

/***************************************************************
*   Set final parameter values.                              *
***************************************************************/
void NumberVectorInfo::setFinalVals(param_init_number_vector& x){
    if (debug) rpt::echo<<"starting NumberVectorInfo::setFinalVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) ppNIs[i-1]->setFinalVal(x(i));
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

/***************************************************************
*   Calculate log prior probability for all parameters.        *
***************************************************************/
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
    os<<"#id init_val phase resample? prior_wgt prior_type prior_params prior_consts  label"<<endl;
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
////////////////////////////NumberVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  BoundedNumberVectorInfo
 *----------------------------------------------------------------------------*/
/**
 * Get a pointer to the ith BoundedNumberInfo element in the vector. 
 * i mustbe in the interval 1<=i<=nNIs.
 * @param i - index (1-based) to BoundedNumberInfo
 * @return pointer to ith BoundedNumberInfo object
 */
BoundedNumberInfo* BoundedNumberVectorInfo::operator[](int i){
    if (ppNIs&&(i<=nNIs)) return (static_cast<BoundedNumberInfo*>(ppNIs[i-1])); 
    return 0;
}
/***************************************************************
*   get lower bounds for parameters as vector                  *
***************************************************************/
dvector BoundedNumberVectorInfo::getLowerBounds(void){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector lbs(1,nNIs);
    lbs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) lbs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getLowerBound();
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::getLowerBounds(void) "<<this<<endl;
    return lbs;
}

/***************************************************************
*   get upper bounds for parameters as vector                  *
***************************************************************/
dvector BoundedNumberVectorInfo::getUpperBounds(void){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nNIs);
    ubs.initialize();
    if (ppNIs) for (int i=1;i<=nNIs;i++) ubs(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->getUpperBound();
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

/***************************************************************
*   Set initial parameter values.                              *
***************************************************************/
void BoundedNumberVectorInfo::setInitVals(param_init_bounded_number_vector& x){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::setInitVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->setInitVal(x(i));
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::setInitVals(x)"<<this<<endl;
}

/***************************************************************
*   Set final parameter values.                              *
***************************************************************/
void BoundedNumberVectorInfo::setFinalVals(param_init_bounded_number_vector& x){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::setFinalVals(x)"<<this<<endl;
    for (int i=1;i<=nNIs;i++) (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->setFinalVal(x(i));
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::setFinalVals(x)"<<this<<endl;
}

/***************************************************************
*   Draw initial parameter values.                             *
***************************************************************/
dvector BoundedNumberVectorInfo::drawInitVals(random_number_generator& rng, double vif){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::drawInitVals(random_number_generator& rng,vif)"<<this<<endl;
    dvector initVals(1,nNIs);
    for (int i=1;i<=nNIs;i++) initVals(i) = (static_cast<BoundedNumberInfo*>(ppNIs[i-1]))->drawInitVal(rng,vif);
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::drawInitVals(random_number_generator& rng,vif)"<<this<<endl;
    return initVals;
}

/***************************************************************
*   Read from stream.                                          *
***************************************************************/
void BoundedNumberVectorInfo::read(cifstream & is){
    if (debug) rpt::echo<<"starting BoundedNumberVectorInfo::read(cifstream & is) "<<name<<endl;
    is>>nNIs;
    if (debug) rpt::echo<<"nNIs ="<<tb<<nNIs<<endl;
    if (ppNIs) deallocate();
    if (nNIs>0) {
        int idx;
        ppNIs = new NumberInfo*[nNIs];
        for (int p=0;p<nNIs;p++) {
            is>>idx;
            if (idx<=nNIs){
                BoundedNumberInfo* pBNI = new BoundedNumberInfo();
                is>>(*pBNI);
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
            for (int p=0;p<nNIs;p++) rpt::echo<<p+1<<tb<<(*ppNIs[p])<<endl;
        }
    }
    if (debug) rpt::echo<<"finished BoundedNumberVectorInfo::read(cifstream & is) "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedNumberVectorInfo::write(ostream & os){
    if (debug) rpt::echo<<"Starting BoundedNumberVectorInfo::write(ostream & os) for "<<name<<endl;
    os<<tb<<nNIs<<"  #number of bounded parameters"<<endl;
    os<<"#id lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts  label"<<endl;
    if (nNIs){
        for (int p=0;p<(nNIs-1);p++) os<<(p+1)<<tb<<(*ppNIs[p])<<endl;
        os<<nNIs<<tb<<(*ppNIs[nNIs-1]);
    }
    if (debug) rpt::echo<<endl<<"Finished BoundedNumberVectorInfo::write(ostream & os)"<<endl;
}

/***************************************************************
*   Write parameters info to stream in R format.               *
***************************************************************/
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
/***************************************************************
*   deallocation                                               *
***************************************************************/
void VectorVectorInfo::deallocate(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::deallocate(void) for "<<name<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) rpt::echo<<"finished VectorVectorInfo::deallocate(void) for "<<name<<endl;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector VectorVectorInfo::getMinIndices(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getMinIndices(void) for "<<name<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) idxs=1;
    if (debug) rpt::echo<<"finished VectorVectorInfo::getMinIndices(void) for "<<name<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector VectorVectorInfo::getMaxIndices(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getMaxIndices(void) for "<<name<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getSize();
    if (debug) rpt::echo<<"finished VectorVectorInfo::getMaxIndices(void) for "<<name<<endl;
    return idxs;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector VectorVectorInfo::getPhases(void){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getPhases(void) for "<<name<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) rpt::echo<<"finished VectorVectorInfo::getPhases(void) for "<<name<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector VectorVectorInfo::getPriorWgts(){
    if (debug) rpt::echo<<"starting VectorVectorInfo::getPriorWgts() for "<<name<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) rpt::echo<<"finished VectorVectorInfo::getPriorWgts() for "<<name<<endl;
    return wts;
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
        if (debug) {
            for (int p=0;p<nVIs;p++) rpt::echo<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) rpt::echo<<"finished VectorVectorInfo::read(cifstream & is) for "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void VectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of parameters"<<endl;
    os<<"#id idx.type  idx.block  read? init_val phase resample? prior_wgt prior_type prior_params prior_consts  label"<<endl;
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
/***************************************************************
*   get lower bounds for parameters as vector                  *
***************************************************************/
dvector BoundedVectorVectorInfo::getLowerBounds(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector lbs(1,nVIs);
    lbs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) lbs(i) = (ppVIs[i-1])->getLowerBound();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    return lbs;
}

/***************************************************************
*   get upper bounds for parameters as vector                  *
***************************************************************/
dvector BoundedVectorVectorInfo::getUpperBounds(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nVIs);
    ubs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) ubs(i) = (ppVIs[i-1])->getUpperBound();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
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
        ppVIs = new BoundedVectorInfo*[nVIs];
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
            BoundedVectorInfo* pVI = ppVIs[p];
            if (pVI->readVals) {
                if (debug) rpt::echo<<"Reading vector values for "<<p+1<<"th bounded vector"<<endl;
                pVI->readInitVals(is);
                if (debug) rpt::echo<<pVI->getInitVals()<<endl;
            }
        }
        if (debug) {
            for (int p=0;p<nVIs;p++) rpt::echo<<p+1<<tb<<(*ppVIs[p])<<endl;
        }
    }
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::read(cifstream & is) "<<name<<endl;
}

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void BoundedVectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of bounded parameters"<<endl;
    os<<"#id   idx.block   read? lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts  label"<<endl;
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
/***************************************************************
*   deallocation                                               *
***************************************************************/
void BoundedVectorVectorInfo::deallocate(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::deallocate(void) "<<this<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::deallocate(void) "<<this<<endl;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector BoundedVectorVectorInfo::getMinIndices(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) idxs=1;
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector BoundedVectorVectorInfo::getMaxIndices(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getSize();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector BoundedVectorVectorInfo::getPhases(void){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector BoundedVectorVectorInfo::getPriorWgts(){
    if (debug) rpt::echo<<"starting BoundedVectorVectorInfo::getPriorWgts()"<<this<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) rpt::echo<<"finished BoundedVectorVectorInfo::getPriorWgts()"<<this<<endl;
    return wts;
}
////////////////////////////BoundedVectorVectorInfo/////////////////////////////////

/*------------------------------------------------------------------------------
 *  DevsVectorVectorInfo
 *----------------------------------------------------------------------------*/
/***************************************************************
*   get lower bounds for parameters as vector                  *
***************************************************************/
dvector DevsVectorVectorInfo::getLowerBounds(void){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    dvector lbs(1,nVIs);
    lbs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) lbs(i) = (ppVIs[i-1])->getLowerBound();
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::getLowerBounds(void) "<<this<<endl;
    return lbs;
}

/***************************************************************
*   get upper bounds for parameters as vector                  *
***************************************************************/
dvector DevsVectorVectorInfo::getUpperBounds(void){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    dvector ubs(1,nVIs);
    ubs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) ubs(i) = (ppVIs[i-1])->getUpperBound();
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::getUpperBounds(void) "<<this<<endl;
    return ubs;
}

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
        ppVIs = new DevsVectorInfo*[nVIs];
        for (int p=0;p<nVIs;p++) {
            is>>idx;
            if ((idx>0)&&(idx<=nVIs)){
                DevsVectorInfo* pDVI = new DevsVectorInfo();
                is>>(*pDVI);
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
            DevsVectorInfo* pVI = ppVIs[p];
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

/***************************************************************
*   Write to stream.                                           *
***************************************************************/
void DevsVectorVectorInfo::write(ostream & os){
    os<<tb<<nVIs<<"  #number of devs vectors"<<endl;
    os<<"#id   idx.block   read? lb ub jitter? init_val phase resample? prior_wgt prior_type prior_params prior_consts  label"<<endl;
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
void DevsVectorVectorInfo::writeToR(ostream& os, adstring nm, int indent){
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
void DevsVectorVectorInfo::writeFinalValsToR(ostream& os){
    if (nVIs){
        os<<"list("<<endl;
        for (int p=1;p<nVIs;p++) {os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<","<<endl;}
        int p=nVIs;               os<<tb<<"'"<<p<<"'="; ppVIs[p-1]->writeFinalValsToR(os); os<<")";
    } else {
        os<<"NULL";
    }
}
/***************************************************************
*   deallocation                                               *
***************************************************************/
void DevsVectorVectorInfo::deallocate(void){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::deallocate(void) "<<this<<endl;
    if (ppVIs) {
        for (int p=0;p<nVIs;p++) if (ppVIs[p]!=0) delete ppVIs[p];
        delete[] ppVIs;
        ppVIs = 0;
    }
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::deallocate(void) "<<this<<endl;
}

/***************************************************************
*   get min indices                                            *
***************************************************************/
ivector DevsVectorVectorInfo::getMinIndices(void){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) idxs=1;
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::getMinIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get max indices                                            *
***************************************************************/
ivector DevsVectorVectorInfo::getMaxIndices(void){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    ivector idxs(1,nVIs);
    idxs.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) idxs(i) = ppVIs[i-1]->getSize();
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::getMaxIndices(void) "<<this<<endl;
    return idxs;
}

/***************************************************************
*   get phase info                                             *
***************************************************************/
ivector DevsVectorVectorInfo::getPhases(void){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::getPhases(void) "<<this<<endl;
    ivector phases(1,nVIs);
    phases.initialize();
    if (ppVIs) for (int i=1;i<=nVIs;i++) phases(i) = ppVIs[i-1]->getPhase();
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::getPhases(void) "<<this<<endl;
    return phases;
}

/***************************************************************
*   get likelihood weights for log prior probabilities         *
***************************************************************/
dvector DevsVectorVectorInfo::getPriorWgts(){
    if (debug) rpt::echo<<"starting DevsVectorVectorInfo::getPriorWgts()"<<this<<endl;
    dvector wts(1,nVIs);
    for (int i=1;i<=nVIs;i++) wts(i) = ppVIs[i-1]->getPriorWgt();
    if (debug) rpt::echo<<"finished DevsVectorVectorInfo::getPriorWgts()"<<this<<endl;
    return wts;
}
////////////////////////////DevsVectorVectorInfo/////////////////////////////////
