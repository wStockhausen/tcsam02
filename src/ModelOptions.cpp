#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelOptions.hpp"

//**********************************************************************
//  Includes
//      EffAvgScenario
//      EffAvgScenarios
//      CapRateAvgScenario
//      CapRateAvgScenarios
//      EffXtrapScenarios
//      ProjectionScenarios
//      ModelOptions
//**********************************************************************
using namespace std;

//--------------------------------------------------------------------------------
//          EmpiricalSelFcn
//--------------------------------------------------------------------------------
int EmpiricalSelFcn::debug = 0;
EmpiricalSelFcn::EmpiricalSelFcn(ModelConfiguration& mc){
    ptrMC = &mc;
}

EmpiricalSelFcn::~EmpiricalSelFcn(){
    ptrMC = 0;
}
        
void EmpiricalSelFcn::read(cifstream& is){
    if (debug) cout<<"starting EmpiricalSelFcn::read(cifstream& is)"<<endl;
    zBs.allocate(1,ptrMC->nZBs);
    esf.allocate(1,ptrMC->nZBs);
    is>>id;
    is>>zBs;
    is>>esf;
    
    if (debug) {
        cout<<"finished EmpiricalSelFcn::read(cifstream& is)"<<endl;
    }
}

void EmpiricalSelFcn::write(std::ostream& os){
    if (debug) cout<<"starting EmpiricalSelFcn::write(ostream& os)"<<endl;
    os<<id<<endl;
    os<<zBs<<endl;
    os<<esf<<endl;
    if (debug) cout<<"finished EmpiricalSelFcn::write(ostream& os)"<<endl;
}

void EmpiricalSelFcn::writeToR(std::ostream& os){
    os<<"list(id="<<id<<cc;
    os<<"esf="; wts::writeToR(os,esf,ptrMC->dimZBsToR); os<<")";
}
//--------------------------------------------------------------------------------
//          EmpiricalSelFcns
//--------------------------------------------------------------------------------
int EmpiricalSelFcns::debug = 0;
EmpiricalSelFcns::EmpiricalSelFcns(ModelConfiguration& mc){
    ptrMC  = &mc;
    nESFs  = 0;
    ppESFs = 0;
}

EmpiricalSelFcns::~EmpiricalSelFcns(){
    ptrMC = 0;
    if (ppESFs) {
        for (int p=0;p<nESFs;p++) delete ppESFs[p];
        delete ppESFs; ppESFs = 0;
    }
}

void EmpiricalSelFcns::read(cifstream& is){
    is>>nESFs;
    if (ppESFs) {
        for (int p=0;p<nESFs;p++) delete ppESFs[p];
        delete ppESFs; ppESFs = 0;
    }
    ppESFs = new EmpiricalSelFcn*[nESFs];
    for (int p=0;p<nESFs;p++) {
        ppESFs[p] = new EmpiricalSelFcn(*ptrMC);
        ppESFs[p]->read(is);
    }
}

void EmpiricalSelFcns::write(std::ostream& os){
    os<<"#------Empirical Selectivity Functions"<<endl;
    os<<nESFs<<tb<<"#number of empirical selectivity functions"<<endl;
    os<<"# id  z's  values"<<endl;
    for (int p=0;p<nESFs;p++) {ppESFs[p]->write(os); os<<endl;}
}

void EmpiricalSelFcns::writeToR(std::ostream& os){
    os<<"list("<<"nESFs="<<nESFs<<cc<<endl;
    for (int n=1;n<nESFs;n++) {os<<"`"<<n<<"`="; ppESFs[n-1]->writeToR(os); os<<","<<endl;}
    os<<"`"<<nESFs<<"`="; ppESFs[nESFs-1]->writeToR(os); os<<endl;
    os<<")";
}
//--------------------------------------------------------------------------------
//          EmpiricalSelFcnPrior
//--------------------------------------------------------------------------------
int EmpiricalSelFcnPrior::debug = 0;
/**
 * Class instantiator
 * 
 * @param mc - pointer to model configuration object
 */
EmpiricalSelFcnPrior::EmpiricalSelFcnPrior(ModelConfiguration& mc){
    ptrMC = &mc;
    pMPI = 0;
}

/**
 * Class destructor.
 */
EmpiricalSelFcnPrior::~EmpiricalSelFcnPrior(){
    ptrMC = 0;
    if (pMPI) delete pMPI; pMPI=0;
}
        
/**
 * Sets the prior type, based on an adstring value.
 * 
 * @param prior - the prior type, as an adstring
 */
void EmpiricalSelFcnPrior::setPriorType(adstring & prior){
    if (debug) rpt::echo<<"starting NumberInfo::setPriorType(adstring prior)"<<this<<endl;
    if (pMPI) delete pMPI;
    priorType=prior;
    pMPI = ModelPDFInfo::getInfo(prior);
    if (pMPI) {
        if (pMPI->getNumParams()!=2) {
            cout<<"Error defining EmpiricalSelFcnPrior."<<endl;
            cout<<"Number of parameters for prior must be two, but "<<prior<<" requires "<<pMPI->getNumParams()<<endl;
            cout<<"Please specify another prior in the ModelOptions file."<<endl;
            exit(-1);
        }
        if (pMPI->getNumConsts()>0)  {
            cout<<"Error defining EmpiricalSelFcnPrior."<<endl;
            cout<<"Number of constants for prior must be two, but "<<prior<<" requires "<<pMPI->getNumConsts()<<endl;
            cout<<"Please specify another prior in the ModelOptions file."<<endl;
            exit(-1);
        }

    }
    if (debug) rpt::echo<<"finished NumberInfo::setPriorType(adstring prior)"<<this<<endl;
}

void EmpiricalSelFcnPrior::read(cifstream& is){
    if (debug) cout<<"starting EmpiricalSelFcnPrior::read(cifstream& is)"<<endl;
    is>>id;
    is>>sel_id;
    is>>priorWgt;
    is>>priorType;
    setPriorType(priorType);
    zBs.allocate(1,ptrMC->nZBs); zBs.initialize();
    is>>zBs;
    p1.allocate(1,ptrMC->nZBs);  p1.initialize();
    is>>p1;
    p2.allocate(1,ptrMC->nZBs);  p2.initialize();
    is>>p2;
    
    if (debug) {
        cout<<"finished EmpiricalSelFcnPrior::read(cifstream& is)"<<endl;
    }
}

void EmpiricalSelFcnPrior::write(std::ostream& os){
    if (debug) cout<<"starting EmpiricalSelFcnPrior::write(ostream& os)"<<endl;
    os<<id       <<tb<<"#  id for selectivity function prior"<<endl;
    os<<sel_id   <<tb<<"#  selectivity function id as defined in Model Parameters Info"<<endl;
    os<<priorWgt <<tb<<"#  multiplicative weight to apply to prior"<<endl;
    os<<priorType<<tb<<"#name of prior to apply"<<endl;
    os<<zBs      <<tb<<"#  sizes at which to prior function can be evaluated"<<endl;
    os<<p1       <<tb<<"#  1st parameter values at which prior can be evaluated"<<endl;
    os<<p2       <<tb<<"#  2nd parameter values at which prior can be evaluated"<<endl;
    if (debug) cout<<"finished EmpiricalSelFcnPrior::write(ostream& os)"<<endl;
}

void EmpiricalSelFcnPrior::writeToR(std::ostream& os){
    os<<"list(id="<<id<<cc<<"sel_id="<<sel_id<<cc<<"priorWgt="<<priorWgt<<cc<<"priorType='"<<priorType<<"',"<<endl;
    os<<"p1="; wts::writeToR(os,p1,ptrMC->dimZBsToR); os<<cc<<endl;
    os<<"p2="; wts::writeToR(os,p2,ptrMC->dimZBsToR); os<<")";
}
//--------------------------------------------------------------------------------
//          EmpiricalSelFcnPriors
//--------------------------------------------------------------------------------
int EmpiricalSelFcnPriors::debug = 0;
EmpiricalSelFcnPriors::EmpiricalSelFcnPriors(ModelConfiguration& mc){
    ptrMC  = &mc;
    nESPs  = 0;
    ppESPs = 0;
}

EmpiricalSelFcnPriors::~EmpiricalSelFcnPriors(){
    ptrMC = 0;
    if (ppESPs) {
        for (int p=0;p<nESPs;p++) delete ppESPs[p];
        delete ppESPs; ppESPs = 0;
    }
    nESPs=0;
}

void EmpiricalSelFcnPriors::read(cifstream& is){
    is>>nESPs;
    if (ppESPs) {
        for (int p=0;p<nESPs;p++) delete ppESPs[p];
        delete ppESPs; ppESPs = 0;
    }
    ppESPs = new EmpiricalSelFcnPrior*[nESPs];
    for (int p=0;p<nESPs;p++) {
        ppESPs[p] = new EmpiricalSelFcnPrior(*ptrMC);
        ppESPs[p]->read(is);
    }
}

void EmpiricalSelFcnPriors::write(std::ostream& os){
    os<<"#------Empirical Selectivity Functions"<<endl;
    os<<nESPs<<tb<<"#number of empirical selectivity functions"<<endl;
    os<<"# id  z's  values"<<endl;
    for (int p=0;p<nESPs;p++) {ppESPs[p]->write(os); os<<endl;}
}

void EmpiricalSelFcnPriors::writeToR(std::ostream& os){
    os<<"list("<<"nESPs="<<nESPs<<cc<<endl;
    for (int n=1;n<nESPs;n++) {os<<"`"<<n<<"`="; ppESPs[n-1]->writeToR(os); os<<","<<endl;}
    os<<"`"<<nESPs<<"`="; ppESPs[nESPs-1]->writeToR(os); os<<endl;
    os<<")";
}
//--------------------------------------------------------------------------------
//          EffAvgScenario
//--------------------------------------------------------------------------------
int EffAvgScenario::debug = 0;
EffAvgScenario::EffAvgScenario(ModelConfiguration& mc){
    ptrMC = &mc;
    ptrIB = new IndexBlock(ptrMC->mnYr,ptrMC->mxYr);
    strVals.allocate(1,2);
}

EffAvgScenario::~EffAvgScenario(){
    if (ptrIB) {delete ptrIB; ptrIB=0;}
    ptrMC = 0;
}
        
/**
 * Sets the max year for the averaging period (for retrospective analysis).
 * 
 * @param mxYr - max year for averaging 
 */
void EffAvgScenario::setMaxYearForAveraging(int mxYr){
    if (ptrIB) {
        ptrIB->checkMaxAndReset(mxYr);
        rpt::echo<<"years = "<<ptrIB->asString()<<endl;
    }
}

void EffAvgScenario::read(cifstream& is){
    if (debug) cout<<"starting EffAvgScenario::read(cifstream& is)"<<endl;
    is>>id;
    for (int i=strVals.indexmin();i<=strVals.indexmax(); i++) is>>strVals(i);
    
    int k = 1;
    f = wts::which(strVals(k++),ptrMC->lblsFsh);
    ptrIB->parse(strVals(k++));
    
    if (debug) {
        int k = 1;
        cout<<"f = "<<f<<tb<<strVals(k++)<<endl;
        cout<<"years = "<<ptrIB->asString()<<endl;
        cout<<"finished EffAvgScenario::read(cifstream& is)"<<endl;
    }
}

void EffAvgScenario::write(std::ostream& os){
    if (debug) cout<<"starting EffAvgScenario::write(ostream& os)"<<endl;
    os<<id<<tb;
    for (int i=strVals.indexmin();i<=strVals.indexmax(); i++) os<<strVals(i)<<tb;
    if (debug) cout<<"finished EffAvgScenario::write(ostream& os)"<<endl;
}

void EffAvgScenario::writeToR(std::ostream& os){
    os<<"list(f="<<f<<cc<<"years="<<getYDimsForR()<<")";
}
//--------------------------------------------------------------------------------
//          EffAvgOptions
//--------------------------------------------------------------------------------
int EffAvgScenarios::debug = 0;
EffAvgScenarios::EffAvgScenarios(ModelConfiguration& mc){
    ptrMC = &mc;
    nAvgs = 0;
    ppEASs = 0;
}

EffAvgScenarios::~EffAvgScenarios(){
    ptrMC = 0;
    if (ppEASs) {
        for (int p=0;p<nAvgs;p++) delete ppEASs[p];
        delete ppEASs; ppEASs = 0;
    }
}

/**
 * Sets the max year for all averaging periods (for retrospective analysis).
 * 
 * @param mxYr - max year for averaging 
 */
void EffAvgScenarios::setMaxYearForAveraging(int mxYr){
    if (ppEASs) {
        for (int p=0;p<nAvgs;p++) ppEASs[p]->setMaxYearForAveraging(mxYr);
    }    
}

void EffAvgScenarios::read(cifstream& is){
    is>>nAvgs;
    if (ppEASs) {
        for (int p=0;p<nAvgs;p++) delete ppEASs[p];
        delete ppEASs; ppEASs = 0;
    }
    ppEASs = new EffAvgScenario*[nAvgs];
    for (int p=0;p<nAvgs;p++) {
        ppEASs[p] = new EffAvgScenario(*ptrMC);
        ppEASs[p]->read(is);
    }
}

void EffAvgScenarios::write(std::ostream& os){
    os<<"#----effort averaging scenarios"<<endl;
    os<<nAvgs<<tb<<"#number of effort averaging scenarios"<<endl;
    os<<"#id  fishery  year_block"<<endl;
    for (int p=0;p<nAvgs;p++) {ppEASs[p]->write(os); os<<endl;}
}

void EffAvgScenarios::writeToR(std::ostream& os){
    os<<"list("<<"nAvgs="<<nAvgs<<cc<<endl;
    for (int n=1;n<nAvgs;n++) {os<<"`"<<n<<"`="; ppEASs[n-1]->writeToR(os); os<<","<<endl;}
    os<<"`"<<nAvgs<<"`="; ppEASs[nAvgs-1]->writeToR(os); os<<endl;
    os<<")";
}
//--------------------------------------------------------------------------------
//          CapRateAvgScenario
//--------------------------------------------------------------------------------
int CapRateAvgScenario::debug = 0;
CapRateAvgScenario::CapRateAvgScenario(ModelConfiguration& mc){
    ptrMC = &mc;
    strVals.allocate(1,4);
}

CapRateAvgScenario::~CapRateAvgScenario(){ptrMC = 0;}

///**
// * Sets the max year for the averaging period (for retrospective analysis).
// * 
// * @param mxYr - max year for averaging 
// */
//void CapRateAvgScenario::setMaxYearForAveraging(int mxYr){
//    
//}
//
void CapRateAvgScenario::read(cifstream& is){
    if (debug) cout<<"starting CapRateAvgScenario::read(cifstream& is)"<<endl;
    is>>id;
    for (int i=strVals.indexmin();i<=strVals.indexmax(); i++) is>>strVals(i);
    is>>idParam;
    is>>idEffAvgInfo;
    is>>optAvg;
    is>>llWgt;
    
    int k = 1;
    f = wts::which(strVals(k++),ptrMC->lblsFsh);
    x = tcsam::getSexType(strVals(k++));
    m = tcsam::getMaturityType(strVals(k++));
    s = tcsam::getShellType(strVals(k++));
    
    if (debug) {
        int k = 1;
        cout<<"f = "<<f<<tb<<strVals(k++)<<endl;
        cout<<"x = "<<x<<tb<<strVals(k++)<<endl;
        cout<<"m = "<<m<<tb<<strVals(k++)<<endl;
        cout<<"s = "<<s<<tb<<strVals(k++)<<endl;
        cout<<"finished CapRateAvgScenario::read(cifstream& is)"<<endl;
    }
}

void CapRateAvgScenario::write(std::ostream& os){
    if (debug) cout<<"starting CapRateAvgScenario::write(ostream& os)"<<endl;
    os<<id<<tb;
    for (int i=strVals.indexmin();i<=strVals.indexmax(); i++) os<<strVals(i)<<tb;
    os<<idParam<<tb;
    os<<idEffAvgInfo<<tb;
    os<<optAvg<<tb;
    os<<llWgt;
    if (debug) cout<<"finished CapRateAvgScenario::write(ostream& os)"<<endl;
}

void CapRateAvgScenario::writeToR(std::ostream& os){
    os<<"list(f="<<f<<cc;
    os<<"x='"<<tcsam::getSexType(x)<<"', ";
    os<<"m='"<<tcsam::getSexType(m)<<"', ";
    os<<"s='"<<tcsam::getSexType(s)<<"', ";
    os<<"paramID="<<idParam<<",effAvgID="<<idEffAvgInfo<<",llWgt="<<llWgt<<",optAvg="<<optAvg;
    os<<")";
}
//--------------------------------------------------------------------------------
//          CapRateAvgScenarios
//--------------------------------------------------------------------------------
int CapRateAvgScenarios::debug = 1;
CapRateAvgScenarios::CapRateAvgScenarios(ModelConfiguration& mc){
    ptrMC = &mc;
    nAvgs = 0;
    ppCRASs = 0;
}

CapRateAvgScenarios::~CapRateAvgScenarios(){
    ptrMC = 0;
    if (ppCRASs) {
        for (int p=0;p<nAvgs;p++) delete ppCRASs[p];
        delete ppCRASs; ppCRASs = 0;
    }
}

///**
// * Sets the max year for all averaging periods (for retrospective analysis).
// * 
// * @param mxYr - max year for averaging 
// */
//void CapRateAvgScenarios::setMaxYearForAveraging(int mxYr){
//    if (ppCRASs) {
//        for (int p=0;p<nAvgs;p++) ppCRASs[p]->setMaxYearForAveraging(mxYr);
//    }    
//}
//
void CapRateAvgScenarios::read(cifstream& is){
    is>>nAvgs;
    if (ppCRASs) {
        for (int p=0;p<nAvgs;p++) delete ppCRASs[p];
        delete ppCRASs; ppCRASs = 0;
    }
    ppCRASs = new CapRateAvgScenario*[nAvgs];
    for (int p=0;p<nAvgs;p++) {
        ppCRASs[p] = new CapRateAvgScenario(*ptrMC);
        ppCRASs[p]->read(is);
    }
}

void CapRateAvgScenarios::write(std::ostream& os){
    os<<"#----fishery capture rate averaging scenarios"<<endl;
    os<<"#-----capture rate averaging options:"<<endl;
    os<<"# 1 - average fully-selected capture rate"<<endl;
    os<<"# 2 - average mean size-specific capture rate"<<endl;
    os<<nAvgs<<tb<<"#number of capture rate averaging scenarios"<<endl;
    os<<"#id  fishery  sex   maturity  shell   param_id  effort_averaging_id   averaging_option   llWgt"<<endl;
    for (int p=0;p<nAvgs;p++) {ppCRASs[p]->write(os); os<<endl;}
}

void CapRateAvgScenarios::writeToR(std::ostream& os){
    os<<"list("<<"nAvgs="<<nAvgs<<cc<<endl;
    for (int n=1;n<nAvgs;n++) {os<<"`"<<n<<"`="; ppCRASs[n-1]->writeToR(os); os<<","<<endl;}
    os<<"`"<<nAvgs<<"`="; ppCRASs[nAvgs-1]->writeToR(os); os<<endl;
    os<<")";
}
//--------------------------------------------------------------------------------
//          EffXtrapOptions
//--------------------------------------------------------------------------------
int EffXtrapScenarios::debug = 1;
EffXtrapScenarios::EffXtrapScenarios(ModelConfiguration& mc){
    ptrMC = &mc;
    ptrEffAvgScenarios = new EffAvgScenarios(*ptrMC);
    ptrCapRateAvgScenarios = new CapRateAvgScenarios(*ptrMC);
}

EffXtrapScenarios::~EffXtrapScenarios(){
    ptrMC = 0;
    if (ptrEffAvgScenarios) {delete ptrEffAvgScenarios; ptrEffAvgScenarios=0;}
    if (ptrCapRateAvgScenarios) {delete ptrCapRateAvgScenarios; ptrCapRateAvgScenarios = 0;}
}
        
/**
 * Sets the max year for all averaging periods (for retrospective analysis).
 * 
 * @param mxYr - max year for averaging 
 */
void EffXtrapScenarios::setMaxYearForAveraging(int mxYr){
    if (ptrEffAvgScenarios)     ptrEffAvgScenarios->setMaxYearForAveraging(mxYr);
//    if (ptrCapRateAvgScenarios) ptrCapRateAvgScenarios->setMaxYearForAveraging(mxYr);
}

ivector EffXtrapScenarios::getTimePeriodForCapRateAveraging(int id){
    CapRateAvgScenario* ptrCRAS = ptrCapRateAvgScenarios->ppCRASs[id-1];
    int idEAS = ptrCRAS->idEffAvgInfo;
    return ptrEffAvgScenarios->ppEASs[idEAS-1]->ptrIB->getFwdIndexVector();
}

void EffXtrapScenarios::read(cifstream& is){
    is>>(*ptrEffAvgScenarios);
    is>>(*ptrCapRateAvgScenarios);
}

void EffXtrapScenarios::write(std::ostream& os){
    os<<"#----Effort Extrapolation Scenarios"<<endl;
    os<<(*ptrEffAvgScenarios); os<<endl;
    os<<(*ptrCapRateAvgScenarios); os<<endl;
}

void EffXtrapScenarios::writeToR(std::ostream& os){
    os<<"list("<<endl;
    os<<"effAvgOpts="; ptrEffAvgScenarios->writeToR(os); os<<","<<endl;
    os<<"capRateAvgOpts="; ptrCapRateAvgScenarios->writeToR(os); os<<endl;
    os<<")"<<endl;
}

//--------------------------------------------------------------------------------
//          ProjectionOptions
//--------------------------------------------------------------------------------
int ProjectionOptions::debug = 1;
ProjectionOptions::ProjectionOptions(){};
void ProjectionOptions::read(cifstream& is){
    is>>nReps;
    is>>nYrs;
    is>>nFs;
    if (nFs>0) {
        Fs.allocate(1,nFs);
        is>>Fs;
    }
    is>>nFMs;
    if (nFMs>0) {
        FMs.allocate(1,nFMs);
        is>>FMs;
    }
}
void ProjectionOptions::write(std::ostream& os){
    os<<nReps<<tb<<"#--number of repetitions"<<endl;
    os<<nYrs <<tb<<"#--number of years to project"<<endl;
    os<<nFs  <<tb<<"#--number of F's"<<endl;
    os<<"#--fishing mortality rates"<<endl;
    if (nFs>0) os<<Fs<<endl;
    os<<nFMs <<tb<<"#--number of Fofl multipliers"<<endl;
    os<<"#--fishing mortality multipliers"<<endl;
    if (nFMs>0) os<<FMs<<endl;
}
void ProjectionOptions::writeToR(std::ostream& os){
    os<<"list(nReps="<<nReps<<cc<<"nYrs="<<nYrs<<cc<<endl;
    os<<tb<<"nFs ="<<nFs <<cc<<"Fs ="; wts::writeToR(os,Fs);  os<<cc<<endl;
    os<<tb<<"nFMs="<<nFMs<<cc<<"FMs="; wts::writeToR(os,FMs); os<<")";
}
//--------------------------------------------------------------------------------
//          CatchDataSimOptions
//--------------------------------------------------------------------------------
int CatchDataSimOptions::debug = 0;
CatchDataSimOptions::CatchDataSimOptions(){
    if (debug) {
        rpt::echo<<"in CatchDataSimOptions constructor"<<endl;
        rpt::echo<<"name = '"<<name<<"'"<<endl;
    }
}
CatchDataSimOptions::~CatchDataSimOptions(){
    if (debug) {
        rpt::echo<<"in CatchDataSimOptions destructor"<<endl;
        rpt::echo<<"name = '"<<name<<"'"<<endl;
    }
}
void CatchDataSimOptions::read(cifstream & is){
    if (debug) {
        rpt::echo<<"Starting CatchDataSimOptions::read()"<<endl;
        rpt::echo<<name<<"  is the name before"<<endl;
        rpt::echo<<is.get_file_name()<<endl;
    }
    is>>name;      if (debug) rpt::echo<<name<<tb;
    is>>rngSeed;   if (debug) rpt::echo<<rngSeed<<tb;
    is>>expFacAbd; if (debug) rpt::echo<<expFacAbd<<tb;
    is>>expFacBio; if (debug) rpt::echo<<expFacBio<<tb;
    is>>expFacZCs; if (debug) rpt::echo<<expFacZCs<<endl;
    if (debug) rpt::echo<<"Finished CatchDataSimOptions::read()"<<endl;
}
void CatchDataSimOptions::write(ostream & os){
    if (debug) rpt::echo<<"starting CatchDataSimOptions::write() "<<endl;
    os<<name<<tb<<rngSeed<<tb<<expFacAbd<<tb<<expFacBio<<tb<<expFacZCs;
}
void CatchDataSimOptions::writeToR(ostream & os){
    if (debug) rpt::echo<<"starting CatchDataSimOptions::writeToR() "<<endl;
    os<<"list(name="<<name<<cc<<"rngSeed="<<rngSeed<<"expFacAbd="<<expFacAbd<<cc<<
            "expFacBio="<<expFacBio<<cc<<"expFacZCs="<<expFacZCs<<")";
}
//--------------------------------------------------------------------------------
//          SimOptions
//--------------------------------------------------------------------------------
int SimOptions::debug = 0;
SimOptions::SimOptions(){
    if (debug) rpt::echo<<"starting SimOptions::SimOptions "<<endl;
    ppRetCatch = 0;
    ppTotCatch = 0;
    ppDscCatch = 0;
    ppIdxCatch = 0;
    if (debug) rpt::echo<<"finished SimOptions::SimOptions "<<endl;
}
SimOptions::~SimOptions(){
    if (debug) rpt::echo<<"starting SimOptions::~SimOptions "<<endl;
    if (ppRetCatch) {                                                        
        for (int f=0;f<nRetCatch;f++) delete ppRetCatch[f];
        delete ppRetCatch; ppRetCatch = 0;
    } 
    if (ppTotCatch) {                                                        
        for (int f=0;f<nTotCatch;f++) delete ppTotCatch[f];
        delete ppTotCatch; ppTotCatch = 0;
    } 
    if (ppDscCatch) {                                                        
        for (int f=0;f<nDscCatch;f++) delete ppDscCatch[f];
        delete ppDscCatch; ppDscCatch = 0;
    } 
    if (ppIdxCatch) {                                                        
        for (int f=0;f<nIdxCatch;f++) delete ppIdxCatch[f];
        delete ppIdxCatch; ppIdxCatch = 0;
    } 
    if (debug) rpt::echo<<"finished SimOptions::~SimOptions "<<endl;
}
void SimOptions::read(cifstream& is){
    if (debug) rpt::echo<<"#--starting SimOptions::read()"<<std::endl;
    is>>nRetCatch; if (debug) rpt::echo<<"nRetCatch = "<<nRetCatch<<endl;
    if (nRetCatch){
        if (debug) rpt::echo<<"#-------Retained Catch SimOptions---------"<<std::endl;
        if (debug) rpt::echo<<"#--Creating retained catch SimOptions pointer array"<<std::endl;        
        ppRetCatch = new CatchDataSimOptions*[nRetCatch];
        if (debug) rpt::echo<<"#--Created retained CatchDataSimOptions pointer array"<<std::endl;        
        for (int i=0;i<nRetCatch;i++) {
            ppRetCatch[i] = new CatchDataSimOptions();
            if (debug) rpt::echo<<"Attempting to read "<<i<<"th retained CatchDataSimOption object "<<ppRetCatch[i]<<endl;
            is>>*(ppRetCatch[i]);
        }
    }
    is>>nTotCatch; if (debug) rpt::echo<<"nTotCatch = "<<nTotCatch<<endl;
    if (nTotCatch){
        if (debug) rpt::echo<<"#-------Total Catch SimOptions---------"<<std::endl;
        ppTotCatch = new CatchDataSimOptions*[nTotCatch];
        for (int i=0;i<nTotCatch;i++) {
            ppTotCatch[i] = new CatchDataSimOptions();
            is>>(*ppTotCatch[i]);
        }
    }
    is>>nDscCatch; if (debug) rpt::echo<<"nDscCatch = "<<nDscCatch<<endl;
    if (nDscCatch){
        if (debug) rpt::echo<<"#-------Discard Catch SimOptions---------"<<std::endl;
        ppDscCatch = new CatchDataSimOptions*[nDscCatch];
        for (int i=0;i<nDscCatch;i++) {
            ppDscCatch[i] = new CatchDataSimOptions();
            is>>(*ppDscCatch[i]);
        }
    }
    is>>nIdxCatch; if (debug) rpt::echo<<"nIdxCatch = "<<nIdxCatch<<endl;
    if (nIdxCatch){
        if (debug) rpt::echo<<"#-------Index Catch SimOptions---------"<<std::endl;
        ppIdxCatch = new CatchDataSimOptions*[nIdxCatch];
        for (int i=0;i<nIdxCatch;i++) {
            ppIdxCatch[i] = new CatchDataSimOptions();
            is>>(*ppIdxCatch[i]);
        }
    }
    is>>grwRngSeed;  if (debug) rpt::echo<<"grwRngSeed = "<<grwRngSeed<<endl;
    is>>grwMultFac;  if (debug) rpt::echo<<"grwMultFac = "<<grwMultFac<<endl;
    is>>modRngSeed;  if (debug) rpt::echo<<"modRngSeed = "<<modRngSeed<<endl;
    is>>modDivFac;   if (debug) rpt::echo<<"modDivFac  = "<<modDivFac<<endl;
}
void SimOptions::write(std::ostream& os){
    if (debug) rpt::echo<<"#--starting SimOptions::write()"<<std::endl;
    os<<"#--Simulation Options"<<endl;
    os<<"#--RNG seed: 0 = use default"<<endl;
    os<<"#--expansion factors:"<<endl;
    os<<"#---- =0: no observation noise added"<<endl;
    os<<"#---- >0: multiply CV"<<endl;
    os<<"#---- <0: replace CV"<<endl;
    os<<"#--retained catch fishery data"<<endl;
    os<<nRetCatch<<tb<<"#--number of fleets";
    if (nRetCatch){
        if (debug) rpt::echo<<"#-------Retained Catch SimOptions---------"<<std::endl;
        for (int i=0;i<nRetCatch;i++) os<<(*ppRetCatch[i])<<endl;
    }
    os<<"#--total catch fishery data"<<endl;
    os<<nTotCatch<<tb<<"#--number of fleets";
    if (nTotCatch){
        if (debug) rpt::echo<<"#-------Total Catch SimOptions---------"<<std::endl;
        for (int i=0;i<nTotCatch;i++) os<<(*ppTotCatch[i])<<endl;
    }
    os<<"#--discard catch fishery data"<<endl;
    os<<nDscCatch<<tb<<"#--number of fleets";
    if (nDscCatch){
        if (debug) rpt::echo<<"#-------Discard Catch SimOptions---------"<<std::endl;
        for (int i=0;i<nDscCatch;i++) os<<(*ppDscCatch[i])<<endl;
    }
    os<<"#--index catch data"<<endl;
    os<<nIdxCatch<<tb<<"#--number of fleets";
    if (nIdxCatch){
        if (debug) rpt::echo<<"#-------Index Catch SimOptions---------"<<std::endl;
        for (int i=0;i<nIdxCatch;i++) os<<(*ppIdxCatch[i])<<endl;
    }
    os<<"#----Growth data"<<endl;
    os<<grwRngSeed<<tb<<"#--RNG seed (0: don't reset seed)"<<endl;
    os<<grwMultFac   <<tb<<"#--observation noise scale factor (CV multiplier)"<<endl;
    os<<"#----Maturity ogive data"<<endl;
    os<<modRngSeed<<tb<<"#--RNG seed (0: don't reset seed)"<<endl;
    os<<modDivFac   <<tb<<"#--observation noise scale factor (reduces input sample sizes)"<<endl;
    if (debug) rpt::echo<<"#--finished SimOptions::write()"<<std::endl;
}
void SimOptions::writeToR(std::ostream& os){
    if (debug) rpt::echo<<"#--starting SimOptions::writeToR()"<<std::endl;
    os<<"list("<<endl;
    os<<tb<<"ret=list("<<endl;
    if (nRetCatch){
        for (int i=0;i<(nRetCatch-1);i++) {
            adstring name = ppRetCatch[i]->name;
            os<<tb<<tb<<name<<"="; ppRetCatch[i]->writeToR(os); os<<endl;
        }
        adstring name = ppRetCatch[nRetCatch-1]->name;
        os<<tb<<tb<<name<<"="; ppRetCatch[nRetCatch-1]->writeToR(os); os<<endl;
        os<<"),"<<endl;
    }
    //
    os<<tb<<"tot=list("<<endl;
    if (nTotCatch){
        for (int i=0;i<(nTotCatch-1);i++) {
            adstring name = ppTotCatch[i]->name;
            os<<tb<<tb<<name<<"="; ppTotCatch[i]->writeToR(os); os<<endl;
        }
        adstring name = ppTotCatch[nTotCatch-1]->name;
        os<<tb<<tb<<name<<"="; ppTotCatch[nTotCatch-1]->writeToR(os); os<<endl;
    }
    os<<"),"<<endl;
    //
    os<<tb<<"dsc=list("<<endl;
    if (nDscCatch){
        for (int i=0;i<(nDscCatch-1);i++) {
            adstring name = ppDscCatch[i]->name;
            os<<tb<<tb<<name<<"="; ppDscCatch[i]->writeToR(os); os<<endl;
        }
        adstring name = ppDscCatch[nDscCatch-1]->name;
        os<<tb<<tb<<name<<"="; ppDscCatch[nDscCatch-1]->writeToR(os); os<<endl;
    }
    os<<"),"<<endl;
    //
    os<<tb<<"idx=list("<<endl;
    if (nIdxCatch){
        for (int i=0;i<(nIdxCatch-1);i++) {
            adstring name = ppIdxCatch[i]->name;
            os<<tb<<tb<<name<<"="; ppIdxCatch[i]->writeToR(os); os<<endl;
        }
        adstring name = ppIdxCatch[nIdxCatch-1]->name;
        os<<tb<<tb<<name<<"="; ppIdxCatch[nIdxCatch-1]->writeToR(os); os<<endl;
    }
    os<<"),"<<endl;
    //
    os<<tb<<"grw=list(rngSeed ="<<grwRngSeed<<cc<<"multFac ="<<grwMultFac<<")"<<cc<<endl;
    os<<tb<<"mod=list(rngSeed ="<<modRngSeed<<cc<<"divFac ="<<modDivFac<<")"<<endl;
    os<<")";
    if (debug) rpt::echo<<"#--finished SimOptions::writeToR()"<<std::endl;
}
//--------------------------------------------------------------------------------
//          ModelOptions
//--------------------------------------------------------------------------------
int ModelOptions::debug = 0;
const adstring ModelOptions::VERSION = "2022.08.05";

ModelOptions::ModelOptions(ModelConfiguration& mc){
    ptrMC=&mc;
    
    //labels for initial n-at-z option 1
    opts1_InitNatZ.allocate(0,3);
    opts1_InitNatZ(0) = "#0 - build up n-at-z from recruitments (like TCSAM2013)"; 
    opts1_InitNatZ(1) = "#1 - calculate initial n-at-z using equilibrium calculations (like Gmacs)";
    opts1_InitNatZ(2) = "#2 - estimate initial n-at-z parameters on log scale      ('UNSCALED' version: pLnBaseInitN is ln-scale abundance in 1st population size bin)";
    opts1_InitNatZ(3) = "#3 - estimate initial n-at-z parameters on logistic scale (  'SCALED' version: pLnBaseInitN is ln-scale total abundance in 1st population size bin)";
    
    //labels for initial n-at-z option 2
    opts2_InitNatZ.allocate(1,2);
    opts2_InitNatZ(1) = "#1 - initialize using input values from MPI or pin file"; 
    opts2_InitNatZ(2) = "#2 - initialize using equilibrium calculations, ignore pin file values";
    
    //labels for initial n-at-z option 3
    opts3_InitNatZ.allocate(1,2);
    opts3_InitNatZ(1) = "#1 - normalize relative abundances to sum to 1         (guarantees sum-to-1)"; 
    opts3_InitNatZ(2) = "#2 - impose sum-to-1 constraint on relative abundances (does not guarantee sum-to-1)";
    
    //options for natural mortality parameterization
    optsParamNM.allocate(0,1);
    optsParamNM(0) = "use log-scale parameterization (default)";
    optsParamNM(1) = "use TCSAM2013 parameterization (arithmetic scale)"; 
    
    //growth parameterization options
    optsGrowthParam.allocate(0,2);
    optsGrowthParam(0) = "TCSAM2013 parameterization (ln-scale intercept, slope)"; 
    optsGrowthParam(1) = "parameterization based on min, max pre-molt sizes";
    optsGrowthParam(2) = "parameterization based on min pre-molt size, ln-scale slope";
    
    //growth pdf options
    optsGrowthPDF.allocate(0,1);
    optsGrowthPDF(0) = "use gamma probability distribution (like TCSAM2013)"; 
    optsGrowthPDF(1) = "use cumulative gamma distribution (like Gmacs)";
    
    //penalty options for prM2M parameters/ogives smoothness
    optsPenSmthPrM2M.allocate(0,1);
    optsPenSmthPrM2M(0) = "evaluate smoothness using parameters";
    optsPenSmthPrM2M(1) = "evaluate smoothness using ogives";    
    //penalty options for prM2M parameters/ogives being non-decreasing w/ size
    optsPenNonDecPrM2M.allocate(0,3);
    optsPenNonDecPrM2M(0) = "use posfun function on parameters";
    optsPenNonDecPrM2M(1) = "use exponential function on parameters";    
    optsPenNonDecPrM2M(2) = "use posfun function on ogives";
    optsPenNonDecPrM2M(3) = "use exponential function on ogives";    
    
    //smoothness penalty options for nonparametric selectivity function parameters/curves
    optsPenSmthNPSel.allocate(0,1);
    optsPenSmthNPSel(0) = "evaluate smoothness using selectivity parameters";
    optsPenSmthNPSel(1) = "evaluate smoothness using selectivity functions";    
    
    //empirical selectivity function options
    ptrEmpiricalSelFcns = new EmpiricalSelFcns(mc);
    
    //empirical selectivity function prior options
    ptrEmpiricalSelFcnPriors = new EmpiricalSelFcnPriors(mc);
    
    //effort extrapolation options
    ptrEffXtrapScenarios = new EffXtrapScenarios(mc);
    
    //options for OFL calculations: capture rate/selectivity function averaging
    optsOFLAvgCapRate.allocate(0,1);
    optsOFLAvgCapRate(0) = "average max capture rates, selectivity functions (like TCSAM2013)";
    optsOFLAvgCapRate(1) = "average size-specific capture rates";
    
    //options for iterative re-weighting
    optsIterativeReweighting.allocate(0,2);
    optsIterativeReweighting(0) = "no iterative re-weighting";
    optsIterativeReweighting(1) = "use harmonic means of McAllister-Ianelli effective N's";
    optsIterativeReweighting(2) = "use Francis weights";
    
    //options for projections
    ptrProjOpts     = new ProjectionOptions();
    ptrProjOptsMCMC = new ProjectionOptions();
    
    //options for simulations
    ptrSimOpts = new SimOptions();
}
/**
 * Read from input stream in ADMB format.
 * 
 * @param is - input stream
 */
void ModelOptions::read(cifstream & is) {
    cout<<"ModelOptions::read(cifstream & is)"<<endl;
    int idx;
    adstring str;
    is>>str;
    if (str!=ModelOptions::VERSION){
        std::cout<<"Reading Model Options file."<<endl;
        std::cout<<"Model Options version does not match!"<<endl;
        std::cout<<"Got '"<<str<<"' but expected '"<<ModelOptions::VERSION<<"'."<<endl;
        std::cout<<"Please update '"<<is.get_file_name()<<"'"<<endl;
        exit(-1);
    }
    cout<<ModelOptions::VERSION<<tb<<"# Model Options version"<<endl;
    
    //initial numbers-at-size options
    if (debug) cout<<"##Initial numbers-at-size Option 1:"<<endl;
    is>>opt1_InitNatZ;
    cout<<opt1_InitNatZ<<tb<<"#"<<opts1_InitNatZ(opt1_InitNatZ)<<endl;
    if (debug) cout<<"##Initial numbers-at-size Option 2:"<<endl;
    is>>opt2_InitNatZ;
    cout<<opt1_InitNatZ<<tb<<"#"<<opts1_InitNatZ(opt2_InitNatZ)<<endl;
    if (debug) cout<<"##Initial numbers-at-size Option 3:"<<endl;
    is>>opt3_InitNatZ;
    cout<<opt1_InitNatZ<<tb<<"#"<<opts1_InitNatZ(opt3_InitNatZ)<<endl;
    if (debug) cout<<"##Initial numbers-at-size sum-to-1 penalty weight:"<<endl;
    is>>wgtSumTo1_InitNatZ;
    cout<<wgtSumTo1_InitNatZ<<endl;
    
    //natural mortality options
    cout << "##Natural mortality options:"<<endl;
    is>>optParamNM;
    cout<<optParamNM<<tb<<"#"<<optsParamNM(optParamNM)<<endl;
    
    //growth parameterization options
    cout<<"##Growth parameterization options:"<<endl;
    is>>optGrowthParam;
    cout<<optGrowthParam<<tb<<"#"<<optsGrowthParam(optGrowthParam)<<endl;
    
    //growth pdf options
    cout<<"##Growth pdf options:"<<endl;
    is>>optGrowthPDF;
    cout<<optGrowthPDF<<tb<<"#"<<optsGrowthPDF(optGrowthPDF)<<endl;
    is>>maxGrowthZBEx;
    cout<<maxGrowthZBEx<<tb<<"#max extent of size bins for growth probabilities"<<endl;
            
    //likelihood penalty options for mean growth approaching negative increments
    cout<<"##Options for likelihood penalties on negative growth increments"<<endl;
    is>>minGrowthCW;
    cout<<minGrowthCW<<tb<<"#minGrowthCW"<<endl;
    is>>maxGrowthCW;
    cout<<maxGrowthCW<<tb<<"#maxGrowthCW"<<endl;
    is>>wgtNegGrowth;
    cout<<wgtNegGrowth<<tb<<"#wgtNegGrowth"<<endl;
    is>>epsNegGrowth;
    cout<<epsNegGrowth<<tb<<"#epsNegGrowth"<<endl;
    
    //terminal molt (prM2M) options
    cout<<"##Terminal molt (prM2M) options:"<<endl;
    {
        int nw; is>>nw;
        cout<<nw<<tb<<"#number of defined prM2M parameter combinations"<<endl;
        //--options for smoothness penalties on prM2M in likelihood
        is>>optPenSmthPrM2M;
        cout<<optPenSmthPrM2M<<tb<<"#option for calculating smoothness penalties on prM2M"<<endl;
        wgtPenSmthPrM2M.allocate(1,nw);
        is>>wgtPenSmthPrM2M;
        cout<<wgtPenSmthPrM2M<<tb<<"#weights for smoothness penalties on prM2M"<<endl;
        //--options for non-decreasing penalties on prM2M in likelihood
        is>>optPenNonDecPrM2M;
        cout<<optPenNonDecPrM2M<<tb<<"#option for calculating non-decreasing penalties on prM2M"<<endl;
        wgtPenNonDecPrM2M.allocate(1,nw);
        is>>wgtPenNonDecPrM2M;
        cout<<wgtPenNonDecPrM2M<<tb<<"#weights for non-decreasing penalties on prM2M"<<endl;
    }
    
    //nonparametric selectivity function (NPSel) options
    cout<<"##Nonparameteric selectivity function (NPSel) options:"<<endl;
    {
        int nw; is>>nw;
        cout<<nw<<tb<<"#number of defined NPSel parameter combinations"<<endl;
        //--options for smoothness penalties on prM2M in likelihood
        is>>optPenSmthNPSel;
        cout<<optPenSmthNPSel<<tb<<"#option for calculating smoothness penalties on NPSels"<<endl;
        wgtPenSmthNPSel.allocate(1,nw);
        is>>wgtPenSmthNPSel;
        cout<<wgtPenSmthNPSel<<tb<<"#weights for smoothness penalties on NPSels"<<endl;
    }
    
    //Empirical selectivity function options
    cout<<"##Empirical Selectivity Options:"<<endl;
    is>>(*ptrEmpiricalSelFcns);
    cout<<(*ptrEmpiricalSelFcns)<<endl;
    
    //Empirical selectivity function priors options
    cout<<"##Empirical Selectivity Function Priors Options:"<<endl;
    is>>(*ptrEmpiricalSelFcnPriors);
    cout<<(*ptrEmpiricalSelFcnPriors)<<endl;
    
    //effort extrapolation options
    cout<<"##Effort extrapolation scenarios:"<<endl;
    is>>(*ptrEffXtrapScenarios);
    cout<<(*ptrEffXtrapScenarios)<<endl;
    
    //likelihood penalties on F-devs
    cout<<"##Likelihood pnalties on F-devs:"<<endl;
    is>>cvFDevsPen;
    cout<<cvFDevsPen<<tb<<"#initial cv for F-devs penalties"<<endl;
    is>>phsDecrFDevsPen;
    cout<<phsDecrFDevsPen<<tb<<"#phase at which to start decreasing the penalties on F-devs"<<endl;
    is>>phsZeroFDevsPen;
    cout<<phsZeroFDevsPen<<tb<<"#phase at which to turn off the penalties on F-devs"<<endl;
    
    //likelihood penalty weights on non-zero sums for devs vectors
    cout<<"##Likelihood penalty weights on on non-zero sums for devs vectors:"<<endl;
    is>>wgtSqSumDevsPen;
    cout<<wgtSqSumDevsPen<<tb<<"#weight for penalties on non-zero sums for devs vectors"<<endl;
    is>>phsSqSumDevsPen;
    cout<<phsSqSumDevsPen<<tb<<"#min phase to apply penalty"<<endl;
    
    //OFL calculation options
    //--averaging options for capture rate/selectivity functions
    cout<<"#OFL capture rate/selectivity functions averaging options"<<endl;
    optOFLAvgCapRate.allocate(1,ptrMC->nFsh);
    optOFLAvgCapRate.initialize();
    for (int f=1;f<=ptrMC->nFsh;f++){
        is>>str; cout<<str<<"# fishery"<<tb;
        idx = wts::which(str,ptrMC->lblsFsh);
        cout<<idx<<tb;
        is>>optOFLAvgCapRate(idx);
        cout<<"= "<<optOFLAvgCapRate(idx)<<endl;
    }
    //--averaging periods
    cout<<"##OFL averaging periods"<<endl;
    oflNumYrsForAvgCapRate.allocate(1,ptrMC->nFsh);
    oflNumYrsForAvgCapRate.initialize();
    for (int f=1;f<=ptrMC->nFsh;f++){
        is>>str; cout<<str<<"# fishery"<<tb;
        idx = wts::which(str,ptrMC->lblsFsh);
        cout<<idx<<tb;
        is>>oflNumYrsForAvgCapRate(idx);
        cout<<"= "<<oflNumYrsForAvgCapRate(idx)<<endl;
    }
    //externally-calculated max capture rates
    cout<<"##Externally-calculated max capture rates for OFL calculations"<<endl;
    oflAvgCapRateInfo.allocate(1,ptrMC->nFsh);
    for (int f=1;f<=ptrMC->nFsh;f++){
        is>>str; cout<<str<<"# fishery"<<tb;
        idx = wts::which(str,ptrMC->lblsFsh);
        cout<<idx<<tb;
        is>>oflAvgCapRateInfo(idx);
        cout<<"= "<<oflAvgCapRateInfo(idx)<<endl;
    }
    
    //Iterative re-weighting options for size compositions
    cout<<"#Option for iterative re-weighting of size compositions"<<endl;
    is>>optIterativeReweighting;
    cout<<optIterativeReweighting<<tb<<"#"<<optsIterativeReweighting(optIterativeReweighting)<<endl;
    is>>phsIterativeReweighting;
    cout<<phsIterativeReweighting<<tb<<"#phase to start iterative re-weighting"<<endl;
    is>>maxIterations;
    cout<<maxIterations<<tb<<"#max iterations"<<endl;
    
    //MSE-related options
    cout<<"#MSE-related options"<<endl;
    is>>opModRecStatsMinYr;
    is>>opModRecStatsMaxYr;
    is>>HCR;
    cout<<HCR<<tb<<"#harvest control rule scenario"<<endl;
    int tst = 1;
    while (tst){
        is>>str;
        if (str=="HCR1"){
            cout<<"#--options for "<<str<<endl;
            is>>HCR1_avgMinYr;
            is>>HCR1_avgMaxYr;
            cout<<HCR1_avgMinYr<<tb<<HCR1_avgMaxYr<<tb<<"#min, max years for averaging"<<endl;
            if (HCR1_avgMaxYr==-1) HCR1_avgMaxYr = ptrMC->mxYr;
            tst=0;
        } else if (str=="HCR2"){
            cout<<"#--options for "<<str<<endl;
            is>>HCR2_rampID;
            cout<<HCR2_rampID<<tb<<tb<<"#ramp id"<<endl;
            tst=0;
        } else {tst=0;}
    }
    cout<<endl;
    
    //Projection options
    cout<<"#Projection-related options"<<endl;
    is>>(*ptrProjOpts);     cout<<"#--non-MCMC options"<<endl<<(*ptrProjOpts)    <<endl;
    is>>(*ptrProjOptsMCMC); cout<<"#--    MCMC options"<<endl<<(*ptrProjOptsMCMC)<<endl;
    
//    if (debug){
//        cout<<"enter 1 to continue : ";
//        cin>>debug;
//        if (debug<0) exit(1);
//    }
    
    //Simulation options
    cout<<"#In ModelOptions, reading simulation-related options"<<endl;
    is>>(*ptrSimOpts); cout<<"#--Simulation options"<<endl<<(*ptrSimOpts)<<endl;
    
    if (debug){
        cout<<"enter 1 to continue : ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    if (debug) cout<<"end ModelOptions::read(cifstream & is)"<<endl;
}

/***************************************************************
*   function to write to file in ADMB format                   *
***************************************************************/
void ModelOptions::write(ostream & os) {
    if (debug) cout<<"#start ModelOptions::write(ostream)"<<endl;
    os<<"#######################################"<<endl;
    os<<"#TCSAM02 Model Options File           #"<<endl;
    os<<"#######################################"<<endl;
    os<<ModelOptions::VERSION<<tb<<"# Model Options version"<<endl;

    //initial n-at-z options
    os<<"#----Initial Numbers-At-Size Options: Option 1"<<endl;
    for (int o=opts1_InitNatZ.indexmin();o<=opts1_InitNatZ.indexmax();o++) {
        os<<opts1_InitNatZ(o)<<endl;
    }
    os<<opt1_InitNatZ<<tb<<"#selected option"<<endl;
    os<<endl;
    os<<"#----Initial Numbers-At-Size Options: Option 2 (applies if Option 1 is 2 or 3)"<<endl;
    for (int o=opts2_InitNatZ.indexmin();o<=opts2_InitNatZ.indexmax();o++) {
        os<<opts2_InitNatZ(o)<<endl;
    }
    os<<opt2_InitNatZ<<tb<<"#selected option"<<endl;
    os<<endl;
    os<<"#----Initial Numbers-At-Size Options: Option 3 (applies if Option 1 is 3)"<<endl;
    for (int o=opts3_InitNatZ.indexmin();o<=opts3_InitNatZ.indexmax();o++) {
        os<<opts3_InitNatZ(o)<<endl;
    }
    os<<opt3_InitNatZ<<tb<<"#selected option"<<endl;
    os<<endl;
    os<<"#----Initial Numbers-At-Size Options: Option 4 (applies if Option 1 is 3 and Option 3 is 2)"<<endl;
    os<<wgtSumTo1_InitNatZ<<tb<<"#--penalty weight on sum-to-1 constraint"<<endl;
    os<<endl;
    
    //natural mortality options
    os<<"#----Options for parameterizing natural mortality"<<endl;
    for (int o=optsParamNM.indexmin();o<=optsParamNM.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsParamNM(o)<<endl;
    }
    os<<optParamNM<<tb<<"#selected option"<<endl;
    os<<endl;
    
    //growth parameterization options
    os<<"#----Growth parameterization options"<<endl;
    for (int o=optsGrowthPDF.indexmin();o<=optsGrowthPDF.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsGrowthPDF(o)<<endl;
    }
    os<<optGrowthPDF<<tb<<"#selected option"<<endl;

    //growth pdf options
    os<<"#----Growth pdf options"<<endl;
    for (int o=optsGrowthPDF.indexmin();o<=optsGrowthPDF.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsGrowthPDF(o)<<endl;
    }
    os<<optGrowthPDF<<tb<<"#selected option"<<endl;
    os<<"#------"<<endl;
    os<<maxGrowthZBEx<<tb<<"#max extent of size bins for growth probabilities"<<endl;

    //likelihood penalty options for mean growth approaching negative increments
    os<<"#----Options for likelihood penalties on negative growth increments"<<endl;
    os<<minGrowthCW <<tb<<"#min pre-molt CW to apply penalty on approaching negative growth increments"<<endl;
    os<<maxGrowthCW <<tb<<"#min pre-molt CW to apply penalty on approaching negative growth increments"<<endl;
    os<<wgtNegGrowth<<tb<<"#likelihood weight for penalty on approaching negative growth increments"<<endl;
    os<<epsNegGrowth<<tb<<"#eps parameter in posfun() for penalty on approaching negative growth increments"<<endl;
    os<<endl;
    
    //prM2M options
    //--smoothness likelihood options
    os<<"#----prM2M Options"<<endl;
    os<<wgtPenSmthPrM2M.size()<<tb<<"#number of prM2M parameter combinations"<<endl;
    os<<"#----Options for penalties on prM2M smoothness"<<endl;
    for (int o=optsPenSmthPrM2M.indexmin();o<=optsPenSmthPrM2M.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsPenSmthPrM2M(o)<<endl;
    }
    os<<optPenSmthPrM2M<<tb<<"#selected option"<<endl;
    os<<wgtPenSmthPrM2M<<tb<<"#weights for prM2M smoothness penalties"<<endl;
    //--non-decreasing likelihood options
    os<<"#----Options for penalties on non-decreasing prM2M"<<endl;
    for (int o=optsPenNonDecPrM2M.indexmin();o<=optsPenNonDecPrM2M.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsPenNonDecPrM2M(o)<<endl;
    }
    os<<optPenNonDecPrM2M<<tb<<"#selected option"<<endl;
    os<<wgtPenNonDecPrM2M<<tb<<"#weights for prM2M non-decreasing penalties"<<endl;
    os<<endl;
    
    //Nonparameteric selectivity function options
    //--smoothness likelihood options
    os<<"#----Nonparameteric selectivity function options"<<endl;
    os<<wgtPenSmthNPSel.size()<<tb<<"#number of parameter combinations"<<endl;
    os<<"#----Options for smoothness penalties on nonparametric selectivity functions"<<endl;
    for (int o=optsPenSmthNPSel.indexmin();o<=optsPenSmthNPSel.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsPenSmthNPSel(o)<<endl;
    }
    os<<optPenSmthNPSel<<tb<<"#selected option"<<endl;
    os<<wgtPenSmthNPSel<<tb<<"#weights for smoothness penalties on nonparametric selectivity functions"<<endl;
    os<<endl;
    
    //empirical selectivity functions options
    os<<"#----Empirical Selectivity Function Options"<<endl;
    os<<(*ptrEmpiricalSelFcns);
    os<<endl;

    //empirical selectivity function priors options
    os<<"#----Empirical Selectivity Function Priors Options"<<endl;
    os<<(*ptrEmpiricalSelFcnPriors);
    os<<endl;

    //effort extrapolation options
    os<<"#----Effort Extrapolation Scenarios"<<endl;
    os<<(*ptrEffXtrapScenarios);
    os<<endl;

    //F-devs penalties
    os<<"#----F-devs penalty options"<<endl;
    os<<cvFDevsPen<<tb<<"#initial cv for F-devs penalties"<<endl;
    os<<phsDecrFDevsPen<<tb<<"#phase at which to start decreasing the penalties on F-devs"<<endl;
    os<<phsZeroFDevsPen<<tb<<"#phase at which to turn off the penalties on F-devs"<<endl;
    os<<endl;
    
    //likelihood penalties on final values of dev vectors
    os<<"#----Likelihood penalty weights on on non-zero sums for devs vectors"<<endl;
    os<<wgtSqSumDevsPen<<tb<<"#weight for penalties on non-zero sums for devs vectors"<<endl;
    os<<phsSqSumDevsPen<<tb<<"#min phase to apply penalty"<<endl;
    os<<endl;
    
    //OFL options
    os<<"#----OFL Calculation Options"<<endl;
    os<<"#-----capture rate/selectivity function options"<<endl;
    for (int o=optsOFLAvgCapRate.indexmin();o<=optsOFLAvgCapRate.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsOFLAvgCapRate(o)<<endl;
    }
    os<<"#  Fishery    Option"<<endl;
    for (int f=1;f<=ptrMC->nFsh;f++){
        os<<tb<<ptrMC->lblsFsh(f)<<tb<<tb<<optOFLAvgCapRate(f)<<endl;
    }
    os<<"#-----averaging period (years)"<<endl;
    os<<"#  Fishery    Option"<<endl;
    for (int f=1;f<=ptrMC->nFsh;f++){
        os<<tb<<ptrMC->lblsFsh(f)<<tb<<tb<<oflNumYrsForAvgCapRate(f)<<endl;
    }
    os<<"#-----externally-calculated average capture rates"<<endl;
    os<<"#  Fishery    Option"<<endl;
    for (int f=1;f<=ptrMC->nFsh;f++){
        os<<tb<<ptrMC->lblsFsh(f)<<tb<<tb<<oflAvgCapRateInfo(f)<<endl;
    }
    os<<endl;
    
    //Iterative re-weighting options
    for (int o=optsIterativeReweighting.indexmin();o<=optsIterativeReweighting.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsIterativeReweighting(o)<<endl;
    }
    os<<optIterativeReweighting<<tb<<"#selected option"<<endl;
    os<<phsIterativeReweighting<<tb<<"#phase to start iterative re-weighting"<<endl;
    os<<maxIterations<<tb<<"#max number of iterations"<<endl;
    os<<endl;
        
    //MSE-related options
    os<<"#----MSE-related options"<<endl;
    os<<"#------Time period for recruitment projection statistics"<<endl;
    os<<"#min, max year for statistics"<<endl;
    os<<opModRecStatsMinYr<<tb<<opModRecStatsMaxYr<<endl;
    
    os<<"#------Harvest Control Rule"<<endl;
    os<<"# 1 - HCR1: "<<endl;
    os<<"# 2 - HCR2: "<<endl;
    os<<"# 3 - HCR3: "<<endl;
    os<<"# 4 - HCR4: "<<endl;
    os<<"# 5 - HCR5: "<<endl;
    os<<"# 6 - HCR6: "<<endl;
    os<<HCR<<tb<<"#selected harvest control rule scenario"<<endl;
    os<<"#---------HCR-specific options (uncomment for selected HCR)"<<endl;
    if (HCR==1){
        os<<"#--options for HCR1"<<endl;
        os<<HCR1_avgMinYr<<tb<<HCR1_avgMaxYr<<tb<<"#min, max years for female biomass averaging"<<endl;
    } else if (HCR==2){
        os<<"#--options for HCR2"<<endl;
        os<<"#Ramp ID: 1, 2 or 3"<<endl;
        os<<HCR2_rampID<<tb<<tb<<"#ramp id"<<endl;
    }
    os<<endl;
    
    //projection options
    os<<"#----Projection Options"<<endl;
    os<<"#------non-MCMC projection options "<<endl;
    os<<(*ptrProjOpts);
    os<<"#------    MCMC projection options "<<endl;
    os<<(*ptrProjOptsMCMC);
    
    //projection options
    os<<"#----Simulation Options"<<endl;
    os<<(*ptrSimOpts);
    
    os<<"#---------------------"<<endl<<endl;
    
    if (debug) cout<<"#end ModelOptions::write(ostream)"<<endl;
}

/**
 * Write ModelOptions info as an R list. TODO: complete this!
 * 
 * @param os - output stream to write to
 * @param nm - name to use for list
 * @param indent - number of tabs to indent
 */
void ModelOptions::writeToR(ostream& os, std::string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list("<<endl;
    indent++;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"initNatZ=list(opt1="<<opt1_InitNatZ<<cc<<"opt2="<<opt2_InitNatZ<<cc<<"opt3="<<opt3_InitNatZ<<cc<<"pen_wgt="<<wgtSumTo1_InitNatZ<<")"<<cc<<endl;
        os<<"natmort="<<optParamNM<<cc<<endl;
        os<<"growth=list(";
            os<<"optGrowthPDF="<<optGrowthPDF<<cc<<endl;
            os<<"maxGrowthZBEx="<<maxGrowthZBEx;
        os<<"),"<<endl;
        os<<"prM2M=list(";
            os<<"pen.smoothing="; wts::writeToR(os,wgtPenSmthPrM2M); os<<cc<<endl;
            os<<"pen.nondecreasing="; wts::writeToR(os,wgtPenNonDecPrM2M); 
        os<<"),"<<endl;
        os<<"npSel=list(";
            os<<"pen.smoothing="; wts::writeToR(os,wgtPenSmthNPSel); os<<"),"<<endl;
        os<<"cvFDevsPen="<<cvFDevsPen<<cc<<"phsDecr="<<phsDecrFDevsPen<<cc<<"phsZero="<<phsZeroFDevsPen<<cc
          <<"wgtLastDevPen="<<wgtSqSumDevsPen<<cc<<"phsLastDevsPen="<<phsSqSumDevsPen<<cc;
        os<<"empSelFcns="; ptrEmpiricalSelFcns->writeToR(os); os<<"),"<<endl;
        os<<"empSelFcnPriors="; ptrEmpiricalSelFcnPriors->writeToR(os); os<<"),"<<endl;
        os<<"effXtrapScenarios="; ptrEffXtrapScenarios->writeToR(os); os<<"),"<<endl;
        os<<"oflOptions=list(";
            os<<"optAvgCapRate="; wts::writeToR(os,optOFLAvgCapRate); os<<cc<<endl;
            os<<"numYears="; wts::writeToR(os,oflNumYrsForAvgCapRate); os<<cc<<endl;
            os<<"rateInfo="; wts::writeToR(os,oflAvgCapRateInfo); 
        os<<")"<<cc<<endl;
        os<<"itRewgtOptions=list(";
            os<<"option="<<optIterativeReweighting<<cc<<endl;
            os<<"phase="<<phsIterativeReweighting<<cc<<endl;
            os<<"maxIts="<<maxIterations;
        os<<")"<<cc<<endl;
        os<<"projOpts=";     ptrProjOpts->writeToR(os);     os<<cc<<endl;
        os<<"projOptsMCMC="; ptrProjOptsMCMC->writeToR(os); os<<cc<<endl;
        os<<"simOpts=";      ptrSimOpts->writeToR(os);      os<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
        os<<")";
}
/////////////////////////////////end ModelOptions/////////////////////////

