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
//      ModelOptions
//**********************************************************************
using namespace std;

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
int EffXtrapScenarios::debug = 0;
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
//          ModelOptions
//--------------------------------------------------------------------------------
int ModelOptions::debug = 0;
const adstring ModelOptions::VERSION = "2017.08.24";

ModelOptions::ModelOptions(ModelConfiguration& mc){
    ptrMC=&mc;
    
    //initial n-at-z options
    optsInitNatZ.allocate(0,1);
    optsInitNatZ(0) = "build up n-at-z from recruitments (like TCSAM2013)"; 
    optsInitNatZ(1) = "calculate initial n-at-z using equilibrium calculations (like Gmacs)";
    
    //options for natural mortality parameterization
    optsParamNM.allocate(0,1);
    optsParamNM(0) = "use log-scale parameterization (default)";
    optsParamNM(1) = "use TCSAM2013 parameterization (arithmetic scale)"; 
    
    //growth parameterization options
    optsGrowthParam.allocate(0,1);
    optsGrowthParam(0) = "TCSAM2013 parameterization (ln-scale intercept, slope)"; 
    optsGrowthParam(1) = "parameterization based on min, max pre-molt sizes";
    
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
    
    //effort extrapolation options
    ptrEffXtrapScenarios = new EffXtrapScenarios(mc);
//    //--estimation
//    optsEffXtrEst.allocate(0,2);
//    optsEffXtrEst(0) = "no extrapolation"; 
//    optsEffXtrEst(1) = "use calculated average capture rate/average effort"; 
//    optsEffXtrEst(2) = "estimate extrapolation parameters via likelihood";
//    //--capture rate averaging
//    optsEffXtrAvgFc.allocate(0,2);
//    optsEffXtrAvgFc(0) = "no averaging"; 
//    optsEffXtrAvgFc(1) = "average fully-selected capture rate";
//    optsEffXtrAvgFc(2) = "average mean size-specific capture rate";
    
    //options for OFL calculations: capture rate/selectivity function averaging
    optsOFLAvgCapRate.allocate(0,1);
    optsOFLAvgCapRate(0) = "average max capture rates, selectivity functions (like TCSAM2013)";
    optsOFLAvgCapRate(1) = "average size-specific capture rates";
}
/***************************************************************
*   function to read from file in ADMB format                  *
***************************************************************/
void ModelOptions::read(cifstream & is) {
    if (debug) cout<<"ModelOptions::read(cifstream & is)"<<endl;
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
    cout<<"##Initial numbers-at-size options:"<<endl;
    is>>optInitNatZ;
    cout<<optInitNatZ<<tb<<"#"<<optsInitNatZ(optInitNatZ)<<endl;
    
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
    
    //terminal molt (prM2M) options
    cout<<"##Terminal molt (prM2M) options:"<<endl;
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
    
    //effort extrapolation options
    cout<<"##Effort extrapolation scenarios:"<<endl;
    is>>(*ptrEffXtrapScenarios);
//    //--effort extrapolation estimation options
//    cout<<"##Effort extrapolation estimation options:"<<endl;
//    optEffXtrEst.allocate(1,ptrMC->nFsh);
//    for (int f=1;f<=ptrMC->nFsh;f++){
//        is>>str; cout<<str<<"fishery"<<tb;
//        idx = wts::which(str,ptrMC->lblsFsh);
//        cout<<idx<<tb;
//        is>>optEffXtrEst(idx);
//        cout<<"= "<<optEffXtrEst(idx)<<endl;
//    }
//    //--effort extrapolation estimation options
//    cout<<"##Effort extrapolation estimation options:"<<endl;
//    optEffXtrAvgFc.allocate(1,ptrMC->nFsh);
//    for (int f=1;f<=ptrMC->nFsh;f++){
//        is>>str; cout<<str<<"fishery"<<tb;
//        idx = wts::which(str,ptrMC->lblsFsh);
//        cout<<idx<<tb;
//        is>>optEffXtrAvgFc(idx);
//        cout<<"= "<<optEffXtrAvgFc(idx)<<endl;
//    }
    cout<<(*ptrEffXtrapScenarios)<<endl;
    
    //likelihood penalties on F-devs
    cout<<"##Likelihood pnalties on F-devs:"<<endl;
    is>>cvFDevsPen;
    cout<<cvFDevsPen<<tb<<"#initial cv for F-devs penalties"<<endl;
    is>>phsDecrFDevsPen;
    cout<<phsDecrFDevsPen<<tb<<"#phase at which to start decreasing the penalties on F-devs"<<endl;
    is>>phsZeroFDevsPen;
    cout<<phsZeroFDevsPen<<tb<<"#phase at which to turn off the penalties on F-devs"<<endl;
    
    //likelihood penalties on last element in devs vectors
    cout<<"##Likelihood penalties on last element in devs vectors:"<<endl;
    is>>wgtLastDevsPen;
    cout<<wgtLastDevsPen<<tb<<"#weight for last-devs penalties"<<endl;
    is>>phsLastDevsPen;
    cout<<phsLastDevsPen<<tb<<"#min phase to apply penalty"<<endl;
    
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
    
    if (debug) cout<<"end ModelOptions::read(cifstream & is)"<<endl;
    if (debug){
        cout<<"enter 1 to continue : ";
        cin>>debug;
        if (debug<0) exit(1);
    }
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
    os<<"#----Initial Numbers-At-Size Options"<<endl;
    for (int o=optsInitNatZ.indexmin();o<=optsInitNatZ.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsInitNatZ(o)<<endl;
    }
    os<<optInitNatZ<<tb<<"#selected option"<<endl;
    
    //natural mortality options
    os<<"#----Options for parameterizing natural mortality"<<endl;
    for (int o=optsParamNM.indexmin();o<=optsParamNM.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsParamNM(o)<<endl;
    }
    os<<optParamNM<<tb<<"#selected option"<<endl;
    
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
    
    //effort extrapolation options
    os<<"#----Effort Extrapolation Scenarios"<<endl;
    os<<(*ptrEffXtrapScenarios);
//    os<<"#------estimation options"<<endl;
//    for (int o=optsEffXtrEst.indexmin();o<=optsEffXtrEst.indexmax();o++) {
//        os<<"#"<<o<<" - "<<optsEffXtrEst(o)<<endl;
//    }
//    os<<"#  Fishery    Option"<<endl;
//    for (int f=1;f<=ptrMC->nFsh;f++){
//        os<<tb<<ptrMC->lblsFsh(f)<<tb<<tb<<optEffXtrEst(f)<<endl;
//    }
//    os<<"#------capture rate averaging options"<<endl;
//    for (int o=optsEffXtrAvgFc.indexmin();o<=optsEffXtrAvgFc.indexmax();o++) {
//        os<<"#"<<o<<" - "<<optsEffXtrAvgFc(o)<<endl;
//    }
//    os<<"#  Fishery    Option"<<endl;
//    for (int f=1;f<=ptrMC->nFsh;f++){
//        os<<tb<<ptrMC->lblsFsh(f)<<tb<<tb<<optEffXtrAvgFc(f)<<endl;
//    }

    //F-devs penalties
    os<<"#----F-devs penalty options"<<endl;
    os<<cvFDevsPen<<tb<<"#initial cv for F-devs penalties"<<endl;
    os<<phsDecrFDevsPen<<tb<<"#phase at which to start decreasing the penalties on F-devs"<<endl;
    os<<phsZeroFDevsPen<<tb<<"#phase at which to turn off the penalties on F-devs"<<endl;
    
    //likelihood penalties on final values of dev vectors
    os<<"#----Last dev penalty options"<<endl;
    os<<wgtLastDevsPen<<tb<<"#weight for last-dev penalties"<<endl;
    os<<phsDecrFDevsPen<<tb<<"#min phase to apply penalty"<<endl;
    
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
        os<<"initNatZ="<<optInitNatZ<<cc<<"natmort="<<optParamNM<<cc<<"growth="<<optGrowthPDF<<cc<<endl;
        os<<"prM2M=list(";
            os<<"wgtSmthLgtPrMat="; wts::writeToR(os,wgtPenSmthPrM2M); os<<cc<<endl;
            os<<"wgtNonDecLgtPrMat="; wts::writeToR(os,wgtPenNonDecPrM2M); os<<")"<<endl;
        os<<"cvFDevsPen="<<cvFDevsPen<<cc<<"phsDecr="<<phsDecrFDevsPen<<cc<<"phsZero="<<phsZeroFDevsPen<<cc
          <<"wgtLastDevPen="<<wgtLastDevsPen<<cc<<"phsLastDevsPen="<<phsLastDevsPen<<cc;
        os<<"effXtrapScenarios="; ptrEffXtrapScenarios->writeToR(os); os<<"),"<<endl;
        os<<"oflOptions=list(";
            os<<"optAvgCapRate="; wts::writeToR(os,optOFLAvgCapRate); os<<cc<<endl;
            os<<"numYears="; wts::writeToR(os,oflNumYrsForAvgCapRate); os<<cc<<endl;
            os<<"rateInfo="; wts::writeToR(os,oflAvgCapRateInfo); os<<")"<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
        os<<")";
}
/////////////////////////////////end ModelOptions/////////////////////////

