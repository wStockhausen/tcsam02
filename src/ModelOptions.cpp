#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelOptions.hpp"

//**********************************************************************
//  Includes
//      ModelConfiguration
//**********************************************************************
using namespace std;

//--------------------------------------------------------------------------------
//          ModelOptions
//--------------------------------------------------------------------------------
int ModelOptions::debug = 0;
const adstring ModelOptions::VERSION = "2017.02.27";

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
    
    //growth function options
    optsGrowth.allocate(0,1);
    optsGrowth(0) = "use gamma probability distribution (like TCSAM2013)"; 
    optsGrowth(1) = "use cumulative gamma distribution (like Gmacs)";
    
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
    //--estimation
    optsEffXtrEst.allocate(0,2);
    optsEffXtrEst(0) = "no extrapolation"; 
    optsEffXtrEst(1) = "use calculated average capture rate/average effort"; 
    optsEffXtrEst(2) = "estimate extrapolation parameters via likelihood";
    //--capture rate averaging
    optsEffXtrAvgFc.allocate(0,2);
    optsEffXtrAvgFc(0) = "no averaging"; 
    optsEffXtrAvgFc(1) = "average fully-selected capture rate";
    optsEffXtrAvgFc(2) = "average mean size-specific capture rate";
    
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
    is>>optInitNatZ;
    cout<<optInitNatZ<<tb<<"#"<<optsInitNatZ(optInitNatZ)<<endl;
    
    //growth options
    is>>optGrowth;
    cout<<optGrowth<<tb<<"#"<<optsGrowth(optGrowth)<<endl;
    
    //terminal molt (prM2M) options
    int nw; is>>nw;
    cout<<nw<<tb<<"#number of defined prM2M parameter combinations"<<endl;
    wgtPenSmthPrM2M.allocate(1,nw);
    wgtPenNonDecPrM2M.allocate(1,nw);
    //--options for smoothness penalties on prM2M in likelihood
    is>>optPenSmthPrM2M;
    cout<<optPenSmthPrM2M<<tb<<"#option for calculating smoothness penalties on prM2M"<<endl;
    is>>wgtPenSmthPrM2M;
    cout<<wgtPenSmthPrM2M<<tb<<"#weights for smoothness penalties on prM2M"<<endl;
    //--options for non-decreasing penalties on prM2M in likelihood
    is>>optPenNonDecPrM2M;
    cout<<optPenNonDecPrM2M<<tb<<"#option for calculating non-decreasing penalties on prM2M"<<endl;
    is>>wgtPenNonDecPrM2M;
    cout<<wgtPenNonDecPrM2M<<tb<<"#weights for non-decreasing penalties on prM2M"<<endl;
    
    //effort extrapolation options
    //--effort extrapolation estimation options
    optEffXtrEst.allocate(1,ptrMC->nFsh);
    for (int f=1;f<=ptrMC->nFsh;f++){
        is>>str; cout<<str<<"# fishery"<<tb;
        idx = wts::which(str,ptrMC->lblsFsh);
        cout<<idx<<tb;
        is>>optEffXtrEst(idx);
        cout<<"= "<<optEffXtrEst(idx)<<endl;
    }
    //--effort extrapolation capture rate averaging options
    optEffXtrAvgFc.allocate(1,ptrMC->nFsh);
    for (int f=1;f<=ptrMC->nFsh;f++){
        is>>str; cout<<str<<"# fishery"<<tb;
        idx = wts::which(str,ptrMC->lblsFsh);
        cout<<idx<<tb;
        is>>optEffXtrAvgFc(idx);
        cout<<"= "<<optEffXtrAvgFc(idx)<<endl;
    }
    
    //likelihood penalties on F-devs
    is>>cvFDevsPen;
    cout<<cvFDevsPen<<tb<<"#initial cv for F-devs penalties"<<endl;
    is>>phsDecrFDevsPen;
    cout<<phsDecrFDevsPen<<tb<<"#phase at which to start decreasing the penalties on F-devs"<<endl;
    is>>phsZeroFDevsPen;
    cout<<phsZeroFDevsPen<<tb<<"#phase at which to turn off the penalties on F-devs"<<endl;
    is>>wgtLastDevsPen;
    
    //likelihood penalties on last element in devs vectors
    cout<<wgtLastDevsPen<<tb<<"#weight for last-devs penalties"<<endl;
    is>>phsLastDevsPen;
    cout<<phsLastDevsPen<<tb<<"#min phase to apply penalty"<<endl;
    
    //OFL calculation options
    //--averaging options for capture rate/selectivity functions
    cout<<"#OFL capture rate/selectivity functions averaging options"<<endl;
    optOFLAvgCapRate.allocate(1,ptrMC->nFsh);
    for (int f=1;f<=ptrMC->nFsh;f++){
        is>>str; cout<<str<<"# fishery"<<tb;
        idx = wts::which(str,ptrMC->lblsFsh);
        cout<<idx<<tb;
        is>>optOFLAvgCapRate(idx);
        cout<<"= "<<optOFLAvgCapRate(idx)<<endl;
    }
    //--averaging periods
    cout<<"#OFL averaging periods"<<endl;
    oflNumYrsForAvgCapRate.allocate(1,ptrMC->nFsh);
    for (int f=1;f<=ptrMC->nFsh;f++){
        is>>str; cout<<str<<"# fishery"<<tb;
        idx = wts::which(str,ptrMC->lblsFsh);
        cout<<idx<<tb;
        is>>optOFLAvgCapRate(idx);
        cout<<"= "<<optOFLAvgCapRate(idx)<<endl;
    }
    //externally-calculated max capture rates
    cout<<"#externally-calculated max capture rates for OFL calculations"<<endl;
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
    
    //growth options
    os<<"#----Growth Function Options"<<endl;
    for (int o=optsGrowth.indexmin();o<=optsGrowth.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsGrowth(o)<<endl;
    }
    os<<optGrowth<<tb<<"#selected option"<<endl;

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
    os<<"#----Effort Extrapolation Options"<<endl;
    os<<"#------estimation options"<<endl;
    for (int o=optsEffXtrEst.indexmin();o<=optsEffXtrEst.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsEffXtrEst(o)<<endl;
    }
    os<<"#  Fishery    Option"<<endl;
    for (int f=1;f<=ptrMC->nFsh;f++){
        os<<tb<<ptrMC->lblsFsh(f)<<tb<<tb<<optEffXtrEst(f)<<endl;
    }
    os<<"#------capture rate averaging options"<<endl;
    for (int o=optsEffXtrAvgFc.indexmin();o<=optsEffXtrAvgFc.indexmax();o++) {
        os<<"#"<<o<<" - "<<optsEffXtrAvgFc(o)<<endl;
    }
    os<<"#  Fishery    Option"<<endl;
    for (int f=1;f<=ptrMC->nFsh;f++){
        os<<tb<<ptrMC->lblsFsh(f)<<tb<<tb<<optEffXtrAvgFc(f)<<endl;
    }

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
        os<<"initNatZ="<<optInitNatZ<<cc<<"natmort="<<optParamNM<<cc<<"growth="<<optGrowth<<cc<<endl;
        os<<"prM2M=list(";
            os<<"wgtSmthLgtPrMat="; wts::writeToR(os,wgtPenSmthPrM2M); os<<cc<<endl;
            os<<"wgtNonDecLgtPrMat="; wts::writeToR(os,wgtPenNonDecPrM2M); os<<")"<<endl;
        os<<"cvFDevsPen="<<cvFDevsPen<<cc<<"phsDecr="<<phsDecrFDevsPen<<cc<<"phsZero="<<phsZeroFDevsPen<<cc
          <<"wgtLastDevPen="<<wgtLastDevsPen<<cc<<"phsLastDevsPen="<<phsLastDevsPen<<cc;
        os<<"effXtr=list(";
            os<<"optEffXtrEst="; wts::writeToR(os,optEffXtrEst); os<<cc<<endl;
            os<<"optEffXtrAvgFc="; wts::writeToR(os,optEffXtrAvgFc); os<<")"<<endl;
        os<<"oflOptions=list(";
            os<<"optAvgCapRate="; wts::writeToR(os,optOFLAvgCapRate); os<<cc<<endl;
            os<<"numYears="; wts::writeToR(os,oflNumYrsForAvgCapRate); os<<cc<<endl;
            os<<"rateInfo="; wts::writeToR(os,oflAvgCapRateInfo); os<<")"<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
        os<<")";
}
/////////////////////////////////end ModelOptions/////////////////////////
