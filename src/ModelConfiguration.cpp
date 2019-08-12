#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"

//**********************************************************************
//  Includes
//      ModelConfiguration
//**********************************************************************
using namespace std;

const adstring ModelConfiguration::VERSION = "2016.11.15";

int ModelConfiguration::debug=0;
//--------------------------------------------------------------------------------
//          ModelConfiguration
//--------------------------------------------------------------------------------
int    ModelConfiguration::mnYr     = -1;//min model year
int    ModelConfiguration::asYr     = -1;//model assessment year
int    ModelConfiguration::mxYr     = -1;//max model year
int    ModelConfiguration::nSrv     = -1;//number of model surveys
int    ModelConfiguration::nFsh     = -1;//number of model fisheries
int    ModelConfiguration::nZBs     = -1;//number of model size bins
int    ModelConfiguration::jitter   = OFF;//flag to jitter initial parameter values
double ModelConfiguration::jitFrac  = 1.0;//fraction to jitter bounded parameter values
int    ModelConfiguration::resample = OFF;//flag to resample initial parameter values
double ModelConfiguration::vif      = 1.0;//variance inflation factor for resampling parameter values
/***************************************************************
*   creation                                                   *
***************************************************************/
ModelConfiguration::ModelConfiguration(){
    runOpMod=fitToPriors=INT_TRUE;//default to TRUE
}
/***************************************************************
*   destruction                                                *
***************************************************************/
ModelConfiguration::~ModelConfiguration(){
    if (debug) cout<<"destroying ModelConfiguration "<<this<<endl;
    if (debug) cout<<"destruction complete "<<endl;
}

/***************************************************************
*   function to read from file in ADMB format                  *
***************************************************************/
void ModelConfiguration::read(const adstring & fn) {
    if (debug) cout<<"ModelConfiguration::read(fn). Reading from '"<<fn<<"'"<<endl;
    cifstream strm(fn);
    read(strm);
    if (debug) cout<<"end ModelConfiguration::read(fn). Read from '"<<fn<<"'"<<endl;
}

/***************************************************************
*   function to write to file in ADMB format                   *
***************************************************************/
void ModelConfiguration::write(const adstring & fn) {
    if (debug) cout<<"#start ModelConfiguration::write(fn). Writing to '"<<fn<<"'"<<endl;
    ofstream strm(fn,ofstream::out|ofstream::trunc);
    write(strm); //write to file
    strm.close();
    if (debug) cout<<"#end ModelConfiguration::write(fn). Wrote to '"<<fn<<"'"<<endl;
}

/***************************************************************
*   function to read from file in ADMB format                  *
***************************************************************/
void ModelConfiguration::read(cifstream & is) {
    if (debug) cout<<"ModelConfiguration::read(cifstream & is)"<<endl;
    cout<<"ModelConfiguration input file name: '"<<is.get_file_name()<<"'"<<endl;
    adstring parent = wts::getParentFolder(is);
    cout<<"parent folder is '"<<parent<<"'"<<endl;
    adstring ver;
    is>>ver;
    if (ver!=ModelConfiguration::VERSION){
        std::cout<<"Reading Model Configuration file."<<endl;
        std::cout<<"Model Configuration version does not match!"<<endl;
        std::cout<<"Got '"<<ver<<"' but expected '"<<ModelConfiguration::VERSION<<"'."<<endl;
        std::cout<<"Please update '"<<is.get_file_name()<<"'"<<endl;
        exit(-1);
    }
    is>>cfgName;
    if (debug) cout<<cfgName<<endl;
    is>>mnYr;      //min model year
    is>>asYr;      //assessment year
    mxYr = asYr-1; //max model year
    is>>nZBs;      //number of model size bins
    if (debug){
        cout<<mnYr <<tb<<"#model min year"<<endl;
        cout<<mxYr <<tb<<"#model max year"<<endl;
        cout<<nZBs<<tb<<"#number of size bins"<<endl;
    }
    zMidPts.allocate(1,nZBs); 
    zCutPts.allocate(1,nZBs+1); 
    onesZMidPts.allocate(1,nZBs); onesZMidPts = 1.0;
    is>>zCutPts;
    for (int z=1;z<=nZBs;z++) zMidPts(z) = 0.5*(zCutPts(z)+zCutPts(z+1));
    if (debug){
        cout<<"#size bins (mm CW)"<<endl;
        cout<<zMidPts<<endl;
        cout<<"#size bin cut points (mm CW)"<<endl;
        cout<<zCutPts <<endl;
        cout<<"enter 1 to continue : ";
        cin>>debug;
        if (debug<0) exit(1);
    }
        
    is>>nFsh; //number of fisheries
    lblsFsh.allocate(1,nFsh);
    for (int i=1;i<=nFsh;i++) is>>lblsFsh(i); //labels for fisheries
    if (debug){
        cout<<nFsh<<tb<<"#number of fisheries"<<endl;
        for (int i=1;i<=nFsh;i++) cout<<lblsFsh(i)<<tb;
        cout<<tb<<"#labels for fisheries"<<endl;
    }
    
    is>>nSrv; //number of surveys
    lblsSrv.allocate(1,nSrv);
    for (int i=1;i<=nSrv;i++) is>>lblsSrv(i);//labels for surveys
    if (debug){
        cout<<nSrv<<tb<<"#number of surveys"<<endl;
        for (int i=1;i<=nSrv;i++) cout<<lblsSrv(i)<<tb;
        cout<<tb<<"#labels for surveys"<<endl;
    }
    
    adstring str1;
    is>>str1; runOpMod    = wts::getBooleanType(str1);//run population model?
    is>>str1; fitToPriors = wts::getBooleanType(str1);//fit priors?
    
    is>>fnMPI;//model parameters information file
    is>>fnMDS;//model datasets file
    is>>fnMOs;//model options file
    
    fnMPI = wts::concatenateFilePaths(parent,fnMPI);
    fnMDS = wts::concatenateFilePaths(parent,fnMDS);
    fnMOs = wts::concatenateFilePaths(parent,fnMOs);
    
    is>>str1; ModelConfiguration::jitter = wts::getOnOffType(str1);
    is>>ModelConfiguration::jitFrac;
    is>>str1; ModelConfiguration::resample = wts::getOnOffType(str1);
    is>>ModelConfiguration::vif;
    
    //convert model quantities to csv strings
    csvYrs  =qt+str(mnYr)+qt; for (int y=(mnYr+1);y<=mxYr;    y++) csvYrs  += cc+qt+str(y)+qt;
    csvYrsP1=qt+str(mnYr)+qt; for (int y=(mnYr+1);y<=(mxYr+1);y++) csvYrsP1 += cc+qt+str(y)+qt;
    csvSXs=qt+tcsam::getSexType(1)     +qt; for (int i=2;i<=tcsam::nSXs;i++) csvSXs += cc+qt+tcsam::getSexType(i)     +qt;
    csvMSs=qt+tcsam::getMaturityType(1)+qt; for (int i=2;i<=tcsam::nMSs;i++) csvMSs += cc+qt+tcsam::getMaturityType(i)+qt;
    csvSCs=qt+tcsam::getShellType(1)   +qt; for (int i=2;i<=tcsam::nSCs;i++) csvSCs += cc+qt+tcsam::getShellType(i)   +qt;
    csvZCs=wts::to_qcsv(zCutPts);
    csvZBs=wts::to_qcsv(zMidPts);
    csvFsh=wts::to_qcsv(lblsFsh);
    csvSrv=wts::to_qcsv(lblsSrv);
    
    dimYrsToR   = "y=c("+csvYrs+")";
    dimYrsP1ToR = "y=c("+csvYrsP1+")";
    dimSXsToR   = tcsamDims::getSXsForR(1,tcsam::nSXs);
    dimMSsToR   = tcsamDims::getMSsForR(1,tcsam::nMSs);
    dimSCsToR   = tcsamDims::getSCsForR(1,tcsam::nSCs);
    dimZCsToR   = "zc=c("+wts::to_qcsv(zCutPts)+")";
    dimZBsToR   = "z=c("+wts::to_qcsv(zMidPts)+")";
    dimZPsToR   = "zp=c("+wts::to_qcsv(zMidPts)+")";
    dimFshToR   = "f=c("+wts::to_qcsv(lblsFsh)+")";
    dimSrvToR   = "v=c("+wts::to_qcsv(lblsSrv)+")";
    
    if (debug){
        cout<<wts::getBooleanType(runOpMod)   <<"   #run operating model?"<<endl;
        cout<<wts::getBooleanType(fitToPriors)<<"   #fit to priors?"<<endl;
        cout<<fnMPI<<"   #model parameters configuration file"<<endl;
        cout<<fnMDS<<"   #model datasets file"<<endl;
        cout<<fnMOs<<"   #model options file"<<endl;
        cout<<wts::getOnOffType(ModelConfiguration::jitter)<<tb<<"#jitter?"<<endl;
        cout<<ModelConfiguration::jitFrac<<tb<<"#jitter fraction"<<endl;
        cout<<wts::getOnOffType(ModelConfiguration::resample)<<tb<<"#resmple?"<<endl;
        cout<<ModelConfiguration::vif<<tb<<"#variance inflation factor"<<endl;
        cout<<"enter 1 to continue : ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    
    if (debug) cout<<"end ModelConfiguration::read(cifstream & is)"<<endl;
}

/**
 * Set max model year (for retrospective model runs).
 * 
 * @param yr - new max model year
 */
void ModelConfiguration::setMaxModelYear(int yr){
    mxYr = yr;
    asYr = mxYr+1;
    csvYrs  =qt+str(mnYr)+qt; for (int y=(mnYr+1);y<=mxYr;    y++) csvYrs  += cc+qt+str(y)+qt;
    csvYrsP1=qt+str(mnYr)+qt; for (int y=(mnYr+1);y<=(mxYr+1);y++) csvYrsP1 += cc+qt+str(y)+qt;
    dimYrsToR   = "y=c("+csvYrs+")";
    dimYrsP1ToR = "y=c("+csvYrsP1+")";
}
/***************************************************************
*   function to write to file in ADMB format                   *
***************************************************************/
void ModelConfiguration::write(ostream & os) {
    if (debug) cout<<"#start ModelConfiguration::write(ostream)"<<endl;
    os<<"#######################################################"<<endl;
    os<<"#TCSAM02 Model Configuration File                     #"<<endl;
    os<<"#######################################################"<<endl;
    os<<ModelConfiguration::VERSION<<tb<<"#ModelConfiguration version"<<endl;
    os<<cfgName<<tb<<"#Model configuration name"<<endl;
    os<<mnYr<<tb<<"#Min model year"<<endl;
    os<<asYr<<tb<<"#Assessment year"<<endl;
    os<<nZBs<<tb<<"#Number of model size classes"<<endl;
    os<<"#size bin cut points"<<endl;
    os<<zCutPts <<endl;
    
    os<<nFsh<<tb<<"#number of fisheries"<<endl;
    for (int i=1;i<=nFsh;i++) cout<<lblsFsh(i)<<tb;
        cout<<tb<<"#labels for fisheries"<<endl;
    os<<nSrv<<tb<<"#number of surveys"<<endl;
    for (int i=1;i<=nSrv;i++) cout<<lblsSrv(i)<<tb;
        cout<<tb<<"#labels for surveys"<<endl;    
        
    os<<wts::getBooleanType(runOpMod)   <<tb<<"#run operating model?"<<endl;
    os<<wts::getBooleanType(fitToPriors)<<tb<<"#fit priors?"<<endl;
    
    os<<fnMPI<<tb<<"#Model parameters info file"<<endl;
    os<<fnMDS<<tb<<"#Model datasets file"<<endl;
    os<<fnMOs<<tb<<"#Model options file"<<endl;
    
    os<<wts::getOnOffType(ModelConfiguration::jitter)<<tb<<"#jitter?"<<endl;
    os<<ModelConfiguration::jitFrac<<tb<<"#jitter fraction"<<endl;
    os<<wts::getOnOffType(ModelConfiguration::resample)<<tb<<"#resmple?"<<endl;
    os<<ModelConfiguration::vif<<tb<<"#variance inflation factor"<<endl;

    if (debug) cout<<"#end ModelConfiguration::write(ostream)"<<endl;
}

/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void ModelConfiguration::writeToR(ostream& os, std::string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list("<<endl;
    indent++;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"version='"<<ModelConfiguration::VERSION<<"'"<<cc<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"configName='"<<cfgName<<"'"<<cc<<endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"dims=list("<<endl;
        indent++;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"y=list(n="<<mxYr-mnYr<<cc<<"mny="<<mnYr<<cc<<"asy="<<asYr<<cc<<"mxy="<<mxYr<<cc<<
                           "nms=c("<<csvYrsP1<<"),vls="<<mnYr<<":"<<asYr<<"),"<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"x=list(n="<<tcsam::nSXs<<",nms=c("<<tcsamDims::formatForR(csvSXs)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"m=list(n="<<tcsam::nMSs<<",nms=c("<<tcsamDims::formatForR(csvMSs)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"s=list(n="<<tcsam::nSCs<<",nms=c("<<tcsamDims::formatForR(csvSCs)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"z=list(n="<<nZBs<<",nms=c("<<csvZBs<<"),vls=c("<<wts::to_csv(zMidPts)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"zc=list(n="<<nZBs+1<<",nms=c("<<csvZCs<<"),vls=c("<<wts::to_csv(zCutPts)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"f=list(n="<<nFsh<<cc<<"nms=c("<<wts::replace('_',' ',csvFsh)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"v=list(n="<<nSrv<<cc<<"nms=c("<<wts::replace('_',' ',csvSrv)<<"))"<<endl;
        indent--;
    for (int n=0;n<indent;n++) os<<tb;
        os<<")"<<cc;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"flags=list(";
        os<<"runOpMod="<<runOpMod<<cc;
        os<<"fitToPriors="<<fitToPriors<<"),";
        os<<endl;
    for (int n=0;n<indent;n++) os<<tb; os<<"fnMPI='"<<wts::replace('\\','/',fnMPI)<<"',"<<endl;
    for (int n=0;n<indent;n++) os<<tb; os<<"fnMDS='"<<wts::replace('\\','/',fnMDS)<<"',"<<endl;
    for (int n=0;n<indent;n++) os<<tb; os<<"fnMOs='"<<wts::replace('\\','/',fnMOs)<<"'"<<endl;
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
        os<<")";
}
/////////////////////////////////end ModelConfiguration/////////////////////////

