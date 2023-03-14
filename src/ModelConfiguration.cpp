#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"

//**********************************************************************
//  Includes
//      ModelConfiguration
//**********************************************************************
using namespace std;

const adstring ModelConfiguration::VERSION = "2023.03.13";

int ModelConfiguration::debug=0;
//--------------------------------------------------------------------------------
//          ModelConfiguration
//--------------------------------------------------------------------------------
int    ModelConfiguration::mnYr       = -1;//min model year
int    ModelConfiguration::asYr       = -1;//model assessment year
int    ModelConfiguration::mxYr       = -1;//max model year
int    ModelConfiguration::mnYrAvgRec = -1;//min year to calculate average recruitment for OFL
int    ModelConfiguration::mxYrOffsetAvgRec = -1;//max year to calculate average recruitment for OFL
int    ModelConfiguration::yRetro     =  0;//number of retrospective years
int    ModelConfiguration::nSrv       = -1;//number of model surveys
int    ModelConfiguration::nFsh       = -1;//number of model fisheries
int    ModelConfiguration::nZBs       = -1;//number of model size bins
double ModelConfiguration::maxZC      = -1;//max size bin cutpoint
double ModelConfiguration::minZC      = -1;//min size bin cutpoint
double ModelConfiguration::delZ       = -1;//size bin size
dvector ModelConfiguration::zMidPts; /* size bin midpoints (CW in mm) */
dvector ModelConfiguration::zCutPts; /* size bin cutpoints (CW in mm) */
int    ModelConfiguration::jitter     = OFF;//flag to jitter initial parameter values
double ModelConfiguration::jitFrac    = 1.0;//fraction to jitter bounded parameter values
int    ModelConfiguration::resample   = OFF;//flag to resample initial parameter values
double ModelConfiguration::vif        = 1.0;//variance inflation factor for resampling parameter values
dvector ModelConfiguration::maxZs(1,tcsam::nSXs);//sex-specific max sizes
ivector ModelConfiguration::maxZBs(1,tcsam::nSXs);//sex-specific max size bins corresponding to maxZs
double ModelConfiguration::maxZRec;//max size possible at recruitment
int ModelConfiguration::maxZBRec;//index to max possible size bin for recruitment
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
    adstring strD; //use to create PRINT2B strings
    strD = "ModelConfiguration input file name: '"+is.get_file_name()+"'"; PRINT2B1(strD)
    adstring parent = wts::getParentFolder(is);
    strD="parent folder is '"+parent+"'"; PRINT2B1(strD);
    adstring ver;
    is>>ver;//model configuration file version
    if (ver!=ModelConfiguration::VERSION){
        std::cout<<"Reading Model Configuration file."<<endl;
        std::cout<<"Model Configuration version does not match!"<<endl;
        std::cout<<"Got '"<<ver<<"' but expected '"<<ModelConfiguration::VERSION<<"'."<<endl;
        std::cout<<"Please update '"<<is.get_file_name()<<"'"<<endl;
        exit(-1);
    }
    is>>cfgName;         //model configuration name
    if (debug) cout<<cfgName<<endl;
    is>>mnYr;            //min model year
    is>>asYr;            //assessment year
    mxYr = asYr-1;       //max model year
    is>>mnYrAvgRec;      //min year to calculate average recruitment for OFL
    is>>mxYrOffsetAvgRec;//max year to calculate average recruitment for OFL
    is>>maxZC;           //maximum size bin cutpoint
    is>>minZC;           //minimum size bin cutpoint
    is>>delZ;            //bin size
    is>>maxZs;           //sex-specific max sizes
    is>>maxZRec;         //max size at recruitment
    if (debug){
        cout<<mnYr <<tb<<"#model min year"<<endl;
        cout<<mxYr <<tb<<"#model max year"<<endl;
        cout<<mnYrAvgRec       <<tb<<"#min year for OFL average recruitment calculation"<<endl;
        cout<<mxYrOffsetAvgRec <<tb<<"#offset from max year for OFL average recruitment calculation"<<endl;
        cout<<maxZC   <<tb<<"#max size bin cutpoint"<<endl;
        cout<<minZC  <<tb<<"#min size bin cutpoint"<<endl;
        cout<<delZ   <<tb<<"#size bin size"<<endl;
        cout<<maxZs  <<tb<<"#max sizes, by sex"<<endl;
        cout<<maxZRec<<tb<<"#max size at recruitment"<<endl;
    }
    nZBs = (maxZC-minZC)/delZ;
    zMidPts.allocate(1,nZBs); 
    zCutPts.allocate(1,nZBs+1); 
    onesZMidPts.allocate(1,nZBs); onesZMidPts = 1.0;
    zCutPts(1) = minZC;
    for (int z=1;z<=nZBs;z++) {
      zCutPts(z+1) = zCutPts(z)+delZ;
      zMidPts(z) = 0.5*(zCutPts(z)+zCutPts(z+1));
    }
    maxZBs = nZBs;
    for (int x=1;x<=tcsam::nSXs;x++){
        for (int z=1;z<=nZBs;z++) 
            if ((zCutPts[z]<=maxZs[x])&&(maxZs[x]<zCutPts[z+1])) maxZBs[x] = z;
    }
    maxZBRec = nZBs;
    for (int z=1;z<=nZBs;z++) 
        if ((zCutPts[z]<=maxZRec)&&(maxZRec<zCutPts[z+1])) maxZBRec = z;
    PRINT2B2("#max sizes, by sex: ",maxZs)
    PRINT2B2("#max size bins, by sex: ",maxZBs)
    PRINT2B2("#max size at recruitment: ",maxZRec)
    PRINT2B2("#max size bin for recruitment: ",maxZBRec)
    PRINT2B2("#size bins (mm CW): ",zMidPts)
    PRINT2B2("#size bin cut points (mm CW): ",zCutPts)
        
    is>>nFsh; //number of fisheries
    if (nFsh){
        lblsFsh.allocate(1,nFsh);
        for (int i=1;i<=nFsh;i++) is>>lblsFsh(i); //labels for fisheries
    }
    PRINT2B2("#number of fisheries: ",nFsh);
    if (nFsh){
        PRINT2B1("#labels for fisheries:");
        for (int i=1;i<=nFsh;i++) {
            adstring tmp = lblsFsh(i);
            PRINT2B1(tmp);
        }
    }
    
    is>>nSrv; //number of surveys
    if (nSrv){
        lblsSrv.allocate(1,nSrv);
        for (int i=1;i<=nSrv;i++) is>>lblsSrv(i);//labels for surveys
    }
    PRINT2B2("#number of surveys: ",nSrv);
    if (nSrv){
        PRINT2B1("#labels for surveys:");
        for (int i=1;i<=nSrv;i++) {
            adstring tmp = lblsSrv(i);
            PRINT2B1(tmp);
        }
    }
    
    adstring str1;
    is>>str1; runOpMod    = wts::getBooleanType(str1);//run population model?
    is>>str1; fitToPriors = wts::getBooleanType(str1);//fit priors?
    PRINT2B2("run operating model only? ",runOpMod)
    PRINT2B2("fit to pripors?           ",fitToPriors)
    
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
    if (nFsh) csvFsh=wts::to_qcsv(lblsFsh);
    if (nSrv) csvSrv=wts::to_qcsv(lblsSrv);
    
    dimYrsToR   = "y=c("+csvYrs+")";
    dimYrsP1ToR = "y=c("+csvYrsP1+")";
    dimSXsToR   = tcsamDims::getSXsForR(1,tcsam::nSXs);
    dimMSsToR   = tcsamDims::getMSsForR(1,tcsam::nMSs);
    dimSCsToR   = tcsamDims::getSCsForR(1,tcsam::nSCs);
    dimZCsToR   = "zc=c("+wts::to_qcsv(zCutPts)+")";
    dimZBsToR   = "z=c("+wts::to_qcsv(zMidPts)+")";
    dimZPsToR   = "zp=c("+wts::to_qcsv(zMidPts)+")";
    if (nFsh) dimFshToR   = "f=c("+wts::to_qcsv(lblsFsh)+")";
    if (nSrv) dimSrvToR   = "v=c("+wts::to_qcsv(lblsSrv)+")";
    
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

///**
// * Set max model year (for retrospective model runs).
// * 
// * @param yr - new max model year
// */
//void ModelConfiguration::setMaxModelYear(int yr){
//    mxYr = yr;
//    asYr = mxYr+1;
//    csvYrs  =qt+str(mnYr)+qt; for (int y=(mnYr+1);y<=mxYr;    y++) csvYrs  += cc+qt+str(y)+qt;
//    csvYrsP1=qt+str(mnYr)+qt; for (int y=(mnYr+1);y<=(mxYr+1);y++) csvYrsP1 += cc+qt+str(y)+qt;
//    dimYrsToR   = "y=c("+csvYrs+")";
//    dimYrsP1ToR = "y=c("+csvYrsP1+")";
//}
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
    os<<mnYrAvgRec<<tb<<"#Min year for OFL average recruitment calculation"<<endl;
    os<<mxYrOffsetAvgRec<<tb<<"#Offset to maxYr for OFL average recruitment calculation"<<endl;
    os<<maxZC  <<tb<<"#maximum size bin cutpoint"<<endl;
    os<<minZC  <<tb<<"#minimum size bin cutpoint"<<endl;
    os<<delZ   <<tb<<"#size bin size"<<endl;
    os<<maxZs  <<tb<<"#max sizes, by sex"<<endl;
    os<<maxZRec<<tb<<"#max size at recruitment"<<endl;
    
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
            os<<"y=list(n="<<mxYr-mnYr<<cc<<
                       "mny="<<mnYr<<cc<<"asy="<<asYr<<cc<<"mxy="<<mxYr<<cc<<
                       "mnyAvgRec="<<mnYrAvgRec<<cc<<"mxyOffsetAvgRec="<<mxYrOffsetAvgRec<<cc<<"yRetro="<<yRetro<<cc<<
                       "nms=c("<<csvYrsP1<<"),vls="<<mnYr<<":"<<asYr<<"),"<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"x=list(n="<<tcsam::nSXs<<",nms=c("<<tcsamDims::formatForR(csvSXs)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"m=list(n="<<tcsam::nMSs<<",nms=c("<<tcsamDims::formatForR(csvMSs)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"s=list(n="<<tcsam::nSCs<<",nms=c("<<tcsamDims::formatForR(csvSCs)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"zLims=list(maxZRec="<<maxZRec<<",maxZs=c("<<wts::to_csv(maxZs)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"z=list(n="<<nZBs<<",nms=c("<<csvZBs<<"),vls=c("<<wts::to_csv(zMidPts)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"zc=list(n="<<nZBs+1<<",nms=c("<<csvZCs<<"),vls=c("<<wts::to_csv(zCutPts)<<"))"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"f=list(n="<<nFsh;
            if (nFsh) os<<cc<<"nms=c("<<wts::replace('_',' ',csvFsh)<<")";
            os<<")"<<cc<<endl;
        for (int n=0;n<indent;n++) os<<tb;
            os<<"v=list(n="<<nSrv;
            if (nSrv) os<<cc<<"nms=c("<<wts::replace('_',' ',csvSrv)<<")";
            os<<")"<<endl;
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

