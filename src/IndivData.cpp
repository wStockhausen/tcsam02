#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelData.hpp"
#include "SummaryFunctions.hpp"

using namespace tcsam;

//**********************************************************************
//  Includes
//      GrowthData
//      ChelaHeightData
//      MaturityOgiveData
//**********************************************************************
/* debug flag for growth data */
int GrowthData::debug  = 0;
/* debug flag for chela height data */
int ChelaHeightData::debug  = 0;
/* debug flag for chela height data */
int MaturityOgiveData::debug  = 0;
//----------------------------------------------------------------------
//          GrowthData
//----------------------------------------------------------------------
/* keyword for chela height data */
const adstring GrowthData::KW_GROWTH_DATA = "GROWTH_DATA";
/**
 * Constructor.
 */
GrowthData::GrowthData(){
    nObs_x.allocate(1,tcsam::nSXs);
    obsYears_xn.allocate(1,tcsam::nSXs);    //note: allocation incomplete
    inpData_xcn.allocate(1,tcsam::nSXs,1,3);//note: allocation incomplete
}
/**
 * Destructor.
 */
GrowthData::~GrowthData(){}
/**
 * Set the maximum year in which to include growth data.
 * 
 * @param mxYr - the max year in which to include growth data
 */
void GrowthData::setMaxYear(int mxYr){
    if (debug) {
        cout     <<"Starting GrowthData::setMaxYear("<<mxYr<<")"<<endl;
        rpt::echo<<"Starting GrowthData::setMaxYear("<<mxYr<<")"<<endl;
    }
    //define temporary arrays
    imatrix newYears_xn;  newYears_xn.allocate(1,tcsam::nSXs);    //note: allocation incomplete
    d3_array newData_xcn; newData_xcn.allocate(1,tcsam::nSXs,1,3);//note: allocation incomplete
    //determine number of observations to keep, by sex
    for (int x=1;x<=nSXs;x++){
        int nx = 0;
        for (int iy=obsYears_xn(x).indexmin(); iy<=obsYears_xn(x).indexmax(); iy++){
            if (obsYears_xn(x)(iy)<=mxYr) nx++;
        }
        //finish allocations
        newYears_xn(x).allocate(1,nx);
        newData_xcn(x,1).allocate(1,nx);
        newData_xcn(x,2).allocate(1,nx);
        newData_xcn(x,3).allocate(1,nx);
        //copy observations
        int iyp = 1;
        for (int iy=obsYears_xn(x).indexmin(); iy<=obsYears_xn(x).indexmax(); iy++){
            if (obsYears_xn(x)(iy)<=mxYr) {
                newYears_xn(x,iyp)   = obsYears_xn(x)(iy);
                newData_xcn(x,1,iyp) = inpData_xcn(x,1,iy);
                newData_xcn(x,2,iyp) = inpData_xcn(x,2,iy);
                newData_xcn(x,3,iyp) = inpData_xcn(x,3,iy);
                iyp++;
            }
        }
        nObs_x(x) = nx;
        obsYears_xn(x).deallocate(); obsYears_xn(x).allocate(1,nx);
        obsYears_xn(x) = newYears_xn(x);
        inpData_xcn(x,1).deallocate(); inpData_xcn(x,1).allocate(1,nx); 
        inpData_xcn(x,2).deallocate(); inpData_xcn(x,2).allocate(1,nx); 
        inpData_xcn(x,3).deallocate(); inpData_xcn(x,3).allocate(1,nx); 
        for (int ic=1;ic<=3;ic++) inpData_xcn(x,ic) = newData_xcn(x,ic);
        if (debug) {
            rpt::echo<<"obsYears_xn("<<x<<") = "<<obsYears_xn(x)<<endl;
            rpt::echo<<"inpData{"<<x<<") = "<<endl<<inpData_xcn(x)<<endl;
        }
    }
    if (debug) {
        cout     <<"Finished GrowthData::setMaxYear("<<mxYr<<")"<<endl;
        rpt::echo<<"Finished GrowthData::setMaxYear("<<mxYr<<")"<<endl;
    }
}
/**
 * Read growth data from input file stream.
 * 
 * @param is - the file stream to read from
 */
void GrowthData::read(cifstream & is){
    if (debug){
        cout<<"start GrowthData::read(...) "<<this<<std::endl;
        cout<<"#------------------------------------------"<<std::endl;
        cout<<"#file name is "<<is.get_file_name()<<std::endl;
        cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        cout<<"Apparent error reading GrowthData."<<std::endl;
        cout<<"#file name is "<<is.get_file_name()<<std::endl;
        cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    if (!(str==KW_GROWTH_DATA)){
        cout<<"#Error reading effort data from "<<is.get_file_name()<<std::endl;
        cout<<"Expected keyword '"<<KW_GROWTH_DATA<<"' but got '"<<str<<"'"<<std::endl;
        cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    is>>name;
    rpt::echo<<name<<tb<<"#dataset name"<<endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>llWgt;
    rpt::echo<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    is>>nSXs;//number of sex categories
    rpt::echo<<nSXs<<tb<<"#number of sex categories"<<std::endl;
    nObs_x.initialize();
    for (int x=1;x<=nSXs;x++){
        is>>str; int xp = tcsam::getSexType(str);
        rpt::echo<<tcsam::getSexType(xp)<<tb<<"#sex"<<std::endl;
        is>>nObs_x(xp);
        rpt::echo<<nObs_x(xp)<<tb<<"#number of observations"<<std::endl;
        dmatrix tmp(1,nObs_x(xp),1,3);
        is>>tmp;
        rpt::echo<<"#year"<<tb<<"pre-molt size"<<tb<<"post-molt size"<<std::endl;
        rpt::echo<<tmp<<std::endl;
        inpData_xcn(xp,1).allocate(1,nObs_x(xp));
        inpData_xcn(xp,2).allocate(1,nObs_x(xp));
        inpData_xcn(xp,3).allocate(1,nObs_x(xp));
        inpData_xcn(xp) = trans(tmp);
        obsYears_xn(xp).allocate(1,nObs_x(xp));
        ivector itmp(inpData_xcn(xp,1));
        obsYears_xn(xp) = itmp;
    }
    
    if (debug) cout<<"end GrowthData::read(...) "<<this<<std::endl;
}
/**
 * Write data to output stream.
 * 
 * @param os - the output stream
 */
void GrowthData::write(ostream & os){
    if (debug) cout<<"start GrowthData::write(...) "<<this<<std::endl;
    os<<KW_GROWTH_DATA<<tb<<"#required keyword"<<std::endl;
    os<<name<<tb<<"#dataset name"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    os<<nSXs<<tb<<"#number of sex categories"<<std::endl;
    for (int x=1;x<=tcsam::nSXs;x++){
        if (nObs_x(x)>0){
            os<<tcsam::getSexType(x)<<tb<<"#sex"<<std::endl;
            os<<nObs_x(x)<<tb<<"#number of observations"<<std::endl;
            os<<"#year"<<tb<<"pre-molt size"<<tb<<"post-molt size"<<std::endl;
            os<<trans(inpData_xcn(x))<<std::endl;
        }
    }
    if (debug) cout<<"end GrowthData::write(...) "<<this<<std::endl;
}
/**
 * Write to an output stream as an R list with name given by name variable. 
 * 
 * @param os - the output stream
 * @param nm - the name of the list (ignored)
 * @param indent - the number of tabs to indent, initially
 */
void GrowthData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"GrowthData::writing to R"<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"`"<<name<<"`"<<"=list("<<std::endl;
        indent++; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"name="<<qt<<name<<qt<<cc;
            os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc;
            os<<"llWgt="<<llWgt<<cc<<std::endl;
            adstring colnames="'y','z','zp'";
            for (int x=1;x<=tcsam::nSXs;x++){
                if (nObs_x(x)>0){
                    for (int n=0;n<indent;n++) os<<tb;
                    os<<tcsam::getSexType(x)<<"="; wts::writeToR(os,trans(inpData_xcn(x)),colnames); os<<cc<<std::endl;
                }
            }//x
            os<<"NULL"<<std::endl;
        indent--;
    for (int n=0;n<indent;n++) os<<tb; os<<")"<<std::endl;
    if (debug) cout<<"GrowthData::done writing to R"<<std::endl;
}
/////////////////////////////////end GrowthData/////////////////////////
//----------------------------------------------------------------------
//          ChelaHeightData
//----------------------------------------------------------------------
/* keyword for chela height data */
const adstring ChelaHeightData::KW_CHELAHEIGHT_DATA = "CHELAHEIGHT_DATA";
/**
 * Constructor for class.
 */
ChelaHeightData::ChelaHeightData(){}
/**
 * Destructor for class.
 */
ChelaHeightData::~ChelaHeightData(){}
/**
 * Set the maximum year in which to include chela height data.
 * 
 * @param mxYr - the max year in which to include chela height data
 */
void ChelaHeightData::setMaxYear(int mxYr, const dvector& zCs){
    if (debug) {
        cout     <<"start ChelaHeightData::setMaxYear("<<mxYr<<")"<<endl;    
        rpt::echo<<"start ChelaHeightData::setMaxYear("<<mxYr<<")"<<endl;    
    }
    //count data with years <= mxYr
    int n = 0;
    for (int i=obsYear_n.indexmin();i<=obsYear_n.indexmax();i++)
        if (obsYear_n(i)<=mxYr) n++;
    nObs = n;
    //revise input data
    n = 1;
    dmatrix newData_nc(1,nObs,1,4);
    for (int i=obsYear_n.indexmin();i<=obsYear_n.indexmax();i++)
        if (obsYear_n(i)<=mxYr) newData_nc(n++) = inpData_nc(i);
    
    inpData_nc.deallocate(); inpData_nc.allocate(1,nObs,1,5);
    for (int i=1;i<=nObs;i++) inpData_nc(i) = newData_nc(i);
    
    obsYear_n.deallocate();         obsYear_n.allocate(1,nObs);
    obsSize_n.deallocate();         obsSize_n.allocate(1,nObs);
    obsSS_n.deallocate();           obsSS_n.allocate(1,nObs);
    obsPrMat_n.deallocate();        obsPrMat_n.allocate(1,nObs);
    obsSizeBinIndex_n.deallocate(); obsSizeBinIndex_n.allocate(1,nObs);
   
    obsYear_n  = wts::to_ivector(column(inpData_nc,1));
    obsSize_n  = column(inpData_nc,2);
    obsSS_n    = column(inpData_nc,4);
    obsPrMat_n = column(inpData_nc,5);
    obsSizeBinIndex_n = 0;//can't fill this in yet, need to use calcSizeBinIndices(zCs)
    calcSizeBinIndices(zCs);
    
    rpt::echo<<"obsYear_n         = "<<obsYear_n <<endl;
    rpt::echo<<"obsSize_n         = "<<obsSize_n <<endl;
    rpt::echo<<"obsSS_n           = "<<obsSS_n   <<endl;
    rpt::echo<<"obsPrMat_n        = "<<obsPrMat_n<<endl;
    rpt::echo<<"obsSizeBinIndex_n = "<<obsSizeBinIndex_n<<endl;
    
    if (debug) {
        cout     <<"end ChelaHeightData::setMaxYear("<<mxYr<<")"<<endl;    
        rpt::echo<<"end ChelaHeightData::setMaxYear("<<mxYr<<")"<<endl;    
    }
}
/**
 * Calculates indices for model size bins corresponding to observed sizes.
 * 
 * @param zCs - model size bin cutpoints
 */
void ChelaHeightData::calcSizeBinIndices(const dvector& zCs){
    rpt::echo<<"Starting ChelaHeightData::calcSizeBinIndices(zCs)"<<endl;
    obsSizeBinIndex_n = wts::assignBinIndices(obsSize_n,zCs,0,0);
    rpt::echo<<"zCs = "<<zCs<<endl;
    for (int n=1;n<=nObs;n++) rpt::echo<<obsSize_n(n)<<tb<<obsSizeBinIndex_n(n)<<tb<<zCs(obsSizeBinIndex_n(n))<<endl;
    rpt::echo<<"Finished ChelaHeightData::calcSizeBinIndices(zCs)"<<endl;
}
/**
 * Read chela height data from input file stream.
 * 
 * @param is - the file stream to read from
 */
void ChelaHeightData::read(cifstream & is){
    if (debug){
        cout<<"start ChelaHeightData::read(...) "<<this<<std::endl;
        cout<<"#------------------------------------------"<<std::endl;
        cout<<"#file name is "<<is.get_file_name()<<std::endl;
        cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        cout<<"Apparent error reading ChelaHeightData."<<std::endl;
        cout<<"#file name is "<<is.get_file_name()<<std::endl;
        cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    if (!(str==KW_CHELAHEIGHT_DATA)){
        cout<<"#Error reading effort data from "<<is.get_file_name()<<std::endl;
        cout<<"Expected keyword '"<<KW_CHELAHEIGHT_DATA<<"' but got '"<<str<<"'"<<std::endl;
        cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    is>>name;
    rpt::echo<<name<<tb<<"#dataset name"<<endl;
    is>>survey;
    rpt::echo<<survey<<tb<<"#survey name"<<endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>llWgt;
    rpt::echo<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    is>>nObs;//number of observations
    rpt::echo<<nObs<<tb<<"#number of observations"<<std::endl;
    inpData_nc.allocate(1,nObs,1,4);
    is>>inpData_nc;
    rpt::echo<<"#year    size    nIndivs     fraction mature"<<std::endl<<inpData_nc<<std::endl;
    
    obsYear_n.allocate(1,nObs);
    obsSize_n.allocate(1,nObs);
    obsSS_n.allocate(1,nObs);
    obsPrMat_n.allocate(1,nObs);
    obsSizeBinIndex_n.allocate(1,nObs);
   
    obsYear_n  = wts::to_ivector(column(inpData_nc,1));
    obsSize_n  = column(inpData_nc,2);
    obsSS_n    = column(inpData_nc,3);
    obsPrMat_n = column(inpData_nc,4);
    obsSizeBinIndex_n = 0;//can't fill this in yet, need to use calcSizeBinIndices(zCs)
    
    rpt::echo<<"obsYear_n  = "<<obsYear_n<<endl;
    rpt::echo<<"obsSize_n  = "<<obsSize_n<<endl;
    rpt::echo<<"obsSS_n    = "<<obsSS_n  <<endl;
    rpt::echo<<"obsPrMat_n = "<<obsPrMat_n<<endl;
    
    if (debug) cout<<"end ChelaHeightData::read(...) "<<this<<std::endl;    
}
/**
 * Write data to output stream.
 * 
 * @param os - the output stream
 */
void ChelaHeightData::write(ostream & os){
    if (debug) cout<<"start ChelaHeightData::write(...) "<<this<<std::endl;
    os<<KW_CHELAHEIGHT_DATA<<tb<<"#required keyword"<<std::endl;
    os<<name<<tb<<"#dataset name"<<std::endl;
    os<<survey<<tb<<"#survey name"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    os<<nObs<<tb<<"#number of observations"<<std::endl;
    os<<"#year    size    nIndivs     fraction mature"<<std::endl;
    os<<inpData_nc<<std::endl;
    if (debug) cout<<"end ChelaHeightData::write(...) "<<this<<std::endl;
}
/**
 * Write to an output stream as an R list with name given by name variable. 
 * 
 * @param os - the output stream
 * @param nm - the name of the list (ignored)
 * @param indent - the number of tabs to indent, initially
 */
void ChelaHeightData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"ChelaHeightData::writing to R"<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"`"<<name<<"`"<<"=list("<<std::endl;
        indent++; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"name="<<qt<<name<<qt<<cc;
            os<<"survey="<<qt<<survey<<qt<<cc;
            os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc;
            os<<"llWgt="<<llWgt<<cc<<std::endl;
            for (int n=0;n<indent;n++) os<<tb;
            adstring colnames="'y','z','nIndivs','fraction mature'";
            os<<"data="; wts::writeToR(os,inpData_nc,colnames); os<<std::endl;
        indent--;
    for (int n=0;n<indent;n++) os<<tb; os<<")"<<std::endl;
    if (debug) cout<<"ChelaHeightData::done writing to R"<<std::endl;
}
/////////////////////////////////end ChelaHeightData/////////////////////////
//----------------------------------------------------------------------
//          MaturityOgiveData
//----------------------------------------------------------------------
/* keyword for maturity ogive data */
const adstring MaturityOgiveData::KW_MATURITYOGIVE_DATA = "MATURITYOGIVE_DATA";
/**
 * Constructor for class.
 */
MaturityOgiveData::MaturityOgiveData(){}
/**
 * Destructor for class.
 */
MaturityOgiveData::~MaturityOgiveData(){}
/**
 * Set the maximum year in which to include maturity ogive data.
 * 
 * @param mxYr - the max year in which to include maturity ogive data
 */
void MaturityOgiveData::setMaxYear(int mxYr){
    if (debug) {
        cout     <<"start MaturityOgiveData::setMaxYear("<<mxYr<<")"<<endl;    
        rpt::echo<<"start MaturityOgiveData::setMaxYear("<<mxYr<<")"<<endl;    
    }
    //count data with years <= mxYr
    int n = 0;
    for (int i=obsYear_n.indexmin();i<=obsYear_n.indexmax();i++)
        if (obsYear_n(i)<=mxYr) n++;
    nObs = n;
    //revise input data
    n = 1;
    dmatrix newData_nc(1,nObs,1,5);
    for (int i=obsYear_n.indexmin();i<=obsYear_n.indexmax();i++)
        if (obsYear_n(i)<=mxYr) newData_nc(n++) = inpData_nc(i);
    
    inpData_nc.deallocate(); inpData_nc.allocate(1,nObs,1,5);
    for (int i=1;i<=nObs;i++) inpData_nc(i) = newData_nc(i);
    
    obsYear_n.deallocate();  obsYear_n.allocate(1,nObs);
    obsSize_n.deallocate();  obsSize_n.allocate(1,nObs);
    obsZBI_n.deallocate();   obsZBI_n.allocate(1,nObs);
    obsSS_n.deallocate();    obsSS_n.allocate(1,nObs);
    obsPrMat_n.deallocate(); obsPrMat_n.allocate(1,nObs);
   
    obsYear_n  = wts::to_ivector(column(inpData_nc,1));
    obsSize_n  = column(inpData_nc,2);
    obsZBI_n   = wts::to_ivector(column(inpData_nc,3));
    obsSS_n    = column(inpData_nc,4);
    obsPrMat_n = column(inpData_nc,5);
    
    if (debug){
        rpt::echo<<"obsYear_n  = "<<obsYear_n <<endl;
        rpt::echo<<"obsSize_n  = "<<obsSize_n <<endl;
        rpt::echo<<"obsZBI_n   = "<<obsZBI_n  <<endl;
        rpt::echo<<"obsSS_n    = "<<obsSS_n   <<endl;
        rpt::echo<<"obsPrMat_n = "<<obsPrMat_n<<endl;
    
        cout     <<"end MaturityOgiveData::setMaxYear("<<mxYr<<")"<<endl;    
        rpt::echo<<"end MaturityOgiveData::setMaxYear("<<mxYr<<")"<<endl;    
    }
}
/**
 * Calculates matrix to re-map model size bins maturity ogive size bins.
 * 
 * @param zCs - model size bin cutpoints
 */
void MaturityOgiveData::calcSizeBinRemapper(const dvector& zCs){
    if (debug) rpt::echo<<"Starting MaturityOgiveData::calcSizeBinRemapper(zCs)"<<endl;
    zbRemapper = tcsam::getRebinMatrix(zCs,cutpts);
    if (debug){
        rpt::echo<<"model zCs = "<<zCs<<endl;
        rpt::echo<<"MO    zCs = "<<cutpts<<endl;
        rpt::echo<<"sizeBinRemapper: "<<endl; rpt::echo<<zbRemapper<<endl;
        rpt::echo<<"Finished MaturityOgiveData::calcSizeBinRemapper(zCs)"<<endl;
    }
}
/**
 * Read maturity ogive data from input file stream.
 * 
 * @param is - the file stream to read from
 */
void MaturityOgiveData::read(cifstream & is){
    if (debug){
        cout<<"start MaturityOgiveData::read(...) "<<this<<std::endl;
        cout<<"#------------------------------------------"<<std::endl;
        cout<<"#file name is "<<is.get_file_name()<<std::endl;
        cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        cout<<"Apparent error reading MaturityOgiveData."<<std::endl;
        cout<<"#file name is "<<is.get_file_name()<<std::endl;
        cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    if (!(str==KW_MATURITYOGIVE_DATA)){
        cout<<"#Error reading effort data from "<<is.get_file_name()<<std::endl;
        cout<<"Expected keyword '"<<KW_MATURITYOGIVE_DATA<<"' but got '"<<str<<"'"<<std::endl;
        cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    is>>name;
    rpt::echo<<name<<tb<<"#dataset name"<<endl;
    is>>survey;
    rpt::echo<<survey<<tb<<"#survey name"<<endl;
    is>>str; sex = tcsam::getSexType(str); 
    rpt::echo<<tcsam::getSexType(sex)<<tb<<"#sex"<<endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>llWgt;
    rpt::echo<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    is>>nZBs;//number of size bins
    rpt::echo<<nZBs<<tb<<"#number of size bins"<<std::endl;
    cutpts.allocate(1,nZBs+1);
    is>>cutpts;
    rpt::echo<<"#size bin cut points"<<std::endl<<cutpts<<std::endl;
    is>>nObs;//number of observations
    rpt::echo<<nObs<<tb<<"#number of observations"<<std::endl;
    inpData_nc.allocate(1,nObs,1,5);
    is>>inpData_nc;
    rpt::echo<<"#year    size    index  nIndivs     fraction mature"<<std::endl<<inpData_nc<<std::endl;
    
    obsYear_n.allocate(1,nObs);
    obsSize_n.allocate(1,nObs);
    obsZBI_n.allocate(1,nObs);
    obsSS_n.allocate(1,nObs);
    obsPrMat_n.allocate(1,nObs);
   
    obsYear_n  = wts::to_ivector(column(inpData_nc,1));
    obsSize_n  = column(inpData_nc,2);
    obsZBI_n   = wts::to_ivector(column(inpData_nc,3));
    obsSS_n    = column(inpData_nc,4);
    obsPrMat_n = column(inpData_nc,5);
    
    rpt::echo<<"obsYear_n  = "<<obsYear_n <<endl;
    rpt::echo<<"obsSize_n  = "<<obsSize_n <<endl;
    rpt::echo<<"obsZBI_n   = "<<obsZBI_n  <<endl;
    rpt::echo<<"obsSS_n    = "<<obsSS_n   <<endl;
    rpt::echo<<"obsPrMat_n = "<<obsPrMat_n<<endl;
    
    if (debug) cout<<"end MaturityOgiveData::read(...) "<<this<<std::endl;    
}
/**
 * Write data to output stream.
 * 
 * @param os - the output stream
 */
void MaturityOgiveData::write(ostream & os){
    if (debug) cout<<"start MaturityOgiveData::write(...) "<<this<<std::endl;
    os<<KW_MATURITYOGIVE_DATA<<tb<<"#required keyword"<<std::endl;
    os<<name<<tb<<"#dataset name"<<std::endl;
    os<<survey<<tb<<"#survey name"<<std::endl;
    os<<tcsam::getSexType(sex)<<tb<<"#sex"<<endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    os<<nZBs<<tb<<"#number of size bins"<<std::endl;
    os<<"#size bin cut points"<<std::endl;
    os<<cutpts<<std::endl;
    os<<nObs<<tb<<"#number of observations"<<std::endl;
    os<<"#year    size    zbIndex  nIndivs     fraction mature"<<std::endl;
    os<<inpData_nc<<std::endl;
    if (debug) cout<<"end MaturityOgiveData::write(...) "<<this<<std::endl;
}
/**
 * Write to an output stream as an R list with name given by name variable. 
 * 
 * @param os - the output stream
 * @param nm - the name of the list (ignored)
 * @param indent - the number of tabs to indent, initially
 */
void MaturityOgiveData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"MaturityOgiveData::writing to R"<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"`"<<name<<"`"<<"=list("<<std::endl;
        indent++; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"name="<<qt<<name<<qt<<cc;
            os<<"survey="<<qt<<survey<<qt<<cc;
            os<<"sex="<<qt<<tcsam::getSexType(sex)<<qt<<cc;
            os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc;
            os<<"llWgt="<<llWgt<<cc<<std::endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"cutpts="; wts::writeToR(os,cutpts);  os<<cc<<std::endl;
            for (int n=0;n<indent;n++) os<<tb;
            adstring colnames="'y','z','zbIndex','nIndivs','fraction mature'";
            os<<"data="; wts::writeToR(os,inpData_nc,colnames); os<<std::endl;
        indent--;
    for (int n=0;n<indent;n++) os<<tb; os<<")"<<std::endl;
    if (debug) cout<<"MaturityOgiveData::done writing to R"<<std::endl;
}
/////////////////////////////////end MaturityOgiveData/////////////////////////
