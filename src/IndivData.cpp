#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelData.hpp"

using namespace tcsam;

//**********************************************************************
//  Includes
//      GrowthData
//      ChelaHeightData
//**********************************************************************
/* debug flag for growth data */
int GrowthData::debug  = 0;
/* debug flag for chela height data */
int ChelaHeightData::debug  = 0;
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
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<name<<tb<<"#dataset name"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    os<<nSXs<<tb<<"#number of sex categories"<<std::endl;
    for (int x=1;x<=tcsam::nSXs;x++){
        if (nObs_x(x)>0){
            os<<tcsam::getSexType(x)<<tb<<"#sex"<<std::endl;
            os<<nObs_x(x)<<tb<<"#number of observations"<<std::endl;
            rpt::echo<<"#year"<<tb<<"pre-molt size"<<tb<<"post-molt size"<<std::endl;
            rpt::echo<<trans(inpData_xcn(x))<<std::endl;
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
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>llWgt;
    rpt::echo<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    is>>nObs;//number of observations
    rpt::echo<<nObs<<tb<<"#number of observations"<<std::endl;
    inpData_nc.allocate(1,nObs,1,4);
    is>>inpData_nc;
    rpt::echo<<"#year    size    nIndivs     fraction mature"<<std::endl<<inpData_nc<<std::endl;
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
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    os<<nObs<<tb<<"#number of observations"<<std::endl;
    os<<"#year    size    nIndivs     fraction mature"<<std::endl<<inpData_nc<<std::endl;
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
