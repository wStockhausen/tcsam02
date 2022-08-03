#include <admodel.h>
#include <qfclib.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelOptions.hpp"
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
 * Replace existing growth data with new values based on randomization.
 * 
 * @param rng - random_number_generator object
 * @param ptrMOs - pointer to ModelOptions object
 * @param grA_xy - matrix of estimated "a" parameters for mean growth, by xy
 * @param grB_xy - matrix of estimated "b" parameters for mean growth, by xy
 * @param grBeta_xy - matrix of estimated scale ("beta") parameters for growth, by xy
 * @param zGrA_xy - dmatrix of pre-molt sizes corresponding to values in grA_xy
 * @apram zGrB_xy - dmatrix of pre-molt sizes corresponding to values in grB_xy
 * @param debug - flag to print debugging info
 * @param cout - output stream for debugging info
 * 
 * @return nothing
 * 
 * @details Modifies column 5 of inpData_xcn (i.e., the observed postmolt sizes)
 */
void GrowthData::replaceGrowthData(random_number_generator& rng,
                                    ModelOptions* ptrMOs,
                                    dvar_matrix& grA_xy, 
                                    dvar_matrix& grB_xy, 
                                    dvar_matrix& grBeta_xy,
                                    dmatrix& zGrA_xy,
                                    dmatrix& zGrB_xy,
                                    int debug, ostream& cout){
    debug=1;
    if (debug) cout<<"Starting GrowthData::replaceGrowthData(...)"<<endl;
    double rfacGrw2 = ptrMOs->ptrSimOpts->grwMultFac*ptrMOs->ptrSimOpts->grwMultFac;//squared value
    if (rfacGrw2>0){
        if (ptrMOs->ptrSimOpts->grwRngSeed>0) rng.reinitialize(ptrMOs->ptrSimOpts->grwRngSeed);
        for (int x=1;x<=nSXs;x++){
            ivector year_n = obsYears_xn(x);
            dvector zpre_n = inpData_xcn(x,2);//--premolt size
            //rgamma(alpha,beta,rng): 
            //--alpha = shape parameter = 1/CV^2
            //--beta  = rate parameter (1/scale) = 1/(mean * CV^2)
            //--mean  = alpha/beta = alpha*scale
            //--var   = alpha/beta^2 = alpha*scale^2
            //--pdf(x,alpha,beta) = (beta^alpha)/Gamma(alpha) * x^(alpha-1) * exp(-beta*x)
            //--line 494 of qfc_sim.cpp
            //to scale std dev (sigma) by factor f (sigma' = f*sigma):
            //--alpha' = alpha/f^2
            //--beta'  = beta/f^2
            dvar_vector mnZ_n;
            if (ptrMOs->optGrowthParam==0){
                mnZ_n = mfexp(grA_xy(x)(year_n)+elem_prod(grB_xy(x)(year_n),log(zpre_n)));
            } else if (ptrMOs->optGrowthParam==1){
                mnZ_n = elem_prod(
                            grA_xy(x)(year_n),
                            mfexp(
                                elem_prod(
                                    elem_div(log(elem_div( grB_xy(x)(year_n), grA_xy(x)(year_n))),
                                             log(elem_div(zGrB_xy(x)(year_n),zGrA_xy(x)(year_n)))),
                                    log(elem_div(zpre_n,zGrA_xy(x)(year_n)))
                                )
                            )
                        );
            } else if (ptrMOs->optGrowthParam==2){
                mnZ_n = elem_prod(
                            grA_xy(x)(year_n),
                            mfexp(
                                elem_prod(
                                    grB_xy(x)(year_n),
                                    log(elem_div(zpre_n,zGrA_xy(x)(year_n)))
                                )
                            )
                        );
            } else {
                //throw error
                PRINT2B1(" ")
                PRINT2B1("#---------------------")
                PRINT2B2("Unknown growth parameterization option",ptrMOs->optGrowthParam)
                PRINT2B1("Terminating model run. Please correct.")
                ad_exit(-1);
            }
            /* multiplicative scale factor, by observation */
            dvar_vector ibeta_n = 1.0/grBeta_xy(x)(year_n);
            /* location factor, by observation */
            dvar_vector alpha_n = elem_prod(mnZ_n-zpre_n,ibeta_n)/rfacGrw2;
            for (int n=1;n<=nObs_x(x);n++){
                double rg_inc = rgamma(value(alpha_n(n)),value(grBeta_xy(x,year_n(n)))/rfacGrw2,rng);
                inpData_xcn(x,3,n) = zpre_n(n) + rg_inc;//--postmolt size
            }
        }//--x
    }//--(rfacGrw2>0)
    if (debug) cout<<"Finished GrowthData::replaceGrowthData(...)"<<endl;
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
        inpData_xcn(xp,1).allocate(1,nObs_x(xp));//--year
        inpData_xcn(xp,2).allocate(1,nObs_x(xp));//--premolt size
        inpData_xcn(xp,3).allocate(1,nObs_x(xp));//--postmolt size
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
 * Simulate the maturity ogive data using randomization
 * 
 * @param rng -  random_number_generator object
 * @param vn_yxmsz - d5_array of model-predicted abundance (numbers) for the survey associated with this data, by y,x,m,s,z
 * @param debug - flag to print debugging info
 * @param cout - output stream for debugging info
 * 
 * @return nothing
 * 
 * @details Modifies column 5 of inpData_xcn and obsPrMat_n (i.e., the observed maturity ogive values).
 * 
 */
void MaturityOgiveData::replaceMaturityOgiveData(random_number_generator& rng,
                                                 SimOptions* ptrSOs,
                                                 d5_array& vn_yxmsz,
                                                 int debug, ostream& cout){
    if (debug) cout<<"Starting MaturityOgiveData::replaceMaturityOgiveData(...)"<<endl;
    double rfacMOD = ptrSOs->modDivFac;
    if (rfacMOD>0){
        if (ptrSOs->modRngSeed>0) rng.reinitialize(ptrSOs->modRngSeed);
        if (debug) cout<<"nObs= "<<nObs<<tb;
        if (nObs>0) {
            const int SEX = sex;
            /* year corresponding to observed fractions */
            const ivector y_n = obsYear_n;
           /* observation size bin corresponding to observation */
            const ivector obsZ_n = obsSize_n;
           /* observation size bin index for size corresponding to observation */
            const ivector obsIZ_n = obsZBI_n;
            /* sample sizes for observed fractions */
            const dvector ss_n = obsSS_n/rfacMOD;
            /* observed fractions of new shell mature crab at size */
            const dvector obsPM_n = obsPrMat_n;

            /* calculate model-predicted ratios by year for 
             * all observed size bins */
            int mny = vn_yxmsz.indexmin();
            int mxy = vn_yxmsz.indexmax();
            if (debug) cout<<"mny = "<<mny<<tb<<"mxy = "<<mxy<<tb<<"nZBs = "<<nZBs<<tb<<"nCutPts = "<<cutpts.size()<<endl;
            dmatrix modPrMat_yzp(mny,mxy,1,nZBs); modPrMat_yzp.initialize();
    //                dvector vmMat_z(1,nZBs);
    //                dvector vmTot_z(1,nZBs);
            for (int y=mny;y<=mxy;y++){
                dvector vmMat_z  = vn_yxmsz(y,SEX,MATURE,NEW_SHELL);
                dvector vmTot_z  = vmMat_z + vn_yxmsz(y,SEX,IMMATURE,NEW_SHELL);
                if (debug) cout<<"indexmin(vmMat_z) = "<<vmMat_z.indexmin()<<tb<<"indexmax(vmMat_z) = "<<vmMat_z.indexmax()<<endl;
                if (debug) cout<<y<<" got here"<<endl;
                if (debug) cout<<"zbRemapper indices: "<<zbRemapper.indexmin()<<tb<<zbRemapper.indexmax()<<endl;
                if (debug) cout<<"bounds(zbRemapper) = "<<wts::getBounds(zbRemapper)<<endl;
                dvector tmp = elem_div(zbRemapper*vmMat_z,
                                           zbRemapper*vmTot_z+1.0e-10);//--predicted maturity ogive (NOTE: small value added)
                if (debug) cout<<"indexmin(tmp) = "<<tmp.indexmin()<<tb<<"indexmax(tmp) = "<<tmp.indexmax()<<endl;
                modPrMat_yzp(y) = tmp;
                if (debug) cout<<y<<tb<<modPrMat_yzp(y)<<endl;
            }

            //extract model predictions corresponding to observations and simulate a random sample
            dvector modPM_n(1,nObs);  modPM_n.initialize();
            if (debug) cout<<"n"<<tb<<"y"<<tb<<"z"<<tb<<"iz"<<tb<<"ss"<<tb<<"obsPM"<<tb<<"modPM"<<endl;
            for (int n=1;n<=nObs;n++){
                if (debug) cout<<n<<tb<<y_n(n)<<tb<<obsZ_n(n)<<tb<<obsIZ_n(n)<<tb<<ss_n(n)<<tb<<obsPM_n(n)<<tb;
                    modPM_n(n) = modPrMat_yzp(y_n(n),obsIZ_n(n));//--predicted maturity ratio for size bin
                    if (debug) cout<<modPM_n(n)<<tb;
                    if ((modPM_n(n)>0.0)&(modPM_n(n)<1.0)){
                        int ss  = min(1,(int)round(ss_n(n)));//sample size
                        int tmp = 0;
                        for (int s=1;s<=ss;s++) tmp += (randu(rng)<=modPM_n(n));
                        inpData_nc(n,5) = ((double)tmp)/((double)ss);
                    } else {
                        inpData_nc(n,5) = modPM_n(n);
                    }
            }//loop over n
        }//nObs>0
    }//(rfacMOD>0)
    obsPrMat_n = column(inpData_nc,5);
    if (debug) cout<<"Finished MaturityOgiveData::replaceMaturityOgiveData(...)"<<endl;
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
