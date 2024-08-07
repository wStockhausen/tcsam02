#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelData.hpp"
#include "ModelOptions.hpp"

using namespace tcsam;

//**********************************************************************
//  Includes
//      EffortData
//      CatchData
//      FleetData
//**********************************************************************
int EffortData::debug = 0;
int CatchData::debug  = 0;
int FleetData::debug  = 0;
ostream& EffortData::os = std::cout;
ostream& CatchData::os  = std::cout;
ostream& FleetData::os  = std::cout;
//----------------------------------------------------------------------
//          EffortData
//----------------------------------------------------------------------
const adstring EffortData::KW_EFFORT_DATA = "EFFORT_DATA";
/**
 * Class destructor.
 */
EffortData::~EffortData(){
    delete ptrAvgIB; ptrAvgIB=0;
}
/**
 * Set the maximum year in which to fit the data.
 * 
 * @param mxYr - the max year to include data
 */
void EffortData::setMaxYear(int mxYr){
    if (debug) {
        cout     <<"Starting EffortData::setMaxYear("<<mxYr<<")"<<endl;
        rpt::echo<<"Starting EffortData::setMaxYear("<<mxYr<<")"<<endl;
    }
    //determine number of years to keep
    int nyp = 0;
    for (int iy=1;iy<=ny;iy++) {if (yrs(iy)<=mxYr) nyp++;}
    
    //re-allocate input array
    inpEff_yc.deallocate();
    inpEff_yc.allocate(1,nyp,1,2);
    int iyp = 0;
    for (int iy=1;iy<=ny;iy++) {
        if (yrs(iy)<=mxYr) {
            iyp++;
            inpEff_yc(iyp,1) = yrs(iy);
            inpEff_yc(iyp,2) = eff_y(yrs(iy));
        }
    }
    if (iyp!=nyp){
        cout<<"Something wrong in EffortData::setMaxYear"<<endl;
        cout<<"iyp != nyp"<<endl;
        exit(0);
    }
    
    //copy back to class members
    ny = nyp;
    yrs.deallocate();
    eff_y.deallocate();
    yrs.allocate(1,ny);
    yrs = (ivector) column(inpEff_yc,1);
    int mny = min(yrs);
    int mxy = max(yrs);
    eff_y.allocate(mny,mxy); eff_y = 0.0;
    for (int iy=1;iy<=ny;iy++) eff_y(yrs(iy)) = inpEff_yc(iy,2);
    if (debug) {
        cout     <<"Finished EffortData::setMaxYear("<<mxYr<<")"<<endl;
        rpt::echo<<"Finished EffortData::setMaxYear("<<mxYr<<")"<<endl;
    }
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void EffortData::read(cifstream & is){
    if (debug){
        cout<<"start EffortData::read(...) "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading EffortData."<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    if (!(str==KW_EFFORT_DATA)){
        cout<<"#Error reading effort data from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword '"<<KW_EFFORT_DATA<<"' but got '"<<str<<"'"<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    ptrAvgIB = new IndexBlock(ModelConfiguration::mnYr,ModelConfiguration::mxYr);
    is>>(*ptrAvgIB);
    rpt::echo<<(*ptrAvgIB)<<tb<<"#intervals over which to average effort/fishing mortality"<<endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>llWgt;
    rpt::echo<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    is>>units;
    rpt::echo<<units<<tb<<"#units"<<endl;
    is>>ny;//number of years of effort data
    rpt::echo<<ny<<tb<<"#number of years"<<endl;
    inpEff_yc.allocate(1,ny,1,2);
    is>>inpEff_yc;
    rpt::echo<<"#year potlifts ("<<units<<")"<<endl<<inpEff_yc<<endl;
    
    yrs.allocate(1,ny);
    yrs = (ivector) column(inpEff_yc,1);
    int mny = min(yrs);
    int mxy = max(yrs);
    eff_y.allocate(mny,mxy); eff_y = 0.0;
    for (int iy=1;iy<=ny;iy++) eff_y(yrs(iy)) = inpEff_yc(iy,2);
    if (debug) cout<<"end EffortData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void EffortData::write(ostream & os){
    if (debug) cout<<"start EffortData::write(...) "<<this<<endl;
    os<<KW_EFFORT_DATA<<tb<<"#required keyword"<<endl;
    os<<(*ptrAvgIB)<<tb<<"#intervals over which to average effort/fishing mortality"<<endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    os<<units<<tb<<"#units"<<endl;
    os<<ny<<tb<<"#number of years"<<endl;
    os<<"#year   potlifts"<<endl<<inpEff_yc<<endl;
    if (debug) cout<<"end EffortData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void EffortData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"EffortData::writing to R"<<endl;
    adstring y  = "year=c("+wts::to_qcsv(yrs)+")";
    for (int n=0;n<indent;n++) os<<tb;
        os<<"effort=list("<<endl;
        indent++; 
            for (int n=0;n<indent;n++) os<<tb;
            os<<"avgRng="; ptrAvgIB->writeToR(os); os<<cc<<endl;
            os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc<<endl;
            os<<"llWgt="<<llWgt<<cc<<endl;
            os<<"units="<<qt<<units<<qt<<cc<<endl;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="; wts::writeToR(os,column(inpEff_yc,2),y); os<<endl;
        indent--;
    for (int n=0;n<indent;n++) os<<tb; os<<")"<<endl;
    if (debug) cout<<"EffortData::done writing to R"<<endl;
}
/////////////////////////////////end EffortData/////////////////////////
//----------------------------------------------------------------------
//          CatchData
//----------------------------------------------------------------------
const adstring CatchData::KW_CATCH_DATA = "CATCH_DATA";
/***************************************************************
*   instantiation.                                             *
***************************************************************/
CatchData::CatchData(adstring& name_){
    name = name_;
    hasN = 0;   ptrN = 0;
    hasB = 0;   ptrB = 0;
    hasZFD = 0; ptrZFD = 0;
}

/***************************************************************
*   destruction.                                               *
***************************************************************/
CatchData::~CatchData(){
    if (ptrN)   delete ptrB;   ptrN = 0;
    if (ptrB)   delete ptrB;   ptrB = 0;
    if (ptrZFD) delete ptrZFD; ptrZFD = 0;
}

/**
 * Set the maximum year in which to fit the data.
 * 
 * @param mxYr - the max year to include data
 */
void CatchData::setMaxYear(int mxYr){
    if (debug) {
        rpt::echo<<"Starting CatchData::setMaxYear("<<mxYr<<") for "<<name<<" data"<<endl;
        AggregateCatchData::debug=1;
        SizeFrequencyData::debug=1;
    }
    
    if (ptrN)   ptrN->setMaxYear(mxYr);
    if (ptrB)   ptrB->setMaxYear(mxYr);
    if (ptrZFD) ptrZFD->setMaxYear(mxYr);
    
    if (debug) {
        AggregateCatchData::debug=0;
        SizeFrequencyData::debug=0;
        rpt::echo<<"Finished CatchData::setMaxYear("<<mxYr<<") for "<<name<<" data"<<endl;
    }
}

/*************************************************\n
 * Replaces catch data based on newNatZ_yxmsz.
 * 
 * @param rng - random number generator object
 * @param ptrCDSOs - pointer to CatchDataSimOptions object
 * @param newNatZ_yxmsz - POINTER to d5_array of catch-at-size by sex/maturity/shell condition/year
 * @param wAtZ_xmz - weight-at-size by sex/maturity
 */
void CatchData::replaceCatchData(random_number_generator& rng,
                                 CatchDataSimOptions* ptrCDSOs,
                                 const d5_array& newNatZ_yxmsz, 
                                 const d3_array& wAtZ_xmz,
                                 int debug,
                                 ostream& cout){
    if (debug) cout<<"***Replacing catch data for "<<name<<endl;
    int mnY = newNatZ_yxmsz.indexmin();
    int mxY = newNatZ_yxmsz.indexmax();
    if (hasN) {
        if (debug) cout<<"replacing abundance data"<<endl;
        d4_array newN_yxms(mnY,mxY,1,nSXs,1,nMSs,1,nSCs); newN_yxms.initialize();
        for (int y=mnY;y<=mxY;y++){
            for (int x=1;x<=nSXs;x++) {
                for (int m=1;m<=nMSs;m++) {
                    for (int s=1;s<=nSCs;s++) {
                        newN_yxms(y,x,m,s) = sum((newNatZ_yxmsz)(y,x,m,s));
                    }
                }
            }                
        }
        ptrN->replaceCatchData(rng,ptrCDSOs->rngFlag,ptrCDSOs->expFacAbd,newN_yxms,debug,cout);
        if (debug) cout<<"replaced abundance data"<<endl;
    }
    if (hasB){
        if (debug) cout<<"replacing biomass data"<<endl;
        d4_array newB_yxms(mnY,mxY,1,nSXs,1,nMSs,1,nSCs); newB_yxms.initialize();
        for (int y=mnY;y<=mxY;y++){
            for (int x=1;x<=nSXs;x++){
                for (int m=1;m<=nMSs;m++){
                    for (int s=1;s<=nSCs;s++) 
                        newB_yxms(y,x,m,s) += newNatZ_yxmsz(y,x,m,s)*wAtZ_xmz(x,m);
                }
            }
        }
        ptrB->replaceCatchData(rng,ptrCDSOs->rngFlag,ptrCDSOs->expFacBio,newB_yxms,debug,cout);
        if (debug) cout<<"replaced biomass data"<<endl;
    }
    if (hasZFD){
        if (debug) cout<<"replacing n-at-size data"<<endl;
        ptrZFD->replaceSizeFrequencyData(rng,ptrCDSOs->rngFlag,ptrCDSOs->expFacZCs,newNatZ_yxmsz,debug,cout);
        if (debug) cout<<"replaced n-at-size data"<<endl;
    }
    if (debug) cout<<"finished replacing catch data for "<<name<<endl;
}

/**
 * Adds catch data based on dvar4_arrray newNatZ_xmsz.
 * 
 * @param y - year to add
 * @param newNatZ_xmsz - dvar4_array of catch numbers-at-size by sex/maturity/shell condition
 * @param wAtZ_xmz - weight-at-size by sex/maturity
 * @param cv - cv for aggregated catch sampling error
 * @param ss - sample size for size frequency sampling error
 * @param rng - random number generator
 * 
 * @return void
 */
void CatchData::addCatchData(int y, 
                            dvar4_array& newNatZ_xmsz, 
                            d3_array& wAtZ_xmz, 
                            double cv, 
                            double ss,
                            random_number_generator& rng){
    if (debug) cout<<"***starting CatchData::addCatchData("<<y<<",dvar4_array)"<<endl;
    //convert newNatZ_xmsz to d4_array and process using associated method.
    d4_array vNatZ_xmsz = wts::value(newNatZ_xmsz);
    addCatchData(y,vNatZ_xmsz,wAtZ_xmz,cv,ss,rng);
    if (debug) cout<<"***finished CatchData::addCatchData("<<y<<",dvar4_array)"<<endl;
}

/**
 * Adds catch data based on d4_arrray newNatZ_xmsz.
 * 
 * @param y - year to add
 * @param newNatZ_xmsz - d4_array of catch numbers-at-size by sex/maturity/shell condition
 * @param wAtZ_xmz - weight-at-size by sex/maturity
 * @param cv - cv for aggregated catch sampling error
 * @param ss - sample size for size frequency sampling error
 * @param rng - random number generator
 * 
 * @return void
 */
void CatchData::addCatchData(int y, 
                            d4_array& newNatZ_xmsz, 
                            d3_array& wAtZ_xmz, 
                            double cv, 
                            double ss,
                            random_number_generator& rng){
    if (debug) cout<<"***starting CatchData::addCatchData("<<y<<",d4_array)"<<endl;
    if (debug) {
        AggregateCatchData::debug=1;
        SizeFrequencyData::debug=1;
    }
    if (hasN) {
        if (debug) cout<<"adding abundance data"<<endl;
        d3_array newN_xms(1,nSXs,1,nMSs,1,nSCs); newN_xms.initialize();
        for (int x=1;x<=nSXs;x++) {
            for (int m=1;m<=nMSs;m++) {
                for (int s=1;s<=nSCs;s++) {
                    if (debug) cout<<"x,m,s = "<<x<<cc<<m<<cc<<s<<endl;
                    newN_xms(x,m,s) = sum(newNatZ_xmsz(x,m,s));
                }
            }
        }                
        ptrN->addCatchData(y,newN_xms);
        if (debug) cout<<"added abundance data"<<endl;
    }
    if (hasB){
        if (debug) cout<<"adding biomass data"<<endl;
        d3_array newB_xms(1,nSXs,1,nMSs,1,nSCs); newB_xms.initialize();
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++) 
                    newB_xms(x,m,s) += newNatZ_xmsz(x,m,s)*wAtZ_xmz(x,m);
            }
        }
        ptrB->addCatchData(y,newB_xms);
        if (debug) cout<<"added biomass data"<<endl;
    }
    if (hasZFD){
        if (debug) cout<<"adding n-at-size data"<<endl;
        ptrZFD->addSizeFrequencyData(y,newNatZ_xmsz);
        if (debug) cout<<"added n-at-size data"<<endl;
    }
    if (debug) {
        AggregateCatchData::debug=0;
        SizeFrequencyData::debug=0;
    }
    if (debug) cout<<"***finished CatchData::addCatchData("<<y<<",d4_array)"<<endl;
}

/***************************************************************
*   read.                                                      *
***************************************************************/
void CatchData::read(cifstream & is){
    if (debug){
        cout<<"start CatchData::read(...) for "<<type<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading CatchData for "<<type<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword"<<endl;
    if (!(str==KW_CATCH_DATA)){
        cout<<"#Error reading catch data from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword '"<<KW_CATCH_DATA<<"' but got '"<<str<<"'"<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    is>>str; hasN = wts::getBooleanType(str);  //has aggregate abundance data?
    is>>str; hasB = wts::getBooleanType(str);  //has aggregate biomass data?
    is>>str; hasZFD = wts::getBooleanType(str);//has size frequency data?

    rpt::echo<<wts::getBooleanType(hasN)<<tb<<"#has aggregate catch abundance (numbers) data?"<<endl;
    rpt::echo<<wts::getBooleanType(hasB)<<tb<<"#has aggregate catch biomass (weight) data?"<<endl;
    rpt::echo<<wts::getBooleanType(hasZFD)<<tb<<"#has size frequency data?"<<endl;
    rpt::echo<<"#-----------AGGREGATE CATCH ABUNDANCE (NUMBERS)---------------#"<<endl;

    
    //ABUNDANCE
    if (hasN){
        ptrN = new AggregateCatchData(name);
        rpt::echo<<"#---Reading abundance data"<<endl;
        //AggregateCatchData::debug=1;
        is>>(*ptrN);
        //AggregateCatchData::debug=0;
        rpt::echo<<"#---Read abundance data"<<endl;
    }
    
    //BIOMASS
    if (hasB){
        ptrB = new AggregateCatchData(name);
        rpt::echo<<"#---Reading biomass data"<<endl;
        //AggregateCatchData::debug=1;
        is>>(*ptrB);
        //AggregateCatchData::debug=0;
        rpt::echo<<"#---Read biomass data"<<endl;
    }
    
    //NUMBERS-AT-SIZE 
    if (hasZFD){
        ptrZFD = new SizeFrequencyData(name);
        rpt::echo<<"#---Reading size frequency data"<<endl;
        is>>(*ptrZFD);
        rpt::echo<<"#---Read size frequency data"<<endl;
    }
    if (debug) cout<<"end CatchData::read(...) "<<this<<endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void CatchData::write(ostream & os){
    if (debug) cout<<"start CatchData::write(...) "<<this<<endl;
    os<<KW_CATCH_DATA<<tb<<"#required keyword"<<endl;
    os<<wts::getBooleanType(hasN)<<tb<<"#has aggregate catch abundance (numbers) data?"<<endl;
    os<<wts::getBooleanType(hasB)<<tb<<"#has aggregate catch biomass (weight) data?"<<endl;
    os<<wts::getBooleanType(hasZFD)<<tb<<"#has size frequency data?"<<endl;
    os<<"#-----------AGGREGATE CATCH ABUNDANCE (NUMBERS)---------------#"<<endl;
    if (hasN) os<<(*ptrN)<<endl;
    os<<"#-----------AGGREGATE CATCH BIOMASS (WEIGHT)------------------#"<<endl;
    if (hasB) os<<(*ptrB)<<endl;
    os<<"#-----------NUMBERS-AT-SIZE-----------------------------------#"<<endl;
    if (hasZFD) os<<(*ptrZFD);
    if (debug) cout<<"end CatchData::write(...) "<<this<<endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void CatchData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"CatchData::writing to R"<<endl;
    for (int n=0;n<indent;n++) os<<tb;
    os<<nm<<"=list(name="<<qt<<name<<qt<<cc<<endl;
    indent++;
        //abundance
        if (hasN) {ptrN->writeToR(os,"abundance",indent); os<<cc<<endl;}
        
        //biomass
        if (hasB) {ptrB->writeToR(os,"biomass",indent); os<<cc<<endl;}
        
        //NatZ
        if (hasZFD) {ptrZFD->writeToR(os,"nAtZ",indent); os<<cc<<endl;}
    indent--;
    os<<"dummy=0)";
    if (debug) cout<<"CatchData::done writing to R"<<endl;
}
/////////////////////////////////end CatchData/////////////////////////
//----------------------------------------------------------------------
//          FleetData
//----------------------------------------------------------------------
const adstring FleetData::KW_FISHERY = "FISHERY";
const adstring FleetData::KW_SURVEY  = "SURVEY";

/**
 * Class constructor.
 */
FleetData::FleetData(){
    hasICD = 0; ptrICD = 0;
    hasRCD = 0; ptrRCD = 0;
    hasDCD = 0; ptrDCD = 0;
    hasTCD = 0; ptrTCD = 0;
    hasEff = 0; ptrEff = 0;
}
/**
 * Class destructor.
 */
FleetData::~FleetData(){
    if (ptrICD) delete ptrICD; ptrICD = 0;
    if (ptrRCD) delete ptrRCD; ptrRCD = 0;
    if (ptrDCD) delete ptrDCD; ptrDCD = 0;
    if (ptrTCD) delete ptrTCD; ptrTCD = 0;
    if (ptrEff) delete ptrEff; ptrEff = 0;
}
/**
 * Set the maximum year in which to fit the data.
 * 
 * @param mxYr - the max year to include data
 */
void FleetData::setMaxYear(int mxYr){
    if (debug) {
        rpt::echo<<"Starting FleetData::setMaxYear("<<mxYr<<") for fleet data"<<endl;
        CatchData::debug=1;
    }
        
    if (hasICD) ptrICD->setMaxYear(mxYr);
    if (hasRCD) ptrRCD->setMaxYear(mxYr);
    if (hasDCD) ptrDCD->setMaxYear(mxYr);
    if (hasTCD) ptrTCD->setMaxYear(mxYr);
    if (hasEff) ptrEff->setMaxYear(mxYr);
    
    if (debug) {
        CatchData::debug=0;
        cout<<"Finished FleetData::setMaxYear("<<mxYr<<") for fleet data"<<endl;
    }
}
/**
 * Replace existing index (survey) catch data with new values.
 * 
* @param rng - random number generator
* @param ptrSOs - pointer to SimOptions object
* @param newNatZ_yxmsz - catch data array
* @param wAtZ_xmz - weight-at-size array
* @param debug - flag to print debugging info
* @param cout - output stream for debugging info
 */
void FleetData::replaceIndexCatchData(random_number_generator& rng,
                                      SimOptions* ptrSOs,
                                      const d5_array& newNatZ_yxmsz, 
                                      const d3_array& wAtZ_xmz,
                                      int debug, 
                                      ostream& cout){
    if (debug) cout<<"starting FleetData::replaceIndexCatchData for "<<this->name<<endl;
    if (hasICD) {
        if (debug) cout<<"replacing index catch data"<<endl;
        int vi = 0;
        for (int v=1;v<=ptrSOs->nIdxCatch;v++) {
          if (debug) cout<<"v = "<<v<<tb<<"ptrSOs->ppIdxCatch[v-1]->name = "<<ptrSOs->ppIdxCatch[v-1]->name<<endl;
          if (name==ptrSOs->ppIdxCatch[v-1]->name) {
            vi=v;
            break;
          }
        }
        if (vi) {
          if (debug) cout<<"Found vi = "<<vi<<tb<<"ptrSOs->ppIdxCatch[vi-1]->name = "<<ptrSOs->ppIdxCatch[vi-1]->name<<endl;
          if (ptrICD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrICD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            if (debug) cout<<"FleetData::replaceIndexCatchData: got here 1"<<endl;
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            if (debug) cout<<"FleetData::replaceIndexCatchData: got here 2"<<endl;
            ivector bnds = wts::getBounds(newNatZ_yxmsz);
            if (debug) cout<<"bounds = "<<bnds<<endl;
            d5_array newNatZD_yxmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],bnds[7],bnds[8],1,nZBDs);
            for (int y=bnds[1];y<=bnds[2];y++){
              for (int x=bnds[3];x<=bnds[4];x++){
                for (int m=bnds[5];m<=bnds[6];m++){
                  for (int s=bnds[7];s<=bnds[8];s++){
                    if (debug) {
                      cout<<"y,x,m,s = "<<y<<tb<<x<<tb<<m<<tb<<s<<endl;
                      cout<<newNatZ_yxmsz(y,x,m,s)<<endl;
                    }
                    newNatZD_yxmsz(y,x,m,s) = mZBtoZBDs * newNatZ_yxmsz(y,x,m,s);
                    if (debug) cout<<"newNatZD_yxmsz(y,x,m,s) = "<<newNatZD_yxmsz(y,x,m,s)<<endl;
                  }
                }
              }
            }
            if (debug) cout<<"FleetData::replaceIndexCatchData: got here 3"<<endl;
            ptrICD->replaceCatchData(rng,ptrSOs->ppIdxCatch[vi-1],newNatZD_yxmsz,wAtZD_xmz,debug,cout);
            if (debug) cout<<"FleetData::replaceIndexCatchData: got here 4"<<endl;
          } else {
            ptrICD->replaceCatchData(rng,ptrSOs->ppIdxCatch[vi-1],newNatZ_yxmsz,wAtZ_xmz,debug,cout);
          }
        }
        if (debug) cout<<"replaced index catch data"<<endl;
    }
    if (debug) cout<<"finished FleetData::replaceIndexCatchData for "<<this->name<<endl;
}
/**
 * Replace existing fishery catch (retained, discarded, total) data with new values.
 * 
 * @param rng - random number generator
 * @param ptrSOs - pointer to SimOptions object
 * @param newCatZ_yxmsz - total catch data array
 * @param newRatZ_yxmsz - retained catch data array
 * @param wAtZ_xmz      - weight-at-size array
 * @param debug - flag to print debugging info
 * @param cout - output stream for debugging info
 */
void FleetData::replaceFisheryCatchData(random_number_generator& rng,
                                        SimOptions* ptrSOs,
                                        const d5_array& newCatZ_yxmsz,
                                        const d5_array& newRatZ_yxmsz,
                                        const d3_array& wAtZ_xmz,
                                        int debug, 
                                        ostream& cout){
    if (debug) cout<<"Starting FleetData::replaceFisheryCatchData for fishery "<<this->name<<endl;
    if (hasTCD) {
        if (debug) cout<<"replacing total catch data"<<endl;
        int vi = 0;
        if (debug) cout<<"Checking "<<ptrSOs->nTotCatch<<" total catch sim options"<<endl;
        for (int v=1;v<=ptrSOs->nTotCatch;v++) {
          if (debug) cout<<"v = "<<v<<tb<<"ptrSOs->ppTotCatch[v-1]->name = "<<ptrSOs->ppTotCatch[v-1]->name<<endl;
          if (name==ptrSOs->ppTotCatch[v-1]->name) {
            vi=v;
            break;
          }
        }
        if (vi) {//found a matching fleet name in ptrSOs->ppTotCatch[v]
          if (debug) cout<<"Found vi = "<<vi<<tb<<"ptrSOs->ppTotCatch[v-1]->name = "<<ptrSOs->ppTotCatch[vi-1]->name<<endl;
          if (ptrTCD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrTCD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 1"<<endl;
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 2"<<endl;
            ivector bnds = wts::getBounds(newCatZ_yxmsz);
            d5_array newCatZD_yxmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],bnds[7],bnds[8],1,nZBDs);
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 3"<<endl;
            for (int y=bnds[1];y<=bnds[2];y++){
              for (int x=bnds[3];x<=bnds[4];x++){
                for (int m=bnds[5];m<=bnds[6];m++){
                  for (int s=bnds[7];s<=bnds[8];s++){
                    if (debug) cout<<"y,x,m,s = "<<y<<tb<<x<<tb<<m<<tb<<s<<endl;
                    newCatZD_yxmsz(y,x,m,s) = mZBtoZBDs * newCatZ_yxmsz(y,x,m,s);
                    if (debug) cout<<"newCatZD_yxmsz(y,x,m,s) = "<<newCatZD_yxmsz(y,x,m,s)<<endl;
                  }
                }
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 4"<<endl;
            ptrTCD->replaceCatchData(rng,ptrSOs->ppTotCatch[vi-1],newCatZD_yxmsz,wAtZD_xmz,debug,cout);
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 5"<<endl;
          } else {
            ptrTCD->replaceCatchData(rng,ptrSOs->ppTotCatch[vi-1],newCatZ_yxmsz,wAtZ_xmz,debug,cout);
          }
        }
        if (debug) cout<<"replaced total catch data"<<endl;
    }//--has TCD
    if (hasRCD) {
        if (debug) cout<<"replacing retained catch data"<<endl;
        int vi = 0;
        for (int v=1;v<=ptrSOs->nRetCatch;v++) {
          if (debug) cout<<"v = "<<v<<tb<<"ptrSOs->ppRetCatch[v-1]->name = "<<ptrSOs->ppRetCatch[v-1]->name<<endl;
          if (name==ptrSOs->ppRetCatch[v-1]->name) {
            vi=v;
            break;
          }
        }
        if (vi) {
          if (debug) cout<<"Found vi = "<<vi<<tb<<"ptrSOs->ppRetCatch[vi-1]->name = "<<ptrSOs->ppRetCatch[vi-1]->name<<endl;
          if (ptrRCD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrRCD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData ret catch got here 1"<<endl;
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData ret catch got here 2"<<endl;
            ivector bnds = wts::getBounds(newRatZ_yxmsz);
            d5_array newRatZD_yxmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],bnds[7],bnds[8],1,nZBDs);
            if (debug) cout<<"FleetData::replaceFisheryCatchData ret catch got here 3"<<endl;
            for (int y=bnds[1];y<=bnds[2];y++){
              for (int x=bnds[3];x<=bnds[4];x++){
                for (int m=bnds[5];m<=bnds[6];m++){
                  for (int s=bnds[7];s<=bnds[8];s++){
                    if (debug) cout<<"y,x,m,s = "<<y<<tb<<x<<tb<<m<<tb<<s<<endl;
                    newRatZD_yxmsz(y,x,m,s) = mZBtoZBDs * newRatZ_yxmsz(y,x,m,s);
                    if (debug) cout<<"newRatZD_yxmsz(y,x,m,s) = "<<newRatZD_yxmsz(y,x,m,s)<<endl;
                  }
                }
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData ret catch got here 4"<<endl;
            ptrRCD->replaceCatchData(rng,ptrSOs->ppRetCatch[vi-1],newRatZD_yxmsz,wAtZD_xmz,debug,cout);
            if (debug) cout<<"FleetData::replaceFisheryCatchData ret catch got here 5"<<endl;
          } else {
            ptrRCD->replaceCatchData(rng,ptrSOs->ppRetCatch[vi-1],newRatZ_yxmsz,wAtZ_xmz,debug,cout);
          }
        }
        if (debug) cout<<"replaced retained catch data"<<endl;
    }//--has RCD
    if (hasDCD) {
        if (debug) cout<<"replacing discard catch data"<<endl;
        ivector bnds = wts::getBounds(newCatZ_yxmsz);
        d5_array newDatZ_yxmsz(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8),bnds(9),bnds(10)); 
        newDatZ_yxmsz.initialize();
        for (int y=bnds[1];y<=bnds[2];y++){
          for (int x=bnds[3];x<=bnds[4];x++){
            for (int m=bnds[5];m<=bnds[6];m++){
              for (int s=bnds[7];s<=bnds[8];s++){
                newDatZ_yxmsz(y,x,m,s) = newCatZ_yxmsz(y,x,m,s)-newRatZ_yxmsz(y,x,m,s);
              }
            }
          }
        }
        if (debug) cout<<"FleetData::replaceFisheryCatchData got here 1"<<endl;
        int vi = 0;
        for (int v=1;v<=ptrSOs->nDscCatch;v++) {
          if (debug) cout<<"v = "<<v<<tb<<"ptrSOs->ppDscCatch[v-1]->name = "<<ptrSOs->ppDscCatch[v-1]->name<<endl;
          if (name==ptrSOs->ppDscCatch[v-1]->name) {
            vi=v;
            break;
          }
        }
        if (debug) cout<<"FleetData::replaceFisheryCatchData got here 2"<<endl;
        if (vi) {
          if (debug) cout<<"Found vi = "<<vi<<tb<<"ptrSOs->ppDscCatch[vi-1]->name = "<<ptrSOs->ppDscCatch[vi-1]->name<<endl;
          if (ptrDCD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrDCD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 3"<<endl;
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 4"<<endl;
            ivector bnds = wts::getBounds(newDatZ_yxmsz);
            d5_array newDatZD_yxmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],bnds[7],bnds[8],1,nZBDs);
            for (int y=bnds[1];y<=bnds[2];y++){
              for (int x=bnds[3];x<=bnds[4];x++){
                for (int m=bnds[5];m<=bnds[6];m++){
                  for (int s=bnds[7];s<=bnds[8];s++){
                    if (debug) cout<<"y,x,m,s = "<<y<<tb<<x<<tb<<m<<tb<<s<<endl;
                    newDatZD_yxmsz(y,x,m,s) = mZBtoZBDs * newDatZ_yxmsz(y,x,m,s);
                    if (debug) cout<<"newDatZD_yxmsz(y,x,m,s) = "<<newDatZD_yxmsz(y,x,m,s)<<endl;
                  }
                }
              }
            }
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 5"<<endl;
            ptrDCD->replaceCatchData(rng,ptrSOs->ppDscCatch[vi-1],newDatZD_yxmsz,wAtZD_xmz,debug,cout);
            if (debug) cout<<"FleetData::replaceFisheryCatchData got here 6"<<endl;
          } else {
            ptrDCD->replaceCatchData(rng,ptrSOs->ppDscCatch[vi-1],newDatZ_yxmsz,wAtZ_xmz,debug,cout);
          }
        }
        if (debug) if (debug) cout<<"replaced discard catch data"<<endl;
    }//--hasDCD
    if (debug) cout<<"Finished FleetData::replaceFisheryCatchData for fishery "<<this->name<<endl;
}
/**
 * Add new year of index (survey) catch data to existing data.
 * 
 * @param y - year
 * @param newNatZ_xmsz - dvar4_array with index catch data
 * @param wAtZ_xmz - weight-at-size array
 * @param cv - cv for aggregated catch sampling error
 * @param ss - sample size for size frequency sampling error
 * @param rng - random number generator
 * 
 * @return void
 */
void FleetData::addIndexCatchData(int y, 
                                  dvar4_array& newNatZ_xmsz, 
                                  d3_array& wAtZ_xmz, 
                                  double cv, 
                                  double ss,
                                  random_number_generator& rng){
    if (debug) cout<<"starting FleetData::addIndexCatchData("<<y<<",...)"<<endl;
    if (debug) CatchData::debug = 1;
    if (hasICD) {
        if (debug) cout<<"adding index catch data for year "<<y<<endl;
          if (ptrICD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrICD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            ivector bnds = wts::getBounds(newNatZ_xmsz);
            d4_array newNatZD_xmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],1,nZBDs);
            for (int x=bnds[1];x<=bnds[2];x++){
              for (int m=bnds[3];m<=bnds[4];m++){
                for (int s=bnds[5];s<=bnds[6];s++){
                  newNatZD_xmsz(x,m,s) = mZBtoZBDs * value(newNatZ_xmsz(x,m,s));
                }
              }
            }
            ptrICD->addCatchData(y, newNatZD_xmsz, wAtZD_xmz, cv, ss, rng);
          } else {
            ptrICD->addCatchData(y, newNatZ_xmsz, wAtZ_xmz, cv, ss, rng);
          }
        if (debug) cout<<"added index catch data"<<endl;
    }
    if (debug) CatchData::debug = 0;
    if (debug) cout<<"finished FleetData::updateIndexCatchData("<<y<<",...)"<<endl;
}
/**
 * Add a new year of fishery catch (retained, discarded, total) data to existing data.
 * 
 * @param y - year
 * @param newCatZ_xmsz - total catch data array
 * @param newRatZ_xmsz - retained catch data array
 * @param wAtZ_xmz      - weight-at-size array
 * @param cv - cv for aggregated catch sampling error
 * @param ss - sample size for size frequency sampling error
 * @param rng - random number generator
 * 
 * @return void
 */
void FleetData::addFisheryCatchData(int y, 
                                    dvar4_array& newCatZ_xmsz, 
                                    dvar4_array& newRatZ_xmsz,
                                    d3_array& wAtZ_xmz, 
                                    double cv,
                                    double ss,
                                    random_number_generator& rng){
    if (debug) cout<<"starting FleetData::addFisheryCatchData("<<y<<",...)"<<endl;
    if (debug) CatchData::debug = 1;
    if (hasTCD) {
        if (debug) cout<<"adding total catch data"<<endl;
        if (ptrTCD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrTCD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            ivector bnds = wts::getBounds(newCatZ_xmsz);
            d4_array newCatZD_xmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],1,nZBDs);
            for (int x=bnds[1];x<=bnds[2];x++){
              for (int m=bnds[3];m<=bnds[4];m++){
                for (int s=bnds[5];s<=bnds[6];s++){
                  newCatZD_xmsz(x,m,s) = mZBtoZBDs * value(newCatZ_xmsz(x,m,s));
                }
              }
            }
            ptrTCD->addCatchData(y, newCatZD_xmsz, wAtZD_xmz, cv, ss, rng);
        } else {
            ptrTCD->addCatchData(y, newCatZ_xmsz, wAtZ_xmz, cv, ss, rng);
        }
        if (debug) cout<<"added total catch data"<<endl;
    }
    if (hasRCD) {
        if (debug) cout<<"adding retained catch data"<<endl;
        if (ptrRCD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrRCD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            ivector bnds = wts::getBounds(newRatZ_xmsz);
            d4_array newRatZD_xmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],1,nZBDs);
            for (int x=bnds[1];x<=bnds[2];x++){
              for (int m=bnds[3];m<=bnds[4];m++){
                for (int s=bnds[5];s<=bnds[6];s++){
                  newRatZD_xmsz(x,m,s) = mZBtoZBDs * value(newRatZ_xmsz(x,m,s));
                }
              }
            }
            ptrRCD->addCatchData(y, newRatZD_xmsz, wAtZD_xmz, cv, ss, rng);
        } else {
            ptrRCD->addCatchData(y, newRatZ_xmsz, wAtZ_xmz, cv, ss, rng);
        }
        if (debug) cout<<"added retained catch data"<<endl;
    }
    if (hasDCD) {
        if (debug) cout<<"adding discard catch data"<<endl;
        ivector bnds = wts::getBounds(newCatZ_xmsz);
        d4_array newDatZ_xmsz(bnds(1),bnds(2),bnds(3),bnds(4),bnds(5),bnds(6),bnds(7),bnds(8)); 
        newDatZ_xmsz.initialize();
        for (int x=1;x<=nSXs;x++){
            for (int m=1;m<=nMSs;m++){
                for (int s=1;s<=nSCs;s++) {
                    newDatZ_xmsz(x,m,s) = value(newCatZ_xmsz(x,m,s)-newRatZ_xmsz(x,m,s));
                }
            }
        }
        if (ptrDCD->hasZFD){
            //need to convert from model bins to data bins
            dvector zCDs = ptrDCD->ptrZFD->zCs;         //data cutpoints
            dvector zCs  = ModelConfiguration::zCutPts; //model cutpoints
            int nZBDs    = zCDs.size()-1;               //number of data size bins 
            int nZBs     = zCs.size()-1;                //number of model size bins
            dmatrix mZBtoZBDs(1,nZBDs,1,nZBs);          //matrix to convert from model to data size bins
            mZBtoZBDs.initialize();
            for (int i=1;i<=nZBs;i++){
              for (int j=1;j<=nZBDs;j++){
                if ((zCDs[j]<=zCs[i])&&(zCs[i+1]<=zCDs[j+1])) mZBtoZBDs(j,i) = 1;
              }
            }
            d3_array wAtZD_xmz(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBDs);
            wAtZD_xmz.initialize();
            for (int x=1;x<=tcsam::nSXs;x++){
              for (int m=1;m<=tcsam::nMSs;m++){
                wAtZD_xmz(x,m) = mZBtoZBDs * wAtZ_xmz(x,m);
              }
            }
            ivector bnds = wts::getBounds(newDatZ_xmsz);
            d4_array newDatZD_xmsz(bnds[1],bnds[2],bnds[3],bnds[4],bnds[5],bnds[6],1,nZBDs);
            for (int x=bnds[1];x<=bnds[2];x++){
              for (int m=bnds[3];m<=bnds[4];m++){
                for (int s=bnds[5];s<=bnds[6];s++){
                  newDatZD_xmsz(x,m,s) = mZBtoZBDs * value(newDatZ_xmsz(x,m,s));
                }
              }
            }
            ptrDCD->addCatchData(y, newDatZD_xmsz, wAtZD_xmz, cv, ss, rng);
        } else {
            ptrDCD->addCatchData(y, newDatZ_xmsz, wAtZ_xmz, cv, ss, rng);
        }
        if (debug) cout<<"added discard catch data"<<endl;
    }
    if (debug) CatchData::debug = 0;
    if (debug) cout<<"finished FleetData::addFisheryCatchData("<<y<<",...)"<<endl;
}

/***************************************************************
*   read.                                                      *
***************************************************************/
void FleetData::read(cifstream & is){
    if (debug) {
        cout<<"start FleetData::read(...) for "<<this<<endl;
        cout<<"#------------------------------------------"<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"#------------------------------------------"<<endl;
    }
    if (!is) {
        cout<<"Apparent error reading FleetData for "<<this<<endl;
        cout<<"#file name is "<<is.get_file_name()<<endl;
        cout<<"File stream is 'bad'--file may not exist!"<<endl;
        cout<<"Terminating!!"<<endl;
        exit(-1);
    }
    
    is>>type; //fleet type
    rpt::echo<<type<<tb<<"#Required keyword"<<endl;
    if ((type!=KW_FISHERY)&&(type!=KW_SURVEY)){
        cout<<"#Error reading fleet data from "<<is.get_file_name()<<endl;
        cout<<"got '"<<type<<"' but expected of the following keywords:"<<endl;
        cout<<"'"<<KW_SURVEY<<"' or '"<<KW_FISHERY<<"'"<<endl; 
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    is>>name; //fleet name
    
    adstring str;
    is>>str; hasICD = wts::getBooleanType(str);//has index (survey) catch data?
    is>>str; hasRCD = wts::getBooleanType(str);//has retained catch data?
    is>>str; hasDCD = wts::getBooleanType(str);//has discard catch data?
    is>>str; hasTCD = wts::getBooleanType(str);//has total catch data?
    is>>str; hasEff = wts::getBooleanType(str);//has effort data?
    
    rpt::echo<<name<<tb<<"#fleet source name"<<endl;
    rpt::echo<<wts::getBooleanType(hasICD)<<tb<<"#has index catch data?"<<endl;
    rpt::echo<<wts::getBooleanType(hasRCD)<<tb<<"#has retained catch data?"<<endl;
    rpt::echo<<wts::getBooleanType(hasDCD)<<tb<<"#has observed discard catch data?"<<endl;
    rpt::echo<<wts::getBooleanType(hasTCD)<<tb<<"#has observed total catch data?"<<endl;
    rpt::echo<<wts::getBooleanType(hasEff)<<tb<<"#has effort data?"<<endl;
    
    //-----------Index Catch--------------------------
    if (hasICD){
        ptrICD = new CatchData(name);
        rpt::echo<<"#---Reading index catch data for "<<name<<endl;
        is>>(*ptrICD);
        rpt::echo<<"#---Read index catch data"<<endl;
    }
    //-----------Retained Catch--------------------------
    if (hasRCD){
        ptrRCD = new CatchData(name);
        rpt::echo<<"#---Reading retained catch data for "<<name<<endl;
        is>>(*ptrRCD);
        rpt::echo<<"#---Read retained catch data"<<endl;
    }
    //-----------Discard Catch--------------------------
    if (hasDCD){
        ptrDCD = new CatchData(name);
        rpt::echo<<"#---Reading discard catch data for "<<name<<endl;
        is>>(*ptrDCD);
        rpt::echo<<"#---Read discard catch data"<<endl;
    }
    //-----------Total catch--------------------------
    if (hasTCD){
        ptrTCD = new CatchData(name);
        rpt::echo<<"#---Reading total catch data for"<<name<<endl;
        is>>(*ptrTCD);
        rpt::echo<<"#---Read total catch data"<<endl;
    }
    //-----------Effort--------------------------
    if (hasEff){
        ptrEff = new EffortData(name);
        rpt::echo<<"#---Reading effort data for "<<name<<endl;
        is>>(*ptrEff);
        rpt::echo<<"#---Read effort data"<<endl;
    }
    if (debug) cout<<"end FleetData::read(...) "<<this<<endl;
}
/**
 * Write fleet data to output stream in input file format.
 * 
 * @param os - output stream to write to
 */
void FleetData::write(ostream & os){
    if (debug) cout<<"start FleetData::write(...) "<<this<<endl;
    os<<type<<tb<<"#required keyword"<<endl;
    os<<name<<tb<<"#fleet source name"<<endl;
    os<<wts::getBooleanType(hasICD)<<tb<<"#has index catch data?"<<endl;
    os<<wts::getBooleanType(hasRCD)<<tb<<"#has retained catch data?"<<endl;
    os<<wts::getBooleanType(hasDCD)<<tb<<"#has observed discard catch data?"<<endl;
    os<<wts::getBooleanType(hasTCD)<<tb<<"#has observed total catch data?"<<endl;
    os<<wts::getBooleanType(hasEff)<<tb<<"#has effort data?"<<endl;
    os<<"#-----------Index Catch Data---------------#"<<endl;
    if (hasICD) os<<(*ptrICD);
    os<<"#-----------Retained Catch Data---------------#"<<endl;
    if (hasRCD) os<<(*ptrRCD);
    os<<"#-----------Observed Discard Catch Data----------------#"<<endl;
    if (hasDCD) os<<(*ptrDCD);
    os<<"#-----------Observed Total Catch Data------------------#"<<endl;
    if (hasTCD) os<<(*ptrTCD);
    os<<"#-----------Effort Data---------------#"<<endl;
    if (hasEff) os<<(*ptrEff);
    if (debug) cout<<"end FleetData::write(...) "<<this<<endl;
}
/**
 * Write fleet data to an output stream as an R list.
 * 
 * @param os - reference to an output stream
 * @param nm - 
 * @param indent
 */
void FleetData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) cout<<"FleetData::writing to R"<<endl;
    for (int n=0;n<indent;n++) os<<tb;
    os<<nm<<"=list(name="<<qt<<wts::replace('_',' ',name)<<qt<<cc<<endl;
    indent++;
        //index catch data
        if (hasICD) {ptrICD->writeToR(os,"index.catch",indent); os<<cc<<endl;}
        
        //retained catch data
        if (hasRCD) {ptrRCD->writeToR(os,"retained.catch",indent); os<<cc<<endl;}
        
        //observed discard catch data
        if (hasDCD) {ptrDCD->writeToR(os,"discard.catch",indent++); os<<cc<<endl;}
        
        //observed total catch data
        if (hasTCD) {ptrTCD->writeToR(os,"total.catch",indent++); os<<cc<<endl;}
    
        //effort
        if (hasEff) {ptrEff->writeToR(os,"effort",indent); os<<cc<<endl;}
    indent--;
    os<<"dummy=0)";
    if (debug) cout<<"FleetData::done writing to R"<<endl;
}
/////////////////////////////////end FleetData/////////////////////////

