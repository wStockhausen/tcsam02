#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelData.hpp"
#include "SummaryFunctions.hpp"

//using namespace tcsam;

//**********************************************************************
//  Includes
//      SizeFrequencyData
//      BioData
//      ModelDatasets
//**********************************************************************
int AggregateCatchData::debug = 0;
int SizeFrequencyData::debug  = 0;
int BioData::debug            = 0;
int ModelDatasets::debug      = 0;
//----------------------------------------------------------------------
//          AggregateCatchData
//----------------------------------------------------------------------
const adstring AggregateCatchData::KW_ABUNDANCE_DATA = "AGGREGATE_ABUNDANCE";
const adstring AggregateCatchData::KW_BIOMASS_DATA   = "AGGREGATE_BIOMASS";
    
/**
 * Aggregate catch data C_xmsy, cv_xmsy, sd_xmsy over summary indices for 
 * fitting in the objective function.
 */
void AggregateCatchData::aggregateData(void){
    if (debug) rpt::echo<<"starting AggregateCatchData::aggregateData()"<<endl;
    //get conversion factor to millions (abundance) or thousands mt (biomass)
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
 //        rpt::echo<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<std::endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
//        rpt::echo<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<std::endl;
    }
    
    int xxmn, xxmx; int mmmn, mmmx; int ssmn, ssmx;
    if (optFit==tcsam::FIT_BY_TOT){
        //aggregate over all indices
        xxmn=xxmx=tcsam::ALL_SXs;
        mmmn=mmmx=tcsam::ALL_MSs;
        ssmn=ssmx=tcsam::ALL_SCs;
    } else if ((optFit==tcsam::FIT_BY_X)||(optFit==tcsam::FIT_BY_XE)){
        //aggregate over maturity, shell condition
        xxmn=1; xxmx=tcsam::ALL_SXs;
        mmmn=mmmx=tcsam::ALL_MSs;
        ssmn=ssmx=tcsam::ALL_SCs;
    } else if ((optFit==tcsam::FIT_BY_XM)||(optFit==tcsam::FIT_BY_X_ME)||(optFit==tcsam::FIT_BY_XME)||(optFit==tcsam::FIT_BY_X_MATONLY)){
        //aggregate over shell condition
        xxmn=1; xxmx=tcsam::ALL_SXs;
        mmmn=1; mmmx=tcsam::ALL_MSs;
        ssmn=ssmx=tcsam::ALL_SCs;
    } else if ((optFit==tcsam::FIT_BY_XS)||(optFit==tcsam::FIT_BY_X_SE)){
        //aggregate over maturity
        xxmn=1; xxmx=tcsam::ALL_SXs;
        mmmn=mmmx=tcsam::ALL_MSs;
        ssmn=1; ssmx=tcsam::ALL_SCs;
    } else if ((optFit==tcsam::FIT_BY_XMS)||(optFit==tcsam::FIT_BY_XM_SE)){
        //don't aggregate
        return;
    }



    //fill in aggregate values, as necessary
    //need to be careful handling values already partially aggregated
    double tC; double vC;
    for (int x=xxmn;x<=xxmx;x++){
        int xmn, xmx;
        xmn = xmx = x; if (x==tcsam::ALL_SXs) {xmn = 1; xmx = tcsam::ALL_SXs;}
        for (int m=mmmn;m<=mmmx;m++){
            int mmn, mmx;
            mmn = mmx = m; if (m==tcsam::ALL_MSs) {mmn = 1; mmx = tcsam::ALL_MSs;}
            for (int s=ssmn;s<=ssmx;s++){
                int smn, smx;
                smn = smx = s; if (s==tcsam::ALL_SCs) {smn = 1; smx = tcsam::ALL_SCs;}
                if ((x==tcsam::ALL_SXs)||(m==tcsam::ALL_MSs)||(s==tcsam::ALL_SCs)){
                    //check to see if aggregated factor combination has been read in 
                    if (debug) rpt::echo<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(m)<<cc<<tcsam::getShellType(s)<<endl;
                    int fmn = factors.indexmin(); int fmx = factors.indexmax();
                    int qf = 0;
                    for (int f=fmn;f<=fmx;f++){
                        int xf = tcsam::getSexType(factors(f,1));
                        int mf = tcsam::getMaturityType(factors(f,2));
                        int sf = tcsam::getShellType(factors(f,3));
                        if ((x==xf)&&(m==mf)&&(s==sf)) qf++;
                    }
                    if (qf){
                        if (debug) rpt::echo<<"Skipping aggregation calculation"<<endl;
                    } else {
                        //calculate aggregate quantities
                        if (debug) rpt::echo<<"Aggregating data"<<endl;
                        for (int y=1;y<=ny;y++){
                            if (debug) rpt::echo<<y<<endl;
                            tC = 0.0;
                            vC = 0.0;
                            for (int xp=xmn;xp<=xmx;xp++){
                                for (int mp=mmn;mp<=mmx;mp++){
                                    for (int sp=smn;sp<=smx;sp++){
                                        if (debug) rpt::echo<<xp<<" "<<mp<<" "<<sp<<" "<<convFac*inpC_xmsyc(xp,mp,sp,y,2)<<" "<<inpC_xmsyc(xp,mp,sp,y,3)<<endl;
                                        tC += convFac*inpC_xmsyc(xp,mp,sp,y,2);//mean
                                        vC += square(convFac*inpC_xmsyc(xp,mp,sp,y,2)*inpC_xmsyc(xp,mp,sp,y,3));//variance = (mean*cv)^2
                                    }//sp
                                }//mp
                            }//xp
                            C_xmsy(x,m,s,y) = tC;
                            if (tC>0.0) {
                                cv_xmsy(x,m,s,y) = sqrt(vC)/tC;
                            }
                        }//y
                        if (llType==tcsam::LL_LOGNORMAL){
                            sd_xmsy(x,m,s) = sqrt(log(1.0+elem_prod(cv_xmsy(x,m,s),cv_xmsy(x,m,s))));
                        } else {
                            sd_xmsy(x,m,s) = elem_prod(cv_xmsy(x,m,s),C_xmsy(x,m,s));
                        }
                    }//qf==0
                    if (debug) {
                        rpt::echo<<"C  = "<< C_xmsy(x,m,s)<<endl;
                        rpt::echo<<"cv = "<<cv_xmsy(x,m,s)<<endl;
                        rpt::echo<<"sd = "<<sd_xmsy(x,m,s)<<endl;
                    }
                }//aggregate?
            }//s
        }//m
    }//x
    if (debug) rpt::echo<<"finished AggregateCatchData::aggregateData()"<<endl;
}

/****************************************************************
 * Replace catch data C_xmsy with new data. Units are MILLIONS for 
 * abundance data and 1000's mt for biomass data.
 * Also modifies inpC_xmsyc to reflect new data, but keeps original units.
 * Error-related quantities remain the same.
 * 
 * @param iSeed - flag to add noise to data (if !=0)
 * @param rng - random number generator
 * @param newC_yxms - d4_array of new catch data
 */
void AggregateCatchData::replaceCatchData(int iSeed,random_number_generator& rng,d4_array& newC_yxms){
    if (debug) rpt::echo<<"starting AggregateCatchData::replaceCatchData(d4_array& newC_yxms)"<<std::endl;
    //get conversion factor to millions (abundance) or thousands mt (biomass)
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
 //        rpt::echo<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<std::endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
//        rpt::echo<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<std::endl;
    }
         
    //year limits on new data
    int mnY = newC_yxms.indexmin();
    int mxY = newC_yxms.indexmax();
    
    //copy old values in temporary variables
    int oldNY = ny;
    ivector oldYrs(1,oldNY); oldYrs = yrs;
    d4_array oldSD_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY);
    d5_array oldInpC_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY,1,3);
    for (int x=1;x<=tcsam::ALL_SXs;x++){
        for (int m=1;m<=tcsam::ALL_MSs;m++){
            for (int s=1;s<=tcsam::ALL_SCs;s++) {
                oldSD_xmsy(x,m,s)    = sd_xmsy(x,m,s);
                oldInpC_xmsyc(x,m,s) = inpC_xmsyc(x,m,s);
            }
        }
    }
    if (debug) rpt::echo<<"Copied old variables"<<endl;
    
    ny = 0;//new ny [need to recalculate in case of retrospective runs]
    for (int y=1;y<=oldNY;y++){//year index for old data
        int yr = oldYrs(y);
        if ((mnY<=yr)&&(yr<=mxY)) ny++;
    }
    if (debug) rpt::echo<<"Calculated new ny"<<endl;
    
    //reallocate yrs for new number of years
    yrs.deallocate(); yrs.allocate(1,ny); yrs.initialize();
    //reallocate inpC_xmsyc for new number of years
    inpC_xmsyc.deallocate(); 
    inpC_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,3);
    inpC_xmsyc.initialize();

    //copy old data from appropriate years & factor combinations
    int yctr = 0;//year index for new data
    for (int y=1;y<=oldNY;y++){//year index for old data
        int yr = oldYrs(y);
        if ((mnY<=yr)&&(yr<=mxY)){
            yctr++;
            yrs(yctr) = yr;
            for (int i=1;i<=factors.indexmax();i++){
                int x = tcsam::getSexType(factors(i,1));
                int m = tcsam::getMaturityType(factors(i,2));
                int s = tcsam::getShellType(factors(i,3));
                double v = tcsam::extractFromYXMS(yr,x,m,s,newC_yxms);
                if (iSeed) {
                    double sd = oldSD_xmsy(x,m,s,y);
                    if (llType==tcsam::LL_LOGNORMAL){
                        v *= exp(wts::drawSampleNormal(rng,0.0,sd)-0.5*sd*sd);
                    } else {
                        v += wts::drawSampleNormal(rng,0.0,sd);
                    }
                }
                inpC_xmsyc(x,m,s,yctr,1) = oldInpC_xmsyc(x,m,s,y,1);//old year
                inpC_xmsyc(x,m,s,yctr,3) = oldInpC_xmsyc(x,m,s,y,3);//old cv
                inpC_xmsyc(x,m,s,yctr,2) = v/convFac;//new value
            }//i loop
        }//if ((mnY<=yr)&&(yr<=mxY))
    }//y loop
    if (debug) rpt::echo<<"Created new inpC_xmsyc"<<endl;
    
    //re-constitute other arrays
    C_xmsy.allocate( 1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);  C_xmsy.initialize();
    cv_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); cv_xmsy.initialize();
    sd_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); sd_xmsy.initialize();    
    int nc = factors.indexmax();
    for (int i=1;i<=nc;i++){
        int x = tcsam::getSexType(factors(i,1));
        int m = tcsam::getMaturityType(factors(i,2));
        int s = tcsam::getShellType(factors(i,3));
        C_xmsy(x,m,s)  = convFac*column(inpC_xmsyc(x,m,s),2);
        cv_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),3);
        if (llType==tcsam::LL_LOGNORMAL){
            sd_xmsy(x,m,s) = sqrt(log(1.0+elem_prod(cv_xmsy(x,m,s),cv_xmsy(x,m,s))));
        } else {
            sd_xmsy(x,m,s) = elem_prod(cv_xmsy(x,m,s),C_xmsy(x,m,s));
        }
        if (debug) {
            rpt::echo<<factors(i,1)<<tb<<factors(i,2)<<tb<<factors(i,3)<<tb<<"#factors"<<std::endl;
            rpt::echo<<"C_xmsy  = "<< C_xmsy(x,m,s)<<endl;
            rpt::echo<<"cv_xmsy = "<<cv_xmsy(x,m,s)<<endl;
            rpt::echo<<"sd_xmsy = "<<sd_xmsy(x,m,s)<<endl;
        }
    }
    
    aggregateData();
    
    if (debug) rpt::echo<<"finished AggregateCatchData::replaceCatchData(d4_array& newC_yxms)"<<std::endl;
}

/***************************************************************
*   read.                                                      *
***************************************************************/
void AggregateCatchData::read(cifstream & is){
    if (debug) {
        std::cout<<"start AggregateCatchData::read(...) "<<this<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        std::cout<<"Apparent error reading AggregateCatchData."<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        std::cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword"<<std::endl;
    if (str==KW_ABUNDANCE_DATA){type = KW_ABUNDANCE_DATA;} else
    if (str==KW_BIOMASS_DATA)  {type = KW_BIOMASS_DATA;}   else
    {   std::cout<<"#Error reading effort data from "<<is.get_file_name()<<std::endl;
        std::cout<<"Expected keyword '"<<KW_ABUNDANCE_DATA<<"' or '"<<KW_BIOMASS_DATA<<"' but got '"<<str<<"'."<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    is>>str; optFit = tcsam::getFitType(str);
    rpt::echo<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood type"<<std::endl;
    is>>ny;//number of years of catch data
    rpt::echo<<ny<<tb<<"#number of years"<<std::endl;
    is>>units;
    rpt::echo<<units<<tb<<"#units"<<std::endl;
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
        rpt::echo<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<std::endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
        rpt::echo<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<std::endl;
    }
    
    inpC_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,3);
    inpC_xmsyc.initialize();
    
    yrs.allocate(1,ny);
    C_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);   C_xmsy.initialize();
    cv_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); cv_xmsy.initialize();
    sd_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); sd_xmsy.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    rpt::echo<<nc<<tb<<"#number of factor combinations to read"<<std::endl;
    factors.allocate(1,nc,1,3);
    for (int i=1;i<=nc;i++){
        is>>factors(i);
        int x = tcsam::getSexType(factors(i,1));
        int m = tcsam::getMaturityType(factors(i,2));
        int s = tcsam::getShellType(factors(i,3));
        rpt::echo<<factors(i,1)<<tb<<factors(i,2)<<tb<<factors(i,3)<<tb<<"#factors"<<std::endl;
        if (x&&m&&s){
            is>>inpC_xmsyc(x,m,s);
            rpt::echo<<"#year    value     cv"<<std::endl;
            rpt::echo<<inpC_xmsyc(x,m,s)<<std::endl;
            yrs = (ivector) column(inpC_xmsyc(x,m,s),1);
            C_xmsy(x,m,s)  = convFac*column(inpC_xmsyc(x,m,s),2);
            cv_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),3);
            if (llType==tcsam::LL_LOGNORMAL){
                sd_xmsy(x,m,s) = sqrt(log(1.0+elem_prod(cv_xmsy(x,m,s),cv_xmsy(x,m,s))));
            } else {
                sd_xmsy(x,m,s) = elem_prod(cv_xmsy(x,m,s),C_xmsy(x,m,s));
            }
            rpt::echo<<"C_xmsy  = "<< C_xmsy(x,m,s)<<endl;
            rpt::echo<<"cv_xmsy = "<<cv_xmsy(x,m,s)<<endl;
            rpt::echo<<"sd_xmsy = "<<sd_xmsy(x,m,s)<<endl;
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"At least one factor for AggregateCatch inpC_xmsyc not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(i,1)<<tb<<factors(i,2)<<tb<<factors(i,3)<<std::endl;
            std::cout<<"Aborting..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }//nc
    
    aggregateData();
    
    if (debug) std::cout<<"end AggregateCatchData::read(...) "<<this<<std::endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void AggregateCatchData::write(ostream & os){
    if (debug) std::cout<<"start AggregateCatchData::write(...) "<<this<<std::endl;
    os<<type<<tb<<"#required keyword"<<std::endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood type"<<std::endl;
    os<<ny<<tb<<"#number of years of catch data"<<std::endl;
    os<<units<<tb<<"#units for catch data"<<std::endl;
    
    int nCs = factors.indexmax();
    os<<nCs<<tb<<"#number of sex x shell x maturity factor combinations to read in"<<std::endl;
    for (int c=1;c<=nCs;c++){
        int x = tcsam::getSexType(factors(c,1));
        int m = tcsam::getMaturityType(factors(c,2));
        int s = tcsam::getShellType(factors(c,3));
        os<<"#-------"<<factors(c,1)<<cc<<factors(c,2)<<cc<<factors(c,3)<<std::endl;
        os<<factors(c,1)<<tb<<factors(c,2)<<tb<<factors(c,3)<<std::endl;
        os<<"#year    value     cv"<<std::endl;
        os<<inpC_xmsyc(x,m,s)<<std::endl;
    }
    if (debug) std::cout<<"end AggregateCatchData::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void AggregateCatchData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) std::cout<<"AggregateCatchData::writing to R"<<std::endl;
    
//    ivector bnds = wts::getBounds(C_xmsy);
//    adstring x = tcsamDims::getSXsForR(bnds[1],bnds[2]);
//    adstring y  = "year=c("+wts::to_qcsv(yrs)+")";
    
    adstring str; adstring unitsp;
    if (type==KW_ABUNDANCE_DATA) {
        str="abundance"; 
        unitsp = tcsam::UNITS_MILLIONS;
    } else {
        str = "biomass";
        unitsp = tcsam::UNITS_KMT;
    }
    
    ivector bnds = wts::getBounds(C_xmsy);
    adstring x = tcsamDims::getSXsForR(bnds(1),bnds(2));
    adstring m = tcsamDims::getMSsForR(bnds(3),bnds(4));
    adstring s = tcsamDims::getSCsForR(bnds(5),bnds(6));
    adstring y = "y=c("+wts::to_qcsv(yrs)+")";
    
    for (int n=0;n<indent;n++) os<<tb;
        os<<str<<"=list(units="<<qt<<units<<qt<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"optFit="<<qt<<tcsam::getFitType(optFit)<<qt<<cc; 
        os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc<<std::endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"y="; wts::writeToR(os,yrs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"data="<<std::endl;
        wts::writeToR(os,C_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"cvs="<<std::endl;
        wts::writeToR(os,cv_xmsy,x,m,s,y); os<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
         os<<")";
    if (debug) std::cout<<"AggregateCatchData::done writing to R"<<std::endl;
}
/////////////////////////////////end AggregateCatchData/////////////////////////
//----------------------------------------------------------------------
//          SizeFrequencyData
//----------------------------------------------------------------------
const adstring SizeFrequencyData::KW_SIZEFREQUENCY_DATA = "SIZE_FREQUENCY_DATA";
/**
 * Normalize the size frequency data to sum to 1 over x,m,s,z.
 */
void SizeFrequencyData::normalize(void){
    dvector nT(1,ny); nT.initialize();
    for (int y=1;y<=ny;y++){
        //calculate total numbers
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++) nT(y) += sum(NatZ_xmsyz(x,m,s,y));//sum over size bins
            }
        }
        //normalize numbers-at-size by total numbers
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++) PatZ_xmsyz(x,m,s,y) = NatZ_xmsyz(x,m,s,y)/nT(y);//norm over size bins
            }
        }
    }
}

/**
 * Replace catch-at-size data NatZ_xmsyz with new data. 
 * Also modifies inpNatZ_xmsyc to reflect new data.
 * Error-related quantities remain the same.
 * 
 * @param iSeed - flag to add noise to data (if !=0)
 * @param rng - random number generator
 * @param newNatZ_yxmsz - d5_array of numbers-at-size by yxms
 */
void SizeFrequencyData::replaceSizeFrequencyData(int iSeed,random_number_generator& rng,d5_array& newNatZ_yxmsz){
    if (debug) std::cout<<"starting SizeFrequencyData::replaceSizeFrequencyData(...) "<<this<<std::endl;
    
    //year limits on new data
    int mnY = newNatZ_yxmsz.indexmin();
    int mxY = newNatZ_yxmsz.indexmax();
    
    //copy old values in temporary variables
    int oldNY = ny;
    ivector oldYrs(1,oldNY); oldYrs = yrs;
    d5_array oldInpNatZ_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY,1,2+(nZCs-1));
    for (int x=1;x<=tcsam::ALL_SXs;x++){
        for (int m=1;m<=tcsam::ALL_MSs;m++){
            for (int s=1;s<=tcsam::ALL_SCs;s++) {
                oldInpNatZ_xmsyc(x,m,s) = inpNatZ_xmsyc(x,m,s);
            }
        }
    }
    
    ny = 0;//new ny
    for (int y=1;y<=oldNY;y++){//year index for old data
        int yr = oldYrs(y);
        if ((mnY<=yr)&&(yr<=mxY)) ny++;
    }
    
    //reallocate yrs for new number of years
    yrs.deallocate(); yrs.allocate(1,ny); yrs.initialize();
    //reallocate inpNatZ_xmsyc for new number of years
    inpNatZ_xmsyc.deallocate(); 
    inpNatZ_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,2+(nZCs-1));
    inpNatZ_xmsyc.initialize();

    //copy old data from appropriate years
    int yctr = 0;//year index for new data
    for (int y=1;y<=oldNY;y++){//year index for old data
        int yr = oldYrs(y);
        if ((mnY<=yr)&&(yr<=mxY)){
            yctr++;
            yrs(yctr) = yr;
            for (int x=1;x<=tcsam::ALL_SXs;x++){
                for (int m=1;m<=tcsam::ALL_MSs;m++){
                    for (int s=1;s<=tcsam::ALL_SCs;s++) {
                        dvector n_z = tcsam::extractFromYXMSZ(yr,x,m,s,newNatZ_yxmsz);
                        if (iSeed){
                            //add stochastic component by resampling using input sample size
                            double n = sum(n_z);
                            if ((n>0)&&(ss_xmsy(x,m,s,y)>0)){
                                cout<<"org n_z = "<<n_z<<endl;
                                dvector p_z = n_z/n;
                                cout<<"p_z = "<<p_z<<endl;
                                dvector pr_z = wts::rmvlogistic(p_z,ss_xmsy(x,m,s,y),rng);//resample
                                cout<<"pr_z = "<<n_z<<endl;
                                n_z = n*pr_z; //re-scale to original size
                                cout<<"new n_z = "<<n_z<<endl;
                            }
                        }
                        inpNatZ_xmsyc(x,m,s,yctr)(1,2) = oldInpNatZ_xmsyc(x,m,s,y)(1,2);//old year, ss
                        inpNatZ_xmsyc(x,m,s,yctr)(3,2+nZCs-1).shift(1) = n_z;//new size frequency
                    }
                }
            }
        }
    }
    
    //re-constitute other arrays
    ss_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    NatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    PatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    
    ss_xmsy.initialize();
    NatZ_xmsyz.initialize();
    PatZ_xmsyz.initialize();
    
    int nc = factors.indexmax();
    for (int i=0;i<nc;i++){
        int x = tcsam::getSexType(factors(i+1,1));
        int m = tcsam::getMaturityType(factors(i+1,2));
        int s = tcsam::getShellType(factors(i+1,3));
            ss_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),2);
            for (int y=1;y<=ny;y++){
                NatZ_xmsyz(x,m,s,y)  = (inpNatZ_xmsyc(x,m,s,y)(3,2+(nZCs-1))).shift(1);
            }
    }
    normalize();
    if (debug) std::cout<<"end SizeFrequencyData::replaceSizeFrequencyData(...) "<<this<<std::endl;
}

/*******************************************************\n
*   read from input stream.\n
******************************************************/
void SizeFrequencyData::read(cifstream & is){
    if (debug) {
        std::cout<<"start SizeFrequencyData::read(...) "<<this<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        std::cout<<"Apparent error reading SizeFrequencyData."<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        std::cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword"<<std::endl;
    if (!(str==KW_SIZEFREQUENCY_DATA)){
        std::cout<<"#Error reading effort data from "<<is.get_file_name()<<std::endl;
        std::cout<<"Expected keyowrd '"<<KW_SIZEFREQUENCY_DATA<<"' but got '"<<str<<"'"<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    
    //NUMBERS-AT-SIZE 
    is>>str; optFit = tcsam::getFitType(str);
    rpt::echo<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>ny;//number of years of numbers-at-size data
    rpt::echo<<ny<<tb<<"#number of years for size frequency data"<<std::endl;
    is>>units;
    rpt::echo<<units<<tb<<"#units for numbers at size data"<<std::endl;
    is>>nZCs;
    rpt::echo<<nZCs<<tb<<"#number of size bin cut points (nZCs)"<<std::endl;
    zCs.allocate(1,nZCs);
    is>>zCs;
    rpt::echo<<zCs<<tb<<"#zCutPts (mm CW)"<<std::endl;
    zBs.allocate(1,nZCs-1);
    zBs = 0.5*(zCs(1,nZCs-1)+(--zCs(2,nZCs)));
    if (debug) std::cout<<zBs<<tb<<"#zBins"<<std::endl;
    
    yrs.allocate(1,ny);
    ss_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    NatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    PatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    inpNatZ_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,2+(nZCs-1));
    
    yrs.initialize();
    ss_xmsy.initialize();
    NatZ_xmsyz.initialize();
    PatZ_xmsyz.initialize();
    inpNatZ_xmsyc.initialize();
    
    int nc; //number of factor combinations to read in data for
    is>>nc;
    rpt::echo<<nc<<tb<<"#number of factor combinations to read"<<std::endl;
    factors.allocate(1,nc,1,3);
    for (int i=0;i<nc;i++){
        is>>factors(i+1);
        int x = tcsam::getSexType(factors(i+1,1));
        int m = tcsam::getMaturityType(factors(i+1,2));
        int s = tcsam::getShellType(factors(i+1,3));
        rpt::echo<<factors(i+1,1)<<tb<<factors(i+1,2)<<tb<<factors(i+1,3)<<tb<<"#factors"<<std::endl;
        if (x&&m&&s){
            is>>inpNatZ_xmsyc(x,m,s);
            rpt::echo<<"#year    ss     "<<zBs<<std::endl;
            rpt::echo<<inpNatZ_xmsyc(x,m,s)<<std::endl;
            yrs = (ivector) column(inpNatZ_xmsyc(x,m,s),1);
            ss_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),2);
            for (int y=1;y<=ny;y++){
                NatZ_xmsyz(x,m,s,y)  = (inpNatZ_xmsyc(x,m,s,y)(3,2+(nZCs-1))).shift(1);
//                std::cout<<"Bounds NatZ_xmsyz(x,m,s,y)   : "<<wts::getBounds(NatZ_xmsyz(x,m,s,y))<<std::endl;
//                std::cout<<"Bounds inpNatZ_xmsyc(x,m,s,y): "<<wts::getBounds(inpNatZ_xmsyc(x,m,s,y))<<std::endl;
            }
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"At least one factor for NatZ not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(i+1,1)<<tb<<factors(i+1,2)<<tb<<factors(i+1,3)<<std::endl;
            std::cout<<"Aborting..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }
    normalize();
    if (debug) std::cout<<"end SizeFrequencyData::read(...) "<<this<<std::endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void SizeFrequencyData::write(ostream & os){
    if (debug) std::cout<<"start SizeFrequencyData::write(...) "<<this<<std::endl;
    os<<KW_SIZEFREQUENCY_DATA<<tb<<"#required keyword"<<std::endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<ny<<tb<<"#number of years of size data"<<std::endl;
    os<<units<<tb<<"#units for numbers at size"<<std::endl;
    os<<nZCs<<tb<<"#number of size bin cutpoints"<<std::endl;
    os<<"#size bin cutpoints (mm CW)"<<std::endl<<zCs<<std::endl;
    
    int nCs = factors.indexmax();
    os<<nCs<<tb<<"#number of sex x shell x maturity factor combinations to read in"<<std::endl;
    for (int c=1;c<=nCs;c++){
        int x = tcsam::getSexType(factors(c,1));
        int m = tcsam::getMaturityType(factors(c,2));
        int s = tcsam::getShellType(factors(c,3));
        os<<"#-------"<<factors(c,1)<<cc<<factors(c,2)<<cc<<factors(c,3)<<std::endl;
        os<<factors(c,1)<<tb<<factors(c,2)<<tb<<factors(c,3)<<std::endl;
        os<<"#year    ss     "<<zBs<<std::endl;
        os<<inpNatZ_xmsyc(x,m,s)<<std::endl;
    }
    if (debug) std::cout<<"end SizeFrequencyData::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void SizeFrequencyData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) std::cout<<"SizeFrequencyData::writing to R"<<std::endl;
    ivector bnds = wts::getBounds(NatZ_xmsyz);
    adstring x = tcsamDims::getSXsForR(bnds(1),bnds(2));
    adstring m = tcsamDims::getMSsForR(bnds(3),bnds(4));
    adstring s = tcsamDims::getSCsForR(bnds(5),bnds(6));
    adstring y = "y=c("+wts::to_qcsv(yrs)+")";
    adstring z = "z=c("+wts::to_qcsv(zBs)+")";        
    for (int n=0;n<indent;n++) os<<tb;
        os<<"nAtZ=list(units="<<qt<<units<<qt<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"optFit="<<qt<<tcsam::getFitType(optFit)<<qt<<cc<<std::endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc<<std::endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"y="; wts::writeToR(os,yrs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"zc="; wts::writeToR(os,zCs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"sample.sizes="; wts::writeToR(os,ss_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"data="<<std::endl;
        wts::writeToR(os,NatZ_xmsyz,x,m,s,y,z); os<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
         os<<")";
    if (debug) std::cout<<"SizeFrequencyData::done writing to R"<<std::endl;
}
/////////////////////////////////end SizeFrequencyData/////////////////////////
//----------------------------------------------------------------------
//          BioData
//----------------------------------------------------------------------
const adstring BioData::KW_BIO_DATA = "BIO_DATA";
/***************************************************************
*   read.                                                      *
***************************************************************/
void BioData::read(cifstream & is){
    if (debug) std::cout<<"start BioData::read(...) "<<this<<std::endl;
    rpt::echo<<"#------------------------------------------"<<std::endl;
    rpt::echo<<"#BioData"<<std::endl;
    rpt::echo<<"#file name is "<<is.get_file_name()<<std::endl;
    rpt::echo<<"#------------------------------------------"<<std::endl;
    if (!is) {
        std::cout<<"Apparent error reading Bio Data."<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"File stream is 'bad'--file may not exist!"<<std::endl;
        std::cout<<"Terminating!!"<<std::endl;
        exit(-1);
    }
    
    adstring str;
    is>>str;
    if (!(str==KW_BIO_DATA)){
        std::cout<<"#Error reading bio data from "<<is.get_file_name()<<std::endl;
        std::cout<<"Expected keyword '"<<KW_BIO_DATA<<"' but got '"<<str<<"'"<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    rpt::echo<<str<<tb<<"#required keyword"<<std::endl;
    
    //RECRUITMENT LAG
    is>>recLag;
    rpt::echo<<recLag<<tb<<"#recLag"<<std::endl;
    
    //SIZE BINS
    is>>nZBins;
    rpt::echo<<nZBins<<tb<<"#nZBins"<<std::endl;
    zBins.allocate(1,nZBins);
    is>>zBins;
    rpt::echo<<zBins<<tb<<"#zBins (mm CW)"<<std::endl;
    
    //WEIGHT-AT-SIZE
    {is>>unitsWatZ;
    rpt::echo<<unitsWatZ<<tb<<"#unitsWatZ"<<std::endl;
    double convToKG = tcsam::getConversionMultiplier(unitsWatZ,tcsam::UNITS_KG);
    rpt::echo<<"#using conversion factor from "<<unitsWatZ<<" to kg: "<<convToKG<<std::endl;
    wAtZ_xmz.allocate(1,tcsam::nSXs,1,tcsam::nMSs,1,nZBins);
    wAtZ_xmz.initialize();
    int nc;
    is>>nc; //number of factor combinations to read in data for
    rpt::echo<<nc<<tb<<"#number of factor combinations"<<std::endl;
    adstring_array factors(1,2);
    for (int i=0;i<nc;i++){
        is>>factors;
        int sex   = tcsam::getSexType(factors(1));
        int mat   = tcsam::getMaturityType(factors(2));
        rpt::echo<<factors(1)<<tb<<factors(2)<<tb<<"#factors"<<std::endl;
        if (sex||mat){
            is>>wAtZ_xmz(sex,mat);
            rpt::echo<<wAtZ_xmz(sex,mat)<<std::endl;
            wAtZ_xmz(sex,mat) *= convToKG;//convert to kg
            if (debug) {
                std::cout<<"#wAtZ_xmz("<<factors(1)<<cc<<factors(2)<<") [kg] ="<<std::endl;
                std::cout<<"#"<<zBins<<std::endl;
                std::cout<<wAtZ_xmz(sex,mat)<<std::endl;
            }
        } else {
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"Reading file name "<<is.get_file_name()<<std::endl;
            std::cout<<"Factors for wAtZ not recognized!"<<std::endl;
            std::cout<<"Factors: "<<factors(1)<<tb<<factors(2)<<std::endl;
            std::cout<<"Terminating..."<<std::endl;
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            exit(-1);
        }
    }}
        
    //MIDPOINTS OF FISHERY SEASONS BY YEAR
    rpt::echo<<"#Timing of fisheries and mating"<<std::endl;
    is>>fshTimingTypical;
    rpt::echo<<fshTimingTypical<<tb<<"#timing of midpoint of typical fishing season"<<std::endl;
    is>>matTimingTypical;
    rpt::echo<<matTimingTypical<<tb<<"#typical mating timing"<<std::endl;
    is>>nyAtypicalT;
    rpt::echo<<nyAtypicalT<<tb<<"#number of years of atypical timing"<<std::endl;
    if (nyAtypicalT){
        timing_yc.allocate(1,nyAtypicalT,1,3);
        is>>timing_yc;
        rpt::echo<<"#atypical timing"<<std::endl<<timing_yc<<std::endl;
    }
    //now construct timing for all model years
    fshTiming_y.allocate(ModelConfiguration::mnYr,ModelConfiguration::mxYr);
    matTiming_y.allocate(ModelConfiguration::mnYr,ModelConfiguration::mxYr);
    fshTiming_y = fshTimingTypical;
    matTiming_y = matTimingTypical;
    for (int i=1;i<=nyAtypicalT;i++){
        //substitute values for atypical years
        fshTiming_y(timing_yc(i,1)) = timing_yc(i,2);
        matTiming_y(timing_yc(i,1)) = timing_yc(i,3);
    }
    
    if (debug) std::cout<<"end BioData::read(...) "<<this<<std::endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void BioData::write(ostream & os){
    if (debug) std::cout<<"start BioData::write(...) "<<this<<std::endl;
    
    os<<KW_BIO_DATA<<std::endl;
    
    os<<"#-----------RECRUITMENT LAG---------------------------#"<<std::endl;
    os<<recLag<<tb<<"#recLag (recruitment lag in years)"<<std::endl;
    
    os<<"#-----------SIZE BINS---------------------------------#"<<std::endl;
    os<<nZBins<<tb<<"#number of size bins"<<std::endl;
    os<<"#size bins (mm CW)"<<std::endl<<zBins<<std::endl;
    
    {os<<"#-----------WEIGHT-AT-SIZE----------------------------#"<<std::endl;
    os<<unitsWatZ<<tb<<"#units for weight-at-size"<<std::endl;
    os<<tcsam::nSXs*tcsam::nMSs<<tb<<"#number of factor combinations (sex x maturity state)"<<std::endl;
    adstring_array factors(1,2);
    for (int sex=1;sex<=tcsam::nSXs;sex++){
        factors(1) = tcsam::getSexType(sex);
        for (int mat=1;mat<=tcsam::nMSs;mat++){
            factors(2) = tcsam::getMaturityType(mat);
            os<<"#-------"<<factors(1)<<cc<<factors(2)<<std::endl;
            os<<factors(1)<<tb<<factors(2)<<std::endl;
            os<<wAtZ_xmz(sex,mat)<<std::endl;
        }
    }}
    
    {os<<"#-----------TIMING------------------------------------#"<<std::endl;
    os<<fshTimingTypical<<tb<<"#timing of midpoint of typical fishing season"<<std::endl;
    os<<matTimingTypical<<tb<<"#typical mating timing"<<std::endl;
    os<<nyAtypicalT<<tb<<"#number of years"<<std::endl;
    os<<"#year   midptFisheries  matingTime"<<std::endl;
    if (nyAtypicalT>0) os<<timing_yc<<std::endl;}
    
    if (debug) std::cout<<"end BioData::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void BioData::writeToR(ostream& os, string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
        os<<nm<<"=list("<<std::endl;
    indent++;
        //size bins
        for (int n=0;n<indent;n++) os<<tb;
        os<<"recLag="<<recLag<<","<<std::endl;
        
        //size bins
        for (int n=0;n<indent;n++) os<<tb;
        os<<"z=";wts::writeToR(os,zBins); os<<","<<std::endl;
        
        //weight-at-size
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"wAtZ=list(units="<<qt<<unitsWatZ<<qt<<cc<<std::endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            os<<"data="<<std::endl;
            ivector bnds = wts::getBounds(wAtZ_xmz);
            adstring x = tcsamDims::getSXsForR(bnds[1],bnds[2]);
            adstring m = tcsamDims::getMSsForR(bnds[3],bnds[4]);
            adstring z = "z=c("+wts::to_qcsv(zBins)+")";
            wts::writeToR(os,wAtZ_xmz,x,m,z); os<<std::endl;
        indent--;}
        for (int n=0;n<indent;n++) os<<tb; os<<"),"<<std::endl;
        
        //timing of fishery season midpoints and mating
        {for (int n=0;n<indent;n++) os<<tb;
        os<<"timing="<<std::endl; 
        indent++;
            for (int n=0;n<indent;n++) os<<tb;
            dmatrix tmp(1,2,ModelConfiguration::mnYr,ModelConfiguration::mxYr);
            tmp(1) = fshTiming_y;
            tmp(2) = matTiming_y;
            adstring cols = "type=c('midptFisheries','matingTime')";
            adstring yrs = "y="+str(ModelConfiguration::mnYr)+":"+str(ModelConfiguration::mxYr);
            wts::writeToR(os,trans(tmp),yrs,cols); os<<std::endl;
        indent--;}
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
}
/////////////////////////////////end BioData/////////////////////////
//----------------------------------------------------------------------
//          ModelDatasets
//----------------------------------------------------------------------
/***************************************************************
*   Instantiation                                              *
***************************************************************/
ModelDatasets::ModelDatasets(ModelConfiguration* ptrMC){
    pMC=ptrMC;
    ptrBio=0;
//    ppFsh=0;
    ppSrv=0;
}
/***************************************************************
*   Destruction                                                *
***************************************************************/
ModelDatasets::~ModelDatasets(){
    pMC=0;
    delete ptrBio;  ptrBio=0;
//    if (ppFsh) {                                                        TODO: uncomment
//        for (int f=0;f<nFsh;f++) delete ppFsh[f];
//        delete ppFsh; ppFsh = 0;
//    } 
    if (ppSrv) {
        for (int s=0;s<nSrv;s++) delete ppSrv[s];
        delete ppSrv; ppSrv = 0;
    } 
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void ModelDatasets::read(cifstream & is){
    cout<<"ModelDatasets input file name: '"<<is.get_file_name()<<"'"<<endl;
    adstring parent = wts::getParentFolder(is);
    cout<<"parent folder is '"<<parent<<endl;
    
    rpt::echo<<"#------------------------------------------"<<std::endl;
    rpt::echo<<"reading Model Datasets file."<<std::endl;
    rpt::echo<<"#file name is '"<<is.get_file_name()<<"'"<<std::endl;
    rpt::echo<<"#------------------------------------------"<<std::endl;
    is>>fnBioData;
    fnBioData = wts::concatenateFilePaths(parent,fnBioData);
    rpt::echo<<fnBioData<<tb<<"#bio data filename"<<std::endl;
    is>>nFsh;
    rpt::echo<<nFsh<<tb<<"#number of fishery datasets to read in"<<std::endl;
    if (nFsh){
        fnsFisheryData.allocate(1,nFsh);
        is>>fnsFisheryData; 
        for (int i=1;i<=nFsh;i++) {
            fnsFisheryData[i] = wts::concatenateFilePaths(parent,fnsFisheryData[i]);
            rpt::echo<<fnsFisheryData[i]<<tb<<"#fishery dataset "<<i<<std::endl;
        }
    }
    is>>nSrv;
    rpt::echo<<nSrv<<tb<<"#number of survey datasets to read in"<<std::endl;
    if (nSrv){
        fnsSurveyData.allocate(1,nSrv);
        is>>fnsSurveyData; 
        for (int i=1;i<=nSrv;i++) {
            fnsSurveyData[i] = wts::concatenateFilePaths(parent,fnsSurveyData[i]);
            rpt::echo<<fnsSurveyData[i]<<tb<<"#survey dataset "<<i<<std::endl;
        }
    }
    
    //          Bio data
    {ptrBio = new BioData();
    cifstream strm(fnBioData,ios::in);
    rpt::echo<<"#------Bio data ------------"<<std::endl;
    strm>>(*ptrBio);
    }  
    
    //          Fishery data 
    if (nFsh) {
        rpt::echo<<"#-------Fishery Datasets---------"<<std::endl;
        ppFsh = new FleetData*[nFsh];
        for (int i=0;i<nFsh;i++) {
            ppFsh[i] = new FleetData();
            cifstream strm(fnsFisheryData(i+1),ios::in);
            rpt::echo<<std::endl<<"#----------Fishery Data "<<i+1<<"-----"<<std::endl;
            strm>>(*ppFsh[i]);
        }
    }
    
    //          Survey data
    if (nSrv) {
        rpt::echo<<"#-------Survey Datasets---------"<<std::endl;
        ppSrv = new FleetData*[nSrv];
        for (int i=0;i<nSrv;i++) {
            ppSrv[i] = new FleetData();
            cifstream strm(fnsSurveyData(i+1),ios::in);
            rpt::echo<<std::endl<<"#----------Survey Data "<<i+1<<"-----"<<std::endl;
            strm>>(*ppSrv[i]);
        }
    }
}
/***************************************************************
*   read.                                                      *
***************************************************************/
void ModelDatasets::write(ostream & os){
    if (debug) std::cout<<"start ModelDatasets::write(...) "<<this<<std::endl;
    os<<fnBioData<<tb<<"#tanner crab biological data file"<<std::endl;
    os<<"#-------fishery data files---------"<<std::endl;
    os<<nFsh<<tb<<"#number of fishery data files"<<std::endl;
    for (int i=1;i<+nFsh;i++) os<<fnsFisheryData[i]<<tb<<"#fishery dataset "<<i<<std::endl;
    os<<"#-------survey data files---------"<<std::endl;
    os<<nSrv<<tb<<"#number of survey data files"<<std::endl;
    for (int i=1;i<+nSrv;i++) os<<fnsSurveyData[i]<<tb<<"#survey dataset "<<i<<std::endl;
    os<<"#-----biological data---------"<<std::endl;
    os<<(*ptrBio)<<std::endl;
    os<<"#-------fishery data ---------"<<std::endl;
    if (nFsh){
        for (int i=1;i<nFsh;i++) os<<"#----fishery dataset "<<i<<std::endl<<(*ppFsh[i-1])<<std::endl;
        os<<"#----fishery dataset "<<nFsh<<std::endl<<(*ppFsh[nFsh-1])<<std::endl;
    }
    os<<"#-------survey data ---------"<<std::endl;
    if (nSrv){
        for (int i=1;i<nSrv;i++) os<<"#----survey dataset "<<i<<std::endl<<(*ppSrv[i-1])<<std::endl;
        os<<"#----survey dataset "<<nSrv<<std::endl<<(*ppSrv[nSrv-1]);
    }
   if (debug) std::cout<<"end ModelDatasets::write(...) "<<this<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void ModelDatasets::writeToR(ostream& os, std::string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
    os<<nm<<"=list("<<std::endl;
    indent++;
        //bio data as list
            ptrBio->writeToR(os,"bio",indent); os<<cc<<std::endl;
        
        //survey data
        adstring srvNm;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"surveys=list("<<std::endl;
        indent++;
            if (ppSrv) {
                for (int i=0;i<(nSrv-1);i++) {
                    srvNm = "`"+wts::replace('_',' ',ppSrv[i]->name)+"`";
                    ppSrv[i]->writeToR(os,(char*)srvNm,indent); os<<cc<<std::endl;
                }
                srvNm = "`"+wts::replace('_',' ',ppSrv[nSrv-1]->name)+"`";
                ppSrv[nSrv-1]->writeToR(os,(char*)srvNm,indent); os<<std::endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<std::endl;            
        
        //fishery data
        adstring fshNm;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"fisheries=list("<<std::endl;
        indent++;
            if (ppFsh) {
                for (int i=0;i<(nFsh-1);i++) {
                    fshNm = "`"+wts::replace('_',' ',ppFsh[i]->name)+"`";
                    ppFsh[i]->writeToR(os,(char*)fshNm,indent); os<<cc<<std::endl;
                }
                fshNm = "`"+wts::replace('_',' ',ppFsh[nFsh-1]->name)+"`";
                ppFsh[nFsh-1]->writeToR(os,(char*)ppFsh[nFsh-1]->name,indent); os<<std::endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<")"<<std::endl;            
    indent--;
    for (int n=0;n<indent;n++) os<<tb;
    os<<")";
}
/////////////////////////////////end ModelDatasets/////////////////////////
