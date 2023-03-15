#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "ModelData.hpp"
#include "SummaryFunctions.hpp"
#include "ModelFunctions.hpp"

//**********************************************************************
//  Includes
//      AggregateCatchData
//      SizeFrequencyData
//      BioData
//      ModelDatasets
//**********************************************************************
int AggregateCatchData::debug = 0;
int SizeFrequencyData::debug  = 0;
int BioData::debug            = 0;
int ModelDatasets::debug      = 0;
ostream& AggregateCatchData::os = std::cout;
ostream& SizeFrequencyData::os  = std::cout;
ostream& BioData::os            = std::cout;
ostream& ModelDatasets::os      = std::cout;
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
    if (debug) os<<"starting AggregateCatchData::aggregateData() on "<<name<<endl;
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
    double tC; //summed catch value
    double vC; //summed catch variance
    double tU; //summed use flags
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
                    if (debug) os<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(m)<<cc<<tcsam::getShellType(s)<<endl;
                    int fmn = factors.indexmin(); int fmx = factors.indexmax();
                    int qf = 0;
                    for (int f=fmn;f<=fmx;f++){
                        int xf = tcsam::getSexType(factors(f,1));
                        int mf = tcsam::getMaturityType(factors(f,2));
                        int sf = tcsam::getShellType(factors(f,3));
                        if ((x==xf)&&(m==mf)&&(s==sf)) qf++;
                    }
                    if (qf){
                        if (debug) os<<"Skipping aggregation calculation"<<endl;
                    } else {
                        //calculate aggregate quantities
                        if (debug) os<<"Aggregating data"<<endl;
                        for (int y=1;y<=ny;y++){
                            if (debug) os<<y<<endl;
                            tC = 0.0;
                            vC = 0.0;
                            tU = 0.0;
                            for (int xp=xmn;xp<=xmx;xp++){
                                for (int mp=mmn;mp<=mmx;mp++){
                                    for (int sp=smn;sp<=smx;sp++){
                                        if (debug) os<<xp<<" "<<mp<<" "<<sp<<" "<<convFac*inpC_xmsyc(xp,mp,sp,y,idEV)<<" "<<inpC_xmsyc(xp,mp,sp,y,idCV)<<endl;
                                        tC += convFac*inpC_xmsyc(xp,mp,sp,y,idEV);//estimated value
                                        vC += square(convFac*inpC_xmsyc(xp,mp,sp,y,idEV)*inpC_xmsyc(xp,mp,sp,y,idCV));//variance = (mean*cv)^2
                                        tU += inpC_xmsyc(xp,mp,sp,y,idUF);//use flag
                                    }//sp
                                }//mp
                            }//xp
                            C_xmsy(x,m,s,y) = tC;
                            if (tC>0.0) cv_xmsy(x,m,s,y) = sqrt(vC)/tC;
                            if (tU>0.0) uf_xmsy(x,m,s,y) = 1.0;//aggregated use flag set to 1 unless all aggregated components had uf = 0 
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
                        rpt::echo<<"uf = "<<uf_xmsy(x,m,s)<<endl;
                    }
                }//aggregate?
            }//s
        }//m
    }//x
    if (debug) os<<"finished AggregateCatchData::aggregateData() on "<<name<<endl;
}

/****************************************************************
 * Replace catch data C_xmsy with new data. Units are MILLIONS for 
 * abundance data and 1000's mt for biomass data.
 * Also modifies inpC_xmsyc to reflect new data, but keeps original units.
 * Error-related quantities remain the same.
 * Use flags remain the same.
 * 
 * @param rng - random number generator
 * @param iSeed - flag to add noise to data (if !=0)
 * @param expFac - error expansion factor (or replacement)
  * @param newC_yxms - d4_array of new catch data
 */
void AggregateCatchData::replaceCatchData(random_number_generator& rng,int iSeed,double expFac,d4_array& newC_yxms){
    if (debug) os<<"starting AggregateCatchData::replaceCatchData(d4_array& newC_yxms) on "<<name<<std::endl;
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
    d4_array oldCV_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY);
    d5_array oldInpC_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY,1,4);
    for (int x=1;x<=tcsam::ALL_SXs;x++){
        for (int m=1;m<=tcsam::ALL_MSs;m++){
            for (int s=1;s<=tcsam::ALL_SCs;s++) {
                oldSD_xmsy(x,m,s)    = sd_xmsy(x,m,s);
                oldCV_xmsy(x,m,s)    = cv_xmsy(x,m,s);
                oldInpC_xmsyc(x,m,s) = inpC_xmsyc(x,m,s);
            }
        }
    }
    if (debug) os<<"Copied old variables"<<endl;
    
    ny = 0;//new ny [need to recalculate in case of retrospective runs]
    for (int y=1;y<=oldNY;y++){//year index for old data
        int yr = oldYrs(y);
        if ((mnY<=yr)&&(yr<=mxY)) ny++;
    }
    if (debug) os<<"Calculated new ny"<<endl;
    
    //reallocate yrs for new number of years
    yrs.deallocate(); yrs.allocate(1,ny); yrs.initialize();
    //reallocate inpC_xmsyc for new number of years
    inpC_xmsyc.deallocate(); 
    inpC_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,4);
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
                double v = tcsam::extractFromYXMS(yr,x,m,s,newC_yxms);//aggregate as necessary
                double cv = oldCV_xmsy(x,m,s,y);
                double sd = oldSD_xmsy(x,m,s,y);
                if (iSeed) {
                    cv = abs(expFac);
                    if (expFac>0) cv = expFac*oldCV_xmsy(x,m,s,y);
                    if (llType==tcsam::LL_LOGNORMAL){
                        sd = sqrt(log(1.0 + square(cv)));
                        v *= exp(wts::drawSampleNormal(rng,0.0,sd)-0.5*sd*sd);
                    } else {
                        sd = cv*v;
                        v += wts::drawSampleNormal(rng,0.0,sd);
                    }
                }
                inpC_xmsyc(x,m,s,yctr,idUF) = oldInpC_xmsyc(x,m,s,y,idUF);//old use flag
                inpC_xmsyc(x,m,s,yctr,idYr) = oldInpC_xmsyc(x,m,s,y,idYr);//old year
                inpC_xmsyc(x,m,s,yctr,idEV) = v/convFac;                  //new value
                inpC_xmsyc(x,m,s,yctr,idCV) = oldInpC_xmsyc(x,m,s,y,idCV);//old cv
            }//i loop
        }//if ((mnY<=yr)&&(yr<=mxY))
    }//y loop
    if (debug) os<<"Created new inpC_xmsyc"<<endl;
    
    //re-constitute other arrays
    C_xmsy.allocate( 1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);  C_xmsy.initialize();
    cv_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); cv_xmsy.initialize();
    sd_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); sd_xmsy.initialize();    
    uf_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); uf_xmsy.initialize();    
    int nc = factors.indexmax();
    for (int i=1;i<=nc;i++){
        int x = tcsam::getSexType(factors(i,1));
        int m = tcsam::getMaturityType(factors(i,2));
        int s = tcsam::getShellType(factors(i,3));
        C_xmsy(x,m,s)  = convFac*column(inpC_xmsyc(x,m,s),idEV);
        cv_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),idCV);
        uf_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),idUF);
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
            rpt::echo<<"uf_xmsy = "<<uf_xmsy(x,m,s)<<endl;
        }
    }
    
    aggregateData();
    
    if (debug) os<<"finished AggregateCatchData::replaceCatchData(d4_array& newC_yxms) on "<<name<<std::endl;
}

/**
 * Update catch data C_xmsy with data for new year "y". 
 * Units are MILLIONS for abundance data and 1000's mt for biomass data.
 * Also modifies inpC_xmsyc to reflect new data, but keeps original units.
 * Use flag for new year is 1.
 * 
 * @param y - year to add
 * @param newC_xms - d3_array of new catch data
 */
void AggregateCatchData::addCatchData(int y, d3_array& newC_xms){
    if (debug) os<<"starting AggregateCatchData::addCatchData(d3_array& newC_xms) on "<<name<<std::endl;
    //get conversion factor to millions (abundance) or thousands mt (biomass)
    double convFac = 1.0;
    if (type==KW_ABUNDANCE_DATA){
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_MILLIONS);
 //        os<<"#conversion factor from "<<units<<" to MILLIONS is "<<convFac<<std::endl;
    } else {
        convFac = tcsam::getConversionMultiplier(units,tcsam::UNITS_KMT);
//        os<<"#conversion factor from "<<units<<" to 1000's MT is "<<convFac<<std::endl;
    }
         
    //copy old values to temporary values
    int oldNY = ny;
    ivector oldYrs(1,oldNY); oldYrs = yrs;
    d5_array oldInpC_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY,1,4);
    for (int x=1;x<=tcsam::ALL_SXs;x++){
        for (int m=1;m<=tcsam::ALL_MSs;m++){
            for (int s=1;s<=tcsam::ALL_SCs;s++) {
                for (int iy=1;iy<=oldNY;iy++) {
                    oldInpC_xmsyc(x,m,s,iy) = inpC_xmsyc(x,m,s,iy);
                }
            }
        }
    }
    if (debug) os<<"Copied old data from old arrays"<<endl;
    
    //reallocate yrs for new number of years
    ny = ny+1;
    yrs.deallocate(); yrs.allocate(1,ny); yrs.initialize();
    yrs(1,oldNY) = oldYrs; yrs(ny) = y;
    
    //reallocate inpC_xmsyc for new number of years
    inpC_xmsyc.deallocate(); 
    inpC_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,4);
    inpC_xmsyc.initialize();

    //copy old data from appropriate years & factor combinations
    for (int x=1;x<=tcsam::ALL_SXs;x++){
        for (int m=1;m<=tcsam::ALL_MSs;m++){
            for (int s=1;s<=tcsam::ALL_SCs;s++) {
                for (int iy=1;iy<=oldNY;iy++) {
                    inpC_xmsyc(x,m,s,iy) = oldInpC_xmsyc(x,m,s,iy);
                }
            }
        }
    }
    if (debug) os<<"Copied old data to updated arrays"<<endl;
    
    //copy new data 
    if (debug) os<<"Adding new year to updated arrays"<<endl;
    for (int i=1;i<=factors.indexmax();i++){
        int x = tcsam::getSexType(factors(i,1));
        int m = tcsam::getMaturityType(factors(i,2));
        int s = tcsam::getShellType(factors(i,3));
        if (debug) os<<"x,m,s = "<<x<<cc<<m<<cc<<s<<endl;
        double v = tcsam::extractFromXMS(x,m,s,newC_xms);//aggregate as necessary
        inpC_xmsyc(x,m,s,ny,idUF) = 1;//use flag set to 1
        inpC_xmsyc(x,m,s,ny,idYr) = y;//new year
        inpC_xmsyc(x,m,s,ny,idEV) = v/convFac;//new value
        inpC_xmsyc(x,m,s,ny,idCV) = oldInpC_xmsyc(x,m,s,oldNY,4);//set new cv = old cv
    }//i loop
    if (debug) os<<"Created new inpC_xmsyc"<<endl;
    
    //re-constitute other arrays
    C_xmsy.allocate( 1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);  C_xmsy.initialize();
    cv_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); cv_xmsy.initialize();
    sd_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); sd_xmsy.initialize();    
    uf_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); uf_xmsy.initialize();    
    int nc = factors.indexmax();
    for (int i=1;i<=nc;i++){
        int x = tcsam::getSexType(factors(i,1));
        int m = tcsam::getMaturityType(factors(i,2));
        int s = tcsam::getShellType(factors(i,3));
        C_xmsy(x,m,s)  = convFac*column(inpC_xmsyc(x,m,s),idEV);
        cv_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),idCV);
        uf_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),idUF);
        if (llType==tcsam::LL_LOGNORMAL){
            sd_xmsy(x,m,s) = sqrt(log(1.0+elem_prod(cv_xmsy(x,m,s),cv_xmsy(x,m,s))));
        } else {
            sd_xmsy(x,m,s) = elem_prod(cv_xmsy(x,m,s),C_xmsy(x,m,s));
        }
        if (debug) {
            os<<factors(i,1)<<tb<<factors(i,2)<<tb<<factors(i,3)<<tb<<"#factors"<<std::endl;
            os<<"C_xmsy  = "<< C_xmsy(x,m,s)<<endl;
            os<<"cv_xmsy = "<<cv_xmsy(x,m,s)<<endl;
            os<<"sd_xmsy = "<<sd_xmsy(x,m,s)<<endl;
            os<<"uf_xmsy = "<<uf_xmsy(x,m,s)<<endl;
        }
    }
    
    aggregateData();
    
    if (debug) os<<"finished AggregateCatchData::addCatchData(d3_array& newC_xms) on "<<name<<std::endl;
}

/**
 * Set the maximum year in which to fit the data.
 * 
 * @param mxYr - the max year to include data
 */
void AggregateCatchData::setMaxYear(int mxYr){
    if (debug) {
        cout     <<"Starting AggregateCatchData::setMaxYear("<<mxYr<<") on "<<name<<endl;
        rpt::echo<<"Starting AggregateCatchData::setMaxYear("<<mxYr<<") on "<<name<<endl;
    }
    //determine number of years to keep
    int nyp = 0;
    for (int iy=1;iy<=ny;iy++) {if (yrs(iy)<=mxYr) nyp++;}
    //allocate temporary arrays
    ivector new_yrs(1,nyp);
    d5_array new_inpC_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp,1,4);
    d4_array new_C_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp);
    d4_array new_cv_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp);
    d4_array new_sd_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp);  
    d4_array new_uf_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp);  
    
    new_yrs.initialize();
    new_inpC_xmsyc.initialize();
    new_C_xmsy.initialize();
    new_cv_xmsy.initialize();
    new_sd_xmsy.initialize();
    new_uf_xmsy.initialize();
    
    int iyp = 0;
    for (int iy=1;iy<=ny;iy++) {
        if (yrs(iy)<=mxYr) {
            iyp++;
            new_yrs(iyp) = yrs(iy);
            for (int x=1;x<=tcsam::ALL_SXs;x++){
                for (int m=1;m<=tcsam::ALL_MSs;m++){
                    for (int s=1;s<=tcsam::ALL_SCs;s++){
                        new_inpC_xmsyc(x,m,s,iyp) = inpC_xmsyc(x,m,s,iy);
                        new_C_xmsy(x,m,s,iyp)     = C_xmsy(x,m,s,iy);
                        new_cv_xmsy(x,m,s,iyp)    = cv_xmsy(x,m,s,iy);
                        new_sd_xmsy(x,m,s,iyp)    = sd_xmsy(x,m,s,iy);
                        new_uf_xmsy(x,m,s,iyp)    = uf_xmsy(x,m,s,iy);
                    }//--s
                }//--m
            }//--x
         }//--if
    }//--iy
    if (iyp!=nyp){
        cout<<"Something wrong in AggregateCatchData::setMaxYear"<<endl;
        cout<<"iyp != nyp"<<endl;
        exit(0);
    }
    
    //copy back to class members
    ny = nyp;
    yrs.deallocate();
    inpC_xmsyc.deallocate();
    C_xmsy.deallocate();
    cv_xmsy.deallocate();
    sd_xmsy.deallocate();
    uf_xmsy.deallocate();
    
    yrs.allocate(new_yrs);
    inpC_xmsyc.allocate(new_inpC_xmsyc);
    C_xmsy.allocate(new_C_xmsy);
    cv_xmsy.allocate(new_cv_xmsy);
    sd_xmsy.allocate(new_sd_xmsy);
    uf_xmsy.allocate(new_sd_xmsy);
    
    for (int iy=1;iy<=ny;iy++) {
        yrs(iy) = new_yrs(iy);
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++){
                    inpC_xmsyc(x,m,s,iy) = new_inpC_xmsyc(x,m,s,iy);
                    C_xmsy(x,m,s,iy)     = new_C_xmsy(x,m,s,iy);
                    cv_xmsy(x,m,s,iy)    = new_cv_xmsy(x,m,s,iy);
                    sd_xmsy(x,m,s,iy)    = new_sd_xmsy(x,m,s,iy);
                    uf_xmsy(x,m,s,iy)    = new_uf_xmsy(x,m,s,iy);
                }//--s
            }//--m
        }//--x
    }//--iy
    
    int old_debug = debug;
    debug = 1;
    aggregateData();
    debug = old_debug;
    
    if (debug) {
        cout     <<"Finished AggregateCatchData::setMaxYear("<<mxYr<<") on "<<name<<endl;
        rpt::echo<<"Finished AggregateCatchData::setMaxYear("<<mxYr<<") on "<<name<<endl;
    }
}

/***************************************************************
*   read.                                                      *
***************************************************************/
void AggregateCatchData::read(cifstream & is){
    if (debug) {
        std::cout<<"start AggregateCatchData::read(...) on "<<name<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
        std::cout<<"#file name is "<<is.get_file_name()<<std::endl;
        std::cout<<"#------------------------------------------"<<std::endl;
    }
    if (!is) {
        std::cout<<"Apparent error reading AggregateCatchData for "<<name<<std::endl;
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
    is>>llWgt;
    rpt::echo<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
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
    
    //last index is use flag, year, catch, sd or cv
    inpC_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,4);
    inpC_xmsyc.initialize();
    
    yrs.allocate(1,ny);
    C_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);   C_xmsy.initialize();
    cv_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); cv_xmsy.initialize();
    sd_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); sd_xmsy.initialize();
    uf_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny); uf_xmsy.initialize();
    
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
            rpt::echo<<"#use? year    value     cv"<<std::endl;
            rpt::echo<<inpC_xmsyc(x,m,s)<<std::endl;
            uf_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),idUF);
            yrs = (ivector) column(inpC_xmsyc(x,m,s),idYr);
            C_xmsy(x,m,s)  = convFac*column(inpC_xmsyc(x,m,s),idEV);
            cv_xmsy(x,m,s) = column(inpC_xmsyc(x,m,s),idCV);
            if (llType==tcsam::LL_LOGNORMAL){
                sd_xmsy(x,m,s) = sqrt(log(1.0+elem_prod(cv_xmsy(x,m,s),cv_xmsy(x,m,s))));
            } else {
                sd_xmsy(x,m,s) = elem_prod(cv_xmsy(x,m,s),C_xmsy(x,m,s));
            }
            rpt::echo<<"C_xmsy  = "<< C_xmsy(x,m,s)<<endl;
            rpt::echo<<"cv_xmsy = "<<cv_xmsy(x,m,s)<<endl;
            rpt::echo<<"sd_xmsy = "<<sd_xmsy(x,m,s)<<endl;
            rpt::echo<<"uf_xmsy = "<<uf_xmsy(x,m,s)<<endl;
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
    
    if (debug) std::cout<<"end AggregateCatchData::read(...) on "<<name<<std::endl;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void AggregateCatchData::write(ostream & os){
    if (debug) std::cout<<"start AggregateCatchData::write(...) on "<<name<<std::endl;
    os<<type<<tb<<"#required keyword"<<std::endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood type"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
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
        os<<"#use? year    value     cv"<<std::endl;
        os<<inpC_xmsyc(x,m,s)<<std::endl;
    }
    if (debug) std::cout<<"end AggregateCatchData::write(...) on "<<name<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void AggregateCatchData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) std::cout<<"AggregateCatchData::writing to R on "<<name<<std::endl;
    
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
        os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc;
        os<<"llWgt="<<llWgt<<cc<<std::endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"y="; wts::writeToR(os,yrs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"data="<<std::endl;
        wts::writeToR(os,C_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"cvs="<<std::endl;
        wts::writeToR(os,cv_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"ufs="<<std::endl;
        wts::writeToR(os,uf_xmsy,x,m,s,y); os<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
         os<<")";
    if (debug) std::cout<<"AggregateCatchData::done writing to R on "<<name<<std::endl;
}
/////////////////////////////////end AggregateCatchData/////////////////////////
//----------------------------------------------------------------------
//          SizeFrequencyData
//----------------------------------------------------------------------
const adstring SizeFrequencyData::KW_SIZEFREQUENCY_DATA = "SIZE_FREQUENCY_DATA";
/**
 * Calculate re-weighting factors for iterative weighting.
 * 
 * @param newF_xms - d3_array by xms with new (incremental) weighting factors
 * @param debug - integer flag to print debugging info
 * @cout - output stream to print debugging info to
 */
void SizeFrequencyData::calcReWeightingFactors(d3_array& newF_xms,int debug,ostream& cout){
    if (debug) cout<<"#--starting SizeFrequencyData::calcReWeightingFactors(...) on "<<name<<endl;
    for (int x=newF_xms.indexmin();x<=newF_xms.indexmax();x++){
        dmatrix newF_ms = newF_xms(x);
        int mnX= x; int mxX=x;
        if (x==tcsam::ALL_SXs) {mnX= tcsam::MALE; mxX=tcsam::ALL_SXs;}
        for (int xp=mnX;xp<=mxX;xp++){
            for (int m=newF_ms.indexmin();m<=newF_ms.indexmax();m++){
                dvector newF_s = newF_ms(m);
                int mnM= m; int mxM=m;
                if (m==tcsam::ALL_MSs) {mnM= tcsam::IMMATURE; mxM=tcsam::ALL_MSs;}
                for (int mp=mnM;mp<=mxM;mp++){
                    for (int s=newF_s.indexmin();s<=newF_s.indexmax();s++){
                        int mnS= s; int mxS=s;
                        if (s==tcsam::ALL_SCs) {mnS= tcsam::NEW_SHELL; mxS=tcsam::ALL_SCs;}
                        for (int sp=mnS;sp<=mxS;sp++){
                            itrF_xms(xp,mp,sp) = newF_s(s);
                            cumF_xms(xp,mp,sp) = newF_s(s)*cumF_xms(xp,mp,sp);
                        }//sp
                    }//s
                }//mp
            }//m
        }//xp
    }//x
    if (debug>0) {
        cout<<"Iterative reweighting factors = " <<endl; wts::print(itrF_xms,cout,0); cout<<endl;
        cout<<"Cumulative reweighting factors = "<<endl; wts::print(cumF_xms,cout,0); cout<<endl;
        cout<<"#--finished SizeFrequencyData::calcReWeightingFactors(...) on "<<name<<endl;
    }
    if (debug<0) {
        ivector bnds = wts::getBounds(itrF_xms);
        adstring x = tcsamDims::getSXsForR(bnds(1),bnds(2));
        adstring m = tcsamDims::getMSsForR(bnds(3),bnds(4));
        adstring s = tcsamDims::getSCsForR(bnds(5),bnds(6));
        cout<<"list("<<endl;
        cout<<"\titFacs="; wts::writeToR(cout,itrF_xms,x,m,s); cout<<","<<endl;
        cout<<"\tcmFacs="; wts::writeToR(cout,cumF_xms,x,m,s); cout<<endl;
        cout<<")"<<endl;
        cout<<"#--finished SizeFrequencyData::calcReWeightingFactors(...) on "<<name<<endl;
    }
}

/**
 * Apply re-weighting factors for iterative weighting to 
 * input sample sizes.
 * 
 */
void SizeFrequencyData::applyReWeightingFactors(){
    if (debug) os<<"Starting SizeFrequencyData::applyReWeightingFactors() on "<<name<<endl;
    for (int x=cumF_xms.indexmin();x<=cumF_xms.indexmax();x++){
        dmatrix f_ms = cumF_xms(x);
        for (int m=f_ms.indexmin();m<=f_ms.indexmax();m++){
            dvector f_s = f_ms(m);
            for (int s=f_s.indexmin();s<=f_s.indexmax();s++){
                ss_xmsy(x,m,s) = f_s(s)*inpSS_xmsy(x,m,s);
            }//s
        }//m
    }//x
    aggregateRawNatZ();//update ss_xmsy
    if (debug) os<<"Finished SizeFrequencyData::applyReWeightingFactorse() on "<<name<<endl;
}

/**
 * Normalize the size frequency data to sum to 1 over all x,m,s,z.
 */
void SizeFrequencyData::normalize(void){
    if (debug) os<<"Starting SizeFrequencyData::normalize() on "<<name<<endl;
    
    aggPatZ_xmsyz.deallocate();
    aggPatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    aggPatZ_xmsyz.initialize();
    
    dvector nT(1,ny); nT.initialize();
    for (int y=1;y<=ny;y++){
        //calculate total numbers
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++) nT(y) += sum(aggNatZ_xmsyz(x,m,s,y));//sum over size bins
            }
        }
        //normalize numbers-at-size by total numbers
        if (nT(y)>0.0){
            for (int x=1;x<=tcsam::ALL_SXs;x++){
                for (int m=1;m<=tcsam::ALL_MSs;m++){
                    for (int s=1;s<=tcsam::ALL_SCs;s++) aggPatZ_xmsyz(x,m,s,y) = aggNatZ_xmsyz(x,m,s,y)/nT(y);//norm over size bins
                }
            }
        }
    }
    if (debug) os<<"Finished SizeFrequencyData::normalize() on "<<name<<endl;
}

/**
 * Create the indices used for tail compression based on aggregated size comps.
 * 
 * Tail compression is calculated using aggregated size comps by xmsy.
 */
void SizeFrequencyData::doTailCompression(void){
    int olddebug=debug;
    debug=1;
    if (debug) rpt::echo<<"starting SizeFrequencyData::doTailCompression() "<<this->name<<std::endl;
    
    tc_xmsyc.deallocate();
    tc_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,2);
    //tc_xmsyc.initialize(); <- body not defined in ADMB's i5arr.cpp
    //work-around:
    for (int x=1;x<=tcsam::ALL_SXs;x++) tc_xmsyc(x).initialize();
    
    tcdNatZ_xmsyz.deallocate();
    tcdNatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    tcdNatZ_xmsyz.initialize();
    
    for (int y=1;y<=ny;y++){
        
        //calculate total numbers
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++) {
                    int mnZ = aggNatZ_xmsyz(x,m,s,y).indexmin();
                    int mxZ = aggNatZ_xmsyz(x,m,s,y).indexmax();
                    tc_xmsyc(x,m,s,y)(1) = mnZ;//--default min value
                    tc_xmsyc(x,m,s,y)(2) = mxZ;//--default max value
                    double nT = sum(aggNatZ_xmsyz(x,m,s,y));//sum over size bins for xmsy
                    if (nT>0.0){
                        if (debug) {
                            rpt::echo<<yrs[y]<<tb<<tcsam::getSexType(x)<<tb<<tcsam::getMaturityType(m)<<tb<<tcsam::getShellType(s)<<tb<<"nT = "<<nT<<endl;
                            rpt::echo<<"aggNatZ_xmsyz = "<<aggNatZ_xmsyz(x,m,s,y)<<endl;
                        }
                        //find tail compression limits
                        dvector cumsum = aggNatZ_xmsyz(x,m,s,y)/nT;
                        int imnFound = 0; int imxFound = 0;
                        if (tc_limits(1)<=cumsum(1)) {tc_xmsyc(x,m,s,y)(1) = 1;   imnFound=1;}
                        if (tc_limits(2) == 0)       {tc_xmsyc(x,m,s,y)(2) = mxZ; imxFound=1;}
                        for (int i=2;i<=mxZ;i++) {
                            cumsum(i) = cumsum(i)+cumsum(i-1);//cumulative sum
                            if ((!imnFound)&&(    tc_limits(1)<=cumsum(i))) {tc_xmsyc(x,m,s,y)(1) = i; imnFound=1;}
                            if ((!imxFound)&&((1-tc_limits(2))<=cumsum(i))) {tc_xmsyc(x,m,s,y)(2) = i; imxFound=1;}
                        }
                        if (!imxFound) tc_xmsyc(x,m,s,y)(2) = mxZ;
                        if (debug) {
                            rpt::echo<<"cumsum = "<<cumsum<<endl;
                            rpt::echo<<"limits = "<<tc_xmsyc(x,m,s,y)(1)<<tb<<tc_xmsyc(x,m,s,y)(2)<<endl;
                        }
                        
                        //apply limits
                        int mnZp = tc_xmsyc(x,m,s,y)(1);
                        int mxZp = tc_xmsyc(x,m,s,y)(2);
                        tcdNatZ_xmsyz(x,m,s,y)(mnZp)          = sum(aggNatZ_xmsyz(x,m,s,y)(1,mnZp));
                        if (mnZp<mxZp) tcdNatZ_xmsyz(x,m,s,y)(mnZp+1,mxZp-1) = aggNatZ_xmsyz(x,m,s,y)(mnZp+1,mxZp-1);
                        tcdNatZ_xmsyz(x,m,s,y)(mxZp)          = sum(aggNatZ_xmsyz(x,m,s,y)(mxZp,nZCs-1));
                        
                        if (debug){
                            rpt::echo<<"yxms      = "<<yrs[y]<<tb<<x<<tb<<m<<tb<<s<<std::endl;
                            rpt::echo<<"tc_limits = "<<tc_limits<<std::endl;
                            rpt::echo<<"tc_xmysc  = "<<tc_xmsyc(x,m,s,y)<<std::endl;
                            rpt::echo<<"cumsum    = "<<cumsum<<std::endl;
                            rpt::echo<<"aggNatZ   = "<<aggNatZ_xmsyz(x,m,s,y)<<std::endl;
                            rpt::echo<<"tcdNatZ   = "<<tcdNatZ_xmsyz(x,m,s,y)<<std::endl;
                        }
                    }//--nT>0
                }//--s
            }//--m
        }//--x
    }//--y
    
    if (debug) rpt::echo<<"finished SizeFrequencyData::doTailCompression() "<<this->name<<std::endl;
    debug=olddebug;
}

/**
 * Set the DM parameter index based on previous value and current value.
 * 
 * @param x
 * @param m
 * @param s
 * @param y - year
 * @param dm - previous value of DM parameter index
 * @param inpDM - current value of DM parameter index
 * 
 * @return possibly updated value for dm.
 */
int SizeFrequencyData::setDM(int x,int m,int s,int y,int dm,int inpDM){
    if (debug>1){
        rpt::echo<<"starting SizeFrequencyData::setDM"<<endl;
        rpt::echo<<"x,m,s,y="<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(m)<<cc<<tcsam::getShellType(s)<<cc<<y<<std::endl;
        rpt::echo<<"previous dm = "<<dm<<". current dm = "<<inpDM<<std::endl;
    }
    if (inpDM) {            //want to set dm
        if (!dm){                //dm has not yet been set
            dm = inpDM;              //set dm
        } else {                 //dm has already been set, so check for consistency
            if (dm!=inpDM){          //if not consistent, terminate
                tcsam::writeMessage("Error in defining indices for Dirichlet-Multinomial parameter.");
                std::cout<<"See EchoData.dat."<<std::endl;
                rpt::echo<<"Indices must be consistent across aggregation level."<<std::endl;
                rpt::echo<<"current x,m,s,y = "<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(m)<<cc<<tcsam::getShellType(s)<<cc<<y<<std::endl;
                rpt::echo<<"previous dm = "<<dm<<". current dm = "<<inpDM<<std::endl;
                std::cout<<"Terminating..."<<std::endl;
                rpt::echo<<"Terminating..."<<std::endl;
                exit(-1);
            }
        }
    }
    if (debug>1){
        rpt::echo<<"final dm = "<<dm<<std::endl;
        rpt::echo<<"finished SizeFrequencyData::setDM"<<std::endl;
    }
    return(dm);//return possibly updated DM parameter index
}

/**
 * Calculate aggregated size compositions from raw size comps prior to tail compression.
 * 
 * Aggregation is done according to value of optFit.
 * 
 * Calculates aggNatZ_xmsyz from rawNatZ_xmsyz.
 * 
 */
void SizeFrequencyData::aggregateRawNatZ(void){
    int olddebug=debug;
    debug=1;
    if (debug) {
      rpt::echo<<"Starting aggregateRawNatZ() for "<<name<<endl;
      std::cout<<"Starting aggregateRawNatZ() for "<<name<<endl;
    }
    
    //working use flags
    uf_xmsy.deallocate();
    uf_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    uf_xmsy.initialize();
    
    //working DM indices
    dm_xmsy.deallocate();
    dm_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    dm_xmsy.initialize();
    
    //working sample sizes
    ss_xmsy.deallocate();
    ss_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    ss_xmsy.initialize();
    
    aggNatZ_xmsyz.deallocate();
    aggNatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    aggNatZ_xmsyz.initialize();
    
    int nSXs = tcsam::nSXs; int ALL_SXs = tcsam::ALL_SXs;
    int nMSs = tcsam::nMSs; int ALL_MSs = tcsam::ALL_MSs;
    int nSCs = tcsam::nSCs; int ALL_SCs = tcsam::ALL_SCs;
    
    if (debug>1){
        rpt::echo<<"aggNatZ_xmsyz.initialize():"<<endl;
        for (int iy=1;iy<=ny;iy++){
            for (int x=1;x<=ALL_SXs;x++){
                for (int m=1;m<=ALL_MSs;m++) {
                    for (int s=1;s<=ALL_SCs;s++) {
                        rpt::echo<<yrs[iy]<<tb<<tcsam::getSexType(x)<<tb<<tcsam::getMaturityType(m)<<tb<<tcsam::getShellType(s)<<aggNatZ_xmsyz(x,m,s,iy)<<endl;
                    }                    
                }
            }
        }
    }
    
    int y;
    int uf; 
    int dm;
    double ss;
    int nZBs = nZCs-1;
    dvector oP_z(1,nZBs);//observed size comp.
    
    if (debug) rpt::echo<<"optFit = "<<tcsam::getFitType(optFit)<<std::endl;
    if (optFit==tcsam::FIT_BY_TOT){
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug) rpt::echo<<"y = "<<iy<<tb<<yrs[iy]<<endl;
            uf = 0;
            dm = 0;
            ss = 0;
            oP_z.initialize();//observed size comp.
            for (int x=1;x<=ALL_SXs;x++){
                for (int m=1;m<=ALL_MSs;m++) {
                    for (int s=1;s<=ALL_SCs;s++) {
                        uf   += inpUF_xmsy(x,m,s,iy);
                        dm    = setDM(x,m,s,y,dm,inpDM_xmsy(x,m,s,iy));
                        ss   += inpSS_xmsy(x,m,s,iy);
//                        oP_z += rawNatZ_xmsyz(x,m,s,iy);
                        if (debug) rpt::echo<<"ModelConfiguration::maxZBs[x] = "<<ModelConfiguration::maxZBs[x]<<endl;
                        double mxZCm = ModelConfiguration::zCutPts[ModelConfiguration::maxZBs[x]+1];//righthand cut line
                        if (debug) rpt::echo<<"mxZCm = "<<mxZCm<<endl;
                        int mxZB = nZCs-1;//--index of largest size bin
                        for (int z=1;z<nZCs-1;z++){
                          if ((zCs[z]<=mxZCm)&&(mxZCm<zCs[z+1])) mxZB = z;
                        }
                        if (debug) rpt::echo<<"mxZB = "<<mxZB<<tb<<"mxZ = "<<zCs[mxZB+1]<<endl;
                        if (debug) rpt::echo<<yrs[iy]<<tb<<tcsam::getSexType(x)<<tb<<tcsam::getMaturityType(m)<<tb<<tcsam::getShellType(s)<<tb<<rawNatZ_xmsyz(x,m,s,iy)<<std::endl;
                        for (int z=1;   z<mxZB; z++) oP_z[z]    += rawNatZ_xmsyz(x,m,s,iy,z);
                        for (int z=mxZB;z<=nZBs;z++) oP_z[mxZB] += rawNatZ_xmsyz(x,m,s,iy,z);
                    }
                }
            }
            uf_xmsy(ALL_SXs,ALL_MSs,ALL_SCs,iy) = uf;
            dm_xmsy(ALL_SXs,ALL_MSs,ALL_SCs,iy) = dm;
            ss_xmsy(ALL_SXs,ALL_MSs,ALL_SCs,iy) = ss;
            aggNatZ_xmsyz(ALL_SXs,ALL_MSs,ALL_SCs,iy) = oP_z;
            if (debug){
                rpt::echo<<"y, x, m, s = "<<yrs[iy]<<cc<<tcsam::getSexType(ALL_SXs)<<cc<<tcsam::getMaturityType(ALL_MSs)<<cc<<tcsam::getShellType(ALL_SCs)<<std::endl;
                rpt::echo<<"uf, dm, ss = "<<uf<<cc<<dm<<cc<<ss<<std::endl;
                rpt::echo<<"oP_z                    = "<<oP_z<<endl;
                rpt::echo<<"aggNatZ_xmsyz(x,m,s,iy) = "<<aggNatZ_xmsyz(ALL_SXs,ALL_MSs,ALL_SCs,iy)<<endl;
            }
        } //loop over iy
        //FIT_BY_TOT
    } else
    if ((optFit==tcsam::FIT_BY_X)||(optFit==tcsam::FIT_BY_XE)){
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug) rpt::echo<<"y = "<<iy<<tb<<yrs[iy]<<endl;
            for (int x=1;x<=nSXs;x++) {
                uf = 0;
                dm = 0;
                ss = 0;
                oP_z.initialize();//observed size comp.
                for (int m=1;m<=ALL_MSs;m++) {
                    for (int s=1;s<=ALL_SCs;s++) {
                        uf   += inpUF_xmsy(x,m,s,iy);
                        dm    = setDM(x,m,s,y,dm,inpDM_xmsy(x,m,s,iy));
                        ss   += inpSS_xmsy(x,m,s,iy);
//                        oP_z += rawNatZ_xmsyz(x,m,s,iy);
                        if (debug) rpt::echo<<"ModelConfiguration::maxZBs[x] = "<<ModelConfiguration::maxZBs[x]<<endl;
                        double mxZCm = ModelConfiguration::zCutPts[ModelConfiguration::maxZBs[x]+1];
                        if (debug) rpt::echo<<"mxZCm = "<<mxZCm<<endl;
                        int mxZB = nZCs-1;//--index of largest size bin
                        for (int z=1;z<nZCs-1;z++){
                          if ((zCs[z]<=mxZCm)&&(mxZCm<zCs[z+1])) mxZB = z;
                        }
                        if (debug) rpt::echo<<"mxZB = "<<mxZB<<tb<<"mxZ = "<<zCs[mxZB+1]<<endl;
                        if (debug) rpt::echo<<yrs[iy]<<tb<<tcsam::getSexType(x)<<tb<<tcsam::getMaturityType(m)<<tb<<tcsam::getShellType(s)<<tb<<rawNatZ_xmsyz(x,m,s,iy)<<std::endl;
                        for (int z=1;   z<mxZB; z++) oP_z[z]    += rawNatZ_xmsyz(x,m,s,iy,z);
                        for (int z=mxZB;z<=nZBs;z++) oP_z[mxZB] += rawNatZ_xmsyz(x,m,s,iy,z);
                    }//--s
                }//--m
                uf_xmsy(x,ALL_MSs,ALL_SCs,iy) = uf;
                dm_xmsy(x,ALL_MSs,ALL_SCs,iy) = dm;
                ss_xmsy(x,ALL_MSs,ALL_SCs,iy) = ss;
                aggNatZ_xmsyz(x,ALL_MSs,ALL_SCs,iy) = oP_z;
                if (debug){
                    rpt::echo<<"y, x, m, s = "<<yrs[iy]<<cc<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(ALL_MSs)<<cc<<tcsam::getShellType(ALL_SCs)<<std::endl;
                    rpt::echo<<"uf, dm, ss = "<<uf<<cc<<dm<<cc<<ss<<std::endl;
                    rpt::echo<<"oP_z                    = "<<oP_z<<endl;
                    rpt::echo<<"aggNatZ_xmsyz(x,m,s,iy) = "<<aggNatZ_xmsyz(x,ALL_MSs,ALL_SCs,iy)<<endl;
                }
            }//x
        } //loop over iy
        //FIT_BY_X
    } else 
    if ((optFit==tcsam::FIT_BY_XM)||(optFit==tcsam::FIT_BY_X_ME)||(optFit==tcsam::FIT_BY_XME)){
        //oP_z.allocate(1,nZBs);
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug) rpt::echo<<"y = "<<iy<<tb<<yrs[iy]<<endl;
            for (int x=1;x<=nSXs;x++) {
                for (int m=1;m<=nMSs;m++){
                    uf = 0;
                    dm = 0;
                    ss = 0;
                    oP_z.initialize();//size comp. aggregated over s
                    for (int s=1;s<=ALL_SCs;s++) {
                        uf   += inpUF_xmsy(x,m,s,iy);
                        dm    = setDM(x,m,s,y,dm,inpDM_xmsy(x,m,s,iy));
                        ss += inpSS_xmsy(x,m,s,iy);
//                        oP_z += rawNatZ_xmsyz(x,m,s,iy);
                        if (debug) rpt::echo<<"ModelConfiguration::maxZBs[x] = "<<ModelConfiguration::maxZBs[x]<<endl;
                        double mxZCm = ModelConfiguration::zCutPts[ModelConfiguration::maxZBs[x]+1];
                        if (debug) rpt::echo<<"mxZCm = "<<mxZCm<<endl;
                        int mxZB = nZCs-1;//--index of largest size bin
                        for (int z=1;z<nZCs-1;z++){
                          if ((zCs[z]<=mxZCm)&&(mxZCm<zCs[z+1])) mxZB = z;
                        }
                        if (debug) rpt::echo<<"mxZB = "<<mxZB<<tb<<"mxZ = "<<zCs[mxZB+1]<<endl;
                        for (int z=1;   z<mxZB; z++) oP_z[z]    += rawNatZ_xmsyz(x,m,s,iy,z);
                        for (int z=mxZB;z<=nZBs;z++) oP_z[mxZB] += rawNatZ_xmsyz(x,m,s,iy,z);
                    }//--s
                    uf_xmsy(x,m,ALL_SCs,iy) = uf;
                    dm_xmsy(x,m,ALL_SCs,iy) = dm;
                    ss_xmsy(x,m,ALL_SCs,iy) = ss;
                    aggNatZ_xmsyz(x,m,ALL_SCs,iy) = oP_z;
                    if (debug){
                        rpt::echo<<"y, x, m, s = "<<yrs[iy]<<cc<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(m)<<cc<<tcsam::getShellType(ALL_SCs)<<std::endl;
                        rpt::echo<<"uf, dm, ss = "<<uf<<cc<<dm<<cc<<ss<<std::endl;
                        rpt::echo<<"oP_z                    = "<<oP_z<<endl;
                        rpt::echo<<"aggNatZ_xmsyz(x,m,s,iy) = "<<aggNatZ_xmsyz(x,m,ALL_SCs,iy)<<endl;
                    }
                }//m
            }//x
        } //loop over iy
        //FIT_BY_XM
    } else 
    if ((optFit==tcsam::FIT_BY_XS)||(optFit==tcsam::FIT_BY_X_SE)){
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug) rpt::echo<<"y = "<<iy<<tb<<yrs[iy]<<endl;
            for (int x=1;x<=nSXs;x++) {
                for (int s=1;s<=nSCs;s++){
                    uf = 0;
                    dm = 0;
                    ss = 0;
                    oP_z.initialize();//size comp. aggregated over m
                    for (int m=1;m<=ALL_MSs;m++) {
                        uf   += inpUF_xmsy(x,m,s,iy);
                        dm    = setDM(x,m,s,y,dm,inpDM_xmsy(x,m,s,iy));
                        ss += inpSS_xmsy(x,m,s,iy);
//                        oP_z += rawNatZ_xmsyz(x,m,s,iy);
                        if (debug) rpt::echo<<"ModelConfiguration::maxZBs[x] = "<<ModelConfiguration::maxZBs[x]<<endl;
                        double mxZCm = ModelConfiguration::zCutPts[ModelConfiguration::maxZBs[x]+1];
                        if (debug) rpt::echo<<"mxZCm = "<<mxZCm<<endl;
                        int mxZB = nZCs-1;//--index of largest size bin
                        for (int z=1;z<nZCs-1;z++){
                          if ((zCs[z]<=mxZCm)&&(mxZCm<zCs[z+1])) mxZB = z;
                        }
                        if (debug) rpt::echo<<"mxZB = "<<mxZB<<tb<<"mxZ = "<<zCs[mxZB+1]<<endl;
                        for (int z=1;   z<mxZB; z++) oP_z[z]    += rawNatZ_xmsyz(x,m,s,iy,z);
                        for (int z=mxZB;z<=nZBs;z++) oP_z[mxZB] += rawNatZ_xmsyz(x,m,s,iy,z);
                    }//--m
                    uf_xmsy(x,ALL_MSs,s,iy) = uf;
                    dm_xmsy(x,ALL_MSs,s,iy) = dm;
                    ss_xmsy(x,ALL_MSs,s,iy) = ss;
                    aggNatZ_xmsyz(x,ALL_MSs,s,iy) = oP_z;
                    if (debug){
                        rpt::echo<<"x, m, s = "<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(ALL_MSs)<<cc<<tcsam::getShellType(s)<<std::endl;
                        rpt::echo<<"uf, dm, ss = "<<uf<<cc<<dm<<cc<<ss<<std::endl;
                        rpt::echo<<"aggNatZ_xmsyz(x,m,s,iy) = "<<oP_z<<endl;
                    }
                }//s
            }//x
        } //loop over iy
        //FIT_BY_XS
    } else 
    if ((optFit==tcsam::FIT_BY_XMS)||(optFit==tcsam::FIT_BY_XM_SE)||(optFit==tcsam::FIT_BY_X_MSE)){
        for (int iy=1;iy<=yrs.size();iy++) {
            y = yrs[iy];
            if (debug) rpt::echo<<"y = "<<iy<<tb<<yrs[iy]<<endl;
                for (int x=1;x<=nSXs;x++) {
                    for (int m=1;m<=nMSs;m++){
                        for (int s=1;s<=nSCs;s++) {
                            uf = 0;
                            dm = 0;
                            ss = 0;
                            oP_z.initialize();//observed size comp.
                            uf   += inpUF_xmsy(x,m,s,iy);
                            dm    = setDM(x,m,s,y,dm,inpDM_xmsy(x,m,s,iy));
                            ss += inpSS_xmsy(x,m,s,iy);
//                            oP_z += rawNatZ_xmsyz(x,m,s,iy);
                            if (debug) rpt::echo<<"ModelConfiguration::maxZBs[x] = "<<ModelConfiguration::maxZBs[x]<<endl;
                            double mxZCm = ModelConfiguration::zCutPts[ModelConfiguration::maxZBs[x]+1];
                            if (debug) rpt::echo<<"mxZCm = "<<mxZCm<<endl;
                            int mxZB = nZCs-1;//--index of largest size bin
                            for (int z=1;z<nZCs-1;z++){
                              if ((zCs[z]<=mxZCm)&&(mxZCm<zCs[z+1])) mxZB = z;
                            }
                            if (debug) rpt::echo<<"mxZB = "<<mxZB<<tb<<"mxZ = "<<zCs[mxZB+1]<<endl;
                            for (int z=1;   z<mxZB; z++) oP_z[z]    += rawNatZ_xmsyz(x,m,s,iy,z);
                            for (int z=mxZB;z<=nZBs;z++) oP_z[mxZB] += rawNatZ_xmsyz(x,m,s,iy,z);
                            uf_xmsy(x,m,s,iy) = uf;
                            dm_xmsy(x,m,s,iy) = dm;
                            ss_xmsy(x,m,s,iy) = ss;
                            aggNatZ_xmsyz(x,m,s,iy) = oP_z;
                            if (debug){
                                rpt::echo<<"x, m, s = "<<tcsam::getSexType(x)<<cc<<tcsam::getMaturityType(m)<<cc<<tcsam::getShellType(s)<<std::endl;
                                rpt::echo<<"uf, dm, ss = "<<uf<<cc<<dm<<cc<<ss<<std::endl;
                                rpt::echo<<"aggNatZ_xmsyz(x,m,s,iy) = "<<oP_z<<endl;
                            }
                        }//s
                    }//m
               }//x
        } //loop over iy
        //FIT_BY_XMS
    } else 
    {
        std::cout<<"Calling calcAggregatedNatZ with invalid fit option."<<endl;
        std::cout<<"Invalid fit option was '"<<tcsam::getFitType(optFit)<<qt<<endl;
        std::cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug) {
      rpt::echo<<"Finished SizeFrequencyData::aggregateRawNatZ() for"<<name<<endl;
      std::cout<<"Finished SizeFrequencyData::aggregateRawNatZ() for"<<name<<endl;
    }
    debug=olddebug;
}

/**
 * Replace catch-at-size data NatZ_xmsyz with new data. 
 * Also modifies inpNatZ_xmsyc to reflect new data.
 * Error-related quantities remain the same.
 * 
 * @param rng - random number generator
 * @param iSeed - flag to add noise to data (if !=0)
 * @param expFac - error expansion factor
 * @param newNatZ_yxmsz - d5_array of numbers-at-size by yxmsz
 */
void SizeFrequencyData::replaceSizeFrequencyData(random_number_generator& rng,int iSeed,double expFac,const d5_array& newNatZ_yxmsz){
    if (debug) std::cout<<"starting SizeFrequencyData::replaceSizeFrequencyData(...) for "<<name<<std::endl;
    
    //year limits on new data
    int mnY = newNatZ_yxmsz.indexmin();
    int mxY = newNatZ_yxmsz.indexmax();
    
    //copy old values in temporary variables
    int oldNY = ny;
    ivector oldYrs(1,oldNY); oldYrs = yrs;
    d5_array oldInpNatZ_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY,1,LAST_COL+(nZCs-1));
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
    inpNatZ_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,LAST_COL+(nZCs-1));
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
                        double ss = oldInpNatZ_xmsyc(x,m,s,y,4);//original sample size
                        if (iSeed){
                            ss = abs(expFac);//fix sample size
                            if (expFac>0) ss = expFac*oldInpNatZ_xmsyc(x,m,s,y,4);//expand input sample size 
                            //add stochastic component by resampling using input sample size
                            double n = sum(n_z);
                            if ((n>0)&&(inpSS_xmsy(x,m,s,y)>0)){
                                cout<<"org n_z = "<<n_z<<endl;
                                dvector p_z = n_z/n;
                                cout<<"p_z = "<<p_z<<endl;
                                dvector pr_z = wts::rmvlogistic(p_z,inpSS_xmsy(x,m,s,y),rng);//resample
                                cout<<"pr_z = "<<n_z<<endl;
                                n_z = n*pr_z; //re-scale to original size
                                cout<<"new n_z = "<<n_z<<endl;
                            }
                        }
                        inpNatZ_xmsyc(x,m,s,yctr)(1,LAST_COL-1) = oldInpNatZ_xmsyc(x,m,s,y)(1,LAST_COL-1);//old use flag, DM, year
                        inpNatZ_xmsyc(x,m,s,yctr)(LAST_COL)     = ss; //new sample size
                        inpNatZ_xmsyc(x,m,s,yctr)(LAST_COL+1,LAST_COL+nZCs-1).shift(1) = n_z;//new size frequency
                    }
                }
            }
        }
    }
    
    //re-constitute other arrays
    inpUF_xmsy.deallocate();
    inpDM_xmsy.deallocate();
    inpSS_xmsy.deallocate();
    rawNatZ_xmsyz.deallocate();
    
    inpUF_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpDM_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpSS_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    rawNatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    
    inpUF_xmsy.initialize();
    inpDM_xmsy.initialize();
    inpSS_xmsy.initialize();
    rawNatZ_xmsyz.initialize();
    
    //reset re-weighting multipliers
    cumF_xms.deallocate();
    cumF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    cumF_xms.initialize();
    cumF_xms = 1.0;
    
    itrF_xms.deallocate();
    itrF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    itrF_xms.initialize();
    
    int nc = factors.indexmax();
    for (int i=0;i<nc;i++){
        int x = tcsam::getSexType(factors(i+1,1));
        int m = tcsam::getMaturityType(factors(i+1,2));
        int s = tcsam::getShellType(factors(i+1,3));
        int k=0;
        inpUF_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),1);
        inpDM_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),2);
        inpSS_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),4);
        for (int y=1;y<=ny;y++){
            rawNatZ_xmsyz(x,m,s,y)  = (inpNatZ_xmsyc(x,m,s,y)(LAST_COL+1,LAST_COL+(nZCs-1))).shift(1);
        }
    }
    aggregateRawNatZ();
    doTailCompression();
    normalize();
    if (debug) std::cout<<"end SizeFrequencyData::replaceSizeFrequencyData(...) for "<<name<<std::endl;
}

/**
 * Update catch-at-size data NatZ_xmsyz with new data for year y. 
 * Also modifies inpNatZ_xmsyc to reflect new data.
 * 
 * @param y - year to add
 * @param newNatZ_xmsz - d4_array of numbers-at-size by xmsz
 */
void SizeFrequencyData::addSizeFrequencyData(int y, d4_array& newNatZ_xmsz){
    if (debug) os<<"starting SizeFrequencyData::addSizeFrequencyData(...) for "<<name<<std::endl;
    
    //copy old values to temporary variables
    if (debug) os<<"Copying old values to temporary storage"<<endl;
    int oldNY = ny;
    ivector oldYrs(1,oldNY); oldYrs = yrs;
    d5_array oldInpNatZ_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,oldNY,1,LAST_COL+(nZCs-1));
    for (int x=1;x<=tcsam::ALL_SXs;x++){
        for (int m=1;m<=tcsam::ALL_MSs;m++){
            for (int s=1;s<=tcsam::ALL_SCs;s++) {
                oldInpNatZ_xmsyc(x,m,s) = inpNatZ_xmsyc(x,m,s);
            }
        }
    }
    if (debug) os<<"Copied old values to temporary storage"<<endl;
        
    //reallocate yrs for new number of years
    if (debug) os<<"Reallocating for new number of years"<<endl;
    ny = ny+1;
    yrs.deallocate(); yrs.allocate(1,ny); yrs.initialize();
    yrs(1,oldNY) = oldYrs; yrs(ny) = y;
    
    //reallocate inpNatZ_xmsyc for new number of years
    inpNatZ_xmsyc.deallocate(); 
    inpNatZ_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,LAST_COL+(nZCs-1));
    inpNatZ_xmsyc.initialize();
    if (debug) os<<"Done reallocating for new number of years"<<endl;

    //copy old data
    if (debug) os<<"Copying old data to new arrays"<<endl;
    for (int iy=1;iy<=oldNY;iy++){//year index for old data
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++) {
                    inpNatZ_xmsyc(x,m,s,iy) = oldInpNatZ_xmsyc(x,m,s,iy);//old use flag, DM, year, ss, size frequency
                }//-s
            }//-m
        }//-x
    }//-iy
    if (debug) os<<"Done copying old data to new arrays"<<endl;
    
    //copy new data
    if (debug) os<<"Copying new data to new arrays"<<endl;
    for (int x=1;x<=tcsam::ALL_SXs;x++){
        for (int m=1;m<=tcsam::ALL_MSs;m++){
            for (int s=1;s<=tcsam::ALL_SCs;s++) {
                if (debug) os<<"x,m,s = "<<x<<cc<<m<<cc<<s<<endl;
                dvector n_z = tcsam::extractFromXMSZ(x,m,s,newNatZ_xmsz);
                if (debug) os<<"n_z = "<<n_z<<endl;
                int k=0;
                inpNatZ_xmsyc(x,m,s,ny)(++k) = 1;                                 //set UF=1 to use new ZC
                inpNatZ_xmsyc(x,m,s,ny)(++k) = oldInpNatZ_xmsyc(x,m,s,oldNY)(k);  //new DM = last year's DM
                inpNatZ_xmsyc(x,m,s,ny)(++k) = y;                                 //new year
                inpNatZ_xmsyc(x,m,s,ny)(++k) = oldInpNatZ_xmsyc(x,m,s,oldNY)(k);  //new ss = last year's ss
                inpNatZ_xmsyc(x,m,s,ny)(LAST_COL+1,LAST_COL+nZCs-1).shift(1) = n_z;//new size frequency
            }//-s
        }//-m
    }//-x
    if (debug) os<<"Done copying new data to new arrays"<<endl;
    
    //re-constitute other arrays
    if (debug) os<<"Reconstituting other arrays"<<endl;
    inpUF_xmsy.deallocate();
    inpDM_xmsy.deallocate();
    inpSS_xmsy.deallocate();
    rawNatZ_xmsyz.deallocate();
    
    inpUF_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpDM_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpSS_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    rawNatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    
    inpUF_xmsy.initialize();
    inpDM_xmsy.initialize();
    inpSS_xmsy.initialize();
    rawNatZ_xmsyz.initialize();
    
    //reset re-weighting multipliers
    cumF_xms.deallocate();
    cumF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    cumF_xms.initialize();
    cumF_xms = 1.0;
    
    itrF_xms.allocate();
    itrF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    itrF_xms.initialize();
    
    int nc = factors.indexmax();
    for (int i=0;i<nc;i++){
        int x = tcsam::getSexType(factors(i+1,1));
        int m = tcsam::getMaturityType(factors(i+1,2));
        int s = tcsam::getShellType(factors(i+1,3));
            inpUF_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),1);
            inpDM_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),2);
            inpSS_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),4);
            for (int iy=1;iy<=ny;iy++){
                rawNatZ_xmsyz(x,m,s,iy)  = (inpNatZ_xmsyc(x,m,s,iy)(LAST_COL+1,LAST_COL+(nZCs-1))).shift(1);
            }
    }
    
    aggregateRawNatZ();
    doTailCompression();
    normalize();
    if (debug) os<<"Done reconstituting other arrays"<<endl;
    if (debug) std::cout<<"end SizeFrequencyData::addSizeFrequencyData(...) for "<<name<<std::endl;
}

/**
 * Set the maximum year in which to fit the data.
 * 
 * @param mxYr - the max year to include data
 */
void SizeFrequencyData::setMaxYear(int mxYr){
    debug=1;
    if (debug) {
        cout     <<"Starting SizeFrequencyData::setMaxYear("<<mxYr<<") for "<<name<<endl;
        rpt::echo<<"Starting SizeFrequencyData::setMaxYear("<<mxYr<<") for "<<name<<endl;
    }
    //determine number of years to keep
    int nyp = 0;
    for (int iy=1;iy<=ny;iy++) {if (yrs(iy)<=mxYr) nyp++;}
    if (debug) rpt::echo<<"ny, nyp = "<<ny<<tb<<nyp<<endl;
    
    //allocate temporary arrays
    ivector new_yrs(1,nyp);
    d4_array new_inpUF_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp);
    i4_array new_idxDM_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp);
    d4_array new_inpSS_xmsy(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp);
    d5_array new_NatZ_xmsyz(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp,1,nZCs-1);
    d5_array new_inpNatZ_xmsyc(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,nyp,1,LAST_COL+(nZCs-1));
    
    new_yrs.initialize();
    new_inpUF_xmsy.initialize();
    new_idxDM_xmsy.initialize();
    new_inpSS_xmsy.initialize();
    new_NatZ_xmsyz.initialize();
    new_inpNatZ_xmsyc.initialize();
    
    //set re-weighting multipliers
    cumF_xms.deallocate();
    cumF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    cumF_xms.initialize();
    cumF_xms = 1.0;
    
    itrF_xms.deallocate();
    itrF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    itrF_xms.initialize();
    
    int iyp = 0;
    for (int iy=1;iy<=ny;iy++) {
        if (yrs(iy)<=mxYr) {
            iyp++;
            new_yrs(iyp) = yrs(iy);
            for (int x=1;x<=tcsam::ALL_SXs;x++){
                for (int m=1;m<=tcsam::ALL_MSs;m++){
                    for (int s=1;s<=tcsam::ALL_SCs;s++){
                        new_inpUF_xmsy(x,m,s,iyp)    = inpUF_xmsy(x,m,s,iy);
                        new_idxDM_xmsy(x,m,s,iyp)    = inpDM_xmsy(x,m,s,iy);
                        new_inpSS_xmsy(x,m,s,iyp)    = inpSS_xmsy(x,m,s,iy);
                        new_NatZ_xmsyz(x,m,s,iyp)    = rawNatZ_xmsyz(x,m,s,iy);
                        new_inpNatZ_xmsyc(x,m,s,iyp) = inpNatZ_xmsyc(x,m,s,iy);
                    }//--s
                }//--m
            }//--x
         }//--if
    }//--iy
    if (iyp!=nyp){
        cout<<"Something wrong in AggregateCatchData::setMaxYear"<<endl;
        cout<<"iyp != nyp"<<endl;
        exit(0);
    }
    
    //copy back to class members
    yrs.deallocate();
    inpUF_xmsy.deallocate();
    inpDM_xmsy.deallocate();
    inpSS_xmsy.deallocate();
    rawNatZ_xmsyz.deallocate();
    inpNatZ_xmsyc.deallocate();
    
    ny = nyp;
    yrs.allocate(1,ny);
    inpUF_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpDM_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpSS_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    rawNatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    inpNatZ_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,LAST_COL+(nZCs-1));
    
    yrs.initialize();
    inpUF_xmsy.initialize();
    inpDM_xmsy.initialize();
    inpSS_xmsy.initialize();
    rawNatZ_xmsyz.initialize();
    inpNatZ_xmsyc.initialize();
    
    for (int iy=1;iy<=ny;iy++) {
        yrs(iy) = new_yrs(iy);
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++){
                    inpUF_xmsy(x,m,s,iy)    = new_inpUF_xmsy(x,m,s,iy);
                    inpDM_xmsy(x,m,s,iy)    = new_idxDM_xmsy(x,m,s,iy);
                    inpSS_xmsy(x,m,s,iy)    = new_inpSS_xmsy(x,m,s,iy);
                    rawNatZ_xmsyz(x,m,s,iy) = new_NatZ_xmsyz(x,m,s,iy);
                    inpNatZ_xmsyc(x,m,s,iy) = new_inpNatZ_xmsyc(x,m,s,iy);
                }//--s
            }//--m
        }//--x
    }//--iy
    
    if (debug){
        for (int x=1;x<=tcsam::ALL_SXs;x++){
            for (int m=1;m<=tcsam::ALL_MSs;m++){
                for (int s=1;s<=tcsam::ALL_SCs;s++){
                    if (sum(inpNatZ_xmsyc(x,m,s))>0){
                        rpt::echo<<tcsam::getSexType(x)<<tb<<tcsam::getMaturityType(m)<<tb<<tcsam::getShellType(s)<<endl;
                        rpt::echo<<inpNatZ_xmsyc(x,m,s)<<endl;
                        rpt::echo<<yrs<<endl;
                        for (int iy=1;iy<=ny;iy++){
                            rpt::echo<<inpUF_xmsy(x,m,s,iy)<<tb<<inpDM_xmsy(x,m,s,iy)<<tb<<yrs(iy)<<tb<<inpSS_xmsy(x,m,s,iy)<<tb<<rawNatZ_xmsyz(x,m,s,iy)<<endl;
                        }//--iy
                    }//--if
                }//--s
            }//--m
        }//--x
    }//--if(debug)
    
    aggregateRawNatZ();
    doTailCompression();
    normalize();
    
    if (debug) {
        cout     <<"Finished SizeFrequencyData::setMaxYear("<<mxYr<<") for "<<name<<endl;
        rpt::echo<<"Finished SizeFrequencyData::setMaxYear("<<mxYr<<") for "<<name<<endl;
    }
}

/*******************************************************\n
*   read from input stream.\n
******************************************************/
void SizeFrequencyData::read(cifstream & is){
  int olddebug=debug;
  debug=1;
    if (debug) {
        std::cout<<"start SizeFrequencyData::read(...) for "<<name<<std::endl;
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
        std::cout<<"Expected keyword '"<<KW_SIZEFREQUENCY_DATA<<"' but got '"<<str<<"'"<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    
    //NUMBERS-AT-SIZE 
    is>>str; optFit = tcsam::getFitType(str);
    rpt::echo<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    is>>str; llType = tcsam::getLikelihoodType(str);
    rpt::echo<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    is>>llWgt;
    rpt::echo<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    tc_limits.allocate(1,2);//tail compression limits (pmin, 1-pmax)
    is>>tc_limits;
    rpt::echo<<tc_limits<<tb<<"#tail compression limits (pmin, 1-pmax)"<<std::endl;
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
    
    qcsvZBs = wts::to_qcsv(zBs);//--size bin centers as comma-separated string of quoted values
    
    yrs.allocate(1,ny);
    inpUF_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpDM_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    inpSS_xmsy.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny);
    rawNatZ_xmsyz.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,nZCs-1);
    inpNatZ_xmsyc.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs,1,ny,1,LAST_COL+(nZCs-1));
    
    yrs.initialize();
    inpUF_xmsy.initialize();
    inpDM_xmsy.initialize();
    inpSS_xmsy.initialize();
    rawNatZ_xmsyz.initialize();
    inpNatZ_xmsyc.initialize();
    
    //set re-weighting multipliers
    cumF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    cumF_xms.initialize();
    cumF_xms = 1.0;
    itrF_xms.allocate(1,tcsam::ALL_SXs,1,tcsam::ALL_MSs,1,tcsam::ALL_SCs);
    itrF_xms.initialize();
    
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
            rpt::echo<<"#use? dm  year  ss     "<<zBs<<std::endl;
            rpt::echo<<inpNatZ_xmsyc(x,m,s)<<std::endl;
            int k=1;
            inpUF_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),k++);
            inpDM_xmsy(x,m,s) = (ivector) column(inpNatZ_xmsyc(x,m,s),k++);
            yrs               = (ivector) column(inpNatZ_xmsyc(x,m,s),k++);
            inpSS_xmsy(x,m,s) = column(inpNatZ_xmsyc(x,m,s),k++);
            for (int y=1;y<=ny;y++){
                rawNatZ_xmsyz(x,m,s,y)  = (inpNatZ_xmsyc(x,m,s,y)(LAST_COL+1,LAST_COL+(nZCs-1))).shift(1);
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

    if (debug) {std::cout<<"starting aggregateRawNatZ"<<endl; rpt::echo<<"starting aggregateRawNatZ"<<endl;}
    aggregateRawNatZ();
    if (debug) {std::cout<<"starting doTailCompression"<<endl; rpt::echo<<"starting doTailCompression"<<endl;}
    doTailCompression();
    if (debug) {std::cout<<"starting normalize"<<endl; rpt::echo<<"starting normalize"<<endl;}
    normalize();    
    
    if (debug) {
      std::cout<<"end SizeFrequencyData::read(...) for "<<name<<std::endl; 
      rpt::echo<<"end SizeFrequencyData::read(...) for "<<name<<std::endl;
    }
    debug=olddebug;
}
/***************************************************************
*   write.                                                     *
***************************************************************/
void SizeFrequencyData::write(ostream & os){
    if (debug) std::cout<<"start SizeFrequencyData::write(...) for "<<name<<std::endl;
    os<<KW_SIZEFREQUENCY_DATA<<tb<<"#required keyword"<<std::endl;
    os<<tcsam::getFitType(optFit)<<tb<<"#objective function fitting option"<<std::endl;
    os<<tcsam::getLikelihoodType(llType)<<tb<<"#likelihood function type"<<std::endl;
    os<<llWgt<<tb<<"#likelihood weight (multiplier)"<<std::endl;
    os<<tc_limits<<tb<<"#tail compression limits (pmin,1-pmax)"<<std::endl;    
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
        os<<"#use? dm  year  ss     "<<zBs<<std::endl;
        os<<inpNatZ_xmsyc(x,m,s)<<std::endl;
    }
    if (debug) std::cout<<"end SizeFrequencyData::write(...) for "<<name<<std::endl;
}
/***************************************************************
*   Function to write object to R list.                        *
***************************************************************/
void SizeFrequencyData::writeToR(ostream& os, std::string nm, int indent) {
    if (debug) tcsam::writeMessage(adstring("SizeFrequencyData::writing to R: ")+adstring(nm.c_str()));
    ivector bnds = wts::getBounds(rawNatZ_xmsyz);
    if (debug) tcsam::writeMessage("got here");
    adstring x = tcsamDims::getSXsForR(bnds(1),bnds(2));
    adstring m = tcsamDims::getMSsForR(bnds(3),bnds(4));
    adstring s = tcsamDims::getSCsForR(bnds(5),bnds(6));
    adstring y = "y=c("+wts::to_qcsv(yrs)+")";
    adstring z = "z=c("+wts::to_qcsv(zBs)+")";        
    for (int n=0;n<indent;n++) os<<tb;
        os<<"nAtZ=list(units="<<qt<<units<<qt<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"optFit="<<qt<<tcsam::getFitType(optFit)<<qt<<cc; 
        os<<"llType="<<qt<<tcsam::getLikelihoodType(llType)<<qt<<cc; 
        os<<"llWgt="<<llWgt<<cc<<std::endl; 
    for (int n=0;n<indent;n++) os<<tb;
        os<<"y="; wts::writeToR(os,yrs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"zc="; wts::writeToR(os,zCs); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"inpUF="; wts::writeToR(os,inpUF_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"inpDM="; wts::writeToR(os,inpDM_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"inpSS="; wts::writeToR(os,inpSS_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"aggUF="; wts::writeToR(os,uf_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"aggDM="; wts::writeToR(os,dm_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"aggSS="; wts::writeToR(os,ss_xmsy,x,m,s,y); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"cumFs="; wts::writeToR(os,cumF_xms,x,m,s); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"itrFs="; wts::writeToR(os,itrF_xms,x,m,s); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"tail.compression="; wts::writeToR(os,tc_xmsyc,x,m,s,y,"c=c(1,2)"); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"data="<<std::endl;
        wts::writeToR(os,rawNatZ_xmsyz,x,m,s,y,z); os<<cc<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
        os<<"aggdata="<<std::endl;
        wts::writeToR(os,aggNatZ_xmsyz,x,m,s,y,z); os<<std::endl;
    for (int n=0;n<indent;n++) os<<tb;
         os<<")";
    if (debug) tcsam::writeMessage("SizeFrequencyData::done writing to R");
}
/////////////////////////////////end SizeFrequencyData/////////////////////////
//----------------------------------------------------------------------
//          BioData
//----------------------------------------------------------------------
const adstring BioData::KW_BIO_DATA = "BIO_DATA";
const adstring BioData::VERSION     = "2023.03.07";
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
    
    is>>str;
    if (!(str==VERSION)){
        std::cout<<"#Error reading bio data from "<<is.get_file_name()<<std::endl;
        std::cout<<"Expected version '"<<VERSION<<"' but got '"<<str<<"'"<<std::endl;
        std::cout<<"Aborting..."<<std::endl;
        exit(-1);
    }
    rpt::echo<<str<<tb<<"#version"<<std::endl;
    
    //RECRUITMENT LAG
    is>>recLag;
    rpt::echo<<recLag<<tb<<"#recLag"<<std::endl;
    
    //WEIGHT-AT-SIZE
    {is>>unitsWatZ;
    rpt::echo<<unitsWatZ<<tb<<"#unitsWatZ"<<std::endl;
    double convToKG = tcsam::getConversionMultiplier(unitsWatZ,tcsam::UNITS_KG);
    rpt::echo<<"#using conversion factor from "<<unitsWatZ<<" to kg: "<<convToKG<<std::endl;
    alpha_xm.allocate(1,tcsam::nSXs,1,tcsam::nMSs);
    alpha_xm.initialize();
    beta_xm.allocate(1,tcsam::nSXs,1,tcsam::nMSs);
    beta_xm.initialize();
    wAtZ_xmz.allocate(1,tcsam::nSXs,1,tcsam::nMSs,1,ModelConfiguration::nZBs);
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
            is>>alpha_xm(sex,mat);
            is>>beta_xm(sex,mat);
            wAtZ_xmz(sex,mat) = alpha_xm(sex,mat)*pow(ModelConfiguration::zMidPts,beta_xm(sex,mat));
            rpt::echo<<wAtZ_xmz(sex,mat)<<std::endl;
            wAtZ_xmz(sex,mat) *= convToKG;//convert to kg
            if (debug) {
                std::cout<<"#"<<factors(1)<<cc<<factors(2)<<": alpha = "<<alpha_xm(sex,mat)<<" beta = "<<beta_xm(sex,mat)<<std::endl;
                std::cout<<"#wAtZ_xmz("<<factors(1)<<cc<<factors(2)<<") [kg] ="<<std::endl;
                std::cout<<"#"<<ModelConfiguration::zMidPts<<std::endl;
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
    
    {os<<"#-----------WEIGHT-AT-SIZE----------------------------#"<<std::endl;
    os<<unitsWatZ<<tb<<"#units for weight-at-size"<<std::endl;
    os<<tcsam::nSXs*tcsam::nMSs<<tb<<"#number of factor combinations (sex x maturity state)"<<std::endl;
    os<<"# sex   maturity   alpha   beta"<<endl;
    adstring_array factors(1,2);
    for (int sex=1;sex<=tcsam::nSXs;sex++){
        factors(1) = tcsam::getSexType(sex);
        for (int mat=1;mat<=tcsam::nMSs;mat++){
            factors(2) = tcsam::getMaturityType(mat);
            os<<factors(1)<<tb<<factors(2)<<tb<<alpha_xm(sex,mat)<<tb<<beta_xm(sex,mat)<<std::endl;
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
/**
 * Write to an output stream as an R list. 
 * 
 * @param os - the output stream
 * @param nm - the name for the list
 * @param indent - the number of tabs to indent, initially
 */
void BioData::writeToR(ostream& os, string nm, int indent) {
    dvector zBins = ModelConfiguration::zMidPts;
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
/**
 * Class constructor.
 * 
 * @param ptrMC - pointer to the global ModelConfiguration object
 */
ModelDatasets::ModelDatasets(ModelConfiguration* ptrMC){
    pMC=ptrMC;
    ptrBio=0;
    ppFsh=0;
    ppSrv=0;
    ppGrw=0;
    ppCHD=0;
    ppMOD=0;
}
/**
 * Class destructor. Cleans up various pointers.
 */
ModelDatasets::~ModelDatasets(){
    pMC=0;
    delete ptrBio;  ptrBio=0;
    if (ppFsh) {                                                        
        for (int f=0;f<nFsh;f++) delete ppFsh[f];
        delete ppFsh; ppFsh = 0;
    } 
    if (ppSrv) {
        for (int s=0;s<nSrv;s++) delete ppSrv[s];
        delete ppSrv; ppSrv = 0;
    } 
    if (ppGrw) {
        for (int s=0;s<nGrw;s++) delete ppGrw[s];
        delete ppGrw; ppGrw = 0;
    } 
    if (ppCHD) {
        for (int s=0;s<nCHD;s++) delete ppCHD[s];
        delete ppCHD; ppCHD = 0;
    } 
    if (ppMOD) {
        for (int s=0;s<nMOD;s++) delete ppMOD[s];
        delete ppMOD; ppMOD = 0;
    } 
}
/**
 * Read ModelDatasets file from input filestream.
 * 
 * @param is - the input filestream
 */
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
    is>>nGrw;
    rpt::echo<<nGrw<<tb<<"#number of growth datasets to read in"<<std::endl;
    if (nGrw){
        fnsGrowthData.allocate(1,nGrw);
        is>>fnsGrowthData; 
        for (int i=1;i<=nGrw;i++) {
            fnsGrowthData[i] = wts::concatenateFilePaths(parent,fnsGrowthData[i]);
            rpt::echo<<fnsGrowthData[i]<<tb<<"#growth dataset "<<i<<std::endl;
        }
    }
    is>>nCHD;
    rpt::echo<<nCHD<<tb<<"#number of chela height datasets to read in"<<std::endl;
    if (nCHD){
        fnsChelaHeightData.allocate(1,nCHD);
        is>>fnsChelaHeightData; 
        for (int i=1;i<=nCHD;i++) {
            fnsChelaHeightData[i] = wts::concatenateFilePaths(parent,fnsChelaHeightData[i]);
            rpt::echo<<fnsChelaHeightData[i]<<tb<<"#chela height dataset "<<i<<std::endl;
        }
    }
    is>>nMOD;
    rpt::echo<<nMOD<<tb<<"#number of maturity ogive datasets to read in"<<std::endl;
    if (nMOD){
        fnsMaturityOgiveData.allocate(1,nMOD);
        is>>fnsMaturityOgiveData; 
        for (int i=1;i<=nMOD;i++) {
            fnsMaturityOgiveData[i] = wts::concatenateFilePaths(parent,fnsMaturityOgiveData[i]);
            rpt::echo<<fnsMaturityOgiveData[i]<<tb<<"#maturity ogive dataset "<<i<<std::endl;
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
    //          Growth data
    if (nGrw) {
        rpt::echo<<"#-------Growth Datasets---------"<<std::endl;
        ppGrw = new GrowthData*[nGrw];
        for (int i=0;i<nGrw;i++) {
            ppGrw[i] = new GrowthData();
            cifstream strm(fnsGrowthData(i+1),ios::in);
            rpt::echo<<std::endl<<"#----------Growth Data "<<i+1<<"-----"<<std::endl;
            strm>>(*ppGrw[i]);
        }
    }
    //          Chela Height Data
    if (nCHD) {
        rpt::echo<<"#-------Chela Height Datasets---------"<<std::endl;
        ppCHD = new ChelaHeightData*[nCHD];
        for (int i=0;i<nCHD;i++) {
            ppCHD[i] = new ChelaHeightData();
            cifstream strm(fnsChelaHeightData(i+1),ios::in);
            rpt::echo<<std::endl<<"#----------Chela Height Data "<<i+1<<"-----"<<std::endl;
            strm>>(*ppCHD[i]);
        }
    }
    //          Maturity Ogive Data
    if (nMOD) {
        rpt::echo<<"#-------Maturity Ogive Datasets---------"<<std::endl;
        ppMOD = new MaturityOgiveData*[nMOD];
        for (int i=0;i<nMOD;i++) {
            ppMOD[i] = new MaturityOgiveData();
            cifstream strm(fnsMaturityOgiveData(i+1),ios::in);
            rpt::echo<<std::endl<<"#----------Maturity Ogive Data "<<i+1<<"-----"<<std::endl;
            strm>>(*ppMOD[i]);
        }
    }
}
/**
 * Write ModelDatasets info to output text file in ADMB format.
 * 
 * @param ofn - the output file name
 */
void ModelDatasets::write(adstring fn){
    ofstream os; os.open(fn, ios::trunc);
    write(os);
    os.close();
}
/**
 * Write ModelDatasets info to output stream in ADMB format.
 * 
 * @param os - the output stream
 */
void ModelDatasets::write(ostream & os){    
    if (debug) std::cout<<"start ModelDatasets::write(...) "<<this<<std::endl;
    os<<fnBioData<<tb<<"#tanner crab biological data file"<<std::endl;
    os<<"#-------fishery data files---------"<<std::endl;
    os<<nFsh<<tb<<"#number of fishery data files"<<std::endl;
    for (int i=1;i<=nFsh;i++) os<<fnsFisheryData[i]<<tb<<"#fishery dataset "<<i<<std::endl;
    os<<"#-------survey data files---------"<<std::endl;
    os<<nSrv<<tb<<"#number of survey data files"<<std::endl;
    for (int i=1;i<=nSrv;i++) os<<fnsSurveyData[i]<<tb<<"#survey dataset "<<i<<std::endl;
    os<<"#-------growth data files---------"<<std::endl;
    os<<nGrw<<tb<<"#number of growth data files"<<std::endl;
    for (int i=1;i<=nGrw;i++) os<<fnsGrowthData[i]<<tb<<"#growth dataset "<<i<<std::endl;
    os<<"#-------chela height data files---------"<<std::endl;
    os<<nCHD<<tb<<"#number of chela height data files"<<std::endl;
    for (int i=1;i<=nCHD;i++) os<<fnsChelaHeightData[i]<<tb<<"#chela height dataset "<<i<<std::endl;
    os<<"#-------maturity ogive data files---------"<<std::endl;
    os<<nMOD<<tb<<"#number of maturity ogive data files"<<std::endl;
    for (int i=1;i<=nMOD;i++) os<<fnsMaturityOgiveData[i]<<tb<<"#maturity ogive dataset "<<i<<std::endl;
   if (debug) std::cout<<"end ModelDatasets::write(...) "<<this<<std::endl;
}
/**
 * Write to an output stream as an R list. 
 * 
 * @param os - the output stream
 * @param nm - the name for the list
 * @param indent - the number of tabs to indent, initially
 */
void ModelDatasets::writeToR(ostream& os, std::string nm, int indent) {
    for (int n=0;n<indent;n++) os<<tb;
    os<<nm<<"=list("<<std::endl;
    indent++;
        //bio data as list
            ptrBio->writeToR(os,"bio",indent); os<<cc<<std::endl;
        
        //growth data
        adstring grwNm;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"growth=list("<<std::endl;
        indent++;
            if (ppGrw) {
                for (int i=0;i<(nGrw-1);i++) {
                    grwNm = "growth_"+str(i+1);
                    ppGrw[i]->writeToR(os,(char*)grwNm,indent); os<<cc<<std::endl;
                }
                grwNm = "growth_"+str(nGrw);
                ppGrw[nGrw-1]->writeToR(os,(char*)grwNm,indent); os<<std::endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<std::endl;            
        
        //chela height data
        adstring chdNm;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"chelaheight=list("<<std::endl;
        indent++;
            if (ppCHD) {
                for (int i=0;i<(nCHD-1);i++) {
                    chdNm = "chelaheight_"+str(i+1);
                    ppCHD[i]->writeToR(os,(char*)chdNm,indent); os<<cc<<std::endl;
                }
                chdNm = "chelaheight_"+str(nCHD);
                ppCHD[nCHD-1]->writeToR(os,(char*)chdNm,indent); os<<std::endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<std::endl;            
        
        //maturity ogive data
        adstring modNm;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"maturityogives=list("<<std::endl;
        indent++;
            if (ppMOD) {
                for (int i=0;i<(nCHD-1);i++) {
                    chdNm = "maturityogives_"+str(i+1);
                    ppMOD[i]->writeToR(os,(char*)chdNm,indent); os<<cc<<std::endl;
                }
                chdNm = "chelaheight_"+str(nCHD);
                ppMOD[nMOD-1]->writeToR(os,(char*)modNm,indent); os<<std::endl;
            }
        indent--;
        for (int n=0;n<indent;n++) os<<tb;
        os<<"),"<<std::endl;            
        
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
