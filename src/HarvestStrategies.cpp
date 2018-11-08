/* 
 * File:   HarvestStrategies.cpp
 * Author: william.stockhausen
 *
 * Created on September 26, 2018, 10:35 AM
 * 2018-9-28 Updated Directories in NetBeans Project
 * 2018-10-1 Meeting with Andre 
 * 2018-11-07 Incorporated into TCSAM02.
 */

//#include <cstdlib>
#include <admodel.h>
#include <cstdlib>
#include "HarvestStrategies.hpp"

using namespace std;

int HarvestStrategies::debug = 1;
const double HarvestStrategies::MMB_AveMilLb =  24; //dummy values for testing
const double HarvestStrategies::MFB_AveMilLb = 12; // dummy values for testing 
const double HarvestStrategies::z = 4; // dummy values for testing, number of ELM size classes 

/**
 * HCR 1: female ramp.
 * 
 * @param MFB - Mature female biomass for the year 
 * @param aveMFB - averaged MFB from 1982-2016
 * @param MMB - Mature Male Biomass for the Year 
 * 
 * @return TAC = exploitation rate time Mature Male Biomass 
 */
double HarvestStrategies::HCR1_FemaleRamp(double MFB, double aveMFB, double MMB){
    double femRatio = MFB/aveMFB;
    double xp;
    double TAC;
    if (femRatio < 0.25){
        xp = 0;
    } else if (femRatio >= 1.0){
        xp = 0.2;
    } else{
        xp = femRatio*0.2;
    }
    TAC = MMB*xp;
    return TAC;
}

/**
 * HCR 2. 3 ramps.
 * 
 * @param MMB Mature Male Biomass in millions of lbs for the year 
 * @param aveMMB constant, MMB average from 1982-2016
 * @param rampID identifies the level of exploitation rate, low, moderate, or aggressive 
 * 
 * @return TAC
 */
double HarvestStrategies::HCR2_MaleRamp(double MMB, double aveMMB, int rampID){
    //Define look up table for ramps
    double slopes[] = {5.0/75.0, 10.0/75.0, 15.0/75.0};
    double intercepts[] = {0.0333, 0.01675, 0.0};
    
    double maleRatio = MMB/aveMMB;
    double xp = 0;
    if (maleRatio < 0.25){
        return 0;
    } else if (maleRatio >= 1.0){
        maleRatio = 1;
    }
    xp = maleRatio*slopes[rampID -1]+intercepts[rampID -1];
    
    double TAC = MMB*xp;
    
    return TAC;
}

/**
 * HCR 3. Takes OFL from operating model, uses pre-determined buffer.
 * 
 * @param OFL from the model
 * @param buffer <- current buffer on OFL (20% currently?)
 * 
 * @return TAC
 */
double HarvestStrategies::HCR3_ABC(double OFL, double buffer){
    //Define look up table for ramps
    double TAC; 
    double ABC = OFL*(1-buffer);
    
    TAC = ABC;
    
    return TAC;
}

/**
 * HCR 4. Scaled with females--> moving.
 * 
 * @param MFB
 * @param aveMFB
 * @param MMB
 * @param aveMMB
 * 
 * @return TAC
 */
double HarvestStrategies::HCR4_FemaleDimmer(double MFB, double aveMFB, double MMB, double aveMMB){
    double femRatio = MFB/aveMFB;
    double maleRatio = MMB/aveMMB;
    double slope;
    double intercept;
    
    double xp;
    double TAC;
    
    if (maleRatio < 0.25){
        return 0;
    } 
    
    if (maleRatio >= 0.25 && femRatio < 0.25){
            xp = 0.05;
    }
    
    else if (maleRatio >= 1.0){
        maleRatio = 1.0;
         if (femRatio >=1.0){
            femRatio = 1.0;
        }
        xp = femRatio/5;
    }
    
    if (femRatio >=1.0){
            femRatio = 1.0;
        }
       
    slope = ((femRatio-0.25)/5)/(maleRatio-0.25);
    intercept = -1*((slope*0.25)-0.05); // CHECK THIS WITH ANDRE y1???
        
    xp = maleRatio*slope + intercept; //I don't think i need to apply an intercept because of how the slope is taken, please confirm
    
    
    TAC = MMB*xp;
    
    return TAC;  
}

/**
 * HCR 5. Scaled with females --> blocks.
 * 
 * @param MFB
 * @param aveMFB
 * @param MMB
 * @param aveMMB
 * @param rampID
 * 
 * @return 
 */

double HarvestStrategies::HCR5_FemaleBlocks(double MFB, double aveMFB, double MMB, double aveMMB){
    double femRatio = MFB/aveMFB;
    double maleRatio = MMB/aveMMB;
    double slopes[] = {0, 5.0/75.0, 10.0/75.0, 15.0/75.0};
    double intercepts[4];
        for( int i = 0; i <=4; ++i){
            intercepts[i]= -1*((slopes[i]*0.25)-0.05);
        }
    // CHECK THIS LINE 
    double xp;
    double TAC;
   
    if (maleRatio< 0.25){
        return 0;
    }
    if ( femRatio < 0.3 && maleRatio >=0.25){
        xp = 0.05;
    }
   if ( femRatio >= 0.3 && femRatio <= 0.5 && maleRatio >=0.25){
       if(maleRatio >= 1.0){
           maleRatio =1.0;
       }
       xp =maleRatio*slopes[1]+intercepts[1];
       
    }
    if ( femRatio >= 0.5 && femRatio <= 0.7 && maleRatio >=0.25){
       if(maleRatio >= 1.0){
           maleRatio =1.0;
       }
       xp =maleRatio*slopes[2]+intercepts[2];  
       
    }
    if ( femRatio >= 0.7 && femRatio <= 1.0 && maleRatio >=0.25){
       if(maleRatio >= 1.0){
           maleRatio =1.0;
       }
       xp =maleRatio*slopes[3]+intercepts[3];   
       
    }
    TAC = MMB*xp;
    
    return TAC;   
}

/**
 * HCR 6. ELM. 
 * 
 * @param propNS proportion new shell
 * @param abunELM
 * @param weights
 * @param sOS
 * @param xpRate 30%, 40%, or 50%
 * 
 * @return 
 */
double HarvestStrategies::HCR6_ELM(double propNS, dvector abunELM, dvector weights, double sOS, int xpRate){
    double xpLM = 0; // exploitable legal males 
    double TAC;
    
    for(int i = abunELM.indexmin(); i<=abunELM.indexmax(); ++i){ // check this indexing format     
        xpLM += propNS*abunELM[i]*weights[i]+sOS*(1-propNS)*abunELM[i]*weights[i];
    }
    TAC=xpLM*xpRate;
    return TAC;
}
