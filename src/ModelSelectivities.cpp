/* 
 * File:   ModelSelectivites.cpp
 * Author: William.Stockhausen
 *
 * Created on August 13, 2012, 5:32 AM
 */

#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelSelectivities.hpp"

using namespace std;

//--------------------------------------------------------------------------------
//  Includes:
//      SelFcns
//--------------------------------------------------------------------------------

int SelFcns::debug = 0;//debug flag

//the following values should all be lower case
const adstring SelFcns::STR_ASCLOGISTIC        ="asclogistic";
const adstring SelFcns::STR_ASCLOGISTICLN50    ="asclogisticln50";
const adstring SelFcns::STR_ASCLOGISTIC5095    ="asclogistic5095";
const adstring SelFcns::STR_ASCLOGISTIC50LN95  ="asclogistic50ln95";
const adstring SelFcns::STR_ASCLOGISTICLN50LN95="asclogisticln50ln95";
const adstring SelFcns::STR_DBLLOGISTIC        ="dbllogistic";
const adstring SelFcns::STR_DBLLOGISTICLND50   ="dbllogisticlnd50";
const adstring SelFcns::STR_DBLLOGISTIC5095    ="dbllogistic5095";
const adstring SelFcns::STR_DBLLOGISTIC50LN95  ="dbllogistic50ln95";
const adstring SelFcns::STR_DBLLOGISTICLN50LN95="dbllogisticln50ln95";
const adstring SelFcns::STR_DBLNORMAL          ="dblnormal";
const adstring SelFcns::STR_ASCLOGISTIC5099    ="asclogistic5099";
const adstring SelFcns::STR_ASCLOGISTIC95LN50  ="asclogistic95ln50";

//--------------------------------------------------------------------------------
//          SelFcns
//  Encapsulates selectivity functions
//--------------------------------------------------------------------------------
/**
 * Class constructor.
 */
SelFcns::SelFcns(){}

/**
 * Class destructor.
 */
SelFcns::~SelFcns(){}

/**
 * Returns the integer id of the requested selectivity function.
 * 
 * @param str - selectivity function name
 * @return - integer id for selectivity function (or 0)
 */
int SelFcns::getSelFcnID(adstring str){
    int id = 0;
    str.to_lower();
    if (str==STR_ASCLOGISTIC)         return ID_ASCLOGISTIC;
    if (str==STR_ASCLOGISTICLN50)     return ID_ASCLOGISTICLN50;
    if (str==STR_ASCLOGISTIC5095)     return ID_ASCLOGISTIC5095;
    if (str==STR_ASCLOGISTIC5099)     return ID_ASCLOGISTIC5099;
    if (str==STR_ASCLOGISTIC50LN95)   return ID_ASCLOGISTIC50LN95;
    if (str==STR_ASCLOGISTICLN50LN95) return ID_ASCLOGISTICLN50LN95;
    if (str==STR_DBLLOGISTIC)         return ID_DBLLOGISTIC;
    if (str==STR_DBLLOGISTICLND50)    return ID_DBLLOGISTICLND50;
    if (str==STR_DBLLOGISTIC5095)     return ID_DBLLOGISTIC5095;
    if (str==STR_DBLLOGISTIC50LN95)   return ID_DBLLOGISTIC50LN95;
    if (str==STR_DBLLOGISTICLN50LN95) return ID_DBLLOGISTICLN50LN95;
    if (str==STR_DBLNORMAL)           return ID_DBLNORMAL;
    if (str==STR_ASCLOGISTIC5099)     return ID_ASCLOGISTIC5099;
    if (str==STR_ASCLOGISTIC95LN50)   return ID_ASCLOGISTIC95LN50;
    cout<<"Error in SelFcns::getSelFcnID(adstring str)"<<endl;
    cout<<"Function name '"<<str<<"' not a valid selectivity function name."<<endl;
    cout<<"Aborting..."<<endl;
    exit(-1);
    return id;
}

/**
 * Returns the selectivity function name based on its integer id.
 * 
 * @param id - integer id for selectivity function
 * @return   - selectivity function name (or "")
 */
adstring SelFcns::getSelFcnID(int id){
    if (debug) cout<<"Starting SelFcns::getSelFcnID("<<id<<")"<<endl;
    adstring str = "";
    switch(id){
        case ID_ASCLOGISTIC: 
            if (debug) cout<<"SelFcn = "<<STR_ASCLOGISTIC<<endl;
            return STR_ASCLOGISTIC;
        case ID_ASCLOGISTICLN50: 
            if (debug) cout<<"SelFcn = "<<STR_ASCLOGISTICLN50<<endl;
            return STR_ASCLOGISTICLN50;
        case ID_ASCLOGISTIC5095: 
            if (debug) cout<<"SelFcn = "<<STR_ASCLOGISTIC5095<<endl;
            return STR_ASCLOGISTIC5095;
        case ID_ASCLOGISTIC50LN95: 
            if (debug) cout<<"SelFcn = "<<STR_ASCLOGISTIC50LN95<<endl;
            return STR_ASCLOGISTIC50LN95;
        case ID_ASCLOGISTICLN50LN95: 
            if (debug) cout<<"SelFcn = "<<STR_ASCLOGISTICLN50LN95<<endl;
            return STR_ASCLOGISTICLN50LN95;
        case ID_DBLLOGISTIC: 
            if (debug) cout<<"SelFcn = "<<STR_DBLLOGISTIC<<endl;
            return STR_DBLLOGISTIC;
        case ID_DBLLOGISTICLND50: 
            if (debug) cout<<"SelFcn = "<<STR_DBLLOGISTICLND50<<endl;
            return STR_DBLLOGISTICLND50;
        case ID_DBLLOGISTIC5095: 
            if (debug) cout<<"SelFcn = "<<STR_DBLLOGISTIC5095<<endl;
            return STR_DBLLOGISTIC5095;
        case ID_DBLLOGISTIC50LN95: 
            if (debug) cout<<"SelFcn = "<<STR_DBLLOGISTIC50LN95<<endl;
            return STR_DBLLOGISTIC50LN95;
        case ID_DBLLOGISTICLN50LN95: 
            if (debug) cout<<"SelFcn = "<<STR_DBLLOGISTICLN50LN95<<endl;
            return STR_DBLLOGISTICLN50LN95;
        case ID_DBLNORMAL: 
            if (debug) cout<<"SelFcn = "<<STR_DBLNORMAL<<endl;
            return STR_DBLNORMAL;
        case ID_ASCLOGISTIC5099: 
            if (debug) cout<<"SelFcn = "<<STR_ASCLOGISTIC5099<<endl;
            return STR_ASCLOGISTIC5099;
        case ID_ASCLOGISTIC95LN50: 
            if (debug) cout<<"SelFcn = "<<STR_ASCLOGISTIC95LN50<<endl;
            return STR_ASCLOGISTIC95LN50;
        default:
        {
            cout<<endl;
            cout<<"Invalid id for SelFcns.getSelFcnID(id)"<<endl;
            cout<<"id was "<<id<<"."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }
    return str;
}

/**
 * Calculates the selectivity function identified by "id".
 * 
 * Inputs:
 * @param id - integer id of selectivity function
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameter values
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::calcSelFcn(int id,dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::calcSelFcns(...) id: "<<id<<endl;
    dvar_vector s(z.indexmin(),z.indexmax());
    s.initialize();
    switch (id){
        case ID_ASCLOGISTIC:         {s=asclogistic(z,params,fsZ);         break;}
        case ID_ASCLOGISTICLN50:     {s=asclogisticLn50(z,params,fsZ);     break;}
        case ID_ASCLOGISTIC5095:     {s=asclogistic5095(z,params,fsZ);     break;}
        case ID_ASCLOGISTIC50LN95:   {s=asclogistic50Ln95(z,params,fsZ);   break;}
        case ID_ASCLOGISTICLN50LN95: {s=asclogisticLn50Ln95(z,params,fsZ); break;}
        case ID_DBLLOGISTIC:         {s=dbllogistic(z,params,fsZ);         break;}
        case ID_DBLLOGISTICLND50:    {s=dbllogisticLnD50(z,params,fsZ);    break;}
        case ID_DBLLOGISTIC5095:     {s=dbllogistic5095(z,params,fsZ);     break;}
        case ID_DBLLOGISTIC50LN95:   {s=dbllogistic50Ln95(z,params,fsZ);   break;}
        case ID_DBLLOGISTICLN50LN95: {s=dbllogisticLn50Ln95(z,params,fsZ); break;}
        case ID_DBLNORMAL:           {s=dblnormal(z,params,fsZ);           break;}
        case ID_ASCLOGISTIC5099:     {s=asclogistic5099(z,params,fsZ);     break;}
        case ID_ASCLOGISTIC95LN50:   {s=asclogistic95Ln50(z,params,fsZ);   break;}
        default:
        {
            cout<<"Invalid id for SelFcns.calcSelFcn(id,...)";
            cout<<"id was "<<id<<"."<<endl;
            cout<<"Aborting..."<<endl;
            exit(-1);
        }
    }
    if (debug) cout<<"Finished SelFcns::calcSelFcns(...) id: "<<id<<endl;
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates ascending logistic function parameterized by 
 *      params[1]: size at 50% selected (z50)
 *      params[2]: slope
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return -
 */
dvar_vector SelFcns::asclogistic(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"Starting SelFcns::asclogistic(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    s = 1.0/(1.0+exp(-params(2)*(z-params(1))));
    if (fsZ>0){
        n = 1.0+exp(-params(2)*(fsZ-params(1)));//normalize so (s(fsZ) = 1
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::asclogistic(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates ascending logistic function parameterized by 
 *      params[1]: ln(size at 50% selected) (z50)
 *      params[2]: slope
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::asclogisticLn50(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) rpt::echo<<"Starting SelFcns::asclogisticLn50(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    s = 1.0/(1.0+exp(-params(2)*(z-exp(params(1)))));
    if (fsZ>0){
        n = 1.0+exp(-params(2)*(fsZ-exp(params(1))));//normalize so (s(fsZ) = 1
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::asclogisticLn50(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}      

/**
 * Calculates ascending logistic function parameterized by 
 *      params[1]: size at 50% selected (z50)
 *      params[2]: increment from z50 to z95
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::asclogistic5095(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::asclogistic5095(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    s = 1.0/(1.0+exp(-log(19.0)*(z-params(1))/params(2)));
    if (fsZ>0){
        n = 1.0+exp(-log(19.0)*(fsZ-params(1))/params(2));//normalization constant
        s *= n;
     } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::asclogistic5095(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates ascending logistic function parameterized by 
 *      params[1]: size at 50% selected (z50)
 *      params[2]: size at 99% selected (z99)
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::asclogistic5099(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::asclogistic5099(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    s = 1.0/(1.0+exp(-log(99.0)*(z-params(1))/(params(2)-params(1))));
    if (fsZ>0){
        n = 1.0+exp(-log(99.0)*(fsZ-params(1))/params(2));//normalization constant
        s *= n;
     } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::asclogistic5099(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates ascending logistic function parameterized by 
 *      params[1]: size at 50% selected (z50)
 *      params[2]: ln-scale increment from z50 to size at 95% selected
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return -
 */
dvar_vector SelFcns::asclogistic50Ln95(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::asclogistic50Ln95(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    dvariable z50    = params(1);
    dvariable dz5095 = exp(params(2));
    s = 1.0/(1.0+exp(-log(19.0)*(z-z50)/dz5095));
    if (fsZ>0){
        n = 1.0+exp(-log(19.0)*(fsZ-z50)/dz5095);//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::asclogistic50Ln95(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates ascending logistic function parameterized by 
 *      params[1]: size at 95% selected (z95)
 *      params[2]: ln-scale increment from z50 to size at 95% selected
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return -
 */
dvar_vector SelFcns::asclogistic95Ln50(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::asclogistic95Ln50(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    dvariable dz5095 = exp(params(2));
    dvariable z50    = params(1) - dz5095;
    s = 1.0/(1.0+exp(-log(19.0)*(z-z50)/dz5095));
    if (fsZ>0){
        n = 1.0+exp(-log(19.0)*(fsZ-z50)/dz5095);//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::asclogistic95Ln50(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates ascending logistic function parameterized by 
 *      params[1]: ln-scale size at 50% selected (ln(z50))
 *      params[2]: ln-scale increment to size at 95% selected
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return -
 */
dvar_vector SelFcns::asclogisticLn50Ln95(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::asclogisticLn50Ln95(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    dvariable z50    = exp(params(1));
    dvariable dz5095 = exp(params(2));
    s = 1.0/(1.0+exp(-log(19.0)*(z-z50)/dz5095));
    if (fsZ>0){
        n = 1.0+exp(-log(19.0)*(fsZ-z50)/dz5095);//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::asclogisticLn50Ln95(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}
            
/**
 * Calculates double logistic function parameterized by 
 *      params[1]: size where ascending limb = 0.5  (z50)
 *      params[2]: ascending limb rate parameter    (slope)
 *      params[3]: size where descending limb = 0.5 (z50)
 *      params[4]: descending limb rate parameter   (slope)
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::dbllogistic(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::dbllogistic(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    s = elem_prod(1.0/(1.0+exp(-params(2)*(z-params(1)))),1.0/(1.0+exp(-params(4)*(params(3)-z))));
    if (fsZ>0){
        n = (1.0+exp(-params(2)*(fsZ-params(1))))*(1.0+exp(-params(4)*(params(3)-fsZ)));//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<params(3)<<tb<<params(4)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::dbllogistic(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates double logistic function parameterized by 
 *      params[1]: size where ascending limb = 0.5  (z50)
 *      params[2]: ascending limb rate parameter    (slope)
 *      params[3]: ln-scale increment from z50 to where descending limb = 0.5 (z50)
 *      params[4]: descending limb rate parameter   (slope)
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::dbllogisticLnD50(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::dbllogisticLnD50(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    s = elem_prod(1.0/(1.0+exp(-params(2)*(z-params(1)))),1.0/(1.0+exp(-params(4)*(params(1)+exp(params(3))-z))));
    if (fsZ>0){
        n = (1.0+exp(-params(2)*(fsZ-params(1))))*(1.0+exp(-params(4)*(params(1)+exp(params(3))-fsZ)));//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<params(3)<<tb<<params(4)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::dbllogisticLnD50(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates double logistic function parameterized by 
 *      params[1]: size where ascending limb = 0.5 (z50) 
 *      params[2]: increment from z50 to z95 on ascending limb
 *      params[3]: size where descending limb = 0.5 (z50) 
 *      params[4]: increment from z95 to z50 on descending limb
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::dbllogistic5095(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::dbllogistic5095(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    s = elem_prod(1.0/(1.0+exp(-log(19.0)*(z-params(1))/params(2))),1.0/(1.0+exp(-log(19.0)*(params(3)-z)/params(4))));
    if (fsZ>0){
        n = (1.0+exp(-log(19.0)*(fsZ-params(1))/params(2)))*(1.0+exp(-log(19.0)*(params(3)-fsZ)/params(4)));//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<params(3)<<tb<<params(4)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::dbllogistic5095(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates double logistic function parameterized by 
 *      params[1]: size where ascending limb = 0.5 (z50) 
 *      params[2]: log-scale increment from z50 to z95 on ascending limb
 *      params[3]: log-scale increment z95 on ascending limb to z95 on descending limb
 *      params[4]: log-scale increment from z95 to z50 on descending limb
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::dbllogistic50Ln95(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::dbllogistic50Ln95(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    dvariable z50a = params(1);
    dvariable dza  = exp(params(2));//increment from z50a to z95a
    dvariable dzz  = exp(params(3));//increment from z95a to z95d
    dvariable dzd  = exp(params(4));//increment from z95d to z50d
    dvariable z50d = z50a+dza+dzz+dzd; 
    s = elem_prod(1.0/(1.0+exp(-log(19.0)*(z-z50a)/dza)),1.0/(1.0+exp(-log(19.0)*(z50d-z)/dzd)));
    if (fsZ>0){
        n = (1.0+exp(-log(19.0)*(fsZ-z50a)/dza))*(1.0+exp(-log(19.0)*(z50d-fsZ)/dzd));//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<params(3)<<tb<<params(4)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::dbllogistic50Ln95(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

/**
 * Calculates double logistic function parameterized by 
 *      params[1]: log-scale size where ascending limb = 0.5 (z50) 
 *      params[2]: log-scale increment from z50 to z95 on ascending limb
 *      params[3]: log-scale increment z95 on ascending limb to z95 on descending limb
 *      params[4]: log-scale increment from z95 to z50 on descending limb
 * Inputs:
 * @param z      - dvector of sizes at which to compute function values
 * @param params - dvar_vector of function parameters
 * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
 * 
 * @return - selectivity function values as dvar_vector
 */
dvar_vector SelFcns::dbllogisticLn50Ln95(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::dbllogisticLn50Ln95(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    dvariable z50a = exp(params(1));
    dvariable dza  = exp(params(2));//increment from z50a to z95a
    dvariable dzz  = exp(params(3));//increment from z95a to z95d
    dvariable dzd  = exp(params(4));//increment from z95d to z50d
    dvariable z50d = z50a+dza+dzz+dzd; 
    s = elem_prod(1.0/(1.0+exp(-log(19.0)*(z-z50a)/dza)),1.0/(1.0+exp(-log(19.0)*(z50d-z)/dzd)));
    if (fsZ>0){
        n = (1.0+exp(-log(19.0)*(fsZ-z50a)/dza))*(1.0+exp(-log(19.0)*(z50d-fsZ)/dzd));//normalization constant
        s *= n;
    } else if (fsZ<0) {
        n = 1.0/max(s);
        s *= n; //normalize by max
    } //otherwise don't normalize it
    if (debug) {
        rpt::echo<<"params, fsZ = "<<params(1)<<tb<<params(2)<<tb<<params(3)<<tb<<params(4)<<tb<<fsZ<<endl;
        rpt::echo<<"n = "<<n<<endl;
        rpt::echo<<"z = "<<z<<endl;
        rpt::echo<<"s = "<<s<<endl;
        cout<<"Finished SelFcns::dbllogisticLn50Ln95(...)"<<endl;
    }
    RETURN_ARRAYS_DECREMENT();
    return s;
}

//TODO: implement this!
dvar_vector SelFcns::dblnormal(dvector& z, dvar_vector& params, double fsZ){
    RETURN_ARRAYS_INCREMENT();
    if (debug) cout<<"Starting SelFcns::dblnormal(...)"<<endl;
    dvariable n; n.initialize();
    dvar_vector s(z.indexmin(),z.indexmax()); s.initialize();
    if (debug) cout<<"Finished SelFcns::dblnormal(...)"<<endl;
    RETURN_ARRAYS_DECREMENT();
    return s;
}


