/* 
 * File:   ModelSelectivities.hpp
 * Author: William.Stockhausen
 *
 * Created on August 13, 2012, 9:22 AM
 * 
 * 2014-12-05: 1. Added asclogisticLn50Ln95 and dbllogisticLn50Ln95 functions
 * 2015-04-17: 1. Added dbllogistic50Ln95 function
 * 2016-11-15: 1. Added asclogisticLn50 and dbllogisticLn50 functions to match TCSAM2013
 * 2017-08-03: 1. Added asclogistic5099 function
 * 2017-08-23: 1. Added asclogistic95Ln50 function
 */

#ifndef MODELSELECTIVITIES_HPP
    #define	MODELSELECTIVITIES_HPP

    #include <admodel.h>

//--------------------------------------------------------------------------------
//          SelFcns
//  Encapsulates selectivity functions
//--------------------------------------------------------------------------------
    class SelFcns {
        public:
            static int debug;
        protected:
            adstring_array stdNames;//names of standard inputs
            adstring_array parNames;//parameter names
            adstring_array xtrNames;//names of extra inputs

        public:
            int nSel;        //number of selectivity functions defined
            adstring_array selNames;
            adstring_array typNames;
            imatrix paramIndexInfo;
            dmatrix xtrInfo;
        public:
            SelFcns();
            ~SelFcns();
            
            int static getSelFcnID(adstring str);
            adstring static getSelFcnID(int i);
            
            dvar_vector static calcSelFcn(int id, dvector& z, dvar_vector& params, double fsZ);
            
        public:
            const static int ID_ASCLOGISTIC        = 1; const static adstring STR_ASCLOGISTIC; 
            const static int ID_ASCLOGISTICLN50    = 2; const static adstring STR_ASCLOGISTICLN50; 
            const static int ID_ASCLOGISTIC5095    = 3; const static adstring STR_ASCLOGISTIC5095; 
            const static int ID_ASCLOGISTIC50LN95  = 4; const static adstring STR_ASCLOGISTIC50LN95; 
            const static int ID_ASCLOGISTICLN50LN95= 5; const static adstring STR_ASCLOGISTICLN50LN95; 
            const static int ID_DBLLOGISTIC        = 6; const static adstring STR_DBLLOGISTIC; 
            const static int ID_DBLLOGISTICLND50   = 7; const static adstring STR_DBLLOGISTICLND50; 
            const static int ID_DBLLOGISTIC5095    = 8; const static adstring STR_DBLLOGISTIC5095; 
            const static int ID_DBLLOGISTIC50LN95  = 9; const static adstring STR_DBLLOGISTIC50LN95; 
            const static int ID_DBLLOGISTICLN50LN95=10; const static adstring STR_DBLLOGISTICLN50LN95; 
            const static int ID_ASCLOGISTIC5099    =11; const static adstring STR_ASCLOGISTIC5099; 
            const static int ID_ASCLOGISTIC95LN50  =12; const static adstring STR_ASCLOGISTIC95LN50; 
            const static int ID_ASCNORMAL          =13; const static adstring STR_ASCNORMAL; 
            const static int ID_DBLNORMAL4         =14; const static adstring STR_DBLNORMAL4; 
            const static int ID_DBLNORMAL6         =15; const static adstring STR_DBLNORMAL6; 
            const static int ID_CONSTANT           =16; const static adstring STR_CONSTANT; 
            const static int ID_NONPARAMETRIC      =17; const static adstring STR_NONPARAMETRIC; 
            const static int ID_CUBICSPLINE        =18; const static adstring STR_CUBICSPLINE; 
            const static int ID_DBLNORMAL4A        =19; const static adstring STR_DBLNORMAL4A; 

            /**
             * Calculates "constant" selectivity function (= 1 at all sizes).
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters (IGNORED)
             * @param idZ    - index at which function = 1 (IGNORED) [int]
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static constant(dvector& z);       
            
            /**
             * Calculates "nonparametric" selectivity function with smoothness imposed
             * on the resulting curve by way of penalties in the objective function.
             * 
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of logit-scale parameters, 1 for each size bin
             * @param idZ    - index at which function = 1 (i.e., fully-selected size) [int]
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static nonparametric(dvector& z, dvar_vector& params, int idZ);       
            
            /**
             * Calculates ascending logistic function parameterized by 
             *      params[1]: size at 50% selected (z50)
             *      params[2]: slope
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static asclogistic(dvector& z, dvar_vector& params, double fsZ);       
            
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
            dvar_vector static asclogisticLn50(dvector& z, dvar_vector& params, double fsZ);       
            
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
            dvar_vector static asclogistic5095(dvector& z, dvar_vector& params, double fsZ);
            
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
            dvar_vector static asclogistic5099(dvector& z, dvar_vector& params, double fsZ);
            
            
            /**
             * Calculates ascending logistic function parameterized by 
             *      params[1]: size at 50% selected (z50)
             *      params[2]: ln-scale increment from z50 to z95
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static asclogistic50Ln95(dvector& z, dvar_vector& params, double fsZ);

            /**
             * Calculates ascending logistic function parameterized by 
             *      params[1]: size at 95% selected (z95)
             *      params[2]: ln-scale increment from z50 to z95
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static asclogistic95Ln50(dvector& z, dvar_vector& params, double fsZ);

            /**
             * Calculates ascending logistic function parameterized by 
             *      params[1]: ln-scale size at 50% selected (ln(z50))
             *      params[2]: ln-scale offset to z95
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static asclogisticLn50Ln95(dvector& z, dvar_vector& params, double fsZ);       
            
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
            dvar_vector static dbllogistic(dvector& z, dvar_vector& params, double fsZ);
            
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
            dvar_vector static dbllogisticLnD50(dvector& z, dvar_vector& params, double fsZ);
            
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
            dvar_vector static dbllogistic5095(dvector& z, dvar_vector& params, double fsZ);
            
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
            dvar_vector static dbllogistic50Ln95(dvector& z, dvar_vector& params, double fsZ);
            
            /**
             * Calculates double logistic function parameterized by 
             *      params[1]: log-scale size where ascending limb = 0.5 (ln(z50))
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
            dvar_vector static dbllogisticLn50Ln95(dvector& z, dvar_vector& params, double fsZ);
            
            /**
             * Calculates ascending normal function parameterized by 
             *      params[1]: size at which ascending limb reaches 1
             *      params[2]: width of ascending limb
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double] NOTE: ignored!
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static ascnormal(dvector& z, dvar_vector& params, double fsZ);
            /**
             * Calculates 4-parameter double normal function parameterized by 
             *      params[1]: size at which ascending limb reaches 1
             *      params[2]: width of ascending limb
             *      params[3]: size at which descending limb departs from 1
             *      params[4]: width of descending limb
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double] NOTE: ignored!
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static dblnormal4(dvector& z, dvar_vector& params, double fsZ);
            /**
             * Calculates 6-parameter double normal function parameterized by 
             *      params[1]: size at which ascending limb reaches 1
             *      params[2]: width of ascending limb
             *      params[3]: size at which descending limb departs from 1
             *      params[4]: width of descending limb
             *      params[5]: floor of ascending limb
             *      params[6]: floor of descending limb
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double] NOTE: ignored!
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static dblnormal6(dvector& z, dvar_vector& params, double fsZ);
            /**
             * Calculates an n-parameter cubic spline function parameterized by 
             *      params[1:n]:    knots
             *      params[n+1:2n]: logit-scale values at knots
             * 
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters (knots + logit-scale values)
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double] NOTE: ignored unless > 0
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static cubic_spline(dvector& z, dvar_vector& params, double fsZ);
            /**
             * Calculates 4-parameter double normal function parameterized by 
             *      params[1]: size at which ascending limb reaches 1
             *      params[2]: width of ascending limb
             *      params[3]: scaled increment to params[1] at which descending limb departs from 1
             *      params[4]: width of descending limb
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - max possible size
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static dblnormal4a(dvector& z, dvar_vector& params, double fsZ);
    };
    
#endif	/* MODELSELECTIVITIES_HPP */

