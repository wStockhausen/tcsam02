/* 
 * File:   ModelSelectivities.hpp
 * Author: William.Stockhausen
 *
 * Created on August 13, 2012, 9:22 AM
 * 
 * 2014-12-05: 1. Added asclogisticLn50Ln95 and dbllogisticLn50Ln95 functions
 * 2015-04-17: 1. Added dbllogistic50Ln95 function
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
            const static int ID_ASCLOGISTIC        =1; const static adstring STR_ASCLOGISTIC; 
            const static int ID_ASCLOGISTIC5095    =2; const static adstring STR_ASCLOGISTIC5095; 
            const static int ID_ASCLOGISTIC50LN95  =3; const static adstring STR_ASCLOGISTIC50LN95; 
            const static int ID_ASCLOGISTICLN50LN95=4; const static adstring STR_ASCLOGISTICLN50LN95; 
            const static int ID_DBLLOGISTIC        =5; const static adstring STR_DBLLOGISTIC; 
            const static int ID_DBLLOGISTIC5095    =6; const static adstring STR_DBLLOGISTIC5095; 
            const static int ID_DBLLOGISTIC50LN95  =7; const static adstring STR_DBLLOGISTIC50LN95; 
            const static int ID_DBLLOGISTICLN50LN95=8; const static adstring STR_DBLLOGISTICLN50LN95; 
            const static int ID_DBLNORMAL          =9; const static adstring STR_DBLNORMAL; 

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
             * Calculates double normal function parameterized by 
             *      params[1]: 
             *      params[2]: 
             *      params[3]:
             *      params[4]:
             *      params[5]:
             *      params[6]:
             * Inputs:
             * @param z      - dvector of sizes at which to compute function values
             * @param params - dvar_vector of function parameters
             * @param fsZ    - size at which function = 1 (i.e., fully-selected size) [double]
             * 
             * @return - selectivity function values as dvar_vector
             */
            dvar_vector static dblnormal(dvector& z, dvar_vector& params, double fsZ);
    };
    
#endif	/* MODELSELECTIVITIES_HPP */

