/* 
 * File:   ModelParameterFunctions.hpp
 * Author: WilliamStockhausen
 *
 * Created on June 5, 2014, 9:30 AM
 */

#ifndef MODELPARAMETERFUNCTIONS_HPP
#define	MODELPARAMETERFUNCTIONS_HPP
#include <admodel.h>
#include "ModelParameterInfoTypes.hpp"
#include "ModelParameterVectorInfoTypes.hpp"

namespace tcsam {
    /******************************************************************************
    * Set initial values for a param_init_bounded_number_vector.
    * 
    * Description: Sets initial values for a parameter vector.
    * 
    * Inputs:
    *  @param pI : pointer to BoundedNumberVectorInfo object
    *  @param p  : reference to a param_init_bounded_number_vector
    *     
    * @return void
     * 
    * @alters - changes initial values of p
     * 
    ******************************************************************************/
    void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, std::ostream& cout);

    /******************************************************************************
    * Set initial values for a param_init_bounded_vector_vector.
    * 
    * Description: Sets initial values for a vector of parameter vectors.
    * 
    * Inputs:
    *  @param pI : pointer to a BoundedVectorVectorInfo object
    *  @param p  : reference to a param_init_bounded_vector_vector
    *     
    * @return void
     * 
    * @alters - changes initial values of p
     * 
    ******************************************************************************/
    void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout);

    /******************************************************************************
    * Set initial devs values for a param_init_bounded_vector_vector.
    * 
    * Description: Sets initial values for a vector of devs vectors.
    * 
    * Inputs:
    *  @param pI : pointer to a DevsVectorVectorInfo object
    *  @param p  : reference to a param_init_bounded_vector_vector
    *     
    * @return void
     * 
    * @alters - changes initial values of p
     * 
    ******************************************************************************/
    void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout);

    /**
     * Set values for a dvar_matrix, intended to function as a vector of devs vectors, 
     * based on values of a param_init_bounded_vector_vector.
     * 
     * @param devs - the dvar_matrix
     * @param pDevs - the param_init_bounded_number_vector
     * @param pI    - pointer to the associated DevsVectorVectorInfo object
     * @param debug - debugging level
     * @param cout  - output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - This alters the values in the devs dvar_matrix.
     * 
     */
    void setDevsVectorVector(dvar_matrix& devs, param_init_bounded_number_vector& pDevs,  DevsVectorVectorInfo* pI, int debug, std::ostream& cout);
    
    /**
     * Set values for a dvar_matrix, intended to function as a vector of bounded parameter vectors, 
     * based on values of a param_init_bounded_number_vector.
     * 
     * @param mat   [output]  - the dvar_matrix
     * @param pDevs [input] - the param_init_bounded_number_vector
     * @param pI    [input] - pointer to the associated DevsVectorVectorInfo object
     * @param debug [input] - debugging level
     * @param os  - [input] output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - This alters the values in the mat dvar_matrix.
     * 
     */
    void setBoundedVectorVector(dvar_matrix& mat, param_init_bounded_number_vector& pVals, BoundedVectorVectorInfo* pI, int debug, std::ostream& os);

    /**
     * Calculate ln-scale priors for a param_init_number_vector, based on its associated NumberVectorInfo, 
     * and add the weighted NLL to the objective function.
     * 
     * @param objFun - the objective function
     * @param ptrVI  - pointer to the NumberVectorInfo object associated with pv
     * @param pv     - the param_init_number_vector
     * @param debug  - debugging level
     * @param cout   - output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - the value of objFun
     */                                      
    void calcPriors(objective_function_value& objFun, NumberVectorInfo* ptrVI,param_init_number_vector& pv, int debug, std::ostream& cout);

    /**
     * Calculate ln-scale priors for a param_init_bounded_number_vector, 
     * based on its associated BoundedNumberVectorInfo, 
     * and add the weighted NLL to the objective function.
     * 
     * @param objFun - the objective function
     * @param ptrVI  - pointer to the BoundedNumberVectorInfo object associated with pv
     * @param pv     - the param_init_bounded_number_vector
     * @param debug  - debugging level
     * @param cout   - output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - the value of objFun
     */                                      
    void calcPriors(objective_function_value& objFun, BoundedNumberVectorInfo* ptrVI,param_init_bounded_number_vector& pv, int debug, std::ostream& cout);

    /**
     * Calculate ln-scale priors for a dvar_matrix acting as a vector_vector, 
     * based on its associated VectorVectorInfo, 
     * and add the weighted NLL to the objective function.
     * 
     * Note that this pertains to VectorVectorInfo and subclasses BoundedVectorVectorInfo and 
     * DevsVectoroVectorInfo.
     * 
     * @param objFun - the objective function
     * @param ptrVI  - pointer to the VectorVectorInfo object associated with pv
     * @param pv     - the dvar_matrix
     * @param debug  - debugging level
     * @param cout   - output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - the value of objFun
     */                                      
    void calcPriors(objective_function_value& objFun, VectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout);

//    /**
//     * Calculate ln-scale priors for a param_init_bounded_vector_vector, 
//     * based on its associated BoundedVectorVectorInfo, 
//     * and add the weighted NLL to the objective function.
//     * 
//     * @param objFun - the objective function
//     * @param ptrVI  - pointer to the BoundedVectorVectorInfo object associated with pv
//     * @param pv     - the param_init_bounded_vector_vector
//     * @param debug  - debugging level
//     * @param cout   - output stream object for debugging info
//     * 
//     * @return - void
//     * 
//     * @alters - the value of objFun
//     */                                      
//    void calcPriors(objective_function_value& objFun, BoundedVectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout);

//    /**
//     * Calculate ln-scale priors for a dvar_matrix acting as a devs_vector_vector (not defined yet), 
//     * based on its associated DevsVectorVectorInfo, 
//     * and add the weighted NLL to the objective function.
//     * 
//     * @param objFun - the objective function
//     * @param ptrVI  - pointer to the DevsVectorVectorInfo object associated with pv
//     * @param pv     - the dvar_matrix
//     * @param debug  - debugging level
//     * @param cout   - output stream object for debugging info
//     * 
//     * @return - void
//     * 
//     * @alters - the value of objFun
//     */                                      
//    void calcPriors(objective_function_value& objFun, DevsVectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout);

    /**
     * Writes the header (column names) to a stream which subsequent writeParameters functions
     * will write to using csv format.
     * 
     * @param os - the output stream
     */
    void writeCSVHeaderForParametersFile(ostream& os);
/**
 * Writes information for a parameter to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter (param_init_number)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated NumberInfo object
 * @param toR - flag to write to R format (=1) or csv (=0)
 * @param willBeActive - flag to write only if parameter will be active
 */
    void writeParameter(ostream& os, 
                        param_init_number& p,             
                        adstring& ctg1, adstring& ctg2, 
                        NumberInfo* i, 
                        int toR, int willBeActive);
/**
 * Writes information for a bounded parameter to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter (param_init_bounded_number)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - BoundedNumberInfo
 * @param toR - flag to write to R format (=1) or csv (=0)
 * @param willBeActive - flag to write only if parameter will be active
 */
    void writeParameter(ostream& os, 
                        param_init_bounded_number& p,     
                        adstring& ctg1, adstring& ctg2, 
                        BoundedNumberInfo* i, 
                        int toR, int willBeActive);
/**
 * Writes information for a parameter vector to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter vector (as dvar_vector)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param lbl - parameter-specific label
 * @param toR - flag to write to R format (=1) or csv (=0)
 * @param willBeActive - flag to write only if parameter vector will be active
 */
    void writeParameters(ostream& os, 
                         dvar_vector& p,             
                         adstring& ctg1, adstring& ctg2, 
                         VectorInfo* i, 
                         int toR, int willBeActive);
/**
 * Writes information for a bounded parameter vector to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter vector (as dvar_vector)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedVectorInfo info object
 * @param toR - flag to write to R format (=1) or csv (=0)
 * @param willBeActive - flag to write only if parameter vector will be active
 */
    void writeParameters(ostream& os, 
                         dvar_vector& p,     
                         adstring& ctg1, adstring& ctg2, 
                         BoundedVectorInfo* i, 
                         int toR, int willBeActive);
/**
 * Writes information for a devs parameter vector to an output stream.
 * 
 * @param os - output stream
 * @param p - dvar_vector with devs parameter vector values
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated DevsVectorInfo info object
 * @param toR - flag to write to R format (=1) or csv (=0)
 * @param willBeActive - flag to write only if parameter vector will be active
 */
    void writeParameters(ostream& os, 
                         dvar_vector& p, 
                         adstring& ctg1, adstring& ctg2, 
                         DevsVectorInfo* i, 
                         int toR, int willBeActive);
/**
 * Writes a vector of parameters (param_init_number_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_number_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated NumberVectorInfo info object
 * @param toR - flag to write to R format (=1) or csv (=0)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
    void writeParameters(ostream& os, 
                         param_init_number_vector& p,        
                         adstring& ctg1, adstring& ctg2, 
                         NumberVectorInfo* pI, 
                         int toR, int willBeActive);
/**
 * Writes a vector of parameters (param_init_bounded_number_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_bounded_number_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedNumberVectorInfo info object
 * @param toR - flag to write to R format (=1) or csv (=0)
* @param willBeActive - flag to write if parameters will be active in some phase
 */
    void writeParameters(ostream& os, 
                         param_init_bounded_number_vector& p,
                         adstring& ctg1, adstring& ctg2, 
                         BoundedNumberVectorInfo* pI, 
                         int toR, int willBeActive);
/**
 * Writes a vector of parameter vectors (as a dvar_matrix) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - dvar_matrix representing the vector_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated VectorVectorInfo info object
 * @param toR - flag to write to R format (=1) or csv (=0)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
    void writeParameters(ostream& os, 
                         dvar_matrix& p,        
                         adstring& ctg1, adstring& ctg2, 
                         VectorVectorInfo* pI, 
                         int toR, int willBeActive);
/**
 * Writes a vector of parameter vectors (as a dvar_matrix) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - dvar_matrix representing the bounded_vector_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedVectorVectorInfo info object
 * @param toR - flag to write to R format (otherwise csv)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
    void writeParameters(ostream& os, 
                         dvar_matrix& p,
                         adstring& ctg1, adstring& ctg2, 
                         BoundedVectorVectorInfo* pI, 
                         int toR, int willBeActive);
    
/**
 * Writes a vector of parameter vectors (as a dvar_matrix) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a dvar_matrix
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated DevsVectorVectorInfo info object
 * @param toR - flag to write to R format (otherwise csv)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
    void writeParameters(ostream& os, 
                         dvar_matrix& p,
                         adstring& ctg1, 
                         adstring& ctg2, 
                         DevsVectorVectorInfo* pI, 
                         int toR, int willBeActive);
    
/**
 * Sets the info for a param_init_number_vector from a NumberVectorInfo object.
 * 
 * @param pNVI - pointer to a NumberVectorInfo instance
 * @param npT [out] - size of vector
 * @param phs [out] - ivector of phases for parameters
 * @param os  - output stream to write processing info to
 */
    void setParameterInfo(NumberVectorInfo* pNVI,                           
                          int& npT,
                          ivector& phs, 
                          ostream& os = std::cout);
/**
 * Sets the info for a param_init_bounded_number_vector from a BoundedNumberVectorInfo object.
 * 
 * @param pBNVI - pointer to a BoundedNumberVectorInfo instance
 * @param npT [out]  - size of vector
 * @param lb  [out]  - dvector of lower bounds on parameter scale
 * @param ub  [out]  - dvector of upper bounds on parameter scale
 * @param phs [out] - ivector of phases for parameters
 * @param os - output stream to write to
 */
    void setParameterInfo(BoundedNumberVectorInfo* pBNVI,
                          int& npT,
                          dvector& lb, 
                          dvector& ub, 
                          ivector& phs, 
                          ostream& os = std::cout);
    /**
        * Sets the info for a param_init_number_vector substituting as a 
        * param_init_vector_vector from a VectorVectorInfo object.
     * 
     * @param pVVI - pointer to a VectorVectorInfo instance
     * @param npV [out] - number of vectors represented
     * @param npT [out] - total size (number of elements)
     * @param mns [out] - ivector with minimum indices for each vector
     * @param mxs [out] - ivector with maximum indices for each vector
     * @param idxs [out] - imatrix of reverse indices
     * @param phs [out] - ivector of phases for parameters
     * @param os - output stream to write to
     */
    void setParameterInfo(VectorVectorInfo* pVVI,   
                          int& npV,
                          int& npT,
                          ivector& mns, 
                          ivector& mxs,
                          imatrix& idxs,
                          ivector& phs, 
                          ostream& os = std::cout);
    /**
    * Sets the info for a param_init_bounded_number_vector substituting as a 
    * param_init_bounded_vector_vector from a BoundedVectorVectorInfo object.
     * 
     * @param pBVVI - pointer to a BoundedVectorVectorInfo instance
     * @param npV [out] - number of devs vectors represented
     * @param npT [out] - total size (number of elements)
     * @param mns [out] - ivector with minimum indices for each vector
     * @param mxs [out] - ivector with maximum indices for each vector
     * @param idxs [out] - imatrix of reverse indices
     * @param lb [out] - dvector of lower bounds
     * @param ub [out] - dvector of upper bounds
     * @param phs [out] - ivector of phases for parameters
     * @param os - output stream to write to
     */
    void setParameterInfo(BoundedVectorVectorInfo* pBVVI,                           
                            int& npV,
                            int& npT,
                            ivector& mns, 
                            ivector& mxs,
                            imatrix& idxs,
                            dvector& lb, 
                            dvector& ub, 
                            ivector& phs, 
                            ostream& os = std::cout);
/**
 * Sets the info for a param_init_bounded_vector_vector acting as a devs_vector_vector
 * from a DevsVectorVectorInfo object.
 * 
 * @param pDVVI - pointer to a DevsVectorVectorInfo instance
 * @param npV [out] - number of devs vectors represented
 * @param npT [out] - total size (number of elements)
 * @param mns [out]  - ivector(1,npV) with minimum indices for each vector
 * @param mxs [out]  - ivector(1,npV) with maximum indices for each vector
 * @param idxs [out] - imatrix of reverse indices
 * @param lb [out]   - dvector(1,npV) of lower bounds
 * @param ub [out]   - dvector(1,npV) of upper bounds
 * @param phs [out]  - ivector(1,npT) of phases for parameters
 * @param os - output stream to write to
 */
    void setParameterInfo(DevsVectorVectorInfo* pDVVI, 
                          int& npV,
                          int& npT,
                          ivector& mns, 
                          ivector& mxs,
                          imatrix& idxs,
                          dvector& lb, 
                          dvector& ub,
                          ivector& phs,
                          ostream& os = std::cout);
} //namespace tcsam

#endif	/* MODELPARAMETERFUNCTIONS_HPP */

