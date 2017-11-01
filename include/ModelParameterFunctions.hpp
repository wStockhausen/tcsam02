/* 
 * File:   ModelParameterFunctions.hpp
 * Author: WilliamStockhausen
 *
 * Created on June 5, 2014, 9:30 AM
 */

#ifndef MODELPARAMETERFUNCTIONS_HPP
#define	MODELPARAMETERFUNCTIONS_HPP
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
     * Set values for a dvar_matrix, intended to function as a vector of devs vector, 
     * based on values of a param_init_bounded_vector_vector.
     * 
     * @param devs - the dvar_matrix
     * @param pDevs - the param_init_bounded_vector_vector
     * @param debug - debugging level
     * @param cout  - output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - This alters the values in the dvar_matrix.
     * 
     */
    void setDevs(dvar_matrix& devs, param_init_bounded_vector_vector& pDevs, int debug, std::ostream& cout);

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
     * Calculate ln-scale priors for a param_init_bounded_vector_vector, 
     * based on its associated BoundedVectorVectorInfo, 
     * and add the weighted NLL to the objective function.
     * 
     * @param objFun - the objective function
     * @param ptrVI  - pointer to the BoundedVectorVectorInfo object associated with pv
     * @param pv     - the param_init_bounded_vector_vector
     * @param debug  - debugging level
     * @param cout   - output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - the value of objFun
     */                                      
    void calcPriors(objective_function_value& objFun, BoundedVectorVectorInfo* ptrVVI, param_init_bounded_vector_vector& pm, int debug, std::ostream& cout);

    /**
     * Calculate ln-scale priors for a dvar_matrix acting as a devs_vector_vector (not defined yet), 
     * based on its associated DevsVectorVectorInfo, 
     * and add the weighted NLL to the objective function.
     * 
     * @param objFun - the objective function
     * @param ptrVI  - pointer to the DevsVectorVectorInfo object associated with pv
     * @param pv     - the dvar_matrix
     * @param debug  - debugging level
     * @param cout   - output stream object for debugging info
     * 
     * @return - void
     * 
     * @alters - the value of objFun
     */                                      
    void calcPriors(objective_function_value& objFun, DevsVectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout);

    void writeParameter(ostream& os, param_init_number& p,             adstring& ctg1, adstring& ctg2, 
                                     NumberInfo* i, int toR, int willBeActive);
    void writeParameter(ostream& os, param_init_bounded_number& p,     adstring& ctg1, adstring& ctg2, 
                                     BoundedNumberInfo* i, int toR, int willBeActive);
    void writeParameters(ostream& os, param_init_vector& p,             adstring& ctg1, adstring& ctg2, 
                                     VectorInfo* i, int toR, int willBeActive);
    void writeParameters(ostream& os, param_init_bounded_vector& p,     adstring& ctg1, adstring& ctg2, 
                                     BoundedVectorInfo* i, int toR, int willBeActive);
    void writeParameters(ostream& os, param_init_bounded_dev_vector& p, adstring& ctg1, adstring& ctg2, 
                                     DevsVectorInfo* i, int toR, int willBeActive);
    void writeParameters(ostream& os, param_init_number_vector& p,        adstring& ctg1, adstring& ctg2, 
                                      NumberVectorInfo* pI, int toR, int willBeActive);
    void writeParameters(ostream& os, param_init_bounded_number_vector& p,adstring& ctg1, adstring& ctg2, 
                                      BoundedNumberVectorInfo* pI, int toR, int willBeActive);
    void writeParameters(ostream& os, param_init_vector_vector& p,        adstring& ctg1, adstring& ctg2, 
                                      VectorVectorInfo* pI, int toR, int willBeActive);
    void writeParameters(ostream& os, param_init_bounded_vector_vector& p,adstring& ctg1, adstring& ctg2, 
                                      BoundedVectorVectorInfo* pI, int toR, int willBeActive);
    
/**
 * Set the parameter info for a NumberVectorInfo object.
 * 
 * @param pNVI - pointer to a NumberVectorInfo instance
 * @param npT - size of vector
 * @param phs - ivector of phases for parameters
 * @param os - output stream to write to
 */
    void setParameterInfo(NumberVectorInfo* pNVI,                           
                          int& npT,
                          ivector& phs, 
                          ostream& os = std::cout);
    /**
     * Set the parameter info for a BoundedNumberVectorInfo object.
     * 
     * @param pBNVI - pointer to a BoundedNumberVectorInfo instance
     * @param npT - size of vector
     * @param lb - dvector of lower bounds
     * @param ub - dvector of upper bounds
     * @param phs - ivector of phases for parameters
     * @param os - output stream to write to
     */
    void setParameterInfo(BoundedNumberVectorInfo* pBNVI,
                          int& npT,
                          dvector& lb, dvector& ub, 
                          ivector& phs, 
                          ostream& os = std::cout);
    /**
     * Set the parameter info for a VectorVectorInfo object.
     * 
     * @param pVVI - pointer to a VectorVectorInfo instance
     * @param npT - size of vector
     * @param lb - dvector of lower bounds
     * @param ub - dvector of upper bounds
     * @param phs - ivector of phases for parameters
     * @param os - output stream to write to
     */
    void setParameterInfo(VectorVectorInfo* pVVI,   
                          int& npT,
                          ivector& mns, ivector& mxs,
                          ivector& phs, 
                          ostream& os = std::cout);
    /**
     * Set the parameter info for a BoundedVectorVectorInfo object.
     * 
     * @param pBVVI - pointer to a BoundedVectorVectorInfo instance
     * @param npT - size of vector
     * @param mns - ivector with minimum indices for each vector
     * @param mxs - ivector with maximum indices for each vector
     * @param idxs - imatrix of reverse indices
     * @param lb - dvector of lower bounds
     * @param ub - dvector of upper bounds
     * @param phs - ivector of phases for parameters
     * @param os - output stream to write to
     */
    void setParameterInfo(BoundedVectorVectorInfo* pBVVI,                           
                          int& npT,
                          ivector& mns, ivector& mxs,
                          imatrix& idxs,
                          dvector& lb, dvector& ub,
                          ivector& phs,
                          ostream& os = std::cout);
    /**
     * Set the parameter info for a DevsVectorVectorInfo object.
     * 
     * @param pDVVI - pointer to a DevsVectorVectorInfo instance
     * @param npT - size of vector
     * @param mns - ivector with minimum indices for each vector
     * @param mxs - ivector with maximum indices for each vector
     * @param idxs - imatrix of reverse indices
     * @param lb - dvector of lower bounds
     * @param ub - dvector of upper bounds
     * @param phs - ivector of phases for parameters
     * @param os - output stream to write to
     */
    void setParameterInfo(DevsVectorVectorInfo* pDVVI,                           
                          int& npT,
                          ivector& mns, ivector& mxs,
                          imatrix& idxs,
                          dvector& lb, dvector& ub,
                          ivector& phs,
                          ostream& os = std::cout);
} //namespace tcsam

#endif	/* MODELPARAMETERFUNCTIONS_HPP */

