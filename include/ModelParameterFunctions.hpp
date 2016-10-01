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

} //namespace tcsam

#endif	/* MODELPARAMETERFUNCTIONS_HPP */

