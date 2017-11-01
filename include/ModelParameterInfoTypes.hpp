/**
 * @file
 * Includes:
 * <ul>
 * <li> NumberInfo
 * <li> BoundedNumberInfo
 * <li> NumberVectorInfo
 * <li> BoundedNumberVectorInfo
 * <li> VectorInfo
 * <li> VectorVectorInfo
 * <li> BoundedVectorInfo
 * <li> BoundedVectorVectorInfo
 * <li> DevsVectorInfo
 * <li> DevsVectorVectorInfo
 * </ul>
 */
#pragma once
#ifndef MODELPARAMETERINFOTYPES_HPP
    #define MODELPARAMETERINFOTYPES_HPP

    #include <admodel.h>
    #include "ModelIndexBlocks.hpp"
    #include "ModelFunctions.hpp"
    #include "ModelConfiguration.hpp"
#include "ModelConstants.hpp"

/**
 * NumberInfo
 * 
 *  Encapsulates characteristics for a param_init_number parameter,
 *  with values assumed to be on an arithmetic scale (i.e., no transformations
 *  are made between the parameter value and the NumberInfo's initVal and finlVal).
 * 
 * Class fields:
 * <ul>
 * <li> scaleType   - parameter scale type (assumed arithmetic for NumberInfo objects)
 * <li> initVal     - initial value for associated parameter
 * <li> finlVal     - final value for associated parameter
 * <li> phase       - phase to start estimation of associated parameter
 * <li> resample    - flag to resample initial value from prior
 * <li> priorWgt    - weight assigned to prior in likelihood
 * <li> priorType   - name identifying pdf for prior
 * <li> priorParams - dvar_vector of parameters for prior
 * <li> priorConsts - dvector of constants for prior
 * <li> name        - name of associated parameter
 * <li> label       - label for associated parameter
 * </ul>
 */
    class NumberInfo {
        public:
            /* flag to turn on debugging output */
            static int debug;
        protected:
            /* flag indicating parameter scale */
            int scaleType; 
            /* arithmetic-scale value used as the initial value for the associated parameter */
            double initVal;
            /* final arithmetic-scale value from the associated parameter (for output to R) */
            double finlVal;     
            /* phase in which to start estimating associated parameter */
            int phase;          
            /* weight to assign to prior probability */
            double priorWgt;
            /* pdf type for prior */
            adstring priorType;
            /* parameters for prior, specified on arithmetic scale */
            dvector priorParams;
            /* constants for prior, specified on arithmetic scale */
            dvector priorConsts;//specified on parameter scale
        public:
            /* the parameter name */
            adstring name;
            /* the label for the associated parameter */
            adstring label;
            /* flag to do resampling of initial values */
            bool resample;
            /* pointer to info for resmapling pdf */
            ModelPDFInfo*  pMPI;
        public:
            /**
             * Class constructor.
             * 
             * Sets name to "" and the pointer pMPI to 0.
             */
            NumberInfo(){this->name="";pMPI=0;scaleType=-1;finlVal=0.0;}
            /**
             * Class constructor.
             * 
             * Also sets the pointer pMPI to 0.
             * 
             * @param name - name (as adstring&) to use for associated parameter
             */
            NumberInfo(adstring& name){this->name=name;pMPI=0;scaleType=-1;finlVal=0.0;}
            /**
             * Class constructor.
             * 
             * Also sets the pointer pMPI to 0.
             * 
             * @param name - name (as const char*) to use for associated parameter
             */
            NumberInfo(const char* name){this->name=name;pMPI=0;scaleType=-1;finlVal=0.0;}
            /**
             * Class destructor.
             * 
             * Deletes pointer pMPI and sets it to 0.
             */
            ~NumberInfo(){delete pMPI;pMPI=0;}
            
            /**
             * Gets the parameter scale type, as an adstring.
             * 
             * @return the scale type
             */
            adstring getScaleType(){return tcsam::getScaleType(scaleType);}
            /**
             * Calculates the arithmetic-scale value corresponding to a given parameter-scale value.
             * 
             * @param x - the parameter-scale value
             * @return - the corresponding arithmetic-scale value
             */
            virtual double calcArithScaleVal(double x);
            /**
             * Calculates the arithmetic-scale value corresponding to given parameter-scale value.
             * 
             * @param x - parameter-scale value
             * @return - arithmetic-scale value
             */
            virtual dvariable calcArithScaleVal(const dvariable& x);
            /**
             * Calculates parameter-scale value corresponding to given bounded arithmetic-scale value.
             * 
             * Use to calculate the initial value for the parameter
             * 
             * @param x - arithmetic-scale value
             * @return - parameter-scale value
             */
            virtual double calcParamScaleVal(double x);
            
            /**
             * Gets the phase in which parameter estimation is started.
             * 
             * @return - the phase
             */
            int getPhase(){return phase;}
            /**
             * Gets the multiplicative weight set on the prior probability in the likelihood.
             * 
             * @return - the weight
             */
            double getPriorWgt(){return priorWgt;}
            /**
             * Gets the prior type, as an adstring.
             * 
             * @return - the prior type, as an adstring object 
             */
            adstring getPriorType(){return priorType;}
            /**
             * Calculates the log-scale prior probability based on an input value.
             * 
             * The input value is presumably from the associated parameter value.
             * 
             * @param x - the value from the associated parameter
             * 
             * @return - the log-scale prior, as a dvariable
             */
            dvariable calcLogPrior(prevariable& x);             
            /**
             * Gets the value to use as the initial value for the associated parameter.
             * 
             * @return - the value, as a double
             */
            double getInitVal(){return initVal;}
            /**
             * Sets the value used as the initial value for the associated parameter.
             * <p>
             * The value given is assumed to be on the arithmetic scale.
             * 
             * @param x - value to set
             */
            virtual void setInitVal(double x){initVal=x;}
            /**
             * Sets the value used as the initial value, from the associated parameter.
             * 
             * Use this method to update \code initVal \endcode in the case a pin file is used to set
             * initial parameter values.
             * 
             * @param x - the associated parameter (an instance of a param_init_number)
             */
            virtual void setInitVal(param_init_number& x){initVal=value(x);}
            /**
             * Draws a random value based on the specified prior probability distribution.
             * 
             * Note that this DOES NOT update initVal.
             * If the estimation phase is \< 0, the value of initVal is returned
             * 
             * @param rng - the random number generator
             * @param vif - the variance inflation factor
             * 
             * @return - the random number, or the value of initVal
             */
            virtual double drawInitVal(random_number_generator& rng, double vif);
            /**
             * Gets the final value from the associated parameter.
             * 
             * @param x - the associated parameter (an instance of a param_init_number)
             */
            virtual double getFinalVal(){return finlVal;}
            /**
             * Sets the final value from the associated parameter.
             * 
             * @param x - the associated parameter (an instance of a param_init_number)
             */
            virtual void setFinalVal(param_init_number& x){finlVal=value(x);}
            /**
             * Reads the parameter info in ADMB format from an input stream.
             * 
             * @param is - the input stream
             */
            virtual void read(cifstream & is);
            /**
             * Writes the parameter info in ADMB format to an output stream.
             * 
             * @param os - the output stream
             */
            virtual void write(std::ostream & os);
            /**
             * Writes the parameter info as part of an R list to an output stream.
             * 
             * @param os - the output stream.
             */
            virtual void writeToR(std::ostream & os);
            
            friend cifstream& operator >>(cifstream & is, NumberInfo & obj){obj.read(is); return is;}
            friend std::ostream& operator <<(std::ostream & os, NumberInfo & obj){obj.write(os); return os;}
        protected:
            /**
             * Writes the NumberInfo part of the parameter info as part of an R list to an output stream.
             * 
             * @param os - the output stream.
             */
            virtual void writeToR1(std::ostream & os);
            /**
             * Sets the prior type, based on an adstring value.
             * 
             * @param prior - the prior type, as an adstring
             */
            virtual void setPriorType(adstring & prior);
    };

/**
 * BoundedNumberInfo : NumberInfo
 * 
 *  Encapsulates characteristics for a param_init_bounded_number parameter on a 
 *  possibly-transformed (logit, lognormal, probit) scale.
 * 
 * Additional class fields:
 * <ul>
 * <li> jitter    - integer flag indicating whether or not to jitter initial values
 * <li> lower     - lower bound on arithmetic scale
 * <li> upper     - upper bound on arithmetic scale
 * </ul>
 * 
 * Note that the priorParams and priorConsts for this class are interpreted on the
 * arithmetic scale, not on the (possibly transformed) scale of the associated parameter.
 */
    class BoundedNumberInfo : public NumberInfo {
        public:
            static int debug;
        protected:
            double lower;  //lower bound
            double upper;  //upper bound
        public:
            int jitter;  //flag to jitter
        public:
            /**
             * Class constructor.
             * 
             * Sets name to "" and the pointer pMPI to 0.
             */
            BoundedNumberInfo():NumberInfo(){}
            /**
             * Class constructor.
             * 
             * Also sets the pointer pMPI to 0.
             * 
             * @param name - name (as adstring&) to use for associated parameter
             */
            BoundedNumberInfo(adstring& name):NumberInfo(name){}
            /**
             * Class constructor.
             * 
             * Also sets the pointer pMPI to 0.
             * 
             * @param name - name (as const char*) to use for associated parameter
             */
            BoundedNumberInfo(const char * name):NumberInfo(name){}
            /**
             * Class destructor.
             * 
             * Deletes pointer pMPI and sets it to 0 via the NumberInfo destructor.
             */
            ~BoundedNumberInfo(){}
            
            /**
             * Returns the parameter-scale lower bound.
             * 
             * @return 
             */            
            double getLowerBound(){return this->calcParamScaleVal(lower);}
            /**
             * Returns the parameter-scale upper bound.
             * 
             * @return 
             */            
            double getUpperBound(){return this->calcParamScaleVal(upper);}           
            /**
             * Calculates the arithmetic-scale value corresponding to a given parameter-scale value.
             * 
             * @param x - the parameter-scale value
             * @return - the corresponding arithmetic-scale value
             */
            virtual double calcArithScaleVal(double x);
            /**
             * Calculates the arithmetic-scale value corresponding to given parameter-scale value.
             * 
             * @param x - parameter-scale value
             * @return - arithmetic-scale value
             */
            virtual dvariable calcArithScaleVal(const dvariable& x);
            /**
             * Calculates parameter-scale value corresponding to given bounded arithmetic-scale value.
             * 
             * Use to calculate the initial value for the parameter
             * 
             * @param x - arithmetic-scale value
             * @return - parameter-scale value
             */
            virtual double calcParamScaleVal(double x);
            /**
             * Sets the initial value to be used for the associated parameter, assuming
             * the value given is on the arithmetic scale.
             * The method ensures the initial value is set within the prescribed bounds.
             * 
             * @param x - the initial value to set
             */
            virtual void setInitVal(double x);           
            /**
             * Sets the initial value on the arithmetic scale, based on a 
             * bounded parameter on a (possibly) non-arithmetic scale.
             * 
             * @param x - the parameter-scale value
             */
            virtual void setInitVal(param_init_bounded_number& x){initVal=calcArithScaleVal(value(x));}
            /**
             * Set the final value on the arithmetic scale, based on a 
             * bounded parameter on a (possibly) non-arithmetic scale.
             * 
             * @param x - the bounded parameter (an instance of a param_init_bounded_number)
             */
            virtual void setFinalVal(param_init_bounded_number& x){finlVal=calcArithScaleVal(value(x));}            
            /**
             * Reads the parameter info in ADMB format from an input stream.
             * 
             * The read order for the parameter info is:
             * <ul>
             *  <li> scaleType (as an adstring)
             *  <li> lower bound, on arithmetic scale
             *  <li> upper bound, on arithmetic scale
             *  <li> initVal
             *  <li> phase
             *  <li> resample
             *  <li> priorWgt
             *  <li> priorType
             *  <li> priorParams
             *  <li> priorConsts
             *  <li> label
             * </ul>
             * 
             * @param is - the input stream
             * 
             */
            virtual void read(cifstream & is);
            /**
             * Writes the parameter info in ADMB format to an output stream.
             * 
             * @param os - the output stream
             */
            virtual void write(std::ostream & os);
            /**
             * Writes the parameter info as part of an R list to an output stream.
             * 
             * @param os - the output stream.
             */
            virtual void writeToR(std::ostream& os);
        protected:
            /**
             * Writes the BoundedNumberInfo part of the parameter info as part of an R list to an output stream.
             * 
             * @param os - the output stream.
             */
            virtual void writeToR1(std::ostream & os);
    };//BoundedNumberInfo

/**
 * VectorInfo : NumberInfo
 * 
 *  Encapsulates characteristics for a param_init_vector parameter,
 *  with values assumed to be on an arithmetic scale (i.e., no transformations
 *  are made between the parameter value and the VectorInfo's initVals and finlVals.
 * 
 * Additional public class fields:
 * <ul>
 * <li> N        -  size of parameter vector
 * <li> readVals -  flag to read vector of values from input stream
 * </ul>
 * Additional protected class fields:
 * <ul>
 * <li> idxType  - adstring with type of indices (e.q., YEAR)
 * <li> ptrIB    - pointer to IndexBlock object defining indices
 * <li> initVals - vector of initial values
 * <li> finlVals - vector of final values
 * </ul>
 * 
 * Note that the priorParams and priorConsts for this class are interpreted on the
 * arithmetic scale.
 */
    class VectorInfo: public NumberInfo {
        public:
            /* flag to print debugging info */
            static int debug;
            /* size of initVals vector*/
            int N;
            /* flag to read initial values */
            int readVals;
        protected:
            /* index type (ie., model dimension, e.g., YEAR) */
            adstring idxType;
            /* pointer to IndexBlock defining vector */
            IndexBlock* ptrIB; 
            /* vector of initial values */
            dvector initVals;
            /* final values (for output to R) */
            dvector finlVals;
        public:
            /**
             * Class constructor.
             * 
             * Sets pointer ptrIB = 0.
             */
            VectorInfo():NumberInfo(){ptrIB=0;}
            /**
             * Class constructor.
             * 
             * @param [in] name - name of associated param_init_vector
             */
            VectorInfo(adstring& name):NumberInfo(name){ptrIB=0;}
            /**
             * Class constructor.
             * 
             * @param [in] name - name of associate param_init_vector
             */
            VectorInfo(const char* name):NumberInfo(name){ptrIB=0;}
            /**
             * Class destructor.
             */
            ~VectorInfo(){if (ptrIB) delete ptrIB;}
            
            /**
             * Gets the size of the vector (indices run 1:N).
             * 
             * @return 
             */
            int getSize(){return N;}
            /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale value.
             * 
             * @param x - parameter-scale value as a double
             * @return - arithmetic-scale value as a double
             */
            virtual double calcArithScaleVal(double x){return NumberInfo::calcArithScaleVal(x);}
            /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
             * 
             * @param x - parameter-scale values as a dvector
             * @return - arithmetic-scale values as a dvector
             */
            dvector calcArithScaleVal(const dvector& x);
            /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
             * 
             * @param x - parameter-scale values as dvar_vector
             * @return - arithmetic-scale values as dvar_vector
             */
            dvar_vector calcArithScaleVal(const dvar_vector& x);
            /**
             * Calculates parameter-scale values corresponding to the input arithmetic-scale values.
             * 
             * @param x - the arithmetic-scale values as dvector
             * @return - the parameter-scale values as dvector
             */
            dvector calcParamScaleVal(dvector& x);
            /**
             * Gets an ivector of model indices corresponding to integer indices 1:N.
             * 
             * @return - ivector iv(1:N) such that iv(i) is the model index value 
             *           corresponding to the ith element.
             */
            ivector getFwdIndices(){return ptrIB->getFwdIndexVector();}
            /**
             * Gets an ivector of integer indices 1:N corresponding to model indices.
             * 
             * @return - ivector iv(modMin:modMax) such that iv(i) is the index into
             *           the block corresponding to model index value i. iv(i) is 0
             *           for model index values i that do not correspond to a block
             *           element.
             */
            ivector getRevIndices(){return ptrIB->getRevIndexVector();}
            /**
             * Sets initial values to input constant.
             * 
             * @param [in] x - the value to set
             */
            void setInitVal(double x){initVals = x;}
            /**
             * Gets the vector of initial values. 
             * Vector indices run 1:N.
             * 
             * @return - dvector of initial values
             */
            dvector getInitVals(){return initVals;}
            /**
             * Sets the vector of initial values to the input dvector element-by-element.
             * 
             * @param [in] x - a dvector of initial values to set. Indices should run 1:N
             */
            virtual void setInitVals(dvector& x){initVals=x;}     
            /**
             * Sets the vector of initial values to the input param_init_vector element-by-element.
             * 
             * @param [in] x - a param_init_vector whose values will be copied as the initial values. Indices should run 1:N
             */
            virtual void setInitVals(param_init_vector & x){initVals=value(x);}     
            /**
             * Gets a copy of the vector of final values.
             * 
             * @return a dvector of final values
             */
            virtual dvector getFinalVals(){return finlVals;}
            /**
             * Sets final values for output to R.
             * 
             * @param [in] x - a param_init_vector whose values will be copied as the final values.
             */
            virtual void setFinalVals(param_init_vector & x){finlVals.allocate(x.indexmin(),x.indexmax()); finlVals=value(x);}     
            /**
             * Calculates a dvar_vector of log prior values based on an input vector.
             * 
             * @return - a dvar_vector
             */
            dvar_vector calcLogPrior(dvar_vector& x);
            /**
             * Draw initial values based on jittering or resampling the prior.
             * 
             * If the parameter estimation phase is \< 0, a copy of the vector of initial 
             * values (initVals) is returned.
             * 
             * @param [in] rng - random number generator object
             * @param [in] vif - variance inflation factor
             * 
             * @return - a dvector of random (or initial) values. Indices run 1:N
             */
            virtual dvector drawInitVals(random_number_generator& rng, double vif);//draw initial values by resampling prior
            /**
             * Reads a vector of initial values from an input filestream.
             * 
             * @param [in] is - the input filestream
             */
            virtual void readInitVals(cifstream& is){is>>initVals;}
            /**
             * Reads the parameter info from an input filestream.
             * The read order is:
             * <ul>
             * <li> idxType - the index type (as an adstring)
             * <li> the index block defining the vector indices (which determines N), as an IndexBlock
             * <li> readVals - an adstring flag to subsequently read a vector of initial values
             *  <li> initVal
             *  <li> phase
             *  <li> resample
             *  <li> priorWgt
             *  <li> priorType
             *  <li> priorParams
             *  <li> priorConsts
             *  <li> label
             * <\ul>
             * 
             * @param [in] is - the filestream to read from
             */
            virtual void read(cifstream & is);
            /**
             * Writes the parameter info to an output stream, in ADMB format.
             * 
             * @param [in] os - the output stream to write to
             */
            virtual void write(std::ostream & os);
            /**
             * Writes the parameter info to an output stream as an R list.
             * 
             * @param [in] os - the output stream to write to
             */
            virtual void writeToR(std::ostream & os);
            /**
             * Writes the final values to an output stream as an R structure.
             * 
             * @param [in] os - the output stream to write to
             */
            virtual void writeFinalValsToR(std::ostream & os);
    };

/**
 * BoundedVectorInfo : BoundedNumberInfo
 * 
 *  Encapsulates characteristics for a param_init_bounded_number_vector 
 *  parameter vector on a possibly-transformed (logit, lognormal, probit) scale.
 *  Basically, a bounded parameter version of VectorInfo.
 * 
 * Additional public class fields:
 * <ul>
 * <li> N        -  size of parameter vector
 * <li> readVals -  flag to read vector of values from input stream
 * </ul>
 * Additional protected class fields:
 * <ul>
 * <li> idxType  - adstring with type of indices (e.q., YEAR)
 * <li> ptrIB    - pointer to IndexBlock object defining indices
 * <li> initVals - vector of initial values
 * <li> finlVals - vector of final values
 * </ul>
 * 
 */
    class BoundedVectorInfo : public BoundedNumberInfo {
        public:
            /* flag to print debugging info */
            static int debug;
            /* size of initVals vector*/
            int N;
            /* flag to read initial values */
            int readVals;
        protected:
            /* index type (ie., model dimension, e.g., YEAR) */
            adstring idxType;
            /* pointer to IndexBlock defining vector */
            IndexBlock* ptrIB; 
            /* vector of initial values */
            dvector initVals;
            /* final values (for output to R) */
            dvector finlVals;
        public:
            /**
             * Class constructor.
             */
            BoundedVectorInfo():BoundedNumberInfo(){}
            /**
             * Class constructor.
             * 
             * @param [in] name - adstring for name of associated param_init_bo
             */
            BoundedVectorInfo(adstring& name):BoundedNumberInfo(name){}
            /**
             * Class constructor.
             * 
             * @param [in] name - adstring for name of associated param_init_bo
             */
            BoundedVectorInfo(const char * name):BoundedNumberInfo(name){}
            /**
             * Class destructor.
             */
            ~BoundedVectorInfo(){}
                                  
            /**
             * Returns size of vector (indices run 1:N)
             * @return 
             */
            int getSize(){return N;}
            /**
             * Return ivector of model indices corresponding to integer indices 1:N.
             * 
             * @return - ivector iv(1:N) such that iv(i) is the model index value 
             *           corresponding to the ith element.
             */
            ivector getFwdIndices(){return ptrIB->getFwdIndexVector();}
            /**
             * Return ivector of integer indices 1:N corresponding to model indices.
             * 
             * @return - ivector iv(modMin:modMax) such that iv(i) is the index into
             *           the block corresponding to model index value i. iv(i) is 0
             *           for model index values i that do not correspond to a block
             *           element.
             */
            ivector getRevIndices(){return ptrIB->getRevIndexVector();}           
            /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale value.
             * 
             * @param x - parameter-scale value as a double
             * @return - arithmetic-scale value as a double
             */
            virtual double calcArithScaleVal(double x){return BoundedNumberInfo::calcArithScaleVal(x);}
            /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
             * 
             * @param x - parameter-scale values as a dvector
             * @return - arithmetic-scale values as a dvector
             */
            virtual dvector calcArithScaleVal(const dvector& x);
            /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
             * 
             * @param x - parameter-scale values as dvar_vector
             * @return - arithmetic-scale values as dvar_vector
             */
            virtual dvar_vector calcArithScaleVal(const dvar_vector& x);
            /**
             * Calculates parameter-scale values corresponding to the input arithmetic-scale values.
             * 
             * @param x - the arithmetic-scale values as dvector
             * @return - the parameter-scale values as dvector
             */
            virtual dvector calcParamScaleVal(dvector& x);
            /**
             * Calculates the (ln-scale) prior probabilities (not NLLS!) corresponding to 
             * the values of the input vector.
             * 
             * @param x - dvar_vector to calculate ln-scale priors for
             * @return - dvar_vector of ln-scale priors (not NLLs)
             */
            dvar_vector calcLogPrior(dvar_vector& x);            
            /**
             * Gets a copy of the vector of initial values.
             * 
             * @return - the initial values.
             */
            virtual dvector getInitVals(){return initVals;}           
            /**
             * Sets initVals to the input value.
             * 
             * @param x - value to set as initial values
             */
            virtual void setInitVal(double x){initVals = x;}            
            /**
             * Sets initVals to the input vector.
             * 
             * @param x - vector to set as initial values
             */
            virtual void setInitVals(dvector& x);
            /**
             * Sets initVals on the arithmetic scale based on the input parameter vector.
             * 
             * @param x - 
             */
            virtual void setInitVals(param_init_bounded_vector & x){initVals=calcArithScaleVal(value(x));} 
            /**
             * Returns a copy of the vector of final values.
             * @return dvector of final values
             */
            virtual dvector getFinalVals(){return finlVals;}
            /**
             * Sets final values for output to R.
             * @param x
             */
            virtual void setFinalVals(param_init_bounded_vector & x){finlVals.allocate(x.indexmin(),x.indexmax()); finlVals=calcArithScaleVal(value(x));} 
            /**
             * Reads the initial values from an input stream.
             * 
             * @param [in] is - the input stream
             */
            virtual void readInitVals(cifstream & is);
            /**
             * Draws the initial values on the arithmetic scale from the 
             * specified prior distribution or by jittering.
             * 
             * If the parameter estimation phase is \<0, the initial values (initVals)
             * are returned instead.
             * 
             * @param [in] rng
             * @param [in] vif
             * 
             * @return a dvector of random values (or the initial values, if the parameter estimation phase is \<0)
             */
            virtual dvector drawInitVals(random_number_generator& rng, double vif);//draw initial values by resampling prior
            /**
             * Reads the parameter info from an input filestream.
             * The read order is:
             * <ul>
             * <li> idxType - the index type (as an adstring)
             * <li> the index block defining the vector indices (which determines N), as an IndexBlock
             * <li> readVals - an adstring flag to subsequently read a vector of initial values
             * <li> lower bound - the lower bound on the arithmetic scale
             * <li> upper bound - the upper bound on the arithmetic scale
             * <li> jitter - adstring flag to "jitter" initial values
             *  <li> initVal
             *  <li> phase
             *  <li> resample
             *  <li> priorWgt
             *  <li> priorType
             *  <li> priorParams
             *  <li> priorConsts
             *  <li> label
             * <\ul>
             * 
             * @param [in] is - the filestream to read from
             */
            virtual void read(cifstream & is);
            /**
             * Writes the parameter info to an output stream, in ADMB format.
             * 
             * @param [in] os - the output stream to write to
             */
            virtual void write(std::ostream & os);
            /**
             * Writes the parameter info to an output stream as an R list.
             * 
             * @param [in] os - the output stream to write to
             */
            virtual void writeToR(std::ostream & os);
            /**
             * Writes the final values to an output stream as an R structure.
             * 
             * @param [in] os - the output stream to write to
             */
            virtual void writeFinalValsToR(std::ostream & os);
    };

/**
 * DevsVectorInfo : BoundedVectorInfo
 * 
  *  Encapsulates characteristics for a param_init_bounded_vector
  *  to be used as the basis for a devs vector, which has the
  *  property that the sum of the elements of the vector should be 0.
  *  This should be implemented by enforcing sum(v[])=0 on the arithmetic
 *   scale in the likelihood.
  *  Otherwise, DevsVectorInfo is identical to BoundedVectorInfo.
 * 
 */
    class DevsVectorInfo : public BoundedVectorInfo {
        public:
            static int debug;
        public:
            /**
             * Class constructor.
             */
            DevsVectorInfo():BoundedVectorInfo(){}
            /**
             * Class constructor.
             * 
             * @param [in] name - adstring for name of associated param_init_bo
             */
            DevsVectorInfo(adstring& name):BoundedVectorInfo(name){}
            /**
             * Class constructor.
             * 
             * @param [in] name - adstring for name of associated param_init_bo
             */
            DevsVectorInfo(const char * name):BoundedVectorInfo(name){}
            /**
             * Class destructor.
             */
            ~DevsVectorInfo(){}
                      
            /**
             * Sets initial values to 0, no matter what @param x is.
             * @param x
             */
            void setInitVal(double x){initVals = 0.0;}
            /**
             * Draws a vector of random values by resampling the prior or jittering.
             * 
             * If the parameter estimation phase is \<0, the vector of initVals
             * is returned. The returned vector is guaranteed to sum to 0.
             * 
             * @param [in] rng - the random number generator
             * @param [in] vif - the variance inflation factor
             * 
             * @return - dvector
             */
            virtual dvector drawInitVals(random_number_generator& rng, double vif);//draw initial values by resampling prior
    };

/**
 * NumberVectorInfo
 * 
 * Encapsulates parameter info for a param_init_number_vector as a 1-d array
 * of NumberInfo objects.
 * 
 * Public fields:
 * <ul>
 * <li> debug - static int flag to print debugging info
 * <li> name - adstring with name of associated param_init_number_vector object
 * <\ul>
 * 
 * Protected fields:
 * <ul>
 * <li> nNIs - number of elements (parameters) in the associated param_init_number_vector
 * <li> ppNIs - pointer to vector of pointers to NumberInfo instances
 * <\ul>
 */
    class NumberVectorInfo {
        public:
            /* int flag to print debugging info */
            static int debug;
        public:
            /* name of associated param_init_number_vector */
            adstring name;
        protected:
            /* number of elements (parameters) in param_init_number_vector */
            int nNIs; 
            /* ptr to vector of ptrs to NumberInfo instances */
            NumberInfo** ppNIs;
        public:
            /**
             * Class constructor.
             * 
             * Sets name to "".
             */
            NumberVectorInfo(){this->name="";nNIs=0;ppNIs=0;}
            /**
             * Class constructor.
             * 
             * @param name - adstring with name of associated param_init_number_vector
             */
            NumberVectorInfo(adstring& name){this->name=name;nNIs=0;ppNIs=0;}
            /**
             * Class constructor.
             * 
             * @param name - adstring with name of associated param_init_number_vector
             */
            NumberVectorInfo(const char * name){this->name=name;nNIs=0;ppNIs=0;}
            /**
             * Class destructor.
             */
            ~NumberVectorInfo(){deallocate();}
            
            /**
             * Operator to extract a pointer to a NumberInfo object.
             * 
             * @param i - the 1-based index of the parameter which to extract the NumberInfo for
             * 
             * @return - the pointer to the requested NumberInfo object
             */
            NumberInfo* operator[](int i){if (ppNIs&&(i<=nNIs)) return ppNIs[i-1]; return 0;}
            /**
             * Deallocates the pointer ppNIs and the associated pointers in the array.
             */
            void deallocate();
            /**
             * Gets the number of elements (parameters) in the associated param_init_number_vector.
             * 
             * @return the size of the associated param_init_number_vector 
             */
            int getSize(void){return nNIs;}
            /**
             * Gets an ivector of the estimation phase for each parameter in the
             * associated param_init_number_vector.
             * 
             * @return - an ivector
             */
            ivector getPhases(void);
            /**
             * Gets a dvector of the likelihood weight for each parameter in the 
             * associated param_init_number_vector.
             * 
             * @return - a dvector
             */
            dvector getPriorWgts(void);
            /**
             * Gets a dvector of the initial values for each parameter in the 
             * associated param_init_number_vector.
             * 
             * @return - a dvector
             */
            dvector getInitVals(void);
            /**
             * Calculates a dvar_vector of the log-scale prior probability for
             * each element in the input vector.
             * The input vector should have the same size and indexing as the 
             * associated param_init_number_vector.
             * 
             * @param [in] pv - a dvar_vector of values to compute log priors for
             * 
             * @return a dvar_vector of the log priors
             */
            dvar_vector calcLogPriors(dvar_vector & pv);    
            
            /** 
             * Gets an adstring_array of labels corresponding to the 
             * individual parameters in the associated param_init_number_vector.
             * 
             * @return an adstring_array
             */
            adstring_array getLabels(void);
            
            /** 
             * Gets an adstring_array of parameer scale types corresponding to the 
             * individual parameters in the associated param_init_number_vector.
             * 
             * @return an adstring_array
             */
            adstring_array getParamScales(void);
            
            /**
             * Draws a dvector of random values to use as
             * initial values for each parameter in the associated param_init_number_vector, 
             * based on the prior defined for each parameter.
             * If the estimation phase for a given parameter is \<0, then the initial value
             * for that parameter will be returned rather than a random value.
             * 
             * @param [in] rng - a random number generator
             * @param [in] vif - a variance inflation factor
             * 
             * @return a dvar_vector of same size and indices as the associated param_init_number_vector
             */
            virtual dvector drawInitVals(random_number_generator& rng, double vif);
            virtual void setInitVals(param_init_number_vector& x);
            virtual void setFinalVals(param_init_number_vector& x);
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
            virtual void writeFinalValsToR(std::ostream& os);
            friend cifstream& operator >>(cifstream & is, NumberVectorInfo & obj){obj.read(is);return is;}
            friend std::ostream& operator <<(std::ostream & os, NumberVectorInfo & obj){obj.write(os);return os;}
    };

/**
 * Encapsulates info for a param_init_bounded_number_vector as a 1-d array of BoundedNumberInfo's.
 * 
 * Subclass of NumberVectorInfo
 */
    class BoundedNumberVectorInfo: public NumberVectorInfo {
        public:
            static int debug;
        public:
            /**
             * Class constructor.
             */
            BoundedNumberVectorInfo():NumberVectorInfo(){}
            /**
             * Class constructor.
             * 
             * @param name - name of associated parameter vector as an adstring
             */
            BoundedNumberVectorInfo(adstring& name):NumberVectorInfo(name){}
            /**
             * Class constructor.
             * 
             * @param name - name of associated parameter vector as a const char*
             */
            BoundedNumberVectorInfo(const char * name):NumberVectorInfo(name){}
            /**
             * Class destructor.
             */
            ~BoundedNumberVectorInfo(){deallocate();}
            
            /**
             * Gets a pointer to the ith BoundedNumberInfo element in the vector. 
             * i must be in the interval 1<=i<=nNIs.
             * 
             * @param i - index (1-based) to BoundedNumberInfo
             * @return pointer to ith BoundedNumberInfo object
             */
            BoundedNumberInfo* operator[](int i);
            
            /**
             * Gets the parameter scale types for the associated aprameters.
             * 
             * @return - an adstring_array of scale types
             */
            adstring_array getScaleTypes();
            
            /**
             * Get the lower bounds for parameters as vector.
             * 
             * @return - a dvector with the lower bound for each parameter
             */
            dvector getLowerBounds();
            /**
             * Get the upper bounds for parameters as vector.
             * 
             * @return - a dvector with the upper bound for each parameter
            */
            dvector getUpperBounds();
              /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
             * 
             * @param x - parameter-scale values as a dvector
             * @return - arithmetic-scale values as a dvector
             */
            dvector calcArithScaleVals(const dvector& x);
            /**
             * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
             * 
             * @param x - parameter-scale values as dvar_vector
             * @return - arithmetic-scale values as dvar_vector
             */
            dvar_vector calcArithScaleVals(const dvar_vector& x);
            /**
             * Calculates parameter-scale values corresponding to the input arithmetic-scale values.
             * 
             * @param x - the arithmetic-scale values as dvector
             * @return - the parameter-scale values as dvector
             */
            dvector calcParamScaleVals(dvector& x);
          
            virtual dvector drawInitVals(random_number_generator& rng, double vif);
            virtual void setInitVals(param_init_bounded_number_vector& x);
            virtual void setFinalVals(param_init_bounded_number_vector& x);
            void read(cifstream & is);
            void write(std::ostream & os);
            void writeToR(std::ostream& os, adstring nm, int indent=0);
    };
    
//--------------------------------------------------------------------------------
//          VectorVectorInfo
//  Encapsulates info for a param_init_vector_vector as a 1-d array of VectorInfo's
//--------------------------------------------------------------------------------
    class VectorVectorInfo {
        public:
            static int debug;
        public:
            adstring name;     //name of vector of parameter vectors
        protected:
            int nVIs;          //number of parameter vectors in vector
            VectorInfo** ppVIs;//ptr to vector of ptrs to VectorInfo instances
        public:
            VectorVectorInfo(){this->name="";nVIs=0;ppVIs=0;}
            VectorVectorInfo(adstring& name){this->name=name;nVIs=0;ppVIs=0;}
            VectorVectorInfo(const char * name){this->name=name;nVIs=0;ppVIs=0;}
            ~VectorVectorInfo(){deallocate();}
            
            VectorInfo* operator[](int i){if ((ppVIs>0)&&(i<=nVIs)) return ppVIs[i-1]; return 0;}
            void deallocate();
            int getSize(void){return nVIs;}
            ivector getMinIndices(void);
            ivector getMaxIndices(void);
            ivector getPhases(void);
            dvector getPriorWgts(void);
            
             /* Function to get array of labels corresponding to parameters */
            adstring_array getLabels(void);
            
           virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
            virtual void writeFinalValsToR(std::ostream& os);
            friend cifstream& operator >>(cifstream & is, VectorVectorInfo & obj){obj.read(is);return is;}
            friend std::ostream& operator <<(std::ostream & os, VectorVectorInfo & obj){obj.write(os);return os;}
    };

//--------------------------------------------------------------------------------
//          BoundedVectorVectorInfo
//  Encapsulates info for a param_init_bounded_vector_vector as a 1-d array of BoundedVectorInfo's
//--------------------------------------------------------------------------------
    class BoundedVectorVectorInfo {
        public:
            static int debug;
        public:
            adstring name;     //name of vector of parameter vectors
        protected:
            int nVIs;                 //number of parameter vectors in vector
            BoundedVectorInfo** ppVIs;//ptr to vector of ptrs to BoundedVectorInfo instances
        public:
            BoundedVectorVectorInfo(){this->name="";nVIs=0;ppVIs=0;}
            BoundedVectorVectorInfo(adstring& name){this->name=name;nVIs=0;ppVIs=0;}
            BoundedVectorVectorInfo(const char * name){this->name=name;nVIs=0;ppVIs=0;}
            ~BoundedVectorVectorInfo(){deallocate();}
            
            BoundedVectorInfo* operator[](int i){if ((ppVIs>0)&&(i<=nVIs)) return ppVIs[i-1]; return 0;}
            void deallocate();
            int getSize(void){return nVIs;}
            ivector getMinIndices(void);
            ivector getMaxIndices(void);
            ivector getPhases(void);
            dvector getPriorWgts(void);
            dvector getLowerBounds(void);
            dvector getUpperBounds(void);
            
            /* Function to get array of labels corresponding to parameters */
            adstring_array getLabels(void);
            
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
            virtual void writeFinalValsToR(std::ostream& os);
            friend cifstream& operator >>(cifstream & is, BoundedVectorVectorInfo & obj){obj.read(is);return is;}
            friend std::ostream& operator <<(std::ostream & os, BoundedVectorVectorInfo & obj){obj.write(os);return os;}
    };

/**
 * DevsVectorVectorInfo
 * 
  *  Encapsulates characteristics for a param_init_bounded_vector_vector
  *  to be used as the basis for a vector of devs vectors. Each devs vector has the
  *  property that the sum of the elements of the vector should be 0.
  *  This should be implemented by enforcing sum(v[])=0 on the arithmetic
 *   scale in the likelihood, for each devs vector v.
  *  Otherwise, DevsVectorVectorInfo is identical to BoundedVectorVectorInfo,
 *   of which it is a subclass.
 * 
 */
    class DevsVectorVectorInfo: public BoundedVectorVectorInfo {
        public:
            static int debug;
        public:
            DevsVectorVectorInfo():BoundedVectorVectorInfo(){}
            DevsVectorVectorInfo(adstring& name):BoundedVectorVectorInfo(name){}
            DevsVectorVectorInfo(const char * name):BoundedVectorVectorInfo(name){}
            
            DevsVectorInfo* operator[](int i){if ((ppVIs>0)&&(i<=nVIs)) return (DevsVectorInfo*) ppVIs[i-1]; return 0;}
            
            virtual void read(cifstream & is);
  //          virtual void write(std::ostream & os);
  //          virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
  //          virtual void writeFinalValsToR(std::ostream& os);
            friend cifstream& operator >>(cifstream & is, DevsVectorVectorInfo & obj){obj.read(is);return is;}
            friend std::ostream& operator <<(std::ostream & os, DevsVectorVectorInfo & obj){obj.write(os);return os;}
    };

#endif //MODELPARAMETERINFOTYPES_HPP
