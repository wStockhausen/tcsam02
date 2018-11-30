/**
 * @file
 * 
 * Header definitions for model parameter info type classes.
 * 
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
 *  This class encapsulates characteristics for a param_init_number parameter on a possibly 
 *  transformed scale (log). Values (initial, final)
 *  are on the "arithmetic" scale, which may be different from the parameter scale. Prior
 *  distributions are also described on the arithmetic scale. The "scaleType"
 *  indicates the transformation from the arithmetic to parameter scale.
 * 
 * Class fields:
 * <ul>
 * <li> scaleType   - parameter scale type
 * <li> initVal     - initial value, on the arithmetic scale, for the associated parameter
 * <li> finlVal     - final value, on the arithmetic scale, for associated parameter
 * <li> phase       - phase to start estimation of associated parameter
 * <li> resample    - flag to resample initial value from prior
 * <li> priorWgt    - weight assigned to prior in likelihood
 * <li> priorType   - name identifying pdf for prior
 * <li> priorParams - dvar_vector of parameters (on the arithmetic scale) for the prior
 * <li> priorConsts - dvector of constants (on the arithmetic scale) for the prior
 * <li> pMPI        - pointer to ModelPDFInfo object for resampling, prior calculations
 * <li> name        - the name of the associated parameter
 * <li> label       - a label for the associated parameter
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
        dvector priorConsts;
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
        virtual dvariable calcArithScaleVal(const prevariable& x);
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
         * Gets the value, on the arithmetic scale, to use as the initial value.
         * Use \code getInitValOnParamScale() \endcode to get the initial value
         * on the correct scale for the associated parameter.
         * 
         * @return - the value, as a double
         */
        double getInitVal(){return initVal;}
        /**
         * Gets the value, on the parameter scale, to assign as the initial value
         * to the associated parameter.
         * 
         * @return - the value, as a double
         */
        virtual double getInitValOnParamScale(){return calcParamScaleVal(initVal);}
        /**
         * Sets the value, on the arithmetic scale, used as the initial value.
         * <p>
         * The value given should be on the arithmetic scale.
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
        virtual void setInitValFromParamVal(const prevariable& x){initVal=calcArithScaleVal(value(x));}
        /**
         * Draws a random value on the arithmetic scale based on the 
         * specified prior probability distribution.
         * 
         * Note that this DOES NOT update initVal.
         * If the estimation phase is \< 0, the value of initVal is returned
         * 
         * @param rng - the random number generator
         * @param vif - the variance inflation factor
         * 
         * @return - the random number on the arithmetic scale, or the value of initVal
         */
        virtual double drawInitVal(random_number_generator& rng, double vif);
        /**
         * Calculates the log-scale prior probability based on an input value.
         * 
         * The input value is presumably from the associated parameter value, 
         * but should be on the "arithmetic" scale.
         * 
         * @param x - the input value on the arithmetic scale, presumably based on the associated parameter
         * 
         * @return - the log-scale prior, as a dvariable
         */
        virtual dvariable calcLogPrior(prevariable& x);             
         /**
         * Gets the "final value", on the arithmetic scale. This value needs to be set
         * using \code setFinalValFromParamVal(...) \endcode
         * 
         * @return x - the final value
         */
        virtual double getFinalVal(){return finlVal;}
        /**
         * Sets the "final value", on the arithmetic scale, from the value for the
         * associated parameter.
         * 
         * @param x - the associated parameter (an instance of a param_init_number)
         */
        virtual void setFinalValFromParamVal(const prevariable& x){finlVal=calcArithScaleVal(value(x));}
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
};//NumberInfo

/**
* BoundedNumberInfo : NumberInfo
* 
*  Encapsulates characteristics for a param_init_bounded_number parameter on a 
*  possibly-transformed (logit, lognormal, probit) scale. Values (initial, final)
*  are on the "arithmetic" scale, which may be different from the parameter scale. Prior
*  distributions are also described on the "arithmetic" scale. The "scaleType"
*  indicates the transformation from "arithmetic" to parameter scale.
* 
* Inherited class fields:
* <ul>
* <li> scaleType   - parameter scale type
* <li> initVal     - initial value, on the arithmetic scale, for the associated parameter
* <li> finlVal     - final value, on the arithmetic scale, for associated parameter
* <li> phase       - phase to start estimation of associated parameter
* <li> resample    - flag to resample initial value from prior
* <li> priorWgt    - weight assigned to prior in likelihood
* <li> priorType   - name identifying pdf for prior
* <li> priorParams - dvar_vector of parameters (on the arithmetic scale) for the prior
* <li> priorConsts - dvector of constants (on the arithmetic scale) for the prior
* <li> name        - the name of the associated parameter
* <li> label       - a label for the associated parameter
* </ul>
* 
* New class fields:
* <ul>
* <li> jitter    - integer flag indicating whether or not to jitter initial values
* <li> lower     - lower bound, on the "arithmetic" scale
* <li> upper     - upper bound, on the "arithmetic" scale
* </ul>
* 
* Note that the priorParams and priorConsts for this class are interpreted on the
* arithmetic scale, not on the (possibly transformed) scale of the associated parameter.
*/
class BoundedNumberInfo : public NumberInfo {
    public:
        static int debug;
    protected:
        /** the lower bound, on the arithmetic scale */
        double lower;
        /** the upper bound, on the arithmetic scale */
        double upper;
    public:
        /** flag indicating whether or not to jitter initial value */
        int jitter;
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
         * Returns the arithmetic-scale lower bound.
         * 
         * @return 
         */            
        double getLowerBound(){return lower;}
        /**
         * Returns the parameter-scale lower bound.
         * 
         * @return 
         */            
        double getLowerBoundOnParamScale(){return this->calcParamScaleVal(lower);}
        /**
         * Returns the arithmetic-scale upper bound.
         * 
         * @return 
         */            
        double getUpperBound(){return upper;}           
        /**
         * Returns the parameter-scale upper bound.
         * 
         * @return 
         */            
        double getUpperBoundOnParamScale(){return this->calcParamScaleVal(upper);}           
        /**
         * Calculates the arithmetic-scale value corresponding to a given parameter-scale value.
         * 
         * @param x - the parameter-scale value
         * @return - the corresponding arithmetic-scale value
         * 
         * @overrides NumberInfo::calcArithScaleVal(double x)
         */
        virtual double calcArithScaleVal(double x);
        /**
         * Calculates the arithmetic-scale value corresponding to given parameter-scale value.
         * 
         * @param x - the parameter-scale value
         * @return - the corresponding arithmetic-scale value, as a dvariable
         * 
         * @overrides NumberInfo::calcArithScaleVal(const prevariable& x)
         */
        virtual dvariable calcArithScaleVal(const prevariable& x);
        /**
         * Calculates the parameter-scale value corresponding to an arithmetic-scale value.
         * 
         * @param x - the arithmetic-scale value
         * @return - the parameter-scale value
         * 
         * @overrides NumberInfo::calcParamScaleVal(double x)
         */
        virtual double calcParamScaleVal(double x);
        /**
         * Gets the value, on the parameter scale, to assign as the initial value
         * to the associated parameter.
         * 
         * @return - the value, as a double
         * 
         * @overrides NumberInfo::getInitValOnParamScale()
        */
        virtual double getInitValOnParamScale(){return calcParamScaleVal(initVal);}
        /**
         * Sets the initial value to be used for the associated parameter, assuming
         * the value given is on the arithmetic scale.
         * The method ensures the initial value is set within the prescribed bounds.
         * 
         * @param x - the initial value, on the arithmetic scale, to set
         * 
         * @overrides NumberInfo::setInitVal(double x)
         */
        virtual void setInitVal(double x);           
        /**
         * Sets the initial value on the arithmetic scale, based on a 
         * bounded parameter on a (possibly) non-arithmetic scale.
         * 
         * @param x - the associated parameter (a param_init_bounded_number)
         * 
         * @overrides NumberInfo::setInitValFromParamVal(prevariable& x)
         */
        virtual void setInitValFromParamVal(const prevariable& x){initVal=calcArithScaleVal(value(x));}
        /**
         * Set the "final value" on the arithmetic scale, based on a 
         * bounded parameter on a (possibly) non-arithmetic scale.
         * 
         * @param x - the associated bounded parameter (an instance of a param_init_bounded_number)
         * 
         * @overrides NumberInfo::setFinalValFromParamVal(prevariable& x)
         */
        virtual void setFinalValFromParamVal(const prevariable& x){finlVal=calcArithScaleVal(value(x));}            
        /**
         * Draws a jittered random value on the arithmetic scale based on the 
         * bounds set and the fraction of the range to jitter (jitFrac).
         * 
         * Note that this DOES NOT update initVal.
         * If the estimation phase is \< 0, the value of initVal is returned
         * 
         * @param rng - the random number generator
         * @param jitFrac - the fraction of the range across which to jitter
         * 
         * @return - the random number on the arithmetic scale, or the value of initVal
         */
        virtual double jitterInitVal(random_number_generator& rng, double jitFrac);
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
         * @overrides NumberInfo::read(cifstream & is)
         */
        virtual void read(cifstream & is);
        /**
         * Writes the parameter info in ADMB format to an output stream.
         * 
         * @param os - the output stream
         * 
         * @overrides NumberInfo::write(std::ostream & os)
         */
        virtual void write(std::ostream & os);
        /**
         * Writes the parameter info as part of an R list to an output stream.
         * 
         * @param os - the output stream.
         * 
         * @overrides NumberInfo::writeToR(std::ostream& os)
         */
        virtual void writeToR(std::ostream& os);
    protected:
        /**
         * Writes the BoundedNumberInfo part of the parameter info as part of an R list to an output stream.
         * 
         * @param os - the output stream.
         * 
         * @overrides NumberInfo::writeToR1(std::ostream & os)
         */
        virtual void writeToR1(std::ostream & os);
};//BoundedNumberInfo

/**
* VectorInfo
* 
*  Encapsulates characteristics for a param_init_vector parameter on a possibly 
*  transformed scale (log). Values (initial, final)
*  are on the "arithmetic" scale, which may be different from the parameter scale. Prior
*  distributions are also described on the "arithmetic" scale. The "scaleType"
*  indicates the transformation from "arithmetic" to parameter scale.
* 
* Class fields:
* <ul>
* <li> N        -  size of parameter vector
* <li> idxType  - adstring with type of indices (e.q., YEAR)
* <li> ptrIB    - pointer to IndexBlock object defining indices
* <li> initVals - vector of initial values, on the arithmetic scale
* <li> finlVals - vector of final values, on th arithmetic scale
* <li> scaleType   - parameter scale type
* <li> phase       - phase to start estimation of associated parameter
* <li> resample    - flag to resample initial value from prior
* <li> priorWgt    - weight assigned to prior in likelihood
* <li> priorType   - name identifying pdf for prior
* <li> priorParams - dvar_vector of parameters (on the arithmetic scale) for the prior
* <li> priorConsts - dvector of constants (on the arithmetic scale) for the prior
* <li> readVals    - flag to read vector of values from input stream
* <li> name        - the name of the associated parameter
* <li> label       - a label for the associated parameter
* </ul>
* 
* Note that the priorParams and priorConsts for this class are interpreted on the
* arithmetic scale.
*/
class VectorInfo {
    public:
        /* flag to turn on debugging output */
        static int debug;
        /* size of initVals vector*/
        int N;
        /* the parameter name */
        adstring name;
        /* the label for the associated parameter */
        adstring label;
        /* flag to do resampling of initial values */
        bool resample;
        /* pointer to info for resmapling pdf */
        ModelPDFInfo*  pMPI;
        /* flag to read initial values */
        int readVals;
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
        dvector priorConsts;
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
        VectorInfo(){this->name="";pMPI=0;scaleType=-1;ptrIB=0;}
        /**
         * Class constructor.
         * 
         * @param [in] name - name of associated param_init_vector
         */
        VectorInfo(adstring& name){this->name=name;pMPI=0;scaleType=-1;ptrIB=0;}
        /**
         * Class constructor.
         * 
         * @param [in] name - name of associated param_init_vector
         */
        VectorInfo(const char* name){this->name=name;pMPI=0;scaleType=-1;ptrIB=0;}
        /**
         * Class destructor.
         */
        ~VectorInfo(){if(ptrIB){delete ptrIB;} if(pMPI){delete pMPI;}}

        /**
         * Gets the size of the vector (indices run 1:N).
         * 
         * @return 
         */
        int getSize(){return N;}
        
        /**
         * Add a value to the end of the init and final vectors.
         * 
         * @param val - value to add
         * @param ibVal - corresponding IndexBlock value to add
         */
        void addValueOnArithmeticScale(double val, int ibVal);
        
        /**
         * Add a value to the end of the init and final vectors.
         * 
         * @param val - value to add
         * @param ibVal - corresponding IndexBlock value to add
         */
        void addValueOnParameterScale(double val, int ibVal);
        
        /**
         * Gets the parameter scale type, as an adstring.
         * 
         * @return the scale type
         */
        adstring getScaleType(){return tcsam::getScaleType(scaleType);}

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
         * Sets the prior type, based on an adstring value.
         * 
         * @param prior - the prior type, as an adstring
         */
        virtual void setPriorType(adstring & prior);
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
         * Gets the vector of initial values on the arithmetic scale. 
         * Vector indices run 1:N.
         * 
         * @return - dvector of initial values
         */
        dvector getInitVals(){return initVals;}
        /**
         * Gets the vector of initial values on the parameter scale. 
         * Vector indices run 1:N.
         * 
         * @return - dvector of initial values
         */
        virtual dvector getInitValsOnParamScale(){return calcParamScaleVals(initVals);}
        /**
         * Gets a copy of the vector of final values on the arithmetic scale.
         * 
         * @return a dvector of final values
         */
        dvector getFinalVals(){return finlVals;}
        /**
         * Calculates arithmetic-scale values corresponding to the input parameter-scale value.
         * 
         * @param x - parameter-scale value as a double
         * @return - arithmetic-scale value as a double
         */
        virtual double calcArithScaleVals(double x);
        /**
         * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
         * 
         * @param x - parameter-scale values as a dvector
         * @return - arithmetic-scale values as a dvector
         */
        virtual dvector calcArithScaleVals(const dvector& x);
        /**
         * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
         * 
         * @param x - parameter-scale values as dvar_vector
         * @return - arithmetic-scale values as dvar_vector
         */
        virtual dvar_vector calcArithScaleVals(const dvar_vector& x);
        /**
         * Calculates parameter-scale value corresponding to the input arithmetic-scale value.
         * 
         * @param x - the arithmetic-scale value as double
         * @return - the parameter-scale value as double
         * 
         * @overrides VectorInfo::calcParamScaleVals(double)
         */
        virtual double calcParamScaleVals(double x);
        /**
         * Calculates parameter-scale values corresponding to the input arithmetic-scale values.
         * 
         * @param x - the arithmetic-scale values as dvector
         * @return - the parameter-scale values as dvector
         */
        virtual dvector calcParamScaleVals(dvector& x);
        /**
         * Sets initial values, on the arithmetic scale, to the input constant.
         * 
         * @param [in] x - the value to set
         */
        virtual void setInitVals(double x){initVals = x;}
        /**
         * Sets the vector of initial values to the input dvector element-by-element.
         * 
         * @param [in] x - a dvector of initial values to set. Indices should run 1:N
         */
        virtual void setInitVals(dvector& x);     
        /**
         * Sets the vector of initial values, on the arithmetic scale, to the 
         * possibly inverse-transformed values of the input vector, 
         * which should be on the parameter scale.
         * 
         * @param [in] x - a dvar_vector on the parameter scale. Indices should run 1:N
         */
        virtual void setInitValsFromParamVals(const dvar_vector& x){initVals=calcArithScaleVals(value(x));}     
        /**
         * Sets the vector of "final values" to the (possibly inverse-transformed to the arithmetic scale)
         * values of the input vector, which should be on the parameter scale.
         * 
         * @param [in] x - a dvar_vector on the parameter scale.
         */
        virtual void setFinalValsFromParamVals(const dvar_vector& x){finlVals=calcArithScaleVals(value(x));}     
        /**
         * Calculates a vector of the log-scale prior probability based on an input vector of values.
         * 
         * The input values are presumably from the associated parameter values, 
         * but should be on the "arithmetic" scale.
         * 
         * @param x - the input dvar_vector on the arithmetic scale, presumably based on the associated parameters
         * 
         * @return - the log-scale prior, as a dvar_vector
         */
        virtual dvar_vector calcLogPrior(dvar_vector& x);
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
         *  <li> idxType - the index type (as an adstring)
         *  <li> the index block defining the vector indices (which determines N), as an IndexBlock
         *  <li> readVals - an adstring flag to subsequently read a vector of initial values
         *  <li> initVal
         *  <li> scaleType
         *  <li> phase
         *  <li> resample
         *  <li> priorWgt
         *  <li> priorType
         *  <li> priorParams
         *  <li> priorConsts
         *  <li> label
         * </ul>
         * 
         * @param [in] is - the filestream to read from
         */
        virtual void read(cifstream & is);
        /**
         * Reads part of the parameter info from an input filestream.
         * The read order is:
         * <ul>
         *  <li> idxType - the index type (as an adstring)
         *  <li> the index block defining the vector indices (which determines N), as an IndexBlock
         *  <li> readVals - an adstring flag to subsequently read a vector of initial values
         * </ul>
         * 
         * @param [in] is - the filestream to read from
         */
        void readPart1(cifstream& is);
        /**
         * Reads part of the parameter info from an input filestream.
         * The read order is:
         * <ul>
         *  <li> initVal
         *  <li> scaleType
         *  <li> phase
         *  <li> resample
         *  <li> priorWgt
         *  <li> priorType
         *  <li> priorParams
         *  <li> priorConsts
         *  <li> label
         * </ul>
         * 
         * @param [in] is - the filestream to read from
         */
        void readPart2(cifstream& is);
        /**
         * Writes the parameter info to an output stream, in ADMB format.
         * 
         * @param [in] os - the output stream to write to
         */
        virtual void write(std::ostream & os);
        void writePart1(std::ostream & os);
        void writePart2(std::ostream & os);
        /**
         * Writes the parameter info to an output stream as an R list.
         * 
         * @param [in] os - the output stream to write to
         */
        virtual void writeToR(std::ostream & os);
        void writeToR1(std::ostream & os);
        /**
         * Writes the final values to an output stream as an R structure.
         * 
         * @param [in] os - the output stream to write to
         */
        virtual void writeFinalValsToR(std::ostream & os);
        friend cifstream& operator >>(cifstream & is, VectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, VectorInfo & obj){obj.write(os);return os;}
};//VectorInfo

/**
* BoundedVectorInfo : VectorInfo
* 
*  Encapsulates characteristics for a param_init_bounded_vector 
*  parameter vector on a possibly-transformed (logit, lognormal, probit) scale.
*  Basically, a bounded parameter version of VectorInfo.
* 
* Inherited class fields:
* <ul>
* <li> N        -  size of parameter vector
* <li> idxType  - adstring with type of indices (e.q., YEAR)
* <li> ptrIB    - pointer to IndexBlock object defining indices
* <li> initVals - vector of initial values, on the arithmetic scale
* <li> finlVals - vector of final values, on th arithmetic scale
* <li> scaleType   - parameter scale type
* <li> phase       - phase to start estimation of associated parameter
* <li> resample    - flag to resample initial value from prior
* <li> priorWgt    - weight assigned to prior in likelihood
* <li> priorType   - name identifying pdf for prior
* <li> priorParams - dvar_vector of parameters (on the arithmetic scale) for the prior
* <li> priorConsts - dvector of constants (on the arithmetic scale) for the prior
* <li> readVals    - flag to read vector of values from input stream
* <li> name        - the name of the associated parameter
* <li> label       - a label for the associated parameter
* </ul>
* 
* New class fields:
* <ul>
* <li> jitter - integer flag indicating whether or not to jitter initial values
* <li> lower  - lower bound, on the "arithmetic" scale
* <li> upper  - upper bound, on the "arithmetic" scale
* </ul>
*/
class BoundedVectorInfo : public VectorInfo {
    public:
        static int debug;
    protected:
        /** the lower bound, on the arithmetic scale */
        double lower;
        /** the upper bound, on the arithmetic scale */
        double upper;
    public:
        /** flag indicating whether or not to jitter initial value */
        int jitter;
    public:
        /**
         * Class constructor.
         */
        BoundedVectorInfo():VectorInfo(){}
        /**
         * Class constructor.
         * 
         * @param [in] name - adstring for name of associated param_init_bo
         */
        BoundedVectorInfo(adstring& name):VectorInfo(name){}
        /**
         * Class constructor.
         * 
         * @param [in] name - adstring for name of associated param_init_bo
         */
        BoundedVectorInfo(const char * name):VectorInfo(name){}
        /**
         * Returns the arithmetic-scale lower bound.
         * 
         * @return 
         */            
        double getLowerBound(){return lower;}
        /**
         * Returns the parameter-scale lower bound.
         * 
         * @return 
         */            
        double getLowerBoundOnParamScale(){return this->calcParamScaleVals(lower);}
        /**
         * Returns the arithmetic-scale upper bound.
         * 
         * @return 
         */            
        double getUpperBound(){return upper;}           
        /**
         * Returns the parameter-scale upper bound.
         * 
         * @return 
         */            
        double getUpperBoundOnParamScale(){return this->calcParamScaleVals(upper);}           
        /**
         * Gets the vector of initial values on the parameter scale. 
         * Vector indices run 1:N.
         * 
         * @return - dvector of initial values
         */
        virtual dvector getInitValsOnParamScale(){return calcParamScaleVals(initVals);}
        /**
         * Calculates arithmetic-scale values corresponding to the input parameter-scale value.
         * 
         * @param x - parameter-scale value as a double
         * @return - arithmetic-scale value as a double
         * 
         * @overrides VectorInfo::calcArithScaleVals(double)
         */
        virtual double calcArithScaleVals(double x);
        /**
         * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
         * 
         * @param x - parameter-scale values as a dvector
         * @return - arithmetic-scale values as a dvector
         * 
         * @overrides VectorInfo::calcArithScaleVals(const dvector&)
         */
        virtual dvector calcArithScaleVals(const dvector& x);
        /**
         * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
         * 
         * @param x - parameter-scale values as dvar_vector
         * @return - arithmetic-scale values as dvar_vector
         * 
         * @overrides VectorInfo::calcArithScaleVals(const dvar_vector&)
         */
        virtual dvar_vector calcArithScaleVals(const dvar_vector& x);
        /**
         * Calculates parameter-scale value corresponding to the input arithmetic-scale value.
         * 
         * @param x - the arithmetic-scale value as double
         * @return - the parameter-scale value as double
         * 
         * @overrides VectorInfo::calcParamScaleVals(double)
         */
        virtual double calcParamScaleVals(double x);
        /**
         * Calculates parameter-scale values corresponding to the input arithmetic-scale values.
         * 
         * @param x - the arithmetic-scale values as dvector
         * @return - the parameter-scale values as dvector
         * 
         * @overrides VectorInfo::calcParamScaleVals(double)
         */
        virtual dvector calcParamScaleVals(dvector& x);
        /**
         * Sets initVals to the input value.
         * 
         * @param x - value to set as initial values
         * 
         * @overrides VectorInfo::setInitVals(double)
         */
        virtual void setInitVals(double x);           
        /**
         * Sets initVals to the input vector.
         * 
         * @param x - vector to set as initial values
         * 
         * @overrides VectorInfo::setInitVals(dvector& x)
         */
        virtual void setInitVals(dvector& x);
        /**
         * Sets initVals on the arithmetic scale based on the input parameter vector.
         * 
         * @param x - the associated parameter vector
         * 
         * @overrides VectorInfo::setInitValsFromParamVals(dvar_vector& x)
         */
        virtual void setInitValsFromParamVals(const dvar_vector & x){initVals=calcArithScaleVals(value(x));} 
        /**
         * Sets final values for output to R. Input vector is assumed to be
         * on the parameter scale, so values are transformed to the arithmetic
         * scale.
         * 
         * @param x  - vector of parameter-scale values to set as final values
         * 
         * @overrides VectorInfo::setFinalValsFromParamVals(dvar_vector& x)
         */
        virtual void setFinalValsFromParamVals(const dvar_vector & x){finlVals=calcArithScaleVals(value(x));} 
        /**
         * Calculates a vector of the log-scale prior probability based on an input vector of values.
         * 
         * The input values are presumably from the associated parameter values, 
         * but should be on the "arithmetic" scale.
         * 
         * @param x - the input dvar_vector on the arithmetic scale, presumably based on the associated parameters
         * 
         * @return - the log-scale prior, as a dvar_vector
         */
        virtual dvar_vector calcLogPrior(dvar_vector& x);
        /**
         * Reads the initial values from an input stream.
         * 
         * @param [in] is - the input stream
         */
        virtual void readInitVals(cifstream & is);
        /**
         * Draws initial values on the arithmetic scale from the 
         * specified prior distribution or by jittering.
         * 
         * If the parameter estimation phase is \<0, the initial values (initVals)
         * are returned instead.
         * 
         * @param [in] rng - the random number generator object
         * @param [in] vif - the desired variance inflation factor
         * 
         * @return a dvector of random values (or the initial values, if the parameter estimation phase is \<0)
         */
        virtual dvector drawInitVals(random_number_generator& rng, double vif);//draw initial values by resampling prior
        /**
         * Draws a jittered random value on the arithmetic scale based on the 
         * bounds set and the fraction of the range to jitter (jitFrac).
         * 
         * Note that this DOES NOT update initVal.
         * If the estimation phase is \< 0, the value of initVal is returned
         * 
         * @param rng - the random number generator
         * @param jitFrac - the fraction of the range across which to jitter
         * 
         * @return - the random number on the arithmetic scale, or the value of initVal
         */
        virtual dvector jitterInitVals(random_number_generator& rng, double jitFrac);
        /**
         * Reads the parameter info from an input filestream.
         * The read order is:
         * <ul>
         *  <li> idxType - the index type (as an adstring)
         *  <li> the index block defining the vector indices (which determines N), as an IndexBlock
         *  <li> readVals - an adstring flag to subsequently read a vector of initial values
         *  <li> lower bound - the lower bound on the arithmetic scale
         *  <li> upper bound - the upper bound on the arithmetic scale
         *  <li> jitter - adstring flag to "jitter" initial values
         *  <li> scaleType
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
         * @param [in] is - the filestream to read from
         * 
         * @overrides VectorInfo::read(cifstream& is)
         */
        virtual void read(cifstream & is);
        /**
         * Writes the parameter info to an output stream, in ADMB format.
         * 
         * @param [in] os - the output stream to write to
         * 
         * @overrides VectorInfo::write(std::ostream & os)
         */
        virtual void write(std::ostream & os);
        /**
         * Writes the parameter info to an output stream as an R list.
         * 
         * @param [in] os - the output stream to write to
         * 
         * @overrides VectorInfo::writeToR(std::ostream & os)
         */
        virtual void writeToR(std::ostream & os);
        /**
         * Writes the final values to an output stream as an R structure.
         * 
         * @param [in] os - the output stream to write to
         * 
         * @overrides VectorInfo::writeFinalValsToR(std::ostream & os)
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
        virtual void setInitVals(double x){initVals = 0.0;}
        /**
         * Sets initVals to the input vector.
         * 
         * @param x - vector to set as initial values
         * 
         * @overrides VectorInfo::setInitVals(dvector& x)
         */
        virtual void setInitVals(dvector& x);
        /**
         * Draws a vector of random values by resampling the prior or jittering.
         * 
         * If the parameter estimation phase is \<0, the vector of initVals
         * is returned. The returned vector is guaranteed to sum to 0.
         * 
         * @param [in] rng - the random number generator
         * @param [in] vif - the variance inflation factor
         * 
         * @return - dvector of values (TODO: on what scale??)
         */
        virtual dvector drawInitVals(random_number_generator& rng, double vif);
        /**
         * Draws a jittered random value on the arithmetic scale based on the 
         * bounds set and the fraction of the range to jitter (jitFrac).
         * 
         * Note that this DOES NOT update initVal.
         * If the estimation phase is \< 0, the value of initVal is returned
         * 
         * @param rng - the random number generator
         * @param jitFrac - the fraction of the range across which to jitter
         * 
         * @return - the random number on the arithmetic scale, or the value of initVal
         */
        virtual dvector jitterInitVals(random_number_generator& rng, double jitFrac);
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
* </ul>
* 
* Protected fields:
* <ul>
* <li> nNIs - number of elements (parameters) in the associated param_init_number_vector
* <li> ppNIs - pointer to vector of pointers to NumberInfo instances
* </ul>
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
         * Gets the parameter scale types for the associated parameters.
         * 
         * @return - an adstring_array of scale types
         */
        adstring_array getScaleTypes();
        /**
         * Gets a dvector of the initial values, on the arithmetic scale, for 
         * each parameter in theassociated param_init_number_vector.
         * 
         * @return - a dvector
         */
        dvector getInitVals(void);
        /**
         * Gets a dvector of the initial values, on the parameter scale, for 
         * each parameter in theassociated param_init_number_vector.
         * 
         * @return - a dvector
         */
        virtual dvector getInitValsOnParamScales(void);
        /**
         * Calculates a vector of the log-scale prior probability based on an input vector of values.
         * 
         * The input values are presumably from the associated parameter values, 
         * but should be on the "arithmetic" scale.
         * 
         * @param x - the input dvar_vector on the arithmetic scale, presumably based on the associated parameters
         * 
         * @return - the log-scale prior, as a dvar_vector
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
         * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
         * 
         * @param x - parameter-scale values as a dvector
         * @return - arithmetic-scale values as a dvector
         */
        virtual dvector calcArithScaleVals(const dvector& x);
        /**
         * Calculates arithmetic-scale values corresponding to the input parameter-scale values.
         * 
         * @param x - parameter-scale values as dvar_vector
         * @return - arithmetic-scale values as dvar_vector
         */
        virtual dvar_vector calcArithScaleVals(const dvar_vector& x);
        /**
         * Calculates parameter-scale values corresponding to the input arithmetic-scale values.
         * 
         * @param x - the arithmetic-scale values as dvector
         * @return - the parameter-scale values as dvector
         */
        virtual dvector calcParamScaleVals(dvector& x);            
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
        virtual void setInitValsFromParamVals(const dvar_vector& x);
        virtual void setFinalValsFromParamVals(const dvar_vector& x);
        /**
         * Read values from stream in TCSAM02 format
         * 
         * @param is - input stream
         */
        virtual void read(cifstream & is);
        /**
         * Write values to stream in TCSAM02 format
         * 
         * @param os - output stream
         */
        virtual void write(std::ostream & os);
        /**
         * Write values to output stream in ADMB pin-file format
         * @param os
         */
        virtual void writeToPin(std::ostream & os);
        /**
         * Write values to output stream in R format.
         * 
         * @param os
         * @param nm - name to assign 
         * @param indent - number of tabs to indent
         */
        virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
        /**
         * Write final values to output stream in R format.
         * 
         * @param os
         */
        virtual void writeFinalValsToR(std::ostream& os);
        friend cifstream& operator >>(cifstream & is, NumberVectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, NumberVectorInfo & obj){obj.write(os);return os;}
};

/**
* BoundedNumberVectorInfo: NumberVectorInfo
* 
* Encapsulates parameter info for a param_init_bounded_number_vector as a 1-d array
* of BoundedNumberInfo objects.
* 
* Inherited class fields:
 * 
* Public fields:
* <ul>
* <li> debug - static int flag to print debugging info
* <li> name - adstring with name of associated param_init_number_vector object
* </ul>
* 
* Protected fields:
* <ul>
* <li> nNIs - number of elements (parameters) in the associated param_init_number_vector
* <li> ppNIs - pointer to vector of pointers to NumberInfo instances
* </ul>
* 
* New class fields:
*  none
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
         * Gets a pointer to the ith BoundedNumberInfo element in the vector. 
         * i must be in the interval 1<=i<=nNIs.
         * 
         * @param i - index (1-based) to BoundedNumberInfo
         * @return pointer to ith BoundedNumberInfo object
         */
        BoundedNumberInfo* operator[](int i);            
        /**
         * Get the lower bounds on the arithmetic scale for parameters as vector.
         * 
         * @return - a dvector with the arithmetic-scale lower bound for each parameter
         */
        dvector getLowerBounds();
        /**
         * Get the lower bounds on the parameter scale for parameters as vector.
         * 
         * @return - a dvector with the lower bound on the parameter scale for each parameter
         */
        dvector getLowerBoundsOnParamScales();
        /**
         * Get the upper bounds for parameters as vector.
         * 
         * @return - a dvector with the upper bound for each parameter
        */
        dvector getUpperBounds();
        /**
         * Get the upper bounds on the parameter scale for parameters as vector.
         * 
         * @return - a dvector with the upper bound on the parameter scale for each parameter
         */
        dvector getUpperBoundsOnParamScales();
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream& os, adstring nm, int indent=0);
};

/**
* VectorVectorInfo
* 
* Encapsulates parameter info for a param_init_vector_vector as a 1-d array
* of VectorInfo objects.
* 
* Public fields:
* <ul>
* <li> debug - static int flag to print debugging info
* <li> name - adstring with name of associated param_init_vector_vector object
* </ul>
* 
* Protected fields:
* <ul>
* <li> nVIs - number of elements (parameters) in the associated param_init_vector_vector
* <li> ppVIs - pointer to vector of pointers to VectorInfo instances
* </ul>
* 
* New class fields:
*  none
*/
class VectorVectorInfo {
    public:
        /** flag to write out debugging info */
        static int debug;
        /** name of vector of associated param_init_vector_vector */
        adstring name;     //
    protected:
        /** number of parameter vectors in vector */
        int nVIs;
        /** ptr to vector of ptrs to VectorInfo instances */
        VectorInfo** ppVIs;
    public:
        VectorVectorInfo(){this->name="";nVIs=0;ppVIs=0;}
        VectorVectorInfo(adstring& name){this->name=name;nVIs=0;ppVIs=0;}
        VectorVectorInfo(const char * name){this->name=name;nVIs=0;ppVIs=0;}
        ~VectorVectorInfo(){deallocate();}

        void deallocate();
        /**
         * Gets the number of parameter vectors represented by the instance.
         * 
         * @return - the number of vectors
         */
        int getSize(void){return nVIs;}
        ivector getMinIndices(void);
        ivector getMaxIndices(void);
        /**
         * Gets the parameter-scale types for each associated parameter vector.
         * 
         * @return - an adstring_array
         */
        adstring_array getScaleTypes(void);
        /**
         * Gets the estimation phase for each associated parameter vector.
         * @return 
         */
        ivector getPhases(void);
        dvector getPriorWgts(void);
        /**
         * Function to get array of labels corresponding to parameter vectors.
         * 
         * @return adstring_array of labels for associated VectorInfo instances
         */
        adstring_array getLabels(void);

        virtual VectorInfo* operator[](int i){if ((ppVIs>0)&&(i<=nVIs)) return ppVIs[i-1]; return 0;}

        virtual void read(cifstream & is);
        virtual void write(std::ostream & os);
        virtual void writeToPin(std::ostream& os);
        virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
        virtual void writeFinalValsToR(std::ostream& os);
        friend cifstream& operator >>(cifstream & is, VectorVectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, VectorVectorInfo & obj){obj.write(os);return os;}
};

/**
* BoundedVectorVectorInfo: VectorVectorInfo
* 
* Encapsulates parameter info for a param_init_bounded_vector_vector as a 1-d array
* of BoundedVectorInfo objects.
* 
* Inherited class fields:\n
* Public:
* <ul>
*   <li> debug - static int flag to print debugging info
*   <li> name - adstring with name of associated param_init_vector_vector object
* </ul>
* Protected:
* <ul>
*   <li> nVIs - number of elements (parameters) in the associated param_init_vector_vector </li>
*   <li> ppVIs - pointer to vector of pointers to VectorInfo instances </li>
* </ul>
* 
* New class fields:
*  none
*/
class BoundedVectorVectorInfo: public VectorVectorInfo {
    public:
        static int debug;
    public:
        BoundedVectorVectorInfo():VectorVectorInfo(){}
        BoundedVectorVectorInfo(adstring& name):VectorVectorInfo(name){}
        BoundedVectorVectorInfo(const char * name):VectorVectorInfo(name){}

        /**
         * Gets a vector of the lower bounds (on the arithmetic scale) for the associated BoundedVectorInfo instances.
         * 
         * @return - a dvector of the lower bounds
         */
        dvector getLowerBounds(void);
        /**
         * Gets a vector of the lower bounds (on the parameter vector scales) for the associated BoundedVectorInfo instances.
         * 
         * @return - a dvector of the lower bounds, on the parameter vector scales
         */
        dvector getLowerBoundsOnParamScales(void);
        /**
         * Gets a vector of the upper bounds (on the arithmetic scale) for the associated BoundedVectorInfo instances.
         * 
         * @return - a dvector of the upper bounds
         */
        dvector getUpperBounds(void);
        /**
         * Gets a vector of the upper bounds (on the parameter vector scales) for the associated BoundedVectorInfo instances.
         * 
         * @return - a dvector of the upper bounds, on the parameter vector scales
         */
        dvector getUpperBoundsOnParamScales(void);
        /**
         * Returns a pointer to the ith (1s-based) BoundedVectorInfo instance
         * encapsulated by the BoundedVectorVectorInfo instance.
         * 
         * @param i
         * @return 
         * 
         * @overrides VectorVectorInfo::operator[](int i)
         */
        virtual BoundedVectorInfo* operator[](int i){if ((ppVIs>0)&&(i<=nVIs)) return static_cast<BoundedVectorInfo*>(ppVIs[i-1]); return 0;}

        virtual void read(cifstream & is);
        virtual void write(std::ostream & os);
        virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
        virtual void writeFinalValsToR(std::ostream& os);
        friend cifstream& operator >>(cifstream & is, BoundedVectorVectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, BoundedVectorVectorInfo & obj){obj.write(os);return os;}
};

/**
* DevsVectorVectorInfo: BoundedVectorVectorInfo
* 
*  Encapsulates characteristics for a param_init_bounded_vector_vector
*  to be used as the basis for a vector of devs vectors. Each devs vector has the
*  property that the sum of the elements of the vector should be 0.
*  This should be implemented by enforcing sum(v[])=0 on the arithmetic
*  scale in the likelihood, for each devs vector v.
*  Otherwise, DevsVectorVectorInfo is identical to BoundedVectorVectorInfo,
*  of which it is a subclass.
* 
* Inherited class fields:
* Public:
* <ul>
* <li> debug - static int flag to print debugging info
* <li> name - adstring with name of associated param_init_vector_vector object
* </ul>
* 
* Protected:
* <ul>
* <li> nVIs - number of elements (parameters) in the associated param_init_vector_vector
* <li> ppVIs - pointer to vector of pointers to VectorInfo instances
* </ul>
* 
* New class fields:
*  none
*/
class DevsVectorVectorInfo: public BoundedVectorVectorInfo {
    public:
        static int debug;
    public:
        DevsVectorVectorInfo():BoundedVectorVectorInfo(){}
        DevsVectorVectorInfo(adstring& name):BoundedVectorVectorInfo(name){}
        DevsVectorVectorInfo(const char * name):BoundedVectorVectorInfo(name){}

        /**
         * Returns a pointer to the ith (1s-based) DevsVectorInfo instance
         * encapsulated by the DevsVectorVectorInfo instance.
         * 
         * @param i
         * @return 
         * 
         * @overrides BoundedVectorVectorInfo::operator[](int i)
         */
        virtual DevsVectorInfo* operator[](int i){if ((ppVIs>0)&&(i<=nVIs)) return (DevsVectorInfo*) ppVIs[i-1]; return 0;}

        virtual void read(cifstream & is);
        friend cifstream& operator >>(cifstream & is, DevsVectorVectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, DevsVectorVectorInfo & obj){obj.write(os);return os;}
};

#endif //MODELPARAMETERINFOTYPES_HPP
