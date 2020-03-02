/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ModelParameterVectorInfoTypes.hpp
 * Author: WilliamStockhausen
 *
 * Created on February 28, 2020, 4:28 PM
 */

#pragma once
#ifndef MODELPARAMETERVECTORINFOTYPES_HPP
#define MODELPARAMETERVECTORINFOTYPES_HPP

#include <admodel.h>
#include "ModelIndexBlocks.hpp"
#include "ModelFunctions.hpp"
#include "ModelConfiguration.hpp"
#include "ModelConstants.hpp"

/**
* VectorInfo
* 
*  Encapsulates characteristics for a vector of parameters on a possibly 
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
* <li> phases   - vector of estimation phases for the corresponding parameters
* <li> initVals - vector of initial values, on the arithmetic scale
* <li> finlVals - vector of final values, on th arithmetic scale
* <li> scaleType   - parameter scale type
* <li> phase       - default phase to start estimation of associated parameter
* <li> resample    - flag to resample initial value from prior
* <li> priorWgt    - weight assigned to prior in likelihood
* <li> priorType   - name identifying pdf for prior
* <li> priorParams - dvar_vector of parameters (on the arithmetic scale) for the prior
* <li> priorConsts - dvector of constants (on the arithmetic scale) for the prior
* <li> readPhases  - flag to read vector of estimation phases from an input stream
* <li> readVals    - flag to read vector of values from an input stream
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
        /* flag to read phases */
        int readPhases;
        /* vector of phases, by vector element */
        ivector phases;
        /* flag to read initial values */
        int readVals;
    protected:
        /* flag indicating parameter scale */
        int scaleType; 
        /* arithmetic-scale value used as the initial value for the associated parameter */
        double initVal;
        /* final arithmetic-scale value from the associated parameter (for output to R) */
        double finlVal;     
        /* default phase in which to start estimating associated parameter vector */
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
        /* vector of initial values, by vector element */
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
         * Gets the phase for each element of the associated param_init_number_vector
         * 
         * @return - estimation phase for each element of the parameter vector
         */
        ivector getPhases(){return phases;}
        
        /**
         * Sets the phases for the associated param_init_vector_vector
         * 
         * @param - ivector of phases
         */
        void setPhases(const ivector& _phases){phases = _phases;}
        
        /**
         * Reads a vector of phases from an input filestream.
         * 
         * @param [in] is - the input filestream
         */
        virtual void readPhaseVector(cifstream& is){is>>phases;}
        
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
         * Sets all phases for the associated parameters
         * 
         * @param - value for all phases
         */
        void setPhase(int _phase){phase = _phase; phases = _phase;}
        
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
         *  <li> readPhases - an adstring flag to subsequently read a vector of estimation phases
         *  <li> readVals - an adstring flag to subsequently read a vector of initial values
         *  <li> initVal
         *  <li> scaleType
         *  <li> phase (default phase)
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
         *  <li> readPhases - an adstring flag to subsequently read a vector of phases
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
////////////////////////////////////////////////////////////////////////////////
/**
* BoundedVectorInfo : VectorInfo
* 
*  Encapsulates characteristics for a bounded ADMB parameter vector
*  on a possibly-transformed (logit, lognormal, probit) scale.
*  Basically, a bounded parameter version of VectorInfo.
* 
* Inherited class fields:
* <ul>
* <li> N        -  size of parameter vector
* <li> idxType  - adstring with type of indices (e.q., YEAR)
* <li> ptrIB    - pointer to IndexBlock object defining indices
* <li> readPhases    - flag to read vector of values from input stream
* <li> phases - vector of parameter estimation phases
* <li> readVals    - flag to read vector of values from input stream
* <li> initVals - vector of initial values, on the arithmetic scale
* <li> finlVals - vector of final values, on th arithmetic scale
* <li> scaleType   - parameter scale type
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
         * @param [in] name - adstring for name of associated ADMB parameter vector
         */
        BoundedVectorInfo(adstring& name):VectorInfo(name){}

        /**
         * Class constructor.
         * 
         * @param [in] name - adstring for name of associated ADMB parameter vector
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
         * 
         * @overrides VectorInfo::getInitValsOnParamScale(dvector& x)
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
         * 
         * @overrides VectorInfo::calcLogPrior(dvar_vector& x)
         */
        virtual dvar_vector calcLogPrior(dvar_vector& x);

        /**
         * Reads a vector initial values from an input stream.
         * Any values outside the bounds are set close to the
         * the nearest bound.
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
         *  <li> readPhases - an adstring flag to subsequently read a vector of estimation phases
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
};//BoundedVectorInfo : VectorInfo
////////////////////////////////////////////////////////////////////////////////
/**
* DevsVectorInfo : BoundedVectorInfo
* 
*  Encapsulates characteristics for a devs parameter vector, which has the
*  property that the sum of the elements of the vector should be 0.
*  This should be implemented by enforcing sum(v[])=0 on the arithmetic
*  scale in the likelihood.
 * 
*  DevsVectorInfo is identical to BoundedVectorInfo other than in how 
 * initial values are set, drawn, or jittered in order to enforce the devs
 * constraint.
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
};//DevsVectorInfo : BoundedVectorInfo
////////////////////////////////////////////////////////////////////////////////
/**
* VectorVectorInfo
* 
* Encapsulates parameter info as a 1-d array
* of VectorInfo objects for a param_init_number_vector  object that 
 * substitutes for a param_init_vector_vector object.
* 
* Public fields:
* <ul>
* <li> debug - static int flag to print debugging info
* <li> name - adstring with name of associated param_init_number_vector object
* </ul>
* 
* Protected fields:
* <ul>
* <li> nVIs - number of parameter vectors represented in the associated param_init_number_vector
* <li> mnIdxs - ivector of minimum index for each parameter vector 
* <li> mxIdxs - ivector of maximum index for each parameter vector 
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
        /** number of parameter vectors represented in the VectorVector */
        int nVIs;
        /** total number of all parameters (i.e., size of the associated param_init_number_vector) */
        int npT; 
        /** ptr to vector of ptrs to VectorInfo instances */
        VectorInfo** ppVIs;
        /** ivector of minimum index for each parameter vector represented in the VectorVector */
        ivector mnIdxs;
        /** ivector of maximum index for each parameter vector represented in the VectorVector */
        ivector mxIdxs;
    public:
        VectorVectorInfo(){this->name="";nVIs=0;ppVIs=0;npT=0;}
        VectorVectorInfo(adstring& name){this->name=name;nVIs=0;ppVIs=0;npT=0;}
        VectorVectorInfo(const char * name){this->name=name;nVIs=0;ppVIs=0;npT=0;}
        ~VectorVectorInfo(){deallocate();}

        void deallocate();
        /**
         * Gets the total number of parameters for all associated devs vectors
         * (i.e., for the associated param_init_number_vector)
         * 
         * @return the total number of parameters
         */
        int getNumParameters(void){return npT;}
        /**
         * Gets the number of parameter vectors represented by the instance.
         * 
         * @return - the number of vectors
         */
        int getSize(void){return nVIs;}
        /**
         * Gets the ivector of phases for the associated param_init_number_vector (a
         * concatenation of the phases ivectors for all associated devs vectors).
         * 
         * @return ivector(1,npT) of phases for each element of the associated param_init_number_vector
         */
        ivector getParameterPhases(void);
        /**
         * Convenience method to set estimation phases for elements of all associated 
         * devs vectors (and thus the associated param_init_number_vector).
         * 
         * @param phases - ivector(1,npT) to use to set phases
         */
        void setParameterPhases(const ivector& phases);
        /**
         * Gets the min index in the param_number_vector for each parameter vector represented by the instance.
         * 
         * @return - the number of vectors
         */
        ivector getMinIndices(void){return mnIdxs;}
        /**
         * Gets the max index in the param_number_vector for each parameter vector represented by the instance.
         * 
         * @return - the number of vectors
         */
        ivector getMaxIndices(void){return mxIdxs;}
        /**
         * Gets the parameter-scale types for each associated parameter vector.
         * 
         * @return - an adstring_array
         */
        adstring_array getScaleTypes(void);
        /**
         * Gets the prior weights for each associated parameter vector.
         * @return 
         */
        dvector getPriorWgts(void);
        /**
         * Function to get array of labels corresponding to parameter vectors.
         * 
         * @return adstring_array of labels for associated VectorInfo instances
         */
        adstring_array getLabels(void);

        /**
         * Get a pointer to the ith (1's-based) VectorInfo object
         * 
         * @param i - 1's based index
         * 
         * @return pointer to the ith VectorInfo object
         */
        virtual VectorInfo* operator[](int i){if ((ppVIs)&&(i<=nVIs)) return ppVIs[i-1]; return 0;}

        /**
         * Read parameter info for VectorVectorInfo instance.
         * 
         * @param is - input stream to read from
         */
        virtual void read(cifstream & is);
        /**
         * Read optional phases and initial values, as required.
         * 
         * @param is - input stream to read from
         */
        void readOptional(cifstream & is);
        /**
         * Allocate mnIdxs, mxIdxs, and determine npT.
         */
        void allocateIndices();

        virtual void write(std::ostream & os);
        virtual void writeToPin(std::ostream& os);
        virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
        virtual void writeFinalValsToR(std::ostream& os);
        friend cifstream& operator >>(cifstream & is, VectorVectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, VectorVectorInfo & obj){obj.write(os);return os;}
};//VectorVectorInfo
////////////////////////////////////////////////////////////////////////////////
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
         * Returns a pointer to the ith (1's-based) BoundedVectorInfo instance
         * encapsulated by the BoundedVectorVectorInfo instance.
         * 
         * @param i - 1's based index
         * @return 
         * 
         * @overrides VectorVectorInfo::operator[](int i)
         */
        virtual BoundedVectorInfo* operator[](int i){if ((ppVIs)&&(i<=nVIs)) return static_cast<BoundedVectorInfo*>(ppVIs[i-1]); return 0;}

        virtual void read(cifstream & is);
        virtual void write(std::ostream & os);
        virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
        virtual void writeFinalValsToR(std::ostream& os);
        friend cifstream& operator >>(cifstream & is, BoundedVectorVectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, BoundedVectorVectorInfo & obj){obj.write(os);return os;}
};//BoundedVectorVectorInfo: VectorVectorInfo
////////////////////////////////////////////////////////////////////////////////
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
        /** flag to print debugging info */
        static int debug;
    public:
        DevsVectorVectorInfo():BoundedVectorVectorInfo(){}
        DevsVectorVectorInfo(adstring& name):BoundedVectorVectorInfo(name){}
        DevsVectorVectorInfo(const char * name):BoundedVectorVectorInfo(name){}

        /**
         * Returns a pointer to the ith (1's-based) DevsVectorInfo instance
         * encapsulated by the DevsVectorVectorInfo instance.
         * 
         * @param i - 1's based index
         * @return 
         * 
         * @overrides BoundedVectorVectorInfo::operator[](int i)
         */
        virtual DevsVectorInfo* operator[](int i){if ((ppVIs)&&(i<=nVIs)) return (DevsVectorInfo*) ppVIs[i-1]; return 0;}

        virtual void read(cifstream & is);
        //virtual void write(ostream & os);//BoundedVectorVectorInfo::write(...) should work
        friend cifstream& operator >>(cifstream & is, DevsVectorVectorInfo & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os, DevsVectorVectorInfo & obj){obj.write(os);return os;}
};//DevsVectorVectorInfo: BoundedVectorVectorInfo
////////////////////////////////////////////////////////////////////////////////

#endif /* MODELPARAMETERVECTORINFOTYPES_HPP */

