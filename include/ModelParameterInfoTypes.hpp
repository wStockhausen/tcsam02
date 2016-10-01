/*------------------------------------------------------------------------------
 * Includes:
 *      NumberInfo
 *      BoundedNumberInfo
 *      NumberVectorInfo
 *      BoundedNumberVectorInfo
 *      VectorInfo
 *      VectorVectorInfo
 *      BoundedVectorInfo
 *      BoundedVectorVectorInfo
 *      DevsVectorInfo
 *      DevsVectorVectorInfo
 *----------------------------------------------------------------------------*/
#pragma once
#ifndef MODELPARAMETERINFOTYPES_HPP
    #define MODELPARAMETERINFOTYPES_HPP

    #include <admodel.h>
    #include "ModelIndexBlocks.hpp"
    #include "ModelFunctions.hpp"
    #include "ModelConfiguration.hpp"

/*------------------------------------------------------------------------------
 * NumberInfo
 *  Encapsulates characteristics for a param_init_number
 *      initVal
 *      phase
 *      resample
 *      priorWgt
 *      priorType
 *      priorParams
 *      priorConsts
 *----------------------------------------------------------------------------*/
    class NumberInfo {
        public:
            static int debug;
        protected:
            double initVal;     //initial value on "natural" scale
            double finlVal;     //final value (for output to R)
            int phase;          //phase in which to turn on parameter
            double priorWgt;    //weight to assign to prior probability
            adstring priorType; //pdf type for prior
            dvector priorParams;//specified on "natural" scale
            dvector priorConsts;//specified on "natural" scale
        public:
            adstring      name;
            bool          resample; //flag to do resampling of initial values
            ModelPDFInfo* pMPI;
        public:
            NumberInfo(){this->name="";pMPI=0;}
            NumberInfo(adstring& name){this->name=name;pMPI=0;}
            NumberInfo(const char* name){this->name=name;pMPI=0;}
            ~NumberInfo(){delete pMPI;pMPI=0;}
            
            int      getPhase(){return phase;}
            double   getPriorWgt(){return priorWgt;}
            adstring getPriorType(){return priorType;}
            double   getInitVal(){return initVal;}
            dvariable calcLogPrior(prevariable& x);             
            
            virtual double drawInitVal(random_number_generator& rng, double vif);//draw initial value by resampling prior
            virtual void   setInitVal(double x){initVal=x;}          
            virtual void   setInitVal(param_init_number& x){initVal=value(x);}
            virtual double getFinalVal(){return finlVal;}
            virtual void   setFinalVal(param_init_number& x){finlVal=value(x);}
            
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream & os);
            
            friend cifstream& operator >>(cifstream & is, NumberInfo & obj){obj.read(is); return is;}
            friend std::ostream& operator <<(std::ostream & os, NumberInfo & obj){obj.write(os); return os;}
        protected:
            void writeToR1(std::ostream & os);
            void setPriorType(adstring & prior);
    };

/*------------------------------------------------------------------------------
 * BoundedNumberInfo : NumberInfo
 *  Encapsulates characteristics for a param_init_bounded_number parameter
 *      jitter
 *      lower
 *      upper
 *---------------------------------------------------------------------------*/
    class BoundedNumberInfo : public NumberInfo {
        public:
            static int debug;
        protected:
            double lower;//lower bound
            double upper;//upper bound
        public:
            int jitter;  //flag to jitter
        public:
            BoundedNumberInfo():NumberInfo(){}
            BoundedNumberInfo(adstring& name):NumberInfo(name){}
            BoundedNumberInfo(const char * name):NumberInfo(name){}
            ~BoundedNumberInfo(){}
            
            double getLowerBound(){return lower;}
            double getUpperBound(){return upper;}
           
            virtual double drawInitVal(random_number_generator& rng,double vif);//draw initial value based on jitter or resampling prior
            virtual void setInitVal(double x);           
            virtual void setInitVal(param_init_number& x) {std::cout<<"cannot use setInitVal(param_init_number) for BoundedNumberInfo. Aborting..."<<endl; exit(-1);}
            virtual void setFinalVal(param_init_number& x){std::cout<<"cannot use setFinalVal(param_init_number) for BoundedNumberInfo. Aborting..."<<endl; exit(-1);}
            virtual void setInitVal(param_init_bounded_number& x){initVal=value(x);}
            virtual void setFinalVal(param_init_bounded_number& x){finlVal=value(x);}
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream& os);
        protected:
            virtual void writeToR1(std::ostream & os);
    };

/*------------------------------------------------------------------------------
 * VectorInfo: NumberInfo
 *  Encapsulates characteristics for a param_init_vector
 *      idxType  - adstring with type of indices (e.q., YEAR)
 *      ptrIB    - pointer to IndexBlock object defining indices
 *      readVals -  flag to read vector of values
 *      initVals - vector of initial values
 *----------------------------------------------------------------------------*/
    class VectorInfo: public NumberInfo {
        public:
            static int debug;
            int N;       //size of initVals vector
            int readVals;//flag to read initial values
        protected:
            adstring idxType;//index type (ie., model dimension, e.g., YEAR)
            IndexBlock* ptrIB;//pointer to index block
            dvector initVals;//initial values
            dvector finlVals;//final values (for output to R)
        public:
            adstring      name;
            ModelPDFInfo* pMPI;
        public:
            VectorInfo():NumberInfo(){ptrIB=0;}
            VectorInfo(adstring& name):NumberInfo(name){ptrIB=0;}
            VectorInfo(const char* name):NumberInfo(name){ptrIB=0;}
            ~VectorInfo(){if (ptrIB) delete ptrIB;}
            
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
             * Sets initial values to input constant.
             * @param x
             */
            void setInitVal(double x){initVals = x;}
            /**
             * Return initial values vector (indices run 1:N).
             * @return 
             */
            dvector getInitVals(){return initVals;}
            /**
             * Sets initial values to input dvector element-by-element..
             * @param x - vector of initial values
             */
            virtual void setInitVals(dvector& x){initVals=x;}     
            /**
             * Sets initial values to input vector element-by-element..
             * @param x
             */
            virtual void setInitVals(param_init_vector & x){initVals=value(x);}     
            /**
             * Returns a copy of the vector of final values.
             * @return dvector of final values
             */
            virtual dvector getFinalVals(){return finlVals;}
            /**
             * Sets final values for output to R.
             * @param x - vector of final values
             */
            virtual void setFinalVals(param_init_vector & x){finlVals.allocate(x.indexmin(),x.indexmax()); finlVals=value(x);}     
            /**
             * Calculates vector of log prior values based on input vector.
             * @return 
             */
            dvar_vector calcLogPrior(dvar_vector& x);
            /**
             * Draw inital values based on jittering or resampling the prior.
             * 
             * @param rng - random number generator object
             * @param vif - variance inflation factor
             * @return 
             */
            virtual dvector drawInitVals(random_number_generator& rng, double vif);//draw initial values by resampling prior
            /**
             * Reads initial values from filestream.
             * @param is
             */
            virtual void readInitVals(cifstream& is){is>>initVals;}
            /**
             * Reads vector info from input filestream.
             * @param is
             */
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream & os);
            virtual void writeFinalValsToR(std::ostream & os);
    };

/*----------------------------------------------------------------------------\n
  * BoundedVectorInfo : BoundedNumberInfo\n
  *  Encapsulates characteristics for a vector of \n
  *  param_init_bounded_number parameters.\n
  *      BoundedNumberInfo characteristics 
 *      idxType  - adstring with type of indices (e.q., YEAR)
 *      ptrIB    - pointer to IndexBlock object defining indices
 *      readVals -  flag to read vector of values
 *      initVals - vector of initial values
*---------------------------------------------------------------------------*/
    class BoundedVectorInfo : public BoundedNumberInfo {
        public:
            static int debug;
            int N;//size of initVals vector
            int readVals;     //flag to read initial values
        protected:
            adstring idxType; //index type (ie., model dimension, e.g., YEAR)
            IndexBlock* ptrIB;//pointer to index block
            dvector initVals; //initial values
            dvector finlVals; //final values (for output to R)
        public:
            BoundedVectorInfo():BoundedNumberInfo(){}
            BoundedVectorInfo(adstring& name):BoundedNumberInfo(name){}
            BoundedVectorInfo(const char * name):BoundedNumberInfo(name){}
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
             * Calculates the (ln-scale) prior probabilities (not NLLS!)vcorresponding to 
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
             * Sets initial values to input value.
             * 
             * @param x - value to set as initial values
             */
            virtual void    setInitVal(double x){initVals = x;}            
            /**
             * Sets initial values to input vector element-by-element.
             * 
             * @param x - vector to set as initial values
             */
            virtual void setInitVals(dvector& x);
            /**
             * Sets initial values to input vector element-by-element..
             * @param x
             */
            virtual void setInitVals(param_init_bounded_vector & x){initVals=value(x);} 
            /**
             * Returns a copy of the vector of final values.
             * @return dvector of final values
             */
            virtual dvector getFinalVals(){return finlVals;}
            /**
             * Sets final values for output to R.
             * @param x
             */
            virtual void setFinalVals(param_init_bounded_vector & x){finlVals.allocate(x.indexmin(),x.indexmax()); finlVals=value(x);} 
            /**
             * Reads initial values from an input stream.
             * @param is
             */
            virtual void readInitVals(cifstream & is);
            /**
             * Draw initial values from prior distribution or by jittering.
             * @param rng
             * @param vif
             * @return 
             */
            virtual dvector drawInitVals(random_number_generator& rng, double vif);//draw initial values by resampling prior
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream & os);
            virtual void writeFinalValsToR(std::ostream & os);
    };

/*----------------------------------------------------------------------------\n
  * DevsVectorInfo : BoundedVectorInfo\n
  *  Encapsulates characteristics for a vector of  param_init_bounded_number  \n
  *  parameters to be used as the basis for a devs vector, which has the      \n
  *  property that the sum of the elements of the vector is identically 0.    \n
  *  This is implemented by setting v(N) = -sum(v[1,(N-1)]).            \n
  *      BoundedNumberInfo characteristics 
  *      idxType  - adstring with type of indices (e.q., YEAR)
  *      ptrIB    - pointer to IndexBlock object defining indices
  *      readVals -  flag to read vector of values
  *      initVals - vector of initial values
*---------------------------------------------------------------------------*/
    class DevsVectorInfo : public BoundedVectorInfo {
        public:
            static int debug;
        public:
            DevsVectorInfo():BoundedVectorInfo(){}
            DevsVectorInfo(adstring& name):BoundedVectorInfo(name){}
            DevsVectorInfo(const char * name):BoundedVectorInfo(name){}
            ~DevsVectorInfo(){}
                      
            /**
             * Sets initial values to 0, no matter what @param x is.
             * @param x
             */
            void setInitVal(double x){initVals = 0.0;}
            
            /**
             * Sets initial values 1:(N-1) to those of the vector x, but
             * sets the value for element N to -sum(initVals(1,N-1)) so
             * the sum over all elements is 0. x may have size N-1.
             * 
             * @param x - dvector of initial values
             */
            virtual void setInitVals(dvector& x);
            /**
             * Sets initial values 1:(N-1) to those of the vector x, but
             * sets the value for element N to -sum(initVals(1,N-1)) so
             * the sum over all elements is 0. x may have size N-1.
             * 
             * @param x - param_init_bounded_vector of initial values
             */
            virtual void setInitVals(param_init_bounded_vector & x);     
            /**
             * Sets final values 1:(N-1) to those of the vector x, but
             * sets the value for element N to -sum(initVals(1,N-1)) so
             * the sum over all elements is 0. x may have size N-1.
             * 
             * @param x - param_init_bounded_vector of final values
             */
            virtual void setFinalVals(param_init_bounded_vector & x);     
            /**
             * Reads initial values 1:N from a file stream and sets the 
             * values for 1:(N-1) to those of the read-in vector x, but sets
             * the value for element N to -sum(initVals(1,N-1)) so the
             * sum over all elements is 0.
             * 
             * @param is - filestream object from which to read initial values
             */
            void    readInitVals(cifstream & is);
            
            virtual dvector drawInitVals(random_number_generator& rng, double vif);//draw initial values by resampling prior
            virtual void read(cifstream & is);
            virtual void writeToR(std::ostream & os);
            virtual void writeFinalValsToR(std::ostream & os);
        protected:
            void calcDevs(void);
    };

//--------------------------------------------------------------------------------
//          NumberVectorInfo
//  Encapsulates info for a param_init_number_vector as a 1-d array of NumberInfo's
//--------------------------------------------------------------------------------
    class NumberVectorInfo {
        public:
            static int debug;
        public:
            adstring name;     //parameter vector name
        protected:
            int nNIs;          //number of parameters in vector
            NumberInfo** ppNIs;//ptr to vector of ptrs to NumberInfo instances
        public:
            NumberVectorInfo(){this->name="";nNIs=0;ppNIs=0;}
            NumberVectorInfo(adstring& name){this->name=name;nNIs=0;ppNIs=0;}
            NumberVectorInfo(const char * name){this->name=name;nNIs=0;ppNIs=0;}
            ~NumberVectorInfo(){deallocate();}
            
            NumberInfo* operator[](int i){if (ppNIs&&(i<=nNIs)) return ppNIs[i-1]; return 0;}
            void deallocate();
            int getSize(void){return nNIs;}
            ivector getPhases(void);
            dvector getPriorWgts(void);
            dvector getInitVals(void);
            dvar_vector calcLogPriors(dvar_vector & pv);         
            
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

//--------------------------------------------------------------------------------
//          BoundedNumberVectorInfo
//  Encapsulates info for a param_init_bounded_number_vector as a 1-d array of BoundedNumberInfo's
//--------------------------------------------------------------------------------
    class BoundedNumberVectorInfo: public NumberVectorInfo {
        public:
            static int debug;
        public:
            BoundedNumberVectorInfo():NumberVectorInfo(){}
            BoundedNumberVectorInfo(adstring& name):NumberVectorInfo(name){}
            BoundedNumberVectorInfo(const char * name):NumberVectorInfo(name){}
            ~BoundedNumberVectorInfo(){deallocate();}
            
            BoundedNumberInfo* operator[](int i);
            dvector getLowerBounds();
            dvector getUpperBounds();
            
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
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
            virtual void writeFinalValsToR(std::ostream& os);
            friend cifstream& operator >>(cifstream & is, BoundedVectorVectorInfo & obj){obj.read(is);return is;}
            friend std::ostream& operator <<(std::ostream & os, BoundedVectorVectorInfo & obj){obj.write(os);return os;}
    };

//--------------------------------------------------------------------------------
//          DevsVectorVectorInfo
//  Encapsulates info for a param_init_bounded_vector_vector as a 1-d array of DevsVectorInfo's
//--------------------------------------------------------------------------------
    class DevsVectorVectorInfo {
        public:
            static int debug;
        public:
            adstring name;         //name of vector of devs vectors
        protected:
            int nVIs;              //number of devs vectors in vector
            DevsVectorInfo** ppVIs;//ptr to vector of ptrs to DevsVectorInfo instances
        public:
            DevsVectorVectorInfo(){this->name="";nVIs=0;ppVIs=0;}
            DevsVectorVectorInfo(adstring& name){this->name=name;nVIs=0;ppVIs=0;}
            DevsVectorVectorInfo(const char * name){this->name=name;nVIs=0;ppVIs=0;}
            ~DevsVectorVectorInfo(){deallocate();}
            
            DevsVectorInfo* operator[](int i){if ((ppVIs>0)&&(i<=nVIs)) return ppVIs[i-1]; return 0;}
            void deallocate();
            int getSize(void){return nVIs;}
            ivector getMinIndices(void);
            ivector getMaxIndices(void);
            ivector getPhases(void);
            dvector getPriorWgts(void);
            dvector getLowerBounds(void);
            dvector getUpperBounds(void);
            virtual void read(cifstream & is);
            virtual void write(std::ostream & os);
            virtual void writeToR(std::ostream& os, adstring nm, int indent=0);
            virtual void writeFinalValsToR(std::ostream& os);
            friend cifstream& operator >>(cifstream & is, DevsVectorVectorInfo & obj){obj.read(is);return is;}
            friend std::ostream& operator <<(std::ostream & os, DevsVectorVectorInfo & obj){obj.write(os);return os;}
    };

#endif //MODELPARAMETERINFOTYPES_HPP
