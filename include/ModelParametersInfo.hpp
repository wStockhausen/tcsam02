/* 
 * File:   ModelParametersInfo.hpp
 * Author: WilliamStockhausen
 *
 * Created on February 12, 2014, 3:04 PM
 */

#ifndef MODELPARAMETERSINFO_HPP
#define	MODELPARAMETERSINFO_HPP

#include <admodel.h>
#include "ModelConfiguration.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelParameterInfoTypes.hpp"

/*------------------------------------------------------------------------------\n
 * ParameterGroupInfo\n
 * Description: Abstract base class for concrete ..ParametersInfo classes \n
 * except ModelParametersInfo.\n
 *----------------------------------------------------------------------------*/
class ParameterGroupInfo{
    public:
        static int debug;
    public:
        adstring name;        //name of ParameterGroup
        int nIVs;             //number of index variables
        adstring_array lblIVs;//names for index variables
        int nPVs;             //number of parameter variables
        adstring_array lblPVs;//names for parameter variables
        adstring_array dscPVs;//descriptions for parameter variables
        int nXIs;             //number of extra indices/values
        adstring_array lblXIs;//names for extra indices/values
        
        int nIBSs;              //number of index variables (IVs) that are blocks
        ivector ibsIdxs;        //indices corresponding to index variables (IVs) that are blocks
        IndexBlockSet** ppIBSs;//pointer to a vector of pointers to IndexBlockSet objects
        
        int nPCs;  //number of rows in parameter combinations matrix
        imatrix in;//input parameter combinations matrix (all integers)
        dmatrix xd;//extra values by parameter combination as doubles
        imatrix** ppIdxs; //pointer to array of pointers to indices matrices for use via getModelIndices(pc)     
    public:
        ParameterGroupInfo();
        ~ParameterGroupInfo();
        
        /* 
         * Returns a pointer to the index block set identified by "type".
         * Inputs:
         *  adstring type:  "type" identifying index block set to return
         * Returns:
         *  pointer to the identified IndexBlockSet
         */
        IndexBlockSet* getIndexBlockSet(adstring type);
        
        /*******************************************
         * get indices for parameter combination.
         * @param pc : id for desired parameter combination
         ******************************************/
        ivector getPCIDs(int pc);
        /*******************************************
         * get model indices for parameter combination.
         * @param pc: id for desired parameter combination
         ******************************************/
        imatrix getModelIndices(int pc);
        
        virtual void read(cifstream & is);
        virtual void write(std::ostream & os);
        virtual void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, ParameterGroupInfo & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, ParameterGroupInfo & obj){obj.write(os); return os;}
    protected:
        void createIndexBlockSets(void);
        void createIndices(void);
        BoundedNumberVectorInfo* read(cifstream& is, adstring& lbl, BoundedNumberVectorInfo* pBNVI);
        BoundedVectorVectorInfo* read(cifstream& is, adstring& lbl, BoundedVectorVectorInfo* pBVVI);
        DevsVectorVectorInfo*    read(cifstream& is, adstring& lbl, DevsVectorVectorInfo* pDVVI);
};

namespace tcsam{
    /**
     * Function to convert parameter combinations to R dimensions 
     * @param pgi
     * @return wts::adstring_matrix
     */
    wts::adstring_matrix convertPCs(ParameterGroupInfo * pgi);
}//namespace tcsam

/*------------------------------------------------------------------------------
 * RecruitmentInfo\n
 * Encapsulates the following recruitment-related parameters:\n
 *  pLnR   : mean ln-scale recruitment\n
 *  pLnRCV : recruitment cv's\n
 *  pLgtRX : logit-scale sex ratio for males\n
 *  pLnRa  : size-at-recruitment parameter\n
 *  pLnRb  : size-at-recruitment parameter\n
 *  pvLnRDevs: ln-scale annual recruitment devs
 * Notes:
 *  1. YEAR_BLOCK is the index variable for the parameters
 *----------------------------------------------------------------------------*/
class RecruitmentInfo: public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"recruitment"
    public:        
        BoundedNumberVectorInfo* pLnR;
        BoundedNumberVectorInfo* pLnRCV;
        BoundedNumberVectorInfo* pLgtRX;
        BoundedNumberVectorInfo* pLnRa;
        BoundedNumberVectorInfo* pLnRb;
        DevsVectorVectorInfo* pDevsLnR; //parameter vectors for annual recruitment devs
        
        RecruitmentInfo();
        ~RecruitmentInfo();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
        
        friend cifstream& operator >>(cifstream & is, RecruitmentInfo & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, RecruitmentInfo & obj){obj.write(os); return os;}
};

/*------------------------------------------------------------------------------
 * NaturalMortalityInfo\n
 * Encapsulates the following recruitment-related parameters:\n
 *   pLnM   : base natural mortality (mature males)
 *   pLnDMT : main ln-scale temporal offsets
 *   pLnDMX : female offsets
 *   pLnDMM : immature offsets
 *   pLnDMXM: female-immature offsets    
 * Notes:
 *  1. index variables for parameters
 *      a. YEAR_BLOCK
 *      b. SEX
 *      c. MATURITY
 *      d. SHELL
*----------------------------------------------------------------------------*/
class NaturalMortalityInfo : public ParameterGroupInfo {
    public:
        /* flag to print debugging info (if >0)*/
        static int debug;
    protected:
        /* parameter group name */
        static adstring NAME;//"natural_mortality"
    public:
        /* reference size for size-specific mortality */
        double zRef;
        /* info for base (immature male) natural mortality rates */
        BoundedNumberVectorInfo* pLnM;
        /* info for ln-scale temporal offset */
        BoundedNumberVectorInfo* pLnDMT;
        /* info for ln-scale female offset */
        BoundedNumberVectorInfo* pLnDMX;
        /* info for ln-scale mature offset */
        BoundedNumberVectorInfo* pLnDMM;
        /* info for ln-scale mature female offset */
        BoundedNumberVectorInfo* pLnDMXM;
        
        /**
         * Class constructor.
         */
        NaturalMortalityInfo();
        /**
         * Class destructor
         */
        ~NaturalMortalityInfo();
        
        /**
         * Read from text file input stream.
         * 
         * @param is - input stream
         */
        void read(cifstream & is);
        /**
         * Write to output stream in ADMB format
         * @param os - output stream
         */
        void write(std::ostream & os);
        /**
         * Write component info to output stream as
         * R-format list.
         * 
         * @param os - output stream
         */
        void writeToR(std::ostream & os);
};

/*------------------------------------------------------------------------------
 * GrowthInfo\n
 * Encapsulates the following recruitment-related parameters:\n
 *   pLnGrA : ln-scale mean growth coefficient "a"
 *   pLnGrB : ln-scale mean growth exponent "b"
 *   pLnGrBeta : scale factor for growth transition matrix
 * Notes:
 *  1. YEAR_BLOCK, SEX are the index variables for the parameters
 *----------------------------------------------------------------------------*/
class GrowthInfo : public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"growth"
    public:
        BoundedNumberVectorInfo* pLnGrA;
        BoundedNumberVectorInfo* pLnGrB;
        BoundedNumberVectorInfo* pLnGrBeta;
        
        GrowthInfo();
        ~GrowthInfo();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/*------------------------------------------------------------------------------
 * MaturityInfo\n
 * Encapsulates the following molt-to-maturity-related parameters:\n
 *  pvLgtM2M: parameter vectors for logit-scale pr(molt-to-maturity|pre-molt size)
 * Notes:
 *  1. YEAR_BLOCK is the 1st index variable for the parameters
 *  1. SEX        is the 2nd index variable for the parameters
*----------------------------------------------------------------------------*/
class Molt2MaturityInfo: public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"molt_to_maturity"
    public:        
        BoundedVectorVectorInfo* pLgtPrM2M; //parameter vectors for logit-scale pr(molt-to-maturity|size)
        
        Molt2MaturityInfo();
        ~Molt2MaturityInfo();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
        
        friend cifstream& operator >>(cifstream & is, Molt2MaturityInfo & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, Molt2MaturityInfo & obj){obj.write(os); return os;}
};

/*------------------------------------------------------------------------------
 * SelectivityInfo\n
 * Encapsulates the following selectivity-related parameters:\n
 *   pS1  : 1st input to selectivity function
 *   pS2  : 2nd input to selectivity function
 *   pS3  : 3rd input to selectivity function
 *   pS4  : 4th input to selectivity function
 *   pS5  : 5th input to selectivity function
 *   pS6  : 6th input to selectivity function
 *   pDevsS1  : devs to 1st input to selectivity function
 *   pDevsS2  : devs to 2nd input to selectivity function
 *   pDevsS3  : devs to 3rd input to selectivity function
 *   pDevsS4  : devs to 4th input to selectivity function
 *   pDevsS5  : devs to 5th input to selectivity function
 *   pDevsS6  : devs to 6th input to selectivity function
 * Notes:
 *  1. no index variables for the parameters
*----------------------------------------------------------------------------*/
class SelectivityInfo : public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"selectivities"
    public:
        BoundedNumberVectorInfo* pS1;
        BoundedNumberVectorInfo* pS2;
        BoundedNumberVectorInfo* pS3;
        BoundedNumberVectorInfo* pS4;
        BoundedNumberVectorInfo* pS5;
        BoundedNumberVectorInfo* pS6;
        DevsVectorVectorInfo* pDevsS1;
        DevsVectorVectorInfo* pDevsS2;
        DevsVectorVectorInfo* pDevsS3;
        DevsVectorVectorInfo* pDevsS4;
        DevsVectorVectorInfo* pDevsS5;
        DevsVectorVectorInfo* pDevsS6;
        
        SelectivityInfo();
        ~SelectivityInfo();
        
        /**
         * Returns the values of the "extra variables" as doubles for
         * the pc-th parameter combination 
         * @param pc
         * @return 
         */
        dvector getPCXDs(int pc){return xd(pc);}
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/*------------------------------------------------------------------------------
 * FisheriesInfo\n
 * Encapsulates the following fishery-related parameters:\n
 *   pHM    : handling mortality (0-1)
 *   pLnC   : ln-scale base mean capture rate (mature males)
 *   pLnCT  : main year_block ln-scale offsets
 *   pLnDCX : ln-scale female offsets
 *   pLnDCM : ln-scale immature offsets
 *   pLnDCXM: ln-scale female-immature offsets    
 *   pLnEffX: ln-scale effort extrapolation parameters
 *
 *   pDevsLnC : annual ln-scale devs w/in year_blocks
 * 
 * Notes:
 *  1. index variables for parameters:
 *      a. SURVEY
 *      b. YEAR_BLOCK
 *      c. SEX
 *      d. MATURITY
 *      e. SHELL
*----------------------------------------------------------------------------*/
class FisheriesInfo : public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"fisheries"
    public:
        int idxUseEX;//column in parameter combinations matrix w/ useEffortRatio flag
        BoundedNumberVectorInfo* pHM;    //handling mortality (0-1)
        BoundedNumberVectorInfo* pLnC;   //ln-scale base mean capture rate (mature males)
        BoundedNumberVectorInfo* pLnDCT; //main year_block ln-scale offsets
        BoundedNumberVectorInfo* pLnDCX; //ln-scale female offsets
        BoundedNumberVectorInfo* pLnDCM; //ln-scale immature offsets
        BoundedNumberVectorInfo* pLnDCXM;//ln-scale female-immature offsets 
        BoundedNumberVectorInfo* pLnEffX;//ln-scale effort extrapolation 
        
        DevsVectorVectorInfo* pDevsLnC;//annual ln-scale devs w/in year_blocks
        
        FisheriesInfo();
        ~FisheriesInfo();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/*------------------------------------------------------------------------------
 * SurveysInfo\n
 * Encapsulates the following recruitment-related parameters:\n
 *   pLnQ   : base q (mature males)
 *   pLnQT  : main temporal offset
 *   pLnDQX : female offsets
 *   pLnDQM : immature offsets
 *   pLnDQXM: female-immature offsets    
 * Notes:
 *  1. index variables for parameters:
 *      a. SURVEY
 *      b. YEAR_BLOCK
 *      c. SEX
 *      d. MATURITY
 *      e. SHELL
*----------------------------------------------------------------------------*/
class SurveysInfo : public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"surveys"
    public:
        BoundedNumberVectorInfo* pLnQ;
        BoundedNumberVectorInfo* pLnDQT;
        BoundedNumberVectorInfo* pLnDQX;
        BoundedNumberVectorInfo* pLnDQM;
        BoundedNumberVectorInfo* pLnDQXM;
        
        SurveysInfo();
        ~SurveysInfo();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/*------------------------------------------------------------------------------
 * ModelParametersInfo
 *----------------------------------------------------------------------------*/
class ModelParametersInfo{
    public:
        static int debug;        
    public:
        ModelConfiguration* ptrMC;
        
        RecruitmentInfo*      ptrRec; //pointer to recruitment info
        NaturalMortalityInfo* ptrNM;  //pointer to natural mortality info
        GrowthInfo*           ptrGrw;  //pointer to growth info
        Molt2MaturityInfo*    ptrM2M; //pointer to molt-to-maturity info
        
        SelectivityInfo*      ptrSel; //pointer to selectivity functions info
        FisheriesInfo*        ptrFsh; //pointer to fisheries info
        SurveysInfo*          ptrSrv; //pointer to surveys info
    public:
        ModelParametersInfo(ModelConfiguration& mc);
        ~ModelParametersInfo();
        
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, ModelParametersInfo & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, ModelParametersInfo & obj){obj.write(os); return os;}
};

namespace tcsam{
    void readError(cifstream & is, const char * expP, adstring gotP);
    
    void setParameterInfo(NumberVectorInfo* pNVI,                           
                          int& npT,
                          ivector& phs, 
                          ostream& os = std::cout);
    void setParameterInfo(BoundedNumberVectorInfo* pBNVI,
                          int& npT,
                          dvector& lb, dvector& ub, 
                          ivector& phs, 
                          ostream& os = std::cout);
    void setParameterInfo(VectorVectorInfo* pVVI,   
                          int& npT,
                          ivector& mns, ivector& mxs,
                          ivector& phs, 
                          ostream& os = std::cout);
    void setParameterInfo(BoundedVectorVectorInfo* pBVVI,                           
                          int& npT,
                          ivector& mns, ivector& mxs,
                          imatrix& idxs,
                          dvector& lb, dvector& ub,
                          ivector& phs,
                          ostream& os = std::cout);
    void setParameterInfo(DevsVectorVectorInfo* pDVVI,                           
                          int& npT,
                          ivector& mns, ivector& mxs,
                          imatrix& idxs,
                          dvector& lb, dvector& ub,
                          ivector& phs,
                          ostream& os = std::cout);
    void writeParameter(ofstream& os, param_init_number& p, int toR, int willBeActive);
    void writeParameter(ofstream& os, param_init_bounded_number& p,int toR, int willBeActive);
    void writeParameter(ofstream& os, param_init_vector& p, int toR, int willBeActive);
    void writeParameter(ofstream& os, param_init_bounded_vector& p, int toR, int willBeActive);
    void writeParameterBounds(ofstream& os, param_init_bounded_dev_vector& p, int toR, int willBeActive);
}
#endif	/* MODELPARAMETERSINFO_HPP */

