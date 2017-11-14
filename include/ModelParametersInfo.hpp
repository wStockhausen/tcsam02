/**
 * @file
 * Header definitions for model classes associated with model parameters
 * information.
 * 
 * Includes:
 * <ul>
 *  <li> ParameterGroupInfo
 *  <li> RecruitmentInfo
 *  <li> NaturalMortalityInfo
 *  <li> GrowthInfo
 *  <li> Molt2MaturityInfo
 *  <li> SelectivityInfo
 *  <li> FisheriesInfo
 *  <li> SurveysInfo
 *  <li> ModelParametersInfo
 * </ul>
 */
#ifndef MODELPARAMETERSINFO_HPP
#define	MODELPARAMETERSINFO_HPP

#include <admodel.h>
#include "ModelConfiguration.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelParameterInfoTypes.hpp"

/*------------------------------------------------------------------------------\n
 * ParameterGroupInfo\n
 * This is the abstract base class for all concrete ..ParametersInfo classes,
 * except ModelParametersInfo.
 *----------------------------------------------------------------------------*/
class ParameterGroupInfo{
    public:
        /** flag to print debugginh info */
        static int debug;
    public:
        /** name of the ParameterGroup */
        adstring name;
        /** number of index variables */
        int nIVs;         
        /** names (labels) for the index variables */
        adstring_array lblIVs;
        /** number of parameter variables */
        int nPVs; 
        /** names (labels) for the parameter variables */
        adstring_array lblPVs;
        /** descriptions for the parameter variables */
        adstring_array dscPVs;
        /** number of "extra" indices or values */
        int nXIs;             //number of extra indices/values
        /** names (labels) for the extra indices/variables */
        adstring_array lblXIs;
        
        /** number of index variables (IVs) that are blocks */
        int nIBSs;
        /** indices corresponding to index variables (IVs) that are blocks */
        ivector ibsIdxs;
        /** pointer to a vector of pointers to IndexBlockSet objects */
        IndexBlockSet** ppIBSs;
        
        /** number of rows in the parameter combinations matrix */
        int nPCs;
        /** input parameter combinations matrix (all integers) */
        imatrix in;
        /** "extra" values by parameter combination as doubles */
        dmatrix xd;
        /** pointer to array of pointers to indices matrices for use via getModelIndices(pc) */
        imatrix** ppIdxs;
        /** array of labels, one for each pc */
        adstring_array pcLabels;
    public:
        /**
         * Constructor for the class.
         */
        ParameterGroupInfo();
        /** 
         * Destructor for the class 
         */
        ~ParameterGroupInfo();
        
        /* 
         * Returns a pointer to the index block set identified by "type".
         * 
         * @param type - adstring "type" identifying index block set to return
         * 
         * @return pointer to the identified IndexBlockSet
         */
        IndexBlockSet* getIndexBlockSet(adstring type);
        
        /**
         * Gets the indices for the parameter combination.
         * 
         * @param pc - the id for the desired parameter combination
         * 
         * @return an ivector of the parameter indices specifying the parameter combination
         */
        ivector getPCIDs(int pc);
        /**
         * Gets the model indices associated with the parameter combination.
         * 
         * @param pc - the id for the desired parameter combination
         * 
         * @return an imatrix of model indices associated with the parameter combination
         */
        imatrix getModelIndices(int pc);
        
        /**
         * Returns the values of the "extra variables" as doubles for
         * the pc-th parameter combination 
         * 
         * @param pc - the index of the desired parameter combination
         * 
         * @return dvector of the values of the "extra" variables associated with the pc
         */
        dvector getPCXDs(int pc);
                
        /**
         * Reads default info in ADMB format for the parameter group from an input filestream.
         * Subclasseses should override this function as appropriate.
         * 
         * @param is - the input filestream
         */
        virtual void read(cifstream & is);
        /**
         * Writes default info in ADMB format for the parameter group to an output stream.
         * Subclasses should override this function as appropriate.
         * 
         * @param os - the output stream
         */
        virtual void write(std::ostream & os);
        /**
         * Writes default info in R format for the parameter group to an output stream.
         * Subclasses should override this function as appropriate.
         * 
         * @param os - the output stream
         */
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

/**
 * @class RecruitmentInfo:ParameterGroupInfo
 * 
 * This class encapsulates the ParameterGroupInfo for
 * the following recruitment-related parameters:
 * <ul>
 *  <li> pLnR   : mean ln-scale recruitment
 *  <li> pRCV   : recruitment cv's
 *  <li> pRX    : fraction males at recruitment to population
 *  <li> pRa    : size-at-recruitment parameter
 *  <li> pRb    : size-at-recruitment parameter
 *  <li> pvLnRDevs: ln-scale annual recruitment devs
 * </ul>
 * Notes:
 * <ol>
 *  <li> YEAR_BLOCK is the index variable for the parameters
 * </ol>
 */
class RecruitmentInfo: public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"recruitment"
    public:        
        /** pointer to a vector of bounded parameters for the ln-scale mean recruitment */
        BoundedNumberVectorInfo* pLnR;
        /** pointer to a vector of bounded parameters for the recruitment cv */
        BoundedNumberVectorInfo* pRCV;
        /** pointer to a vector of bounded parameters for the fraction of males at recruitment */
        BoundedNumberVectorInfo* pRX;
        /** pointer to a vector of bounded parameters for the recruitment distribution scale parameter */
        BoundedNumberVectorInfo* pRa;
        /** pointer to a vector of bounded parameters for the recruitment distribution shape parameter */
        BoundedNumberVectorInfo* pRb;
        /** pointer to info for a vector of parameter vectors for annual recruitment devs */
        DevsVectorVectorInfo* pDevsLnR; 
        
        /**
         * Class constructor.
         */
        RecruitmentInfo();
        /**
         * Class destructor.
         */
        ~RecruitmentInfo();
        
        /**
         * Reads the ParameterGroupInfo for recruitment from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
        void read(cifstream & is);
        /**
         * Writes the ParameterGroupInfo for recruitment to an output stream in ADMB format.
         * 
         * @param is - the output stream
         */
        void write(std::ostream & os);
        /**
         * Writes the ParameterGroupInfo for recruitment to an output stream in R format.
         * 
         * @param is - the output stream
         */
        void writeToR(std::ostream & os);
        
        friend cifstream& operator >>(cifstream & is, RecruitmentInfo & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, RecruitmentInfo & obj){obj.write(os); return os;}
};

/**
 * @class NaturalMortalityInfo:ParameterGroupInfo
 * 
 * This class encapsulates the following natural mortality-related parameters:
 * <ul>
 *   <li> pM   : base arithmetic-scale natural mortality
 *   <li> pDM1 : ln-scale offsets
 *   <li> pDM2 : ln-scale offsets
 *   <li> pDM3 : ln-scale offsets
 *   <li> pDM4 : ln-scale offsets    
 * </ul>
 * Notes:
 * <ol type="1">
 *  <li> index variables for parameters
 *    <ol type="a">
 *      <li> YEAR_BLOCK
 *      <li> SEX
 *      <li> MATURITY
 *      <li> SHELL
 *    </ol>
 *  </li>
 * </ol>
*/
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
        BoundedNumberVectorInfo* pM;
        /* info for offset 1 */
        BoundedNumberVectorInfo* pDM1;
        /* info for offset 2 */
        BoundedNumberVectorInfo* pDM2;
        /* info for offset 3 */
        BoundedNumberVectorInfo* pDM3;
        /* info for offset 4 */
        BoundedNumberVectorInfo* pDM4;
        
        /**
         * Class constructor.
         */
        NaturalMortalityInfo();
        /**
         * Class destructor
         */
        ~NaturalMortalityInfo();
        
        /**
         * Reads the ParameterGroupInfo for natural mortality from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
        void read(cifstream & is);
        /**
         * Writes to an output stream in ADMB format.
         * 
         * @param os - output stream
         */
        void write(std::ostream & os);
        /**
         * Writes component info to an output stream as an
         * R-format list.
         * 
         * @param os - output stream
         */
        void writeToR(std::ostream & os);
};

/**
 * @class GrowthInfo:ParameterGroupInfo
 * 
 * This class encapsulates the following recruitment-related parameters:
 * <ul>
 *   <li> pGrA : mean growth coefficient "a"
 *   <li> pGrB : mean growth exponent "b"
 *   <li> pGrBeta : scale factor for growth transition matrix
 * </ul>
 * Notes:
 * <ol>
 *  <li> the index variables for the parameters are
 *      <ol type="a">
 *          <li> YEAR_BLOCK
 *          <li> SEX
 *      </ol>
 *  </li>
 * </ol>
 */
class GrowthInfo : public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        /** name for this parameter group */
        static adstring NAME;
    public:
        /** pointer to info on mean growth coefficient "a" parameters */
        BoundedNumberVectorInfo* pGrA;
        /** pointer to info on mean growth coefficient "b" parameters */
        BoundedNumberVectorInfo* pGrB;
        /** pointer to info on the scale factor for the growth transition matrix */
        BoundedNumberVectorInfo* pGrBeta;
        
        /**
         * Class constructor.
         */
        GrowthInfo();
        /**
         * Class destructor.
         */
        ~GrowthInfo();
        
        /**
         * Reads the ParameterGroupInfo for growth from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/**
 * @class MaturityInfo:ParameterGroupInfo
 * 
 * This class encapsulates the following molt-to-maturity-related parameters:
 * <ul>
 *  <li> pvLgtM2M - parameter vectors for logit-scale pr(molt-to-maturity|pre-molt size)
 * </ul>
 * <ol>
 *  <li> the index variables for the parameters are
 *      <ol type="a">
 *          <li> YEAR_BLOCK
 *          <li> SEX
 *      </ol>
 *  </li>
 * </ol>
*/
class Molt2MaturityInfo: public ParameterGroupInfo {
    public:
        static int debug;
    protected:
        static adstring NAME;//"molt_to_maturity"
    public:        
        /** pointer to info for parameter vectors describing logit-scale pr(molt-to-maturity|size) */
        BoundedVectorVectorInfo* pLgtPrM2M; 
        
        /**
         * Class constructor.
         */
        Molt2MaturityInfo();
        /**
         * Class destructor.
         */
        ~Molt2MaturityInfo();
        
        /**
         * Reads the ParameterGroupInfo for the molt-to-maturity from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
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
         * Reads the ParameterGroupInfo for selectivity from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/*------------------------------------------------------------------------------
 * FisheriesInfo\n
 * Encapsulates the following fishery-related parameters:\n
 *   pHM    : handling mortality (0-1)
 *   pLnC   : ln-scale base mean capture rate (mature males)
 *   pDC1   : ln-scale offset 1 (e.g., main year_block ln-scale offsets)
 *   pDC2   : ln-scale offset 2 (e.g., female offsets)
 *   pDC3   : ln-scale offset 3 (e.g., immature offsets)
 *   pDC4   : ln-scale offset 4 (e.g., female-immature offsets)   
 *
 *   pDevsLnC : annual ln-scale devs w/in year_blocks
 * 
 *   pLnEffX: ln-scale effort extrapolation parameters
 *   pLgtRet: logit-scale max retention parameter
 * 
 * Notes:
 *  1. index variables for parameters:
 *      a. FISHERY
 *      b. YEAR_BLOCK
 *      c. SEX
 *      d. MATURITY
 *      e. SHELL
*----------------------------------------------------------------------------*/
class FisheriesInfo : public ParameterGroupInfo {
    public:
        static int debug;
        static int idxHM;  //column in parameter combinations matrix with parameter index for column in parameter combinations matrix indicating handling mortality parameters
        static int idxLnC; //column in parameter combinations matrix with parameter index for ln-scale base mean capture rate (mature males)
        static int idxDC1; //column in parameter combinations matrix with parameter index for ln-scale offsets pDC1
        static int idxDC2; //column in parameter combinations matrix with parameter index for ln-scale offsets pDC2
        static int idxDC3; //column in parameter combinations matrix with parameter index for ln-scale offsets pDC3
        static int idxDC4; //column in parameter combinations matrix with parameter index for ln-scale offsets pDC4
        static int idxLnDevs;//column in parameter combinations matrix with parameter index for annual ln-scale devs w/in year_blocks
        static int idxLnEffX;//column in parameter combinations matrix with parameter index for ln-scale effort extrapolation 
        static int idxLgtRet;//column in parameter combinations matrix with parameter index for logit-scale retained fraction (for old shell crab)
        static int idxSelFcn;//column in parameter combinations matrix indicating selectivity function index
        static int idxRetFcn;//column in parameter combinations matrix indicating retention function index
        static int idxUseEX; //column in parameter combinations matrix indicating effort extrapolation use
    protected:
        static adstring NAME;//"fisheries"
    public:
        BoundedNumberVectorInfo* pHM;  //handling mortality (0-1)
        BoundedNumberVectorInfo* pLnC; //ln-scale base mean capture rate (mature males)
        BoundedNumberVectorInfo* pDC1; //ln-scale offsets
        BoundedNumberVectorInfo* pDC2; //ln-scale offsets
        BoundedNumberVectorInfo* pDC3; //ln-scale offsets
        BoundedNumberVectorInfo* pDC4; //ln-scale offsets 
        BoundedNumberVectorInfo* pLnEffX;//ln-scale effort extrapolation 
        BoundedNumberVectorInfo* pLgtRet;//logit-scale retained fraction (for old shell crab)
        
        DevsVectorVectorInfo* pDevsLnC;//annual ln-scale devs w/in year_blocks
        
        FisheriesInfo();
        ~FisheriesInfo();
        
        /**
         * Reads the ParameterGroupInfo for the fisheries from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/*------------------------------------------------------------------------------
 * SurveysInfo\n
 * Encapsulates the following recruitment-related parameters:\n
 *   pQ : base q (mature males)
 *   pDQ1 : ln-scale offset 1 (e.g., main temporal offset)
 *   pDQ2 : ln-scale offset 1 (e.g., female offsets)
 *   pDQ3 : ln-scale offset 1 (e.g., immature offsets)
 *   pDQ4 : ln-scale offset 1 (e.g., female-immature offsets)
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
        /** base Q */
        BoundedNumberVectorInfo* pQ;
        /** ln-scale offset 1 */
        BoundedNumberVectorInfo* pDQ1;
        /** ln-scale offset 2 */
        BoundedNumberVectorInfo* pDQ2;
        /** ln-scale offset 3 */
        BoundedNumberVectorInfo* pDQ3;
        /** ln-scale offset 4 */
        BoundedNumberVectorInfo* pDQ4;
        
        SurveysInfo();
        ~SurveysInfo();
        
        /**
         * Reads the ParameterGroupInfo for the surveys from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);
};

/**
 * @class ModelParametersInfo
 * 
 * This class encapsulates information for all the model parameters by parameter group and
 * provides read/write methods for ModelParametersInfo files.
 * 
 * The following "parameter groups" are incorporated here:
 * <ul>
 *  <li> recruitment       (ptrRec)
 *  <li> natural mortality (ptrNM)
 *  <li> growth            (ptrGrw)
 *  <li> molt-to-maturity  (ptrM2M)
 *  <li> selectivity       (ptrSel)
 *  <li> fisheries         (ptrFsh)
 *  <li> surveys           (ptrSrv)
 * </ul>
 */
class ModelParametersInfo{
    public:
        /** a flag to print debugging info */
        static int debug;      
        /** an adstring with the ModelParametersInfo version */
        static const adstring version;
    public:
        /** pointer to the ModelConfiguration object */
        ModelConfiguration* ptrMC;
        
        RecruitmentInfo*      ptrRec; //pointer to recruitment info
        NaturalMortalityInfo* ptrNM;  //pointer to natural mortality info
        GrowthInfo*           ptrGrw; //pointer to growth info
        Molt2MaturityInfo*    ptrM2M; //pointer to molt-to-maturity info
        
        SelectivityInfo*      ptrSel; //pointer to selectivity functions info
        FisheriesInfo*        ptrFsh; //pointer to fisheries info
        SurveysInfo*          ptrSrv; //pointer to surveys info
    public:
        /**
         * Class constructor.
         * 
         * @param mc - a reference to the ModelConsfiguration object
         */
        ModelParametersInfo(ModelConfiguration& mc);
        /**
         * Class destructor.
         */
        ~ModelParametersInfo();
        
        /**
         * Reads the model parameters info for all parameter groups from an input filestream in ADMB format.
         * 
         * @param is - the input filestream
         */
        void read(cifstream & is);
        void write(std::ostream & os);
        void writeToR(std::ostream & os);

        friend cifstream& operator >>(cifstream & is, ModelParametersInfo & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os, ModelParametersInfo & obj){obj.write(os); return os;}
};

#endif	/* MODELPARAMETERSINFO_HPP */

