/* 
 * File:   ModelData.hpp
 * Author: william.stockhausen
 *
 * Created on 2014-02-11.
 * 
 * 20140930: 1. Added recLag to BioData.
 * 20150209: 1. Removed prMature_xz, frMature_xsz, cvMnMxZ_xc from BioData.
 * 20150121: 1. Combined FisheryData, SurveyData classes into a FleetData class
 *           2. Updated some documentation.
 * 20161109: 1. added likelihood multiplier (llWgt) to Data classes
 *           2. added likelihood type to EffortData class
 *           3. changed effort averaging interval from IndexRange to IndexBlock
 * 20161109: 1. added GrowthData, ChelaHeightData classes
 *           2. revised ModelDatasets to incorporate new data classes
 * 20190530: 1. added MaturityOgiveData class
 *           2. revised ModelDatasets to incorporate new data class
 * 20201127: 1. Completed adding tail compression for size composition data.
 *           2. Added index to the Dirichlet-Multinomial parameter used for
 *                size composition data.
*/

#ifndef MODELDATA_HPP
#define MODELDATA_HPP

//**********************************************************************
//  Includes
//      AggregateCatchData
//      EffortData
//      SizeFrequencyData
//      CatchData
//      BioData
//      FleetData
//      GrowthData
//      ChelaHeightData
//      MaturtyOgiveData
//      ModelDatasets
//**********************************************************************
class ModelConfiguration; //forward definition
class IndexBlock;

//--------------------------------------------------------------------------------
//          AggregateCatchData
// Change notes:
//  2014-12-08: 1. changed input data row format for so columns are now
//                  year, sex, maturity, shell_condition, value, cv
//              2. changed inpC_yc dimensions from yc to yxmsc to match new input format
//              3. changed C_xy, cv_xy, sd_xy dimensions from xy to xmsy (and renamed them accordingly)
//
//--------------------------------------------------------------------------------
    /**
     * Class encapsulating aggregate (numbers, biomass) catch data from a data file.
     */
    class AggregateCatchData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static ostream& os;
        /* keyword indicating abundance (numbers) data */
        const static adstring KW_ABUNDANCE_DATA;
        /* keyword indicating biomass (weight) data */
        const static adstring KW_BIOMASS_DATA;  
        /* index into "c" for use flags in inpC_xmsyc */
        const static int idUF = 1;
        /* index into "c" for years in inpC_xmsyc */
        const static int idYr = 2;
        /* index into "c" for estimated values in inpC_xmsyc */
        const static int idEV = 3;
        /* index into "c" for cv's in inpC_xmsyc */
        const static int idCV = 4;
    protected:
        /* input aggregate catch data (c: use flag, value, cv by xmsy) */
        d5_array inpC_xmsyc; 
    public:
        /* type (abundance, biomass) of data */
        adstring type;  
        /* fleet name */
        adstring name;  
        /* objective function fitting option */
        int optFit;     
        /* likelihood function type */
        int llType; 
        /* likelihood weight (i.e., multiplier) */
        double llWgt;   
        int ny;         //number of data rows for aggregate catch
        adstring units; //units for aggregate catch data
        ivector yrs;    //years for aggregate catch data
        wts::adstring_matrix factors;//factor combinations for input numbers-at-size
        d4_array C_xmsy;   //aggregate catch by sex, maturity, shell_condition, year (converted from units to THOUSANDS of crab or MT))
        d4_array cv_xmsy;  //aggregate catch cv's by sex, maturity, shell_condition, year
        d4_array sd_xmsy;  //aggregate catch stdv's by sex, maturity, shell_condition, year
        d4_array uf_xmsy;  //use flags, by sex, maturity, shell_condition, year
        
    public:
        AggregateCatchData(adstring& name_){name = name_;}
        ~AggregateCatchData(){}
        /**
         * Set the maximum year in which to fit the data.
         * 
         * @param mxYr - the max year to include data
         */
        void setMaxYear(int mxYr);
        /**
         * Replace catch data C_xmsy with new data. 
         * Also modifies inpC_xmsyc to reflect new data.
         * Error-related quantities remain the same.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newC_xmsy - d4_array with new catch data
         */
        void replaceCatchData(int iSeed,random_number_generator& rng,d4_array& newC_xmsy);
        /**
         * Update catch data C_xmsy with data for new year 
         * Also modifies inpC_xmsyc to reflect new data.
         * 
         * @param y - year to add
         * @param newC_xms - d3_array with new catch data
         */
        void addCatchData(int y ,d3_array& newC_xms);
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        friend cifstream& operator >>(cifstream & is, AggregateCatchData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   AggregateCatchData & obj){obj.write(os); return os;}
    
    protected:
        /**
         * Aggregate catch data C_xmsy, cv_xmsy, sd_xmsy over summary indices.
         */
        void aggregateData(void);

};

//--------------------------------------------------------------------------------
//          EffortData
//--------------------------------------------------------------------------------
    /**
     * Class encapsulating effort data from a data file.
     */
    class EffortData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static ostream& os;
        /* keyword indicating effort data */
        const static adstring KW_EFFORT_DATA;
    public:
        /** fleet name */
        adstring name;
        /* likelihood function type */
        int llType; 
        /* likelihood weight (i.e., multiplier) */
        double llWgt;   
        /* pointer to intervals (IndexBlock) over which to average effort/fishing mortality */
        IndexBlock* ptrAvgIB; 
        /* units for potlifts */
        adstring units;   
        /* number of years of effort data */
        int ny;
        /* input effort data (year)x(year,potlifts) */
        dmatrix  inpEff_yc;   
        /* vector of years w/ effort data */
        dvector yrs;
        /* effort data vector corresponding to \code{yrs} */
        dvector eff_y;        
    public:
        /**
         * Constructor.
         */
        EffortData(adstring& name_){name=name_; ptrAvgIB=0;}
        /**
         * Destructor.
         */
        ~EffortData();
        /**
         * Set the maximum year in which to fit the data.
         * 
         * @param mxYr - the max year to include data
         */
        void setMaxYear(int mxYr);
        /**
         * Read input data in ADMB format from a file stream
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write data to an output stream in ADMB format
         * 
         * @param os output stream
         */
        void write(std::ostream & os); //write object to file in ADMB format
        /**
         * Write data to an output stream as an R-formatted list object
         * 
         * @param os output stream
         */
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into an EffortData object.
         */
        friend cifstream& operator >>(cifstream & is, EffortData & obj){obj.read(is); return is;}
        /**
         * Operators to write data to an output stream in ADMB format from an EffortData object.
         */
        friend std::ostream&   operator <<(std::ostream & os,   EffortData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          SizeFrequencyData
//--------------------------------------------------------------------------------
    /**
     * Class encapsulating size frequency data from a data file.
     */
    class SizeFrequencyData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static std::ostream& os;
        /* keyword indicating size frequency data */
        const static adstring KW_SIZEFREQUENCY_DATA;
        /* integer indicating last column in inpNatZ_xmsyc BEFORE nAtZ data */
        const static int LAST_COL=4;//use flag, dm index, year, ss
    private:
        /* factor combinations for input numbers-at-size */
        wts::adstring_matrix factors;
        /* input numbers-at-size data (sex,maturity state,shell condition,year,use+idxDM+year+sample_size+nAtZ) */
        d5_array inpNatZ_xmsyc;       
    public:
        /** fleet name */
        adstring name;
        /* objective function fitting option */
        int optFit; 
        /* likelihood function type */
        int llType; 
        /* likelihood weight (i.e., multiplier) */
        double llWgt;   
        /* number of size bin cut pts */
        int nZCs;    
        /* vector of cut points for size bins (1:nZCs) */
        dvector zCs; 
        /* vector of size bin centers (1:nZCs-1) */
        dvector zBs; 
        
        /* units for numbers-at-size data */
        adstring units; 
        
        /* number of years of size frequency data */
        int ny;         
        /* years of size frequency data */
        ivector  yrs;       
        /* use flags for size frequency data */
        d4_array inpUF_xmsy;   
        /* input indices for Dirichlet-Multinomial parameter */
        d4_array inpDM_xmsy;
        /* input sample sizes for size frequency data */
        d4_array inpSS_xmsy;   
        /* raw size frequency data */
        d5_array rawNatZ_xmsyz;
        /* size frequency data aggregated to fit option */
        d5_array aggNatZ_xmsyz;
        /* tail-compressed aggregated size frequency data */
        d5_array tcdNatZ_xmsyz;
        /* normalized aggregated size frequency data (sums to 1 over xmsz for each y) */
        d5_array aggPatZ_xmsyz;
        
        /* tail compression limits (min, 1-max) */
        dvector tc_limits;
        /* tail compression indices */
        i5_array tc_xmsyc;   
        
        
        /* cumulative iterative re-weighting factors */
        d3_array cumF_xms;
        /* last set of iterative re-weighting factors */
        d3_array itrF_xms;
        
        /* working use flags for size frequency data */
        i4_array uf_xmsy;   
        /* working Dirichlet-Multinomial parameter indices for size frequency data */
        i4_array dm_xmsy;   
        /* working (possibly re-weighted) sample sizes for size frequency data */
        d4_array ss_xmsy;   
        
    public:
        /**
         * Constructor.
         */
        SizeFrequencyData(adstring& name_){name=name_;}
        /**
         * Destructor.
         */
        ~SizeFrequencyData(){}
        /**
         * Set the maximum year in which to fit the data.
         * 
         * @param mxYr - the max year to include data
         */
        void setMaxYear(int mxYr);
        /**
         * Replace catch-at-size data NatZ_xmsyz with new data. 
         * Also modifies inpNatZ_xmsyc to reflect new data.
         * Error-related quantities remain the same.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newNatZ_yxmsz - d5_array with new numbers-at-size data
         */
        void replaceSizeFrequencyData(int iSeed,random_number_generator& rng,d5_array& newNatZ_xmsyz);
        /**
         * Update catch-at-size data NatZ_xmsyz with new data for year y. 
         * Also modifies inpNatZ_xmsyc to reflect new data.
         * 
         * @param y - year
         * @param newNatZ_xmsz - d4_array with new numbers-at-size data
         */
        void addSizeFrequencyData(int y, d4_array& newNatZ_xmsz);
        /**
         * Save the negative log-likelihoods from a model fit (values only).
         * 
         * @param nlls
         */
        void saveNLLs(dvar_matrix& nlls);
        /**
         * Calculate re-weighting factors for iterative weighting.
         * 
         * @param newF_xms - d3_array by xms with new (incremental) weighting factors
         * @param debug - integer flag to print debugging info
         * @cout - output stream to print debugging info to
         */
        void calcReWeightingFactors(d3_array& newF_xms,int debug,ostream& cout);
        /**
         * Apply re-weighting factors for iterative weighting to 
         * input sample sizes.
         * 
         */
        void applyReWeightingFactors();
        /**
         * Read input data in ADMB format from a file stream
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write data to an output stream in ADMB format
         * 
         * @param os output stream
         */
        void write(std::ostream & os); //write object to file in ADMB format
        /**
         * Write data to an output stream as an R-formatted list object
         * 
         * @param os output stream
         */
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into an SizeFrequencyData object.
         */
        friend cifstream& operator >>(cifstream & is, SizeFrequencyData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   SizeFrequencyData & obj){obj.write(os); return os;}
    public:
        /**
         * Calculate normalized size compositions PatZ_xmsyz based on NatZ_xmsyz.
         * 
         * Normalization is by year across xmsz.
         */
        void normalize(void);
        
        /**
         * Create the indices used for tail compression based on aggregated size comps.
         * 
         * Tail compression is specified for aggregated size comps by xmsy.
         */
        void doTailCompression(void);
        
        /**
         * Calculate aggregated size compositions from raw size comps prior to tail compression.
         * 
         * Aggregation is done according to value of optFit.
         * 
         * Calculates aggNatZ_xmsyz from rawNatZ_xmsyz.
         * 
         */
        void aggregateRawNatZ(void);
    private:
        /**
         * Set the DM parameter index based on previous value and current value.
         * 
         * @param x
         * @param m
         * @param s
         * @param y - year
         * @param dm - previous value of DM parameter index
         * @param inpDM - current value of DM parameter index
         * 
         * @return possibly updated value for dm.
         */
        int setDM(int x,int m,int s,int y,int dm,int inpDM);
    };

//--------------------------------------------------------------------------------
//          BioData
//--------------------------------------------------------------------------------
    /**
     * Class encapsulating biological information (weight-at-age, etc.)
     */
    class BioData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static std::ostream& os;
        const static adstring KW_BIO_DATA;
    public:
        int nZBins;           //number of size bin cut pts
        dvector zBins;        //size bins
        adstring unitsWatZ;   //units for weight-at-size
        d3_array wAtZ_xmz;    //weight at size (in kg)
        
        int recLag;         //recruitment lag (in years)
        double fshTimingTypical;//typical midpoint of fishery seasons
        double matTimingTypical;//typical timing of mating
        int nyAtypicalT;        //number of years for atypical fishery season midpoints, time-at-mating
        dvector fshTiming_y;    //timing of midpoint of fisheries seasons
        dvector matTiming_y;    //timing of mating
    protected:
        dmatrix timing_yc;//midpoint of fishery season, time-at-mating by year
    public:
        BioData(){}
        ~BioData(){}
        void read(cifstream & is);//read file in ADMB format
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a BioData object.
         */
        friend cifstream& operator >>(cifstream & is, BioData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   BioData & obj){obj.write(os); return os;}
    };

//------------------------------------------------------------------------------
    /**
     * CatchData class.
     */
    class CatchData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static std::ostream& os;
        const static adstring KW_CATCH_DATA;
    protected:
        dmatrix inpN_yc;       //input catch abundance data (year,female abundance,female cv,male abundance, male cv, total abundance, cv_total)
        dmatrix inpB_yc;       //input catch biomass   data (year,female biomass,female cv,male biomass, male cv, total biomass, cv_total)
        //TODO: need sample sizes (tows/potlifts, non-zero tows/potlifts, num individ.s by xms, calculated ss)
        d5_array inpNatZ_xmsyc;//input numbers-at-size data (sex,maturity,shell,year) x (year,sample_size,nAtZ)
    public:
        adstring type;    //data type (e.g. survey or fishery)
        adstring name;    //data source name
        int hasN;         //abundance data flag
        AggregateCatchData* ptrN;//pointer to aggregate abundance data
        int hasB;         //biomass data flag
        AggregateCatchData* ptrB;//pointer to aggregate biomass data
        int hasZFD;       //numbers-at-size data flag
        SizeFrequencyData* ptrZFD;//pointer to numbers-at-size data
        
    public:
        /**
         * Constructor.
         */
        CatchData(adstring& name_);
        /**
         * Destructor
         */
        ~CatchData();
        /**
         * Set the maximum year in which to fit the data.
         * 
         * @param mxYr - the max year to include data
         */
        void setMaxYear(int mxYr);
        /**
         * Replace catch data based on newNatZ_yxmsz.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newNatZ_yxmsz - d5_array of catch-at-size by sex/maturity/shell condition/year
         * @param wAtZ_xmz - d3_array of weight-at-size by sex/maturity
         */
        virtual void replaceCatchData(int iSeed,random_number_generator& rng,d5_array& newNatZ_yxmsz, d3_array& wAtZ_xmz);
        /**
         * Adds a new year of catch data based on dvar4_array newNatZ_xmsz to existing data.
         * 
         * @param y - year to add
         * @param newNatZ_xmsz - dvar4_array of catch numbers-at-size by sex/maturity/shell condition
         * @param wAtZ_xmz - weight-at-size by sex/maturity
         * @param cv - cv for aggregated catch sampling error
         * @param ss - sample size for size frequency sampling error
         * @param rng - random number generator
         * 
         * @return void
         */
        virtual void addCatchData(int y, 
                                  dvar4_array& newNatZ_xmsz, 
                                  d3_array& wAtZ_xmz, 
                                  double cv, 
                                  double ss,
                                  random_number_generator& rng);
        /**
         * Adds a new year of catch data based on d4_array newNatZ_xmsz to existing data.
         * 
         * @param y - year to add
         * @param newNatZ_xmsz - d4_array of catch numbers-at-size by sex/maturity/shell condition
         * @param wAtZ_xmz - weight-at-size by sex/maturity
         * @param cv - cv for aggregated catch sampling error
         * @param ss - sample size for size frequency sampling error
         * @param rng - random number generator
         * 
         * @return void
         */
        virtual void addCatchData(int y, 
                                  d4_array& newNatZ_xmsz, 
                                  d3_array& wAtZ_xmz, 
                                  double cv, 
                                  double ss,
                                  random_number_generator& rng);
        /**
         * Read a data file in ADMB format.
         * 
         * @param is - input filestream to read from
         */
        virtual void read(cifstream & is);//read file in ADMB format
        /**
         * Write a data file in ADMB format.
         * 
         * @param os
         * @param nm
         * @param indent
         */
        virtual void write(std::ostream & os); //write object to file in ADMB format
        /**
         * Write a data file in R format.
         * 
         * @param os
         * @param nm
         * @param indent
         */
        virtual void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a CatchData object.
         */
        friend cifstream& operator >>(cifstream & is, CatchData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   CatchData & obj){obj.write(os); return os;}
    };

//------------------------------------------------------------------------------    
    /**
     * Fleet data class.
     */
    class FleetData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static ostream& os;
        /** key word indicating fishery data */
        const static adstring KW_FISHERY;
        /** key word indicating survey data */
        const static adstring KW_SURVEY;
    public:
        /** fleet name */
        adstring name;     
        /** fleet type (fishery or survey) */
        adstring type;     
        /** flag indicating observed index (survey) catch data */
        int hasICD;        
        /**pointer to CatchData object for index survey) catch data */
        CatchData* ptrICD; 
        /** flag indicating retained catch data */
        int hasRCD;        
        /** pointer to CatchData object for retained catch */
        CatchData* ptrRCD; 
        /** flag indicating observed discard catch data */
        int hasDCD;        
        /** pointer to CatchData object for observed discards */
        CatchData* ptrDCD; 
        /** flag indicating observed total catch data */
        int hasTCD;        
        /** pointer to CatchData object for observed total catch */
        CatchData* ptrTCD; 
        /** flag indicating effort data */
        int hasEff;        
        /** pointer to effort data */
        EffortData* ptrEff;
    public:
        /**
         * Constructor.
         */
        FleetData();
        /**
         * Destructor.
         */
        ~FleetData();
        /**
         * Set the maximum year in which to fit the data.
         * 
         * @param mxYr - the max year to include data
         */
        void setMaxYear(int mxYr);
        /**
         * Replace existing index catch data with new values.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newNatZ_yxmsz - catch data array
         * @param wAtZ_xmz - weight-at-size array
         */
        void replaceIndexCatchData(int iSeed,random_number_generator& rng,d5_array& newNatZ_yxmsz, d3_array& wAtZ_xmz);
        /**
         * Replace existing fishery catch (retained, discarded, total) data with new values.
         * 
         * @param iSeed - flag (!=0) to add random noise
         * @param rng - random number generator
         * @param newCatZ_yxmsz - total catch data array
         * @param newRatZ_yxmsz - retained catch data array
         * @param wAtZ_xmz      - weight-at-size array
         */
        void replaceFisheryCatchData(int iSeed,random_number_generator& rng,d5_array& newCatZ_yxmsz,d5_array& newRatZ_yxmsz,d3_array& wAtZ_xmz);
        /**
         * Add a new year of index catch data to existing data.
         * 
         * @param y - year to add
         * @param newNatZ_xmsz - dvar4_array of catch data
         * @param wAtZ_xmz - weight-at-size array
         * @param cv - cv for aggregated catch sampling error
         * @param ss - sample size for size frequency sampling error
         * @param rng - random number generator
         * 
         * @return void
         */
        void addIndexCatchData(int y, 
                                dvar4_array& newNatZ_xmsz, 
                                d3_array& wAtZ_xmz, 
                                double cv, 
                                double ss,
                                random_number_generator& rng);
        /**
         * Add a new year of fishery catch (retained, discarded, total) data to existing.
         * 
         * @param y - year to add
         * @param newCatZ_xmsz - dvar4_array of total catch data
         * @param newRatZ_xmsz - dvar4_array of retained catch data
         * @param wAtZ_xmz - weight-at-size array
         * @param cv - cv for aggregated catch sampling error
         * @param ss - sample size for size frequency sampling error
         * @param rng - random number generator
         * 
         * @return void
         */
        void addFisheryCatchData(int y, 
                                dvar4_array& newCatZ_xmsz, 
                                dvar4_array& newRatZ_xmsz, 
                                d3_array& wAtZ_xmz,
                                double cv, 
                                double ss,
                                random_number_generator& rng);
        /**
         * Read input file in ADMB format from input stream.
         * 
         * @param is - input stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write the fleet data in ADMB format (i.e., as an input file)
         * 
         * @param os - output stream to write to
         */
        void write(std::ostream & os); //write object to file in ADMB format
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a FleetData object.
         */
        friend cifstream& operator >>(cifstream & is, FleetData & obj){obj.read(is); return is;}
        friend std::ostream&   operator <<(std::ostream & os,   FleetData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//          GrowthData
//--------------------------------------------------------------------------------
    /**
     * Class encapsulating growth data.
     */
    class GrowthData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static ostream& os;
        /** keyword indicating effort data */
        const static adstring KW_GROWTH_DATA;
    public:
        /** dataset name */
        adstring name;
        /** likelihood function type */
        int llType; 
        /** likelihood weight (i.e., multiplier) */
        double llWgt;   
        /** number of sex classes */
        int nSXs;
        /** number of observations, by sex */
        ivector nObs_x;
        /** years associated with observations, by sex */
        imatrix obsYears_xn;
        /** input data, by sex, c: year,pre-molt size, post-molt size, n: observations */
        d3_array inpData_xcn;  
    public:
        /**
         * Constructor.
         */
        GrowthData();
        /**
         * Destructor.
         */
        ~GrowthData();
        /**
         * Set the maximum year in which to fit the data.
         * 
         * @param mxYr - the max year to include data
         */
        void setMaxYear(int mxYr);
        /**
         * Read input data in ADMB format from a file stream
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write data to an output stream in ADMB format
         * 
         * @param os output stream
         */
        void write(std::ostream & os); //write object to file in ADMB format
        /**
         * Write data to an output stream as an R-formatted list object
         * 
         * @param os output stream
         */
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into an GrowthData object.
         */
        friend cifstream& operator >>(cifstream & is, GrowthData & obj){obj.read(is); return is;}
        /**
         * Operators to write data to an output stream in ADMB format from a GrowthData object.
         */
        friend std::ostream&   operator <<(std::ostream & os,   GrowthData & obj){obj.write(os); return os;}
    };

    /**
     * Class encapsulating a dataset with male maturity ogives based on chela height data
     */
    class ChelaHeightData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static ostream& os;
        /** keyword indicating effort data */
        const static adstring KW_CHELAHEIGHT_DATA;
    public:
        /** dataset name */
        adstring name;
        /** survey name */
        adstring survey;
        /** likelihood function type */
        int llType; 
        /** likelihood weight (i.e., multiplier) */
        double llWgt;   
        /** number of observations */
        int nObs;
        /** input data (columns: year,size,N,fraction mature) */
        dmatrix  inpData_nc;  
        /** ivector of year corresponding to observations */
        ivector obsYear_n;
        /** dvector of observed sizes */
        dvector obsSize_n;
        /** dvector for sample sizes corresponding to observations */
        dvector obsSS_n;
        /** dvector for observed fraction mature */
        dvector obsPrMat_n;
        /** ivector of indices to model size bins corresponding to observed sizes */
        ivector obsSizeBinIndex_n;
    public:
        /**
         * Constructor.
         */
        ChelaHeightData();
        /**
         * Destructor.
         */
        ~ChelaHeightData();
        /**
         * Set the maximum year in which to include chela height.
         * 
         * @param mxYr - the max year in which to include chela height data
         */
        void setMaxYear(int mxYr, const dvector& zCs);    
        /**
         * Calculates indices for model size bins corresponding to observed sizes.
         * 
         * @param zCs - model size bin cutpoints
         */
        void calcSizeBinIndices(const dvector& zCs);
        /**
         * Read input data in ADMB format from a file stream
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write data to an output stream in ADMB format
         * 
         * @param os output stream
         */
        void write(std::ostream & os); //write object to file in ADMB format
        /**
         * Write data to an output stream as an R-formatted list object
         * 
         * @param os output stream
         */
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a ChelaHeightData object.
         */
        friend cifstream& operator >>(cifstream & is, ChelaHeightData & obj){obj.read(is); return is;}
        /**
         * Operators to write data to an output stream in ADMB format from an ChelaHeightData object.
         */
        friend std::ostream& operator <<(std::ostream & os, ChelaHeightData & obj){obj.write(os); return os;}
    };

    /**
     * Class encapsulating a dataset reflecting annual maturity ogives (male or female)
     */
    class MaturityOgiveData {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static ostream& os;
        /** keyword indicating effort data */
        const static adstring KW_MATURITYOGIVE_DATA;
    private:
        int mxYr; //max year to which data will be fitted
    public:
        /** dataset name */
        adstring name;
        /** survey name */
        adstring survey;
        /** sex (as int) */
        int sex;
        /** likelihood function type */
        int llType; 
        /** likelihood weight (i.e., multiplier) */
        double llWgt;   
        /**number of size bins used */
        int nZBs;
        /** cutpoints for the maturity ogives */
        dvector cutpts;
        /** number of observations */
        int nObs;
        /** input data (columns: year,size,N,fraction mature) */
        dmatrix  inpData_nc;  
        /** ivector of year corresponding to observations */
        ivector obsYear_n;
        /** dvector of observed sizes */
        dvector obsSize_n;
        /** dvector for sample sizes corresponding to observations */
        dvector obsSS_n;
        /** dvector for observed fraction mature */
        dvector obsPrMat_n;
        /** ivector of indices to model size bins corresponding to observed sizes */
        ivector obsZBI_n;
        /** matrix to re-map model size bins to maturity ogive size bins */
        dmatrix zbRemapper;
    public:
        /**
         * Constructor.
         */
        MaturityOgiveData();
        /**
         * Destructor.
         */
        ~MaturityOgiveData();
        /**
         * Set the maximum year in which to include maturity ogive.
         * 
         * @param mxYr - the max year in which to include maturity ogive
         */
        void setMaxYear(int mxYr);    
        /**
         * Calculates matrix to re-map model size bins maturity ogive size bins.
         * 
         * @param zCs - model size bin cutpoints
         */
        void calcSizeBinRemapper(const dvector& zCs);
        /**
         * Read input data in ADMB format from a file stream
         * 
         * @param is - input file stream
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write data to an output stream in ADMB format
         * 
         * @param os output stream
         */
        void write(std::ostream & os); //write object to file in ADMB format
        /**
         * Write data to an output stream as an R-formatted list object
         * 
         * @param os output stream
         */
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a MaturityOgiveData object.
         */
        friend cifstream& operator >>(cifstream & is, MaturityOgiveData & obj){obj.read(is); return is;}
        /**
         * Operator to write data to an output stream in ADMB format from a MaturityOgiveData object.
         */
        friend std::ostream& operator <<(std::ostream & os, MaturityOgiveData & obj){obj.write(os); return os;}
    };

//--------------------------------------------------------------------------------
//         ModelDatasets
//--------------------------------------------------------------------------------
    /**
     * Class encapsulating model dataset objects.
     */
    class ModelDatasets {
    public:
        /** flag to print debugging info */
        static int debug;
        /** stream to write debugging info */
        static ostream& os;
    public:
        /** pointer to ModelConfiguration object */
        ModelConfiguration* pMC;
        /** bio data file name */
        adstring fnBioData;
        /** pointer to bio dataset object */
        BioData* ptrBio;   
        
        /** number of fishery datasets to read */
        int nFsh;
        /** fishery data file names */
        adstring_array fnsFisheryData;      
        /** pointer to array of pointers to fishery dataset objects */
        FleetData**    ppFsh;         
        
        /** number of survey datasets to read */
        int nSrv;
        /** survey data files names */
        adstring_array fnsSurveyData;
        /** pointer to array of pointers to survey dataset objects */
        FleetData**    ppSrv;        
        
        /** number of growth datasets to read */
        int nGrw;
        /** growth data files names */
        adstring_array fnsGrowthData;
        /** pointer to array of pointers to growth dataset objects */
        GrowthData**    ppGrw;        
        
        /** number of chela height datasets to read */
        int nCHD;
        /** chela height data files names */
        adstring_array fnsChelaHeightData;
        /** pointer to array of pointers to chela height dataset objects */
        ChelaHeightData**    ppCHD;        
        
        /** number of maturity ogive datasets to read */
        int nMOD;
        /** maturity ogive data files names */
        adstring_array fnsMaturityOgiveData;
        /** pointer to array of pointers to maturity ogive dataset objects */
        MaturityOgiveData**    ppMOD;        
    public:
        ModelDatasets(ModelConfiguration* ptrMC);
        ~ModelDatasets();
        /**
         * Read from ModelDatasets file in ADMB format.
         * 
         * @param is - input filestream to read from
         */
        void read(cifstream & is);//read file in ADMB format
        /**
         * Write to ModelDatasets file in ADMB format.
         * 
         * @param fn - file name
         */
        void write(adstring fn); 
        /**
         * Write to ModelDatasets file in ADMB format.
         * 
         * @param os - output stream
         */
        void write(std::ostream & os); //write object to file in ADMB format
        /**
         * Write to MosdelDatasets file in R format.
         * 
         * @param os - output stream
         * @param nm - name of object to write
         * @param indent - number of tabs to indent by
         */
        void writeToR(std::ostream& os, std::string nm, int indent=0);//write object to R file as list
        /**
         * Operator to read ADMB-formatted data from an input stream into a ModelDatasets object.
         */
        friend cifstream&    operator >>(cifstream & is,   ModelDatasets & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os,ModelDatasets & obj){obj.write(os); return os;}
    };
#endif  /* MODELDATA_HPP */

