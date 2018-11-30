/* 
 * File:   ModelIndexBlocks.hpp
 * Author: WilliamStockhausen
 * 
 * Includes:
 *  class IndexRange
 *  class IndexBlock
 *  class IndexBlockSet
 *  class IndexBlockSets
 *
 * Created on March 11, 2014, 7:04 AM
 */

#ifndef MODELINDEXBLOCKS_HPP
#define	MODELINDEXBLOCKS_HPP

#include <admodel.h>

/*
 * This class is intended to simplify the use of index ranges
 * in the definition of index blocks. It can interpret single numbers
 * or strings of the form "min:max" as a range of indices. 
 */
class IndexRange{
    public:
        static int debug;
    protected:
        int modMin;//minimum possible value in model for index
        int modMax;//maximum possible value in model for index
        int mn;
        int mx;
        ivector iv;
    public:
        /*
         * Constructor with default range values.
         */
        IndexRange(int modMn, int modMx);
        ~IndexRange(){}
        
        int getMax(){return mx;}
        int getMin(){return mn;}
        /* 
         * Creates an ivector with elements that span the range.
         * iv[i+1] = min+i for i=0:[max-min]
         */
        void createRangeVector(int min, int max);    
        /* 
         * Returns an ivector with elements that span the range.
         * iv[i+1] = min+i for i=0:[max-min]
         */
        ivector getRangeVector(){return iv;}
        /**
         * Parse a range string ("x:y" or "y") to obtain actual min, max for range.
         * If result finds x (y) < 0, then x (y) will be substituted using
         *  if x<0, x = modMin+1+x (so x=-1->x=modMin, x=-2->x=modMin-1, etc)
         *  if y<0, y = modMax-1-y (so y=-1->y=modMax, y=-2->y=modMax+1, etc)
         * @param str
         */
        void parse(adstring str);
        /* 
         * Reads a range in ADMB format (single number or range formatted as "min:max") 
         * from an input stream.
         */
        void read(cifstream & is);      //read object in ADMB format
        /*
         * Writes the range in ADMB format to an output stream.
         */
        void write(std::ostream & os);  //write object to file in ADMB format
        /*
         * Writes the range as an un-named R list.
         */
        void writeToR(std::ostream& os);//write object to R file as list
        /*
         * Returns the range as an adstring object.
         */
        adstring asString();
        
        friend cifstream&    operator >>(cifstream & is,   IndexRange & obj){obj.read(is);return is;}
        friend std::ostream& operator <<(std::ostream & os,IndexRange & obj){obj.write(os);return os;}
};

/*-------------------------------------------------------------\n
 * Class that defines an index block (e.g., a time block).
 *-------------------------------------------------------------*/
class IndexBlock{
    public:
        static int debug;//debug flag (0=OFF, 1=ON)
    public:
        int modMin;//minimum possible value in model for index
        int modMax;//maximum possible value in model for index
        int nRCs;  //number of individual range components defining block
        int nIDs;  //number of block indices
        ivector ivFwd;//forward index vector from block indices to model indices (1:nIDs)
        ivector ivRev;//reverse index vector from model indices to block indices (modMin:modMax)
        IndexRange** ppIRs;//pointer to vector of IndexRange pointers
    private:
        adstring rDim;//as an R array dimension
    public:
        IndexBlock(int modMn, int modMx){nRCs=0;nIDs=0;ppIRs=0;modMin=modMn;modMax=modMx;}
        ~IndexBlock();
        
        /**
         * Returns size (max index=nIDs) of associated index vector.
         * @return 
         */
        int getSize(void){return nIDs;}    
        /**
         * Get the minimum index across the block.
         * 
         * @return the min index
         */
        int getMin(void);
        /**
         * Get the maximum index across the block.
         * 
         * @return the max index
         */
        int getMax(void);
        
        /**
         * Add an element to the end of the block.
         * 
         * @param e - value of element to add
         */
        void addElement(int e);
        
        /** 
         * Returns an ivector that maps block indices (1:nIDs) to model indices.\n
         * Thus, if iv is the ivector result of the function, iv(i) is the model\n
         * index value corresponding to the ith block element.
         * Indices of the resulting ivector run 1:nIDs.\n
         * 
         * @return forward index vector
         */
        ivector getFwdIndexVector(){return ivFwd;}        
        /** 
         * Returns an ivector that maps model indices to block indices (1:nIDs).\n
         * Thus, if iv is the ivector result of the function, iv(i) is the block\n
         * index value corresponding to the model index i.
         * Indices of the resulting ivector run modMin:modMax.\n
         * 
         * @return - reverse index vector
         */
        ivector getRevIndexVector(){return ivRev;}        
        /**
         * Parse the string str1 as an index block. 
         * String str1 must start with "[" and end with "]".
         * Individual ranges are separated using a semi-colon (";").
         * Individual ranges have the form "x:y" or "x".
         * 
         * An example is "[1962:2000;2005;-1:1959]". In this example, the 
         * "-1" would be replaced by the specified minimum range for the block.
         * 
         * @param str1 - adstring to parse
         */
        void parse(adstring str1);        
        /* 
         * Reads an adstring from an input stream and parses it
         * to define the IndexBlock.
         */
        void read(cifstream & is);      //read object in ADMB format
        /*
         * Writes an IndexBlock to an output stream in ADMB format.\n
         */
        void write(std::ostream & os);  //write object to file in ADMB format
        /*
         * Writes an IndexBlock as an un-named R list.\n
         */
        void writeToR(std::ostream& os);//write object to R file as list
        /**
         * Creates a representation of the IndexBlock as an R array dimension
         */
        void createRDim();
        /*
         * Returns IndexBlock as an R array dimension.\n
         */
        adstring getAsRDim(){return rDim;}
        /*
         * Returns an IndexBlock as an adstring object.
         */
        adstring asString();
    protected:
        /**
         * Creates vectors of indices that map from block to model (forward)\n
         * and from model to block (reverse).\n
         */
        void createIndexVectors(void);
        
    public:    
        friend cifstream&    operator >>(cifstream & is,   IndexBlock & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os,IndexBlock & obj){obj.write(os);return os;}
};

/*-------------------------------------------------------------\n
 * Class that defines an index block set (e.g., a set of 
 * time blocks).
 *-------------------------------------------------------------*/
class IndexBlockSet{
    public:
        static int debug;//debug flag (0=OFF, 1=ON)
    protected:
        adstring type;//index type (dimension name)
        int modMin;//minimum possible value in model for index
        int modMax;//maximum possible value in model for index
        int nIBs;//number of index blocks
        IndexBlock** ppIBs;//pointer to a vector of pointers to IndexBlock objects
    public:
        IndexBlockSet(){type=""; nIBs=0;ppIBs=0;modMin=-1;modMax=-1;}
        IndexBlockSet(adstring idxType, int modMn, int modMx){type=idxType; nIBs=0;ppIBs=0;modMin=modMn;modMax=modMx;}
        ~IndexBlockSet();
        
        /**
         * Allocate n IndexBlocks for this IndexBlockSet.
         * @param n
         */
        void allocate(int n);
        
        /**
         * Get the number of IndexBlocks in this set.
         * 
         * @return the number of IndexBlocks.
         */
        int getSize(void){return nIBs;}
        
        /**
         * Gets the ith IndexBlock in this set.
         * @param i - index in interval 1:nIBs 
         * @return 
         */
        IndexBlock* getIndexBlock(int i){return ppIBs[i-1];}
        
        /*
         * Return the index type (dimension name) associated with this IndexBlockSet.
         */
        adstring getType(void){return type;}
        
        /**
         * Sets the dimension type associated with this IndexBlockSet.
         * @param theType
         */
        void setType(adstring theType);        
        
        /** 
         * Returns the forward index vector corresponding to the ith index block.
         * 
         * @param int i:  index (starting at 1) of the index block to access
         * 
         * @return - ivector iv, where iv(i) is the model index value \n
         *           corresponding to the ith block element. \n
         *           Limits are 1:nIDs. \n
         */
        ivector getFwdIndexVector(int i){return ppIBs[i-1]->getFwdIndexVector();}
        
        /** 
         * Returns the reverse index vector corresponding to the ith index block.
         * 
         * @param int i:  index (starting at 1) of the index block to access
         * 
         * @return - ivector iv, where iv(i) is the block index \n
         *           corresponding to the model index value i. \n
         *           If iv(i) is 0, there is no corresponding block element.\n
         *           Limits are modMin:modMax. \n
         */
        ivector getRevIndexVector(int i){return ppIBs[i-1]->getRevIndexVector();}
        
        /* 
         * Reads an IndexBlockSet in ADMB format from an input stream.\n
         */
        void read(cifstream & is);      //read object in ADMB format
        /*
         * Writes an IndexBlockSet in ADMB format to an output stream.\n
         */
        void write(std::ostream & os);  //write object to file in ADMB format
        /*
         * Writes an IndexBlockSet as an un-named R list.\n
         */
        void writeToR(std::ostream& os);//write object to R file as list
        
    public:    
        friend cifstream&    operator >>(cifstream & is,   IndexBlockSet & obj){obj.read(is); return is;}
        friend std::ostream& operator <<(std::ostream & os,IndexBlockSet & obj){obj.write(os);return os;}
};

namespace tcsam {
    void getIndexLimits(adstring& idxType,int& mn,int& mx);
}
#endif	/* MODELINDEXBLOCKS_HPP */

