/* 
 * File:   ModelIndexFunctions.hpp
 * Author: William Stockhausen
 * 
 * Created on February 20, 2014, 4:55 PM
 */

#ifndef MODELINDEXFUNCTIONS_HPP
#define	MODELINDEXFUNCTIONS_HPP

#include <admodel.h>

/*-------------------------------------------------------\n
 * Abstract base class in a hierarchy of sub-classes that 
 * provide the ability to compute an index into a vector 
 * based on its equivalency to a multidimensional array.
 *-------------------------------------------------------*/
class IndexFunction{
    protected:
        int nDims;          //number of dimensions
        adstring_array dims;//dimension names
        ivector idx;        //vector of indices
    public:
        void allocate(void);
        int getDimensionality(){return nDims;}
        adstring getDimName(int i){return dims[i];}
        void setDimName(int i,adstring dim){dims[i]=dim;}
        void setDimNames(adstring_array& names){for (int i=1;i<=nDims;i++) dims[i]=names(i);}
        void defineIndexMapping(imatrix m);
        
        virtual int getSize()=0;
        virtual int calcOffset(ivector i)=0; //calc offset into idx vector
        virtual int getIndex(ivector i)=0;  //get value of index corresponding to indices vector
};

/*-------------------------------------------------------\n
 * Class that provides the ability to compute an index 
 * into a vector based on its equivalency to another vector.
 *-------------------------------------------------------*/
class IndexFunction1:public IndexFunction{
    public:
        int mn;
        int mx;
    public:
        IndexFunction1(int imn, int imx);
        ~IndexFunction1(){}
        
        virtual int getSize(){return mx-mn+1;}
        virtual int calcOffset(ivector i);
        virtual int getIndex(ivector i);
        int getIndex(int i);
};

/*-------------------------------------------------------\n
 * Class that provides the ability to compute an index 
 * into a vector based on its equivalency to a matrix.
 *-------------------------------------------------------*/
class IndexFunction2: public IndexFunction1{
    public:
        int mn;
        int mx;
    public:
        IndexFunction2(int imn, int imx, int jmn, int jmx);
        ~IndexFunction2(){}
        
        virtual int getSize(){return (mx-mn+1)*IndexFunction1::getSize();}
        virtual int calcOffset(ivector i);
        virtual int getIndex(ivector i);
        int getIndex(int i, int j);
};

/*-------------------------------------------------------\n
 * Class that provides the ability to compute an index 
 * into a vector based on its equivalency to a 3d array.
 *-------------------------------------------------------*/
class IndexFunction3: public IndexFunction2{
    public:
        int mn;
        int mx;
    public:
        IndexFunction3(int imn, int imx, int jmn, int jmx, int kmn, int kmx);
        ~IndexFunction3(){}
        
        virtual int getSize(){return (mx-mn+1)*IndexFunction2::getSize();}
        virtual int calcOffset(ivector i);
        virtual int getIndex(ivector i);
        int getIndex(int i,int j, int k);
};

/*-------------------------------------------------------\n
 * Class that provides the ability to compute an index 
 * into a vector based on its equivalency to a 4d array.
 *-------------------------------------------------------*/
class IndexFunction4: public IndexFunction3{
    public:
        int mn;
        int mx;
    public:
        IndexFunction4(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx);
        ~IndexFunction4(){}
        
        virtual int getSize(){return (mx-mn+1)*IndexFunction3::getSize();}
        virtual int calcOffset(ivector i);
        virtual int getIndex(ivector i);
        int getIndex(int i,int j, int k, int l);
};

/*-------------------------------------------------------\n
 * Class that provides the ability to compute an index 
 * into a vector based on its equivalency to a 5d array.
 *-------------------------------------------------------*/
class IndexFunction5: public IndexFunction4{
    public:
        int mn;
        int mx;
    public:
        IndexFunction5(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx, int mmn, int mmx);
        ~IndexFunction5(){}
        
        virtual int getSize(){return (mx-mn+1)*IndexFunction4::getSize();}
        virtual int calcOffset(ivector i);
        virtual int getIndex(ivector i);
        int getIndex(int i, int j, int k, int l, int m);
};

/*-------------------------------------------------------\n
 * Class that provides the ability to compute an index 
 * into a vector based on its equivalency to a 6d array.
 *-------------------------------------------------------*/
class IndexFunction6: public IndexFunction5{
    public:
        int mn;
        int mx;
    public:
        IndexFunction6(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx, int mmn, int mmx, int nmn, int nmx);
        ~IndexFunction6(){}
        
        virtual int getSize(){return (mx-mn+1)*IndexFunction5::getSize();}
        virtual int calcOffset(ivector i);
        virtual int getIndex(ivector i);
        int getIndex(int i, int j, int k, int l, int m, int n);
};

/*-------------------------------------------------------\n
 * Class that provides the ability to compute an index 
 * into a vector based on its equivalency to a 7d array.
 *-------------------------------------------------------*/
class IndexFunction7: public IndexFunction6{
    public:
        int mn;
        int mx;
    public:
        IndexFunction7(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx, int mmn, int mmx, int nmn, int nmx, int omn, int omx);
        ~IndexFunction7(){}
        
        virtual int getSize(){return (mx-mn+1)*IndexFunction6::getSize();}
        virtual int calcOffset(ivector i);
        virtual int getIndex(ivector i);
        int getIndex(int i, int j, int k, int l, int m, int n, int o);
};

#endif	/* MODELINDEXFUNCTIONS_HPP */

