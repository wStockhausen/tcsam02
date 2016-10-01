#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelIndexBlocks.hpp"
#include "ModelConfiguration.hpp"
#include "ModelData.hpp"


int IndexRange::debug = 0;
int IndexBlock::debug = 0;
int IndexBlockSet::debug = 0;
/*----------------------------------------------------------------------------*/
/**
 * Construct an IndexRange that will substitute the given 
 * min and/or max limits when defaults (<0) are specified on input. 
 * @param modMn
 * @param modMx
 */
IndexRange::IndexRange(int modMn,int modMx){
    modMin = modMn; modMax = modMx;
}
/**
 * Construct index vector covering actual range from min to max.
 * @param min
 * @param max
 */
void IndexRange::createRangeVector(int min, int max){
    mn = min; mx = max;
    int n = mx-mn+1;
    iv.deallocate(); iv.allocate(1,n);
    for (int i=0;i<n;i++) iv(i+1)=mn+i;
}
/**
 * Parse a range string ("x:y" or "y") to obtain actual min, max for range.
 * If result finds x (y) < 0, then x (y) will be substituted using
 *  if x<0, x = modMin+1+x (so x=-1->x=modMin, x=-2->x=modMin-1, etc)
 *  if y<0, y = modMax-1-y (so y=-1->y=modMax, y=-2->y=modMax+1, etc)
 * @param str
 */
void IndexRange::parse(adstring str){
    if (debug) cout<<"IndexRange parse("<<str<<")"<<endl;
    int i = str.pos(":");
    if (debug) cout<<"':' located at position "<<i<<endl;
    if (i){
        mn = ::atoi(str(1,i-1));          
        if (mn<0){mn=modMin+1+mn;} else
        if (mn<modMin){mn=modMin;} else
        if (mn>modMax){mn=modMax;}
        mx = ::atoi(str(i+1,str.size())); 
        if (mx<0){mx=modMax-1-mx;} else
        if (mx<modMin){mx=modMin;} else
        if (mx>modMax){mx=modMax;}
    } else {
        mx = ::atoi(str);
        if (mx<0){mx=modMax-1-mx;} else
        if (mx<modMin){mx=modMin;} else
        if (mx>modMax){mx=modMax;}
        mn=mx;
    }
    if (debug) cout<<"mn,mx = "<<mn<<cc<<mx<<endl;
    createRangeVector(mn,mx);
}
/**
 * Read string ("x:y" or "x") from the filestream object and parse it
 * to obtain actual min, max for range.
 * @param is
 */
void IndexRange::read(cifstream& is){
    if (debug) cout<<"starting void IndexRange::read(is)"<<endl;
    adstring str; 
    is>>str;
    parse(str);
    if (debug) {
        cout<<"input = '"<<str<<"'"<<endl;
        cout<<"mn = "<<mn<<tb<<"mx = "<<mx<<endl;
        cout<<"iv = "<<iv<<endl;
        cout<<"finished void IndexRange::read(is)"<<endl;
    }
}
/**
 * Writes range to output stream in parse-able format.
 * @param os
 */
void IndexRange::write(std::ostream& os){
    os<<asString();
}
/**
 * Writes range to file as R vector min:max.
 * @param os
 */
void IndexRange::writeToR(std::ostream& os){
    os<<mn<<":"<<mx;
}
/*
 * Returns a range as an adstring object.
 */
adstring IndexRange::asString(){
    adstring a;
    if (iv.size()>1) a = str(mn)+":"+str(mx);
    else a = str(mn);
    return(a);
}

/**
 * Destructor for class.
 */
IndexBlock::~IndexBlock(){
    if (ppIRs) {
        for (int i=0;i<nRCs;i++) delete ppIRs[i]; 
        delete ppIRs;
    }
}

/**
 * Creates vectors of indices that map from block to model (forward)\n
 * and from model to block (reverse).\n
 */
void IndexBlock::createIndexVectors(){
    int mnM = modMin; 
    int mxM = modMax;
    nIDs = 0;
    for (int c=0;c<nRCs;c++) {
        mnM = min(mnM,ppIRs[c]->getMin());//need to do these checks
        mxM = max(mxM,ppIRs[c]->getMax());//need to do these checks
        nIDs += ppIRs[c]->getMax()-ppIRs[c]->getMin()+1;
    }
    ivFwd.allocate(1,nIDs);
    ivRev.allocate(mnM,mxM);
    ivRev = 0;//model elements that DON'T map to block will have value "0"
    nIDs = 0;//reset nIDs to use as counter
    for (int c=0;c<nRCs;c++) {
        int mx = ppIRs[c]->getMax();
        int mn = ppIRs[c]->getMin();
        for (int i=mn;i<=mx;i++) {
            ivFwd[++nIDs] = i;
            ivRev[i] = nIDs;
        }
    }
    //nIDs should now be same as before
}

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
void IndexBlock::parse(adstring str1){
    if (debug) cout<<"Parsing '"<<str1<<"'"<<endl;
    if (!(str1(1,1)=="[")) {
        cout<<"Error parsing index block '"<<str1<<"'."<<endl;
        cout<<"Expected 1st character to be '[' but got '"<<str1(1,1)<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    adstring str2 = str1(2,str1.size()-1);
    int k = 0; 
    int p=str2.pos(";");
    if (debug) cout<<"k,p,str = "<<k<<tb<<p<<tb<<str2<<endl;
    while(p){ 
        k++;
        str2 = str2(p+1,str2.size());
        p = str2.pos(";");
        if (debug) cout<<"k,p,str = "<<k<<tb<<p<<tb<<str2<<endl;
    }
    k++;
    if (debug) cout<<"Found "<<k<<" intervals"<<endl;
    
    nRCs = k;
    str2 = str1(2,str1.size()-1);
    ppIRs = new IndexRange*[nRCs];
    for (k=1;k<nRCs;k++){
        p = str2.pos(";");
        ppIRs[k-1] = new IndexRange(modMin,modMax);
        ppIRs[k-1]->parse(str2(1,p-1));
        str2 = str2(p+1,str2.size());
    }
    ppIRs[nRCs-1] = new IndexRange(modMin,modMax);
    ppIRs[nRCs-1]->parse(str2);
    if (debug) {
        cout<<"IndexBlock = "<<(*this)<<endl;
        cout<<"Enter 1 to contiinue > ";
        cin>>debug;
    }
    createIndexVectors();
}
/* 
 * Reads an adstring from an input stream and parses it
 * to define the IndexBlock.
 */
void IndexBlock::read(cifstream & is){
    adstring str1; 
    is>>str1;
    parse(str1);
}
/*
 * Writes an IndexBlock in ADMB format to an output stream.
 */
void IndexBlock::write(std::ostream & os){
    os<<asString();
}
/*
 * Writes an IndexBlock as an unnamed character string.
 */
void IndexBlock::writeToR(std::ostream& os){
    os<<"'"<<asString()<<os<<"'";
}
/*
 * Returns an IndexBlock as an adstring object.
 */
adstring IndexBlock::asString(){
    adstring a = "[";
    for (int i=1;i<nRCs;i++) a += ppIRs[i-1]->asString()+";";
    a += ppIRs[nRCs-1]->asString()+"]";
    return a;
}

/*----------------------------------------------------------------------------*/
IndexBlockSet::~IndexBlockSet(){
    if (ppIBs) {
        for (int i=0;i<nIBs;i++) delete ppIBs[i]; 
        delete ppIBs;
    }
}

/**
 * Allocate n IndexBlocks for this IndexBlockSet.
 * @param n
 */
void IndexBlockSet::allocate(int n){
    nIBs = n;
    ppIBs = new IndexBlock*[nIBs];
    for (int i=1;i<=nIBs;i++) ppIBs[i-1] = new IndexBlock(modMin,modMax);
}

/**
 * Sets the dimension type for this IndexBlockSet.
 * @param theType
 */
void IndexBlockSet::setType(adstring theType){
    type = theType;
    int p = theType.pos("_");
    if (p)  type = theType(1,p-1);//strip off type
    if (type==tcsam::STR_YEAR){
        modMin = ModelConfiguration::mnYr;
        modMax = ModelConfiguration::mxYr;
    } else
    if (type==tcsam::STR_SEX){
        modMin = 1;
        modMax = tcsam::nSXs;
    } else
    if (type==tcsam::STR_MATURITY_STATE){
        modMin = 1;
        modMax = tcsam::nMSs;
    } else
    if (type==tcsam::STR_SHELL_CONDITION){
        modMin = 1;
        modMax = tcsam::nSCs;
    } else 
    if (type==tcsam::STR_SIZE){
        modMin = 1;
        modMax = ModelConfiguration::nZBs;
    } else 
    if (type==tcsam::STR_FISHERY){
        modMin = 1;
        modMax = ModelConfiguration::nFsh;
    } else 
    if (type==tcsam::STR_SURVEY){
        modMin = 1;
        modMax = ModelConfiguration::nSrv;
    } else 
    {
        cout<<"WARNING...WARNING...WARNING..."<<endl;
        cout<<"Defining non-standard index type '"<<type<<"'."<<endl;
        cout<<"Make sure this is what you want."<<endl;
        cout<<"WARNING...WARNING...WARNING..."<<endl;
    }
    if (debug) cout<<"modMin = "<<modMin<<tb<<"modMax = "<<modMax<<endl;
}

/* 
 * Reads an IndexBlockSet in ADMB format from an input stream.
 */
void IndexBlockSet::read(cifstream & is){
    if (debug) cout<<"Starting IndexBlockSet::read(cifstream & is)"<<endl;
    if (type=="") {
        if (debug) cout<<"reading type"<<endl;
        is>>type;
        setType(type);
    }
    is>>nIBs;//number of index blocks
    ppIBs = new IndexBlock*[nIBs];
    int id;
    for (int i=0;i<nIBs;i++) ppIBs[i] = new IndexBlock(modMin,modMax);
    for (int i=0;i<nIBs;i++){
        is>>id;
        is>>(*ppIBs[id-1]);
        if (debug) cout<<id<<tb<<(*ppIBs[id-1])<<endl;
    }
    if (debug) cout<<"Finished IndexBlockSet::read(cifstream & is)"<<endl;
}
/*
 * Writes an IndexBlockSet in ADMB format to an output stream.
 */
void IndexBlockSet::write(std::ostream & os){
    os<<type<<tb<<"#index type (dimension name)"<<endl;
    os<<nIBs<<tb<<"#number of index blocks defined"<<endl;
    os<<"#id  Blocks"<<endl;
    for (int i=1;i<nIBs;i++) os<<i<<tb<<(*ppIBs[i-1])<<endl;
    os<<nIBs<<tb<<(*ppIBs[nIBs-1]);
}
/*
 * Writes an IndexBlockSet as an un-named R list.
 */
void IndexBlockSet::writeToR(std::ostream& os){
    os<<"list(type="<<type<<cc<<"nIBs="<<nIBs<<cc;
    os<<"IBs=c("; 
    for (int i=1;i<nIBs;i++) {ppIBs[i-1]->writeToR(os); os<<cc;}
    ppIBs[nIBs-1]->writeToR(os); os<<")"<<")";
}



void tcsam::getIndexLimits(adstring& idxType,int& mn,int& mx){
    if (idxType==tcsam::STR_YEAR) {
        mn = ModelConfiguration::mnYr; 
        mx = ModelConfiguration::mxYr;
    }else
    if (idxType==tcsam::STR_SIZE) {
        mn = 1;
        mx = ModelConfiguration::nZBs;
    } else
    if (idxType==tcsam::STR_SEX)  {
        mn = 1;
        mx = tcsam::nSXs;
    } else
    if (idxType==tcsam::STR_MATURITY_STATE)  {
        mn = 1;
        mx = tcsam::nMSs;
    } else
    if (idxType==tcsam::STR_SHELL_CONDITION) {
        mn = 1;
        mx = tcsam::nSCs;
    } else
    {
        cout<<"Error in tcsam::getIndexLimits("<<idxType<<")"<<endl;
        cout<<"Unrecognized idxType '"<<idxType<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
}