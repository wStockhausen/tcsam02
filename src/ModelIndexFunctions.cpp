#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelIndexFunctions.hpp"


/*----------------------------------------------------------------------------*/
void IndexFunction::allocate(void){
    dims.allocate(1,nDims); idx.allocate(1,getSize()); idx = 0;
}
void IndexFunction::defineIndexMapping(imatrix m){
    int ic = nDims+1;
    ivector iv(1,nDims);
    int o = 0;
    for (int r=m.rowmin();r<=m.rowmax();r++){
        iv = 0;
        for (int c=1;c<=nDims;c++) iv(c) = m(r,c);
        o = calcOffset(iv);    
        if (o) idx(o) = m(r,ic);
    }
}

IndexFunction1::IndexFunction1(int imn, int imx){
    nDims=1;mn=imn;mx=imx;
}
int IndexFunction1::calcOffset(ivector i){
    if (i[1]>mx) return 0; if (i[1]<mn) return 0; return i[1]-mn+1;
}
int IndexFunction1::getIndex(ivector i){
    int o = calcOffset(i);
    if (o) return idx(o);
    return o;
}
int IndexFunction1::getIndex(int i){
    ivector iv(1,nDims);
    iv(1) = i;
    return getIndex(iv);
}

IndexFunction2::IndexFunction2(int imn, int imx, int jmn, int jmx):IndexFunction1(jmn,jmx){
    nDims=2;mn=imn;mx=imx;
}
int IndexFunction2::calcOffset(ivector i){
    if (i[1]>mx) return 0; 
    if (i[1]<mn) return 0; 
    int idx = IndexFunction1::calcOffset(i(2,2));
    if (!idx) return 0;
    return (i[1]-mn)*IndexFunction1::getSize()+idx;
}
int IndexFunction2::getIndex(ivector i){
    int o = calcOffset(i);
    if (o) return idx(o);
    return o;
}
int IndexFunction2::getIndex(int i, int j){
    ivector iv(1,nDims);
    iv(1) = i; iv(2) = j;
    return getIndex(iv);
}

IndexFunction3::IndexFunction3(int imn, int imx, int jmn, int jmx, int kmn, int kmx):IndexFunction2(jmn,jmx,kmn,kmx){
    nDims=3;mn=imn;mx=imx;
}
int IndexFunction3::calcOffset(ivector i){
    if (i[1]>mx) return 0; 
    if (i[1]<mn) return 0; 
    int idx = IndexFunction2::calcOffset(i(2,3));
    if (!idx) return 0;
    return (i[1]-mn)*IndexFunction2::getSize()+idx;
}
int IndexFunction3::getIndex(ivector i){
    int o = calcOffset(i);
    if (o) return idx(o);
    return o;
}
int IndexFunction3::getIndex(int i, int j, int k){
    ivector iv(1,nDims);
    iv(1) = i; iv(2) = j; iv(3) = k;
    return getIndex(iv);
}

IndexFunction4::IndexFunction4(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx):IndexFunction3(jmn,jmx,kmn,kmx,lmn,lmx){
    nDims=4;mn=imn;mx=imx;
}
int IndexFunction4::calcOffset(ivector i){
    if (i[1]>mx) return 0; 
    if (i[1]<mn) return 0; 
    int idx = IndexFunction3::calcOffset(i(2,4));
    if (!idx) return 0;
    return (i[1]-mn)*IndexFunction3::getSize()+idx;
}
int IndexFunction4::getIndex(ivector i){
    int o = calcOffset(i);
    if (o) return idx(o);
    return o;
}
int IndexFunction4::getIndex(int i, int j, int k, int l){
    ivector iv(1,nDims);
    iv(1) = i; iv(2) = j; iv(3) = k; iv(4) = l;
    return getIndex(iv);
}

IndexFunction5::IndexFunction5(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx, int mmn, int mmx):IndexFunction4(jmn,jmx,kmn,kmx,lmn,lmx,mmn,mmx){
    nDims=5;mn=imn;mx=imx;
}
int IndexFunction5::calcOffset(ivector i){
    if (i[1]>mx) return 0; 
    if (i[1]<mn) return 0; 
    int idx = IndexFunction4::calcOffset(i(2,5));
    if (!idx) return 0;
    return (i[1]-mn)*IndexFunction4::getSize()+idx;
}
int IndexFunction5::getIndex(ivector i){
    int o = calcOffset(i);
    if (o) return idx(o);
    return o;
}
int IndexFunction5::getIndex(int i, int j, int k, int l, int m){
    ivector iv(1,nDims);
    iv(1) = i; iv(2) = j; iv(3) = k; iv(4) = l; iv(5) = m;
    return getIndex(iv);
}

IndexFunction6::IndexFunction6(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx, int mmn, int mmx, int nmn, int nmx):IndexFunction5(jmn,jmx,kmn,kmx,lmn,lmx,mmn,mmx,nmn,nmx){
    nDims=6;mn=imn;mx=imx;
}
int IndexFunction6::calcOffset(ivector i){
    if (i[1]>mx) return 0; 
    if (i[1]<mn) return 0; 
    int idx = IndexFunction5::calcOffset(i(2,6));
    if (!idx) return 0;
    return (i[1]-mn)*IndexFunction6::getSize()+idx;
}
int IndexFunction6::getIndex(ivector i){
    int o = calcOffset(i);
    if (o) return idx(o);
    return o;
}
int IndexFunction6::getIndex(int i, int j, int k, int l, int m, int n){
    ivector iv(1,nDims);
    iv(1) = i; iv(2) = j; iv(3) = k; iv(4) = l; iv(5) = m; iv(6) = n;
    return getIndex(iv);
}


IndexFunction7::IndexFunction7(int imn, int imx, int jmn, int jmx, int kmn, int kmx,int lmn, int lmx, int mmn, int mmx, int nmn, int nmx, int omn, int omx):IndexFunction6(jmn,jmx,kmn,kmx,lmn,lmx,mmn,mmx,nmn,nmx,omn,omx){
    nDims=7;mn=imn;mx=imx;
}
int IndexFunction7::calcOffset(ivector i){
    if (i[1]>mx) return 0; 
    if (i[1]<mn) return 0; 
    int idx = IndexFunction6::calcOffset(i(2,7));
    if (!idx) return 0;
    return (i[1]-mn)*IndexFunction6::getSize()+idx;
}
int IndexFunction7::getIndex(ivector i){
    int o = calcOffset(i);
    if (o) return idx(o);
    return o;
}
int IndexFunction7::getIndex(int i, int j, int k, int l, int m, int n, int o){
    ivector iv(1,nDims);
    iv(1) = i; iv(2) = j; iv(3) = k; iv(4) = l; iv(5) = m; iv(6) = n; iv(7) = o;
    return getIndex(iv);
}

