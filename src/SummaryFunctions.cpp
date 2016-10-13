#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "SummaryFunctions.hpp"
    
    
/**
 * Calculate biomass (1000's t) associated with abundance array.
 * 
 * @param n_xmsz - abundance array (millions)
 * @param w_xmz - individual weight-at-size by sex, maturity state (kg)
 * @return associated biomass as d4_array (1000's t)
 */
d4_array tcsam::calcBiomass(const d4_array& n_xmsz, const d3_array& w_xmz){
    d4_array b_xmsz = n_xmsz;
    for (int x=n_xmsz.indexmin();x<=n_xmsz.indexmax();x++){
        for (int m=n_xmsz(x).indexmin();m<=n_xmsz(x).indexmax();m++){
            for (int s=n_xmsz(x,m).indexmin();s<=n_xmsz(x,m).indexmax();s++) 
                b_xmsz(x,m,s) = elem_prod(w_xmz(x,m),n_xmsz(x,m,s));
        }
    }
    return b_xmsz;
}

/**
 * Calculate biomass (1000's t) associated with abundance array.
 * 
 * @param n_yxmsz - abundance array (millions)
 * @param w_xmz - individual weight-at-size by sex, maturity state (kg)
 * @return associated biomass as d5_array (1000's t)
 */
d5_array tcsam::calcBiomass(const d5_array& n_yxmsz, const d3_array& w_xmz){
    d5_array b_yxmsz = n_yxmsz;
    for (int y=n_yxmsz.indexmin();y<=n_yxmsz.indexmax();y++) b_yxmsz(y) = calcBiomass(n_yxmsz(y),w_xmz);
    return b_yxmsz;
}

/**
 * Calculate biomass (1000's t) associated with abundance array.
 * 
 * @param n_fyxmsz - abundance array (millions)
 * @param w_xmz - individual weight-at-size by sex, maturity state (kg)
 * @return associated biomass as d6_array (1000's t)
 */
d6_array tcsam::calcBiomass(const d6_array& n_fyxmsz, const d3_array& w_xmz){
    d6_array b_fyxmsz = n_fyxmsz;
    for (int f=n_fyxmsz.indexmin();f<=n_fyxmsz.indexmax();f++) b_fyxmsz(f) = calcBiomass(n_fyxmsz(f),w_xmz);
    return b_fyxmsz;
}

/**
 * Compute marginal weighted sums over last dimension.
 * 
 * @param  n_ij - dmatrix
 * @param  w    - weighting vector
 * @return n_i  - dvector
 */
dvector tcsam::sumOverLastDim(const dmatrix& n_ij, dvector& w){
    ivector bnds = wts::getBounds(n_ij);
    int d = 1;
    int mni = bnds(d++); int mxi = bnds(d++);
    dvector n(mni,mxi);
    for (int i=mni;i<=mxi;i++){
        n(i) = n_ij(i)*w;//dot product sum over j index
    }
    return n;
}

/**
 * Compute marginal weighted sums over last dimension.
 * 
 * @param  n_ijk - d3_array
 * @param  w     - weighting vector
 * @return n_ij  - dmatrix
 */
dmatrix tcsam::sumOverLastDim(const d3_array& n_ijk, dvector& w){
    ivector bnds = wts::getBounds(n_ijk);
    int d = 1;
    int mni = bnds(d++); int mxi = bnds(d++);
    int mnj = bnds(d++); int mxj = bnds(d++);
    dmatrix n(mni,mxi,mnj,mxj);
    for (int i=mni;i<=mxi;i++){
        n(i) = sumOverLastDim(n_ijk(i),w);
    }
    return n;
}

/**
 * Compute marginal weighted sums over last dimension.
 * 
 * @param  n_ijkl - d4_array
 * @param  w      - weighting vector
 * @return n_ijk  - d3_array
 */
d3_array tcsam::sumOverLastDim(const d4_array& n_ijkl, dvector& w){
    ivector bnds = wts::getBounds(n_ijkl);
    int d = 1;
    int mni = bnds(d++); int mxi = bnds(d++);
    int mnj = bnds(d++); int mxj = bnds(d++);
    int mnk = bnds(d++); int mxk = bnds(d++);
    d3_array n(mni,mxi,mnj,mxj,mnk,mxk);
    for (int i=mni;i<=mxi;i++){
        n(i) = sumOverLastDim(n_ijkl(i),w);
    }
    return n;
}

/**
 * Compute marginal weighted sums over last dimension.
 * 
 * @param  n_ijklm - d5_array
 * @param  w    - weighting vector
 * @return n_ijkl  - d4_array
 */
d4_array tcsam::sumOverLastDim(const d5_array& n_ijklm, dvector& w){
    ivector bnds = wts::getBounds(n_ijklm);
    int d = 1;
    int mni = bnds(d++); int mxi = bnds(d++);
    int mnj = bnds(d++); int mxj = bnds(d++);
    int mnk = bnds(d++); int mxk = bnds(d++);
    int mnl = bnds(d++); int mxl = bnds(d++);
    d4_array n(mni,mxi,mnj,mxj,mnk,mxk,mnl,mxl);
    for (int i=mni;i<=mxi;i++){
        n(i) = sumOverLastDim(n_ijklm(i),w);
    }
    return n;
}

/**
 * Compute marginal weighted sums over last dimension.
 * 
 * @param  n_ijklmn - d6_array
 * @param  w    - weighting vector
 * @return n_ijklm  - d4_array
 */
d5_array tcsam::sumOverLastDim(const d6_array& n_ijklmn, dvector& w){
    ivector bnds = wts::getBounds(n_ijklmn);
    int d = 1;
    int mni = bnds(d++); int mxi = bnds(d++);
    int mnj = bnds(d++); int mxj = bnds(d++);
    int mnk = bnds(d++); int mxk = bnds(d++);
    int mnl = bnds(d++); int mxl = bnds(d++);
    int mnm = bnds(d++); int mxm = bnds(d++);
    d5_array n(mni,mxi,mnj,mxj,mnk,mxk,mnl,mxl,mnm,mxm);
    for (int i=mni;i<=mxi;i++){
        n(i) = sumOverLastDim(n_ijklmn(i),w);
    }
    return n;
}

/**
 * Compute marginal weighted sums over last dimension.
 * 
 * @param  n_ijklmno - d7_array
 * @param  w         - weighting vector
 * @return n_ijklmn  - d6_array
 */
d6_array tcsam::sumOverLastDim(const d7_array& n_ijklmno, dvector& w){
    ivector bnds = wts::getBounds(n_ijklmno);
    int d = 1;
    int mni = bnds(d++); int mxi = bnds(d++);
    int mnj = bnds(d++); int mxj = bnds(d++);
    int mnk = bnds(d++); int mxk = bnds(d++);
    int mnl = bnds(d++); int mxl = bnds(d++);
    int mnm = bnds(d++); int mxm = bnds(d++);
    int mno = bnds(d++); int mxo = bnds(d++);
    d6_array n(mni,mxi,mnj,mxj,mnk,mxk,mnl,mxl,mnm,mxm,mno,mxo);
    for (int i=mni;i<=mxi;i++){
        n(i) = sumOverLastDim(n_ijklmno(i),w);
    }
    return n;
}
    
/**
 * Calculate marginal sums over msz by iyx.
 * @param n_iyxmsz - array to calculate sums on
 * @return d3_array with dims iyx.
 */
d3_array tcsam::calcIYXfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d3_array n_iyx(mni,mxi,mny,mxy,mnx,mxx);
    n_iyx.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                n_iyx(i,y,x) = sum(n_iyxmsz(i,y,x));
            }
        }
    }
    return n_iyx;
}

/**
 * Calculate marginal sums over msz by iyx, weighted by w_xmz.
 * @param n_iyxmsz - array to calculate sums on
 * @param w_xmz - weight-at-xmz array
 * @return d3_array with dims iyx.
 */
d3_array tcsam::calcIYXfromIYXMSZ(d6_array& n_iyxmsz, d3_array w_xmz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d3_array b_iyx(mni,mxi,mny,mxy,mnx,mxx);
    b_iyx.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int m=mnm;m<=mxm;m++){
                    for (int s=mns;s<=mxs;s++) b_iyx(i,y,x) += n_iyxmsz(i,y,x,m,s)*w_xmz(x,m);//dot product here
                }
            }
        }
    }
    return b_iyx;
}

/**
 * Calculate marginal sums over z by iyxms, weighted by w_xmz.
 * @param n_iyxmsz - array to calculate sums on
 * @param w_xmz - weight-at-xmz array
 * @return d5_array with dims iyxms.
 */
d5_array tcsam::calcIYXMSfromIYXMSZ(d6_array& n_iyxmsz, d3_array w_xmz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int d = 1;
    int mni = bnds(d++); int mxi = bnds(d++);
    int mny = bnds(d++); int mxy = bnds(d++);
    int mnx = bnds(d++); int mxx = bnds(d++);
    int mnm = bnds(d++); int mxm = bnds(d++);
    int mns = bnds(d++); int mxs = bnds(d++);
    int mnz = bnds(d++); int mxz = bnds(d++);
    d5_array b_iyxms(mni,mxi,mny,mxy,mnx,mxx,mnm,mxm,mns,mxs);
    b_iyxms.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int m=mnm;m<=mxm;m++){
                    for (int s=mns;s<=mxs;s++) b_iyxms(i,y,x,m,s) = n_iyxmsz(i,y,x,m,s)*w_xmz(x,m);//dot product here
                }
            }
        }
    }
    return b_iyxms;
}

/**
 * Calculate marginal sums over z by yxms, weighted by w_xmz.
 * @param n_yxmsz - array to calculate sums on
 * @param w_xmz - weight-at-xmz array
 * @return d4_array with dims yxms.
 */
d4_array tcsam::calcYXMSfromYXMSZ(d5_array& n_yxmsz, d3_array w_xmz){
    ivector bnds = wts::getBounds(n_yxmsz);
    int d = 1;
    int mny = bnds(d++); int mxy = bnds(d++);
    int mnx = bnds(d++); int mxx = bnds(d++);
    int mnm = bnds(d++); int mxm = bnds(d++);
    int mns = bnds(d++); int mxs = bnds(d++);
    int mnz = bnds(d++); int mxz = bnds(d++);
    d4_array b_yxms(mny,mxy,mnx,mxx,mnm,mxm,mns,mxs);
    b_yxms.initialize();
    for (int x=mnx;x<=mxx;x++){
        for (int y=mny;y<=mxy;y++){
            for (int m=mnm;m<=mxm;m++){
                for (int s=mns;s<=mxs;s++) b_yxms(y,x,m,s) = n_yxmsz(y,x,m,s)*w_xmz(x,m);//dot product here
            }
        }
    }
    return b_yxms;
}

/**
 * Calculate marginal sums over msz by ixy.
 * @param n_iyxmsz - array to calculate sums on
 * @return d3_array with dims ixy.
 */
d3_array tcsam::calcIXYfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d3_array n_ixy(mni,mxi,mnx,mxx,mny,mxy);
    n_ixy.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                n_ixy(i,x,y) = sum(n_iyxmsz(i,y,x));
            }
        }
    }
    return n_ixy;
}

/**
 * Calculate marginal sums over msz by ixy, weighted by w_xmz.
 * @param n_iyxmsz - array to calculate sums on
 * @param w_xmz - weight-at-xmz array
 * @return d3_array with dims ixy.
 */
d3_array tcsam::calcIXYfromIYXMSZ(d6_array& n_iyxmsz, d3_array w_xmz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d3_array b_ixy(mni,mxi,mnx,mxx,mny,mxy);
    b_ixy.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int m=mnm;m<=mxm;m++){
                    for (int s=mns;s<=mxs;s++) b_ixy(i,x,y) += n_iyxmsz(i,y,x,m,s)*w_xmz(x,m);//dot product here
                }
            }
        }
    }
    return b_ixy;
}

/**
 * Calculate marginal sums over ms by ixyz.
 * @param n_iyxmsz - array to calculate sums on
 * @return d4_array with dims ixyz.
 */
d4_array tcsam::calcIXYZfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d4_array n_ixyz(mni,mxi,mnx,mxx,mny,mxy,mnz,mxz);
    n_ixyz.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int z=mnz;z<=mxz;z++) {
                    for (int m=mnm;m<=mxm;m++) {
                        for (int s=mns;s<=mxs;s++) n_ixyz(i,x,y,z) += n_iyxmsz(i,y,x,m,s,z);
                    }
                }
            }
        }
    }
    return n_ixyz;
}

/**
 * Calculate marginal sums over s by ixmyz.
 * @param n_iyxmsz - array to calculate sums on
 * @return d5_array with dims ixmyz.
 */
d5_array tcsam::calcIXMYZfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d5_array n_ixmyz(mni,mxi,mnx,mxx,mnm,mxm,mny,mxy,mnz,mxz);
    n_ixmyz.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int z=mnz;z<=mxz;z++) {
                    for (int m=mnm;m<=mxm;m++) {
                        for (int s=mns;s<=mxs;s++) n_ixmyz(i,x,m,y,z) += n_iyxmsz(i,y,x,m,s,z);
                    }
                }
            }
        }
    }
    return n_ixmyz;
}

/**
 * Sum iyxmsz array over maturity states and rearrange indices to ixsyz.
 * 
 * @param n_iyxmsz - d6_array w/ indices iyxmsz to sum
 * 
 * @return d5_array with indices ixsyz
 */
d5_array tcsam::calcIXSYZfromIYXMSZ(d6_array& n_iyxmsz){
    ivector bnds = wts::getBounds(n_iyxmsz);
    int mni = bnds( 1); int mxi = bnds( 2);
    int mny = bnds( 3); int mxy = bnds( 4);
    int mnx = bnds( 5); int mxx = bnds( 6);
    int mnm = bnds( 7); int mxm = bnds( 8);
    int mns = bnds( 9); int mxs = bnds(10);
    int mnz = bnds(11); int mxz = bnds(12);
    d5_array n_ixsyz(mni,mxi,mnx,mxx,mns,mxs,mny,mxy,mnz,mxz);
    n_ixsyz.initialize();
    for (int i=mni;i<=mxi;i++){
        for (int x=mnx;x<=mxx;x++){
            for (int y=mny;y<=mxy;y++){
                for (int z=mnz;z<=mxz;z++) {
                    for (int m=mnm;m<=mxm;m++) {
                        for (int s=mns;s<=mxs;s++) n_ixsyz(i,x,s,y,z) += n_iyxmsz(i,y,x,m,s,z);
                    }
                }
            }
        }
    }
    return n_ixsyz;
}

d5_array tcsam::rearrangeYXMSZtoXMSYZ(d5_array& n_yxmsz){
//    cout<<"In rearrangeYXMSZtoXMSYZ(...)"<<endl;
    ivector perm(1,5);//{4,1,2,3,5};
    perm[1]=4;perm[2]=1;perm[3]=2;perm[4]=3;perm[5]=5;
    d5_array n_xmsyz = wts::permuteDims(perm,n_yxmsz);
//    cout<<"Finished rearrangeYXMSZtoXMSYZ(...)"<<endl;
    return n_xmsyz;
}

/**
 * Rearrange indices for a d5_array from IYXMS to IXMSY.
 * 
 * @param n_iyxms - d5_array with indices iyxms
 * 
 * @return d5_array with indices ixmsy
 */
d5_array tcsam::rearrangeIYXMStoIXMSY(d5_array& n_iyxms){
//    cout<<"starting rearrangeIYXMStoIXMSY(...)"<<endl;
    ivector perm(1,5);//{1,5,2,3,4};
    perm[1]=1;perm[2]=5;perm[3]=2;perm[4]=3;perm[5]=4;
    d5_array n_ixmsy = wts::permuteDims(perm,n_iyxms);
//    cout<<"finished rearrangeIYXMStoIXMSY(...)"<<endl;
    return n_ixmsy;
}

/**
 * Rearrange indices for a d6_array from IYXMSZ to IXMSYZ.
 * 
 * @param n_iyxmsz - d6_array with indices iyxmsz
 * 
 * @return d6_array with indices ixmsyz
 */
d6_array tcsam::rearrangeIYXMSZtoIXMSYZ(d6_array& n_iyxmsz){
//    cout<<"starting rearrangeIYXMSZtoIXMSYZ(...)"<<endl;
    ivector perm(1,6);//{1,5,2,3,4,6};
    perm[1]=1;perm[2]=5;perm[3]=2;perm[4]=3;perm[5]=4;perm[6]=6;
    d6_array n_ixmsyz = wts::permuteDims(perm,n_iyxmsz);
//    cout<<"finished rearrangeIYXMSZtoIXMSYZ(...)"<<endl;
    return n_ixmsyz;
}

/**
 * Extract (possibly summary) vector from 5d array.
 * 
 * @param x - sex index             (tcsam::ALL_SXs yields sum over sex)
 * @param m - maturity index        (tcsam::ALL_MSs yields sum over maturity state)
 * @param s - shell condition index (tcsam::ALL_SCs yields sum over shell condition)
 * @param y - year index
 * @param n_xmsyz - d5_array from which to extract vector at z's
 * 
 * @return extracted vector (indices consistent with z's)
 */
dvector tcsam::extractFromXMSYZ(int x, int m, int s, int y, d5_array& n_xmsyz){
    ivector bnds = wts::getBounds(n_xmsyz);
    dvector n_z(bnds(9,10));//dimension for z index
    int xmn, xmx;
    xmn = xmx = x; if (x==tcsam::ALL_SXs) {xmn = 1; xmx = tcsam::nSXs;}
    int mmn, mmx;
    mmn = mmx = m; if (m==tcsam::ALL_MSs) {mmn = 1; mmx = tcsam::nMSs;}
    int smn, smx;
    smn = smx = s; if (s==tcsam::ALL_SCs) {smn = 1; smx = tcsam::nSCs;}
    n_z.initialize();
    for (int xp=xmn;xp<=xmx;xp++){
        for (int mp=mmn;mp<=mmx;mp++){
            for (int sp=smn;sp<=smx;sp++){
                n_z += n_xmsyz(xp,mp,sp,y);
            }
        }
    }
    return(n_z);
}

/**
 * Extract (possibly summary) value from 4d array w/ indices
 * y,x,m,s.
 * 
 * @param y - year index
 * @param x - sex index             (tcsam::ALL_SXs yields sum over sex)
 * @param m - maturity index        (tcsam::ALL_MSs yields sum over maturity state)
 * @param s - shell condition index (tcsam::ALL_SCs yields sum over shell condition)
 * @param n_yxms - d4_array from which to extract value
 * 
 * @return extracted double value
 */
double tcsam::extractFromYXMS(int y, int x, int m, int s, d4_array& n_yxms){
//    rpt::echo<<"in extractFromYXMS"<<endl;
//    rpt::echo<<y<<" "<<x<<" "<<m<<" "<<" "<<s<<endl;
    int xmn, xmx;
    xmn = xmx = x; if (x==tcsam::ALL_SXs) {xmn = 1; xmx = tcsam::nSXs;}
    int mmn, mmx;
    mmn = mmx = m; if (m==tcsam::ALL_MSs) {mmn = 1; mmx = tcsam::nMSs;}
    int smn, smx;
    smn = smx = s; if (s==tcsam::ALL_SCs) {smn = 1; smx = tcsam::nSCs;}
    double v = 0.0;
    for (int xp=xmn;xp<=xmx;xp++){
        for (int mp=mmn;mp<=mmx;mp++){
            for (int sp=smn;sp<=smx;sp++){
//                rpt::echo<<xp<<" "<<mp<<" "<<sp<<endl;
                v += n_yxms(y,xp,mp,sp);
            }
        }
    }
//    rpt::echo<<"finished extractFromYXMS"<<endl;
    return(v);
}

/**
 * Extract (possibly summary) vector from 5d array w/ indices
 * y,x,m,s,z.
 * 
 * @param y - year index
 * @param x - sex index             (tcsam::ALL_SXs yields sum over sex)
 * @param m - maturity index        (tcsam::ALL_MSs yields sum over maturity state)
 * @param s - shell condition index (tcsam::ALL_SCs yields sum over shell condition)
 * @param n_yxmsz - d5_array from which to extract vector at z's
 * 
 * @return extracted vector (indices consistent with z's)
 */
dvector tcsam::extractFromYXMSZ(int y, int x, int m, int s, d5_array& n_yxmsz){
    ivector bnds = wts::getBounds(n_yxmsz);
//    cout<<"in extractFromYXMSZ"<<endl;
    dvector n_z(bnds(9),bnds(10));//dimension for z index
//    cout<<n_z.indexmin()<<" "<<n_z.indexmax()<<endl;
//    cout<<y<<" "<<x<<" "<<m<<" "<<" "<<s<<endl;
    int xmn, xmx;
    xmn = xmx = x; if (x==tcsam::ALL_SXs) {xmn = 1; xmx = tcsam::nSXs;}
    int mmn, mmx;
    mmn = mmx = m; if (m==tcsam::ALL_MSs) {mmn = 1; mmx = tcsam::nMSs;}
    int smn, smx;
    smn = smx = s; if (s==tcsam::ALL_SCs) {smn = 1; smx = tcsam::nSCs;}
    n_z.initialize();
    for (int xp=xmn;xp<=xmx;xp++){
        for (int mp=mmn;mp<=mmx;mp++){
            for (int sp=smn;sp<=smx;sp++){
//                cout<<xp<<" "<<mp<<" "<<sp<<endl;
                n_z += n_yxmsz(y,xp,mp,sp);
            }
        }
    }
//    cout<<"finished extractFromYXMSZ"<<endl;
    return(n_z);
}
