/*
 * MSEClasses.cpp
 */

#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelConfiguration.hpp"
#include "MSEClasses.hpp"

int MSE_OpModInfo::debug = 0;

MSE_OpModInfo::MSE_OpModInfo(ModelConfiguration* ptrMC){
    nSXs = tcsam::nSXs;
    nMSs = tcsam::nMSs;
    nSCs = tcsam::nSCs;
    nZBs = ptrMC->nZBs;
    nFsh = ptrMC->nFsh;
    nSrv = ptrMC->nSrv;
}
void MSE_OpModInfo::read(cifstream& is){
    is>>mnYr; //start year
    is>>mxYr; //current year
    is>>dtF;  //dtF
    is>>dtM;  //dtM
    allocate();
    is>>wAtZ_xmz; 
    is>>R_y(mnYr,mxYr);
    is>>R_x;
    is>>R_z;
    is>>M_xmsz;
    is>>prGr_xszz;
    is>>prM2M_xz;
    is>>hmF_f;
    is>>cpF_fxmsz;
    is>>ret_fxmsz;
    is>>sel_fxmsz;
    is>>q_vxmsz;
    is>>n_xmsz;
}

void MSE_OpModInfo::allocate(){
    wAtZ_xmz.allocate(1,nSXs,1,nMSs,1,nZBs);
    R_y.allocate(mnYr,mxYr);
    R_x.allocate(1,nSXs);
    R_z.allocate(1,nZBs);
    M_xmsz.allocate(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    prGr_xszz.allocate(1,nSXs,1,nSCs,1,nZBs,1,nZBs);
    prM2M_xz.allocate(1,nSXs,1,nZBs);
    hmF_f.allocate(1,nFsh);
    cpF_fxmsz.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    ret_fxmsz.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    sel_fxmsz.allocate(1,nFsh,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    q_vxmsz.allocate(1,nSrv,1,nSXs,1,nMSs,1,nSCs,1,nZBs);
    n_xmsz.allocate(1,nSXs,1,nMSs,1,nSCs,1,nZBs);
}