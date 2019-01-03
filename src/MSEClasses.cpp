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
/**
 * Read from an input filestream in ADMB OpModMode format
 * @param is - input filestream
 */
void MSE_OpModInfo::read(cifstream& is){
    is>>mnYr; //start year
    is>>mxYr; //current year
    is>>dtF;  //dtF
    is>>dtM;  //dtM
    allocate();
    for (int x=1;x<=nSXs;x++)
        for (int m=1;m<=nMSs;m++)
            is>>wAtZ_xmz(x,m); 
    is>>R_y;
    is>>R_x;
    is>>R_z;
    
    for (int x=1;x<=nSXs;x++)
        for (int m=1;m<=nMSs;m++)
            for (int s=1;s<=nSCs;s++)
                is>>M_xmsz(x,m,s);
        
    for (int x=1;x<=nSXs;x++)
        for (int s=1;s<=nSCs;s++)
            for (int z=1;z<=nZBs;z++)
                is>>prGr_xszz(x,s,z);
        
    for (int x=1;x<=nSXs;x++)
        is>>prM2M_xz(x);
        
    is>>hmF_f;
    
    for (int f=1;f<=nFsh;f++)
        for (int x=1;x<=nSXs;x++)
            for (int m=1;m<=nMSs;m++)
                for (int s=1;s<=nSCs;s++)
                    is>>cpF_fxmsz(f,x,m,s);
            
    for (int f=1;f<=nFsh;f++)
        for (int x=1;x<=nSXs;x++)
            for (int m=1;m<=nMSs;m++)
                for (int s=1;s<=nSCs;s++)
                    is>>ret_fxmsz(f,x,m,s);
            
    for (int f=1;f<=nFsh;f++)
        for (int x=1;x<=nSXs;x++)
            for (int m=1;m<=nMSs;m++)
                for (int s=1;s<=nSCs;s++)
                    is>>sel_fxmsz(f,x,m,s);
            
    for (int v=1;v<=nSrv;v++)
        for (int x=1;x<=nSXs;x++)
            for (int m=1;m<=nMSs;m++)
                for (int s=1;s<=nSCs;s++)
                    is>>q_vxmsz(v,x,m,s);
            
    for (int x=1;x<=nSXs;x++)
        for (int m=1;m<=nMSs;m++)
            for (int s=1;s<=nSCs;s++)
                is>>n_xmsz(x,m,s);
}

/**
 * Write to an output filestream in ADMB OpModMode format
 * @param os - output filestream
 */
void MSE_OpModInfo::write(ostream& os){
    os<<mnYr<<tb<<"#mnYr"<<endl; //start year
    os<<mxYr<<tb<<"#mxYr"<<endl; //current year
    os<<dtF<<tb<<"#dtF"<<endl;  //dtF
    os<<dtM<<tb<<"#dtM"<<endl;  //dtM
    os<<endl;
    
    os<<"#wAtZ_xmz:"<<endl; wts::print(wAtZ_xmz,os,0); os<<endl;
    
    os<<"#R_y:"<<endl; os<<R_y<<endl;
    os<<"#R_x:"<<endl; os<<R_x<<endl;
    os<<"#R_z:"<<endl; os<<R_z<<endl;
    os<<endl;
    
    os<<"#M_xmsz:"   <<endl; wts::print(M_xmsz,os,0);    os<<endl;
    os<<"#prGr_xszz:"<<endl; wts::print(prGr_xszz,os,0); os<<endl;
    os<<"#prM2M_xz:" <<endl; wts::print(prM2M_xz,os,0);  os<<endl;
    
    os<<"#hmF_f:"<<endl; os<<hmF_f<<endl;
    os<<endl;
    
    os<<"#cpF_fxmsz:"<<endl; wts::print(cpF_fxmsz,os,0); os<<endl;
    os<<"#ret_fxmsz:"<<endl; wts::print(ret_fxmsz,os,0); os<<endl;
    os<<"#sel_fxmsz:"<<endl; wts::print(sel_fxmsz,os,0); os<<endl;
    os<<"#q_vxmsz:"  <<endl; wts::print(q_vxmsz,os,0);   os<<endl;
    os<<"#n_xmsz:"   <<endl; wts::print(n_xmsz,os,0);    os<<endl;
}

/**
 * Write to an output filestream in R OpModMode format
 * @param os - output filestream
 */
void MSE_OpModInfo::writeToR(ostream& os){
    os<<mnYr<<tb<<"#mnYr"<<endl; //start year
    os<<mxYr<<tb<<"#mxYr"<<endl; //current year
    os<<dtF<<tb<<"#dtF"<<endl;  //dtF
    os<<dtM<<tb<<"#dtM"<<endl;  //dtM
    os<<"#--wAtZ_xmz:"<<endl;
    wts::print(wAtZ_xmz,os,0); 
    os<<"#--R_y:"<<endl;
    os<<R_y<<endl;
    os<<"#--R_x:"<<endl;
    os<<R_x<<endl;
    os<<"#--R_z:"<<endl;
    os<<R_z<<endl;
    os<<"#--M_xmsz:"<<endl;
    wts::print(M_xmsz,os,0);
    os<<"#--prGr_xszz:"<<endl;
    wts::print(prGr_xszz,os,0);
    os<<"#--prM2M_xz:"<<endl;
    wts::print(prM2M_xz,os,0);
    os<<"#--hmF_f:"<<endl;
    os<<hmF_f<<endl;
    os<<"#--cpF_fxmsz:"<<endl;
    wts::print(cpF_fxmsz,os,0);
    os<<"#--ret_fxmsz:"<<endl;
    wts::print(ret_fxmsz,os,0);
    os<<"#--sel_fxmsz:"<<endl;
    wts::print(sel_fxmsz,os,0);
    os<<"#--q_vxmsz:"<<endl;
    wts::print(q_vxmsz,os,0);
    os<<"#--n_xmsz:"<<endl;
    wts::print(n_xmsz,os,0);
}

void MSE_OpModInfo::allocate(){
    wAtZ_xmz.allocate(1,nSXs,1,nMSs,1,nZBs);
    R_y.allocate(mnYr,mxYr-1);
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