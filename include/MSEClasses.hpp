/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MSEClasses.hpp
 * Author: WilliamStockhausen
 *
 * Created on November 19, 2018, 9:22 AM
 */

#ifndef MSECLASSES_HPP
#define MSECLASSES_HPP

class MSE_OpModInfo {
public:
    static int debug;
public:
    ModelConfiguration* ptrMC;
    int nSXs;
    int nMSs;
    int nSCs;
    int nZBs;
    int nFsh;
    int nSrv;
    int mnYr; //start year
    int mxYr; //current year
    double dtF;  //dtF
    double dtM;  //dtM
    d3_array wAtZ_xmz; 
    dvector R_y;
    dvector R_x;
    dvector R_z;
    d4_array M_xmsz;
    d4_array prGr_xszz;
    dmatrix prM2M_xz;
    dvector hmF_f;
    d5_array cpF_fxmsz;
    d5_array ret_fxmsz;
    d5_array sel_fxmsz;
    d5_array q_vxmsz;
    d4_array n_xmsz;
public:
    MSE_OpModInfo(ModelConfiguration* ptrMC);
    ~MSE_OpModInfo(){ptrMC=0;}
    /**
     * Read from an input filestream in ADMB OpModMode format
     * @param is - input filestream
     */
    void read(cifstream& is);
    /**
     * Write to an output filestream in ADMB OpModMode format
     * @param is - output filestream
     */
    void write(ostream& os);
    /**
     * Write to an output filestream in R format
     * @param is - output filestream
     */
    void writeToR(ostream& os);
    friend cifstream& operator >>(cifstream & is, MSE_OpModInfo & obj){obj.read(is);  return is;}
    friend ostream&   operator <<(ostream & os,   MSE_OpModInfo & obj){obj.write(os); return os;}
protected:
    void allocate();
};

#endif /* MSECLASSES_HPP */

