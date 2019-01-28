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
    int nTACs;      //number of years for TACs/OFLs
    ivector yrsTAC; //years for TACs/OFLs
    dvector TAC_y;  //TACs, indexed by year
    dvector OFL_y;  //OFLs, indexed by year
    d3_array wAtZ_xmz;  //weight at size, by sex and maturity state
    dvector R_y;        //recruitment, by year
    dvector R_x;        //recruitment proportions, indexed by sex
    dvector R_z;        //recruitment proportions, indexed by size
    d4_array M_xmsz;    //natural mortality by sex, maturity state, shell condition, size
    d4_array prGr_xszz; //growth probabilities, by sex, shell condition, postmolt size, and premolt size
    dmatrix prM2M_xz;   //terminal olt probabilities, by sex and size
    dvector hmF_f;      //handling mortality, by fishery
    d5_array cpF_fxmsz; //capture rates, by fishery, sex, maturity state, shell condition and size
    d5_array ret_fxmsz; //retention rates, by fishery, sex, maturity state, shell condition and size
    d5_array sel_fxmsz; //selectivity rates, by fishery, sex, maturity state, shell condition and size
    d5_array q_vxmsz;   //survey catchability, by survey, sex, maturity state, shell condition and size
    d4_array n_xmsz;    //initial population abundance, by sex, maturity state, shell condition and size
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
    /**
     * Add one year to recruitment, TAC, and OFL time series.
     * @param prjRec
     * @param TAC
     * @param OFL
     */
    void addOneYear(double prjRec, double TAC, double OFL);
    friend cifstream& operator >>(cifstream & is, MSE_OpModInfo & obj){obj.read(is);  return is;}
    friend ostream&   operator <<(ostream & os,   MSE_OpModInfo & obj){obj.write(os); return os;}
protected:
    void allocate();
};

#endif /* MSECLASSES_HPP */

