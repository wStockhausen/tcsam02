//ModelParameterFunctions.cpp

#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelParameterInfoTypes.hpp"
#include "ModelParameterFunctions.hpp"



namespace tcsam{

/******************************************************************************
* Set initial values for a param_init_bounded_number_vector.
* 
* Description: Sets initial values for a parameter vector.
* 
* Inputs:
*  @param pI : pointer to BoundedNumberVectorInfo object
*  @param p  : reference to a param_init_bounded_number_vector
*     
* @return void
 * 
* @alters - changes initial values of p
 * 
******************************************************************************/
void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, std::ostream& cout){
    debug=tcsam::dbgAll;
    if (debug>=tcsam::dbgAll) std::cout<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitVals();
        p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).get_name()<<": "<<p<<std::endl;
        if (debug>=tcsam::dbgAll) {
            std::cout<<"vls = "<<vls<<std::endl;
            std::cout<<"p   = "<<p<<std::endl;
            double_index_guts* ptr = new dvector_index(vls);
            std::cout<<ptr<<tb<<ptr->indexmin()<<tb<<ptr->indexmax()<<endl;
            for (int i=ptr->indexmin();i<=ptr->indexmax();i++) std::cout<<(*((*ptr)[i]))<<tb; std::cout<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=tcsam::dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

/******************************************************************************
* Set initial values for a param_init_bounded_vector_vector.
* 
* Description: Sets initial values for a vector of parameter vectors.
* 
* Inputs:
*  @param pI : pointer to a BoundedVectorVectorInfo object
*  @param p  : reference to a param_init_bounded_vector_vector
*     
* @return void
 * 
* @alters - changes initial values of p
 * 
******************************************************************************/
void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout){
    debug=tcsam::dbgAll;
    if (debug>=tcsam::dbgAll) std::cout<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            dvector vls = (*pI)[i]->getInitVals();
            if (debug>=tcsam::dbgAll) cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<std::endl;
            for (int j=vls.indexmin();j<=vls.indexmax();j++) p(i,j) = vls(j);
            if (debug>=tcsam::dbgAll) {
                std::cout<<"vls  = "<<vls<<std::endl;
                std::cout<<"p(i) = "<<p(i)<<std::endl;
            }
        }
        for (int i=1;i<=np;i++) rpt::echo<<"InitVals for "<<p(i).get_name()<<":"<<tb<<p(i)<<std::endl;
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=tcsam::dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

/******************************************************************************
* Set initial devs values for a param_init_bounded_vector_vector.
* 
* Description: Sets initial values for a vector of devs vectors.
* 
* Inputs:
*  @param pI : pointer to a DevsVectorVectorInfo object
*  @param p  : reference to a param_init_bounded_vector_vector
*     
* @return void
 * 
* @alters - changes initial values of p
******************************************************************************/
void setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& cout){
    debug=tcsam::dbgAll;
    if (debug>=tcsam::dbgAll) std::cout<<"Starting setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            dvector vls = (*pI)[i]->getInitVals();
            if (debug>=tcsam::dbgAll) cout<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<std::endl;
            for (int j=vls.indexmin();j<=(vls.indexmax()-1);j++) p(i,j)=vls(j);
            if (debug>=tcsam::dbgAll) {
                std::cout<<"vls  = "<<vls<<std::endl;
                std::cout<<"p(i) = "<<p(i)<<std::endl;
            }
        }
        for (int i=1;i<=np;i++) rpt::echo<<"InitVals for "<<p(i).get_name()<<":"<<tb<<p(i)<<std::endl;
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=tcsam::dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(DevsVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

/**
 * Set values for a dvar_matrix, intended to function as a vector of devs vector, 
 * based on values of a param_init_bounded_vector_vector.
 * 
 * @param devs  - the dvar_matrix
 * @param pDevs - the param_init_bounded_vector_vector
 * @param pI    - pointer to the associated DevsVectorVectorInfo object
 * @param debug - debugging level
 * @param cout  - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - This alters the values in the devs dvar_matrix.
 * 
 */
void setDevs(dvar_matrix& devs, param_init_bounded_vector_vector& pDevs, DevsVectorVectorInfo* pI, int debug, std::ostream& cout){
    if (debug>=tcsam::dbgAll) cout<<"starting setDevs(devs,pDevs,pI)"<<std::endl;
    devs.initialize();
    if (pI->getSize()){
        int nv = pDevs.indexmax();//number of devs vectors defined
        int mni; int mxi;
        for (int v=1;v<=nv;v++){
            mni = pDevs(v).indexmin();
            mxi = pDevs(v).indexmax();
            devs(v)(mni,mxi) = (*pI)[v]->calcArithScaleVal(pDevs(v));
            if (debug>=(tcsam::dbgAll)) cout<<v<<":  "<<devs(v)<<endl;
        }
        if (debug>=(tcsam::dbgAll+1)) {
            std::cout<<"Enter 1 to continue >>";
            std::cin>>nv;
            if (nv<0) exit(-1);
            cout<<"finished setDevs(devs,pDevs,pI)"<<std::endl;
        }
    } else {
        if (debug>=(tcsam::dbgAll)) cout<<"size=0, so no devs to set"<<endl;
    }
}

/**
 * Calculate ln-scale priors for a param_init_number_vector, based on its associated NumberVectorInfo, 
 * and add the weighted NLL to the objective function.
 * 
 * @param objFun - the objective function
 * @param ptrVI  - pointer to the NumberVectorInfo object associated with pv
 * @param pv     - the param_init_number_vector
 * @param debug  - debugging level
 * @param cout   - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - the value of objFun
 */                                      
void calcPriors(objective_function_value& objFun, NumberVectorInfo* ptrVI,param_init_number_vector& pv, int debug, std::ostream& cout){
    if (ptrVI->getSize()){
        dvar_vector tmp = pv*1.0;
        dvar_vector pri = ptrVI->calcLogPriors(tmp);//ln-scale prior (NOT NLL!)
        if (debug>=tcsam::dbgPriors) cout<<"priors("<<ptrVI->name<<") = "<<pri<<std::endl;
        objFun += -ptrVI->getPriorWgts()*pri;
        if (debug<0){
            dvector wts = ptrVI->getPriorWgts();
            cout<<ptrVI->name<<"=list("<<endl;
            for (int i=pri.indexmin();i<pri.indexmax();i++){
                adstring type = (*ptrVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-pri(i)<<cc<<"objfun="<<-wts(i)*pri(i)<<"),"<<endl;
            }
            {
                int i=pri.indexmax();
                adstring type = (*ptrVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-pri(i)<<cc<<"objfun="<<-wts(i)*pri(i)<<")"<<endl;
            }
            cout<<")";
        }
    } else {
        if (debug<0) cout<<ptrVI->name<<"=NULL";
    }
}

/**
 * Calculate ln-scale priors for a param_init_bounded_number_vector, 
 * based on its associated BoundedNumberVectorInfo, 
 * and add the weighted NLL to the objective function.
 * 
 * @param objFun - the objective function
 * @param ptrVI  - pointer to the BoundedNumberVectorInfo object associated with pv
 * @param pv     - the param_init_bounded_number_vector
 * @param debug  - debugging level
 * @param cout   - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - the value of objFun
 */                                      
void calcPriors(objective_function_value& objFun, BoundedNumberVectorInfo* ptrVI,param_init_bounded_number_vector& pv, int debug, std::ostream& cout){
    if (ptrVI->getSize()){
        dvar_vector tmp = pv*1.0;
        dvar_vector pri = ptrVI->calcLogPriors(tmp);//ln-scale prior (NOT NLL!)
        if (isnan(value(sum(pri)))){
            std::cout<<"Found NAN for priors("<<ptrVI->name<<") = "<<pri<<std::endl;
            std::cout<<"param values = "<<pv<<endl;
            ptrVI->write(std::cout);
            std::cout<<"Aborting!!"<<endl;
            exit(-1);
        }
        if (debug>=tcsam::dbgPriors) cout<<"priors("<<ptrVI->name<<") = "<<pri<<std::endl;
        objFun += -ptrVI->getPriorWgts()*pri;
        if (debug<0){
            dvector wts = ptrVI->getPriorWgts();
            cout<<ptrVI->name<<"=list("<<endl;
            for (int i=pri.indexmin();i<pri.indexmax();i++){
                adstring type = (*ptrVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-pri(i)<<cc<<"objfun="<<-wts(i)*pri(i)<<"),"<<endl;
            }
            {
                int i=pri.indexmax();
                adstring type = (*ptrVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-pri(i)<<cc<<"objfun="<<-wts(i)*pri(i)<<")"<<endl;
            }
            cout<<")";
        }
    } else {
        if (debug<0) cout<<ptrVI->name<<"=NULL";
    }
}

/**
 * Calculate ln-scale priors for a param_init_bounded_vector_vector, 
 * based on its associated BoundedVectorVectorInfo, 
 * and add the weighted NLL to the objective function.
 * 
 * @param objFun - the objective function
 * @param ptrVI  - pointer to the BoundedVectorVectorInfo object associated with pv
 * @param pv     - the param_init_bounded_vector_vector
 * @param debug  - debugging level
 * @param cout   - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - the value of objFun
 */                                      
void calcPriors(objective_function_value& objFun, BoundedVectorVectorInfo* ptrVVI, param_init_bounded_vector_vector& pm, int debug, std::ostream& cout){
    if (ptrVVI->getSize()){
        if (debug<0) cout<<ptrVVI->name<<"=list("<<endl;        
        dvector wts = ptrVVI->getPriorWgts();
        for (int i=pm.indexmin();i<pm.indexmax();i++) {
            dvar_vector tmp = 1.0*pm(i);
            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
            if (isnan(value(sum(pri)))){
                std::cout<<"Found NAN for priors("<<ptrVVI->name<<"["<<i<<"]) = "<<pri<<std::endl;
                std::cout<<"param values = "<<tmp<<endl;
                ptrVVI->write(std::cout);
                std::cout<<"Aborting!!"<<endl;
                exit(-1);
            }
            if (debug>=tcsam::dbgPriors){
                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
            }
            objFun += -wts(i)*sum(pri);
            if (debug<0) {
                adstring type = (*ptrVVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-sum(pri)<<cc<<"objfun="<<-wts(i)*sum(pri)<<"),"<<endl;
            }
        }
        {
            int i = pm.indexmax();
            dvar_vector tmp = 1.0*pm(i);
            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
            if (isnan(value(sum(pri)))){
                std::cout<<"Found NAN for priors("<<ptrVVI->name<<"["<<i<<"]) = "<<pri<<std::endl;
                std::cout<<"param values = "<<tmp<<endl;
                ptrVVI->write(std::cout);
                std::cout<<"Aborting!!"<<endl;
                exit(-1);
            }
            if (debug>=tcsam::dbgPriors){
                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
            }
            objFun += -wts(i)*sum(pri);
            if (debug<0) {
                adstring type = (*ptrVVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-sum(pri)<<cc<<"objfun="<<-wts(i)*sum(pri)<<")"<<endl;
            }
        }
        if (debug<0) cout<<")";
    } else {
        if (debug<0) cout<<ptrVVI->name<<"=NULL";
    }
}

/**
 * Calculate ln-scale priors for a dvar_matrix acting as a devs_vector_vector (not defined yet), 
 * based on its associated DevsVectorVectorInfo, 
 * and add the weighted NLL to the objective function.
 * 
 * @param objFun - the objective function
 * @param ptrVI  - pointer to the DevsVectorVectorInfo object associated with pv
 * @param pv     - the dvar_matrix
 * @param debug  - debugging level
 * @param cout   - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - the value of objFun
 */                                      
void calcPriors(objective_function_value& objFun, DevsVectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout){
    if (ptrVVI->getSize()){
        if (debug<0) cout<<ptrVVI->name<<"=list("<<endl;        
        dvector wts = ptrVVI->getPriorWgts();
        for (int i=1;i<pm.indexmax();i++) {
            dvar_vector tmp = 1.0*pm(i);
            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
            if (isnan(value(sum(pri)))){
                std::cout<<"Found NAN for priors("<<ptrVVI->name<<"["<<i<<"]) = "<<pri<<std::endl;
                std::cout<<"param values = "<<tmp<<endl;
                ptrVVI->write(std::cout);
                std::cout<<"Aborting!!"<<endl;
                exit(-1);
            }
            if (debug>=tcsam::dbgPriors){
                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
            }
            objFun += -wts(i)*sum(pri);
            if (debug<0) {
                adstring type = (*ptrVVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-sum(pri)<<cc<<"objfun="<<-wts(i)*sum(pri)<<"),"<<endl;
            }
        }
        {
            int i = pm.indexmax();
            dvar_vector tmp = 1.0*pm(i);
            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
             if (isnan(value(sum(pri)))){
                std::cout<<"Found NAN for priors("<<ptrVVI->name<<"["<<i<<"]) = "<<pri<<std::endl;
                std::cout<<"param values = "<<tmp<<endl;
                ptrVVI->write(std::cout);
                std::cout<<"Aborting!!"<<endl;
                exit(-1);
            }
           if (debug>=tcsam::dbgPriors){
                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
            }
            objFun += -wts(i)*sum(pri);
            if (debug<0) {
                adstring type = (*ptrVVI)[i]->getPriorType();
                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-sum(pri)<<cc<<"objfun="<<-wts(i)*sum(pri)<<")"<<endl;
            }
        }
        if (debug<0) cout<<")";
    } else {
        if (debug<0) cout<<ptrVVI->name<<"=NULL";
    }
}

/**
 * Function to write parameter information to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter (param_init_number)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated NumberInfo object
 * @param toR - flag to write to R format
 * @param willBeActive - flag to write only if parameter will be active
 */
void writeParameter(ostream& os, param_init_number& p, adstring& ctg1, adstring& ctg2, 
                   NumberInfo* pI,int toR, int willBeActive){
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_number'"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"value="<<pI->calcArithScaleVal(value(p))
                                <<"),";
        } else {
            os<<1<<cc<<p.get_phase_start()<<cc<<1<<cc<<1<<cc<<"-Inf"<<cc<<"Inf"<<cc
              <<pI->calcArithScaleVal(p)<<cc<<p.get_name()<<cc<<"\"param_init_number\",\""<<ctg1<<"\",\""<<ctg2<<"\",\""
              <<pI->label<<"\""<<endl;
        }
    }
}    
/**
 * Function to write parameter information to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter (param_init_bounded_number)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - BoundedNumberInfo
 * @param toR - flag to write to R format
 * @param willBeActive - flag to write only if parameter will be active
 */
void writeParameter(ostream& os, param_init_bounded_number& p,adstring& ctg1, adstring& ctg2, 
                    BoundedNumberInfo* pI, int toR, int willBeActive){
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        std::cout<<p.get_name()<<tb<<p.get_phase_start()<<tb<<p.get_minb()<<tb<<p.get_maxb()<<tb<<p<<endl;
        std::cout<<pI->label<<tb<<pI->getLowerBound()<<tb<<pI->getUpperBound()<<tb<<pI->getInitVal()<<tb<<pI->calcParamScaleVal(pI->getInitVal())<<endl;
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_bounded_number'"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"phase="<<p.phase_start<<cc
                                <<"bounds=c("<<pI->calcArithScaleVal(p.get_maxb())<<cc
                                             <<pI->calcArithScaleVal(p.get_minb())<<")"<<cc
                                <<"value="<<pI->calcArithScaleVal(p)
                                <<"),";
        } else {
            os<<1<<cc<<p.get_phase_start()<<cc<<1<<cc<<1<<cc
              <<pI->calcArithScaleVal(p.get_minb())<<cc<<pI->calcArithScaleVal(p.get_maxb())<<cc
              <<pI->calcArithScaleVal(p)<<cc<<p.get_name()<<cc<<"\"param_init_bounded_number\",\""<<ctg1<<"\",\""<<ctg2<<"\",\""
              <<pI->label<<"\""<<endl;;
        }
    }
}    
/**
 * Function to write parameter vector information to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter vector (param_init_vector)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param lbl - parameter-specific label
 * @param toR - flag to write to R format
 * @param willBeActive - flag to write only if parameter vector will be active
 */
void writeParameter(ostream& os, param_init_vector& p,adstring& ctg1, adstring& ctg2, 
                    VectorInfo* pI, int toR, int willBeActive){
    int mn = p.indexmin();
    int mx = p.indexmax();
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        dvector vals = pI->calcArithScaleVal(value(p));
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_vector'"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"value=c(";for (int i=mn;i<mx;i++) {os<<vals(i)<<cc;} os<<vals(mx)<<")";
            os<<"),";
        } else {        
            for (int i=mn;i<=mx;i++) os<< i<<cc<<p.get_phase_start()<<cc
                    <<mn<<cc<<mx<<cc<<"-Inf"<<cc<<"Inf"<<cc
                    <<vals(i)<<cc<<p.get_name()<<cc<<"\"param_init_vector\",\""<<ctg1<<"\",\""<<ctg2<<"\",\""
                    <<pI->label<<"\""<<endl;
        }
    }
}       

/**
 * Function to write parameter vector information to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter vector (param_init_bounded_vector)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedVectorInfo info object
 * @param toR - flag to write to R format
 * @param willBeActive - flag to write only if parameter vector will be active
 */
void writeParameter(ostream& os, param_init_bounded_vector& p,adstring& ctg1, adstring& ctg2, 
                    BoundedVectorInfo* pI, int toR, int willBeActive){
    int mn = p.indexmin();
    int mx = p.indexmax();
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        dvector vals = pI->calcArithScaleVal(value(p));
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_bounded_vector'"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"bounds=c("<<pI->calcArithScaleVal(p.get_minb())<<cc
                                             <<pI->calcArithScaleVal(p.get_maxb())<<")"<<cc
                                <<"value=c("; for (int i=mn;i<mx;i++) {os<<vals(i)<<cc;} os<<vals(mx)<<")";
            os<<"),";
        } else {
            for (int i=mn;i<=mx;i++) os<<i<<cc<<p.get_phase_start()<<cc<<mn<<cc<<mx<<cc
                                       <<pI->calcArithScaleVal(p.get_minb())<<cc
                                       <<pI->calcArithScaleVal(p.get_maxb())<<cc
                                       <<vals(i)<<cc<<p.get_name()<<cc<<"\"param_init_bounded_vector\",\""
                                       <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
        }
    }
}

/**
 * Function to write parameter vector information to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter vector (param_init_bounded_dev_vector)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated DevsVectorInfo info object
 * @param toR - flag to write to R format
 * @param willBeActive - flag to write only if parameter vector will be active
 */
void writeParameter(ostream& os, param_init_bounded_dev_vector& p,adstring& ctg1, adstring& ctg2, 
                    DevsVectorInfo* pI, int toR, int willBeActive){
    int mn = p.indexmin();
    int mx = p.indexmax();
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        if (toR){
            os<<p.get_name()<<"=list("<<"type=param_init_bounded_dev_vector"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"bounds=c("<<pI->calcArithScaleVal(p.get_minb())<<cc
                                             <<pI->calcArithScaleVal(p.get_maxb())<<")"<<cc
                                <<"value=c("; 
                for (int i=mn;i<mx;i++) os<<pI->calcArithScaleVal(value(p(i)))<<cc; 
                os<<pI->calcArithScaleVal(value(p(mx)))<<")";
           os<<"),";
        } else {
            for (int i=mn;i<=mx;i++) os<< i<<cc<<p.get_phase_start()<<cc<<mn<<cc<<mx<<cc
                                       <<pI->calcArithScaleVal(p.get_minb())<<cc
                                       <<pI->calcArithScaleVal(p.get_maxb())<<cc
                                       <<pI->calcArithScaleVal(value(p(i)))<<cc<<p.get_name()<<cc
                                       <<"\"param_init_bounded_dev_vector\",\""<<ctg1<<"\",\""<<ctg2<<"\",\""
                                       <<pI->label<<"\""<<endl;
        }
    }
}    

/**
 * Function to write a vector of parameters (param_init_number_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_number_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated NumberVectorInfo info object
 * @param toR - flag to write to R format (otherwise csv)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
void writeParameters(ostream& os, param_init_number_vector& p, adstring& ctg1, adstring& ctg2, 
                     NumberVectorInfo* pI,int toR, int willBeActive){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toR,willBeActive);
        }
    }
}

/**
 * Function to write a vector of parameters (param_init_bounded_number_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_bounded_number_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedNumberVectorInfo info object
 * @param toR - flag to write to R format (otherwise csv)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
void writeParameters(ostream& os, param_init_bounded_number_vector& p, adstring& ctg1, adstring& ctg2, 
                     BoundedNumberVectorInfo* pI,int toR, int willBeActive){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
//            rpt::echo<<"writeParameters(BNVI) for "<<i<<p[i].get_name()<<endl;
//            rpt::echo<<"typeof((*pI)[i]) = "<<typeid((*pI)[i]).name()<<endl;
//            BoundedNumberInfo* ptrI = (*pI)[i];
//            rpt::echo<<"typeof(ptrI) = "<<typeid(ptrI).name()<<endl;
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toR,willBeActive);
        }
    }
}

/**
 * Function to write a vector of parameters (param_init_vector_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_vector_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated VectorVectorInfo info object
 * @param toR - flag to write to R format (otherwise csv)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
void writeParameters(ostream& os, param_init_vector_vector& p, adstring& ctg1, adstring& ctg2, 
                     VectorVectorInfo* pI,int toR, int willBeActive){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toR,willBeActive);
        }
    }
}

/**
 * Function to write a vector of parameters (param_init_bounded_vector_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_bounded_vector_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedVectorVectorInfo info object
 * @param toR - flag to write to R format (otherwise csv)
 * @param willBeActive - flag to write if parameters will be active in some phase
 */
void writeParameters(ostream& os, param_init_bounded_vector_vector& p, adstring& ctg1, adstring& ctg2, 
                     BoundedVectorVectorInfo* pI,int toR, int willBeActive){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toR,willBeActive);
        }
    }
}

/**
 * Sets the parameter info for a NumberVectorInfo object.
 * 
 * 
 * @param [in]      pNVI - pointer to a NumberVectorInfo instance
 * @param [in][out] npT - size of vector
 * @param [in][out] phs - ivector of phases for parameters
 * @param [in]      os  - output stream to write processing info to
 */
void setParameterInfo(NumberVectorInfo* pNVI,                           
                             int& npT,
                             ivector& phs, 
                             ostream& os){
    int np = pNVI->getSize();
    if (np){npT = np;} else {npT = 1;}
    phs.allocate(1,npT);
    if (np){
        phs = pNVI->getPhases();
        os<<"parameter "<<pNVI->name<<":"<<endl;
        os<<"#phase"<<endl;
        for (int n=1;n<=np;n++) os<<n<<tb<<phs(n)<<endl;
    } else {
        phs = -1;
        os<<"number vector parameter "<<pNVI->name<<" has no parameter values"<<endl;
    }
}
  
/**
 * Set the parameter info for a BoundedNumberVectorInfo object.
 * 
 * @param [in]      pBNVI - pointer to a BoundedNumberVectorInfo instance
 * @param [in][out] npT   - size of vector
 * @param [in][out] lb    - dvector of lower bounds on parameter scale
 * @param [in][out] ub    - dvector of upper bounds on parameter scale
 * @param [in][out] phs   - ivector of phases for parameters
 * @param [in]      os    - output stream to write to
 */
void setParameterInfo(BoundedNumberVectorInfo* pBNVI,
                             int& npT,
                             dvector& lb, dvector& ub, 
                             ivector& phs, 
                             ostream& os){
    int np = pBNVI->getSize();
    if (np){npT = np;} else {npT = 1;}
    phs.allocate(1,npT);
    lb.allocate(1,npT);
    ub.allocate(1,npT);
    adstring_array scales(1,npT);
    if (np){
        phs = pBNVI->getPhases();
        lb  = pBNVI->getLowerBounds();
        ub  = pBNVI->getUpperBounds();
        scales = pBNVI->getScaleTypes();
        os<<"parameter "<<pBNVI->name<<":"<<endl;
        os<<"#lower  upper  phase  scale"<<endl;
        for (int n=1;n<=np;n++) os<<n<<tb<<lb(n)<<tb<<ub(n)<<tb<<phs(n)<<tb<<scales(n)<<endl;
    } else {
        phs = -1;
        lb  = -1.0;
        ub  =  1.0;
        os<<"bounded number parameter vector "<<pBNVI->name<<" has no parameter values"<<endl;
    }
}
  
/**
 * Set the parameter info for a VectorVectorInfo object.
 * 
 * @param pVVI - pointer to a VectorVectorInfo instance
 * @param npT - size of vector
 * @param mns - ivector with minimum indices for each vector
 * @param mxs - ivector with maximum indices for each vector
 * @param phs - ivector of phases for parameters
 * @param os - output stream to write to
 */
void setParameterInfo(VectorVectorInfo* pVVI,
                             int& npT,
                             ivector& mns, ivector& mxs,
                             ivector& phs, 
                             ostream& os){
    int np = pVVI->getSize();
    if (np){npT = np;} else {npT = 1;}
    mns.allocate(1,npT);
    mxs.allocate(1,npT);
    phs.allocate(1,npT);
    if (np){
        phs = pVVI->getPhases();
        mns = pVVI->getMinIndices();
        mxs = pVVI->getMaxIndices();
        os<<"parameter vector"<<pVVI->name<<":"<<endl;
        os<<"#mnIdx  mxIdx  phase"<<endl;
        for (int n=1;n<=np;n++) os<<n<<tb<<mns(n)<<tb<<mxs(n)<<tb<<phs(n)<<endl;
    } else {
        mns =  0;
        mxs =  0;
        phs = -1;
        os<<"vector vector "<<pVVI->name<<" has no parameter values"<<endl;
    }
}
  
/**
 * Set the parameter info for a BoundedVectorVectorInfo object.
 * 
 * @param pBVVI - pointer to a BoundedVectorVectorInfo instance
 * @param npT - size of vector
 * @param mns - ivector with minimum indices for each vector
 * @param mxs - ivector with maximum indices for each vector
 * @param idxs - imatrix of reverse indices
 * @param lb - dvector of lower bounds
 * @param ub - dvector of upper bounds
 * @param phs - ivector of phases for parameters
 * @param os - output stream to write to
 */
void setParameterInfo(BoundedVectorVectorInfo* pBVVI,                           
                             int& npT,
                             ivector& mns, ivector& mxs,
                             imatrix& idxs,
                             dvector& lb, dvector& ub, 
                             ivector& phs, 
                             ostream& os){
    int np = pBVVI->getSize();
    if (np){npT = np;} else {npT = 1;}
    mns.allocate(1,npT);
    mxs.allocate(1,npT);
    lb.allocate(1,npT);
    ub.allocate(1,npT);
    phs.allocate(1,npT);
    if (np){
        phs = pBVVI->getPhases();
        mns = pBVVI->getMinIndices();
        mxs = pBVVI->getMaxIndices();
        lb  = pBVVI->getLowerBounds();
        ub  = pBVVI->getUpperBounds();
        os<<"parameter vector "<<pBVVI->name<<":"<<endl;
        os<<"#mnIdx  mxIdx  lower  upper  phase"<<endl;
        for (int n=1;n<=np;n++) os<<n<<tb<<mns(n)<<tb<<mxs(n)<<tb<<lb(n)<<tb<<ub(n)<<tb<<phs(n)<<endl;
        idxs.allocate(1,np);
        for (int n=1;n<=np;n++) idxs(n) = (*pBVVI)[n]->getRevIndices();
        os<<"Reverse indices:"<<endl;
        int mnc = idxs(1).indexmin(); int mxc = idxs(1).indexmax();
        for (int n=2;n<=np;n++) {
            mnc = min(mnc,idxs(n).indexmin());
            mxc = max(mxc,idxs(n).indexmax());
        }
        os<<"mnc = "<<mnc<<tb<<"mxc = "<<mxc<<endl;
        imatrix idxps(mnc,mxc,1,np); idxps = -1;
        for (int c=mnc;c<=mxc;c++){
            for (int n=1;n<=np;n++){
                if ((idxs(n).indexmin()<=c)&&(c<=idxs(n).indexmax())) idxps(c,n) = idxs(n,c);
            }
        }
        for (int c=mnc;c<=mxc;c++) os<<c<<tb<<idxps(c)<<endl;
    } else {
        mns =  0;
        mxs =  0;
        phs = -1;
        lb  = -1.0;
        ub  =  1.0;
        idxs.allocate(1,1,1,1);//dummy allocation
        idxs(1,1) = 0;
        os<<"bounded vector vector "<<pBVVI->name<<" has no parameter values"<<endl;
    }
}

/**
 * Set the parameter info for a DevsVectorVectorInfo object.
 * 
 * @param pDVVI - pointer to a DevsVectorVectorInfo instance
 * @param npT - size of vector
 * @param mns - ivector with minimum indices for each vector
 * @param mxs - ivector with maximum indices for each vector
 * @param idxs - imatrix of reverse indices
 * @param lb - dvector of lower bounds
 * @param ub - dvector of upper bounds
 * @param phs - ivector of phases for parameters
 * @param os - output stream to write to
 */
void setParameterInfo(DevsVectorVectorInfo* pDVVI,                           
                             int& npT,
                             ivector& mns, ivector& mxs,
                             imatrix& idxs,
                             dvector& lb, dvector& ub, 
                             ivector& phs, 
                             ostream& os){
    int np = pDVVI->getSize();
    if (np){npT = np;} else {npT = 1;}
    mns.allocate(1,npT);
    mxs.allocate(1,npT);
    lb.allocate(1,npT);
    ub.allocate(1,npT);
    phs.allocate(1,npT);
    if (np){
        phs = pDVVI->getPhases();
        mns = pDVVI->getMinIndices();
        mxs = pDVVI->getMaxIndices();
        lb  = pDVVI->getLowerBounds();
        ub  = pDVVI->getUpperBounds();
        os<<"parameter vector"<<pDVVI->name<<":"<<endl;
        os<<"#mnIdx  mxIdx  lower  upper  phase"<<endl;
        for (int n=1;n<=np;n++) os<<n<<tb<<mns(n)<<tb<<mxs(n)<<tb<<lb(n)<<tb<<ub(n)<<tb<<phs(n)<<endl;
        idxs.allocate(1,np);
        for (int n=1;n<=np;n++) idxs(n) = (*pDVVI)[n]->getRevIndices();
        os<<"Reverse indices:"<<endl;
        int mnc = idxs(1).indexmin(); int mxc = idxs(1).indexmax();
        for (int n=2;n<=np;n++) {
            mnc = min(mnc,idxs(n).indexmin());
            mxc = max(mxc,idxs(n).indexmax());
        }
        os<<"mnc = "<<mnc<<tb<<"mxc = "<<mxc<<endl;
        imatrix idxps(mnc,mxc,1,np); idxps = -1;
        for (int c=mnc;c<=mxc;c++){
            for (int n=1;n<=np;n++){
                if ((idxs(n).indexmin()<=c)&&(c<=idxs(n).indexmax())) idxps(c,n) = idxs(n,c);
            }
        }
        for (int c=mnc;c<=mxc;c++) os<<c<<tb<<idxps(c)<<endl;
    } else {
        mns =  0;
        mxs =  0;
        phs = -1;
        lb  = -1.0;
        ub  =  1.0;
        idxs.allocate(1,1,1,1);//dummy allocation
        idxs(1,1) = 0;
        os<<"devs vector vector "<<pDVVI->name<<" has no parameter values"<<endl;
    }
}



} //namespace tcsam