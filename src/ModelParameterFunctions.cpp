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
 * @param devs - the dvar_matrix
 * @param pDevs - the param_init_bounded_vector_vector
 * @param debug - debugging level
 * @param cout  - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - This alters the values in the dvar_matrix.
 * 
 */
void setDevs(dvar_matrix& devs, param_init_bounded_vector_vector& pDevs, int debug, std::ostream& cout){
    if (debug>=tcsam::dbgAll) cout<<"starting setDevs(devs,pDevs)"<<std::endl;
    int nv = pDevs.indexmax();//number of devs vectors defined
    int mni; int mxi;
    for (int v=1;v<=nv;v++){
        mni = pDevs(v).indexmin();
        mxi = pDevs(v).indexmax();
        devs(v)(mni,mxi) = pDevs(v);
        devs(v,mxi+1) = -sum(devs(v)(mni,mxi));
        if (debug>=(tcsam::dbgAll+1)) cout<<v<<":  "<<devs(v)<<endl;
    }
    if (debug>=(tcsam::dbgAll+1)) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>nv;
        if (nv<0) exit(-1);
        cout<<"finished setDevs(devs,pDevs)"<<std::endl;
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

} //namespace tcsam