//ModelParameterFunctions.cpp

#include <admodel.h>
#include "wtsADMB.hpp"
#include "ModelConstants.hpp"
#include "ModelParameterInfoTypes.hpp"
#include "ModelParameterFunctions.hpp"
#include "ModelParametersInfo.hpp"
#include "ModelOptions.hpp"



namespace tcsam{

/**
 *Set initial values for a param_init_bounded_number_vector from an associated 
 * BoundedNumberVectorInfo instance.
 * 
 * @param pI : pointer to BoundedNumberVectorInfo object
 * @param p  : reference to a param_init_bounded_number_vector
 * @param debug: flag to print debugging info
 * @param os: stream to print debugging info to
 *     
 * @return void
 * 
 * @alters - changes initial values of p
 * 
 */
void setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p, int debug, std::ostream& os){
    debug=dbgAll;
    if (debug>=dbgAll) os<<"Starting setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        dvector vls = pI->getInitValsOnParamScales();
        p.set_initial_value(vls);
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" on arith scale: "<<pI->getInitVals()<<std::endl;
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" on param scale: "<<p<<std::endl;
        if (debug>=dbgAll) {
            os<<"vls = "<<vls<<std::endl;
            os<<"p   = "<<p<<std::endl;
            double_index_guts* ptr = new dvector_index(vls);
            os<<ptr<<tb<<ptr->indexmin()<<tb<<ptr->indexmax()<<endl;
            for (int i=ptr->indexmin();i<=ptr->indexmax();i++) os<<(*((*ptr)[i]))<<tb; std::cout<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedNumberVectorInfo* pI, param_init_bounded_number_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

/**
* Set initial values for a param_init_bounded_vector_vector from its associated
* BoundedVectorVectorInfo instance.
* 
* Inputs:
*  @param pI : pointer to a BoundedVectorVectorInfo object
*  @param p  : reference to a param_init_bounded_vector_vector
 * @param debug: flag to print debugging info
 * @param os: stream to print debugging info to
*     
* @return void
* 
* @alters - changes initial values of p
* 
*/
void setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p, int debug, std::ostream& os){
    debug=dbgAll;
    if (debug>=dbgAll) std::cout<<"Starting setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    int np = pI->getSize();
    if (np){
        for (int i=1;i<=np;i++) {
            dvector vls = (*pI)[i]->getInitValsOnParamScale();
            if (debug>=dbgAll) os<<"pc "<<i<<" :"<<tb<<p(i).indexmin()<<tb<<p(i).indexmax()<<tb<<vls.indexmin()<<tb<<vls.indexmax()<<std::endl;
            for (int j=vls.indexmin();j<=vls.indexmax();j++) p(i,j) = vls(j);
            if (debug>=dbgAll) {
                os<<"vls  = "<<vls<<std::endl;
                os<<"p(i) = "<<p(i)<<std::endl;
            }
        }
        for (int i=1;i<=np;i++) {
            rpt::echo<<"InitVals for "<<p(i).get_name()<<" on arith scale:"<<tb<<(*pI)[i]->getInitVals()<<std::endl;
            rpt::echo<<"InitVals for "<<p(i).get_name()<<" on param scale:"<<tb<<p(i)<<std::endl;
        }
    } else {
        rpt::echo<<"InitVals for "<<p(1).get_name()<<" not defined because np = "<<np<<std::endl;
    }
    
    if (debug>=dbgAll) {
        std::cout<<"Enter 1 to continue >>";
        std::cin>>np;
        if (np<0) exit(-1);
        std::cout<<"Finished setInitVals(BoundedVectorVectorInfo* pI, param_init_bounded_vector_vector& p) for "<<p(1).get_name()<<std::endl; 
    }
}

/**
 * Set values for a dvar_matrix, intended to function as a vector of devs parameter vectors, 
 * based on values of a param_init_bounded_number_vector.
 * 
 * @param devs  - the dvar_matrix
 * @param pDevs - the param_init_bounded_vector_vector
 * @param pI    - pointer to the associated DevsVectorVectorInfo object
 * @param debug - debugging level
 * @param os  - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - This alters the values in the devs dvar_matrix.
 * 
 */
void setDevsVectorVector(dvar_matrix& devs, param_init_bounded_number_vector& pDevs, DevsVectorVectorInfo* pI, int debug, std::ostream& os){
    if (debug>=dbgAll) os<<"starting setDevs(devs,pDevs,pI)"<<std::endl;
    devs.initialize();
    if (pI->getSize()){
        int nv = pI->getSize();//number of devs vectors defined
        ivector mniv = pI->getMinIndices();
        ivector mxiv = pI->getMaxIndices();
        int ctr = 1;
        for (int v=1;v<=nv;v++){
            dvar_vector ps(mniv(v),mxiv(v));
            for (int j=mniv(v);j<=mxiv(v);j++) ps(j) = pDevs(ctr++);
            devs(v) = (*pI)[v]->calcArithScaleVals(ps);
            if (debug>=(dbgAll)) os<<v<<":  "<<devs(v)<<endl;
        }
        if (debug>=(dbgAll+1)) {
            std::cout<<"Enter 1 to continue >>";
            std::cin>>nv;
            if (nv<0) exit(-1);
            os<<"finished setDevs(devs,pDevs,pI)"<<std::endl;
        }
    } else {
        if (debug>=(dbgAll)) os<<"size=0, so no devs to set"<<endl;
    }
}

/**
 * Set values for a dvar_matrix, intended to function as a vector of parameter vectors, 
 * based on values of a param_init_bounded_number_vector.
 * 
 * @param mat   [output]  - the dvar_matrix
 * @param pDevs [input] - the param_init_bounded_vector_vector
 * @param pI    [input] - pointer to the associated DevsVectorVectorInfo object
 * @param debug [input] - debugging level
 * @param os  - [input] output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - This alters the values in the mat dvar_matrix.
 * 
 */
void setBoundedVectorVector(dvar_matrix& mat, param_init_bounded_number_vector& pVals, BoundedVectorVectorInfo* pI, int debug, std::ostream& os){
    if (debug>=dbgAll) os<<"starting setVectorVector(mat,pVals,pI)"<<std::endl;
    mat.initialize();
    if (pI->getSize()){
        int nv = pI->getSize();//number of devs vectors defined
        ivector mniv = pI->getMinIndices();
        ivector mxiv = pI->getMaxIndices();
        int ctr = 1;
        for (int v=1;v<=nv;v++){
            dvar_vector ps(mniv(v),mxiv(v));
            for (int j=mniv(v);j<=mxiv(v);j++) ps(j) = pVals(ctr++);
            mat(v) = (*pI)[v]->calcArithScaleVals(ps);
            if (debug>=(dbgAll)) os<<v<<":  "<<mat(v)<<endl;
        }
        if (debug>=(dbgAll+1)) {
            std::cout<<"Enter 1 to continue >>";
            std::cin>>nv;
            if (nv<0) exit(-1);
            os<<"finished setVectorVector(mat,pVals,pI)"<<std::endl;
        }
    } else {
        if (debug>=(dbgAll)) os<<"size=0, so no values to set"<<endl;
    }
}

/**
 * Calculate ln-scale prior likelihood values for selectivity functions using empirically-determined
 * selectivity function priors.
 * 
 * @param objFun - the objective function
 * @param ptrESPs  - pointer to the EmpiricalSelFcnPriors object from ModelOptions
 * @param sel_cz  - dvar_matrix with parameterized selectivity functions by parameter combination "c"
 * @param debug  - debugging level
 * @param cout   - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - the value of objFun
 */                                      
void calcPriors(objective_function_value& objFun, EmpiricalSelFcnPriors* ptrESPs, dvar_matrix& sel_cz, int debug, std::ostream& cout){
    int nESPs = ptrESPs->nESPs;
    if (nESPs){
        if (debug<0) cout<<"ESPs=list("<<endl;
        dvar_vector params(1,2);
        dvector consts;//leave unallocated
        ivector sel_ids(1,nESPs);
        dvar_vector nlls(1,nESPs); nlls.initialize();
        dvector wgts(1,nESPs);
        adstring_array types(1,nESPs);
        for (int i=1;i<=nESPs;i++){
            EmpiricalSelFcnPrior* ptrESP = ptrESPs->ppESPs[i-1];
            sel_ids[i] = ptrESP->sel_id;
            wgts[i]    = ptrESP->priorWgt;
            types[i]   = ptrESP->priorType;
            for (int j=ptrESP->p1.indexmin();j<=ptrESP->p1.indexmax();j++) {
                params[1] = ptrESP->p1[j];
                params[2] = ptrESP->p2[j];
                dvariable val = sel_cz(sel_ids[i],j);
                nlls[i]  -= ptrESP->pMPI->calcLogPDF(val,params,consts);
            }//--j
            if (debug>=dbgPriors) cout<<"ESP priors("<<sel_ids[i]<<") = "<<nlls[i]<<std::endl;
            objFun += wgts[i]*nlls[i];
        }//--i
        if (debug<0){
            for (int i=1;i<nESPs;i++){
                cout<<tb<<"'"<<sel_ids[i]<<"'=list(sel_id="<<sel_ids[i]<<cc<<"type='"<<types(i)<<"'"<<cc<<"wgt="<<wgts(i)<<cc<<"nll="<<nlls(i)<<cc<<"objfun="<<wgts(i)*nlls(i)<<")"<<endl;
            }
            {
                int i=nESPs;
                cout<<tb<<"'"<<sel_ids[i]<<"'=list(sel_id="<<sel_ids[i]<<cc<<"type='"<<types(i)<<"'"<<cc<<"wgt="<<wgts(i)<<cc<<"nll="<<nlls(i)<<cc<<"objfun="<<wgts(i)*nlls(i)<<")"<<endl;
            }
            cout<<")";
        }
    } else {
        if (debug<0) cout<<"ESPs=NULL";
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
void calcPriors(objective_function_value& objFun, BoundedNumberVectorInfo* ptrVI, param_init_bounded_number_vector& pv, int debug, std::ostream& cout){
    if (ptrVI->getSize()){
        dvar_vector tmp = pv*1.0;
        tmp = ptrVI->calcArithScaleVals(tmp);
        dvar_vector pri = ptrVI->calcLogPriors(tmp);//ln-scale prior (NOT NLL!)
        if (isnan(value(sum(pri)))){
            std::cout<<"Found NAN for priors("<<ptrVI->name<<") = "<<pri<<std::endl;
            std::cout<<"param values = "<<pv<<endl;
            ptrVI->write(std::cout);
            std::cout<<"Aborting!!"<<endl;
            exit(-1);
        }
        if (debug>=dbgPriors) cout<<"priors("<<ptrVI->name<<") = "<<pri<<std::endl;
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

///**
// * Calculate ln-scale priors for a param_init_bounded_vector_vector, 
// * based on its associated BoundedVectorVectorInfo, 
// * and add the weighted NLL to the objective function.
// * 
// * @param objFun - the objective function
// * @param ptrVI  - pointer to the BoundedVectorVectorInfo object associated with pv
// * @param pv     - the param_init_bounded_vector_vector
// * @param debug  - debugging level
// * @param cout   - output stream object for debugging info
// * 
// * @return - void
// * 
// * @alters - the value of objFun
// */                                      
//void calcPriors(objective_function_value& objFun, BoundedVectorVectorInfo* ptrVVI, param_init_bounded_vector_vector& pm, int debug, std::ostream& cout){
//    if (ptrVVI->getSize()){
//        if (debug<0) cout<<ptrVVI->name<<"=list("<<endl;        
//        dvector wts = ptrVVI->getPriorWgts();
//        for (int i=pm.indexmin();i<pm.indexmax();i++) {
//            dvar_vector tmp = 1.0*pm(i);
//            tmp = (*ptrVVI)[i]->calcArithScaleVals(tmp);
//            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
//            if (isnan(value(sum(pri)))){
//                std::cout<<"Found NAN for priors("<<ptrVVI->name<<"["<<i<<"]) = "<<pri<<std::endl;
//                std::cout<<"param values = "<<tmp<<endl;
//                ptrVVI->write(std::cout);
//                std::cout<<"Aborting!!"<<endl;
//                exit(-1);
//            }
//            if (debug>=dbgPriors){
//                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
//                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
//            }
//            objFun += -wts(i)*sum(pri);
//            if (debug<0) {
//                adstring type = (*ptrVVI)[i]->getPriorType();
//                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-sum(pri)<<cc<<"objfun="<<-wts(i)*sum(pri)<<"),"<<endl;
//            }
//        }
//        {
//            int i = pm.indexmax();
//            dvar_vector tmp = 1.0*pm(i);
//            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
//            if (isnan(value(sum(pri)))){
//                std::cout<<"Found NAN for priors("<<ptrVVI->name<<"["<<i<<"]) = "<<pri<<std::endl;
//                std::cout<<"param values = "<<tmp<<endl;
//                ptrVVI->write(std::cout);
//                std::cout<<"Aborting!!"<<endl;
//                exit(-1);
//            }
//            if (debug>=dbgPriors){
//                cout<<"wts["<<ptrVVI->name<<"]("<<i<<") = "<<wts(i)<<std::endl;
//                cout<<"priors["<<ptrVVI->name<<"]("<<i<<") = "<<pri<<std::endl;
//            }
//            objFun += -wts(i)*sum(pri);
//            if (debug<0) {
//                adstring type = (*ptrVVI)[i]->getPriorType();
//                cout<<tb<<"'"<<i<<"'=list(type='"<<type<<"'"<<cc<<"wgt="<<wts(i)<<cc<<"nll="<<-sum(pri)<<cc<<"objfun="<<-wts(i)*sum(pri)<<")"<<endl;
//            }
//        }
//        if (debug<0) cout<<")";
//    } else {
//        if (debug<0) cout<<ptrVVI->name<<"=NULL";
//    }
//}

/**
 * Calculate ln-scale priors for a dvar_matrix acting as a vector_vector, 
 * based on its associated VectorVectorInfo, 
 * and add the weighted NLL to the objective function.
 * 
 * Note that this pertains to VectorVectorInfo and subclasses BoundedVectorVectorInfo and 
 * DevsVectoroVectorInfo.
 * 
 * @param objFun - the objective function
 * @param ptrVI  - pointer to the VectorVectorInfo object associated with pv
 * @param pv     - the dvar_matrix
 * @param debug  - debugging level
 * @param cout   - output stream object for debugging info
 * 
 * @return - void
 * 
 * @alters - the value of objFun
 */                                      
void calcPriors(objective_function_value& objFun, VectorVectorInfo* ptrVVI, dvar_matrix& pm, int debug, std::ostream& cout){
    if (ptrVVI->getSize()){
        if (debug<0) cout<<ptrVVI->name<<"=list("<<endl;        
        dvector wts = ptrVVI->getPriorWgts();
        for (int i=1;i<pm.indexmax();i++) {
            dvar_vector tmp = 1.0*pm(i);
            tmp = (*ptrVVI)[i]->calcArithScaleVals(tmp);
            dvar_vector pri = (*ptrVVI)[i]->calcLogPrior(tmp);//ln-scale prior (NOT NLL!)
            if (isnan(value(sum(pri)))){
                std::cout<<"Found NAN for priors("<<ptrVVI->name<<"["<<i<<"]) = "<<pri<<std::endl;
                std::cout<<"param values = "<<tmp<<endl;
                ptrVVI->write(std::cout);
                std::cout<<"Aborting!!"<<endl;
                exit(-1);
            }
            if (debug>=dbgPriors){
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
           if (debug>=dbgPriors){
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
 * Writes the header (column names) to a stream which subsequent writeParameters functions
 * will write to using csv format.
 * 
 * @param os - the output stream
 */
void writeCSVHeaderForParametersFile(ostream& os){
    os<<"index, phase, min_index, max_index, parameter_scale"<<cc;
    os<<"min_arith, max_arith, value_arith"<<cc;
    os<<"min_param, max_param, value_param"<<cc;
    os<<"name, type, category, process, label"<<std::endl;
}

/**
 * Writes information for a parameter to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter (param_init_number)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated NumberInfo object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write only if parameter will be active
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameter(ostream& os, 
                    param_init_number& p, 
                    adstring& ctg1, 
                    adstring& ctg2, 
                    NumberInfo* pI,
                    int toF, 
                    int willBeActive,
                    int phase){
    if (!willBeActive||((p.get_phase_start()>0)&&((phase==0)||(p.get_phase_start()<=phase)))){
        if (toF>0){
            //write to R file
            os<<p.get_name()<<"=list("<<"type='param_init_number'"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"scale="<<qt<<pI->getScaleType()<<qt<<cc
                                <<"avalue="<<pI->calcArithScaleVal(value(p))<<cc
                                <<"pvalue="<<value(p)
                                <<"),";
        } else if (toF<0){
            os<<"# "<<p.get_name()<<cc<<ctg1<<cc<<ctg2<<cc<<pI->label<<endl;
            os<<value(p)<<endl;
        } else {
            //write to csv file
            os<<1<<cc<<p.get_phase_start()<<cc<<1<<cc<<1<<cc<<pI->getScaleType()<<cc
              <<"-Inf"<<cc<<"Inf"<<cc<<pI->calcArithScaleVal(p)<<cc
              <<"-Inf"<<cc<<"Inf"<<cc<<value(p)                <<cc
              <<p.get_name()<<cc<<"\"param_init_number\",\""
              <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
        }
    }
}    
/**
 * Writes information for a bounded parameter to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter (param_init_bounded_number)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - BoundedNumberInfo
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write only if parameter will be active
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameter(ostream& os, 
                    param_init_bounded_number& p,
                    adstring& ctg1, 
                    adstring& ctg2, 
                    BoundedNumberInfo* pI, 
                    int toF, 
                    int willBeActive,
                    int phase){
    if (!willBeActive||((p.get_phase_start()>0)&&((phase==0)||(p.get_phase_start()<=phase)))){
        if (toF>0){
            //write to R file
            os<<p.get_name()<<"=list("<<"type='param_init_bounded_number'"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"phase="<<p.phase_start<<cc
                                <<"scale="<<qt<<pI->getScaleType()<<qt<<cc
                                <<"bounds=c("<<pI->getLowerBound()<<cc
                                             <<pI->getUpperBound()<<")"<<cc
                                <<"avalue="<<pI->calcArithScaleVal(p)<<cc
                                <<"pvalue="<<value(p)
                                <<"),";
        } else if (toF<0){
            os<<"# "<<p.get_name()<<cc<<ctg1<<cc<<ctg2<<cc<<pI->label<<endl;
            os<<value(p)<<endl;
        } else {
            //write to csv file
            os<<1<<cc<<p.get_phase_start()<<cc<<1<<cc<<1<<cc<<pI->getScaleType()<<cc
              <<pI->getLowerBound()<<cc<<pI->getUpperBound()<<cc<<pI->calcArithScaleVal(p)<<cc
              <<p.get_minb()       <<cc<<p.get_maxb()       <<cc<<value(p)                <<cc
              <<p.get_name()<<cc<<"\"param_init_bounded_number\",\""
              <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;;
        }
    }
}    
/**
 * Writes information for a parameter vector to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter vector (as a dvar_vector)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param lbl - parameter-specific label
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write only if parameter vector will be active
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameter(ostream& os, 
                    dvar_vector& p,
                    adstring& ctg1, 
                    adstring& ctg2, 
                    VectorInfo* pI, 
                    int toF, 
                    int willBeActive,
                    int phase){
    int mn = p.indexmin();
    int mx = p.indexmax();
    int est_flag = 0;
    ivector phases = pI->getPhases();
    for (int i=mn;i<=mx;i++) est_flag += 1.0*(phases[i]>0);
     if (!willBeActive||est_flag){
        dvector vals = pI->calcArithScaleVals(value(p));
        if (toF>0){
            //writing to R file
            os<<pI->name<<"=list("<<"type=param_init_vector"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<max(phases)<<cc
                                <<"scale="<<qt<<pI->getScaleType()<<qt<<cc<<endl;
            os<<"phases=c("; for (int i=mn;i<mx;i++) {os<<phases(i)<<cc;} os<<phases(mx)<<")"<<cc<<endl;
            os<<"avalue=c("; for (int i=mn;i<mx;i++) {os<<vals(i)<<cc;} os<<vals(mx)<<")"<<cc<<endl;
            os<<"pvalue=c("; for (int i=mn;i<mx;i++) {os<<p(i)   <<cc;} os<<p(mx)   <<")";
            os<<"),";
        } else if (toF<0){
            os<<"# "<<pI->name<<cc<<ctg1<<cc<<ctg2<<cc<<pI->label<<endl;
            os<<"# "; for (int i=mn;i<mx;i++) {os<<  i <<" ";} os<<  mx <<endl;
            os<<"  "; for (int i=mn;i<mx;i++) {os<<p(i)<<" ";} os<<p(mx)<<endl;
        } else {
            //writing to csv file
            for (int i=mn;i<=mx;i++) {
                if ((!willBeActive)||((phases[i]>0)&&(phases[i]<=phase))){
                    os<<i<<cc<<phases(i)<<cc<<mn<<cc<<mx<<cc<<pI->getScaleType()<<cc
                                           <<"-Inf"<<cc<<"Inf"<<cc<<vals(i)<<cc
                                           <<"-Inf"<<cc<<"Inf"<<cc<<p(i)<<cc
                                           <<pI->name<<cc<<"\"param_init_vector\",\""
                                           <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
                } else if((phase==0)&&willBeActive){
                    //writing only parameters that will be active (i.e., phases[i]>0)
                    if (phases(i)>0){
                        os<<i<<cc<<phases(i)<<cc<<mn<<cc<<mx<<cc<<pI->getScaleType()<<cc
                                           <<"-Inf"<<cc<<"Inf"<<cc<<vals(i)<<cc
                                           <<"-Inf"<<cc<<"Inf"<<cc<<p(i)<<cc
                                           <<pI->name<<cc<<"\"param_init_bounded_dev_vector\",\""
                                           <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
                    }
                }
            }
        }
    }
}       

/**
 * Writes information for a bounded parameter vector to an output stream.
 * 
 * @param os - output stream
 * @param p - parameter vector (as a dvar_vector)
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedVectorInfo info object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write only if parameter vector will be active
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameter(ostream& os, 
                    dvar_vector& p,
                    adstring& ctg1, 
                    adstring& ctg2, 
                    BoundedVectorInfo* pI, 
                    int toF, 
                    int willBeActive,
                    int phase){
    int mn = p.indexmin();
    int mx = p.indexmax();
    int est_flag = 0;
    ivector phases = pI->getPhases();
    for (int i=mn;i<=mx;i++) est_flag += 1.0*(phases[i]>0);
    if (!willBeActive||est_flag){
        dvector vals = pI->calcArithScaleVals(value(p));
        if (toF>0){
            //writing to R file
            os<<pI->name<<"=list("<<"type=param_init_bounded_vector"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<max(phases)<<cc
                                <<"scale="<<qt<<pI->getScaleType()<<qt<<cc
                                <<"bounds=c("<<pI->getLowerBound()<<cc
                                             <<pI->getUpperBound()<<")"<<cc<<endl;
            os<<"phases=c("; for (int i=mn;i<mx;i++) {os<<phases(i)<<cc;} os<<phases(mx)<<")"<<cc<<endl;
            os<<"avalue=c("; for (int i=mn;i<mx;i++) {os<<vals(i)<<cc;} os<<vals(mx)<<")"<<cc<<endl;
            os<<"pvalue=c("; for (int i=mn;i<mx;i++) {os<<p(i)   <<cc;} os<<p(mx)   <<")";
               os<<"),";
        } else if (toF<0){
            os<<"# "<<pI->name<<cc<<ctg1<<cc<<ctg2<<cc<<pI->label<<endl;
            os<<"# "; for (int i=mn;i<mx;i++) {os<<  i <<" ";} os<<  mx <<endl;
            os<<"  "; for (int i=mn;i<mx;i++) {os<<p(i)<<" ";} os<<p(mx)<<endl;
        } else {
            //writing to csv file
            for (int i=mn;i<=mx;i++) {
                if ((!willBeActive)||((phases[i]>0)&&(phases[i]<=phase))){
                    os<<i<<cc<<phases(i)<<cc<<mn<<cc<<mx<<cc<<pI->getScaleType()<<cc
                                           <<pI->getLowerBound()<<cc<<pI->getUpperBound()<<cc<<vals(i)<<cc
                                           <<pI->calcParamScaleVals(pI->getLowerBound())<<cc
                                           <<pI->calcParamScaleVals(pI->getUpperBound())<<cc
                                           <<p(i)<<cc
                                           <<pI->name<<cc<<"\"param_init_bounded_vector\",\""
                                           <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
                } else if((phase==0)&&willBeActive){
                    //writing only parameters that will be active (i.e., phases[i]>0)
                    if (phases(i)>0){
                        os<<i<<cc<<phases(i)<<cc<<mn<<cc<<mx<<cc<<pI->getScaleType()<<cc
                                           <<pI->getLowerBound()<<cc<<pI->getUpperBound()<<cc<<vals(i)<<cc
                                           <<pI->calcParamScaleVals(pI->getLowerBound())<<cc
                                           <<pI->calcParamScaleVals(pI->getUpperBound())<<cc
                                           <<p(i)<<cc
                                           <<pI->name<<cc<<"\"param_init_bounded_dev_vector\",\""
                                           <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
                    }
                }
            }
        }
    }
}

/**
 * Writes information for a devs parameter vector to an output stream.
 * 
 * @param os - output stream
 * @param p - dvar_vector of devs parameter values
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated DevsVectorInfo info object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write only if parameter vector will be active
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameter(ostream& os, 
                    dvar_vector& p,
                    adstring& ctg1, 
                    adstring& ctg2, 
                    DevsVectorInfo* pI, 
                    int toF, 
                    int willBeActive,
                    int phase){
    int mn = p.indexmin();
    int mx = p.indexmax();
    int est_flag = 0;
    ivector phases = pI->getPhases();
    for (int i=mn;i<=mx;i++) est_flag += 1.0*(phases[i]>0);
    if (!willBeActive||est_flag){
        dvector vals = pI->calcArithScaleVals(value(p));
        if (toF>0){
            //writing to R file
            os<<pI->name<<"=list("<<"type=param_init_bounded_dev_vector"<<cc
                                <<"ctg1="<<qt<<ctg1<<qt<<cc
                                <<"ctg2="<<qt<<ctg2<<qt<<cc
                                <<"lbl="<<qt<<pI->label<<qt<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<max(phases)<<cc
                                <<"scale="<<qt<<pI->getScaleType()<<qt<<cc
                                <<"bounds=c("<<pI->getLowerBound()<<cc
                                             <<pI->getUpperBound()<<")"<<cc<<endl;
            os<<"phases=c("; for (int i=mn;i<mx;i++) {os<<phases(i)<<cc;} os<<phases(mx)<<")"<<cc<<endl;
            os<<"avalue=c("; for (int i=mn;i<mx;i++) {os<<vals(i)<<cc;} os<<vals(mx)<<")"<<cc<<endl;
            os<<"pvalue=c("; for (int i=mn;i<mx;i++) {os<<p(i)   <<cc;} os<<p(mx)   <<")";
               os<<"),";
        } else if (toF<0){
            os<<"# "<<pI->name<<cc<<ctg1<<cc<<ctg2<<cc<<pI->label<<endl;
            os<<"# "; for (int i=mn;i<mx;i++) {os<<  i <<" ";} os<<  mx <<endl;
            os<<"  "; for (int i=mn;i<mx;i++) {os<<p(i)<<" ";} os<<p(mx)<<endl;
        } else {
            //writing to csv file
            for (int i=mn;i<=mx;i++) {
                if ((!willBeActive)||((phases[i]>0)&&(phases[i]<=phase))){
                    //writing all parameters OR currently active parameters
                    os<<i<<cc<<phases(i)<<cc<<mn<<cc<<mx<<cc<<pI->getScaleType()<<cc
                                       <<pI->getLowerBound()<<cc<<pI->getUpperBound()<<cc<<vals(i)<<cc
                                       <<pI->calcParamScaleVals(pI->getLowerBound())<<cc
                                       <<pI->calcParamScaleVals(pI->getUpperBound())<<cc
                                       <<p(i)<<cc
                                       <<pI->name<<cc<<"\"param_init_bounded_dev_vector\",\""
                                       <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
                } else if((phase==0)&&willBeActive){
                    //writing only parameters that will be active (i.e., phases[i]>0)
                    if (phases(i)>0){
                        os<<i<<cc<<phases(i)<<cc<<mn<<cc<<mx<<cc<<pI->getScaleType()<<cc
                                           <<pI->getLowerBound()<<cc<<pI->getUpperBound()<<cc<<vals(i)<<cc
                                           <<pI->calcParamScaleVals(pI->getLowerBound())<<cc
                                           <<pI->calcParamScaleVals(pI->getUpperBound())<<cc
                                           <<p(i)<<cc
                                           <<pI->name<<cc<<"\"param_init_bounded_dev_vector\",\""
                                           <<ctg1<<"\",\""<<ctg2<<"\",\""<<pI->label<<"\""<<endl;
                    }
                }
            }
        }
    }
}    

/**
 * Writes a vector of parameters (param_init_number_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_number_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated NumberVectorInfo info object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write if parameters will be active in some phase
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameters(ostream& os, 
                     param_init_number_vector& p, 
                     adstring& ctg1, 
                     adstring& ctg2, 
                     NumberVectorInfo* pI,
                     int toF, 
                     int willBeActive,
                     int phase){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toF,willBeActive,phase);
        }
    } else if (toF<0) {
        os<<"# "<<pI->name<<endl;
        os<<0<<endl;
    }
}

/**
 * Writes a vector of parameters (param_init_bounded_number_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a param_init_bounded_number_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedNumberVectorInfo info object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write if parameters will be active in some phase
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameters(ostream& os, 
                     param_init_bounded_number_vector& p, 
                     adstring& ctg1, 
                     adstring& ctg2, 
                     BoundedNumberVectorInfo* pI,
                     int toF, 
                     int willBeActive,
                     int phase){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toF,willBeActive,phase);
        }
    } else if (toF<0) {
        os<<"# "<<pI->name<<endl;
        os<<0<<endl;
    }
}

/**
 * Writes a vector of parameter vectors (param_init_vector_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a dvar_matrix representing a param_init_vector_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated VectorVectorInfo info object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write if parameters will be active in some phase
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameters(ostream& os, 
                     dvar_matrix& p, 
                     adstring& ctg1, 
                     adstring& ctg2, 
                     VectorVectorInfo* pI,
                     int toF, 
                     int willBeActive,
                     int phase){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toF,willBeActive,phase);
        }
    } else if (toF<0) {
        os<<"# "<<pI->name<<endl;
        os<<0<<endl;
    }
}

/**
 * Writes a vector of bounded parameter vectors (param_init_bounded_vector_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a dvar_matrix representing a param_init_bounded_vector_vector
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated BoundedVectorVectorInfo info object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write if parameters will be active in some phase
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameters(ostream& os, 
                     dvar_matrix& p, 
                     adstring& ctg1, 
                     adstring& ctg2, 
                     BoundedVectorVectorInfo* pI,
                     int toF, 
                     int willBeActive,
                     int phase){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toF,willBeActive,phase);
        }
    } else if (toF<0) {
        os<<"# "<<pI->name<<endl;
        os<<0<<endl;
    }
}

/**
 * Writes a vector of devs parameter vectors (representing a param_init_devs_vector_vector) to R or csv.
 * 
 * @param os - output stream to write to
 * @param p - a dvar_matrix representing a "param_init_devs_vector_vector"
 * @param ctg1 - category 1 label
 * @param ctg2 - category 2 label
 * @param pI - pointer to associated DevsVectorVectorInfo info object
 * @param toF - format flag for output: csv (=0), R list (=1), or par file (=-1) 
 * @param willBeActive - flag to write if parameters will be active in some phase
 * @param phase - current phase (or 0 to ignore)
 */
void writeParameters(ostream& os, 
                     dvar_matrix& p, 
                     adstring& ctg1, 
                     adstring& ctg2, 
                     DevsVectorVectorInfo* pI,
                     int toF, 
                     int willBeActive,
                     int phase){
    if (pI->getSize()){
        for (int i=p.indexmin();i<=p.indexmax();i++) {
            tcsam::writeParameter(os,p[i],ctg1,ctg2,(*pI)[i],toF,willBeActive,phase);
        }
    } else if (toF<0) {
        os<<"# "<<pI->name<<endl;
        os<<0<<endl;
    }
}

/**
 * Sets the info for a param_init_number_vector from a NumberVectorInfo object.
 * 
 * @param pNVI - pointer to a NumberVectorInfo instance
 * @param npT [out] - size of vector
 * @param phs [out] - ivector of phases for parameters
 * @param os  - output stream to write processing info to
 */
void setParameterInfo(NumberVectorInfo* pNVI,                           
                     int& npT,
                     ivector& phs, 
                     ostream& os){
    int np = pNVI->getSize();
    if (np){npT = np;} else {npT = 1;}
    phs.allocate(1,npT);
    adstring_array scales(1,npT);
    if (np){
        phs = pNVI->getPhases();
        scales = pNVI->getScaleTypes();
        os<<"parameter "<<pNVI->name<<":"<<endl;
        os<<"#phase   scale"<<endl;
        for (int n=1;n<=np;n++) os<<n<<tb<<phs(n)<<tb<<scales(n)<<endl;
    } else {
        phs = -1;
        os<<"number vector parameter "<<pNVI->name<<" has no parameter values"<<endl;
    }
}
  
/**
 * Sets the info for a param_init_bounded_number_vector from a BoundedNumberVectorInfo object.
 * 
 * @param pBNVI - pointer to a BoundedNumberVectorInfo instance
 * @param npT [out]  - size of vector
 * @param lb  [out]  - dvector of lower bounds on parameter scale
 * @param ub  [out]  - dvector of upper bounds on parameter scale
 * @param phs [out] - ivector of phases for parameters
 * @param os - output stream to write to
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
        lb  = pBNVI->getLowerBoundsOnParamScales();
        ub  = pBNVI->getUpperBoundsOnParamScales();
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
 * Sets the info for a param_init_number_vector substituting as a
 * param_init_vector_vector from a VectorVectorInfo object.
 * 
 * @param pVVI - pointer to a VectorVectorInfo instance
 * @param npV [out] - number of vectors represented
 * @param npT [out] - total size (number of elements)
 * @param mns [out] - ivector with minimum indices for each vector
 * @param mxs [out] - ivector with maximum indices for each vector
 * @param phs [out] - ivector of phases for parameters
 * @param os - output stream to write to
 */
void setParameterInfo(VectorVectorInfo* pVVI,
                      int& npV,
                      int& npT,
                      ivector& mns, 
                      ivector& mxs,
                      imatrix& idxs,
                      ivector& phs, 
                      ostream& os){
    os<<"tcsam::setParameterInfo(VectorVectorInfo*,...):"<<endl;
    int np = pVVI->getSize();
    if (np){npV = np; npT = pVVI->getNumParameters();} else {npV = 1; npT=1;}
    mns.deallocate(); mns.allocate(1,npV);
    mxs.deallocate(); mxs.allocate(1,npV);
    phs.deallocate(); phs.allocate(1,npT);
    adstring_array scales(1,npV);
    if (np){
        os<<"npV = "<<npV<<endl;
        os<<pVVI->getMinIndices()<<endl;
        os<<pVVI->getMaxIndices()<<endl;
        mns = pVVI->getMinIndices();
        mxs = pVVI->getMaxIndices();
        os<<"npT = "<<npT<<endl;
        os<<"parameter vector "<<pVVI->name<<":"<<endl;
        os<<"#mnIdx  mxIdx  phase  scale"<<endl;
        phs = pVVI->getParameterPhases();
        scales = pVVI->getScaleTypes();
        int ctr = 1;
        for (int ip=1;ip<=np;ip++) {
            for (int j=mns(ip);j<=mxs(ip);j++) {
                os<<ip<<tb<<mns(ip)<<tb<<mxs(ip)<<tb<<phs(ctr)<<tb<<scales(ip)<<endl;
                ctr++;
            }
        }
        
        //create reverse indices
        idxs.allocate(1,np);//allocating first index of imatrix idxs
        for (int n=1;n<=np;n++) idxs(n) = (*pVVI)[n]->getRevIndices();
        os<<"Reverse indices:"<<endl;
        int mnc = idxs(1).indexmin(); int mxc = idxs(1).indexmax();
        for (int n=2;n<=np;n++) {
            mnc = min(mnc,idxs(n).indexmin());
            mxc = max(mxc,idxs(n).indexmax());
        }
        //check for correctness
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
        os<<"vector vector "<<pVVI->name<<" has no parameter values"<<endl;
    }
}

/**
 * Sets the info for a param_init_bounded_number_vector substituting as a 
 * param_init_bounded_vector_vector from a BoundedVectorVectorInfo object.
 * 
 * @param pBVVI - pointer to a BoundedVectorVectorInfo instance
 * @param npV [out] - number of vectors represented
 * @param npT [out] - total size (number of elements)
 * @param mxs [out] - ivector with maximum indices for each vector
 * @param idxs [out] - imatrix of reverse indices
 * @param lb [out] - dvector of lower bounds
 * @param ub [out] - dvector of upper bounds
 * @param phs [out] - ivector of phases for parameters
 * @param os - output stream to write to
 */
void setParameterInfo(BoundedVectorVectorInfo* pBVVI,                           
                      int& npV,
                      int& npT,
                      ivector& mns, 
                      ivector& mxs,
                      imatrix& idxs,
                      dvector& lb, 
                      dvector& ub, 
                      ivector& phs, 
                      ostream& os){
    os<<"tcsam::setParameterInfo(BoundedVectorVectorInfo*,...):"<<endl;
    int np = pBVVI->getSize();
    if (np){npV = np; npT = pBVVI->getNumParameters();} else {npV = 1; npT=1;}
    mns.deallocate(); mns.allocate(1,npV);
    mxs.deallocate(); mxs.allocate(1,npV);
    lb.deallocate();  lb.allocate(1,npT);
    ub.deallocate();  ub.allocate(1,npT);
    phs.deallocate(); phs.allocate(1,npT);
    adstring_array scales(1,npV);
    if (np){
        os<<"npV = "<<npV<<endl;
        os<<pBVVI->getMinIndices()<<endl;
        os<<pBVVI->getMaxIndices()<<endl;
        mns = pBVVI->getMinIndices();//vector of minimum indices
        mxs = pBVVI->getMaxIndices();//vector of maximum indices
        os<<"npT = "<<npT<<endl;
        os<<"parameter vector "<<pBVVI->name<<":"<<endl;
        os<<"#mnIdx  mxIdx  lower  upper  phase  scale"<<endl;
        int ctr = 1;
        for (int ip=1;ip<=np;ip++) {
            double lbv  = (*pBVVI)[ip]->getLowerBoundOnParamScale();
            double ubv  = (*pBVVI)[ip]->getUpperBoundOnParamScale();
            for (int j=mns(ip);j<=mxs(ip);j++) {
                lb(ctr)    = lbv;
                ub(ctr)    = ubv;
                ctr++;
            }
        }//--ip loop
        phs = pBVVI->getParameterPhases();
        scales = pBVVI->getScaleTypes();
        ctr = 1;
        for (int ip=1;ip<=np;ip++) {
            for (int j=mns(ip);j<=mxs(ip);j++) {
                os<<ip<<tb<<mns(ip)<<tb<<mxs(ip)<<tb<<lb(ctr)<<tb<<ub(ctr)<<tb<<phs(ctr)<<tb<<scales(ip)<<endl;
                ctr++;
            }
        }
        
        //create reverse indices
        idxs.allocate(1,np);//allocating first index of imatrix idxs
        for (int n=1;n<=np;n++) idxs(n) = (*pBVVI)[n]->getRevIndices();
        os<<"Reverse indices:"<<endl;
        int mnc = idxs(1).indexmin(); int mxc = idxs(1).indexmax();
        for (int n=2;n<=np;n++) {
            mnc = min(mnc,idxs(n).indexmin());
            mxc = max(mxc,idxs(n).indexmax());
        }
        //check for correctness
        imatrix idxps(mnc,mxc,1,np); idxps = -1;
        for (int c=mnc;c<=mxc;c++){
            for (int n=1;n<=np;n++){
                if ((idxs(n).indexmin()<=c)&&(c<=idxs(n).indexmax())) idxps(c,n) = idxs(n,c);
            }
        }
        for (int c=mnc;c<=mxc;c++) os<<c<<tb<<idxps(c)<<endl;
    } else {
        mns =  1;
        mxs =  1;
        phs = -1;
        lb  = -1.0;
        ub  =  1.0;
        idxs.allocate(1,1,1,1);//dummy allocation
        idxs(1,1) = 0;
        os<<"BoundedVectorVector "<<pBVVI->name<<" has no parameter values"<<endl;
    }
}

/**
 * Sets the info for a param_init_bounded_number_vector acting as a devs_vector_vector
 * from a DevsVectorVectorInfo object.
 * 
 * @param pDVVI - pointer to a DevsVectorVectorInfo instance
 * @param npV [out] - number of devs vectors represented
 * @param npT [out] - total size (number of elements)
 * @param mns [out] - ivector with minimum indices for each devs vector
 * @param mxs [out] - ivector with maximum indices for each devs vector
 * @param idxs [out] - imatrix of reverse indices
 * @param lb [out] - dvector of lower bounds
 * @param ub [out] - dvector of upper bounds
 * @param phs [out] - ivector of phases for parameters
 * @param os - output stream to write to
 */
void setParameterInfo(DevsVectorVectorInfo* pDVVI,     
                      int& npV,
                      int& npT,
                      ivector& mns, 
                      ivector& mxs,
                      imatrix& idxs,
                      dvector& lb, 
                      dvector& ub, 
                      ivector& phs, 
                      ostream& os){
    os<<"tcsam::setParameterInfo(DevsVectorVectorInfo*,...):"<<endl;
    int np = pDVVI->getSize();//number of devs vectors
    if (np){npV = np; npT = pDVVI->getNumParameters();} else {npV = 1; npT=1;}
    mns.deallocate(); mns.allocate(1,npV);
    mxs.deallocate(); mxs.allocate(1,npV);
    lb.deallocate();  lb.allocate(1,npT);
    ub.deallocate();  ub.allocate(1,npT);
    phs.deallocate(); phs.allocate(1,npT);
    adstring_array scales(1,npV);
    if (np){
        os<<"npV = "<<npV<<endl;
        os<<pDVVI->getMinIndices()<<endl;
        os<<pDVVI->getMaxIndices()<<endl;
        mns = pDVVI->getMinIndices();//vector of minimum indices
        mxs = pDVVI->getMaxIndices();//vector of maximum indices
        os<<"npT = "<<npT<<endl;
        os<<"parameter vector "<<pDVVI->name<<":"<<endl;
        os<<"#mnIdx  mxIdx  lower  upper  phase  scale"<<endl;
        int ctr = 1;
        for (int ip=1;ip<=np;ip++) {
            double lbv  = (*pDVVI)[ip]->getLowerBoundOnParamScale();
            double ubv  = (*pDVVI)[ip]->getUpperBoundOnParamScale();
            for (int j=mns(ip);j<=mxs(ip);j++) {
                lb(ctr)    = lbv;
                ub(ctr)    = ubv;
                ctr++;
            }
        }//--ip loop
        phs = pDVVI->getParameterPhases();
        scales = pDVVI->getScaleTypes();
        ctr = 1;
        for (int ip=1;ip<=np;ip++) {
            for (int j=mns(ip);j<=mxs(ip);j++) {
                os<<ip<<tb<<mns(ip)<<tb<<mxs(ip)<<tb<<lb(ctr)<<tb<<ub(ctr)<<tb<<phs(ctr)<<tb<<scales(ip)<<endl;
                ctr++;
            }
        }
        
        //create reverse indices
        idxs.allocate(1,np);//allocating first index of imatrix idxs
        for (int n=1;n<=np;n++) idxs(n) = (*pDVVI)[n]->getRevIndices();
        os<<"Reverse indices:"<<endl;
        int mnc = idxs(1).indexmin(); int mxc = idxs(1).indexmax();
        for (int n=2;n<=np;n++) {
            mnc = min(mnc,idxs(n).indexmin());
            mxc = max(mxc,idxs(n).indexmax());
        }
        //check for correctness
        imatrix idxps(mnc,mxc,1,np); idxps = -1;
        for (int c=mnc;c<=mxc;c++){
            for (int n=1;n<=np;n++){
                if ((idxs(n).indexmin()<=c)&&(c<=idxs(n).indexmax())) idxps(c,n) = idxs(n,c);
            }
        }
        for (int c=mnc;c<=mxc;c++) os<<c<<tb<<idxps(c)<<endl;
    } else {
        mns =  1;
        mxs =  1;
        phs = -1;
        lb  = -1.0;
        ub  =  1.0;
        idxs.allocate(1,1,1,1);//dummy allocation
        idxs(1,1) = 0;
        os<<"DevsVectorVector "<<pDVVI->name<<" has no parameter values"<<endl;
    }
}

} //namespace tcsam