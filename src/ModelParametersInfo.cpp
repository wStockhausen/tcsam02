#include <admodel.h>
#include <wtsADMB.hpp>
#include "ModelConstants.hpp"
#include "ModelFunctions.hpp"
#include "ModelIndexFunctions.hpp"
#include "ModelParameterInfoTypes.hpp"
#include "ModelParametersInfo.hpp"
#include "ModelSelectivities.hpp"

/*----------------------------------------------------------------------------\n
**  Includes:
**      ParametersGroupInfo\n
**      RecruitmentInfo
**      NaturalMortalityInfo\n
**      GrowthInfo\n
**      MaturityInfo\n
**      SelectivityInfo\n
**      FisheriesInfo\n
**      SurveysInfo\n
**      ModelParametersInfo\n
**----------------------------------------------------------------------------*/
int ParameterGroupInfo::debug   = 0;
int RecruitmentInfo::debug      = 0;
int NaturalMortalityInfo::debug = 0;
int GrowthInfo::debug           = 0;
int Molt2MaturityInfo::debug    = 0;
int SelectivityInfo::debug      = 0;
int FisheriesInfo::debug        = 0;
int SurveysInfo::debug          = 0;
int MSE_Info::debug             = 0;
int ModelParametersInfo::debug  = 0;
const adstring ModelParametersInfo::version = "2018.10.29";
    
/*----------------------------------------------------------------------------*/
/**
 * Function to convert parameter combinations to R dimensions 
 * @param pgi
 * @return wts::adstring_matrix
 */
wts::adstring_matrix tcsam::convertPCs(ParameterGroupInfo * pgi){
    //cout<<"starting tcsam::convertPCs for "<<pgi->name<<endl;
    int nIVs = pgi->nIVs;
    int nPCs = pgi->nPCs;
    //cout<<"nIVs = "<<nIVs<<"; nPCs = "<<nPCs<<endl;
    //wts::adstring_matrix::debug = 1;
    wts::adstring_matrix a(1,nIVs,0,nPCs);
    //wts::adstring_matrix::debug = 0;
    //cout<<"allocated a"<<endl;
    int ibsIdx = 1;
    for (int i=1;i<=nIVs;i++){
        if (pgi->lblIVs(i)==tcsam::STR_SEX)             {a(i,0)=pgi->lblIVs(i);} else
        if (pgi->lblIVs(i)==tcsam::STR_MATURITY_STATE)  {a(i,0)=pgi->lblIVs(i);} else
        if (pgi->lblIVs(i)==tcsam::STR_SHELL_CONDITION) {a(i,0)=pgi->lblIVs(i);} else 
        if ((pgi->ibsIdxs.allocated())&&(i==pgi->ibsIdxs(ibsIdx))){
            a(i,0)=pgi->ppIBSs[ibsIdx-1]->getType();
            if (ibsIdx<pgi->nIBSs) ibsIdx++;//increment to next
        }
    }
    for (int r=1;r<=nPCs;r++){//loop over rows
        //cout<<"r = "<<r<<endl;
        int ibsIdx = 1;
        for (int i=1;i<=nIVs;i++){//loop over index variables
            //cout<<tb<<"i = "<<i<<" "<<pgi->lblIVs(i)<<tb<<pgi->in(r,i)<<endl;
            if (pgi->lblIVs(i)==tcsam::STR_SEX)             {a(i,r)=tcsam::getSexType(pgi->in(r,i));}      else
            if (pgi->lblIVs(i)==tcsam::STR_MATURITY_STATE)  {a(i,r)=tcsam::getMaturityType(pgi->in(r,i));} else
            if (pgi->lblIVs(i)==tcsam::STR_SHELL_CONDITION) {a(i,r)=tcsam::getShellType(pgi->in(r,i));}    else 
            if ((pgi->ibsIdxs.allocated())&&(i==pgi->ibsIdxs(ibsIdx))){
                a(i,r)=pgi->ppIBSs[ibsIdx-1]->getIndexBlock(r)->asString();
                if (ibsIdx<pgi->nIBSs) ibsIdx++;//increment to next
            } else {a(i,r)=str(pgi->in(r,i));}
        }
    }
    //cout<<a<<endl;
    //cout<<"finished tcsam::convertPCs for "<<pgi->name<<endl;
    return a;
}
/*----------------------------------------------------------------------------\n
 * ParameterGroupInfo\n
 -----------------------------------------------------------------------------*/
/**
 * Constructor.
 */
ParameterGroupInfo::ParameterGroupInfo(){
    nIVs=0; nPVs=0; nXIs=0; nIBSs=0; nPCs=0;
    ppIBSs=0; ppIdxs=0;
}
/**
 * Destructor.
 */
ParameterGroupInfo::~ParameterGroupInfo(){
    if (ppIBSs) {
        for (int i=0;i<nIBSs;i++) delete ppIBSs[i]; 
        delete ppIBSs;
        nIBSs=0;
    }
    if (ppIdxs) {
        for (int p=0;p<nPCs;p++) delete ppIdxs[p]; 
        delete ppIdxs;
        ppIdxs=0;
    }
}
/**
 * Gets indices for parameter combination pc.
 * @param pc - id for desired parameter combination
 * @return 
 */
ivector ParameterGroupInfo::getPCIDs(int pc){
    if (debug) cout<<"starting ParameterGroupInfo::getPCIDs(int pc)"<<endl;
    if (debug) {
        cout<<"pcids for "<<pc<<"th parameter combination: "<<in(pc)<<endl;
        cout<<"finished ParameterGroupInfo::getPCIDs(int pc)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    return in(pc);
}

/**
 * Gets "extra" values for parameter combination pc.
 * @param pc - id for desired parameter combination
 * @return 
 */
dvector ParameterGroupInfo::getPCXDs(int pc){
    if (debug) cout<<"starting ParameterGroupInfo::getPCXDs(int pc)"<<endl;
    if (debug) {
        cout<<"pcids for "<<pc<<"th parameter combination: "<<xd(pc)<<endl;
        cout<<"finished ParameterGroupInfo::getPCXDs(int pc)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    return xd(pc);
}

/**
 * Gets a matrix of model indices corresponding 
 * to the given parameter combination pc.
 * 
 * @param pc
 * @return - imatrix of model indices
 */
imatrix ParameterGroupInfo::getModelIndices(int pc){
    if (debug) cout<<"starting ParameterGroupInfo::getModelIndices(int pc)"<<endl;
    if (pc>nPCs) {
        cout<<"Error in ParameterGroupInfo::getModelIndices(int pc) for "<<name<<endl;
        cout<<"pc was "<<pc<<" but max defined is "<<nPCs<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    imatrix idxs = (*ppIdxs[pc-1]);
    if (debug) {
        cout<<"imatrix idxs has rows "<<idxs.indexmin()<<":"<<idxs.indexmax()<<endl;
        cout<<"finished ParameterGroupInfo::getModelIndices(int pc)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
    return idxs;
}

/**
 * Creates the pointer to the IndexBlockSets object 
 * and the individual IndexBlockSet objects corresponding
 * to the blocks defined for model indices for this 
 * ParameterGroup object.
 */
void ParameterGroupInfo::createIndexBlockSets(){
    if (debug) cout<<"starting ParameterGroupInfo::createIndexBlockSets() "<<nIBSs<<endl;
    if (nIBSs){
        ppIBSs = new IndexBlockSet*[nIBSs];
        for (int i=1;i<=nIBSs;i++){
            ppIBSs[i-1] = new IndexBlockSet();
        }
    }
    if (debug) cout<<"finished ParameterGroupInfo::createIndexBlockSets()"<<endl;
}

/* 
 * Returns a pointer to the index block set identified by "type".
 * Inputs:
 *  adstring type:  "type" identifying index block set to return
 * Returns:
 *  pointer to the identified IndexBlockSet
 */
IndexBlockSet* ParameterGroupInfo::getIndexBlockSet(adstring type){
    if (debug) cout<<"starting  IndexBlockSets::getIndexBlockSet("<<type<<")"<<endl;
    IndexBlockSet* p = 0;
    int s=0;
    while(s<nIBSs){
        p = ppIBSs[s];
        if (p->getType()==type) break;
        s++;
    }
    if (debug) cout<<"finished  IndexBlockSets::getIndexBlockSet("<<type<<")"<<endl;
    return p;
}

/**
 * Creates the set of imatrix's representing the model indices 
 * corresponding to each parameter combination.
 */
void ParameterGroupInfo::createIndices(void){
    if (debug) cout<<"starting void ParameterGroupInfo::createIndices(void)"<<endl;
    //create pointer array
    ppIdxs = new imatrix*[nPCs];
    //loop over parameter combinations and create an indices imatrix for each
    for (int p=1;p<=nPCs;p++){
        int nc=1;//number of rows for indices imatrix
        imatrix tmp(1,nIVs);//matrix of values of indices
        for (int i=1;i<=nIVs;i++){//loop over index variables
            if (lblIVs(i)==tcsam::STR_SEX){
                tmp(i).allocate(1,1);
                tmp(i,1) = in(p,i);
                nc *= 1;//just for clarity
            } else
            if (lblIVs(i)==tcsam::STR_MATURITY_STATE){
                tmp(i).allocate(1,1);
                tmp(i,1) = in(p,i);
                nc *= 1;//just for clarity
            } else
            if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {
                tmp(i).allocate(1,1);
                tmp(i,1) = in(p,i);
                nc *= 1;//just for clarity
            } else
            if (lblIVs(i)==tcsam::STR_FISHERY) {
                tmp(i).allocate(1,1);
                tmp(i,1) = in(p,i);
                nc *= 1;//just for clarity
            } else
            if (lblIVs(i)==tcsam::STR_SURVEY) {
                tmp(i).allocate(1,1);
                tmp(i,1) = in(p,i);
                nc *= 1;//just for clarity
            } else {
                if (debug) cout<<"Parsing label '"<<lblIVs(i)<<"' for IndexBlockSet dim type."<<endl;
                //parse label (should be of form 'type_BLOCK')
                int n = lblIVs(i).pos("_");
                if (n) { //can parse label correctly
                    adstring type = lblIVs(i)(1,n-1);
                    int id = in(p,i);//block index in IndexBlockSet
                    IndexBlockSet* pIBS = 0;
                    if (debug) cout<<"Checking local IndexBlockSets using type '"<<type<<"'."<<endl;
//                    pIBS = ptrIBSs->getIndexBlockSet(type);
                    pIBS = ppIBSs[in(p,i)-1];
                    if (pIBS){
                        ivector idxs = pIBS->getFwdIndexVector(p);
                        tmp(i).allocate(idxs.indexmin(),idxs.indexmax());
                        tmp(i) = idxs;
                        nc *= idxs.size();
                        if (debug) cout<<"indices for block = "<<tmp(i)<<endl;
                    } else {
                        cout<<"Error in ParameterGroupInfo::createIndices(void) for "<<name<<endl;
                        cout<<"Could not find IndexBlockSet for dim type '"<<type<<endl;
                        cout<<"Dim label was '"<<lblIVs(i)<<"'"<<endl;
                        cout<<"Aborting..."<<endl;
                        exit(-1);
                    }
                } else {
                    cout<<"Error in ParameterGroupInfo::createIndices(void) for "<<name<<endl;
                    cout<<"Could not parse index variable label '"<<lblIVs(i)<<"' correctly"<<endl;
                    cout<<"Should have form 'type_BLOCK', where 'type' is dimension type (e.g., 'YEAR')."<<endl;
                    cout<<"Aborting..."<<endl;
                    exit(-1);
                }
            }            
        }
        ppIdxs[p-1] = new imatrix();
        (*ppIdxs[p-1]).allocate(1,nc,1,nIVs);
        imatrix m(1,nIVs,1,nc);
        for (int i=1;i<=nIVs;i++){//loop over index variables
            ivector itmp = tmp(i);
            int sz = itmp.size();
            int cmn = 1;
            while ((cmn+sz-1)<=nc){
                m(i)(cmn,cmn+sz-1) = itmp.shift(cmn);
                cmn += sz;
            }
        }
        (*ppIdxs[p-1]) = trans(m);
        if (debug) {
            cout<<"Index matrix for pc "<<p<<endl;
            cout<<"#";
            for (int i=1;i<=nIVs;i++) cout<<lblIVs(i)<<tb; cout<<endl;
            cout<<(*ppIdxs[p-1])<<endl;  
        }
    }
    if (debug) {
        cout<<"finished void ParameterGroupInfo::createIndices(void)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Reads info for a BoundedNumberVector object from an input filestream.
 * @param is
 * @param lbl
 * @param pBNVI
 * @return 
 */
BoundedNumberVectorInfo* ParameterGroupInfo::read(cifstream& is, adstring& lbl, BoundedNumberVectorInfo* pBNVI){
    adstring param;
    is>>param;
    if (param==lbl){
        if (debug) cout<<"Reading "<<lbl<<endl;
        pBNVI = new BoundedNumberVectorInfo(lbl);
        is>>(*pBNVI);
        if (debug) cout<<"ptr to "<<lbl<<": "<<pBNVI<<endl<<(*pBNVI)<<endl;
    } else {
        tcsam::readError(is,lbl,param);
        exit(-1);
    }
    return pBNVI;
}

/**
 * 
 * @param is
 * @param lbl
 * @param pBVVI
 * @return 
 */
BoundedVectorVectorInfo* ParameterGroupInfo::read(cifstream& is, adstring& lbl, BoundedVectorVectorInfo* pBVVI){
    if (debug) cout<<"Starting BVVI* PGI::read(is,lbl,BVVI*) for "<<lbl<<endl;
    adstring param;
    is>>param;
    if (param==lbl){
        if (debug) cout<<"Reading "<<lbl<<endl;
        pBVVI = new BoundedVectorVectorInfo(lbl);
        is>>(*pBVVI);
        if (debug) cout<<"ptr to "<<lbl<<": "<<pBVVI<<endl<<(*pBVVI)<<endl;
    } else {
        tcsam::readError(is,lbl,param);
        exit(-1);
    }
    if (debug) cout<<"Finished BVVI* PGI::read(is,lbl,BVVI*) for "<<lbl<<endl;
    return pBVVI;
}

/**
 * 
 * @param is
 * @param lbl
 * @param pDVVI
 * @return 
 */
DevsVectorVectorInfo* ParameterGroupInfo::read(cifstream& is, adstring& lbl, DevsVectorVectorInfo* pDVVI){
    adstring param;
    is>>param;
    if (param==lbl){
        if (debug) cout<<"Reading "<<lbl<<endl;
        pDVVI = new DevsVectorVectorInfo(lbl);
        is>>(*pDVVI);
        if (debug) cout<<"ptr to "<<lbl<<": "<<pDVVI<<endl<<(*pDVVI)<<endl;
    } else {
        tcsam::readError(is,lbl,param);
        exit(-1);
    }
    return pDVVI;
}

/**
 * Reads info for the parameter group from an input filestream.
 * @param is
 */
void ParameterGroupInfo::read(cifstream& is){
    if (debug) cout<<"starting void ParameterGroupInfo::read(is)"<<endl;
    adstring str; 
    //read parameter combinations
    is>>str;
    if (debug) cout<<str<<tb<<"#Required keyword (PARAMETER_COMBINATIONS)"<<endl;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETER_COMBINATIONS)"<<endl;
    if (str=="PARAMETER_COMBINATIONS"){
        is>>nPCs;
        if (debug) cout<<nPCs<<tb<<"#number of parameter combinations"<<endl;
        if (debug) cout<<"#id  "; 
        rpt::echo<<nPCs<<tb<<"#number of parameter combinations"<<endl;
        rpt::echo<<"#id  "; 
        if (debug){
            for (int i=1;i<=nIVs;i++) cout<<lblIVs(i)<<tb; 
            for (int i=1;i<=nPVs;i++) cout<<lblPVs(i)<<tb; 
            for (int i=1;i<=nXIs;i++) cout<<lblXIs(i)<<tb; 
        }
        for (int i=1;i<=nIVs;i++) rpt::echo<<lblIVs(i)<<tb; 
        for (int i=1;i<=nPVs;i++) rpt::echo<<lblPVs(i)<<tb; 
        for (int i=1;i<=nXIs;i++) rpt::echo<<lblXIs(i)<<tb; 
        if (debug) cout<<tb<<"label"<<tb<<endl;
        rpt::echo<<tb<<"label"<<tb<<endl;
        if (debug) cout<<"About to allocate ppIBSs"<<endl;
        for (int i=0;i<nIBSs;i++) ppIBSs[i]->allocate(nPCs);
        if (debug) cout<<"Finished allocating ppIBSs"<<endl;
        //read parameters combinations definition matrix
        int ibsIdx=1;
        in.allocate(1,nPCs,1,nIVs+nPVs+nXIs);
        if (nXIs) xd.allocate(1,nPCs,1,nXIs);
        pcLabels.allocate(1,nPCs);
        for (int r=1;r<=nPCs;r++){//loop over rows
            is>>str; //read id
            rpt::echo<<str<<tb;
            if (debug) cout<<"pc = "<<r<<tb<<"looping over index variables"<<endl;
            for (int i=1;i<=nIVs;i++){//loop over index variables
                is>>str;
                if (debug) cout<<tb<<"index variable"<<tb<<i<<tb<<str<<endl;
                rpt::echo<<str<<tb;
                if (lblIVs(i)==tcsam::STR_SEX)             {in(r,i) = tcsam::getSexType(str);}      else
                if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {in(r,i) = tcsam::getMaturityType(str);} else
                if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {in(r,i) = tcsam::getShellType(str);}    else
                if ((ibsIdxs.allocated())&&(i==ibsIdxs(ibsIdx))){//variable is a block
                    ppIBSs[ibsIdx-1]->getIndexBlock(r)->parse(str);
                    in(r,i)=ibsIdx;//index to associated IndexBlockSet
                    if (ibsIdx<nIBSs) ibsIdx++;//increment to next IBS
                } else {in(r,i)=::atoi(str);}
            }
            if (debug) cout<<"looping over parameter indices"<<endl;
            for (int p=1;p<=nPVs;p++) {is>>in(r,nIVs+p); rpt::echo<<in(r,nIVs+p)<<tb;}  
            if (debug) cout<<"looping over extra variables"<<endl;
            for (int x=1;x<=nXIs;x++) {//loop over "extra" variables
                is>>str;
                if (debug) cout<<str<<endl;
                rpt::echo<<str<<tb;
                if (lblXIs(x)==tcsam::STR_SELFCN){
                    //identify selectivity function and return function index
                    in(r,nIVs+nPVs+x) = SelFcns::getSelFcnID(str);
                } else {
                    in(r,nIVs+nPVs+x)=::atoi(str);
                    xd(r,x)=::atof(str);
                }
                //rpt::echo<<"x = "<<x<<tb<<"str = "<<str<<tb<<"in() = "<<in(r,nIVs+nPVs+x)<<tb<<"xd() = "<<xd(r,x)<<endl;
            }
            is>>pcLabels(r);
            if (debug) cout<<pcLabels(r)<<endl;
            rpt::echo<<pcLabels(r)<<tb;
            rpt::echo<<endl;
            if (debug) {
                cout<<"pc row "<<r<<": "<<in(r);
                if (nXIs) cout<<tb<<"xd: "<<xd(r);
                cout<<endl;
            }
        }
        //revert reading back to sub-class to read values for parameters
    } else {
        cout<<"Reading ParameterGroupInfo for '"<<name<<"' from file "<<is.get_file_name()<<endl;
        cout<<"Expected 'PARAMETER_COMBINATIONS' keyword but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug) ParameterGroupInfo::write(cout);        
    
    if (nIVs) createIndices();
    
    //now back to sub-class to read values for parameters
    if (debug) cout<<"finished void ParameterGroupInfo::read(is)"<<endl;
}

//void ParameterGroupInfo::setToWriteVectorInitialValues(bool flag){
//    //do nothing [not sure why have to provide a concrete
//}

void ParameterGroupInfo::write(std::ostream& os){
    if (debug) cout<<"starting ParameterGroupInfo::write(std::ostream& os)"<<endl;
//    os<<(*ptrIBSs)<<endl;
    os<<"PARAMETER_COMBINATIONS"<<endl;
    os<<nPCs<<tb<<"#number of parameter combinations"<<endl;
    os<<"#id  "; 
    for (int i=1;i<=nIVs;i++) os<<lblIVs(i)<<tb; 
    for (int i=1;i<=nPVs;i++) os<<lblPVs(i)<<tb; 
    for (int i=1;i<=nXIs;i++) os<<lblXIs(i)<<tb; 
    os<<"label"<<tb;
    os<<endl;
    for (int r=1;r<=nPCs;r++){//loop over rows
        os<<r<<tb;
        int ibsIdx = 1;
        for (int i=1;i<=nIVs;i++){//loop over index variables           
            if (lblIVs(i)==tcsam::STR_SEX)             {os<<tcsam::getSexType(in(r,i))<<tb;} else
            if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {os<<tcsam::getMaturityType(in(r,i))<<tb;} else
            if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {os<<tcsam::getShellType(in(r,i))<<tb;} else 
            if ((ibsIdxs.allocated())&&(i==ibsIdxs(ibsIdx))){
                os<<(*ppIBSs[ibsIdx-1]->getIndexBlock(r))<<tb;
                if (ibsIdx<nIBSs) ibsIdx++;//increment to next
            } else {os<<in(r,i)<<tb;}
        }
        for (int p=1;p<=nPVs;p++) os<<in(r,nIVs+p)<<tb;      //loop over parameter variables
        for (int x=1;x<=nXIs;x++) { //loop over extra indices
            if (lblXIs(x)==tcsam::STR_SELFCN) os<<SelFcns::getSelFcnID(in(r,nIVs+nPVs+x))<<tb; else
            os<<xd(r,x)<<tb;
    //        os<<in(r,nIVs+nPVs+x)<<tb;
            
        }
        os<<pcLabels(r)<<tb;
        os<<endl;
    }//r
    if (debug) cout<<"finished ParameterGroupInfo::write(std::ostream& os)"<<endl;
}
/**
 * This should be called from sub-classes to write generic information to R.
 * @param os
 */
void ParameterGroupInfo::writeToR(std::ostream& os){
    adstring lbls = "";
//    if (nIVs) lbls += wts::to_qcsv(lblIVs);
//    if (nPVs) {if (lbls=="") lbls += wts::to_qcsv(lblPVs); else lbls += cc+wts::to_qcsv(lblPVs);}
    if (nPVs) lbls += wts::to_qcsv(lblPVs);
    if (nXIs) {if (lbls=="") lbls += wts::to_qcsv(lblXIs); else lbls += cc+wts::to_qcsv(lblXIs);}
    os<<"pgi=list(name="<<qt<<name<<qt<<cc<<endl;
        os<<"pcs=list("<<endl;
            for (int p=1;p<=nPCs;p++){
                os<<qt<<p<<qt<<"=list(";
                os<<"label='"<<pcLabels(p)<<"',"<<endl;
                int ibsIdx = 1;
                for (int i=1;i<=nIVs;i++){//loop over index variables
                    os<<lblIVs(i)<<"='";
                    if (lblIVs(i)==tcsam::STR_SEX)             {os<<tcsam::getSexType(in(p,i))     <<"',";} else
                    if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {os<<tcsam::getMaturityType(in(p,i))<<"',";} else
                    if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {os<<tcsam::getShellType(in(p,i))   <<"',";} else 
                    if ((ibsIdxs.allocated())&&(i==ibsIdxs(ibsIdx))){
                        os<<(*ppIBSs[ibsIdx-1]->getIndexBlock(p))<<"',";
                        if (ibsIdx<nIBSs) ibsIdx++;//increment to next
                    } else {os<<in(p,i)<<"',";}
                }
                ivector iv = getPCIDs(p);
                imatrix im = getModelIndices(p);
                os<<"ids.PC="; wts::writeToR(os,iv(nIVs+1,nIVs+nPVs+nXIs),lbls); os<<cc<<endl;
                adstring ids = "index=c(1:"+str(im.indexmax())+")";
                adstring tps = "type=c("+wts::to_qcsv(lblIVs)+")";
                os<<"ids.Mod="; wts::writeToR(os,im,ids,tps); os<<"),"<<endl;
            }
        os<<"NULL)"<<endl;
    os<<")";
}
/*------------------------------------------------------------------------------
 * RecruitmentInfo
 -----------------------------------------------------------------------------*/
adstring RecruitmentInfo::NAME = "recruitment";
RecruitmentInfo::RecruitmentInfo(){
    if (debug) cout<<"starting RecruitmentInfo::RecruitmentInfo()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 1;
    lblIVs.allocate(1,nIVs);
    lblIVs(k=1) = "YEAR_BLOCK";
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs = 6;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pLnR";    dscPVs(k++) = "ln-scale mean recruitment";
    lblPVs(k) = "pRCV";    dscPVs(k++) = "recruitment cv's";
    lblPVs(k) = "pRX";     dscPVs(k++) = "fraction males at recruitment";
    lblPVs(k) = "pRa";     dscPVs(k++) = "size-at-recruitment parameter a";
    lblPVs(k) = "pRb";     dscPVs(k++) = "size-at-recruitment parameter b";    
    lblPVs(k) = "pDevsLnR";dscPVs(k++) = "ln-scale recruitment devs";    
    pLnR = 0;
    pRCV = 0;
    pRX  = 0;
    pRa  = 0;
    pRb  = 0;
    pDevsLnR = 0;
    
    nXIs=1;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "nllWgt";
}

RecruitmentInfo::~RecruitmentInfo(){
    if (pLnR) delete pLnR;  pLnR = 0;
    if (pRCV) delete pRCV;  pRCV = 0;
    if (pRX)  delete pRX;   pRX  = 0;
    if (pRa)  delete pRa;   pRa  = 0;
    if (pRb)  delete pRb;   pRb  = 0;
    if (pDevsLnR) delete pDevsLnR; pDevsLnR = 0;
}

void RecruitmentInfo::read(cifstream & is){
    if (debug) cout<<"starting void RecruitmentInfo::read(cifstream & is)"<<endl;    
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in RecruitmentInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        pLnR   = ParameterGroupInfo::read(is,lblPVs(k),pLnR);   
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnR)<<endl;   k++;
        pRCV = ParameterGroupInfo::read(is,lblPVs(k),pRCV); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRCV)<<endl; k++;
        pRX = ParameterGroupInfo::read(is,lblPVs(k),pRX); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRX)<<endl; k++;
        pRa  = ParameterGroupInfo::read(is,lblPVs(k),pRa);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRa)<<endl;  k++;
        pRb  = ParameterGroupInfo::read(is,lblPVs(k),pRb);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pRb)<<endl;  k++;
        pDevsLnR = ParameterGroupInfo::read(is,lblPVs(k),pDevsLnR); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsLnR)<<endl; 
    } else {
        cout<<"Error reading RecruitmentInfo from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    if (debug){
        cout<<"RecruitmentInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished RecruitmentInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void RecruitmentInfo::setToWriteVectorInitialValues(bool flag){
    if (pDevsLnR){
        for (int i=1;i<=pDevsLnR->getSize();i++){
            DevsVectorInfo* dvi = (*pDevsLnR)[i];
            if (flag) dvi->readVals = INT_TRUE; else dvi->readVals = INT_FALSE;        
        }
    }
}

/**
 * Write RecruitmentInfo to output stream in ADMB format
 * 
 * @param os - output stream
 */
void RecruitmentInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnR)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pRCV)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pRX)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pRa)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pRb)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsLnR)<<endl;
}

/**
 * Write RecruitmentInfo to output stream in R format
 * 
 * @param os - output stream
 */
void RecruitmentInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"rec=list("<<endl;
        ParameterGroupInfo::writeToR(os);   os<<cc<<endl;
        pLnR->writeToR(os,"pLnR",indent+1); os<<cc<<endl;
        pRCV->writeToR(os,"pRCV",indent+1); os<<cc<<endl;
        pRX->writeToR(os, "pRX", indent+1); os<<cc<<endl;
        pRa->writeToR(os, "pRa", indent+1); os<<cc<<endl;
        pRb->writeToR(os, "pRb", indent+1); os<<cc<<endl;
        pDevsLnR->writeToR(os,"pDevsLnR",indent+1); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * NaturalMortalityInfo
 -----------------------------------------------------------------------------*/
adstring NaturalMortalityInfo::NAME = "natural_mortality";
NaturalMortalityInfo::NaturalMortalityInfo(){
    if (debug) cout<<"starting NaturalMortalityInfo::NaturalMortalityInfo()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    if (debug) cout<<1<<endl;
    
    int k;
    //define index variables for parameters
    nIVs = 4;
    lblIVs.allocate(1,nIVs);
    k = 1;
    lblIVs(k++) = "YEAR_BLOCK";
    lblIVs(k++) = tcsam::STR_SEX;
    lblIVs(k++) = tcsam::STR_MATURITY_STATE;
    lblIVs(k++) = tcsam::STR_SHELL_CONDITION;
    if (debug) cout<<2<<endl;
    //define index block sets for index variables
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    if (debug) cout<<3<<endl;
    
    //define parameters
    nPVs=5;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pM";   dscPVs(k++) = "base natural mortality rate";
    lblPVs(k) = "pDM1"; dscPVs(k++) = "offset 1";
    lblPVs(k) = "pDM2"; dscPVs(k++) = "offset 2";
    lblPVs(k) = "pDM3"; dscPVs(k++) = "offset 3";
    lblPVs(k) = "pDM4"; dscPVs(k++) = "offset 4";
    if (debug) cout<<3<<endl;
    pM  = 0;
    pDM1  = 0;
    pDM2  = 0;
    pDM3  = 0;
    pDM4 = 0;
    if (debug) cout<<4<<endl;
    
    //define "extra" indices
    nXIs=1;
    lblXIs.allocate(1,nXIs);
    lblXIs(k=1) = "zScaling";    
    if (debug) cout<<5<<endl;
    if (debug) cout<<"finished NaturalMortalityInfo::NaturalMortalityInfo()"<<endl;
}

NaturalMortalityInfo::~NaturalMortalityInfo(){
    if (pM)   delete pM;    pM =0;
    if (pDM1) delete pDM1;  pDM1 =0;
    if (pDM2) delete pDM2;  pDM2 =0;
    if (pDM3) delete pDM3;  pDM3 =0;
    if (pDM4) delete pDM4;  pDM4 =0;
}

void NaturalMortalityInfo::read(cifstream & is){
    if (debug) cout<<"starting NaturalMortalityInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in NaturalMortalityInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    is>>zRef;    
    rpt::echo<<"zRef = "<<zRef<<endl;
    if (debug) ParameterGroupInfo::debug=1;
    ParameterGroupInfo::read(is);
    if (debug) ParameterGroupInfo::debug=0;
    
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        pM    = ParameterGroupInfo::read(is,lblPVs(k),pM);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pM)<<endl;    k++;
        pDM1  = ParameterGroupInfo::read(is,lblPVs(k),pDM1);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM1)<<endl;  k++;
        pDM2  = ParameterGroupInfo::read(is,lblPVs(k),pDM2);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM2)<<endl;  k++;
        pDM3  = ParameterGroupInfo::read(is,lblPVs(k),pDM3);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM3)<<endl;  k++;
        pDM4 = ParameterGroupInfo::read(is,lblPVs(k),pDM4); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDM4)<<endl; k++;
    } else {
        cout<<"Error reading NaturalMortalityInfo from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    if (debug){
        cout<<"NMInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished NaturalMortalityInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void NaturalMortalityInfo::setToWriteVectorInitialValues(bool flag){
    //does nothing: no parameter vectors
}

void NaturalMortalityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    os<<zRef<<tb<<"#reference size for scaling"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pM)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDM1)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDM2)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDM3)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDM4)<<endl;
 }

void NaturalMortalityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"nm=list("<<endl;
        ParameterGroupInfo::writeToR(os);    os<<cc<<endl;
        pM  ->writeToR(os,"pM", indent+1);   os<<cc<<endl;
        pDM1->writeToR(os,"pDM1", indent+1); os<<cc<<endl;
        pDM2->writeToR(os,"pDM2", indent+1); os<<cc<<endl;
        pDM3->writeToR(os,"pDM3", indent+1); os<<cc<<endl;
        pDM4->writeToR(os,"pDM4", indent+1); os<<endl;
    os<<")";
}
        
/*------------------------------------------------------------------------------
 * GrowthInfo
 -----------------------------------------------------------------------------*/
adstring GrowthInfo::NAME = "growth";
GrowthInfo::GrowthInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 2;
    lblIVs.allocate(1,nIVs);
    lblIVs(1) = "YEAR_BLOCK";
    lblIVs(2) = tcsam::STR_SEX;
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs=3;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pGrA";    dscPVs(k++) = "mean growth coefficient 'a'";
    lblPVs(k) = "pGrB";    dscPVs(k++) = "mean growth coefficient 'b'";
    lblPVs(k) = "pGrBeta"; dscPVs(k++) = "growth transition matrix scale factor";
    pGrA    = 0;
    pGrB    = 0;
    pGrBeta = 0;
    
    //define "extra" indices
    nXIs=2;
    lblXIs.allocate(1,nXIs);
    lblXIs(k=1) = "zScaleGrA"; //pre-molt size corresponding to pGrA in alt parameterization
    lblXIs(k=2) = "zScaleGrB"; //pre-molt size corresponding to pGrA in alt parameterization
    
}

GrowthInfo::~GrowthInfo(){
    if (pGrA)    delete pGrA;    pGrA    =0;
    if (pGrB)    delete pGrB;    pGrB    =0;
    if (pGrBeta) delete pGrBeta; pGrBeta =0;
}

void GrowthInfo::read(cifstream & is){
    if (debug) cout<<"starting GrowthInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in GrowthInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        pGrA     = ParameterGroupInfo::read(is,lblPVs(k),pGrA);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pGrA)<<endl;    k++;
        pGrB     = ParameterGroupInfo::read(is,lblPVs(k),pGrB);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pGrB)<<endl;    k++;
        pGrBeta  = ParameterGroupInfo::read(is,lblPVs(k),pGrBeta); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pGrBeta)<<endl; k++;
     } else {
        cout<<"Error reading GrowthInfo from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    if (debug){
        cout<<"GrowthInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished GrowthInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void GrowthInfo::setToWriteVectorInitialValues(bool flag){
    //does nothing: no parameter vectors
}

void GrowthInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
     ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pGrA)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pGrB)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pGrBeta)<<endl;
 }

void GrowthInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"grw=list("<<endl;
        ParameterGroupInfo::writeToR(os);         os<<cc<<endl;
        pGrA->writeToR(os,   "pGrA",   indent+1); os<<cc<<endl;
        pGrB->writeToR(os,   "pGrB",   indent+1); os<<cc<<endl;
        pGrBeta->writeToR(os,"pGrBeta",indent+1); os<<endl;
    os<<")";
}
        
/*------------------------------------------------------------------------------
 * MaturityInfo
 -----------------------------------------------------------------------------*/
adstring Molt2MaturityInfo::NAME = "molt_to_maturity";
Molt2MaturityInfo::Molt2MaturityInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 2;
    lblIVs.allocate(1,nIVs);
    lblIVs(1) = "YEAR_BLOCK";
    lblIVs(2) = tcsam::STR_SEX;
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs = 1;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pvLgtPrM2M";dscPVs(k++) = "logit-scale parameter vectors for Pr(molt-to-maturity|size)";    
    pvLgtPrM2M = 0;
    
    nXIs = 0;    
}

Molt2MaturityInfo::~Molt2MaturityInfo(){
    if (pvLgtPrM2M) delete pvLgtPrM2M; pvLgtPrM2M = 0;
}

void Molt2MaturityInfo::read(cifstream & is){
    if (debug) cout<<"starting void MaturityInfo::read(cifstream & is)"<<endl;    
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in MaturityInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        pvLgtPrM2M = ParameterGroupInfo::read(is,lblPVs(k),pvLgtPrM2M);
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pvLgtPrM2M)<<endl;  
    } else {
        cout<<"Error reading MaturityInfo from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    if (debug){
        cout<<"MaturityInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished MaturityInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void Molt2MaturityInfo::setToWriteVectorInitialValues(bool flag){
    if (pvLgtPrM2M){
        for (int i=1;i<=pvLgtPrM2M->getSize();i++){
            BoundedVectorInfo* vi = (*pvLgtPrM2M)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
}

void Molt2MaturityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pvLgtPrM2M)<<endl;
}

void Molt2MaturityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"mat=list("<<endl;
        ParameterGroupInfo::writeToR(os);             os<<cc<<endl;
        pvLgtPrM2M->writeToR(os,"pvLgtPrM2M",indent+1); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * SelectivityInfo
 -----------------------------------------------------------------------------*/
adstring SelectivityInfo::NAME = "selectivities";
SelectivityInfo::SelectivityInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 1;
    lblIVs.allocate(1,nIVs);
    lblIVs(1) = "YEAR_BLOCK";
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 1;
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs=13;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pS1"; dscPVs(k++) = "1st input to selectivity function";
    lblPVs(k) = "pS2"; dscPVs(k++) = "2nd input to selectivity function";
    lblPVs(k) = "pS3"; dscPVs(k++) = "3rd input to selectivity function";
    lblPVs(k) = "pS4"; dscPVs(k++) = "4th input to selectivity function";
    lblPVs(k) = "pS5"; dscPVs(k++) = "5th input to selectivity function";
    lblPVs(k) = "pS6"; dscPVs(k++) = "6th input to selectivity function";
    lblPVs(k) = "pDevsS1"; dscPVs(k++) = "devs to 1st input to selectivity function";
    lblPVs(k) = "pDevsS2"; dscPVs(k++) = "devs to 2nd input to selectivity function";
    lblPVs(k) = "pDevsS3"; dscPVs(k++) = "devs to 3rd input to selectivity function";
    lblPVs(k) = "pDevsS4"; dscPVs(k++) = "devs to 4th input to selectivity function";
    lblPVs(k) = "pDevsS5"; dscPVs(k++) = "devs to 5th input to selectivity function";
    lblPVs(k) = "pDevsS6"; dscPVs(k++) = "devs to 6th input to selectivity function";
    lblPVs(k) = "pvNPSel"; dscPVs(k++) = "non-parametric selectivity functions";
    
    nXIs=2;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "fsZ";
    lblXIs(k++) = "selFcn";
    
}

SelectivityInfo::~SelectivityInfo(){
    if (pS1) delete pS1; pS1=0;
    if (pS2) delete pS2; pS2=0;
    if (pS3) delete pS3; pS3=0;
    if (pS4) delete pS4; pS4=0;
    if (pS5) delete pS5; pS5=0;
    if (pS6) delete pS6; pS6=0;
    if (pDevsS1) delete pDevsS1; pDevsS1=0;
    if (pDevsS2) delete pDevsS2; pDevsS2=0;
    if (pDevsS3) delete pDevsS3; pDevsS3=0;
    if (pDevsS4) delete pDevsS4; pDevsS4=0;
    if (pDevsS5) delete pDevsS5; pDevsS5=0;
    if (pDevsS6) delete pDevsS6; pDevsS6=0;
    if (pvNPSel) delete pvNPSel; pvNPSel=0;
}

void SelectivityInfo::read(cifstream & is){
    if (debug) cout<<"starting SelectivityInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in SelectivityInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    
    int debugPGI = ParameterGroupInfo::debug;
    if (debug) ParameterGroupInfo::debug=1;
    ParameterGroupInfo::read(is);
    if (debug) ParameterGroupInfo::debug=debugPGI;
    
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        pS1 = ParameterGroupInfo::read(is,lblPVs(k),pS1); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS1)<<endl;  k++;
        pS2 = ParameterGroupInfo::read(is,lblPVs(k),pS2); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS2)<<endl;  k++;
        pS3 = ParameterGroupInfo::read(is,lblPVs(k),pS3); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS3)<<endl;  k++;
        pS4 = ParameterGroupInfo::read(is,lblPVs(k),pS4); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS4)<<endl;  k++;
        pS5 = ParameterGroupInfo::read(is,lblPVs(k),pS5); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS5)<<endl;  k++;
        pS6 = ParameterGroupInfo::read(is,lblPVs(k),pS6); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pS6)<<endl;  k++;
        pDevsS1 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS1); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS1)<<endl;  k++;
        pDevsS2 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS2); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS2)<<endl;  k++;
        pDevsS3 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS3); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS3)<<endl;  k++;
        pDevsS4 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS4); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS4)<<endl;  k++;
        pDevsS5 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS5); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS5)<<endl;  k++;
        pDevsS6 = ParameterGroupInfo::read(is,lblPVs(k),pDevsS6); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsS6)<<endl;  k++;
        pvNPSel = ParameterGroupInfo::read(is,lblPVs(k),pvNPSel); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pvNPSel)<<endl;  k++;
    } else {
        cout<<"Error reading SelectivityInfo from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
   
    if (debug){
        cout<<"SelectivityInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished SelectivityInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void SelectivityInfo::setToWriteVectorInitialValues(bool flag){
    if (pDevsS1){
        for (int i=1;i<=pDevsS1->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS1)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
    if (pDevsS2){
        for (int i=1;i<=pDevsS2->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS2)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
    if (pDevsS3){
        for (int i=1;i<=pDevsS3->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS3)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
    if (pDevsS4){
        for (int i=1;i<=pDevsS4->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS4)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
    if (pDevsS5){
        for (int i=1;i<=pDevsS5->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS5)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
    if (pDevsS6){
        for (int i=1;i<=pDevsS6->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsS6)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
    if (pvNPSel){
        for (int i=1;i<=pvNPSel->getSize();i++){
            BoundedVectorInfo* vi = (*pvNPSel)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
}

void SelectivityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pS1)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pS2)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pS3)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pS4)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pS5)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pS6)<<endl;
    
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsS1)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsS2)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsS3)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsS4)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsS5)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsS6)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pvNPSel)<<endl;
 }

void SelectivityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"sel=list("<<endl;
        ParameterGroupInfo::writeToR(os); os<<cc<<endl;
        pS1->writeToR(os,"pS1",indent+1); os<<cc<<endl;
        pS2->writeToR(os,"pS2",indent+1); os<<cc<<endl;
        pS3->writeToR(os,"pS3",indent+1); os<<cc<<endl;
        pS4->writeToR(os,"pS4",indent+1); os<<cc<<endl;
        pS5->writeToR(os,"pS5",indent+1); os<<cc<<endl;
        pS6->writeToR(os,"pS6",indent+1); os<<cc<<endl;
        pDevsS1->writeToR(os,"pDevsS1",indent+1); os<<cc<<endl;
        pDevsS2->writeToR(os,"pDevsS2",indent+1); os<<cc<<endl;
        pDevsS3->writeToR(os,"pDevsS3",indent+1); os<<cc<<endl;
        pDevsS4->writeToR(os,"pDevsS4",indent+1); os<<cc<<endl;
        pDevsS5->writeToR(os,"pDevsS5",indent+1); os<<cc<<endl;
        pDevsS6->writeToR(os,"pDevsS6",indent+1); os<<cc<<endl;
        pvNPSel->writeToR(os,"pvNPSel",indent+1); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * FisheriesInfo
 -----------------------------------------------------------------------------*/
adstring FisheriesInfo::NAME = "fisheries";
int FisheriesInfo::idxHM     = 5+1;//column in parameter combinations matrix with parameter index for column in parameter combinations matrix indicating handling mortality parameters
int FisheriesInfo::idxLnC    = 5+2;//column in parameter combinations matrix with parameter index for ln-scale base mean capture rate (mature males)
int FisheriesInfo::idxDC1    = 5+3;//column in parameter combinations matrix with parameter index for main year_block ln-scale offsets
int FisheriesInfo::idxDC2    = 5+4;//column in parameter combinations matrix with parameter index for ln-scale female offsets
int FisheriesInfo::idxDC3    = 5+5;//column in parameter combinations matrix with parameter index for ln-scale immature offsets
int FisheriesInfo::idxDC4    = 5+6;//column in parameter combinations matrix with parameter index for ln-scale female-immature offsets 
int FisheriesInfo::idxLnDevs = 5+7;//column in parameter combinations matrix with parameter index for annual ln-scale devs w/in year_blocks
int FisheriesInfo::idxLnEffX = 5+8;//column in parameter combinations matrix with parameter index for ln-scale effort extrapolation 
int FisheriesInfo::idxLgtRet = 5+9;//column in parameter combinations matrix with parameter index for logit-scale retained fraction (for old shell crab)
int FisheriesInfo::idxSelFcn = 5+9+1;//column in parameter combinations matrix indicating selectivity function index
int FisheriesInfo::idxRetFcn = 5+9+2;//column in parameter combinations matrix indicating retention function index
int FisheriesInfo::idxUseEX  = 5+9+3;//column in parameter combinations matrix indicating effort extrapolation use
/*------------------------------------------------------------------------------
 * FisheriesInfo\n
 * Encapsulates the following fishery-related parameters:\n
 *   pHM    : handling mortality (0-1)
 *   pLnC   : ln-scale base mean capture rate (mature males)
 *   pDC1   : ln-scale offset 1 (e.g., main year_block ln-scale offsets)
 *   pDC2   : ln-scale offset 2 (e.g., female offsets)
 *   pDC3   : ln-scale offset 3 (e.g., immature offsets)
 *   pDC4   : ln-scale offset 4 (e.g., female-immature offsets)   
 *
 *   pDevsLnC : annual ln-scale devs w/in year_blocks
 * 
 *   pLnEffX: ln-scale effort extrapolation parameters
 *   pLgtRet: logit-scale max retention parameter
 *
 * Notes:
 *  1. FISHERY, YEAR_BLOCK, SEX, MATURITY_STATE, SHELL_CONDITION are the index variables for the parameters
*----------------------------------------------------------------------------*/
FisheriesInfo::FisheriesInfo(){
    if (debug) cout<<"starting FisheriesInfo::FisheriesInfo()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    //create "independent variables" for parameter group assignment
    nIVs = 5;//number of independent variables
    lblIVs.allocate(1,nIVs);
    k=1;
    lblIVs(k++) = tcsam::STR_FISHERY;
    lblIVs(k++) = "YEAR_BLOCK";
    lblIVs(k++) = tcsam::STR_SEX;
    lblIVs(k++) = tcsam::STR_MATURITY_STATE;
    lblIVs(k++) = tcsam::STR_SHELL_CONDITION;
    //create index block sets for "BLOCKS" (e.g. year blocks)
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 2; //YEAR_BLOCK
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs=9;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pHM";      dscPVs(k++) = "handling mortality (0-1)";
    lblPVs(k) = "pLnC";     dscPVs(k++) = "ln-scale base mean capture rate (mature male crab)";
    lblPVs(k) = "pDC1";     dscPVs(k++) = "ln-scale capture rate offset 1";
    lblPVs(k) = "pDC2";     dscPVs(k++) = "ln-scale capture rate offset 2";
    lblPVs(k) = "pDC3";     dscPVs(k++) = "ln-scale capture rate offset 3";
    lblPVs(k) = "pDC4";     dscPVs(k++) = "ln-scale capture rate offset 4";
    lblPVs(k) = "pDevsLnC"; dscPVs(k++) = "ln-scale annual capture rate devs";
    lblPVs(k) = "pLnEffX";  dscPVs(k++) = "ln-scale effort extrapolation parameters";
    lblPVs(k) = "pLgtRet";  dscPVs(k++) = "logit-scale retained fraction parameters";
    pHM   = 0;
    pLnC  = 0;
    pDC1  = 0;
    pDC2  = 0;
    pDC3  = 0;
    pDC4  = 0;
    pDevsLnC = 0;
    pLnEffX  = 0;
    pLgtRet  = 0;
    
    //create "extra indices"
    nXIs=3;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "idx.SelFcn";
    lblXIs(k++) = "idx.RetFcn";
    lblXIs(k++) = "useEffX";
    
    if (debug) cout<<"finished FisheriesInfo::FisheriesInfo()"<<endl;
}

FisheriesInfo::~FisheriesInfo(){
    if (pHM)  delete pHM;   pHM=0;
    if (pLnC) delete pLnC;  pLnC=0;
    if (pDC1) delete pDC1;  pDC1=0;
    if (pDC2) delete pDC2;  pDC2=0;
    if (pDC3) delete pDC3;  pDC3=0;
    if (pDC4) delete pDC4;  pDC4=0;
    if (pDC4) delete pDC4;  pDC4=0;
    if (pDevsLnC) delete pDevsLnC; pDevsLnC=0;
    if (pLnEffX)  delete pLnEffX;  pLnEffX=0;
    if (pLgtRet)  delete pLgtRet;  pLgtRet=0;
}

void FisheriesInfo::read(cifstream & is){
    if (debug) cout<<"starting FisheriesInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in FisheriesInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        pHM     = ParameterGroupInfo::read(is,lblPVs(k),pHM);     
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pHM)<<endl;  k++;
        pLnC    = ParameterGroupInfo::read(is,lblPVs(k),pLnC);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnC)<<endl;  k++;
        pDC1  = ParameterGroupInfo::read(is,lblPVs(k),pDC1);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC1)<<endl;  k++;
        pDC2  = ParameterGroupInfo::read(is,lblPVs(k),pDC2);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC2)<<endl;  k++;
        pDC3  = ParameterGroupInfo::read(is,lblPVs(k),pDC3);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC3)<<endl;  k++;
        pDC4 = ParameterGroupInfo::read(is,lblPVs(k),pDC4); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDC4)<<endl;  k++;
        pDevsLnC = ParameterGroupInfo::read(is,lblPVs(k),pDevsLnC); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsLnC)<<endl;  k++;
        pLnEffX = ParameterGroupInfo::read(is,lblPVs(k),pLnEffX); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnEffX)<<endl;  k++;
        pLgtRet = ParameterGroupInfo::read(is,lblPVs(k),pLgtRet); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLgtRet)<<endl;  k++;
    } else {
        cout<<"Error reading FisheriesInfo from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    if (debug){
        cout<<"FisheriesInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished FisheriesInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void FisheriesInfo::setToWriteVectorInitialValues(bool flag){
    if (pDevsLnC){
        for (int i=1;i<=pDevsLnC->getSize();i++){
            BoundedVectorInfo* vi = (*pDevsLnC)[i];
            if (flag) vi->readVals = INT_TRUE; else vi->readVals = INT_FALSE;        
        }
    }
}
void FisheriesInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pHM)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnC)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDC1)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDC2)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDC3)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDC4)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsLnC)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnEffX)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLgtRet)<<endl;
 }

void FisheriesInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"fsh=list("<<endl;
        ParameterGroupInfo::writeToR(os);           os<<cc<<endl;
        pLnC->writeToR(os, "pLnC", indent++); os<<cc<<endl;
        pDC1->writeToR(os, "pDC1", indent++); os<<cc<<endl;
        pDC2->writeToR(os, "pDC2", indent++); os<<cc<<endl;
        pDC3->writeToR(os, "pDC3", indent++); os<<cc<<endl;
        pDC4->writeToR(os, "pDC4", indent++); os<<cc<<endl;
        pLnEffX->writeToR(os, "pLnEffX", indent++); os<<cc<<endl;
        pLgtRet->writeToR(os, "pLgtRet", indent++); os<<cc<<endl;
        pDevsLnC->writeToR(os,"pDevsLnC",indent++); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * SurveysInfo\n
 * Encapsulates the following recruitment-related parameters:\n
 *   pQ   : base q (mature males)
 *   pDQ1 : ln-scale offset 1 (e.g., main temporal offset)
 *   pDQ2 : ln-scale offset 1 (e.g., female offsets)
 *   pDQ3 : ln-scale offset 1 (e.g., immature offsets)
 *   pDQ4 : ln-scale offset 1 (e.g., female-immature offsets)
 * Notes:
 *  1. index variables for parameters:
 *      a. SURVEY
 *      b. YEAR_BLOCK
 *      c. SEX
 *      d. MATURITY
 *      e. SHELL
*----------------------------------------------------------------------------*/
adstring SurveysInfo::NAME = "surveys";
SurveysInfo::SurveysInfo(){
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    nIVs = 5;
    lblIVs.allocate(1,nIVs);
    k=1;
    lblIVs(k++) = tcsam::STR_SURVEY;
    lblIVs(k++) = "YEAR_BLOCK";
    lblIVs(k++) = tcsam::STR_SEX;
    lblIVs(k++) = tcsam::STR_MATURITY_STATE;
    lblIVs(k++) = tcsam::STR_SHELL_CONDITION;
    
    nIBSs = 1;
    ibsIdxs.allocate(1,nIBSs);
    ibsIdxs(1) = 2;//index for YEAR_BLOCK
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs=5;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pQ";   dscPVs(k++) = "base catchability (e.g. mature male crab)";
    lblPVs(k) = "pDQ1"; dscPVs(k++) = "offset 1 for ln-scale catchability";
    lblPVs(k) = "pDQ2"; dscPVs(k++) = "offset 2 for ln-scale catchability";
    lblPVs(k) = "pDQ3"; dscPVs(k++) = "offset 3 for ln-scale catchability";
    lblPVs(k) = "pDQ4"; dscPVs(k++) = "offset 4 for ln-scale catchability";
    pQ   = 0;
    pDQ1 = 0;
    pDQ2 = 0;
    pDQ3 = 0;
    pDQ4 = 0;
    
    nXIs=2;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "idx.AvlFcn";
    lblXIs(k++) = "idx.SelFcn";
    
}

SurveysInfo::~SurveysInfo(){
    if (pQ)   delete pQ;   pQ   =0;
    if (pDQ1) delete pDQ1; pDQ1 =0;
    if (pDQ2) delete pDQ2; pDQ2 =0;
    if (pDQ3) delete pDQ3; pDQ3 =0;
    if (pDQ4) delete pDQ4; pDQ4 =0;
}

void SurveysInfo::read(cifstream & is){
    if (debug) cout<<"starting SurveysInfo::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in SurveysInfo::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    ParameterGroupInfo::read(is);
    
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        pQ = ParameterGroupInfo::read(is,lblPVs(k),pQ);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pQ)<<endl;  k++;
        pDQ1 = ParameterGroupInfo::read(is,lblPVs(k),pDQ1);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ1)<<endl;  k++;
        pDQ2 = ParameterGroupInfo::read(is,lblPVs(k),pDQ2);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ2)<<endl;  k++;
        pDQ3 = ParameterGroupInfo::read(is,lblPVs(k),pDQ3);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ3)<<endl;  k++;
        pDQ4 = ParameterGroupInfo::read(is,lblPVs(k),pDQ4); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDQ4)<<endl;  k++;
    } else {
        cout<<"Error reading SurveysInfo from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    
    if (debug){
        cout<<"SurveysInfo object "<<this<<endl<<(*this)<<endl;
        cout<<"finished SurveysInfo::read(cifstream & is)"<<endl;
        cout<<"Enter 1 to continue: ";
        cin>>debug;
        if (debug<0) exit(1);
    }
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void SurveysInfo::setToWriteVectorInitialValues(bool flag){
    //do nothing: no parameter vectors defined
}

void SurveysInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pQ)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDQ1)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDQ2)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDQ3)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDQ4)<<endl;
 }

void SurveysInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"srv=list("<<endl;
        ParameterGroupInfo::writeToR(os);   os<<cc<<endl;
        pQ  ->writeToR(os,"pQ",indent++);   os<<cc<<endl;
        pDQ1->writeToR(os,"pDQ1",indent++); os<<cc<<endl;
        pDQ2->writeToR(os,"pDQ2",indent++); os<<cc<<endl;
        pDQ3->writeToR(os,"pDQ3",indent++); os<<cc<<endl;
        pDQ4->writeToR(os,"pDQ4",indent++); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * MSE_Info
 -----------------------------------------------------------------------------*/
adstring MSE_Info::NAME = "MSE";
int MSE_Info::idxMSE_LnC  = 1+1;//column in parameter combinations matrix with parameter index for capture rate for MSE op mod 
/*------------------------------------------------------------------------------
 * MSE_Info\n
 * Encapsulates the following MSE-related parameters:\n
 *   pMSE_LnC : ln-scale capture rate for MSE op mod   
 *
 * Notes:
 *  1. FISHERY is the only index variable for the parameters
*----------------------------------------------------------------------------*/
MSE_Info::MSE_Info(){
    if (debug) cout<<"starting MSE_Info::MSE_Info()"<<endl;
    name = NAME;//assign static NAME to ParameterGroupInfo::name
    
    int k;
    //create "independent variables" for parameter group assignment
    nIVs = 1;//number of independent variables
    lblIVs.allocate(1,nIVs);
    k=1;
    lblIVs(k++) = tcsam::STR_FISHERY;
//    lblIVs(k++) = "YEAR_BLOCK";
//    lblIVs(k++) = tcsam::STR_SEX;
//    lblIVs(k++) = tcsam::STR_MATURITY_STATE;
//    lblIVs(k++) = tcsam::STR_SHELL_CONDITION;
    //create index block sets for "BLOCKS" (e.g. year blocks)
    nIBSs = 0;
//    ibsIdxs.allocate(1,nIBSs);
//    ibsIdxs(1) = 2; //YEAR_BLOCK
    ParameterGroupInfo::createIndexBlockSets();
    for (int i=1;i<=nIBSs;i++) ppIBSs[i-1]->setType(lblIVs(ibsIdxs(i)));
    
    nPVs=1;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pMSE_LnC"; dscPVs(k++) = "ln-scale capture rate for directed fishery in MSE op mod";
    pMSE_LnC = 0;
    
    //create "extra indices"
    nXIs=0;
//    lblXIs.allocate(1,nXIs);//
//    k=1;
//    lblXIs(k++) = "idx.SelFcn";
//    lblXIs(k++) = "idx.RetFcn";
    
    if (debug) cout<<"finished MSE_Info::MSE_Info()"<<endl;
}

MSE_Info::~MSE_Info(){
    if (pMSE_LnC) delete pMSE_LnC;  pMSE_LnC=0;
}

void MSE_Info::read(cifstream & is){
    if (debug) cout<<"starting MSE_Info::read(cifstream & is)"<<endl;
    
    adstring str;
    is>>str;
    rpt::echo<<str<<tb<<"#Required keyword ("<<NAME<<")"<<endl;
    if (!(str==NAME)){
        cout<<"Error in MSE_Info::read(cifstream & is)"<<endl;
        tcsam::readError(is,NAME,str);
        exit(-1);
    }
    if (debug) {
        cout<<"MSE_Info starting ParameterGroupInfo::read(is)"<<endl;
        ParameterGroupInfo::debug=1;
    }
    ParameterGroupInfo::read(is);
    if (debug) {
        ParameterGroupInfo::debug=0;
        cout<<"MSE_Info finished ParameterGroupInfo::read(is)"<<endl;
    }
    
    is>>str;
    if (debug) cout<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETERS)"<<endl;
    if (str=="PARAMETERS"){
        int k=1;
        if (debug) cout<<"MSE_Info::read--starting PGI::read(is,lblPVs(k),pMSE_LnC"<<endl;
        pMSE_LnC = ParameterGroupInfo::read(is,lblPVs(k),pMSE_LnC);    
        if (debug) {
            cout<<"MSE_Info::read--finished PGI::read(is,lblPVs(k),pMSE_LnC"<<endl;
            cout<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; cout<<(*pMSE_LnC)<<endl;
        }
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pMSE_LnC)<<endl;  k++;
    } else {
        cout<<"Error reading MSE_Info from "<<is.get_file_name()<<endl;
        cout<<"Expected keyword 'PARAMETERS' but got '"<<str<<"'."<<endl;
        cout<<"Aborting..."<<endl;
        exit(-1);
    }
    if (debug) cout<<"finished MSE_Info::read(is)"<<endl;
}

/**
 * Sets the flags to write initial values for vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void MSE_Info::setToWriteVectorInitialValues(bool flag){
    //do nothing: no parameter vectors defined
}

void MSE_Info::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pMSE_LnC)<<endl;
 }

void MSE_Info::writeToR(std::ostream & os){
    int indent=0;
    os<<"mse=list("<<endl;
        ParameterGroupInfo::writeToR(os);         os<<cc<<endl;
        pMSE_LnC->writeToR(os, "pMSE_LnC", indent++); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * ModelParametersInfo
 -----------------------------------------------------------------------------*/
ModelParametersInfo::ModelParametersInfo(ModelConfiguration& mc){
    ptrMC=&mc;
    ptrRec=0;
    ptrNM =0;
    ptrGrw=0;
    ptrM2M=0;
    ptrSel=0;
    ptrFsh=0;
    ptrSrv=0;
    ptrMSE=0;
}

ModelParametersInfo::~ModelParametersInfo(){
    if (ptrRec) {delete ptrRec;     ptrRec  = 0;}
    if (ptrNM)  {delete ptrNM;      ptrNM   = 0;}
    if (ptrGrw) {delete ptrGrw;     ptrGrw  = 0;}
    if (ptrM2M) {delete ptrM2M;     ptrM2M  = 0;}
    if (ptrSel) {delete ptrSel;     ptrSel  = 0;}
    if (ptrFsh) {delete ptrFsh;     ptrFsh  = 0;}
    if (ptrSrv) {delete ptrSrv;     ptrSrv  = 0;}
    if (ptrMSE) {delete ptrMSE;     ptrMSE  = 0;}
}

/**
 * Sets the flags to write initial values for all vector parameters to file 
 * when writing parameter info to file.
 * 
 * @param flag - true/false to set to write initial values to file
 */
void ModelParametersInfo::setToWriteVectorInitialValues(bool flag){
    if (ptrRec) ptrRec->setToWriteVectorInitialValues(flag);
    if (ptrNM)  ptrNM ->setToWriteVectorInitialValues(flag);
    if (ptrGrw) ptrGrw->setToWriteVectorInitialValues(flag);
    if (ptrM2M) ptrM2M->setToWriteVectorInitialValues(flag);
    if (ptrSel) ptrSel->setToWriteVectorInitialValues(flag);
    if (ptrFsh) ptrFsh->setToWriteVectorInitialValues(flag);
    if (ptrSrv) ptrSrv->setToWriteVectorInitialValues(flag); //not needed for SurveysInfo
    if (ptrMSE) ptrMSE->setToWriteVectorInitialValues(flag); //not needed for MSE_Info
}

void ModelParametersInfo::read(cifstream & is){
    if (debug) cout<<"starting void ModelParametersInfo::read(cifstream & is)"<<endl;
    
    rpt::echo<<"#------Reading Model Parameters Info file-----"<<endl;
    rpt::echo<<"file is '"<<is.get_file_name()<<"'."<<endl;
    
    adstring str;
    is>>str;
    if (str!=ModelParametersInfo::version){
        cout<<"-----------------------------------------------------------"<<endl;
        cout<<"Model Parameters Info file '"<<is.get_file_name()<<"'"<<endl;
        cout<<"has incorrect version number!!!"<<endl;
        cout<<"Expected '"<<ModelParametersInfo::version<<"'. Got '"<<str<<"'."<<endl;
        cout<<"Terminating run..."<<endl;
        exit(-1);
    }
    
    //read recruitment parameters
    rpt::echo<<"#---reading  Recruitment Info"<<endl;
    ptrRec = new RecruitmentInfo();
    is>>(*ptrRec);
    rpt::echo<<"#---created  RecruitmentInfo object"<<endl;
    
    //read natural mortality parameters
    rpt::echo<<"#---reading Natural Mortality Info"<<endl;
    ptrNM = new NaturalMortalityInfo();
    is>>(*ptrNM);
    if (debug) cout<<"created NaturalMortalityInfo object"<<endl;
    
    //read growth parameters
    rpt::echo<<"#---reading Growth Info"<<endl;
    ptrGrw = new GrowthInfo();
    is>>(*ptrGrw);
    if (debug) cout<<"created GrowthInfo object"<<endl;
    
    //read maturity parameters
    rpt::echo<<"#---reading Molt2Maturity Info"<<endl;
    ptrM2M = new Molt2MaturityInfo();
    is>>(*ptrM2M);
    if (debug) cout<<"created Molt2MaturityInfo object"<<endl;
    
    //read selectivity function parameters
    rpt::echo<<"#---reading Selectivity Info"<<endl;
    ptrSel = new SelectivityInfo();
    is>>(*ptrSel);
    if (debug) cout<<"created SelectivityInfo object"<<endl;
    
    //read fisheries parameters
    rpt::echo<<"#---reading Fisheries Info"<<endl;
    ptrFsh = new FisheriesInfo();
    is>>(*ptrFsh);
    if (debug) cout<<"created FisheriesInfo object"<<endl;
    
    //read surveys parameters
    rpt::echo<<"#---reading Surveys Info"<<endl;
    ptrSrv = new SurveysInfo();
    is>>(*ptrSrv);
    rpt::echo<<"#---read Surveys Info"<<endl;
    
    //read MSE-related parameters
    rpt::echo<<"#---reading MSE Info"<<endl;
    ptrMSE = new MSE_Info();
    is>>(*ptrMSE);
    rpt::echo<<"#---read MSE Info"<<endl;
    
    if (debug) cout<<"finished void ModelParametersInfo::read(cifstream & is)"<<endl;
}
    
void ModelParametersInfo::write(std::ostream & os){
    os<<version<<tb<<"#model parameters info version"<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Recruitment parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrRec)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Natural mortality parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrNM)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Growth parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrGrw)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Molt-to-maturity parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrM2M)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Selectivity parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrSel)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Fisheries parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrFsh)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# Surveys parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrSrv)<<endl;
    
    os<<"#-------------------------------"<<endl;
    os<<"# MSE-related parameters  "<<endl;
    os<<"#-------------------------------"<<endl;
    os<<(*ptrMSE)<<endl;
}

void ModelParametersInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"mpi=list(version='"<<version<<"'"<<cc<<endl;
    ptrRec->writeToR(os); os<<cc<<endl;
    ptrNM->writeToR(os);  os<<cc<<endl;
    ptrGrw->writeToR(os); os<<cc<<endl;
    ptrM2M->writeToR(os); os<<cc<<endl;
    ptrSel->writeToR(os); os<<cc<<endl;
    ptrFsh->writeToR(os); os<<cc<<endl;
    ptrSrv->writeToR(os); os<<endl;
    ptrMSE->writeToR(os); os<<endl;
    os<<")";
}
        
