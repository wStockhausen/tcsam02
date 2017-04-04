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
int ModelParametersInfo::debug  = 0;
const adstring ModelParametersInfo::version = "2017.04.03";
    
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
        if (i==pgi->ibsIdxs(ibsIdx)){
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
            if (i==pgi->ibsIdxs(ibsIdx)){
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
/********************************************\n
 * Gets indices for parameter combination pc.\n
 * pc : id for desired parameter combination.\n
 *******************************************/
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
 * Reads info for a BoundedNumberVector object from an input  filestream.
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
    rpt::echo<<str<<tb<<"#Required keyword (PARAMETER_COMBINATIONS)"<<endl;
    if (str=="PARAMETER_COMBINATIONS"){
        is>>nPCs;
        rpt::echo<<nPCs<<tb<<"#number of parameter combinations"<<endl;
        rpt::echo<<"#id  "; 
        for (int i=1;i<=nIVs;i++) rpt::echo<<lblIVs(i)<<tb; 
        for (int i=1;i<=nPVs;i++) rpt::echo<<lblPVs(i)<<tb; 
        for (int i=1;i<=nXIs;i++) rpt::echo<<lblXIs(i)<<tb; 
        rpt::echo<<endl;
        for (int i=0;i<nIBSs;i++) ppIBSs[i]->allocate(nPCs);
        //read parameters combinations definition matrix
        int ibsIdx=1;
        in.allocate(1,nPCs,1,nIVs+nPVs+nXIs);
        if (nXIs) xd.allocate(1,nPCs,1,nXIs);
        for (int r=1;r<=nPCs;r++){//loop over rows
            is>>str; //read id
            rpt::echo<<str<<tb;
            for (int i=1;i<=nIVs;i++){//loop over index variables
                is>>str;
                rpt::echo<<str<<tb;
                if (lblIVs(i)==tcsam::STR_SEX)             {in(r,i) = tcsam::getSexType(str);}      else
                if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {in(r,i) = tcsam::getMaturityType(str);} else
                if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {in(r,i) = tcsam::getShellType(str);}    else
                if (i==ibsIdxs(ibsIdx)){//variable is a block
                    ppIBSs[ibsIdx-1]->getIndexBlock(r)->parse(str);
                    in(r,i)=ibsIdx;//index to associated IndexBlockSet
                    if (ibsIdx<nIBSs) ibsIdx++;//increment to next IBS
                } else {in(r,i)=::atoi(str);}
            }
            for (int p=1;p<=nPVs;p++) {is>>in(r,nIVs+p); rpt::echo<<in(r,nIVs+p)<<tb;}  
            //rpt::echo<<"looping over extra variables"<<endl;
            for (int x=1;x<=nXIs;x++) {//loop over "extra" variables
                is>>str;
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
            rpt::echo<<endl;
            if (debug) cout<<"pc row "<<r<<": "<<in(r)<<tb<<"xd: "<<xd(r)<<endl;
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

void ParameterGroupInfo::write(std::ostream& os){
//    os<<(*ptrIBSs)<<endl;
    os<<"PARAMETER_COMBINATIONS"<<endl;
    os<<nPCs<<tb<<"#number of parameter combinations"<<endl;
    os<<"#id  "; 
    for (int i=1;i<=nIVs;i++) os<<lblIVs(i)<<tb; 
    for (int i=1;i<=nPVs;i++) os<<lblPVs(i)<<tb; 
    for (int i=1;i<=nXIs;i++) os<<lblXIs(i)<<tb; 
    os<<endl;
    for (int r=1;r<=nPCs;r++){//loop over rows
        os<<r<<tb;
        int ibsIdx = 1;
        for (int i=1;i<=nIVs;i++){//loop over index variables           
            if (lblIVs(i)==tcsam::STR_SEX)             {os<<tcsam::getSexType(in(r,i))<<tb;} else
            if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {os<<tcsam::getMaturityType(in(r,i))<<tb;} else
            if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {os<<tcsam::getShellType(in(r,i))<<tb;} else 
            if (i==ibsIdxs(ibsIdx)){
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
        os<<endl;
    }
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
//        if (nIVs) {os<<"IVs=c("<<wts::to_qcsv(lblIVs)<<")"<<cc;} else {os<<"IVs=NULL"<<cc;}
//        if (nPVs) {os<<"PVs=c("<<wts::to_qcsv(lblPVs)<<")"<<cc;} else {os<<"PVs=NULL"<<cc;}
//        if (nXIs) {os<<"XIs=c("<<wts::to_qcsv(lblXIs)<<")"<<cc;} else {os<<"XIs=NULL"<<cc;}
//        os<<"nPVs="<<nPVs<<cc;
//        os<<"nXIs="<<nXIs<<cc;
//        os<<"nPCs="<<nPCs<<cc;
//        os<<endl;
        os<<"pcs=list("<<endl;
            for (int p=1;p<=nPCs;p++){
                os<<qt<<p<<qt<<"=list(";
                int ibsIdx = 1;
                for (int i=1;i<=nIVs;i++){//loop over index variables
                    os<<lblIVs(i)<<"='";
                    if (lblIVs(i)==tcsam::STR_SEX)             {os<<tcsam::getSexType(in(p,i))<<"',";} else
                    if (lblIVs(i)==tcsam::STR_MATURITY_STATE)  {os<<tcsam::getMaturityType(in(p,i))<<"',";} else
                    if (lblIVs(i)==tcsam::STR_SHELL_CONDITION) {os<<tcsam::getShellType(in(p,i))<<"',";} else 
                    if (i==ibsIdxs(ibsIdx)){
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
    lblPVs(k) = "pLnRCV";  dscPVs(k++) = "recruitment cv's";
    lblPVs(k) = "pLgtRX";  dscPVs(k++) = "logit-scale male sex ratio";
    lblPVs(k) = "pLnRa";   dscPVs(k++) = "size-at-recruitment parameter";
    lblPVs(k) = "pLnRb";   dscPVs(k++) = "size-at-recruitment parameter";    
    lblPVs(k) = "pDevsLnR";dscPVs(k++) = "ln-scale recruitment devs";    
    pLnR     = 0;
    pLnRCV   = 0;
    pLgtRX   = 0;
    pLnRa    = 0;
    pLnRb    = 0;
    pDevsLnR = 0;
    
    nXIs = 0;    
}

RecruitmentInfo::~RecruitmentInfo(){
    if (pLnR)     delete pLnR;     pLnR     = 0;
    if (pLnRCV)   delete pLnRCV;   pLnRCV   = 0;
    if (pLgtRX)   delete pLgtRX;   pLgtRX   = 0;
    if (pLnRa)    delete pLnRa;    pLnRa    = 0;
    if (pLnRb)    delete pLnRb;    pLnRb    = 0;
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
        pLnRCV = ParameterGroupInfo::read(is,lblPVs(k),pLnRCV); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnRCV)<<endl; k++;
        pLgtRX = ParameterGroupInfo::read(is,lblPVs(k),pLgtRX); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLgtRX)<<endl; k++;
        pLnRa  = ParameterGroupInfo::read(is,lblPVs(k),pLnRa);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnRa)<<endl;  k++;
        pLnRb  = ParameterGroupInfo::read(is,lblPVs(k),pLnRb);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnRb)<<endl;  k++;
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

void RecruitmentInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnR)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnRCV)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLgtRX)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnRa)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnRb)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsLnR)<<endl;
}

void RecruitmentInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"rec=list("<<endl;
        ParameterGroupInfo::writeToR(os);           os<<cc<<endl;
        pLnR->writeToR(os,  "pLnR",  indent+1);     os<<cc<<endl;
        pLnRCV->writeToR(os,"pLnRCV",indent+1);     os<<cc<<endl;
        pLgtRX->writeToR(os,"pLgtRX",indent+1);     os<<cc<<endl;
        pLnRa->writeToR(os, "pLnRa", indent+1);     os<<cc<<endl;
        pLnRb->writeToR(os, "pLnRb", indent+1);     os<<cc<<endl;
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
    lblPVs(k) = "pLnM";      dscPVs(k++) = "ln-scale base natural mortality rate (immature male crab)";
    lblPVs(k) = "pLnDMT";    dscPVs(k++) = "main temporal ln-scale natural mortality offsets";
    lblPVs(k) = "pLnDMX";    dscPVs(k++) = "ln-scale natural mortality offset for female crabs";
    lblPVs(k) = "pLnDMM";    dscPVs(k++) = "ln-scale natural mortality offset for mature crabs";
    lblPVs(k) = "pLnDMXM";   dscPVs(k++) = "ln-scale natural mortality offset for mature female crabs";
    if (debug) cout<<3<<endl;
    pLnM    = 0;
    pLnDMT  = 0;
    pLnDMX  = 0;
    pLnDMM  = 0;
    pLnDMXM = 0;
    if (debug) cout<<4<<endl;
    
    //define "extra" indices
    nXIs=1;
    lblXIs.allocate(1,nXIs);
    lblXIs(k=1) = "zScaling";    
    if (debug) cout<<5<<endl;
    if (debug) cout<<"finished NaturalMortalityInfo::NaturalMortalityInfo()"<<endl;
}

NaturalMortalityInfo::~NaturalMortalityInfo(){
    if (pLnM)    delete pLnM;     pLnM   =0;
    if (pLnDMT)  delete pLnDMT;   pLnDMT =0;
    if (pLnDMX)  delete pLnDMX;   pLnDMX =0;
    if (pLnDMM)  delete pLnDMM;   pLnDMM =0;
    if (pLnDMXM) delete pLnDMXM;  pLnDMXM=0;
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
        pLnM    = ParameterGroupInfo::read(is,lblPVs(k),pLnM);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnM)<<endl;    k++;
        pLnDMT  = ParameterGroupInfo::read(is,lblPVs(k),pLnDMT);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDMT)<<endl;  k++;
        pLnDMX  = ParameterGroupInfo::read(is,lblPVs(k),pLnDMX);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDMX)<<endl;  k++;
        pLnDMM  = ParameterGroupInfo::read(is,lblPVs(k),pLnDMM);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDMM)<<endl;  k++;
        pLnDMXM = ParameterGroupInfo::read(is,lblPVs(k),pLnDMXM); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDMXM)<<endl; k++;
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

void NaturalMortalityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    os<<zRef<<tb<<"#reference size for scaling"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnM)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDMT)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDMX)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDMM)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDMXM)<<endl;
 }

void NaturalMortalityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"nm=list("<<endl;
        ParameterGroupInfo::writeToR(os);         os<<cc<<endl;
        pLnM->writeToR(os,   "pLnM",   indent+1); os<<cc<<endl;
        pLnDMT->writeToR(os, "pLnDXT", indent+1); os<<cc<<endl;
        pLnDMX->writeToR(os, "pLnDMX", indent+1); os<<cc<<endl;
        pLnDMM->writeToR(os, "pLnDMM", indent+1); os<<cc<<endl;
        pLnDMXM->writeToR(os,"pLnDMXM",indent+1); os<<endl;
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
    lblPVs(k) = "pLnGrA";    dscPVs(k++) = "ln-scale mean growth coefficient 'a'";
    lblPVs(k) = "pLnGrB";    dscPVs(k++) = "ln-scale mean growth coefficient 'b'";
    lblPVs(k) = "pLnGrBeta"; dscPVs(k++) = "ln-scale growth transition matrix scale factor";
    pLnGrA    = 0;
    pLnGrB    = 0;
    pLnGrBeta = 0;
    
    nXIs=0;
//    lblXIs.allocate(1,nXIs);
    
}

GrowthInfo::~GrowthInfo(){
    if (pLnGrA)    delete pLnGrA;    pLnGrA    =0;
    if (pLnGrB)    delete pLnGrB;    pLnGrB    =0;
    if (pLnGrBeta) delete pLnGrBeta; pLnGrBeta =0;
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
        pLnGrA     = ParameterGroupInfo::read(is,lblPVs(k),pLnGrA);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnGrA)<<endl;    k++;
        pLnGrB     = ParameterGroupInfo::read(is,lblPVs(k),pLnGrB);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnGrB)<<endl;    k++;
        pLnGrBeta  = ParameterGroupInfo::read(is,lblPVs(k),pLnGrBeta); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnGrBeta)<<endl; k++;
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

void GrowthInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
     ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnGrA)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnGrB)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnGrBeta)<<endl;
 }

void GrowthInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"grw=list("<<endl;
        ParameterGroupInfo::writeToR(os);             os<<cc<<endl;
        pLnGrA->writeToR(os,   "pLnGrA",   indent+1); os<<cc<<endl;
        pLnGrB->writeToR(os,   "pLnGrB",   indent+1); os<<cc<<endl;
        pLnGrBeta->writeToR(os,"pLnGrBeta",indent+1); os<<endl;
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
    pLgtPrM2M = 0;
    
    nXIs = 0;    
}

Molt2MaturityInfo::~Molt2MaturityInfo(){
    if (pLgtPrM2M) delete pLgtPrM2M; pLgtPrM2M = 0;
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
        pLgtPrM2M = ParameterGroupInfo::read(is,lblPVs(k),pLgtPrM2M);
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLgtPrM2M)<<endl;  
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

void Molt2MaturityInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLgtPrM2M)<<endl;
}

void Molt2MaturityInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"mat=list("<<endl;
        ParameterGroupInfo::writeToR(os);             os<<cc<<endl;
        pLgtPrM2M->writeToR(os,"pLgtPrM2M",indent+1); os<<endl;
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
    
    nPVs=12;
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
        pDevsS6->writeToR(os,"pDevsS6",indent+1); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * FisheriesInfo
 -----------------------------------------------------------------------------*/
adstring FisheriesInfo::NAME = "fisheries";
/*------------------------------------------------------------------------------
 * FisheriesInfo\n
 * Encapsulates the following fishery-related parameters:\n
 *   pHM    : handling mortality (0-1)
 *   pLnC   : ln-scale base mean capture rate (mature males)
 *   pLnCT  : main year_block offset
 *   pLnDCX : female offsets
 *   pLnDCM : immature offsets
 *   pLnDCXM: female-immature offsets    
 *
 *   pDevsLnC : annual ln-scale devs w/in year_blocks
 * 
 *   pLnEffX: ln-scale effort extrapolation parameters
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
    
    nPVs=8;
    lblPVs.allocate(1,nPVs); dscPVs.allocate(1,nPVs);
    k=1;
    lblPVs(k) = "pHM";       dscPVs(k++) = "handling mortality (0-1)";
    lblPVs(k) = "pLnC";      dscPVs(k++) = "ln-scale base mean capture rate (mature male crab)";
    lblPVs(k) = "pLnDCT";    dscPVs(k++) = "main year_block offset for ln-scale capture rate";
    lblPVs(k) = "pLnDCX";    dscPVs(k++) = "ln-scale capture rate offset for female crabs";
    lblPVs(k) = "pLnDCM";    dscPVs(k++) = "ln-scale capture rate offset for immature crabs";
    lblPVs(k) = "pLnDCXM";   dscPVs(k++) = "ln-scale capture rate offset for immature female crabs";
    lblPVs(k) = "pDevsLnC";  dscPVs(k++) = "ln-scale annual capture rate devs";
    lblPVs(k) = "pLnEffX";   dscPVs(k++) = "ln-scale effort extrapolation parameters";
    pHM     = 0;
    pLnC    = 0;
    pLnDCT  = 0;
    pLnDCX  = 0;
    pLnDCM  = 0;
    pLnDCXM = 0;
    pDevsLnC = 0;
    pLnEffX  = 0;
    
    //create "extra indices"
    nXIs=3;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "idx.SelFcn";
    lblXIs(k++) = "idx.RetFcn";
    lblXIs(k++) = "useEffX";
    
    idxHM    = nIVs+1;     //1st PV column
    idxLnEX  = nIVs+nPVs;  //last PV column
    idxUseEX = nIVs+nPVs+3;//column in parameter combinations matrix for useEffortExtrapolation flag
    
    if (debug) cout<<"finished FisheriesInfo::FisheriesInfo()"<<endl;
}

FisheriesInfo::~FisheriesInfo(){
    if (pHM)     delete pHM;      pHM=0;
    if (pLnC)    delete pLnC;     pLnC   =0;
    if (pLnDCT)  delete pLnDCT;   pLnDCT =0;
    if (pLnDCX)  delete pLnDCX;   pLnDCX =0;
    if (pLnDCM)  delete pLnDCM;   pLnDCM =0;
    if (pLnDCXM) delete pLnDCXM;  pLnDCXM=0;
    if (pLnDCXM) delete pLnDCXM;  pLnDCXM=0;
    if (pLnEffX) delete pLnEffX;  pLnEffX=0;
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
        pLnDCT  = ParameterGroupInfo::read(is,lblPVs(k),pLnDCT);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDCT)<<endl;  k++;
        pLnDCX  = ParameterGroupInfo::read(is,lblPVs(k),pLnDCX);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDCX)<<endl;  k++;
        pLnDCM  = ParameterGroupInfo::read(is,lblPVs(k),pLnDCM);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDCM)<<endl;  k++;
        pLnDCXM = ParameterGroupInfo::read(is,lblPVs(k),pLnDCXM); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDCXM)<<endl;  k++;
        pDevsLnC = ParameterGroupInfo::read(is,lblPVs(k),pDevsLnC); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pDevsLnC)<<endl;  k++;
        pLnEffX = ParameterGroupInfo::read(is,lblPVs(k),pLnEffX); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnEffX)<<endl;  k++;
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
    os<<(*pLnDCT)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDCX)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDCM)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDCXM)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pDevsLnC)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnEffX)<<endl;
 }

void FisheriesInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"fsh=list("<<endl;
        ParameterGroupInfo::writeToR(os);           os<<cc<<endl;
        pLnC->writeToR(os,    "pLnC",    indent++); os<<cc<<endl;
        pLnDCT->writeToR(os,  "pLnDXT",  indent++); os<<cc<<endl;
        pLnDCX->writeToR(os,  "pLnDCX",  indent++); os<<cc<<endl;
        pLnDCM->writeToR(os,  "pLnDCM",  indent++); os<<cc<<endl;
        pLnDCXM->writeToR(os, "pLnDCXM", indent++); os<<cc<<endl;
        pLnEffX->writeToR(os, "pLnEffX", indent++); os<<cc<<endl;
        pDevsLnC->writeToR(os,"pDevsLnC",indent++); os<<endl;
    os<<")";
}

/*------------------------------------------------------------------------------
 * SurveysInfo
 -----------------------------------------------------------------------------*/
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
    lblPVs(k) = "pLnQ";      dscPVs(k++) = "ln-scale base catchability (mature male crab)";
    lblPVs(k) = "pLnDQT";    dscPVs(k++) = "main temporal offset for ln-scale catchability";
    lblPVs(k) = "pLnDQX";    dscPVs(k++) = "ln-scale catchability offset for female crabs";
    lblPVs(k) = "pLnDQM";    dscPVs(k++) = "ln-scale catchability offset for immature crabs";
    lblPVs(k) = "pLnDQXM";   dscPVs(k++) = "ln-scale catchability offset for immature female crabs";
    pLnQ    = 0;
    pLnDQT  = 0;
    pLnDQX  = 0;
    pLnDQM  = 0;
    pLnDQXM = 0;
    
    nXIs=2;
    lblXIs.allocate(1,nXIs);
    k=1;
    lblXIs(k++) = "idx.AvlFcn";
    lblXIs(k++) = "idx.SelFcn";
    
}

SurveysInfo::~SurveysInfo(){
    if (pLnQ)    delete pLnQ;     pLnQ   =0;
    if (pLnDQT)  delete pLnDQT;   pLnDQT =0;
    if (pLnDQX)  delete pLnDQX;   pLnDQX =0;
    if (pLnDQM)  delete pLnDQM;   pLnDQM =0;
    if (pLnDQXM) delete pLnDQXM;  pLnDQXM=0;
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
        pLnQ    = ParameterGroupInfo::read(is,lblPVs(k),pLnQ);    
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnQ)<<endl;  k++;
        pLnDQT  = ParameterGroupInfo::read(is,lblPVs(k),pLnDQT);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDQT)<<endl;  k++;
        pLnDQX  = ParameterGroupInfo::read(is,lblPVs(k),pLnDQX);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDQX)<<endl;  k++;
        pLnDQM  = ParameterGroupInfo::read(is,lblPVs(k),pLnDQM);  
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDQM)<<endl;  k++;
        pLnDQXM = ParameterGroupInfo::read(is,lblPVs(k),pLnDQXM); 
        rpt::echo<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; rpt::echo<<(*pLnDQXM)<<endl;  k++;
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

void SurveysInfo::write(std::ostream & os){
    os<<NAME<<tb<<"#process name"<<endl;
    ParameterGroupInfo::write(os);
    
    os<<"PARAMETERS"<<endl;
    
    int k=1;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnQ)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDQT)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDQX)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDQM)<<endl;
    os<<lblPVs(k)<<tb<<"#"<<dscPVs(k)<<endl; k++;
    os<<(*pLnDQXM)<<endl;
 }

void SurveysInfo::writeToR(std::ostream & os){
    int indent=0;
    os<<"srv=list("<<endl;
        ParameterGroupInfo::writeToR(os);         os<<cc<<endl;
        pLnQ->writeToR(os,"pLnQ",indent++);       os<<cc<<endl;
        pLnDQT->writeToR(os,"pLnDXT",indent++);   os<<cc<<endl;
        pLnDQX->writeToR(os,"pLnDQX",indent++);   os<<cc<<endl;
        pLnDQM->writeToR(os,"pLnDQM",indent++);   os<<cc<<endl;
        pLnDQXM->writeToR(os,"pLnDQXM",indent++); os<<endl;
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
}

ModelParametersInfo::~ModelParametersInfo(){
    if (ptrRec) {delete ptrRec;     ptrRec  = 0;}
    if (ptrNM)  {delete ptrNM;      ptrNM   = 0;}
    if (ptrGrw) {delete ptrGrw;     ptrGrw  = 0;}
    if (ptrM2M) {delete ptrM2M;     ptrM2M  = 0;}
    if (ptrSel) {delete ptrSel;     ptrSel  = 0;}
    if (ptrFsh) {delete ptrFsh;     ptrFsh  = 0;}
    if (ptrSrv) {delete ptrSrv;     ptrSrv  = 0;}
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
    rpt::echo<<"#---read  Recruitment Info"<<endl;
    
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
    os<<")";
}
        
void tcsam::readError(cifstream & is, const char * expP, adstring gotP){
    cout<<"Reading "<<is.get_file_name()<<endl;
    cout<<"Expected parameter name '"<<expP<<"' but got '"<<gotP<<"'."<<endl;
    cout<<"Aborting..."<<endl;
}
        
void tcsam::setParameterInfo(NumberVectorInfo* pNVI,                           
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
  
void tcsam::setParameterInfo(BoundedNumberVectorInfo* pBNVI,
                             int& npT,
                             dvector& lb, dvector& ub, 
                             ivector& phs, 
                             ostream& os){
    int np = pBNVI->getSize();
    if (np){npT = np;} else {npT = 1;}
    phs.allocate(1,npT);
    lb.allocate(1,npT);
    ub.allocate(1,npT);
    if (np){
        phs = pBNVI->getPhases();
        lb  = pBNVI->getLowerBounds();
        ub  = pBNVI->getUpperBounds();
        os<<"parameter "<<pBNVI->name<<":"<<endl;
        os<<"#lower  upper  phase  "<<endl;
        for (int n=1;n<=np;n++) os<<n<<tb<<lb(n)<<tb<<ub(n)<<tb<<phs(n)<<endl;
    } else {
        phs = -1;
        lb  = -1;
        ub  = -1;
        os<<"bounded number parameter vector "<<pBNVI->name<<" has no parameter values"<<endl;
    }
}
  
void tcsam::setParameterInfo(VectorVectorInfo* pVVI,
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
  
void tcsam::setParameterInfo(BoundedVectorVectorInfo* pBVVI,                           
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

void tcsam::setParameterInfo(DevsVectorVectorInfo* pDVVI,                           
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
        mxs = pDVVI->getMaxIndices()-1;//max index for parameters is one less than max indices
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

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
void tcsam::writeParameter(ofstream& os, param_init_number& p, int toR, int willBeActive){
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_number'"<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"value="<<value(p)
                                <<"),";
        } else {
            os<<1<<cc<<p.get_phase_start()<<cc<<1<<cc<<1<<cc<<"-Inf"<<cc<<"Inf"<<cc<<p<<cc<<p.get_name()<<cc<<"'param_init_number'"<<endl;
        }
    }
}    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
void tcsam::writeParameter(ofstream& os, param_init_bounded_number& p,int toR, int willBeActive){
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_bounded_number'"<<cc
                                <<"phase="<<p.phase_start<<cc
                                <<"bounds=c("<<p.get_minb()<<cc<<p.get_maxb()<<")"<<cc
                                <<"value="<<value(p)
                                <<"),";
        } else {
            os<<1<<cc<<p.get_phase_start()<<cc<<1<<cc<<1<<cc<<p.get_minb()<<cc<<p.get_maxb()<<cc<<p<<cc<<p.get_name()<<cc<<"'param_init_bounded_number'"<<endl;
        }
    }
}    
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
void tcsam::writeParameter(ofstream& os, param_init_vector& p, int toR, int willBeActive){
    int mn = p.indexmin();
    int mx = p.indexmax();
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_vector'"<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"value=c("; {for (int i=mn;i<mx;i++) os<<value(p(i))<<cc;} os<<value(p(mx))<<")";
            os<<"),";
        } else {        
            for (int i=mn;i<=mx;i++) os<< i<<cc<<p.get_phase_start()<<cc<<mn<<cc<<mx<<cc<<"-Inf"<<cc<<"Inf"<<cc<<p(i)<<cc<<p.get_name()<<cc<<"'param_init_vector'"<<endl;
        }
    }
}       

void tcsam::writeParameter(ofstream& os, param_init_bounded_vector& p, int toR, int willBeActive){
    int mn = p.indexmin();
    int mx = p.indexmax();
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        if (toR){
            os<<p.get_name()<<"=list("<<"type='param_init_bounded_vector'"<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"bounds=c("<<p.get_minb()<<cc<<p.get_maxb()<<")"<<cc
                                <<"value=c("; {for (int i=mn;i<mx;i++) os<<value(p(i))<<cc;} os<<value(p(mx))<<")";
           os<<"),";
        } else {
            for (int i=mn;i<=mx;i++) os<< i<<cc<<p.get_phase_start()<<cc<<mn<<cc<<mx<<cc<<p.get_minb()<<cc<<p.get_maxb()<<cc<<p(i)<<cc<<p.get_name()<<cc<<"'param_init_bounded_vector'"<<endl;
        }
    }
}

void tcsam::writeParameterBounds(ofstream& os, param_init_bounded_dev_vector& p, int toR, int willBeActive){
    int mn = p.indexmin();
    int mx = p.indexmax();
    if (!willBeActive||(willBeActive&&(p.get_phase_start()>0))){
        if (toR){
            os<<p.get_name()<<"=list("<<"type=param_init_bounded_dev_vector"<<cc
                                <<"dims=c("<<mn<<cc<<mx<<")"<<cc
                                <<"phase="<<p.get_phase_start()<<cc
                                <<"bounds=c("<<p.get_minb()<<cc<<p.get_maxb()<<")"<<cc
                                <<"value=c("; {for (int i=mn;i<mx;i++) os<<value(p(i))<<cc;} os<<value(p(mx))<<")";
           os<<"),";
        } else {
            for (int i=mn;i<=mx;i++) os<< i<<cc<<p.get_phase_start()<<cc<<mn<<cc<<mx<<cc<<p.get_minb()<<cc<<p.get_maxb()<<cc<<p(i)<<cc<<p.get_name()<<cc<<"'param_init_bounded_dev_vector'"<<endl;
        }
    }
}    
