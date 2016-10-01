//--------------------------------------------------------------------------------
//  Includes:
//      ModelParametersInfo
//--------------------------------------------------------------------------------
#pragma once
#ifndef MODELPARAMETERSCONFIGURATION_HPP
    #define MODELPARAMETERSCONFIGURATION_HPP

    #include <admodel.h>

#include "ModelProcesses.hpp"

    //forward declarations
    class ModelConfiguration;
    class SexRatio;
    class Recruitment;
    class Maturity;
    class NaturalMortality;

//--------------------------------------------------------------------------------
//          ModelParametersInfo
//  Encapsulates definitions for all model parameters
//--------------------------------------------------------------------------------
    class ModelParametersConfiguration {
        public:
            static int debug;
        protected:
            int nSrv;               //number of surveys defined
            int incFshTnr;          //include directed Tanner fishery?
            int incFshOpi;          //include opilio fishery bycatch?
            int incFshRKC;          //include red king crab fishery bycatch?
            int incFshGrf;          //include groundfish trawl fishery bycatch?

        public:
            adstring cfgName; //model parameters configuration name

        public:
            ModelConfiguration*  ptrMC;  ///ptr to model configuration

            SexRatio*            ptrSxR;//sex ratio
            Recruitment*         ptrRec;//recruitment            
            Maturity*            ptrMat;//maturity
            NaturalMortality*    ptrNM; //natural mortality           
            
        public:
            ModelParametersConfiguration(ModelConfiguration & modCon);
            ~ModelParametersConfiguration();
            void read(void);
            void write(const adstring& fn);
            void write(ostream & os);
            void writeInfoToR(ostream& os, adstring nm, int indent=0);
            void writeValuesToR(ostream& os, adstring nm, int indent=0);
            void writeResultsToR(ostream& os, adstring nm, int indent=0);
            friend ostream& operator <<(ostream & os, ModelParametersConfiguration & obj){obj.write(os); return os;}

        private:
            void read(cifstream & is);
    };
#endif //MODELPARAMETERSINFO_HPP
