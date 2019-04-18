/********************************************************************************
    HypoViscosity.hpp
    nairn-mpm-fea
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		IsotropicMat.hpp (Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef _HYPOVISCOSITY_

#define _HYPOVISCOSITY_

#define HYPOVISCOSITY 110

// plastic law properties
typedef struct {
	double G;
	double K;
    double mu_s;
    double mu_2;
    double xi;
    int Is_visco;
    double volMax;
    double pMin;
    int Is_GeState;
    int Is_KeState;
} HV_Properties;

class HardeningLawBase;

#include "Materials/IsotropicMat.hpp"

class HypoViscosity : public IsotropicMat
{
    public:
        // constructors and destructors
		HypoViscosity();
		HypoViscosity(char *matName);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void SetInitialParticleState(MPMBase *,int,int) const;
	
		// history data
        virtual int SizeOfHistoryData(void) const;
        virtual char *InitHistoryData(char *,MPMBase *);
        virtual double GetHistory(int,char *) const;
 	
		// const methods
        virtual void PrintMechanicalProperties(void) const;
		
		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
		virtual void ConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;

		// custom methods: Find yield function and solve for lambda
        virtual double DoubleInner(Tensor *,Tensor *,int) const;
		virtual void GetDfDsigma(double,Tensor *,int,Tensor *) const;
		
		// accessors
        virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
        const char *MaterialType(void) const;
		virtual int AltStrainContains(void) const;
		
    protected:
		HV_Properties pr;
		double G0;
};

#endif

