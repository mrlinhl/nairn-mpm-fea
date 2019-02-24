/********************************************************************************
    DP_Plasticity.hpp
    nairn-mpm-fea
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		IsotropicMat.hpp (Elastic.hpp, MaterialBase.hpp)
********************************************************************************/

#ifndef _DP_PLASTICITY_

#define _DP_PLASTICITY_

#define DP_PLASTICITY 109

// plastic law properties
typedef struct {
	double G;
	double K;
    double aphi;
    double apsi;
    double cphi;
} DP_Properties;

class HardeningLawBase;

#include "Materials/IsotropicMat.hpp"

class DP_Plasticity : public IsotropicMat
{
    public:
        // constructors and destructors
		DP_Plasticity();
		DP_Plasticity(char *matName);
		
		// initialize
        virtual char *InputMaterialProperty(char *,int &,double &);
        virtual const char *VerifyAndLoadProperties(int);
        virtual void SetInitialParticleState(MPMBase *,int,int) const;
	
		// history data
        //virtual int SizeOfHistoryData(void) const;
        //virtual char *InitHistoryData(char *,MPMBase *);
        //virtual double GetHistory(int,char *) const;
 	
		// const methods
        virtual void PrintMechanicalProperties(void) const;
		
		// methods
        virtual int SizeOfMechanicalProperties(int &) const;
		virtual void *GetCopyOfMechanicalProps(MPMBase *,int,void *,void *,int) const;
		virtual void MPMConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *,int) const;
		virtual void ConstitutiveLaw(MPMBase *,Matrix3,double,int,void *,ResidualStrains *) const;
        //virtual void PlasticityConstLaw(MPMBase *,double,double,double,double,double,int,
        //                                double,double,DP_Properties *,ResidualStrains *,Matrix3 *) const;

		// custom methods: Find yield function and solve for lambda
		//virtual void UpdatePressure(MPMBase *,double,int,DP_Properties *,ResidualStrains *,double,double &,double &) const;
        virtual double GetMagnitudeSFromDev(Tensor *,int) const;
		virtual void GetDfDsigma(double,Tensor *,int,Tensor *) const;
		
		// accessors
        virtual Tensor GetStress(Tensor *,double,MPMBase *) const;
        const char *MaterialType(void) const;
		virtual int AltStrainContains(void) const;
		
    protected:
		DP_Properties pr;
		double G0;
        double phi;
        double psi;
        double c;
        int par_phi;
        //HardeningLawBase *plasticLaw;

};

#endif

