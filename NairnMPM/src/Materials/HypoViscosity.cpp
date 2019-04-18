/********************************************************************************
    HypoViscosity.cpp
    nairn-mpm-fea
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "stdafx.h"
#include "HypoViscosity.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Materials/LinearHardening.hpp"
#include "System/UnitsController.hpp"

#pragma mark HypoViscosity::Constructors and Destructors

// Constructors
HypoViscosity::HypoViscosity() {}

// Constructors
// throws std::bad_alloc
HypoViscosity::HypoViscosity(char *matName) : IsotropicMat(matName)
{

}

#pragma mark HypoViscosity::Initialization

// Read material properties
char *HypoViscosity::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
   
    if (strcmp(xName, "mu_s") == 0)
    {
        input = DOUBLE_NUM;
        return((char *)&pr.mu_s);
    }

    else if (strcmp(xName, "mu_2") == 0)
    {
        input = DOUBLE_NUM;
        return((char *)&pr.mu_2);
    }

    else if (strcmp(xName, "xi") == 0)
    {
        input = DOUBLE_NUM;
        return((char *)&pr.xi);
    }

    else if (strcmp(xName, "Is_visco") == 0)
    {
        input = INT_NUM;
        return((char *)&pr.Is_visco);
    }

    else if (strcmp(xName, "volMax") == 0)
    {
        input = DOUBLE_NUM;
        return((char *)&pr.volMax);
    }

    else if (strcmp(xName, "pMin") == 0)
    {
        input = DOUBLE_NUM;
        return((char *)&pr.pMin);
    }

    else if (strcmp(xName, "Is_GeState") == 0)
    {
        input = INT_NUM;
        return((char *)&pr.Is_GeState);
    }

    else if (strcmp(xName, "Is_KeState") == 0)
    {
        input = INT_NUM;
        return((char *)&pr.Is_KeState);
    }
    // otherwise get material properties
    return(IsotropicMat::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *HypoViscosity::VerifyAndLoadProperties(int np)
{

    //Start Verify elastic parameters
	// check in superclass (along with its initialization)
    const char *ptr = IsotropicMat::VerifyAndLoadProperties(np);
    
	// reduced prooperties
	G0 = C66;
	pr.G = G0;
    // from C33 = lambda + 2G = K + 4G/3
	pr.K = C33 - 4.*G0/3.;							
	return ptr;
}

// print mechanical properties to the results
void HypoViscosity::PrintMechanicalProperties(void) const
{	
    IsotropicMat::PrintMechanicalProperties();
    PrintProperty("mu_s", pr.mu_s, "");
    PrintProperty("mu_2", pr.mu_2, "");
    PrintProperty("xi", pr.xi, "");
    PrintProperty("Is_visco", pr.Is_visco, "");
    PrintProperty("volMax", pr.volMax, "");
    PrintProperty("pMin", pr.pMin, "");
    PrintProperty("Is_GeState", pr.Is_GeState, "");
    PrintProperty("Is_KeState", pr.Is_KeState, "");
    cout << endl;
}

// Set intial particle Left Cauchy strain tensor to identity
void HypoViscosity::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
{
    double p0 = 0.;
	Tensor *sp=mptr->GetStressTensor();
    p0 = -(sp->xx + sp->yy + sp->zz) / 3;
    mptr->IncrementPressure(p0);
    sp->xx += p0;
    sp->yy += p0;
    sp->zz += p0;

	MaterialBase::SetInitialParticleState(mptr,np,offset);
}

#pragma mark HypoViscosity::History Data Methods

// return number of bytes needed for history data
int HypoViscosity::SizeOfHistoryData(void) const { return 4*sizeof(double); }

char *HypoViscosity::InitHistoryData(char *pchr,MPMBase *mptr)
{
	double *p = CreateAndZeroDoubles(pchr,4);
    Tensor *sp=mptr->GetStressTensor();
    double p0 = -(sp->xx + sp->yy + sp->zz) / 3;
    double J0 = 1.;
    if ( p0 > 0 ) {
        J0 = 1 - p0 * rho / pr.K;
        if (J0 < 0.01)
            J0 = 0.01;
        Matrix3 F(pow(J0,1./3),0.,0.,pow(J0,1./3),pow(J0,1./3));
        mptr->SetDeformationGradientMatrix(F);      
    }
    
    p[0]=J0;				   // J
    p[2]=pr.mu_s;			   // mu_s
	return (char *)p;
}

// archive material data for this material type when requested.
double HypoViscosity::GetHistory(int num,char *historyPtr) const
{   double history=0.;
    if(num>0 && num<=4)
    {	double *J=(double *)historyPtr;
        history = J[num-1];
    }
    return history;
}

#pragma mark HypoViscosity::Methods

/* Take increments in strain and calculate new Particle: strains, rotation strain, plastic strain,
		stresses, strain energy, plastic energy, dissipated energy, angle
	dvij are (gradient rates X time increment) to give deformation gradient change
	For Axisymmetry: x->R, y->Z, z->theta, np==AXISYMMEtRIC_MPM, otherwise dvzz=0
	This is general analysis for isotropic plastic material. Subclass must define
		GetYield() and GetKPrime() and optionally can override more. Those methods
		require history dependent variables and rates (e.g. cum. plastic strain (alpint) and
		plastic strain rate (dalpha/delTime) to be set before they are called.
	This material tracks pressure and stores deviatoric stress only in particle stress
		tensor
*/
void HypoViscosity::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{   
    ConstitutiveLaw(mptr, du, delTime, np, properties, res, historyOffset);
}

#pragma mark HypoViscosity::Methods (Large Rotation)

// Entry point for constitutive law
void HypoViscosity::ConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{
	// get incremental deformation gradient
	const Matrix3 dF = du.Exponential(1);
     Matrix3 dR;
     Matrix3 dV = dF.LeftDecompose(&dR,NULL);
	
	// current deformation gradient
	Matrix3 Fn = mptr->GetDeformationGradientMatrix();
    // without implement GetCurrentRelativeVolume, always J =1
    double Jn = 1.;
    Matrix3 Fk = dF*Fn;
    double detFk = Fk.determinant();
    if ( pr.volMax >= 1. && detFk > pr.volMax ) {
        Fk.Scale(pow(pr.volMax/dF.determinant()/detFk,1./3));
    }
    mptr->SetDeformationGradientMatrix( Fk );
    


    // without implement GetCurrentRelativeVolume, always J =1
    double Jk = 1.;
    
    double Vn = mptr->GetRelativeVolume();
    mptr->SetHistoryDble(0,Vn,historyOffset);

    // Deviatoric D
    double delV = du.trace();
    double thirdDelV = delV/3.;
    Tensor DdevDt;
	DdevDt.xx = du(0,0)-thirdDelV;
	DdevDt.yy = du(1,1)-thirdDelV;
    DdevDt.zz = du(2,2)-thirdDelV;
	DdevDt.xy = 0.5*(du(0,1)+du(1,0));
    double DdevDt_M = sqrt(DoubleInner(&DdevDt,&DdevDt,np))/ SQRT_TWO ;
    // ------- velocity (Legacy units mm/sec)
    mptr->SetHistoryDble(1,0.001*sqrt(mptr->vel.x*mptr->vel.x+mptr->vel.y*mptr->vel.y),historyOffset);
    mptr->SetTemperature(DdevDt_M/delTime, 0);

    // End of Update Deformation
    // Start to Check Stress
    HV_Properties *p = (HV_Properties *)properties;
    Tensor *sp=mptr->GetStressTensor();

    // Update Pressure
    double P0 = rho/Jn*mptr->GetPressure();
    double dP = -p->K* du.trace();
    double Pk = P0 + dP;
    if ( Pk < pr.pMin  ) {
      mptr->SetPressure( Jk/rho*pr.pMin );
      sp->xx = 0.;
      sp->xy = 0.;
      sp->yy = 0.;
      sp->zz = 0.;
      // pressure = pMin means it's rigid
      mptr->SetHistoryDble(2,pr.mu_s,historyOffset);
      mptr->SetHistoryDble(3,0.,historyOffset);
      return;
    }

    mptr->SetPressure( Jk/rho*Pk );

    // If is GeSate, use pressure to update Ge 
    if (pr.Is_GeState == 1) {
        p->G = pr.mu_s / 0.2 * Pk;
    }

    // Deviatoric stress
    Matrix3 Sdev_n(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
    Sdev_n.Scale( rho/Jn );
    double Sdev_n_M = rho/Jn*sqrt(DoubleInner(sp,sp,np))/ SQRT_TWO;

    // rotation of stress
	Matrix3 Sdev_R = Sdev_n.RMRT(dR);

    // End of All Preparation
    // Start to Update Stress
    
    // Calculate mu_k
    double mu_k = pr.mu_s;
    if (DdevDt_M > 1e-5 * pr.xi * sqrt(Pk) ) {
        mu_k += ( pr.mu_2 - pr.mu_s ) * DdevDt_M / (  pr.xi * sqrt(Pk)*delTime + DdevDt_M );
    }
    //if (DdevDt_M > 1e-5 * pr.xi ) {
    //    mu_k += ( pr.mu_2 - pr.mu_s ) * DdevDt_M / (  pr.xi * delTime + DdevDt_M );
    //}
    mptr->SetHistoryDble(2,mu_k,historyOffset);

    // Calculate Deviatoric stress increment
    Tensor dSdev_E;
	dSdev_E.xx = 2.*p->G*DdevDt.xx;
	dSdev_E.yy = 2.*p->G*DdevDt.yy;
    dSdev_E.zz = 2.*p->G*DdevDt.zz;
	dSdev_E.xy = 2.*p->G*DdevDt.xy;

	// incremental rotate of plastic strain
	Tensor *eplast=mptr->GetAltStrainTensor();
	Matrix3 etn(eplast->xx,eplast->xy,eplast->xy,eplast->yy,eplast->zz);
	Matrix3 etr = etn.RMRT(dR);
	eplast->xx = etr(0,0);
	eplast->yy = etr(1,1);
	eplast->xy = etr(0,1);

    // Calculate alpha
    double alpha = 0.;
    double dSdev_E_M = 2 * pr.G * DdevDt_M;
    if (  pr.Is_visco == 1 
          //&& 1.e-10 * Sdev_n_M < dSdev_E_M  &&  dSdev_E_M < 1.e10 * Sdev_n_M
          && Sdev_n_M > 0. && DdevDt_M > 0.
        ) 
    {
        alpha = 0.5*rho/Jn*DoubleInner(sp,&DdevDt,np)/Sdev_n_M/DdevDt_M;
        alpha = (alpha > 0.) ? alpha : 0.;
        alpha = (alpha < 1.) ? alpha : 1.;
    }


    // Calculate trial deviatoric stress
    Tensor Sdev_tr;
        Sdev_tr.xx = alpha * Sdev_n(0,0) + (1. - alpha ) * Sdev_R(0,0) + dSdev_E.xx ;
        Sdev_tr.yy = alpha * Sdev_n(1,1) + (1. - alpha ) * Sdev_R(1,1) + dSdev_E.yy;
	    Sdev_tr.zz = alpha * Sdev_n(2,2) + (1. - alpha ) * Sdev_R(2,2) + dSdev_E.zz;
        Sdev_tr.xy = alpha * Sdev_n(0,1) + (1. - alpha ) * Sdev_R(0,1) + dSdev_E.xy;
	
    // Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
    // ||s|| = sqrt(s.s) = sqrt(2J2) where J2 = (1/2)s.s
    // strial = sqrt(2 J2)
    double Sdev_tr_M = sqrt(DoubleInner(&Sdev_tr,&Sdev_tr,np))/ SQRT_TWO;
    // MPM P155 (5.111)
    double ftrial = Sdev_tr_M - mu_k*Pk;
    
    
    if (ftrial <= 0.) {
        // Elasticity
        mptr->SetHistoryDble(3,0.,historyOffset);

        // Calculate S(n+1) Elasticity
        sp->xx = Jk/rho*( Sdev_R(0,0) + dSdev_E.xx );
        sp->xy = Jk/rho*( Sdev_R(0,1) + dSdev_E.xy );
        sp->yy = Jk/rho*( Sdev_R(1,1) + dSdev_E.yy );
        sp->zz = Jk/rho*( Sdev_R(2,2) + dSdev_E.zz );

    }
    else{
        // Viscosity update plasticity strain
        mptr->SetHistoryDble(3,alpha,historyOffset);

	    // Find  lambda for this plastic state
        // Plastic strain increments on particle    
        double lambdak = (Sdev_tr_M - mu_k*Pk) / ( p->G );
        double dsnpDxx = lambdak * Sdev_tr.xx / 2. / Sdev_tr_M;
        double dsnpDyy = lambdak * Sdev_tr.yy / 2. / Sdev_tr_M;
        double dsnpDzz = lambdak * Sdev_tr.zz / 2. / Sdev_tr_M;
        double dsnpDxy = lambdak * Sdev_tr.xy / 2. / Sdev_tr_M;
    	eplast->xx += dsnpDxx;
        eplast->yy += dsnpDyy;
    	eplast->zz += dsnpDzz;
        eplast->xy += dsnpDxy;

        // Calculate S(n+1) Viscosity
        Sdev_tr.xx *= mu_k*Pk / Sdev_tr_M;
        Sdev_tr.yy *= mu_k*Pk / Sdev_tr_M;
	    Sdev_tr.zz *= mu_k*Pk / Sdev_tr_M;
        Sdev_tr.xy *= mu_k*Pk / Sdev_tr_M;
        sp->xx = Jk/rho*Sdev_tr.xx;
        sp->xy = Jk/rho*Sdev_tr.xy;
        sp->yy = Jk/rho*Sdev_tr.yy;
        sp->zz = Jk/rho*Sdev_tr.zz; 
    }
   
    return;
    
}


#pragma mark HypoViscosity::Custom Methods

double HypoViscosity::DoubleInner(Tensor *st1,Tensor *st2,int np) const
{
	double s,t;
	
	switch(np)
    {   case THREED_MPM:
            s = st1->xx*st2->xx + st1->yy*st2->yy + st1->zz*st2->zz;
            t = st1->xy*st2->xy + st1->xz*st2->xz + st1->yz*st2->yz;
            break;
            
		default:
			s = st1->xx*st2->xx + st1->yy*st2->yy + st1->zz*st2->zz;
			t = st1->xy*st2->xy;
			break;
	}
	return (s+t+t);
}

// return derivatives of the yield function wrt to components of deviatoric stress
// which for isotropic material with f = ||s|| - sqrt(2/3)*sy
// reduces to s/||s|| (written for tensorial plastic strains)
// But do not call for plane stress, which must be found in special case
void HypoViscosity::GetDfDsigma(double smag,Tensor *st0,int np,Tensor *dfds) const
{
    // s/||s|| = n
    dfds->xx = st0->xx/smag;
    dfds->yy = st0->yy/smag;
    dfds->zz = st0->zz/smag;
    dfds->xy = st0->xy/smag;		// tensorial shear strain
    if(np==THREED_MPM)
    {	dfds->xz = st0->xz/smag;
        dfds->yz = st0->yz/smag;
    }
}

#pragma mark HypoViscosity::Accessors

// store plastic strain in alt strain
int HypoViscosity::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}

// buffer size for mechanical properties
int HypoViscosity::SizeOfMechanicalProperties(int &altBufferSize) const
{
    altBufferSize = 0;
    return sizeof(HV_Properties);
}

// Isotropic material can use read-only initial properties
void *HypoViscosity::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	HV_Properties *p = (HV_Properties *)matBuffer;
	*p = pr;
    return p;

}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor HypoViscosity::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// return unique, short name for this material
const char *HypoViscosity::MaterialType(void) const { return "HypoViscosity"; }



