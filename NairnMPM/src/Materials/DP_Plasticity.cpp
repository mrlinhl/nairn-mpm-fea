/********************************************************************************
    DP_Plasticity.cpp
    nairn-mpm-fea
    
    Created by John Nairn, June 16, 2008.
    Copyright (c) 2008 John A. Nairn, All rights reserved.

	Dependencies
		Orthotropic.hpp (TranIsotropic.hpp, MaterialBase.hpp)
********************************************************************************/

#include "stdafx.h"
#include "DP_Plasticity.hpp"
#include "MPM_Classes/MPMBase.hpp"
#include "Custom_Tasks/ConductionTask.hpp"
#include "Custom_Tasks/DiffusionTask.hpp"
#include "Exceptions/CommonException.hpp"
#include "Materials/HardeningLawBase.hpp"
#include "Materials/LinearHardening.hpp"
#include "System/UnitsController.hpp"

#pragma mark DP_Plasticity::Constructors and Destructors

// Constructors
DP_Plasticity::DP_Plasticity() {}

// Constructors
// throws std::bad_alloc
DP_Plasticity::DP_Plasticity(char *matName) : IsotropicMat(matName)
{

}

#pragma mark DP_Plasticity::Initialization

// Read material properties
char *DP_Plasticity::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
   
    if (strcmp(xName, "phi") == 0)
    {
        input = DOUBLE_NUM;
        return((char *)&phi);
    }

    else if (strcmp(xName, "psi") == 0)
    {
        input = DOUBLE_NUM;
        return((char *)&psi);
    }

    else if (strcmp(xName, "c") == 0)
    {
        input = DOUBLE_NUM;
        return UnitsController::ScaledPtr((char *)&c, gScaling, 1.e6);
    }

    else if (strcmp(xName, "par_phi") == 0)
    {
        input = INT_NUM;
        return((char *)&par_phi);
    }


    // otherwise get material properties
    return(IsotropicMat::InputMaterialProperty(xName,input,gScaling));
}

// verify settings and some initial calculations
const char *DP_Plasticity::VerifyAndLoadProperties(int np)
{
    //Start Verify plastic parameters
    double phi_arc = phi*PI_CONSTANT / 180.;
    double psi_arc = psi*PI_CONSTANT / 180.;
    if (c < 0.) return "The cohesion is negative";
    
    switch (par_phi)
    {
    case 1:
        // MPM P156 (5.112)
        pr.aphi = 6.*sin(phi_arc) / (sqrt(3.)*(3-sin(phi_arc))) ;
        pr.apsi = 6.*sin(psi_arc) / (sqrt(3.)*(3 - sin(psi_arc)));
        pr.cphi = 6.*c*cos(phi_arc) / (sqrt(3.)*(3 - sin(phi_arc)));
        break;

    case 2:
        pr.aphi = 6.*sin(phi_arc) / (sqrt(3.)*(3 + sin(phi_arc)));
        pr.apsi = 6.*sin(psi_arc) / (sqrt(3.)*(3 + sin(psi_arc)));
        pr.cphi = 6.*c*cos(phi_arc) / (sqrt(3.)*(3 + sin(phi_arc)));
        break;
    case 3:
        pr.aphi = 3.*tan(phi_arc) / sqrt(9.+12.*tan(phi_arc)*tan(phi_arc));
        pr.apsi = 3.*tan(psi_arc) / sqrt(9. + 12.*tan(psi_arc)*tan(psi_arc));
        pr.cphi = 3.*c / sqrt(9. + 12.*tan(phi_arc)*tan(phi_arc));
        break;
    default:
        return "par_phi must be 1 or 2 or 3";
    }

    //End verify plastic parameter

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
void DP_Plasticity::PrintMechanicalProperties(void) const
{	
    IsotropicMat::PrintMechanicalProperties();
    PrintProperty("phi", phi, "");
    PrintProperty("psi", psi, "");
    PrintProperty("c", c*UnitsController::Scaling(1.e-6), "");
    PrintProperty("par_phi", par_phi, "");
    PrintProperty("aphi", pr.aphi, "");
    PrintProperty("apsi", pr.apsi, "");
    PrintProperty("cphi", pr.cphi*UnitsController::Scaling(1.e-6), "");
    cout << endl;
}

// Set intial particle Left Cauchy strain tensor to identity
void DP_Plasticity::SetInitialParticleState(MPMBase *mptr,int np,int offset) const
{
    double p0 = 0.;
	Tensor *sp=mptr->GetStressTensor();
    p0 = -(sp->xx + sp->yy + sp->zz) / 3;
    mptr->IncrementPressure(p0);
    sp->xx += p0;
    sp->yy += p0;
    sp->zz += p0;

    // current deformation gradient
	Matrix3 F = mptr->GetDeformationGradientMatrix();
	
    mptr->SetDeformationGradientMatrix(F);



	MaterialBase::SetInitialParticleState(mptr,np,offset);
}

#pragma mark DP_Plasticity::History Data Methods
//
//// return number of bytes needed for history data
//int DP_Plasticity::SizeOfHistoryData(void) const { return 1*sizeof(double); }
//
//// Store J, which is calculated incrementally, and available for archiving
//char *DP_Plasticity::InitHistoryData(char *pchr,MPMBase *mptr)
//{
//	double *p = CreateAndZeroDoubles(pchr,1);
//	p[0] = 1.;
//	return (char *)p;
//}
//
//// archive material data for this material type when requested.
//double DP_Plasticity::GetHistory(int num,char *historyPtr) const
//{   double history=0.;
//    if(num>0 && num<=2)
//    {	double *J=(double *)historyPtr;
//        history = J[num-1];
//    }
//    return history;
//}

#pragma mark DP_Plasticity::Methods

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
void DP_Plasticity::MPMConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res,int historyOffset) const
{   
    ConstitutiveLaw(mptr, du, delTime, np, properties, res);
}

#pragma mark DP_Plasticity::Methods (Large Rotation)

// Entry point for constitutive law
void DP_Plasticity::ConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// get incremental deformation gradient
	const Matrix3 dF = du.Exponential(1);
     Matrix3 dR;
     Matrix3 dV = dF.LeftDecompose(&dR,NULL);
	
	// current deformation gradient
	Matrix3 pF = mptr->GetDeformationGradientMatrix();
	const Matrix3 F = dF*pF;	

    // without implement GetCurrentRelativeVolume
    // so Tracked stress is (Cauchy stress)/rho = (Cauchy stress)/rho0
    double Jn = 1.;
    mptr->SetDeformationGradientMatrix(F);
    double Jk = 1.;

    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	//  (note: using unreduced terms in CTE3 and CME3)
	double eres=CTE3*res->dT;
	if(DiffusionTask::active)
		eres+=CME3*res->dC;

	DP_Properties *p = (DP_Properties *)properties;
    // Trial update assuming elastic response
    double delV = du.trace() - 3.*eres;

//    // 2D Only
//    PlasticityConstLaw(mptr,du(0,0),du(1,1),0.5*(du(0,1)+du(1,0)),du(2,2),delTime,np,delV,eres,p,res,&dR);
//}
////For 2D
//// deformation gradient (the dvij), volume change (delV)
//// handled first by the subclass. This method then finishes the constitutive law
//void DP_Plasticity::PlasticityConstLaw(MPMBase *mptr,double dexx,double deyy,double dexy,
//									   double dezz,double delTime,int np,double delV,double eres,
//									   DP_Properties *p,ResidualStrains *res,Matrix3 *dR) const
//{
    // Update current relative volume change by incremental method

    double dexx = du(0,0);
    double deyy = du(1,1);
    double dexy = 0.5*(du(0,1)+du(1,0));
    double dezz = du(2,2);

    // here deij is total strain increment, dsnxx is relative strain by subtracting off eres
    double dsnxx = dexx-eres;
    double dsnyy = deyy-eres;
	double dsnzz = dezz-eres;			

    // MPM P142 (5.54) ; P = - sigma_m  = -ssm
    double P0 = rho/Jn*mptr->GetPressure();
    double dPtrial = -p->K*delV;
    double Ptrial = P0 + dPtrial;

    // Deviatoric stress
	Tensor *sp=mptr->GetStressTensor();
    Matrix3 ssD_n(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
    ssD_n.Scale( rho/Jn );

    // Deviatoric stress increment
    Tensor dssD_J,ssD_k;
    double thirdDelV = delV/3.;
	dssD_J.xx = 2.*p->G*(dsnxx-thirdDelV);
	dssD_J.yy = 2.*p->G*(dsnyy-thirdDelV);
    dssD_J.zz = 2.*p->G*(dsnzz-thirdDelV);
	dssD_J.xy = 2.*p->G*dexy;
	
	// incremental rotate of plastic strain
	Tensor *eplast=mptr->GetAltStrainTensor();
	Matrix3 etn(eplast->xx,eplast->xy,eplast->xy,eplast->yy,eplast->zz);
	Matrix3 etr = etn.RMRT(dR);
	eplast->xx = etr(0,0);
	eplast->yy = etr(1,1);
	eplast->xy = etr(0,1);
	
	// incremental rotation of stress
    // MPM P136 (5.8)
	Matrix3 ssD_R = ssD_n.RMRT(dR);
	
    // trial deviatoric stress
    // MPM P138 (5.29)
    ssD_k.xx = ssD_R(0,0) + dssD_J.xx;
    ssD_k.yy = ssD_R(1,1) + dssD_J.yy;
	ssD_k.zz = ssD_R(2,2) + dssD_J.zz;
    ssD_k.xy = ssD_R(0,1) + dssD_J.xy;
	
    // Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
    // ||s|| = sqrt(s.s) = sqrt(2J2) where J2 = (1/2)s.s
    // strial = sqrt(2 J2)
    double sijsij = GetMagnitudeSFromDev(&ssD_k,np);
    double sqrtJ2trial = sijsij / SQRT_TWO;

    // MPM P155 (5.111)
    double ftrial = sqrtJ2trial - pr.aphi*Ptrial - pr.cphi;
    
    // tensile strength; stress of tensile
    double sst_dp = pr.cphi / pr.aphi;

    //  tensile failure
    if (-Ptrial > sst_dp) {
        // MPM P160 (5.144)
        double lambdat = (-Ptrial - sst_dp) / pr.K;

        // MPM P160 (5.145)
        eplast->xx += lambdat / 3.;
        eplast->yy += lambdat / 3.;
        eplast->zz += lambdat / 3.;

        mptr->IncrementPressure(Jk/rho*(-Ptrial-sst_dp) );
        sp->xx = Jk/rho*ssD_k.xx;
        sp->xy = Jk/rho*ssD_k.xy;
        sp->yy = Jk/rho*ssD_k.yy;
        sp->zz = Jk/rho*ssD_k.zz;
        return;
    }
    
    // elastic, update stress and strain energy as usual
	if(ftrial<=0.)
	{	
        // update pressure
        mptr->IncrementPressure( Jk/rho*dPtrial );
		// set in plane deviatoric stress
		sp->xx = Jk/rho*ssD_k.xx;
		sp->xy = Jk/rho*ssD_k.xy;
		sp->yy = Jk/rho*ssD_k.yy;
        sp->zz = Jk/rho*ssD_k.zz;
        
        //UpdatePressure(mptr, delV, np, p, res, eres, dTq0, dispEnergy);
		// work energy increment per unit mass (dU/(rho0 V0)) (by midpoint rule) (nJ/g)
        // energy units are also Pa mm^3/g, i.e., same as stress units
        //mptr->AddWorkEnergy(sp->xx*dexx + sp->yy*deyy + sp->xy*dgxy);		
		// give subclass material chance to update history variables that change in elastic updates
		//plasticLaw->ElasticUpdateFinished(mptr,np,delTime,0);
		// heat energy from pressure and any pressure method items
		//IncrementHeatEnergy(mptr,res->dT,dTq0,dispEnergy);
	
		return;
    }
    
	// Find  lambda for this plastic state
	// Base class finds it numerically, subclass can override if solvable by more efficient methods
    
    // MPM P159 (5.132) 
    double lambdak = ftrial / ( p->G + p->K*pr.aphi*pr.apsi );
    
    // Now have lambda, finish update on this particle

    // Plastic strain increments on particle    
    // MPM P159 (5.133)
    double dsnpDxx = lambdak * ssD_k.xx / 2. / sqrtJ2trial;
    double dsnpDyy = lambdak * ssD_k.yy / 2. / sqrtJ2trial;
    double dsnpDzz = lambdak * ssD_k.zz / 2. / sqrtJ2trial;
    double dsnpDxy = lambdak * ssD_k.xy / 2. / sqrtJ2trial;

    // MPM P159 (5.134)
    double dsnpkk = lambdak*pr.apsi;

	/*==================*/
    /* Start Update Data*/
    /*==================*/
    // MPM P160 (5.133) (5.134)
	eplast->xx += dsnpDxx + dsnpkk / 3.;
    eplast->yy += dsnpDyy + dsnpkk / 3.;
	eplast->zz += dsnpDzz + dsnpkk / 3.;
    eplast->xy += dsnpDxy;
    
    // MPM P142 (5.52) P160 (5.137)
    mptr->IncrementPressure( Jk/rho*(dPtrial + pr.K*dsnpkk));

	// increment particle deviatoric stresses (plane stress increments found above)
	// MPM P142 (5.50) (5.51)
    dssD_J.xx -= 2.*p->G*dsnpDxx;
    dssD_J.yy -= 2.*p->G*dsnpDyy;
	dssD_J.zz -= 2.*p->G*dsnpDzz;
    dssD_J.xy -= 2.*p->G*dsnpDxy;
		
	// update in-plane stressees
    // MPM P142 (5.51) P160 (5.136)
	sp->xx = Jk/rho*(ssD_R(0,0) + dssD_J.xx);
	sp->yy = Jk/rho*(ssD_R(1,1) + dssD_J.yy);
	sp->xy = Jk/rho*(ssD_R(0,1) + dssD_J.xy);
	//Plane Strain Condition
    sp->zz = Jk/rho*(ssD_R(2,2) + dssD_J.zz);
    
    // Elastic work increment per unit mass (dU/(rho0 V0)) (nJ/g)
	//double workEnergy = sp->xx*dexx + sp->yy*deyy + sp->xy*dgxy;
	//double workEnergy = 0.5*((st0.xx+sp->xx)*dexx + (st0.yy+sp->yy)*deyy + (st0.xy+sp->xy)*dgxy);  
    // total work
    //mptr->AddWorkEnergy(workEnergy);
    // plastic strain work
    //double plastEnergy = lambdak*(sp->xx*dfds.xx + sp->yy*dfds.yy + sp->zz*dfds.zz + 2.*sp->xy*dfds.xy);
    // dand subtract q dalpha to get isispated energy per unit mass (dPhi/(rho0 V0)) (nJ/g)
    //double qdalphaTerm = lambdak*SQRT_TWOTHIRDS*plasticLaw->GetYieldIncrement(mptr,np,delTime,&alpha,p->hardProps);
    //double qdalphaTerm = 0.;
    //dispEnergy += plastEnergy - qdalphaTerm;
	// The cumulative dissipated energy is tracked in plastic energy
    //mptr->AddPlastEnergy(dispEnergy); 
    // heat energy is Cv(dT-dTq0) - dPhi
    // The dPhi is subtracted here because it will show up in next
    //		time step within Cv dT (if adibatic heating occurs)
    //IncrementHeatEnergy(mptr,res->dT,dTq0,dispEnergy);
    
}


#pragma mark DP_Plasticity::Methods (Small Rotation)



#pragma mark DP_Plasticity::Custom Methods

// This method handles the pressure equation of state. Its tasks are
// 1. Calculate the new pressure
// 2. Update particle pressure
// 3. Increment the particle energy
// 4. Call plasticLaw to see if it wants to change the shear modulus
// 5. Optionally change delV (which is passed by reference)
// Notes:
//  delV is incremental volume change on this step.
//void DP_Plasticity::UpdatePressure(MPMBase *mptr,double delV,int np,DP_Properties *p,ResidualStrains *res,double eres,
//								   double &dTq0,double &dispEnergy) const
//{   // pressure change
//    double dP = -p->K*delV;
//    mptr->IncrementPressure(dP);
//    
//	// get total dV
//	//double dVoverV;
//	//if(np==PLANE_STRESS_MPM)
//	//	dVoverV = delV + 2.*p->psRed*eres;
//	//else
//	//	dVoverV = delV + 3.*eres;
//    //dVoverV = delV + 3.*eres;
//    // work energy is dU = -P dVtot + s.de(total)
//	// Here do hydrostatic term
//    // Work energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
//    //double avgP = mptr->GetPressure()-0.5*dP;
//    //mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);	
//	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
//	//dTq0 = -gamma0*mptr->pPreviousTemperature*dVoverV;
//}

// Get magnitude of the deviatoric stress tensor when input is a deviatoric stress
// ||s|| = sqrt(s.s) = sqrt(2J2) where J2 = (1/2)s.s
// In 2D 2J2 = sx^2 + sy^2 + sz^2 + 2*txy^2
// In 3D 2J2 = sx^2 + sy^2 + sz^2 + 2*txy^2 + 2*txz^2 + 2*tyz^2
double DP_Plasticity::GetMagnitudeSFromDev(Tensor *st,int np) const
{
	double s,t;
	
	switch(np)
    {   case THREED_MPM:
            s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
            t = st->xy*st->xy + st->xz*st->xz + st->yz*st->yz;
            break;
            
		default:
			s = st->xx*st->xx + st->yy*st->yy + st->zz*st->zz;
			t = st->xy*st->xy;
			break;
	}
	return sqrt(s+t+t);
}

// return derivatives of the yield function wrt to components of deviatoric stress
// which for isotropic material with f = ||s|| - sqrt(2/3)*sy
// reduces to s/||s|| (written for tensorial plastic strains)
// But do not call for plane stress, which must be found in special case
void DP_Plasticity::GetDfDsigma(double smag,Tensor *st0,int np,Tensor *dfds) const
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

#pragma mark DP_Plasticity::Accessors

// store plastic strain in alt strain
int DP_Plasticity::AltStrainContains(void) const
{	return ENG_BIOT_PLASTIC_STRAIN;
}

// buffer size for mechanical properties
int DP_Plasticity::SizeOfMechanicalProperties(int &altBufferSize) const
{
    altBufferSize = 0;
    return sizeof(DP_Properties);
}

// Isotropic material can use read-only initial properties
void *DP_Plasticity::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	DP_Properties *p = (DP_Properties *)matBuffer;
	*p = pr;
    return p;

	//p->Gred = G0red*Gratio;    
	//if(np==PLANE_STRESS_MPM)
	//{	// these are terms for plane stress calculations only
	//	p->psRed = 1./(p->Kred/(2.*p->Gred) + 2./3.);					// (1-2nu)/(1-nu) for plane stress
	//	p->psLr2G = (p->Kred/(2.*p->Gred) - 1./3.)*p->psRed;			// nu/(1-nu) to find ezz
	//	p->psKred = p->Kred*p->psRed;									// E/(3(1-v)) to find lambda
	//}
}

// Copy stress to a read-only tensor variable
// Subclass material can override, such as to combine pressure and deviatory stress into full stress
Tensor DP_Plasticity::GetStress(Tensor *sp,double pressure,MPMBase *mptr) const
{   Tensor stress = *sp;
    stress.xx -= pressure;
    stress.yy -= pressure;
    stress.zz -= pressure;
    return stress;
}

// return unique, short name for this material
const char *DP_Plasticity::MaterialType(void) const { return "DP_Plastic"; }



