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
    //plasticLaw = new LinearHardening(this);
}

#pragma mark DP_Plasticity::Initialization

// Read material properties
char *DP_Plasticity::InputMaterialProperty(char *xName,int &input,double &gScaling)
{
    //// look for different plastic law
    //if(strcmp(xName,"Hardening")==0)
    //{	input = HARDENING_LAW_SELECTION;
    //    return (char *)this;
    //}
    //// check plastic law
    //char *ptr = plasticLaw->InputMaterialProperty(xName,input,gScaling);
    //if(ptr != NULL) return ptr;
    
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
        return "par_phi must be 1-3";
    }
    pr.cphired = pr.cphi/this->GetRho(NULL);
    //End verify plastic parameter

    //Start Verify elastic parameters
	// check in superclass (along with its initialization)
    const char *ptr = IsotropicMat::VerifyAndLoadProperties(np);
    
	// reduced prooperties
	G0red = C66/rho;
	pr.Gred = G0red;
    // from C33 = lambda + 2G = K + 4G/3
	pr.Kred = C33/rho - 4.*G0red/3.;							
	return ptr;
}

// Allows any hardening law
//bool DP_Plasticity::AcceptHardeningLaw(HardeningLawBase *pLaw,int lawID)
//{   delete plasticLaw;
//    plasticLaw = pLaw;
//    return TRUE;
//}

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

#pragma mark DP_Plasticity::History Data Methods

// The DP_Plasticity has not history data, but its plasticity law might
//char *DP_Plasticity::InitHistoryData(char *pchr,MPMBase *mptr)
//{
////	int num = plasticLaw->HistoryDoublesNeeded();
//    int num = 0;
//	if(num==0) return NULL;
//	//double *p = CreateAndZeroDoubles(pchr,num);
//	//plasticLaw->InitPlasticHistoryData(p);
//    //return (char *)p;
//}
//
//// DP_Plasticity has no history data, by the hardening law might
//double DP_Plasticity::GetHistory(int num,char *historyPtr) const
//{	//return plasticLaw->GetHistory(num,historyPtr);
//    return (double) 0.;
//
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
    LRConstitutiveLaw(mptr, du, delTime, np, properties, res);
}

#pragma mark DP_Plasticity::Methods (Large Rotation)

// Entry point for large rotation
void DP_Plasticity::LRConstitutiveLaw(MPMBase *mptr,Matrix3 du,double delTime,int np,void *properties,ResidualStrains *res) const
{
	// current previous deformation gradient and stretch
	Matrix3 pFnm1 = mptr->GetDeformationGradientMatrix();
	
    // get incremental deformation gradient and decompose it
	const Matrix3 dF = du.Exponential(incrementalDefGradTerms);
    Matrix3 dR;
    Matrix3 dV = dF.LeftDecompose(&dR,NULL);
	
	// decompose to get previous stretch
	Matrix3 Vnm1 = pFnm1.LeftDecompose(NULL,NULL);
	
	// get strain increments de = (dV-I) dR Vnm1 dRT
	dV(0,0) -= 1.;
	dV(1,1) -= 1.;
	dV(2,2) -= 1.;
	Matrix3 de = dV*Vnm1.RMRT(dR);
	
	// Update total deformation gradient
	Matrix3 pF = dF*pFnm1;
    // Update Total Strain
	mptr->SetDeformationGradientMatrix(pF);
	
    // Effective strain by deducting thermal strain (no shear thermal strain because isotropic)
	//  (note: using unreduced terms in CTE3 and CME3)
	double eres=CTE3*res->dT;
	if(DiffusionTask::active)
		eres+=CME3*res->dC;
	
	// Trial update assuming elastic response
	double delV;
    
    // 3D or 2D
	DP_Properties *p = (DP_Properties *)properties;
    delV = de.trace() - 3.*eres;
    LRPlasticityConstLaw(mptr,de(0,0),de(1,1),2.*de(0,1),de(2,2),delTime,np,delV,eres,p,res,&dR);
}
//For 2D
// To allow some subclasses to support large deformations, the initial calculation for incremental
// deformation gradient (the dvij), volume change (delV)
// handled first by the subclass. This method then finishes the constitutive law
void DP_Plasticity::LRPlasticityConstLaw(MPMBase *mptr,double dexx,double deyy,double dgxy,
									   double dezz,double delTime,int np,double delV,double eres,
									   DP_Properties *p,ResidualStrains *res,Matrix3 *dR) const
{
	// here deij is total strain increment, dexxr is relative strain by subtracting off eres
    double dexxr = dexx-eres;
    double deyyr = deyy-eres;
	double dezzr = dezz-eres;			
	
	// allow arbitrary equation of state for pressure
    //double P0 = mptr->GetPressure();
	//double dTq0 = 0.,dispEnergy = 0.;
	//UpdatePressure(mptr,delV,np,p,res,eres,dTq0,dispEnergy);
    //double Pfinal = mptr->GetPressure();
    //double Pfinal = P0 - p->Kred*delV;

    // MPM P142 (5.54) ; P = - sigma_m
    double P0 = mptr->GetPressure();
    double dPtrial = -p->Kred*delV;    
    double Ptrial = P0 + dPtrial;

    // Deviatoric stress increment
	Tensor *sp=mptr->GetStressTensor();
    Tensor dels,stk;
	//Tensor st0=*sp;
    double thirdDelV = delV/3.;
	dels.xx = 2.*p->Gred*(dexxr-thirdDelV);
	dels.yy = 2.*p->Gred*(deyyr-thirdDelV);
    dels.zz = 2.*p->Gred*(dezzr - thirdDelV);
    // dgxy = 2.*de(0,1)
	dels.xy = p->Gred*dgxy;
	
	// incremental rotate of plastic strain
	Tensor *eplast=mptr->GetAltStrainTensor();
	Matrix3 etn(eplast->xx,0.5*eplast->xy,0.5*eplast->xy,eplast->yy,eplast->zz);
	Matrix3 etr = etn.RMRT(*dR);
	eplast->xx = etr(0,0);
	eplast->yy = etr(1,1);
	eplast->xy = 2.*etr(0,1);
	
	// incremental rotation of stress
    // MPM P136 (5.8)
	Matrix3 stn(sp->xx,sp->xy,sp->xy,sp->yy,sp->zz);
	Matrix3 str = stn.RMRT(*dR);
	
    // trial deviatoric stress
    // MPM P138 (5.29)
    stk.xx = str(0,0) + dels.xx;
    stk.yy = str(1,1) + dels.yy;
	stk.zz = str(2,2) + dels.zz;
    stk.xy = str(0,1) + dels.xy;
	
    // Calculate plastic potential f = ||s|| - sqrt(2/3)*sy(alpha,rate,...)
	//HardeningAlpha alpha;
	//plasticLaw->UpdateTrialAlpha(mptr,np,&alpha,0);			// initialize to last value and zero plastic strain rate
    //||s|| = sqrt(s.s) = sqrt(2J2) where J2 = (1/2)s.s
	double strial = GetMagnitudeSFromDev(&stk,np);
    double sqrtJ2trial = strial / SQRT_TWO;

    //Origin Yield Function
    //double ftrial = strial - SQRT_TWOTHIRDS*plasticLaw->GetYield(mptr,np,delTime,&alpha,p->hardProps);   
    //Linhl  Add DP Yield Criterion
    //double qphi_dp = 2.*SQRT_THREE/5.;  //phi=30deg, circumcircle, triaxial compression 
    //double qphi_dp = 0.; //phi=20deg, inscribed,MPM P156(5.113)
    //double qphi_dp = 0.4803844614; //phi=30deg, inscribed,MPM P156(5.113)
    //double kphi_dp = 0.75e6/SQRT_THREE/this->GetRho(NULL);    
    //double kphi_dp = 1.e6 / SQRT_THREE / this->GetRho(NULL);    
    //double ftrial = strial/SQRT_TWO - qphi_dp*Pfinal - kphi_dp;

    // strial = sqrt(2 J2)
    // MPM P155 (5.111)
    double ftrial = sqrtJ2trial - pr.aphi*Ptrial - pr.cphired;
    
    //Division of plastic flow area
    // MPM P158 (5.121)
    //double ab_dp = sqrt(1+pr.aphi*pr.aphi)-pr.aphi;
    //double htrial = strial / SQRT_TWO + ab_dp*(Pfinal+pr.cphired/pr.aphi);

    //tensile strength; stress of tensile
    double sst_dp = pr.cphired / pr.aphi;
    /*if ( -Pfinal > sst_dp  && htrial <=0 ) {
        UpdatePressure(mptr, 0, np, p, res, eres, dTq0, dispEnergy);
        return;   
    }*/

    //  tensile failure
    if (-Ptrial > sst_dp) {
        //UpdatePressure(mptr, 0, np, p, res, eres, dTq0, dispEnergy);
        // MPM P160 (5.144)
        double lambdat = (-Ptrial - sst_dp) / pr.Kred;

        // MPM P160 (5.145)
        eplast->xx += lambdat / 3.;
        eplast->yy += lambdat / 3.;
        eplast->zz += lambdat / 3.;

        mptr->IncrementPressure(-Ptrial-sst_dp );
        sp->xx = stk.xx;
        sp->xy = stk.xy;
        sp->yy = stk.yy;
        sp->zz = stk.zz;
        return;
    }
    
    // elastic, update stress and strain energy as usual
	if(ftrial<0.)
	{	
		//Linhl, update pressure
        // update pressure
        mptr->IncrementPressure(dPtrial);
		// set in plane deviatoric stress
		sp->xx = stk.xx;
		sp->xy = stk.xy;
		sp->yy = stk.yy;
        sp->zz = stk.zz;
        
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
    
    // MPM P156 (5.132) 
    double lambdak = ftrial / ( p->Gred + p->Kred*pr.aphi*pr.apsi );
	
    //Linhl, update pressure

    
    // Now have lambda, finish update on this particle

    // MPM P159 (5.129) dfds = sij/sqrt(sij:sij) = sij/sqrt(2 J2) = sij/sqrt(2)/sqrt(J2)
	//Tensor dfds;
    //GetDfDsigma(strial, &stk, np, &dfds);
	
    // Plastic strain increments on particle    
    //double dexxp = lambdak*dfds.xx;
    //double deyyp = lambdak*dfds.yy;
	//double dezzp = lambdak*dfds.zz;
    //double dgxyp = 2.*lambdak*dfds.xy;     // 2 for engineering plastic shear strain

    // Plastic strain increments on particle    
    // MPM P159 (5.133)
    double dexxp = lambdak*stk.xx / 2. / sqrtJ2trial;
    double deyyp = lambdak*stk.yy / 2. / sqrtJ2trial;
    double dezzp = lambdak*stk.zz / 2. / sqrtJ2trial;
    double dgxyp = 2.*lambdak*stk.xy / 2. / sqrtJ2trial;     // 2 for engineering plastic shear strain
    // MPM P159 (5.134)
    double dekkp = lambdak*pr.apsi;

	/*==================*/
    /* Start Update Data*/
    /*==================*/
    // MPM P160 (5.133) (5.134)
	eplast->xx += dexxp + dekkp / 3.;
    eplast->yy += deyyp + dekkp / 3.;
	eplast->zz += dezzp + dekkp / 3.;
    eplast->xy += dgxyp;
    
    // MPM P142 (5.52) P160 (5.137)
    mptr->IncrementPressure(dPtrial + pr.Kred*dekkp);

	// increment particle deviatoric stresses (plane stress increments found above)
	// MPM P142 (5.50) (5.51)
    dels.xx -= 2.*p->Gred*dexxp;
    dels.yy -= 2.*p->Gred*deyyp;
	dels.zz -= 2.*p->Gred*dezzp;
    dels.xy -= p->Gred*dgxyp;
		
	// update in-plane stressees
    // MPM P142 (5.51) P160 (5.136)
	sp->xx = str(0,0) + dels.xx;
	sp->yy = str(1,1) + dels.yy;
	sp->xy = str(0,1) + dels.xy;
	//Plane Strain Condition
    sp->zz += dels.zz;
    //UpdatePressure(mptr, delV, np, p, res, eres, dTq0, dispEnergy);
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
void DP_Plasticity::UpdatePressure(MPMBase *mptr,double delV,int np,DP_Properties *p,ResidualStrains *res,double eres,
								   double &dTq0,double &dispEnergy) const
{   // pressure change
    double dP = -p->Kred*delV;
    mptr->IncrementPressure(dP);
    
	// get total dV
	//double dVoverV;
	//if(np==PLANE_STRESS_MPM)
	//	dVoverV = delV + 2.*p->psRed*eres;
	//else
	//	dVoverV = delV + 3.*eres;
    //dVoverV = delV + 3.*eres;
    // work energy is dU = -P dVtot + s.de(total)
	// Here do hydrostatic term
    // Work energy increment per unit mass (dU/(rho0 V0)) (nJ/g)
    //double avgP = mptr->GetPressure()-0.5*dP;
    //mptr->AddWorkEnergyAndResidualEnergy(-avgP*delV,-3.*avgP*eres);	
	// Isoentropic temperature rise = -(K 3 alpha T)/(rho Cv) (dV/V) = - gamma0 T (dV/V)
	//dTq0 = -gamma0*mptr->pPreviousTemperature*dVoverV;
}

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
{   //altBufferSize = plasticLaw->SizeOfHardeningProps();
    altBufferSize = 0;
    return sizeof(DP_Properties);
}

// Isotropic material can use read-only initial properties
void *DP_Plasticity::GetCopyOfMechanicalProps(MPMBase *mptr,int np,void *matBuffer,void *altBuffer,int offset) const
{
	DP_Properties *p = (DP_Properties *)matBuffer;
	*p = pr;
    return p;
 	//p->hardProps = plasticLaw->GetCopyOfHardeningProps(mptr,np,altBuffer,offset);
	//double Gratio = plasticLaw->GetShearRatio(mptr,mptr->GetPressure(),1.,p->hardProps,offset);
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



