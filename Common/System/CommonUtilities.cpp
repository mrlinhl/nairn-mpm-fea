/********************************************************************************
    CommonUtilities.hpp
    nairn-mpm-fea
    
    Created by John Nairn on 1/12/06.
    Copyright (c) 2006 John A. Nairn, All rights reserved.
********************************************************************************/

#include "stdafx.h"

#define T_SQRT3    1.73205080756887729352744634151   // sqrt(3)

// local globals
static int section=1;

#pragma mark Miscellaneous Functions

/*********************************************************************
    Start next section in results window
*********************************************************************/

void PrintSection(const char *text)
{
    char sline[200];
    
    sprintf(sline,"*****%3d. %s\n\n",section,text);
    cout << sline;
	section++;
}

/************************************************************
	See if two double precision numbers are sufficiently
		equal (currently 100 ppm)
	If both are less than 1e-12, they are assumed equal
************************************************************/

#define MAX_DIFF 1.0e-16
#define MAX_REL_DIFF 1.0e-7

bool DbleEqual(double A, double B)
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	double diff = fabs(A - B);
	if (diff <= MAX_DIFF)
		return true;

	// Check relative difference
	A = fabs(A);
	B = fabs(B);
	double largest = (B > A) ? B : A;
	if (diff <= largest * MAX_REL_DIFF)
		return true;

	// assuming unequal
	return false;
}

/************************************************************
    Reverse block of bytes to convert between Intel
	processor and Mac processor
************************************************************/

int Reverse(char *nptr,int nlen)
{
    int i,imax;
    char *rptr,temp;
    
    rptr=nptr+nlen-1;
    imax=nlen>>1;
    for(i=0;i<imax;i++)
    {	temp=*nptr;
        *nptr++=*rptr;
        *rptr--=temp;
    }
    return nlen;
}

/************************************************************
	Factor integer num into prime factors (not including 1)
************************************************************/

void PrimeFactors(int num,vector<int> &factors)
{	int maxFactor = (int)(sqrt((double)num)+.001);
	for(int i=2;i<=maxFactor;i++)
	{	// does i go into num evenly?
		if((num % i) == 0)
		{	factors.push_back(i);
			PrimeFactors(num/i, factors);
			return;
		}
	}
	// if did not find any factors, it is a prime
	factors.push_back(num);
}

#pragma mark Vector Functions

// return vector from components
Vector MakeVector(double x,double y,double z)
{	Vector v;
	v.x=x;
	v.y=y;
	v.z=z;
	return v;
}

// zero a vector and return pointer to zeroed v
Vector *ZeroVector(Vector *v)
{	v->x=v->y=v->z=0.;
	return v;
}

// copy vector and return pointer to new one
Vector *CopyVector(Vector *v,const Vector *oldv)
{	v->x=oldv->x;
	v->y=oldv->y;
	v->z=oldv->z;
	return v;
}

// scale a vector and return pointer to changed vector
Vector *ScaleVector(Vector *v,const double scale)
{	v->x*=scale;
	v->y*=scale;
	v->z*=scale;
	return v;
}

// scale a vector and copy to new one, then return pointer to new v
Vector *CopyScaleVector(Vector *v,const Vector *oldv,const double scale)
{	v->x=scale*oldv->x;
	v->y=scale*oldv->y;
	v->z=scale*oldv->z;
	return v;
}

// add vector toadd to vector v and return pointer to changed v
Vector *AddVector(Vector *v,const Vector *toadd)
{	v->x+=toadd->x;
	v->y+=toadd->y;
	v->z+=toadd->z;
	return v;
}

// sub vector toadd from vector v and return pointer to changed v
Vector *SubVector(Vector *v,const Vector *toadd)
{	v->x-=toadd->x;
	v->y-=toadd->y;
	v->z-=toadd->z;
	return v;
}

// add vector toadd scaled by scale to vector v and return pointer to changed v
Vector *AddScaledVector(Vector *v,const Vector *toadd,const double scale)
{	v->x+=toadd->x*scale;
	v->y+=toadd->y*scale;
	v->z+=toadd->z*scale;
	return v;
}

// Dot product of two vectors (if used for 2D make sure z's are zero)
double DotVectors(const Vector *v1,const Vector *v2)
{	return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z ;
}

// Dot prodiuct of two 2D vectors
double DotVectors2D(const Vector *v1,const Vector *v2)
{	return v1->x*v2->x + v1->y*v2->y ;
}

// Vector cross product of two vectors in a third vector
Vector *CrossProduct(Vector *cp,const Vector *a,const Vector *b)
{	cp->x = a->y*b->z - b->y*a->z ;
	cp->y = b->x*a->z - a->x*b->z ;
	cp->z = a->x*b->y - b->x*a->y ;
	return cp;
}

// Z component of vector cross product of two vectors ion x-y plane
double CrossProduct2D(const Vector *a,const Vector *b)
{	return a->x*b->y - b->x*a->y ;
}

// Print vector to cout when debugging
void PrintVector(const char *label,const Vector *v)
{	if(strlen(label)>0) cout << label;
	cout << "(" << v->x << ", " << v->y << ", " << v->z << ") ";
}

#pragma mark Tensor Functions

// return tensor from components
Tensor MakeTensor(double xx,double yy,double zz,double yz,double xz,double xy)
{	Tensor t;
	t.xx = xx;
	t.yy = yy;
	t.zz = zz;
#ifdef MPM_CODE
	t.yz = yz;
	t.xz = xz;
#endif
	t.xy = xy;
	return t;
}

// return tensor from components
Tensor MakeTensor2D(double xx,double yy,double zz,double xy)
{	Tensor t;
	t.xx = xx;
	t.yy = yy;
	t.zz = zz;
#ifdef MPM_CODE
	t.yz = 0.;
	t.xz = 0.;
#endif
	t.xy = xy;
	return t;
}

// zero a vector and return pointer to zeroed v
Tensor *ZeroTensor(Tensor *t)
{	t->xx=0.;
	t->yy=0.;
	t->zz=0.;
#ifdef MPM_CODE
	t->yz=0.;
	t->xz=0.;
#endif
	t->xy=0.;
	return t;
}

// add tensor toadd to tensor t and return pointer to changed t
Tensor *AddTensor(Tensor *t,Tensor *toadd)
{	t->xx+=toadd->xx;
	t->yy+=toadd->yy;
	t->zz+=toadd->zz;
#ifdef MPM_CODE
	t->yz+=toadd->yz;
	t->xz+=toadd->xz;
#endif
	t->xy+=toadd->xy;
	return t;
}

// scale tensor t and return pointer to changed t
Tensor *ScaleTensor(Tensor *t,double scale)
{	t->xx*=scale;
	t->yy*=scale;
	t->zz*=scale;
#ifdef MPM_CODE
	t->yz*=scale;
	t->xz*=scale;
#endif
	t->xy*=scale;
	return t;
}

// Dot product of two vectors (if used for 2D make sure z's are zero)
double DotTensors2D(const Tensor *t1,const Tensor *t2)
{	return t1->xx*t2->xx + t1->yy*t2->yy + t1->zz*t2->zz + t1->xy*t2->xy;
}


#ifdef MPM_CODE
// Dot product of two vectors (if used for 2D make sure z's are zero)
double DotTensors(const Tensor *t1,const Tensor *t2)
{	return t1->xx*t2->xx + t1->yy*t2->yy + t1->zz*t2->zz
	+ t1->yz*t2->yz + t1->xz*t2->xz + t1->xy*t2->xy;
}

// Find Eigenvalues of symmetric tensor in Voight form
// stress is true for stress vector or false for strain vector
// is2D is true for 2D tensor or false for 3D
Vector TensorEigenvalues(Tensor *t,bool stress,bool is2D)
{
	Vector lam;
	
	if(is2D)
	{   // solving x^2 + bx + c = 0 (see Numerical Recipes in C, page 156)
		double b = -(t->xx+t->yy);
		double c = stress ? t->xx*t->yy - t->xy*t->xy : t->xx*t->yy - 0.25*t->xy*t->xy ;
		double arg = b*b-4.*c;
		if(arg<0.)
		{	// assuming here all matrices are positive definite, which means
			// a negative value should be interpreted as zero
			lam.x = -0.5*b;
		}
		else
		{	arg = sqrt(arg);
			lam.x = b>0 ? -0.5*(b+arg) : -0.5*(b-arg) ;
		}
		lam.y = lam.x==0. ? 0. : c/lam.x;
		lam.z = t->zz;
	}
	else
	{	double mm, c1, c0;
		
		// Determine coefficients of characteristic poynomial. We write
		//       | a   d   f  |
		//  m =  | d   b   e  |
		//       | f   e   c  |
		double fde,dd,ee,ff;		// products of elements
		if(stress)
		{	fde = t->xy*t->yz*t->xz;
			dd = t->xy*t->xy;
			ee = t->yz*t->yz;
			ff = t->xz*t->xz;
		}
		else
		{	fde = 0.125*t->xy*t->yz*t->xz;
			dd = 0.25*t->xy*t->xy;
			ee = 0.25*t->yz*t->yz;
			ff = 0.25*t->xz*t->xz;
		}
		mm  = t->xx + t->yy + t->zz;
		c1 = (t->xx*t->yy + t->xx*t->zz + t->yy*t->zz)
					- (dd + ee + ff);				// a*b + a*c + b*c - d^2 - e^2 - f^2
		c0 = t->zz*dd + t->xx*ee + t->yy*ff - t->xx*t->yy*t->zz
					- 2.0*fde;						// c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
		
		double p, sqrt_p, q, c, s, phi;
		p = mm*mm - 3.0*c1;
		q = mm*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
		sqrt_p = sqrt(fabs(p));
		
		phi = 27.0 * ( 0.25*c1*c1*(p - c1) + c0*(q + 27.0/4.0*c0));
		phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
		
		c = sqrt_p*cos(phi);
		s = (1.0/T_SQRT3)*sqrt_p*sin(phi);
		
		lam.y  = (1.0/3.0)*(mm - c);
		lam.z  = lam.y + s;
		lam.x  = lam.y + c;
		lam.y -= s;
	}
	
	return lam;
}

// convert tensor to matrix (stress is true for straight convertion of false to assume strain tensor)
Matrix3 TensorToMatrix(Tensor *t,bool stress)
{	return stress ?
			Matrix3(t->xx,t->xy,t->xz,t->xy,t->yy,t->yz,t->xz,t->yz,t->zz) :
			Matrix3(t->xx,0.5*t->xy,0.5*t->xz,0.5*t->xy,t->yy,0.5*t->yz,0.5*t->xz,0.5*t->yz,t->zz) ;
}

// convert tensor to matrix (stress is true for straight convertion of false to assume strain tensor)
Matrix3 TensorToMatrix2D(Tensor *t,bool stress)
{	return stress ?
			Matrix3(t->xx,t->xy,t->xy,t->yy,t->zz) :
			Matrix3(t->xx,0.5*t->xy,0.5*t->xy,t->yy,t->zz) ;
}

#endif

// get element of contracted tensor by ID
double Tensor_i(Tensor *t,int contractedID)
{
	switch(contractedID)
	{	case XX:
			return t->xx;
		case YY:
			return t->yy;
		case ZZ:
			return t->zz;
#ifdef MPM_CODE
		case YZ:
			return t->yz;
		case XZ:
			return t->xz;
#endif
		case XY:
			return t->xy;
		default:
			return 0.;
	}
}

// get element of tensor by row and colume
// For FEA can only get xx or 1,1, xy or 1,2 or 2,1, yy or 2,2, anbd zz or 3,3
double Tensor_ij(Tensor *t,int row,int col)
{
	switch(row)
	{	case 1:
			switch(col)
			{	case 1:
					return t->xx;
				case 2:
					return t->xy;
#ifdef MPM_CODE
				case 3:
					return t->xz;
#endif
				default:
					return 0.;
			}
		case 2:
			switch(col)
			{	case 1:
					return t->xy;
				case 2:
					return t->yy;
#ifdef MPM_CODE
				case 3:
					return t->yz;
#endif
				default:
					return 0.;
			}
		case 3:
			switch(col)
			{
#ifdef MPM_CODE
				case 1:
					return t->xz;
				case 2:
					return t->yz;
#endif
				case 3:
					return t->zz;
				default:
					return 0.;
			}
		default:
			return 0.;
	}
}

// Print tensor to cout when debuggingin form of a matrix
void PrintTensor(const char *label,Tensor *t)
{
	int i,lead=(int)strlen(label);
	for(i=0;i<lead;i++) cout << " ";
#ifdef MPM_CODE
	cout << "(" << t->xx << ", " << t->xy << "," << t->xz << ")" << endl;
	cout << label << "(" << t->xy << ", " << t->yy << "," << t->yz << ")" << endl;
#else
	cout << "(" << t->xx << ", " << t->xy << ", 0.0)" << endl;
	cout << label << "(" << t->xy << ", " << t->yy << ", 0.0)" << endl;
#endif
	for(i=0;i<lead;i++) cout << " ";
#ifdef MPM_CODE
	cout << "(" << t->xz << ", " << t->yz << "," << t->zz << ")" << endl;
#else
	cout << "( 0.0, 0.0," << t->zz << ")" << endl;
#endif
}

#ifdef MPM_CODE

TensorAntisym *ZeroTensorAntisym(TensorAntisym *a)
{	a->xy=a->xz=a->yz=0.;
	return a;
}

#endif

#pragma mark String Functions

// convert character to lower case
unsigned char ConvertToLower(unsigned char ch)
{	if (ch >= 'A' && ch <= 'Z')
		ch = 'a' + (ch - 'A');
	return ch;
}

// convert path with / to DOS path with back slashes
// return pointer to input string
char *MakeDOSPath(char *unixPath)
{	for (unsigned int i = 0; i < strlen(unixPath); i++)
	{	if (unixPath[i] == '/')
			unixPath[i] = '\\';
	}
	return unixPath;
}

// case insensitive compare
// return 0 if match, return <0 or >0 if first case insentive mismatch character
//      has lower value in s1 or in s2, respectively
int CIstrcmp(const char *s1, const char *s2)
{
	unsigned char *us1 = (unsigned char *)s1;
	unsigned char *us2 = (unsigned char *)s2;
	
	while(ConvertToLower(*us1) == ConvertToLower(*us2++))
	{	// they are equal, so still matching
		
		// if matching at the end, then found a match
		if(*us1 == '\0') return (0);
		
		// otherwise on to next character
		us1++;
	}
	
	// found mismatch
	us2--;
	return (int)(ConvertToLower(*us1) - ConvertToLower(*us2));
}

// file extension of last component in the path (or empty string if none)
void GetFileExtension(const char *fileName,char *ext,int maxLength)
{
	int endLoc = (int)strlen(fileName)-1;
	int dotLoc = endLoc;
	while(dotLoc>=0 && fileName[dotLoc]!='.' && fileName[dotLoc]!='/') dotLoc--;
	
	// found an extension
	ext[0] = 0;
	if(dotLoc>=0)
	{	// empty if hit folder divider
		if(fileName[dotLoc]=='/') return;
		
		// first character after the period
		dotLoc++;
		int ei = 0;
		
		// last one gets the zero terminator
		while(ei<maxLength && dotLoc<=endLoc)
			ext[ei++]=fileName[dotLoc++];
		ext[ei] = 0;
	}
}

#pragma mark Other Function

/****************************************************************
 *  Functions for an intersect of plane and  unit cube
 ****************************************************************/
// Returns a vector with coordinates of vertex
Vector MakeCubeCorner(int vertex){
	Vector corner;
	switch(vertex){
		case(0): 
			corner.x = 0.5;
			corner.y = 0.5;
			corner.z = 0.5;			
			break;
		case(1): 			
			corner.x = 0.5;
			corner.y = -0.5;
			corner.z = 0.5;			
			break;
		case(2): 
			corner.x = 0.5;
			corner.y = 0.5;
			corner.z = -0.5;			
			break;
		case(3): 
			corner.x = -0.5;
			corner.y = 0.5;
			corner.z = 0.5;			
			break;
		case(4): 
			corner.x = -0.5;
			corner.y = -0.5;
			corner.z = 0.5;			
			break;
		case(5): 			
			corner.x = 0.5;
			corner.y = -0.5;
			corner.z = -0.5;			
			break;
		case(6): 			
			corner.x = -0.5;
			corner.y = 0.5;
			corner.z = -0.5;			
			break;
		case(7): 
			corner.x = -0.5;
			corner.y = -0.5;
			corner.z = -0.5;			
			break;
	}
	return corner;
}
// Do these edges intersect?
// Inspired by paper: "A Vertex Program for Efficient Box-Plane Intersection"
// by Christof Rezk Salama and Andreas Kolb, 2005
bool IntersectEdges(int i,int j,Vector *Point,Vector *n)
{
// Vertices of edge ij
Vector Vertex_i = MakeCubeCorner(i);  // vertex i
Vector Edge_ij = MakeCubeCorner(j);   // vertex j
Edge_ij.x -=Vertex_i.x;    // difference of coordinates
Edge_ij.y -=Vertex_i.y;
Edge_ij.z -=Vertex_i.z;
// Start calculating lambda
double lambda = -DotVectors(&Vertex_i,n);  // numerator of lambda
double denom = DotVectors(&Edge_ij,n); // denominator of lambda

// Don't divide by zero
if(DbleEqual(denom,0.0))
{
	return false;  // colinear
}

// Calculate lambda
lambda = lambda/denom;
// if not in [0,1] then it doesn't intersect
if(lambda>1.0 || lambda<0.0){
	return false;  // doesn't intersect
}
// Find point
Point->x = Vertex_i.x +lambda*Edge_ij.x;
Point->y = Vertex_i.y +lambda*Edge_ij.y;
Point->z = Vertex_i.z +lambda*Edge_ij.z;
return true; // it does intersect
}


// get the area of  the  polygon
double AreaOverVolume3D(Vector *norm,double dx,double dy,double dz)
{
	// define stuff
	Vector Polygon[6];
	double area =0.0;
	int j,k,v;
	bool intersect;
	// Get the polygon from the intersection
	// Also inspired by paper: "A Vertex Program for Efficient Box-Plane Intersection"

	// Get the first point P0
	v=0;
	if(!IntersectEdges(0,1,&Polygon[v],norm)){
		if(!IntersectEdges(1,4,&Polygon[v],norm)){
			intersect = IntersectEdges(4,7,&Polygon[v],norm);
			
		}
	}
	v++;
	// P1 Second point (possibily)
	intersect = IntersectEdges(1,5,&Polygon[v],norm);
	if(intersect) v++;  // if it does intersect go the next point

	// Get P2
	if(!IntersectEdges(0,2,&Polygon[v],norm)){
		if(!IntersectEdges(2,5,&Polygon[v],norm)){
			intersect = IntersectEdges(5,7,&Polygon[v],norm);
		}
	}
	v++;

	//Get P3 (if it is there)
	intersect = IntersectEdges(2,6,&Polygon[v],norm);
	if(intersect) v++;  // if it does intersect go the next point

	// Get P4
	if(!IntersectEdges(0,3,&Polygon[v],norm)){
		if(!IntersectEdges(3,6,&Polygon[v],norm)){
			intersect = IntersectEdges(6,7,&Polygon[v],norm);
		}
	}
	v++;

	//Get P5 (if it is there)
	intersect = IntersectEdges(3,4,&Polygon[v],norm);
	if(intersect) v++;  // if it does intersect go the next point
	
	// Find which coordinate project out
	double nx = fabs(norm->x);
	double ny = fabs(norm->y);
	double nz = fabs(norm->z);
	double norm_n = sqrt(nx*nx+ny*ny+nz*nz);
	// Loop through and find area of polygon
	if(nx>ny&&nx>nz){  // project out x
		 // area of polygon projected onto z-y
		for(int i =0;i<v;i++){
			j = (i+1)%v;
			k = (i+2)%v;
			area += Polygon[j].y*(Polygon[k].z-Polygon[i].z);
		}
		area = fabs(area*norm_n/(2.0 * nx )); // correction factor for projection
		area /= dx; // Convert from unit cube and divide by volume  
	}else if(ny>nx&&ny>nz){// project out y
		// area of polygon projected onto z-x
		for(int i =0;i<v;i++){
			j = (i+1)%v;
			k = (i+2)%v;
			area += Polygon[j].z*(Polygon[k].x-Polygon[i].x);
		}
		area = fabs(area*norm_n/(2.0 * ny )); // correction factor for projection
		area /= dy; // Convert from unit cube and divide by volume 
	}else{
		// area of polygon projected onto x-y
		for(int i =0;i<v;i++){
			j = (i+1)%v;
			k = (i+2)%v;
			area += Polygon[j].x*(Polygon[k].y-Polygon[i].y);
		}
		area = fabs(area*norm_n/(2.0 * nz)); // correction factor for projection
		area /= dz; // Convert from unit cube and divide by volume 
	}
	return area;

}



