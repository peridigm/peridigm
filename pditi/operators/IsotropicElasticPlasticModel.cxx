/*
 * IsotropicElasticPlasticModel.cxx
 *
 *  Created on: Jul 19, 2010
 *      Author: jamitch
 */

#include "IsotropicElasticPlasticModel.h"
#include "PdMaterialUtilities.h"
#include "PdITI_Utilities.h"
#include <iostream>
using std::cout;
using std::endl;
using PdImp::IsotropicHookeSpec;

namespace PdITI {

vector<FieldSpec> IsotropicElasticPlasticModel::registerTemporalBondVariables() const {
	return fieldSpecs;
}

IsotropicElasticPlasticModel::IsotropicElasticPlasticModel(IsotropicElasticPlasticSpec spec)
:
ConstitutiveModel(),
fieldSpecs(2,Field_NS::DEVIATORIC_PLASTIC_EXTENSION),
matSpec(spec)
{
/*
 * Initialized above fieldSpecs with DEVIATORIC_PLASTIC_EXTENSION
 * Need to set second field spec to 'Lambda' -- consistency parameter
 */
	fieldSpecs[1]=Field_NS::LAMBDA;
}


void IsotropicElasticPlasticModel::computeInternalForce
(
   const double* xOverlapPtr,
   const double* yOverlapPtr,
   const double* mOwned,
   const double* volumeOverlapPtr,
   const double* dilatationOwned,
   const double* bondDamage,
   const double* dsf,
   double* fInternalOverlapPtr,
   const int*  localNeighborList,
   int numOwnedPoints,
   vector< TemporalField<double> > temporalBondFields
){

	/*
	 * COMPUTE INTERNAL FORCE
	 */
	IsotropicHookeSpec hookSpec = matSpec.getHookeSpec();
	double K  = hookSpec.bulkModulus();
	double MU = hookSpec.shearModulus();
	double HORIZON=matSpec.materialHorizon();
	double YIELD_STRESS=matSpec.yieldStress();
	Field<double> fieldN   = temporalBondFields[0].getField(Field_NS::FieldSpec::STEP_N);
	Field<double> fieldNp1 = temporalBondFields[0].getField(Field_NS::FieldSpec::STEP_NP1);
	const double* deviatoricPlasticExtensionStateN   = fieldN.getArray().get();
	double* deviatoricPlasticExtensionStateNp1 = fieldNp1.getArray().get();
	Field<double> lambdaFieldN   = temporalBondFields[1].getField(Field_NS::FieldSpec::STEP_N);
	Field<double> lambdaFieldNP1 = temporalBondFields[1].getField(Field_NS::FieldSpec::STEP_NP1);
	const double* lambdaN   = lambdaFieldN.getArray().get();
	double* lambdaNP1 = lambdaFieldNP1.getArray().get();


	PdMaterialUtilities::computeInternalForceIsotropicElasticPlastic
	(
			xOverlapPtr, //
			yOverlapPtr, //
			mOwned,      //
			volumeOverlapPtr, //
			dilatationOwned,  //
			bondDamage,       //
			deviatoricPlasticExtensionStateN, //
			deviatoricPlasticExtensionStateNp1, //
			lambdaN, //
			lambdaNP1, //
			fInternalOverlapPtr, //
			localNeighborList, //
			numOwnedPoints, //
			K,  //
			MU, //
			HORIZON,
			YIELD_STRESS
	);

}

/*
 **Linearized tangent stiffness about the displacement field 'u'
 * @param I, P, Q are local indices
 * @param x -- overlap pointer to original mesh coordinates
 * @param u -- overlap pointer to displacement field
 * @param volume -- overlap pointer to cell volumes
 * @param m -- overlap pointer to weighted volumes
 * @param dilatation -- overlap pointer to pre-computed dilatation for displacement u
 * @param k 3x3 stiffness associated with these three points
 * @return k
 */
double* IsotropicElasticPlasticModel::
kIPQ3x3(
		int I,
		int P,
		int Q,
		const double* x,
		const double* u,
		const double* volume,
		const double *m,
		const double *dilatation,
		double *k,
		double horizon,
		double dsf_I
) {

	IsotropicHookeSpec hookSpec = matSpec.getHookeSpec();
	double K  = hookSpec.bulkModulus();
	double MU = hookSpec.shearModulus();

	/*
	 * Initialize stiffness
	 */
	SET(k,k+9,0.0);

	/*
	 * Construct bonds
	 */
	// Unit vectors along deformation states
	double MIP[3];
	double MIQ[3];
	// magnitudes of bonds and deformation states
	double  IP = 0.0,  IQ = 0.0;
	double YIP = 0.0, YIQ = 0.0;
	{
		/*
		 * Get points I, P, Q
		 */
		const int Ix3 = 3*I;
		const int Px3 = 3*P;
		const int Qx3 = 3*Q;
		const double *xI = x+Ix3;
		const double *xP = x+Px3;
		const double *xQ = x+Qx3;
		double bondIP[3]; BOND(xI,xP,bondIP);
		double bondIQ[3]; BOND(xI,xQ,bondIQ);
		IP = MAGNITUDE(bondIP);
		IQ = MAGNITUDE(bondIQ);
		if(IP > horizon || IQ > horizon) return k;

		/*
		 * Construct geometry for I, P, Q
		 */

		const double *uI = u+Ix3;
		const double *uP = u+Px3;
		const double *uQ = u+Qx3;
		double yI[3], yP[3], yQ[3];
		UPDATE_GEOMETRY(xI,uI,yI);
		UPDATE_GEOMETRY(xP,uP,yP);
		UPDATE_GEOMETRY(xQ,uQ,yQ);

		/*
		 * Normalized Deformation State for I and P;
		 * Normalized Deformation State for I and Q
		 */
		BOND(yI,yP,MIP);
		BOND(yI,yQ,MIQ);
		YIP = NORMALIZE(MIP);
		YIQ = NORMALIZE(MIQ);
	}

	double OMEGA_IP=1.0;
	double OMEGA_IQ=1.0;

	/*
	 * Volume associated with P and Q
	 */
	double volP = *(volume+P);
	double volQ = *(volume+Q);

	/*
	 * Weighted volume for I
	 */
	double mI = *(m+I);


	/*
	 * Outer product of M(Y)<IP> X M(Y)<IQ>
	 * ** 3x3 matrix
	 * ** Elements in a column are contiguous (column major)
	 * ** Interface with BLAS means that LDA=3 -- This is the stride between consecutive elements of the same row
	 */
	TENSOR_PRODUCT(MIP,MIQ,k);
	double c = 9*K-15*MU;
	c *= OMEGA_IP * OMEGA_IQ * IP * IQ * volP * volQ / mI / mI;
	SCALE_BY_VALUE(k,k+9,c);
	if(P!=Q) return k;

	/*
	 * P==Q from here on out
	 */
	c = (15 * MU / mI) * OMEGA_IQ * volQ;
	double k2[9];
	TENSOR_PRODUCT(MIQ,MIQ,k2);
	SCALE_BY_VALUE(k2,k2+9,c);
	SUMINTO(k2,k2+9,k);

	/*
	 * Geometric Stiffness;
	 * ONLY for a single bond (IP == IQ)
	 */

	/*
	 * Geometric stiffness is associated with the force state
	 * in a particular bond
	 */
	{
		double tI(0.0);
		double MIQxMIQ[9];
		TENSOR_PRODUCT(MIQ,MIQ,MIQxMIQ);
		double thetaI = *(dilatation+I);
		double c = OMEGA_IQ * (3.0 * K - 5.0 * MU ) / mI;
		tI = (c * thetaI * IP + OMEGA_IQ * (15.0 * MU / mI) * (YIQ - IQ));
		//		std::cout << "tI = " << tI << std::endl;
		//		std::cout << "theta = " << theta << std::endl;
		//		std::cout << "m = " << mI << std::endl;
		//		std::cout << "Bulk Modulus K = " << K << std::endl;
		//		std::cout << "tI = 9*K*(YIP-IP)/m = " << 9.0*K*(YIP-IP)/mI << std::endl;
		//		std::cout << "theta = 3*(Y/X - 1) = " << 3.0*(YIP/IP-1) << std::endl;
		//		std::cout << "m = Vol * X * X = " << volP*IP*IP << std::endl;
		SCALE_BY_VALUE(MIQxMIQ,MIQxMIQ+9,-1);
		MIQxMIQ[0]+=1.0;
		MIQxMIQ[4]+=1.0;
		MIQxMIQ[8]+=1.0;
		SCALE_BY_VALUE(MIQxMIQ,MIQxMIQ+9,tI * volP / YIP);
		SUMINTO(MIQxMIQ,MIQxMIQ+9,k);
	}


	return k;
}


}
