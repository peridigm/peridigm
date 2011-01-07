/*
 * PimpIsotropicHookeSpec.h
 *
 *  Created on: Feb 16, 2010
 *      Author: jamitch
 */

#ifndef PDIMP_MATERIALS_H_
#define PDIMP_MATERIALS_H_

#include "Property.h"

namespace PdImp {

class YoungsModulus : public PositiveNumber<double> {
public:
	explicit YoungsModulus(double e) : PositiveNumber<double>(e, "YoungsModulus"){}
};

class YieldStress : public PositiveNumber<double> {
public:
	explicit YieldStress(double eY) : PositiveNumber<double>(eY, "YieldStress"){}
};

class MaterialHorizon : public PositiveNumber<double> {
public:
	explicit MaterialHorizon(double radius) : PositiveNumber<double>(radius, "MaterialHorizon"){}
};

class PoissonsRatio : public DomainNumber<double> {
public:
	explicit PoissonsRatio(double nu) : DomainNumber<double>(Inclusive<double>(0.0), Exclusive<double>(0.5), nu, "PoissonsRatio"){
		if(nu>=.5){
			std::string message("PoissonsRatio");
			message += " must be < .5; Input value = ";
			message += Property<double>::toString(nu);
			throw std::domain_error(message);
		}
	}
};

class BulkModulus : public PositiveNumber<double> {
public:
	explicit BulkModulus(double k) : PositiveNumber<double>(k, "BulkModulus"){}
};

class ShearModulus : public PositiveNumber<double> {
public:
	explicit ShearModulus(double g) : PositiveNumber<double>(g, "ShearModulus"){}
};

class MassDensity : public PositiveNumber<double> {
public:
	explicit MassDensity(double rho) : PositiveNumber<double>(rho, "MassDensity"){}
};

class IsotropicHookeSpec {

public:
	static YoungsModulus youngsModulus(double e) {return YoungsModulus(e);}
	static PoissonsRatio poissonsRatio(double nu) {return PoissonsRatio(nu);}
	static BulkModulus bulkModulus(double k) {return BulkModulus(k);}
	static ShearModulus shearModulus(double g) {return ShearModulus(g);}

	explicit IsotropicHookeSpec(BulkModulus k, PoissonsRatio poissonsRatio)
	: E(3.0*k()*(1-2.0*poissonsRatio())),
	  nu(poissonsRatio) {}

	explicit IsotropicHookeSpec(YoungsModulus youngsModulus, PoissonsRatio poissonsRatio)
	: E(youngsModulus),
	  nu(poissonsRatio) {}

	const double youngsModulus() const { return E.getValue(); }
	double poissonsRatio() const { return nu.getValue(); }
	double shearModulus()  const { return E.getValue()/(2.*(1.0 + nu.getValue())); }
	double lameConstant()  const {
		double e = E.getValue();
		double n = nu.getValue();
		return -((e*n)/((1 + n)*(-1 + 2*n)));
	}
	double bulkModulus() const {
		double lambda = lameConstant();
		double G = shearModulus();
		return lambda+2.0*G/3.0;
	}

private:
	const YoungsModulus E;
	const PoissonsRatio nu;
};

class IsotropicElasticPlasticSpec {
public:
	static YieldStress yieldStress(double Y) {return YieldStress(Y);}
	static MaterialHorizon materialHorizon(double radius) {return MaterialHorizon(radius);}
	explicit IsotropicElasticPlasticSpec(YieldStress Y, MaterialHorizon r, IsotropicHookeSpec H)
	: yield(Y),
	  delta(r),
	  hookeSpec(H) {}
	const IsotropicHookeSpec& getHookeSpec() const { return hookeSpec; }
	double materialHorizon() const { return delta.getValue(); }
	double yieldStress() const { return yield.getValue(); }
private:
	const YieldStress yield;
	const MaterialHorizon delta;
	const IsotropicHookeSpec hookeSpec;

};

}

#endif /* PDIMP_MATERIALS_H_ */
