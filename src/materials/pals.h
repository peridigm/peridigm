#include <math.h>
#include <vector>

#include "Peridigm_InfluenceFunction.hpp"



namespace MATERIAL_EVALUATION {

namespace PALS {

const int NUM_LAGRANGE_MULTIPLIERS=6;

typedef typename PeridigmNS::InfluenceFunction::functionPointer FunctionPointer;

/**
 *
 * Following 'lambda' functions are a much cleaner implementation of the influence functions,
 * but this approach requires c++11.

// #include <functional>
std::function<double (double*)> get_gaussian(double horizon) {
   return [horizon](double *xi) {
      double dx=*(xi+0);
      double dy=*(xi+1);
      double dz=*(xi+2);
      double h2=horizon*horizon;
      double xi2=dx*dx+dy*dy+dz*dz;
      return exp(-xi2/h2);
   };
}
std::function<double (double*)> get_constant(double c) {
   return [c](double *xi) {  return c; };
}
   std::function<double (double*)> g=get_gaussian(1.0);
   std::function<double (double*)> c=get_constant(2.0);
   double zero[]={0.0,0.0,0.0};
   cout << std::scientific << std::setw(15) << g(zero) << endl;
   cout << std::scientific << std::setw(15) << c(zero) << endl;

*/

struct dilatation_influence {
   static double eval(const double *bond, const double *weights){
      const double *w=weights;
      double a=bond[0], b=bond[1], c=bond[2];
      return w[0]*a*a+w[1]*b*b+w[2]*c*c+2.0*(w[3]*a*b+w[4]*a*c+w[5]*b*c);
   }
};

struct deviatoric_influence {
	static double eval(const double *bond, const double *weights){
		const double *w=weights;
		double a=bond[0], b=bond[1], c=bond[2];
		double a2=a*a, b2=b*b, c2=c*c, ab=a*b, ac=a*c, bc=b*c;
		double x2=a2+b2+c2;
		double r=std::sqrt(x2);
		double avg=r/3.0;
		double ek[]={a2/r-avg,b2/r-avg,c2/r-avg,2*ab/r,2*ac/r,2*bc/r};
		return w[0]*ek[0]*ek[0]+
			   w[1]*ek[1]*ek[1]+
			   w[2]*ek[2]*ek[2]+
			   w[3]*ek[3]*ek[3]+
			   w[4]*ek[4]*ek[4]+
			   w[5]*ek[5]*ek[5];
   }
};

template<typename T>
struct pals_influence {

   pals_influence(FunctionPointer p, double _c, const double *w):
      I0(p), normalize(_c), weights(w) {}

   double operator()(const double *bond, double horizon){
      double a=bond[0], b=bond[1], c=bond[2];
      double x=std::sqrt(a*a+b*b+c*c);
      /*
       * NOTE: 'normalize is distributed across entire influence
       * function -- both the 'initial' influence function and
       * the Lagrange multipliers piece.
       */
      return normalize*(I0(x,horizon)+T::eval(bond,weights));
   }

   // Initial influence function
   const FunctionPointer I0;
   // Normalizing scalar for I0
   const double normalize;
   // Lagrange multipliers
   const double *weights;
};


double 
compute_normalizing_constant_point
(
 struct pals_influence<dilatation_influence>& OMEGA,
 const double *X,
 const double *xOverlap,
 const double *volumeOverlap,
 const int *neigh,
 double horizon
);

double 
compute_normalizing_constant_point
(
 struct pals_influence<deviatoric_influence>& SIGMA,
 const double *X,
 const double *xOverlap,
 const double *volumeOverlap,
 const int *neigh,
 double horizon
);


void
compute_lagrange_multipliers
(
	const double *xOverlap,
	const double *volumeOverlap,
	int num_owned_points,
	const int *localNeighborList,
	double horizon,
	std::vector<double *>& sig_owned,
	double *beta_sig_owned,
	std::vector<double *>& tau_owned,
	double *beta_tau_owned,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0
);

void
compute_lagrange_multipliers_point
(
	const double *X,
	const double *xOverlap,
	const double *volumeOverlap,
	const int* neigh,
	double h,
	double *omega_multipliers,
	double *omega_constant,
	double *sigma_multipliers,
	double *sigma_constant,
	const FunctionPointer dilation_influence_function,
	const FunctionPointer deviatoric_influence_function
);

/*
 * Computes a weighted volume but with the pals deviatoric influence function
 */
void computeWeightedVolume
(
	const double *xOverlap,
	const double *volumeOverlap,
	const std::vector<const double *>& _sigma_multipliers,
	const double *sigma_constant,
	double *weighted_volume,
	int myNumPoints,
	const int* localNeighborList,
	double horizon,
	const FunctionPointer SIGMA_0=PeridigmNS::InfluenceFunction::self().getInfluenceFunction()
);

/*
 * Computes 'normalized' weighted volume
 * Sanity check: all values should be == 3.0
 */
void computeNormalizedWeightedVolume
(
	const double *xOverlap,
	const double *volumeOverlap,
	const double *omega_constant,
	const double *bondDamage,
	double *weighted_volume,
	int myNumPoints,
	const int* localNeighborList,
	double horizon,
	const FunctionPointer SIGMA_0=PeridigmNS::InfluenceFunction::self().getInfluenceFunction()
);

void computeDilatation
(
	const double *xOverlap,
	const double *yOverlap,
	const double *volumeOverlap,
	const std::vector<const double *>& _omega_multipliers,
	const double *omega_constant,
	double *dilatation,
	const int *localNeighborList,
	int numOwnedPoints,
	double horizon,
	const FunctionPointer OMEGA_0
);

void computeDilatationAndPalsPressure
(
	const double *xOverlap,
	const double *yOverlap,
	const double *volumeOverlap,
	const std::vector<const double *>& _omega_multipliers,
	const double *omega_constant,
	const std::vector<const double *>& _sigma_multipliers,
	const double *sigma_constant,
	const double *weighted_volume,
	double *dilatation,
	double *pals_pressure,
	const int *localNeighborList,
	int numOwnedPoints,
	double BULK_MODULUS,
	double SHEAR_MODULUS,
	double horizon,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0
);

void computeInternalForcePals
(
	const double *xOverlap,
	const double *yOverlap,
	const double *volumeOverlap,
	const std::vector<const double *>& _omega_multipliers,
	const double *omega_constant,
	const std::vector<const double *>& _sigma_multipliers,
	const double *sigma_constant,
	const double *dilatation,
	const double *pals_pressure,
	double *fInternalOverlap,
	const int *localNeighborList,
	int numOwnedPoints,
	double BULK_MODULUS,
	double SHEAR_MODULUS,
	double horizon,
	const FunctionPointer OMEGA_0,
	const FunctionPointer SIGMA_0
);

}

}


