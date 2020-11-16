#include "pals.h"
#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

namespace MATERIAL_EVALUATION {

namespace PALS {

//Lapack linear equations
//http://www.netlib.org/lapack/lug/node38.html

// Factorize general matrix
// http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f_source.html

// Solve general matrix
// http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f_source.html

//Cholesky factorize (dpotrf)
//http://www.netlib.org/lapack/explore-html/df/da8/VARIANTS_2cholesky_2RL_2dpotrf_8f_source.html

//Solve using factorize (dpotrs)
//http://www.math.utah.edu/software/lapack/lapack-d/dpotrs.html

extern "C" {
  /*
   * Lapack
   * Cholesky factorize dense matrix
   * Cholesky backsolve dense matrix
   */
  void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
  void dpotrs_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);

  /*
   * Lapack
   * Factorize general 'dense' matrix
   */
  void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
  void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
}

using std::setw;
using std::vector;


void print_point_3d(std::ostream& out, int local_id, const double *x){

  out << std::setprecision(5);
  out <<  std::scientific << "\t" << setw(5) << local_id  << ": "
  << setw(13) << x[0] << ", "
  << setw(13) << x[1] << ", "
  << setw(13) << x[2] << ", " << "\n";

}

void print_vector3d(std::ostream& out, const std::string &label, const double *rhs){

  out << label << "\n";
  out <<  std::scientific << "\t" << std::setw(13) << rhs[0] << ", ";
  out <<  std::scientific << "\t" << std::setw(13) << rhs[1] << ", ";
  out <<  std::scientific << "\t" << std::setw(13) << rhs[2] << "\n";

}

void print_N_vector(std::ostream& out, const std::string &label, int N, const double *rhs){

  out << label << "\n";
  for(int n=0;n<N;n++)
    out << std::scientific << std::setw(13) << rhs[n] << " ";
  out << "\n";
}

void print_symmetrix_3x3(std::ostream& out, const std::string &label, const double *k){

  /*
   * Print symmetric matrix
    double k11=k[0];
    double k22=k[1];
    double k33=k[2];
    double k12=k[3];
    double k23=k[4];
    double k13=k[5];
   */
  out << label << "\n";
  out <<  std::scientific << "\t" << setw(13) << k[0] << " " << setw(13) << k[3] << " " << setw(13) << k[5] << "\n";
  out <<  std::scientific << "\t" << setw(13) << k[3] << " " << setw(13) << k[1] << " " << setw(13) << k[4] << "\n";
  out <<  std::scientific << "\t" << setw(13) << k[5] << " " << setw(13) << k[4] << " " << setw(13) << k[2] << "\n";

}

void print_symmetrix_6x6(std::ostream& out, const std::string &label, const double *k){

  out << label << "\n";
  for(int r=0,p=0;r<NUM_LAGRANGE_MULTIPLIERS;r++){
    for(int c=0;c<NUM_LAGRANGE_MULTIPLIERS;c++,p++){
      out << std::scientific;
      out << setw(13) << k[p] << " ";
    }
    out << "\n";
  }
}

void
solve_symmetric_3x3(const double *k, const double *rhs, double *solution)
{

  double *w=solution;

  /*
   * Symmetric matrix
   */
  double k11=k[0];
  double k22=k[1];
  double k33=k[2];
  double k12=k[3];
  double k23=k[4];
  double k13=k[5];

  /*
   * Determinant
   */
  double det = -k13 * k13 * k22 + 2 * k12 * k13 * k23 - k11 * k23 * k23 - k12 * k12 * k33 + k11 * k22 * k33;

  /*
   * Each row of inverse
   */
  vector<double> row1(3), row2(3), row3(3);
  row1[0] = -k23 * k23 + k22 * k33; row1[1] =  k13 * k23 - k12 * k33; row1[2] = -k13 * k22 + k12 * k23;
  row2[0] =  k13 * k23 - k12 * k33; row2[1] = -k13 * k13 + k11 * k33; row2[2] =  k12 * k13 - k11 * k23;
  row3[0] = -k13 * k22 + k12 * k23; row3[1] =  k12 * k13 - k11 * k23; row3[2] = -k12 * k12 + k11 * k22;

  /*
   * Solution
   */
  w[0]=(row1[0]*rhs[0]+row1[1]*rhs[1]+row1[2]*rhs[2])/det;
  w[1]=(row2[0]*rhs[0]+row2[1]*rhs[1]+row2[2]*rhs[2])/det;
  w[2]=(row3[0]*rhs[0]+row3[1]*rhs[1]+row3[2]*rhs[2])/det;

}

void solve_linear_problem(double *k, double *rhs){

  /*
   * Factorize k
   */
  int N=NUM_LAGRANGE_MULTIPLIERS;
  int INFO, M=N, LDA=N;
  int IPIV[NUM_LAGRANGE_MULTIPLIERS];
  dgetrf_(&M,&N,k,&LDA,IPIV,&INFO);
  if(0!=INFO){
    std::string message="ERROR Pals model: compute_lagrange_multipliers \n\t";
    message+="Matrix factorization error -- bad value or factorization will yield division by zero during solve.";
    throw std::runtime_error(message);
  }
  /*
   * Solve
   */
  char TRANS='N';
  int NRHS=1, LDB=N;
  dgetrs_(&TRANS,&N,&NRHS,k,&LDA,IPIV,rhs,&LDB,&INFO);
  if(0!=INFO){
    std::string message="ERROR Pals model: compute_lagrange_multipliers \n\t";
    message+="Matrix solve error -- bad value.";
    throw std::runtime_error(message);
  }
}

void
compute_lagrange_multipliers
(
  const double *xOverlap,
  const double *volumeOverlap,
  int num_owned_points,
  const int *localNeighborList,
  double horizon,
  vector<double *>& omega_multipliers,
  double *omega_constants,
  vector<double *>& sigma_multipliers,
  double *sigma_constants,
  const FunctionPointer OMEGA_0,
  const FunctionPointer SIGMA_0
)
{

  const double *X=xOverlap;
  const int *neigh_X=localNeighborList;
  int num_neigh_X=*neigh_X;
  double *oc_X=omega_constants;
  double *sc_X=sigma_constants;
  double omega_X[NUM_LAGRANGE_MULTIPLIERS], sigma_X[NUM_LAGRANGE_MULTIPLIERS];
  for(int p=0;p<num_owned_points;p++,X+=3,oc_X++,sc_X++){

    // evaluate lagrange multipliers and normalizing constants for point X
    //print_point_3d(std::cout,p,X);
    compute_lagrange_multipliers_point(X,xOverlap,volumeOverlap,neigh_X,horizon,omega_X,oc_X,sigma_X,sc_X,OMEGA_0,SIGMA_0);

    // Collect computed Lagrange multipliers for this point
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      omega_multipliers[i][p]=omega_X[i];
      sigma_multipliers[i][p]=sigma_X[i];
    }
    //print_N_vector(std::cout,"lagrange_multipliers dilatation",6,sig_X);
    //print_N_vector(std::cout,"lagrange_multipliers deviatoric",6,tau_X);

    // move to next point
    num_neigh_X=*neigh_X;
    neigh_X+=(num_neigh_X+1);

  }

}


void
compute_lagrange_multipliers_point
(
  const double *X,
  const double *xOverlap,
  const double *volumeOverlap,
  const int *neighborhood,
  double h,
  double *omega_multipliers,
  double *omega_constant,
  double *sigma_multipliers,
  double *sigma_constant,
  const FunctionPointer OMEGA_0,
  const FunctionPointer SIGMA_0
)
{
  /*
   * Compute Lagrange multipliers for point X: 6 each (omega_multipliers, sigma_multipliers)
   * INPUT
   * X: 3D point
   * q: neighbors of X
   * vol_q: volume of each neighbor q
   * num_neigh: number of neighbors @ X
   * h: horizon
   * OUTPUT
   * omega_multipliers[6]: Lagrange multipliers for dilatation
   * sigma_multipliers[6]: Lagrange multipliers for deviatoric
   * omega_constant[1]: normalizing constant for trial dilatation influence function
   * sigma_constant[1]: normalizing constant for trial deviatoric influence function
   */

  /*
   * RHS side vectors of linear problems
   */
  double rhs_dil[NUM_LAGRANGE_MULTIPLIERS];
  double rhs_dev[NUM_LAGRANGE_MULTIPLIERS];

  // Matrix at point
  double k[NUM_LAGRANGE_MULTIPLIERS*NUM_LAGRANGE_MULTIPLIERS];
  double k_dev[NUM_LAGRANGE_MULTIPLIERS*NUM_LAGRANGE_MULTIPLIERS];
  // Initialize RHS vectors and Matrices
  for(int i=0,p=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
    rhs_dil[i]=0.0; rhs_dev[i]=0.0;
    for(int j=0;j<NUM_LAGRANGE_MULTIPLIERS;j++,p++){
      k[p]=0.0;
      k_dev[p]=0.0;
    }
  }

  /*
   * Normalize Gaussian influence function at point
   */
   const int *neigh=neighborhood;
  int num_neigh=*neigh; neigh++;
  double m(0.0),sum_dev(0.0);
  const double x=*X, y=*(X+1), z=*(X+2);
  for(int iq=0;iq<num_neigh;iq++){

    int neigh_local_id=*neigh; neigh++;
    const double *q=xOverlap+3*neigh_local_id;
    double vol=volumeOverlap[neigh_local_id];
    double bond[3]={*(q+0)-x,*(q+1)-y,*(q+2)-z};
    double a=bond[0], b=bond[1], c=bond[2];
    double a2=a*a, b2=b*b, c2=c*c;
    double ab=a*b, ac=a*c, bc=b*c;
    /*
     * Components of extension state
     * NOTE
     * Each of these extension state has a bond in the denominator;
     * But the dilatation matrix is formed by the product of the bond
     * with the extension state and hence the bond cancels.
     */
    double e_k[]={a2,b2,c2,2*ab,2*ac,2*bc};

    /*
     * Symmetric dilatation matrix;
     */
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      double e_i=e_k[i];
      int col=i*NUM_LAGRANGE_MULTIPLIERS;
      for(int j=0;j<NUM_LAGRANGE_MULTIPLIERS;j++){
        double e_j=e_k[j];
        double kij=e_i*e_j*vol;
        k[col+j]+=kij;
      }
    }

    /*
     * Components of deviatoric extension state
     */
    double xi2=a2+b2+c2;
    double r=std::sqrt(xi2);
    double avg=r/3.0;
    double epsilon_k[]={a2/r-avg,b2/r-avg,c2/r-avg,2*ab/r,2*ac/r,2*bc/r};

    /*
     * Symmetric deviatoric matrix;
     */
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      double epsilon_i=epsilon_k[i];
      double epsilon_i2=epsilon_i*epsilon_i;
      int col=i*NUM_LAGRANGE_MULTIPLIERS;
      for(int j=0;j<NUM_LAGRANGE_MULTIPLIERS;j++){
        double epsilon_j=epsilon_k[j];
        double epsilon_j2=epsilon_j*epsilon_j;
        double kij=epsilon_i2*epsilon_j2*vol;
        k_dev[col+j]+=kij;
      }
    }

    /*
     * Reference influence function evaluation
     */
    double omega_0=OMEGA_0(r,h);
    double sigma_0=SIGMA_0(r,h);

    // sums for normalization of OMEGA_0 and SIGMA_0
    m+=omega_0*xi2*vol;
    double epsilon_0=2*(ab+bc+ac)/r;
    sum_dev+=sigma_0*epsilon_0*epsilon_0*vol;

    /*
     * RHS vectors
     */
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      rhs_dil[i]+=omega_0*e_k[i]*vol;
      rhs_dev[i]+=sigma_0*epsilon_k[i]*epsilon_k[i]*vol;
    }
  }

   // normalizing scalars
   // const double D=3.0;
   // const double norm_omega_0=D/m;
   // const double norm_sigma_0=6.0/sum_dev;

  // finalize RHS vectors
  const double c_1=1.0, c_2=1.0;
  rhs_dil[0]=1.0-c_1*rhs_dil[0];
  rhs_dil[1]=1.0-c_1*rhs_dil[1];
  rhs_dil[2]=1.0-c_1*rhs_dil[2];
  rhs_dil[3]=0.0-c_1*rhs_dil[3];
  rhs_dil[4]=0.0-c_1*rhs_dil[4];
  rhs_dil[5]=0.0-c_1*rhs_dil[5];

  const double c23=2.0/3.0;
  rhs_dev[0]=c23-c_2*rhs_dev[0];
  rhs_dev[1]=c23-c_2*rhs_dev[1];
  rhs_dev[2]=c23-c_2*rhs_dev[2];
  rhs_dev[3]=2.0-c_2*rhs_dev[3];
  rhs_dev[4]=2.0-c_2*rhs_dev[4];
  rhs_dev[5]=2.0-c_2*rhs_dev[5];

  /*
   * ISSUE is that dx, dy, or dz may be zero;
   * For 2D meshes in x-y plane, this occurs because all of the
   * z-coordinate values have the same constant value and hence dz=0.
   * The linear problem in these cases is ill-defined; Following code fixes.
   */
  int i_diag[]={0,1*6+2-1,2*6+3-1,3*6+4-1,4*6+5-1,5*6+6-1};
  double trace_k=0;
  for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++)
    trace_k+=k[i_diag[i]];
  double small=1.0e-15;
  for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
    int d=i_diag[i];
    if(k[d]/trace_k < small) {
      /*
       * Zero out row/col for column i
       * Remember that k is column major
       * 'r' tracks across  row 'i'
       * 'c' tracks down column 'i'
       */
      int N=NUM_LAGRANGE_MULTIPLIERS;
      for(int j=0,r=i,c=i*N;j<N;j++,r+=N,c++){
        k[r]=0.0;
        k[c]=0.0;
        k_dev[r]=0.0;
        k_dev[c]=0.0;
      }
      // Set diagonal to 1.0
      k[d]=1.0;
      k_dev[d]=1.0;
      // Set rhs vectors to zero
      rhs_dil[i]=0.0;
      rhs_dev[i]=0.0;
    }
  }

  //print_symmetrix_6x6(std::cout,"K_DIL",k);
  //print_N_vector(std::cout,"RHS_DIL",NUM_LAGRANGE_MULTIPLIERS,rhs_dil);
  solve_linear_problem(k,rhs_dil);
  //print_N_vector(std::cout,"RHS_DIL",NUM_LAGRANGE_MULTIPLIERS,rhs_dil);

//  print_symmetrix_6x6(std::cout,"K_DEV",k_dev);
//  print_N_vector(std::cout,"RHS_DEV",NUM_LAGRANGE_MULTIPLIERS,rhs_dev);
    solve_linear_problem(k_dev,rhs_dev);
//  print_N_vector(std::cout,"RHS_DEV",NUM_LAGRANGE_MULTIPLIERS,rhs_dev);

  // Copy solution into output vectors
  for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
    omega_multipliers[i]=rhs_dil[i];
    sigma_multipliers[i]=rhs_dev[i];
  }

  /*
   * This code computes scaling factors for the composite influence
   * functions.  Generally, it turns out the scaling for 'OMEGA' is 1.0
   * but scaling for SIGMA is a non-trivial number.
   * ALSO -- note that the 'pals_influence' function definition can
   * be adjusted to scale on OMEGA_0 and SIGMA_0 only.  Its a bit
   * subtle here -- the functions below are computed on the assumption
   * that both parts of the influence function are scaled, ie, these
   * are computed using the 'composite/whole' function.  Assuming that
   * there are no bugs, it turns out that scaling 'SIGMA' does
   * not always work;  Turning scaling off has no effect on 'OMEGA' since 
    * the scaling is implicitly correct with the choice of 'matching 
    * deformations' used here.  However, scaling on the deviatoric 
    * piece does not generally work but it does work when scaling is 
    * not used.
   */
   //const double one(1.0);
   //pals_influence<dilatation_influence> OMEGA(OMEGA_0,one,rhs_dil);
   //pals_influence<deviatoric_influence> SIGMA(SIGMA_0,one,rhs_dev);
   //const double norm_omega=compute_normalizing_constant_point(OMEGA,X,xOverlap,volumeOverlap,neighborhood,h);
   //const double norm_sigma=compute_normalizing_constant_point(SIGMA,X,xOverlap,volumeOverlap,neighborhood,h);
  //*omega_constant=norm_omega;
  //*sigma_constant=norm_sigma;

   /*
    * Turn scaling off;
    */
  *omega_constant=1.0;
  *sigma_constant=1.0;

}


double compute_normalizing_constant_point
(
    struct pals_influence<dilatation_influence>& OMEGA,
    const double *X,
    const double *xOverlap,
    const double *volumeOverlap,
    const int *neigh,
    double horizon
) {
   double m=0.0;
   const double D=3.0;
  int num_neigh=*neigh; neigh++;
  const double x=*X, y=*(X+1), z=*(X+2);
  for(int iq=0;iq<num_neigh;iq++){
    int neigh_local_id=*neigh; neigh++;
    const double *q=xOverlap+3*neigh_local_id;
    double vol=volumeOverlap[neigh_local_id];
    double bond[3]={*(q+0)-x,*(q+1)-y,*(q+2)-z};
    double omega = OMEGA(bond,horizon);
    double a = bond[0];
    double b = bond[1];
    double c = bond[2];
    double xi2 = a*a+b*b+c*c;
    m+=omega*xi2*vol;
   }
   return D/m;
}

double 
compute_normalizing_constant_point
(
    struct pals_influence<deviatoric_influence>& SIGMA,
    const double *X,
    const double *xOverlap,
    const double *volumeOverlap,
    const int *neigh,
    double horizon
) {
  double norm=0.0;
  const double D=6.0;
  int num_neigh=*neigh; neigh++;
  const double x=*X, y=*(X+1), z=*(X+2);
  for(int iq=0;iq<num_neigh;iq++){
    int neigh_local_id=*neigh; neigh++;
    const double *q=xOverlap+3*neigh_local_id;
    double vol=volumeOverlap[neigh_local_id];
    double bond[3]={*(q+0)-x,*(q+1)-y,*(q+2)-z};
    double sigma =SIGMA(bond,horizon);
    double a = bond[0];
    double b = bond[1];
    double c = bond[2];
    double ab=a*b, ac=a*c, bc=b*c;
    double xi2 = a*a+b*b+c*c;
    double r=std::sqrt(xi2);
    double epsilon=2*(ab+bc+ac)/r;
    norm+=sigma*epsilon*epsilon*vol;
   }
   return D/norm;
}


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
    int numOwnedPoints,
    const int* localNeighborList,
    double horizon,
    const FunctionPointer OMEGA_0
)
{
  double bond[3];
  const double *xOwned = xOverlap;
  const double *oc=omega_constant;
  double *m=weighted_volume;
  const int *neighPtr = localNeighborList;
  double cellVolume;
  for(int q=0; q<numOwnedPoints;q++, xOwned+=3, oc++, m++){
    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    *m=0.0;

    for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
      int localId = *neighPtr;
      cellVolume = volumeOverlap[localId];
      const double *XP = &xOverlap[3*localId];
      bond[0]=XP[0]-X[0];
      bond[1]=XP[1]-X[1];
      bond[2]=XP[2]-X[2];
      double a = bond[0];
      double b = bond[1];
      double c = bond[2];
      double xi2 = a*a+b*b+c*c;
      double omega=(*oc)*OMEGA_0(std::sqrt(xi2),horizon);
      *m += omega*(1.0-*bondDamage)*xi2*cellVolume;
    }
  }
}
/*
 * Computes a weighted volume but with the pals deviatoric influence function
 */
void computeWeightedVolume
(
    const double *xOverlap,
    const double *volumeOverlap,
    const vector<const double *>& sigma_multipliers,
    const double *sigma_constant,
    double *weighted_volume,
    int numOwnedPoints,
    const int* localNeighborList,
    double horizon,
    const FunctionPointer SIGMA_0
)
{
  double bond[3];
  const double *xOwned = xOverlap;
  double tau_X[NUM_LAGRANGE_MULTIPLIERS];
  const double *sc=sigma_constant;
  double *m=weighted_volume;
  const int *neighPtr = localNeighborList;
  double cellVolume;
  for(int q=0; q<numOwnedPoints;q++, xOwned+=3, sc++, m++){
    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    *m=0.0;

    // Collect computed Lagrange multipliers for this point
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      tau_X[i]=sigma_multipliers[i][q];
    }

    for(int n=0;n<numNeigh;n++,neighPtr++){
      int localId = *neighPtr;
      cellVolume = volumeOverlap[localId];
      const double *XP = &xOverlap[3*localId];
      bond[0]=XP[0]-X[0];
      bond[1]=XP[1]-X[1];
      bond[2]=XP[2]-X[2];
      double a = bond[0];
      double b = bond[1];
      double c = bond[2];
      double xi2 = a*a+b*b+c*c;
      pals_influence<deviatoric_influence> SIGMA(SIGMA_0,*sc,tau_X);
      double sigma = SIGMA(bond,horizon);
      *m += sigma*xi2*cellVolume;
    }
  }
}

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
){

  double bond[3];
  const double *xOwned = xOverlap;
  const double *yOwned = yOverlap;
  double lambda_X[NUM_LAGRANGE_MULTIPLIERS];
  const double *oc=omega_constant;
  double *theta = dilatation;
  double cellVolume;
  const int *neighPtr = localNeighborList;
  for(int q=0; q<numOwnedPoints;q++, xOwned+=3, yOwned+=3, oc++, theta++){
    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const double *Y = yOwned;
    // Collect computed Lagrange multipliers for this point
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      lambda_X[i]=_omega_multipliers[i][q];
    }
    *theta = double(0.0);
    for(int n=0;n<numNeigh;n++,neighPtr++){
      int localId = *neighPtr;
      cellVolume = volumeOverlap[localId];
      const double *XP = &xOverlap[3*localId];
      const double *YP = &yOverlap[3*localId];
      bond[0]=XP[0]-X[0];
      bond[1]=XP[1]-X[1];
      bond[2]=XP[2]-X[2];
      double a = bond[0];
      double b = bond[1];
      double c = bond[2];
      double xi2 = a*a+b*b+c*c;
      a = YP[0]-Y[0];
      b = YP[1]-Y[1];
      c = YP[2]-Y[2];
      double dY = a*a+b*b+c*c;
      double x = sqrt(xi2);
      double e = sqrt(dY)-x;
      pals_influence<dilatation_influence> OMEGA(OMEGA_0,*oc,lambda_X);
      double omega = OMEGA(bond,horizon);
      double omega_x=omega*x;
      *theta+=omega_x*e*cellVolume;
    }
  }
}

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
){

  double K = BULK_MODULUS;
  double TWO_MU = 2.0 * SHEAR_MODULUS;
  double bond[3];
  const double *xOwned = xOverlap;
  const double *yOwned = yOverlap;
  double lambda_X[NUM_LAGRANGE_MULTIPLIERS];
  const double *oc=omega_constant;
  double tau_X[NUM_LAGRANGE_MULTIPLIERS];
  const double *sc=sigma_constant;
  const double *m=weighted_volume;
  double *theta = dilatation;
  double *p = pals_pressure;
  double cellVolume;
  const int *neighPtr = localNeighborList;
  for(int q=0; q<numOwnedPoints;q++, xOwned+=3, yOwned+=3, oc++, sc++, m++, theta++, p++){
    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const double *Y = yOwned;
    // Collect computed Lagrange multipliers for this point
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      lambda_X[i]=_omega_multipliers[i][q];
      tau_X[i]=_sigma_multipliers[i][q];
    }
    *theta = double(0.0);
    *p = double(0.0);
    for(int n=0;n<numNeigh;n++,neighPtr++){
      int localId = *neighPtr;
      cellVolume = volumeOverlap[localId];
      const double *XP = &xOverlap[3*localId];
      const double *YP = &yOverlap[3*localId];
      bond[0]=XP[0]-X[0];
      bond[1]=XP[1]-X[1];
      bond[2]=XP[2]-X[2];
      double a = bond[0];
      double b = bond[1];
      double c = bond[2];
      double xi2 = a*a+b*b+c*c;
      a = YP[0]-Y[0];
      b = YP[1]-Y[1];
      c = YP[2]-Y[2];
      double dY = a*a+b*b+c*c;
      double x = sqrt(xi2);
      double e = sqrt(dY)-x;
      pals_influence<dilatation_influence> OMEGA(OMEGA_0,*oc,lambda_X);
      pals_influence<deviatoric_influence> SIGMA(SIGMA_0,*sc,tau_X);
      double omega = OMEGA(bond,horizon);
      double sigma = SIGMA(bond,horizon);
      double omega_x=omega*x;
      double sigma_x=sigma*x;
      *theta+=omega_x*e*cellVolume;
      *p+=-(TWO_MU*sigma_x/3.0)*e*cellVolume;
    }
    // Final piece of pals_pressure requires dilatation that is only ready here
    *p+=(K+TWO_MU*(*m)/9.0)*(*theta);
  }
}

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
)
{
  /*
   * Compute processor local contribution to internal force
   */
  //double K = BULK_MODULUS;
  double TWO_MU = 2.0 * SHEAR_MODULUS;

  const double *xOwned = xOverlap;
  const double *yOwned = yOverlap;
  double lambda_X[NUM_LAGRANGE_MULTIPLIERS];
  const double *oc=omega_constant;
  double tau_X[NUM_LAGRANGE_MULTIPLIERS];
  const double *sc=sigma_constant;
  double *fOwned = fInternalOverlap;

  double bond[3];
  double a, b, c;
  double xi;
  double Y_dx, Y_dy, Y_dz;
  double dY;
  double e;
  double eps;
  double omega, sigma;
  double t;
  double fx, fy, fz;
  double cell_volume;
  const double *theta = dilatation;
  const double *p = pals_pressure;
  const int *neighPtr = localNeighborList;
  for(int q=0; q<numOwnedPoints;q++, xOwned+=3, yOwned+=3, fOwned+=3, oc++, sc++, theta++, p++){

    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const double *Y = yOwned;
    // Collect computed Lagrange multipliers for this point
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      lambda_X[i]=_omega_multipliers[i][q];
      tau_X[i]=_sigma_multipliers[i][q];
    }

    double self_cell_volume = volumeOverlap[q];
    for(int n=0;n<numNeigh;n++,neighPtr++){
      int localId = *neighPtr;
      const double *XP = &xOverlap[3*localId];
      const double *YP = &yOverlap[3*localId];
      bond[0] = XP[0]-X[0];
      bond[1] = XP[1]-X[1];
      bond[2] = XP[2]-X[2];
      a = bond[0];
      b = bond[1];
      c = bond[2];
      xi = sqrt(a*a+b*b+c*c);
      Y_dx = YP[0]-Y[0];
      Y_dy = YP[1]-Y[1];
      Y_dz = YP[2]-Y[2];
      dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
      e = dY - xi;
      eps=e-(*theta)*xi/3.0;
      pals_influence<dilatation_influence> OMEGA(OMEGA_0,*oc,lambda_X);
      pals_influence<deviatoric_influence> SIGMA(SIGMA_0,*sc,tau_X);
      omega = OMEGA(bond,horizon);
      sigma = SIGMA(bond,horizon);


      t=(*p)*omega*xi+TWO_MU*sigma*eps;
      fx = t * Y_dx / dY;
      fy = t * Y_dy / dY;
      fz = t * Y_dz / dY;
      cell_volume= volumeOverlap[localId];
      *(fOwned+0) += fx*cell_volume;
      *(fOwned+1) += fy*cell_volume;
      *(fOwned+2) += fz*cell_volume;
      fInternalOverlap[3*localId+0] -= fx*self_cell_volume;
      fInternalOverlap[3*localId+1] -= fy*self_cell_volume;
      fInternalOverlap[3*localId+2] -= fz*self_cell_volume;
    }

  }
}
}
}
