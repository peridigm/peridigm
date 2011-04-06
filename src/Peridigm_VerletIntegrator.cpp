/*! \file Peridigm_VerletIntegrator.cpp */
//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include "Peridigm_VerletIntegrator.hpp"

PeridigmNS::VerletIntegrator::VerletIntegrator(Teuchos::RCP<EpetraExt::ModelEvaluator> model_, Teuchos::RCP<Teuchos::ParameterList> PL) {
  model = model_;
  solverParams = PL;
  t_initial = PL->get("Initial Time", 1.0); 
  t_current = t_initial;
  dt        = PL->get("Fixed dt", 1.0); 
}

int PeridigmNS::VerletIntegrator::setIntegrationObserver(Teuchos::RCP<PeridigmNS::VerletObserver> observer_) {observer = observer_; return 0;}

int PeridigmNS::VerletIntegrator::integrate(double t_final) {

  // Get initial condition from model 
  Teuchos::RCP<Epetra_Vector> x_in = Teuchos::rcp(new Epetra_Vector(*(model->get_x_init())));
  Teuchos::RCP<Epetra_Vector> x_out = Teuchos::rcp(new Epetra_Vector(x_in->Map()));

  EpetraExt::ModelEvaluator::InArgs   inArgs = model->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs outArgs = model->createOutArgs();

  inArgs.set_x(x_in);
  outArgs.set_f(x_out);

  // Get one-dimensional map from the model evaluator; This gives GID info for each node, and # nodes on my processor
  Teuchos::RCP<PeridigmNS::ModelEvaluator> PDmodel = Teuchos::rcp_dynamic_cast<PeridigmNS::ModelEvaluator>(model);
  Teuchos::RCP<const Epetra_Map> dataMap = PDmodel->getOneDimensionalMap();
  int length = dataMap->NumMyElements();

  // Initial call to model evaluator to construct x_out (needed by velocity Verlet)
  model->evalModel(inArgs,outArgs);

  // Use BLAS for local-only vector updates (BLAS-1)
  Epetra_BLAS blas;

  // Pointer index into sub-vectors for use with BLAS
  double *x0,*x1,*x2;
  x_in->ExtractView( &x0 );
  x1 = &(x0[1]);
  x2 = &(x0[2]);
  double *v0,*v1,*v2;
  v0 = &(x0[3]);
  v1 = &(x0[4]);
  v2 = &(x0[5]);
  double *a0,*a1,*a2;
  x_out->ExtractView( &a2 );
  a0 = &(a2[3]);
  a1 = &(a2[4]);
  a2 = &(a2[5]);

  // Integrate from t_current up to t_final with timestep dt
  double dt2 = dt/2.0;
  double tmp = t_current;
  int step = 0;
  while (tmp < t_final) {

    // Do one step of velocity-Verlet
  
    // V^{n+1/2} = V^{n} + (dt/2)*A^{n}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const  
    blas.AXPY(length, dt2, a0, v0, 6, 6);
    blas.AXPY(length, dt2, a1, v1, 6, 6);
    blas.AXPY(length, dt2, a2, v2, 6, 6);

    // X^{n+1}   = X^{n} + (dt)*V^{n+1/2}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const  
    blas.AXPY(length, dt, v0, x0, 6, 6);
    blas.AXPY(length, dt, v1, x1, 6, 6);
    blas.AXPY(length, dt, v2, x2, 6, 6);

    // Update forces based on new positions
    model->evalModel(inArgs,outArgs);

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const  
    blas.AXPY(length, dt2, a0, v0, 6, 6);
    blas.AXPY(length, dt2, a1, v1, 6, 6);
    blas.AXPY(length, dt2, a2, v2, 6, 6);

    step = step + 1;
    tmp = t_current + step*dt;

    // Call the observer
    observer->observeCompletedTimeStep(x_in,tmp);

  }

  t_current = tmp;

  // Return success
  return 0;

}
