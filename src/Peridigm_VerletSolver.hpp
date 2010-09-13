#ifndef PERIDIGM_VERLETSOLVER_H
#define PERIDIGM_VERLETSOLVER_H

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <EpetraExt_ModelEvaluator.h>
#include <Teuchos_ParameterList.hpp>
#include <Peridigm_VerletIntegrator.hpp>
#include <Peridigm_VerletObserver.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include <Epetra_SerialComm.h>
#endif

namespace PeridigmNS {

class VerletSolver: public EpetraExt::ModelEvaluator {

  public:

  VerletSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
               Teuchos::RCP<EpetraExt::ModelEvaluator> model,
               Teuchos::RCP<PeridigmNS::VerletObserver > observer = Teuchos::null
               );

  ~VerletSolver();

  EpetraExt::ModelEvaluator::InArgs createInArgs() const;

  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;

  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  private:

  Teuchos::RCP<const Epetra_Map> get_x_map() const;

  Teuchos::RCP<const Epetra_Map> get_f_map() const;

  private:

  //These are set in the constructor and used in evalModel
  mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
  Teuchos::RCP<EpetraExt::ModelEvaluator> model;

  Teuchos::RCP<Teuchos::FancyOStream> out;
  Teuchos::EVerbosityLevel solnVerbLevel;
  Teuchos::RCP<PeridigmNS::VerletIntegrator> integrator;
  double t_final;

};
}
#endif
