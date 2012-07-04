#include "Peridigm_VerletSolver.hpp"
#include <Thyra_EpetraModelEvaluator.hpp>

PeridigmNS::VerletSolver::VerletSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
                                     Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
                                     Teuchos::RCP<PeridigmNS::VerletObserver> observer ) :
  appParams(appParams_),
  model(model_) {

  Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::rcp(&(appParams_->sublist("Solver")),false);

  bool verbose = solverParams->get("Verbose", false);

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if(verbose)
    *out << "\nA) Get the base parameter list ...\n";

  Teuchos::RCP<Teuchos::ParameterList> verletPL = sublist(solverParams, "Verlet", true);

  if(verbose){
        std::cout << "Verlet parameter list:" << std::endl;
        verletPL->print(std::cout, 2, true, true) ;
  }

  // options:  VERB_DEFAULT, VERB_NONE, VERB_LOW, VERB_MEDIUM, VERB_HIGH, VERB_EXTREME
  solnVerbLevel = Teuchos::VERB_HIGH;

  // should get this from input deck
  t_final   = verletPL->get("Final Time", 1.0);

  if(verbose)
        *out << "\nB) Create the integrator for the forward problem ...\n";

  integrator = Teuchos::rcp(new PeridigmNS::VerletIntegrator(model,verletPL));

  if (observer != Teuchos::null)
    integrator->setIntegrationObserver(observer);

}

PeridigmNS::VerletSolver::~VerletSolver() {}

Teuchos::RCP<const Epetra_Map> PeridigmNS::VerletSolver::get_x_map() const {
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> PeridigmNS::VerletSolver::get_f_map() const {
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

EpetraExt::ModelEvaluator::InArgs PeridigmNS::VerletSolver::createInArgs() const {
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
//   inArgs.set_Np(1);
//  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs PeridigmNS::VerletSolver::createOutArgs() const {
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
//   outArgs.set_Np_Ng(1, 1);
//  outArgs.setSupports(OUT_ARG_DgDp, 0, 0, DerivativeSupport(DERIV_MV_BY_COL));
  return outArgs;
}

void PeridigmNS::VerletSolver::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const{

  bool verbose = appParams->get("Verbose", false);

  // Parse InArgs

  // Parse OutArgs

  Teuchos::RCP<const Epetra_Vector> finalSolution;

  if(verbose)
        *out << "\nC) Solve the forward problem ...\n";

  // the call to getNominalValues() returns initial_conditions which is filled using
  // model.get_x_init(), model.get_x_dot_init(), etc., for various data that the model supports
  // also makes calls to model.get_x_lower_bounds(), model.get_x_upper_bounds(), etc.
//  Thyra::ModelEvaluatorBase::InArgs<Scalar> initial_conditions = fwdStateModel->getNominalValues();

//HERE  *out << "Initial x:\n" << Teuchos::describe(*initial_conditions.get_x(), solnVerbLevel);

//  initial_conditions.set_p(0,Thyra::create_Vector(p_in, fwdStateModel->get_p_space(0)));

  // Integration occurs here, with the call to getFwdPoints()

/*
  Teuchos::Array<RCP<const Thyra::VectorBase<Scalar> > > x_final_array;
  fwdStateIntegrator->getFwdPoints(Teuchos::tuple<Scalar>(t_final), &x_final_array, NULL, NULL);
  const RCP<const Thyra::VectorBase<Scalar> > x_final = x_final_array[0];
*/

  integrator->integrate(t_final);

//   *out << "\nFinal x:\n" << Teuchos::describe(*x_final, solnVerbLevel);

//   finalSolution = Thyra::get_Epetra_Vector(*model->get_x_map(), x_final);

//   cout << "Final Solution\n" << *finalSolution << std::endl;

}
