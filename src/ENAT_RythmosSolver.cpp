#include "ENAT_RythmosSolver.hpp"
#include <Rythmos_ForwardEulerStepper.hpp>
//#include <Rythmos_ExplicitRKStepper.hpp>
#include <Rythmos_SimpleIntegrationControlStrategy.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_EpetraModelEvaluator.hpp>

ENAT::RythmosSolver::RythmosSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
                          Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
                          Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer ) :
  appParams(appParams_),
  model(model_)//,
//   p_init(model_->get_p_init(0)),
//   g_map(model_->get_g_map(0))
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::rcp(&(appParams_->sublist("Solver")),false);

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  bool verbose = solverParams->get("Verbose", false);

  if(verbose)
	*out << "\nA) Get the base parameter list ...\n";

  RCP<Teuchos::ParameterList> rythmosPL = sublist(solverParams, "Rythmos", true);

  if(verbose){
	std::cout << "Rythmos parameter list:" << std::endl;
	rythmosPL->print(std::cout, 2, true, true) ;
  }

  // options:  VERB_DEFAULT, VERB_NONE, VERB_LOW, VERB_MEDIUM, VERB_HIGH, VERB_EXTREME
  solnVerbLevel = Teuchos::VERB_HIGH;

//   const int numTimeSteps = rythmosPL->get("Num Time Steps", 10);
  const Scalar t_init = 0.0;
  t_final = rythmosPL->get("Final Time", 0.1);
      
  const Rythmos::TimeRange<Scalar> fwdTimeRange(t_init, t_final); // THIS IS NEVER USED??

  //
  // This is the linear solve strategy that will be used to solve for the
  // linear system with the W.
  //
  if(verbose)
	*out << "\nB) Create the Stratimikos linear solver factory ...\n";
      
  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  linearSolverBuilder.setParameterList(sublist(rythmosPL, "Stratimikos", true));
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory = createLinearSolveStrategy(linearSolverBuilder);
      
  if(verbose)
	*out << "\nC) Create and initalize the forward model ...\n";
      
  // C.1) Create the underlying EpetraExt::ModelEvaluator
      
  // already constructed as "model"
      
  // C.2) Create the Thyra-wrapped ModelEvaluator
      
  fwdStateModel = epetraModelEvaluator(model, W_factory);
      
  if(verbose)
	*out << "\nD) Create the stepper and integrator for the forward problem ...\n";

  // Forward Euler integrator
  fwdStateStepper = Rythmos::forwardEulerStepper<double>(fwdStateModel);

  fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));
  {
    RCP<Teuchos::ParameterList>
      integrationControlPL = sublist(rythmosPL, "Rythmos Integration Control", true);
    
    // RCP<Rythmos::IntegratorBase<Scalar> >
    RCP<Rythmos::DefaultIntegrator<Scalar> >
      defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(Rythmos::simpleIntegrationControlStrategy<Scalar>(integrationControlPL));
    fwdStateIntegrator = defaultIntegrator;
	if(verbose)
	  cout << "integrationControlPL:\n" << *integrationControlPL << endl ;
  }
  fwdStateIntegrator->setParameterList(sublist(rythmosPL, "Rythmos Integrator", true));

  if (observer != Teuchos::null) 
    fwdStateIntegrator->setIntegrationObserver(observer);
  
}

ENAT::RythmosSolver::~RythmosSolver()
{
}

Teuchos::RCP<const Epetra_Map> ENAT::RythmosSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> ENAT::RythmosSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

// Teuchos::RCP<const Epetra_Map> ENAT::RythmosSolver::get_p_map(int l) const
// {
//   TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
//                      std::endl <<
//                      "Error!  App::ModelEval::get_p_map() only " <<
//                      " supports 1 parameter vector.  Supplied index l = " <<
//                      l << std::endl);

//   return model->get_p_map(l);
// }

// Teuchos::RCP<const Epetra_Map> ENAT::RythmosSolver::get_g_map(int j) const
// {
//   TEST_FOR_EXCEPTION(j != 0, Teuchos::Exceptions::InvalidParameter,
//                      std::endl <<
//                      "Error!  ENAT::RythmosSolver::get_g_map() only " <<
//                      " supports 1 parameter vector.  Supplied index l = " <<
//                      j << std::endl);

//   return g_map;
// }

// Teuchos::RCP<const Epetra_Vector> ENAT::RythmosSolver::get_x_init() const
// {
//   Teuchos::RCP<const Epetra_Vector> neverused;
//   return neverused;
// }

// Teuchos::RCP<const Epetra_Vector> ENAT::RythmosSolver::get_p_init(int l) const
// {
//   TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
//                      std::endl <<
//                      "Error!  App::ModelEval::get_p_map() only " <<
//                      " supports 1 parameter vector.  Supplied index l = " <<
//                      l << std::endl);

//   return p_init;
// }

EpetraExt::ModelEvaluator::InArgs ENAT::RythmosSolver::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
//   inArgs.set_Np(1);
//  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs ENAT::RythmosSolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
//   outArgs.set_Np_Ng(1, 1);
//  outArgs.setSupports(OUT_ARG_DgDp, 0, 0, DerivativeSupport(DERIV_MV_BY_COL));
  return outArgs;
}

void ENAT::RythmosSolver::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool verbose = appParams->sublist("Problem").get("Verbose", false);

  // Parse InArgs
//   RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
//   if (!p_in.get()) cout << "ERROR: ENAT::RythmosSolver requires p as inargs" << std::endl;

//   // Parse OutArgs

//   RCP<Epetra_Vector> g_out = outArgs.get_g(0); 

//   // Parse out-args for sensitivity calculation
//   RCP<Epetra_MultiVector> dgdp_out;
//   RCP<Epetra_MultiVector> dxdp;

  RCP<const Epetra_Vector> finalSolution;

  if(verbose)
 	*out << "\nE) Solve the forward problem ...\n";

  // the call to getNominalValues() returns initial_conditions which is filled using
  // model.get_x_init(), model.get_x_dot_init(), etc., for various data that the model supports
  // also makes calls to model.get_x_lower_bounds(), model.get_x_upper_bounds(), etc.
  Thyra::ModelEvaluatorBase::InArgs<Scalar> initial_conditions = fwdStateModel->getNominalValues();

//HERE  *out << "Initial x:\n" << Teuchos::describe(*initial_conditions.get_x(), solnVerbLevel);

//  initial_conditions.set_p(0,Thyra::create_Vector(p_in, fwdStateModel->get_p_space(0)));

//   // THIS IS NOT WORKING FOR EXPLICIT RK STEPPER
  fwdStateStepper->setInitialCondition(initial_conditions);

//   //fwdStateStepper->setModel(fwdStateModel); // HERE, redundant, also blows away settings?  never get to Peridigm::ModelEvaluator::evalModel() if this is set
  fwdStateIntegrator->setStepper(fwdStateStepper, t_final, true);

  // Integration occurs here, with the call to getFwdPoints()
  
  Teuchos::Array<RCP<const Thyra::VectorBase<Scalar> > > x_final_array;
  fwdStateIntegrator->getFwdPoints(Teuchos::tuple<Scalar>(t_final), &x_final_array, NULL, NULL);
  const RCP<const Thyra::VectorBase<Scalar> > x_final = x_final_array[0];

//   *out << "\nFinal x:\n" << Teuchos::describe(*x_final, solnVerbLevel);
    
//   finalSolution = Thyra::get_Epetra_Vector(*model->get_x_map(), x_final);

//   cout << "Final Solution\n" << *finalSolution << std::endl;

//   // As post-processing step, calc responses at final solution
//   if (g_out != Teuchos::null) {
//     EpetraExt::ModelEvaluator::InArgs model_inargs = model->createInArgs();
//     EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
//     model_inargs.set_x(finalSolution);
//     if(model_inargs.Np() > 0)
//       model_inargs.set_p(0, p_in);
//     g_out->PutScalar(0.0);
//     model_outargs.set_g(0, g_out);
//     model->evalModel(model_inargs, model_outargs);
//   }

}
