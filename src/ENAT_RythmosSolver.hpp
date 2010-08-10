#ifndef ENAT_RYTHMOSSOLVER_H
#define ENAT_RYTHMOSSOLVER_H

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include <Epetra_SerialComm.h>
#endif

#include <EpetraExt_ModelEvaluator.h>
#include <Rythmos_DefaultIntegrator.hpp>
#include <Rythmos_IntegrationObserverBase.hpp>

/** 
 * \brief This set of classes will support a wide number of different types of abstract
 * problem types that will allow NOX, LOCA, Rythmos, Aristos, and MOOCHO to
 * solve different types of problems with Charon.
 * 
 * \todo Finish documentation for ENAT.
 */

namespace ENAT {

class RythmosSolver
    : public EpetraExt::ModelEvaluator
{

  typedef double Scalar;

  public:

  /*! \name Constructors/initializers */
  //@{

  /*! \brief Takes the number of elements in the discretization . */
  RythmosSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
                Teuchos::RCP<EpetraExt::ModelEvaluator> model,
                Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer = Teuchos::null
                );

  //@}

  ~RythmosSolver();

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /*! \brief Response function, for example the objective function in
   * an optimization problem.  Not required for the case of simple time
   * integration.
   */
//   Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;

  /*! \brief Parameter space, treated as independent variables for 
   *  optimization and sensitivity analysis.
   */
//   Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

  /*! \brief Jacobian, written in the form
   *  \f$ \mathbf{W} = \alpha \frac{\partial \mathbf{f}}{\partial \dot{x}}
   *  + \beta \frac{\partial \mathbf{f}}{\partial x} \f$. */
  //Teuchos::RCP<Epetra_Operator> create_W() const;

  /*! \brief . */
   EpetraExt::ModelEvaluator::InArgs createInArgs() const;

  /*! \brief . */
  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;

  /*! \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  private:

  /*! \brief Solution space. */
   Teuchos::RCP<const Epetra_Map> get_x_map() const;

  /*! \brief Governing equation.  For implicit time integration,
   *  \f$ f \f$ is defined by moving all terms to one side of the equation
   * (setting it equal to zero).  For explicit time integration, 
   * \f$ f = \frac{\partial \mathbf{f}}{\partial x} \f$.
   *
   */
   Teuchos::RCP<const Epetra_Map> get_f_map() const;

  /*! \brief Solution space initial values. */
//   Teuchos::RCP<const Epetra_Vector> get_x_init() const;

  /*! \brief Parameter space. */
//   Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

  /*! \brief Sets default values for problem parameters. */
//   void setProblemParamDefaults(Teuchos::ParameterList* appParams_);

  /** \brief Sets default values for solver parameters. */
//   void setSolverParamDefaults(Teuchos::ParameterList* appParams_);

  //@}

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;

//    Teuchos::RCP<const Epetra_Vector> p_init;
//    Teuchos::RCP<const Epetra_Map> g_map;

    Teuchos::RCP<Rythmos::StepperBase<Scalar> > fwdStateStepper;
    Teuchos::RCP<Teuchos::FancyOStream> out;
    Teuchos::EVerbosityLevel solnVerbLevel;
    Teuchos::RCP<Rythmos::DefaultIntegrator<Scalar> > fwdStateIntegrator;
    Teuchos::RCP<Thyra::ModelEvaluator<double> > fwdStateModel;
    Scalar t_final;

};
}
#endif
