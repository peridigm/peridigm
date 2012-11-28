/*! \file Peridigm.hpp */

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

#ifndef PERIDIGM_HPP
#define PERIDIGM_HPP

#include <vector>
#include <set>

#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_FECrsMatrix.h>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>

#include "Peridigm_Block.hpp"
#include "Peridigm_AbstractDiscretization.hpp"
#include "Peridigm_ModelEvaluator.hpp"
#include "Peridigm_DataManager.hpp"
#include "Peridigm_SerialMatrix.hpp"
#include "Peridigm_OutputManagerContainer.hpp"
#include "Peridigm_ComputeManager.hpp"
#include "Peridigm_BoundaryAndInitialConditionManager.hpp"
#include "materials/Peridigm_Material.hpp"
#include "damage/Peridigm_DamageModel.hpp"
#include "contact/Peridigm_ContactModel.hpp"
#include "muParser/muParserDef.h"

#include "mesh_input/quick_grid/QuickGridData.h"

namespace PeridigmNS {

  class Peridigm : public NOX::Epetra::Interface::Required, public NOX::Epetra::Interface::Jacobian, public NOX::Epetra::Interface::Preconditioner {

  public:

    //! Constructor
    Peridigm(const Teuchos::RCP<const Epetra_Comm>& comm,
             const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    ~Peridigm(){};

    //! Instantiate material models
    void instantiateMaterials();

    //! Instantiate damage models
    void instantiateDamageModels();

    //! Initialize discretization and maps
    void initializeDiscretization(Teuchos::RCP<AbstractDiscretization> peridigmDisc);

    //! Load node sets from input deck and/or input mesh into nodeSets container
    void initializeNodeSets(Teuchos::RCP<Teuchos::ParameterList>& bcParams,
                            Teuchos::RCP<AbstractDiscretization> peridigmDisc);

    //! Initialize contact
    void initializeContact();

    //! Initialize the workset
    void initializeWorkset();

    //! Instantiate the compute manager
    void instantiateComputeManager();

    //! Initialize the output manager
    void initializeOutputManager();

    //! Initialize blocks
    void initializeBlocks(Teuchos::RCP<AbstractDiscretization> disc);

    //! Main routine to drive time integration
    void execute();

    //! Compute the residual vector (pure virtual method in NOX::Epetra::Interface::Required).
    bool computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType fillType = Residual);

    //! Compute the Jacobian (pure virtual method in NOX::Epetra::Interface::Jacobian).
    bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

    //! Number of nonlinear iterations between Jacobian updates for NOX Nonlinear CG solves.
    int m_noxTriggerJacobianUpdate;

    //! Counter for updating the Jacobian for NOX Nonlinear CG solves.
    int m_noxJacobianUpdateCounter;

    //! Compute the preconditioner (pure virtual method in NOX::Epetra::Interface::Preconditioner).
    bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& Prec, Teuchos::ParameterList* precParams = 0);

    //Return a reference to the Epetra_CrsMatrix withh the Jacobian that is
    //generated by the interface
    Teuchos::RCP<Epetra_CrsMatrix> getJacobian();
    
    //! Residual and Jacobian matrix fills for NOX interface
    virtual bool evaluateNOX(FillType f, const Epetra_Vector *solnVector, Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix);

    //! Perform diagnostics on Jacobian and print results to screen.
    void jacobianDiagnostics(Teuchos::RCP<NOX::Epetra::Group> noxGroup);

    void executeExplicit();

    //! Main routine to drive problem solution for quasistatics
    void executeQuasiStatic();

    //! Main routine to drive problem solution for quasistatics using NOX
    void executeNOXQuasiStatic(); 

    //! Set the preconditioner for the global linear system
    void quasiStaticsSetPreconditioner(Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>& linearProblem);

    //! Damp the tangent matrix by scaling the diagonal and adding a small value to each entry in the diagonal
    void quasiStaticsDampTangent(double dampedNewtonDiagonalScaleFactor,
                                 double dampedNewtonDiagonalShiftFactor);

    //! Solve the global linear system
    Belos::ReturnType quasiStaticsSolveSystem(Teuchos::RCP<Epetra_Vector> residual,
					      Teuchos::RCP<Epetra_Vector> lhs,
					      Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>& linearProblem,
					      Teuchos::RCP< Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> >& belosSolver);

    //! Perform line search
    double quasiStaticsLineSearch(Teuchos::RCP<Epetra_Vector> residual,
				  Teuchos::RCP<Epetra_Vector> lhs);

    //! Main routine to drive problem solution with implicit time integration
    void executeImplicit();

    //! Allocate memory for non-zeros in global Jacobian
    void allocateJacobian();

    //! Compute the Jacobian for implicit dynamics
    void computeImplicitJacobian(double beta);

    //! Compute the residual for quasi-statics
    double computeQuasiStaticResidual(Teuchos::RCP<Epetra_Vector> residual);

    //! Synchronize data in DataManagers across processes (needed before call to OutputManager::write() )
    void synchDataManagers();

    //! Rebalance the mesh
    void rebalance();

    //! Compute a parallel decomposion based on the current configuration
    QUICKGRID::Data currentConfigurationDecomp();

    //! Create a rebalanced bond map
    Teuchos::RCP<Epetra_BlockMap> createRebalancedBondMap(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                          Teuchos::RCP<const Epetra_Import> oneDimensionalMapToRebalancedOneDimensionalMapImporter);

    //! Create a global ID neighbor list in a rebalanced partitioning
    Teuchos::RCP<Epetra_Vector> createRebalancedNeighborGlobalIDList(Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap,
                                                                     Teuchos::RCP<const Epetra_Import> bondMapToRebalancedBondMapImporter);

    //! Create a rebalanced NeighborhoodData object
    Teuchos::RCP<PeridigmNS::NeighborhoodData> createRebalancedNeighborhoodData(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                                Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalOverlapMap,
                                                                                Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap,
                                                                                Teuchos::RCP<Epetra_Vector> rebalancedNeighborGlobalIDs);

    //! Fill the contact neighbor information in rebalancedDecomp and populate contactNeighborsGlobalIDs and offProcesorContactIDs
    void contactSearch(Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap,
                       Teuchos::RCP<const Epetra_BlockMap> rebalancedBondMap,
                       Teuchos::RCP<const Epetra_Vector> rebalancedNeighborGlobalIDs,
                       QUICKGRID::Data& rebalancedDecomp,
                       Teuchos::RCP< std::map<int, std::vector<int> > > contactNeighborGlobalIDs,
                       Teuchos::RCP< std::set<int> > offProcessorContactIDs);

    //! Create a rebalanced NeighborhoodData object for contact
    Teuchos::RCP<PeridigmNS::NeighborhoodData> createRebalancedContactNeighborhoodData(Teuchos::RCP<std::map<int, std::vector<int> > > contactNeighborGlobalIDs,
                                                                                       Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                                       Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalOverlapMap);

    //! Accessor for comm object
    Teuchos::RCP<const Epetra_Comm> getEpetraComm(){ return peridigmComm; }

    //! @name Accessors for maps 
    //@{ 
    Teuchos::RCP<const Epetra_BlockMap> getOneDimensionalMap() { return oneDimensionalMap; }
    Teuchos::RCP<const Epetra_BlockMap> getThreeDimensionalMap() { return threeDimensionalMap; }
    Teuchos::RCP<const Epetra_BlockMap> getBondMap() { return bondMap; }
    Teuchos::RCP<const Epetra_BlockMap> getOneDimensionalOverlapMap() { return oneDimensionalOverlapMap; }
    //@}

    //! @name Accessors for main solver-level vectors
    //@{ 
    Teuchos::RCP<const Epetra_Vector> getBlockIDs() { return blockIDs; }
    Teuchos::RCP<const Epetra_Vector> getX() { return x; }
    Teuchos::RCP<const Epetra_Vector> getU() { return u; }
    Teuchos::RCP<const Epetra_Vector> getY() { return y; }
    Teuchos::RCP<const Epetra_Vector> getV() { return v; }
    Teuchos::RCP<const Epetra_Vector> getA() { return a; }
    Teuchos::RCP<const Epetra_Vector> getForce() { return force; }
    Teuchos::RCP<const Epetra_Vector> getContactForce() { return contactForce; }
    Teuchos::RCP<const Epetra_Vector> getDeltaU() { return deltaU; }
    //@}

    //! @name Accessors for neighborhood data
    //@{
    Teuchos::RCP<const PeridigmNS::NeighborhoodData> getGlobalNeighborhoodData() { return globalNeighborhoodData; }
    Teuchos::RCP<const PeridigmNS::NeighborhoodData> getGlobalContactNeighborhoodData() { return globalContactNeighborhoodData; }
    //@}

    //! Accessor for vector of Blocks
    Teuchos::RCP< std::vector<PeridigmNS::Block> > getBlocks() { return blocks; }

    //! Accessor for individual Blocks
    Teuchos::RCP<PeridigmNS::Block> getBlock(unsigned int blockNumber) {
      return Teuchos::rcpFromRef(blocks->at(blockNumber));
    }

    //! Accessor for Material Models
    std::map< std::string, Teuchos::RCP<const PeridigmNS::Material> > getMaterialModels() {
      return materialModels;
    }

    //! Accessor for compute manager
    Teuchos::RCP< PeridigmNS::ComputeManager > getComputeManager() { return computeManager; }

    //! Display a progress bar
    void displayProgress(string title, double percentComplete);

  private:

    //! Friend classes
    friend class OutputManager_ExodusII;

    //! Parameterlist of entire input deck
    Teuchos::RCP<Teuchos::ParameterList> peridigmParams;

    //! Epetra communicator established by Peridigm_Factory
    Teuchos::RCP<const Epetra_Comm> peridigmComm;

    //! Output stream
    Teuchos::RCP<Teuchos::FancyOStream> out;

    //! Maps for scalar, vector, and bond data
    Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<const Epetra_BlockMap> threeDimensionalMap;
    Teuchos::RCP<const Epetra_BlockMap> bondMap;
    Teuchos::RCP<const Epetra_BlockMap> oneDimensionalOverlapMap;

    //! Current simulation time
    double timeCurrent;

    //! Number of blocks in the model
    int numBlocks;

    //! Rebalance flag
    bool analysisHasRebalance;

    //! Rebalance frequency
    int rebalanceFrequency;

    //! Contact flag
    bool analysisHasContact;

    //! Contact search frequency
    int contactRebalanceFrequency;

    //! Contact search radius
    double contactSearchRadius;

    //! Partial volume flag
    bool analysisHasPartialVolumes;

    //! Material models
    std::map< std::string, Teuchos::RCP<const PeridigmNS::Material> > materialModels;

    //! Damage models
    std::map< std::string, Teuchos::RCP<const PeridigmNS::DamageModel> > damageModels;

    //! Contact models
    std::map< std::string, Teuchos::RCP<const PeridigmNS::ContactModel> > contactModels;

    //! Boundary and initial condition manager
    Teuchos::RCP<PeridigmNS::BoundaryAndInitialConditionManager> boundaryAndInitialConditionManager;

    //! Compute manager
    Teuchos::RCP<PeridigmNS::ComputeManager> computeManager;

    //! Mothership multivector that contains all the three-dimensional global vectors (x, u, y, v, a, force, etc.)
    Teuchos::RCP<Epetra_MultiVector> threeDimensionalMothership;

    //! Mothership multivector that contains all the one-dimensional global vectors (blockID, volume)
    Teuchos::RCP<Epetra_MultiVector> oneDimensionalMothership;

    //! Blocks
    Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks;

    //! Block iterator, for convenience
    std::vector<PeridigmNS::Block>::iterator blockIt;

    //! Node sets
    Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets;

    //! Overlap Jacobian; filled by each processor and then assembled into the mothership Jacobian;
    Teuchos::RCP<PeridigmNS::SerialMatrix> overlapJacobian;

    //! Global vector for initial positions
    Teuchos::RCP<Epetra_Vector> x;

    //! Global vector for displacement
    Teuchos::RCP<Epetra_Vector> u;

    //! Global vector for current position
    Teuchos::RCP<Epetra_Vector> y;

    //! Global vector for velocity
    Teuchos::RCP<Epetra_Vector> v;

    //! Global vector for acceleration
    Teuchos::RCP<Epetra_Vector> a;

    //! Global vector for force
    Teuchos::RCP<Epetra_Vector> force;

    //! Global vector for contact force (used only in simulations with contact)
    Teuchos::RCP<Epetra_Vector> contactForce;

    //! Global vector for delta u (used only in implicit time integration)
    Teuchos::RCP<Epetra_Vector> deltaU;

    //! Global scratch space vector
    Teuchos::RCP<Epetra_Vector> scratch;

    //! Global vector for block ID 
    Teuchos::RCP<Epetra_Vector> blockIDs;

    //! Global vector for cell volume 
    Teuchos::RCP<Epetra_Vector> volume;

    //! Global vector for cell density
    Teuchos::RCP<Epetra_Vector> density;

    //! Map for global tangent matrix (note, must be an Epetra_Map, not an Epetra_BlockMap)
    Teuchos::RCP<Epetra_Map> tangentMap;

    //! Global tangent matrix
    Teuchos::RCP<Epetra_FECrsMatrix> tangent;

    //! List of neighbors for all locally-owned nodes
    Teuchos::RCP<PeridigmNS::NeighborhoodData> globalNeighborhoodData;

    //! List of potential contact neighbors for all locally-owned nodes
    Teuchos::RCP<PeridigmNS::NeighborhoodData> globalContactNeighborhoodData;

    //! Workset that is passed to the modelEvaluator
    Teuchos::RCP<PHAL::Workset> workset;

    //! The peridigm model evaluator
    Teuchos::RCP<PeridigmNS::ModelEvaluator> modelEvaluator;

    //! The peridigm output manager
    Teuchos::RCP<PeridigmNS::OutputManagerContainer> outputManager;

    //! BLAS for local-only vector updates (BLAS-1)
    Epetra_BLAS blas;

    // field ids for all relevant data
    int blockIdFieldId;
    int volumeFieldId;
    int modelCoordinatesFieldId;
    int coordinatesFieldId;
    int displacementFieldId;
    int velocityFieldId;
    int accelerationFieldId;
    int forceDensityFieldId;
    int contactForceDensityFieldId;
    int partialVolumeFieldId;
    int tangentReferenceCoordinatesFieldId;
  };
}

#endif // PERIDIGM_HPP
