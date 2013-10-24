/*! \file utPeridigm_State.cpp  without boost library: a test case  */

//@HEADER
// ************************************************************************
//                            
// ************************************************************************
//@HEADER 


#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>
#include "Peridigm_State.hpp"
#include <vector>
#include <string>
#include <math.h>
#include <iostream>


using namespace Teuchos;
using namespace PeridigmNS;
using namespace std;

// TestResult class

class TestResult
{
public:
    TestResult (): failureCount (0), testCount(0){}

    virtual ~TestResult() {};

	virtual void testWasRun (){ testCount++;}
	virtual void startTests (){}
	virtual void addFailure(std::string condition, std::string testName, std::string fileName, long lineNumber){
		
		failureCount++;
		std::cout << fileName.c_str () << "(" << lineNumber << "):" << "Failure: \"" << condition.c_str () << "\" " <<  std::endl;
    }
	virtual void endTests (){

	        std::cerr << "Running " << testCount << " test cases" << std::endl;
    
		if (failureCount > 0)
                    std::cerr << "*** " << failureCount << " failures detected.";
                else
                    std::cerr << "*** No errors detected.";
	}

       int getFailureCount() const { return failureCount; }

protected:
	int failureCount;
	int testCount;
       
};


// TestSetup class

class TestSetup
{
public:
    virtual void setup() = 0;
    virtual void teardown() = 0;
};

//Test class

class Test : public TestSetup
{
public:
	Test (const std::string& testName);

	virtual void	run (TestResult& result);
	virtual void	runTest (TestResult& result) = 0;

protected:
	std::string		name;

};


//TestRegistry class

class TestRegistry
{
public:
	static void addTest (Test *test){
	       instance ().add (test);
	}
	
        static void runAllTests (TestResult& result){
		   instance ().run (result);
	}
      
      

private:

	static TestRegistry&	instance (){
		   static TestRegistry registry;
	       return registry;
	}
	void add (Test *test){
		tests.push_back (test);
	}
	
        void run (TestResult& result){
		 
		result.startTests ();
	    
	       for (std::vector<Test *>::iterator it = tests.begin (); it != tests.end (); ++it)
		    (*it)->run (result);
	       result.endTests ();
	}
       
	std::vector<Test *> tests;
};


// Test class methods

Test::Test (const std::string& testName) : name (testName) 
{

	TestRegistry::addTest (this);
}



void Test::run (TestResult& result) 
{
#ifndef DONT_CATCH_EXCEPTIONS
    try
    {
#endif
                setup();
	        runTest (result);
		
#ifndef DONT_CATCH_EXCEPTIONS
    }
    catch (...) {
                result.addFailure ("Unhandled exception", name, "", 0);
    }
#endif
                teardown();
	        result.testWasRun();
}

#define TESTWITHSETUP(name)\
	class name##Test : public Test, name##Setup\
	{ \
		public: \
			name##Test () : Test (#name) {} \
                        void setup() {name##Setup::setup();} \
                        void teardown() {name##Setup::teardown();} \
			void runTest (TestResult& result_); \
	} name##Instance; \
	void name##Test::runTest (TestResult& result_)


// Here is a collection of testing macros that can be used in the  bodies of tests.  

#ifndef DONT_CATCH_EXCEPTIONS

#define CHECK(condition) \
    try { \
    if (!(condition)) \
       result_.addFailure (#condition, name, __FILE__, __LINE__); \
    } catch(...) { \
        result_.addFailure ("Unhandled exception", name, __FILE__, __LINE__); \
    }


#define CHECK_LONGS_EQUAL(expected,actual)\
{\
    try { \
        long _expected = (expected);\
        long _actual = (actual);\
        if (_expected != _actual) {\
            char message [80];\
            sprintf (message, "expected %ld but was: %ld", _expected, _actual);\
            result_.addFailure (message, name, __FILE__, __LINE__);\
        }\
    } catch(...) { \
    result_.addFailure ("Unhandled exception", name, __FILE__, __LINE__); \
    }\
}


#define CHECK_DOUBLES_EQUAL(expected,actual)\
{\
    try { \
        double _expected = (expected);\
        double _actual = (actual);\
        if (fabs ((_expected)-(_actual)) > 0.001) {\
            char message [80];\
            sprintf (message, "expected %lf but was: %lf", (_expected), (_actual));\
            result_.addFailure (message, name, __FILE__, __LINE__);\
        }\
    } catch(...) { \
    result_.addFailure ("Unhandled exception", name, __FILE__, __LINE__); \
    }\
}

#define CHECK_DOUBLES_CLOSE(expected,actual, diff)\
{\
    try { \
        double _expected = (expected);\
        double _actual = (actual);\
        if (fabs ((_expected)-(_actual)) > diff) {\
            char message [80];\
            sprintf (message, "expected %lf but was: %lf", (_expected), (_actual));\
            result_.addFailure (message, name, __FILE__, __LINE__);\
        }\
    } catch(...) { \
    result_.addFailure ("Unhandled exception", name, __FILE__, __LINE__); \
    }\
}

#else

#define CHECK(condition) \
	if (!(condition)) \
		result_.addFailure (#condition, name, __FILE__, __LINE__);


#define CHECK_LONGS_EQUAL(expected,actual)\
{\
	long _expected = (expected);\
	long _actual = (actual);\
	if (_expected != _actual) {\
		char message [80];\
		sprintf (message, "expected %ld but was: %ld", _expected, _actual);\
		result_.addFailure (message, name, __FILE__, __LINE__);\
	}\
}


#define CHECK_DOUBLES_EQUAL(expected,actual)\
{\
	double _expected = (expected);\
	double _actual = (actual);\
	if (fabs ((_expected)-(_actual)) > 0.001) {\
		char message [80];\
		sprintf (message, "expected %lf but was: %lf", (_expected), (_actual));\
		result_.addFailure (message, name, __FILE__, __LINE__);\
	}\
}

#define CHECK_DOUBLES_CLOSE(expected,actual, diff)\
{\
        double _expected = (expected);\
        double _actual = (actual);\
        if (fabs ((_expected)-(_actual)) > diff) {\
            char message [80];\
            sprintf (message, "expected %lf but was: %lf", (_expected), (_actual));\
            result_.addFailure (message, name, __FILE__, __LINE__);\
        }\
}


#endif // DONT_CATCH_EXCEPTIONS


// TwoPointProblem test class

class TwoPointProblemSetup:public TestSetup{


public:
      PeridigmNS::State createTwoPointProblem(){

      Teuchos::RCP<Epetra_Comm> comm;
           #ifdef HAVE_MPI
           comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
           #else
           comm = rcp(new Epetra_SerialComm);
           #endif
         
          
           // set up overlap maps, which include ghosted nodes

           int numGlobalElements(2), numMyElements(2), elementSize(1), indexBase(0);
           std::vector<int> myGlobalElements(numMyElements);
           for(int i=0; i<numMyElements ; ++i)   myGlobalElements[i] = i;
           
           // overlapScalarPointMap
           // used for cell volumes and scalar constitutive data

           overlapScalarPointMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));   

          // overlapVectorPointMap
          // used for positions, displacements, velocities and vector constitutive data
          elementSize = 3;
          overlapVectorPointMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));
         
          // ownedScalarBondMap
          // used for bond damage and bond constitutive data

          std::vector<int> bondElementSize(numMyElements, 1);

          ownedScalarBondMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], &bondElementSize[0], indexBase, *comm));
          PeridigmNS::State state;
                             
          FieldManager& fm = FieldManager::self();

           // create a list of scalar field specs and allocate the data

          

          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Element_Id") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Proc_Num") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Number_Of_Neighbors") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Lambda") );
          scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Surface_Correction_Factor") );

          state.allocateScalarPointData(scalarPointFieldIds, overlapScalarPointMap);


          // create a list of vector field specs and allocate the data
         

         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Contact_Force_Density") );
      
         state.allocateVectorPointData(vectorPointFieldIds, overlapVectorPointMap);

         
         // create a list of bond field specs and allocate the data

         bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage") );
         bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deviatoric_Plastic_Extension") );

         state.allocateScalarBondData(bondFieldIds, ownedScalarBondMap);

         return state;

     }
      void setup(){

          state = createTwoPointProblem();
 
          FieldManager& fm = FieldManager::self();
      
          coordinatesFieldId = fm.getFieldId("Coordinates");
          elementIdFieldId = fm.getFieldId("Element_Id");
          forceDensityFieldId = fm.getFieldId("Force_Density");
          damageFieldId = fm.getFieldId("Damage");

          // Initialization of coordinates that needs to be checked later

          coordinates = state.getData(coordinatesFieldId);

         
          
          // Set scalar data
   
          Epetra_Vector& ids = *(state.getData(elementIdFieldId));

          for(int i=0 ; i<ids.MyLength() ; ++i)  ids[i] = i;
          
          // Set vector data

          Epetra_Vector& force = *(state.getData(forceDensityFieldId));
          
          
          for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
              force[i*3] = i*3;
              force[i*3+1] = i*3+1;
              force[i*3+2] = i*3+2;
          }

          // Set bond data

        Epetra_Vector& bondDamage = *(state.getData(damageFieldId));

          for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){

              int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
             
              for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
                  bondDamage[firstPointInElement+j] = firstPointInElement+j;
          }
   

  }

void teardown(){}
      
protected:

     PeridigmNS::State state;
     
     
     Teuchos::RCP<Epetra_BlockMap> overlapScalarPointMap;
     Teuchos::RCP<Epetra_BlockMap> overlapVectorPointMap;
     Teuchos::RCP<Epetra_BlockMap> ownedScalarBondMap;
     Teuchos::RCP<Epetra_Vector> coordinates;

     vector<int> scalarPointFieldIds;
     vector<int> vectorPointFieldIds;
     vector<int> bondFieldIds;

     int coordinatesFieldId;
     int elementIdFieldId;
     int forceDensityFieldId;
     int damageFieldId;

     
   
};


TESTWITHSETUP (TwoPointProblem)
{
    
    CHECK_LONGS_EQUAL(state.getScalarPointMultiVector()->NumVectors(), (int)scalarPointFieldIds.size());
    CHECK_LONGS_EQUAL(state.getScalarPointMultiVector()->MyLength(), overlapScalarPointMap->NumMyPoints());
    CHECK( state.getScalarPointMultiVector()->Map().SameAs( *overlapScalarPointMap ) );

    CHECK_LONGS_EQUAL( state.getVectorPointMultiVector()->NumVectors(), (int)vectorPointFieldIds.size() );
    CHECK_LONGS_EQUAL( state.getVectorPointMultiVector()->MyLength(), overlapVectorPointMap->NumMyPoints() );
    CHECK( state.getVectorPointMultiVector()->Map().SameAs( *overlapVectorPointMap ) );

    CHECK_LONGS_EQUAL( state.getScalarBondMultiVector()->NumVectors(), (int)bondFieldIds.size() );
    CHECK_LONGS_EQUAL( state.getScalarBondMultiVector()->MyLength(), ownedScalarBondMap->NumMyPoints() );
    CHECK( state.getScalarBondMultiVector()->Map().SameAs( *ownedScalarBondMap ) );

    for(int i=0 ; i<coordinates->MyLength() ; ++i)
        CHECK_DOUBLES_EQUAL( (*coordinates)[i], 0.0 );

    // Check scalar data
    Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    for(int i=0 ; i<ids.MyLength() ; ++i)
        CHECK_DOUBLES_CLOSE(ids[i], (double)(i), 1.0e-14);

    // Check vector data
    Epetra_Vector& force = *(state.getData(forceDensityFieldId));

    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      
        CHECK_DOUBLES_CLOSE(force[i*3], (double)(i*3), 1.0e-14);
        CHECK_DOUBLES_CLOSE(force[i*3+1], (double)(i*3+1), 1.0e-14);
        CHECK_DOUBLES_CLOSE(force[i*3+2], (double)(i*3+2), 1.0e-14);
    }

    // Check bond data
    Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
          CHECK_DOUBLES_CLOSE(bondDamage[firstPointInElement+j], (double)(firstPointInElement+j), 1.0e-14);
    }
}


// ThreePointProblem test class

class ThreePointProblemSetup:public TestSetup{


public:
      PeridigmNS::State createThreePointProblem(){

           Teuchos::RCP<Epetra_Comm> comm;
           #ifdef HAVE_MPI
           comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
           #else
           comm = rcp(new Epetra_SerialComm);
           #endif
           int numProcs = comm->NumProc();
           int myPID = comm->MyPID();
  

          // set up overlap maps, which include ghosted nodes

          int numGlobalElements(3), numMyElements(0), elementSize(1), indexBase(0);
          std::vector<int> myGlobalElements;

          if(numProcs == 1){
             numMyElements = 3;
             myGlobalElements.resize(numMyElements);
             myGlobalElements[0] = 0;
             myGlobalElements[1] = 1;
             myGlobalElements[2] = 2;
          }
          else if(numProcs == 2){
                  if(myPID == 0){
                     numMyElements = 2;
                     myGlobalElements.resize(numMyElements);
                     myGlobalElements[0] = 0;
                     myGlobalElements[1] = 1;
                  }
          else if(myPID == 1){
                  numMyElements = 2;
                  myGlobalElements.resize(numMyElements);
                  myGlobalElements[0] = 2;
                  myGlobalElements[1] = 1;
              }
          }
         
          // overlapScalarPointMap
          // used for cell volumes and scalar constitutive data
         
          overlapScalarPointMap = Teuchos::rcp(new Epetra_BlockMap(-1, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));

          // overlapVectorPointMap
          // used for positions, displacements, velocities and vector constitutive data

          elementSize = 3;
     
          overlapVectorPointMap = Teuchos::rcp(new Epetra_BlockMap(-1, numMyElements, &myGlobalElements[0], elementSize, indexBase, *comm));
         
          // ownedScalarBondMap
          // used for bond damage and bond constitutive data

          std::vector<int> bondElementSize;

          if(numProcs == 1){
             numMyElements = 3;
             myGlobalElements.resize(numMyElements);
             myGlobalElements[0] = 0;
             myGlobalElements[1] = 1;
             myGlobalElements[2] = 2;
             bondElementSize.resize(numMyElements);
             bondElementSize[0] = 1;
             bondElementSize[1] = 2;
             bondElementSize[2] = 1;
          }
          else if(numProcs == 2){
                  if(myPID == 0){
                     numMyElements = 2;
                     myGlobalElements.resize(numMyElements);
                     myGlobalElements[0] = 0;
                     myGlobalElements[1] = 1;
                     bondElementSize.resize(numMyElements);
                     bondElementSize[0] = 1;
                     bondElementSize[1] = 2;
                  }
                  else if(myPID == 1){
                          numMyElements = 1;
                          myGlobalElements.resize(numMyElements);
                          myGlobalElements[0] = 2;
                          bondElementSize.resize(numMyElements);
                          bondElementSize[0] = 1;
                  }
          }

          ownedScalarBondMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &myGlobalElements[0], &bondElementSize[0], indexBase, *comm));
          PeridigmNS::State state;
                             
          FieldManager& fm = FieldManager::self();

          // create a list of scalar field specs and allocate the data

        

         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Element_Id") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Proc_Num") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Number_Of_Neighbors") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Lambda") );
         scalarPointFieldIds.push_back( fm.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Surface_Correction_Factor") );
  
         state.allocateScalarPointData(scalarPointFieldIds, overlapScalarPointMap);


          // create a list of vector field specs and allocate the data

         
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density") );
         vectorPointFieldIds.push_back( fm.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Contact_Force_Density") );
         
         state.allocateVectorPointData(vectorPointFieldIds, overlapVectorPointMap);
      
         // create a list of bond field specs and allocate the data

         bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage") );
         bondFieldIds.push_back( fm.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deviatoric_Plastic_Extension") );

         state.allocateScalarBondData(bondFieldIds, ownedScalarBondMap);

         return state;

     }
      void setup(){

         state = createThreePointProblem();

         FieldManager& fm = FieldManager::self();
         coordinatesFieldId = fm.getFieldId("Coordinates");
         elementIdFieldId = fm.getFieldId("Element_Id");
         forceDensityFieldId = fm.getFieldId("Force_Density");
         damageFieldId = fm.getFieldId("Damage");

         // Initialization of coordinates that needs to be checked later

         coordinates = state.getData(coordinatesFieldId);

         
         // Set scalar data

        Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    
        for(int i=0 ; i<ids.MyLength() ; ++i) ids[i] = i;

        // Set vector data

        Epetra_Vector& force = *(state.getData(forceDensityFieldId));

        for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
            force[i*3] = i*3;
            force[i*3+1] = i*3+1;
            force[i*3+2] = i*3+2;
        }

        // bond data

       Epetra_Vector& bondDamage = *(state.getData(damageFieldId));

       for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
           int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
           for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
               bondDamage[firstPointInElement+j] = firstPointInElement+j;
       }

     }

     void teardown(){}
      
protected:

     PeridigmNS::State state;
     
     Teuchos::RCP<Epetra_BlockMap> overlapScalarPointMap;
     Teuchos::RCP<Epetra_BlockMap> overlapVectorPointMap;
     Teuchos::RCP<Epetra_BlockMap> ownedScalarBondMap;
     Teuchos::RCP<Epetra_Vector> coordinates;

     vector<int> scalarPointFieldIds;
     vector<int> vectorPointFieldIds;
     vector<int> bondFieldIds;
     

     int coordinatesFieldId;
     int elementIdFieldId;
     int forceDensityFieldId;
     int damageFieldId; 
};

TESTWITHSETUP (ThreePointProblem)
{
    
    CHECK_LONGS_EQUAL(state.getScalarPointMultiVector()->NumVectors(), (int)scalarPointFieldIds.size());
    CHECK_LONGS_EQUAL(state.getScalarPointMultiVector()->MyLength(), overlapScalarPointMap->NumMyPoints());
    CHECK( state.getScalarPointMultiVector()->Map().SameAs( *overlapScalarPointMap ) );

    CHECK_LONGS_EQUAL( state.getVectorPointMultiVector()->NumVectors(), (int)vectorPointFieldIds.size() );
    CHECK_LONGS_EQUAL( state.getVectorPointMultiVector()->MyLength(), overlapVectorPointMap->NumMyPoints() );
    CHECK( state.getVectorPointMultiVector()->Map().SameAs( *overlapVectorPointMap ) );

    CHECK_LONGS_EQUAL( state.getScalarBondMultiVector()->NumVectors(), (int)bondFieldIds.size() );
    CHECK_LONGS_EQUAL( state.getScalarBondMultiVector()->MyLength(), ownedScalarBondMap->NumMyPoints() );
    CHECK( state.getScalarBondMultiVector()->Map().SameAs( *ownedScalarBondMap ) );

    for(int i=0 ; i<coordinates->MyLength() ; ++i)
        CHECK_DOUBLES_EQUAL( (*coordinates)[i], 0.0 );

    // Check scalar data

    Epetra_Vector& ids = *(state.getData(elementIdFieldId));
    for(int i=0 ; i<ids.MyLength() ; ++i)
        CHECK_DOUBLES_CLOSE(ids[i], (double)(i), 1.0e-14);

    // Check vector data
    Epetra_Vector& force = *(state.getData(forceDensityFieldId));

    for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
      
        CHECK_DOUBLES_CLOSE(force[i*3], (double)(i*3), 1.0e-14);
        CHECK_DOUBLES_CLOSE(force[i*3+1], (double)(i*3+1), 1.0e-14);
        CHECK_DOUBLES_CLOSE(force[i*3+2], (double)(i*3+2), 1.0e-14);
    }

    // Check bond data
    Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
    for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
      int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
      for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j)
          CHECK_DOUBLES_CLOSE(bondDamage[firstPointInElement+j], (double)(firstPointInElement+j), 1.0e-14);
    }
}


//  start copyFrom

class CopyFromSetup:public TestSetup{

public:
      
      void setup(){

         state = threePointProblem.createThreePointProblem();

         FieldManager& fm = FieldManager::self();

         elementIdFieldId = fm.getFieldId("Element_Id");
         forceDensityFieldId = fm.getFieldId("Force_Density");
         damageFieldId = fm.getFieldId("Damage");

         // set some scalar data

         Epetra_Vector& ids = *(state.getData(elementIdFieldId));

         for(int i=0 ; i<ids.MyLength() ; ++i) ids[i] = i;

         // set some vector data

         Epetra_Vector& force = *(state.getData(forceDensityFieldId));

         for(int i=0 ; i<force.Map().NumMyElements() ; ++i){
             force[i*3] = i*3;
             force[i*3+1] = i*3+1;
             force[i*3+2] = i*3+2;
         }
    
        // set some bond data

        Epetra_Vector& bondDamage = *(state.getData(damageFieldId));
             
        for(int i=0 ; i<bondDamage.Map().NumMyElements() ; ++i){
            int firstPointInElement = bondDamage.Map().FirstPointInElement(i);
            for(int j=0 ; j<bondDamage.Map().ElementSize(i); ++j){
                bondDamage[firstPointInElement+j] = firstPointInElement+j;
            }
        }

        // create a temporary State
        // this mimics what is done behind the scenes in the finite-difference Jacobian routine
        // there is one owned point, processor 0 has GID 0, and processor 1 (if it exists) has GID 2
        // in either case there is one ghosted point, its neighbor, GID = 1

        int procID = state.getScalarPointMultiVector()->Comm().MyPID();
        Epetra_SerialComm serialComm;
        int numOwnedIDs = 1;
        std::vector<int> tempMyGlobalIDs(2); // includes ghost

        // main owned point

        if(procID == 0)
           tempMyGlobalIDs[0] = 0;
        else if(procID == 1)
           tempMyGlobalIDs[0] = 2;

        // neighbor

        tempMyGlobalIDs[1] = 1;
        std::vector<int> bondElementSize(1);
        bondElementSize[0] = 1;

        tempOverlapScalarPointMap = Teuchos::rcp(new Epetra_BlockMap(tempMyGlobalIDs.size(), tempMyGlobalIDs.size(), &tempMyGlobalIDs[0], 1, 0, serialComm));
        tempOverlapVectorPointMap = Teuchos::rcp(new Epetra_BlockMap(tempMyGlobalIDs.size(), tempMyGlobalIDs.size(), &tempMyGlobalIDs[0], 3, 0, serialComm));
        tempOwnedScalarBondMap = Teuchos::rcp(new Epetra_BlockMap(numOwnedIDs, numOwnedIDs, &tempMyGlobalIDs[0], &bondElementSize[0], 0, serialComm));

        PeridigmField::Relation relation;
        PeridigmField::Length length;

        // allocate data for the field specs, which are taken from the initial state object

        // scalar point data

        relation = PeridigmField::ELEMENT;
        length = PeridigmField::SCALAR;
        vector<int> scalarPointFieldIds = state.getFieldIds(relation, length);
        tempState.allocateScalarPointData(scalarPointFieldIds, tempOverlapScalarPointMap);

        // vector point data

        relation = PeridigmField::NODE;
        length = PeridigmField::VECTOR;
        vector<int> vectorFieldIds = state.getFieldIds(relation, length);
        tempState.allocateVectorPointData(vectorFieldIds, tempOverlapVectorPointMap);

        // scalar bond data

        relation = PeridigmField::BOND;
        length = PeridigmField::SCALAR;
        vector<int> bondFieldIds = state.getFieldIds(relation, length);
        tempState.allocateScalarBondData(bondFieldIds, tempOwnedScalarBondMap);

     }

     void teardown(){}
      
protected:

     ThreePointProblemSetup threePointProblem;

     PeridigmNS::State state;
     int elementIdFieldId;
     int forceDensityFieldId;
     int damageFieldId;


     Teuchos::RCP<Epetra_BlockMap> tempOverlapScalarPointMap;
     Teuchos::RCP<Epetra_BlockMap> tempOverlapVectorPointMap;
     Teuchos::RCP<Epetra_BlockMap> tempOwnedScalarBondMap;
     
     PeridigmNS::State tempState;
        
};


TESTWITHSETUP (CopyFrom)
{
    CHECK_LONGS_EQUAL( tempState.getVectorPointMultiVector()->NumVectors(), state.getVectorPointMultiVector()->NumVectors());
    CHECK_LONGS_EQUAL( tempState.getVectorPointMultiVector()->MyLength(), tempOverlapVectorPointMap->NumMyPoints() );
    CHECK( tempState.getVectorPointMultiVector()->Map().SameAs( *tempOverlapVectorPointMap ) );

    CHECK_LONGS_EQUAL( tempState.getVectorPointMultiVector()->NumVectors(), state.getVectorPointMultiVector()->NumVectors());
    CHECK_LONGS_EQUAL( tempState.getVectorPointMultiVector()->MyLength(), tempOverlapVectorPointMap->NumMyPoints() );
    CHECK( tempState.getVectorPointMultiVector()->Map().SameAs( *tempOverlapVectorPointMap ) );

    CHECK_LONGS_EQUAL( tempState.getScalarBondMultiVector()->NumVectors(), state.getScalarBondMultiVector()->NumVectors());
    CHECK_LONGS_EQUAL( tempState.getScalarBondMultiVector()->MyLength(), tempOwnedScalarBondMap->NumMyPoints() );
    CHECK( tempState.getScalarBondMultiVector()->Map().SameAs( *tempOwnedScalarBondMap ) );

    // copy the data from the initial State to the smaller, temporary state
  tempState.copyLocallyOwnedDataFromState(Teuchos::RCP<PeridigmNS::State>(&state, false));

  // check the data
  vector<int> tempFieldIds = tempState.getFieldIds();
  for(unsigned int iId=0 ; iId<tempFieldIds.size() ; ++iId){
    Teuchos::RCP<Epetra_Vector> data = state.getData( tempFieldIds[iId] );
    Teuchos::RCP<Epetra_Vector> tempData = tempState.getData( tempFieldIds[iId] );
    for(int iLID=0 ; iLID<tempData->Map().NumMyElements() ; ++iLID){
      int globalID = tempData->Map().GID(iLID);
      CHECK(globalID != -1);
      int stateLID = data->Map().LID(globalID);
      CHECK(stateLID != -1);
      int tempElementSize = tempData->Map().ElementSize(iLID);
      int stateElementSize = data->Map().ElementSize(stateLID);
      CHECK_LONGS_EQUAL( tempElementSize, stateElementSize );
      int tempFirstPointInElement = tempData->Map().FirstPointInElement(iLID);
      int stateFirstPointInElement = data->Map().FirstPointInElement(stateLID);
      for(int i=0 ; i<tempElementSize ; ++i){
         CHECK_LONGS_EQUAL( (*tempData)[tempFirstPointInElement + i], (*data)[stateFirstPointInElement + i] );
      }
    }
  }

}

// end copyFrom

int main(int argc, char* argv[])
{   
    int numProcs = 1;

    #ifdef HAVE_MPI
         MPI_Init(&argc,&argv);
         MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    #endif

    int returnCode = -1;
   
    TestResult result;
    TestRegistry::runAllTests(result);  
     
   
    if(numProcs == 1 || numProcs == 2){
       returnCode = result.getFailureCount();
    }
    else{
       std::cerr << "Unit test runtime ERROR: utPeridigm_State only makes sense on 1 or 2 processors." << std::endl;
    }

    #ifdef HAVE_MPI
       MPI_Finalize();
    #endif    
             
    return returnCode;  
}





