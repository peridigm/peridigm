#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>

using Teuchos::ArrayRCP;
using Teuchos::RCP;

using namespace boost::unit_test;

class A {
public:
	A() {}
	~A() { std::cout << "~A" << std::endl; }
};

void smartPointer()
{
	std::size_t size=10;
	double fill(0.0);
	ArrayRCP<double> pPtr(size,fill);
	RCP<A> aPtr(new A());
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utTeuchosSmartPointer" );
	proc->add(BOOST_TEST_CASE( &smartPointer ));
	framework::master_test_suite().add( proc );
	return success;
}


bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(
		int argc,
		char* argv[]
)
{


	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
