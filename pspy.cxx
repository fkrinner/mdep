#include<boost/python.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>


#ifdef ADOL_ON
#define ADOL_IS_ON
#undef ADOL_ON // Does currently not work with adouble on.
#endif//ADOL_ON

#include "phaseSpace.h"

#ifdef ADOL_IS_ON
#undef ADOL_IS_ON
#define ADOL_ON
#endif//ADOL_IS_ON

namespace bp = boost::python;

double phaseSpace2(double m, int i){
	return phaseSpace(m,i);
};

BOOST_PYTHON_MODULE(libpspy){
	bp::def("phaseSpace",	phaseSpace	);
	bp::def("phaseSpace",	phaseSpace2	);
};
