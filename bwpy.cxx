#include<boost/python.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>

#ifdef ADOL_ON
#define ADOL_IS_ON
#undef ADOL_ON // Does currently not work with adouble on.
#endif//ADOL_ON

#include "breitWigners.h"

#ifdef ADOL_IS_ON
#undef ADOL_IS_ON
#define ADOL_ON
#endif//ADOL_IS_ON


namespace bp = boost::python;

template< typename T >
std::vector<T> to_std_vector(const bp::object& iterable){
    return std::vector<T>(bp::stl_input_iterator<T>(iterable), bp::stl_input_iterator<T>( ));
};

std::complex<double> bwpy(double m, bp::list par, int model){
	std::vector<double> xxx =  to_std_vector<double>(par);
	return bw<double>(m,xxx,model);
//	return std::complex<double>(1.,0.);
};



BOOST_PYTHON_MODULE(libbwpy){
	bp::def("bw"		,bwpy		);
	bp::def("getNpars"	,getNpars	);
	bp::def("getModelName"	,getModelName	);
};
