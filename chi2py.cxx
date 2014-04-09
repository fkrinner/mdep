#include<boost/python.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>

#define IS_RHO_CLASS_QWERT // Enabels methods specific for fit on Spin-Density-Matrices
#include "chi_to_py.h"
#include "chi2.h"


struct chi2py : public chi_to_py<chi2>, bp::wrapper<chi_to_py<chi2> >{
	chi2py(const chi_to_py<chi2> chi2in):
		chi_to_py<chi2>(chi2in),
		bp::wrapper<chi_to_py<chi2> >(){};
	chi2py(const int nWaves):
		chi_to_py<chi2>(nWaves),
		bp::wrapper<chi_to_py<chi2> >(){};
};

BOOST_PYTHON_MODULE(libchi2py){
	export_chi2<chi2py>("chi2");
};
