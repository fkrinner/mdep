#include<boost/python.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>

#define IS_AMP_CLASS_KUL // Enabels methods specific for fit on Amplitudes
#include "chi_to_py.h"
#include "chi2amp.h"


struct chi2amppy : public chi_to_py<chi2amp>, bp::wrapper<chi_to_py<chi2amp> >{
	chi2amppy(const chi_to_py<chi2amp> chi2ampin):
		chi_to_py<chi2amp>(chi2ampin),
		bp::wrapper<chi_to_py<chi2amp> >(){};
	chi2amppy(const int nWaves):
		chi_to_py<chi2amp>(nWaves),
		bp::wrapper<chi_to_py<chi2amp> >(){};
};
BOOST_PYTHON_MODULE(libchi2amppy){
	export_chi2<chi2amppy>("chi2amp");
};
