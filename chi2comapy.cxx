#include<boost/python.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>

#define IS_COMA_CLASS_FRRR // Enables methods specific for covariance matrices
#define IS_RHO_CLASS_QWERT // Enabels methods specific for fit on Spin-Density-Matrices
#include "chi_to_py.h"
#include "chi2coma.h"


struct chi2comapy : public chi_to_py<chi2coma>, bp::wrapper<chi_to_py<chi2coma> > {
	chi2comapy(const chi_to_py<chi2coma> chi2comain):
		chi_to_py<chi2coma>(chi2comain),
		bp::wrapper<chi_to_py<chi2coma> >(){};
	chi2comapy(const int nWaves):
		chi_to_py<chi2coma>(nWaves),
		bp::wrapper<chi_to_py<chi2coma> >(){};
};

BOOST_PYTHON_MODULE(libchi2comapy){
	export_chi2<chi2comapy>("chi2coma");
};
