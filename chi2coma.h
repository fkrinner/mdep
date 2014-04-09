#ifndef CHI2COMA_MY_GLOB
#define CHI2COMA_MY_GLOB

#include "chi2.h"
#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON

class chi2coma : public chi2 {
	public:
		void vf();

		chi2coma(const int nWaves);

		template<typename xdouble>
		xdouble EvalComa(std::vector<xdouble> &par);

#ifdef ADOL_ON
		adouble EvalAdouble(std::vector<adouble> par);
#endif//ADOL_ON

		bool CheckConsistency();

                // Indexing is different from chi2, but will be handeled internally:
                // i1,j1, i2, j2 
		void SetComa(int i1, int j1, int i2, int j2, std::vector<double> bins);
		void SetComa2(int ii, int jj, std::vector<double> bins);
		std::vector<double> GetComa(int i1, int j1, int i2, int j2);
		std::vector<double> GetComa2(int ii, int jj);
		void GetComaFromErrors(); // Gets a diagonal covariance matrix from the errors.

//		void vf(){std::cout << "cold from koma" << std::endl;};

	protected:
		std::vector<std::vector<std::vector<double> > > _coma; // INVERTED covariance matrix
};




#endif//CHI2COMA_MY_GLOB
