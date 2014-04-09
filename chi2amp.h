#ifndef CHI2_TOP_KEK
#define CHI2_TOP_KEK

#include "waveset.h"
#include <string>
#include <complex>

#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON

class chi2amp : public waveset{
	public:
		chi2amp(const int nWaves);

		template<typename xdouble>//Evaluate Chi2 with all parameters given in 'par'
		xdouble Eval(std::vector<xdouble> &par);

//		template<typename xdouble>//Evaluate Chi2 with all parameters given in 'par'
//		xdouble Eval2(std::vector<xdouble> &par);

		template<typename xdouble>//Evaluate Chi2 with only released parametrers given in 'par'. (In '_waves[i]._release_parameter[j])
		xdouble EvalRel(std::vector<xdouble> &par);

		template<typename xdouble>//Evaluate Chi2 with only released parameters from 0. to 1. given in par. They will be rescaled by the interval given through '_waves[i]._upper/lowerLim[j]'
		xdouble EvalNorm(std::vector<xdouble> &par);

		template<typename xdouble>//GetDeltas
		std::vector<xdouble> Deltas(int bin,std::vector<xdouble> &par);


#ifdef ADOL_ON
		virtual adouble EvalAdouble(std::vector<adouble> par);
#endif//ADOL_ON

		std::string ClassName();

		void SetData(int i, std::vector<double> data);
		void SetComa(int i, int j, std::vector<double> errors);
		std::vector<double> GetData(int i);
		std::vector<double> GetComa(int i, int j);

		bool CheckConsistency();
		
		std::vector<std::complex<double> > GetPoints(int i,int nPoints = 1000);
		std::vector<std::complex<double> > GetPointsComponent(int i, std::vector<int> comp, int nPoints = 1000);

		void SetRandomCouplings();
		void SetCouplingSize();


	protected:
		int _anchorWave;						// Number of anchor wave
		std::vector<std::vector<double > > _points;			// Data points to fit	// _points[i][j] = real part for i>j, intensity for i==j, imag part for i>j
		std::vector<std::vector<std::vector<double> > > _coma;	// Corresponding errors	// _points[i][j] = real part for i>j, intensity for i==j, imag part for i>j


};


#endif //CHI2_TOP_KEK
