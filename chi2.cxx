#include "chi2.h"
#include "waveset.h"
#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include <fstream>

#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON


chi2::chi2(const int nWaves): waveset(nWaves) ,_anchorWave(-1) {
	for (int i=0; i < _nWaves;i++){
		_points.push_back(std::vector<std::vector<double> >());
		_errors.push_back(std::vector<std::vector<double> >());
	};
	for (int i=0; i<_nWaves;i++){
		for(int j=0; j<_nWaves;j++){
			_points[i].push_back(std::vector<double>());
			_errors[i].push_back(std::vector<double>());
		};
	};
};


template<typename xdouble> // Using Deltas(...) to calculate the chi2.
xdouble chi2::Eval(std::vector<xdouble> &par){
	xdouble chi2_tot = 0.;
	for (int bin=0; bin<_nBins; bin++){			// Loop over mass-bins
		std::vector<xdouble> deltas = Deltas<xdouble>(bin,par);
		for (int i=0; i<_nWaves*_nWaves;i++){
			// Do not use values, with delta == 0. (bin is excluded or would not contribute anyway)
			int ii = i/_nWaves;
			int jj = i - _nWaves*ii;
			if (_errors[ii][jj][bin] != 0.){			
				chi2_tot+= pow(deltas[i]/_errors[ii][jj][bin],2);
			};
		};
	};
	if(_nCalls%1000 == 0){
		std::cout << "Function call #"<<_nCalls<<": "<<chi2_tot<<std::endl;
	};
	return chi2_tot;	
};
template double chi2::Eval(std::vector<double> &par);

template<typename xdouble> // Gives the differences between points and fit-curve.
std::vector<xdouble> chi2::Deltas(int bin, std::vector<xdouble> &par){
	double min = _binning[bin];
	double max = _binning[bin+1];
	double m   = (min+max)/2;
	std::vector<std::complex<xdouble> >amps = Amps(bin,par);
	if(_anchorWave > -1){	
		std::complex<xdouble> abs = std::complex<xdouble>(std::norm(amps[_anchorWave]),0.);
		std::complex<xdouble> phaseFactor = amps[_anchorWave]* pow(abs,-.5);//Still wrong: norm = re^2 + im^2
		for(int ww=0;ww<_nWaves;ww++){
			amps[ww]/=phaseFactor;
		};
	};
	std::vector<xdouble> deltas(_nWaves*_nWaves);
	for (int i=0; i<_nWaves; i++){
		xdouble norm = std::norm(amps[i]);
		double iMax = _waves[i].GetMaxLimit();
		double iMin = _waves[i].GetMinLimit();
		if (m > iMin and m < iMax){
			deltas[i*_nWaves+i] = std::norm(amps[i]) - _points[i][i][bin];
//			std::cout << "    enter: m="<< m << " mmin="<< iMin << " mmax=" << iMax <<std::endl;
		}else{
			deltas[i*_nWaves+i] = 0.;	
//			std::cout << "not enter: m="<< m << " mmin="<< iMin << " mmax=" << iMax <<std::endl;
		};
		for (int j=0;j<i;j++){
			double jMax = _waves[j].GetMaxLimit();
			double jMin = _waves[j].GetMinLimit();
			std::complex<xdouble> inter = std::conj(amps[i])*amps[j];
			if (m > jMin and m < jMax and m > iMin and m < iMax){
				deltas[i*_nWaves+j] = std::real(inter)-_points[i][j][bin];
				deltas[j*_nWaves+i] = std::imag(inter)-_points[j][i][bin];
			}else{
				deltas[i*_nWaves+j] = 0.;
				deltas[j*_nWaves+i] = 0.;
			};
		};
	};
	return deltas; // Exclude bins with delta == 0.
};
template std::vector<double> chi2::Deltas(int bin, std::vector<double> &par);

template<typename xdouble> //Take a vector of only relreased parameters as argument. Is equal to Eval, if all parameters are released.
xdouble chi2::EvalRel(std::vector<xdouble> &par){
	int actRel=0;
	std::vector<xdouble> xRel;
	for (int i=0;i<_nWaves;i++){
		for (int j=0;j<_waves[i].GetNpars();j++){
			if(_waves[i].GetParStat(j)){
				xRel.push_back(par[actRel]);
				actRel+=1;
			}else{
				xdouble act=(xdouble)(_waves[i].GetParameter(j));
				xRel.push_back(act);
			};
		};
	};
	return Eval<xdouble>(xRel);
};
template double chi2::EvalRel(std::vector<double> &par);

template<typename xdouble>// Equal to EvalRel(GetParFromNorm(par))
xdouble chi2::EvalNorm(std::vector<xdouble> &par){
	int actRel=0;
	std::vector<xdouble> xRel;
	for (int i=0;i<_nWaves;i++){
		for (int j=0;j<_waves[i].GetNpars();j++){
			if(_waves[i].GetParStat(j)){
				xdouble scaledPar = _waves[i].GetLowerLim(j)*(1-par[actRel]) + _waves[i].GetUpperLim(j)*par[actRel]; 
				xRel.push_back(scaledPar);
				actRel+=1;
			}else{
				xdouble act=(xdouble)(_waves[i].GetParameter(j));
				xRel.push_back(act);
			};
		};
	};
	return Eval<xdouble>(xRel);
};
template double chi2::EvalNorm(std::vector<double> &par);

#ifdef ADOL_ON
adouble chi2::EvalAdouble(std::vector<adouble> par){
	return Eval(par);
};
#endif//ADOL_ON

std::string chi2::ClassName(){
	return "chi2";
};


void chi2::SetData(int i, int j, std::vector<double> data){
//	if (data.size() != _nBins){
//		std::cout << "Warning: Size of new data does not match number of bins." << std::endl;
//	};
	_points[i][j] = data;
};
void chi2::SetError(int i, int j, std::vector<double> error){
//	if (error.size() != _nBins){
//		std::cout << "Warning: Size of new errors does not match number of bins." << std::endl;
//	};
	_errors[i][j] = error;
};

std::vector<double> chi2::GetData(int i, int j){
	return _points[i][j];
};
std::vector<double> chi2::GetError(int i, int j){
	return _errors[i][j];
};
void chi2::Conjugate(int i, int j){
	for (unsigned int bin=0;bin<_points[i][j].size();bin++){
		_points[i][j][bin]*=-1;
	};
};

bool chi2::CheckConsistency(){
	int nErr=0;
	if (not waveset::CheckConsistency()){
		nErr+=1;
	};
	for (int i=0;i<_nWaves;i++){
		for (int j=0;j<_nWaves;j++){
			if (_points[i][j].size() != _nBins){
				std::cout << "Inconsistency found: Data points in _points["<<i<<"]["<<j<<"] does not match _nBins" << std::endl;
				nErr+=1;
			};
			if (_errors[i][j].size() != _nBins){
				std::cout << "Inconsistency found: Errors in _errors["<<i<<"]["<<j<<"] does not match _nBins" << std::endl;
				nErr+=1;
			};

		};
	};
	for (int i =0; i<_nWaves;i++){
		double upLimI = _waves[i].GetMaxLimit();
		double lowLimI= _waves[i].GetMinLimit();
		for (int j=0;j<_nWaves;j++){
			double upLimJ = _waves[j].GetMaxLimit();
			double lowLimJ= _waves[j].GetMinLimit();
			for (int bin =0; bin<_nBins;bin++){
				double mmm = (_binning[bin]+_binning[bin+1])/2;
				if (_errors[i][j][bin] == 0. and mmm > lowLimI and mmm > lowLimJ and mmm < upLimI and mmm < upLimJ){
					std::cout << "Inconsisntency found: Zero error at: _errors["<<i<<"]["<<j<<"]["<<bin<<"]"<<std::endl;
					nErr+=1;
				};
			};			
		};
	};
	if (nErr >0){
		return false;
	};
	return true;
};

std::vector<std::complex<double> > chi2::GetPoints(int i, int nPoints){
	std::vector<std::complex<double> > waveResult;
	double step = (_mMax - _mMin)/nPoints;
	for (int j=0;j<nPoints;j++){
		std::vector<double> pars = _waves[i].GetParameters();
		waveResult.push_back(_waves[i].Ampl<double>(_mMin + (j+0.5)*step,pars));
	};
	return waveResult;
};

std::vector<std::complex<double> > chi2::GetPointsComponent(int i, std::vector<int> comp, int nPoints){
	std::vector<bool> cpl= _waves[i].Couplings();
	std::vector<double> par_val;
	int nPars = _waves[i].GetNpars();
	for(int j=0;j<nPars; j++){ 							// Step 1: Store parameters
		par_val.push_back(_waves[i].GetParameter(j));
	};
	int parametrization = 0;
	for(int j=0;j<nPars;j++){							// Step 2: Set non used parameters to zero
		if(cpl[j]){
			bool zero = true;
			for(unsigned int k=0;k<comp.size();k++){
				if(comp[k] == parametrization){
					zero = false;
				};
			};
			if (zero){
				_waves[i].SetParameter(j,0.);
				_waves[i].SetParameter(j+1,0.);
			};
			j+=1;
			parametrization+=1;
		};
	};
	std::vector<std::complex<double> > result = GetPoints(i,nPoints);		// Step 3: Evaluate with the values
	for(int j=0;j<nPars;j++){							// Step 4: Set Parameters back to original value
		_waves[i].SetParameter(j,par_val[j]);
	};	
	return result;
};

void chi2::SetRandomCouplings(){
	_waves[0].SetRandomCouplings(true);
	for (int i=1; i<_nWaves;i++){
		_waves[i].SetRandomCouplings(false);
	};
};
void chi2::SetCouplingSize(){
	int expMin=-4;
	int expMax=8;
	double basis=10.;
	std::vector<std::vector<bool> > stats;
	for (int i=0;i<_nWaves;i++){
		stats.push_back(_waves[i].Couplings());
	};
	std::vector<double> pars;
	for (int i=0;i<_nWaves;i++){
		for (int j=0;j<_waves[i].GetNpars();j++){
			if (stats[i][j]){
				pars.push_back(0.);
			}else{
				pars.push_back(_waves[i].GetParameter(j));
			};
		};
	};
	int atPar=0;
	for (int i=0;i<_nWaves;i++){
		bool switched = false;
		if (_waves[i].GetRphi()){
			_waves[i].SetRphi(false);
			switched = true;
		};
		for (int j=0;j<_waves[i].GetNpars();j++){
			if (stats[i][j]){
				pars[atPar] = pow(basis,expMax);
				double best = Eval<double>(pars);
				double act;
				int bestExp = expMax;
				for (int exp = expMin; exp<expMax; exp++){
					pars[atPar] = pow(basis,exp);
					act = Eval<double>(pars);
//					std::cout << exp <<"::"<<act<<std::endl;
					if (act < best){
						best = act;
						bestExp = exp;
					};
				};
//				std::cout << i <<"::"<<j<<"::bestExp::"<<bestExp<<std::endl;
				_waves[i].SetUpperLim(j,pow(basis,bestExp+2));
				_waves[i].SetLowerLim(j,-pow(basis,bestExp+2));
				pars[atPar]=0.;
			};
			atPar+=1;
		};
		if (switched){
			_waves[i].SetRphi(true);
		};
	};
};


void chi2::SetAnchorWave(int i){
	_anchorWave=i;
};







