#include "chi2amp.h"
#include "waveset.h"
#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include <fstream>

#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON

chi2amp::chi2amp(const int nWaves): waveset(nWaves) ,_anchorWave(0) {
	_points = std::vector<std::vector<double> >(2*_nWaves-1,std::vector<double>());
	_coma   = std::vector<std::vector< std::vector<double> > >(2*_nWaves-1,std::vector<std::vector<double> >(2*_nWaves-1,std::vector<double>()));
};


template<typename xdouble> // Using Deltas(...) to calculate the chi2amp.
xdouble chi2amp::Eval(std::vector<xdouble> &par){
	xdouble chi2amp_tot = 0.;
	for (int bin=0; bin<_nBins; bin++){			// Loop over mass-bins
		std::vector<xdouble> deltas = Deltas<xdouble>(bin,par);
		for (int i=0; i<2*_nWaves-1;i++){
			for (int j=0;j<2*_nWaves-1;j++){
				chi2amp_tot+=deltas[i]*_coma[i][j][bin]*deltas[j];
				if (deltas[i]*_coma[i][j][bin]*deltas[j]>1000.){
//					std::cout<< "----------------------------------------" << std::endl;
//					std::cout<< "delta[" << i << "][" << j << "][" << bin << "] =" << deltas[i]*_coma[i][j][bin]*deltas[j] <<std::endl;
//					std::cout<<"_coma["<<i<<"]["<<j<<"]["<<bin<<"] = "<<_coma[i][j][bin]<<std::endl;
//					std::cout<<"deltas["<<i<<"] ="<<deltas[i]<<std::endl;
//					std::cout<<"deltas["<<j<<"] ="<<deltas[j]<<std::endl;

				};
			};
		};
	};
	if(_nCalls%1000 == 0){
		std::cout << "Function call #"<<_nCalls<<": "<<chi2amp_tot<<std::endl;
	};
	return chi2amp_tot;	
};
template double chi2amp::Eval(std::vector<double> &par);

template<typename xdouble> // Gives the differences between points and fit-curve.
std::vector<xdouble> chi2amp::Deltas(int bin, std::vector<xdouble> &par){
	double min = _binning[bin];
	double max = _binning[bin+1];
	double m   = (min+max)/2;
	std::vector<std::complex<xdouble> >amps = Amps(bin,par);
	std::complex<xdouble> anchorNorm = std::complex<xdouble>(std::norm(amps[_anchorWave]),0.);
	std::vector<xdouble> deltas(2*_nWaves-1,0.);
	double aMax = _waves[_anchorWave].GetMaxLimit();
	double aMin = _waves[_anchorWave].GetMinLimit();
	if (m > aMin and m< aMax){
		std::complex<xdouble> phaseFactor = amps[_anchorWave] * pow(anchorNorm,-.5);
		for (int i=0;i<_nWaves-1;i++){
			amps[i]/=phaseFactor;
		};
		deltas[0] = std::real(amps[_anchorWave])-_points[0][bin]; //Assumes _anchorWave is _waves[0]. Does not work otherwise
		for (int i=1;i<_nWaves;i++){
			double iMax = _waves[i].GetMaxLimit();
			double iMin = _waves[i].GetMinLimit();
			if (iMin < m and iMax > m){
				deltas[2*i-1] = std::real(amps[i]) - _points[2*i-1][bin];
				deltas[2*i]   = std::imag(amps[i]) - _points[2*i  ][bin];
			};
		};
	};
	return deltas; // Exclude bins with delta == 0.
};
template std::vector<double> chi2amp::Deltas(int bin, std::vector<double> &par);

template<typename xdouble> //Take a vector of only relreased parameters as argument. Is equal to Eval, if all parameters are released.
xdouble chi2amp::EvalRel(std::vector<xdouble> &par){
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
template double chi2amp::EvalRel(std::vector<double> &par);

template<typename xdouble>// Equal to EvalRel(GetParFromNorm(par))
xdouble chi2amp::EvalNorm(std::vector<xdouble> &par){
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
template double chi2amp::EvalNorm(std::vector<double> &par);



#ifdef ADOL_ON
adouble chi2amp::EvalAdouble(std::vector<adouble> par){
	return Eval(par);
};
#endif//ADOL_ON

std::string chi2amp::ClassName(){
	return "chi2amp";
};


void chi2amp::SetData(int i, std::vector<double> data){
	_points[i] = data;
	// i in [0,2*_nWaves-2] : 
	// _points[0] = norm(T_0)
	// _points[1] = real(T_1) (rotated that T_0 is real)
	// _points[2] = imag(T_1) (rotated that T_0 is real)
	// _points[3] = real(T_2) (rotated that T_0 is real)
	// .
	// .
	// _points[2_n-2] = imag(T_{_nWaves-1}) (rotated that T_0 is real)
};
void chi2amp::SetComa(int i, int j, std::vector<double> error){
//	if (error.size() != _nBins){
//		std::cout << "Warning: Size of new errors does not match number of bins." << std::endl;
//	};
	_coma[i][j] = error;
};

std::vector<double> chi2amp::GetData(int i){
	return _points[i];
};
std::vector<double> chi2amp::GetComa(int i, int j){
	return _coma[i][j];
};

bool chi2amp::CheckConsistency(){
	int nErr=0;
	if (not waveset::CheckConsistency()){
		nErr+=1;
	};
	for (int i=0;i<_nWaves;i++){
		if (_points[i].size() != _nBins){
			std::cout << "Inconsistency found: Data points in _points["<<i<<"] does not match _nBins" << std::endl;
			nErr+=1;
		};
		for (int j=0;j<_nWaves;j++){
			if (_coma[i][j].size() != _nBins){
				std::cout << "Inconsistency found: Errors in _coma["<<i<<"]["<<j<<"] does not match _nBins" << std::endl;
				nErr+=1;
			};

		};
	};
	double minA =_waves[0].GetMinLimit();
	double maxA =_waves[0].GetMaxLimit();
	for (int i=1;i<_nWaves;i++){
		double minn =_waves[i].GetMinLimit();
		double maxx =_waves[i].GetMaxLimit();
		if (minn < minA){
			std::cout << "Inconsistency found: _waves["<<i<<"].GetMinLimit() smaller than lower limit of the anchor wave" << std::endl;
			nErr+=1;
		};
		if (maxx > maxA){
			std::cout << "Inconsistency found: _waves["<<i<<"].GetMaxLimit() bigger than upper limit of the anchor wave" << std::endl;
			nErr+=1;
		};
	};
	if (nErr >0){
		return false;
	};
	return true;
};

std::vector<std::complex<double> > chi2amp::GetPoints(int i, int nPoints){
	std::vector<std::complex<double> > waveResult;
	double step = (_mMax - _mMin)/nPoints;
	for (int j=0;j<nPoints;j++){
		std::vector<double> pars = _waves[i].GetParameters();
		waveResult.push_back(_waves[i].Ampl<double>(_mMin + (j+0.5)*step,pars));
	};
	return waveResult;
};

std::vector<std::complex<double> > chi2amp::GetPointsComponent(int i, std::vector<int> comp, int nPoints){
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

void chi2amp::SetRandomCouplings(){
	_waves[0].SetRandomCouplings(true);
	for (int i=1; i<_nWaves;i++){
		_waves[i].SetRandomCouplings(false);
	};
};

void chi2amp::SetCouplingSize(){
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









