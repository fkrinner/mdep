#include "waveset.h"
#include "wave.h"
#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include <fstream>

#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON


waveset::waveset(const int nWaves): _nWaves(nWaves),_nCalls(0),_nPars(0),_mMin(1.),_mMax(2.),_binWidth(0.02),_isSampleAmps(false){
	for (int i=0; i < _nWaves;i++){
		_waves.push_back(wave());
		_waveNpars.push_back(0);
	};
	_nBins = (_mMax-_mMin)/_binWidth;
};

template<typename xdouble>
std::vector<std::complex<xdouble> > waveset::Amps(int bin, std::vector<xdouble> &par){
	std::vector<std::complex<xdouble> > amps;
	double min = _binning[bin];
	double max = _binning[bin+1];
	double m   = (min+max)/2;
	if (not _isSampleAmps){
		for (int wave=0;wave<_nWaves; wave++){
			double iMax = _waves[wave].GetMaxLimit();
			double iMin = _waves[wave].GetMinLimit();
			if(m < iMax and m > iMin){
				int nSkip=0;
				for (int j=0; j<wave; j++){
					nSkip+=_waveNpars[j];
				};
				std::vector<xdouble> actPars = std::vector<xdouble>();
				for (int j=nSkip; j< nSkip+_waveNpars[wave]; j++){
					actPars.push_back(par[j]);
				};
				amps.push_back(WaveAmp<xdouble>(wave,min,max,actPars));
			} else {
				amps.push_back(std::complex<xdouble>(0.,0.));
			};
		};
		if (bin==0){
			_nCalls+=1;
		};
	}else{
		for (int i=0;i<_nWaves;i++){
			amps.push_back(_sampleAmps[bin][i]);
		};
	};
	return amps;
};
template std::vector<std::complex<double> > waveset::Amps(int bin, std::vector<double> &par);


template <typename xdouble>
std::complex<xdouble> waveset::WaveAmp(int i, double min,double max, std::vector<xdouble> &par) const{
	double m = (max+min)/2; 		// To be able to implement some integration in the waves
	return _waves[i].Ampl<xdouble>(m,par);
};
template std::complex<double> waveset::WaveAmp(int i, double min,double max, std::vector<double> &par) const;

#ifdef ADOL_ON
template std::vector<std::complex<adouble> > waveset::Amps(int bin, std::vector<adouble> &par);

template std::complex<adouble> waveset::WaveAmp(int i, double min,double max, std::vector<adouble> &par) const;

adouble waveset::EvalAdouble(std::vector<adouble> par){
	std::cout <<"Warning: waveset.cxx: EvalAdouble(...) not oveloaded, Gradient(...) not valid." << std::endl;
	return 0.;
};

std::vector<double> waveset::Gradient(std::vector<double> par){
	int nTape = 0;
	const int parSize = par.size();
	trace_on(nTape);
	double x[parSize];
	vector<adouble> ax;
	adouble a;
	for (int i=0; i < parSize; i++){
		ax.push_back(a);
		x[i] = par[i];
		ax[i] <<= x[i];
	};
	double y;
	adouble ay;
	ay = EvalAdouble(ax);
	ay >>= y;
	trace_off();
	double grad[parSize];
	gradient(nTape,parSize,x,grad);	
	vector<double> gradient;
	for (int i=0;i<parSize;i++){
		gradient.push_back(grad[i]);
	};
	return gradient;
};
#endif//ADOL_ON

void waveset::SetWaveModel(int i, std::vector<int> model){
	_waves[i].SetModel(model);
	_waveNpars[i]=_waves[i].GetNpars();
	_nPars=0;
	for (int j=0; j<_nWaves;j++){
		_nPars+=_waveNpars[j];
	};
};

std::vector<int> waveset::GetWaveModel(int i){
	return _waves[i].GetModel();
};

void waveset::SetRelParameters(std::vector<double> par){
	int actPar=0;
	for (int i=0; i<_nWaves; i++){
		std::vector<int> actRels = _waves[i].GetRelParNumbers();
		for (unsigned int j=0; j< actRels.size(); j++){
			_waves[i].SetParameter(actRels[j],par[actPar]);
			actPar+=1;
		};
	};
};
void waveset::SetRelNormParameters(std::vector<double> par){
	int actPar=0;
	for (int i=0; i<_nWaves; i++){
		std::vector<int> actRels = _waves[i].GetRelParNumbers();
		for (unsigned int j=0; j< actRels.size(); j++){
			double paramet = (1-par[actPar])*_waves[i].GetLowerLim(actRels[j]) + par[actPar]*_waves[i].GetUpperLim(actRels[j]);
			_waves[i].SetParameter(actRels[j],paramet);
			actPar+=1;
		};
	};
};
std::vector<double> waveset::GetRelParameters(){
	std::vector<double> rels;
	for (int i=0; i<_nWaves; i++){
		std::vector<int> actRels = _waves[i].GetRelParNumbers();
		for (unsigned int j=0; j< actRels.size(); j++){
			rels.push_back(_waves[i].GetParameter(actRels[j]));
		};
	};
	return rels;
};
std::vector<double> waveset::GetRelNormParameters(){
	std::vector<double> rels;
	for (int i=0; i<_nWaves; i++){
		std::vector<int> actRels = _waves[i].GetRelParNumbers();
		for (unsigned int j=0; j< actRels.size(); j++){
			double up = _waves[i].GetUpperLim(actRels[j]);
			double low= _waves[i].GetLowerLim(actRels[j]);
			rels.push_back((_waves[i].GetParameter(actRels[j])-low)/(up-low));
		};
	};
	return rels;
};

void waveset::SetWaveParameter(int i, int n, double param){
	_waves[i].SetParameter(n,param);
};
void waveset::SetWaveParameters(int i, std::vector<double> pars){
	_waves[i].SetParameters(pars);
};

double waveset::GetWaveParameter(int i, int n){
	return _waves[i].GetParameter(n);
};
std::vector<double> waveset::GetParameters(){
	std::vector<double> pars;
	for(int i=0;i<_nWaves;i++){
		int nPar = _waves[i].GetNpars();
		for (int j=0;j<nPar;j++){
			pars.push_back(_waves[i].GetParameter(j));
		};
	};
	return pars;
};
std::vector<double> waveset::GetWaveParameters(int i){
	return _waves[i].GetParameters();
};

//////////////////////////////// NormParameters are in the range from 0. to 1. and correspond to the range upperLim to lowerLim for each parameter:
//   par = (1-normPar)*lowerLim + normPar*upperLim
//   normPar = (par - lowerLim)/(upperLim - lowerLim)


void waveset::SetWaveNormParameter(int i, int n, double param){
	_waves[i].SetParameter(n, (1-param)*_waves[i].GetLowerLim(n) + param*_waves[i].GetUpperLim(n));
};
double waveset::GetWaveNormParameter(int i, int n){
	return (_waves[i].GetParameter(n)-_waves[i].GetLowerLim(n))/(_waves[i].GetUpperLim(n)-_waves[i].GetLowerLim(n));
};

double waveset::GetNormFromPar(int i, int n, double par){
	return (par - _waves[i].GetLowerLim(n))/(_waves[i].GetUpperLim(n)- _waves[i].GetLowerLim(n));
};

double waveset::GetParFromNorm(int i, int n, double norm){
	return (1. - norm)*_waves[i].GetLowerLim(n) + norm*_waves[i].GetUpperLim(n);
};

std::string waveset::ClassName(){
	return "waveset";
};

void waveset::SetWaveParameterName(int i, int n, std::string name){
	_waves[i].SetParName(n,name);
};

std::string waveset::GetWaveParameterName(int i, int n){
	return _waves[i].GetParName(n);
};
std::vector<std::string> waveset::GetParNames(){
	std::vector<std::string> names;
	for(int i=0;i<_nWaves;i++){
		int nPar =_waves[i].GetNpars();
		for (int j=0;j<nPar;j++){
			names.push_back(_waves[i].GetParName(j));
		};
	};
	return names;
};
void waveset::SetWaveLimits(int i, double lower, double upper){
	_waves[i].SetLimits(lower,upper);
};
double waveset::GetWaveMaxLimit(int i){
	return _waves[i].GetMaxLimit();
};
double waveset::GetWaveMinLimit(int i){
	return _waves[i].GetMinLimit();
};

int waveset::GetNwaves(){
	return _nWaves;
};
int waveset::GetNpars(){
	return _nPars;
};
int waveset::GetWaveNpars(int i){
	return _waveNpars[i];
};

void waveset::PrintParameters(){
	for (int i=0;i<_nWaves;i++){
		std::cout << "-------------------------------------------------------" << std::endl;
		std::cout << "_waves["<< i <<"]: " << _waves[i].GetName() << " r/phi: " <<_waves[i].GetRphi() << std::endl;
		std::cout << "-------------------------------------------------------" << std::endl;
		_waves[i].PrintParameters();
		std::cout <<std::endl;
	};

};
void waveset::PrintWaveParameters(int i){
	_waves[i].PrintParameters();
};
void waveset::PrintWaveModel(int i){
	_waves[i].PrintModelInfo();
};

void waveset::SetWaveName(int i,std::string name){
	_waves[i].SetName(name);
};

std::string waveset::GetWaveName(int i){
	return _waves[i].GetName();
};

void waveset::SetMmax(double m){
	_mMax=m;
	_nBins = (_mMax-_mMin)/_binWidth+0.001;
	_binning = std::vector<double>();
	_binning.push_back(_mMin);
	for (int i=0; i<_nBins; i++){
		_binning.push_back(_mMin + (i+1)*_binWidth);
	};
//	for (int i=0;i<_nWaves;i++){
//		for (int j=0;j<_nWaves;j++){
//			if (_points[i][j].size() != _nBins){
//				std::cout << "Warning: Number of bins does not match data-size for _points[" << i << "][" << j << "]" << std::endl;
//			};
//			if (_errors[i][j].size() != _nBins){
//				std::cout << "Warning: Number of errors does not match data-size for _errors[" << i << "][" << j << "]" << std::endl;
//			};
//		};
//	};
};
void waveset::SetMmin(double m){
	_mMin=m;
	_nBins = (_mMax-_mMin)/_binWidth+0.001;
	_binning = std::vector<double>();
	_binning.push_back(_mMin);
	for (int i=0; i<_nBins; i++){
		_binning.push_back(_mMin + (i+1)*_binWidth);
	};
//	for (int i=0;i<_nWaves;i++){
//		for (int j=0;j<_nWaves;j++){
//			if (_points[i][j].size() != _nBins){
//				std::cout << "Warning: Number of bins does not match data-size for _points[" << i << "][" << j << "]" << std::endl;
//			};
//			if (_errors[i][j].size() != _nBins){
//				std::cout << "Warning: Number of errors does not match data-size for _errors[" << i << "][" << j << "]" << std::endl;
//			};
//		};
//	};
};
void waveset::SetBinWidth(double m){
	_binWidth=m;
	_nBins = (_mMax-_mMin)/_binWidth+0.001;
	_binning = std::vector<double>();
	_binning.push_back(_mMin);
	for (int i=0; i<_nBins; i++){
		_binning.push_back(_mMin + (i+1)*_binWidth);
	};
//	for (int i=0;i<_nWaves;i++){
//		for (int j=0;j<_nWaves;j++){
//			if (_points[i][j].size() != _nBins){
//				std::cout << "Warning: Number of bins does not match data-size for _points[" << i << "][" << j << "]" << std::endl;
//			};
//			if (_errors[i][j].size() != _nBins){
//				std::cout << "Warning: Number of errors does not match data-size for _errors[" << i << "][" << j << "]" << std::endl;
//			};
//		};
//	};
};
double waveset::GetMmax(){
	return _mMax;
};
double waveset::GetMmin(){
	return _mMin;
};
double waveset::GetBinWidth(){
	return _binWidth;
};
int waveset::GetNbins(){
	return _nBins;
};
std::vector<bool> waveset::GetParStat(){
	std::vector<bool> stats;
	for(int i=0;i<_nWaves;i++){
		int nPar=_waves[i].GetNpars();
		for (int j=0;j<nPar;j++){
			stats.push_back(_waves[i].GetParStat(j));
		};
	};
	return stats;
};

bool waveset::GetWaveParStat(int i, int n){
	return _waves[i].GetParStat(n);
};
void waveset::FixWaveParameter(int i, int n){
	_waves[i].FixParameter(n);
};
void waveset::RelWaveParameter(int i, int n){
	_waves[i].RelParameter(n);
};
int waveset::GetNparsReleased(){
	int nRel=0;
	for (int i=0;i<_nWaves;i++){
		nRel+=_waves[i].GetNparsReleased();
	};
	return nRel;
};
int waveset::GetWaveNparsReleased(int i){
	return _waves[i].GetNparsReleased();
};
std::vector<std::string> waveset::GetRelParNames(){
	std::vector<std::string> names;
	for (int i=0;i<_nWaves;i++){
		std::vector<std::string> act = _waves[i].GetRelParNames();
		names.insert(names.end(),act.begin(),act.end());
	};
	return names;
};
std::vector<int> waveset::GetWaveRelParNumbers(int i){
	return _waves[i].GetRelParNumbers();
};

std::vector<double> waveset::GetWaveUpperLims(int i){
	return _waves[i].GetUpperLims();
};
std::vector<double> waveset::GetWaveLowerLims(int i){
	return _waves[i].GetLowerLims();
};
double waveset::GetWaveUpperLim(int i, int n){
	return _waves[i].GetUpperLim(n);
};	
double waveset::GetWaveLowerLim(int i, int n){
	return _waves[i].GetLowerLim(n);
};
void waveset::SetWaveUpperLims(int i, std::vector<double> lim){
	_waves[i].SetUpperLims(lim);
};
void waveset::SetWaveLowerLims(int i, std::vector<double> lim){
	_waves[i].SetLowerLims(lim);
};
void waveset::SetWaveUpperLim(int i, int n, double lim){
	_waves[i].SetUpperLim(n,lim);
};
void waveset::SetWaveLowerLim(int i, int n, double lim){
	_waves[i].SetLowerLim(n,lim);
};

void waveset::SetWaveSpin(int i,int L){
	_waves[i].SetSpin(L);
};
int waveset::GetWaveSpin(int i){
	return _waves[i].GetSpin();
};
void waveset::SetWaveIsobarMass(int i, double m){
	_waves[i].SetIsobarMass(m);
};
double waveset::GetWaveIsobarMass(int i){
	return _waves[i].GetIsobarMass();
};

void waveset::SetWavePhaseSpace(int i, std::vector<int> ps){
	_waves[i].SetPhaseSpace(ps);
};
std::vector<int> waveset::GetWavePhaseSpace(int i){
	return _waves[i].GetPhaseSpace();
};
void waveset::SetPhaseSpace(std::vector<int> ps){
	for (int i=0; i<_nWaves;i++){
		_waves[i].SetPhaseSpace(ps);
	};
};

void waveset::SetBinning(std::vector<double> binning){
	_binning = binning;
	_nBins = binning.size()-1;
	_mMin = binning[0];
	_mMax = binning[_nBins];
	_binWidth = (_mMax-_mMin)/_nBins;
};

std::vector<double> waveset::GetBinning(){
	return _binning;
};

void waveset::SetWaveRphi(int i, bool flag){
	_waves[i].SetRphi(flag);
};
bool waveset::GetWaveRphi(int i){
	return _waves[i].GetRphi();
};
void waveset::SetRphi(bool flag){
	for (int i=0;i<_nWaves;i++){
		_waves[i].SetRphi(flag);
	};
};


int waveset::GetNcalls(){
	return _nCalls;
};

bool waveset::CheckConsistency(){
	int nErr=0;
	if (_nBins != (int)((_mMax-_mMin)/_binWidth+0.1)){
		std::cout << "Inconsistecy found: _nBins does not compute with _mMin, _mMax, _binWidth" << std::endl;
		nErr+=1;
	};
	for (int i=0;i<_nWaves; i++){
		if (_waves[i].GetModel().size() ==0){
			std::cout << "Inconsistecy found: Model for _wave["<<i<<"] not set"<<std::endl;
			nErr+=1;
		};
	};
	for (int i=0;i<_nWaves; i++){
		if(_waves[i].GetUpperLims().size() != _waves[i].GetParameters().size()){
			std::cout << "Inconsistency found: _waves["<<i<<"]: _parameters.size() not equal to _upperLim.size()." << std::endl;
			nErr+=1;
		};
		if(_waves[i].GetLowerLims().size() != _waves[i].GetParameters().size()){
			std::cout << "Inconsistency found: _waves["<<i<<"]: _parameters.size() not equal to _lowerLim.size()." << std::endl;
			nErr+=1;
		};
	};
	if (_binning.size() != _nBins+1){
		std::cout << "Inconsistency found: _nBins does not coincide with _binning.size()" << std::endl;
		nErr+=1;
	}else{
		if (_binning[0] != _mMin){
			std::cout << "Inconsisntency found: _binning[0] does not coincide with _mMin" << std::endl;
			nErr+=1;
		};
		if (_binning[_nBins] != _mMax){
			std::cout << "Inconsistency found: _binning[_nBins] does not coincide with _mMax" << std::endl;
			nErr+=1;
		};
	};
	if (_isSampleAmps){
		std::cout << "Warning: Sample mode active. Amplitudes will not be evaluated, but taken from '_sampleAmps[...]'"<<std::endl;
	};
	if (nErr >0){
		return false;
	};
	return true;
};

void waveset::IsSampleAmps(bool flag){
	if (flag and _sampleAmps.size() != _nBins){
		_sampleAmps = std::vector<std::vector<std::complex<double> > >();
		for (int i=0;i<_nBins;i++){
			_sampleAmps.push_back(std::vector<std::complex<double> >());
		};	
	};
	_isSampleAmps = flag;
};
void waveset::SetSampleAmps(int bin, std::vector<std::complex<double> > in){
	_sampleAmps[bin] = in;
};






