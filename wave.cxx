#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h> 

#include "wave.h"
#include "breitWigners.h"
#include "phaseSpace.h"

#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON




///////////////////////// CONSTRUCTORS OF WAVE

wave::wave():_name("unnamed"), _rPhi(false){};
	

wave::wave(std::vector<int> model):_model(model),_minLim(0.), _maxLim(0.), _name("unnamed"), _L(0), _rPhi(false){  // KREATOR
	for (unsigned int i=0;i<_model.size();i++){
		_nPars.push_back(getNpars(_model[i]));
	};
	_nPar= 2*_model.size();
	for (unsigned int i=0; i < _model.size();i++){ 
		_nPar+=_nPars[i];	
	};
	for (int i =0; i< _nPar;i++){
		_parameters.push_back(0.);
	};
	std::string prefix="unnamed";
	std::stringstream sufstream;
	for (int i=0; i < _nPar;i++){ 
		sufstream << i;
		_parameterNames.push_back(prefix+sufstream.str());
		sufstream.str("");
	};
	_release_parameter.clear();	
	for (unsigned int i=0;i<_model.size();i++){
		_release_parameter.push_back(true);
		_release_parameter.push_back(true);
		for (int j=0;j<_nPars[i];j++){
			_release_parameter.push_back(false);
		};
	};	
};

/////////////// NONTRIVIAL METHODS OF WAVE

template<typename xdouble>
std::complex<xdouble> wave::Ampl(double m,std::vector<xdouble> &x) const{
	if (m >= _minLim and m <= _maxLim){
		std::complex<xdouble> amp = std::complex<xdouble>(0,0);
		std::vector<xdouble> param;	// Will contain the parameters for a single parametrization
		int skipParam;
		for (unsigned int i=0;i<_model.size();i++){
			param.clear();
			skipParam = 2;	// Skip the coupling of the first amplitude
			for (unsigned int j=0; j < i; j++){	// 
				skipParam+=_nPars[j];	// Skip the parameters belonging to previous parametrizations
				skipParam+=2;		// Skip also the couplings of the provious parametrizations
			};
	//				std::cout << "Skipping "<< skipParam<< " parameters for parametrization #" << i+1 <<std::endl; 
			for (int j=skipParam;j<skipParam+_nPars[i];j++){
	//					std::cout << "Setting par "<<j<<"  "<<x[j]<<std::endl;
				param.push_back(x[j]);	// Add parameters for current parametrizations to the std::vector
			};
			if (_rPhi){
				amp+=std::complex<xdouble>(x[skipParam-2],0.)*std::complex<xdouble>(cos(x[skipParam-1]),sin(x[skipParam-1]))*bw<xdouble>(m, param, _model[i],_L);
//				std::cout<<_name<<":: " << std::complex<xdouble>(x[skipParam-2],0.)*std::complex<xdouble>(cos(x[skipParam-1]),sin(x[skipParam-1])) << std::endl;
			}else{
				amp+=std::complex<xdouble>(x[skipParam-2],x[skipParam-1])*bw<xdouble>(m, param, _model[i],_L); 	// Evaluate --  Add to the amplitude. (The skipped parameters from before are here used as couplings)
			};
	//				std::cout<< "amp = " << amp <<" = " <<  std::complex<double>(x[2*i],x[2*i+1]) << " * "<< bw<double>(m, param, model[i])<<std::endl;
		};
		if (_phaseSpace.size() > 0){
			double factor = 1.;
			for (unsigned int ii = 0; ii < _phaseSpace.size();ii++){
				factor*=phaseSpace(m,_phaseSpace[ii],_L,_isobar_mass);
			};
			amp*=factor;
		};
		return amp;
	};
	return std::complex<xdouble>(0.,0.);
};
template std::complex<double> wave::Ampl(double m,std::vector<double> &x) const;

template<typename xdouble>
std::complex<xdouble> wave::AmplRel(double m,std::vector<xdouble> &x) const{
	int irel=0;
	std::vector<xdouble> xrel;
	for (int i=0;i<_nPar;i++){
		if (_release_parameter[i]){
			xrel.push_back(x[irel]);
			irel+=1;
		}else{
			xrel.push_back(_parameters[i]);
		};
	};
	return Ampl<xdouble>(m,xrel);
};
template std::complex<double> wave::AmplRel(double m,std::vector<double> &x) const;

#ifdef ADOL_ON
template std::complex<adouble> wave::Ampl(double m,std::vector<adouble> &x) const;
template std::complex<adouble> wave::AmplRel(double m,std::vector<adouble> &x) const;
#endif//ADOL_ON

double wave::RealPart(std::vector<double> m,std::vector<double> param){
			std::complex<double> amplitude = Ampl(m[0],param);
			return std::real(amplitude);
		};

double wave::ImagPart(std::vector<double> m,std::vector<double> param){
			std::complex<double> amplitude = Ampl(m[0],param);
			return std::imag(amplitude);
		};

void wave::PrintModelInfo(){
	for (unsigned int i=0;i<_model.size();i++){
		std::cout << "Function #" << _model[i] << " using " << _nPars[i] << " parameters." << std::endl;
	};
};

void wave::PrintParameters(){
	for(unsigned int i=0; i<_parameters.size();i++){
		std::cout << "_parameters["<<i<<"] = "<< _parameters[i]<<  " with Name: '"<< _parameterNames[i]<<"'";
		if (_release_parameter[i]){
			std::cout << " free" << std::endl;
		}else{
			std::cout << " fixed" << std::endl;
		};
	};
};


//////////////////////// SETTERS AND GETTERS FOR WAVE

std::string wave::GetName(){
	return _name;
};

void wave::SetName(std::string name){
	_name=name;
};

std::string wave::ClassName(){
	return "wave";
};

double wave::GetParameter(int i){
	return _parameters[i];
};

void wave::SetParameter(int i, double par){
	if (i < _nPar){
		_parameters[i]=par;
	}else{
		std::cerr << "wave.cxx: Error: Specified parameter number does not exist: _parameters[" << i << "]" << std::endl;
	};
};

std::vector<double> wave::GetParameters(){
	return _parameters;
};
void wave::SetParameters(std::vector<double> pars){
	_parameters = pars;
};


int wave::GetNpars(){
	return _nPar;
};

void wave::SetLimits(double min, double max){
	_minLim=min;
	_maxLim=max;
};

double wave::GetMaxLimit(){
	return _maxLim;
};

double wave::GetMinLimit(){
	return _minLim;
};

std::string wave::GetParName(int i){
	return _parameterNames[i];
//			return defaultNames[i];
};

void wave::SetParName(int i,std::string name){
	if (i<_nPar) {
		_parameterNames[i]=name;	
	}else{
		std::cerr <<  "wave.cxx: Error: Specified parameter number does not exist: _parameterNames[" << i << "]" << std::endl;
	};
};

void wave::SetModel(std::vector<int> model){
	_model=model;
	_nPars.clear();
	for (unsigned int i=0;i<_model.size();i++){
		_nPars.push_back(getNpars(_model[i]));
	};
	_nPar= 2*_model.size();
	for (unsigned int i=0; i < _model.size();i++){ 
		_nPar+=_nPars[i];
	};
	_parameters.clear();
	for (int i =0; i< _nPar;i++){
		_parameters.push_back(0.);
	};	
	std::string prefix="unnamed";
	std::stringstream sufstream;
	for (int i=0; i < _nPar;i++){ 
		sufstream << i+1;
		_parameterNames.push_back(prefix+sufstream.str());
		sufstream.str("");
	};
	_release_parameter.clear();
	for (unsigned int i=0;i<_model.size();i++){
		_release_parameter.push_back(true);
		_release_parameter.push_back(true);
		for (int j=0;j<_nPars[i];j++){
			_release_parameter.push_back(false);
		};
	};			
};		

std::vector<int> wave::GetModel(){
	return _model;
};

bool wave::GetParStat(int i){
	return _release_parameter[i];
};
void wave::FixParameter(int i){
	_release_parameter[i]=false;
};
void wave::RelParameter(int i){
	_release_parameter[i]=true;
};
int wave::GetNparsReleased(){
	int nRel=0;
	if (_nPar !=  _release_parameter.size()){
		std::cerr << "wave.cxx: Error: Parameter status not equal to number of parameters" << std::endl;
	};
	for (unsigned int i=0; i<_nPar;i++){
		if (_release_parameter[i]) {
			nRel+=1;
		};
	};
	return nRel;
};
std::vector<std::string> wave::GetRelParNames(){
	std::vector<std::string> names;
	for (int i=0;i<_nPar;i++){
		if (_release_parameter[i]){
			names.push_back(_parameterNames[i]);
		};
	};	
	return names;
};
std::vector<double> wave::GetRelParValues(){
	std::vector<double> pars;
	for (int i=0;i<_nPar;i++){
		if (_release_parameter[i]){
			pars.push_back(_parameters[i]);
		};
	};	
	return pars;
};
std::vector<int> wave::GetRelParNumbers(){
	std::vector<int> pars;
	for (int i=0;i<_nPar;i++){
		if (_release_parameter[i]){
			pars.push_back(i);
		};
	};	
	return pars;
};

std::vector<double> wave::GetUpperLims(){
	return _upperLim;
};
std::vector<double> wave::GetLowerLims(){
	return _lowerLim;
};
double wave::GetUpperLim(int n){
	return _upperLim[n];
};
double wave::GetLowerLim(int n){
	return _lowerLim[n];
};
void wave::SetUpperLims(std::vector<double> lim){
	_upperLim = lim;
};
void wave::SetLowerLims(std::vector<double> lim){
	_lowerLim = lim;
};
void wave::SetUpperLim(int n, double lim){
	_upperLim[n]=lim;
};
void wave::SetLowerLim(int n, double lim){
	_lowerLim[n]=lim;
};

void wave::SetSpin(int L){
	_L=L;
};
int wave::GetSpin(){
	return _L;
};
void wave::SetIsobarMass(double m){
	_isobar_mass = m;
};
double wave::GetIsobarMass(){
	return _isobar_mass;
};

void wave::SetPhaseSpace(std::vector<int> ps){
	_phaseSpace=ps;
};
std::vector<int> wave::GetPhaseSpace(){
	return _phaseSpace;
};

void wave::SetRandomCouplings(bool real){
	std::vector<bool> toSet;
	for (unsigned int i =0; i< _nPars.size(); i++){
		toSet.push_back(true);
		toSet.push_back(true);
		for (int j=0; j<_nPars[i];j++){
			toSet.push_back(false);
		};
	};
	for (int i=0; i<_nPar;i++){
		if (toSet[i]){
			double upper = _upperLim[i];
			double lower = _lowerLim[i];
			double factor = double(rand())/RAND_MAX;
			SetParameter(i, lower * (1-factor) + upper * factor);
		};
	};
	if (real and _nPar > 1){	// Works for re/im and r/phi im = 0. <= phi = 0.
		SetParameter(1,0.);
	};
};
std::vector<bool> wave::Couplings(){
	std::vector<bool> toSet;
	for (unsigned int i =0; i< _nPars.size(); i++){
		toSet.push_back(true);
		toSet.push_back(true);
		for (int j=0; j<_nPars[i];j++){
			toSet.push_back(false);
		};
	};
	return toSet;
};

void wave::SetRphi(bool flag){
	if (flag != _rPhi){
		_rPhi = flag;
		std::vector<bool> cpls = Couplings();
		for (int i=0; i< _nPar;i++){
			if (cpls[i]){	
				double first  = _parameters[i];
				double secon  = _parameters[i+1];	
				double uFirst = _upperLim[i];
				double lFirst = _lowerLim[i];
				double uSecon = _upperLim[i+1];
				double lSecon = _lowerLim[i+1];	
				if (flag){ // Convert from re/im -> r/phi
					_parameters[i] = pow(first*first+secon*secon,0.5);
					_parameters[i+1]= atan2(secon,first);
					_upperLim[i]  = std::max(abs(uFirst),abs(uSecon));
					_lowerLim[i]  = 0.;
					_upperLim[i+1]= PI;
					_lowerLim[i+1]=-PI;					
				}else{     // Convert from r/phi -> re/im
					_parameters[i] = first *cos(secon);
					_parameters[i+1]= first*sin(secon);
					_upperLim[i]  = uFirst;	
					_lowerLim[i]  =-uFirst;
					_upperLim[i+1]= uFirst;
					_lowerLim[i+1]=-uFirst;
				};
			i+=1;	
			};
		};
	};
};
bool wave::GetRphi(){
	return _rPhi;
};


