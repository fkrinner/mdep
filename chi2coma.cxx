#include "chi2coma.h"
#include <iostream>

#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON


chi2coma::chi2coma(const int nWaves): chi2(nWaves){
	std::vector<double> empty;
	std::vector<std::vector<double> > line;
	for (int i=0;i<_nWaves*_nWaves;i++){
		line.push_back(empty);
	};
	for (int i=0;i<_nWaves*_nWaves;i++){
		_coma.push_back(line);
	};
};


template<typename xdouble>
xdouble chi2coma::EvalComa(std::vector<xdouble> &par){
	xdouble chi2_tot=0.;
	for(int bin = 0.;bin<_nBins;bin++){
		std::vector<xdouble> deltas = Deltas<xdouble>(bin,par);
		for(int i=0;i<_nWaves*_nWaves;i++){
			for(int j=0;j<_nWaves*_nWaves;j++){
				chi2_tot+= deltas[i]*deltas[j]*_coma[i][j][bin];
//				if(chi2_tot != chi2_tot and not getting_unreal){
//					std::cout << "deltas["<<i<<"]: "<<deltas[i] << std::endl;
//					std::cout << "deltas["<<j<<"]: "<<deltas[j] << std::endl;
//					std::cout << "_coma["<<i<<"]["<<j<<"]["<<bin<<"]: "<<_coma[i][j][bin]<<std::endl;
//					std::cout << "del*coma*del: " << deltas[i]*deltas[j]*_coma[i][j][bin]<<std::endl;
//				};
			};		
		};
	};
	if(_nCalls%1000 ==0){
		std::cout << "Function call #"<<_nCalls<<": "<<chi2_tot<<std::endl;
	};
	return chi2_tot;
};
template double chi2coma::EvalComa(std::vector<double> &par);

#ifdef ADOL_ON
adouble chi2coma::EvalAdouble(std::vector<adouble> par){
	return EvalComa<adouble>(par);
};
#endif//ADOL_ON

void chi2coma::SetComa(int i1, int j1, int i2, int j2, std::vector<double> coma){
	if (i1 == i2 and j1 == j2){
		std::vector<double> errors;
		for (unsigned int i=0;i<coma.size();i++){
			if (coma[i] != 0.){
				errors.push_back(1/pow(coma[i],.5));
			}else{
				errors.push_back(0.);
			};
		};
		chi2::SetError(i1,j1,errors);
	};
	// i1, j2, interference of wave i1 and j2, correlated with the interference if wave i2 and j2, over all mass bins
	int index1 = i1*_nWaves + j1;
	int index2 = i2*_nWaves + j2; 
	_coma[index1][index2]=coma;
};
void chi2coma::SetComa2(int ii, int jj, std::vector<double> coma){ // This is a little detour, because if takes ii,jj -> i1,j1,i2,j2 and then in SetComa(...) i1,j1,i2,j2 -> ii,jj
	int i1 = ii/_nWaves;					   // But this is not ciritcal, and kept for now for consistency reasons
	int j1 = ii - _nWaves*i1;
	int i2 = jj/_nWaves;
	int j2 = jj - _nWaves*i2;
	SetComa(i1,j1,i2,j2, coma);
};
std::vector<double> chi2coma::GetComa2(int ii, int jj){ // This is a little detour, because if takes ii,jj -> i1,j1,i2,j2 and then in SetComa(...) i1,j1,i2,j2 -> ii,jj
	int i1 = ii/_nWaves;					   // But this is not ciritcal, and kept for now for consistency reasons
	int j1 = ii - _nWaves*i1;
	int i2 = jj/_nWaves;
	int j2 = jj - _nWaves*i2;
	return GetComa(i1,j1,i2,j2);
};



std::vector<double> chi2coma::GetComa(int i1, int j1, int i2, int j2){
	// i1, j2, interference of wave i1 and j2, correlated with the interference if wave i2 and j2, over all mass bins
	int index1 = i1*_nWaves+j1;
	int index2 = i2*_nWaves+j2;
	return _coma[index1][index2];
};

void chi2coma::GetComaFromErrors(){
	std::vector<double> zero;
	for(int bin =0;bin<_nBins;bin++){
		zero.push_back(0.);
	};
	for(int i=0;i<_nWaves*_nWaves;i++){
		for(int j=0;j<_nWaves*_nWaves;j++){
			if (i==j){
				int ii = i/_nWaves;
				int jj = i - ii*_nWaves;
				std::vector<double> err = GetError(ii,jj);
				std::vector<double> comaDiag ;
				for (int bin =0; bin < _nBins; bin++){
					if(err[bin] !=0){
						comaDiag.push_back(1/err[bin]/err[bin]);
					}else{
						comaDiag.push_back(0.);
					};
				};
				_coma[i][j] = comaDiag;
			}else{
				_coma[i][j] = zero;
			};
		};
	};
};


bool chi2coma::CheckConsistency(){
	int nErr=0;	
	if (not chi2::CheckConsistency()){
		nErr+=1;
	};
	for(int i=0;i<_nWaves*_nWaves;i++){
		for(int j=0; j<_nWaves*_nWaves;j++){
			if (_coma[i][j].size() != _nBins){
				int i1 = i/_nWaves;
				int j1 = i - i1*_nWaves;
				int i2 = j/_nWaves;
				int j2 = j - i2*_nWaves;
				std::cout << "Inconsistency found: _coma["<<i<<"]["<<j<<"].size() does not conicide with _nBins. Wave numbers: ["<<i1<<"]["<<j1<<"]#["<<i2<<"]["<<j2<<"]"<<std::endl;
				nErr+=1;
			}else{ // If the number of bins coincides, check symmetry
				for(int bin=0; bin<_nBins;bin++){
					if (pow(_coma[i][j][bin] - _coma[j][i][bin],2)>10.E-15){	//Numerical limit HARDCODED
						std::cout << "Inconsistecy found: Covariance matrix not symmetric: _coma["<<i<<"]["<<j<<"]["<<bin<<"] != _coma["<<j<<"]["<<i<<"]["<<bin<<"]" <<std::endl;
						nErr+=1;
					};
				};
			};
		};
	};
	if (nErr==0){
		return true;
	};
	return false;
};
