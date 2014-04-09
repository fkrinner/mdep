#ifndef WAVESET_COWABUNGA
#define WAVESET_COWABUNGA

#include "wave.h"
#include <string>
#include <complex>

#ifdef ADOL_ON
#include <adolc/adolc.h>  
#endif//ADOL_ON

class waveset {
	public:

		waveset(const int nWaves);

		template<typename xdouble>
		std::vector<std::complex<xdouble> > Amps(int bin, std::vector<xdouble> &par);

		template<typename xdouble>
		std::complex<xdouble> WaveAmp(int i, double min, double max , std::vector<xdouble> &par) const;

#ifdef ADOL_ON
		virtual adouble EvalAdouble(std::vector<adouble> par);
		std::vector<double> Gradient(std::vector<double> per);
#endif//ADOL_ON

		void SetRelParameters(std::vector<double> par);
		void SetRelNormParameters(std::vector<double> par);
		std::vector<double> GetRelParameters();
		std::vector<double> GetRelNormParameters();

		void SetWaveModel(int i, std::vector<int> model);
		std::vector<int> GetWaveModel(int i);
		
		std::vector<double> GetParameters(); ///////////////////////////
		void SetWaveParameter(int i, int n, double par);
		double GetWaveParameter(int i, int n);

		void SetWaveParameters(int i, std::vector<double> pars);
		std::vector<double> GetWaveParameters(int i);

		void SetWaveNormParameter(int i, int n, double par);
		double GetWaveNormParameter(int i , int n);

		double GetParFromNorm(int i, int n, double norm);
		double GetNormFromPar(int i, int n, double par);

		void SetWaveParameterName(int i, int n, std::string name);
		std::string GetWaveParameterName(int i, int n);
		std::vector<std::string> GetParNames();/////////////////////////

		std::string ClassName();

		void SetWaveName(int i, std::string name);
		std::string GetWaveName(int i);

		void SetWaveLimits(int i, double lower, double upper);	
		double GetWaveMaxLimit(int i);
		double GetWaveMinLimit(int i);

		int GetNwaves();
		int GetNpars();
		int GetWaveNpars(int i);

		void PrintParameters();
		void PrintWaveParameters(int i);
		void PrintWaveModel(int i);

		void SetMmax(double m);
		void SetMmin(double m);
		void SetBinWidth(double m);
		double GetMmax();
		double GetMmin();
		double GetBinWidth();
		int GetNbins();
		
		bool GetWaveParStat(int i, int n);
		std::vector<bool> GetParStat();
		void FixWaveParameter(int i, int n);
		void RelWaveParameter(int i, int n);
		int GetNparsReleased();
		int GetWaveNparsReleased(int i);
		std::vector<std::string> GetRelParNames();
		std::vector<int> GetWaveRelParNumbers(int i);

		std::vector<double> GetWaveUpperLims(int i);
		std::vector<double> GetWaveLowerLims(int i);
		double GetWaveUpperLim(int i,int n);
		double GetWaveLowerLim(int i, int n);
		void SetWaveUpperLims(int i, std::vector<double> lim);
		void SetWaveLowerLims(int i, std::vector<double> lim);
		void SetWaveUpperLim(int i, int n, double lim);
		void SetWaveLowerLim(int i, int n, double lim);

		void SetWaveSpin(int i, int L);
		int GetWaveSpin(int i);
		void SetWaveIsobarMass(int i, double m);
		double GetWaveIsobarMass(int i);

		void SetWavePhaseSpace(int i, std::vector<int> ps);
		std::vector<int> GetWavePhaseSpace(int i);
		void SetPhaseSpace(std::vector<int> ps);

		void SetBinning(std::vector<double>);
		std::vector<double> GetBinning();

		void SetRphi(bool flag);
		void SetWaveRphi(int i,bool flag);
		bool GetWaveRphi(int i);

		int GetNcalls();

		bool CheckConsistency();

		void IsSampleAmps(bool flag);
		void SetSampleAmps(int bin, std::vector<std::complex<double> > in);

	protected:
		const int _nWaves;						// Number of waves
		int _nPars;							// Total number of parameters
		std::vector<int> _waveNpars;					// Number of parameters for single waves
		std::vector<wave> _waves;					// The waves
		int _nBins;
		bool _isSampleAmps;						// If true, 'Amps(...)' returns '_sampleAmps[...]' instead of acutal amplitudes
		std::vector<std::vector<std::complex<double> > > _sampleAmps;	//
		double _mMax;							// Mass maximum
		double _mMin;							// Mass minimum
		double _binWidth;						// Width of the binning
		std::vector<double> _binning;					// Binning
		int _nCalls;
};
#endif //WAVESET_COWABUNGA
