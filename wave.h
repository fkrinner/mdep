#ifndef WAVE_ZERGRUSH_OLOLOLOLO
#define WAVE_ZERGRUSH_OLOLOLOLO
#include <vector>
#include <complex>
class wave {
	public:
//		template<typename xdouble>
		wave();
		wave(std::vector<int> model);

		template<typename xdouble>
		std::complex<xdouble> Ampl(double m,std::vector<xdouble> &x) const; 
		template<typename xdouble>
		std::complex<xdouble> AmplRel(double m,std::vector<xdouble> &x) const;
		double RealPart(std::vector<double> m,std::vector<double> param);
		double ImagPart(std::vector<double> m,std::vector<double> param);

		void PrintModelInfo();
		void PrintParameters();

		double GetParameter(int i);
		void SetParameter(int i, double par);
		std::vector<double> GetParameters();
		void SetParameters(std::vector<double> pars);

		std::string GetName();
		void SetName(std::string name);
		int GetNpars();
		void SetLimits(double min, double max);
		double GetMaxLimit();
		double GetMinLimit();
		std::string GetParName(int i);
		std::string ClassName();
		void SetParName(int i,std::string name);
		void SetModel(std::vector<int> model);
		std::vector<int> GetModel();

		bool GetParStat(int i);
		void FixParameter(int i);
		void RelParameter(int i);
		int GetNparsReleased();
		std::vector<std::string> GetRelParNames();
		std::vector<double> GetRelParValues();
		std::vector<int> GetRelParNumbers();

		std::vector<double> GetUpperLims();
		std::vector<double> GetLowerLims();
		double GetUpperLim(int n);
		double GetLowerLim(int n);
		void SetUpperLims(std::vector<double> lim);
		void SetLowerLims(std::vector<double> lim);
		void SetUpperLim(int n, double lim);
		void SetLowerLim(int n, double lim);

		void SetPhaseSpace(std::vector<int> ps);
		std::vector<int> GetPhaseSpace();

		void SetSpin(int L);
		int GetSpin();
		void SetIsobarMass(double m);
		double GetIsobarMass();

		void SetRandomCouplings(bool real);
		std::vector<bool> Couplings();

		void SetRphi(bool flag);
		bool GetRphi();

	protected:
		std::string _name;
		std::vector<int> _model;			// Vector containing the numbers of different parametrizations
		std::vector<int> _nPars;			// Vector containing the numbers of parameters of the different parametrizations
		std::vector<int> _phaseSpace;			// Vector containing the numbers of phase space prefactors
		int _nPar;					// Number of parameters
		std::vector<double> _parameters;		// Vector containing the actual parameters used for print
		std::vector<double> _upperLim;
		std::vector<double> _lowerLim;
		std::vector<bool> _release_parameter;		// Tells, if a parameter is relreased or not
		std::vector<std::string> _parameterNames;	// Parameter names
		double _minLim;					// Fit range 
		double _maxLim;					// Fit range
		int _L;						// Spin of the wave
		double _isobar_mass;				// Mass os the isobar
		bool _rPhi;					// Switch from re/im to r/phi in the couplings
};

#endif //WAVE_ZERGRUSH_OLOLOLOLO
