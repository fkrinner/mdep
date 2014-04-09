#ifndef CHI_TO_PY_COWABUNGA
#define CHI_TO_PY_COWABUNGA

#include<string>

namespace bp = boost::python;

template<class T>
bp::list std_vector_to_py_list(const std::vector<T>& v){
	bp::object get_iter = bp::iterator<std::vector<T> >();
	bp::object iter = get_iter(v);
	bp::list l(iter);
	return l;
};


template< typename T >
std::vector<T> to_std_vector(const bp::object& iterable){
    return std::vector<T>(bp::stl_input_iterator<T>(iterable), bp::stl_input_iterator<T>( ));
};

template< typename T >
std::vector<std::vector<T> > to_std_vector_vector(const bp::object& iterable){
	std::vector<bp::object> intermediate = to_std_vector<bp::object>(iterable);
	std::vector<T> line;
	std::vector<std::vector<T> > final;
	for (unsigned int i=0;i<intermediate.size();i++){
		line=to_std_vector<T>(intermediate[i]);
		final.push_back(line);
	};
	return final;	
};

// Overwrites methods, that take or return a std::vector<...> with versions that take or return bp::objects/bp::lists
// Depending on preprocessor-flags, special methods are added

template<class T>
struct chi_to_py : public T{

	chi_to_py(const int nWaves) : T(nWaves){};


	double Eval(bp::object& param){
		std::vector<double>  xxx = to_std_vector<double>(param);
		return T::Eval(xxx);
	};
	double EvalRel(bp::object& param){
		std::vector<double>  xxx = to_std_vector<double>(param);
		return T::EvalRel(xxx);
	};

	double EvalNorm(bp::object& param){
		std::vector<double> xxx = to_std_vector<double>(param);
		return T::EvalNorm(xxx);
	};

	std::complex<double> WaveAmp(int i,double min, double max, bp::object& param){
		std::vector<double>  xxx = to_std_vector<double>(param);
		return T::WaveAmp(i,min, max,xxx);
	};

	void SetWaveParameters(int i, bp::object pars){
		std::vector<double> xxx = to_std_vector<double>(pars);
		T::SetWaveParameters(i,xxx);
	};
	void SetRelParameters(bp::object pars){
		std::vector<double> xxx = to_std_vector<double>(pars);
		T::SetRelParameters(xxx);
	};
	void SetRelNormParameters(bp::object pars){
		std::vector<double> xxx = to_std_vector<double>(pars);
		T::SetRelNormParameters(xxx);
	};
	bp::list GetParameters(){
		std::vector<double> xxx = T::GetParameters();
		return std_vector_to_py_list(xxx);
	};
	bp::list GetWavePars(int i){
		std::vector<double> pars = T::GetWaveParameters(i);
		return std_vector_to_py_list(pars);
	};
	bp::list GetRelParameters(){
		std::vector<double> pars = T::GetRelParameters();
		return std_vector_to_py_list(pars);
	};
	bp::list GetRelNormParameters(){
		std::vector<double> pars = T::GetRelNormParameters();
		return std_vector_to_py_list(pars);
	};
	void SetWaveModel(int i, bp::object &model){
		std::vector<int> xxx = to_std_vector<int>(model);
		T::SetWaveModel(i,xxx);
	};

	bp::list GetWaveModel(int i){
		std::vector<int> model = T::GetWaveModel(i);
		return std_vector_to_py_list(model);
	};
	bp::list GetParStat(){
		std::vector<bool> xxx = T::GetParStat();
		return std_vector_to_py_list(xxx);
	};
	bp::list GetParNames(){
		std::vector<std::string> xxx = T::GetParNames();
		return std_vector_to_py_list(xxx);
	};
	bp::list GetRelParNames(){
		std::vector<std::string> names = T::GetRelParNames();
		return std_vector_to_py_list(names);
	};
	bp::list GetWaveRelParNumbers(int i){
		std::vector<int> pars = T::GetWaveRelParNumbers(i);
		return std_vector_to_py_list(pars);
	};
	bp::list GetWaveUpperLims(int i){
		std::vector<double> lim = T::GetWaveUpperLims(i);
		return std_vector_to_py_list(lim);
	};
	bp::list GetWaveLowerLims(int i){
		std::vector<double> lim = T::GetWaveLowerLims(i);
		return std_vector_to_py_list(lim);
	};
	void SetWaveUpperLims(int i, bp::object lim){
		std::vector<double> xxx = to_std_vector<double>(lim);
		T::SetWaveUpperLims(i,xxx);
	};
	void SetWaveLowerLims(int i, bp::object lim){
		std::vector<double> xxx = to_std_vector<double>(lim);
		T::SetWaveLowerLims(i,xxx);
	};
	void SetWavePhaseSpace(int i, bp::object ps){
		std::vector<int> xxx = to_std_vector<int>(ps);
		T::SetWavePhaseSpace(i, xxx);
	};
	void SetPhaseSpace(bp::object ps){
		std::vector<int> xxx = to_std_vector<int>(ps);
		T::SetPhaseSpace(xxx);
	};
	bp::list GetWavePhaseSpace(int i){
		std::vector<int> xxx = T::GetWavePhaseSpace(i);
		return std_vector_to_py_list(xxx);
	};

	bp::list GetBinning(){
		std::vector<double> xxx = T::GetBinning();
		return std_vector_to_py_list(xxx);
	};
	void SetBinning(bp::object binning){
		std::vector<double> xxx = to_std_vector<double>(binning);
		T::SetBinning(xxx);
	};
	bp::list GetPoints(int i, int nPoints){
		std::vector<std::complex<double> > xxx = T::GetPoints(i,nPoints);
		return std_vector_to_py_list(xxx);
	};
	bp::list GetPointsComponent(int i, bp::object comps, int nPoints){
		std::vector<int> xxx1 = to_std_vector<int>(comps);
		std::vector<std::complex<double> > xxx2 = T::GetPointsComponent(i,xxx1,nPoints);
		return std_vector_to_py_list(xxx2);
	};
	void SetSampleAmps(int bin, bp::object amps){
		std::vector<double> xxx = to_std_vector<double>(amps);
		std::vector<std::complex<double> > yyy;
		int nWaves = T::GetNwaves();
		for (int i =0;i<nWaves;i++){
			yyy.push_back(std::complex<double>(xxx[2*i],xxx[2*i+1]));
		};
		T::SetSampleAmps(bin,yyy);
	};
#ifdef IS_RHO_CLASS_QWERT
	void SetData(int i, int j, bp::object data){
		std::vector<double> xxx = to_std_vector<double>(data);
		T::SetData(i,j,xxx);
	};
	void SetError(int i, int j, bp::object error){
		std::vector<double> xxx = to_std_vector<double>(error);
		T::SetError(i,j,xxx);
	};
	bp::list GetData(int i, int j){
		std::vector<double> data = T::GetData(i,j);
		return std_vector_to_py_list(data);
	};
	bp::list GetError(int i, int j){
		std::vector<double> error = T::GetError(i,j);
		return std_vector_to_py_list(error);
	};
#endif//IS_RHO_CLASS_QWERT

#ifdef ADOL_ON
	bp::list Gradient(bp::object par){
		std::vector<double> xxx = to_std_vector<double>(par);
		std::vector<double> grad = T::Gradient(xxx);
		return std_vector_to_py_list(grad);
	};
#endif//ADOL_ON

#ifdef IS_COMA_CLASS_FRRR 
	double EvalComa(bp::object& param){ 
		std::vector<double>  xxx = to_std_vector<double>(param);
		return T::EvalComa(xxx);
	};
	bp::list GetComa(int i1, int j1, int i2, int j2){
		std::vector<double> xxx = T::GetComa(i1,j1,i2,j2);
		return std_vector_to_py_list(xxx);
	};
	void SetComa(int i1, int j1, int i2, int j2, bp::object in){
		std::vector<double> xxx = to_std_vector<double>(in);
		T::SetComa(i1,j1,i2,j2, xxx);
	};
	bp::list GetComa2(int ii, int jj){
		std::vector<double> xxx = T::GetComa2(ii,jj);
		return std_vector_to_py_list(xxx);
	};
	void SetComa2(int ii, int jj, bp::object in){
		std::vector<double> xxx = to_std_vector<double>(in);
		T::SetComa2(ii,jj, xxx);
	};
#endif//IS_COMA_CLASS_FRRR

#ifdef IS_AMP_CLASS_KUL
	void SetData(int i, bp::object data){
		std::vector<double> xxx = to_std_vector<double>(data);
		T::SetData(i,xxx);
	};
	bp::list GetData(int i){
		std::vector<double> data = T::GetData(i);
		return std_vector_to_py_list(data);
	};
	bp::list GetComa(int i, int j){
		std::vector<double> xxx = T::GetComa(i,j);
		return std_vector_to_py_list(xxx);
	};
	void SetComa(int i, int j, bp::object in){
		std::vector<double> xxx = to_std_vector<double>(in);
		T::SetComa(i,j, xxx);
	};
#endif//IS_AMP_CLASS_KUL
};


template<class T>
void export_chi2(std::string name) { 
	// Expose the std::vectors<...> to python
	
	bp::class_<std::vector<double> >("vector_double")
		.def(bp::vector_indexing_suite<std::vector<double> >())
	;
	bp::class_<std::vector<bool> >("vector_bool")
		.def(bp::vector_indexing_suite<std::vector<bool> >())
	;
	bp::class_<std::vector<std::string> >("vector_string")
		.def(bp::vector_indexing_suite<std::vector<std::string> >())
	;
	bp::class_<std::vector<int> >("vector_int")
		.def(bp::vector_indexing_suite<std::vector<int> >())
	;

	bp::class_<std::vector<std::complex<double> > >("vector_complex")
		.def(bp::vector_indexing_suite<std::vector<std::complex<double> > >())
	;

	bp::class_<std::vector<std::vector<double> > >("vector_array_double")
		.def(bp::vector_indexing_suite<std::vector<std::vector<double> > >())
	;


	bp::class_<std::vector<std::vector<std::vector<double> > > >("vector_array_array_double") // To be able to have a vector of covariance matrices.
		.def(bp::vector_indexing_suite<std::vector<std::vector<std::vector<double> > > >())
	;

    bp::class_<T>(name.c_str(),bp::init<const int>())
		//Expose methods
		.def("EvalRel",		&T::EvalRel								)
		.def("Eval",		&T::Eval								)
		.def("EvalNorm",	&T::EvalNorm								)
		.def("WaveAmp",		&T::WaveAmp								)
		.def("PrintWaveModel", 	&T::PrintWaveModel							)
		.def("PrintParameters", &T::PrintParameters							)
		.def("PrintWaveParameters",&T::PrintWaveParameters						)
//		//Expose setters and getters
		.def("ClassName",	&T::ClassName								)
		.def("GetWaveName",	&T::GetWaveName								)
		.def("SetWaveName",	&T::SetWaveName								)
		.def("GetParameters",	&T::GetParameters							)
		.def("GetWaveParameter",&T::GetWaveParameter							)
		.def("GetWaveParameters",&T::GetWaveParameters							)
		.def("SetWaveParameter",&T::SetWaveParameter							)
		.def("SetWaveParameters",&T::SetWaveParameters							)
		.def("GetNwaves", 	&T::GetNwaves								)
		.def("GetNpars", 	&T::GetNpars								)
		.def("GetWaveNpars", 	&T::GetWaveNpars							)
		.def("SetWaveLimits", 	&T::SetWaveLimits							)
		.def("GetWaveMaxLimit",	&T::GetWaveMaxLimit							)
		.def("GetWaveMinLimit",	&T::GetWaveMinLimit							)
		.def("GetParNames",	&T::GetParNames								)
		.def("GetWaveParameterName",&T::GetWaveParameterName						)
		.def("SetWaveParameterName",&T::SetWaveParameterName						)
		.def("SetWaveModel", 	&T::SetWaveModel							)		
		.def("GetWaveModel", 	&T::GetWaveModel							)
		.def("SetData",		&T::SetData								)
		.def("GetData",		&T::GetData								)
		.def("SetMmax",		&T::SetMmax								)
		.def("SetMmin",		&T::SetMmin								)
		.def("SetBinWidth",	&T::SetBinWidth								)
		.def("GetMmax",		&T::GetMmax								)
		.def("GetMmin",		&T::GetMmin								)
		.def("GetBinWidth",	&T::GetBinWidth								)
		.def("GetNbins",	&T::GetNbins								)
		.def("CheckConsistency",&T::CheckConsistency							)
		.def("GetParStat",	&T::GetParStat								)
		.def("GetWaveParStat",	&T::GetWaveParStat							)
		.def("FixWaveParameter",&T::FixWaveParameter							)
		.def("RelWaveParameter",&T::RelWaveParameter							)
		.def("GetNparsReleased",&T::GetNparsReleased							)
		.def("GetWaveNparsReleased",&T::GetWaveNparsReleased						)
		.def("GetRelParNames",	&T::GetRelParNames							)
		.def("GetWaveRelParNumbers",&T::GetWaveRelParNumbers						)
		.def("GetWaveUpperLims",&T::GetWaveUpperLims							)
		.def("GetWaveLowerLims",&T::GetWaveLowerLims							)
		.def("SetWaveUpperLims",&T::SetWaveUpperLims							)
		.def("SetWaveLowerLims",&T::SetWaveLowerLims							)
		.def("GetWaveUpperLim",&T::GetWaveUpperLim							)
		.def("GetWaveLowerLim",&T::GetWaveLowerLim							)
		.def("SetWaveUpperLim",&T::SetWaveUpperLim							)
		.def("SetWaveLowerLim",&T::SetWaveLowerLim							)
		.def("GetWaveNormParameter",&T::GetWaveNormParameter						)
		.def("SetWaveNormParameter",&T::SetWaveNormParameter						)
		.def("GetParFromNorm"	,&T::GetParFromNorm							)
		.def("GetNormFromPar"	,&T::GetNormFromPar							)
		.def("GetRelParameters"	,&T::GetRelParameters							)
		.def("GetRelNormParameters",&T::GetRelNormParameters						)
		.def("SetRelParameters"	,&T::SetRelParameters							)
		.def("SetRelNormParameters",&T::SetRelNormParameters						)
		.def("SetWaveSpin"	,&T::SetWaveSpin							)
		.def("GetWaveSpin"	,&T::GetWaveSpin							)
		.def("GetWaveIsobarMass",&T::GetWaveIsobarMass							)
		.def("SetWaveIsobarMass",&T::SetWaveIsobarMass							)
		.def("SetPhaseSpace"	,&T::SetPhaseSpace							)
		.def("SetWavePhaseSpace",&T::SetWavePhaseSpace							)
		.def("GetWavePhaseSpace",&T::GetWavePhaseSpace							)
		.def("SetBinning"	,&T::SetBinning								)
		.def("GetBinning"	,&T::GetBinning								)
		.def("SetRandomCouplings",&T::SetRandomCouplings						)
		.def("SetCouplingSize"	,&T::SetCouplingSize							)
		.def("GetNcalls"	,&T::GetNcalls								)
		.def("SetRphi"		,&T::SetRphi								)
		.def("SetWaveRphi"	,&T::SetWaveRphi							)
		.def("GetWaveRphi"	,&T::GetWaveRphi							)
		.def("GetPoints"	,&T::GetPoints								)
		.def("GetPointsComponent",&T::GetPointsComponent						)
		.def("IsSampleAmps"	,&T::IsSampleAmps							)
		.def("SetSampleAmps"	,&T::SetSampleAmps							)
#ifdef IS_AMP_CLASS_KUL
		.def("GetComa",		&T::GetComa								)
		.def("SetComa",		&T::SetComa								)
#endif//IS_AMP_CLASS_KUL
#ifdef IS_RHO_CLASS_QWERT
		.def("SetError",	&T::SetError								)
		.def("GetError",	&T::GetError								)
		.def("Conjugate",	&T::Conjugate								)
		.def("SetAnchorWave"	,&T::SetAnchorWave							)
#endif//IS_RHO_CLASS_QWERT
#ifdef ADOL_ON
		.def("Gradient",	&T::Gradient								)
#endif//ADOL_ON

#ifdef IS_COMA_CLASS_FRRR
		.def("EvalComa",	&T::EvalComa								)
		.def("GetComaFromErrors",&T::GetComaFromErrors							)
		.def("GetComa",		&T::GetComa								)
		.def("SetComa",		&T::SetComa								)
		.def("GetComa2",	&T::GetComa2								)
		.def("SetComa2",	&T::SetComa2								)
#endif//IS_COMA_CLASS_FRRR
 	;
};


#endif//CHI_TO_PY_COWABUNGA
