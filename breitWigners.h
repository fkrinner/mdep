// This file contains a number of definitions for amplitude parameterizations. 
// Header only, because of templates.
#ifndef BREITWIGNERS_CILLY_BO
#define BREITWIGNERS_CILLY_BO

#include <cmath> 
#include <vector>
#include <complex>
#include <iostream>
#include <string>
// double mPi=1.3957018;///\pm0.00035MeV // Particle Data Booklet 2012

#ifdef ADOL_ON // Some function on std::complex<adouble>, needed for automatic differentiation.
#include <adolc/adolc.h>  
std::complex<adouble> log(std::complex<adouble> z){
	adouble re = std::real(z);
	adouble im = std::imag(z);

	adouble newReal = pow(re*re+im+im,.5);
	adouble newImag = atan2(im,re);

	return std::complex<adouble>(newReal,newImag);
};
std::complex<adouble> sqrt(std::complex<adouble> z){
	return pow(z,0.5);
};
adouble abs(std::complex<adouble> z){
	adouble squared = std::real(z*std::conj(z));
	return pow(squared,0.5);
};
#endif//ADOL_ON

const double PION_MASS 	= 0.139;
const double PI  	= 3.141592653589793238463;

//////////////////////////  SOME COMMON DEFINITIONS  //////////////////////////////////////////////////////////////////////
template< typename xdouble> xdouble breakupMomentumReal(xdouble M2, xdouble m12, xdouble m22){ // Real breakup momentum: sqrt(lambda(M,m1,m2))/(2M) or 0.
	xdouble lambda= M2*M2 + m12*m12 + m22*m22 - 2*M2*m12 -2*M2*m22 - 2*m12*m22;
	if ( lambda >= 0. ){
		return sqrt(lambda)/(2*M2);
	}else{
		std::cerr << "breitWigners.h: Error: Found 0 > q^2("<<M2<<","<<m12<<","<<m22<<") = "<<lambda/(4*M2*M2)<<". Sub-threshold decay."<<std::endl;
		return 0.;
	};
};

template<typename xdouble> xdouble barrierFactor(xdouble q, int L){
	double pr = 0.1973;
	xdouble z=q*q/pr/pr;
	xdouble res;
	if (L == 0){
		res=1.;
	}else if (L==1){
		res=pow(2*z/(z+1),0.5);
	}else if (L==2){
		res=pow(13*z*z/((z-3)*(z-3)+9*z),0.5);
	}else if (L==3){
		res=pow(277*z*z*z/(z*(z-15)*(z-15)+9*pow(2*z-5,2)),0.5);
	}else if (L==4){
		res=pow(12746*pow(z,4)/(pow(z*z-45*z+105,2)+25*z*pow(2*z-21,2)),0.5);
	} else {
		std::cerr << "breitWigners.h: Error: Barrier factors not defined for L =" << L <<std::endl;
		res =0.;
	};
	return res;
};

template<typename xdouble>  // Some different definition for Barrier factors... used by Dimas program.
xdouble fdl(xdouble P, xdouble R, int L){
	xdouble X = P*R;
	if (L==0){
		return 1.;
	}else if(L==1){
		return 1. + X*X;
	}else if(L==2){
		return 9. + 3.*X*X + X*X*X*X;
	}else if(L==3){
		return 225. + 45.*X*X + 6.*X*X*X*X * X*X*X*X*X*X;	
	}else if(L==4){
		return 11025. + 1575.*X*X + 135.*X*X*X*X + 10*X*X*X*X*X*X + X*X*X*X*X*X*X*X;
	}else{
		std::cerr<<"breitWigners.h: Error: 'fdl(...)' not defined for 4 < L = "<< L << std::endl;
		return 1.;
	};
};

template<typename xdouble>
xdouble psl(xdouble m, xdouble m1, xdouble m2, xdouble R, int L){
	xdouble ampor = m1+m2;
	if (m > ampor){
		xdouble E = (m*m + m1*m1 - m2 * m2)/(2*m);
		xdouble P = pow((E*E-m1*m1)*(E*E-m1*m1),.25);
		xdouble f = fdl<xdouble>(P,R,L);
		return pow(P,2*L+1)/f;
	}else{
		return 0.;
	};
};


template<typename xdouble>
xdouble bowler_integral_table(xdouble m){
	double points[]={
0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.1313E-09, 0.8940E-07, 0.6175E-06, 0.2020E-05, 0.4786E-05, 0.9458E-05, 0.1664E-04, 0.2701E-04, 0.4131E-04, 0.6037E-04, 0.8510E-04, 0.1165E-03, 0.1558E-03, 0.2040E-03, 0.2628E-03, 0.3335E-03, 0.4179E-03, 0.5178E-03, 0.6355E-03, 0.7732E-03, 0.9337E-03, 0.1120E-02, 0.1335E-02, 0.1583E-02, 0.1867E-02, 0.2194E-02, 0.2568E-02, 0.2995E-02, 0.3483E-02, 0.4039E-02, 0.4673E-02, 0.5396E-02, 0.6220E-02, 0.7160E-02, 0.8233E-02, 0.9458E-02, 0.1086E-01, 0.1246E-01, 0.1430E-01, 0.1641E-01, 0.1884E-01, 0.2163E-01, 0.2484E-01, 0.2853E-01, 0.3277E-01, 0.3759E-01, 0.4306E-01, 0.4917E-01, 0.5591E-01, 0.6322E-01, 0.7100E-01, 0.7913E-01, 0.8752E-01, 0.9604E-01, 0.1046E+00, 0.1132E+00, 0.1218E+00, 0.1302E+00, 0.1386E+00, 0.1469E+00, 0.1551E+00, 0.1631E+00, 0.1711E+00, 0.1790E+00, 0.1867E+00, 0.1944E+00, 0.2020E+00, 0.2095E+00, 0.2169E+00, 0.2243E+00, 0.2315E+00, 0.2387E+00, 0.2458E+00, 0.2529E+00, 0.2599E+00, 0.2668E+00, 0.2737E+00, 0.2805E+00, 0.2873E+00, 0.2940E+00, 0.3007E+00, 0.3073E+00, 0.3138E+00, 0.3204E+00, 0.3269E+00, 0.3333E+00, 0.3397E+00, 0.3461E+00, 0.3525E+00, 0.3587E+00, 0.3650E+00, 0.3713E+00, 0.3775E+00, 0.3837E+00, 0.3898E+00, 0.3959E+00, 0.4020E+00, 0.4081E+00, 0.4141E+00, 0.4201E+00, 0.4261E+00, 0.4320E+00, 0.4380E+00, 0.4439E+00, 0.4498E+00, 0.4556E+00, 0.4615E+00, 0.4673E+00, 0.4731E+00, 0.4790E+00, 0.4847E+00, 0.4905E+00, 0.4962E+00, 0.5019E+00, 0.5076E+00, 0.5134E+00, 0.5189E+00, 0.5246E+00, 0.5303E+00, 0.5359E+00, 0.5415E+00, 0.5471E+00, 0.5526E+00, 0.5582E+00, 0.5638E+00, 0.5693E+00, 0.5749E+00, 0.5804E+00, 0.5858E+00, 0.5914E+00, 0.5968E+00, 0.6023E+00, 0.6077E+00, 0.6132E+00, 0.6186E+00, 0.6241E+00, 0.6294E+00, 0.6348E+00, 0.6403E+00, 0.6456E+00, 0.6510E+00, 0.6563E+00, 0.6617E+00, 0.6671E+00, 0.6724E+00, 0.6777E+00, 0.6830E+00, 0.6882E+00, 0.6936E+00, 0.6990E+00, 0.7041E+00, 0.7095E+00, 0.7149E+00, 0.7199E+00, 0.7252E+00, 0.7305E+00, 0.7356E+00, 0.7410E+00, 0.7462E+00, 0.7514E+00, 0.7567E+00, 0.7619E+00, 0.7668E+00, 0.7723E+00, 0.7774E+00, 0.7826E+00, 0.7878E+00, 0.7930E+00, 0.7982E+00, 0.8033E+00, 0.8084E+00, 0.8135E+00, 0.8188E+00, 0.8239E+00, 0.8291E+00, 0.8340E+00, 0.8393E+00, 0.8444E+00, 0.8493E+00, 0.8547E+00, 0.8597E+00, 0.8649E+00, 0.8700E+00, 0.8750E+00, 0.8800E+00, 0.8851E+00, 0.8903E+00, 0.8953E+00, 0.9005E+00, 0.9054E+00, 0.9105E+00, 0.9156E+00, 0.9205E+00, 0.9256E+00, 0.9308E+00, 0.9358E+00, 0.9408E+00, 0.9458E+00, 0.9507E+00, 0.9560E+00, 0.9609E+00, 0.9659E+00, 0.9711E+00, 0.9760E+00, 0.9808E+00, 0.9860E+00, 0.9909E+00, 0.9960E+00, 0.1001E+01, 0.1006E+01, 0.1011E+01, 0.1016E+01, 0.1021E+01, 0.1026E+01, 0.1031E+01, 0.1036E+01, 0.1041E+01, 0.1046E+01, 0.1051E+01, 0.1056E+01, 0.1061E+01, 0.1066E+01, 0.1071E+01, 0.1076E+01, 0.1081E+01, 0.1085E+01, 0.1090E+01, 0.1096E+01, 0.1100E+01, 0.1105E+01, 0.1110E+01, 0.1115E+01, 0.1120E+01, 0.1125E+01, 0.1130E+01, 0.1135E+01, 0.1140E+01, 0.1145E+01, 0.1150E+01, 0.1154E+01, 0.1160E+01, 0.1164E+01, 0.1169E+01, 0.1174E+01, 0.1179E+01, 0.1184E+01, 0.1189E+01, 0.1194E+01, 0.1199E+01, 0.1204E+01, 0.1208E+01, 0.1214E+01, 0.1218E+01, 0.1223E+01, 0.1228E+01, 0.1233E+01, 0.1238E+01, 0.1243E+01, 0.1248E+01, 0.1253E+01, 0.1257E+01, 0.1262E+01, 0.1267E+01, 0.1272E+01, 0.1277E+01, 0.1282E+01, 0.1287E+01, 0.1292E+01, 0.1296E+01, 0.1301E+01, 0.1306E+01, 0.1311E+01, 0.1316E+01, 0.1321E+01, 0.1326E+01, 0.1330E+01, 0.1336E+01, 0.1340E+01, 0.1345E+01, 0.1350E+01, 0.1355E+01, 0.1359E+01, 0.1365E+01, 0.1369E+01, 0.1374E+01, 0.1379E+01, 0.1384E+01, 0.1389E+01, 0.1394E+01, 0.1398E+01, 0.1404E+01, 0.1408E+01, 0.1412E+01, 0.1418E+01, 0.1422E+01, 0.1427E+01, 0.1432E+01, 0.1437E+01, 0.1442E+01, 0.1447E+01, 0.1451E+01, 0.1457E+01, 0.1461E+01, 0.1466E+01, 0.1472E+01, 0.1475E+01, 0.1480E+01, 0.1486E+01, 0.1490E+01, 0.1495E+01, 0.1500E+01, 0.1504E+01, 0.1510E+01, 0.1514E+01, 0.1518E+01, 0.1524E+01, 0.1529E+01, 0.1534E+01, 0.1538E+01, 0.1542E+01, 0.1549E+01, 0.1552E+01, 0.1557E+01, 0.1562E+01, 0.1567E+01, 0.1573E+01, 0.1577E+01, 0.1581E+01, 0.1586E+01, 0.1592E+01, 0.1595E+01, 0.1601E+01, 0.1605E+01, 0.1610E+01, 0.1616E+01, 0.1619E+01, 0.1625E+01, 0.1630E+01, 0.1634E+01, 0.1639E+01, 0.1644E+01, 0.1648E+01, 0.1654E+01, 0.1658E+01, 0.1663E+01, 0.1668E+01, 0.1672E+01, 0.1678E+01, 0.1682E+01, 0.1687E+01, 0.1692E+01, 0.1696E+01, 0.1701E+01, 0.1708E+01, 0.1710E+01, 0.1716E+01, 0.1721E+01, 0.1724E+01, 0.1726E+01};
	double xmin = 0.;
	double xmax = 4.;
	double step = 0.01;
	if (m > xmin and m < xmax){
		int nStep = 0;
		while (xmin+(nStep+1)*step < m){ // nStep is the last integer, where the mass is smaller than the input mass.
			nStep+=1;
		};
		double upper = points[nStep+1];
		double lower = points[nStep];
		double mUpper = xmin + (nStep+1)*step;
		double mLower = xmin + nStep*step;
		xdouble x = (m - mLower)/(mUpper-mLower);
		return (1-x)*lower + x*upper;
	}else{
		return 0.;
	};
};

//////////////////////////////////////  BREIT WIGNER DEFINITIONS  //////////////////////////////////////////////////////////////


template< typename xdouble> std::complex<xdouble> bw(double m, std::vector<xdouble> param, int model, int L=0){ 	// Declare parametrization also in 'getNpars()'. Otherwise, the chi2 doesn't know 
	if (model==-1) {											// the number of parameters. It will then work with 20. Any function, with more than 20
		return std::complex<xdouble>(1.,0.);								// parameters which was not declared in 'getNpars()' will crash.
	};
	if (model==0) { //Simple Breit-Wigner// Reproduces Dimas program
		xdouble m0 = param[0];
		xdouble G0 = param[1];
//		std::cout << "Params are "<< m0 << ' ' << G0 << endl;
		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G0);
		std::complex<xdouble> value = std::complex<xdouble>(m0*G0)/denominator;
//		std::cout<<"EVAL SIMPLE BW("<< m <<", "<< m0 <<", "<< G0 <<"): "<< value <<endl;
		return value;
	};
	if (model==1){ //Mass Dependent Breit-Wigner (one decay channel)
		xdouble m0 = param[0];
		xdouble G0 = param[1];
		xdouble mPi= param[2];
		xdouble mIso=param[3];

		xdouble q0 = breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso*mIso);
		xdouble q  = breakupMomentumReal<xdouble>(m* m ,mPi*mPi,mIso*mIso);
		xdouble Fl = barrierFactor<xdouble>(q,L);
		xdouble Fl0= barrierFactor<xdouble>(q0,L);

		xdouble G  = G0* m0/m * q*Fl*Fl/q0/Fl0/Fl0; //G0 * m0/m q*Fl^2/(q0*Fl0^2)
		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);
		return std::complex<xdouble>(m0*G0,0.)/denominator;	
	};
	if (model==2){ //Mass dependent Breit-Wigner (two decay channels) // Copy Dimas code, do not use Stephan's slides// Reproduces Dimas program
		xdouble m0 = param[0];
		xdouble G0 = param[1];
		xdouble mPi= param[2];
		xdouble mIso1 = param[3];
		xdouble mIso2 = param[4];
		xdouble X  = param[5];

		xdouble R = 5.;

		xdouble psl1 = psl<xdouble>(m, mPi, mIso1, R, L);
		xdouble psl2 = psl<xdouble>(m, mPi, mIso2, R, L);

		xdouble psl10= psl<xdouble>(m0, mPi, mIso1, R, L);
		xdouble psl20= psl<xdouble>(m0, mPi, mIso2, R, L);

		xdouble G = G0 * m0/m * ((1-X) * psl1/psl10 + X * psl2/psl20);

		std::complex<xdouble> denominator  = std::complex<xdouble>(m0*m0-m*m,-m0*G);
		return std::complex<xdouble>(m0*G0,0)/denominator;
	};
	if (model==3){ // Vandermeulen Phase Space// Reproduces Dimas program

		xdouble alpha = param[0];
		xdouble mPi   = param[1];
		xdouble mIso  = param[2];

		xdouble ampor = mPi + mIso;
		std::complex<xdouble> value;

		if ( m > ampor){
			double S = m*m;
			xdouble E = (S + mPi * mPi - mIso*mIso)/(2*m);
			xdouble PSQ = E*E - mPi*mPi;
			value = std::complex<xdouble>(exp(alpha*PSQ),0.);
		}else{
			value = std::complex<xdouble>(1.,0.);			
		};
		return value;
	};
	if (model==4){ // Valera, Dorofeev Background
		xdouble m0 = param[0];
		xdouble alpha = param[1];
		xdouble beta  = param[2];
		
		return std::complex<xdouble>(pow((m-m0)/0.5,alpha)*exp(-beta*(m-m0-0.5)),0);
	};
	if (model==5){ // Bowler parametrization// Reproduces Dimas program
		xdouble m0 = param[0];
		xdouble G0 = param[1];

		xdouble G = G0*	bowler_integral_table<double>(m)/bowler_integral_table<xdouble>(m0) * m0/m;
		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);

		return std::complex<xdouble>(pow(m0*G0,.5),0)/denominator;
	};
	if (model==6){ // Flatte (as in the 3 charged pion release-note from Sep2013
		xdouble m0 = param[0];
		xdouble g1 = param[1];
		xdouble g2 = param[2];
		xdouble mPi= param[3];
		xdouble mK = param[4];
		
		xdouble qpp= breakupMomentumReal<xdouble>(m*m,mPi*mPi,mPi*mPi);
		xdouble qKK= breakupMomentumReal<xdouble>(m*m,mK*mK,mK*mK);

		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-(g1*qpp*qpp + g2*qKK*qKK));
		return std::complex<xdouble>(1,0)/denominator;
	};
	if (model==9){	//Simple gaus for test		
		xdouble m0 = param[0];
		xdouble sig= param[1];

		return std::complex<xdouble>(exp(-(m-m0)*(m-m0)/2/sig/sig),0);

	};

	if (model==10){ //Polynomial c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4
		xdouble c0 = param[0];
		xdouble c1 = param[1];
		xdouble c2 = param[2];
		xdouble c3 = param[3];
		xdouble c4 = param[4];

		xdouble ret = c4*m*m*m*m + c3*m*m*m + c2*m*m + c1*m + c0;
		return std::complex<xdouble>(ret,0);
	};
	if (model==22){ //Mass dependent Breit-Wigner (two decay channels) // From Stephan Schmeing's slides
		xdouble m0 = param[0];
		xdouble G0 = param[1];
		xdouble mPi= param[2];
		xdouble mIso1 = param[3];
		xdouble mIso2 = param[4];
		xdouble X  = param[5];

		xdouble q1 = breakupMomentumReal<xdouble>(m*m,mPi*mPi,mIso1*mIso1);
		xdouble q10= breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso1*mIso1);
		xdouble q2 = breakupMomentumReal<xdouble>(m*m,mPi*mPi,mIso2*mIso2);
		xdouble q20= breakupMomentumReal<xdouble>(m0*m0,mPi*mPi,mIso2*mIso2);
		xdouble Fl1= barrierFactor<xdouble>(q1,L);
		xdouble Fl10=barrierFactor<xdouble>(q10,L);
		xdouble Fl2= barrierFactor<xdouble>(q2,L);
		xdouble Fl20=barrierFactor<xdouble>(q20,L);	

		xdouble G  = G0 * m0/m* ((1-X)*q1*Fl1*Fl1/q10/Fl10/Fl10 + X* q2*Fl2*Fl2/q20/Fl20/Fl20);
		std::complex<xdouble> denominator = std::complex<xdouble>(m0*m0-m*m,-m0*G);
		return std::complex<xdouble>(m0*G0,0.)/denominator;	
	};
	if (model==101){ //t'-dependent background // Reproduces Dimas program
		xdouble b     = param[0];
		xdouble c0    = param[1];
		xdouble c1    = param[2];
		xdouble c2    = param[3];
		xdouble m0    = param[4];
		xdouble mPi   = param[5];
		xdouble mIso  = param[6];
		xdouble tPrime= param[7];


		xdouble PSQ = 0.;
		xdouble mpor = mPi + mIso;		
		if (m > mpor){
			xdouble E = (m*m +mPi*mPi - mIso*mIso)/(2*m);
			PSQ = E*E - mPi*mPi;
		};
		return std::complex<xdouble>(pow(m-0.5,b)*exp(PSQ*(c0+c1*tPrime+c2*tPrime*tPrime)),0.);
	};
	if((model >= 510 and model <= 512) or (model >=610 and model <=612)){ // p-vector-formalism given by Michael Pennington (converted to C++)
		xdouble E = m;		// Mass
		
		xdouble m_Pi = param[0];	// Pion-mass
		xdouble m_K  = param[1];	// Kaon-mass

		xdouble Pi = 3.141592653592;	// Pi

		xdouble S1   = 0.26261;            	// Parameters which MAY NOW be touched.
		xdouble F11  = 0.38949;      
		xdouble F21  = 0.24150;          
		xdouble S2   = 1.0811;            
		xdouble F12  = 0.33961;           
		xdouble F22  =-0.78538;            
		xdouble C110 = 0.14760;           
		xdouble C111 = 0.62181E-01;       
		xdouble C112 = 0.29465E-01;      
		xdouble C120 = 0.10914;           
		xdouble C121 =-0.17912;           
		xdouble C122 = 0.10758;           
		xdouble C220 =-0.27253;           
		xdouble C221 = 0.79442;           
		xdouble C222 =-0.49529; 
		if (model >= 610 and model <= 612){
			S1   = param[2];            	// Parameters which may NOT (!!!) be touched.
			F11  = param[3];      
			F21  = param[4];          
			S2   = param[5];            
			F12  = param[6];           
			F22  = param[7];            
			C110 = param[8];           
			C111 = param[9];       
			C112 = param[10];      
			C120 = param[11];           
			C121 = param[12];           
			C122 = param[13];           
			C220 = param[14];           
			C221 = param[15];           
			C222 = param[16];   
		};
		xdouble S0   = 0.41000*m_Pi*m_Pi;	// Adler zero at (S-S0)
		xdouble eps  = 1E-5;

		xdouble S=E*E;
		std::complex<xdouble> S_c=std::complex<xdouble>(S,eps);
		xdouble rho1=sqrt(1-(4*m_Pi*m_Pi/S));
		std::complex<xdouble> rho1_c = sqrt(std::complex<xdouble>(1.,0.)-std::complex<xdouble>(4*m_Pi*m_Pi,0.)/S_c);
		std::complex<xdouble> rho2_c = sqrt(std::complex<xdouble>(1.,0.)-std::complex<xdouble>(4*m_K*m_K,0.)/S_c);
		std::complex<xdouble> F1_c   = rho1_c/Pi * log((rho1_c+std::complex<xdouble>(1.,0.))/(rho1_c-std::complex<xdouble>(1.,0.)));
		std::complex<xdouble> F2_c   = rho2_c/Pi * log((rho2_c+std::complex<xdouble>(1.,0.))/(rho2_c-std::complex<xdouble>(1.,0.)));

		xdouble rho2=0.;
		if (S>4*m_K*m_K){
			rho2 = abs(rho2_c);
		};

		xdouble Q2 = S/(4*m_K*m_K) - 1;
		xdouble Q4 = Q2*Q2;

		xdouble APL11 = F11*F11/((S1-S)*(S1-S0))+F12*F12/((S2-S)*(S2-S0));
		xdouble AQL11 = C110+C111*Q2+C112*Q4;
		xdouble APL12 = F11*F21/((S1-S)*(S1-S0))+F12*F22/((S2-S)*(S2-S0));
		xdouble AQL12 = C120+C121*Q2+C122*Q4;
		xdouble APL22 = F21*F21/((S1-S)*(S1-S0))+F22*F22/((S2-S)*(S2-S0));
		xdouble AQL22 = C220+C221*Q2+C222*Q4;

		xdouble ALL11 = (APL11+AQL11)*(S-S0)/(4*m_K*m_K);
		xdouble ALL12 = (APL12+AQL12)*(S-S0)/(4*m_K*m_K);
		xdouble ALL22 = (APL22+AQL22)*(S-S0)/(4*m_K*m_K);

		xdouble DET = ALL11*ALL22 - ALL12*ALL12;
		
		std::complex<xdouble> DEN_c = std::complex<xdouble>(1.,0.) + F1_c*ALL11 + F2_c*ALL22 + F1_c*F2_c*DET;

		std::complex<xdouble> T11_c = (ALL11+F2_c*DET)/DEN_c;
		std::complex<xdouble> T12_c = ALL12/DEN_c;
		std::complex<xdouble> T22_c = (ALL22+F1_c*DET)/DEN_c;

		std::complex<xdouble> T11_red_c = T11_c/std::complex<xdouble>((S-S0),0); // Remove Adler zero
		std::complex<xdouble> T12_red_c = T12_c/std::complex<xdouble>((S-S0),0); // Remove Adler zero
		std::complex<xdouble> T22_red_c = T22_c/std::complex<xdouble>((S-S0),0); // Remove Adler zero
		std::complex<xdouble>bw_c;
		if (model==510 or model==610){			// Pi Pi --> Pi Pi
			bw_c=T11_red_c;
		};
		if (model==511 or model==611){			// K  K  --> Pi Pi  (and Pi Pi --> K  K )
			bw_c=T12_red_c;
		};
		if (model==512 or model==612){			// K  K   --> K  K  
			bw_c=T22_red_c;
		};
		return(bw_c);
	};

	std::cerr << "breitWigners.h: Error: Invalid number '" << model << "' for fit-function, use vanishing amplitude: 0. + 0. i."<<std::endl;
	return std::complex<xdouble>(0.,0.);
};

//////////////////////////////////////////  NUMBER OF PARAMETERS NEEDED BY THE DEFINITIONS ABOVE  ///////////////////////////

int getNpars(int model){ 	// Currently has to match the function 'bw(...,model)'.
	if (model==-1){			// Constant function.
		return 0;
	};
	if (model==0){			// Simple Breit-Wigner
		return 2;
	};
	if (model==1){			// Mass dependent Breit-Wigner (one channel )
		return 4;		
	};
	if (model==2){			// Mass dependent Breit-Wigner (two channels)
		return 6;
	};
	if (model==3){			// e^{-\alpha q^2} non-resonant component //Vandermeulen Phase Space
		return 3;
	};
	if (model==4){			// Valera, Dorofeev Background
		return 3;
	};
	if (model==5){			// Bowler parametrization
		return 2;
	};
	if (model==6){			// Flatte
		return 5;
	};
	if (model==9){			// Simple gaus
		return 2;
	};
	if (model==10){			// Polynomial up to 4th order
		return 5;
	};
	if (model==22){	
		return 6.;
	};
	if (model==101){		// t' dependent background
		return 8;
	};
	if (model>=510 and model <= 512){// Different matrix elements of a p-vector approach (Given by Michael Pennington).
		return 2;
	};
	if (model>=610 and model <= 612){// Different matrix elements of a p-vector approach (Given by Michael Pennington). (Parameters are touched !!!!!)
		return 17;
	};
	std::cout << "breitWigners.h: Warning: Invalid number '"<< model <<"' for fit-function. Use 20 parameters to be sure." << std::endl;
	return 20;
};

///////////////////////////////////// NAMES OF THE PARAMETRIZATIONS ///////////////////////////////////////////////////

std::string getModelName(int model){ 	// Currently has to match the function 'bw(...,model)'.
	if (model==-1){
		return "constant";
	};
	if (model==0){
		return "breit-wigner";
	};
	if (model==1){
		return "mass-dependent-BW";
	};
	if (model==2){
		return "mass-dependent-BW-2channels-dima";
	};
	if (model==3){
		return "vandermeulen-phase-space";
	};
	if (model==4){
		return "valera-dorofeev";
	};
	if (model==5){
		return "bowler-parametrization";	
	};
	if (model==6){
		return "flatte";
	};
	if (model==9){
		return "gaus";
	};
	if (model==10){
		return "polynomial";
	};
	if (model==22){
		return "mass-dependent-BW-2channels-dima";
	};
	if (model==101){
		return "t'-dependent";
	};
	if (model>=510 and model <= 512){
		return "pennington-p-vector";
	};
	if (model>=610 and model <= 612){
		return "pennington-p-vector-released";
	};
	return "unknown";
};

#endif
