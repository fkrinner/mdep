# Simple script to set up a Chi2() for playing around

# Import chi2
import libchi2py
from libchi2py import chi2

#######################################
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#

Chi2=chi2(6)

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/#
#######################################

### Set Up Model ###


Chi2.SetMmin(0.5)
Chi2.SetMmax(2.5)
Chi2.SetBinWidth(0.04)


################################################################[Define_waves]################################################################################
#
#	Some wave definitions
#

Chi2.SetWaveName(0,'1++0+RhoPiS')
Chi2.SetWaveSpin(0,1)
Chi2.SetWaveIsobarMass(0,0.770)
Chi2.SetWavePhaseSpace(0,[1002,10])
Chi2.SetWaveLimits(0,0.9,2.3)
Chi2.SetWaveModel(0,[5,0,101])
Chi2.SetWaveParameterName(0,0,'r_a11260')
Chi2.SetWaveParameterName(0,1,'i_a11260')
Chi2.SetWaveParameterName(0,2,'m_a11260')
Chi2.SetWaveParameterName(0,3,'G_a11260')
Chi2.SetWaveParameterName(0,4,"r_a1")
Chi2.SetWaveParameterName(0,5,"i_a1")
Chi2.SetWaveParameterName(0,6,"m_a1")
Chi2.SetWaveParameterName(0,7,"G_a1")
Chi2.SetWaveParameterName(0,8,'r_NR_1pp')
Chi2.SetWaveParameterName(0,9,'i_NR_1pp')
Chi2.SetWaveParameterName(0,10,'b_NR_1pp')
Chi2.SetWaveParameterName(0,11,'c0_NR_1pp')
Chi2.SetWaveParameterName(0,12,'c1_NR_1pp')
Chi2.SetWaveParameterName(0,13,'c2_NR_1pp')
Chi2.SetWaveParameterName(0,14,'m0_NR_1pp')
Chi2.SetWaveParameterName(0,15,'mPi_NR_1pp')
Chi2.SetWaveParameterName(0,16,'mR_NR_1pp')
Chi2.SetWaveParameterName(0,17,"t_NR_1pp")
Chi2.SetWaveUpperLims(0, [ 10., 0., 1.5   , 1.    , 10., 0., 2.5  , 1.   , 10., 0., 10.   , 10.     , 10.    , 10.    , 0.5  , 0.139, 0.770, 0.10339])

Chi2.SetWaveParameters(0,[ 	1.5117, 0., 1.2683, 0.3726,
				0.72702E-01 , 0.80465E-01,	1.916, 0.236, 
				1.4313, 0.22655,	0.2185,-1.98419 , 1.57851,-2.18077, 0.5  , 0.139, 0.770, 0.10339])

Chi2.SetWaveLowerLims(0, [ 0., 0., 1.    ,0.1    , 0., 0., 1.5  , 0.05 , 0., 0.,-10.   ,-10.     ,-10.    ,-10.    , 0.5  , 0.139, 0.770, 0.10339])


Chi2.SetWaveName(1,'0-+0+f0PiS')
Chi2.SetWaveSpin(1,0)
Chi2.SetWaveIsobarMass(1,0.980)
Chi2.SetWavePhaseSpace(1,[1000,10])
Chi2.SetWaveLimits(1,1.3,2.3)
Chi2.SetWaveModel(1,[0,3])
Chi2.SetWaveParameterName(1,0,'r_pi1800')
Chi2.SetWaveParameterName(1,1,'i_pi1800')
Chi2.SetWaveParameterName(1,2,'m_pi1800')
Chi2.SetWaveParameterName(1,3,'G_pi1800')
Chi2.SetWaveParameterName(1,4,'r_NR_1mp')
Chi2.SetWaveParameterName(1,5,'i_NR_1mp')
Chi2.SetWaveParameterName(1,6,'a_NR_1mp')
Chi2.SetWaveParameterName(1,7,'mPi_NR_1mp')
Chi2.SetWaveParameterName(1,8,'mf0_NR_1mp')
Chi2.SetWaveUpperLims(1, [ 10     , 10     , 2.2  , 0.5  , 10     , 10     , 10 , 0.139, 0.980])

Chi2.SetWaveParameters(1,[ 	-0.44733     , 0.18768     ,1.810, 0.200, 
				-0.42805     , 0.11115     ,-4.5, 0.139, 1.000])

Chi2.SetWaveLowerLims(1, [ 0.     , 0.     , 1.3  , 0.02 , 0.,-1000000,-10 , 0.139, 1.000])


Chi2.SetWaveName(2,'1++0+f0PiP')
Chi2.SetWaveSpin(2,1)
Chi2.SetWaveIsobarMass(2,0.980)
Chi2.SetWavePhaseSpace(2,[1003,10])
Chi2.SetWaveLimits(2,1.3,1.6)
Chi2.SetWaveModel(2,[0,3])
Chi2.SetWaveParameterName(2,0,'r_a11420')
Chi2.SetWaveParameterName(2,1,'i_a11420')
Chi2.SetWaveParameterName(2,2,'m_a11420')
Chi2.SetWaveParameterName(2,3,'G_a11420')
Chi2.SetWaveParameterName(2,4,'r_NR_1pp_f')
Chi2.SetWaveParameterName(2,5,'i_NR_1pp_f')
Chi2.SetWaveParameterName(2,6,'aNR1pp_f0')
Chi2.SetWaveParameterName(2,7,'mPiNR1ppf0')
Chi2.SetWaveParameterName(2,8,'mf0NR1ppf0')
Chi2.SetWaveUpperLims(2, [ 10     , 10     , 2.   , 0.5  , 10     , 10     , 10  , 0.139, 0.980])

Chi2.SetWaveParameters(2,[ 	0.17946     , -0.41955E-01     ,1.418, 0.136, 
				-0.23176     , -0.12469     ,-4.94, 0.139, 0.980])

Chi2.SetWaveLowerLims(2, [ 0.     , 0.     , 1.   , 0.02 , 0.     , 0.     ,-10  , 0.139, 0.980])


Chi2.SetWaveName(3,'2++1+RhoPiD')
Chi2.SetWaveSpin(3,2)
Chi2.SetWaveIsobarMass(3,0.770)
Chi2.SetWavePhaseSpace(3,[1005,10])
Chi2.SetWaveLimits(3,0.9,2.0)
Chi2.SetWaveModel(3,[2,0,101])
Chi2.SetWaveParameterName(3,0,'r_a21320')
Chi2.SetWaveParameterName(3,1,'i_a21320')
Chi2.SetWaveParameterName(3,2,'m_a21320')
Chi2.SetWaveParameterName(3,3,'G_a21320')
Chi2.SetWaveParameterName(3,4,'m_Pi_1')
Chi2.SetWaveParameterName(3,5,'m_Rho_1')
Chi2.SetWaveParameterName(3,6,'m_2_1')
Chi2.SetWaveParameterName(3,7,'BRpiK_1')
Chi2.SetWaveParameterName(3,8,"r_a2")
Chi2.SetWaveParameterName(3,9,"i_a2")
Chi2.SetWaveParameterName(3,10,"m_a2")
Chi2.SetWaveParameterName(3,11,"G_a2")
Chi2.SetWaveParameterName(3,12,'r_NR_2pp')
Chi2.SetWaveParameterName(3,13,'i_NR_2pp')
Chi2.SetWaveParameterName(3,14,'b_NR_2pp')
Chi2.SetWaveParameterName(3,15,'c0_NR_2pp')
Chi2.SetWaveParameterName(3,16,'c1_NR_2pp')
Chi2.SetWaveParameterName(3,17,'c2_NR_2pp')
Chi2.SetWaveParameterName(3,18,'m0_NR_2pp')
Chi2.SetWaveParameterName(3,19,'mPi_NR_2pp')
Chi2.SetWaveParameterName(3,20,'mRhoNR2pp')
Chi2.SetWaveParameterName(3,21,"t_NR_2pp")
Chi2.SetWaveUpperLims(3, [ 10.     , 10., 1.6    , 0.5  , 0.139, 0.770, 0.547, 0.2, 10. , 0.  , 2.2  , 1.   , 10  , 0.  , 10   , 10   , 0., 0., 0.5, 0.139, 0.770, 0.10339])

Chi2.SetWaveParameters(3,[ 	0.72692     , -0.38851E-01 	, 1.313  , 0.110, 0.139, 0.770, 0.547, 0.2, 
				0.12200  	, -0.64628E-01  	, 1.881, 0.442, 
				0.73039E-01  	, -0.12532   	, 0.937,-0.476, 0., 0., 0.5, 0.139, 0.770, 0.10339])

Chi2.SetWaveLowerLims(3, [ 0.      , 0. , 1.0    , 0.05 , 0.139, 0.770, 0.547, 0.2, 0.0 , 0.  , 1.5  , 0.1  , 0.0 , 0.  ,-10   ,-10   , 0., 0., 0.5, 0.139, 0.770, 0.10339])


Chi2.SetWaveName(4,'2-+0+f2PiS')
Chi2.SetWaveSpin(4,2)
Chi2.SetWaveIsobarMass(4,1.270)
Chi2.SetWavePhaseSpace(4,[1004,10])
Chi2.SetWaveLimits(4,1.50,2.3)
Chi2.SetWaveModel(4,[0,0,101])
Chi2.SetWaveParameterName(4,0,'r_pi21660')
Chi2.SetWaveParameterName(4,1,'i_pi21660')
Chi2.SetWaveParameterName(4,2,'m_pi21660')
Chi2.SetWaveParameterName(4,3,'G_pi21660')
Chi2.SetWaveParameterName(4,4,'r_pi21880')
Chi2.SetWaveParameterName(4,5,'i_pi21880')
Chi2.SetWaveParameterName(4,6,'m_pi21880')
Chi2.SetWaveParameterName(4,7,'G_pi21880')
Chi2.SetWaveParameterName(4,8,'r_NR_2mp')
Chi2.SetWaveParameterName(4,9,'i_NR_2mp')
Chi2.SetWaveParameterName(4,10,'b_NR_2mp')
Chi2.SetWaveParameterName(4,11,'c0_NR_2mp')
Chi2.SetWaveParameterName(4,12,'c1_NR_2mp')
Chi2.SetWaveParameterName(4,13,'c2_NR_2mp')
Chi2.SetWaveParameterName(4,14,'m0_NR_2mp')
Chi2.SetWaveParameterName(4,15,'mPi_NR_2mp')
Chi2.SetWaveParameterName(4,16,'mf2_NR_2mp')
Chi2.SetWaveParameterName(4,17,"t_NR_2mp")
Chi2.SetWaveUpperLims(4, [ 10.    , 0.     , 2.0  , 1    , 10.    , 0.     , 2.5  , 1.   , 10.    , 0.     , 10  , 10  , 20   , 0., 0.5, 0.139, 1.275, 0.10339])

Chi2.SetWaveParameters(4,[ 	0.97652     , -0.35428     , 1.654, 0.273, 
				0.82846E-01     , -0.11728     , 1.968, 0.209, 
				0.29964E-01     , -0.37549E-01    , 2.91,-4.12, 11.50, 0., 0.5, 0.139, 1.275, 0.10339])

Chi2.SetWaveLowerLims(4, [ 0.     , 0.     , 1.2  , 0.05 , 0.     , 0.     , 1.5  , 0.04 , 0.     , 0.     ,-10  ,-10  ,-20   , 0., 0.5, 0.139, 1.275, 0.10339])


Chi2.SetWaveName(5,'4++1+RhoPiG')
Chi2.SetWaveSpin(5,4)
Chi2.SetWaveIsobarMass(5,0.770)
Chi2.SetWavePhaseSpace(5,[1001,10])
Chi2.SetWaveLimits(5,1.30,2.3)
Chi2.SetWaveModel(5,[0,3])
Chi2.SetWaveParameterName(5,0,'r_a42040')
Chi2.SetWaveParameterName(5,1,'i_a42040')
Chi2.SetWaveParameterName(5,2,'m_a42040')
Chi2.SetWaveParameterName(5,3,'G_a42040')
Chi2.SetWaveParameterName(5,4,'r_NR_4pp')
Chi2.SetWaveParameterName(5,5,'i_NR_4pp')
Chi2.SetWaveParameterName(5,6,'a_NR_4pp')
Chi2.SetWaveParameterName(5,7,'mPi_NR_4pp')
Chi2.SetWaveParameterName(5,8,'mRhoNR4pp')
Chi2.SetWaveUpperLims(5, [ 10.    , 0.     , 2.5  , 1.5  , 10.    , 0.     , 10.   , 0.139, 0.770])

Chi2.SetWaveParameters(5,[ 	0.25746      , 0.14865     , 1.952, 0.363, 
				0.17766     , -0.91230E-01     ,-1.0434, 0.139, 0.770])

Chi2.SetWaveLowerLims(5, [ 0.     , 0.     , 1.5  , 0.06 , 0.     , 0.     ,-10    , 0.139, 0.770])

################################################################[Load_data]################################################################################
#
#	Load the data_file
#
from wave_data import wave_data
from wave_data import wave_errors

for i in range(6):
	for j in range(6):
		Chi2.SetData(i,j,wave_data[i][j])
		Chi2.SetError(i,j,wave_errors[i][j])


################################################################[Ready_to_go]##############################################################################

if not Chi2.CheckConsistency():
	print "Something is wrong"
else:
	print "Ready to go"

Chi2.PrintParameters()

par = Chi2.GetParameters()

print "Eval([...]) = "+str(Chi2.Eval(par))






