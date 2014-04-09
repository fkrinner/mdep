import libchi2amppy
import sys
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
from libchi2amppy import chi2amp
from types import MethodType
from convertTextOutput import getComaAmp

"""Some class method definitions for the chi2amp class, that are useful, but not in the C++ code"""



def releaseParameterByName(self,name):
	"""Releases all parmeters of a chi2 with 'name'"""
	nWave = self.GetNwaves()
	for i in range(nWave):
		nPar = self.GetWaveNpars(i)
		for j in range(nPar):
			if self.GetWaveParameterName(i,j) == name:
				self.RelWaveParameter(i,j)

def fixParameterByName(self,name):
	"""Fixes all parmeters of a chi2 with 'name'"""
	nWave = self.GetNwaves()
	for i in range(nWave):
		nPar = self.GetWaveNpars(i)
		for j in range(nPar):
			if self.GetWaveParameterName(i,j) == name:
				self.FixWaveParameter(i,j)

def setParameterByName(self,name,val):
	"""Sets the parameter with 'name' to val"""
	nWave = self.GetNwaves()
	for i in range(nWave):
		nPar = self.GetWaveNpars(i)
		for j in range(nPar):
			if self.GetWaveParameterName(i,j) == name:
				self.SetWaveParameter(i,j,val)


def getPlot(self,i,nBins=0):
	"""
	Returns the plots for data and fit for wave i interfering with j
	The format is: [[mass,binwidths,data,errors],[mass,fit]]
	"""
	if nBins ==0:
		nBins = self.GetNbins()
	binning = self.GetBinning()
	width = (binning[-1]-binning[0])/nBins
	anc = self.GetPoints(0,nBins)
	amp = self.GetPoints(int((i+1)/2),nBins)
	for j in range(len(amp)):
		try:
			amp[j]*=abs(anc[j])/anc[j]
		except:
			amp[j] = 0.+0.j
	fit=[]
	mFit=[]
	mDat=[]
	dat=self.GetData(i)
	err=self.GetComa(i,i)
	for j in range(len(err)):
		err[j]**=-.5
	errm = []
	for bin in range(nBins):
		mFit.append(binning[0] + (bin+0.5)*width)
		mDat.append((binning[bin] + binning[bin+1])/2)
		errm.append((binning[bin+1]-binning[bin])/2)
		if i==0:
			fit.append(abs(amp[bin]))
		elif i%2 == 1:
			fit.append(amp[bin].real)
#			print "darilpat"
			
		else:
			fit.append(amp[bin].imag)
#			print "imagpat"
	return [[mDat,errm,dat,err],[mFit,fit]]


def pyMinuitFunction(self, name='Chi2'):
	"""
	Generates a function, that can be used by PyMinuit (every parameter defined explicitely)
	The Chi2 instance has to be called 'name' to work
	The function to be called by PyMinuit is then 'Chi2PM'
	"""
	parNames = self.GetParNames()
	parString= ', '.join(parNames)
	command = "def Chi2PM("+parString+"): return "+name+".Eval(["+parString+"])"
	return command


def loadFit(self,waves,direct):
	"""Loads the waves in 'waves' from the fit in 'direct'. The fit has to have text output"""
	if not len(waves) == self.GetNwaves():
		print "Wrong number of waves"
		return
	nWaves = self.GetNwaves()
	data = getComaAmp(waves,direct)
	pointsSorted=[]
	for i in range(nWaves):
		self.SetWaveName(i,waves[i].strip())
	for i in range(2*nWaves-1):
		pointsSorted.append([])
	comaSorted=[]
	for i in range(2*nWaves-1):
		comaSortedLine=[]
		for j in range(2*nWaves-1):
			comaSortedLine.append([])
		comaSorted.append(comaSortedLine)
	for bin in range(len(data[0])):
		for i in range(2*nWaves-1):
			if i ==0:
				ii=0
			else:
				ii=i+1
			if ii%2 ==0:
				pointsSorted[i].append(data[0][bin][ii])
			else:
				pointsSorted[i].append(-data[0][bin][ii])
			for j in range(2*nWaves-1):
				if j == 0:
					jj = 0
				else:
					jj = j+1				
				comaSorted[i][j].append((data[1][bin][ii][jj]).real) # Should be real anyway
	for i in range(2*nWaves-1):
		self.SetData(i,pointsSorted[i])
	for ii in range(2*nWaves-1):
		for jj in range(2*nWaves-1):
			self.SetComa(ii,jj,comaSorted[ii][jj])

data = None
def sampleFits(self, waves, direct, level = 1., verbose = False):
	from sample_coma import sample_coma
	import numpy as np
	from numpy import linalg as la
	global data
	if data is None:	
		data = getComaAmp(waves,direct) # Calling again with another set does not work.
	self.IsSampleAmps(True)
	toot_deel=0.
	for bin in range(len(data[0])):
		coma = la.inv(np.matrix(data[1][bin])).tolist()
		setter = data[0][bin][:]
		for i in range(len(setter)):
			if i%2==1:
				setter[i]*=-1
		sample = sample_coma(setter,coma, level).tolist()
		if verbose:
			cco = np.matrix(data[1][bin])
			ddl = np.matrix(sample) - np.matrix(setter)
			deel = (ddl*cco*np.transpose(ddl)).tolist()[0][0]
			print "DELTA["+str(bin)+"]: "+str(deel)
			toot_deel+=deel
		self.SetSampleAmps(bin,sample)
	ev = self.Eval([])
	self.IsSampleAmps(False)
	if verbose:
		print "TOTAL_DELTA: "+str(toot_deel)
	return ev

def sampleParameter(self,npar,waves,direct):
	import random
	from random import gauss
	chi2s=[]
	check=self.SampleAmp(waves,direct,0.,False) # Sets all SampleAmps to the mean values ('level' == 0.)
	if check > 1.E-10:
		print "Chi2 does not give zero at zero deviation. Seems like wrong parameters."
	global data
	if data is None:	
		data = getComaAmp(waves,direct) # Calling again with another set does not work.
	self.IsSampleAmps(True)
	for bin in range(len(data[0])):
		setter = data[0][bin][:]
		for i in range(len(setter)):
			if i%2==1:
				setter[i]*=-1
		sig = data[1][bin][npar][npar]**.5
		val = gauss(setter[npar],sig)
		sample=setter[:]
		sample[npar] = val
		self.SetSampleAmps(bin,sample)
		ev = self.Eval([])
#		print ev
		chi2s.append(ev)
		self.SetSampleAmps(bin,setter)	
	return chi2s

# Add methods tho chi2amp
chi2amp.RelParameterName 	= MethodType(	releaseParameterByName,		None, chi2amp)
chi2amp.FixParameterName 	= MethodType(	fixParameterByName, 		None, chi2amp)
chi2amp.SetParameterByName	= MethodType(	setParameterByName, 		None, chi2amp)
chi2amp.GetPlot 		= MethodType(	getPlot, 			None, chi2amp)
chi2amp.PyMinuitFunction	= MethodType(	pyMinuitFunction,		None, chi2amp)
chi2amp.LoadFit			= MethodType(	loadFit,			None, chi2amp)
chi2amp.SampleAmp		= MethodType(	sampleFits,			None, chi2amp)
chi2amp.SampleParameter		= MethodType(	sampleParameter,		None, chi2amp)
