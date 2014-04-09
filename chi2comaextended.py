import sys
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
import libchi2comapy
from libchi2comapy import chi2coma
from types import MethodType
from convertTextOutput import getComaData


"""Some class method definitions for the chi2 class, that are useful, but not in the C++ code"""


def loadFit(self,waves,direct, flag='PINV',eps=1.E-3):
	"""Loads the waves in 'waves' from the fit in 'direct'. The fit has to have text output"""
	if not len(waves) == self.GetNwaves():
		print "Wrong number of waves"
		return
	nWaves = self.GetNwaves()
	data = getComaData(waves,direct, flagg=flag,eps=1.E-3)
	pointsSorted=[]
	for i in range(nWaves):
		self.SetWaveName(i,waves[i].strip())
		pointsLine=[]
		for j in range(nWaves):
			pointsLine.append([])
		pointsSorted.append(pointsLine)
	comaSorted=[]
	for i in range(nWaves**2):
		comaSortedLine=[]
		for j in range(nWaves**2):
			comaSortedLine.append([])
		comaSorted.append(comaSortedLine)
	for bin in range(len(data[0])):
		for i in range(nWaves):
			for j in range(nWaves):
				ii = i*nWaves+j
				pointsSorted[i][j].append((data[0][bin][ii]).real)
				for k in range(nWaves):
					for l in range(nWaves):
						jj = k*nWaves + l
						comaSorted[ii][jj].append((data[1][bin][ii][jj]).real)
	for i in range(nWaves):
		for j in range(nWaves):
			self.SetData(i,j,pointsSorted[i][j])
			self.SetError(i,j,comaSorted[i*nWaves+j][i*nWaves+j])
	for ii in range(nWaves**2):
		i1 = int(ii/nWaves)
		j1 = ii - i1*nWaves
		for jj in range(nWaves**2):
			i2 = int(jj/nWaves)
			j2 = jj - i2*nWaves
			self.SetComa(i1,j1,i2,j2,comaSorted[ii][jj])

			
	


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


def getPlot(self,i,j,nBins=0):
	"""
	Returns the plots for data and fit for wave i interfering with j
	The format is: [[mass,binwidths,data,errors],[mass,fit]]
	"""
	if nBins ==0:
		nBins = self.GetNbins()
	binning = self.GetBinning()
	width = (binning[-1]-binning[0])/nBins
	amp1 = self.GetPoints(i,nBins)
	amp2 = self.GetPoints(j,nBins)
	fit=[]
	mFit=[]
	mDat=[]
	dat=self.GetData(i,j)
	err=self.GetError(i,j)
	errm = []
	for bin in range(nBins):
		mFit.append(binning[0] + (bin+0.5)*width)
		mDat.append((binning[bin] + binning[bin+1])/2)
		errm.append((binning[bin+1]-binning[bin])/2)
		if i==j:
			fit.append(abs(amp1[bin])**2)
		if i<j:
			fit.append((amp1[bin]*amp2[bin].conjugate()).imag)
		if j<i:
			fit.append((amp1[bin]*amp2[bin].conjugate()).real)
	return [[mDat,errm,dat,err],[mFit,fit]]



def pyMinuitFunctionComa(self, name='Chi2'):
	"""
	Generates a function, that can be used by PyMinuit (every parameter defined explicitely)
	The Chi2 instance has to be called 'name' to work
	The function to be called by PyMinuit is then 'Chi2PM'
	"""
	parNames = self.GetParNames()
	parString= ', '.join(parNames)
	command = "def Chi2PM("+parString+"): return "+name+".EvalComa(["+parString+"])"
	return command

# Add methods to chi2coma
chi2coma.RelParameterName 	= MethodType(	releaseParameterByName,		None, chi2coma)
chi2coma.FixParameterName 	= MethodType(	fixParameterByName, 		None, chi2coma)
chi2coma.SetParameterByName	= MethodType(	setParameterByName, 		None, chi2coma)
chi2coma.GetPlot 		= MethodType(	getPlot, 			None, chi2coma)
chi2coma.PyMinuitFunctionComa	= MethodType(	pyMinuitFunctionComa,		None, chi2coma)
chi2coma.LoadFit		= MethodType(	loadFit,			None, chi2coma)


