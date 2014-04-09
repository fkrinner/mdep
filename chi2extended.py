import libchi2py
from libchi2py import chi2
from types import MethodType


"""Some class method definitions for the chi2 class, that are useful, but not in the C++ code"""


def importDataFile(self,fileName,waves={}):
	""" 
	Imports a data file to the chi2. The mapping 'waves' does:	
	'waves = {1:2,3:4}' => wave 1 from the file will be written to wave 2 in chi2, wave 3 from the file will be written to wave 4 in chi2
	If no such mapping is given => waves = {1:1, 2:2, 3:3, ... }
	"""
	if waves == {}:				
		for i in range(self.GetNwaves()):
			waves[i]=i
	inData = open(fileName,'r')
	for line in inData.readlines():
		chunks = line.split('/')
		vals=[]
		nWaves=self.GetNwaves()
		for i in range(3,len(chunks)):
			vals.append(float(chunks[i]))
		try:
			ii = waves[int(chunks[0])]
			jj = waves[int(chunks[1])]
			md = chunks[2]
			if md == 'v' and ii < nWaves and jj < nWaves:
#				print 'Set data '+str(ii)+'/'+str(jj)+' with length '+str(len(vals))
				self.SetData(ii,jj,vals)
			if md == 'e'and ii < nWaves and jj < nWaves:
#				print 'Set error '+str(ii)+'/'+str(jj)+' with length '+str(len(vals))
				self.SetError(ii,jj,vals)
			if md == 'b':
#				print 'Set binning length: '+str(len(vals))
				self.SetBinning(vals)
		except:
			pass


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


def waveIsInTitle(title, waves):
	for wave in waves:
		if wave.strip() in title:
			return wave.strip()
	return "NOT"

def nDig(inString, length):
	while len(inString) < length:
		inString = '0'+inString
	return inString

def getVals(hist, fac=1):
	nBins = hist.GetNbinsX()
	dat=[]
	err=[]
	for bin in range(1,nBins+1):
		dat.append(hist.GetBinContent(bin)*fac)
		err.append(hist.GetBinError(bin))
	return [dat,err]

def load_waves_from_root(self,waves, ROOTfile, nnn=3):
	import ROOT
	from ROOT import TH1D
	import rootpy
	rootpy.log.basic_config_colorized()
	from rootpy.io import root_open
	data=root_open(ROOTfile,"READ")
	if not len(waves) == self.GetNwaves():
		print "Number of waves does not match."
		return
	for i in range(len(waves)):
		namm = waves[i].strip()
		self.SetWaveName(i,namm)
	for name in  data.walk():
		Anames = name
	histNames = Anames[2]
	waveNumbers={}
	for name in histNames:
		if len(name) < 4:
			hi = data.Get(name)
			title = hi.GetTitle()
			wav = waveIsInTitle(title,waves)
			if not wav =="NOT":
				waveNumbers[wav] = nDig(name[1:], nnn)
	for i in range(len(waves)):
		wave=waves[i]
		waveNumber = waveNumbers[wave.strip()]
		histName = 'h'+str(int(waveNumber)) # Name of the intensity histogram
		hi = data.Get(histName)
		vals = getVals(hi)	
		self.SetData(i,i,vals[0])
		self.SetError(i,i,vals[1])
		for j in range(i):
			wave2=waves[j]
			waveNumber2 = waveNumbers[wave2.strip()]	
			reName1 = 'h1'+waveNumber+waveNumber2
			imName1 = 'h2'+waveNumber+waveNumber2
			hi = data.Get(reName1)
			vals = getVals(hi)
			self.SetData(i,j,vals[0])
			self.SetError(i,j,vals[1])
			hi = data.Get(imName1)
			vals = getVals(hi)
			self.SetData(j,i,vals[0])
			self.SetError(j,i,vals[1])
	data.Close()



# Add methods tho chi2
chi2.ImportDataFile 		= MethodType(	importDataFile, 		None, chi2)
chi2.RelParameterName 		= MethodType(	releaseParameterByName,		None, chi2)
chi2.FixParameterName 		= MethodType(	fixParameterByName, 		None, chi2)
chi2.SetParameterByName		= MethodType(	setParameterByName, 		None, chi2)
chi2.GetPlot 			= MethodType(	getPlot, 			None, chi2)
chi2.PyMinuitFunction		= MethodType(	pyMinuitFunction,		None, chi2)
chi2.LoadHfit			= MethodType(	load_waves_from_root,		None, chi2)

