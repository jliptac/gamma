import numpy as np
import pandas as pd
import math
import itertools

#parameters to iterate over (system, laser, linac)
gammaE = [5.0,10.0,20.0] # output gamma energy in MeV
gammaN = [1e12,1e13,1e14] # number of output gammas
gammaBW = [1e-3] # fractional bandwidth
repRate = [1000.0] # system reprate in Hz
laserWavelength = [343,515] # interaction laser wavelength in nm
laserPulseEnergy = [1.0, 2.0] # interaction laser pulse energy in J
accShuntZ = [160, 160*2.2] # linac shunt impedance in Mohms m
accGrad = [100,150] # linac gradient in MV/m

#linac params
nMicroBunch = 1000 #assume 1000 microbunches
f = 11.424e9 #xband (Hz)
beamPulseWidth = nMicroBunch/f*1e6 #(us)
rfRiseFallTime = 0.5 #(us)
pulseWidth = beamPulseWidth + rfRiseFallTime #total pulse width (us)

#interaction region params
thomsonXS = 0.63e-28 #(m2)
spotSizeDia = 30 #(um)
spotSizeArea = math.pi*(spotSizeDia * 1.0e-6 /2.0) ** 2 #(m2)

x=np.zeros(len(gammaE)*len(gammaN)*len(gammaBW)*len(repRate)*len(laserWavelength)*len(laserPulseEnergy)*len(accShuntZ)*len(accGrad))
df = pd.DataFrame({
	'gammaE':x,'gammaN':x,'gammaBW':x,'repRate':x,'laserWavelength':x,'laserPulseEnergy':x,'accShuntZ':x,'accGrad':x,'laserAvgPower':x,'laserHarmonicEff':x,'laserN':x,'accElectronN':x,'accAvgCurrent':x,'accElectronE':x,'accLength':x,'accAvgPower':x,'accPeakPower':x})
df=df[['gammaE','gammaN','gammaBW','repRate','laserWavelength','laserPulseEnergy','accShuntZ','accGrad','laserAvgPower','laserHarmonicEff','laserN','accElectronN','accAvgCurrent','accElectronE','accLength','accAvgPower','accPeakPower']]
index = 0
for p in itertools.product(gammaE,gammaN,gammaBW,repRate,laserWavelength,laserPulseEnergy,accShuntZ,accGrad):
	df.loc[index,:8] = p
	# interaction laser calculations
	laserPhotonEnergy = 4.1357e-15 * 3.0e8 * 1.0e9 / df.at[index,'laserWavelength'] # laser photon energy (eV) with wavelength (nm)
	if df.at[index,'laserWavelength'] == 343:
		df.at[index,'laserHarmonicEff'] = 0.4
	elif df.at[index,'laserWavelength'] == 515:
		df.at[index,'laserHarmonicEff'] = 0.8
	else:
		df.at[index,'laserHarmonicEff'] = 0.5
	df.at[index,'laserAvgPower'] = df.at[index,'laserPulseEnergy'] * df.at[index,'repRate'] # W from J * Hz
	df.at[index,'laserN'] = df.at[index,'laserAvgPower'] * 6.242e18 * df.at[index,'laserHarmonicEff'] / laserPhotonEnergy # number of interaction photons
	#calc electron params needed for single pass system
	df.at[index,'accElectronN'] = df.at[index,'gammaN'] * spotSizeArea / df.at[index,'laserN'] / thomsonXS / df.at[index,'gammaBW']
	df.at[index,'accAvgCurrent'] = df.at[index,'accElectronN']*1.602e-19*1e6 # uA
	df.at[index,'accElectronE'] = 0.511/2*np.sqrt(df.at[index,'gammaE']/laserPhotonEnergy*1.e6) # approax max electron in MeV
	df.at[index,'accLength'] = df.at[index,'accElectronE']/df.at[index,'accGrad']
	dutyFactor = pulseWidth * 1e-6 * df.at[index,'repRate']
	chargePerRep = df.at[index,'accAvgCurrent'] * 1e-6 * 1e12 / df.at[index,'repRate']
	peakStructurePower = df.at[index,'accElectronE']**2/df.at[index,'accLength']/df.at[index,'accShuntZ']
	peakBeamCurrent = chargePerRep / beamPulseWidth / 1e3 #(mA)
	peakBeamPower = df.at[index,'accElectronE']*peakBeamCurrent/1e3
	rfPowerMarginFactor = 1.1
	df.at[index,'accPeakPower'] = (peakBeamPower + peakStructurePower)*rfPowerMarginFactor #MW
	df.at[index,'accAvgPower'] = df.at[index,'accPeakPower'] * dutyFactor *1e3 #kW
	index = index+1
df=df[['gammaE','gammaN','gammaBW','repRate','laserWavelength','laserN','laserPulseEnergy','laserAvgPower','laserHarmonicEff','accElectronE','accAvgCurrent','accLength','accAvgPower','accPeakPower','accElectronN','accGrad','accShuntZ']]
df.to_csv('out.csv')
print(df.head())

