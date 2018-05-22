def calcSysParams(eG):
	print('Gamma Energy,', eG, 'MeV')
	#,BW,Ng,eE,grad,shuntZ,repRate,laserWavelength,eLaser
	return  'yo'

def calcIntLaser(wavelength, pulseEnergy, repRate, efficiency):
	laserPhotonEnergy = 4.1357e-15 * 3.0e8 * 1.0e9 / wavelength # laser photon energy (eV) with wavelength (nm)
	avgLaserPower = pulseEnergy * repRate # W from J * Hz
	nLaser = avgLaserPower * 6.242e18 * efficiency / laserPhotonEnergy
	print('-- Interaction Laser Parameters --')
	print('Pulse Energy (J):', pulseEnergy)
	print('Repetition Rate (Hz):', repRate)
	print('Wavelength (nm):', wavelength)
	print('Average Laser Power (W):',avgLaserPower)
	print('Number of Laser Photons (#/s):', '{:.2e}'.format(nLaser))
	#df = pd.DataFrame(
	#	{"parameter" : ['Pulse Energy (J)','Repetition Rate (Hz)', 'Wavelength (nm)'],
	#	"value" : [pulseEnergy, repRate, wavelength]})

def calcXband(grad,eMax,repRate,shuntZ,averageBeamCurrent):
	nMicroBunch = 1000 #assume 1000 microbunches
	f = 11.424e9 #xband (Hz)
	beamPulseWidth = nMicroBunch/f*1e6 #(us)
	rfRiseFallTime = 0.5 #(us)
	pulseWidth = beamPulseWidth + rfRiseFallTime #total pulse width (us)
	chargePerRepScan = averageBeamCurrent * 1e-6 * 1e12 / repRate #pc #np.linspace(1,100*1000,1000) #(pC)
	chargePerMicrobunch = chargePerRepScan / 1000 #pC
	dutyFactor = pulseWidth * 1e-6 * repRate
	peakBeamCurrent = chargePerRepScan / beamPulseWidth /1e3 #(mA)
	#averageBeamCurrent = chargePerRepScan * 1e-12 * repRate * 1e6 #(uA)
	eg,ig = np.meshgrid(electronEnergyScan,peakBeamCurrent)
	#shuntImpedance = 160 #(Mohms m) *2.2 LN2 cooling, *4.6 He.
	peakBeamPower = eMax*peakBeamCurrent/1e3
	systemLength = eMax/grad
	peakStructurePower = eMax**2/systemLength/shuntZ
	rfPowerMarginFactor = 1.1
	totalPeakRFPower = (peakBeamPower + peakStructurePower)*rfPowerMarginFactor
	totalAverageRFPower = totalPeakRFPower * dutyFactor *1e3
	nElectrons = averageBeamCurrent*1e-6/1.602e-19
	print('-- X-Band Accelerator Parameters --')
	print('Electron Energy (MeV):', eMax)
	print('Gradient (MV/m):', grad)
	print('Repetition Rate (Hz):', repRate)
	print('Average Current (uA):', averageBeamCurrent)
	print('Number of Electrons (#/s):', '{:.2e}'.format(nElectrons))
	print('Average RF Power (kW):',totalAverageRFPower)
	print('Peak RF Power (MW):',totalPeakRFPower)
	print('Accelerator Length (m):', systemLength)
	print('chargePerMicrobunch (pC):', chargePerMicrobunch)