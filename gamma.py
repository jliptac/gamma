import numpy as np
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

#calc inverse compton scattering parameters and generate relevant plots
#created by jliptac on 05.10.2018 

thomsonCrossSection = 0.63e-28 #(m2)
spotSizeDia = 30 #(um)
spotSizeArea = math.pi*(spotSizeDia * 1.0e-6 /2.0) ** 2 #(m2)

#interaction laser parameters
laserWavelength = np.array([1030/3, 1030/2, 1030]) #laser wavelength (nm)
laserPhotonEnergy = 4.1357e-15 * 3.0e8 * 1.0e9 / laserWavelength # laser photon energy (eV)
a0 = 0.03

#geometry
thetaL = 0.0 
scatteringAngle = np.linspace(-0.005, 0.005, 10000)
gammaE = np.array([200/0.511, 325/0.511, 425/0.511])
approx = 1/gammaE


thetaRMS = np.linspace(1.0e-6, 5.0e-3, 100000)
BWscan = gammaE[2]**2 * thetaRMS**2
Ne = 10e-3/1.602e-19
Nl = 7e18
Nph = Ne*Nl*thomsonCrossSection/spotSizeArea*BWscan
BW = np.zeros([len(gammaE),len(thetaRMS)])
nGammaBW = np.zeros([len(gammaE),len(thetaRMS)])

#print('gammaE',gammaE)
#gammaPhotonEnergy = 2 * gammaE** 2 * (1+np.cos(thetaL)) * laserPhotonEnergy / (1 + (gammaE * scatteringAngle)**2 + a0**2 + 4*gammaE*laserPhotonEnergy/0.511e6)
#gammaPhotonEnergy = 2 * gammaE** 2 * (1+np.cos(thetaL)) * laserPhotonEnergy / (1 + np.power(np.multiply(gammaE,scatteringAngle),2) + a0**2 + 4*gammaE*laserPhotonEnergy/0.511e6)

sns.set(color_codes=True)
colors = (plt.rcParams['axes.prop_cycle'].by_key()['color'])

plt.figure(1)
gammaPhotonEnergy = np.zeros([len(gammaE),len(scatteringAngle)])
labelNames = ['250 MeV','325 MeV','425 MeV']
for i in range(len(gammaE)):
	gammaPhotonEnergy[i,:] = 2 * gammaE[i]** 2 * (1+np.cos(thetaL)) * laserPhotonEnergy[0] / (1 + np.power(np.multiply(gammaE[i],scatteringAngle),2) + a0**2 + 4*gammaE[i]*laserPhotonEnergy[0]/0.511e6)
	plt.plot(scatteringAngle*1000,gammaPhotonEnergy[i]/1e6, label = labelNames[i])
	#plt.axvline(x=1/gammaE[i]*1000*.0316, ls='--', alpha=0.5, color=colors[i])
	#plt.axvline(x=-1/gammaE[i]*1000*.0316, ls='--', alpha=0.5, color=colors[i])
plt.xlabel('Scattering Angle (mrad)')
plt.ylabel('Photon Energy (MeV)')
plt.title('Photon Energy vs Scattering Angle ($\lambda_L$ = 343 nm)')
plt.legend()
plt.xlim([-5,5])
plt.ylim([0,10])

fig,axPrime =plt.subplots()
electronEnergy = np.linspace(0.5, 800, 1000)
gammaE = electronEnergy/0.511
gammaPhotonEnergyV = np.zeros([len(laserPhotonEnergy),len(gammaE)])
for i in range(len(laserPhotonEnergy)):
	gammaPhotonEnergyV[i,:] = 2 * np.power(gammaE,2) * (1+np.cos(thetaL)) * laserPhotonEnergy[i] / (1 + np.power(np.multiply(gammaE,0),2) + a0**2 + 4*gammaE*laserPhotonEnergy[i]/0.511e6)
	plt.plot(electronEnergy,gammaPhotonEnergyV[i]/1e6)
plt.xlabel('Electron Energy (MeV)')
plt.ylabel('Photon Energy (MeV)')
plt.title('Photon Energy vs Electron Energy')
plt.legend(['$\lambda_L$ = 343 nm','$\lambda_L$ = 515 nm','$\lambda_L$ = 1030 nm'])
plt.axhline(y=10, ls='--', alpha=0.5, color='black')
plt.axhline(y=5, ls='--', alpha=0.5, color='black')
plt.xlim([0,800])
plt.ylim([0,20])
plt.savefig('gammEvsEe.png', bbox_inches='tight')
plt.savefig('gammEvsEe.pdf', bbox_inches='tight')

plt.figure(3)
laserPowerScan = np.linspace(0.001,10000,10000)
laserWavelength = 1030/3 #laser wavelength (nm)
laserPhotonEnergy = 4.1357e-15 * 3.0e8 * 1.0e9 / laserWavelength # laser photon energy (eV)
laserFreqMultEff = 0.4 #3w estimate
nLaserPhotonScan = laserPowerScan * 6.242e18 *laserFreqMultEff / laserPhotonEnergy
BW = 1e-3 #fractional bandwidth
gammaOutput = np.geomspace(1.0e12,1.0e14,num=3)
print(gammaOutput)
for i in range(len(gammaOutput)):
	nElectrons = gammaOutput[i]*spotSizeArea/thomsonCrossSection/nLaserPhotonScan/BW*1.602e-19
	plt.loglog(nElectrons, laserPowerScan)
plt.loglog(nLaserPhotonScan*1.602e-19, laserPowerScan, color='black', alpha=0.5, ls='--')
plt.xlim([1e-6,1])
plt.ylim([1e-2,1e4])
plt.xlabel('Average Electron Current (A) $\propto N_e$')
plt.ylabel('Average Laser Power (W) $\propto N_L$')
plt.legend(['N$_\gamma = 10^{12}$ (1/s)','N$_\gamma = 10^{13}$ (1/s)','N$_\gamma = 10^{14}$ (1/s)'])
plt.title(r'Number of Output Photons (N$_\gamma$), with $\lambda_L$ = 343 nm')
plt.savefig('gammaOutput.png', bbox_inches='tight')
plt.savefig('gammaOutput.pdf', bbox_inches='tight')

#rf power calculation based on xband system
gradient = np.array([50, 100, 150, 80])
electronEnergyScan = np.linspace(1, 1000, 1000) #MeV
beamPulseWidth = 1000/11.424e9*1e6 #1000 microbunches at x-band (us)
rfRiseFallTime = 0.5 # (us)
pulseWidth = beamPulseWidth + rfRiseFallTime #total pulse width (us)
repRate = 1000 # (Hz)
chargePerRepScan = np.linspace(1,100*1000,1000) #(pC)
dutyFactor = pulseWidth * 1e-6 * repRate
peakBeamCurrent = chargePerRepScan / beamPulseWidth /1e3 #(mA)
averageBeamCurrent = chargePerRepScan *1e-12 * repRate *1e6 #(uA)
eg,ig = np.meshgrid(electronEnergyScan,peakBeamCurrent)
shuntImpedance = 160 #(Mohms m) *2.2 LN2 cooling, *4.6 He.
peakBeamPower = eg*ig/1e3
systemLength = eg/gradient[2]
peakStructurePower = eg*eg/systemLength/shuntImpedance
rfPowerMarginFactor = 1.1
totalPeakRFPower = (peakBeamPower + peakStructurePower)*rfPowerMarginFactor
totalAverageRFPower = totalPeakRFPower * dutyFactor

#uncomment to generate rf power plot with interactive contour labeling
# plt.ion()
# fig,axPrime =plt.subplots()
# im=axPrime.contourf(electronEnergyScan,averageBeamCurrent,totalAverageRFPower, cmap='YlOrRd', alpha=.6,levels=np.arange(0.0,1.5,0.05))
# #for c in im.collections:
# #    c.set_edgecolor('red')
# cbar = fig.colorbar(im, ax=axPrime)
# cbar.set_label('Average RF Power (MW)')
# powerLevels = np.arange(0.1,1.0,.1)
# contours = axPrime.contour(electronEnergyScan,averageBeamCurrent,totalAverageRFPower, colors='black',levels=powerLevels, alpha=0.75)

# plt.clabel(contours, fontsize=10, colors='black', fmt=r'%1.1f', manual=True)
# axPrime.set_xlabel('Electon Energy (MeV)')
# axPrime.set_ylabel('Average Electron Current ($\mu$A) $\propto N_e$')

# axTop= axPrime.twiny()
# xPrimeTicks = axPrime.get_xticks()
# xTopTicks = xPrimeTicks/gradient[2]
# axTop.set_xticks(xTopTicks)
# axTop.xaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
# axTop.grid(False)
# axTop.set_xlabel(r'Accelerator Length (m) with $\nabla$ = 100 MeV/m')
# axPrime.grid(True, color='#95a5a6')
# axPrime.set_facecolor('#ffffff')
# axPrime.axvline(x=525, ls='--', alpha=0.5, color='black')
# axPrime.axvline(x=300, ls='--', alpha=0.5, color='black')
# plt.ioff()
# plt.savefig('rfPower.png', bbox_inches='tight')
# plt.savefig('rfPower.png', bbox_inches='tight')


#plt.show()



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
	#laserParams = [avgLaserPower,nLaser,laserPhotonEnergy]
	return nLaser

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

#calculate and print system params

print('-- ICS System Parameters --')
eG = 5 #max gamma energy (Mev)
Ng = 1e14 #number of output gammas per second
BW = 1e-3 #fractional bandwidth
grad = 100
repRate = 1000
shuntZ = 160
#prescribe gamma and laser params to calc electrons

#maybe need lookup unless solve analytic
calcXband(grad,300,repRate,shuntZ,10)

x = calcIntLaser(343, 1, 1000, 0.5)

print(x)
calcSysParams(5)