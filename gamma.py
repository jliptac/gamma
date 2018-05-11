import numpy as np
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
from matplotlib.ticker import FormatStrFormatter

thomsonCrossSection = 0.63e-28 #(m2)
spotSizeDia = 30 #(um)
spotSizeArea = math.pi*(spotSizeDia * 1.0e-6 /2.0) ** 2 #(m2)

#interaction laser parameters
laserWavelength = np.array([1030, 1030/2, 1030/3, 10]) #laser wavelength (nm)
laserPhotonEnergy = 4.1357e-15 * 3.0e8 * 1.0e9 / laserWavelength # laser photon energy (eV)

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

#plt.figure(5)
#plt.loglog(BWscan,Nph)

a0 = 0.03
#ocheck = 4*gammaE*laserPhotonEnergy/0.511e6
#print('ochec',ocheck)
print('gammaE',gammaE)
#gammaPhotonEnergy = 2 * gammaE** 2 * (1+np.cos(thetaL)) * laserPhotonEnergy / (1 + (gammaE * scatteringAngle)**2 + a0**2 + 4*gammaE*laserPhotonEnergy/0.511e6)
#gammaPhotonEnergy = 2 * gammaE** 2 * (1+np.cos(thetaL)) * laserPhotonEnergy / (1 + np.power(np.multiply(gammaE,scatteringAngle),2) + a0**2 + 4*gammaE*laserPhotonEnergy/0.511e6)

sns.set(color_codes=True)
colors = (plt.rcParams['axes.prop_cycle'].by_key()['color'])
#print(colors[1])

gammaPhotonEnergy = np.zeros([len(gammaE),len(scatteringAngle)])
BW = np.zeros([len(gammaE),len(thetaRMS)])
nGammaBW = np.zeros([len(gammaE),len(thetaRMS)])
labelNames = ['250 MeV','400 MeV','525 MeV']
plt.figure(1)
for i in range(len(gammaE)):
	gammaPhotonEnergy[i,:] = 2 * gammaE[i]** 2 * (1+np.cos(thetaL)) * laserPhotonEnergy[2] / (1 + np.power(np.multiply(gammaE[i],scatteringAngle),2) + a0**2 + 4*gammaE[i]*laserPhotonEnergy[2]/0.511e6)
	plt.plot(scatteringAngle*1000,gammaPhotonEnergy[i]/1e6, label = labelNames[i])
	#plt.axvline(x=1/gammaE[i]*1000*.0316, ls='--', alpha=0.5, color=colors[i])
	#plt.axvline(x=-1/gammaE[i]*1000*.0316, ls='--', alpha=0.5, color=colors[i])
plt.xlabel('Scattering Angle (mrad)')
plt.ylabel('Photon Energy (MeV)')
plt.title('Photon Energy vs Scattering Angle ($\lambda_L$ = 343 nm)')
plt.legend()
plt.xlim([-5,5])
plt.ylim([0,10])


electronEnergy = np.linspace(0.5, 800, 1000)
gammaE = electronEnergy/0.511
gammaPhotonEnergyV = np.zeros([len(laserPhotonEnergy),len(gammaE)])
plt.figure(2)
fig,axPrime =plt.subplots()
for i in range(len(laserPhotonEnergy)):
	gammaPhotonEnergyV[i,:] = 2 * np.power(gammaE,2) * (1+np.cos(thetaL)) * laserPhotonEnergy[i] / (1 + np.power(np.multiply(gammaE,0),2) + a0**2 + 4*gammaE*laserPhotonEnergy[i]/0.511e6)
plt.plot(electronEnergy,gammaPhotonEnergyV[0]/1e6, color=colors[2])
plt.plot(electronEnergy,gammaPhotonEnergyV[1]/1e6, color=colors[1])
plt.plot(electronEnergy,gammaPhotonEnergyV[2]/1e6, color=colors[0])
plt.plot(electronEnergy,gammaPhotonEnergyV[3]/1e6, color=colors[3], ls='--')
plt.xlabel('Electron Energy (MeV)')
plt.ylabel('Photon Energy (MeV)')
plt.title('Photon Energy vs Electron Energy')
plt.legend(['$\lambda_L$ = 1030 nm','$\lambda_L$ = 515 nm','$\lambda_L$ = 343 nm','$\lambda_L$ = 10 nm'])
plt.axhline(y=10, ls='--', alpha=0.75, color='black')
plt.xlim([0,800])
plt.ylim([0,20])


laserPowerScan = np.linspace(0.001,10000,10000)
laserWavelength = 1030/2 #laser wavelength (nm)
laserPhotonEnergy = 4.1357e-15 * 3.0e8 * 1.0e9 / laserWavelength # laser photon energy (eV)
laserFreqMultEff = 0.5
nLaserPhotonScan = laserPowerScan * 6.242e18 *laserFreqMultEff / laserPhotonEnergy
plt.figure(4)
fig,axPrime2 =plt.subplots()
axTopNe= axPrime2.twiny()
axRightNl= axPrime2.twinx()
#xPrimeTicks = axPrime2.get_xticks()
#xTopTicks = xPrimeTicks/100
#axTopNe.set_xticks(xTopTicks)

BW = 1e-3
neGamma1e12 = 1.0e12*spotSizeArea/thomsonCrossSection/nLaserPhotonScan/BW*1.602e-19
neGamma1e13 = 1.0e13*spotSizeArea/thomsonCrossSection/nLaserPhotonScan/BW*1.602e-19
neGamma1e14 = 1.0e14*spotSizeArea/thomsonCrossSection/nLaserPhotonScan/BW*1.602e-19
axPrime2.loglog(neGamma1e12, laserPowerScan)
axPrime2.loglog(neGamma1e13, laserPowerScan)
axPrime2.loglog(neGamma1e14, laserPowerScan)
axPrime2.loglog(nLaserPhotonScan*1.602e-19, laserPowerScan, color='black', alpha=0.5, ls='--')
axPrime2.set_ylim([1e-2,1e4])
axPrime2.set_xlim([1e-6,1])
axPrime2.set_xlabel('Average Electron Current (A) $\propto N_e$')
axPrime2.set_ylabel('Average Laser Power (W) $\propto N_L$')
#axPrime2.title('Idealized Photon Output ($\lambda_L$ = 515 nm)')
axPrime2.legend(['N$_\gamma = 10^{12}$ (1/s)','N$_\gamma = 10^{13}$ (1/s)','N$_\gamma = 10^{14}$ (1/s)'])

axTopNe.set_xlabel(r'Number of Electrons ($N_e$)')
axTopNe.set_xscale('log')
#axTopNe.set_xlim(axPrime2.get_xlim())
#axPrimelocs,axPrimeLabels = axPrime2.get_xticks()
# xTopTicks = xPrimeTicks/1.602e-19
# print(xPrimeTicks)
# #axTopNe.set_xticks(axPrimelocs,axPrimeLabels)
# #axTopNe.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# axTopNe.grid(False)

# axRightNl.set_ylabel(r'Number of Laser Photons ($N_L$)')
# axRightNl.grid(False)
# axPrime2.set_xlim([1e-5,1])




#axTopNe.set_xlim([1e-5,1])


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
shuntImpedance = 160 #(Mohms m)
peakBeamPower = eg*ig/1e3
systemLength = eg/gradient[2]
peakStructurePower = eg*eg/systemLength/shuntImpedance
rfPowerMarginFactor = 1.1
totalPeakRFPower = (peakBeamPower + peakStructurePower)*rfPowerMarginFactor
totalAverageRFPower = totalPeakRFPower * dutyFactor


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
# axTop.grid(False)
# axTop.set_xlabel(r'Accelerator Length (m) with $\nabla$= 100 MeV/m')
# axPrime.grid(True, color='#95a5a6')
# axPrime.set_facecolor('#ffffff')
# axPrime.axvline(x=525, ls='--', alpha=0.8, color=colors[0])
# axPrime.axvline(x=300, ls='--', alpha=0.8, color=colors[0])
# plt.ioff()

#calculate system params
#Ee = 
#laserWavelength = 
#electronCurremt = 

#test git
plt.show()