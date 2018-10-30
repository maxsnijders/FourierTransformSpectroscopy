from __future__ import division

import numpy as np
import NIST
from matplotlib import pyplot as plt
from collections import OrderedDict
import colors as col
import astronomy as astr
from decimal import *

def is_array(variable):
	try:
		length = len(variable);
		return True, length;
	except TypeError:
		pass;
		
	try:
		length = variable.shape;
			
		if(len(length)):
			return True, length;
			
		else:
			return False, 0;
	except AttributeError:
		return False, 0;

class Spectrum:
	CalibratedWavelength 	= 0;
	LaserRelMag 			= 0;
	SourceRelMag 			= 0;
	NIST					= 0;
	LaserMask				= 0;
	Component_Frac_Min		= 0;
	
	def plot_range(self, outFile, title = "Spectrum", printPeaks = True, printComponents = True, printPlanck = False):
	
		fig = plt.figure();
		ax = fig.add_subplot(111);
		#ax.plot(self.CalibratedWavelength[Mask], self.LaserRelMag[Mask] /  np.max(self.LaserRelMag[Mask]),  label = "Laser");
		#ax.plot(self.CalibratedWavelength[Mask], self.SourceRelMag[Mask] / np.max(self.SourceRelMag[Mask]), label = "Source");
		
		SourceMinusLaser = self.SourceRelMag / np.max(self.SourceRelMag) - self.LaserRelMag / np.max(self.LaserRelMag);
		MaskedSourceMinusLaser = SourceMinusLaser[self.LaserMask];
		ax.fill_between(self.CalibratedWavelength, 0, SourceMinusLaser / np.max(MaskedSourceMinusLaser), label = "Source - Laser", zorder=0);
		
		ax.set_xlabel("Wavelength [nm]");
		ax.set_xlim([np.min(self.CalibratedWavelength), np.max(self.CalibratedWavelength)]);	
		ax.set_ylim([0, 1.1]);
		ax.set_ylabel("Relative Intensity");
		#ax.set_yscale('log')
		ax.grid(b = True, which='minor');
		ax.grid(b = True, which='major');
		ax.minorticks_on();		
		
		#Print PLANCK
		if(printPlanck):
			print "Plancking..."
			PeakWL = self.CalibratedWavelength[np.where(SourceMinusLaser == np.max(SourceMinusLaser))[0][0]];
			Temp = 2.8977*10**(6) / PeakWL; 
			ax.axvline(PeakWL, label="$\lambda_{max} \Rightarrow T="+str(round(Temp,0))+"K$", c = 'red');
			RelIntsPlanck = self.planckspectrum(Temp);
			ax.plot(self.CalibratedWavelength, RelIntsPlanck / np.max(RelIntsPlanck), label = "Planck curve", zorder = 0.5);
			#title = title + " Temp = " + str(Temp) + "K";
		
			#Now calculate a sensitivity spectrum						
			Sensitivities = (SourceMinusLaser * np.max(RelIntsPlanck)) / (np.max(MaskedSourceMinusLaser) / RelIntsPlanck);
			
			fig2 = plt.figure();
			ax2 = fig2.add_subplot(111);
			ax2.plot(self.CalibratedWavelength, Sensitivities / np.max(Sensitivities), label = "Derived Sensitivity");
			#ax2.plot(self.CalibratedWavelength, 1 / Sensitivities, label = "Derived Correction");
			ax2.set_xlabel("Wavelength [nm]");
			ax2.set_ylabel("Relative Sensitivity");
			ax2.set_xlim([np.min(self.CalibratedWavelength), np.max(self.CalibratedWavelength)]);
			ax2.set_ylim([0, 1.1]);
			ax2.grid(b = True, which='minor');
			ax2.grid(b = True, which='major');
			ax2.minorticks_on();		
			ax2.legend(loc=0);
			
			fig2.savefig(outFile[:-4] + "-sensitivity.png");
			
		ax.set_title(title);
		
		#Print the peaks?
		if(printPeaks):
			PeakIndices = self.NIST.Peaks(self.CalibratedWavelength, SourceMinusLaser / np.max(SourceMinusLaser[self.LaserMask]));
			Heights = np.zeros(len(PeakIndices));
			Peaks = np.zeros(len(PeakIndices));
			
			SMLaserNorm = SourceMinusLaser / np.max(SourceMinusLaser[self.LaserMask]);
			
			for i in np.arange(0, len(PeakIndices), 1):
				Heights[i] = SourceMinusLaser[PeakIndices[i]] / np.max(SourceMinusLaser[self.LaserMask]);
				Peaks[i] = self.CalibratedWavelength[PeakIndices[i]];
								
			ax.scatter(Peaks, Heights, label = "Peaks", c="green", zorder=1);
		if(printComponents):
			Components = self.analyse_spectrum();
			
			for Component in Components.items():
				Name = Component[0];
				Fraction = Component[1][1];
				
				if(Fraction < self.Component_Frac_Min):
					continue; #Not enough present
				#Now get the DATA
				C_WLs 		= self.NIST.Candidate_Components_Data[Name]["WLs"];
			
				C_WLs_Mask = np.logical_and(C_WLs > np.min(self.CalibratedWavelength), C_WLs < np.max(self.CalibratedWavelength));
				C_WLs = C_WLs[C_WLs_Mask];
			
				C_RelInt 	= self.NIST.Candidate_Components_Data[Name]["RelInt"][C_WLs_Mask];
			
				if(np.max(C_RelInt) == 0):
					break; #Nope, invalid data
				
				C_RelInt /= np.max(C_RelInt);
				
				
				MaxWLIndex = np.where(C_RelInt == np.max(C_RelInt));
				MaxWL = C_WLs[MaxWLIndex[0][0]];
				
				r,g,b = col.compute_rgb(MaxWL);
				c = col.rgb_to_hex(r * 255, g * 255, b * 255);
			
				percentage = Component[1][1] * 100;
				ax.scatter(C_WLs, C_RelInt, label = "NIST " + Name, c = c, zorder=2);
		
		ax.legend(loc=0);
		fig.savefig(outFile);
		
	def __init__(self, UncalibratedFrequency, LaserRelMag, SourceRelMag, lowLimit = 300, highLimit = 900, withoutLaser = False,
				laserBuffer = 0.2, Component_Frac_Min = 0.2, NM_Peak_Allowance = 0.5, Peak_Min_RelMag=0.2):
				
		#Calibrate the wavelength...
		MaxLaserRelMagVal = np.max(LaserRelMag[100:]);
		MaxLaserRelMagInd = np.where(LaserRelMag == MaxLaserRelMagVal);
		LaserPeakUncalFreq = UncalibratedFrequency[MaxLaserRelMagInd[0][0]];
		CalibratedWavelength = 632.82 * LaserPeakUncalFreq / UncalibratedFrequency;

		#Filter wavelengths!
		WLMask = np.logical_and(CalibratedWavelength >= lowLimit, CalibratedWavelength <= highLimit);
				
		laser_high = 632.82 + laserBuffer; 
		laser_low =  632.82  - laserBuffer;
		
		if(withoutLaser):
			#Filter out the laser as well.
			#WLMask = np.logical_and(WLMask, np.logical_or(CalibratedWavelength >= laser_high, CalibratedWavelength <= laser_low));
			self.LaserMask = np.logical_or(CalibratedWavelength[WLMask] >= laser_high, CalibratedWavelength[WLMask] <= laser_low)
		else:
			self.LaserMask = np.zeros(len(CalibratedWavelength[WLMask]));
			self.LaserMask.fill(1);
		
		self.CalibratedWavelength 	= CalibratedWavelength[WLMask];
		self.LaserRelMag 			= LaserRelMag[WLMask];
		self.SourceRelMag 			= SourceRelMag[WLMask];
		
		self.Component_Frac_Min = Component_Frac_Min;
		
		self.NIST = NIST.NIST(verbose = False, lowLimit = lowLimit, highLimit = highLimit, Peak_Min_RelMag=Peak_Min_RelMag, 
								NM_Peak_Allowance=NM_Peak_Allowance);
		
	def analyse_spectrum(self):
		#Now try to look for components using this MASK.
		SourceMinusLaser = self.SourceRelMag / np.max(self.SourceRelMag[self.LaserMask]) - self.LaserRelMag / np.max(self.LaserRelMag[self.LaserMask]);
		Components = self.NIST.ComponentsInSpectrum(self.CalibratedWavelength, SourceMinusLaser / np.max(SourceMinusLaser[self.LaserMask]));
		
		#Sort them
		Components = OrderedDict(sorted(Components.items(), key=lambda i:i[1], reverse = True));
		return Components;
		
	def planckspectrum(self, Temp):
		#return 2 * astr.h * astr.c**2 / (self.CalibratedWavelength * (np.exp(astr.h * astr.c / (self.CalibratedWavelength * astr.kb * Temp)) - 1));
		specrads = np.zeros(len(self.CalibratedWavelength));
		
		for i in np.arange(0, len(self.CalibratedWavelength), 1):
			item_wl = self.CalibratedWavelength[i];
			
			firstfrac 	= Decimal(2) * Decimal(astr.h) * Decimal(astr.c)**2  / (Decimal(item_wl)**5);
			exponent  	= Decimal(astr.h) * Decimal(astr.c) / (Decimal(item_wl * 10**(-9)) * Decimal(astr.kb) * Decimal(Temp));
			secfrac	  	= Decimal(1) / (Decimal(astr.e) ** exponent - 1);
			specrad		= firstfrac * secfrac;
		
			specrads[i] = float(specrad);
			
		return specrads/np.max(specrads);
		
		