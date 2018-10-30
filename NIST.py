from __future__ import division

import numpy as np
import requests
from matplotlib import pyplot as plt
import os.path
import os
import warnings
import collections

class NIST:
	Candidate_Components_In = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al" , "Si", "P", "S", "Cl", "Ar", "K", "Ca",
		"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", 
		"Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", 
		"Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "CO2", "CO",
		"O2", "N2", "H2"];
	'''Candidate_Components_In = ["Hg","Na","H","N","Ne","Ar","C","Xe"];'''
	Candidate_Components = [];
	
	Candidate_Components_Data = {};
	verbose = False;
	NM_Peak_Allowance = 0.5;
	Peak_Min_RelMag = 0.1;
	
	def log(self, str):
		if(self.verbose):
			print str;
	
	@staticmethod
	def toFloat(str):
		StringStrip = ''.join(c for c in str if (c.isdigit() or c == '.'))
		if len(StringStrip) == 0:
			return 0;
		return float(StringStrip);
	
	def download_nist_array(self, Type, Boundary_Low, Boundary_High):
		#Check if we have it on file
		if(not os.path.isfile("Nist/"+Type+".dat")): 

			Response = requests.get("http://physics.nist.gov/cgi-bin/ASD/lines1.pl", params={"spectra" : Type, "low_wl" : str(Boundary_Low), "upp_wl" : str(Boundary_High), "unit" : "1", "line_out" : "0", "remove_js" : "on", "order_out" : "0", "output" : "0", "show_obs_" : "1", "show_av" : "2", "page_size" : "15", "show_obs_wl" : "1", "show_av" : "2", "tsb_value" : "0", "A_out" : "0", "intens_out": "on", "allowed_out" : "1", "forbid_out" : "1", "conf_out" : "on", "term_out" : "on", "enrg_out" : "on", "de" : "0", "format" : "1", "en_unit" : "1"})
			Raw_Response = Response.text;
			Table_Text = Raw_Response[Raw_Response.find("<pre>") + 5:Raw_Response.find("</pre>")];
			Table_Text_NoEnd = Table_Text[:Table_Text.rfind("\n", 10)];
			Table_Text_NoEnd = Table_Text_NoEnd[:Table_Text_NoEnd.rfind("\n", 10)];
			
			#Table_Text_NoEnd = os.linesep.join([s for s in Table_Text_NoEnd.splitlines() if s])
			if(len(Table_Text_NoEnd) < 5):
				return [], [];
			tmp_file = open("Nist/" + Type + ".dat", "w+");
			tmp_file.write(Table_Text_NoEnd);
			tmp_file.close();
			self.log("Downloaded data for candidate: " + Type);
		try:
			with warnings.catch_warnings():
				warnings.simplefilter("ignore");
				return np.genfromtxt("Nist/" + Type + ".dat", unpack = True, autostrip=True, loose=True, invalid_raise=False, delimiter = "|", filling_values=0, skip_header=7, usecols=np.array([0,1]), converters={0: NIST.toFloat, 1: NIST.toFloat});
		except IndexError:
			self.log( "Error loading.");
			return [], [];

	def __init__(self, verbose = False,lowLimit=300,highLimit=900,NM_Peak_Allowance=0.5,Peak_Min_RelMag=0.2):
	
		self.Peak_Min_RelMag = Peak_Min_RelMag;
		self.NM_Peak_Allowance = NM_Peak_Allowance;
		self.verbose = verbose;
		
		self.log("Downloading data");
		
		Numerics = ["","I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"];
		
		for Component in self.Candidate_Components_In:
			for Numeric in Numerics:
				Name = Component + " " + Numeric;
				self.Candidate_Components.append(Name);
		
		#Download data for all candidate components
		for Candidate_Component in self.Candidate_Components:
			try:
				WLs, RelInt = self.download_nist_array(Candidate_Component, lowLimit, highLimit);
			except ValueError:
				self.log("Candidate invalid: " + Candidate_Component);
				continue;

			if(isinstance(WLs, collections.Iterable)):
				try:
					if(len(WLs) > 0):
						#Store it!
						self.Candidate_Components_Data[Candidate_Component] = {};
						self.Candidate_Components_Data[Candidate_Component]["WLs"] = WLs;
						self.Candidate_Components_Data[Candidate_Component]["RelInt"] = RelInt;
				
						self.log("Verified data for candidate: " + Candidate_Component);
					else:
						self.log("Candidate invalid: " + Candidate_Component);
				except TypeError:
					self.log("Candidate invalid: " + Candidate_Component);
			else:
				self.log("Candidate invalid: " + Candidate_Component);
					
		self.log("Download done...");
		
	def Peaks(self, Wavelengths, RelNormMag):
		#Find peaks in this data
		peaksWLs = np.array([]); #Wavelengths at which peaks occur

		for i in np.arange(1, len(Wavelengths) - 1, 1):
			WL = Wavelengths[i];
			#Is there a PEAK at this wavelength?
			if(RelNormMag[i] > self.Peak_Min_RelMag):
				#if(RelNormMag[i] > RelNormMag[i-1] and RelNormMag[i] > RelNormMag[i+1]):
				#	peaksWLs = np.append(peaksWLs, WL);
				
				#Walk left until we walk non-horizontal
				CurPos = i;
				LeftDown = False;
				RightDown = False;
				NM_Per_Step = abs(Wavelengths[1] - Wavelengths[0]);
				Steps_Per_NM = 1 / NM_Per_Step;
				Steps_Peak_Allowance = int(Steps_Per_NM * self.NM_Peak_Allowance);
				while(True):
					CurPos -= 1;
					if(RelNormMag[CurPos] > RelNormMag[i]):
						break;
					if(RelNormMag[CurPos] < RelNormMag[i]):
						#Does it not go up right after?
						Nope = False;
						for step in np.arange(1, Steps_Peak_Allowance, 1):
							if(CurPos - step < 0):
								break;
							if(RelNormMag[CurPos - step] > RelNormMag[i]):
								Nope = True;
								break;
							
						if(Nope):
							break;
						LeftDown = True;
						break;
					if(CurPos == 0):
						break;
						
				while(True):
					CurPos += 1;
					if(RelNormMag[CurPos] > RelNormMag[i]):
						break;
					if(RelNormMag[CurPos] < RelNormMag[i]):
						Nope = False;
						for step in np.arange(1, Steps_Peak_Allowance, 1):
							if(CurPos + step >= len(Wavelengths)):
								break;
							if(RelNormMag[CurPos + step] > RelNormMag[i]):
								Nope = True;
								break;
							
						if Nope:
							break;
						RightDown = True;
						break;
					if(CurPos == len(Wavelengths) - 1):
						break;
				
				if(LeftDown and RightDown):
					peaksWLs = np.append(peaksWLs, i);
								
		return peaksWLs;
		
	def ComponentsInSpectrum(self, Wavelengths, RelMag):
		#Now try to match peaks.
		
		#Find peaks
		Peaks = self.Peaks(Wavelengths, RelMag); #Find the wavelengths at which peaks occur
							
		ComponentsInSpectrum = [];
		PeaksCount = [];
		
		#Now look through known data if there are peaks at the same wavelengths.
		for Candidate in self.Candidate_Components_Data.keys():
			C_WLs 		= self.Candidate_Components_Data[Candidate]["WLs"];
			
			C_WLs_Mask = np.logical_and(C_WLs >= np.min(Wavelengths), C_WLs <= np.max(Wavelengths));
			C_WLs = C_WLs[C_WLs_Mask];
						
			C_RelInt 	= self.Candidate_Components_Data[Candidate]["RelInt"][C_WLs_Mask];
						
			if(not isinstance(C_RelInt, collections.Iterable) or len(C_RelInt) == 0 or np.max(C_RelInt) == 0):
				continue; #Nope, invalid data
				
			C_RelInt /= np.max(C_RelInt);
			
			#How many peaks match?
			PeakCounter = 0;
			for PeakIndex in Peaks:
				PeakWL = Wavelengths[PeakIndex];
				Mask = np.logical_and(np.logical_and(PeakWL > C_WLs - self.NM_Peak_Allowance, PeakWL < C_WLs + self.NM_Peak_Allowance), C_RelInt >= self.Peak_Min_RelMag);
			
				if(len(C_WLs[Mask]) > 0):
					PeakCounter += 1;
			
			#Now, how many peaks match?
			if(PeakCounter > 0):
				self.log("Peaks match for " + Candidate + ": " + str(PeakCounter));
				ComponentsInSpectrum.append(Candidate);
				PeaksCount.append(PeakCounter);
		
		ReturnDict = {};
		
		for i in np.arange(0, len(ComponentsInSpectrum)):
			Comp = ComponentsInSpectrum[i];
			Matches = PeaksCount[i];
			
			#We need to check if this candidate is viable.
			#How many peaks are there in the NIST data for this candidate?
			C_WLs 		= self.Candidate_Components_Data[Comp]["WLs"];
			
			C_WLs_Mask = np.logical_and(C_WLs > np.min(Wavelengths), C_WLs < np.max(Wavelengths));
			C_WLs = C_WLs[C_WLs_Mask];
			
			C_RelInt 	= self.Candidate_Components_Data[Comp]["RelInt"][C_WLs_Mask];
			
			if(np.max(C_RelInt) == 0):
				break; #Nope, invalid data
				
			C_RelInt /= np.max(C_RelInt);
			
			C_RelInt_Mask = C_RelInt > self.Peak_Min_RelMag;
			
			Amount_Peaks = len(C_RelInt[C_RelInt_Mask]);
			
			Percentage_Peaks_In_Measurement = Matches / Amount_Peaks;
			
			ReturnDict[Comp] = (Matches, Percentage_Peaks_In_Measurement);
			
		return ReturnDict;