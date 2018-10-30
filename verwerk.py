import numpy as np
import sys
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import requests

BVisMask = (500, 900);
BDoubletMask = (588, 591);

if(len(sys.argv) < 2):
    print "You need to pass at least one argument: the file to read. Optional: nist-ref-file";
    sys.exit(0);
    
File = sys.argv[1];
if len(File) is 0:
	print "Filename invalid";
	sys.exit(0);
	
def firstToFloat(data):
	#StringStrip = data.strip().replace("b","").replace("*","").replace("s","").replace("r","").replace("l","").replace("g","").replace("a","").replace(;
	StringStrip = ''.join(c for c in data if (c.isdigit() or c == '.'))
	if len(StringStrip) == 0:
		return 0;
	return float(StringStrip);

def Nist_Download(Type, Boundary_Low = 300, Boundary_High = 900):
	Response = requests.get("http://physics.nist.gov/cgi-bin/ASD/lines1.pl", params={"spectra" : Type, "low_wl" : str(Boundary_Low), "upp_wl" : str(Boundary_High), "unit" : "1", "line_out" : "0", "remove_js" : "on", "order_out" : "0", "output" : "0", "show_obs_" : "1", "show_av" : "2", "page_size" : "15", "show_obs_wl" : "1", "show_av" : "2", "tsb_value" : "0", "A_out" : "0", "intens_out": "on", "allowed_out" : "1", "forbid_out" : "1", "conf_out" : "on", "term_out" : "on", "enrg_out" : "on", "de" : "0", "format" : "1", "en_unit" : "1"})
	Raw_Response = Response.text;
	Table_Text = Raw_Response[Raw_Response.find("<pre>") + 5:Raw_Response.find("</pre>")];
	Table_Text_NoEnd = Table_Text[:Table_Text.rfind("\n", 10)];
	Table_Text_NoEnd = Table_Text_NoEnd[:Table_Text_NoEnd.rfind("\n", 10)];
	return Table_Text_NoEnd;
	
def Nist_Array(Type, Boundary_Low = 300, Boundary_High = 900):
	Str = Nist_Download(Type, Boundary_Low, Boundary_High);
	tmp_file = open("tmp.txt","w+");
	tmp_file.write(Str);
	tmp_file.close();
	
	return np.loadtxt("tmp.txt", unpack = True, delimiter = "|", skiprows=7, usecols=np.array([1,2]), converters={1: firstToFloat, 2:firstToFloat});

if len(sys.argv) > 2:	
	NIST_WL, NIST_Rel_Int = Nist_Array(sys.argv[2], BVisMask[0], BVisMask[1]);
	
	NIST_VisMask = np.logical_and(NIST_WL > BVisMask[0], NIST_WL < BVisMask[1]);
	NIST_DoubletMask = np.logical_and(NIST_WL > BDoubletMask[0], NIST_WL < BDoubletMask[1]);	
	
UncalFreq, LaserRelMag, SourceRelMag = np.loadtxt(File, unpack = True, skiprows=1);

MaxLaserRelMagVal = np.max(LaserRelMag[100:]);
MaxLaserRelMagInd = np.where(LaserRelMag == MaxLaserRelMagVal);

LaserPeakUncalFreq = UncalFreq[MaxLaserRelMagInd[0][0]];
CalWavelength = 632.82 * LaserPeakUncalFreq / UncalFreq;

VisMask = np.logical_and(CalWavelength > BVisMask[0], CalWavelength < BVisMask[1]);
DoubletMask = np.logical_and(CalWavelength > BDoubletMask[0], CalWavelength < BDoubletMask[1]);

#Now plot the resulting FFT for the SOURCE

NormalisedLaserRelMag  = LaserRelMag  / LaserRelMag[MaxLaserRelMagInd]; 
NormalisedSourceRelMag = SourceRelMag / SourceRelMag[MaxLaserRelMagInd];
SourceMinusLaserRelMag = NormalisedSourceRelMag - NormalisedLaserRelMag;

fig = plt.figure();
ax = fig.add_subplot(111);
#ax.plot(CalWavelength[VisMask], NormalisedLaserRelMag,  label="Source");
#ax.plot(CalWavelength[VisMask], NormalisedSourceRelMag, label="Laser");
ax.plot(CalWavelength[VisMask], SourceMinusLaserRelMag[VisMask] / np.max(SourceMinusLaserRelMag[VisMask]), label = "Source - Laser");
if len(sys.argv) > 2:
	ax.scatter(NIST_WL[NIST_VisMask], NIST_Rel_Int[NIST_VisMask] / np.max(NIST_Rel_Int[NIST_VisMask]), label = "NIST Data");
ax.set_xlabel("Calibrated Wavelength [nm]");
ax.set_ylabel("Relative Intensity");
ax.set_ylim([0, 1.1]);
ax.legend(loc=1);
ax.grid(b = True, which='minor');
ax.grid(b = True, which='major');
ax.minorticks_on();
fig.savefig(File[:-4]+"-plot.png");

'''
fig2 = plt.figure();
ax2 = fig2.add_subplot(111);
ax2.plot(CalWavelength[DoubletMask], SourceMinusLaserRelMag[DoubletMask] / np.max(SourceMinusLaserRelMag[DoubletMask]), label = "Source - Laser");
if len(sys.argv) > 2:
	if(len(NIST_Rel_Int[NIST_DoubletMask]) > 5):
		ax2.plot(NIST_WL[NIST_DoubletMask], NIST_Rel_Int[NIST_DoubletMask] / np.max(NIST_Rel_Int[NIST_DoubletMask]), label = "NIST Data");
ax2.set_xlabel("Calibrated Wavelength [nm]");
ax2.set_ylabel("Relative Intensity");
ax2.set_ylim([0,1.1]);
ax2.legend(loc=1);
ax2.grid(b = True, which='minor');
ax2.grid(b = True, which='major');
ax2.minorticks_on();
fig2.savefig(File[:-4]+"-plot-doublet.png");
'''

def findPeaks(Wavelengths, RelMags):
	#Normalise...
	RelMags = RelMags / np.max(RelMags);
	
	PeakMask = np.logical_and(RelMags > 0.1, Wavelengths > 0);
	return Wavelengths[PeakMask], RelMags[PeakMask];
	
PWavelen, PRelMags = findPeaks(CalWavelength[VisMask], SourceMinusLaserRelMag[VisMask] / np.max(SourceMinusLaserRelMag[VisMask]));

print "Searching for candidate components";
Materials = ("Na","K", "O", "He", "Ne", "Ar", "Kr", "Xe", "Rn", "N", "Br", "C", "Li", "Cl", "B", "H", "Mg", "Ca", "Be", "Rb", "Sr", "Cs", "Ba", "Fr", "Ra", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Al", "Si", "P", "S", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Lu", "Ir", "Pt", "Au", "Hg", "Ti", "Pb", "Bi", "Po", "At");
Materials_Found = [];

for i in np.arange(0, len(PWavelen)):
	WL = PWavelen[i];
	RM = PRelMags[i];
	
	#Find the material to which this peak belongs...	
	for Material in Materials:
		if Material in Materials_Found:
			break;
		WLs, RMs = Nist_Array(Material, BVisMask[0], BVisMask[1]);
		if(len(RMs) == 0):
			break;
		if(np.max(RMs) == 0):
			break;
		RMs = RMs / np.max(RMs); #Normalise
		Mask = np.logical_and(np.logical_and(WLs > WL * 0.99, WLs < WL * 1.01), RMs > 0.1);
		if(len(WLs[Mask]) > 1):
			#Peak!
			#print "Material found: " + str(Material);
			Found = False;
			for MF in Materials_Found:
				if(MF is Material):
					Found = True;
			if not Found:
				Materials_Found.append(Material);
				print "Component Found: " + Material;

print "Found: " + str(Materials_Found);