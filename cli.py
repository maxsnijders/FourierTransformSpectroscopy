from __future__ import division
import sys
from Spectrum import *
import numpy as np

if(len(sys.argv) < 3):
	print "You need to pass the file to analyse and the plot out name.";
	sys.exit(0);
	
Data_File = sys.argv[1];
Plot_Out = sys.argv[2];

#Load it!
UncalFreq, LaserRelMag, SourceRelMag = np.loadtxt(Data_File, unpack = True, skiprows=1);

#Now process it
S = Spectrum(UncalFreq, LaserRelMag, SourceRelMag,lowLimit=300,highLimit=1500,withoutLaser=True,laserBuffer=1,
	Component_Frac_Min=100,NM_Peak_Allowance=0.2,Peak_Min_RelMag=50);
	
S.plot_range(Plot_Out, printPeaks = False, printPlanck = False);
Components = S.analyse_spectrum();

for Component in Components.items():
	Name = Component[0];
	Count = Component[1][0];
	Percent = Component[1][1];
	
	print "Component " + Name + " matched " + str(Count) + " peaks (" + str(round(Percent * 100, 1)) + "%)."