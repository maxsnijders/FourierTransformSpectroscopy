#Taken from http://www.physics.sfasu.edu/astro/color/spectra.html
from __future__ import division
import numpy as np
import string
import bisect

#The wavelengths in this file should be in in NANOMETERS

temps, rgb = np.loadtxt("blackbody_colors.txt", usecols=(0, 11), dtype='string', unpack=True);

def blackbody_col_temp(temperature):
	#Look it up in temp;
	for i in np.arange(0, len(temps), 1):
		temp = float(temps[i]);
		if(temp > temperature):
			break;
			
	rgbhex =  "#{rgb}".format(rgb = rgb[i]);
	return rgbhex;
	
	
def color_bminusv(bv):
	col = "#000000";
	if(bv >= -0.33):
		col = "#9bb0ff"; 
	if(bv > -0.17):
		col = "#aabfff";
	if(bv >= 0.15):
		col = "#cad7ff";
	if(bv >= 0.44):
		col = "#f8f7ff";
	if(bv >= 0.68):
		col = "#fff4ea";
	if(bv >= 1.15):
		col = "#ffd2a1";
	if(bv >= 1.64):
		col = "#ffcc6f";
	if(bv >= 2.05):
		col = "#000000";
	return col;

def compute_rgb(wl):
	r = 0.0;
	g = 0.0;
	b = 0.0;
	
	if(wl >= 380 and wl <= 440):
		r = -1.0 * (wl - 440.0) / (440.0 - 380.0);
		g = 0;
		b = 1;
	if(wl >= 440 and wl <= 490):
		r = 0;
		g = (wl - 440.0)/(490.0 - 440.0);
		b = 1;
	if(wl >= 490 and wl <= 510):
		r = 0;
		g = 1;
		b = -1 * (wl - 510.0)/(510.0 - 490.0);
	if(wl >= 510 and wl <= 580):
		r = (wl - 510.0)/(580.0 - 510.0);
		g = 1;
		b = 0;
	if(wl >= 580 and wl <= 645):
		r = 1;
		g = -1 * (wl - 645.0)/(645.0 - 580.0);
		b = 0;
	if(wl >= 645 and wl <= 780):
		r = 1;
		g = 0;
		b = 0;
	
	#Fade off near the edges
	a = 1;
	if(wl >= 700):
		a = 0.3 + 0.7 * (780.0 - wl)/(780.0 - 700.0);
	elif(wl <= 420):
		a = 0.3 + 0.7 * (wl - 380.0)/(420.0 - 380.0);
	
	r *= a;
	g *= a;
	b *= a;
	
	return r,g,b;

def rgb_to_hex(r, g, b):
	return '#%02x%02x%02x' % (r, g, b)

def wl_get_hexcolor(wl):
	r,g,b = compute_rgb(wl);
	return rgb_to_hex(r * 255,g * 255,b * 255);