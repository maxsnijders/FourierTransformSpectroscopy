import os  
import re

for fn in os.listdir('In/'):
	#Verwerken!
	Path = "In/" + fn;
	Out = "Out/" + fn[:-4] + ".png"
	
	os.system("python cli.py " + re.escape(Path) + " " + re.escape(Out));