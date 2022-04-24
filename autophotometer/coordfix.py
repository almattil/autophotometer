#! /usr/bin/env python3
from astropy.io import fits
import subprocess
import ccdproc
import os # Vitaly 20211108


def runscamp2mass():
#Vitaly 20211108
	ph1cat_file='phase1.cat'
	scampconf_file='scamp.conf'
	if not os.path.exists(ph1cat_file):
		ph1cat_file=os.path.expanduser('~/.autophot/phase1.cat')
	if not os.path.exists(scampconf_file):
		scampconf_file=os.path.expanduser('~/.autophot/scamp.conf')
#Vitaly 20211108
	subprocess.call(['scamp', "-c", scampconf_file, ph1cat_file, '-ASTREF_CATALOG', '2MASS'])


def runscamp():
#Vitaly 20211108
	ph1cat_file='phase1.cat'
	scampconf_file='scamp.conf'
	if not os.path.exists(ph1cat_file):
		ph1cat_file=os.path.expanduser('~/.autophot/phase1.cat')
	if not os.path.exists(scampconf_file):
		scampconf_file=os.path.expanduser('~/.autophot/scamp.conf')
#Vitaly 20211108
	subprocess.call(['scamp', "-c", scampconf_file, ph1cat_file])

def runsex(fits_file):
# Vitaly 20211108
	run1_file="run1.param"
	sexdef_file='default.sex'
	if not os.path.exists(sexdef_file):
		sexdef_file=os.path.expanduser('~/.autophot/default.sex')
	if not os.path.exists(run1_file):
		run1_file=os.path.expanduser('~/.autophot/run1.param')		
	subprocess.call(['sex', "-c", sexdef_file, fits_file, "-PARAMETERS_NAME", run1_file])
# Vitaly 20211108

	
def read_header(filename): #function that reads .head files created by SCAMP
	headerlines = []
	with open(filename) as f:
		for line in f.readlines()[3:-1]:
			line = line.replace('/', '=').strip()	#removing all extra bits from the lines
			line = line.split('=')
			value = line[1].replace("'",'').strip()
			try:	#convert the values into floats if possible
				line[1] = float(value)
			except ValueError:
				line[1] = value
			headerlines.append(line)	#put all elements in a single list
	return headerlines

def fits_backup(fits_file):
	fits_inf = fits.open(fits_file, 'update', save_backup=True)
	#fits_inf.info()
	fits_inf.close()


def addheads(fits_head, fits_file):
	heads = read_header(fits_head)
	for head in heads:
		try:
			fits.delval(fits_file, head[0])
		except KeyError:
			continue
		fits.setval(fits_file, head[0], value=head[1])			

def cosmics(fits_file):
	hdu = fits.open(fits_file, mode = 'update')
	image = hdu[0].data
	header = hdu[0].header
	try:
		gain = header['GAIN']
	except:
		try:
			gain = header['GAIN1'] #some images seem to have multiple gain keywords.
		except:
			gain = 1 #if no gain found, set it to 1 and hope for the best.
	hdu[0].data = ccdproc.cosmicray_lacosmic(image, gain = gain)[0]
	hdu.close()
	
def coordcorr(fits_file, filt):
	filename = fits_file.split('.')
	fits_head = 'phase1.head'
	
	print('Do you want to make changes to', fits_file,'before continuing?\n[y]es, [n]o, [q]uit.')
	while True:
		query = str(input('y/n/q: ')).lower()
		if query == 'y':
			print('Do you want to make a backup of the file?')
			bak_query = str(input('y/n: ')).lower()
			
			print('Do you want to remove cosmic rays?')
			cosmic_query = str(input('y/n: ')).lower()
			
			print('Do you want to do the coordinate correction?')
			corr_query = str(input('y/n: ')).lower()
			
			if bak_query == 'y':
				fits_backup(fits_file)
			
			if cosmic_query == 'y':
				cosmics(fits_file)
			
			if corr_query == 'y':
				runsex(fits_file)
				if filt in ['J', 'H', 'K']:
					runscamp2mass()
				else:
					runscamp()
				addheads(fits_head, fits_file)
				print('----Done!----')
			break
		elif query == 'n':
			break
		elif query == 'q':
			exit()
		
		else:
			print('Unknown input.')
