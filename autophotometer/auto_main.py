#! /usr/bin/env python3 

import sys
import os # new
import warnings
try:
    import autophotometer.coordfix as coordfix
    import autophotometer.starseparator as ss
    import autophotometer.db_extract as db_extract
    isModuleInstalled=1
except ModuleNotFoundError:
    isModuleInstalled=0
    import coordfix
    import starseparator as ss
    import db_extract

if not isModuleInstalled:
	sys.stderr.write("\n*******************************************************************************\n\n")
	sys.stderr.write(" The AutoPhotometer package is not installed!     However, it can be used by running \n")
	sys.stderr.write(" 'python autophot.py' from the folder where all the package files are located.\n")

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.visualization.stretch import LinearStretch
import numpy as np
import photutils
from astropy.stats import sigma_clip
import configparser
config = configparser.ConfigParser()

config.read(['conf_autophot.ini', os.path.expanduser('~/.autophotometer/conf_autophot.ini')],)
#config.read('conf_autophot.ini')  #old

sigma_mult = float(config['DEFAULT']['sigma_mult'])
aperture_rad = float(config['DEFAULT']['aperture_rad'])
annulus_min = float(config['DEFAULT']['annulus_min'])
annulus_max = float(config['DEFAULT']['annulus_max'])

def onpick(event, filt1, filt2):
	ind = event.ind[0]
	try:
		print('{:>5} {:11.6f} {:11.6f}  {:8.3f}  {:8.3f}  {:8.3f}   {:8.3f} {:8.3f}     {:8.3f}{:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}   {:8.3f} {:8.3f}    {}'.format(
			ind, 
			obs_array1[ind][4], 
			abs(obs_array1[ind][5]), 
			obs_array1[ind][6], 
			obs_array1[ind][10], 
			abs(obs_array1[ind][8]), 
			obs_array1[ind][7],
			obs_array1[ind][11],
			obs_array2[ind][4], 
			abs(obs_array2[ind][5]), 
			obs_array2[ind][6], 
			obs_array2[ind][10], 
			abs(obs_array2[ind][8]), 
			obs_array2[ind][7],
			obs_array2[ind][11], 
			det_types1[ind])
			)
		log_selected.append('{:>5} {:11.6f} {:11.6f}{:8.3f}{:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}   {:8.3f} {:8.3f}     {:8.3f}{:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}   {:8.3f} {:8.3f}    {}'.format(
			ind, 
			obs_array1[ind][0],
			obs_array1[ind][1],
			obs_array1[ind][4], 
			abs(obs_array1[ind][5]), 
			obs_array1[ind][6], 
			obs_array1[ind][10], 
			abs(obs_array1[ind][8]),
			obs_array1[ind][7],
			obs_array1[ind][11],
			obs_array2[ind][4], 
			abs(obs_array2[ind][5]), 
			obs_array2[ind][6], 
			obs_array2[ind][10], 
			abs(obs_array2[ind][8]), 
			obs_array2[ind][7],
			obs_array2[ind][11], 
			det_types1[ind])
			)
	except:
		print('{:>5} {:11.6f} {:11.6f} {:8.3f} {:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}   {:8.3f} {:8.3f}    {}'.format(
			ind, 
			obs_array1[ind][0],
			obs_array1[ind][1],
			obs_array1[ind][4], 
			abs(obs_array1[ind][5]), 
			obs_array1[ind][6], 
			obs_array1[ind][10], 
			abs(obs_array1[ind][8]), 
			obs_array1[ind][7],
			obs_array1[ind][11], 
			det_types1[ind])
			)
		log_selected.append('{:>5} {:11.6f} {:11.6f}{:8.3f}{:8.3f}  {:8.3f}  {:8.3f}  {:8.3f}   {:8.3f} {:8.3f}    {}'.format(
			ind, obs_array1[ind][0],
			obs_array1[ind][1],
			obs_array1[ind][4], 
			abs(obs_array1[ind][5]),
			obs_array1[ind][6],
			obs_array1[ind][10], 
			abs(obs_array1[ind][8]),
			obs_array1[ind][7],
			obs_array1[ind][11],
			det_types1[ind])
			)
	
def writefiles(obs_array1,obs_array2,filt1,filt2):
	
	sigma_mask1 = sigma_clip(obs_array1[:,3]-obs_array1[:,4], sigma = 3, masked = True)
	sigma_mask2 = sigma_clip(obs_array2[:,3]-obs_array2[:,4], sigma = 3, masked = True)
	semag1 = []
	phmag1 = []
	psmag1 = []
	semagerr1 = []
	for i in range(len(obs_array1)):
		if sigma_mask2[i] != True:
			semag1.append(obs_array1[i][3])
			phmag1.append(obs_array1[i][9])
			psmag1.append(obs_array1[i][4])
			semagerr1.append(obs_array1[i][8])
	semag2 = []
	phmag2 = []
	psmag2 = []
	semagerr2 = []
	for i in range(len(obs_array2)):
		if sigma_mask2[i] != True:
			semag2.append(obs_array2[i][3])
			phmag2.append(obs_array2[i][9])
			psmag2.append(obs_array2[i][4])
			semagerr2.append(obs_array2[i][8])
	print(len(semag1),len(semag2),len(sigma_mask1),len(sigma_mask2))
	uncorrected_array1 = np.stack((semag1,semagerr1,semag2,semagerr2,phmag1,phmag2,psmag1,psmag2), axis=1)
	hdr = 'SExtractor '+filt1+'\tSEx error '+filt1+'\tSExtractor '+filt2+'\tSEx error '+filt2+'\tPhotutils '+filt1+'\tPhotutils '+filt2+'\tPS1 '+filt1+'\t\t\tPS1 '+filt2
	np.savetxt('mags_'+filt1+'-'+filt2+'.txt', uncorrected_array1, header = hdr, fmt = '%.6f', delimiter = '\t\t')

def determine_filter(hdr):
	'''
	header keywords that contain filter names:
		ALFOSC_FASU:
			ALFLTNM
			FAFLTNM
		NOTCAM:
			NCFLTNM1
			NCFLTNM2
	'''
	
	if hdr['INSTRUME'] == 'ALFOSC_FASU':		#checking if the instrument is ALFOSC_FASU, and checking which filter wheel is used.
		if hdr['ALFLTNM'].strip() != 'Open':	#stripping whitespace from header value, as they seemed to be inconsistent between images.
			filt_in = hdr['ALFLTNM']
		elif hdr['FAFLTNM'].strip() != 'Open':
			filt_in = hdr['FAFLTNM']
		else:
			print('Filter in ', hdr['INSTRUME'], ' not found')
	
	elif hdr['INSTRUME'] == 'NOTCAM':
		if hdr['NCFLTNM1'].strip() != 'Open':
			filt_in = hdr['NCFLTNM1']
		elif hdr['NCFLTNM2'].strip() != 'Open':
			filt_in = hdr['NCFLTNM2']
		else:
			print('Filter in ', hdr['INSTRUME'], ' not found')
	
	elif hdr['INSTRUME'] == 'StanCam':
		filt_raw = hdr['STFLTNM'].strip()
		filt_in = filt_raw[0:5]
			# filt_in = hdr['NCFLTNM1']
		# elif hdr['STFLTNM'].strip() != 'Open':
			# filt_in = hdr['NCFLTNM2']
		# else:
			# print('Filter in ', hdr['INSTRUME'], ' not found')
	
	else:
		print(hdr['INSTRUME'], 'instrument not found')
	
	filt_0 = filt_in.split(' ')[0]
	#splitting filter name because the number of spaces in header value inconsistent. Also, many different *'_SDSS filters, and I'm, not sure which ones are needed.
	filters = {					
		'g\'_SDSS'	:	'G_PS',
		'r\'_SDSS'	:	'R_PS',
		'i\'_SDSS'	:	'I_PS',
		'z\'_SDSS'	:	'Z_PS',
		#'y'			:	'Y_PS',
		'B_Bes'		:	'B',
		'V_Bes'		:	'V',
		'R_Bes'		:	'R',
		'i_int'		:	'I',
		'J'			:	'J',
		'H'			:	'H',
		'Ks'		:	'K',
		}
			#should I have the complete names of the filters, instead of just the first part? 
			#If so, I need a way to make the filter names in the header uniform with names in dictionary (from NOT's documentation). So I would have to reduce the number of spaces in the name.
					#That could be done with filt_good = ' '.join(filt_in.split())

	filt_out = filters[filt_0]	
	return(filt_out)	
	
def colorterm_airmass(filt, airmass, ps_bvri, own_mags, dmag_mean3, dmag_StdDev3):  #Vitaly
	'''
	Filters:
		PS1			G R I Z Y (_PS)
		2MASS		J H K_S
		Own images	B V R I   #Vitaly  can also be PS_I PS_Z
	'''

	ext_coeff = {
		'G_PS'	:	0,
		'R_PS'	:	0,
		'I_PS'	:	0,
		'Z_PS'	:	0,
		'Y_PS'	:	0,
		'B'		:	0.22,
		'V'		:	0.12,
		'R'		:	0.08,
		'I'		:	0.04,
		'J'		:	0,
		'H'		:	0,
		'K'		:	0,
		}
	color_term = {
		'G_PS'	:	0,
		'R_PS'	:	0,
		'I_PS'	:	0,
		'Z_PS'	:	0,
		'Y_PS'	:	0,
		'B'		:	0.06,
		'V'		:	-0.06,
		'R'		:	-0.07,
		'I'		:	-0.03,
		'J'		:	0,
		'H'		:	0,
		'K'		:	0,
		}
	# try:
		# filt = determine_filter(hdr)
	# except:
		# filt = 
	if filt in 'BVRI':
		B = ps_bvri[:,0]
		V = ps_bvri[:,1]
		R = ps_bvri[:,2]
		I = ps_bvri[:,3]
		color = {
			'B'		:	B-V,
			'V'		:	B-V,
			'R'		:	V-R,
			'I'		:	V-I,
			}
	if filt in 'JHK':
		J = ps_bvri[:,0]
		H = ps_bvri[:,1]
		K = ps_bvri[:,2]
		color = {
			'J'		:	J-H,
			'H'		:	J-H,
			'K'		:	J-K,
			}
	if filt in ['G_PS', 'R_PS', 'I_PS', 'Z_PS', 'Y_PS']:
		G_PS = ps_bvri[:,0]
		R_PS = ps_bvri[:,1]
		I_PS = ps_bvri[:,2]
		Z_PS = ps_bvri[:,3]
		Y_PS = ps_bvri[:,4]
		color = {
			'G_PS'	:	G_PS-R_PS,
			'R_PS'	:	G_PS-R_PS,
			'I_PS'	:	R_PS-I_PS,
			'Z_PS'	:	I_PS-Z_PS,
			'Y_PS'	:	Z_PS-Y_PS,
			}
	try:
		C = color[filt]
		X = []
		Y = []
		# print(len(C))  #Vitaly
		# print(len(own_mags))  # Vitaly
		for i in range(len(C)):
			if (not np.isnan(C[i])) and (not np.isnan(own_mags[i])):
			# if np.isnan(C[i]) == False and np.isnan(own_mags[i]) == False:
				# print('Color', C[i], own_mags[i] - V[i])      #Vitaly
				# X.append(C[i])
				if filt == 'B':
					dmag = own_mags[i]-B[i]
				elif filt == 'V':
					dmag = own_mags[i]-V[i]
				elif filt == 'R':
					dmag = own_mags[i]-R[i]
				elif filt == 'I':
					dmag = own_mags[i]-I[i]
				elif filt == 'J':
					dmag = own_mags[i]-J[i]
				elif filt == 'H':
					dmag = own_mags[i]-H[i]
				elif filt == 'K':
					dmag = own_mags[i]-K[i]
				elif filt == 'G_PS':
					dmag = own_mags[i]-G_PS[i]
				elif filt == 'R_PS':
					dmag = own_mags[i]-R_PS[i]
				elif filt == 'I_PS':
					dmag = own_mags[i]-I_PS[i]
				elif filt == 'Z_PS':
					dmag = own_mags[i]-Z_PS[i]
				elif filt == 'Y_PS':
					dmag = own_mags[i]-Y_PS[i]
				else:
					print('Unknown filter! Stop!!')
					exit()
				if abs(dmag - dmag_mean3) < dmag_StdDev3 * sigma_mult:  #Vitaly
					X.append(C[i])
					Y.append(dmag)
		# print(len(C),len(X))  # Vitaly
		# print('Color',X,Y)    #Vitaly
		#slope = np.polyfit(X, Y, 1)[0]
		slope1 = np.polyfit(X,Y,1)  #Vitaly
		print('Fit=',slope1)    #Vitaly
		p1 = np.poly1d(slope1)    #Vitaly
		# xp = np.linspace(min(X), max(X), 2)    #Vitaly
		#y2 = np.zeros_like(Y)  # Vitaly
		#print('Y2=',y2)
		y2 = p1(X)  # Vitaly
		dy = Y - y2  # Vitaly
		# print('Y2=',y2)
		#plt.plot(X, Y, '.', X,p1(X), '--')    #Vitaly
		plt.plot(X, Y, '.',color = 'r')    #Vitaly
		plt.plot(X, y2, '--',color = 'b',label="Slope="+'%.3f' % slope1[0] + "   Intersection="+'%.3f' % slope1[1]) #Vitaly
		plt.legend()
		plt.title("Colour-term and Zero-point")
		plt.show()
		plt.close('all')
		#output = color_term[filt] * C + ext_coeff[filt] * airmass
		#output = slope * C + ext_coeff[filt] * airmass
		output = slope1[0] * C + slope1[1]  #Vitaly
	except: #if filter errors or something. Should probably also print some error message in this case.
		print("Something wrong with color terms.")  #Vitaly
		output = 0
	# print('Color', X, Y)  # Vitaly
	return(output)

def photometry_calculations(obs_array, fits_file, xx, yy, fwhm_pix, magerr, ps_bvri, filt, output_data=False):
	# removing escapers and calculating standard deviation 
	magerr_threshold = float(config['DEFAULT']['magerr_threshold'])
	ownmags = obs_array[:,3]
	psmags = obs_array[:,4]
	mag1 = [] # own data
	mag2 = [] # PANSTARRS data
	for i, value in enumerate(ownmags):
		if np.isnan(psmags[i]) == False and obs_array[i][2] == 1 and abs(obs_array[i][5]) < magerr_threshold:
			mag1.append(value)
			mag2.append(psmags[i])
	while len(mag1) < 5:
		if magerr_threshold > 1:
			break
		magerr_threshold += 0.05
		print('Number of low-error detections too small! New magerr_threshold =', format(magerr_threshold, '.3f'))
		mag1 = [] # own data
		mag2 = [] # PANSTARRS data
		for i, value in enumerate(ownmags):
			if np.isnan(psmags[i]) == False and obs_array[i][2] == 1 and abs(obs_array[i][5]) < magerr_threshold:
				mag1.append(value)
				mag2.append(psmags[i])

	mag1 = np.array(mag1)
	mag2 = np.array(mag2)
	dmag = mag1-mag2
	dmag_mean3, dmag_median3, dmag_StdDev3 = sigma_clipped_stats(dmag, sigma=5, maxiters=None)
	
	mag1_corr = mag1 + abs(dmag_mean3)
	
	# getting fits data 
	hdu = fits.open(fits_file)[0]
	header = hdu.header
	image = hdu.data
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore", category=FITSFixedWarning)
		wcs = WCS(header)

	# photutils photometry
	ap_pos = []
	for i in range(len(xx)):
		pos = [xx[i],yy[i]]
		ap_pos.append(pos)

	r_ap = np.mean(fwhm_pix) * aperture_rad
	r_an_min = np.mean(fwhm_pix) * annulus_min
	r_an_max = np.mean(fwhm_pix) * annulus_max
	
	apertures = photutils.CircularAperture(ap_pos, r = r_ap)
	annulus = photutils.CircularAnnulus(ap_pos, r_in = r_an_min, r_out = r_an_max)
	masks = annulus.to_mask(method='center')
	bkg_median = []
	for mask in masks:
		annulus_data = mask.multiply(image)
		annulus_data_1d = annulus_data[mask.data > 0]
		mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(annulus_data_1d)
		bkg_median.append(median_sigclip)
	bkg_median = np.array(bkg_median)
		
	phot = photutils.aperture_photometry(image, apertures)
	phot['annulus_median'] = bkg_median
	phot['aper_bkg'] = bkg_median * apertures.area
	phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
	for col in phot.colnames:
		phot[col].info.format = '%.8g'  # for consistent table output
	ph_mags = []
	ph_dmag = []
	ph_magerrs = []
	for i in range(len(phot)):
		
		flux = phot['aper_sum_bkgsub'][i]
		with warnings.catch_warnings():
			warnings.filterwarnings("ignore", category=RuntimeWarning)
			ph_mag = -2.5*np.log10(flux)
		try:
			d2_mag = ph_mag - psmags[i]
			ph_dmag.append(abs(d2_mag))
			ph_mags.append(ph_mag)
		except ValueError: #if there is no value in database.
			ph_mags.append(np.nan)
			ph_dmag.append(np.nan)
	ph_mags_temp = []
	ph_dmag_temp = []
	for i, value in enumerate(ph_mags):
		if np.isnan(ph_mags[i]) == False and np.isnan(ph_dmag[i]) == False and obs_array[i][2] == 1 and abs(obs_array[i][5]) < magerr_threshold:
			ph_mags_temp.append(value)
			ph_dmag_temp.append(ph_dmag[i])
	
	# Removing escapers and calculating standard deviation for photutils magnitudes.
	ph_dmag_filtered = sigma_clip(ph_dmag_temp, sigma=5, maxiters=None, masked=False, copy=False, cenfunc='median')
	
	ph_dmag_mean3, ph_dmag_median3, ph_dmag_sd3 = sigma_clipped_stats(ph_dmag_temp, sigma=5, maxiters=None)
	airmass = float(header['AIRMASS'])	
	color_corr = colorterm_airmass(filt, airmass, ps_bvri, ownmags, dmag_mean3, dmag_StdDev3)

	ph_truemags = np.array(ph_mags) - color_corr
	truemags = ownmags - color_corr
	
	obs_array = np.insert(obs_array, 6, truemags, axis=1)
	obs_array = np.insert(obs_array, 7, ph_truemags, axis=1)
	obs_array = np.insert(obs_array, 8, magerr, axis=1)

	n_nan = 0
	for x in ph_truemags:
		if np.isnan(x) == True:
			n_nan += 1
	print('photutil magnitude errors: ',n_nan)
	
	p_mag = []
	p_mag_corr = []
	for i, value in enumerate(ownmags):
		if np.isnan(psmags[i]) == False and 1 in obs_array[i] and abs(obs_array[i][5]) < magerr_threshold:
			p_mag.append(ph_mags[i])
			p_mag_corr.append(ph_truemags[i])
	

	obs_array = np.insert(obs_array, 9, ph_mags, axis=1)
	obs_array = np.insert(obs_array, 10, ownmags - dmag_mean3, axis=1)    #Vitaly
	obs_array = np.insert(obs_array, 11, np.array(ph_mags) + ph_dmag_mean3, axis=1)    #Vitaly


	outfile = open('All_Mags.dat', "w")      #Vitaly
	outfile.write('  PS_mag    SEx_mag   SEx_noCT     Ph_mag    Ph_noCT\n')
	for i in range(len(ownmags)):      #Vitaly
		outfile.write('%8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n' % (obs_array[i][4],obs_array[i][6],obs_array[i][10],obs_array[i][7],obs_array[i][11]))
	outfile.close()      #Vitaly

	if output_data == True:
		return(obs_array,image, wcs)
	else:
		return(obs_array)

############################################################################################################################################
def MainProject(fits_file):

	hdu = fits.open(fits_file)[0]
	header = hdu.header
	try:
		filt = determine_filter(header)
		print('Filter 1: ', filt)
	except:
		print()
		print("Available filters")
		print("PANSTARRS, input: G_PS, R_PS, I_PS, Z_PS, Y_PS")
		print("Standard, input:  B, V, R, I")
		print("2MASS, input: J, H, K")
		filterlist = ['G_PS', 'R_PS', 'I_PS', 'Z_PS', 'Y_PS', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
		while True:
			filt = input("Which filter is the image using? ").upper()
			if filt in filterlist:
				break
			else:
				print("Filter not found.")


# coordinate correction 
	coordfix.coordcorr(fits_file, filt)

# star separation 
	global mag, dmag_avg, ph_mags, ph_dmag_avg, fwhm, det_types1    #  Vitaly
	ra,dec,mag,fwhm,xx,yy,meanFWHM,fwhm_pix,det_types1,magerr,ra_all,dec_all,mag_all = ss.findstars(fits_file)     #  Vitaly

	# choose filter
	
			
# PANSTARRS extraction 
	global obs_array1
	print('\nBeginning database extraction...')
	ps_a, ps_d, obs_array1, ps_bvri = db_extract.database_extraction(ra_all, dec_all, mag_all, filt)
	
	obs_array1,image,wcs = photometry_calculations(obs_array1, fits_file, xx, yy, fwhm_pix, magerr, ps_bvri, filt, output_data = True)
	
	stretch = LinearStretch(slope=0.5, intercept=0.5)
	norm = ImageNormalize(image, stretch=LinearStretch(), interval=ZScaleInterval())
	
	fig = plt.figure()
	ax = plt.subplot(projection=wcs)
	ax.imshow(image, origin='lower', norm=norm, cmap='gray', aspect='equal')

	plt.xlabel(r'RA')
	plt.ylabel(r'Dec')
	overlay = ax.get_coords_overlay('icrs')
	overlay.grid(color='white', ls='dotted')
	
	nonps_dec = []
	nonps_ra = []

	for i in range(len(obs_array1)):
		#Vitaly: if not present in the db catalog, and it's not a star, add coords.
		if obs_array1[i][2] == 0 and det_types1[i]=='Likely    ':
			nonps_ra.append(obs_array1[i][0])
			nonps_dec.append(obs_array1[i][1])
			
		else: #if present, add nan to keep the correct length.
		#elif obs_array1[i][2] == 1: #if present, add nan to keep the correct length.
			nonps_ra.append(np.nan)
			nonps_dec.append(np.nan)
	
# drawing fits image and markers
	global ra_sat, dec_sat, ra_nonstar, dec_nonstar, ra_err, dec_err, ra_susp, dec_susp

	plt.scatter(ra_all, dec_all, 10, marker='.', facecolors='w', alpha = 0, transform=ax.get_transform('fk5'), picker=True)
	print('plot:',len(ra_all))
	
	if len(ss.ra_sat)>0:
		plt.scatter(ss.ra_sat, ss.dec_sat, 300, marker='o', facecolors='none', edgecolors='g', transform=ax.get_transform('fk5'), label = "saturated")
	plt.scatter(ss.ra_nonstar, ss.dec_nonstar, marker='X', facecolors='none', edgecolors='c', transform=ax.get_transform('fk5'), label = "nonstar")
	plt.scatter(ss.ra_err, ss.dec_err, 100, marker='h', facecolors='none', edgecolors='c', transform=ax.get_transform('fk5'), label = "large error")
	plt.scatter(ss.ra_susp, ss.dec_susp, 100, marker='p', facecolors='none', edgecolors='m', transform=ax.get_transform('fk5'), label = "large/small fwhm")
	
	plt.scatter(ra, dec, 100, marker='o', facecolors='none', edgecolors='r', transform=ax.get_transform('fk5'), label = "stars")
	print('plot:',np.count_nonzero(~np.isnan(nonps_ra)))
	plt.scatter(nonps_ra,nonps_dec, marker='D', facecolors='none', edgecolors='b', transform=ax.get_transform('fk5'), label = "Transient?")

	#plt.scatter(own_a, own_d, 50, marker='s', facecolors='none', edgecolors='r', transform=ax.get_transform('fk5'), picker=True)
	#plt.scatter(match_a, match_d, 50, marker='D', facecolors='none', edgecolors='b', transform=ax.get_transform('fk5'))
	#plt.scatter(ps_a,ps_d, 100, marker='D', facecolors='none', edgecolors='r', transform=ax.get_transform('fk5'))

	#plt.scatter(obs_array1[:,0],obs_array1[:,1], marker='x', facecolors='b', transform=ax.get_transform('fk5'), picker=True)
	print('plot:',len(obs_array1[:,0]))
	#plt.scatter(obs_array1[:,0],obs_array1[:,1], marker='.', facecolors='w', transform=ax.get_transform('fk5'), picker=True)
	
	filt2 = 'nan'
	global log_selected
	log_selected = []
	
	# calling the function that allows us to click a target and get information.
	fig.canvas.mpl_connect('pick_event', lambda event: onpick(event, filt, filt2))
	plt.legend(loc='best', bbox_to_anchor=(1.0, 1.0))	
	print('Index        RA         DEC     DB_mag    DB_err   SEx_mag   SEx_noCT   SEx_err   Ph_mag   Ph_noCT   Star?')
	plt.show()
	plt.close('all')	
	
	fitsname = fits_file.split('.')[0]
	outfile = open('Selected_Mags_'+fitsname+'.dat', "w")
	outfile.write('Index        RA        DEC     DB_mag   DB_err   SEx_mag   SEx_noCT   SEx_err   Ph_mag   Ph_noCT   Star?\n')
	for i in range(len(log_selected)):
		outfile.write(log_selected[i]+'\n')
	outfile.close()
	
	
############################################################################################################################################
def MainProject_dual(fits_file1, fits_file2):
	# choose filter
	hdu1 = fits.open(fits_file1)[0]
	header1 = hdu1.header
	hdu2 = fits.open(fits_file2)[0]
	header2 = hdu2.header
	try:
		filt1 = determine_filter(header1)
		print('Filter 1: ', filt1)
		filt2 = determine_filter(header2)
		print('Filter 2: ', filt2)
	except:
		print()
		print("Available filters")
		print("PANSTARRS, input: G_PS, R_PS, I_PS, Z_PS, Y_PS")
		print("Standard, input:  B, V, R, I")
		print("2MASS, input: J, H, K")
		filterlist = ['G_PS', 'R_PS', 'I_PS', 'Z_PS', 'Y_PS', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
		while True:
			print("Which filter is the image",fits_file1,"using? ")
			filt1 = input().upper()
			if filt1 in filterlist:
				break
			else:
				print("Filter not found.")

		while True:
			print("Which filter is the image",fits_file2,"using? ")
			filt2 = input().upper()
			if filt2 in filterlist:
				break
			else:
				print("Filter not found.")	
				
# coordinate correction 
	coordfix.coordcorr(fits_file1, filt1)
	coordfix.coordcorr(fits_file2, filt2)

# star separation 
	global mag1, dmag_avg1, ph_mags1, ph_dmag_avg1, fwhm1, det_types1
	global mag2, dmag_avg2, ph_mags2, ph_dmag_avg2, fwhm2, det_types2
	
# PANSTARRS extraction and photometry
	
	
	global obs_array1
	global obs_array2

	print('\nBeginning database extraction...')
	
	ra1,dec1,mag1,fwhm1,xx1,yy1,meanFWHM1,fwhm_pix1,det_types1,magerr1,ra1_all,dec1_all,mag1_all = ss.findstars(fits_file1)
	ra2,dec2,mag2,fwhm2,xx2,yy2,meanFWHM2,fwhm_pix2,det_types2,magerr2,ra2_all,dec2_all,mag2_all = ss.findstars(fits_file2)
	
	# matching the catalog indices
	img1_cat = SkyCoord(ra = ra1_all*u.degree, dec=dec1_all*u.degree)
	img2_cat = SkyCoord(ra = ra2_all*u.degree, dec=dec2_all*u.degree)
	idx, d2d, d3d = match_coordinates_sky(img1_cat, img2_cat)
	max_sep = 1.0*u.arcsec
	sep_constraint = d2d < max_sep
	ra2_all_fix = np.array(ra2_all)[idx[sep_constraint]]
	dec2_all_fix = np.array(dec2_all)[idx[sep_constraint]]
	mag2_all_fix = np.array(mag2_all)[idx[sep_constraint]]
	magerr2_fix = np.array(magerr2)[idx[sep_constraint]]
	xx2_fix = np.array(xx2)[idx[sep_constraint]]
	yy2_fix = np.array(yy2)[idx[sep_constraint]]
	fwhm_pix2_fix = np.array(fwhm_pix2)[idx[sep_constraint]]

	ra1_all = np.array(ra1_all)[sep_constraint]
	dec1_all = np.array(dec1_all)[sep_constraint]
	mag1_all = np.array(mag1_all)[sep_constraint]
	xx1 = np.array(xx1)[sep_constraint]
	yy1 = np.array(yy1)[sep_constraint]
	fwhm_pix1 = np.array(fwhm_pix1)[sep_constraint]
	magerr1 = np.array(magerr1)[sep_constraint]

	ps_a, ps_d, obs_array1, ps_bvri1, = db_extract.database_extraction(ra1_all, dec1_all, mag1_all, filt1)	
	obs_array1,image1,wcs1 = photometry_calculations(obs_array1, fits_file1, xx1, yy1, fwhm_pix1, magerr1, ps_bvri1, output_data = True)
	
	
	ps_a, ps_d, obs_array2, ps_bvri2 = db_extract.database_extraction(ra2_all_fix, dec2_all_fix, mag2_all_fix, filt2)
	obs_array2,image2,wcs2 = photometry_calculations(obs_array2, fits_file2, xx2_fix, yy2_fix, fwhm_pix2_fix, magerr2_fix, ps_bvri2, output_data = True)

	
	stretch = LinearStretch(slope=0.5, intercept=0.5)
	norm1 = ImageNormalize(image1, stretch=LinearStretch(), interval=ZScaleInterval())
	norm2 = ImageNormalize(image2, stretch=LinearStretch(), interval=ZScaleInterval())

	nonps_dec1 = []
	nonps_ra1 = []
	nonps_dec2 = []
	nonps_ra2 = []
	for i in range(len(obs_array1)):
		if obs_array1[i][2] == 0: #if not present in db catalog, add coords.
			nonps_ra1.append(obs_array1[i][0])
			nonps_dec1.append(obs_array1[i][1])
			
		elif obs_array1[i][2] == 1: #if present, add nan.
			nonps_ra1.append(np.nan)
			nonps_dec1.append(np.nan)
		
		if obs_array2[i][2] == 0:
			nonps_ra2.append(obs_array2[i][0])
			nonps_dec2.append(obs_array2[i][1])
		elif obs_array2[i][2] == 1:
			nonps_ra2.append(np.nan)
			nonps_dec2.append(np.nan)

	fig = plt.figure()
	# Plot 1
	ax1 = plt.subplot(1,2,1, projection=wcs1)
	ax1.imshow(image1, origin='lower', norm=norm1, cmap='gray', aspect='equal')
	overlay = ax1.get_coords_overlay('icrs')
	overlay.grid(color='white', ls='dotted')
	ax1.set(xlabel = r'RA', ylabel = r'Dec')

	plt.scatter(ra1, dec1, 100, marker='s', facecolors='none', edgecolors='r', transform=ax1.get_transform('fk5'), label = "stars")
	plt.scatter(obs_array1[:,0],obs_array1[:,1], marker='x', facecolors='b', alpha = 0, transform=ax1.get_transform('fk5'), picker=True)
	plt.scatter(nonps_ra1,nonps_dec1, marker='D', facecolors='none', edgecolors='c', transform=ax1.get_transform('fk5'))

	# Plot 2 
	ax2 = plt.subplot(1,2,2, projection=wcs2)
	ax2.imshow(image2, origin='lower', norm=norm2, cmap='gray', aspect='equal')
	overlay = ax2.get_coords_overlay('icrs')
	overlay.grid(color='white', ls='dotted')
	ax2.set(xlabel = r'RA', ylabel = r'Dec')

	plt.scatter(ra2, dec2, 100, marker='s', facecolors='none', edgecolors='r', transform=ax2.get_transform('fk5'), label = "stars")
	plt.scatter(obs_array2[:,0],obs_array2[:,1], marker='x', facecolors='b', alpha = 0, transform=ax2.get_transform('fk5'), picker=True)
	plt.scatter(nonps_ra2,nonps_dec2, marker='D',  facecolors='none', edgecolors='c', transform=ax2.get_transform('fk5'))


	global log_selected
	log_selected = []
	fig.canvas.mpl_connect('pick_event', lambda event: onpick(event, filt1, filt2))
	print('Filter {}                                                                       Filter {}'.format(filt1, filt2))
	print('Index   DB_mag   DB_err   SEx_mag   SEx_noCT   SEx_err   Ph_mag   Ph_noCT      DB_mag   DB_err   SEx_mag   SEx_noCT   SEx_err   Ph_mag   Ph_noCT   Star?')
	plt.show()
	plt.close('all')	


	outfile = open('Selected_Mags.dat', "w")
	outfile.write('Index   DB_mag_{} DB_err_{} SEx_mag_{} SEx_noCT_{} SEx_err_{} Ph_mag_{} Ph_noCT_{}    DB_mag_{} DB_err_{} SEx_mag_{} SEx_noCT_{} SEx_err_{} Ph_mag_{} Ph_noCT_{} Star?\n'.format(
		filt1,filt1,filt1,filt1,filt1,filt1,filt1,filt2,filt2,filt2,filt2,filt2,filt2,filt2))
	for i in range(len(log_selected)):
		outfile.write(log_selected[i]+'\n')
	outfile.close()
	
######################################################################################################
	
#if __name__ == "__main__":

	#try:
		#fits_file1 = sys.argv[1]
	#except IndexError:
		#print("put the name of the fits file as argument")
		#quit()
	#MainProject(fits_file1)
