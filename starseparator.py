import subprocess
import numpy as np
import configparser
config = configparser.ConfigParser()
config.read('conf_autophot.ini')


def read_file(file):
# SExtractor catalog columns
#   1 NUMBER                 Running object number                                     
#   2 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
#   3 FLUX_ISO               Isophotal flux                                             [count]
#   4 FLUX_MAX               Peak flux above background                                 [count]
#   5 XPEAK_IMAGE            x-coordinate of the brightest pixel                        [pixel]
#   6 YPEAK_IMAGE            y-coordinate of the brightest pixel                        [pixel]
#   7 ALPHAPEAK_J2000        Right ascension of brightest pix (J2000)                   [deg]
#   8 DELTAPEAK_J2000        Declination of brightest pix (J2000)                       [deg]
#   9 CLASS_STAR             S/G classifier output                                     
#  10 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
#  11 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
#  12 FLAGS                  Extraction flags              
	data =[]
	with open(file) as f:
		for line in f.readlines():
			if line.startswith("#") == True:
				continue
			else:
				#print(line)
				line = line.rstrip('\n')
				line_list = [s for s in line.split(' ') if s]
				line = line.strip().split(" ")		
				filtered = list(filter(None, line))
				data.append(filtered)	
	f.close()
	return(data)

#Vitaly
def get_powers(n):
    """Get positive powers of two which add up to n.

    Parameters
    ----------
    n : positive integer

    Returns
    -------
    list of integers which are powers of 2

    Examples
    --------
    >>> get_powers(14)
    [2, 4, 8]

    >>> get_powers(5)
    [1, 4]
    """
    #get_bin = lambda x: x >= 0 and str(bin(x))[2:] or "-" + str(bin(x))[3:]
    n = int(n)
    get_bin = lambda x: x >= 0 and str(bin(x))[2:] or "-" + str(bin(x))[3:]
    bin_str = get_bin(n)  # convert n to a binary string
    powers = []
    for i, digit in enumerate(bin_str[::-1]):
        if digit == '1':
            powers.append(2**i)
    return powers
#Vitaly

# extract specific elements from data
def e_fwhm(lst): 
    return [float(item[1]) for item in lst]
def e_fiso(lst): 
    return [float(item[2]) for item in lst]
def e_fmax(lst): 
    return [float(item[3]) for item in lst]
def e_x(lst): 
    return [float(item[4]) for item in lst]
def e_y(lst): 
    return [float(item[5]) for item in lst]
def e_a(lst): 
    return [float(item[6]) for item in lst]
def e_d(lst): 
    return [float(item[7]) for item in lst]
def e_sg(lst): 
    return [float(item[8]) for item in lst]
def e_mag(lst):
	return [float(item[9]) for item in lst]
def e_magerr(lst):
	return [float(item[10]) for item in lst]
def e_flags(lst):
	return [float(item[11]) for item in lst]
def e_fwhm_pix(lst):
	return [float(item[12]) for item in lst]


def runsex(cmd):
	subprocess.call(cmd)

def findstars(fits_file):
	filename = fits_file.split('.')
	name = filename[0]	
	cat2 = "cat_{}.cat".format(name)
	# command = ["sex", fits_file, "-c", "ph2conf.sex", "-CATALOG_NAME", cat2]
	command = ["sex", fits_file, "-c", "ph2conf.sex", "-CATALOG_NAME", cat2, "-PARAMETERS_NAME", "run2.param"]
	runsex(command)

	data = read_file(cat2)

# Vitaly 20210330 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
	from astropy.stats import sigma_clipped_stats
	from astropy.stats import sigma_clip
	from scipy import stats
	import plotille

	sg_thresh = float(config['DEFAULT']['sg_threshold'])
	MagErrLimit = float(config['DEFAULT']['magerr_threshold'])

	fwhm1 = e_fwhm(data)
	sg1   = e_sg(data)
	magerr1 = e_magerr(data)
	flags1 = e_flags(data)

	FWHMtemp = []
	for i in range(len(fwhm1)):
		if ((magerr1[i] <= 2*MagErrLimit) and (sg1[i] > sg_thresh)) and not(4 in get_powers(flags1[i])):
			FWHMtemp.append(fwhm1[i]*3600)
	if FWHMtemp == []:
		print('Bad config parameters MagErrLimit and sg_thresh!')
		quit()
	avgfwhm = np.median(FWHMtemp)

	command2 = ["sex", fits_file, "-c", "ph2conf.sex", "-CATALOG_NAME", cat2, "-PARAMETERS_NAME", "run2.param", "-SEEING_FWHM", str(avgfwhm)]
	runsex(command2)

	data2 = read_file(cat2)
	sg2 = e_sg(data2)

	fiso = e_fiso(data2)
	fmax = e_fmax(data2)
	x	 = e_x(data2)
	y	 = e_y(data2)
	a	 = e_a(data2)
	d	 = e_d(data2)
	fwhm = e_fwhm(data2)
	sg 	 = e_sg(data2)
	mag	 = e_mag(data2)
	magerr = e_magerr(data2)
	flags = e_flags(data2)
	fwhm_pix = e_fwhm_pix(data2)

#	print("\n\n*******\nFWHM distribution\n")

	HistBins=int((max(FWHMtemp)-min(FWHMtemp))/0.1)
#	print(plotille.hist(FWHMtemp))

#	print("\n\n*******\nClipped FWHM distribution\n")

	filtered_data = sigma_clip(FWHMtemp, sigma=3, maxiters=None, masked=False, copy=False, cenfunc='median')
	HistBins=int((max(filtered_data)-min(filtered_data))/0.05)
#	print(plotille.hist(filtered_data,bins=HistBins))

# scipy.stats.tmean  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.tmean.html#scipy.stats.tmean
	Mean3, Median3, StdDev3 = sigma_clipped_stats(FWHMtemp, sigma=3, maxiters=None)
	print("\n3 sigma Clipped Mean=", Mean3, "=/-", StdDev3)
	print("Number of star-candidates:        ", len(FWHMtemp))
	print("Number of clipped star-candidates:", len(filtered_data))
	print("\n")


	# print("\n\n*******\nSorting the detections\n")
	FWHMtemp = []
	for i in range(len(fwhm)):
		if ((magerr[i] <= MagErrLimit) and (sg[i] > sg_thresh)) and not(4 in get_powers(flags[i])):
			FWHMtemp.append(fwhm[i]*3600)
	filtered_data = sigma_clip(FWHMtemp, sigma=3, maxiters=None, masked=False, copy=False, cenfunc='median')

	global xx_all, yy_all, mag_all, fwhm_all
#	global ra_all, dec_all, ra_sat, dec_sat, ra_nonstar, dec_nonstar, ra_err, dec_err, ra_susp, dec_susp
	global ra_sat, dec_sat, ra_nonstar, dec_nonstar, ra_err, dec_err, ra_susp, dec_susp

	n_sat   = 0      # Number of saturated objects
	n_nonstars = 0   # Number of non-stars
	n_err   = 0      # Large photometric errors
	n_susp  = 0      # Number of suspicious objects with Large FWHMs
	n_stars = 0      # Number of stars

	xx_all   = x
	yy_all   = y
	ra_all   = a
	dec_all  = d
	mag_all  = mag
	fwhm_all = 3600 * fwhm

	ra_sat = []
	dec_sat = []
	ra_nonstar = []
	dec_nonstar = []
	ra_err = []
	dec_err = []
	ra_susp = []
	dec_susp = []

	xx_sat   = []
	yy_sat   = []
	ra_sat   = []
	dec_sat  = []
	mag_sat  = []
	fwhm_sat = []

	xx_nonstar   = []
	yy_nonstar   = []
	ra_nonstar   = []
	dec_nonstar  = []
	mag_nonstar  = []
	fwhm_nonstar = []

	xx_err   = []
	yy_err   = []
	ra_err   = []
	dec_err  = []
	mag_err  = []
	fwhm_err = []

	xx_susp   = []
	yy_susp   = []
	ra_susp   = []
	dec_susp  = []
	mag_susp  = []
	fwhm_susp = []

	xx=[]
	yy=[]
	ra=[]
	dec=[]
	mags=[]
	fwhms = []
	
	det_type = []
	
	for i in range(len(fiso)):
		if (sg[i] < sg_thresh):       # Probably not a star
			n_nonstars += 1
			xx_nonstar.append(x[i])
			yy_nonstar.append(y[i])
			ra_nonstar.append(a[i])
			dec_nonstar.append(d[i])
			mag_nonstar.append(mag[i])
			fwhm_nonstar.append(3600 * fwhm[i])
			det_type.append('Nonstar   ')

		elif 4 in get_powers(flags[i]):   # Saturated?
			# nerr += 1
			n_sat += 1
			xx_sat.append(x[i])
			yy_sat.append(y[i])
			ra_sat.append(a[i])
			dec_sat.append(d[i])
			mag_sat.append(mag[i])
			fwhm_sat.append(3600*fwhm[i])
			det_type.append('Saturated ')

		elif (magerr[i]>=MagErrLimit):  # Large photometric error
			n_err += 1
			xx_err.append(x[i])
			yy_err.append(y[i])
			ra_err.append(a[i])
			dec_err.append(d[i])
			mag_err.append(mag[i])
			fwhm_err.append(3600*fwhm[i])
			det_type.append('Large_err ')
			
		elif (3600*fwhm[i]>max(filtered_data) or 3600*fwhm[i]<min(filtered_data)):        # Suspicious?
			n_susp += 1
			xx_susp.append(x[i])
			yy_susp.append(y[i])
			ra_susp.append(a[i])
			dec_susp.append(d[i])
			mag_susp.append(mag[i])
			fwhm_susp.append(3600*fwhm[i])
			det_type.append('Bad FWHM  ')
			
		else:															#else, it's probably a star.
			n_stars += 1
			xx.append(x[i])
			yy.append(y[i])
			ra.append(a[i])
			dec.append(d[i])
			mags.append(mag[i])
			fwhms.append(3600*fwhm[i])
			det_type.append('Likely    ')


	print("\n\n*******\nSorting the detections:\n")

	print("Number of all extracted objects:", len(fwhm))
	print("Number of saturated detections:", n_sat)
	print("Number of non-stars:", n_nonstars)
	print("Number of stars with large photometric errors: ", n_err)
	print("Number of suspicious objects (large or small FWHM): ", n_susp)
	print("Number of stars: ", n_stars)

	if n_stars < 5:
		print('WARNING!!! A very small number of stars. Consider changing the input parameters.')

	return ra,dec,mags,fwhms,x,y,Mean3,fwhm_pix,det_type,magerr,ra_all,dec_all,mag_all
