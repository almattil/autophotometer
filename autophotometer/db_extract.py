#! /usr/bin/env python3
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import json
import requests
from scipy import spatial
from urllib.parse import quote as urlencode
from urllib.request import urlretrieve
import http.client as httplib 
import matplotlib.pyplot as plt
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
import configparser
import os


config = configparser.ConfigParser()
config.read(['conf_autophot.ini', os.path.expanduser('~/.autophotometer/conf_autophot.ini')],)


def ps1cone(ra,dec,radius,table="mean",release="dr1",format="csv",columns=None,
		   baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
		   **kw):
	"""Do a cone search of the PS1 catalog
	
	Parameters
	----------
	ra (float): (degrees) J2000 Right Ascension
	dec (float): (degrees) J2000 Declination
	radius (float): (degrees) Search radius (<= 0.5 degrees)
	table (string): mean, stack, or detection
	release (string): dr1 or dr2
	format: csv, votable, json
	columns: list of column names to include (None means use defaults)
	baseurl: base URL for the request
	verbose: print info about request
	**kw: other parameters (e.g., 'nDetections.min':2)
	"""
	
	data = kw.copy()
	data['ra'] = ra
	data['dec'] = dec
	data['radius'] = radius
	return ps1search(table=table,release=release,format=format,columns=columns,
					baseurl=baseurl, verbose=verbose, **data)


def ps1search(table="mean",release="dr1",format="csv",columns=None,
		   baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
		   **kw):
	"""Do a general search of the PS1 catalog (possibly without ra/dec/radius)
	
	Parameters
	----------
	table (string): mean, stack, or detection
	release (string): dr1 or dr2
	format: csv, votable, json
	columns: list of column names to include (None means use defaults)
	baseurl: base URL for the request
	verbose: print info about request
	**kw: other parameters (e.g., 'nDetections.min':2).  Note this is required!
	"""
	
	data = kw.copy()
	url = "{baseurl}/{release}/{table}.{format}".format(**locals())
	if columns:
		# check that column values are legal
		# create a dictionary to speed this up
		dcols = {}
		for col in ps1metadata(table,release)['name']:
			dcols[col.lower()] = 1
		badcols = []
		for col in columns:
			if col.lower().strip() not in dcols:
				badcols.append(col)
		if badcols:
			raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
		data['columns'] = '[{}]'.format(','.join(columns))

	r = requests.get(url, params=data)
	if verbose:
		print(r.url)
	r.raise_for_status()
	if format == "json":
		return r.json()
	else:
		return r.text


def ps1metadata(table="mean",release="dr1",
		   baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
	"""Return metadata for the specified catalog and table
	
	Parameters
	----------
	table (string): mean, stack, or detection
	release (string): dr1 or dr2
	baseurl: base URL for the request
	
	Returns an astropy table with columns name, type, description
	"""
	
	url = "{baseurl}/{release}/{table}/metadata".format(**locals())
	r = requests.get(url)
	r.raise_for_status()
	v = r.json()
	# convert to astropy table
	tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
			   names=('name','type','description'))
	return tab


def twomass_query(ra, dec, radius, verbose = True, **kw):
	url = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?&catalog=fp_psc"
	
	data = kw.copy()
	data['spatial'] = 'cone'
	data['objstr'] = str(ra)+' '+str(dec)
	data['radius'] = radius
	data['radunits'] = 'deg'
	data['outfmt'] = '1'
	
	r = requests.get(url, params=data)
	if verbose:
		print(r.url)
	r.raise_for_status()
	
	tab = ascii.read(r.text)
	return(tab)


def read_file(file):
	data =[]
	with open(file) as f:
		for line in f.readlines():
			if line.startswith("#") == True:
				continue
			else:
				line = line.rstrip('\n')
				line_list = [s for s in line.split(' ') if s]
				line = line.strip().split(" ")		
				filtered = list(filter(None, line))
				data.append(filtered)
	f.close()
	return(data)

#extract from PANSTARRS data
def e_a(lst): 
	return [float(item[1]) for item in lst]
def e_d(lst): 
	return [float(item[2]) for item in lst]
def e_g(lst): 
	return [float(item[9]) for item in lst]
def e_r(lst): 
	return [float(item[10]) for item in lst]
def e_i(lst): 
	return [float(item[11]) for item in lst]
def e_z(lst): 
	return [float(item[12]) for item in lst]
def e_y(lst): 
	return [float(item[13]) for item in lst]
def e_g_err(lst): 
	return [float(item[14]) for item in lst]
def e_r_err(lst): 
	return [float(item[15]) for item in lst]
def e_i_err(lst): 
	return [float(item[16]) for item in lst]
def e_z_err(lst): 
	return [float(item[17]) for item in lst]
def e_y_err(lst): 
	return [float(item[18]) for item in lst]
	
#extract from own data
def e2_a(lst): 
	return [float(item[1]) for item in lst]
def e2_d(lst): 
	return [float(item[2]) for item in lst]
def e2_mag(lst):
	return [float(item[3]) for item in lst]
	
	
def gri2bvri(g,r,i):	
	B1 = g + 0.561*(g - r) + 0.194
	B2 = g + 0.540*(g - r)  + 0.016*(g - r)**2 + 0.199

	V1 = g - 0.508*(g - r) - 0.017
	V2 = r + 0.492*(g - r) - 0.017

	R1 = r - 0.166*(g - r) - 0.142
	R2 = r - 0.275*(r - i) - 0.166

	I1 = i - 0.167*(g - r) - 0.376
	I2 = i - 0.214*(r - i) - 0.416
	return(B2, V2, R2, I2)

def catalog_matching(cat1, cat2):	
		idx, d2d, d3d = match_coordinates_sky(cat1, cat2) #idx - indices of cat2 that match cat1, d2d - on-sky separation between closest match
		max_sep = 1.0*u.arcsec
		sep_constraint = d2d < max_sep
		sep_const_miss = d2d > max_sep
		common_detections = cat1[sep_constraint]
		catalog_matches = cat2[idx[sep_constraint]]
		return(common_detections, catalog_matches, idx)

def findlargest(ar1, ar2):
	ar_min = np.full(len(ar1), 0.05)
	arcomb = np.stack((ar1,ar2,ar_min),axis = 1)
	ar3 = np.amax(arcomb, axis=1)
	return(ar3)


def database_extraction(own_a,own_d,own_mag,filt):
	ra = np.median(own_a)
	dec = np.median(own_d)
	cone_radius = str(config['DEFAULT']['cone_radius'])
	radius = cone_radius
	constraints = {'nDetections.gt':1}
	# strip blanks and weed out blank and commented-out values
	columns = """objID,raMean,decMean,nDetections,ng,nr,ni,nz,ny,
		gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag,gMeanPSFMagErr,rMeanPSFMagErr,iMeanPSFMagErr,zMeanPSFMagErr,yMeanPSFMagErr""".split(',')
	columns = [x.strip() for x in columns]
	columns = [x for x in columns if x and not x.startswith('#')]
	
	filterlist1 = ['J', 'H', 'K']
	filterlist2 = ['B', 'V', 'R', 'I']
	filterlist3 = ['G_PS', 'R_PS', 'I_PS', 'Z_PS', 'Y_PS']
	if filt in filterlist1:
		tab = twomass_query(ra,dec,radius,verbose = False, **constraints)
		ps_a = np.array(tab['ra'])
		ps_d = np.array(tab['dec'])
	elif filt in filterlist2 or filt in filterlist3:
		results = ps1cone(ra,dec,radius,release='dr2',columns=columns,verbose=False,**constraints)
		tab = ascii.read(results)
		for filter in 'grizy':
			col = filter+'MeanPSFMag'
			try:
				tab[col].format = ".4f"
				tab[col][tab[col] == -999.0] = np.nan
			except KeyError:
				print("{} not found".format(col))
		ps_a = e_a(tab)
		ps_d = e_d(tab)
		B, V, R, I = gri2bvri(np.array(e_g(tab)), np.array(e_r(tab)), np.array(e_i(tab)))
		B_err = findlargest(np.array(e_g_err(tab)), np.array(e_r_err(tab)))
		V_err = findlargest(np.array(e_g_err(tab)), np.array(e_r_err(tab)))
		R_err = findlargest(np.array(e_r_err(tab)), np.array(e_i_err(tab)))
		I_err = findlargest(np.array(e_r_err(tab)), np.array(e_i_err(tab)))
	print('Data extracted!')

	if filt == "G_PS":
		mags = e_g(tab)
		magerr = e_g_err(tab)
	elif filt == "R_PS":
		mags = e_r(tab)
		magerr = e_r_err(tab)
	elif filt == "I_PS":
		mags = e_i(tab)
		magerr = e_i_err(tab)
	elif filt == "Z_PS":
		mags = e_z(tab)
		magerr = e_z_err(tab)
	elif filt == "Y_PS":
		mags = e_y(tab)
		magerr = e_y_err(tab)
	elif filt == "B":
		mags = B
		magerr = B_err
	elif filt == "V":
		mags = V
		magerr = V_err
	elif filt == "R":
		mags = R
		magerr = R_err
	elif filt == "I":
		mags = I
		magerr = I_err
	elif filt == 'J':
		mags = tab['j_m']
		magerr = tab['j_msigcom']
	elif filt == 'H':
		mags = tab['h_m']
		magerr = tab['h_msigcom']
	elif filt == 'K':
		mags = tab['k_m']
		magerr = tab['k_msigcom']	#I could do other filters like this too, instead of making functions...
	
	own_cat = SkyCoord(ra = own_a*u.degree, dec=own_d*u.degree)
	ps_cat = SkyCoord(ra = ps_a*u.degree, dec=ps_d*u.degree)
	
	c_matches, cat2_matches, idx = catalog_matching(own_cat, ps_cat)
	obs_array = np.stack((own_cat.ra.value, own_cat.dec.value), axis=1)
	presence = [] 
	matches_set = np.stack((c_matches.ra.value, c_matches.dec.value), axis=1) # changing to array from skycoords, presence determination loop works a lot faster this way.
	for obj in obs_array:
		if obj in matches_set:
			presence.append(1)
		else:
			presence.append(0)
			
	obs_array = np.insert(obs_array, 2, presence, axis=1)

	if filt in filterlist1:
		J = np.array(tab['j_m'])
		H = np.array(tab['h_m'])
		K = np.array(tab['k_m'])
	f1_sort = []
	f2_sort = []
	f3_sort = []
	f4_sort = []
	f5_sort = []
	obs_mags = []
	ps_mags = []
	ps_magerr = []
	for i in range(len(idx)):
		obs_mags.append(own_mag[i])
		if presence[i] == 1:
			ps_mags.append(mags[idx[i]])
			ps_magerr.append(magerr[idx[i]])
		elif presence[i] == 0:
			ps_mags.append(np.nan)
			ps_magerr.append(np.nan)
			
		if filt in filterlist2:
			f1_sort.append(B[idx[i]])
			f2_sort.append(V[idx[i]])
			f3_sort.append(R[idx[i]])
			f4_sort.append(I[idx[i]])
			db_mags_full = np.stack((f1_sort, f2_sort, f3_sort, f4_sort), axis=1)
		elif filt in filterlist1:
			f1_sort.append(J[idx[i]])
			f2_sort.append(H[idx[i]])
			f3_sort.append(K[idx[i]])
			db_mags_full = np.stack((f1_sort, f2_sort, f3_sort), axis=1)
		elif filt in filterlist3:
			f1_sort.append(e_g(tab)[idx[i]])
			f2_sort.append(e_r(tab)[idx[i]])
			f3_sort.append(e_i(tab)[idx[i]])
			f4_sort.append(e_z(tab)[idx[i]])
			f5_sort.append(e_y(tab)[idx[i]])
			db_mags_full = np.stack((f1_sort, f2_sort, f3_sort, f4_sort, f5_sort), axis=1)
	obs_array = np.insert(obs_array, 3, obs_mags, axis=1)
	obs_array = np.insert(obs_array, 4, ps_mags, axis=1)
	obs_array = np.insert(obs_array, 5, ps_magerr, axis=1)

	return(ps_a, ps_d, obs_array, db_mags_full)
