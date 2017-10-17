import os
import astropy.io.fits as fits
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from pydl.pydlutils.yanny import yanny
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib
import subprocess
import re
import datetime
import ConfigParser 
from scipy.optimize import curve_fit
import scipy.special as sp
from scipy.stats import chisquare

config = ConfigParser.ConfigParser()
config.read("ELG_config.ini")
stilts_path	= config.get('files', 'stilts_path')

STILTSCMD='java -jar '+stilts_path+'/stilts.jar '

pid      = str(os.getpid())


def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('-', s) ]


def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    
def my_arange(a,b,dr,decimals=6):

    res = [a]

    k = 1
    while res[-1] < b:
    	
    	tmp = round(a + k*dr,decimals)
    	
    	if tmp > b:
        	break	
        
        res.append(tmp)
        k+=1

    return np.asarray(res)
   
   
def skew_gauss(x, mu, sigmag, alpha, a, c1, c2, c3):
    
    x_mu = x-mu
    
    normpdf = np.exp(-0.5*x_mu*x_mu/sigmag/sigmag)
    normcdf = (1. + sp.erf(alpha*(x_mu/sigmag)/np.sqrt(2.)))
    
    return a*normpdf*normcdf + c1*x*x + c2*x + c3
    
    
def fit_redshift_dist(z,nz,init):
	""" Fit the redshift distribution with a skew gaussian dist. 
	"""
	p   = init
	res = []
	
	chi2_old  = 1.e6
	chi2_diff = 1.e6
	k = 0
	
	while ( (chi2_diff > 1.e-4) & (k < 2) ):
	
		try:
			popt, pcov = curve_fit(skew_gauss, z, nz, p)
			#bounds=([min(z),0.0,-np.inf,0.0,-np.inf,-np.inf,-np.inf],[max(z),2.*np.std(z),np.inf,np.inf,np.inf,np.inf,np.inf]))
		except:
			pass
			sys.exit('Bad initial condition for n(z) fit.')

		chi2      = chisquare(nz, skew_gauss(z, *popt))[0]
		chi2_diff = abs(chi2 - chi2_old)
		chi2_old  = chi2
		
		res.append((chi2, popt))
		p = popt
		
		k += 1 # do at least 2 iterations
	
	chi2, popt = min(res, key=lambda x:x[0])
	
	# return the paramter of the fit
	return popt
    
    
# function to read the tile properties
def read_tiles(parfile):
    a  = yanny(parfile)
    tile= np.array(a['STRUCT1']['tile']).astype('str')
    ra = np.array(a['STRUCT1']['racen'])
    tmp= (ra>270.)
    ra[tmp] = ra[tmp] - 360.
    dec= np.array(a['STRUCT1']['deccen'])
    return [tile,ra,dec]


# function to attribute the tiling completeness to a set of (ra,dec)
def Tiling_completeness(specra, specdec, mjd, chunk, write_file=True):
	# inputs:
	# - specra :  spectroscopic ra
	# - specdec:  spectroscopic dec
	# - chunk  :  chunk eboss 21,22 or 23
	# output:
	# - speccompl: np.array with the tiling completeness (0<completeness<1)
    
    # settings

	output_path = 'output/'+chunk+'/'+mjd

	if chunk == 'eboss21':

		photcat = 'data/elg_240_sgc.masked.ForTiling.fits'
		outroot = output_path+'/ELG_SGC.eboss21.Tcompleteness'
	
	elif chunk == 'eboss22':
		
		photcat = 'data/elg_240_sgc.masked.ForTiling.fits'
		outroot = output_path+'/ELG_SGC.eboss22.Tcompleteness'

	elif chunk == 'eboss23':
		
		photcat = 'data/elg_190_ngc.masked.ForTiling.fits'
		outroot = output_path+'/ELG_NGC.eboss23.Tcompleteness'
    
    # reading the tiles
	tile    = np.zeros(0,dtype='str')
	tilera  = np.zeros(0)
	tiledec = np.zeros(0)
	
	tmp		= read_tiles('data/tiles-'+chunk+'.par')
	tile	= np.append(tile, tmp[0])
	tilera	= np.append(tilera, tmp[1])
	tiledec = np.append(tiledec,tmp[2])
	ntile 	= len(tilera)
	
	# reading the photcat
	hdu         = fits.open(photcat)
	photra      = hdu[1].data['ra']
	tmp         = (photra>270.)
	photra[tmp] = photra[tmp]-360.
	photdec     = hdu[1].data['dec']
	
	if chunk == 'eboss21':
		tmp		= (photra < 0.)
		photra	= photra[tmp]
		photdec = photdec[tmp]
	
	elif chunk == 'eboss22':
		tmp		= (photra > 0.)
		photra	= photra[tmp]
		photdec = photdec[tmp]

	photSkyCoord = SkyCoord(ra=photra*u.deg,dec=photdec*u.deg,frame='fk5')
	specSkyCoord = SkyCoord(ra=specra*u.deg,dec=specdec*u.deg,frame='fk5')
	nphotobj     = len(photra)
	nspecobj     = len(specra)
	
	tile_rad = 1.49 # tile radius in degrees
	
	# listing the tiles in which each object lies
	phottile = np.array(['' for x in np.arange(nphotobj)],dtype=object)
	spectile = np.array(['' for x in np.arange(nspecobj)],dtype=object)
	for i in xrange(ntile):
		tileSkyCoord = SkyCoord(ra=tilera[i]*u.deg,dec=tiledec[i]*u.deg,frame='fk5')
		tmpphot = (photSkyCoord.separation(tileSkyCoord).deg <= tile_rad)
		tmpspec = (specSkyCoord.separation(tileSkyCoord).deg <= tile_rad)
		phottile[tmpphot] += tile[i]+'-'
		spectile[tmpspec] += tile[i]+'-'
	
	# unique combinations of tiles
	tilecomb_uniq = np.unique(phottile[(phottile!='')]) # includes all np.unique(spectile), by construction
	ntilecomb_uniq= len(tilecomb_uniq)
	
	# loop on unique combinations of tiles
	completeness = np.zeros(ntilecomb_uniq)
	speccompl    = np.zeros(nspecobj)
	
	for i in xrange(ntilecomb_uniq):
	    tmpphot = (phottile==tilecomb_uniq[i])
	    tmpspec = (spectile==tilecomb_uniq[i])
	    percent = float(len(spectile[tmpspec])) / float(len(phottile[tmpphot]))
	    if (percent>0):
			completeness[i] = percent
			speccompl[tmpspec] = percent
	
	# writing to file
	if write_file:
		f = open(outroot+'.asc', 'w')
		f.write('# Created at GMT '+datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")+'\n')
		f.write('# Chunks: '+chunk+'\n')
		f.write('# COMPLETENESS[%]\tTILE_COMBINATION\n')
		for i in xrange(ntilecomb_uniq):
		    f.write(str(completeness[i])+'\t'+tilecomb_uniq[i]+'\n')
		f.close()

	return speccompl


# function to attribute the redshift completeness to a set of (ra,dec)
def Redshift_completeness(goodra, gooddec, obsra, obsdec, mjd, chunk, write_file=True):
	# inputs:
	# - goodra :  ra position from good ELG spectra with reliable redshift
	# - obsra  :  ra position of ELG spectroscopic targets
	# - chunk  :  chunk eboss 21,22 or 23
	# output:
	# - Zcompl : np.array with the redshift completeness (0<completeness<1)
    
    # settings

	output_path = 'output/'+chunk+'/'+mjd

	if chunk == 'eboss21':

		outroot = output_path+'/ELG_SGC.eboss21.Zcompleteness'
	
	elif chunk == 'eboss22':
		
		outroot = output_path+'/ELG_SGC.eboss22.Zcompleteness'

	elif chunk == 'eboss23':
		
		outroot = output_path+'/ELG_NGC.eboss23.Zcompleteness'
    
    # reading the tiles
	tile    = np.zeros(0,dtype='str')
	tilera  = np.zeros(0)
	tiledec = np.zeros(0)
	
	tmp		= read_tiles('data/tiles-'+chunk+'.par')
	tile	= np.append(tile, tmp[0])
	tilera	= np.append(tilera, tmp[1])
	tiledec = np.append(tiledec,tmp[2])
	ntile 	= len(tilera)

	obsSkyCoord  = SkyCoord(ra=obsra*u.deg,dec=obsdec*u.deg,frame='fk5')
	goodSkyCoord = SkyCoord(ra=goodra*u.deg,dec=gooddec*u.deg,frame='fk5')
	nobsobj      = len(obsra)
	ngoodobj     = len(goodra)
	
	tile_rad = 1.49 # tile radius in degrees
	
	# listing the tiles in which each object lies
	obstile = np.array(['' for x in np.arange(nobsobj)],dtype=object)
	goodtile = np.array(['' for x in np.arange(ngoodobj)],dtype=object)
	for i in xrange(ntile):
		tileSkyCoord = SkyCoord(ra=tilera[i]*u.deg,dec=tiledec[i]*u.deg,frame='fk5')
		tmpobs  = (obsSkyCoord.separation(tileSkyCoord).deg <= tile_rad)
		tmpgood = (goodSkyCoord.separation(tileSkyCoord).deg <= tile_rad)
		obstile[tmpobs]   += tile[i]+'-'
		goodtile[tmpgood] += tile[i]+'-'
	
	# unique combinations of tiles
	tilecomb_uniq = np.unique(obstile[(obstile!='')]) 
	ntilecomb_uniq= len(tilecomb_uniq)
	
	# loop on unique combinations of tiles
	completeness = np.zeros(ntilecomb_uniq)
	Zcompl       = np.zeros(ngoodobj)
	
	for i in xrange(ntilecomb_uniq):
	    tmpobs = (obstile==tilecomb_uniq[i])
	    tmpgood = (goodtile==tilecomb_uniq[i])
	    percent = float(len(goodtile[tmpgood])) / float(len(obstile[tmpobs]))
	    if (percent>0):
			completeness[i] = percent
			Zcompl[tmpgood] = percent
	
	# writing to file
	if write_file:
		f = open(outroot+'.asc', 'w')
		f.write('# Created at GMT '+datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")+'\n')
		f.write('# Chunks: '+chunk+'\n')
		f.write('# COMPLETENESS[%]\tTILE_COMBINATION\n')
		for i in xrange(ntilecomb_uniq):
		    f.write(str(completeness[i])+'\t'+tilecomb_uniq[i]+'\n')
		f.close()

	return Zcompl


# function to remove duplicates, if any
def rmv_dups(specfits,run_info,output_path):
# inputs:
# - specfits: fits catalog of data galaxies
# output:
# - keepind: np.array of True/False
# identifying duplicates
	tmpstr = (STILTSCMD+' tmatch1 action=identify '+
				'matcher=sky values="PLUG_RA PLUG_DEC" params=0.1 '+
				'in='+specfits+' ifmt=fits '+
				'out=tmp.fits_'+pid+' ofmt=fits')
	subprocess.call(tmpstr+' &> '+output_path+'/stilts.log', shell=True)
	# keeping non duplicates and, for duplicates, the one with the best redshift
	hdu     = fits.open('tmp.fits_'+pid)
	groupid = hdu[1].data['GroupID']
	groupsize=hdu[1].data['GroupSize']
	zwarning= hdu[1].data['ZWARNING']
	rchi2   = hdu[1].data['RCHI2']
	
	## non-duplicates
	keepind = (groupsize<0) # empty -> -2147483648
	## for duplicates, the one with the best redshift
	groupid_dups = np.unique(groupid[(groupsize>=2)])
	info_line = ("%0d duplicates found\n" % len(groupid_dups))
	print info_line
	run_info.write(info_line)
	for gid in groupid_dups:
		dupsind = np.arange(len(rchi2))[(groupid==gid)]                 # indexes of the dups from this group
	
	dupsind   = dupsind[np.argsort(zwarning[dupsind])]              # sorting by increasing zwarning
	dupsind   = dupsind[(zwarning[dupsind]==zwarning[dupsind][0])]  # selecting dups with smallest zwarning
	dupsind   = dupsind[np.argsort(np.abs(rchi2[dupsind]-1.))]      # sorting by abs(rchi2-1)
	keepind[dupsind[0]] = True                                      # adding the dups with lowest abs(rchi2-1)
	
	## cleaning
	subprocess.call('rm tmp.fits_'+pid, shell=True)

	return keepind
