import numpy as np
from astropy.io import ascii, fits
import matplotlib
import matplotlib.pyplot as plt
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Ellipse
from astropy.cosmology import wCDM
from ELGlib import read_tiles, fit_redshift_dist, skew_gauss, my_arange


def rand_compl(compl_type, ra, dec, mjd, chunk, outfig=True):

	output_path = 'output/'+chunk+'/'+mjd
	
	if compl_type == 'tiling':
		tmp = 'Tcompleteness'
	elif compl_type == 'redshift':
		tmp = 'Zcompleteness'
	else:
		sys.exit('Error: The completeness type in rand.py must be either "tiling" or "redshift".')
		
	if chunk == 'eboss21':
		outroot = output_path+'/ELG_SGC.eboss21.'+tmp
		completeness_file = output_path+'/ELG_SGC.eboss21.'+tmp+'.asc'
		footprint_ra  = np.array([-43.0, 0.0,0.0,-43.0,-43.0])
		footprint_dec = np.array([-2.0 ,-2.0,2.0, 2.0 , -2.0])
		xlim = [-46,3]
		ylim = [-3,3]
		
	elif chunk == 'eboss22':
		outroot = output_path+'/ELG_SGC.eboss22.'+tmp
		completeness_file = output_path+'/ELG_SGC.eboss22.'+tmp+'.asc'
		footprint_ra  = np.array([ 0.0,45.0,45.0,0.0, 0.0])
		footprint_dec = np.array([-5.0,-5.0, 5.0,5.0,-5.0])
		xlim = [-3,48]
		ylim = [-6,6]
		
	elif chunk == 'eboss23':
		outroot = output_path+'/ELG_NGC.eboss23.'+tmp
		completeness_file = output_path+'/ELG_NGC.eboss23.'+tmp+'.asc'
		footprint_ra  = np.array([126,126,142.5,142.5,157,157,136.5,136.5,126])
		footprint_dec = np.array([16,29,29,27,27,13.8,13.8,16,16])
		xlim = [125,158]
		ylim = [13,30]
		
	# reading the completeness
	completeness  = np.zeros(0)
	tilecomb_uniq = np.zeros(0,dtype='str') 
	with open(completeness_file) as f:
		for line in f:
			if (line[0]!='#'):
				completeness  = np.append(completeness, float(line.split()[0]))
				tilecomb_uniq = np.append(tilecomb_uniq,line.split()[1])
				
	# reading the tiles
	tile   = np.zeros(0,dtype='str')
	tilera = np.zeros(0)
	tiledec= np.zeros(0)
	
	tmp    		= read_tiles('data/tiles-'+chunk+'.par')
	tile   		= np.append(tile, tmp[0])
	tilera 		= np.append(tilera, tmp[1])
	tiledec		= np.append(tiledec,tmp[2])
	ntile 		= len(tilera)
	
	# number of spectroscopic objects
	nspecobj    = len(ra)
	objSkyCoord = SkyCoord(ra=ra*u.deg,dec=dec*u.deg,frame='fk5')
	tile_rad    = 1.49 # tile radius in degrees
	
	# listing the tiles in which each object lies
	spectile = np.array(['' for x in np.arange(nspecobj)],dtype=object)
	for i in xrange(ntile):
		tileSkyCoord = SkyCoord(ra=tilera[i]*u.deg,dec=tiledec[i]*u.deg,frame='fk5')
		tmpspec  = (objSkyCoord.separation(tileSkyCoord).deg <= tile_rad)
		spectile[tmpspec] += tile[i]+'-'
	
	# loop on unique combinations of tiles
	speccompl = np.zeros(nspecobj)
	for i in xrange(len(tilecomb_uniq)):
		tmpspec            = (spectile==tilecomb_uniq[i])
		speccompl[tmpspec] = completeness[i]
		
	if outfig:
		# building ellipses
		ells = [Ellipse(xy=np.array([x,y]),width=3.0/np.cos(y/180.*np.pi),height=3.0,angle=0)
		        for x,y in zip(tilera,tiledec)]
		
		if chunk in ['eboss21','eboss22']:
			fig  = plt.figure(figsize=(20,8.8))
		elif chunk == 'eboss23':
			fig  = plt.figure(figsize=(16,8))
		
		ax = fig.add_subplot(111)
		
		# eBOSS/ELG footprint
		ax.plot(footprint_ra,footprint_dec,c='k',zorder=10,linewidth=2)
		if compl_type == 'redshift':
			clim = [50.,100.*max(speccompl)]
		else:
			clim = [100.*min(speccompl),100.*max(speccompl)]
		SC   = ax.scatter(ra,dec,s=2,c=100.*speccompl,edgecolor='none',\
		        vmin=clim[0],vmax=clim[1],cmap=matplotlib.cm.rainbow)
		for e in ells:
		    e.set_color('k')
		    e.set_facecolor('none')
		    e.set_linewidth(1)
		    e.set_alpha(0.3)
		    e.set_clip_box(ax.bbox)
		    ax.add_artist(e)
		# colorbar
		cbar = plt.colorbar(SC, pad = 0.01)
		if compl_type == 'tiling':
			cbar.set_label(r'Tiling completeness [\%]',size=16)
		elif compl_type == 'redshift':
			cbar.set_label(r'Redshift completeness [\%]',size=16)
		cbar.set_clim(clim)
		# plot stuff
		ax.grid(True,linestyle='-',alpha=0.3)
		ax.set_xlabel(r'R.A. [deg.]',size=16)
		ax.set_ylabel(r'DEC. [deg.]',size=16)
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		plt.tight_layout()
		plt.savefig(outroot+'.png',bbox_inches='tight',format='png')
	   
	return speccompl
	
	


		
###########################################################################


def generate_elg_random(mjd, chunk, area, output_path, random_path, out_random, out_catalog, out_dndz, 
						z_min, z_max, binwidth_z, weights, tiling_compl_cut, redshift_compl_cut, 
						down_sample_fact, H0, omega_m, omega_lam, w_DE):
	
	
	run_info = open(output_path+'/run.info','a')


	# Read plate centers
	tilelist  = ascii.read('data/ebosstiles.dat')
	tiles     = tilelist["TILE"]
	racen     = tilelist["racen"]
	deccen    = tilelist["deccen"]
	obs_tiles = np.loadtxt(output_path+'/tiles.txt')
	
	rac  = []
	decc = []

	i = 0
	for tile in tiles:
		if tile in obs_tiles:
			rac.append(racen[i])
			decc.append(deccen[i])
		i+=1
	
	rac     = np.asarray(rac)
	rac     = np.where(rac>270.,rac-360.,rac)
	decc    = np.asarray(decc)
	nplates = len(rac)
	#----------------------


	# open random fits catalog
	hdu     = fits.open(random_path)
	data    = hdu[1].data
	tmp     = data['ra']
	randra  = np.where(tmp>270.,tmp-360.,tmp)
	randdec = data['dec']
	#----------------------
	
	
	# Modify the random number density
	idx = np.random.choice(len(randra),size=int(len(randra)/down_sample_fact))
	
	new_density = (10000./down_sample_fact)
	info_line   = ("Number density in the random catalog: %.2f deg^-2\n\n" % new_density)
	print info_line
	run_info.write(info_line)

	randra  = randra[idx]
	randdec = randdec[idx]
	nrand   = len(randra)
	
	randSkyCoord = SkyCoord(ra=randra*u.deg,dec=randdec*u.deg,frame='fk5')
	#--------------------


	# open data ascii file created by function 'generate_elg_catalog'
	tmp         = ascii.read(out_catalog+'.dat')
	datara      = tmp['ra']
	datadec     = tmp['dec']
	dataz       = tmp['z']
	dataw_compl = tmp['w_compl']
	dataw_sys   = tmp['w_sys']
	ndata 		= len(datara)
	#--------------------
	
	
	# selecting randoms within the plate
	inplate  = np.zeros(nrand, dtype='bool')
	tile_rad = 1.49

	for k in range(nplates):
		cenSkyCoord = SkyCoord(ra=rac[k]*u.deg,dec=decc[k]*u.deg,frame='fk5')
		tmp =  (randSkyCoord.separation(cenSkyCoord).deg <= tile_rad)
		inplate[tmp] = True

	randra  = np.asarray(randra[inplate])
	randdec = np.asarray(randdec[inplate])
	#-------------------------------------
	
	
	# Remove randoms with a completeness below a certain treshold
	randc = rand_compl('tiling', randra, randdec, mjd, chunk)
	
	if tiling_compl_cut > 0.:
		
		tmp = (randc >= tiling_compl_cut)
		randra  = randra[tmp]
		randdec = randdec[tmp]
	
	randc = rand_compl('redshift', randra, randdec, mjd, chunk)
		
	if redshift_compl_cut > 0.:
		
		tmp = (randc >= redshift_compl_cut)
		randra  = randra[tmp]
		randdec = randdec[tmp]

	nrand  = len(randra)
	
	
	# assign random redshift to random objects, 
	# following the same redshift distribution as in data
	dataw_tot = dataw_compl*dataw_sys
	
	randz = np.random.choice(dataz, size = nrand, p=dataw_tot/sum(dataw_tot))
	#------------------------------------
	
	
	mbins  = my_arange(z_min, z_max , binwidth_z)
	
	
	###---------------###
	### TOTAL WEIGHTS ###
	###---------------###
	# To be assigned to data and random
	
	randw 	  = np.ones(nrand)
	dataw_FKP = np.ones(ndata)
	
	
	if 3 in weights: # The FKP weights
	
		def weight_FKP(x):
			return 1./(1. + x*10000.)
		
		# open the data redshift distribution
		tmp = ascii.read(out_dndz+'.dat')
	
		bincenters = tmp['z']
		nz         = tmp['nz']
		
		# Calculation of the comoving volume
		cosmo  = wCDM(H0=H0, Om0=omega_m, Ode0=omega_lam, w0=w_DE) 
			

		# sky surface = 4pi sr = 129600/pi deg^2
		tmp  = cosmo.comoving_volume(mbins).value * (np.pi / 129600.) * cosmo.h**3. # h^-3 * Mpc^3 / deg^2
		dV   = tmp[1:] - tmp[:-1]
		
		nz   = nz / area / dV
		#-----------------------------------
		
		
		# Fitting of the redshift distribution with a skew gaussian
		# mu, sigmag, alpha, a, c1, c2, c3
		skew_init = [np.mean(dataz), np.std(dataz), 0., max(nz), 0., 0., 0.] # initial guess
		popt      = fit_redshift_dist(bincenters, nz, skew_init) # optimal parameters of the fit
		
		info_line  = ("Initial guess for the fit of n(z):\n")
		info_line += ("mu    = %.2f\n" % skew_init[0])
		info_line += ("sigma = %.2f\n" % skew_init[1])
		info_line += ("alpha = %.2f\n" % skew_init[2])
		info_line += ("amp   = %.2e\n" % skew_init[3])
		info_line += ("c1    = %.2e\n" % skew_init[4])
		info_line += ("c2    = %.2e\n" % skew_init[5])
		info_line += ("c3    = %.2e\n\n" % skew_init[6])	
		
		print info_line
		run_info.write(info_line)
		
		info_line  = ("Output parameters for the fit of n(z):\n")
		info_line += ("mu    = %.2f\n" % popt[0])
		info_line += ("sigma = %.2f\n" % popt[1])
		info_line += ("alpha = %.2f\n" % popt[2])
		info_line += ("amp   = %.2e\n" % popt[3])
		info_line += ("c1    = %.2e\n" % popt[4])
		info_line += ("c2    = %.2e\n" % popt[5])
		info_line += ("c3    = %.2e\n\n" % popt[6])	
		
		print info_line
		run_info.write(info_line)
		#-------------------------------------------------------
		
		
		# Calculate the FKP for data and randoms
		nz_fit    = skew_gauss(randz,*popt)
		randw  	  = weight_FKP(nz_fit)
		
		nz_fit    = skew_gauss(dataz,*popt)
		dataw_FKP = weight_FKP(nz_fit)
		#---------------------------------------
		
		# print n(z)
		plt.figure()
		tmp = np.linspace(z_min,z_max,1000)
		plt.plot(bincenters,nz*1.e4,'bo', label='data')
		plt.plot(tmp,skew_gauss(tmp,*skew_init)*1.e4,'g-', label='inital guess')
		plt.plot(tmp,skew_gauss(tmp,*popt)*1.e4,'r-', label='final fit')
		plt.legend(loc='upper right')	
		plt.xlabel(r'$z$',size=18)
		plt.ylabel(r'$n(z)$ [$10^{-4}h^3$Mpc$^{-3}$]',size=18)
		plt.savefig(output_path+'/nz_fit.pdf',format='pdf')
	
	
	
	
	
	# Calculate the total weights to assign to data
	dataw_tot = np.ones(ndata)
	
	if 1 in weights:
		dataw_tot *= dataw_compl
		
		info_line = ('Using completeness weights\n')
		print info_line
		run_info.write(info_line)
		
	if 2 in weights:
		dataw_tot *= dataw_sys
		
		info_line = ('Using systematic weights\n')
		print info_line
		run_info.write(info_line)
	
	if 3 in weights:
		dataw_tot *= dataw_FKP	
		
		info_line = ('Using FKP weights\n')
		print info_line
		run_info.write(info_line)
	
	if len(weights) == 0:
		
		info_line = ('No weights are used\n')
		print info_line
		run_info.write(info_line)
	
	
	# The effective redshift
	tmp       = np.average(dataz, weights=dataw_tot)
	info_line = ("\n")
	info_line += ("Effective redshift: %.2f\n" % tmp)
	
	# The effective number of data
	tmp        = sum(dataw_tot)
	info_line += ("Effective number of data: %.6f\n" % tmp)
	
	# The effective number of randoms
	tmp        = sum(randw)
	info_line += ("Effective number of randoms: %.6f\n\n" % tmp)
	
	print info_line
	run_info.write(info_line)
	#----------------------------------------------
	
	
		
	
	# writing the catalog for data	
	tmp = np.transpose(np.vstack((datara, datadec, dataz, dataw_tot)))
	ascii.write(tmp,out_catalog+'.dat', format='no_header', 
				names = ['ra','dec','z','w'], 
				formats={'ra':'%.6f', 'dec':'%.6f', 'z':'%.6f', 'w':'%.6f'})
				
	col1 = fits.Column(name='ra',  format='D', array=datara)
	col2 = fits.Column(name='dec', format='D', array=datadec)
	col3 = fits.Column(name='z',   format='D', array=dataz)
	
	col4 = fits.Column(name='w_FKP',   format='D', array=dataw_FKP)
	col5 = fits.Column(name='w_compl', format='D', array=dataw_compl)
	col6 = fits.Column(name='w_sys',   format='D', array=dataw_sys)
	cols = fits.ColDefs([col1, col2, col3, col4, col5, col6])

	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(out_catalog+'.fits', clobber=True)
	
	
	# writing the random
	tmp = np.transpose(np.array([randra,randdec,randz,randw]))
	ascii.write(tmp,out_random+'.dat', format='no_header', 
				names = ['randra','randdec','randz','randw'], 
				formats={'randra':'%.6f', 'randdec':'%.6f', 'randz':'%.6f', 'randw':'%.6f'})

	col1 = fits.Column(name='randra',    format='D', array=randra)
	col2 = fits.Column(name='randdec',   format='D', array=randdec)
	col3 = fits.Column(name='randz',     format='D', array=randz)
	col4 = fits.Column(name='randw', format='D', array=randw)
	cols = fits.ColDefs([col1, col2, col3, col4])

	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(out_random+'.fits', clobber=True)
	
	#---- Figures ----	

	if chunk == 'eboss21':
		plt.figure(figsize=(20,8.8))
		footprint_ra  = np.array([-43.0, 0.0,0.0,-43.0,-43.0])
		footprint_dec = np.array([-2.0 ,-2.0,2.0, 2.0 , -2.0])
		xlim = [-48.5,5]
		ylim = [-2.5,2.5]
	
	elif chunk == 'eboss22':
		
		plt.figure(figsize=(20,8.8))
		footprint_ra  = np.array([ 0.0,45.0,45.0,0.0, 0.0])
		footprint_dec = np.array([-5.0,-5.0, 5.0,5.0,-5.0])
		xlim = [-3,48]
		ylim = [-6,6]

	elif chunk == 'eboss23':
		
		plt.figure(figsize=(16,8))
		footprint_ra  = np.array([126,126,142.5,142.5,157,157,136.5,136.5,126])
		footprint_dec = np.array([16,29,29,27,27,13.8,13.8,16,16])
		xlim = [125,158]
		ylim = [13,30]
	
	# ra, dec
	plt.scatter(datara,datadec,s=1,color='b',alpha=0.4,zorder=2,label='Data : '+str(ndata)+' objects')
	plt.scatter(randra,randdec,s=0.2,color='r',zorder=1,label='Randoms : '+str(nrand)+' objects')
	plt.plot(footprint_ra,footprint_dec,c='k',zorder=10,linewidth=2)
	plt.xlabel(r'R.A. [deg]', size=16)
	plt.ylabel(r'DEC. [deg]', size=16)
	plt.xlim(xlim)
	plt.ylim(ylim)
	plt.tight_layout()
	plt.grid(True,linestyle='-',alpha=0.3)
	plt.savefig(output_path+'/footprint.png',format='png')

	# z histogram of randoms

	plt.figure()
	plt.hist(dataz,bins=mbins,color='b',edgecolor='none',normed=True,histtype='bar',label='Data : '+str(ndata)+' galaxies')
	plt.hist(randz,bins=mbins,color='r',linewidth = 2,normed=True,histtype='step',label='Randoms : '+str(nrand)+' objects')
	plt.xlabel(r'$z$', size=18)
	plt.ylabel(r'$n(z)$', size=18)
	plt.xlim(z_min,z_max)
	plt.legend()
	plt.grid(True,linestyle='-',alpha=0.3)
	plt.savefig(output_path+'/hist_rand.pdf',format='pdf')

	run_info.close()
