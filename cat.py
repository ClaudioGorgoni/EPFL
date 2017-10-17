import numpy as np 
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import sys
from ELGlib import rmv_dups, Tiling_completeness, Redshift_completeness, my_arange



def generate_elg_catalog(mjd, chunk, output_path, catalog_path, out_catalog, out_dndz, z_min, z_max, 
						 binwidth_z, weights, zwarning, plates_to_rmv, template_noqso, add_qso, comparat, 
						 tiling_compl_cut, redshift_compl_cut):
	
	#output_path = 'output/'+chunk+'/'+mjd
	
	run_info = open(output_path+'/run.info','a')
	
	if template_noqso:
		add_qso = False
	
	info_line  = ("Using chunk : "+chunk+"\n")
	info_line += ("Tiling completeness cut   : %.2f%%\n" % (100.*tiling_compl_cut) )
	info_line += ("Redshift completeness cut : %.2f%%\n" % (100.*redshift_compl_cut) )
	print info_line
	run_info.write(info_line)
	
	# Remove the duplicates
	dupl_mask = rmv_dups(catalog_path, run_info, output_path)	
	
	#----------------------------------------
	# OPEN THE DATA FILE AND READ THE COLUMNS
	#----------------------------------------
	
	hdu		  = fits.open(catalog_path)
	data      = hdu[1].data
	
	eboss_target = data["EBOSS_TARGET1"]
	schunk		 = data["chunk"].strip()
	plates       = data["PLATE"]
	
	# Remove plates from list 'plates_to_rmv'
	plates_to_keep = np.array([pp not in plates_to_rmv for pp in plates], dtype=bool)

	if chunk in ['eboss21','eboss22']:
		elg_num   = 2**43
		surf_dens = 240.
	
	elif chunk == 'eboss23':
		elg_num   = 2**44
		surf_dens = 204.
	
	# Selection criterion for ELG target
	tmp = ( (eboss_target == elg_num) & plates_to_keep & dupl_mask )
	num_candidate = float(sum(tmp))
	
	ra     = data["PLUG_RA"][tmp]
	dec    = data["PLUG_DEC"][tmp]
	w_sys  = data["syst_weight"][tmp]
	
	# Save the used tiles in a file
	tiles = np.sort(list(set(list(data["TILE"][tmp]))))
	np.savetxt(output_path+'/tiles.txt', tiles, fmt='%0d')
	
	info_line  = ( "Number of plates: %0d\n" % len(tiles) )
	print info_line
	run_info.write(info_line)
	
	if template_noqso:
		sclass 	  = data["CLASS_NOQSO"].strip()[tmp]
		z         = data["Z_NOQSO"][tmp]
		zW        = data["ZWARNING_NOQSO"][tmp]
		#rchi2     = data["RCHI2_NOQSO"][tmp]
		#z_err     = data["Z_ERR_NOQSO"][tmp]
		#rchi2diff = data["RCHI2DIFF_NOQSO"][tmp]
	else:
		sclass    = data["CLASS"].strip()[tmp]
		z         = data["Z"][tmp]
		zW        = data["ZWARNING"][tmp]
		#rchi2     = data["RCHI2"][tmp]
		#z_err     = data["Z_ERR"][tmp]
		#rchi2diff = data["RCHI2DIFF"][tmp]
	
	if comparat:
	
		if template_noqso:
			zq   				   = data["Z_NOQSO_zQ"][tmp]
			zcont 				   = data["Z_NOQSO_zCont"][tmp]
			zq[np.isnan(zq)]       = -99.0
			zcont[np.isnan(zcont)] = -99.0
		else:
			zq   				   = data["Z_zQ"][tmp]
			zcont 				   = data["Z_zCont"][tmp]
			zq[np.isnan(zq)]       = -99.0
			zcont[np.isnan(zcont)] = -99.0

	old_z = z
	
	ra    = np.where(ra>270.0, ra-360.0,ra)
	
	Tcompl = Tiling_completeness(ra, dec, mjd, chunk)
	
	
	# Remove object with a tiling completeness below a certain treshold
	if tiling_compl_cut > 0.0:
		
		tmp = (Tcompl >= tiling_compl_cut)
		
		sclass    = sclass[tmp]
		ra        = ra[tmp]
		dec       = dec[tmp]
		z         = z[tmp]
		zq        = zq[tmp]
		zcont     = zcont[tmp]
		zW        = zW[tmp]
		w_sys     = w_sys[tmp]
		Tcompl    = Tcompl[tmp]
		#z_err     = z_err[tmp]
		#rchi2diff = rchi2diff[tmp]
		#rchi2     = rchi2[tmp]
	#---------------------------------------	
	
	area  = sum(1./Tcompl)/surf_dens # normalized to target density
	
	
	
	
	
	#-------------------
	# FILTER THE OBJECTS
	#-------------------
	
	gal     = np.array([ra,dec,z,Tcompl,w_sys])
	z_range = ( (z_min <= z) & (z <= z_max) )
	
	
	# STEP 1
	if template_noqso:
		info_line = "Using templates with NOQSO\n"
		
	else:
		info_line = "Using templates with QSO\n"
		
	in_gal_catalog = ( (sclass == "GALAXY") & ([zz in zwarning for zz in zW]) )
	#------------------------------------------------------------------------------



	# STEP 2
	if comparat:
		info_line += "Using flags zq and zcont of J. Comparat\n"
		
		comparat_criterion = ( (zq >= 2.) | ((zq >= 1.) & (zcont > 0.)) | ((zq >= 0.) & (zcont >= 2.5)) )
		
		gal = [ g[(in_gal_catalog & z_range & comparat_criterion)] for g in gal ]
		
	else:
		info_line += "NOT using flags zq and zcont of J. Comparat\n"
		
		gal = [ g[(in_gal_catalog & z_range)] for g in gal ]
	#------------------------------------------------------------------------------
	
	
	
	# STEP 3
	if add_qso:
		info_line += "QSOs added to the catalog\n"
		
		qso = [ra,dec,z,Tcompl,w_sys]
		in_qso_catalog = ( (sclass == "QSO")    & ([zz in zwarning for zz in zW]) )
		
		if comparat:
			qso = [ q[(in_qso_catalog & z_range & comparat_criterion)] for q in qso ]
			
		else:
			qso = [ q[(in_qso_catalog & z_range)] for q in qso ]
			
		for i in range(len(gal)):
			gal[i] = np.append(gal[i],qso[i])
	#------------------------------------------------------------------------------
			
			
			
	print info_line
	run_info.write(info_line)
	

	tmp    = (sclass != 'STAR')
	Zcompl = Redshift_completeness(gal[0], gal[1], ra[tmp], dec[tmp], mjd, chunk)
	
	# Remove object with a redshift completeness below a certain treshold
	if redshift_compl_cut > 0.0:
		
		tmp = (Zcompl >= redshift_compl_cut)
		
		ra     = gal[0][tmp]
		dec    = gal[1][tmp]
		z      = gal[2][tmp]
		Tcompl = gal[3][tmp]
		w_sys  = gal[4][tmp]
	#--------------------------------------------------------------------
	
	
	num_obj = len(z)

	w_compl = 1./Tcompl 

	# writing the catalog
	tmp = np.transpose(np.vstack((ra, dec, z, w_compl, w_sys)))
	ascii.write(tmp, out_catalog+'.dat', names = ['ra','dec','z','w_compl','w_sys'], 
				formats={'ra':'%.6f', 'dec':'%.6f', 'z':'%.6f', 'w_compl':'%.6f', 'w_sys':'%.6f'})
	#---------------------------------------------------------------------------------------------
	
	
	
	# Redshift histogram
	bins       = my_arange(z_min,z_max,binwidth_z)
	bincenters = bins[0:-1]+binwidth_z/2.
	
	w_tot = np.ones(num_obj)
	
	if 1 in weights:
		w_tot *= w_compl
		
	if 2 in weights:
		w_tot *= w_sys
		
	tmp = np.histogram(z, bins = bins, weights = w_tot, density = False )[0]
	
	z_dist = np.transpose(np.vstack((bincenters[bincenters <= z_max], tmp[bincenters <= z_max])))
	ascii.write(z_dist,out_dndz+'.dat', names=['z','nz'], formats={'z':'%.6f', 'nz':'%.6f'})
	#---------------------------------------------------------------------------------------
		
		
	info_line = "\n"
	info_line += ("Number of objects in the catalog in range %.2f < z < %.2f: %.0f/%.0f (%.2f%%)\n"
		% (z_min, z_max, num_obj, num_candidate, 100.*num_obj/num_candidate) )
	print info_line
	run_info.write(info_line)
	
	
	# figures
	f, ax = plt.subplots()
	
	ax.hist(old_z,bins=my_arange(0.,2.,binwidth_z),color='grey',histtype='stepfilled', alpha=0.8, label=r'ELG targets : %0d' % num_candidate)
	ax.hist(z,bins=bins,color='r',linewidth = 2,histtype='step', label=(r'Kept objects : %0d (%.2f' % (num_obj,100.*num_obj/num_candidate))+(r'\%)') )
	
	line1 = (r'Using chunk : '+chunk)
	line2 = (r'T-compl $> %.1f$' % (100.*tiling_compl_cut))+(r'\%')
	line3 = (r'Z-compl $> %.1f$' % (100.*redshift_compl_cut))+(r'\%')
	line4 = (r'$%.2f < z < %.2f$' % (z_min, z_max))
	ax.text(0.61,0.55, s=(line1+'\n'+line2+'\n'+line3+'\n'+line4), transform=ax.transAxes, fontsize=16, bbox={'facecolor':'white', 'pad':8})
	
	ax.set_xlabel(r'$z$', size=18)
	ax.set_ylabel(r'$N(z)$', size=18)
	ax.set_xlim(0.1,2.)
	ax.legend()
	ax.grid(True,linestyle='-',alpha=0.3)
	f.savefig(output_path+'/hist_z.pdf',format='pdf')
	
	# OUTPUT
	
	return area
	
	run_info.close()

