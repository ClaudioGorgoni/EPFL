import numpy as np
import sys
import matplotlib.pyplot as plt
from cat import generate_elg_catalog
from rand import generate_elg_random
from cute import launch_CUTE
import ConfigParser
import time
from os.path import isfile, exists
from os import makedirs
import subprocess
import shutil
from astropy.time import Time
from plotter import angular, _3D_ps, _3D_rm


###INITIALIZATION###

start = time.time()

t   = Time(Time.now(), scale='utc')
mjd = "%.5f" % t.mjd

config = ConfigParser.ConfigParser()
config.read('ELG_config.ini')

chunk 			 = config.get('files', 'chunk')
generate_files   = config.getboolean('generate', 'generate_files')

if generate_files:

	if chunk not in ['eboss21','eboss22','eboss23']:

		info_line = 'ERROR : The chunk you selected is not recognizable.\n'
		info_line += 'Must be one between : eboss21, eboss22 or eboss23.\n'
		print info_line
		sys.exit(0)
		
	if chunk == 'eboss21':
		catalog_path = 'data/eboss21.v5_10_4.latest.fits'
		random_path  = 'data/eboss21.rands10000.unmsk.fits'

	elif chunk == 'eboss22':
		catalog_path = 'data/eboss22.v5_10_4.latest.fits'
		random_path  = 'data/eboss22.rands10000.unmsk.fits'

	elif chunk == 'eboss23':
		catalog_path = 'data/eboss23.v5_10_4.latest.fits'
		random_path  = 'data/eboss23.rands10000.test.unmsk.depthivarcut.fits'


	# store files in output/chunk/mjd
	output_path = 'output/'+chunk+'/'+mjd

	if not exists(output_path):
		makedirs(output_path)

	out_random  = output_path+'/random.'+chunk
	out_catalog = output_path+'/catalog.'+chunk
	out_dndz    = output_path+'/dndz.'+chunk
#---------------------------------------


# IMPORT PARAMETERS FROM CONFIG FILE
tmp = config.get('files', 'plates_to_rmv')

if len(tmp)>0:
	plates_to_rmv = [ int(l.strip()) for l in tmp.split(',') ]
else:
	plates_to_rmv = []

#-------------------------------------------------------------------
down_sample_fact   = config.getfloat(  'files', 'down_sample_fact')
if down_sample_fact < 1.0:
	down_sample_fact = 1.
	print "down_sample_fact is set to 1, because input parameter was < 1.\n"
	
show_fig           = config.getboolean('files', 'show_figures')
#-------------------------------------------------------------------
z_min              = config.getfloat(  'catalog_param', 'z_min')
z_max              = config.getfloat(  'catalog_param', 'z_max')
binwidth_z         = config.getfloat(  'catalog_param', 'binwidth_z')
comparat           = config.getboolean('catalog_param', 'comparat')
tmp          	   = config.get(		 'catalog_param', 'zwarning')
zwarning           = [ int(l.strip()) for l in tmp.split(',') ]
template_noqso     = config.getboolean('catalog_param', 'template_noqso')
add_qso 	       = config.getboolean('catalog_param', 'add_qso')
tiling_compl_cut   = config.getfloat(  'catalog_param', 'tiling_compl_cut')
redshift_compl_cut = config.getfloat(  'catalog_param', 'redshift_compl_cut')
tmp          	   = config.get(		 'catalog_param', 'weights')
if len(tmp)>0:
	weights        = [ int(l.strip()) for l in tmp.split(',') ]
else:
	weights 	   = []
#-------------------------------------------------------------------
H0                 = config.getfloat(  'cosmology', 'H0')
omega_m            = config.getfloat(  'cosmology', 'omega_m')
omega_lam          = config.getfloat(  'cosmology', 'omega_lam')
w                  = config.getfloat(  'cosmology', 'w')
#-------------------------------------------------------------------
run_CUTE		   = config.getboolean('CUTE_param', 'launch_CUTE')
yaxis_scale		   = config.get(		 'CUTE_param', 'yaxis_scale')
plot_label		   = config.get(		 'CUTE_param', 'plot_label')
compare_SHAM	   = config.getboolean('CUTE_param', 'compare_SHAM')
tmp				   = config.get(		 'CUTE_param', 'correlation_type')
if len(tmp)>0 and run_CUTE:
	correlation_type = [ int(l.strip()) for l in tmp.split(',') ]
elif len(tmp) == 0 and run_CUTE:
	sys.exit('Must choose a correlation type')
else:
	correlation_type = []
#-------------------------------------------------------------------



#------------------------------------------------
#GENERATE CATALOG OF GALAXIES AND RANDOM OBJECTS:
#------------------------------------------------

if generate_files:

	run_info = open(output_path+'/run.info','a')
	run_info.write('Created at MJD: '+mjd+'\n\n')
	run_info.close()
	
	print "\n"
	print "Generating ELG catalogs..."
	
	area = generate_elg_catalog(mjd, chunk, output_path, catalog_path, out_catalog, out_dndz, z_min, z_max, 
								binwidth_z, weights, zwarning, plates_to_rmv, template_noqso, add_qso, 
								comparat, tiling_compl_cut, redshift_compl_cut)
	
	generate_elg_random(mjd, chunk, area, output_path, random_path, out_random, out_catalog, out_dndz, 
						z_min, z_max, binwidth_z, weights, tiling_compl_cut, redshift_compl_cut, 
						down_sample_fact, H0, omega_m, omega_lam, w)
	
	print "ELG catalogs generated !\n"
	
else:

	out_catalog = config.get('generate', 'your_catalog')
	out_random  = config.get('generate', 'your_random')
	
	output_path  = config.get('CUTE_param', 'out_folder_cute')
	if not exists(output_path):
		makedirs(output_path)
		#sys.exit('The folder -out_folder_cute- you specified does not exists.')

	run_info = open(output_path+'/run.info','a')
	run_info.write('Running at MJD: '+mjd+'\n\n')
	run_info.close()



#-------------------------
# PLOT CORRELATION RESULTS
#-------------------------


if run_CUTE:


	CUTE_list = [out_catalog, out_random, output_path, correlation_type, \
				 omega_m, omega_lam, w, 0, 21, 200., 50, 200., 50, z_min, z_max, 50]
				 
				 
	run_info = open(output_path+'/run.info','a')

	info_line  = ('Cosmology:\n')
	info_line += ('H0      = %.4f\n' % H0)
	info_line += ('Omega_M = %.4f\n' % omega_m)
	info_line += ('Omega_L = %.4f\n' % omega_lam)
	info_line += ('w       = %.4f\n\n' % w)
	print info_line
	run_info.write(info_line)
	run_info.close()
	#------------------------
	
	
	if 1 in correlation_type:
	
		CUTE_list[3] = 'angular'
		
		CUTE_list[7]  = config.getint(	 'angular', 'log_bin')
		CUTE_list[8]  = config.getint(	 'angular', 'n_logint')
		CUTE_list[9]  = config.getfloat( 'angular', 'dim1_max')
		CUTE_list[10] = config.getint(	 'angular', 'dim1_nbin')
		
		
		x, xi, err_x = launch_CUTE(*CUTE_list)
		angular(x, xi, err_x, output_path, plot_label, yaxis_scale, show_fig)

	if 2 in correlation_type:
		
		if 1 in correlation_type:
			subprocess.call('rm '+output_path+'/output_CUTE.log', shell=True)
		
		CUTE_list[3] = '3D_ps'
		CUTE_list[7] = 0
		
		CUTE_list[9]  = config.getfloat( '3D_ps', 'dim1_max')
		CUTE_list[10] = config.getint(	 '3D_ps', 'dim1_nbin')
		CUTE_list[11] = config.getfloat( '3D_ps', 'dim2_max')
		CUTE_list[12] = config.getint(	 '3D_ps', 'dim2_nbin')
		
		x1, x2, xi, err_x = launch_CUTE(*CUTE_list)
		_3D_ps(x1, x2, xi, err_x, output_path, plot_label, show_fig)
			
	if 3 in correlation_type:
		
		if (1 in correlation_type) or (2 in correlation_type):
			subprocess.call('rm '+output_path+'/output_CUTE.log', shell=True)
		
		CUTE_list[3] = '3D_rm'
		
		CUTE_list[7]  = config.getint(	 '3D_rm', 'log_bin')
		CUTE_list[8]  = config.getint(	 '3D_rm', 'n_logint')
		CUTE_list[9]  = config.getfloat( '3D_rm', 'dim1_max')
		CUTE_list[10] = config.getint(	 '3D_rm', 'dim1_nbin')
		CUTE_list[11] = config.getfloat( '3D_rm', 'dim2_max')
		CUTE_list[12] = config.getint(	 '3D_rm', 'dim2_nbin')
		
		x1, x2, xi, err_x = launch_CUTE(*CUTE_list)
		_3D_rm(x1, x2, xi, err_x, chunk, tiling_compl_cut, redshift_compl_cut, z_min, z_max, weights, yaxis_scale, compare_SHAM, 
		 	  output_path, plot_label, show_fig)
		 	  
	end = time.time()
	print "Elapsed time : %.3f sec" % (end-start)
		 	  
	if show_fig:
		plt.show()


###### END #######


