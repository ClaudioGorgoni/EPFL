import matplotlib.pyplot as plt 
import numpy as np
import subprocess
import os
import sys
import ConfigParser

config = ConfigParser.ConfigParser()
config.read("ELG_config.ini")

path_CUTE      = config.get('CUTE_param', 'cute_path')

def launch_CUTE(out_catalog, out_random, output_path, correlation_type,
				omega_m, omega_lam, w, log_bin, n_logint, dim1_max,
			 	dim1_nbin, dim2_max, dim2_nbin, dim3_min, dim3_max,
			 	dim3_nbin):
	
	out_cute = output_path+'/'+correlation_type

	correlation_estim = 'LS' # w=(DD-2*DR+RR)/RR (Landy & Szalay estimator)
		
	content= ('# input-output files and parameters\n'
			'data_filename= '+out_catalog+'.dat \n'
			'random_filename= '+out_random+'.dat \n'
			'input_format= 2\n'
			'mask_filename= none\n'
			'z_dist_filename= none\n'
			'output_filename= '+out_cute+'.dat \n'
			'num_lines= all\n\n'
			'# estimation parameters\n'
			'corr_type= '+correlation_type+'\n'
			'corr_estimator= '+correlation_estim+'\n'
			'np_rand_fact= 10\n\n'
			'# cosmological parameters\n'
			'omega_M= '+str(omega_m)+'\n'
			'omega_L= '+str(omega_lam)+'\n'
			'w= '+str(w)+'\n\n'
			'# binning\n'
			'log_bin= '+str(log_bin)+'\n'
			'n_logint= '+str(n_logint)+'\n'
			'dim1_max= '+str(dim1_max)+'\n'
			'dim1_nbin= '+str(dim1_nbin)+'\n'
			'dim2_max= '+str(dim2_max)+'\n'
			'dim2_nbin= '+str(dim2_nbin)+'\n'
			'dim3_min= '+str(dim3_min)+'\n'
			'dim3_max= '+str(dim3_max)+'\n'
			'dim3_nbin= '+str(dim3_nbin)+'\n\n'
			'# pixels for radial correlation\n'
			'radial_aperture= 7.0\n\n'
			'# pm parameters\n'
			'use_pm= 0\n'
			'n_pix_sph= 2048'
			)

	param = open('param.ini', 'w')
	param.write(content)
	param.close()
	
	print correlation_type+" correlation launched..."
	
	subprocess.call(path_CUTE+' param.ini &> '+output_path+'/output_CUTE.log', shell=True)
	
	print "done !\n"

	if correlation_type in ['angular', 'radial', 'monopole']:
		x, xi, err_x = np.loadtxt(out_cute+'.dat', delimiter=' ', usecols = (0,1,2), unpack = True)
		return x, xi, err_x
	elif correlation_type in ['3D_ps', '3D_rm']:
		x1, x2, xi, err_x = np.loadtxt(out_cute+'.dat', delimiter=' ', usecols = (0,1,2,3), unpack = True)
		return x1, x2, xi, err_x
	elif correlation_type in ['angular_cross', 'full']:
		x1, x2, x3, xi, err_x = np.loadtxt(out_cute+'.dat', delimiter=' ', usecols = (0,1,2,3,4), unpack = True)
		return x1, x2, x3, xi, err_x
