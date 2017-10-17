
The command to execute for launching the program :
	
	python ELG_correlation.py

#############################################
	
CUTE parameters:
---------------

- For ELG clustering the only relevant correlation types are:

	#    "angular"        -> Angular correlation function
	#    "3D_ps"          -> 3D correlation binning in sigma-pi
	#    "3D_rm"          -> 3D correlation binning in r-mu (gives the monopole and the quadrupole)

- correlation_estim = 'LS' # w=(DD-2*DR+RR)/RR (Landy & Szalay estimator)


- log_bin: Set to 1 to use logarithmic binning. Only relevant for the angular 2PCF,
		   and 3D 2PCF (r-mu). 0 otherwise.


- n_logint: Number of bins in r or theta per decade when logaritmic binning is activated.


- dim1_max: Maximum scale to which the 2PCF will be computed in the first dimension. This corresponds to:
	      - The maximum angular separation in the case of the angular 2PC angular cross-correlations and the "full" correlation function.
		  - The maximum value of r for the monopole and 3D correlation function (in r-mu).
	      - The maximum redshift separation in the case of the radial correlation function.
		  - The maximum radial separation (pi) in the case of the 3D correlation function (in pi-sigma)


- dim1_nbin: Number of bins in the first dimension


- dim2_max: Maximum scale to which the 2PCF will be computed in the second dimension. This corresponds to:
		  - The maximum transverse separation (sigma) in the case of the 3D correlation function (in pi-sigma).
		  	Note that in the case of r-mu binning, mu always goes from 0 to 1.
          - The maximum redshift separation in the case of the "full" correlation function.


##########################################

CUTE output:
------------

The output file contains 6 columns for the radial, angular and monopole correlation functions with
             x   xi(x)   error(x)   DD(x)   DR(x)   RR(x) 
where x is Dz, theta or r. The third column is the Poisson error calculated from DD, DR and RR.

For the 3-D correlation functions the output file has 7 columns with
   x1   x2   xi(x1,x2)   error(x1,x2)   DD(x1,x2)   DR(x1,x2)   RR(x1,x2)
where (x1,x2) is either (pi,sigma) or (mu,r).
				 
