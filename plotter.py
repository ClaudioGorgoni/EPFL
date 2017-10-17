import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.io import ascii
from matplotlib import rc
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('xtick', labelsize=16) 
rc('ytick', labelsize=16) 





# THE ANGULAR CORRELATION FUNCTION
def angular(x, xi, err_x, output_path, plot_label, yaxis_scale, show_fig):
			
	plt.figure()
	
	if yaxis_scale == 'r2':
		r2     = x*x
		xi	  *= r2
		err_x *= r2
		s      = (r'$\theta^2 w(\theta)$')
		
		plt.ylim(-0.15,1.5)
		
	elif yaxis_scale == 'r':
		xi	  *= x
		err_x *= x
		s      = (r'$\theta w(\theta)$')
		
		plt.ylim(-0.025,0.25)
		
	elif yaxis_scale == 'lin':
		s	    = (r'$w(\theta)$')
		
		plt.ylim(-0.005,0.05)
		
	else:
		print '\n'
		print 'ERROR : yaxis_scale input is not recognizable.'
		print 'Must be one between : r2, r or lin.\n'
		sys.exit(0)
	
	plt.errorbar(x, xi, yerr = err_x, fmt='bo',label=plot_label)
	plt.xlabel(r'$\theta$ [deg]',size=18)
	plt.ylabel(s,size=18)
	plt.xscale('log')
	plt.grid(True,linestyle='-',alpha=0.3)
	plt.legend(loc='upper right')
	
	plt.savefig(output_path+'/angular.pdf', bbox_inches='tight', format='pdf')

#--- END DEF ---
		
		
		
		
		
# 3D CORR. FUNC. IN (pi,sigma) BINNING
def _3D_ps(x1, x2, xi, err_x, output_path, plot_label, show_fig):

	y = np.unique(x1) # pi
	x = np.unique(x2) # sigma
	Z = xi.reshape((len(x),len(y)))
	
	y = np.hstack((-y[::-1],y))		    # (-pi,pi)
	x = np.hstack((-x[::-1],x))		    # (-sigma,sigma)

	X,Y = np.meshgrid(x,y)

	Z = np.rot90(Z,1)
	Z = np.hstack((np.fliplr(Z),Z))
	Z = np.vstack((Z,np.flipud(Z)))
	
	plt.figure()
	Z  = gaussian_filter(Z,2.)
	
	if max(y) > 50.:
	
		s2  = X*X + Y*Y
		Z   = s2*Z
	
		Z  = np.where(Z>=150.,150.,Z)
		Z  = np.where(Z<=-100.,-100.,Z)
		
		levels = np.linspace(-100.,150.,11)
		cf = plt.contourf(X,Y,Z, levels = levels ,cmap='jet', origin='lower')
		bar = plt.colorbar(cf)
		bar.set_label(r'$r^2\xi(r_{\perp},r_{\parallel})$',size=18)
		
	else:
	
		Z  = np.where(Z>=30.,30.,Z)
		Z  = np.where(Z<=0.,0.01,Z)
		
		levels = [0.008] + [i for i in np.logspace(-2,1,10)] + [30.]
		
		cf = plt.contourf(X,Y,Z, levels=levels, norm = LogNorm() ,cmap='jet', origin='lower')
		
		#plt.contour(cf, levels = [0.25,0.5,1,2], colors='k',linewidth=2.)
		
		bar = plt.colorbar(cf)
		bar.set_label(r'$\xi(r_{\perp},r_{\parallel})$',size=18)
	
	
	plt.xlabel(r'$r_{\perp}$ [$h^{-1}$Mpc]',size=18) 
	plt.ylabel(r'$r_{\parallel}$ [$h^{-1}$Mpc]',size=18)
	
	
	plt.savefig(output_path+'/3D_ps.'+plot_label+'.pdf', bbox_inches='tight', format='pdf')
				
#--- END DEF ---
		
		
		
		
		
# 3D CORR. FUNC. IN (r,mu) BINNING.
# OUTPUT THE MONOPOLE AND THE QUADRUPOLE
def _3D_rm(x1, x2, xi, err_x, chunk, tiling_compl_cut, redshift_compl_cut, z_min, z_max, weights, yaxis_scale, compare_SHAM, 
		  output_path, plot_label, show_fig):

	r  = np.unique(x2)
	Nr = len(r)

	mu  = np.unique(x1)
	Nm  = len(mu)
	dmu = mu[1]-mu[0]
	
	# The Legendre polynomial of degree 2
	P2_mu = 0.5*(3.*mu*mu-1.)
	
	xi0     = np.zeros(Nr)
	err_xi0 = np.zeros(Nr)
	xi2     = np.zeros(Nr) 
	err_xi2 = np.zeros(Nr)
	
	for i in range(Nr):
		
		Nmi = Nm*i
		
		for k in range(Nm):
		
			xi_ik  = xi[ Nmi + k ]
			err_ik = err_x[ Nmi + k ]
			P2_k   = P2_mu[k]
			
			xi0[i] 	   += dmu*xi_ik
			err_xi0[i] += err_ik*err_ik
			
			xi2[i] 	   += 5.*dmu*xi_ik*P2_k
			err_xi2[i] += (err_ik*P2_k)*(err_ik*P2_k)
		
		
		err_xi0[i] = dmu*np.sqrt(err_xi0[i])
		err_xi2[i] = 5.*dmu*np.sqrt(err_xi2[i])
	

	tmp = np.transpose(np.vstack((r, xi0, xi2, err_xi0, err_xi2)))
	ascii.write(tmp,output_path+'/3D_rm.dat', format='no_header', 
			names = ['r','xi0','xi2','err_xi0','err_xi2'],
			formats={'r':'%.6f', 'xi0':'%.6f', 'xi2':'%.6f', 'err_xi0':'%.6f', 'err_xi2':'%.6f'})
			
			
	x_m, xi_m, err_x_m = np.loadtxt('data/mock.dat',delimiter=' ', usecols = (0,1,2), unpack = True)
	
	f0,ax0 = plt.subplots()
	f2,ax2 = plt.subplots()
	
	if yaxis_scale == 'r2':
		r2      = r*r
		xi0	    *= r2
		xi2		*= r2
		err_xi0 *= r2
		err_xi2 *= r2
		xi_m    *= x_m*x_m
		err_x_m *= x_m*x_m
		s0	    = (r'$s^2\xi_0(s)$')
		s2	    = (r'$s^2\xi_2(s)$')
		
		ax0.set_ylim(-100,150)
		ax2.set_ylim(-100,150)
		
	elif yaxis_scale == 'r':
		xi0	    *= r
		xi2		*= r
		err_xi0 *= r
		err_xi2 *= r
		xi_m    *= x_m
		err_x_m *= x_m
		s0	    = (r'$s\xi_0(s)$')
		s2	    = (r'$s\xi_2(s)$')
		
	elif yaxis_scale == 'lin':
		s0	    = (r'$\xi_0(s)$')
		s2	    = (r'$\xi_2(s)$')
		
		ax0.set_yscale('symlog',linthreshy= 0.01)
		ax0.set_ylim(-1.e-2, 5.)
		ax0.set_yticks([-1.e-2, -5.e-3, 0., 5e-3, 1e-2, 1e-1, 1e0])
		
		ax2.set_yscale('symlog',linthreshy= 0.01)
		
	else:
		print '\n'
		print 'ERROR : yaxis_scale input is not recognizable.'
		print 'Must be one between : r2, r or lin.\n'
		sys.exit(0)
	
	
	if compare_SHAM:
		ax0.plot(x_m, xi_m,'m-',linewidth=2, label = r'SHAM with $z=0.8$ $b=1.4$')

	
	ax0.errorbar(r, xi0, yerr = err_xi0, fmt='bo', label=plot_label)
	ax2.errorbar(r, xi2, yerr = err_xi2, fmt='ro', label=plot_label)	
	
	line1 = (r'Using chunk : '+chunk)
	line2 = (r'T-compl $> %.1f$' % (100.*tiling_compl_cut))+(r'\%')
	line3 = (r'Z-compl $> %.1f$' % (100.*redshift_compl_cut))+(r'\%')
	line4 = (r'$%.2f < z < %.2f$' % (z_min, z_max))	
	
	if len(weights) == 0:
		line5 = (r'Using no weights')
	
	else:
		line5 = r'Using weights : '
		if 1 in weights:
		
			line5 += r'compl'
			if len(weights) > 1:
				line5 += r', '

		if 2 in weights:
			
			if len(weights) == 1:
				line5 += r'sys'
			elif len(weights) == 2:
				if 1 in weights:
					line5 += r'sys'
				else:
					line5 += r'sys, '
			elif len(weights) == 3:
				line5 += r'sys, '
				
		if 3 in weights:
			line5 += r'FKP'
		
	
	ax0.text(0.52,0.6, s=(line1+'\n'+line2+'\n'+line3+'\n'+line4+'\n'+line5), transform=ax0.transAxes, 
			fontsize=14, bbox={'facecolor':'white', 'pad':8})
	ax2.text(0.52,0.06, s=(line1+'\n'+line2+'\n'+line3+'\n'+line4+'\n'+line5), transform=ax2.transAxes, 
			fontsize=14, bbox={'facecolor':'white', 'pad':8})
	
	ax0.set_xlabel(r'$s$ [$h^{-1}$Mpc]',size=18)
	ax2.set_xlabel(r'$s$ [$h^{-1}$Mpc]',size=18)
	ax0.set_ylabel(s0,size=18)
	ax2.set_ylabel(s2,size=18)
	
	ax0.grid(True,linestyle='-',alpha=0.3)
	ax2.grid(True,linestyle='-',alpha=0.3)
	
	ax0.legend(loc='upper right',fontsize=14)
	ax2.legend(loc='upper right',fontsize=14)

	f0.savefig(output_path+'/3D_rm.x0.pdf', bbox_inches='tight', format='pdf')
	f2.savefig(output_path+'/3D_rm.x2.pdf', bbox_inches='tight', format='pdf')

#--- END DEF ---



