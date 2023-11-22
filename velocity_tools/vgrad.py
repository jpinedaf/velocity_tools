import numpy as np
import sys
import astropy.units as u


def vfit(x, y, v, ev, distance=300.0*u.pc, resolution, boot=0):
"""
	It calculates the velocity gradient to a group of velocity 
	measuments.

	INPUTS:
		x:         	array of RA off-set positions in degrees.
		y:         	array of Dec off-set positions in degrees.
		v:         	Velocity at the position in km/s
		ev:   	  	Uncertainty in velocity (km/s)
	KEYWORD PARAMETERS:
		Width:     	FWHM of the weight function, in degrees.
		Boot:      	Set the number of bootstrap calculations
		distance:  	Distance to the region in pc.
	OUTPUTS:
		Vc:       	Mean centroid velocity in km/s
		Vc_err:   	Uncertainty of the mean centroid velocity in km/s
		Grad:     	Velocity gradient in km/s/pc
		Grad_Err: 	Uncertainty associated to the velocity gradient (km/s/pc)
		PosAng:   	Position angle of the fitted gradient, in degrees
		PAErr:    	Uncertainty of the position angle (degrees)
		ChiSq:    	Chi^2
"""

	Xold = np.copy(X)
	Yold = np.copy(Y)
	npts = X.shape
	wt = 1.0 / ev**2
	dtor = np.pi/180.*1000
	if (distance >= 0):
		conv = distance*dtor #X and Y are in degrees => pc
		resolution *= conv
		X = X*conv
		Y = Y*conv
	else:
    	print 'distance is negative!!'
    	return

	# Obtain total weight, and average (x,y,v) to create new variables (dx,dy,dv)
	# which provide a lower uncertainty in the fit.
	#
  	suma = np.sum(wt)
  	x_mean = np.sum(x * wt) / suma
  	y_mean = np.sum(y * wt) / suma
  	v_mean = np.sum(v * wt) / suma
  	# relative units
  	dx = (X - x_mean)
  	dy = (Y - y_mean)
  	dv = (V - v_mean)
  	M = [[np.sum(wt), np.sum(dx * wt), np.sum(dy * wt)], 
     	[np.sum(dx * wt), np.sum(dx**2 * wt), np.sum(dx * dy * wt)], 
     	[np.sum(dy * wt), np.sum(dx * dy * wt), np.sum(dy**2 * wt)]]
  	#print M
  	try:
    	covar = np.linalg.inv(M)
  	except IOError:
	    import sys
    	sys.exit('Singular matrix: no solution returned')
  	coeffs = dot(covar,[[sum(dv * wt)], [sum(dx * dv * wt)],[sum(dy * dv * wt)]])
  	#
  	errx= sqrt(covar[1, 1])
  	erry= sqrt(covar[2, 2])
  	#
  	gx = coeffs[1][0]
  	gy = coeffs[2][0]
  	# best fit mean V_lsr and predicted velocities
  	vc = coeffs[0] + v_mean
  	vp = coeffs[0] + coeffs[1] * dx + coeffs[2] * dy
  	grad = np.sqrt(coeffs[1]**2 + coeffs[2]**2)
  	posang = np.degrees(np.arctan2(gy, -gx))# * 180 / pi
  	# chi^2_red
  	red_chisq = np.sum((dv - vp)**2 * wt) / (len(dv) - 3.)

  	if boot > 0:
		# Calculate uncertainty by randomly changing the velocities 
  		# and recalculate the fit
    	grad_list = np.zeros(boot + 1)
    	posang_list = np.zeros(boot + 1)
    	vc_list = np.zeros(boot + 1)
    	for i in range(boot):
      		v_noise = np.random.normal(0., 1., npts) * ev
      	# Call again this function
      		grad_aux,grad_err_aux,posang_aux,paerr_aux,chisq_aux,vc_aux,vc_err_aux = \
        	VFit_OverAll(Xold,Yold,V+v_noise,Sigma_v,distance,resolution/conv,boot=0)
      		grad_list[i] = grad_aux
      		posang_list[i] = posang_aux
      		vc_list[i] = vc_aux
    	vc_err = np.std(vc_list)
    	grad_err = np.std(grad_list)
    	paerr = np.std(posang_list)
  	else:
    	vc_err = 0.
    	grad_err = sqrt((gx*errx)**2+(gy*erry)**2)/grad
	    grad_err = sqrt((gx*errx)**2+(gy*erry)**2+2*gx*gy*covar[2,1])/grad
    	paerr = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry**2+
                         (gy/(gx**2+gy**2))**2*errx**2)
    	paerr = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry**2+
                         (gy/(gx**2+gy**2))**2*errx**2-2*gx*gy/(gx**2+gy**2)**2*covar[2,1])
	chisq = red_chisq
  	vp += v_mean
	# Restore X and Y
  	X=copy(Xold)
  	Y=copy(Yold)

  	return grad[0], grad_err[0], posang, paerr, chisq, vc, vc_err