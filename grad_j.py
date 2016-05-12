import pyfits
from pylab import *
from numpy import *
import numpy.random  

def fitstoarrays(ffile,fmask):

  fitsfile = pyfits.open(ffile)
  data = fitsfile[0].data

  header = pyfits.getheader(ffile)
  naxis1 = header['naxis1']
  naxis2 = header['naxis2']
  cdelt1 = header['cdelt1']
  cdelt2 = header['cdelt2']
  crpix1 = header['crpix1']
  crpix2 = header['crpix2']
  crval1 = header['crval1']
  crval2 = header['crval2']

  X = zeros(data.shape)
  Y = zeros(data.shape)

  for j in range(data.shape[0]):
    for i in range(data.shape[1]):
      X[j,i] = (1+i)*cdelt1
      Y[j,i] = (1+j)*cdelt2

  maskfile = pyfits.open(fmask)
  datam = maskfile[0].data

  mask = datam!=0
  #Z = (X**2+Y**2)

  return X[mask],Y[mask],data[mask]

def vfit_j(X,Y,V,Sigma_v, distance, resolution,width=2., nmin=7):

# INPUTS:
#      X:         Off-Set in degrees.
#      Y:         Off-Set in degrees.
#      V:         Velocity at the position in km/s
###      IntV:      Integral intensity at the position
#      Sigma_V:   Uncertainty in velocity (km/s)
#      distance:  Distance to the region in pc.
#      resolution:Pixel resolution in degrees.
# KEYWORD PARAMETERS:
#      Width:     FWHM of the weight function, in resolution units. Default=2.
#      NMin:      Minimum number of neighbors to perfom fit. Default 7.
#
# OUTPUTS:
#       Grad:     Velocity gradient in km/s/pc
#       Grad_Err: Uncertainty associated to the velocity gradient (km/s/pc)
#       PosAng:   Position angle of the fitted gradient, in degrees
#       PAErr:    Uncertainty of the position angle (degrees)
#       ChiSq:    Chi^2 at each position

  Xold = copy(X)
  Yold = copy(Y)
  nterms = 2.
  #mask = (IntV!=0)
  grad    = X*0.0
  grad_err= X*0.0
  posang  = X*0.0
  paerr   = X*0.0
  chisq   = X*0.0
  #npts = X.shape
  npts = len(X)
  count = zeros(npts)
  erry = zeros(npts)
  errx = zeros(npts)
  weight = 1/(Sigma_v**2)
  dtor = pi/180.*1000
  if (distance >= 0):
    conv = distance*dtor #X and Y are in degrees => pc
    resolution *= conv
    X = X*conv
    Y = Y*conv
  else:
    print 'distance is negative!!'
    return

#	Choose the points to pass to the the fitting procedure
#       Pick a point, and then decide if there are at least Nmin points
#       adjacent to it (within the grid spacing)

  gs = abs(2.01*resolution)	# 2 x grid spacing
  line_num = 0
  #print 'nmin',nmin
  for i in range(npts):#[0]):
    #for j in range(npts[1]):
    counter = 0.
    dx = abs(X-X[i])
    dy = abs(Y-Y[i])
    counter = len( X[(dx<=gs) * (dy<=gs)] )
    if counter >= nmin:
      count[i] = counter
    else:
      count[i] = 0

#	For each point with at least 6 neighbors, 
#	create arrays with all their relevant parameters

  for i in range(npts):#[0]):
    #for j in range(npts[1]):
      if (count[i] > 0):
      #wt      = zeros(npts)	; The new weights
        dx      = (X[i] - X)/resolution
        #dx      = (X[i] - X)/0.00633/conv#/resolution
        dy      = (Y[i] - Y)/resolution
        dz      = sqrt(dx**2 + dy**2)
        coeff   = exp(-(dz**2)/(2*(width/2.354)**2))
        index   = (coeff >= 1e-4) * isfinite(weight)
# New weights are old weights times a Gaussian to weight nearby points
# higher than more distant ones
        wt      = weight[index]*coeff[index]	
        uniq_x  = X[index]
        uniq_y  = Y[index]
        vcen    = V[index]
        M = matrix([[sum(wt), sum(uniq_x*wt), sum(uniq_y*wt)], 
           [sum(uniq_x*wt), sum(uniq_x**2*wt), sum(uniq_x*uniq_y*wt)], 
           [sum(uniq_y*wt), sum(uniq_x*uniq_y*wt), sum(uniq_y**2*wt)]])
        from scipy import linalg
	try:
          covar = linalg.inv(M)
	except LinAlgError:
#if (q gt 0) then begin
          print 'Singular matrix: no solution returned'
          grad[i] = nan
          continue
        coeffs = dot(covar,[[sum(vcen*wt)], [sum(uniq_x*vcen*wt)],[sum(uniq_y*vcen*wt)]])
        errx[i]= sqrt(covar[1, 1])
        erry[i]= sqrt(covar[2, 2])
        gx = coeffs[1]
        gy = coeffs[2]
        red_chisq = sum((vcen-(coeffs[0]+coeffs[1]*uniq_x+ \
                             coeffs[2]*uniq_y))**2* \
                      wt)/(len(vcen)-3)
        grad[i]     = sqrt(coeffs[1]**2+coeffs[2]**2)
        grad_err[i] = sqrt(errx[i]**2+erry[i]**2)
        posang[i]   = arctan2(gy,-gx)*180/pi
        paerr[i]    = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry[i]**2+ \
                              (gy/(gx**2+gy**2))**2*errx[i]**2)
        chisq[i]    = red_chisq
      else:
        grad[i]     = nan
        grad_err[i] = nan
        posang[i]   = nan
        paerr[i]    = nan
        chisq[i]    = nan

#Restore X and Y
  X=copy(Xold)
  Y=copy(Yold)
  return grad,grad_err,posang,paerr,chisq



def VFit_OverAll( X, Y, V, Sigma_v, distance, resolution, boot=0):

# INPUTS:
#      X:         Off-Set in degrees.
#      Y:         Off-Set in degrees.
#      V:         Velocity at the position in km/s
#      IntV:      Integral intensity at the position
#      Sigma_V:   Uncertainty in velocity (km/s)
#      distance:  Distance to the region in pc.

# KEYWORD PARAMETERS:
#      Width:     FWHM of the weight function, in degrees.
#      Boot:      Set the number of bootstrap calculations

# OUTPUTS:
#       Vc:       Mean centroid velocity in km/s
#       Vc_err:   Uncertainty of the mean centroid velocity in km/s
#       Grad:     Velocity gradient in km/s/pc
#       Grad_Err: Uncertainty associated to the velocity gradient (km/s/pc)
#       PosAng:   Position angle of the fitted gradient, in degrees
#       PAErr:    Uncertainty of the position angle (degrees)
#       ChiSq:    Chi^2

  Xold = copy(X)
  Yold = copy(Y)
  npts = X.shape
  #mask = (IntV!=0)
  wt = 1/(Sigma_v**2)
  dtor = pi/180.*1000
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
  suma  = sum(wt)
  x_mean=sum(X*wt)/suma
  y_mean=sum(Y*wt)/suma
  v_mean=sum(V*wt)/suma
  dx      = (X-x_mean)#[mask]  # remove mean value from inputs 
  dy      = (Y-y_mean)#[mask]  # to reduce fit uncertainties
  dv      = (V-v_mean)#[mask]  #
  M = [[sum(wt),    sum(dx*wt),    sum(dy*wt)], 
     [sum(dx*wt), sum(dx**2*wt),  sum(dx*dy*wt)], 
     [sum(dy*wt), sum(dx*dy*wt), sum(dy**2*wt)]]
  #print M
  try:
    covar = linalg.inv(M)
  except IOError:
    import sys
    sys.exit('Singular matrix: no solution returned')
  coeffs = dot(covar,[[sum(dv*wt)], [sum(dx*dv*wt)],[sum(dy*dv*wt)]])
  #print covar
  errx= sqrt(covar[1, 1])
  erry= sqrt(covar[2, 2])
  #print coeffs
  gx = coeffs[1][0]
  gy = coeffs[2][0]
  #print gx,gy
  vc = coeffs[0]+v_mean
  vp = coeffs[0]+coeffs[1]*dx+coeffs[2]*dy
  grad     = sqrt(coeffs[1]**2+coeffs[2]**2)
  posang   = arctan2(gy, -gx)*180/pi
  #print -gx,gy,posang
  red_chisq = sum( (dv-vp)**2*wt)/(len(dv)-3.)

  if boot > 0:
    grad_list   = zeros(boot+1)
    posang_list = zeros(boot+1)
    vc_list     = zeros(boot+1)
    for i in range(boot):
      v_noise=numpy.random.normal(0.,1.,npts)*Sigma_v
      grad_aux,grad_err_aux,posang_aux,paerr_aux,chisq_aux,vc_aux,vc_err_aux = \
        VFit_OverAll(Xold,Yold,V+v_noise,Sigma_v,distance,resolution/conv,boot=0)
      grad_list[i]   =grad_aux
      posang_list[i] =posang_aux
      vc_list[i]     =vc_aux
    vc_err   = std(vc_list)
    grad_err = std(grad_list)
    paerr    = std(posang_list)
  else:
    vc_err   = 0.
    grad_err = sqrt((gx*errx)**2+(gy*erry)**2)/grad
    #print 'grad_err1', grad_err
    grad_err = sqrt((gx*errx)**2+(gy*erry)**2+2*gx*gy*covar[2,1])/grad
    #print 'grad_err2', grad_err
    paerr    = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry**2+
                         (gy/(gx**2+gy**2))**2*errx**2)
    #print 'paerr1',paerr
    paerr    = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry**2+
                         (gy/(gx**2+gy**2))**2*errx**2-2*gx*gy/(gx**2+gy**2)**2*covar[2,1])
    #print 'paerr2',paerr
  chisq    = red_chisq
  vp += v_mean
# Restore X and Y
  X=copy(Xold)
  Y=copy(Yold)

  return grad[0],grad_err[0],posang,paerr,chisq,vc,vc_err



#-------------


def vfit_j2(X,Y,V,Sigma_v, distance, resolution,width=2., nmin=7):

# INPUTS:
#      X:         Off-Set in degrees.
#      Y:         Off-Set in degrees.
#      V:         Velocity at the position in km/s
###      IntV:      Integral intensity at the position
#      Sigma_V:   Uncertainty in velocity (km/s)
#      distance:  Distance to the region in pc.
#      resolution:Pixel resolution in degrees.
# KEYWORD PARAMETERS:
#      Width:     FWHM of the weight function, in resolution units. Default=2.
#      NMin:      Minimum number of neighbors to perfom fit. Default 7.
#
# OUTPUTS:
#       Grad:     Velocity gradient in km/s/pc
#       Grad_Err: Uncertainty associated to the velocity gradient (km/s/pc)
#       PosAng:   Position angle of the fitted gradient, in degrees
#       PAErr:    Uncertainty of the position angle (degrees)
#       ChiSq:    Chi^2 at each position

  Xold = copy(X)
  Yold = copy(Y)
  nterms = 2.
  #mask = (IntV!=0)
  grad    = X*0.0
  grad_err= X*0.0
  posang  = X*0.0
  paerr   = X*0.0
  chisq   = X*0.0
  #npts = X.shape
  npts = len(X)
  count = zeros(npts)
  erry = zeros(npts)
  errx = zeros(npts)
  GRy = zeros(npts)
  GRx = zeros(npts)
  weight = 1/(Sigma_v**2)
  dtor = pi/180.*1000
  if (distance >= 0):
    conv = distance*dtor #X and Y are in degrees => pc
    resolution *= conv
    X = X*conv
    Y = Y*conv
  else:
    print 'distance is negative!!'
    return

#	Choose the points to pass to the the fitting procedure
#       Pick a point, and then decide if there are at least Nmin points
#       adjacent to it (within the grid spacing)

  gs = abs(2.01*resolution)	# 2 x grid spacing
  line_num = 0
  #print 'nmin',nmin
  for i in range(npts):#[0]):
    #for j in range(npts[1]):
    counter = 0.
    dx = abs(X-X[i])
    dy = abs(Y-Y[i])
    counter = len( X[(dx<=gs) * (dy<=gs)] )
    if counter >= nmin:
      count[i] = counter
    else:
      count[i] = 0

#	For each point with at least 6 neighbors, 
#	create arrays with all their relevant parameters

  for i in range(npts):#[0]):
    #for j in range(npts[1]):
      if (count[i] > 0):
      #wt      = zeros(npts)	; The new weights
        dx      = (X[i] - X)/resolution
        #dx      = (X[i] - X)/0.00633/conv#/resolution
        dy      = (Y[i] - Y)/resolution
        dz      = sqrt(dx**2 + dy**2)
        coeff   = exp(-(dz**2)/(2*(width/2.354)**2))
        index   = (coeff >= 1e-4) * isfinite(weight)
# New weights are old weights times a Gaussian to weight nearby points
# higher than more distant ones
        wt      = weight[index]*coeff[index]
        wt2     = weight[index]*coeff[index]**2
        uniq_x  = X[index]
        uniq_y  = Y[index]
        vcen    = V[index]
        M = matrix([[sum(wt), sum(uniq_x*wt), sum(uniq_y*wt)], 
           [sum(uniq_x*wt), sum(uniq_x**2*wt), sum(uniq_x*uniq_y*wt)], 
           [sum(uniq_y*wt), sum(uniq_x*uniq_y*wt), sum(uniq_y**2*wt)]])
        M2 = matrix([[sum(wt2), sum(uniq_x*wt2), sum(uniq_y*wt2)], 
           [sum(uniq_x*wt2), sum(uniq_x**2*wt2), sum(uniq_x*uniq_y*wt2)], 
           [sum(uniq_y*wt2), sum(uniq_x*uniq_y*wt2), sum(uniq_y**2*wt2)]])
        from scipy import linalg
	try:
          Minv = linalg.inv(M)
	except LinAlgError:
#if (q gt 0) then begin
          print 'Singular matrix: no solution returned'
          grad[i] = nan
          continue

        coeffs = dot(Minv,[[sum(vcen*wt)], [sum(uniq_x*vcen*wt)],[sum(uniq_y*vcen*wt)]])
	covar = dot(dot(Minv,M2),Minv)

        errx[i]= sqrt(covar[1, 1])
        erry[i]= sqrt(covar[2, 2])
        gx = coeffs[1]
        gy = coeffs[2]
	GRx[i] = gx
	GRy[i] = gy
        red_chisq = sum((vcen-(coeffs[0]+coeffs[1]*uniq_x+ \
                             coeffs[2]*uniq_y))**2* \
                      wt)/(len(vcen)-3)
        grad[i]     = sqrt(coeffs[1]**2+coeffs[2]**2)
        grad_err[i] = sqrt(errx[i]**2+erry[i]**2)
        posang[i]   = arctan2(gy,-gx)*180/pi
        paerr[i]    = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry[i]**2+ \
                              (gy/(gx**2+gy**2))**2*errx[i]**2)
        chisq[i]    = red_chisq
      else:
        grad[i]     = nan
        grad_err[i] = nan
        posang[i]   = nan
        paerr[i]    = nan
        chisq[i]    = nan

#Restore X and Y
  X=copy(Xold)
  Y=copy(Yold)
  return grad,grad_err,posang,paerr,chisq,GRx,GRy,errx,erry

#-------------Full with second derivates


def vfit_secondd(X,Y,V,Sigma_v, distance, resolution,width=2., nmin=7):

# INPUTS:
#      X:         Off-Set in degrees.
#      Y:         Off-Set in degrees.
#      V:         Velocity at the position in km/s
###      IntV:      Integral intensity at the position
#      Sigma_V:   Uncertainty in velocity (km/s)
#      distance:  Distance to the region in pc.
#      resolution:Pixel resolution in degrees.
# KEYWORD PARAMETERS:
#      Width:     FWHM of the weight function, in resolution units. Default=2.
#      NMin:      Minimum number of neighbors to perfom fit. Default 7.
#
# OUTPUTS:
#       Grad:     Velocity gradient in km/s/pc
#       Grad_Err: Uncertainty associated to the velocity gradient (km/s/pc)
#       PosAng:   Position angle of the fitted gradient, in degrees
#       PAErr:    Uncertainty of the position angle (degrees)
#       ChiSq:    Chi^2 at each position

  Xold = copy(X)
  Yold = copy(Y)
  nterms = 2.
  #mask = (IntV!=0)
  grad    = X*0.0
  grad_err= X*0.0
  posang  = X*0.0
  paerr   = X*0.0
  chisq   = X*0.0
  #npts = X.shape
  npts = len(X)
  count = zeros(npts)
  erry = zeros(npts)
  errx = zeros(npts)
  weight = 1/(Sigma_v**2)
  dtor = pi/180.*1000
  if (distance >= 0):
    conv = distance*dtor #X and Y are in degrees => pc
    resolution *= conv
    X = X*conv
    Y = Y*conv
  else:
    print 'distance is negative!!'
    return

#	Choose the points to pass to the the fitting procedure
#       Pick a point, and then decide if there are at least Nmin points
#       adjacent to it (within the grid spacing)

  gs = abs(2.01*resolution)	# 2 x grid spacing
  line_num = 0
  #print 'nmin',nmin
  for i in range(npts):#[0]):
    #for j in range(npts[1]):
    counter = 0.
    dx = abs(X-X[i])
    dy = abs(Y-Y[i])
    counter = len( X[(dx<=gs) * (dy<=gs)] )
    if counter >= nmin:
      count[i] = counter
    else:
      count[i] = 0

#	For each point with at least 6 neighbors, 
#	create arrays with all their relevant parameters

  for i in range(npts):#[0]):
    #for j in range(npts[1]):
      if (count[i] > 0):
      #wt      = zeros(npts)	; The new weights
        dx      = (X[i] - X)/resolution
        #dx      = (X[i] - X)/0.00633/conv#/resolution
        dy      = (Y[i] - Y)/resolution
        dz      = sqrt(dx**2 + dy**2)
        coeff   = exp(-(dz**2)/(2*(width/2.354)**2))
        index   = (coeff >= 1e-4) * isfinite(weight)
# New weights are old weights times a Gaussian to weight nearby points
# higher than more distant ones
        wt      = weight[index]*coeff[index]
        wt2     = weight[index]*coeff[index]**2

        uniq_x  = X[index]
        uniq_y  = Y[index]
        vcen    = V[index]
        M = matrix([[sum(wt), sum(uniq_x*wt), sum(uniq_y*wt)], 
           [sum(uniq_x*wt), sum(uniq_x**2*wt), sum(uniq_x*uniq_y*wt)], 
           [sum(uniq_y*wt), sum(uniq_x*uniq_y*wt), sum(uniq_y**2*wt)]])
        M2 = matrix([[sum(wt2), sum(uniq_x*wt2), sum(uniq_y*wt2)], 
           [sum(uniq_x*wt2), sum(uniq_x**2*wt2), sum(uniq_x*uniq_y*wt2)], 
           [sum(uniq_y*wt2), sum(uniq_x*uniq_y*wt2), sum(uniq_y**2*wt2)]])
        from scipy import linalg
	try:
          Minv = linalg.inv(M)
	except LinAlgError:
#if (q gt 0) then begin
          print 'Singular matrix: no solution returned'
          grad[i] = nan
          continue

        coeffs = dot(Minv,[[sum(vcen*wt)], [sum(uniq_x*vcen*wt)],[sum(uniq_y*vcen*wt)]])
	covar = dot(dot(Minv,M2),Minv)

        errx[i]= sqrt(covar[1, 1])
        erry[i]= sqrt(covar[2, 2])
        gx = coeffs[1]
        gy = coeffs[2]
        red_chisq = sum((vcen-(coeffs[0]+coeffs[1]*uniq_x+ \
                             coeffs[2]*uniq_y))**2* \
                      wt)/(len(vcen)-3)
        grad[i]     = sqrt(coeffs[1]**2+coeffs[2]**2)
        grad_err[i] = sqrt(errx[i]**2+erry[i]**2)
        posang[i]   = arctan2(gy,-gx)*180/pi
        paerr[i]    = 180/pi*sqrt((gx/(gx**2+gy**2))**2*erry[i]**2+ \
                              (gy/(gx**2+gy**2))**2*errx[i]**2)
        chisq[i]    = red_chisq
      else:
        grad[i]     = nan
        grad_err[i] = nan
        posang[i]   = nan
        paerr[i]    = nan
        chisq[i]    = nan

#Restore X and Y
  X=copy(Xold)
  Y=copy(Yold)
  return grad,grad_err,posang,paerr,chisq



