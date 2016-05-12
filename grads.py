import pyfits
from pylab import *
from numpy import *
from grad_j import *
import aplpy as apl
from matplotlib import pylab
from matplotlib import rc
rc('text',usetex=True)


import os
import argparse
import pywcs

# create the parser
parser = argparse.ArgumentParser(description="Command line arguments.")

parser.add_argument('-coords', nargs=2)
parser.add_argument('-background', required=True)
parser.add_argument('-v0origfile', default='V0F_orig.fits')
parser.add_argument('-v0finalfile', default='V0F.fits')
parser.add_argument('-v0errfile', default='V0_err.fits')
parser.add_argument('-intfile', default = 'IntF.fits')
parser.add_argument('-sigmafile', default='SigmaF.fits')
parser.add_argument('-out', required=True)
parser.add_argument('-dist', type=float)
parser.add_argument('-size')
parser.add_argument('-sizey')
parser.add_argument('-dir', default=".")
parser.add_argument('-ncut')
#options
parser.add_argument('-cont',action='store_true')
parser.add_argument('-cont2')
parser.add_argument('-vedit', action='store_true', help='match colorscale min and max')
parser.add_argument('-itfield', default='fielditems.txt') #ie default = False
parser.add_argument('-cmin', default=0.05, const=0.05,nargs='?', type=float)
parser.add_argument('-cmax', default=1.0,const=1.0,nargs='?', type=float)
parser.add_argument('-cdelt', default=0.1,const=0.1,nargs='?', type=float)
parser.add_argument('-contdef', action='store_true')
parser.add_argument('-colorm', default='gray')
parser.add_argument('-bunit')
parser.add_argument('-reset', action='store_true')
parser.add_argument('-notnorth', action='store_true', help='update north')
parser.add_argument('-beam', action='store_true')
parser.add_argument('-shape', action='store_true')
parser.add_argument('-shift', action='store_true')
parser.add_argument('-items')

parser.add_argument('-peak', action='store_true', help='plot cross on peak position of fits, or of intensity map')

parser.add_argument('-sfile', default=False)

args = parser.parse_args()

if args.dir:
  os.system("mkdir "+args.dir)

def fitsshift(infile,outfile,intfile):
  sfile = pyfits.open(intfile)
  sdata = sfile[0].data
  sdata[isnan(sdata)]=0
  smax=sdata.max()
  #print 'smax',smax
  peakx = sdata.argmax()/(sdata.shape[1])
  peaky = sdata.argmax()%(sdata.shape[1])
  wcs = pywcs.WCS(sfile[0].header)
  sky = wcs.wcs_pix2sky([[peaky,peakx,0]], 0)
  #print sky
  afile = pyfits.open(infile, mode='update')
  adata = afile[0].data
  adata[isnan(adata)]=0
  wcs2 = pywcs.WCS(afile[0].header)
  pix = (wcs2.wcs_sky2pix(sky,0))[0][0:2]
  #print 'pixs', pix
  vatmax = adata[int(round(pix[1],0)),int(round(pix[0],0))]
  #print 
  adata[adata==0.0] = vatmax
  adata = adata - vatmax
  header = pyfits.getheader(infile)
  pyfits.writeto(outfile,adata,header)
  afile.close()

#create subimages
if args.coords and args.size:
  size = args.size
  if args.sizey:
    size += " "+args.sizey
  os.system("mSubimage "+args.v0origfile+" faux1a.fits "+args.coords[0]+" "+args.coords[1]+" "+size)
  f1file = 'faux1.fits'
  os.system("mSubimage "+args.v0errfile+" faux2.fits "+args.coords[0]+" "+args.coords[1]+" "+size)
  f2file = 'faux2.fits'
  os.system("mSubimage "+args.intfile+" faux3.fits "+args.coords[0]+" "+args.coords[1]+" "+size)
  f3file = 'faux3.fits'
  os.system("mSubimage "+args.v0finalfile+" faux4a.fits "+args.coords[0]+" "+args.coords[1]+" "+size)
  f4file = 'faux4.fits'
  os.system("mSubimage "+args.sigmafile+" faux5.fits "+args.coords[0]+" "+args.coords[1]+" "+size)
  f5file = 'faux5.fits'

  if args.shift:
    os.system("mSubimage "+args.background+" faux0a.fits "+args.coords[0]+" "+args.coords[1]+" "+size)
    f0file = 'faux0.fits'
    fitsshift("faux0a.fits",f0file,f3file)
  else:
    os.system("mSubimage "+args.background+" faux0.fits "+args.coords[0]+" "+args.coords[1]+" "+size)
    f0file = 'faux0.fits'
  
  fitsshift("faux1a.fits",f1file,f3file)
  #fitsshift("faux4a.fits",f4file,f3file)

  #if args.cont:
    #os.system("mSubimage "+args.cont+" faux2.fits "+args.coords[0]+" "+args.coords[1]+" "+args.size)
else:
  f0file = args.background
  f1file = args.v0origfile
  f2file = args.v0errfile
  f3file = args.intfile
  f4file = args.v0finalfile
  f5file = args.sigmafile

#deg1, deg2, vel = fitstoarrays('V0F_orig.fits','IntF.fits')
#deg1b, deg2b, vel_err = fitstoarrays('V0_err.fits','IntF.fits')
#header = pyfits.getheader('V0F_orig.fits')
deg1, deg2, vel = fitstoarrays(f1file,f3file)
#vel *= -1.
deg1b, deg2b, vel_err = fitstoarrays(f2file,f3file)
deg1c, deg2c, sigma = fitstoarrays(f5file,f3file)
deg1c, deg2c, intf = fitstoarrays(f3file,f3file)
deg1c, deg2c, dvel = fitstoarrays(f4file,f3file)

header = pyfits.getheader(f1file)

res = header['cdelt1']
crval1 = header['crval1']
crval2 = header['crval2']
crpix1 = header['crpix1']
crpix2 = header['crpix2']
cdelt1 = header['cdelt1']
cdelt2 = header['cdelt2']
naxis1 = header['naxis1']
naxis2 = header['naxis2']
#res2 = 0.00633

dist = 3.2 #in kpc
if args.dist:
  dist = args.dist

#deg1 *= abs(res2/res)
#deg2 *= res2/res

#Grad,Grad_Err,PosAng,PAErr,ChiSq = vfit_j(deg1,deg2,vel-mean(vel),vel_err,3.2,abs(res))
#Grad2,Grad_Err2,PosAng2,PAErr2,ChiSq2 = vfit_j2(deg1,deg2,vel-mean(vel),vel_err,3.2,abs(res),width=6.)
Grad,Grad_Err,PosAng,PAErr,ChiSq,GRx,GRy,errx,erry = vfit_j2(deg1,deg2,vel-mean(vel),vel_err,dist,abs(res),width=6.)

Grad_A,Grad_err_A,PosAng_A,PAerr_A,Chisq_A,Vc_A,Vc_err_A = VFit_OverAll( deg1, deg2, vel, vel_err, dist,abs(res), boot=0)

boolA = isnan(Grad)
Grad[boolA] = 0.001
Grad_Err[boolA] = 5.0
PosAng[boolA] = 0.0
#Grad = Grad[boolA]
#Grad_Err = Grad_Err[boolA]
#PosAng = PosAng[boolA]
#PAErr = PAErr[boolA]
#ChiSq = ChiSq[boolA]
#GRx = GRx[boolA]
#GRy = GRy[boolA]
#errx = errx[boolA]
#erry = erry[boolA]
#deg1 = deg1[boolA]
#deg2 = deg2[boolA]

gx = Grad*cos(PosAng*pi/180)
gy = Grad*sin(PosAng*pi/180)

#Grad_gx,Grad_Err_gx,PosAng_gx,PAErr_gx,ChiSq_gx,GRx_gx,GRy_gx,errx_gx,erry_gx = vfit_j2(deg1,deg2,abs(GRx-mean(GRx)),errx,3.2,abs(res),width=6.)

#Grad_gy,Grad_Err_gy,PosAng_gy,PAErr_gy,ChiSq_gy,GRx_gy,GRy_gy,errx_gy,erry_gy = vfit_j2(deg1,deg2,abs(GRy-mean(GRy)),erry,3.2,abs(res),width=6.)

Grad_gr,Grad_Err_gr,PosAng_gr,PAErr_gr,ChiSq_gr,GRx_gr,GRy_gr,errx_gr,erry_gr = \
 vfit_j2(deg1,deg2,Grad,Grad_Err,dist,abs(res),width=6.)
Grad_A_gr,Grad_err_A_gr,PosAng_A_gr,PAerr_A_gr,Chisq_A_gr,Vc_A_gr,Vc_err_A_gr = \
   VFit_OverAll( deg1, deg2, Grad, Grad_Err, dist,abs(res), boot=0)
#Grad_gr,Grad_Err_gr,PosAng_gr,PAErr_gr,ChiSq_gr,GRx_gr,GRy_gr,errx_gr,erry_gr = \
# vfit_j2(deg1[boolA],deg2[boolA],Grad[boolA],Grad_Err[boolA],dist,abs(res),width=6.)
#Grad_A_gr,Grad_err_A_gr,PosAng_A_gr,PAerr_A_gr,Chisq_A_gr,Vc_A_gr,Vc_err_A_gr = \
#   VFit_OverAll( deg1[boolA], deg2[boolA], Grad[boolA], Grad_Err[boolA], dist,abs(res), boot=0)


#fitsfile = pyfits.open('V0F.fits')
#plotgrad(f4file,Grad,PosAng)

def plotgrad(fifile,gradient,posangle,gradT,posangT,output,gunit=True):

  fitsfile = pyfits.open(fifile)
  fitsdata = fitsfile[0].data
  fitsdata[isnan(fitsdata)]=0
  vm = max(abs(fitsdata.min()),fitsdata.max())
  #print vm
  #f = apl.FITSFigure('V0F.fits')
  if args.colorm == 'RdBu':
    cmaps = cm.RdBu
  elif args.colorm =='hot':
    cmaps = cm.hot
  else:
    cmaps = cm.gray

  if args.reset:
    fitsdata[fitsdata==0] = fitsdata.max()*100
    os.system("rm faux.fits")
    pyfits.writeto('faux.fits',fitsdata,header)
    f = apl.FITSFigure('faux.fits')
  else:
    f = apl.FITSFigure(fifile)

  if args.dist:
    h1=naxis2*0.90
    delta = abs(cdelt1*args.dist*1000*pi/180)
    #print 'delt',delta
    plot((2.,2.+1./delta),(h1,h1),'k-',linewidth=2)
    dh = 0.2
    plot((2.,2.),(h1+dh,h1-dh),'k-',linewidth=2)
    plot((2.+1./delta,2.+1./delta),(h1+dh,h1-dh),'k-',linewidth=2)
    #print naxis1-2.-1./delta,naxis1-2.
    #plot((1,1+1./delta),(h1,h1),'k-',linewidth=3)
    text(2.,h1+0.5,'1 PC')
    #f.add_label(0.05,h1+0.05,"1 PC",relative=True)

  if args.items:
    import asciitable
    fielditems = asciitable.read(args.items, delimiter="\t",guess=False)
    FIraS = fielditems['ICRS_(J2000)_RA']
    FIdecS = fielditems['ICRS_(J2000)_DEC']
    l = len(FIraS)
    ra=zeros(l)
    dec=zeros(l)
    import coords as C
    for i,raA in enumerate(FIraS):
      ra1=float(raA[:2])
      ra2=float(raA[3:5])
      ra3=float(raA[6:])
      ra[i]=(ra1+ra2/60+ra3/3600)/24*360

    for i,decA in enumerate(FIdecS):
      dec1=float(decA[1:3])
      dec2=float(decA[4:6])
      if(decA[7:]):
        dec3=float(decA[7:])
      else:
        dec3=0.0
      if(decA[0]=='+'):
        dec[i]=dec1+dec2/60+dec3/3600
      else:
        dec[i]=(dec1+dec2/60+dec3/3600)*-1.0

    rapix,decpix=f.world2pixel(ra,dec)
#estrellas
    mssize=12
    maux= fielditems['Otype']=='*'
    plot(rapix[maux],decpix[maux],'y*',ms=mssize,label="Stars")
#UV source
    maux= fielditems['Otype']=='UV'
    plot(rapix[maux],decpix[maux],'m.',ms=mssize,label="UV source")
#mm radio source
    maux= fielditems['Otype']=='mm'
    plot(rapix[maux],decpix[maux],'c.',ms=mssize,label="mm radio source")
#mm radio source
    maux= fielditems['Otype']=='smm'
    plot(rapix[maux],decpix[maux],'cs',ms=8,label="smm radio source")
#ir source
    maux= fielditems['Otype']=='IR'
    plot(rapix[maux],decpix[maux],'r.',ms=mssize,label="infrared source")
#hii region
    maux= fielditems['Otype']=='HII'
    plot(rapix[maux],decpix[maux],'rs',ms=mssize,label="HII region")
#pulsar
    maux= fielditems['Otype']=='Psr'
    plot(rapix[maux],decpix[maux],'b^',ms=mssize,label="Pulsar")
#dark cloud nebula
    maux= fielditems['Otype']=='DNe'
    plot(rapix[maux],decpix[maux],'g.',ms=mssize,label="Dark CLoud Nebula")
#radio source
    maux= fielditems['Otype']=='Rad'
    plot(rapix[maux],decpix[maux],'c^',ms=mssize,label="Radio source")

    params0 = {'legend.fontsize': 6}
    pylab.rcParams.update (params0)
    legend(loc=0,numpoints=1,markerscale=0.5)#loc=(0.85,0.7)

  #f = apl.FITSFigure(f4file)
  if args.vedit:
    f.show_colorscale(pmin=1,stretch='linear', vmin=-vm, vmax=vm,cmap=cmaps)
  else:
    f.show_colorscale(pmin=1,stretch='linear',cmap=cmaps, vmax=vm)
  f.add_colorbar()

  arry = crval2+deg2-cdelt2*(crpix2)
  #print arry
  #arrx = crval1+deg1+0.00633*crpix1
  rescale = 1./cos(arry*pi/180)
  for i in range(1,len(rescale)):
    rescale[i] = mean(rescale[0:i])

  #print res/cos(crval2*pi/180)
  arrx = crval1+(deg1-cdelt1*(crpix1))*rescale
  #print arrx
#scale = max(abs(Grad))*res*2.
#scaley = Grad_A/abs(res)/1.2
#scalex = Grad_A/res2/1.2
  scalex = max(abs(gradient))/cdelt1/2.
  scaley = max(abs(gradient))/cdelt2/2.
#scalex = scaley
  arrdx = gradient*cos(posangle*pi/180)
#arrdx = arrdx/max(abs(Grad))*cdelt1*4.
  arrdx = arrdx/scalex
  arrdy = gradient*sin(posangle*pi/180)
#arrdy = arrdy/max(abs(Grad))*cdelt2*4.
  arrdy = arrdy/scaley

  arrow_len = sqrt(arrdx**2+arrdy**2)*.2/abs(cdelt1)
  #print len(arrow_len)
  f.show_arrows(arrx,arry,arrdx,arrdy)
#for i,l in enumerate(arrow_len):
#  f.show_arrows(arrx[i],arry[i],arrdx[i],arrdy[i],head_width=l,head_length=l)
#  if i%50 ==0:
#    print i

  arr_l = sqrt( (gradT*cos(posangT*pi/180)/scalex)**2 + (gradT*sin(posangT*pi/180)/scaley)**2 )
  f.show_arrows(f.pixel2world(naxis1*0.9,naxis2*0.1)[0],f.pixel2world(naxis1*0.9,naxis2*0.1)[1], \
    gradT*cos(posangT*pi/180)/scalex,gradT*sin(posangT*pi/180)/scaley,\
    head_width=.2*arr_l/abs(cdelt1),head_length=.2*arr_l/abs(cdelt1))

  if args.cont:
    contfile = pyfits.open(f3file)
    contdata = contfile[0].data
    contdata[isnan(contdata)]=0
    vmax=contdata.max()
    cont1=arange(args.cmin+0.1,args.cmax,args.cdelt)*vmax
    if args.contdef:
      f.show_contour(f3file,colors='yellow',layer='cont1')
    else:
      f.show_contour(f3file,colors='yellow',layer='cont1',levels=cont1)

  if args.cont2:
    contfile = pyfits.open(args.cont2)
    contdata = contfile[0].data
    contdata[isnan(contdata)]=0
    vmax=contdata.max()
    cont2=arange(args.cmin+0.1,args.cmax,args.cdelt)*vmax
    f.show_contour(args.cont2,colors='blue',layer='cont2')#,levels=cont2)
  alt = 0.12
  if (gradT*sin(posangT*pi/180)/scaley ) > 0:
    alt = 0.06
  if gunit:
    f.add_label(0.9,alt,"%.2f km/s/pc"%(gradT),relative=True)
  else:
    f.add_label(0.9,alt,"%.2f km/s/pc$^2$"%(gradT),relative=True)
  f.add_beam()
  f.beam.show()
  f.beam.set_facecolor('white')
  f.beam.set_edgecolor('black')
  f.frame.set_color('Black')
  f.ticks.set_color('Black')
  f.save(output)
  f.close()

print 'Grad_A',Grad_A
print 'PosAng',PosAng_A


Grad2x = Grad*cos(PosAng*pi/180)-Grad_A*cos(PosAng_A*pi/180)
Grad2y = Grad*sin(PosAng*pi/180)-Grad_A*sin(PosAng_A*pi/180)
Grad2 = sqrt(Grad2x**2+Grad2y**2)
PosAng2 = arctan2(Grad2y,Grad2x)*180/pi

arry = crval2+deg2-cdelt2*(crpix2)
rescale = 1./cos(arry*pi/180)
for i in range(1,len(rescale)):
  rescale[i] = mean(rescale[0:i])
arrx = crval1+(deg1-cdelt1*(crpix1))*rescale

with file('grad.txt', 'w') as outfile:
    outfile.write("RA DEC Grad PA\n")
    np.savetxt(outfile,transpose([arrx,arry,Grad,PosAng]),fmt=('%.5f','%.5f','%.3f', '%.3f' ) )

plotgrad(f0file,Grad,PosAng,Grad_A,PosAng_A,args.dir+'/'+args.out)
#plotgrad(f4file,Grad_gx,PosAng_gx,(args.out).replace(".png","_gx.png"))
#plotgrad(f4file,Grad_gy,PosAng_gy,(args.out).replace(".png","_gy.png"))
plotgrad(f0file,Grad_gr,PosAng_gr,Grad_A_gr,PosAng_A_gr,(args.dir+'/'+args.out).replace(".png","_gr.png"),gunit=False)

plotgrad(f0file,Grad2,PosAng2,Grad_A,PosAng_A,(args.dir+'/'+args.out).replace(".png","_1m.png"))

os.system("rm faux*.fits")

PosAng2 = copy(PosAng)
PosAng2[PosAng2<0.] += 360
plot(Grad,PosAng2,'.')
xlabel('$|$ Grad $|$')
ylabel('PA')
savefig(args.dir+'/gradvspa.png')
close()

plot(Grad,sigma,'.')
xlabel('$|$ Grad $|$')
ylabel('$\sigma$')
savefig(args.dir+'/gradvssigma.png')
close()

plot(Grad,intf,'.')
xlabel('$|$ Grad $|$')
ylabel('Integrated Intensity')
savefig(args.dir+'/gradvsint.png')
close()

plot(Grad,dvel,'.')
xlabel('$|$ Grad $|$')
ylabel('Shifted velocity')
savefig(args.dir+'/gradvsdvel.png')
close()

from histplot import histOutline
def histo(bins,data,out,xlab="",ylab=""):
  (binT, dataT) = histOutline(data,bins)
  plot(binT,dataT,'k-')
  xlabel(xlab)
  ylabel(ylab)
  if axis()[3] == dataT.max():
    ax1 = array(axis())*[1,1,1,1.1]
    axis(ax1)
  savefig(out)
  close()


histo(arange(25)*15-180,PosAng,'h_pa.png','PA','counts')

histo(arange(0.0,1.51,0.1),Grad,'h_grad.png','Gradient','counts')

bins = [0.0,0.25,0.5,0.75,1.,1.25,1.5,1.75,2,2.25,2.5]
histo(bins,sigma,'h_sigma.png','$\sigma$','counts')

histo(range(9),intf,'h_intf.png','I Intensity','counts')

#fil1 = open("values.txt", 'w')
#fil1.write("dx dy vcen verr\n")
#for i in range(len(deg1)):
#  fil1.write("%.3e %.3e %.3e %.3e\n"%(deg1[i],deg2[i],vel[i],vel_err[i]))

