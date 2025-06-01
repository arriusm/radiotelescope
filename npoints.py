#!/Users/arri/miniconda3/bin/python

VERSION='V250531'

# import numpy as np
# from math import *
# from numpy import *
# from pylab import *
# from random import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import matplotlib.mlab as mlab
from scipy import optimize as opt
import scipy.ndimage as ndimage
# from mpl_toolkits.mplot3d import axes3d
import mpl_toolkits.mplot3d.axes3d as ax3d


def read(table) :
  print(table)
  data = np.genfromtxt(table,usecols=(5, 7, 11, 13, 19),comments='#')
  AZ    = data[:,0]
  EL    = data[:,1]
  AZoff = data[:,2]
  ELoff = data[:,3]
  P     = data[:,4]
  return (AZ,EL,AZoff,ELoff,P)


def click(event):
  global ix,iy
  ix = event.xdata
  iy = event.ydata
  print(' click  %7.4f %7.4f' % (ix,iy))
  return


def gauss(r,A,sig,C,r0) :
  return C+A*exp(-(r-r0)**2/(2.*sig**2))


#def gauss2d(x,y,A,sig,C,x0,y0) :
#  r2=(x-x0)**2+(y-y0)**2
#  return C+A*np.exp(-r2/(2.*sig**2))

def f_gauss2d(C, A, x0, y0, sigx, sigy) :
    return lambda x,y: C+A*np.exp(-((x-x0)**2/(2.*sigx**2)+(y-y0)**2/(2.*sigy**2)))

def fit_gauss2d(x,y,z):
    C = np.min(z)/2.
    C = 0.
    A = np.max(z)-C
    #print('A=',A)
    #print('C=',C)
    (ix,iy) = np.unravel_index(np.argmax(z),z.shape)
    print('image = ')
    print(z)
    print('index_max =',ix,iy)
    xmax = x[ix,iy]
    ymax = y[ix,iy]
    print('xmax ymax =',xmax,ymax)
    sigx = 2.43
    sigy = 2.43
    # p_in = (C,A,x0,y0,sigx,sigy)
    p_in = [C,A,xmax,ymax,sigx,sigy]
    # print('p_in = ',p_in)
    xy = np.stack((x,y))
    # print(xy)
    errfunc = lambda p: np.ravel( f_gauss2d(*p)(*xy) - z )
    p_fit, success = opt.leastsq(errfunc, p_in)
    #print('p_fit = ',p_fit)
    return p_fit

################################################################################


import argparse

parser = argparse.ArgumentParser(description='npoints.py by Arno Riffeser '+VERSION)
parser.add_argument('-v', '--verbose', type=int,  dest='verbose',  default='0', help='verbose level')
parser.add_argument('FILES',           nargs='+', help='files')
args = parser.parse_args()


# table = ['171114_141108_Sun_npoint_.tab']
table=args.FILES
verbose=args.verbose


global ix, iy

outtable = open('npoints.tab','w')

for tab in table :

  (AZ,EL,AZoff,ELoff,P) = read(tab)
  
  EL2d = np.zeros((5,5))
  AZ2d = np.zeros((5,5))
  P2d  = np.zeros((5,5))
  
  z=0
  for i in range(0,5) :
    for j in range(0,5) :
      AZ2d[i,j] = AZoff[z]*np.cos((EL[z])/180.*np.pi)
      EL2d[i,j] = ELoff[z]
      P2d[i,j]  = P[z]
      z+=1

  # P2d[0,0] = P2d[1,0] ## Trick
      
  ################################################################################
  # FIT
  
  # # random walk
  # # 1.63427616205 1111602.536 2.43812476416 286460.783732 1.56207115957 0.00875307509684
  # A        = 150.
  # sig_meas = 2.43812476416
  # C        = 0.
  # x0       = 0.
  # y0       = 0.
  # AT   = A
  # CT   = C
  # sigT = sig_meas
  # x0T  = x0
  # y0T  = y0
  # chi2best=1e10
  # ranfak=0.01
  # ranmin=1.0-ranfak/2.
  # ranlen=1
  # for k in range(0,ranlen) :
  #   chi2=0.
  #   for j in range(0,5) :
  #     for i in range(0,5) :
  #       GG  = gauss2d(AZ2d[i,j],EL2d[i,j],AT,sigT,CT,x0T,y0T)
  #       #EE2 = 4e4 * P2d[i,j]
  #       #chi2 += (GG-P2d[i,j])**2/EE2
  #       chi2 += (GG-P2d[i,j])**2/1e10
  #   if chi2<chi2best :
  #     print(k,chi2,AT,sigT,CT,x0T,y0T)
  #     A=AT
  #     C=CT
  #     sig_meas=sigT
  #     x0=x0T
  #     y0=y0T
  #     chi2best=chi2
  #   AT   =       A*(ranmin+random()*ranfak)
  #   CT   =       C*(ranmin+random()*ranfak)
  #   sigT = sig_meas*(ranmin+random()*ranfak)
  #   x0T  =      x0+random()*ranfak
  #   y0T  =      y0+random()*ranfak
  
  
  #print(AZ2d)
  #print(EL2d)
  p_fit = fit_gauss2d(AZ2d,EL2d,P2d)
  C = p_fit[0]
  A = p_fit[1]
  x0 = p_fit[2]
  y0 = p_fit[3]
  sig_meas = np.sqrt(p_fit[4]*p_fit[5])
  
  
  def gauss2d(x,y,A,sig,C,x0,y0) :
      p = [C,A,x0,y0,sig,sig]
      xy = np.stack((x,y))
      return f_gauss2d(*p)(*xy)
  
  # print(gauss2d(11,11,50.,3.,100,10.,10.))
  # print(gauss2d([11,12],[11,12],50.,3.,100,10.,10.))
  # print(gauss2d([[11,12],[13,14]],[[11,12],[12,15]],50.,3.,100,10.,10.))
  
  
  ################################################################################
  
  # Radius Sonne 1.4e6/2./150e6/pi*180. deg
  sig_obj=0.27 # Sonne
  sig_beam=np.sqrt(sig_meas**2-sig_obj**2) # Eq. 8.15 in Tools_of_Radio_Astronomy.pdf
  
  # exp(-(HPBW/2.)**2/(2.*sig**2) = 1/2
  # HPBW**2 = 8.*sig**2 ln 2
  # HPBW = sig * sqrt(8.* ln 2)
  HPBW_beam = sig_beam*np.sqrt(8.*np.log(2.))
  #OmegaHPBW = 2. * pi * (1.-cos(HPBW_beam/180.*pi / 2.) )
  OmegaHPBW = 1.13309 * (HPBW_beam/180.*np.pi)**2  # Eq. 8.13 in Tools_of_Radio_Astronomy.pdf
  
  
  # OmegaA = Int_4pi Pn(theta,phi) dOmega
  # Integral (auf ebener Flaeche nicht auf Kugel) ueber 2d Gauss mit A=1 ist: 
  IntPn = 2.*np.pi*sig_beam**2  # [deg^2]
  
  # Umrechnung von deg^2 in sr
  # Kugeloberflaeche in [deg^2]: 4pir^2 = 4pi(U/2pi)^2 = 4pi(360./2pi)^2 = 360.^2/pi = 41253 deg^2 =!= 4pi
  
  # OmegaA = IntPn / 41253. * 4*pi
  # OmegaA = IntPn * 4.*pi / (360.**2/np.pi)
  OmegaA = IntPn * 4.*np.pi**2 / 360.**2
  
  lam = 0.21 # m
  Aeff = lam**2/OmegaA
  Deff = 2.*np.sqrt(Aeff/np.pi)
  Dgeom = 2.42
  Ageom = (Dgeom/2.)**2*np.pi

  print('  sig_meas   [deg] = ',sig_meas)
  print('  sig_obj    [deg] = ',sig_obj)
  print('  sig_beam   [deg] = ',sig_beam)
  
  if verbose > 0 :
    print('HPBW_beam  [deg] = ',HPBW_beam)
    print('Omega_HPBW  [sr] = ',OmegaHPBW)
    print('Omega_A     [sr] = ',OmegaA)
    print('A_eff      [m^2] = ',Aeff)
    print('D_eff        [m] = ',Deff)
    print('A_geom     [m^2] = ',Ageom)
    print('D_geom       [m] = ',Dgeom)
  
  ################################################################################
  
  # G2d  = np.zeros((5,5))
  # for j in range(0,5) :
  #   for i in range(0,5) :
  #     G2d[i,j]  = gauss2d(AZ2d[i,j],EL2d[i,j],A,sig_meas,C,x0,y0)
  # 
  # G2d_AZfix = np.zeros((5,101))
  # G2d_ELfix = np.zeros((101,5))
  # xAZ  = np.zeros(101)
  # xEL  = np.zeros(101)
  # for i in range(0,101) :
  #   for j in range(0,5) :
  #     xEL[i] = i/100.*16.-8.
  #     xAZ[i] = i/100.*16.-8.
  #     G2d_AZfix[j,i]  = gauss2d(AZ2d[2,j],xEL[i],A,sig_meas,C,x0,y0)
  #     G2d_ELfix[i,j]  = gauss2d(xAZ[i],EL2d[j,2],A,sig_meas,C,x0,y0)
  
  lenfull = 21
  Gfull  = np.zeros((lenfull,lenfull))
  xfull  = np.zeros((lenfull,lenfull))
  yfull  = np.zeros((lenfull,lenfull))
  for j in range(0,lenfull) :
    for i in range(0,lenfull) :
      xfull[i,j] = i/(lenfull-1.)*16.-8.
      yfull[i,j] = j/(lenfull-1.)*16.-8.
      Gfull[i,j] = gauss2d(xfull[i,j],yfull[i,j],A,sig_meas,C,x0,y0)
  
  ################################################################################
  
#   fig1=plt.figure(figsize=(7,7),facecolor='w')
#   plt.xlabel('AZimut * cos(ELe) [deg]')
#   plt.ylabel('Power 1')
#   for i in range(0,5) :
#     x=AZ2d[i,0:5]
#     y=P2d[i,0:5]
#     plt.scatter(x,y,zorder=2,color='black',marker='+',s=200)
#     plt.plot(x,y,zorder=1,color='blue',lw=1)
#     #g=G2d[i,0:5]
#     #plt.scatter(x,g,zorder=2,color='red',marker='x',s=200)
#     plt.plot(xAZ,G2d_ELfix[0:101,i],zorder=0,color='red',lw=1)
#   plt.axis([-8.,8.,C*0.75,(A+C)*1.1])
#   #plt.title('Sonne 01.07.2015 11:14-11:20')
#   #plt.savefig('npoints_AZ.pdf')  
#   plt.savefig('npoints_AZ_'+tab+'.pdf')  
#   fig1.canvas.mpl_connect('button_press_event',click)
#   
#      
#   fig2=plt.figure(figsize=(7,7),facecolor='w')
#   plt.xlabel('ELevation [deg]')
#   plt.ylabel('Power 1')
#   for i in range(0,5) :
#     x=EL2d[0:5,i]
#     y=P2d[0:5,i]
#     plt.scatter(x,y,zorder=2,color='black',marker='+',s=200)
#     plt.plot(x,y,zorder=1,color='green',lw=1)
#     #g=G2d[0:5,i]
#     #plt.scatter(x,g,zorder=2,color='red',marker='x',s=200)
#     plt.plot(xEL,G2d_AZfix[i,0:101],zorder=0,color='red',lw=1)
#   plt.axis([-8.,8.,C*0.75,(A+C)*1.1])
#   #plt.title('Sonne 01.07.2015 11:14-11:20')
#   #plt.savefig('npoints_EL.pdf')  
#   plt.savefig('npoints_EL_'+tab+'.pdf')  
#   fig2.canvas.mpl_connect('button_press_event',click)
#   
#   
#   fig3=plt.figure(figsize=(7,7),facecolor='w')
#   ax3=fig3.add_subplot(111, projection='3d')
#   ax3.plot_wireframe(xfull,yfull,Gfull,zorder=0,rstride=1,cstride=1,color='red')
#   #ax3.plot_wireframe(AZ2d,EL2d,G2d,zorder=1,rstride=1,cstride=1,color='blue',linewidth=4)
#   ax3.plot_wireframe(AZ2d,EL2d,P2d,zorder=1,rstride=1,cstride=1,color='blue',linewidth=4)
#   ax3.set_xlabel('AZimut [deg]')
#   ax3.set_ylabel('ELevation [deg]')
#   ax3.set_zlabel('Power')
#   #title('Sonne 01.07.2015 11:14-11:20')
#   #plt.savefig('npoints_3d.pdf')  
#   plt.savefig('npoints_3d_'+tab+'.pdf')  
  
  
  fig4=plt.figure(figsize=(7,7),facecolor='w')
  norm=Gfull.max()
  norm=P2d.max()
  xdata = ndimage.zoom(AZ2d,5)
  ydata = ndimage.zoom(EL2d,5)
  # data  = ndimage.zoom(G2d/norm,5)
  data  = ndimage.zoom(P2d/norm,5)
  # http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contour
  # http://matplotlib.org/examples/color/colormaps_reference.html
  flevels = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
  flevels = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1]
  plt.contourf(xdata,ydata,data,flevels,cmap=plt.cm.Greys,origin='lower',extend='both')
  con1=plt.contour(xdata,ydata,data,flevels,colors='green')
  plt.clabel(con1,inline=1,fontsize=10)
  levels = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
  con2=plt.contour(xfull,yfull,Gfull/norm,levels,colors='red')
  plt.clabel(con2,inline=1,fontsize=10)
  plt.axis([-7.,7.,-7.,7.])
  #plt.title('Sonne 01.07.2015 11:14-11:20')
  plt.xlabel('Azimut [deg]')
  plt.ylabel('Elevation [deg]')
  # plt.savefig('npoints_contour.pdf')  
  plt.savefig('npoints_contour_'+tab+'.pdf')  
  ix=x0
  iy=y0
  bb=fig4.canvas.mpl_connect('button_press_event',click)
  plt.show()
  #print('center = %5.4f  %5.4f' % (ix,iy))
  AZEL = '{:7.2f}'.format(AZ[12])+' '+'{:7.2f}'.format(float(EL[12]))
  x0y0 = '{:7.4f}'.format(ix)+'  '+'{:7.4f}'.format(iy)
  print('x0 y0 =',x0y0,'       AZ EL =',AZEL)
  s = tab + '     ' + x0y0 + '    ' + AZEL + '\n'
  outtable.write(s)

outtable.close()
