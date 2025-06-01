#!/Users/arri/miniconda3/bin/python

VERSION='V250601'

import sys, os, math, random, time

#import matplotlib
#matplotlib.use(u'MacOSX')
#matplotlib.use(u'TkAgg')
#matplotlib.use(u'pdf')

import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
#from pylab import *
from math import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, ScalarFormatter

import astropy.units as u
from astropy.coordinates import SkyCoord,FK5
import argparse

parser = argparse.ArgumentParser(description='radio_spectra.py by Arno Riffeser '+VERSION)

parser.add_argument("-d", "--dir",                             dest="dir",       default="", help="directory, default is current directory")
parser.add_argument("-m", "--mask",      action='store_true',  dest="domask",                help="masking")
parser.add_argument("-n", "--no-display",action='store_false', dest="display",               help="turn off showing on display")
parser.add_argument("-f", "--makeflat",  action='store_true',  dest="makeflat",              help="applies a flat division")
parser.add_argument("-o", "--plotflat",  action='store_true',  dest="plotflat",              help="plots the flat")
parser.add_argument("-p", "--pdf",       action='store_true',  dest="pdf",                   help="stores a pdf")
parser.add_argument("-g", "--png",       action='store_true',  dest="png",                   help="stores a png")
parser.add_argument("-T", "--plotTemp",  action='store_true',  dest="plotTemp",              help="plot Temp on y axis")
parser.add_argument("-x", "--xplot",                           dest="xplot",    default="r", help="xaxis: r - register, f - frequency, v - velocity")
parser.add_argument("-v", "--verbose",   type=int,             dest="verbose",  default="2", help="verbose level")
parser.add_argument('FILES',             nargs='+',                                          help="files")

args = parser.parse_args()

### dir = "/Users/arri/ASTRO/Radioteleskop/DATA/"
### dir = "/home/radioadmin/Work/srtnver3/data/2016/"
### dir = "/home/radioadmin/Work/srtnver6/data_170721"
### dir = "/Users/arri/ASTRO/Radioteleskop/DATA/171024"
### dir = "/home/radioadmin/Work/srtnver6/DATA"
### dir = "/Users/arri/ASTRO/Radioteleskop/DATA/171024"
dir =      args.dir
display =  args.display
makeflat = args.makeflat
plotflat = args.plotflat
pdf =      args.pdf
png =      args.png
plotTemp = args.plotTemp
xplot =    args.xplot
verbose =  args.verbose

############## mask

#####  mask version 20250601


newpeak = []

#startpeak=list(range(0,87+1))
startpeak=list(range(0,156+1))

brightpeak=list(range(1629,1647+1))

triplepeak=[ # immer 3 Linien im Abstand 13
    104,105,106,107,108,109,110,111,112,113,114,115,  122,  131,132,133,134,135,136,137,138,139,140,141,
    158,159,160,161,162,163,164,165,166,167,  185,186,187,188,189,190,191,192,193,194,195,196,
    231,232,233,  244,245,246,  257,258,259,
    316,317,318,  329,330,331,  343,344,
    402,403,  415,416,  428,429,
    461,462,  486,487,488,489,  499,500,501,502,  513,514,515,
    572,573,574,  585,586,587,  598,599,600,
    657,658,659,660,  669,670,671,672,673,675,  682,683,684,685,
    742,743,744,745,  755,756,757,758,759,  768,769,770,771,
    801,802,803,804,805,  827,828,829,830,  840,841,842,843,
    853,854,855,856,  880,881,882,  886,887,888,
    907,908,909,910,911,912,913,914,  923,924,925,926,927,928,  938,939,940,941,
    995,996,997,998,999,1000,  1011,1012,1013,  1024,1025,1026,
    1081,1082,1083,1084,1085,  1096,1097,1098,  1109,1110,1111,
    1167,1168,1169,1170,  1181,1182,1183,  1194,1195,1196,
    1236,1249,  1266,1267,1268,  1283,1284,  
    1338,1339,1340,  1351,1352,1353,  1364,1365,1366,
    1423,1424,1425,  1436,1437,1438,  1450,1451,
    1482,  1495,  1509,1510,1511,1512,
    1522,1523,1524,  1535,1536,1537,
    1593,1594,1595,1596,1597,  1606,1607,1608,1609,1610,  1619,1620,1621,1622,1623,
    1646,1647,1648,1649,  1674,  1679,1680,1681,
    1692,1693,1694,1696,1697,  1705,1706,1707,
    1763,1764,1765,1766,  1776,1777,1778,1779,1780,1781,1782, 1790,1791,1792,1793,1797, 1810,
    1849,1850,1851,  1862,1863,1864,  1875,1876,1877,
    1934,1935,1936,  1947,1948,1949,  1961,1962,1963,
    2019,2020,2021,2022,  2032,2033,2034,2035,  2046,2047,2048,
    2105,2106,2107,  2118,2119,2120,  2129,2130,2131,2132,2133,
    2190,2191,2192,  2203,2204,2205,  2216,2217,2218,
    2288,2289,  2290,2291,  2301,2302,2303,2304,  2316,
    2362,  2375, 2388,
    2444,2445,2446,2447,  2457,2458,2459,2460,  2470,2471,2472,2473,2474,
    2453,  2466,  2479,2480,2481,2482,2483,2484,2485,   2492,2493,   2498,2499,
    2530,2531,2532,  2544,2545,2546,  2557,2558,2559,
    2616,2617,2618,  2629,2630,2631,  2642,2643,2644,
    2701,2702,2703,  2714,2715,2716,2717,  2727,2728,2729,
    2786,2787,2788,  2799,2800,2801,  2812,2813,2814,
    2871,2872,2873,  2884,2885,2886,  2898,2899,
    2956,2957,2958,  2969,2970,2971,  2983,2984,
    3042,3043,  3055,3056,3057,  3068,3069,3070,  3083,3085, 
    3127,3128,3129,  3140,3141,3142,  3153,3154,3155,
    3193,  3205,3206,  3219,
    3212,3213,3214,  3225,3226,3227,3228,  3238,3239,3240,
    3297,3298,3299,3300,  3310,3311,3312,3313,  3323,3324,3325,
    3465,3466,3467,3468,3469,3470,  3478,3479,3480,3481,3482,3483,  3490,3491,3492,3493,3494,3495,3496,
    3554,  3565,3566,3567,  3579,3580,
    3638,3639,  3651,3652,3653,  3664,3665,3666,
    3721,3722,3723,3724,3725,  3736,3737,3738,  3749,3750,3751,
    3808,3809,3810,  3822,3823,  3835,3836,
    3683,3684,  3895,  3908,  3921,3922,3923,3924,
    3978,3979,3980,  3992,3993,  4006,
    4015,4016,4017,  4028,4029,4030,4031,  4041,4042,4043,
    4064,4065,  4077,  4090] 

doublepeak=[ # oft 2er Paare im Abstand 7
    88,95,  539,540,546,547,  563,564,565,566,  624,625,626,631,632,  709,710,711,  
    957,  965,966,972,973,  1050,1051,1052,1057,1058,1059,  1128,1136,1143,1145,
    1469,1477,  1562,1568,1569,1574,  1817,1824,  1884,1892,  1987,1988,1994,1995,2003,2004,2005, 
    2152,2157,2158, 2419,2420,2421,2426,  2963,2964,2965,  3266,3272,3273,
    3441,3442,3443,3444,3452 ]

singlepeak=[179, 444, 497, 615,616,  650,651,  736,737, 796,  872,
            957,  1163,1164, 1223,  1253,  1327, 1380, 1403,1404,1405,
            1543,  1671, 1769,  1713,  1738,1739,  1748,  1818,  1834,  1843,1844, 1855, 
            1902,1903,  1920, 1928,1929, 1941,1942, 1953, 1990, 
            2002,  2006,  2027,  2039,  
            2072,2073,2074,2075,2076,  2102, 2123,2124,2125, 2185,2186,  2213,  2224,
            2249,2250,  2275,2276,2277,  2329,    2591, 2608,  2663,2664,
            2761,  2840,  2843,2844,2845,  3084,  3187,  3246,  3351,
            3521,  3547, 3563,3564,  3685,3686,3687,3688,  3804,  4078,  4091] 

mask = startpeak + brightpeak + triplepeak + doublepeak + singlepeak + newpeak

#for i in mask :
#    print(i)


##### mouse click

def onclick(event):
  #print('  %9.4f MHz   %8.2f K   %10.0f ' % (event.xdata, event.ydata, event.ydata/factor))
  print('  %9.4f    %10.0f ' % (event.xdata, event.ydata))
  #ax.plot([event.xdata,event.xdata], [0,1e6], linestyle='-', color='red',  linewidth=1)
  #show()
  return

def onclick_MHz(event):
  #print('  %9.4f MHz   %8.2f K   %10.0f ' % (event.xdata, event.ydata, event.ydata/factor))
  print('  %9.4f MHz   %8.2f K   %10.0f ' % (event.xdata, event.ydata, event.ydata))
  #ax.plot([event.xdata,event.xdata], [0,1e6], linestyle='-', color='red',  linewidth=1)
  #show()
  return

def onclick_kms(event):
  #print('  %9.4f km/s  %8.2f K   %10.0f ' % (event.xdata, event.ydata, event.ydata/factor))
  print('  %9.4f km/s  %8.2f K   %10.0f ' % (event.xdata, event.ydata, event.ydata))
  #ax.plot([event.xdata,event.xdata], [0,1e6], linestyle='-', color='red',  linewidth=1)
  #show()
  return


# void GalactictoRadec(double glat, double glon, double *ra, double *dec)
#   /* galactic to radec  2000 epoch pole at 12h51.4 27.1 */
# {
#     double a, xg, yg, zg, xr, yr, zr, d0, dp, r0, rp;
#     d0 = -(28.0 + 56.0 / 60.0) * PI / 180.0;
#     r0 = (17.0 + 45.5 / 60.0) * PI / 12.0;
#     dp = 27.1 * PI / 180.0;
#     rp = (12.0 + 51.4 / 60.0) * PI / 12.0;
#     zr = sin(d0);
#     xr = cos(r0 - rp) * cos(d0);
#     yr = sin(r0 - rp) * cos(d0);
#     xg = xr * sin(dp) - zr * cos(dp);
#     yg = yr;
#     a = atan2(yg, xg);
#     xg = cos((glon * PI / 180.0) + a) * cos(glat * PI / 180.0);
#     yg = sin((glon * PI / 180.0) + a) * cos(glat * PI / 180.0);
#     zg = sin(glat * PI / 180.0);
#     xr = xg * sin(dp) + zg * cos(dp);
#     yr = yg;
#     zr = zg * sin(dp) - xg * cos(dp);
#     *dec = atan2(zr, sqrt(xr * xr + yr * yr));
#     *ra = atan2(yr, xr) + rp;
#     if (*ra < 0)
#         *ra += TWOPI;
# }



# calculate velocity of the local standard of rest
def calc_vlsr(time, ra_h, de_d) :
  #ra = ra_h*pi/12.
  ra = ra_h*15.*pi/180.
  de = de_d*pi/180.
  ############ velocity sun around LSR #############
  # http://coursewiki.astro.cornell.edu/Astro4410/RadialVelocities
  # sun 20km/s towards ra=18h          de=30d          (equinox B1900.0)
  #                    ra=18h03m50.24s de=30d00m16.8s  (equinox J2000.0)
  ra_v = 18. # h
  de_v = 30. # deg
  # ra_v,de_v = transEpoch( 18., 30., 1900., 2015., 1., 1.)
  x0 = 20. * cos(ra_v*pi/12.) * cos(de_v*pi/180.)
  y0 = 20. * sin(ra_v*pi/12.) * cos(de_v*pi/180.)
  z0 = 20. *                    sin(de_v*pi/180.) 
  # The algorithm for calculating VLSR was written by M. A. Gordon and 
  # can be found in Methods of Experimental Physics,
  #   Volume 12, Part C: Radio Observations , Ed. M. L. Meeks.
  # The IAU Standard Solar Motion assumes:
  #     V = 20 km/s toward
  #     RA = 18 h (Equinox B1900.0)
  #     DEC = 30 deg (Equinox B1900.0) 
  #     transformed Apex vector (Galactic): 
  #         l=   56.15740696    b=   22.76488475    V=   20.000000    
  # The Mihalas and Binney Solar Motion assumes:
  #     V = 16.5 km s-1 toward
  #     l = 53 deg
  #     b = 25 deg
  # consistent to:  http://ned.ipac.caltech.edu/forms/vel_correction.html
  vsun = -x0*cos(ra)*cos(de) - y0*sin(ra)*cos(de) - z0*sin(de)
  #print("vsun=",vsun)
  ############### velocity earth around the sun ############
  x0 = cos(ra) * cos(de)
  y0 = sin(ra) * cos(de)
  z0 =           sin(de)
  earth_axis = 23.5*pi/180.
  x = x0
  y = y0*cos(earth_axis) + z0*sin(earth_axis)
  z = z0*cos(earth_axis) - y0*sin(earth_axis)
  soulat  = atan2(z, sqrt(x*x + y*y))
  soulong = atan2(y, x)
  # long=280 day 1 :
  sunlong = (time / (365.2422 * 86400.) * 360. + 280.) / 180. * np.pi
  #print("sunlong = ",sunlong)
  vearth = -30. * cos(soulat) * np.sin(sunlong - soulong)
  #print("sunlong = ",(sunlong - soulong) % np.pi)
  # if verbose>=1:
  #   print("  vsun          = ",vsun," km/s")
  #   print("  vearth        = ",vearth," km/s")
  #   print("  vsun + vearth = ",(vsun + vearth)," km/s")
  return (vsun + vearth)

# time to calculate vlsr
# for i in 150*S8.tab; do echo $i `awk '(NR==2)' $i`; done
# 150415_1655_S8.tab    1429109873 
# 150424_1549_z0940_S8.tab    1429883362 
# 150508_1459_z0940_S8.tab    1431090003 
# 150518_1343_z0940_S8.tab    1431951303 

# Test for vlsr
# https://www.astro.virginia.edu/~emm8x/utils/vlsr.html
#   RA = 054721.30
#   DE = -014018
#                  Lat = 48.1453
#                  Long = 348.393
#                  Elev = 600
#           geo   ... topo   ... here
#   15.4.   41.52 ... 41.77  ... 41.84
#   24.4    39.29 ... 39.51  ... 39.68
#    8.5.   34.89 ... 35.05  ... 35.29
#   18.5.   31.17 ... 31.28  ... 31.52


# double vlsr(double time, double ra, double dec)
#   // calculate velocity of the local standard of rest
# {
#     double decc, rac, vsun, vearth, x0, y0, z0, sunlong, soulong, soulat, x, y,
#         z, dp, rp, gwest, grad, gpole, lon0, ggwest;
#     ############ velocity sun around LSR #############
#     x0 = 20.0 * cos(18.0 * PI / 12.0) * cos(30.0 * PI / 180.0);
#     y0 = 20.0 * sin(18.0 * PI / 12.0) * cos(30.0 * PI / 180.0);
#     z0 = 20.0 * sin(30.0 * PI / 180.0); /* sun 20km/s towards ra=18h dec=30.0 */
#     vsun = -x0 * cos(ra) * cos(dec) - y0 * sin(ra) * cos(dec)
#         - z0 * sin(dec);
#     ############### velocity earth around the sun ############
#     x0 = cos(ra) * cos(dec);
#     y0 = sin(ra) * cos(dec);
#     z0 = sin(dec);
#     x = x0;
#     y = y0 * cos(23.5 * PI / 180.0) + z0 * sin(23.5 * PI / 180.0);
#     z = z0 * cos(23.5 * PI / 180.0) - y0 * sin(23.5 * PI / 180.0);
#     soulat = atan2(z, sqrt(x * x + y * y));
#     soulong = atan2(y, x);
#     sunlong = (time * 360.0 / (365.2422 * 86400.0) + 280.0) * PI / 180.0; /* long=280 day 1 */
#     vearth = -30.0 * cos(soulat) * sin(sunlong - soulong);
#     ############### calculation glat glon ################
#     dp = 27.1 * PI / 180.0;
#     rp = (12.0 + 51.4 / 60.0) * PI / 12.0;
#     decc = -(28.0 + 56.0 / 60.0) * PI / 180.0;
#     rac = (17.0 + 45.5 / 60.0) * PI / 12.0;
#     gwest = cos(decc) * cos(rp - rac);
#     grad = cos(decc) * sin(rp - rac);
#     ggwest = gwest * sin(dp) - sin(decc) * cos(dp);
#     gpole = gwest * cos(dp) + sin(decc) * sin(dp);
#     lon0 = (atan2(ggwest, grad)) * 180.0 / PI;
#     gwest = cos(dec) * cos(rp - ra);
#     grad = cos(dec) * sin(rp - ra);
#     ggwest = gwest * sin(dp) - sin(dec) * cos(dp);
#     gpole = gwest * cos(dp) + sin(dec) * sin(dp);
#     d1.glat = (atan2(gpole, sqrt(ggwest * ggwest + grad * grad))) * 180.0 / PI;
#     d1.glon = (atan2(ggwest, grad)) * 180.0 / PI - lon0;
#     if (d1.glon < 0.0)
#         d1.glon += 360.0;
#     if (d1.glon > 360.0)
#         d1.glon += -360.0;
#     return vsun + vearth;
# }




#################################################################################################


##### constants
freqHI = 1420.406
regHI = 1804 
# T0     = 273. # K
T0     = 273.15 # K
c      = 299792.458 # km/s


##### frequency grid
# theoretical from ADC specs with 20 MHz clocking and 1416 MHz Mixing
fstart  = 1416.+1e-6/(8192./20e6)
fend    = 1416.+10. # 10 MHz = 20 MHz / 2.
istart  = 1
iend    = 4096
freqsep = (fend-fstart)/(iend-istart)
# print('freqsep = %15.10f kHz' % (freqsep*1000.))
# print('1001: %15.10f MHz' % (1416.+freqsep*1001.))
# print('2000: %15.10f MHz' % (1416.+freqsep*2000.))
# print('2001: %15.10f MHz' % (1416.+freqsep*2001.))
# print('3000: %15.10f MHz' % (1416.+freqsep*3000.))
freq = np.zeros(4095)
reg  = np.zeros(4095)
for i in range(0,4095) :
  freq[i] = fstart+(i+1-istart)*freqsep
  reg[i]  = i


if verbose>=3:
  print("register HI",reg[regHI],freq[regHI],freqHI)
  print("freq_min = %10.5f MHz" % (freq[0]))
  print("freq_max = %10.5f MHz" % (freq[4094]))

##### Temp calibfactor by Lisa Uhlirsch with 150701_1718_Absorber_T30.tab
if verbose>=10:
  Tvanecold  = T0+6. # Kuehlschrank
  Tvanehot   = T0+30. # Sommertag
  Pvanecold  = 399000.
  # Pvanehot angepasst, sodass factor = 0.000845 :
  Pvanehot = Pvanecold + (Tvanehot-Tvanecold)/0.000845
  # Pvanehot = 427400.
  # print("Pvanehot = ",Pvanehot)
  factor2015 = (Tvanehot-Tvanecold)/(Pvanehot-Pvanecold)
  print("calibfactor 2015 (Lisa Uhlirsch) =  %10.5g " % factor2015)


##### Temp calibfactor 161129 white
if verbose>=10:
  #                161129_1556_absorber_white_-1C.tab           583610           1962       363.9       91.9      272.0
  #        161129_1612_absorber_white_21C_contact.tab           618888           1957       385.9       91.9      294.0
  Tvanecold  = T0-1.  # WintertagKuehlschrank
  Tvanehot   = T0+21. # Zimmertemp
  Pvanecold  = 583610.  
  Pvanehot   = 618888.  
  factor2016 = (Tvanehot-Tvanecold)/(Pvanehot-Pvanecold)
  print("calibfactor 2016 (flat, astrolab) = %10.5g" % factor2016)

##### diode
#                161129_1658_Kalt1_151124_diode.tab           435645         277010       271.7      271.7        0.0
#                      161129_1646_Kalt1_151124.tab           293190         257849       182.8      182.8        0.0
if verbose>=10:
  Pcold  = 293190.  
  Pdiode = 435645.
  Tdiode = (435645.-293190.)*factor2016
  print("Tdiode = %10.3f" % Tdiode)


# creating measung.spec
meastab = open(os.path.join(dir,"measurements.tab"),'w')
meastab.write('{0:>50s}   {1:>14s} {2:>12s}   {3:>12s}   {4:>10s} {5:>10s} {6:>10s} {7:>10s}\n'.format("filename","ADUtot","sigADUtot","calibfactor","Ttot","Tobj","Tsys","Tobj"))


for filename in args.FILES:

  title=""
  pngfile=""

  #dfilename = os.path.join(dir,tables[j])
  #filename = tables[j]

  dfilename = os.path.join(dir,filename)

  if verbose>=1:
    print(filename)



  ##### preparing graphical window
  #fig = plt.figure(num=1, figsize=(10,6), dpi=100, facecolor='w', edgecolor='k')
  fig = plt.figure(num=1, figsize=(12,6), facecolor='w', edgecolor='k')
  ax = plt.subplot2grid( (1,1),(0,0) )

  ##### reading spectra 
  in_file  = open( str( dfilename ), "r" )

#     ADUtot         0.00000
#       Tobj         0.00000
#     factor     6.00000e-04
#       ADU0      -232.49410
#       ADU1      2629.72490
#    ADUtick         0.00000
#      Temp0        -0.13950
#      Temp1         1.57783
#   Temptick         0.00000
#      freq0      1416.00244
#      freq1      1426.00000
#   freqtick         1.00000
#         v0       -50.00000
#         v1        50.00000
#      vtick         5.00000
#   flatname flat_161129_1556_absorber_white_-1C.tab_50.spec

  line = in_file.readline();      time=float(line.split()[2])
  line = in_file.readline();        ra=float(line.split()[2])
  line = in_file.readline();        de=float(line.split()[2])
  line = in_file.readline(); ADUtot_in=float(line.split()[2])
  line = in_file.readline();      Tobj=float(line.split()[2])
  line = in_file.readline();    factor=float(line.split()[2])
  line = in_file.readline();      ADU0=float(line.split()[2])
  line = in_file.readline();      ADU1=float(line.split()[2])
  line = in_file.readline();   ADUtick=float(line.split()[2])
  line = in_file.readline();     Temp0=float(line.split()[2])
  line = in_file.readline();     Temp1=float(line.split()[2])
  line = in_file.readline();  Temptick=float(line.split()[2])
  line = in_file.readline();     freq0=float(line.split()[2])
  line = in_file.readline();     freq1=float(line.split()[2])
  line = in_file.readline();  freqtick=float(line.split()[2])
  line = in_file.readline();        v0=float(line.split()[2])
  line = in_file.readline();        v1=float(line.split()[2])
  line = in_file.readline();     vtick=float(line.split()[2])
  line = in_file.readline();  flatname=      line.split()[2] 

  if verbose>=1:
    print("  calibfactor = ",factor)
    
  F = np.zeros(4095)
  Fminmax = np.zeros((2,4095)) 
  Fminmax[0,]=1e20 
  Fminmax[1,]=-1e20
  #Fmin = np.zeros_like(freq)+1e20
  #Fmax = np.zeros_like(freq)-1e20
  #Fsingle = np.zeros_like(freq)
  nr = 0
  line=' '
  time_m=0.
  RA_m=0.
  DE_m=0.
  # while line != '' and nr<10 : # read only 10 lines
  while line != '':
    line = in_file.readline()
    #print('.')
    #print(len(line))
    #print(">",line[len(line)-5:len(line)],"<")
    if line!='' and line[:1]!='#' and line[:1]!=' ':
      nr+=1
      v = [s for s in line.split()]
      #v = [float(s) for s in line.split()]
      #for i in range(199,4294):
      # print(v[27],v[27+4095])
      time_m+=float(v[3])
      RA_m+=float(v[15])
      DE_m+=float(v[17])
      for i in range(0,4095):
        value=float(v[27+i])
        #Fsingle[i]=value
        F[i]+=value
        #if value<Fmin[i]:
        #  Fmin[i]=value
        #if value>Fmax[i]:
        #  Fmax[i]=value            
        if value<Fminmax[0,i]:
          Fminmax[0,i]=value
        if value>Fminmax[1,i]:
          Fminmax[1,i]=value            
      #if xplot=='v':
      #  ax.scatter(x, Fsingle, s=1, color='yellow', marker='o')
  #print(" ")

  if args.domask :
    F[mask]=np.nan
  
  if verbose>=1:
    print("  nr of averaged spectra = ",nr)
  # mean
  # for i in range(0,4095):
  #   if F[i]!=np.nan :
  F/=nr 
  time_m/=nr 
  RA_m/=nr 
  DE_m/=nr 

  print("  <time> = ",time_m)
  print("  <RA>   = ",RA_m)
  print("  <DE>   = ",DE_m)

  sc = SkyCoord(ra=RA_m*15., dec=DE_m, unit='deg', frame=FK5, equinox='J2000.0')
  gc = SkyCoord(l=10*u.degree, b=10*u.degree, frame='galactic',unit='deg')
  gc_new = sc.transform_to(gc)
  lon=gc_new.l.degree
  lat=gc_new.b.degree

  print("  galactic longitude <l>  = ",lon)
  print("  galactic latitude  <b>  = ",lat)

  vsun0 = 220.
  vstarRmin = vsun0
  vrmax = vstarRmin - vsun0*np.sin(lon/180.*np.pi)

  print("  vrmax = %10.2f km/s" % (vrmax))
  
  ################ calculate velocity grid ############
  # see also https://ned.ipac.caltech.edu/forms/calculator.html

  vlsr = 0.
  if ra!=0. and de!=0. :
    vlsr = calc_vlsr(time,ra,de)
    time_all = np.linspace(time-365.25*24.*3600., time, 1000, endpoint=True)
    vlsr_all = calc_vlsr(time_all,ra,de)
    vlsr_min = np.min(vlsr_all)
    vlsr_max = np.max(vlsr_all)
  else :
    vlsr = calc_vlsr(time_m,RA_m,DE_m)
    time_all = np.linspace(time_m-365.25*24.*3600., time_m, 1000, endpoint=True)
    vlsr_all = calc_vlsr(time_all,RA_m,DE_m)
    vlsr_min = np.min(vlsr_all)
    vlsr_max = np.max(vlsr_all)

  if verbose>=1:
    print('  vlsr     = %6.2f km/s' % (vlsr))
    print('  vlsr_min = %6.2f km/s' % (vlsr_min))
    print('  vlsr_max = %6.2f km/s' % (vlsr_max))

  veloc = (1.-freq/freqHI)*c - vlsr
  # print("  veloc = ",veloc)

  #print(' %s    vlsr= %5.2f km/s' % (filename,vlsr))
  freqHI_lsr = (1.-vlsr/c) * freqHI 
  freqHI_min = (1.-vlsr_min/c) * freqHI
  freqHI_max = (1.-vlsr_max/c) * freqHI 
  if verbose>=1:
    print('  freqHI     = %10.3f MHz' % (freqHI))
    print('  freqHI_lsr = %10.3f MHz' % (freqHI_lsr))
    print('  freqHI_min = %10.3f MHz' % (freqHI_min))
    print('  freqHI_max = %10.3f MHz' % (freqHI_max))
  
  ##### flatting spectra with flat frame
  norm = np.ones_like(F) # 1 for no flat correction

  if makeflat :
    if verbose>=1:
      print("  makeflat: ")
      print("      ",flatname)
    #flatname=os.path.join(dirflat,"flat.spec")
    #[mF]=read_columns( flatname, [2])
    # dataframe_flat = pd.read_csv(flatname,header=-1,delim_whitespace=True)
    # mF = np.array(dataframe_flat[1])
    flatdata = np.genfromtxt(flatname, names=['n','mF'])
    mF = flatdata['mF']
    # norm = mF/np.max(mF)
    # Warum [1704:1905]? Zentraler flacher Bereich um HI!! 
    #  a=freq[1704:1905]
    #  print(a)
    #  print("len=",len(freq[1704:1905]),a[0],a[100],a[200])
    norm = mF/np.nanmedian(mF[1704:1905])
    if verbose>=1:
      print("       flat median [1704:1905] = ",np.nanmedian(norm[1704:1905]))
      if verbose>=3:
        print("       flat median [1001:2001] = ",np.nanmedian(norm[1001:2000]))
        print("       flat median [2001:3000] = ",np.nanmedian(norm[2001:3000]))
    for i in range(0,4095) :
      if norm[i]>0. :
        F[i]    /= norm[i]
        #Fmin[i] /= norm[i]
        #Fmax[i] /= norm[i]
        Fminmax[0,i] /= norm[i]
        Fminmax[1,i] /= norm[i]
      else :
        F[i]    = np.nan
        #Fmin[i] = 0.
        #Fmax[i] = 0.
        Fminmax[0,i] = np.nan
        Fminmax[1,i] = np.nan

  # write measuremnts
  Fmed_1001_2000 = np.nanmedian(F[1001:2000])
  smed_1001_2000 =    np.nanstd(F[1001:2000])
  Fmed_2001_3000 = np.nanmedian(F[2001:3000])
  smed_2001_3000 =    np.nanstd(F[2001:3000])
  if verbose>=3:
    print("  Test: spec median [1001:2000] = %8.0f +- %5.0f " % (Fmed_1001_2000,smed_1001_2000))
    print("  Test: spec median [2001:3000] = %8.0f +- %5.0f " % (Fmed_2001_3000,smed_2001_3000))

  ADUtotcalc     = np.nanmedian(F[604:3005])
  sigADUtotcalc  = np.nanmedian([ np.nanstd(F[1204:1404]),np.nanstd(F[1404:1604]),np.nanstd(F[1604:1805]),np.nanstd(F[1805:2005]),np.nanstd(F[2005:2205]) ])
  if verbose>=2:
    print("  ADUtotcalc [604:3005] = %8.0f +- %5.0f " % (ADUtotcalc,sigADUtotcalc))

  #print("  spec aver [763:768] = ",np.average(F[763:768]))

  if ADUtot_in!=0. : 
    ADUtot   = ADUtot_in
    if verbose>=2:
      print("  ADUtot_in             = %8.0f" % (ADUtot_in))
  else :
    ADUtot   = np.nanmedian(F[604:3005])
  sigADUtot = sigADUtotcalc

  Ttot =  ADUtot*factor
  Tsys =  Ttot-Tobj
  Tobj_calc = Ttot-Tsys
  ADUsys = Tsys/factor # ADUtot - Tobj/factor
  ADUobj_calc = ADUtot-ADUsys

  if verbose>=1:
    print("  ADUtot = %7.0f +/- %6.0f     Ttot = %7.2f K = %7.2f C" % (ADUtot,sigADUtot,Ttot,Ttot-T0))
    print("  ADUobj = %7.0f                Tobj = %7.2f K = %7.2f C" % (ADUobj_calc,Tobj_calc,Tobj_calc-T0))
    print("  ADUsys = %7.0f                Tsys = %7.2f K = %7.2f C" % (ADUsys,Tsys,Tsys-T0))

  meastab.write('{0:>50s}   {1:14.0f} {2:12.0f}   {3:12.6g}   {4:10.2f} {5:10.2f} {6:10.2f} {7:10.2f}\n'.format(filename,ADUtot,sigADUtot,factor,Ttot,Tobj,Tsys,Tobj_calc))
                
  #################### plotting spectra #####################

  title=filename  
  ax.set_title(title)

  if plotTemp==1 :
    plt.ylabel('Temp [K]') 
  else:
    plt.ylabel('Analog Digital Units [ ]') 
  
  if xplot=='r': # register nr
    plt.xlabel('register nr [ ]')
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%5.0f'))
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    x_vec = reg
    x_min = -55
    x_max = 4095+55
  
  if xplot=='f':
    plt.xlabel('frequency [MHz]')  
    if freqtick!=0 : ax.xaxis.set_major_locator(MultipleLocator(freqtick))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%5.1f'))
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    x_vec = freq
    #* x_min = freq_min
    #* x_max = freq_max    
    x_min = freq0
    x_max = freq1
    
  if xplot=='v':
    # AUTO
    # v_min = (1.-freq_min/freqHI)*c - vlsr
    # v_max = (1.-freq_max/freqHI)*c - vlsr
    # print('  v_min = %6.2f [km/s]' % (v_min))
    # print('  v_max = %6.2f [km/s]' % (v_max))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # ax.plot( [0.,0.], [0.,1.], linestyle='-', color='green', linewidth=1)
    # ax.plot( [-vlsr,-vlsr], [0.,1.], linestyle='-', color='red', linewidth=1)
    plt.xlabel('velocity with respect to LSR: $v_{rad} - v_{lsr}$ [km/s]')
    if vtick!=0 : ax.xaxis.set_major_locator(MultipleLocator(vtick))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%5.0f'))
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    x_vec = veloc
    #*x_min = v_min
    #*x_max = v_max
    x_min = v1
    x_max = v0
  # freq_min = ax.get_xaxis().get_view_interval()[0]
  # freq_max = ax.get_xaxis().get_view_interval()[1]
  # ax.get_yaxis().set_ticks([])
  # plot min/max range as curve
  # ax.plot(freq,  Fmin, linestyle='-', color='yellow', linewidth=1)
  # ax.plot(freq,  Fmax, linestyle='-', color='yellow', linewidth=1)

  xx_vec = np.zeros((2,4095)) 
  for i in range(0,4095) :
    xx_vec[0,i] = x_vec[i]
    xx_vec[1,i] = x_vec[i]

  #Ftotmin = np.min(Fmin[205:3072])
  #Ftotmax = np.max(Fmax[205:3072])
  Ftotmin = np.min(Fminmax[0,205:3072])
  Ftotmax = np.max(Fminmax[1,205:3072])
  DF = Ftotmax - Ftotmin
  F00 = Ftotmin-0.1*DF
  F11 = Ftotmax+0.1*DF
  #   if F00<=0. :
  #     F00=Ftotmin/2. 
  # F00 = ADUtot-10.*sigADUtot
  # F11 = ADUtot+50.*sigADUtot

  ## plot spectrum
  if plotTemp==1 :
    if Temptick!=0 :
      ax.yaxis.set_major_locator(MultipleLocator(Temptick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%5.0f'))
    if Temp0==Temp1==0. :
      y_min = F00*factor
      y_max = F11*factor
    else :
      y_min = Temp0
      y_max = Temp1
    print("  plotting...")
    #for i in range(0,4095) :
    #  ax.plot([x_vec[i],x_vec[i]],  [Fmin[i]*factor,Fmax[i]*factor], linestyle='-', color='yellow', linewidth=1)
    ax.plot(xx_vec, Fminmax*factor, linestyle='-', color='yellow', linewidth=1)
    ax.plot( x_vec,       F*factor, linestyle='-', color='black',  linewidth=1)
  else :
    if ADUtick!=0 :
      ax.yaxis.set_major_locator(MultipleLocator(ADUtick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%5.0f'))
    if ADU0==ADU1==0. :
      y_min = F00
      y_max = F11
    else :
      y_min = ADU0
      y_max = ADU1
    print("  plotting...")
    # for i in range(0,4095) :
    #   ax.plot([x_vec[i],x_vec[i]],  [Fmin[i],Fmax[i]], linestyle='-', color='yellow', linewidth=1)    
    ax.plot(xx_vec, Fminmax, linestyle='-', color='yellow', linewidth=1)    
    ax.plot( x_vec, F,       linestyle='-', color='black',  linewidth=1)
  
  if plotTemp==1 :
    ax.plot( [x_min,x_max], [Ttot,Ttot], linestyle='-', color='red', alpha=0.5, linewidth=1, label='$T_{tot}$')
    if (Ttot!=Tsys) :
      ax.plot( [x_min,x_max], [Tsys,Tsys], linestyle='-', color='blue', alpha=0.5,  linewidth=1, label='$T_{sys}$')
  else : 
    ax.plot( [x_min,x_max], [ADUtot,ADUtot], linestyle='-', color='red', alpha=0.5, linewidth=1, label='$ADU_{tot}$')
    if (Ttot!=Tsys) :
      ax.plot( [x_min,x_max], [ADUsys,ADUsys], linestyle='-', color='blue', alpha=0.5,  linewidth=1, label='$ADU_{sys}$')
   
  # elif filename.find("S8")!=-1 : 
  # 
  #   # Tgauss = 51.1142007530  # Integral = lit integral
  #   Tgauss = 56.17669676653278 # Integral = 856. km/s
  #   Tsys = 0.
  #   correcttofit = 1.
  # 
  #   if filename=="150415_1655_S8.tab" :
  #     Tsys = 156122*factor
  #     correcttofit =  0.8293
  # 
  #   if filename=="150415_1655_S8.tab" :
  #     Tsys = 156122*factor
  #     correcttofit =  0.8293
  #  
  #   print("  Tsys = ",Tsys,"K"," (",Tsys/factor,")")
  #   print("  correcttofit = ",correcttofit)
  # 
  #   # lesen und plotten der Literaturwerte
  #   [v,w]=read_columns( os.path.join(dirflat+"1982A+A___106__190Kalberla_S8_1x1.tab"), [1,2])
  #   ffreq   = (1. - (     np.array(v)+vlsr)/c) * freqHI  # km/s -> HZ
  #   TBlit = np.array(w)
  #   # area_lit  = np.trapz(TAlit,v)  # not exact enough
  #   # print("  integral of measured Gaussian curve (green) = ",area_lit,"[K*km/s ]" # not exact enough)
  #   print("  literature value for T_B =                     856 [K*km/s ]")
  #   temp = Tsys + TBlit - 0.28   # manual correction by 0.28 K
  #   ax.plot(ffreq,temp,linestyle='-', color='red',   linewidth=1)
  # 
  #   # plotten von eigenem Gauss
  #   vlim0 =  -5.1
  #   vlim1 =  22.3
  #   #S8_vcen = (5.199+6.422)/2.  # km/s 
  #   S8_vcen = 6.7  # km/s 
  #   S8_fcen = (1. - (S8_vcen+vlsr)/c) * freqHI  # km/s -> HZ
  #   print( "  S8 peak v = ",S8_vcen," km/s")
  #   print( "  S8 peak f = ",S8_fcen," MHz")
  # 
  #   xg = np.zeros(1001)
  #   vg = np.zeros_like(xg)
  #   yg = np.zeros_like(xg)
  #   yvg = np.zeros_like(xg)
  #   cc = np.zeros_like(xg)
  # 
  #   N1 =  0.8; N2 = 1.-N1
  #   for i in range(0,1001) :
  #     xg[i] = ffreq[-1] + i*(ffreq[0]-ffreq[-1])/1000.
  #     # yg[i] = Tsys + correcttofit * Tgauss * exp(-pow(xg[i]-S8_fcen,2.)/0.0017) # only one Gaussian
  #     yg[i] = Tsys + correcttofit * Tgauss * ( N1*exp(-pow(xg[i]-S8_fcen,2.)/0.0013) + N2*exp(-pow(xg[i]-S8_fcen,2.)/0.008) )
  #     vg[i] = vlim0     + i*(vlim1-vlim0)/1000.
  #     fvg_i = (1. - ( vg[i]+vlsr)/c) * freqHI
  #     yvg[i] = Tsys + correcttofit * Tgauss * ( N1*exp(-pow(fvg_i-S8_fcen,2.)/0.0013) + N2*exp(-pow(fvg_i-S8_fcen,2.)/0.008) )
  #     cc[i] = Tsys 
  #   # area_meas = np.trapz(  yg-Tsys,   xg) 
  #   area_meas = np.trapz(  yvg-Tsys, vg) 
  #   print('  integral of measured Gaussian curve (green) = %6.2f [K*km/s ]' % (area_meas))
  #   # print('  ratio area_meas/area_lit = %6.2f percent' % (100.*area_meas/area_lit) # not exact enough)
  #   print('  ratio area_meas/area_lit = %6.2f percent' % (100.*area_meas/856.))
  #   print('  T_B / T_A(uncorrected) = %6.3f ' % (856./area_meas))
  # 
  #   if plotveloc==0:
  #     ax.plot(   xg,  yg,linestyle='-', color='green', linewidth=1)
  #     ax.plot(   xg,  cc,linestyle='-', color='grey', linewidth=1)
  #     ax.axis([1420.213-0.5,1420.213+0.5,Tsys-80.,Tsys+200.])
  #     ax.axis([1420.213-0.5,1420.213+0.5,Tsys-10.,Tsys+80.])      
  #     ax.axis([1420.213-0.5,1420.213+0.5,Tsys-50.,Tsys+150.])
  #     #ax.xaxis.set_minor_locator(MultipleLocator(0.02))
  #     #ax.xaxis.set_major_locator(MultipleLocator(0.2))
  #     ax.yaxis.set_minor_locator(MultipleLocator(1))
  #     ax.yaxis.set_major_locator(MultipleLocator(10))
  #
  # else :
  #   print(" object not recognized....")
  #   Ftotmin = np.min(Fmin[205:3072])
  #   Ftotmax = np.max(Fmax[205:3072])
  #   DF = Ftotmax - Ftotmin
  #   F00 = Ftotmin-0.1*DF
  #   F11 = Ftotmax+0.1*DF
  #   print("  ",F00,F11)
  #   if plotveloc==0:
  #     # ax.axis([1419.99,1420.01,F00,F11])  # test radio frequency 1420.0
  #     ax.axis([1416.5,1423.5,F00,F11])  
  #     #ax.set_yscale('log')
  #     #ax.yaxis.set_minor_locator(MultipleLocator(100))
  #     #ax.yaxis.set_major_locator(MultipleLocator(500))

  # if filename.find("_G")!=-1 : # Galaxy
  #   dx = (x_max-x_min)
  #   # G: ax.axis([x_min+0.2*dx,x_max-0.2*dx, F00,F11])
  #   #ax.yaxis.set_minor_locator(MultipleLocator(10))
  #   #ax.yaxis.set_major_locator(MultipleLocator(100))
  # else :
  #   print(" object not recognized....")
  #   Ftotmin = np.min(Fmin[205:3072])
  #   Ftotmax = np.max(Fmax[205:3072])
  #   DF = Ftotmax - Ftotmin
  #   F00 = Ftotmin-0.1*DF
  #   F11 = Ftotmax+0.1*DF
  #   print("  ",F00,F11)
  #   if F00<=0. :
  #     F00=Ftotmin/2. 
  #   if plotveloc==0:
  #     # ax.axis([1419.99,1420.01,F00,F11])  # test radio frequency 1420.0
  #     #ax.axis([1416.5,1423.5,F00,F11])  
  #     ax.axis([freq0,freq1,F00,F11])  
  #     #ax.set_yscale('log')
  #     #ax.yaxis.set_minor_locator(MultipleLocator(100))
  #     #ax.yaxis.set_major_locator(MultipleLocator(500))

  ax.axis([x_min,x_max,y_min,y_max])  

  if plotflat==1 :
    y_min = ax.get_yaxis().get_view_interval()[0]
    y_max = ax.get_yaxis().get_view_interval()[1]
    normscaled = norm/np.max(norm)*(y_max-y_min)+y_min
    #if plotveloc==0:
    ax.plot(x_vec, normscaled , linestyle='-', color='cyan', alpha=0.5, label="flat*.spec",  linewidth=1) 
  
  ###################################
      
  ax.yaxis.set_minor_locator(AutoMinorLocator(10))
  ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
  if xplot=='r':
    ax.plot( [regHI,regHI],            [y_min,y_max], linestyle='-', color='green', alpha=0.5,  linewidth=1, zorder=2, label='$reg_{HI}$')
  if xplot=='f':
    ax.plot( [freqHI,freqHI],          [y_min,y_max], linestyle='-', color='green', alpha=0.5,   linewidth=1, zorder=4, label='$f_{HI}$')
    ax.plot( [freqHI_lsr,freqHI_lsr],  [y_min,y_max], linestyle='-', color='magenta', alpha=0.5, linewidth=1, zorder=4, label='$f_{HI,lsr}$')
    if verbose>=3:   
      ax.plot( [freqHI_min,freqHI_min],  [y_min,y_max], linestyle='-', color='grey', alpha=0.5,    linewidth=3, zorder=3, label='$f_{HI,min}$')
      ax.plot( [freqHI_max,freqHI_max],  [y_min,y_max], linestyle='-', color='grey', alpha=0.5,    linewidth=3, zorder=3, label='$f_{HI,max}$')
  if xplot=='v':
    # ax.plot( [-vlsr,-vlsr], [y_min,y_max],               linestyle='-', color='green',   linewidth=1, zorder=2, label='$-v_{lsr}$')
    ax.plot( [0.,0.],       [y_min,y_max],                 linestyle='-', color='magenta', alpha=0.5, linewidth=1, zorder=4, label='$v_{rad} = v_{lsr}$')
    ax.plot( [vrmax,vrmax], [y_min,y_max], linestyle='-', color='blue', alpha=0.5,    linewidth=1, zorder=3, label='$v_{*}(R_{min})=220$ km/s')
    if verbose>=3:   
      ax.plot( [vlsr_min-vlsr,vlsr_min-vlsr], [y_min,y_max], linestyle='-', color='grey', alpha=0.5,    linewidth=3, zorder=3, label='$v_{lsr,min}-v_{lsr}$')
      ax.plot( [vlsr_max-vlsr,vlsr_max-vlsr], [y_min,y_max], linestyle='-', color='grey', alpha=0.5,    linewidth=3, zorder=3, label='$v_{lsr,max}-v_{lsr}$')

  #** if plotveloc==1: # in extrawindow
  #**   freq_min = ax.get_xaxis().get_view_interval()[0]
  #**   freq_max = ax.get_xaxis().get_view_interval()[1]
  #**   #print(freq_min,freq_max )
  #**   v_min = (1.-freq_min/freqHI)*c - vlsr
  #**   v_max = (1.-freq_max/freqHI)*c - vlsr
  #**   #print(v_min,v_max )
  #**   #plot = plt.subplot2grid( (15,1), (14,0) )
  #**   plot = plt.subplot2grid( (1,1),(0,0) )
  #**   ax.xaxis.set_minor_locator(MultipleLocator(5))
  #**   ax.xaxis.set_major_locator(MultipleLocator(50))
  #**   ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
  #**   ax.get_yaxis().set_ticks([])
  #**   ax.axis([v_min,v_max,0,1])
  #**   ax.plot( [0.,0.], [0.,1.], linestyle='-', color='green', linewidth=1)
  #**   ax.plot( [-vlsr,-vlsr], [0.,1.], linestyle='-', color='red', linewidth=1)
  #**   xlabel('espect to LSR [km/s]')  

  

  # plt.legend(loc='best')
  plt.legend(loc='upper right')
    
  #### output file
  if pdf==1 :
    plt.savefig(filename+"_"+xplot+".pdf")
  if png==1 :
    plt.savefig(filename+"_"+xplot+".png")

  if display :
  #   if plotveloc==0:
  #     buf = fig.canvas.mpl_connect('button_press_event', onclick_MHz)
  #   if plotveloc==1:
  #  buf = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
  else :
    plt.clf()    # clear the current figure

  # k = input("weiter mit <enter>")

meastab.close()
