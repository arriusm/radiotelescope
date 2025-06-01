#!/Users/arri/miniconda3/bin/python

VERSION='V250531'

import sys, os, math, random, time
import matplotlib.pyplot as plt
import numpy as np
#from numpy import *
from pylab import *
from math import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import argparse

parser = argparse.ArgumentParser(description='radio_spectra.py by Arno Riffeser '+VERSION)
parser.add_argument('-r', '--medrad', type=int,   dest='medrad',   default='50',                 help='median radis')
parser.add_argument('FILES',             nargs='+',                                            help='files')
args = parser.parse_args()

# Steuerung
# dir='/Users/arri/ASTRO/Radioteleskop/DATA/'
# dir='/home/radioadmin/Work/srtnver6/DATA/171024'
# dir='/Users/arri/ASTRO/Radioteleskop/DATA/171024'
# table='150227_1332_z0945_Absorber_T13.6.tab'; medrad = 10  # flat.spec flat1.spec
# table='150508_1439_z0940_Absorber_T28.0.tab'; medrad = 10  # flat2.spec
# table='150701_1718_Absorber_T30.tab';         medrad = 10  # flat3.spec
# table='151124_1449_absorbercold.tab'; medrad = 10   # 151124_astroprak/15-11-24/flat.spec
# table='161129_1556_absorber_white_-1C.tab'; medrad = 5
# table='171024_135747_stow_blackside1.tab'; medrad = 10  

dir=''
table=args.FILES[0];
medrad = args.medrad

create=1
display=1
pdf=1
png=1


# theoretical from ADC specs with 20 MHz clocking and 1416 MHz Mixing
fstart  = 1416.+1e-6/(8192./2e7)
fend    = 1416.+10. # 10 = 20 MHz / 2.
istart  = 1
iend    = 4096
# freqsep = (fend-fstart)/(iend-istart)
# print('freq_sep = %10.5f kHz' % (freqsep*1000.))
# print('1001: %15.10f MHz' % (1416.+freqsep*1001.))
# print('2000: %15.10f MHz' % (1416.+freqsep*2000.))
# print('2001: %15.10f MHz' % (1416.+freqsep*2001.))
# print('3000: %15.10f MHz' % (1416.+freqsep*3000.))

freq = np.zeros(4095)
regnr = np.zeros(4095)
for i in range(0,4095) :
  #freq[i] = fstart+(i+1-istart)*freqsep
  regnr[i] = i
  #print(i,freq[i])

#print('freq_min = %10.5f MHz' % (freq[0]))
#print('freq_max = %10.5f MHz' % (freq[4094]))


print('flat: ',table)
dtable = os.path.join(dir,table)
in_file  = open( str( dtable ), 'r' )

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

F = np.zeros_like(freq)
#Fmin = np.zeros_like(freq)+1e20
#Fmax = np.zeros_like(freq)-1e20
Fminmax = np.zeros((2,4095)) 
Fminmax[0,]=1e20 
Fminmax[1,]=-1e20
#Fsingle = np.zeros_like(freq)
 
factor=1.

#fig, ax = plt.subplots()
fig = plt.figure(num=1, figsize=(13, 7), dpi=100, facecolor='w', edgecolor='k')
plot = plt.subplot2grid( (13, 1), (0, 0), rowspan=13)

# reading spectra 
nr = 0
line=' '
while line != '':
  line = in_file.readline()
  if line!='' and line[:1]!='#' and line[:1]!=' ':
    nr+=1
    v = [s for s in line.split()]
    for i in range(0,4095):
      value=float(v[27+i])*factor
      F[i]+=value
      #Fsingle[i]=value
      if value<Fminmax[0,i]:
        Fminmax[0,i]=value
      if value>Fminmax[1,i]:
        Fminmax[1,i]=value            
    #plot.scatter(freq, Fsingle, s=1, color='yellow', marker='o')

print('  lines = ',nr)

F/=nr # mean

print('median smoothing...')
mF=np.zeros_like(F)
for i in range(0+medrad,4095-medrad) :
  mF[i] = np.median(F[i-medrad:i+1+medrad])

if create==1:
  print('creating flat*.spec...')
  flatspecname = 'flat_'+table+'_'+str(medrad)+'.spec'
  print(flatspecname)
  flatspec = open(os.path.join(dir,flatspecname), 'w')
  print('writing: ',flatspecname)
  idx=np.array(range(0,4095))
  for i in range(0,4095) :
    flatspec.write('{0:10.0f} {1:10.4f}\n'.format(idx[i],mF[i]))
  flatspec.close()
  
# norm = mF/np.max(mF)
# print('  aver([1001:2000]) = ',np.average(norm[1001:2000]))
# print('  aver([2001:3000]) = ',np.average(norm[2001:3000]))

# plot min/max range as curve
# plot.plot(freq,  Fmin, linestyle='-', color='yellow', linewidth=1)
# plot.plot(freq,  Fmax, linestyle='-', color='yellow', linewidth=1)    

print('plotting...')

xx_vec = np.zeros((2,4095)) 
for i in range(0,4095) :
  xx_vec[0,i] = regnr[i]
  xx_vec[1,i] = regnr[i]

# plot min/max range as errorbars
#for i in range(0,4095) :
#plot.plot([freq[i],freq[i]], [Fmin[i],Fmax[i]], linestyle='-', color='yellow', linewidth=1)
#plot.plot([regnr[i],regnr[i]], [Fmin[i],Fmax[i]], linestyle='-', color='yellow', linewidth=1)
plot.plot(xx_vec, Fminmax, linestyle='-', color='yellow', linewidth=1)

# plot spectrum
#plot.plot(freq, F, linestyle='-', color='black', linewidth=1)
plot.plot(regnr, F, linestyle='-', color='black', linewidth=1)

title=flatspecname
plot.set_title(title)

# xlabel('frequency [MHz]')  
xlabel('register nr [ ]')  
ylabel('Analog Digital Units [ ]')

Fmax=np.max(F[1000:3000])
print(Fmax)
# plot.axis([1416.5,1423.5,100000.,700000.])
#plot.axis([0.,4096.,0.,700000.])
plot.axis([0.,4096.,0.,Fmax*1.1])

#plot.xaxis.set_minor_locator(MultipleLocator(0.05))
#plot.xaxis.set_major_locator(MultipleLocator(0.5))
#plot.xaxis.set_major_formatter(FormatStrFormatter('%5.1f'))
plot.yaxis.set_minor_locator(MultipleLocator(10000.))
plot.yaxis.set_major_locator(MultipleLocator(50000.))
plot.yaxis.set_major_formatter(FormatStrFormatter('%d'))

# freqHI=1420.406
# plot.plot([freqHI,freqHI],[0.,1000000.], linestyle='-', color='grey', linewidth=1)

#plot.plot(freq,mF,linestyle='-',color='red',linewidth=1,label='smoothed') 
plot.plot(regnr,mF,linestyle='-',color='red',linewidth=1,label='smoothed') 

plt.legend(loc='upper right')

plotfile = os.path.join(dir,flatspecname)

if pdf==1 :
  savefig(plotfile+'.pdf')
if png==1 :
  savefig(plotfile+'.png')

# mouse click
def onclick(event):
  print('  %9.4f MHz     %10.0f ' % (event.xdata, event.ydata))
  return

if display :
  #buf = fig.canvas.mpl_connect('button_press_event', onclick)
  show()
else :
  clf()    # clear the current figure

