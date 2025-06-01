#!/Users/arri/miniconda3/bin/python

VERSION='V250531'

import argparse, os

parser = argparse.ArgumentParser(description='radio_header.py by Arno Riffeser '+VERSION)

parser.add_argument('-k', '--key',   dest='key',   default='', help='header keyword name')
parser.add_argument('-v', '--value', dest='value', default='', help='header keyword value')
parser.add_argument('FILES',         nargs='+',                help='files')

args = parser.parse_args()

key   = args.key
value = args.value

for FILE in args.FILES:
  if key=='' :     
    print(str( FILE ))
  infile   = open( str( FILE ), 'r' )
  outfile  = open( 'bla', 'w' )
  line=' '
  while line != '':
    line = infile.readline()
    if value=='' :
      if line[:1]=='#' :
        v = [s for s in line.split()]
        if key=='' :      
          print('    %10s  %14s' % (v[1],v[2]))
        else :
          if key == v[1] :
            print('%45s  %10s  %14s' % (FILE,v[1],v[2]))
    else :
      if line[:1]=='#' :
        v = [s for s in line.split()]
        if key == v[1] :
          print('%45s  %10s  %14s     ---> %14s' % (FILE,v[1],v[2],value))
          outfile.write('# {0:>10s} {1:>15s}\n'.format(key,value))
        else :
          outfile.write(line)   
      else :
        outfile.write(line)   
  if key=='' :     
    print('') 
  infile.close()
  outfile.close()
  if value!='' :
    os.rename(FILE, FILE+'.bak')
    os.rename('bla', FILE)
