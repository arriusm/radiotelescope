radio_spectra.py by Arno Riffeser V250601

usage:
'''
radio_spectra.py [-h] [-d DIR] [-m] [-n] [-f] [-o] [-p] [-g] [-T]
                        [-x XPLOT] [-v VERBOSE]
                        FILES [FILES ...]
'''

options:
'''
  -h, --help            show this help message and exit
  -d DIR, --dir DIR     directory, default is current directory
  -m, --mask            masking
  -n, --no-display      turn off showing on display
  -f, --makeflat        applies a flat division
  -o, --plotflat        plots the flat
  -p, --pdf             stores a pdf
  -g, --png             stores a png
  -T, --plotTemp        plot Temp on y axis
  -x XPLOT, --xplot XPLOT
                        xaxis: r - register, f - frequency, v - velocity
  -v VERBOSE, --verbose VERBOSE
                        verbose level
'''