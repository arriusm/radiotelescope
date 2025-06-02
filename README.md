# radio_spectra.py
by Arno Riffeser V250601

## usage:
```
radio_spectra.py [-h] [-d DIR] [-m] [-n] [-f] [-o] [-p] [-g] [-T]
                 [-x XPLOT] [-v VERBOSE]
                 FILES [FILES ...]
```

## positional arguments:
```
  FILES                 files
```

## options:
```
  -h, --help                     show this help message and exit
  -d DIR, --dir DIR              directory, default is current directory
  -m, --mask                     masking
  -n, --no-display               turn off showing on display
  -f, --makeflat                 applies a flat division
  -o, --plotflat                 plots the flat
  -p, --pdf                      stores a pdf
  -g, --png                      stores a png
  -T, --plotTemp                 plot Temp on y axis
  -x XPLOT, --xplot XPLOT        xaxis: r - register, f - frequency, v - velocity
  -v VERBOSE, --verbose VERBOSE  verbose level
```


# radio_flat.py
by Arno Riffeser V250602

## usage:
```
radio_flat.py [-h] [-r MEDRAD] FILES [FILES ...]
```

## positional arguments:
```
  FILES                 files
```

## options:
```
  -h, --help                    show this help message and exit
  -r MEDRAD, --medrad MEDRAD    median radis
```

# radio_header.py
by Arno Riffeser V250531

## usage:
```
radio_header.py [-h] [-k KEY] [-v VALUE] FILES [FILES ...]
```

## positional arguments:
```
  FILES                 files
```

## options:
```
  -h, --help                 show this help message and exit
  -k KEY, --key KEY          header keyword name
  -v VALUE, --value VALUE    header keyword value
```


# npoints.py
by Arno Riffeser V250531

## usage:
```
npoints.py [-h] [-v VERBOSE] FILES [FILES ...]
```


## positional arguments:
```
  FILES                 files
```

## options:
```
  -h, --help                     show this help message and exit
  -v VERBOSE, --verbose VERBOSE  verbose level
```
