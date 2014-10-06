# 1. Introduction

[HiTIME website](https://github.com/bjpop/HiTIME)

# 2. License

HiTIME is open source software, released under the [3-clause BSD license](https://github.com/bjpop/HiTIME/blob/master/LICENSE.txt).

# 3. Usage

    usage: hitime [-h] [--format FORMAT] [--intensityRatio R] [--rtWidth W]
                  [--rtSigma RTSIGMA] [--ppm P] [--mzWidth F] [--mzSigma MZSIGMA]
                  [--logFile FILENAME] [--mzDelta D] [--removeLow REMOVELOW]
                  [--outDir DIRECTORY] [--noScore] [--minSample MINSAMPLE]
                  inputFile outputFile
    
    Filter mass spec data for isotope doublets.
    
    positional arguments:
      inputFile             mass spec input data file
      outputFile            file name to save text data to
    
    optional arguments:
      -h, --help            show this help message and exit
      --format FORMAT       file format used for input mass spec data, options
                            are: mzml; mzdata
      --intensityRatio R    ratio of intensities for a doublet (isotope
                            amount/parent amount)
      --rtWidth W           Retention Time full width at half maximum in number of
                            scans
      --rtSigma RTSIGMA     Boundary for retention time width in standard
                            deviations
      --ppm P               m/z tolerance in parts per million
      --mzWidth F           m/z full width at half maximum in parts per million
      --mzSigma MZSIGMA     Boundary for mz window in standard deviations
      --logFile FILENAME    log progress in FILENAME
      --mzDelta D           m/z difference for doublets
      --removeLow REMOVELOW
                            Remove intensity values below the given signal level
      --outDir DIRECTORY    save output in DIRECTORY, if it does not exist it will
                            be created
      --noScore             process without scoring. Use for data exploration
      --minSample MINSAMPLE
                            minimum number of data points required in each sample
                            region
