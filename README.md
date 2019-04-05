# sst-gate
Software to perform ray-tracing simulations for the SST-GATE telescope

The following code can be used for ray-tracing studies of the SST-GATE
telescope. This code was written in c++ and has been tested using the
clang++ compiler: Apple LLVM version 8.1.0 (clang-802.0.42) on an
Apple Macbook Pro running macOS Sierra version 10.12.6

# Prerequisites:

In order to use this software properly a working copy of the **CERN ROOT** analysis framework is required.
Currently this software has only been tested on the production version ROOT v6.10.06 (see https://root.cern.ch)

In addition, for this software to work a copy of the **ROBAST** software version 2.4.4 must be installed.

You must define the the ROBAST environment variable so CMake can find your ROBAST installation.

i.e. export ROBAST=path/to/your/ROBAST/installation

For more information about ROBAST please see http://robast.github.io/

The latest ROBAST releases are available from https://github.com/ROBAST/ROBAST/releases

In addition you need to have **CMake** installed.
This software was tested with CMake version 3.7.2
The minimum version supported is version 2.8

# Installation:

To install the SST-GATE simulation software:

 ```git clone sst-gate-v2-0-0-source```

 ```cd sst-gate-v2-0-0-source```
 
 ```git checkout -b v2-0-0 v2-0-0```

 ```cd ../```
 
  ```mkdir sst-gate-v2-0-0-build sst-gate-v2-0-0-install```
  
  ```cd sst-gate-v2-0-0-build```
  
  ```cmake -DCMAKE_INSTALL_PREFIX=/full/path/to/sst-gate-v2-0-0-install ../sst-gate-v2-0-0-source/```
  
The cmake build system supports parallel builds. You can specify the option -jN where N is the number of cores you have available.
  
  ```cmake --build . -- -j2```
  
  ```cmake --build . --target install```
  
  ```export PATH=${PATH}:/path/to/sst-gate-v2-0-0-install/bin```

You are now ready to test the SST-GATE software, you can do this from the command line by running:

 ```./sst_gate_test 2.0 test-run square 0.0 0.0```

ZEMAX data files are contained within the zemax_data_files directory located within the extras directory.
Also included in the zemax_data_files directory is an internal technical note on the optical design written by Jürgen Schmoll (University of Durham).

# License:
If you use or fork this software you must credit the original authors in all your works.
In particular you should always include the following reference:

```
Simulating the optical performance of a small-sized telescope with secondary optics
for the Cherenkov Telescope Array

by Cameron Rulten, Andreas Zech, Akira Okumura, Philippe Laporte and Jürgen Schmoll

Astroparticle Physics
Volume 82, September 2016, Pages 36-48
https://doi.org/10.1016/j.astropartphys.2016.05.002
```

You can access the publication under a creative commomns license: http://www.sciencedirect.com/science/article/pii/S0927650516300652


Cameron Rulten 2017
