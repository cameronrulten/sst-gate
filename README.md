# sst-gate
Software to perform ray-tracing simulations for the SST-GATE telescope

The following code can be used for ray-tracing studies of the SST-GATE
telescope. This code was written in c++ and has been tested using both the g++
and clang++ compilers on a Macbook Pro OSX 10.9

In in order to use this software properly a working copy of the CERN ROOT analysis framework is required.
Currently this software has only been tested with ROOT v5-34 (see https://root.cern.ch)

In addition, for this software to work a copy of the ROBAST software must be installed.
(see http://sourceforge.net/projects/robast)

For convenience a sample compile script (compile.sh) has been provided to demonstrate how the software can be compiled linking to the ROOT and ROBAST libraries.

The compile script is simple to run from the command line as follows:

./compile.sh <input cpp file> <output executable file>

e.g. ./compile sst_gate_test.cpp sst_gate_test

I will construct a Makefile to make things easier... please be patient :)
