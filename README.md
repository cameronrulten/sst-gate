# sst-gate
Software to perform ray-tracing simulations for the SST-GATE telescope

The following code can be used for ray-tracing studies of the SST-GATE
telescope. This code was written in c++ and has been tested using the
clang++ compiler: Apple LLVM version 8.1.0 (clang-802.0.42) on an
Apple Macbook Pro running macOS Sierra version 10.12.6

In in order to use this software properly a working copy of the CERN ROOT analysis framework is required.
Currently this software has only been tested on the production version ROOT v6.10.06 (see https://root.cern.ch)

In addition, for this software to work a copy of the ROBAST software version 2.4.4 must be installed.
(see http://robast.github.io/)
The latest ROBAST releases are available from https://github.com/ROBAST/ROBAST/releases/tag/v2.4.4

In addition you need to have CMake installed. This was tested with CMake version 3.7.2
although the minimum version support is version 2.8

To install the SST-GATE simulation software:

git clone sst-gate-v2-0-0-source
cd sst-gate-v2-0-0-source
git checkout -b v2-0-0 v2-0-0
cd ../
mkdir sst-gate-v2-0-0-build sst-gate-v2-0-0-install
cd sst-gate-v2-0-0-build
cmake -DCMAKE_INSTALL_PREFIX=/full/path/to/sst-gate-v2-0-0-install ../sst-gate-v2-0-0-source/

The cmake build system supports parallel builds. You can specify the option -jN where N is the number of cores you have available.
cmake --build . -- -j2

cmake --build . --target install

export PATH=${PATH}:/path/to/sst-gate-v2-0-0-install/bin

You are now ready to test the SST-GATE software, you can do this from the command line by running:
./sst_gate_test 2.0 test-run square 0.0 0.0

ZEMAX data files are contained within the zemax_data_files directory located within the extras directory.
Also included in the zemax_data_files directory is an internal technical note on the optical design written by Jürgen Schmoll (University of Durham).
