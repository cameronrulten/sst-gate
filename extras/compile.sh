#!/bin/bash
infile=${1}
outfile=${2}

clang++ -ggdb ${infile} -o ${outfile} `root-config --cflags --libs --glibs` -L/$ROOTLIB -lXMLIO -lMLP -lTreePlayer -lMinuit -lTMVA -lGeom -L/$ROBAST -lROBAST -I/$ROBAST/include -I${UTILITIES}/include -I${BOOST_ROOT}/
#-I${SST_INCLUDE_DIR}/include
