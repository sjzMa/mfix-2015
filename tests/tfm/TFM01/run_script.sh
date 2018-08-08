#!/bin/bash -exl

# set case directory
export CASE_DIR=`pwd`

# load modules
module load autoconf/2.69
module load gnu/4.6.4 
module load openmpi/1.5.5_gnu4.6

# compile MFIX in ./src/
echo "******** Compiling MFIX..."
cd $CASE_DIR
../../../configure_mfix --enable-dmp FC=mpif90 FCFLAGS="-O0 -g"
make
#../../../model/make_mfix.old --dmp --opt=O3 --compiler=gcc --exe=mfix.exe -j

cd $CASE_DIR

# Run mesh_32 (i.e., 32x32 for 2D, 32x32x32 for 3D)
echo "******** Running mesh_32..."
./mfix imax=32 jmax=32 > out.log
rm -f $CASE_DIR/{TFM01.*,out.log}
#rm -f $CASE_DIR/de_norms.dat

rm -f $CASE_DIR/mfix

echo "******** Done."


# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
