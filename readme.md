# To build Linux 
Please see mfix_usr_guide.pdf


# Cross-Compile 

## Ubuntu or WSL(Windows10)

### Cross-Compile tool-clain

sudo apt install -y --no-install-recommends mingw-w64 build-essentail

git clone https://github.com/sjzMa/mfix-2015.git

cd mfix-2015

mkdir trybuild && cd mkdir trybuild

####  computing with parallel
../configure_mfix --host=x86_64-w64-mingw32 FCFLAGS="-O2" FFLAGS="-O2" --enable-smp
####  computing without parallel
../configure_mfix --host=x86_64-w64-mingw32 FCFLAGS="-O2" FFLAGS="-O2"

make -j4 LDFLAGS="-static -Wl,--stack,16777216"

Last, you will get standalone EXE Program