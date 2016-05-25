#! /bin/bash
set +e
clear
echo
echo "Fast-FC Setup and Compilation Script"
echo
echo "Building the Fast-FC Libraries"
cd libSrc
./Allwclean
./Allwmake
echo
echo "Done"
echo
echo 
echo "Compiling the Fast-FC Solvers"
cd ../appSrc
./Allwclean
./Allwmake
cd ..
echo
echo "Done"
echo
