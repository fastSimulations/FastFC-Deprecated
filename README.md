# FastFC

Source the fastfBASHRC file:

open $HOME/.bashrc in an editor of your choice and add the following line:

source $WM_PROJECT_USER_DIR/FastFC/etc/fastfcBASHRC

save and close the editor and type the following in the terminal:

source $HOME/.bashrc

Compile the source code:

./buildFastFC

Extract the case files to your run directory using the following steps:

type "fastfc" in the terminal and this will take you to the fastfc source directory

type ./extractCaseFiles in the terminal and this will copy the case files to $FOAM_RUN

type "run" in the terminal and this will take you to your FOAM and FastFC case file directory

You will see a directory called "defaultCaseFiles", go into this directory:

cd defaultCaseFiles/mixedScalarMEA.v0.0.3.master

type the following and insert your desired casefile name in "newCaseName" 

fcCopyCase -d NewCaseName

type "run" in the terminal

cd newCaseName

type the following to setup and run the case file:

fcCleanCase

fcCaseSetup

fcSerialRun

Post processing is done via paraview and polarization data is located in polarizationData.txt
