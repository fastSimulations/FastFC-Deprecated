// Initialise electron field pointer lists
#include "PMEMFieldPointerLists.H"

// Populate PMEM field pointer lists
forAll(PMEMRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << PMEMRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include "PMEMReadComponentDictionary.H"
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include "PMEMDimScalarAssign.H"
	
    // Create Field Objects for solved variables
    #include "PMEMSolvedFieldObjects.H"

    // Create Field Objects for transport coefficient variables
    #include "PMEMCoefficientFieldObjects.H"

    // Create Field Objects for source term variables
    #include "PMEMSourceTermFieldObjects.H"

	// Create Field Objects for the calculated variables
	#include "PMEMCalculatedFieldObjects.H"

}
