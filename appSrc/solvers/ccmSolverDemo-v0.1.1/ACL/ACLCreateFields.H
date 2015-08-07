// Initialise electron field pointer lists
#include "ACLFieldPointerLists.H"

// Populate ACL field pointer lists
forAll(ACLRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << ACLRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include "ACLReadComponentDictionary.H"
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include "ACLDimScalarAssign.H"
	
    // Create Field Objects for solved variables
    #include "ACLSolvedFieldObjects.H"

    // Create Field Objects for transport coefficient variables
    #include "ACLCoefficientFieldObjects.H"

    // Create Field Objects for source term variables
    #include "ACLSourceTermFieldObjects.H"

	// Create Field Objects for the calculated variables
	#include "ACLCalculatedFieldObjects.H"
}
