// Overpotential
actPot = (potElectron - potProton);

// BV 
reactionRateBVAnodic.internalField() = exchCurrBV_ACL.value()*(Foam::exp(2.*0.5*actPot.internalField()*96485./(8.314*353.)));
reactionRateBVAnodic.boundaryField() = exchCurrBV_ACL.value()*(Foam::exp(2.*0.5*actPot.boundaryField()*96485./(8.314*353.)));

reactionRateBVCathodic.internalField() = exchCurrBV_ACL.value()*(-Foam::exp(-2.*0.5*96485.*actPot.internalField()/(8.314*353.)));
reactionRateBVCathodic.boundaryField() = exchCurrBV_ACL.value()*(-Foam::exp(-2.*0.5*96485.*actPot.boundaryField()/(8.314*353.)));

// Form the Sum
reactionRateBVSum = (reactionRateBVAnodic + reactionRateBVCathodic);
