// Overpotential
actPot = ((potElectron - potProton) - OCV);

// Reaction Rates
reactionRateBVAnodic.internalField() = exchCurrBV_CCL.value()*(Foam::exp(2.*0.5*96485.*actPot.internalField()/(8.314*353.)));
reactionRateBVAnodic.boundaryField() = exchCurrBV_CCL.value()*(Foam::exp(2.*0.5*96485.*actPot.boundaryField()/(8.314*353.)));

reactionRateBVCathodic.internalField() = exchCurrBV_CCL.value()*(-Foam::exp(-2.*0.5*96485.*actPot.internalField()/(8.314*353.)));
reactionRateBVCathodic.boundaryField() = exchCurrBV_CCL.value()*(-Foam::exp(-2.*0.5*96485.*actPot.boundaryField()/(8.314*353.)));

// Form the Sum
reactionRateBVSum = reactionRateBVAnodic + reactionRateBVCathodic;
