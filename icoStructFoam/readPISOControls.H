//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/readPISOControls.H 943 2007-01-25T21:38:27.500728Z bgschaid  $ 
    const dictionary& piso = mesh1.solutionDict().subDict("PISO");

    int nCorr(readInt(piso.lookup("nCorrectors")));

    int nNonOrthCorr = 0;
    if (piso.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(piso.lookup("nNonOrthogonalCorrectors"));
    }

    bool momentumPredictor = true;
    if (piso.found("momentumPredictor"))
    {
        momentumPredictor = Switch(piso.lookup("momentumPredictor"));
    }

    bool fluxGradp = false;
    if (piso.found("fluxGradp"))
    {
        fluxGradp = Switch(piso.lookup("fluxGradp"));
    }

    bool transSonic = false;
    if (piso.found("transSonic"))
    {
        transSonic = Switch(piso.lookup("transSonic"));
    }

    int nOuterCorr = 1;
    if (piso.found("nOuterCorrectors"))
    {
        nOuterCorr = readInt(piso.lookup("nOuterCorrectors"));
    }
