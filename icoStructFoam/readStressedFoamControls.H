//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/readStressedFoamControls.H 1280 2007-04-27T13:55:25.796383Z bgschaid  $ 
   
const dictionary& stressControl = mesh2.solutionDict().subDict("stressAnalysis");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("D")));
Switch compactNormalStress(stressControl.lookup("compactNormalStress"));
