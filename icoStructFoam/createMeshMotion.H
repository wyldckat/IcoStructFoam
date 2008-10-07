//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/createMeshMotion.H 1906 2007-08-28T16:16:19.392553Z bgschaid  $ 

PrimitivePatchInterpolation<polyPatch> inter(solidMesh);

if(solidMesh.meshPoints().size()!=fluidMesh.meshPoints().size()) {
  FatalErrorIn("Coupled Solver ") << "Sizes do not match" << exit(FatalError);
 }	

autoPtr<Foam::motionSolver> motionPtr = motionSolver::New(mesh1);

bool displacementMotionSolver=false;

if(mesh1.objectRegistry::foundObject<volVectorField>("cellDisplacement")) {
    displacementMotionSolver=true;
    Info << " Motionsolver is displacement based\n" << endl;
} else if(mesh1.objectRegistry::foundObject<volVectorField>("cellMotionU")) {
    displacementMotionSolver=false;
    Info << " Motionsolver is motionU based\n" << endl;
} else {
    FatalErrorIn("createMeshMotion.H")
        << " Can't deal with this kind of motion solver\n" << endl
            << abort(FatalError);
}

tractionDisplacementFvPatchVectorField& displace= 
 refCast<tractionDisplacementFvPatchVectorField>(D.boundaryField()[solidSideI]);
