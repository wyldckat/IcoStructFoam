//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/readCoupling.H 1906 2007-08-28T16:16:19.392553Z bgschaid  $ 

Info << "Reading coupling parameters\n" << endl;

IOdictionary coupling
(
 IOobject
 (
  "couplingParameters",
  runTime.constant(),
  mesh2,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  )
 );

word fluidSide=coupling.lookup("fluidSide");
word solidSide=coupling.lookup("solidSide");

dimensionedScalar startMeshMotion(coupling.lookup("startMeshMotion"));
dimensionedScalar motionRelaxation(coupling.lookup("motionRelaxation"));

label fluidSideI = mesh1.boundaryMesh().findPatchID(fluidSide);

if(fluidSideI<0) {
  FatalErrorIn("Coupled Solver ") << "Patch " << fluidSide << " does not exist in mesh1" << exit(FatalError);
 }

label solidSideI = mesh2.boundaryMesh().findPatchID(solidSide);

if(solidSideI<0) {
  FatalErrorIn("Coupled Solver ") << "Patch " << solidSide << " does not exist in mesh2" << exit(FatalError);
 }
Info << "Coupling boundary " <<  fluidSide << "("<< fluidSideI << ") of fluid with " <<
  "boundary " <<  solidSide << "("<< solidSideI << ") of solid\n" << endl;

const polyPatch &fluidMesh=mesh1.boundaryMesh()[fluidSideI];
const polyPatch &solidMesh=mesh2.boundaryMesh()[solidSideI];

// Sanity checks for the coupling

if(solidMesh.size()!=fluidMesh.size()) {
  FatalErrorIn("Coupled Solver ") << "Sizes do not match" << exit(FatalError);
}

labelList exchange(fluidMesh.size());
List<vector> oldPoints(fluidMesh.size());
List<vector> dispVals(fluidMesh.size());

{
  const scalar tolerance=1e-10;

  forAll(fluidMesh,fI) {
    const vector &fV=fluidMesh.faceCentres()[fI];
    scalar minDist=1e20;
    label minInd=-1;

    forAll(solidMesh,sI) {
      const vector &sV=solidMesh.faceCentres()[sI];
      scalar dist=mag(fV-sV);

      if(dist<minDist) {
	minDist=dist;
	minInd=sI;
      }
    }
    if(minDist>tolerance) {
      FatalErrorIn("Coupled Solver ") << "Face centres do not match" << exit(FatalError);
    }

    exchange[fI]=minInd;
    oldPoints[fI]=fV;
  }
}							       
