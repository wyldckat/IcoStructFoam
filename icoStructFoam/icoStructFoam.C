/*---------------------------------------------------------------------------*\
 ##   ####  ######     | 
 ##  ##     ##         | Copyright: ICE Stroemungsfoschungs GmbH
 ##  ##     ####       |
 ##  ##     ##         | http://www.ice-sf.at
 ##   ####  ######     |
-------------------------------------------------------------------------------
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    icoStructFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids 
    coupled with structural mechanics

    Based on icoFoam and solidDisplacementFoam
    Coupling after an idea by conjugateFoam

 ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/icoStructFoam.C 1906 2007-08-28T16:16:19.392553Z bgschaid  $ 
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Switch.H"
#include "fixedValueFvPatchFields.H"
#include "tractionDisplacement/tractionDisplacementFvPatchVectorField.H"
#include "fvMesh.H"
#include "PrimitivePatchInterpolation.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMeshes.H"

#   include "readMechanicalProperties.H"
#   include "readStressedFoamControls.H"
#   include "readThermalProperties.H"

#   include "createIcoFields.H"
#   include "createStructureFields.H"

#   include "readCoupling.H"

#   include "createMeshMotion.H"

#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

	if(runTime.value()>=startMeshMotion.value()) {
	  Info << "\nMoving mesh\n" << endl;

          // Make the fluxes absolute
          //          fvc::makeAbsolute(phi, U);

          //	  dispVals=inter.faceToPointInterpolate(displace);
          dispVals=displace;

	  scalar maxDist=-1e10;
	  scalar maxAway=-1e10;

          if(displacementMotionSolver) {
              volVectorField::GeometricBoundaryField &meshDisplacement=
                  const_cast<volVectorField&>(mesh1.objectRegistry::lookupObject<volVectorField>("cellDisplacement")).boundaryField();

              if(typeid(meshDisplacement[fluidSideI])!=typeid(fixedValueFvPatchField<vector>)) {
                  FatalErrorIn("Coupled Solver ") << "Fluid side not movable" << exit(FatalError);
              }
              
              vectorField &mDisp=refCast<vectorField>(meshDisplacement[fluidSideI]);
              
              scalar factor=1/motionRelaxation.value();

              forAll(fluidMesh,fluidI) {
                  label solidI=exchange[fluidI];

                  //	    vector here=mesh1.points()[fluidPoints[fluidI]];
                  vector here=fluidMesh.faceCentres()[fluidI];
                  vector old =oldPoints[fluidI];
                  vector disp=dispVals[solidI];
                  vector neu=disp*factor+(here-old)*(1-factor);
                  vector move=neu;

                  mDisp[fluidI]=move;

                  if(mag(here-old)>maxDist) {
                      maxDist=mag(here-old);
                  }
                  if(mag(neu-here)>maxAway) {
                      maxAway=mag(here-neu);
                  }
              }

              
          } else {
              volVectorField::GeometricBoundaryField &motionU=
                  const_cast<volVectorField&>(mesh1.objectRegistry::lookupObject<volVectorField>("cellMotionU")).boundaryField();
              
              if(typeid(motionU[fluidSideI])!=typeid(fixedValueFvPatchField<vector>)) {
                  FatalErrorIn("Coupled Solver ") << "Fluid side not movable" << exit(FatalError);
              }
              
              vectorField &patchU=refCast<vectorField>(motionU[fluidSideI]);

              //	  tetPointVectorField neu=inter.faceToPointInterpolate(patchU);
              scalar factor=1/(runTime.deltaT().value()*motionRelaxation.value());
              //          const labelList& fluidPoints=fluidPatch.meshPoints();

              //	  forAll(fluidPoints,fluidI) {
              forAll(fluidMesh,fluidI) {
                  label solidI=exchange[fluidI];

                  //	    vector here=mesh1.points()[fluidPoints[fluidI]];
                  vector here=fluidMesh.faceCentres()[fluidI];
                  vector old =oldPoints[fluidI];
                  vector disp=dispVals[solidI];
                  vector neu=disp+old;
                  vector move=factor*(neu-here);

                  patchU[fluidI]=move;

                  if(mag(here-old)>maxDist) {
                      maxDist=mag(here-old);
                  }
                  if(mag(neu-here)>maxAway) {
                      maxAway=mag(here-neu);
                  }
              }
          }
	  mesh1.movePoints(motionPtr->newPoints());
          //          U.correctBoundaryConditions();

	  Info << "\nBiggest movement: " << maxDist << " Bigges divergence " << maxAway <<  endl;

          // Make the fluxes relative
          //          fvc::makeRelative(phi, U);
	}

	Info << "Solving flow in mesh1\n" << endl;
	{
#         include "readPISOControls.H"
#         include "CourantNo.H"

	  fvVectorMatrix UEqn
	    (
	     fvm::ddt(U)
	     + fvm::div(phi, U)
	     - fvm::laplacian(nu, U)
	     );

	  solve(UEqn == -fvc::grad(p));

	  // --- PISO loop

	  for (int corr=0; corr<nCorr; corr++)
	    {
	      volScalarField rUA = 1.0/UEqn.A();

	      U = rUA*UEqn.H();
	      phi = (fvc::interpolate(U) & mesh1.Sf()) 
                + fvc::ddtPhiCorr(rUA, U, phi);

	      adjustPhi(phi, U, p);

	      for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
		{
		  fvScalarMatrix pEqn
		    (
		     fvm::laplacian(rUA, p) == fvc::div(phi)
		     );

		  pEqn.setReference(pRefCell, pRefValue);
		  pEqn.solve();

		  if (nonOrth == nNonOrthCorr)
		    {
		      phi -= pEqn.flux();
		    }
		}

#             include "continuityErrs.H"

	      U -= rUA*fvc::grad(p);
	      U.correctBoundaryConditions();
	    }

	}

	Info << "\nSolving structure in mesh2\n" << endl;

	{
#         include "readStressedFoamControls.H"
	  int iCorr = 0;
	  scalar initialResidual = 0;

	  do
	    {
                if (thermalStress)
		{
                    volScalarField& T = Tptr();
                    solve
                        (
                            fvm::ddt(T) == fvm::laplacian(DT, T)
                        );
		}

                {
                    fvVectorMatrix DEqn
                        (
                            fvm::d2dt2(D)
                            ==
                            fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
                            + divSigmaExp
                        );

                    if (thermalStress)
                    {
                        const volScalarField& T = Tptr();
                        DEqn += threeKalpha*fvc::grad(T);
                    }

                    //UEqn.setComponentReference(1, 0, vector::X, 0);
                    //UEqn.setComponentReference(1, 0, vector::Z, 0);

                    initialResidual = DEqn.solve().initialResidual();

                    if (!compactNormalStress)
                    {
                        divSigmaExp = fvc::div(DEqn.flux());
                    }
                }
              {
                  volTensorField gradD = fvc::grad(D);
                  sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);
                  
                  if (compactNormalStress)
                  {
                      divSigmaExp = 
                          fvc::div(sigmaD - (2*mu + lambda)*gradD, "div(sigmaD)");
                  }
                  else
                  {
                      divSigmaExp += fvc::div(sigmaD);
                  }
              }
              
	    } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

#       include "calculateStress.H"
 
	  Info << "\nMaximum Displacement: " << max(mag(D)).value() << endl;

	  
	} 

	Info << "\nCoupling the solutions\n" << endl;

	scalarField & fluidP = p.boundaryField()[fluidSideI];
	
 	scalarField &solidP = displace.pressure();
	
	forAll(fluidP,fI) {
	  solidP[exchange[fI]]=-fluidP[fI];
	}

#       include "write.H"
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
