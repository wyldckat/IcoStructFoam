//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/createIcoFields.H 943 2007-01-25T21:38:27.500728Z bgschaid  $ 
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh1,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );


    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh1,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh1
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh1,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh1
    );


#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh1.solutionDict().subDict("PISO"), pRefCell, pRefValue);
