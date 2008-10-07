//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/createMeshes.H 1280 2007-04-27T13:55:25.796383Z bgschaid  $ 
    Info<< "region1: Create mesh for time = "
        << runTime.value() << nl << endl;

    fvMesh mesh1
    (
        IOobject
        (
            "region1",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

//    mesh1.constructAndClear();

    Info << "Mesh 1 with " << mesh1.nCells() << " cells \n" << endl;

    Info<< "region2: Create mesh for time = "
        << runTime.value() << nl << endl;

    fvMesh mesh2
    (
        IOobject
        (
            "region2",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

//    mesh2.constructAndClear();

    Info << "Mesh 2 with " << mesh2.nCells() << " cells \n" << endl;
