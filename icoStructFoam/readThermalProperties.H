//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/readThermalProperties.H 943 2007-01-25T21:38:27.500728Z bgschaid  $ 
    Info<< "Reading thermal properties\n" << endl;

    IOdictionary thermalProperties
    (
        IOobject
        (
            "thermalProperties",
            runTime.constant(),
            mesh2,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Switch thermalStress(thermalProperties.lookup("thermalStress"));

    dimensionedScalar threeKalpha
    (
        "threeKalpha",
        dimensionSet(0, 2, -2 , -1, 0),
        0
    );

    dimensionedScalar DT
    (
        "DT",
        dimensionSet(0, 2, -1 , 0, 0),
        0
    );

    if (thermalStress)
    {
        dimensionedScalar C(thermalProperties.lookup("C"));
        dimensionedScalar rhoK(thermalProperties.lookup("k"));
        dimensionedScalar alpha(thermalProperties.lookup("alpha"));

        Info<< "Normalising k : k/rho\n" << endl;
        dimensionedScalar k = rhoK/rho;

        Info<< "Calculating thermal coefficients\n" << endl;

        threeKalpha = threeK*alpha;
        DT.value() = (k/C).value();

        Info<< "threeKalpha = " << threeKalpha.value() << " Pa/rho\n";
    }
