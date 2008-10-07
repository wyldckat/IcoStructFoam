//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/readMechanicalProperties.H 1280 2007-04-27T13:55:25.796383Z bgschaid  $ 
    Info<< "Reading mechanical properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh2,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho(mechanicalProperties.lookup("rho"));
    dimensionedScalar rhoE(mechanicalProperties.lookup("E"));
    dimensionedScalar nuMech(mechanicalProperties.lookup("nu"));

    Info<< "Normalising E : E/rho\n" << endl;
    dimensionedScalar E = rhoE/rho;

    Info<< "Calculating Lame's coefficients\n" << endl;

    dimensionedScalar mu = E/(2.0*(1.0 + nuMech));
    dimensionedScalar lambda = nuMech*E/((1.0 + nuMech)*(1.0 - 2.0*nuMech));
    dimensionedScalar threeK = E/(1.0 - 2.0*nuMech);

    Switch planeStress(mechanicalProperties.lookup("planeStress"));

    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        //- change lambda and threeK for plane stress
        lambda = nuMech*E/((1.0 + nuMech)*(1.0 - nuMech));
        threeK = E/(1.0 - nuMech);
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }

    Info<< "mu = " << mu.value() << " Pa/rho\n";
    Info<< "lambda = " << lambda.value() << " Pa/rho\n";
    Info<< "threeK = " << threeK.value() << " Pa/rho\n";
