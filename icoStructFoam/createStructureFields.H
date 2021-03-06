//  ICE Revision: $Id: /local/openfoam/branches/WikiVersions/FluidStructCoupling/icoStructFoam/createStructureFields.H 1280 2007-04-27T13:55:25.796383Z bgschaid  $ 

// mesh2.weights();

Info<< "Reading field D (displacemnent)\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh2,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh2
);


autoPtr<volScalarField> Tptr(NULL);

if (thermalStress)
{
    Info<< "Reading field T\n" << endl;
    Tptr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh2,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh2
        )
    );
}


Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh2,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(fvc::grad(D)) + lambda*(I*tr(fvc::grad(D)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh2,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}
