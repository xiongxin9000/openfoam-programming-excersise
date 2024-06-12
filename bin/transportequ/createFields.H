// create T and reaf from initial and boundary conditions
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// -------从U文件中读取速度场---------
Info << "read velocity" << endl;
// 定义一个向量场U
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//reading transportProperties
Info << "Reading transportProperties" << endl;
// 读取transportProperties文件
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// 读取热扩散系数
Info << "Reading DT" << endl;
dimensionedScalar DT (transportProperties.lookupOrDefault<dimensionedScalar>("DT", 0.01));
// dimensionedScalar DT (transportProperties.lookup("DT"));


//reading/creating phi
Info << "Reading phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf() //& is the overloaded dot product
);



