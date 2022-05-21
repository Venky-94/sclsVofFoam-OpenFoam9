#include "createRDeltaT.H"

//Start of sclsvof
    Info<< "Reading field psi0\n" << endl;
  

    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ); 

volScalarField psi0
    (
        IOobject
        (
            "psi0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        psi//mesh,
        //dimensionedScalar("psi0",dimless, 0.0)
    );

    Info<< "Reading field delta\n" << endl;
    volScalarField delta
    (
        IOobject
        (
            "delta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("delta",dimless, 0.0)
    );

    Info<< "Reading field H\n" << endl;
    volScalarField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("H",dimless, 0.0)
    );
    H.oldTime();

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("C",dimless/dimLength, 0.0)
    );
//End of sclsvof

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
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

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
userimmiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<usertwoPhaseChangeModel> phaseChangePtr
(
    usertwoPhaseChangeModel::New(mixture)
);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const volScalarField& nuMod1 = mixture.nuM1();
const volScalarField& nuMod2 = mixture.nuM2();

const volScalarField limAlpha1
    (
        "limAlpha1",
        min(max(alpha1, scalar(0)), scalar(1))
    );

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
    ),
	//rho2 + (rho1 - rho2)*(Foam::exp(25.00*limAlpha1)-1.00)/Foam::exp(25.00)
	//rho2 + (rho1 - rho2)*(Foam::exp(scalar(25)*limAlpha1)-scalar(1))/Foam::exp(scalar(25))
	//rho2 + (rho1 - rho2)*pow(limAlpha1, 1.5*(nuMod1/nuMod2) + 0.75)
	alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Getting dynamic viscocity for postprocessing
volScalarField muproc
(
	IOobject
	(
		"muproc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mixture.mu()
);

// Getting kinematic viscocity for postprocessing
volScalarField nuproc
(
	IOobject
	(
		"nuproc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mixture.nu()
);

// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

//Start of sclscvof
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
	
	dimensionedScalar deltaX
    (
        transportProperties.lookup("deltaX")
    );

    dimensionedScalar gamma
    (
        dimensionedScalar(deltaX*0.75)
    );

    dimensionedScalar epsilon
    (  //        dimensionedScalar(deltaX*1.5)
	 transportProperties.lookup("epsilon")
    );

    dimensionedScalar deltaTau
    (
        dimensionedScalar(deltaX*0.1)
    );

    dimensionedScalar dimChange
    (
        dimensionedScalar("dimChange",dimLength, 1.0)
    );

    dimensionedScalar sigma
    (
        transportProperties.lookup("sigma")
    );
//End of sclsvof