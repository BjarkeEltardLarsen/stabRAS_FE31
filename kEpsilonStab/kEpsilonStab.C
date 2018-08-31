/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kEpsilonStab.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilonStab, 0);
addToRunTimeSelectionTable(RASModel, kEpsilonStab, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilonStab::kEpsilonStab
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
 alphaBS_ //Coefficient for the buoyancy production
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaBS",
            coeffDict_,
            1.36
        )
    ),
    lambda2_ //New limiter suggested by Larsen and Furhman (2018)
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lambda2",
            coeffDict_,
            0.05
        )
    ),
    pOmegaSmall_("pOmegaSmall", dimless/(dimTime*dimTime), SMALL),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),


    gField_ //Needed for the buoyancy production term
(
IOobject
(
"g",
this->db().time().constant(),
this->db(),
IOobject::MUST_READ,
IOobject::NO_WRITE
)
 ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
  nut_ = Cmu_*sqr(k_)/max(epsilon_ + epsilonSmall_,lambda2_*C2_/C1_*2.0*magSqr(symm(fvc::grad(U_)))/(2.0*magSqr(skew(fvc::grad(U_)))+pOmegaSmall_)*epsilon_); //nut  as suggested by Larsen and Fuhrman (2018)
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kEpsilonStab::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kEpsilonStab::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kEpsilonStab::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool kEpsilonStab::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
alphaBS_.readIfPresent(coeffDict());
 lambda2_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsilonStab::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(epsilon_, epsilon0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 // Include density and buoyancy production term. Implemented by Bjarke Eltard Larsen 18-07-2018

    // Load rho and rhoPhi
   const surfaceScalarField& rhoPhi_ = U_.db().objectRegistry::lookupObject<surfaceScalarField>("rho*phi");
    const volScalarField& rho_ = U_.db().objectRegistry::lookupObject<volScalarField>("rho");

    // Calculate the Brunt-Vaisala frequency
 volScalarField N2 = gField_&fvc::grad(rho_)/rho_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate eddy viscosity using the stabilizing approach

 volScalarField p0 = 2.0*magSqr(symm(fvc::grad(U_))); // 2S_ij S_ij
 volScalarField pOmega= 2.0*magSqr(skew(fvc::grad(U_))); // 2 Omega_ij Omega_ij
 volScalarField epsilonTilde = max(epsilon_ + epsilonSmall_,lambda2_*C2_/C1_*p0*epsilon_/(pOmega+pOmegaSmall_));

 nut_ = Cmu_*sqr(k_)/epsilonTilde; //nut as suggested by Larsen and Fuhrman (2018)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    volScalarField G("RASModel::G", nut_*2.0*magSqr(symm(fvc::grad(U_))));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
     fvm::ddt(rho_,epsilon_)
      + fvm::div(rhoPhi_, epsilon_)
      + fvm::SuSp(-fvc::div(rhoPhi_), epsilon_)
      - fvm::laplacian(rho_*DepsilonEff(), epsilon_)
     ==
        rho_*C1_*Cmu_*k_*p0 //Changed according to formulation in Larsen and Fuhrman (2018)
      - fvm::Sp(rho_*C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    // No longer needed: matrix completes at the point of solution
    // HJ, 17/Apr/2012
//     epsEqn().completeAssembly();

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
     fvm::ddt(rho_,k_)
      + fvm::div(rhoPhi_, k_)
      + fvm::SuSp(-fvc::div(rhoPhi_), k_)
      - fvm::laplacian(rho_*DkEff(), k_)
     ==
        G*rho_
      - fvm::Sp(rho_*epsilon_/k_, k_)
  - fvm::Sp(nut_*rho_*alphaBS_*N2/max(k_,k0_),k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
epsilonTilde = max(epsilon_ + epsilonSmall_,lambda2_*C2_/C1_*p0*epsilon_/(pOmega+pOmegaSmall_));

 nut_ = Cmu_*sqr(k_)/epsilonTilde;
  
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
