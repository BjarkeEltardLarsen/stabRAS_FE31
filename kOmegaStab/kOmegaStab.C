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

#include "kOmegaStab.H"
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

defineTypeNameAndDebug(kOmegaStab, 0);
addToRunTimeSelectionTable(RASModel, kOmegaStab, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaStab::kOmegaStab
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
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            coeffDict_,
            0.072
        )
    ),
    alpha_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha",
            coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            coeffDict_,
            0.5
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
    lambda2_ //New limiter suggested by Larsen and Fuhrman 2018
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
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
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
  nut_ = k_/max(omega_ + omegaSmall_,lambda2_*beta_/(Cmu_*alpha_)*2.0*magSqr(symm(fvc::grad(U_)))/(2.0*magSqr(skew(fvc::grad(U_)))+pOmegaSmall_)*omega_); //nut as suggested by Larsen and Fuhrman 2018
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaStab::R() const
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


tmp<volSymmTensorField> kOmegaStab::devReff() const
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


tmp<fvVectorMatrix> kOmegaStab::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool kOmegaStab::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        beta_.readIfPresent(coeffDict());
        alphaK_.readIfPresent(coeffDict());
        alphaOmega_.readIfPresent(coeffDict());
 alphaBS_.readIfPresent(coeffDict());
 lambda2_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaStab::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(omega_, omega0_);
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
// Calculate eddy viscosity using the stabilizing approach described in Larsen and Fuhrman (2018)

 volScalarField p0 = 2.0*magSqr(symm(fvc::grad(U_))); // 2S_ij S_ij
 volScalarField pOmega= 2.0*magSqr(skew(fvc::grad(U_))); // 2 Omega_ij Omega_ij
 volScalarField omegaTilde2 = max (omega_+omegaSmall_, lambda2_*(beta_/(Cmu_*alpha_))*p0*omega_/(pOmega+pOmegaSmall_));
  nut_ = k_/omegaTilde2;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

   volScalarField G("RASModel::G", nut_*2.0*magSqr(symm(fvc::grad(U_))));

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
     fvm::ddt(rho_,omega_)
      + fvm::div(rhoPhi_, omega_)
      + fvm::SuSp(-fvc::div(rhoPhi_), omega_)
      - fvm::laplacian(rho_*DomegaEff(), omega_)
     ==
     rho_*alpha_*p0 //Changed according to formulation in Larsen and Fuhrman (2018)
     - fvm::Sp(beta_*omega_*rho_, omega_)
    );

    omegaEqn().relax();

    // No longer needed: matrix completes at the point of solution
    // HJ, 17/Apr/2012
//     omegaEqn().completeAssembly();

    solve(omegaEqn);
    bound(omega_, omega0_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
     fvm::ddt(rho_,k_)
      + fvm::div(rhoPhi_, k_)
      + fvm::SuSp(-fvc::div(rhoPhi_), k_)
      - fvm::laplacian(rho_*DkEff(), k_)
     ==
    G*rho_
      - fvm::Sp(Cmu_*omega_*rho_, k_)
  - fvm::Sp(nut_*rho_*alphaBS_*N2/max(k_,k0_),k_) // Buoyancy production term
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
 omegaTilde2 = max (omega_+omegaSmall_, lambda2_*(beta_/(Cmu_*alpha_))*p0*omega_/(pOmega+pOmegaSmall_));
  nut_ = k_/omegaTilde2;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
