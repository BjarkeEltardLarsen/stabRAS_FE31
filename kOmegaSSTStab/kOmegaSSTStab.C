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

#include "kOmegaSSTStab.H"
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

defineTypeNameAndDebug(kOmegaSSTStab, 0);
addToRunTimeSelectionTable(RASModel, kOmegaSSTStab, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kOmegaSSTStab::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    volScalarField arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

tmp<volScalarField> kOmegaSSTStab::F2() const
{
    volScalarField arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSSTStab::kOmegaSSTStab
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
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

    y_(mesh_),

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
        autoCreateK("k", mesh_, U_.db())
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
        autoCreateOmega("omega", mesh_, U_.db())
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
        autoCreateNut("nut", mesh_, U_.db())
    )
{
    bound(omega_, omega0_);

    nut_ = a1_*k_/max(a1_*lambda2_*beta1_/(betaStar_*gamma1_)*2.0*magSqr(symm(fvc::grad(U_)))/(2.0*magSqr(skew(fvc::grad(U_)))+pOmegaSmall_)*omega_,max(a1_*omega_, F2()*sqrt(2.0)*mag(symm(fvc::grad(U_))))); //nut as suggested by Larsen and Fuhrman 2018
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaSSTStab::R() const
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


tmp<volSymmTensorField> kOmegaSSTStab::devReff() const
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


tmp<fvVectorMatrix> kOmegaSSTStab::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool kOmegaSSTStab::read()
{
    if (RASModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
	alphaBS_.readIfPresent(coeffDict());
	lambda2_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaSSTStab::correct()
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

    if (mesh_.changing())
    {
        y_.correct();
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
 volScalarField omegaTilde2 = max(omega_,lambda2_*beta1_/(betaStar_*gamma1_)*p0*omega_/(pOmega+pOmegaSmall_));

nut_ = a1_*k_/max(a1_*omegaTilde2,max(a1_*omega_, F2()*sqrt(2.0)*mag(symm(fvc::grad(U_))))); //nut as suggested by Larsen and Fuhrman 2018



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


volScalarField S2 = magSqr(symm(fvc::grad(U_))); 
  volScalarField G("RASModel::G", nut_*2.0*S2);

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    volScalarField CDkOmega =
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_;

    volScalarField F1 = this->F1(CDkOmega);

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
     fvm::ddt(rho_,omega_)
      + fvm::div(rhoPhi_, omega_)
      + fvm::SuSp(-fvc::div(rhoPhi_), omega_)
      - fvm::laplacian(rho_*DomegaEff(F1), omega_)
     ==
        rho_*gamma(F1)*p0 
      - fvm::Sp(beta(F1)*omega_*rho_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_*rho_,
            omega_
        )
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
      - fvm::laplacian(rho_*DkEff(F1), k_)
     ==
        min(G*rho_, rho_*c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_*rho_, k_)
 - fvm::Sp(nut_*rho_*alphaBS_*N2/max(k_,k0_),k_) // Buoyancy production term
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
omegaTilde2 = max(omega_,lambda2_*beta1_/(betaStar_*gamma1_)*p0*omega_/(pOmega+pOmegaSmall_));

nut_ = a1_*k_/max(a1_*omegaTilde2,max(a1_*omega_, F2()*sqrt(2.0)*mag(symm(fvc::grad(U_))))); //nut calculated based as suggested by Larsen and Fuhrman (2018)
 
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
