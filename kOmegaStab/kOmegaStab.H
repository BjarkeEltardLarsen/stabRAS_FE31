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

Class
    Foam::incompressible::RASModels::kOmegaStab

Description
    Stabilized k-omega turbulence model for simulation of free-surface waves. Builds on the turbulence model by Wilcox (1988), includes a buoyancy production term see e.g. Umlauf et al. (2003)or Burchard (2002) and is stabilized for usage in "potential flow" regions according to Larsen and (Fuhrman 2018)

    References:
    @verbatim
        Wilcox, D. C.,
	"Turbulence Modeling for CFD"
        DCW Industries, Inc., La Canada,
        California, 1988.
	
	Burchard, H.
	"Applied Turbulence modelling for Marine Waters"
	Springer
	2002

     Umlauf, L., Burchard, H. and Hutter, K.
	"Extending the k-omega turbulence model towards oceanic applications"
	Ocean Modelling. 5(3), 195-218
	2003

	Larsen, B.E. and Fuhrman, D.R.
	On the over-production of turbulence beneath surface waves in RANS models
	J. Fluid Mech.
	2018 (accepted)
	
        See also:
        http://www.cfd-online.com/Wiki/Wilcox's_k-omega_model
    @endverbatim


    The default model coefficients correspond to the following:
    @verbatim
        kOmegaStabCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            beta        0.072;
            alphak      0.5;
            alphaOmega  0.5;
	    alphaBS     1.36;
	    lambda2     0.05;
	    
        }
    @endverbatim

SourceFiles
    kOmegaStab.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaStab_H
#define kOmegaStab_H

#include "RASModel.H"
#include <map>
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmegaStab Declaration
\*---------------------------------------------------------------------------*/

class kOmegaStab
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar beta_;
            dimensionedScalar alpha_;
            dimensionedScalar alphaK_;
            dimensionedScalar alphaOmega_;
            dimensionedScalar alphaBS_;
            dimensionedScalar lambda2_;
            dimensionedScalar pOmegaSmall_;

        // Fields

            volScalarField k_;
            volScalarField omega_;
            uniformDimensionedVectorField gField_;
            volScalarField nut_;


public:

    //- Runtime type information
    TypeName("kOmegaStab");

    // Constructors

        //- Construct from components
        kOmegaStab
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor
    virtual ~kOmegaStab()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaK_*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaOmega_*nut_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
