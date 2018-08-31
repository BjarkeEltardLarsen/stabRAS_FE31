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
    Foam::incompressible::RASModels::RNGkEpsilonStab

Description
 Stabilized renomalisation group k-epsilon turbulence model for simulation of free-surface waves. Builds on the turbulence model by Yakhot et al. (1991), includes a buoyancy production term see e.g. Umlauf et al. (2003) or Burchard (2002) and is stabilized for usage in "potential flow" regions according to Larsen and Fuhrman (2018).
References:
    @verbatim
        Yakhot, V., Thangham, S., Speziale, C., Orszag, S.A and Gatski, T. B.
	"Development of turbulence models for shear flows by a double expansion technique"
	Inst. Computer Apllication Sci. Eng. NAS1-18605
	1991
	
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
	
     
    @endverbatim

    The default model coefficients correspond to the following:
    @verbatim
        RNGkEpsilonStabCoeffs
        {
            Cmu         0.0845;
            C1          1.42;
            C2          1.68;
            sigmak      0.71942;
            sigmaEps    0.71942;
            eta0        4.38;
            beta        0.012;
alphaBS     1.36;
	    lambda2     0.05;
        }
    @endverbatim

SourceFiles
    RNGkEpsilonStab.C

\*---------------------------------------------------------------------------*/

#ifndef RNGkEpsilonStab_H
#define RNGkEpsilonStab_H

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
                           Class RNGkEpsilonStab Declaration
\*---------------------------------------------------------------------------*/

class RNGkEpsilonStab
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar eta0_;
            dimensionedScalar beta_;
 dimensionedScalar alphaBS_;
            dimensionedScalar lambda2_;
            dimensionedScalar pOmegaSmall_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
  uniformDimensionedVectorField gField_;
            volScalarField nut_;


public:

    //- Runtime type information
    TypeName("RNGkEpsilonStab");

    // Constructors

        //- Construct from components
        RNGkEpsilonStab
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~RNGkEpsilonStab()
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
                new volScalarField("DkEff", nut_/sigmak_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
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
