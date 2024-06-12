/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "xinNewmark.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{
    defineTypeNameAndDebug(xinNewmark, 0);
    addToRunTimeSelectionTable(sixDoFSolver, xinNewmark, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::xinNewmark::xinNewmark
(
    const dictionary& dict,
    sixDoFRigidBodyMotion& body
)
:
    sixDoFSolver(dict, body),
    gamma_(dict.getOrDefault<scalar>("gamma", 0.5)),
    beta_
    (
        max
        (
            0.25*sqr(gamma_ + 0.5),
            dict.getOrDefault<scalar>("beta", 0.25)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::xinNewmark::~xinNewmark()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDoFSolvers::xinNewmark::solve
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    //pi: angular momentum, tau: torque, tconstraints: translational constraints, rconstraints: rotational constraints, aDamp: acceleration damping
    //v: velocity, a: acceleration, centreOfRotation: centre of rotation, Q: orientation, 
    Info<< "xinNewmark::solve print start ------------------------------------" << endl;
    Info << "gamma = " << gamma_ << ", beta = " << beta_ << endl;
    Info << "deltaT = " << deltaT << ", deltaT0 = " << deltaT0 << endl;
    Info << "fGlobal = " << fGlobal << ", tauGlobal = " << tauGlobal << endl;
    Info << "v0 = " << v0() << ", a0 = " << a0() << endl;
    Info <<"v = " << v() << ", a = " << a() << endl;
    Info << "pi0 = " << pi0() << ", tau0 = " << tau0() << endl;
    Info << "pi = " << pi() << ", tau = " << tau() << endl;
    Info << "centreOfRotation0 = " << centreOfRotation0() << endl;
    Info << "centreOfRotation = " << centreOfRotation() << endl;
    Info << "tConstraints = " << tConstraints() << endl;
    Info << "rConstraints = " << rConstraints() << endl;
    Info << "aDamp = " << aDamp() << endl;
    Info << "Q0 = " << Q0() << ", Q = " << Q() << endl;
    //access mass 
    Info<<"mass = "<<body_.mass()<<endl;
    //access stiffness
    dynamicFvMesh()
    //access damping
    Info<< "xinNewmark::solve print end ------------------------------------" << endl;

    // Update the linear acceleration and torque
    updateAcceleration(fGlobal, tauGlobal);
    // Update the constraints to the object
    updateConstraints();

    // Correct linear velocity
    v() =
        tConstraints()
      & (v0() + aDamp()*deltaT*(gamma_*a() + (1 - gamma_)*a0()));

    // Correct angular momentum
    pi() =
        rConstraints()
      & (pi0() + aDamp()*deltaT*(gamma_*tau() + (1 - gamma_)*tau0()));

    // Correct position
    centreOfRotation() =
        centreOfRotation0()
      + (
            tConstraints()
          & (
                deltaT*v0()
              + aDamp()*sqr(deltaT)*(beta_*a() + (0.5 - beta_)*a0())
            )
        );

    // Correct orientation
    vector piDeltaT =
        rConstraints()
      & (
            deltaT*pi0()
          + aDamp()*sqr(deltaT)*(beta_*tau() + (0.5 - beta_)*tau0())
        );
    Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
    Q() = Qpi.first();

    Info<< "xinNewmark::solve print start2 ------------------------------------" << endl;
    Info << "piDeltaT = " << piDeltaT << endl;
    Info << "Qpi.first = " << Qpi.first() <<endl;
    Info<< "xinNewmark::solve print end2 ------------------------------------" << endl;
}

// ************************************************************************* //
