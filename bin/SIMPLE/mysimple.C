/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
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

Application
    mysimple

Description
This is light way of simple algorithm.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar alpha;
    fvSolution.lookup("alpha") >> alpha;
    scalar pRefCell;
    fvSolution.lookup("pRefCell") >> pRefCell;
    scalar pRefValue;
    fvSolution.lookup("pRefValue") >> pRefValue;

    Info<<"read the following parameters: "<<endl;
    Info<<"alpha = "<<alpha<<endl;
    Info<<"pRefCell = "<<pRefCell<<endl;
    Info<<"pRefValue = "<<pRefValue<<endl;

    while (runTime.loop())
    {
        Info<<nl<<"time = "<<runTime.timeName()<<endl;

        fvVectorMatrix Ueqn
        (
          fvm::div(phi, U)- fvm::laplacian(nu, U)==-fvc::grad(p)
        );


        //1. solve the momentum equation
        Ueqn.solve();

        volScalarField A=Ueqn.A();
        volVectorField H=Ueqn.H();

        volScalarField A_inv=1.0/A;
        surfaceScalarField A_inv_f= fvc::interpolate(A_inv);

        volVectorField HbyA=A_inv*H;

        fvScalarMatrix pEqn
        (
            
            //AU=H-grad(p)
            //U=H/A-(1/A)*nabla(p)
            //nabla(U)=0

            //pressure equation
            //nabla((1/A)*nabla(p))=nabla(H/A)
            fvm::laplacian(A_inv_f, p)==fvc::div(HbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
        //2. solve the pressure equation
        pEqn.solve();

        p=alpha*p+(1.0-alpha)*p_old;
        U=(A_inv*H)-(A_inv*fvc::grad(p));

        phi= fvc::interpolate(U) & mesh.Sf();

        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        p_old=p;

        runTime.write();

    }
    


    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
