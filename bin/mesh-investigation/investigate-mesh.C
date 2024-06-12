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
    investigate-mesh

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // #include "labelFwd.H"
    // #define WM_LABEL_SIZE 32
    // #include "label.H"

    Info << "The most recent time folder is " << runTime.timeName() << nl
    //mesh.C() and mesh.Cf() are the cells and faces center coordinates
    <<" The mesh has " <<mesh.C().size() << " cells and " << mesh.Cf().size()
    << " internal faces." << endl;

    Info <<"#############################################################################"<<endl;
    word patchname;
    forAll(mesh.boundary(),patchI)
    {
        scalar surfaceArea = 0.0;
        forAll(mesh.boundary()[patchI].Cf(),myfaceI)
        {
            surfaceArea += mag(mesh.boundary()[patchI].Sf()[myfaceI]);
        }
        patchname=mesh.boundary()[patchI].name();
        Info<<"The surface area of the patch "<<patchname<<" "<<patchI<<" is "<<surfaceArea<<endl;
    }

    // Info<<"Found patch : "<<patchname <<" at the ID "<<patchID<<nl<<endl;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
