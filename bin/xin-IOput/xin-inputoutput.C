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
    xin-inputoutput

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    //The creates the fvmesh system for us(instances called mesh)
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    dictionary mydict;
    const word mydictname("mydictproperties");

    IOobject mydictIO
    (
        mydictname,
        mesh.time().constant(),//path to the file
        mesh,//reference to the  mesh it needed.
        IOobject::MUST_READ//read the file
    );

    if(!mydictIO.typeHeaderOk<dictionary>(true))
    {
        FatalErrorIn(args.executable())<<"cannot open the dictionary file "<<mydictname<<exit(FatalError);
    }

    mydict=IOdictionary(mydictIO);
    word someUsername;
    mydict.lookup("username")>>someUsername;
    scalar offsetvalue(mydict.lookupOrDefault<scalar>("offsetvalue",0.0));
    bool incompressibleflag(mydict.lookupOrDefault<Switch>("incompressibleFlag",true));
    word incom;
    if (incompressibleflag)  
    {
        incom="incompressible";
    }
    else
    {
        incom="compressible";
    } 
    List <scalar>inputvalues (mydict.lookup("inputvalues")); 
    HashTable<float,word>mytable;
    mytable=mydict.lookup("mytable");

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    fileName myoutputdirectory=mesh.time().path()/"postProcessing";
    mkDir(myoutputdirectory);

    autoPtr<OFstream> myoutputfileptr;
    myoutputfileptr.reset(new OFstream(myoutputdirectory/"myoutputfile.dat"));

    myoutputfileptr() << "run by xin" << endl;
    mytable.insert("force",2.0);
    myoutputfileptr() << mytable << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    Info<<"used by "<<someUsername<<nl
    <<"offsetvalue: "<<offsetvalue<<nl
    <<"you are running "<<incom<<" simulation."<<nl
    <<"with inputvalues: "<<inputvalues<<nl
    <<"with mytable: "<<mytable<<endl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
