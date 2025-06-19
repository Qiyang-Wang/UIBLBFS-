/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "mydynamicOversetFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "label.H"
#include "vectorIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mydynamicOversetFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        mydynamicOversetFvMesh,
        IOobject
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mydynamicOversetFvMesh::mydynamicOversetFvMesh(const IOobject& io)
:
    dynamicMotionSolverFvMesh(io),
    faceType_
    (
        IOobject
        (
            "faceType",
            io.time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("myZero",dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    ),
    faceDonor_
    (
        IOobject
        (
            "faceDonor",
            io.time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("myZero",dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    ),
    cellType_
    (
        IOobject
        (
            "cellType",
            io.time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("myZero",dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    cellDonor_
    (
        IOobject
        (
            "cellDonor",
            io.time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("myZero",dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    volZoneID_
    (
        IOobject
        (
            "zoneID",
            io.time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    ),
    zoneID_
    (
        IOobject
        (
            "zoneIDList",
            facesInstance(),
            polyMesh::meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nCells()
    ),
    nZones_(gMax(volZoneID_)+1),
    meshParts_(nZones_)
//    cellMaps_(nZones_)
{
//    forAll(meshParts_, zonei)
//    {
//        fvMeshSubset meshI(*this);
//        meshI.setLargeCellSubset(static_cast<labelList>(zoneID_.internalField()), zonei);
//        //const labelList& cm = meshI.cellMap();
//        cellMaps_.set
//        (
//            zonei,
//            &(meshI.cellMap())
//        );
//        //Info << cm << endl;
//        //cin.get();
//
//        meshParts_.set
//        (
//            zonei,
//            &(meshI.subMesh())
//        );
//         
//        //const fvBoundaryMesh& fvm = meshParts_[zonei].boundary();
//        //forAll(fvm, patchI)
//        //{
//        //    Info << fvm[patchI].name() << endl;
//        //    Info << fvm[patchI].size() << endl;
//        //}
//        //cin.get();
//    }

    forAll(volZoneID_, cellI)
    {
        zoneID_[cellI] = label(volZoneID_[cellI]);
    }

    forAll(meshParts_, zonei)
    {
        meshParts_.set
        (
            zonei,
            new fvMeshSubset(*this)
        );
        meshParts_[zonei].setLargeCellSubset(zoneID_, zonei);
        //meshParts_[zonei].setLargeCellSubset(static_cast<labelList>(zoneID_.internalField()), zonei);
    }

    DynamicList<label> interpolationCells;
    DynamicList<label> holeCells;

    const fvBoundaryMesh& fvm = boundary();
    // *********************Mark the interpolation cells near the oversetPatch**************************//
    // 'overset' patches
    forAll(fvm, patchI)
    {
        //if (isA<zeroGradientFvPatchField>(fvm[patchI]))
        if (fvm[patchI].name()=="oversetPatch" || fvm[patchI].name()=="oversetpatch")
        {
            //Info << "hahahaha111111111111111111111111111111111111111111111111"<< fvm[patchI].name() << endl;
            const labelList& fc = fvm[patchI].faceCells();
            forAll(fc,fci)
            {
                interpolationCells.append(fc[fci]);
            }
        }
    }
    
    // *********************Mark the interpolation circle around the object**************************//
    //Determine zone that contains oversetPatch
    oversetZoneID_ = -1; 
    // 'overset' patches
    forAll(fvm, patchI)
    {
        //if (isA<oversetFvPatch>(fvm[patchI]))
        if (fvm[patchI].name()=="oversetPatch" || fvm[patchI].name()=="oversetpatch")
        {
            const labelList& fc = fvm[patchI].faceCells();
            label cell0 = fc[0];
            oversetZoneID_ = zoneID_[cell0]; 
        }
    }

    //define a radius 
    //if some points of a cell are wihtin the radius, and some are not, then the cell defines the front
    vectorIOList test = time().lookupObject<vectorIOList>("oversetMotion");
    const point center(test[0].x(), test[0].y(), test[0].z());
    // oscillation cylinder
    //scalar t = mesh_.time().value();
    //point center(0, 0.14*sin(1.0472*t), 0);
    //scalar r_ = 1.0;
    
    // swming fish
    rxMin_ = -0.75;
    rxMax_ = 0.75;
    ryMin_ = -0.75;
    ryMax_ = 0.75;
    
    r_ = 0.75;
    
   // sedimentation
   // rxMin_ = -1.4;
   // rxMax_ = 1.4;
   // ryMin_ = -1.4;
   // ryMax_ = 5.;

    const labelListList& cellPts = cellPoints();
    const vectorList& pts = points();
    //const pointList& pts = points();
    for(label cellI = 0; cellI < nCells(); cellI++)
    {
        labelList pointSet = cellPts[cellI];
        point point0 = pts[pointSet[0]];
        vector r0 = point0 - center;
        bool inOut0 = mag(r0)<r_;
        // bool inOut0 = (r0.x()<rxMax_)
		//     &&(r0.x()>rxMin_)
		//     &&(r0.y()<ryMax_)
		//     &&(r0.y()>ryMin_);
        label n = 1;
        for(n = 1; n < pointSet.size(); n++)
        {
            point pointN = pts[pointSet[n]];
            vector rN = pointN-center;
            bool inOutN = mag(rN)<r_;
            // bool inOutN = (rN.x()<rxMax_)
		    //     &&(rN.x()>rxMin_)
		    //     &&(rN.y()<ryMax_)
		    //     &&(rN.y()>ryMin_);
            if(inOut0!=inOutN)
            {
                if(zoneID_[cellI]!=oversetZoneID_)
                {
                    interpolationCells.append(cellI);
                    break;
                }
            }
        }
        if((inOut0)&&(n == pointSet.size()-1))
        {
            if(zoneID_[cellI]!=oversetZoneID_)
            {
                holeCells.append(cellI);
            }
        }
    }

    interpolationCells_.transfer(interpolationCells);
    holeCells_.transfer(holeCells);

    //reset interpolation cell 
    for(label cellI = 0; cellI < nCells(); cellI++)
    {
        cellType_[cellI] = 0;
    }

    forAll(interpolationCells_, i)
    {
        cellType_[interpolationCells_[i]] = 1;
    }

    forAll(holeCells_, i)
    {
        cellType_[holeCells_[i]] = 2;
    }
    
    walkFront();
    findDonor();
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mydynamicOversetFvMesh::~mydynamicOversetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mydynamicOversetFvMesh::update()
{
    if(Foam::dynamicMotionSolverFvMesh::update())
    {
        forAll(meshParts_, zonei)
        {
            meshParts_.set
            (
                zonei,
                new fvMeshSubset(*this)
            );
            meshParts_[zonei].setLargeCellSubset(zoneID_, zonei);
        }

        updateInterpolationCell();
        //interpolateField();
    }

    return true;
}

bool Foam::mydynamicOversetFvMesh::updateInterpolationCell()
{
    DynamicList<label> interpolationCells;
    DynamicList<label> holeCells;

    const fvBoundaryMesh& fvm = boundary();
    // *********************Mark the interpolation cells near the oversetPatch**************************//
    // 'overset' patches
    forAll(fvm, patchI)
    {
        //if (isA<zeroGradientFvPatchField>(fvm[patchI]))
        if (fvm[patchI].name()=="oversetPatch" || fvm[patchI].name()=="oversetpatch")
        {
            //Info << "hahahaha111111111111111111111111111111111111111111111111"<< fvm[patchI].name() << endl;
            const labelList& fc = fvm[patchI].faceCells();
            forAll(fc,fci)
            {
                interpolationCells.append(fc[fci]);
            }
        }
    }
    
    // *********************Mark the interpolation circle around the object**************************//
    //define a radius 
    //if some points of a cell are wihtin the radius, and some are not, then the cell defines the front
    //scalar t = time().value();
    //point center(0, 0.14*sin(1.0472*t), 0);
    vectorIOList test = time().lookupObject<vectorIOList>("oversetMotion");
    const point center(test[0].x(), test[0].y(), test[0].z());
    //point center(0, 0, 0);

    const labelListList& cellPts = cellPoints();
    const vectorList& pts = points();
    //const pointList& pts = points();
    for(label cellI = 0; cellI < nCells(); cellI++)
    {
        labelList pointSet = cellPts[cellI];
        point point0 = pts[pointSet[0]];
        vector r0 = point0 - center;
        bool inOut0 = mag(r0)<r_;
        // bool inOut0 = (r0.x()<rxMax_)
		//     &&(r0.x()>rxMin_)
		//     &&(r0.y()<ryMax_)
		//     &&(r0.y()>ryMin_);
        label n = 1;
        for(n = 1; n < pointSet.size(); n++)
        {
            point pointN = pts[pointSet[n]];
            vector rN = pointN-center;
            bool inOutN = mag(rN)<r_;
            // bool inOutN = (rN.x()<rxMax_)
		    //     &&(rN.x()>rxMin_)
		    //     &&(rN.y()<ryMax_)
		    //     &&(rN.y()>ryMin_);
            if(inOut0!=inOutN)
            {
                if(zoneID_[cellI]!=oversetZoneID_)
                {
                    interpolationCells.append(cellI);
                    break;
                }
            }
            if((inOut0)&&(n == pointSet.size()-1))
            {
                if(zoneID_[cellI]!=oversetZoneID_)
                {
                    interpolationCells.append(cellI);
                    //holeCells.append(cellI);
                }
            }
        }
    }

    interpolationCells_.clear();
    interpolationCells_.transfer(interpolationCells);
    holeCells_.clear();
    holeCells_.transfer(holeCells);
    
    //reset interpolation cell 
    for(label cellI = 0; cellI < nCells(); cellI++)
    {
        cellType_[cellI] = 0;
    }

    forAll(interpolationCells_, i)
    {
        //Info << interpolationCells_[i] << endl;
        cellType_[interpolationCells_[i]] = 1;
    }
    forAll(holeCells_, i)
    {
        //Info << interpolationCells_[i] << endl;
        cellType_[holeCells_[i]] = 2;
    }

    walkFront();
    findDonor();

    return true;
}

bool Foam::mydynamicOversetFvMesh::interpolateField()
{
    
    const volScalarField &rho = lookupObject<volScalarField>("rho");
    const volVectorField &rhoU = lookupObject<volVectorField>("rhoU");
    const volScalarField &rhoE = lookupObject<volScalarField>("rhoE");
     
    //forAll(interpolationCells_, i)
    //{
    //    //Info << interpolationCells_[i] << endl;
    //    rho[interpolationCells_[i]] = rho[donorCells_[i]];
    //    rhoU[interpolationCells_[i]] = rhoU[donorCells_[i]];
    //    rhoE[interpolationCells_[i]] = rhoE[donorCells_[i]];
    //}

    //Dearl with holes

    return true;
}

bool Foam::mydynamicOversetFvMesh::findDonor()
{
    //DynamicList<label> donorCells;

    //forAll(interpolationCells_, i)
    //{
    //    point location = C()[interpolationCells_[i]];
    //    label zonei = zoneID_[interpolationCells_[i]];
    //    if(zonei == 0)
    //    {
    //        zonei=1;        
    //    }else if(zonei ==1) {
    //        zonei=0;        
    //    }

    //    label donor = findNearestCell(location, 0, zonei);
    //    donorCells.append(meshParts_[zonei].cellMap()[donor]);
    //}

    //donorCells_.clear();
    //donorCells_.transfer(donorCells);
    
    for(label cellI = 0; cellI < nCells(); cellI++)
    {
        if(cellType_[cellI] == 1)
        {
            point location = C()[cellI];
            label zonei = zoneID_[cellI];
            if(zonei == 0)
            {
                zonei=1;        
            }else if(zonei == 1) {
                zonei=0;        
            }
            label donor = findNearestCell(location, 0, zonei);
	    cellDonor_[cellI] = meshParts_[zonei].cellMap()[donor];
	} else {
	    cellDonor_[cellI] = -1;
	}
    }

   /////////////////////////////////////// 
    label zonei = 1;
    for(label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        if (faceType_[faceI] == 1) 
        {
            label donor = findNearestCell(Cf()[faceI], 0, zonei);
            faceDonor_[faceI] = meshParts_[zonei].cellMap()[donor];
            cellType_[faceDonor_[faceI]] = 3;
        } else {
            faceDonor_[faceI] = -1;
        }
    }

    const fvBoundaryMesh& fvm = boundary();
    zonei = 0;
    forAll(fvm, patchI)
    {
        forAll(fvm[patchI],fci)
        {
           if(faceType_.boundaryFieldRef()[patchI][fci] == 1)
           {
               label donor = findNearestCell(fvm[patchI].Cf()[fci], 0, zonei);
               faceDonor_.boundaryFieldRef()[patchI][fci] = meshParts_[zonei].cellMap()[donor];
               cellType_[faceDonor_.boundaryFieldRef()[patchI][fci]] = 3;
           } else {
               faceDonor_.boundaryFieldRef()[patchI][fci] = -1;
           }
        }
    }

    return true;
}

Foam::label Foam::mydynamicOversetFvMesh::findNearestCell
(
    const point& location,
    const label seedCelli,
    const label zonei
) const
{
//    const fvMesh& mesh = *this;
    const fvMesh& mesh = meshParts_[zonei].subMesh();
    label curCelli = seedCelli;
    scalar distanceSqr = magSqr(mesh.cellCentres()[curCelli] - location);
 
    bool closer;
 
    do
    {
        // Try neighbours of curCelli
        closer = findNearer
        (
            location,
            mesh.cellCentres(),
            mesh.cellCells()[curCelli],
            curCelli,
            distanceSqr
        );
    } while (closer);
 
    return curCelli;
}

bool Foam::mydynamicOversetFvMesh::findNearer
(
    const point& sample,
    const pointField& points,
    const labelList& indices,
    label& nearestI,
    scalar& nearestDistSqr
) const
{
    bool nearer = false;

    forAll(indices, i)
    {
        label pointi = indices[i];

        scalar distSqr = magSqr(points[pointi] - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            nearestI = pointi;
            nearer = true;
        }
    }

    return nearer;
}

bool Foam::mydynamicOversetFvMesh::walkFront()
{
    DynamicList<label> frontList;

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    for (label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        label ownType = cellType_[own[faceI]];
        label neiType = cellType_[nei[faceI]];

        if
        (
             (ownType == 2 && neiType != 2)
          || (ownType != 2 && neiType == 2)
        )
        {
            //Pout<< "Front at face:" << faceI
            //    << " at:" << mesh_.faceCentres()[faceI] << endl;
            frontList.append(faceI);
        }
    }
    
    frontList_.clear();
    frontList_.transfer(frontList);
    
    //reset interpolation cell 
    for(label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        faceType_[faceI] = 0;
    }

    forAll(frontList_, i)
    {
        //Info << frontList_[i] << endl;
        faceType_[frontList_[i]] = 1;
    }
    
    const fvBoundaryMesh& fvm = boundary();
    forAll(fvm, patchI)
    {
        //if (isA<zeroGradientFvPatchField>(fvm[patchI]))
        if (fvm[patchI].name()=="oversetPatch" || fvm[patchI].name()=="oversetpatch")
        {
            forAll(fvm[patchI],fci)
            {
                faceType_.boundaryFieldRef()[patchI][fci] = 1;
            }
        }else {
            forAll(fvm[patchI],fci)
            {
                faceType_.boundaryFieldRef()[patchI][fci] = 0;
            }
        }
    }
 
    return true;
}


// ************************************************************************* //
