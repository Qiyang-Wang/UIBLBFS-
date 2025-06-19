/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 1991-2008 OpenCFD Ltd.

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gksFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "fvcSurfaceReconstruct.H"
#include "cellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(gksFluxScheme, 0);
addToRunTimeSelectionTable(fluxScheme, gksFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gksFluxScheme::gksFluxScheme
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const compressibleTurbulenceModel& turbulence
)
:
    fluxScheme(typeName, dict),
    mesh_(U.mesh()),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    rhoU_(rhoU),
    rhoE_(rhoE),
    dict_(dict),
    turbulence_(turbulence)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

gksFluxScheme::~gksFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::gksFluxScheme::calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{
/*
        volScalarField rhoOld
        (
	    IOobject
	    (
	      "rhoOld",
	      "0",
	      mesh_,
	      IOobject::NO_READ,
	      IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimensionedScalar("rhoOld",dimensionSet(1,-3,0,0,0,0,0),scalar(1.225))
        );
*/
	surfaceScalarField slen
	(
	    IOobject
	    (
		"slen",
	        "0",	
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimensionedScalar( "slen", dimensionSet(0,1,0,0,0,0,0), 0.0 )
	);

    calcStreamingDistant(slen);
    //dimensionedScalar Cspeed("Cspeed", dimensionSet (0 ,1, -1, 0, 0, 0, 0), Foam::sqrt(3.0*287*288));
    scalar Cspeed(Foam::sqrt(3.0*287*288.15));

    volScalarField nuEff(turbulence_.muEff()/rho_);

    surfaceScalarField  muave  = fvc::interpolate(nuEff);//mu at cell faces

//    Info << muave << endl;
//    Info <<"write muave" << endl;

    volVectorField gradrho = fvc::grad(rho_);
    volTensorField gradU   = fvc::grad(U_);

    //Info <<"builing fluxes...."<<endl;
    // Variables definition
    label  id_L, id_R ;                      // Mesh faces, cells connectivity
    vector CG_L, CG_R, Cf, n_i_;               //
    scalar rho_L_,  rho_R_ ;         // Temporary conservative variables L, R
    vector m_L_, m_R_;                                       //                                  //
    vector m_B_;                                             //
    scalar Frho_;                               // Temporary fluxes
    scalar Fe_;                               // Temporary fluxes
    vector Fm_;                                              //
    vector U_L_, U_R_;                                       //
    tensor gradU_L,gradU_R;
    vector gradrho_L,gradrho_R;
    scalar Sf;
    //==========================================================================
    // A) Loop on internal faces [index (i)]
    //==========================================================================

    const faceList& fs = mesh_.faces();
    const pointField& p = mesh_.points();

    forAll( mesh_.Sf(), i )
    {
        id_L        = mesh_.faceOwner()[i];
        id_R        = mesh_.faceNeighbour()[i];
        rho_L_      = rho_[id_L];
        rho_R_      = rho_[id_R];
        m_L_        = rhoU_[id_L];
        m_R_        = rhoU_[id_R];
        U_L_ = m_L_/rho_L_;
        U_R_ = m_R_/rho_R_;

        const labelList& f = fs[i];
        vector t1 = (p[f[1]] - p[f[0]]);
        // n_i_        = faceNormal[i];
        n_i_        = mesh_.Sf()[i]/mesh_.magSf()[i];

        vector t2 = n_i_^t1;
        vector t1_i_ = t1/mag(t1);
        vector t2_i_ = t2/mag(t2);

        // Cells CGs
        CG_L  = mesh_.C()[id_L];
        CG_R  = mesh_.C()[id_R];
        Cf    = mesh_.Cf()[i];

        gradrho_L = gradrho[id_L];
        gradrho_R = gradrho[id_R];
        gradU_L   =   gradU[id_L];
        gradU_R   =   gradU[id_R];
/*
        Info <<" rho_L_   " << endl;
        Info << rho_L_      << endl;
        Info <<" rho_R_   " << endl;
        Info << rho_R_      << endl;
        Info <<" U_L_     " << endl;
        Info << U_L_        << endl;
        Info <<" U_R_     " << endl;
        Info << U_R_        << endl;
        Info <<" gradU_L  " << endl;
        Info << gradU_L     << endl;
        Info <<" gradU_R  " << endl;
        Info << gradU_R     << endl;
        Info <<" gradrho_L" << endl;
        Info << gradrho_L   << endl;
        Info <<" gradrho_R" << endl;
        Info << gradrho_R   << endl;
        Info <<" CG_L     " << endl;
        Info << CG_L        << endl;
        Info <<" CG_R     " << endl;
        Info << CG_R        << endl;
        Info <<" Cf       " << endl;
        Info << Cf          << endl;
        Info <<" n_i_     " << endl;
        Info << n_i_        << endl;
        Info <<" t1_i_    " << endl;
        Info << t1_i_       << endl;
        Info <<" t2_i_    " << endl;
        Info << t2_i_       << endl;
        Info <<" muave[i] " << endl;
        Info << muave[i]    << endl;
        Info <<" slen[i]  " << endl;
        Info << slen[i]     << endl;
        Info <<" Cspeed   " << endl;
        Info << Cspeed      << endl;
        Info <<" Frho_    " << endl;
        Info << Frho_       << endl;
        Info <<" Fm_      " << endl;
        Info << Fm_        << endl;

        cin.get();
*/
        icoGKSFlux(  rho_L_ ,
                     rho_R_ ,
                     U_L_ ,
                     U_R_ ,
                     gradU_L,
                     gradU_R,
                     gradrho_L ,
                     gradrho_R,
                     CG_L,
                     CG_R,
                     Cf,
                     n_i_,
                     t1_i_,
                     t2_i_,
                     muave[i],
                     slen[i] ,
                     Cspeed,
                     Frho_,
                     Fm_,
                     Fe_
                  );

        // Output
        phi[i]       = Frho_*mesh_.magSf()[i];
        phiUp[i]         = Fm_*mesh_.magSf()[i];
        
        phiEp[i] = Fe_*mesh_.magSf()[i];
    }
    //==========================================================================
    // B) Loop on boundary faces and BCs [local index(ii), global index(i)]
    //==========================================================================
     forAll( mesh_.boundaryMesh(), iPatch )
    { // Patch BCtype
        word BCTypePhysical = mesh_.boundaryMesh().physicalTypes()[iPatch];
        word BCType         = mesh_.boundaryMesh().types()[iPatch];
        word BCName         = mesh_.boundaryMesh().names()[iPatch];
        //======================================================================
        // Loop on patch BCtype faces
        //======================================================================
        // *** BCs ( Boundary type must not be empty ) ***
       // if((BCType== "wall")&&(FFrho.boundaryField()[iPatch].size() > 0))

        if(phi.boundaryField()[iPatch].size() > 0)
        {
            if (mesh_.boundaryMesh()[iPatch].coupled()) {
                forAll(mesh_.boundaryMesh()[iPatch].faceAreas(), ii )
                {
                    scalar rhoBN = (rho_.boundaryField()[iPatch].patchNeighbourField())()[ii];
//                    Info << "rhoBN:" << rhoBN << endl;
                    scalar rhoBI = (rho_.boundaryField()[iPatch].patchInternalField())()[ii];
//                    Info << "rhoBI:" << rhoBI << endl;
                    vector rhoGradBN = (gradrho.boundaryField()[iPatch].patchNeighbourField())()[ii];
//                    Info << "rhoGradBN:" << rhoGradBN << endl;
                    vector rhoGradBI = (gradrho.boundaryField()[iPatch].patchInternalField())()[ii];
//                    Info << "rhoGradBI:" << rhoGradBI << endl;
                    vector rhoUBN = (rhoU_.boundaryField()[iPatch].patchNeighbourField())()[ii];
//                    Info << "rhoUBN:" << rhoUBN << endl;
                    vector rhoUBI = (rhoU_.boundaryField()[iPatch].patchInternalField())()[ii];
//                    Info << "rhoUBI:" << rhoUBI << endl;
                    tensor uGradBN = (gradU.boundaryField()[iPatch].patchNeighbourField())()[ii];
//                    Info << "uGradBN:" << uGradBN << endl;
                    tensor uGradBI = (gradU.boundaryField()[iPatch].patchInternalField())()[ii];
//                    Info << "uGradBI:" << uGradBI << endl;
                    vector CBN = (mesh_.C().boundaryField()[iPatch].patchNeighbourField())()[ii];
//                    Info << "CBN:" << CBN << endl;
                    vector CBI = (mesh_.C().boundaryField()[iPatch].patchInternalField())()[ii];
//                    Info << "CBI:" << CBI << endl;
                    vector CBf = mesh_.Cf().boundaryField()[iPatch][ii];
//                    Info << "CBf:" << CBf << endl;

                    label faceLabel = mesh_.boundaryMesh()[iPatch].start() + ii;
                    const labelList& fB = fs[faceLabel];

                    vector t1B = (p[fB[1]] - p[fB[0]]);
                    // n_i_        = faceNormal[i];
//                  n_i_B        = mesh.Sf().boundaryField()[iPatch][ii]/mesh.magSf().boundaryField()[iPatch][ii];
                    vector n_i_B        = mesh_.boundaryMesh()[iPatch].faceAreas()[ii]/mag(mesh_.boundaryMesh()[iPatch]
.faceAreas()[ii]);

                    vector t2B = n_i_B^t1B;
                    vector t1_i_B = t1B/mag(t1B);
                    vector t2_i_B = t2B/mag(t2B);

                    icoGKSFlux(  rhoBI,
                                 rhoBN,
                                 rhoUBI/rhoBI,
                                 rhoUBN/rhoBN,
                                 uGradBI,
                                 uGradBN,
                                 rhoGradBI,
                                 rhoGradBN,
                                 CBI,
                                 CBN,
                                 CBf,
                                 n_i_B,
                                 t1_i_B,
                                 t2_i_B,
                                 nuEff.boundaryField()[iPatch][ii],
                                 slen.boundaryField()[iPatch][ii],
                                 Cspeed,
                                 Frho_,
                                 Fm_,
                                 Fe_
                              );
                    // Output
                    phi.boundaryFieldRef()[iPatch][ii] = Frho_*mag(mesh_.boundaryMesh()[iPatch].faceAreas()[ii]);
                    phiUp.boundaryFieldRef()[iPatch][ii]   = Fm_*mag(mesh_.boundaryMesh()[iPatch].faceAreas()[ii]);
                    phiEp.boundaryFieldRef()[iPatch][ii]   = Fe_*mag(mesh_.boundaryMesh()[iPatch].faceAreas()[ii]);
                }
            } else {
	       // Info << "treating boundary conditions :"<<BCType<< endl;
		forAll( mesh_.boundaryMesh()[iPatch].faceAreas(), ii )
		{
		    id_L        = mesh_.boundaryMesh()[iPatch].faceCells()[ii];
		    // i           = mesh.boundaryMesh()[iPatch].whichFace(ii);

		    rho_L_      = rho_.boundaryField()[iPatch][ii];
		    U_L_        = U_.boundaryField()[iPatch][ii];
		    scalar nuFace = nuEff.boundaryField()[iPatch][ii];
		    scalar p_L_ =   rho_L_*Cspeed*Cspeed/3.0;

		    Sf   =  mag(mesh_.boundaryMesh()[iPatch].faceAreas()[ii] );

		    n_i_ = mesh_.boundaryMesh()[iPatch].faceAreas()[ii];
		    n_i_ = n_i_/mag(n_i_);

		    tensor  gradU_b = gradU.boundaryField()[iPatch][ii];

		    viscousFlux( n_i_, rho_L_, p_L_,gradU_b, nuFace ,Frho_, Fm_ );

		    phi.boundaryFieldRef()[iPatch][ii]      = ( rho_L_* (U_L_&n_i_)   +   Frho_)*Sf;
		    phiUp.boundaryFieldRef()[iPatch][ii]        = ( rho_L_* U_L_*(U_L_&n_i_) + Fm_)*Sf;

		    scalar H = 3.5*Cspeed*Cspeed/3.0 + magSqr(U_L_)/2;
		    phiEp.boundaryFieldRef()[iPatch][ii]        = ( rho_L_*(U_L_&n_i_)*H)*Sf;
                }
            }
        }
     }
}

void Foam::gksFluxScheme::viscousFlux
(
    vector n, 
    scalar rho, 
    scalar P,
    tensor gradU,
    scalar nu,
    scalar& Grho,
    vector& Gm
)
{
    // Variables definition
    vector S;

    // Strain rate and divergence of velocity field
    S.x() = gradU.xx()*n.x() + gradU.xy()*n.y() + gradU.xz()*n.z() +
            gradU.xx()*n.x() + gradU.yx()*n.y() + gradU.zx()*n.z();
    S.y() = gradU.yx()*n.x() + gradU.yy()*n.y() + gradU.yz()*n.z() +
            gradU.xy()*n.x() + gradU.yy()*n.y() + gradU.zy()*n.z();
    S.z() = gradU.zx()*n.x() + gradU.zy()*n.y() + gradU.zz()*n.z() +
            gradU.xz()*n.x() + gradU.yz()*n.y() + gradU.zz()*n.z();



    //viscous fluxes
    Grho = 0.0;
    Gm   = P*n-rho*nu*S ;
}

void Foam::gksFluxScheme::icoGKSFlux
(
            const scalar& rhoLeft,
            const scalar& rhoRight,
            const vector& ULeft,
            const vector& URight,
            const tensor& gradULeft,
            const tensor& gradURight,
            const vector& gradRhoLeft,
            const vector& gradRhoRight,
            const vector& leftCenter,
            const vector& rightCenter,
            const vector& faceCenter,
            const vector& normalVector,
            const vector& t1,
            const vector& t2,
            const scalar& nuEff,
            const scalar& streamingLength,
            const scalar& Cspeed,
            scalar& rhoFlux,
            vector& rhoUFlux,
            scalar& rhoEFlux
)
{
    //Info <<"here icoGKSFlux" <<endl;
    const scalar n1x= normalVector.x();
    const scalar n1y= normalVector.y();
    const scalar n1z= normalVector.z();

    const scalar n2x= t1.x();
    const scalar n2y= t1.y();
    const scalar n2z= t1.z();

    const scalar n3x= t2.x();
    const scalar n3y= t2.y();
    const scalar n3z= t2.z();

    scalar DenL,DenR;
    vector UL,UR;
    scalar delta_rho1,delta_rho2;
    vector delta_v1,delta_v2;

    delta_rho1 = (faceCenter- leftCenter)&gradRhoLeft;
    delta_rho2 = (faceCenter- rightCenter)&gradRhoRight;
    DenL = rhoLeft + delta_rho1 ;
    DenR = rhoRight + delta_rho2 ;
    delta_v1 = (faceCenter- leftCenter)&gradULeft;
    delta_v2 = (faceCenter- rightCenter)&gradURight;
    UL = ULeft + delta_v1    ;
    UR = URight + delta_v2    ;

    scalar wL1 = DenL;
    scalar wL2 = UL.x()*n1x+UL.y()*n1y+UL.z()*n1z;
    scalar wL3 = UL.x()*n2x+UL.y()*n2y+UL.z()*n2z;
    scalar wL4 = UL.x()*n3x+UL.y()*n3y+UL.z()*n3z;

    scalar wR1 = DenR;
    scalar wR2 = UR.x()*n1x+UR.y()*n1y+UR.z()*n1z;
    scalar wR3 = UR.x()*n2x+UR.y()*n2y+UR.z()*n2z;
    scalar wR4 = UR.x()*n3x+UR.y()*n3y+UR.z()*n3z;

    scalar tp0 = 1./(4.*3.141592654);

    //Left
    scalar tpu1 = gradULeft.xx()* n1x + gradULeft.yx()*n1y + gradULeft.zx()*n1z;
    scalar tpu2 = gradULeft.xy()* n1x + gradULeft.yy()*n1y + gradULeft.zy()*n1z;
    scalar tpu3 = gradULeft.xz()* n1x + gradULeft.yz()*n1y + gradULeft.zz()*n1z;
    scalar tpv1 = gradULeft.xx()* n2x + gradULeft.yx()*n2y + gradULeft.zx()*n2z;
    scalar tpv2 = gradULeft.xy()* n2x + gradULeft.yy()*n2y + gradULeft.zy()*n2z;
    scalar tpv3 = gradULeft.xz()* n2x + gradULeft.yz()*n2y + gradULeft.zz()*n2z;
    scalar tpw1 = gradULeft.xx()* n3x + gradULeft.yx()*n3y + gradULeft.zx()*n3z;
    scalar tpw2 = gradULeft.xy()* n3x + gradULeft.yy()*n3y + gradULeft.zy()*n3z;
    scalar tpw3 = gradULeft.xz()* n3x + gradULeft.yz()*n3y + gradULeft.zz()*n3z;

    scalar gduL1 = tpu1*n1x+tpu2*n1y+tpu3*n1z; //d_u1/d_x1
    scalar gduL2 = tpu1*n2x+tpu2*n2y+tpu3*n2z; //d_u1/d_x2
    scalar gduL3 = tpu1*n3x+tpu2*n3y+tpu3*n3z; //d_u1/d_x3
    scalar gdvL1 = tpv1*n1x+tpv2*n1y+tpv3*n1z; //d_u2/d_x1
    scalar gdvL2 = tpv1*n2x+tpv2*n2y+tpv3*n2z; //d_u2/d_x2
    scalar gdvL3 = tpv1*n3x+tpv2*n3y+tpv3*n3z; //d_u2/d_x3
    scalar gdwL1 = tpw1*n1x+tpw2*n1y+tpw3*n1z; //d_u3/d_x1
    scalar gdwL2 = tpw1*n2x+tpw2*n2y+tpw3*n2z; //d_u3/d_x2
    scalar gdwL3 = tpw1*n3x+tpw2*n3y+tpw3*n3z; //d_u3/d_x3
    scalar gdgL1 = (gradRhoLeft&normalVector)*tp0;
    scalar gdgL2 = (gradRhoLeft&t1)*tp0;
    scalar gdgL3 = (gradRhoLeft&t2)*tp0;

    //Right
    tpu1 = gradURight.xx()* n1x + gradURight.yx()*n1y + gradURight.zx()*n1z;
    tpu2 = gradURight.xy()* n1x + gradURight.yy()*n1y + gradURight.zy()*n1z;
    tpu3 = gradURight.xz()* n1x + gradURight.yz()*n1y + gradURight.zz()*n1z;
    tpv1 = gradURight.xx()* n2x + gradURight.yx()*n2y + gradURight.zx()*n2z;
    tpv2 = gradURight.xy()* n2x + gradURight.yy()*n2y + gradURight.zy()*n2z;
    tpv3 = gradURight.xz()* n2x + gradURight.yz()*n2y + gradURight.zz()*n2z;
    tpw1 = gradURight.xx()* n3x + gradURight.yx()*n3y + gradURight.zx()*n3z;
    tpw2 = gradURight.xy()* n3x + gradURight.yy()*n3y + gradURight.zy()*n3z;
    tpw3 = gradURight.xz()* n3x + gradURight.yz()*n3y + gradURight.zz()*n3z;

    scalar gduR1 = tpu1*n1x+tpu2*n1y+tpu3*n1z; //d_u1/d_x1
    scalar gduR2 = tpu1*n2x+tpu2*n2y+tpu3*n2z; //d_u1/d_x2
    scalar gduR3 = tpu1*n3x+tpu2*n3y+tpu3*n3z; //d_u1/d_x3
    scalar gdvR1 = tpv1*n1x+tpv2*n1y+tpv3*n1z; //d_u2/d_x1
    scalar gdvR2 = tpv1*n2x+tpv2*n2y+tpv3*n2z; //d_u2/d_x2
    scalar gdvR3 = tpv1*n3x+tpv2*n3y+tpv3*n3z; //d_u2/d_x3
    scalar gdwR1 = tpw1*n1x+tpw2*n1y+tpw3*n1z; //d_u3/d_x1
    scalar gdwR2 = tpw1*n2x+tpw2*n2y+tpw3*n2z; //d_u3/d_x2
    scalar gdwR3 = tpw1*n3x+tpw2*n3y+tpw3*n3z; //d_u3/d_x3
    scalar gdgR1 = (gradRhoRight&normalVector)*tp0;
    scalar gdgR2 = (gradRhoRight&t1)*tp0;
    scalar gdgR3 = (gradRhoRight&t2)*tp0;

// Step 2: compute Roe averged quantities for face:

    if(DenL<0.0||DenR<0.0)
    {
        Info<<"DenL = "<< DenL <<"  DenR = "<<DenR<<endl;
        Info<<"rhoLeft = "<< rhoLeft <<"  rhoRight = "<<rhoRight<<endl;
    }
    // Some temporary variables:
    const scalar rhoLeftSqrt = Foam::sqrt(DenL);
    const scalar rhoRightSqrt =  Foam::sqrt(DenR);


    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1.0 - wLeft;

    const scalar wM1 = wL1*wLeft + wR1*wRight;
    const scalar wM2 = wL2*wLeft + wR2*wRight;
    const scalar wM3 = wL3*wLeft + wR3*wRight;
    const scalar wM4 = wL4*wLeft + wR4*wRight;
    const scalar cf =  Cspeed;

    // const scalar LeftDistance =  mag(faceCenter- leftCenter);
    // const scalar RightDistance = mag( faceCenter- rightCenter)a

    scalar dt  = 0.05*streamingLength/( max(max(mag(wM2),mag(wM3)),mag(wM4)) +Cspeed );


    // Info <<"here icoGKSFlux" <<deltaT<< endl;
    scalar gL = wL1*tp0;
    scalar gR = wR1*tp0;

    //coefficients
    scalar a0L = wL2-(gduL1*wM2+gduL2*wM3+gduL3*wM4)*dt;
    scalar a1L = cf-gduL1*cf*dt;
    scalar a2L = -gduL2*cf*dt;
    scalar a3L = -gduL3*cf*dt;
    scalar b0L = wL3-(gdvL1*wM2+gdvL2*wM3+gdvL3*wM4)*dt;
    scalar b1L = -gdvL1*cf*dt;
    scalar b2L = cf-gdvL2*cf*dt;
    scalar b3L = -gdvL3*cf*dt;
    scalar c0L = wL4-(gdwL1*wM2+gdwL2*wM3+gdwL3*wM4)*dt;
    scalar c1L = -gdwL1*cf*dt;
    scalar c2L = -gdwL2*cf*dt;
    scalar c3L = cf-gdwL3*cf*dt;
    scalar g0L = gL-(gdgL1*wM2+gdgL2*wM3+gdgL3*wM4)*dt;
    scalar g1L = -gdgL1*cf*dt;
    scalar g2L = -gdgL2*cf*dt;
    scalar g3L = -gdgL3*cf*dt;
    scalar a0R = wR2-(gduR1*wM2+gduR2*wM3+gduR3*wM4)*dt;
    scalar a1R = cf-gduR1*cf*dt;
    scalar a2R = -gduR2*cf*dt;
    scalar a3R = -gduR3*cf*dt;
    scalar b0R = wR3-(gdvR1*wM2+gdvR2*wM3+gdvR3*wM4)*dt;
    scalar b1R = -gdvR1*cf*dt;
    scalar b2R = cf-gdvR2*cf*dt;
    scalar b3R = -gdvR3*cf*dt;
    scalar c0R = wR4-(gdwR1*wM2+gdwR2*wM3+gdwR3*wM4)*dt;
    scalar c1R = -gdwR1*cf*dt;
    scalar c2R = -gdwR2*cf*dt;
    scalar c3R = cf-gdwR3*cf*dt;
    scalar g0R = gR-(gdgR1*wM2+gdgR2*wM3+gdgR3*wM4)*dt;
    scalar g1R = -gdgR1*cf*dt;
    scalar g2R = -gdgR2*cf*dt;
    scalar g3R = -gdgR3*cf*dt;

    tp0 = 2./3.;
    scalar tp1 = 2.*tp0;
    scalar tp2 = 1./4.;
    scalar pi = 3.141592654;

    //conservative variables at the cell interface
    scalar W01 = ((2.*g0L+g1L) + (2.*g0R-g1R))*pi;
    scalar W02 =
                  pi*
                  (
                      (2.*a0L*g0L+a0L*g1L+a1L*g0L+tp0*(a1L*g1L+a2L*g2L+a3L*g3L))
                    + (2.*a0R*g0R-a0R*g1R-a1R*g0R+tp0*(a1R*g1R+a2R*g2R+a3R*g3R))
                  );
    scalar W03 =
                  pi*
                  (
                      (2.*b0L*g0L+b0L*g1L+b1L*g0L+tp0*(b1L*g1L+b2L*g2L+b3L*g3L))
                    + (2.*b0R*g0R-b0R*g1R-b1R*g0R+tp0*(b1R*g1R+b2R*g2R+b3R*g3R))
                  );
    scalar W04 =
                  pi*
                  (
                      (2.*c0L*g0L+c0L*g1L+c1L*g0L+tp0*(c1L*g1L+c2L*g2L+c3L*g3L))
                    + (2.*c0R*g0R-c0R*g1R-c1R*g0R+tp0*(c1R*g1R+c2R*g2R+c3R*g3R))
                  );
    // constant temperature (rho*E+p)*U==(E+p/rho)*rhoU==H*rhoU
    scalar u1 = W02/W01;
    scalar u2 = W03/W01;
    scalar u3 = W04/W01;

    scalar W05 = 2.5*cf*cf/3 + (u1*u1+u2*u2+u3*u3)/2;
    scalar H = 3.5*cf*cf/3 + (u1*u1+u2*u2+u3*u3)/2;

    //flux attributed to g_face
    scalar pp  = W01*cf*cf/3.0;//Dim=2
    scalar F11 = W02;
    scalar F12 = W02*W02/W01+pp;
    scalar F13 = W02*W03/W01;
    scalar F14 = W02*W04/W01;
    scalar F15 = W02*H;

    //flux attributed to g_sphere
    scalar F21 = W02;
    scalar F22 =
                 pi*
                 (
                      (
                          tp0*g0L*(a1L*a1L+a2L*a2L+a3L*a3L)
                        + (tp1*a0L+0.5*a1L)*(a1L*g1L+a2L*g2L+a3L*g3L)
                        + 2.*a0L*g0L*(a0L+a1L)
                        + tp2*g1L*(4.*a0L*a0L+a2L*a2L+a3L*a3L)
                      )
                      +
                      (
                           tp0*g0R*(a1R*a1R+a2R*a2R+a3R*a3R)
                         + (tp1*a0R-0.5*a1R)*(a1R*g1R+a2R*g2R+a3R*g3R)
                         + 2.*a0R*g0R*(a0R-a1R)
                         - tp2*g1R*(4.*a0R*a0R+a2R*a2R+a3R*a3R)
                      )
                 );
    scalar F23 =
                 pi*
                 (
                      (
                          tp0*g0L*(3.*a0L*b0L+a1L*b1L+a2L*b2L+a3L*b3L)
                        + tp0*a0L*(b1L*g1L+b2L*g2L+b3L*g3L)
                        + tp0*b0L*(a1L*g1L+a2L*g2L+a3L*g3L)
                        + g0L*(a0L*b1L+a1L*b0L)
                        + tp2*g3L*(a1L*b3L+a3L*b1L)
                        + tp2*g2L*(a2L*b1L+a1L*b2L)+tp2*g1L*(4.*a0L*b0L+2.*a1L*b1L+a2L*b2L+a3L*b3L)
                      )
                      +
                      (
                          tp0*g0R*(3.*a0R*b0R+a1R*b1R+a2R*b2R+a3R*b3R)
                        + tp0*a0R*(b1R*g1R+b2R*g2R+b3R*g3R)
                        + tp0*b0R*(a1R*g1R+a2R*g2R+a3R*g3R)
                        - g0R*(a0R*b1R+a1R*b0R)-tp2*g3R*(a1R*b3R+a3R*b1R)
                        - tp2*g2R*(a2R*b1R+a1R*b2R)-tp2*g1R*(4.*a0R*b0R+2.*a1R*b1R+a2R*b2R+a3R*b3R)
                      )
                 );

    scalar F24 =
                 pi*
                (
                      (
                          tp0*g0L*(3.*a0L*c0L+a1L*c1L+a2L*c2L+a3L*c3L)
                        + tp0*a0L*(c1L*g1L+c2L*g2L+c3L*g3L)
                        + tp0*c0L*(a1L*g1L+a2L*g2L+a3L*g3L)
                        + g0L*(a0L*c1L+a1L*c0L)
                        + tp2*g3L*(a1L*c3L+a3L*c1L)
                        + tp2*g2L*(a2L*c1L+a1L*c2L)
                        + tp2*g1L*(4.*a0L*c0L+2.*a1L*c1L+a2L*c2L+a3L*c3L)
                      )
                      +
                      (
                          tp0*g0R*(3.*a0R*c0R+a1R*c1R+a2R*c2R+a3R*c3R)
                        + tp0*a0R*(c1R*g1R+c2R*g2R+c3R*g3R)
                        + tp0*c0R*(a1R*g1R+a2R*g2R+a3R*g3R)
                        - g0R*(a0R*c1R+a1R*c0R)
                        - tp2*g3R*(a1R*c3R+a3R*c1R)
                        - tp2*g2R*(a2R*c1R+a1R*c2R)
                        - tp2*g1R*(4.*a0R*c0R+2.*a1R*c1R+a2R*c2R+a3R*c3R)
                      )
                 );

    //total flux
     scalar taoStar= nuEff/pp/dt;
     //scalar taoStar= nuEff/3.0/magSqr(Cspeed)/deltaT;

     scalar rhoTemp;
     vector rhoUTemp;
     scalar rhoETemp;

     rhoTemp        = (1.0-taoStar)*F11+taoStar*F21;
     rhoUTemp.x()   = (1.0-taoStar)*F12+taoStar*F22;
     rhoUTemp.y()   = (1.0-taoStar)*F13+taoStar*F23;
     rhoUTemp.z()   = (1.0-taoStar)*F14+taoStar*F24;
     rhoETemp       =  F15;
  // Compute conservative variables
     rhoFlux        = rhoTemp;
     rhoUFlux.x()   = rhoUTemp.x()*n1x+ rhoUTemp.y()*n2x+ rhoUTemp.z()*n3x;
     rhoUFlux.y()   = rhoUTemp.x()*n1y+ rhoUTemp.y()*n2y+ rhoUTemp.z()*n3y;
     rhoUFlux.z()   = rhoUTemp.x()*n1z+ rhoUTemp.y()*n2z+ rhoUTemp.z()*n3z;
     rhoEFlux        = rhoETemp;
}

void Foam::gksFluxScheme::calcStreamingDistant(surfaceScalarField& slen)
{
    const volVectorField &center = mesh_.C();
    const surfaceVectorField  &surfacecenter= mesh_.Cf();


    const labelListList & faceEdges = mesh_.faceEdges();
    //determine streaming length for each edge
    forAll (mesh_.owner() ,iface)
    {
	scalar length= GREAT;

	    const labelList& edgenumberfacei = faceEdges[iface];
	    //const  vector & x1= mesh.Sf()[iface];
	    forAll(edgenumberfacei,edgei)
	    {
		const label& edgenumber = edgenumberfacei[edgei];
		const label& startPoint = mesh_.edges()[edgenumber][0];
		const label& endPoint = mesh_.edges()[edgenumber][1];
		vector x2= mesh_.points()[endPoint]- mesh_.points()[startPoint];
		length = Foam::min(length, mag(x2));
	    }
       //Info << length << endl;
       label leftCell , rightCell;

      // Get the left and right cell index    // scalar niu=nu.value();
      leftCell  = mesh_.owner()[iface];
      rightCell = mesh_.neighbour()[iface];

      // Get face center
      point face_center = surfacecenter[iface];

      //Get cell center
      point rightCenter = center[rightCell];
      point leftCenter  = center[leftCell];
      //scalar min;
      scalar a,b;

      a= mag(leftCenter-face_center);
      b= mag(rightCenter-face_center);

    /*            if (iface==288389)
		{
		    Info << a << endl;
		    Info << b << endl;
		}*/
      slen[iface]= Foam::min(0.5*length , Foam::min(a,b));
    }

    forAll(mesh_.boundaryMesh(), iPatch)
    {

	if (mesh_.boundaryMesh()[iPatch].coupled())
	{
		forAll(mesh_.boundaryMesh()[iPatch], fi)
		{
		    label faceLabel = mesh_.boundaryMesh()[iPatch].start() + fi;
		    scalar length= GREAT;

		    const labelList& edgenumberfacei = faceEdges[faceLabel];
		    //const  vector & x1= mesh.Sf()[iface];
		    forAll(edgenumberfacei,edgei)
		    {
			const label& edgenumber = edgenumberfacei[edgei];
			const label& startPoint = mesh_.edges()[edgenumber][0];
			const label& endPoint = mesh_.edges()[edgenumber][1];
			vector x2= mesh_.points()[endPoint]- mesh_.points()[startPoint];
			length = Foam::min(length, mag(x2));
		    }

		    // Get face center
		    point face_center = surfacecenter.boundaryField()[iPatch][fi];

		    //Get cell center
		    point rightCenter = (center.boundaryField()[iPatch].patchNeighbourField())()[fi];
		    point leftCenter  = (center.boundaryField()[iPatch].patchInternalField())()[fi];
		    //scalar min;
		    scalar a,b;

		    a= mag(leftCenter-face_center);
		    b= mag(rightCenter-face_center);
		    slen.boundaryFieldRef()[iPatch][fi]= Foam::min(0.5*length , Foam::min(a,b));
	        }
	 }
    }
}

/*
void Foam::gksFluxScheme::icoGKSFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{

}
*/

/*
void Foam::gksFluxScheme::calcGKSCord(int i)
{
    const faceList& fs = mesh.faces();
    const pointField& p = mesh.points();

    const labelList& f = fs[i];
    vector t1 = (p[f[1]] - p[f[0]]);
    // n_i_        = faceNormal[i];
    n_i_        = mesh.Sf()[i]/mesh.magSf()[i];

    vector t2 = n_i_^t1;
    vector t1_i_ = t1/mag(t1);
    vector t2_i_ = t2/mag(t2);
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
