/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2016 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2016 Oliver Oxtoby - CSIR, South Africa
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

Class
    Foam::lbfsScheme

Description
    AUSM+Up flux splitting scheme as described by Liou 2006, A sequal to
    AUSM, Part II: AUSM+-up for all speeds

SourceFiles
    lbfsScheme.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include <fstream>
#include "fvMeshSubset.H"

inline double Foam::hisaModule::deltaFunc(double r)
{
    double deltaFunc = 0;
    double absR = fabs(r);
    if(absR <= 1)
    {
        deltaFunc = (3-2*absR+sqrt(1+4*absR-4*absR*absR))/8.;
    } else if(absR <= 2)
    {
        deltaFunc = (5-2*absR-sqrt(-7+12*absR-4*absR*absR))/8.;
    } else
    {
        deltaFunc = 0;
    }

    return deltaFunc;
}

inline void Foam::hisaModule::readLagragianPoints()
{
    const IOdictionary dictIB
    (
        IOobject
        (
            "immersedBoundaryProperties",
            mesh_().time().caseConstant(),
            mesh_(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false // Do not register
        )
    );

    Info << "Reading Lagrangian points" <<endl;
    IOList<vector> test
    (
         IOobject
         (
             "LagrangianPoints_I",
             mesh_().time().caseConstant(),
             mesh_(),
             IOobject::MUST_READ,
             IOobject::NO_WRITE
         )
    );
    Info << "Lagrangian points Loaded" <<endl;

    //double refD = 1.0;
    //double deltaTheta;
    //double PI = 3.14159265359;

    //h_ = 1./100;//2.5/1000;
    //N_L= 225/5.;
    h_ = dictIB.lookupOrDefault<scalar>("h", 1.0e-2);
    
    Info << "hahahahaha " << h_ <<endl;

    N_L=test.size();

    solidObject_ = List<SolidObject>(N_L);

    for(int i = 0; i<N_L; i++)
    {   
        solidObject_[i].x = test[i].x();
        solidObject_[i].y = test[i].y();
//        solidObject_[i].z = test[i].z();
    }
    for(int i = 0; i<N_L; i++)
    {
	solidObject_[i].x0 = solidObject_[i].x;
	solidObject_[i].y0 = solidObject_[i].y;
//	solidObject_[i].z0 = solidObject_[i].z;
    }

    //*****************************start overset****************************// 
    const fvMesh& mesh_whole = mesh_();
    zoneID_.set(&(mesh_whole.lookupObject<labelIOList>("zoneIDList")));
    //zoneID_.set(&(mesh_whole.lookupObject<volScalarField>("zoneID")));
    const labelIOList& zoneID = zoneID_(); 
   
    //find the moving zone 
    oversetZoneID_ = 1;
    //const fvBoundaryMesh& fvm = mesh_whole.boundary();
    //// 'overset' patches
    //forAll(fvm, patchI)
    //{
    //    if (isA<oversetFvPatch>(fvm[patchI]))
    //    {
    //         const labelList& fc = fvm[patchI].faceCells();
    //         label cell0 = fc[0];
    //         oversetZoneID_ = zoneID[cell0];
    //     }
    //}

   
    //label nZones = gMax(zoneID)+1;
    //PtrList<fvMeshSubset> meshParts(1);
    meshParts_.setSize(1);
    meshParts_.set
    (
        0,
        new fvMeshSubset(mesh_whole)
    );
    meshParts_[0].setLargeCellSubset(static_cast<labelList>(zoneID_), oversetZoneID_);
//    meshParts_[0].setLargeCellSubset(static_cast<labelList>(zoneID_.internalField()), oversetZoneID_); 
    const fvMesh& mesh = meshParts_[0].subMesh();

//    Info<< meshParts_[0].cellMap() <<endl;
//    Info<< oversetZoneID_ <<endl;
//    cin.get();
//    fvMesh& mesh = mesh_();
    //*****************************end overset****************************// 
    D_ = List<DynamicList<double>>(N_L);
    DCellID_ = List<DynamicList<int>>(N_L);

    for(int i = 0; i<N_L; i++)
    {
        D_[i].clear();
        DCellID_[i].clear();
    }
    for(int i = 0; i<N_L; i++)
    {
        for(int j = 0; j<mesh.C().size(); j++)
        {
            double rx = fabs(mesh.C()[j].x()-solidObject_[i].x)/h_;
            double ry = fabs(mesh.C()[j].y()-solidObject_[i].y)/h_;
            //double rz = fabs(mesh.C()[j].z()-solidObject_[i].z)/h_;
            if ( (rx<2) && (ry<2))
            //if ( (rx<2) && (ry<2) && (rz<2))
            {
                D_[i].append(deltaFunc(rx)*deltaFunc(ry));
                //D_[i].append(deltaFunc(rx)*deltaFunc(ry)*deltaFunc(rz));
                DCellID_[i].append(j);
            }
        }
        // Info << DCellID_[i].size() << endl;
    }

    
    hostCell_ = List<int>(N_L);
    for(int i = 0; i<N_L; i++)
    {
        point location(0,0,0);
        location.x() = solidObject_[i].x;
        location.y() = solidObject_[i].y;
        hostCell_[i] = findNearestCell(location, 0);
    }

    neighbourCells_ = List<NeighbourCells>(mesh.C().size());
    for(int j = 0; j<mesh.C().size(); j++)
    {
        neighbourCells_[j].w_=-1;
        neighbourCells_[j].e_=-1;
        neighbourCells_[j].n_=-1;
        neighbourCells_[j].s_=-1;

        // only neighbouring cells
        const labelList& cells = mesh.cellCells()[j];
        forAll(cells, idx)
        {
            label cell = cells[idx];
            vector direction = mesh.cellCentres()[cell] - mesh.cellCentres()[j];
            direction = direction/mag(direction);
            if(direction.x()>0.5)
            {
        	neighbourCells_[j].e_=cell;
            } else if(direction.x()<-0.5) {
        	neighbourCells_[j].w_=cell;
            } else if(direction.y()>0.5) {
        	neighbourCells_[j].n_=cell;
            } else if(direction.y()<-0.5) {
        	neighbourCells_[j].s_=cell;
            }
        }
    }

    //D_ = List<List<double>>(N_L);
    //DCellID_ = List<List<int>>(N_L);
    //for(int i = 0; i<N_L; i++)
    //{
    //    int N_nonzero = 0;
    //    for(int j = 0; j<mesh.C().size(); j++)
    //    {
    //        double rx = fabs(mesh.C()[j].x()-solidObject_[i].x)/h_;
    //        double ry = fabs(mesh.C()[j].y()-solidObject_[i].y)/h_;
    //        if ( (rx<2) && (ry<2))
    //        {
    //            N_nonzero++;
    //        }
    //    }
    //    D_[i] = List<double>(N_nonzero);
    //    DCellID_[i] = List<int>(N_nonzero);
    //}

    //for(int i = 0; i<N_L; i++)
    //{
    //    int N_nonzero = 0;
    //    for(int j = 0; j<mesh.C().size(); j++)
    //    {
    //        double rx = fabs(mesh.C()[j].x()-solidObject_[i].x)/h_;
    //        double ry = fabs(mesh.C()[j].y()-solidObject_[i].y)/h_;
    //        if ( (rx<2) && (ry<2))
    //        {
    //            D_[i][N_nonzero] = deltaFunc(rx)*deltaFunc(ry);
    //            DCellID_[i][N_nonzero] = j; 
    //            N_nonzero++;
    //        }
    //    }
    //}

    deltaU_ = vectorField(mesh_->C().size(), vector(0, 0, 0));

    // I_ = PI*refD*refD/64.;
    // M_ = 1.25*PI*(refD*refD/4.);
    // rhoS_ = 1.25;
    // rhoF_ = 1;

    // g_ = 9.8*0.4823;

    xc_ =0;
    yc_ =0;

    ut_ = 0;
    vt_ = 0;
    vt_o_ = 0;
    omega_ = 0;

    oversetMotion_[0].x() = xc_;
    oversetMotion_[0].y() = yc_;
    oversetMotion_[0].z() = 0;

    desiredProperties();

    Info<<"IBM initialization done." << endl;
}

inline void Foam::hisaModule::desiredProperties()
{
    Ub_ = vectorField(N_L, vector(0, 0, 0));
    // Tb_ = scalarField(N_L,.25);
}

inline void Foam::hisaModule::calcB(int cmpt, const labelList& cellMap)
{
    volVectorField& U = U_();
	
    for(int i = 0; i<N_L; i++)
    {
        double rowSum = 0;
        for(int j = 0; j<D_[i].size(); j++)
        {
            double u = U.internalField()[cellMap[DCellID_[i][j]]].component(cmpt);
            rowSum += D_[i][j]*u;
        }
        B_[i] = Ub_[i].component(cmpt) - rowSum;
    }
}

inline void Foam::hisaModule::Amul(scalarField& wA, scalarField& psi)
{
    const fvMesh& mesh = meshParts_[0].subMesh();
    scalarField x_star(mesh.C().size(), 0);

    for(int i = 0; i<N_L; i++)
    {
        for(int j = 0; j<D_[i].size(); j++)
        {
            double rowSum = D_[i][j]*psi[i];
            x_star[DCellID_[i][j]] += rowSum;
        }
    }

    for(int i = 0; i<N_L; i++)
    {
        double rowSum = 0;
        for(int j = 0; j<D_[i].size(); j++)
        {
            rowSum += D_[i][j]*x_star[DCellID_[i][j]];
        }
        wA[i] = rowSum;
    }
}

inline void Foam::hisaModule::Tmul(scalarField& wA, scalarField& psi)
{
    const fvMesh& mesh = meshParts_[0].subMesh();
    scalarField x_star(mesh.C().size(), 0);

    for(int i = 0; i<N_L; i++)
    {
        for(int j = 0; j<D_[i].size(); j++)
        {
            double rowSum = D_[i][j]*psi[i];
            x_star[DCellID_[i][j]] += rowSum;
        }
    }

    for(int i = 0; i<N_L; i++)
    {
        double rowSum = 0;
        for(int j = 0; j<D_[i].size(); j++)
        {
            rowSum += D_[i][j]*x_star[DCellID_[i][j]];
        }
        wA[i] = rowSum;
    }
}

//inline void Foam::hisaModule::Amul(scalarField& wA, scalarField& psi)
//{
//    for(int i = 0; i<N_L; i++)
//    {
//        double rowSum = 0;
//        for(int j = 0; j<A_[i].size(); j++)
//        {
//            rowSum += A_[i][j]*psi[APointID_[i][j]];
//        }
//        wA[i] = rowSum;
//    }
//}

//inline void Foam::hisaModule::Tmul(scalarField& wT, scalarField& psi)
//{
//    for(int i = 0; i<N_L; i++)
//    {
//        double rowSum = 0;
//        for(int j = 0; j<A_[i].size(); j++)
//        {
//            rowSum += A_[i][j]*psi[APointID_[i][j]];
//        }
//        wT[i] = rowSum;
//    }
//}


inline void Foam::hisaModule::solveX()
{
    //solve X using CG.
    //PBiCGsolve(X_, B_);
    explicitSolve(X_, B_);
}

inline scalar Foam::hisaModule::calcNormFactor
(
    const scalarField& psi,
    const scalarField& source,
    const scalarField& Apsi,
    scalarField& tmpField
)
{
    // --- Calculate A dot reference value of psi
    for (int i=0; i<N_L; i++)
    {
        tmpField[i] += 0;
        for (int j=0; j<D_[i].size(); j++)
        {
            tmpField[i] += D_[i][j];
        }
    }

    tmpField *= gAverage(psi);

    return
        gSum
        (
            (mag(Apsi - tmpField) + mag(source - tmpField))()
        )
      + 1e-20;
}

inline void Foam::hisaModule::explicitSolve
(
    scalarField& psi,
    const scalarField& source
)
{
    scalarField V(N_L, 0);
    const fvMesh& mesh = meshParts_[0].subMesh();
    scalarField x_star(mesh.C().size(), 0);

    for(int i = 0; i<N_L; i++)
    {
        for(int j = 0; j<D_[i].size(); j++)
        {
            double rowSum = D_[i][j];
            x_star[DCellID_[i][j]] += rowSum;
        }
    }

    for(int i = 0; i<N_L; i++)
    {
        double rowSum = 0;
        for(int j = 0; j<D_[i].size(); j++)
        {
            rowSum += D_[i][j]*x_star[DCellID_[i][j]];
        }
        V[i] = rowSum;
    }

    for(int i = 0; i<N_L; i++)
    {
        psi[i] = source[i]/V[i];
    }
}

inline void Foam::hisaModule::PBiCGsolve
(
    scalarField& psi,
    const scalarField& source
)
{
    Info << "solving IBM using PBiCG... "<< endl;

    scalar initialResidual;
    scalar finalResidual;
    scalar nIterations = 0;

//    const label nCells = psi.size();
    scalar* __restrict__ psiPtr = psi.begin();

    scalarField pA(N_L);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField wA(N_L);
    scalar* __restrict__ wAPtr = wA.begin();

    // --- Calculate A*X
    Amul(wA, psi);

    // --- Calculate initial residual field
    scalarField rA(B_ - wA);
    scalar* __restrict__ rAPtr = rA.begin();

    // --- Calculate normalisation factor
    const scalar normFactor = calcNormFactor(psi, source, wA, pA);

    // --- Calculate normalised residual norm
    initialResidual =
        gSumMag(rA)
       /normFactor;
    finalResidual = initialResidual;

    /*if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )*/
    // --- Check convergence, solve if not converged
    if (initialResidual>1e-8)
    {
        scalarField pT(N_L, 0);
        scalar* __restrict__ pTPtr = pT.begin();

        scalarField wT(N_L);
        scalar* __restrict__ wTPtr = wT.begin();
	
        Tmul(wT, psi);

	scalarField rT(source - wT);
        scalar* __restrict__ rTPtr = rT.begin();

        // --- Initial value not used
        scalar wArT = 0;

	do
	{
	    // --- Store previous wArT
            const scalar wArTold = wArT;


	    // --- Precondition residuals
	    wA = rA;
	    wT = rT;

	    // --- Update search directions:
            wArT = gSumProd(wA, rT);

	    if (nIterations == 0)
            {
                for (label i=0; i<N_L; i++)
                {
                    pAPtr[i] = wAPtr[i];
                    pTPtr[i] = wTPtr[i];
                }
            }
            else
            {
                const scalar beta = wArT/wArTold;

                for (label i=0; i<N_L; i++)
                {
                    pAPtr[i] = wAPtr[i] + beta*pAPtr[i];
                    pTPtr[i] = wTPtr[i] + beta*pTPtr[i];
                }
            }

	    // --- Update preconditioned residuals
            Amul(wA, pA);
            Tmul(wT, pT);

	    const scalar wApT = gSumProd(wA, pT);

	    // --- Test for singularity
            if ((wApT/normFactor)<1e-20)
            {
/*		Info << "singularity in solving X" << endl;
		Info <<"pA: " << pA << endl;
		Info <<"pT: " << pT << endl;
		Info <<"wApT: " << wApT << endl;
		Info <<"normFactor: " << normFactor << endl;
		Info <<"psi: " << psi << endl;
		Info <<"source: " << source << endl;
*/                break;
            }

	    // --- Update solution and residual:

            const scalar alpha = wArT/wApT;

            for (label i=0; i<N_L; i++)
            {
                psiPtr[i] += alpha*pAPtr[i];
                rAPtr[i] -= alpha*wAPtr[i];
                rTPtr[i] -= alpha*wTPtr[i];
            }

            finalResidual = gSumMag(rA)/normFactor;

	    nIterations++;
	} while (finalResidual > 1e-8);

//	Info << "X_: "<< endl;
//	Info << X_ << endl;
//	Info << "B_: "<< endl;
//	Info << B_ << endl;
//	Info << "Residuals:" << endl;
//	Info << initialResidual << endl;
//	Info << finalResidual << endl;
//	cin.get();
    }
}

inline void Foam::hisaModule::calcCorrection(label cmpt, const labelList& cellMap)
{
    scalar dt = time_.deltaT().value();
    scalar t = time_.value();
    
    volScalarField& rho = scalarVars_[0];

    scalar Cf = 0;
    scalar Cm = 0;
    for(int i = 0; i<N_L; i++)
    {
        for(int j = 0; j<D_[i].size(); j++)
        {
            double rowSum = D_[i][j]*X_[i];
            deltaU_[cellMap[DCellID_[i][j]]].component(cmpt) += rowSum;
	    Cf += 2*rowSum*rho[cellMap[DCellID_[i][j]]]*h_*h_/dt;
            if (cmpt == 0)
            {
		scalar Ly = solidObject_[i].y - yc_;
	        Cm = Cf*Ly;
            }
            
	    if (cmpt == 1)
            {
		scalar Lx = solidObject_[i].x - xc_;
	        Cm = Cf*Lx;
            }
        }
    }

    if (cmpt == 0)
    {
	Cd_ = -Cf;
	Cm_ = -Cm;
        std::ofstream outfile;
        outfile.open("Cd.txt", std::ios_base::app);
        outfile << t <<"    " <<  Cd_ << "\n";
    }

    if (cmpt == 1)
    {
	Cl_ = -Cf;
	Cm_ += -Cm;
        std::ofstream outfile;
        outfile.open("Cl.txt", std::ios_base::app);
        outfile << t <<"    " <<  Cl_ << "\n";
        
	std::ofstream outfileCm;
        outfileCm.open("Cm.txt", std::ios_base::app);
        outfileCm << t <<"    " << Cm_ << "\n";
    }

    Info << cmpt << ": " << -Cf << endl;
}

inline void Foam::hisaModule::ibm()
{
//    fvMesh& mesh = mesh_();
    const labelList& cellMap = meshParts_[0].cellMap();

    deltaU_ = vectorField(mesh_->C().size(), vector(0, 0, 0));
    B_ = scalarField(N_L,0);
    X_ = scalarField(N_L,0);

    calcB(0, cellMap);
    solveX();
    calcCorrection(0, cellMap);

    calcB(1, cellMap);
    solveX();
    calcCorrection(1, cellMap);
}

inline vectorField Foam::hisaModule::getDeltaU()
{
    Info <<"IBM starts..." << endl;	
    updateLagrangian();
    ibm();
    return deltaU_;
}

// oscillation cylinder
//inline void Foam::hisaModule::updateLagrangian()
//{
//    const Time& runTime = time_;
//    scalar t = runTime.value();
//    for(int i = 0; i<N_L; i++)
//    {
//        // update velocity
//        Ub_[i].component(0) = 0;
//        Ub_[i].component(1) = 0.14*1.0472*cos(1.0472*t);
//    }
//}

inline void Foam::hisaModule::updateLagrangian()
{
    calcMotion();
    //locateCell();   // no need to update D and A matrices since there is no relative motion
}

// stationary cylinder
inline void Foam::hisaModule::calcMotion()
{

    const Time& runTime = time_;
    scalar t = runTime.value();
    scalar dt = runTime.deltaT().value();

    for(int i = 0; i<N_L; i++)
    {
	// update velocity
	Ub_[i].component(0) = 0;
	Ub_[i].component(1) = 0;
    }

    oversetMotion_[0].x() = 0;
    oversetMotion_[0].y() = 0;
    //oversetMotion_[0].y() = 0.5*sin(2*3.141592653*t/1);
    oversetMotion_[0].z() = 0;
}


// Swiming fish
// inline void Foam::hisaModule::calcMotion()
// {
//     const scalar f_swim = 2;
//     const scalar lambda = 1;

//     const Time& runTime = time_;
//     scalar t = runTime.value();
//     scalar dt = runTime.deltaT().value();

//     ut_ += Cd_*0.5*dt/0.12;
//     vt_ = 0;

//     scalar dxc, dyc;
//     dxc = ut_*dt;
//     dyc = vt_*dt;

//     xc_ += dxc;
//     yc_ += dyc;

//     std::ofstream outfile;
//     outfile.open("motion.txt", std::ios_base::app);
//     outfile << t <<"    " 
// 	    <<  xc_ <<"    "
// 	    <<  yc_ <<"    "
// 	    <<  ut_ <<"    "
// 	    <<  vt_ <<"    "
// 	    << "\n";

//     for(int i = 0; i<N_L; i++)
//     {
// 	scalar xLoc = solidObject_[i].x0;
// 	scalar A = 0.02 - 0.0825*xLoc + 0.1625*xLoc*xLoc;
// 	scalar h = A*cos(2*3.141592653*(xLoc/lambda-f_swim*t));
// 	// update location
// 	solidObject_[i].x = solidObject_[i].x0 ;
// 	solidObject_[i].y = solidObject_[i].y0  + h;
// 	// update velocity
// 	Ub_[i].component(0) = ut_;
// 	Ub_[i].component(1) = 2*3.141592653*f_swim*A*sin(2*3.141592653*(xLoc/lambda-f_swim*t));
// 	//// relative to overset
// 	//solidObject_[i].x = solidObject_[i].x0;
// 	//solidObject_[i].y = solidObject_[i].y0 + h;
// 	//// update velocity
// 	//Ub_[i].component(0) = 0;//ut_;
// 	//Ub_[i].component(1) = 2*3.141592653*f_swim*A*sin(2*3.141592653*(xLoc/lambda-f_swim*t));
//     }
//     oversetMotion_[0].x() = xc_;
//     oversetMotion_[0].y() = yc_;
//     oversetMotion_[0].z() = 0;
// }


// Sedimentation
//inline void Foam::hisaModule::calcMotion()
//{
//    const Time& runTime = time_;
//    scalar t = runTime.value();
//    scalar dt = runTime.deltaT().value();
//
//    scalar ax = (vt_-vt_o_)/dt;
//    vt_o_ = vt_;
//
//    ut_ += Cd_*0.5*dt/M_;
//    vt_ += (Cl_*0.5 - M_*g_*(1-rhoF_/rhoS_) + M_/1.25*ax)*dt/M_;
//    omega_ += 0;//Cm_*0.5*dt/I_;
//
//    scalar dxc, dyc, dTheta;
//    dxc = ut_*dt;
//    dyc = vt_*dt;
//    dTheta = omega_*dt;
//
//    xc_ += dxc;
//    yc_ += dyc;
//
//    scalar wall = -16.+0.5+3./2.*h_ ;
//   // if (yc_<=wall)
//   // {
//   //     //vt_ = -vt_;
//   //     //yc_ = yc_-dyc+vt_*dt; 
//   //     vt_ = 0;
//   //     yc_ = wall; 
//   // }
//
//    if (fabs(yc_-wall)>1e-4)
//    {
//        if (yc_<=wall)
//        {
//            xc_ -= dxc;
//            yc_ -= dyc;
//
//            dt = (wall-yc_)/vt_; 
//            dxc = ut_*dt;
//            dyc = vt_*dt;
//            dTheta = omega_*dt;
//
//            xc_ += dxc;
//            yc_ += dyc;
//        }
//    } else {
//        vt_ = 0;
//    }
//
//    std::ofstream outfile;
//    outfile.open("motion.txt", std::ios_base::app);
//    outfile << t <<"    " 
//	    <<  xc_ <<"    "
//	    <<  yc_ <<"    "
//	    <<  ut_ <<"    "
//	    <<  vt_ <<"    "
//	    << "\n";
//
//    for(int i = 0; i<N_L; i++)
//    {
//	// update location
//	solidObject_[i].x = solidObject_[i].x + dxc;
//	solidObject_[i].y = solidObject_[i].y + dyc;
//	// update velocity
//	Ub_[i].component(0) = ut_ + omega_*(solidObject_[i].y-yc_);
//	Ub_[i].component(1) = vt_ - omega_*(solidObject_[i].x-xc_);
//    }
//    oversetMotion_[0].x() = xc_;
//    oversetMotion_[0].y() = yc_;
//    oversetMotion_[0].z() = 0;
//}

//inline void Foam::hisaModule::locateCell()
//{
//    const fvMesh& mesh = meshParts_[0].subMesh();
//
//    for(int i = 0; i < N_L;i++)
//    {
//        D_[i].clear();
//        DCellID_[i].clear();
//        for(int j = 0; j<mesh.C().size(); j++)
//        {
//            double rx = fabs(mesh.C()[j].x()-solidObject_[i].x)/h_;
//            double ry = fabs(mesh.C()[j].y()-solidObject_[i].y)/h_;
//            if ( (rx<2) && (ry<2))
//            {
//                D_[i].append(deltaFunc(rx)*deltaFunc(ry));
//                DCellID_[i].append(j);
//            }
//        }
//    }
//
//}
inline void Foam::hisaModule::locateCell()
{
    const fvMesh& mesh = meshParts_[0].subMesh();

    for(int i = 0; i<N_L; i++)
    {
        D_[i].clear();
        DCellID_[i].clear();
	point location(0,0,0);
	location.x() = solidObject_[i].x;
	location.y() = solidObject_[i].y;
    
        //update host cells
	hostCell_[i] = findNearestCell(location,hostCell_[i]);
        //update matrices D and A
	findEulerian(hostCell_[i], i);
    }
}

inline label Foam::hisaModule::findNearestCell(const point& location, const label seedCelli)
{
//    fvMesh& mesh = mesh_();
    const fvMesh& mesh = meshParts_[0].subMesh();
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

inline bool Foam::hisaModule::findNearer
(
    const point& sample,
    const pointField& points,
    const labelList& indices,
    label& nearestI,
    scalar& nearestDistSqr
)
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

inline void Foam::hisaModule::findEulerian(const label hostCell, const label i)
{
    // host cell
    int count = 0;
    isNearby(hostCell, i, count);
    // west cells
    label cell1 = neighbourCells_[hostCell].w_;
    isNearby(cell1, i, count);
    label cell2 = neighbourCells_[cell1].w_;
    isNearby(cell2, i, count);
    // east cells
    label cell3 = neighbourCells_[hostCell].e_;
    isNearby(cell3, i, count);
    label cell4 = neighbourCells_[cell3].e_;
    isNearby(cell4, i, count);
    // north cells
    label cell5 = neighbourCells_[hostCell].n_;
    isNearby(cell5, i, count);
    label cell6 = neighbourCells_[cell5].n_;
    isNearby(cell6, i, count);
    // south cells
    label cell7 = neighbourCells_[hostCell].s_;
    isNearby(cell7, i, count);
    label cell8 = neighbourCells_[cell7].s_;
    isNearby(cell8, i, count);
    // northwest cells
    label cell9 = neighbourCells_[cell1].n_;
    isNearby(cell9, i, count);
    label cell10 = neighbourCells_[cell9].n_;
    isNearby(cell10, i, count);
    label cell11 = neighbourCells_[cell2].n_;
    isNearby(cell11, i, count);
    label cell12 = neighbourCells_[cell11].n_;
    isNearby(cell12, i, count);
    // southwest cells
    label cell13 = neighbourCells_[cell1].s_;
    isNearby(cell13, i, count);
    label cell14 = neighbourCells_[cell13].s_;
    isNearby(cell14, i, count);
    label cell15 = neighbourCells_[cell2].s_;
    isNearby(cell15, i, count);
    label cell16 = neighbourCells_[cell15].s_;
    isNearby(cell16, i, count);
    // northeast cells
    label cell17 = neighbourCells_[cell3].n_;
    isNearby(cell17, i, count);
    label cell18 = neighbourCells_[cell17].n_;
    isNearby(cell18, i, count);
    label cell19 = neighbourCells_[cell4].n_;
    isNearby(cell19, i, count);
    label cell20 = neighbourCells_[cell19].n_;
    isNearby(cell20, i, count);
    // southeast cells
    label cell21 = neighbourCells_[cell3].s_;
    isNearby(cell21, i, count);
    label cell22 = neighbourCells_[cell21].s_;
    isNearby(cell22, i, count);
    label cell23 = neighbourCells_[cell4].s_;
    isNearby(cell23, i, count);
    label cell24 = neighbourCells_[cell23].s_;
    isNearby(cell24, i, count);

//    Info << "number of eulerian cells" <<endl;
//    Info << count <<endl;
//    cin.get();
}

inline void Foam::hisaModule::isNearby(const label cell, const label i, int& count)
{
    const fvMesh& mesh = meshParts_[0].subMesh();
    const pointField& cellCentres = mesh.cellCentres();
    const point particle(solidObject_[i].x, solidObject_[i].y, 0);

    point dist = (cellCentres[cell] - particle)/h_;
    
    if ( (mag(dist.x())<2) && (mag(dist.y())<2) )
    {
        //D_[i][count] = deltaFunc(dist.x())*deltaFunc(dist.y());
        //DCellID_[i][count] = cell; 
        D_[i].append(deltaFunc(dist.x())*deltaFunc(dist.y()));
        DCellID_[i].append(cell);
        count++;
    }
}

/*inline scalar Foam::hisaModule::searchNear(const label cell)
{
    scalar minDist = GREAT;
    labelList& cellCells = mesh.cellCells(cell);

    forAll(cellCells, idx)
    {
        label cell = cellCells[idx];
        scalar dist = mag(points[cell] - );
	if (dist)
	{
	}
        scalar minDist = min(minDist,dist);
    }
    if (dist>2)
    {
        return 1;
    } else {
        forAll(cellCells, idx)
        {
            label cell = cellCells[idx];
            minDist = searchNear(cell);
        }
    }
}
*/

inline void Foam::hisaModule::postProcess(scalar& Cl_p, scalar& Cd_p, scalar& Cl_v, scalar& Cd_v)
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

