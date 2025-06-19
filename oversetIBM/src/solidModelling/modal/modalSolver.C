/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2017 Ridhwaan Suliman - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2011-2013 OpenFOAM Foundation

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

#include "modalSolver.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wordReList.H"

#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

#include "PrimitivePatchInterpolation.H"
#include "surfMesh.H"

#include "modalMotionPointPatchVectorField.H"

#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(modalSolver, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        modalSolver,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::modalSolver::writeFileHeader(const label i)
{
    file() << "# Time" << tab;
    for (int i = 0; i < nModes_; i++)
    {
        file() << "modalCoordinate" << i << tab;
    }
    for (int i = 0; i < nModes_; i++)
    {
        file() << "modalVelocity" << i << tab;
    }
    for (int i = 0; i < nModes_; i++)
    {
        file() << "generalisedForce" << i << tab;
    }
    for (int i = 0; i < nModes_; i++)
    {
        file() << "generalisedForceLin" << i << tab;
    }
    for (int i = 0; i < nModes_; i++)
    {
        for (int j = i; j < nModes_; j++)
        {
            file() << "generalisedForceNonlinCoeff" << i << j << tab;
        }
    }

    file()<< endl;
}


Foam::tmp<Foam::volSymmTensorField> Foam::functionObjects::modalSolver::devRhoReff(const objectRegistry& obr) const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (obr.foundObject<cmpTurbModel>(cmpTurbModel::typeName))
    {
        const cmpTurbModel& turb =
            obr.lookupObject<cmpTurbModel>(cmpTurbModel::typeName);

        return turb.devRhoReff();
    }
    else if (obr.foundObject<icoTurbModel>(icoTurbModel::typeName))
    {
        const incompressible::turbulenceModel& turb =
            obr.lookupObject<icoTurbModel>(icoTurbModel::typeName);

        return rho(obr)*turb.devReff();
    }
    else if (obr.foundObject<fluidThermo>(fluidThermo::typeName))
    {
        const fluidThermo& thermo =
            obr.lookupObject<fluidThermo>(fluidThermo::typeName);

        const volVectorField& U = obr.lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if
    (
        obr.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            obr.lookupObject<transportModel>("transportProperties");

        const volVectorField& U = obr.lookupObject<volVectorField>(UName_);

        return -rho(obr)*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (obr.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             obr.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

        const volVectorField& U = obr.lookupObject<volVectorField>(UName_);

        return -rho(obr)*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorIn("forces::devRhoReff()")
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}

Foam::tmp<Foam::volScalarField> Foam::functionObjects::modalSolver::rho(const objectRegistry& obr) const
{
    if (rhoName_ == "rhoInf")
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr);

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(obr.lookupObject<volScalarField>(rhoName_));
    }
}


Foam::scalar Foam::functionObjects::modalSolver::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorIn("modalSolver::rho(const volScalarField& p)")
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::modalSolver::modalSolver
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    writeFiles(name, runTime, dict, name),
    name_(name),
    time_(time),
    dict_(dict),
    outputControl_(time, dict, name),
    active_(true),
    log_(false),
    nModes_(0),
    steadyState_(false),
    massMatrix_(0),
    dampingMatrix_(0),
    stiffnessMatrix_(0),
    u_(0),
    uprev_(0),
    v_(0),
    vprev_(0),
    modalVars_(IOobject
                (
                    "modalVariables",
                    time_.timeName(),
                    "uniform",
                    time_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                )
              ),
    lastTimeIndex_(0),
    calledForPatch_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(VGREAT),
    pRef_(0),
    genForcesConst_(0),
    genForcesLinCoeffs_(0),
    genForcesConst0_(0),
    genForcesLinCoeffs0_(0),
    forcesRelax_(0),
    forcesLinCoeffsRelax_(0),
    firstIteration_(true),
    secondIteration_(false),
    forcesIncrementPrevIter_(0),
    forcesLinCoeffsIncrementPrevIter_(0),
    forcesRelaxCurr_(0),
    forcesLinCoeffsRelaxCurr_(0)
{
    readDict(dict);
    resetNames(this->typeName());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::modalSolver::~modalSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::modalSolver::on()
{
    active_ = true;
}

void Foam::functionObjects::modalSolver::off()
{
    active_ = false;
}

bool Foam::functionObjects::modalSolver::start()
{
    readDict(dict_);
    return true;
}

bool Foam::functionObjects::modalSolver::read(const dictionary& dict)
{
    if (dict != dict_)
    {
        dict_ = dict;
        outputControl_.read(dict);
        return start();
    }
    else
    {
        return false;
    }
}

void Foam::functionObjects::modalSolver::readDict(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);

        steadyState_ = dict.lookupOrDefault<Switch>("steadyState", false);

        if (!dict.readIfPresent("stiffnessMatrix", stiffnessMatrix_))
        {
            FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                << "'stiffnessMatrix' must be specified." << exit(FatalIOError);
        }
        nModes_ = round(sqrt(scalar(stiffnessMatrix_.size())));
        if (nModes_*nModes_ != stiffnessMatrix_.size())
        {
            FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                << "Stiffness matrix must be square." << exit(FatalIOError);
        }

        if (!steadyState_)
        {
            if (!dict.readIfPresent("massMatrix", massMatrix_))
            {
                FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                    << "'massMatrix' must be specified." << exit(FatalIOError);
            }
            if (nModes_*nModes_ != massMatrix_.size())
            {
                FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                    << "Mass matrix must be square and same dimensions as stiffness matrix." << exit(FatalIOError);

            }
            if (!dict.readIfPresent("dampingMatrix", dampingMatrix_))
            {
                FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                    << "'dampingMatrix' must be specified." << exit(FatalIOError);
            }
            if (nModes_*nModes_ != dampingMatrix_.size())
            {
                FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                    << "Damping matrix must be square and same dimensions as stiffness matrix." << exit(FatalIOError);
            }
        }

        u_ = modalVars_.lookupOrAddDefault<scalarField>("coordinates", scalarField(nModes_, 0.0));
        v_ = modalVars_.lookupOrAddDefault<scalarField>("velocities", scalarField(nModes_, 0.0));
        if (u_.size() != nModes_ || v_.size() != nModes_)
        {
            FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                << "Number of initial modal coordinates and/or velocities "
                << "not compatible with dimensions of stiffness matrix." << exit(FatalIOError);
        }
        uprev_ = u_;
        vprev_ = v_;

        directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

        patchSet_.clear();
        regions_.clear();
        const entry* regionsEntryPtr = dict.lookupEntryPtr("regions", false, false);
        if (regionsEntryPtr && regionsEntryPtr->isDict())
        {
            const dictionary& subDicts = regionsEntryPtr->dict();
            forAllConstIter(dictionary, subDicts, iter)
            {
                if (!iter().isDict())
                {
                    FatalIOErrorIn("void Foam::functionObjects::modalSolver::read(const dictionary& dict)",dict)
                        << "Unexpected entry inside 'region' dictionary." << exit(FatalIOError);
                }
                const word& region = iter().keyword();
                const dictionary& subdict = iter().dict();

                const fvMesh& mesh = time_.lookupObject<fvMesh>(region);
                const polyBoundaryMesh& pbm = mesh.boundaryMesh();

                patchSet_.append(pbm.patchSet(wordReList(subdict.lookup("patches"))));
                regions_.append(region);
            }
        }
        else
        {
            const fvMesh& mesh = time_.lookupObject<fvMesh>(polyMesh::defaultRegion);
            const polyBoundaryMesh& pbm = mesh.boundaryMesh();
            patchSet_.append(pbm.patchSet(wordReList(dict.lookup("patches"))));
            regions_.append(polyMesh::defaultRegion);
        }

        if (directForceDensity_)
        {
            // Optional entry for fDName
            fDName_ = dict.lookupOrDefault<word>("fDName", "fD");
        }
        else
        {
            // Optional entries U and p
            pName_ = dict.lookupOrDefault<word>("pName", "p");
            UName_ = dict.lookupOrDefault<word>("UName", "U");
            rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");

            if (rhoName_ == "rhoInf")
            {
                // Reference density needed since no rho found
                rhoRef_ = readScalar(dict.lookup("rhoInf"));
            }

            // Reference pressure, 0 by default
            pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
        }

        calledForPatch_.resize(patchSet_.size());
        forAll(patchSet_, iRegion)
        {
            forAllConstIter(labelHashSet, patchSet_[iRegion], iter)
            {
                calledForPatch_[iRegion].insert(iter.key(), true); // Must recalculate
            }
            calledForPatch_[iRegion].shrink();
        }

        genForcesConst_.resize(nModes_, 0.0);
        genForcesLinCoeffs_.resize((nModes_+1)*nModes_/2, 0.0);
        genForcesConst0_ = modalVars_.lookupOrAddDefault<scalarField>("oldGenForcesConst", genForcesConst_);
        genForcesLinCoeffs0_ = modalVars_.lookupOrAddDefault<scalarField>("oldGenForcesLinCoeffs", genForcesLinCoeffs_);

        ode_.initialise(nModes_, steadyState_, &massMatrix_, &dampingMatrix_, &stiffnessMatrix_, &genForcesConst_, &genForcesLinCoeffs_, &genForcesConst0_, &genForcesLinCoeffs0_, &time_);
        if (!steadyState_)
        {
            dictionary dict;
            dict.add("solver", "RKCK45");
            odeSolver_ = ODESolver::New(ode_, dict); // Make Runge-Kutta solver
        }

        forcesRelax_.resize(nModes_, dict.lookupOrDefault<scalar>("initialForceRelaxation", 1.0));
        forcesLinCoeffsRelax_.resize((nModes_+1)*nModes_/2, dict.lookupOrDefault<scalar>("initialForceRelaxation", 1.0));
        forcesRelaxCurr_ = forcesRelax_;
        forcesLinCoeffsRelaxCurr_ = forcesLinCoeffsRelax_;

        overRelaxLimit_ = dict.lookupOrDefault<scalar>("overRelaxLimit", 10.0);
        negOverRelaxLimit_ = dict.lookupOrDefault<scalar>("negOverRelaxLimit", -1.0);

    }
}


bool Foam::functionObjects::modalSolver::execute(const bool forceWrite)
{
    // Called at end of time step

    if (active_)
    {
        if (forceWrite || outputControl_.output())
        {
            write();
        }
    }
    return true;
}


bool Foam::functionObjects::modalSolver::end()
{
    // Called at end of time loop

    if (active_)
    {
        if (outputControl_.output())
        {
            write();
        }
    }

    return true;
}

bool Foam::functionObjects::modalSolver::timeSet()
{
    // Called at beginning of timestep

    return true;
}


void Foam::functionObjects::modalSolver::write()
{
    if (!active_)
    {
        return;
    }

    if (Pstream::master())
    {
        scalarField genForces(genForcesConst_);
        for (int i = 0; i < nModes_; i++)
        {
            for (int j = 0; j < nModes_; j++)
            {
                const label k = (j >= i ? (2*nModes_ - i + 1)*i/2 + j - i : (2*nModes_ - j + 1)*j/2 + i - j);
                genForces[i] += genForcesLinCoeffs_[k]*u_[j];
            }
        }

        functionObjectFile::write();

        if (log_)
        {
            Info<< type() << " output:" << nl
                << "    modal co-ordinates = (";
        }
        file() << time_.value() << tab;

        for (int i = 0; i < nModes_; i++)
        {
            if (log_) Info << u_[i];
            file() << u_[i] << tab;
            if (log_ && i < nModes_-1)
            {
                Info << " ";
            }
        }

        if (log_)
        {
            Info<< ")" << nl
                << "    modal velocities = (";
        }
        for (int i = 0; i < nModes_; i++)
        {
            if (log_) Info << v_[i];
            file() << v_[i] << tab;
            if (log_ && i < nModes_-1)
            {
                Info << " ";
            }
        }

        if(log_)
        {
            Info<< ")" << nl
                << "    generalised forces = (";
        }
        for (int i = 0; i < nModes_; i++)
        {
            if (log_) Info << genForces[i];
            file() << genForces[i] << tab;
            if (log_ && i < nModes_-1)
            {
                Info << " ";
            }
        }
        for (int i = 0; i < nModes_; i++)
        {
            file() << genForcesConst_[i] << tab;
        }
        for (int i = 0; i < nModes_; i++)
        {
            for (int j = i; j < nModes_; j++)
            {
                const label k = (j >= i ? (2*nModes_ - i + 1)*i/2 + j - i : (2*nModes_ - j + 1)*j/2 + i - j);
                file() << genForcesLinCoeffs_[k] << tab;
            }
        }

        if(log_) Info<< ")" << nl << endl;
        file() << endl;

    }
}

void Foam::functionObjects::modalSolver::updateGenForces()
{
    // Save current values for Aitken relaxation
    scalarField forcesPrevIter;
    scalarField forcesLinCoeffsPrevIter;
    if (firstIteration_)
    {
        // Current values already estimated based on forward-projection
        forcesPrevIter = genForcesConst0_;
        forcesLinCoeffsPrevIter = genForcesLinCoeffs0_;
    }
    else
    {
        forcesPrevIter = genForcesConst_;
        forcesLinCoeffsPrevIter = genForcesLinCoeffs_;
    }

    genForcesConst_ = 0.0;
    genForcesLinCoeffs_ = 0.0;

    forAll(patchSet_, iRegion)
    {
        const fvMesh& mesh = time_.lookupObject<fvMesh>(regions_[iRegion]);

        forAllConstIter(labelHashSet, patchSet_[iRegion], iter)
        {
            label patchI = iter.key();

            vectorField force;

            if (directForceDensity_)
            {
                const volVectorField& fD = mesh.lookupObject<volVectorField>(fDName_);

                const surfaceVectorField::GeometricBoundaryField& Sfb =
                    mesh.Sf().boundaryField();

                scalarField sA(mag(Sfb[patchI]));

                force = sA*fD.boundaryField()[patchI];
            }
            else
            {
                const volScalarField& p = mesh.lookupObject<volScalarField>(pName_);

                const surfaceVectorField::GeometricBoundaryField& Sfb =
                    mesh.Sf().boundaryField();

                tmp<volSymmTensorField> tdevRhoReff = devRhoReff(mesh);
                const volSymmTensorField::GeometricBoundaryField& devRhoReffb
                    = tdevRhoReff().boundaryField();

                // Scale pRef by density for incompressible simulations
                scalar pRef = pRef_/rho(p);

                force =
                    rho(p)*Sfb[patchI]*(p.boundaryField()[patchI] - pRef)
                    + (Sfb[patchI] & devRhoReffb[patchI]);

            }

            // Modes are point fields while forces are face quantities.
            // Interpolate modes to face positions before dotting.
            PrimitivePatchInterpolation<primitivePatch> interp(mesh.C().boundaryField()[patchI].patch().patch());
            const pointPatchField<vector>& ppf = mesh.lookupObject<pointVectorField>("pointDisplacement").boundaryField()[patchI];
            try
            {
                const modalMotionPointPatchVectorField& mmPatch = dynamic_cast<const modalMotionPointPatchVectorField&>(ppf);
                const List<vectorField>& modeShapes = mmPatch.getModeShapes();
                if (modeShapes.size() != nModes_ + (nModes_+1)*nModes_/2)
                {
                    FatalErrorIn("void Foam::functionObjects::modalSolver::updateCoeffs()")
                        << "Number of modes in function object '" << this->name_ << "'"
                        << " does not match number of mode shapes in "
                        << "'" << modalMotionPointPatchVectorField::typeName_() << "' patch '"
                        << mesh.boundaryMesh()[patchI].name() << "'." << exit(FatalError);
                }

                for(label iMode = 0; iMode < nModes_; iMode++) // Only linear modes
                {
                    vectorField faceMode = interp.pointToFaceInterpolate(modeShapes[iMode]);
                    genForcesConst_[iMode] += gSum(faceMode & force);
                }
                for(label iMode = 0; iMode < nModes_; iMode++) // Quadratic
                {
                    for (label jMode = iMode; jMode < nModes_; jMode++) // Only doing half as symmetric.
                    {
                        const int N = nModes_;
                        vectorField faceMode = interp.pointToFaceInterpolate(modeShapes[N+(2*N-iMode+1)*iMode/2+jMode-iMode]);
                        genForcesLinCoeffs_[(2*N-iMode+1)*iMode/2+jMode-iMode] += gSum(faceMode & force);
                    }
                }
            }
            catch(const std::bad_cast&)
            {
                FatalErrorIn("void Foam::functionObjects::modalSolver::updateCoeffs()")
                    << "pointDisplacement patch '" << mesh.boundaryMesh()[patchI].name() << "' must be of type '"
                    << modalMotionPointPatchVectorField::typeName_() << "'." << exit(FatalError);
            }

        }
    }

    // Relax accelerations
    if (firstIteration_ || secondIteration_)
    {
        forcesRelaxCurr_ = forcesRelax_;
        forcesLinCoeffsRelaxCurr_ = forcesLinCoeffsRelax_;
        if (firstIteration_)
        {
            firstIteration_ = false;
            secondIteration_ = true;
        }
        else if (secondIteration_)
        {
            secondIteration_ = false;
        }
    }
    else
    {

        // Aitken acceleration following the method of Kuttler and Wall,
        // Comput. Mech. (2008) 43:61-72
        // Note this is the scalar version
        forcesRelaxCurr_ *= forcesIncrementPrevIter_/stabilise(forcesIncrementPrevIter_-(genForcesConst_-forcesPrevIter), VSMALL);
        forcesLinCoeffsRelaxCurr_ *= forcesLinCoeffsIncrementPrevIter_/stabilise(forcesLinCoeffsIncrementPrevIter_-(genForcesLinCoeffs_-forcesLinCoeffsPrevIter), VSMALL);

        // Limit over-relaxation; allow negative values for removal of linear growth
        forcesRelaxCurr_ = min(forcesRelaxCurr_, overRelaxLimit_);
        forcesRelaxCurr_ = max(forcesRelaxCurr_, negOverRelaxLimit_);
        forcesLinCoeffsRelaxCurr_ = min(forcesLinCoeffsRelaxCurr_, overRelaxLimit_);
        forcesLinCoeffsRelaxCurr_ = max(forcesLinCoeffsRelaxCurr_, negOverRelaxLimit_);

    }

    forcesIncrementPrevIter_ = genForcesConst_ - forcesPrevIter;
    forcesLinCoeffsIncrementPrevIter_ = genForcesLinCoeffs_ - forcesLinCoeffsPrevIter;

    if (log_)
    {
        Info << "Aitken relaxation parameters: " << forcesRelaxCurr_ << " " << forcesLinCoeffsRelaxCurr_ << endl;
    }

    genForcesConst_ = forcesRelaxCurr_*genForcesConst_ + (1 - forcesRelaxCurr_)*forcesPrevIter;
    genForcesLinCoeffs_ = forcesLinCoeffsRelaxCurr_*genForcesLinCoeffs_ + (1 - forcesLinCoeffsRelaxCurr_)*forcesLinCoeffsPrevIter;

}

void Foam::functionObjects::modalSolver::updateCoeffs()
{

    updateGenForces();

    // Is this a new timestep?
    if (time_.timeIndex() != lastTimeIndex_)
    {
        lastTimeIndex_ = time_.timeIndex();
        uprev_ = u_;
        vprev_ = v_;
        //These are set here (after updateGenForces) because this function
        //is called at beginning of a timestep (mesh is moved before fluid solve)
        scalarField fc0 = genForcesConst0_;
        scalarField fl0 = genForcesLinCoeffs0_;
        genForcesConst0_ = genForcesConst_;
        genForcesLinCoeffs0_ = genForcesLinCoeffs_;

        //Forward project to get current values initial estimate
        scalar deltaT0 = time_.deltaT0Value();
        scalar deltaT = time_.deltaTValue();
        genForcesConst_ += (genForcesConst_-fc0)/deltaT0*deltaT;
        genForcesLinCoeffs_ += (genForcesLinCoeffs_-fl0)/deltaT0*deltaT;

        firstIteration_ = true;
        secondIteration_ = false;

    }

    // Solve modal ODEs and use quadratic mode shapes to get updated positions

    scalarField y(vprev_);
    y.append(uprev_);
    if(!steadyState_)
    {
        scalarField yerr(ode_.nEqns());
        scalar hEst = time_.deltaTValue();
        odeSolver_->solve(time_.timeOutputValue()-time_.deltaTValue(), time_.timeOutputValue(), y, hEst);
    }
    else
    {
        // In steady state case, just assume all derivatives are zero and invert modal equations directly
        scalarField yin(y);
        ode_.derivatives(time_.timeOutputValue(), yin, y);
    }
    for (int i = 0; i < nModes_; i++)
    {
        v_[i] = y[i];
        u_[i] = y[nModes_+i];
    }

    // Update dictionary entries for writing
    modalVars_.set("coordinates", u_);
    modalVars_.set("velocities", v_);

}

void Foam::functionObjects::modalSolver::getPatchDispl(const word& region, const label& patchId, const List<vectorField>& modeShapes, vectorField& patchDispl)
{
    forAll(regions_, iRegion)
    {
        if (regions_[iRegion] == region || region == "")
        {
            if (calledForPatch_[iRegion][patchId]) // We've already been called for this patch; can't use cached data, must update.
            {
                updateCoeffs();
                forAll(regions_, jRegion)
                {
                    forAllConstIter(labelHashSet, patchSet_[jRegion], iter)
                    {
                        calledForPatch_[jRegion][iter.key()] = false; // Can use cached data for all other patches
                    }
                }
            }
            calledForPatch_[iRegion][patchId] = true;
        }
    }

    if (modeShapes.size() != nModes_+(nModes_+1)*nModes_/2)
    {
        FatalErrorIn("void Foam::functionObjects::modalSolver::getPatchDispl(...)")
            << "Number of modes in function object '" << this->name_ << "'"
            << " does not match number of mode shapes in patch "
            << patchId << abort(FatalError);
    }


    // Calculate displacement of patch using quadratic mode shapes
    patchDispl = vector(0,0,0);
    for(label iMode = 0; iMode < nModes_; iMode++) // Linear modes
    {
        patchDispl += u_[iMode]*modeShapes[iMode];
    }
    for(label iMode = 0; iMode < nModes_; iMode++) // Quadratic
    {
        for (label jMode = iMode; jMode < nModes_; jMode++) // Only doing half to avoid repetition but adding conditional factor of 2 below.
        {
            const int N = nModes_;
            patchDispl += 0.5*(iMode == jMode ? 1 : 2)*u_[iMode]*u_[jMode]*modeShapes[N+(2*N-iMode+1)*iMode/2+jMode-iMode];
        }
    }

}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //defineNamedTemplateTypeNameAndDebug(modalSolverFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        modalSolver,
        dictionary
    );
}

// ************************************************************************* //
