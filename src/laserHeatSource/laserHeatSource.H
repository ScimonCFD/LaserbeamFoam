/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    laserHeatSource

Description
    Calculate the deposition heat source field from a laser using ray-tracing.

    Laser parameters are defined via the provided dictionary.

    It is assumed that the laser enters the domain at the boundary with an
    outward-facing unit normal of (0 1 0).

    Laser rays are written as VTKs for visualisation.

Author
    Thomas Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    laserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef laserHeatSource_H
#define laserHeatSource_H

#include "volFields.H"
#include "dictionary.H"
#include "interpolationTable.H"
#include "OFstream.H"
#include "compactRay.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                        Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/


class laserHeatSource
:
    public IOdictionary
{
    // Private data

        // Heat source field
        volScalarField deposition_;

        //- Field indicating the boundary faces where the laser enters the
        //  domain
        volScalarField laserBoundary_;

        //- Field for tracking error codes
        volScalarField errorTrack_;

        //- Field for tracking the latest ray to vist each cell
        volScalarField rayNumber_;

        //- Field for tracking the energy deposited in each cell
        volScalarField rayQ_;

        //- y dimension field
        volScalarField yDim_;

        //- Powder sim flag
        //  Needed for cases with spherical particles or spurious currents will
        //  smear out the particle interface un-physically
        const Switch powderSim_;

        //- new laser discretisation flag
        //  Tell the heat source if you want to discreetise the traditional laserbeamfoam way or in radial and polar coordinates
        const Switch radialPolarHeatSource_;

        //- List of the laser names
        wordList laserNames_;

        //- Copy of the laser dicts
        PtrList<dictionary> laserDicts_;

        //- List of the time vs laser position interpolation series
        PtrList<interpolationTable<vector>> timeVsLaserPosition_;

        //- List of the time vs laser power interpolation series
        PtrList<interpolationTable<scalar>> timeVsLaserPower_;

        //- Global bounding box of the mesh
        boundBox globalBB_;


    // Private Member Functions

        //- Create the initial rays
        //  This involves creating a list of points (starting position of the
        //  laser rays) and their associated power (J/s)
        //  The rayCoords and rayPowers fields should be the same on all
        //  processors when returned from this function
        void createInitialRays
        (
            pointField& rayCoords,
            scalarField& rayPowers,
            const vectorField& CI,
            const vector& currentLaserPosition,
            const scalar laserRadius,
            const label N_sub_divisions,
            const label nRadial,
            const label nAngular,
            const vector& V_incident,
            const scalar Radius_Flavour,
            const scalar Q_cond,
            const scalar beam_radius
        ) const;

        //- Disallow default bitwise copy construct
        laserHeatSource(const laserHeatSource&);

        //- Disallow default bitwise assignment
        void operator=(const laserHeatSource&);

public:

    //- Runtime type information
    TypeName("laserHeatSource");

    // Constructors

        //- Construct from the mesh
        laserHeatSource
        (
            const fvMesh& mesh
        );

    // Destructor

        virtual ~laserHeatSource()
        {}


    // Member Functions

        //- Const access to the deposition field
        const volScalarField& deposition() const
        {
            return deposition_;
        }

        //- Non-const access to the yDim field
        volScalarField& yDim()
        {
            return yDim_;
        }

        //- Const access to the powderSim flag
        const Switch& powderSim() const
        {
            return powderSim_;
        }

        //- Const access to the HS discretisation flag
        const Switch& radialPolarHeatSource() const
        {
            return radialPolarHeatSource_;
        }

        //- Update the deposition field for all lasers
        //  The filtered volume frame field (alphaFiltered) and the filtered
        //  interface normal field are required to define the interface
        //  The electrical resistivity field is passed as an argument, allowing
        //  a general laser heat source for multi-component and temperature
        //  dependent resistivity solvers
        //  This resets the deposition field to zero before updating it
        void updateDeposition
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in
        );

        //- Update the deposition for the given laser parameters
        //  This adds to the existing deposition field and does not reset it to
        //  zero first
        void updateDeposition
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in,
            const word& laserName,
            const vector& currentLaserPosition,
            const scalar currentLaserPower,
            const scalar laserRadius,
            const label N_sub_divisions,
            const label nRadial,
            const label nAngular,
            const vector& V_incident,
            const scalar wavelength,
            const scalar e_num_density,
            const scalar dep_cutoff,
            const scalar Radius_Flavour,
            const Switch useLocalSearch,
            const label maxLocalSearch,
            const boundBox& globalBB  // FOR ACCESS TO BOUNDING BOX GLOBALLY
        );

        //- Const access to the global bounding box
        const boundBox& globalBB() const
        {
            return globalBB_;
        }

        // bool pointOnFace(const face& f, const point& p, const pointField& meshPoints);



// void writeLinesToVTK(const pointField& points, const word& filename)
// {
//     OFstream file(filename);
    
//     // VTK header
//     file << "# vtk DataFile Version 3.0" << nl;
//     file << "Line data" << nl;
//     file << "ASCII" << nl;
//     file << "DATASET POLYDATA" << nl;
    
//     // Points section
//     file << "POINTS " << points.size() << " float" << nl;
//     forAll(points, i)
//     {
//         const point& pt = points[i];
//         file << pt.x() << " " << pt.y() << " " << pt.z() << nl;
//     }
    
//     // Lines section (connects consecutive points)
//     label nLines = points.size() - 1;
//     file << "LINES " << nLines << " " << (nLines * 3) << nl;
    
//     for (label i = 0; i < nLines; i++)
//     {
//         file << "2 " << i << " " << (i + 1) << nl;
//     }
// }

// void writeLinesToVTK(const DynamicList<point>& points, const fileName& filename)
// {
//     OFstream file(filename);
    
//     if (!file.good())
//     {
//         FatalErrorInFunction
//             << "Cannot open file " << filename
//             << exit(FatalError);
//     }
    
//     // VTK header
//     file << "# vtk DataFile Version 3.0" << nl;
//     file << "Line data" << nl;
//     file << "ASCII" << nl;
//     file << "DATASET POLYDATA" << nl;
    
//     // Points section
//     file << "POINTS " << points.size() << " float" << nl;
//     forAll(points, i)
//     {
//         const point& pt = points[i];
//         file << pt.x() << " " << pt.y() << " " << pt.z() << nl;
//     }
    
//     // Lines section
//     label nLines = points.size() - 1;
//     file << "LINES " << nLines << " " << (nLines * 3) << nl;
    
//     for (label i = 0; i < nLines; i++)
//     {
//         file << "2 " << i << " " << (i + 1) << nl;
//     }
// }


void writeMultipleRaysToVTK(
    const DynamicList<DynamicList<point>>& rays,
    const fileName& filename)//const fileName& filename
{

    
    OFstream file(filename);
    
    if (!file.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << filename
            << exit(FatalError);
    }
    
    // Calculate total number of points and lines
    label totalPoints = 0;
    label totalLines = 0;
    
    forAll(rays, rayI)
    {
        totalPoints += rays[rayI].size();
        if (rays[rayI].size() > 1)
        {
            totalLines += (rays[rayI].size() - 1); // n-1 line segments per ray
        }
    }
    
    // VTK header
    file << "# vtk DataFile Version 3.0" << nl;
    file << "Multiple ray data" << nl;
    file << "ASCII" << nl;
    file << "DATASET POLYDATA" << nl;
    
    // Write all points
    file << "POINTS " << totalPoints << " float" << nl;
    forAll(rays, rayI)
    {
        const DynamicList<point>& ray = rays[rayI];
        forAll(ray, pointI)
        {
            const point& pt = ray[pointI];
            file << pt.x() << " " << pt.y() << " " << pt.z() << nl;
        }
    }
    
    // Write line connectivity
    // Each line segment is defined by 2 points
    file << "LINES " << totalLines << " " << (totalLines * 3) << nl;
    
    label pointOffset = 0;
    forAll(rays, rayI)
    {
        const DynamicList<point>& ray = rays[rayI];
        
        // Connect consecutive points within this ray only
        for (label i = 0; i < ray.size() - 1; i++)
        {
            file << "2 " << (pointOffset + i) << " " << (pointOffset + i + 1) << nl;
        }
        
        pointOffset += ray.size(); // Update offset for next ray
    }
    
    Info << "Written VTK file with " << rays.size() << " rays: " << filename << endl;
}


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
