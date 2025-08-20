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
    compactRay

Description
    Class for storing information about individual rays propagating through
    a polMesh in a ray tracing algorithm.

Author
    Thomas Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    compactRay.C

\*---------------------------------------------------------------------------*/

#ifndef compactRay_H
#define compactRay_H

#include "volFields.H"
#include "dictionary.H"
#include "interpolationTable.H"

// #include <vector>
// #include "globalIndex.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class compactRay Declaration
\*---------------------------------------------------------------------------*/


class compactRay
{
    public:

        // Public data

        point origin_;
        vector direction_;
        scalar power_;
        scalar remainingPower_;
        label currentCell_;
        label bounceCount_;
        label global_Ray_number_;
        bool active_;
        DynamicList<point> path_;

    // Public functions

        //- Empty constructor
        compactRay()
        :
            power_(0),
            remainingPower_(0),
            currentCell_(-1), 
            bounceCount_(0),
            global_Ray_number_(-1),
            active_(false)
        {}

        //- Construct from position, direction and power
        compactRay(const point& origin, const vector& dir, scalar power)
        :
            origin_(origin),
            direction_(dir),
            power_(power), 
            remainingPower_(power),
            currentCell_(-1),
            bounceCount_(0),
            global_Ray_number_(-1),
            active_(true)
        {
            direction_ /= mag(direction_) + VSMALL;
            path_.clear();
        }

        // Destructor

            virtual ~compactRay()
            {}


        //- Reset the position, direction and power
        void reset(const point& origin, const vector& dir, scalar power)
        {
            origin_ = origin;
            direction_ = dir / (mag(dir) + VSMALL);
            power_ = power;
            remainingPower_ = power;
            currentCell_ = -1;
            bounceCount_ = 0;
            global_Ray_number_ = -1;
            active_ = true;
        }

        //- Copy
        bool operator==(const compactRay& other) const
        {
            return 
            origin_ == other.origin_ && 
            direction_ == other.direction_ && 
            power_ == other.power_ &&
            global_Ray_number_ == other.global_Ray_number_;
        }

        //- Compare
        bool operator!=(const compactRay& other) const
        {
            return !(*this == other);
        }

        // Compare
        bool operator<(const compactRay& other) const
        {
            if (origin_ != other.origin_)
            {
                return origin_ < other.origin_;
            }

            if (direction_ != other.direction_)
            {
                return direction_ < other.direction_;
            }

            return power_ < other.power_;
        }

        //- Output stream operator
        friend Ostream& operator<<(Ostream& os, const compactRay& ray)
        {
            return os  << ray.origin_
                       << ray.direction_
                       << ray.power_
                       << ray.global_Ray_number_;
        }
    
        //- Input stream operator
        friend Istream& operator>>(Istream& is, compactRay& ray)
        {
            return is  >> ray.origin_
                       >> ray.direction_
                       >> ray.power_
                       >> ray.global_Ray_number_;
        }
};


/*---------------------------------------------------------------------------*\
                        Helper functions and structs
\*---------------------------------------------------------------------------*/

struct combineRayLists
{
    void operator()(DynamicList<compactRay>& x, const DynamicList<compactRay>& y) const
    {
        // Create a hash set of existing items for O(1) lookup
        HashSet<compactRay> existing;
        forAll(x, i)
        {
            existing.insert(x[i]);
        }
        
        // Add only new items
        forAll(y, i)
        {
            if (!existing.found(y[i]))
            {
                x.append(y[i]);
            }
        }
    }
};

struct combineLists
{
    void operator()(DynamicList<point>& x, const DynamicList<point>& y) const
    {
        x.append(y);
    }
};

struct combineRayPaths
{
    void operator()(DynamicList<DynamicList<point>>& x, const DynamicList<DynamicList<point>>& y) const
    {
        forAll(y, i)
        {
            x.append(y[i]);
        }
    }
};

    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
