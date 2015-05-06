/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#ifndef __PARTITION_REFINER_HPP__
#define __PARTITION_REFINER_HPP__

#include <vector>
#include "LookUpTableUtils.hpp"

class PartitionRefiner {
public:
        PartitionRefiner(LUTTools::KmerID numElem);
        void splitWith(std::vector<LUTTools::KmerID> splitter);
        void relabel();
        const std::vector<LUTTools::KmerID>& partitionMembership();

private:
        LUTTools::KmerID numElem_;
        std::vector<LUTTools::KmerID> membership_;
        LUTTools::KmerID maxSetIdx_;
};

#endif // __PARTITION_REFINER_HPP__
