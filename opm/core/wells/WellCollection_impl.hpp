/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/core/wells/WellCollection.hpp>


namespace Opm
{


    template <class WellState>
    void WellCollection::updateWellTargets(const WellState& well_state)
    {
        // TODO: currently, we only handle the level of the well groups for the moment, i.e. the level just above wells
        // We believe the relations between groups are similar to the relations between different wells inside the same group.
        // While there will be somre more complication invloved for sure.
        for (size_t i = 0; i < leaf_nodes_.size(); ++i) {
            // find a node needs to update targets, then update targets for all the wellls inside the group.
            // if (leaf_nodes_[i]->shouldUpdateWellTargets() && !leaf_nodes_[i]->individualControl()) {
            if (!leaf_nodes_[i]->individualControl()) {
                // TODO: will remove dynamic_cast with interface revision.
                WellsGroup* parent_node = dynamic_cast<Opm::WellsGroup *>(leaf_nodes_[i]->getParent());
                // update the target within this group.
                parent_node->updateWellTargets(well_state);
            }
        }

        setJustUpdateWellTargets(true);
    }

}
