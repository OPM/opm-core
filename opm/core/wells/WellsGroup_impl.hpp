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

#include <opm/core/wells/WellsGroup.hpp>
#include <iostream>


namespace Opm
{


    template <class WellState>
    void WellsGroup::updateWellProductionTargets(const WellState& well_state)
    {
        // TODO: currently, we only handle the level of the well groups for the moment, i.e. the level just above wells
        // We believe the relations between groups are similar to the relations between different wells inside the same group.
        // While there will be somre more complication invloved for sure.
        // Basically, we need to update the target rates for the wells still under group control.

        ProductionSpecification::ControlMode control_mode = prodSpec().control_mode_;
        double target_rate = prodSpec().liquid_max_rate_;

        if (ProductionSpecification::toString(control_mode) == "FLD") {
            auto* parent_node = getParent();
            control_mode = parent_node->prodSpec().control_mode_;
            target_rate = parent_node->prodSpec().liquid_max_rate_;
        }

        double rate_individual_control = 0;

        for (size_t i = 0; i < children_.size(); ++i) {
            if (children_[i]->individualControl() && std::dynamic_pointer_cast<Opm::WellNode>(children_[i])->isProducer()) {
                // get the rate here.
                const std::string well_name = children_[i]->name();
                typedef typename WellState::WellMapType::const_iterator const_iterator;
                const_iterator it = well_state.wellMap().find(well_name);
                const int well_index = (*it).second[0];
                const int np = well_state.numPhases();

                const double oil_rate = well_state.wellRates()[np * well_index + 1];
                const double water_rate = well_state.wellRates()[np * well_index];
                rate_individual_control += std::abs(oil_rate + water_rate);
            }
        }

        const double rate_for_group_control = target_rate - rate_individual_control;

        const double my_guide_rate = productionGuideRate(true);

        for (size_t i = 0; i < children_.size(); ++i) {
            // if (children_[i]->shouldUpdateWellTargets() && !children_[i]->individualControl()) {
            if (!children_[i]->individualControl() && std::dynamic_pointer_cast<Opm::WellNode>(children_[i])->isProducer()) {
                const double children_guide_rate = children_[i]->productionGuideRate(true);
                // children_[i]->applyProdGroupControl(control_mode, (children_guide_rate/my_guide_rate) * rate_for_group_control, false);
                children_[i]->applyProdGroupControl(control_mode, (children_guide_rate/my_guide_rate) * rate_for_group_control, true);
                children_[i]->setShouldUpdateWellTargets(false);
            }
        }
    }






    template <class WellState>
    void WellsGroup::updateWellInjectionTargets(const WellState&) {
        // NOT doing anything yet.
        // Will finish it when having an examples with more than one injection wells within same injection group.
        for (size_t i = 0; i < children_.size(); ++i) {
            if (!children_[i]->individualControl() && std::dynamic_pointer_cast<Opm::WellNode>(children_[i])->isInjector()) {
                children_[i]->setShouldUpdateWellTargets(false);
            }
        }

    }

}
