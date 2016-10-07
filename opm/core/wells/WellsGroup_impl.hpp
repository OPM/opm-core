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
    void WellsGroup::updateWellTargets(const WellState& well_state)
    {
        // TODO: currently, we only handle the level of the well groups for the moment, i.e. the level just above wells
        // We believe the relations between groups are similar to the relations between different wells inside the same group.
        // While there will be somre more complication invloved for sure.
        // Basically, we need to update the target rates for the wells still under group control.

        // What facility need to be put in?
        // update the production control
        // What we need there? 1. The target. 2. The control mode. 3. The share produced by the wells under individual control
        //                     4. The guide rate.
        // control mode
        ProductionSpecification::ControlMode control_mode = prodSpec().control_mode_;
        double target_rate = prodSpec().liquid_max_rate_;
        std::cout << " control mode is " << ProductionSpecification::toString(control_mode);
        std::cout << " rate is " << target_rate << std::endl;

        if (ProductionSpecification::toString(control_mode) == "FLD") {
            auto* parent_node = getParent();
            control_mode = parent_node->prodSpec().control_mode_;
            target_rate = parent_node->prodSpec().liquid_max_rate_;
            std::cout << " control mode is " << ProductionSpecification::toString(control_mode);
            std::cout << " rate is " << target_rate << std::endl;
        }

        double rate_individual_control = 0;

        for (size_t i = 0; i < children_.size(); ++i) {
            if (children_[i]->individualControl()) {
                // get the rate here.
                const std::string well_name = children_[i]->name();
                std::cout << "well_name " << well_name;
                typedef typename WellState::WellMapType::const_iterator const_iterator;
                const_iterator it = well_state.wellMap().find(well_name);
                const int well_index = (*it).second[0];
                const int np = well_state.numPhases();

                const double oil_rate = well_state.wellRates()[np * well_index + 1];
                const double water_rate = well_state.wellRates()[np * well_index];
                std::cout << " oil_rate is " << oil_rate << " water_rate is " << water_rate << " liquid rate is " << oil_rate + water_rate << std::endl;
                rate_individual_control += std::abs(oil_rate + water_rate);
            }
        }

        const double rate_for_group_control = target_rate - rate_individual_control;
        std::cout << " rate_for_group_control " << rate_for_group_control << std::endl;

        const double my_guide_rate = productionGuideRate(true);
        std::cout << " my_guide_rate is " << my_guide_rate << " when updateWellTargets " << std::endl;
        if (my_guide_rate == 0) {
            std::cout << " something wrong with the my_guide_rate, need to check, it should have come here " << std::endl;
            std::cin.ignore();
        }

        for (size_t i = 0; i < children_.size(); ++i) {
            // if (children_[i]->shouldUpdateWellTargets() && !children_[i]->individualControl()) {
            if (!children_[i]->individualControl()) {
                const double children_guide_rate = children_[i]->productionGuideRate(true);
                std::cout << " well_name " << children_[i]->name() << " children_guide_rate " << children_guide_rate << " my_guide_rate " << my_guide_rate << std::endl;
                std::cout << " new target rate is " << (children_guide_rate/my_guide_rate) * rate_for_group_control * 86400/0.1590 << std::endl;
                // children_[i]->applyProdGroupControl(control_mode, (children_guide_rate/my_guide_rate) * rate_for_group_control, false);
                children_[i]->applyProdGroupControl(control_mode, (children_guide_rate/my_guide_rate) * rate_for_group_control, true);
                children_[i]->setShouldUpdateWellTargets(false);
            }
        }
        std::cin.ignore();
        // liquid_max_rate

        // update the injectioin control
    }

}
