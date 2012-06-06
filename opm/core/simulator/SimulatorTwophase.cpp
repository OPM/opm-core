/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/simulator/SimulatorTwophase.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

namespace Opm
{

    class SimulatorTwophase::Impl
    {
    public:
        Impl() {}
        void init(const parameter::ParameterGroup& param);
        void run(const SimulatorTimer& timer, const Wells& wells, TwophaseState& state, WellState& well_state);
    private:
#if 0
        bool output_;
        std::string output_dir_;
        int output_interval;
        int num_transport_substeps_;
        double gravity_[3];
        bool use_gravity_;
        bool use_segregation_split = false; // 
        int nl_pressure_maxiter = 0;
        double nl_pressure_tolerance = 0.0;
        const double nl_tolerance = param.getDefault("nl_tolerance", 1e-9);
        const int nl_maxiter = param.getDefault("nl_maxiter", 30);

        // boost::scoped_ptr<Opm::GridManager> grid_;
        // boost::scoped_ptr<Opm::IncompPropertiesInterface> props_;
        // boost::scoped_ptr<Opm::WellsManager> wells_;
        // boost::scoped_ptr<Opm::RockCompressibility> rock_comp_;
        // Opm::SimulatorTimer simtimer_;
        // Opm::TwophaseState state_;
        // std::vector<double> src(num_cells, 0.0);
        // Opm::FlowBCManager bcs;
        // Opm::LinearSolverFactory linsolver(param);




        std::vector<double> totmob;
        std::vector<double> omega; // Will remain empty if no gravity.
        std::vector<double> rc; // Will remain empty if no rock compressibility.
        std::vector<double> porevol;
        double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);
        std::vector<double> transport_src;

        Opm::IncompTpfa psolver(*grid->c_grid(), props->permeability(), grav, linsolver);

        Opm::TransportModelTwophase tsolver(*grid->c_grid(), *props, nl_tolerance, nl_maxiter);

        typedef std::pair<std::vector<int>, std::vector<std::vector<int> > > ColMap;
        ColMap columns;

        std::vector<int> allcells(num_cells);
#endif
    };

    SimulatorTwophase::SimulatorTwophase()
    {
        pimpl_.reset(new Impl);
    }



    void SimulatorTwophase::init(const parameter::ParameterGroup& param)
    {
        pimpl_->init(param);
    }



    void SimulatorTwophase::run(const SimulatorTimer& timer, const Wells& wells, TwophaseState& state, WellState& well_state)
    {
        pimpl_->run(timer, wells, state, well_state);
    }



    void SimulatorTwophase::Impl::init(const parameter::ParameterGroup& /*param*/)
    {
    }



    void SimulatorTwophase::Impl::run(const SimulatorTimer& timer, const Wells& wells, TwophaseState& state, WellState& well_state)
    {
        (void) timer;
        (void) wells;
        (void) state;
        (void) well_state;
    }

} // namespace Opm
