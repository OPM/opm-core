/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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


#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/EclipseGridParser.hpp>
#include <dune/porsol/blackoil/BlackoilFluid.hpp>


int main(int argc, char** argv)
{
    // Parameters.
    Dune::parameter::ParameterGroup param(argc, argv);

    // Parser.
    std::string ecl_file = param.get<std::string>("filename");
    Dune::EclipseGridParser parser(ecl_file);

    // Look at the BlackoilFluid behaviour
    Opm::BlackoilFluid fluid;
    fluid.init(parser);
    Opm::BlackoilFluid::CompVec z(0.0);
    z[Opm::BlackoilFluid::Water] = param.getDefault("z_w", 0.0);
    z[Opm::BlackoilFluid::Oil] = param.getDefault("z_o", 1.0);
    z[Opm::BlackoilFluid::Gas] = param.getDefault("z_g", 0.0);
    int num_pts = param.getDefault("num_pts", 41);
    double min_press = param.getDefault("min_press", 1e7);
    double max_press = param.getDefault("max_press", 3e7);
    bool print_compr = param.getDefault("print_compr", true);
    for (int i = 0; i < num_pts; ++i) {
        double factor = double(i)/double(num_pts - 1);
        double p = (1.0 - factor)*min_press + factor*max_press;
        Opm::BlackoilFluid::FluidState state = fluid.computeState(Opm::BlackoilFluid::PhaseVec(p), z);
        double totmob = 0.0;
        for (int phase = 0; phase < Opm::BlackoilFluid::numPhases; ++phase) {
            totmob += state.mobility_[phase];
        }
        if (print_compr) {
            std::cout.precision(6);
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << p << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << state.total_compressibility_ << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << totmob << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << state.total_phase_volume_density_ << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << state.experimental_term_ << '\n';
        } else {
            std::cout.precision(6);
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << p << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << state.saturation_[0] << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << state.saturation_[1] << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << state.saturation_[2] << "  ";
            std::cout.width(15);
            std::cout.fill(' ');
            std::cout << state.total_phase_volume_density_ << '\n';
        }
    }
}
