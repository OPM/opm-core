/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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


#include "config.h"

#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/EclipseGridParser.hpp>
#include <dune/porsol/blackoil/BlackoilFluid.hpp>


int main(int argc, char** argv)
{
    std::cout << "%{\n";

    // Parameters.
    Dune::parameter::ParameterGroup param(argc, argv);

    // Parser.
    std::string ecl_file = param.get<std::string>("filename");
    Dune::EclipseGridParser parser(ecl_file);
    // Look at the BlackoilFluid behaviour
    Opm::BlackoilFluid fluid;
    fluid.init(parser);
    Opm::BlackoilFluid::CompVec z0(0.0);
    z0[Opm::BlackoilFluid::Water] = param.getDefault("z_w", 0.0);
    z0[Opm::BlackoilFluid::Oil] = param.getDefault("z_o", 1.0);
    z0[Opm::BlackoilFluid::Gas] = param.getDefault("z_g", 0.0);
    int num_pts_p = param.getDefault("num_pts_p", 41);
    int num_pts_z = param.getDefault("num_pts_z", 51);
    double min_press = param.getDefault("min_press", 1e7);
    double max_press = param.getDefault("max_press", 3e7);
    int changing_component = param.getDefault("changing_component", int(Opm::BlackoilFluid::Gas));
    double min_z = param.getDefault("min_z", 0.0);
    double max_z = param.getDefault("max_z", 500.0);
    int variable = param.getDefault("variable", 0);
    Opm::BlackoilFluid::CompVec z = z0;
    std::cout << "%}\n"
              << "data = [\n";

    for (int i = 0; i < num_pts_p; ++i) {
        double pfactor = num_pts_p < 2 ? 0.0 : double(i)/double(num_pts_p - 1);
        double p = (1.0 - pfactor)*min_press + pfactor*max_press;
        for (int j = 0; j < num_pts_z; ++j) {
            double zfactor = num_pts_z < 2 ? 0.0 : double(j)/double(num_pts_z - 1);
            z[changing_component] = (1.0 - zfactor)*min_z + zfactor*max_z;
//             std::cout << p << ' ' << z << '\n';
            Opm::BlackoilFluid::FluidState state = fluid.computeState(Opm::BlackoilFluid::PhaseVec(p), z);
            std::cout.precision(6);
            std::cout.width(15);
            std::cout.fill(' ');
            double var = 0.0;
            switch (variable) {
            case 0:
                var = state.total_compressibility_;
                break;
            case 1:
                var = state.experimental_term_;
                break;
            case 2:
                var = state.saturation_[0];
                break;
            case 3:
                var = state.saturation_[1];
                break;
            case 4:
                var = state.saturation_[2];
                break;
            case 5:
                var = state.formation_volume_factor_[0];
                break;
            case 6:
                var = state.formation_volume_factor_[1];
                break;
            case 7:
                var = state.formation_volume_factor_[2];
                break;
            case 8:
                var = state.solution_factor_[0];
                break;
            case 9:
                var = state.solution_factor_[1];
                break;
            case 10:
                var = state.solution_factor_[2];
                break;
            default:
                THROW("Unknown varable specification: " << variable);
                break;
            }
            std::cout << var << ' ';
        }
        std::cout << '\n';
    }
    std::cout << "];\n\n"
              << "paxis = [\n";
    for (int i = 0; i < num_pts_p; ++i) {
        double pfactor = double(i)/double(num_pts_p - 1);
        double p = (1.0 - pfactor)*min_press + pfactor*max_press;
        std::cout << p << '\n';
    }
    std::cout << "];\n\n"
              << "zaxis = [\n";
    for (int j = 0; j < num_pts_z; ++j) {
        double zfactor = double(j)/double(num_pts_z - 1);
        std::cout << (1.0 - zfactor)*min_z + zfactor*max_z << '\n';
    }
    std::cout << "];\n";
}

