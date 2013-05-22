#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/pvt/SinglePvtConstCompr.hpp>
#include <opm/core/props/pvt/SinglePvtDead.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
using namespace Opm;
using namespace std;



int main () {
    // read parameters from command-line
    const string filename = "../tests/SPE9small.DATA";
    cout << "Reading deck: " << filename << endl;
    const EclipseGridParser deck (filename);
    std::string mu_output = "mu_output";

    PhaseUsage phase_usage_;

    std::vector<std::tr1::shared_ptr<SinglePvtInterface> > props_;

    phase_usage_ = phaseUsageFromDeck(deck);
    enum PhaseIndex { Aqua = 0, Liquid = 1, Vapour = 2 };

    int samples = 0;

    std::fstream muos(mu_output.c_str(), std::fstream::out | std::fstream::trunc);
    if(!(muos.good())){
        std::cout << "Could not open"<< mu_output << std::endl;
        exit(3);
    }

    // Set the properties.
        props_.resize(phase_usage_.num_phases);
        // Water PVT
        if (phase_usage_.phase_used[Aqua]) {
            if (deck.hasField("PVTW")) {
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(deck.getPVTW().pvtw_));
            } else {
                // Eclipse 100 default.
                props_[phase_usage_.phase_pos[Aqua]].reset(new SinglePvtConstCompr(0.5*Opm::prefix::centi*Opm::unit::Poise));
            }
        }
        // Oil PVT
        if (phase_usage_.phase_used[Liquid]) {
            if (deck.hasField("PVDO")) {
                if (samples > 0) {
                    //props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDeadSpline(deck.getPVDO().pvdo_, samples));
                } else {
                    props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDead(deck.getPVDO().pvdo_));
                }
            } else if (deck.hasField("PVTO")) {
                //props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtLiveOil(deck.getPVTO().pvto_));
            } else if (deck.hasField("PVCDO")) {
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtConstCompr(deck.getPVCDO().pvcdo_));
            } else {
                THROW("Input is missing PVDO or PVTO\n");
            }
        }
        // Gas PVT
        if (phase_usage_.phase_used[Vapour]) {
            if (deck.hasField("PVDG")) {
                if (samples > 0) {
                    //props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDeadSpline(deck.getPVDG().pvdg_, samples));
                } else {
                    props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDead(deck.getPVDG().pvdg_));
                }
            } else if (deck.hasField("PVTG")) {
                //props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtLiveGas(deck.getPVTG().pvtg_));
            } else {
                THROW("Input is missing PVDG or PVTG\n");
            }
        }


        int n = 1;
        int np = 1; //phase_usage_.num_phases;
        double p[n];
        double r[n];
        double z[n];

        double mu[n];
        double dmudp[n];
        double dmudr[n];
        double mu_new[n];

        p[0] = 10000;

        // not in use yet
        r[0] = 0;
        z[0] = 0;


        for (int phase = 0; phase < np; ++phase) {
            props_[phase]->mu(n, p, r, mu_new,dmudp,dmudr);
            props_[phase]->mu(n, z, r, mu);

        }
        std::copy(mu,mu + np*n, std::ostream_iterator<double>(muos, " "));
        std::copy(mu_new,mu_new + np*n, std::ostream_iterator<double>(muos, " "));



}
