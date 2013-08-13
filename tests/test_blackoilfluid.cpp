#include <config.h>

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/pvt/SinglePvtConstCompr.hpp>
#include <opm/core/props/pvt/SinglePvtDead.hpp>
#include <opm/core/props/pvt/SinglePvtDeadSpline.hpp>
#include <opm/core/props/pvt/SinglePvtLiveOil.hpp>
#include <opm/core/props/pvt/SinglePvtLiveGas.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing
#define BOOST_TEST_MODULE BlackoilFluidTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <memory>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>

using namespace Opm;
using namespace std;

BOOST_AUTO_TEST_CASE(test_blackoilfluid)
{


    // read eclipse deck
    const string filename = "testFluid.DATA";
    cout << "Reading deck: " << filename << endl;
    const EclipseGridParser deck (filename);

    // setup pvt interface
    std::vector<std::shared_ptr<SinglePvtInterface> > props_;
    PhaseUsage phase_usage_ = phaseUsageFromDeck(deck);
    enum PhaseIndex { Aqua = 0, Liquid = 1, Vapour = 2 };
    int samples = 0;


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
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDeadSpline(deck.getPVDO().pvdo_, samples));
            } else {
                props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtDead(deck.getPVDO().pvdo_));
            }
        } else if (deck.hasField("PVTO")) {

            props_[phase_usage_.phase_pos[Liquid]].reset(new SinglePvtLiveOil(deck.getPVTO().pvto_));
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
                props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDeadSpline(deck.getPVDG().pvdg_, samples));
            } else {
                props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtDead(deck.getPVDG().pvdg_));
            }
        } else if (deck.hasField("PVTG")) {
            props_[phase_usage_.phase_pos[Vapour]].reset(new SinglePvtLiveGas(deck.getPVTG().pvtg_));
        } else {
            THROW("Input is missing PVDG or PVTG\n");
        }
    }

    // setup a test case. We will check 6 [p,r] pairs and compare them to both the [p,z] interface and a finite difference
    // approximation of the derivatives.
    const int n = 6;
    const int np = phase_usage_.num_phases;

    // the tolerance for acceptable difference in values
    const double reltol = 1e-9;

    std::vector<double> p(n);
    std::vector<double> r(n);
    std::vector<double> z(n * np);

    std::vector<double> mu(n);
    std::vector<double> dmudp(n);
    std::vector<double> dmudr(n);
    std::vector<double> mu_new(n);
    double dmudp_diff;
    double dmudr_diff;
    double dmudp_diff_u;
    double dmudr_diff_u;


    // Used for forward difference calculations
    const double h_p = 100000;
    const double h_r = 1;

    // saturated
    p[0] = 10000000;
    p[1] = p[0] + h_p;
    p[2] = p[0];

    r[0] = 200;
    r[1] = 200;
    r[2] = 200 + h_r;

    // undersaturated
    p[3] = p[0];
    p[4] = p[1];
    p[5] = p[2];

    r[3] = 50;
    r[4] = 50;
    r[5] = 50 +h_r;


    // Corresponing z factors, used to compare with the [p,z] interface
    for (int i = 0; i < n; ++i) {
        z[0+i*np] = 0; z[1+i*np] = 1;
        z[2+i*np] = r[i];

    }

    // test mu
    for (int phase = 1; phase < 2; ++phase) {
        props_[phase]->mu(n, &p[0], &r[0], &mu_new[0], &dmudp[0], &dmudr[0]);
        props_[phase]->mu(n, &p[0], &z[0], &mu[0]);
        dmudp_diff = (mu_new[1]-mu_new[0])/h_p;
        dmudr_diff = (mu_new[2]-mu_new[0])/h_r;
        dmudp_diff_u = (mu_new[4]-mu_new[3])/h_p;
        dmudr_diff_u = (mu_new[5]-mu_new[3])/h_r;

        for (int i = 0; i < n; ++i){
            BOOST_CHECK_CLOSE(mu_new[i],mu[i],reltol);
        }
        // saturated case
        BOOST_CHECK_CLOSE(dmudp_diff,dmudp[0],reltol);
        BOOST_CHECK_CLOSE(dmudr_diff,dmudr[0],reltol);

        // unsaturated case
        BOOST_CHECK_CLOSE(dmudp_diff_u,dmudp[3] , reltol);
        BOOST_CHECK_CLOSE(dmudr_diff_u,dmudr[3] , reltol);

    }

    // test b
    std::vector<double> b(n);
    std::vector<double> B(n);
    std::vector<double> invB(n);
    std::vector<double> dinvBdp(n);
    std::vector<double> dBdp(n);
    std::vector<double> dbdr(n);
    std::vector<double> dbdp(n);
    double dbdp_diff;
    double dbdr_diff;
    double dbdp_diff_u;
    double dbdr_diff_u;

    for (int phase = 1; phase < 2; ++phase) {
        props_[phase]->b(n, &p[0], &r[0], &b[0], &dbdp[0], &dbdr[0]);
        //props_[phase]->B(n, p, z, B);
        props_[phase]->dBdp(n, &p[0], &z[0], &B[0], &dBdp[0]);
        dbdp_diff = (b[1]-b[0])/h_p;
        dbdr_diff = (b[2]-b[0])/h_r;
        dbdp_diff_u = (b[4]-b[3])/h_p;
        dbdr_diff_u = (b[5]-b[3])/h_r;
        for (int i = 0; i < n; ++i){
            invB[i] = 1/B[i];
            dinvBdp[i] = -1/pow(B[i],2) * dBdp[i];
        }

        for (int i = 0; i < n; ++i){
            BOOST_CHECK_CLOSE(invB[i],b[i] , reltol);
            BOOST_CHECK_CLOSE(dinvBdp[i],dbdp[i] , reltol);

        }
        // saturated case
        BOOST_CHECK_CLOSE(dbdp_diff,dbdp[0], reltol);
        BOOST_CHECK_CLOSE(dbdr_diff,dbdr[0], reltol);

        // unsaturated case
        BOOST_CHECK_CLOSE(dbdp_diff_u,dbdp[3], reltol);
        BOOST_CHECK_CLOSE(dbdr_diff_u,dbdr[3], reltol);
    }

    // test bublepoint pressure
    std::vector<double> rbub(n);
    std::vector<double> drbubdp(n);
    double drbubdp_diff;
    double drbubdp_diff_u;

    for (int phase = 1; phase < 2; ++phase) {
        props_[phase] ->rbub(n, &p[0], &rbub[0], &drbubdp[0]);

        drbubdp_diff = (rbub[1]-rbub[0])/h_p;
        drbubdp_diff_u = (rbub[4]-rbub[3])/h_p;

        // saturated case
        BOOST_CHECK_CLOSE(drbubdp_diff,drbubdp[0], reltol);

        // unsaturad case
        BOOST_CHECK_CLOSE(drbubdp_diff_u,drbubdp[3], reltol);

    }





}
