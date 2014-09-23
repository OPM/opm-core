#include <config.h>

#include <opm/core/props/pvt/PvtConstCompr.hpp>
#include <opm/core/props/pvt/PvtDead.hpp>
#include <opm/core/props/pvt/PvtDeadSpline.hpp>
#include <opm/core/props/pvt/PvtLiveOil.hpp>
#include <opm/core/props/pvt/PvtLiveGas.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/pvt/BlackoilPvtProperties.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

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


std::vector<std::shared_ptr<PvtInterface> > getProps(Opm::DeckConstPtr deck, Opm::EclipseStateConstPtr eclipseState, PhaseUsage phase_usage_){
    enum PhaseIndex { Aqua = 0, Liquid = 1, Vapour = 2 };
    int samples = 0;

    std::vector<std::shared_ptr<PvtInterface> > props_;
    // Set the properties.
    props_.resize(phase_usage_.num_phases);

    // Water PVT
    if (phase_usage_.phase_used[Aqua]) {
        if (deck->hasKeyword("PVTW")) {
            std::shared_ptr<PvtConstCompr> pvtw(new PvtConstCompr);
            pvtw->initFromWater(deck->getKeyword("PVTW"));

            props_[phase_usage_.phase_pos[Aqua]] = pvtw;
        } else {
            // Eclipse 100 default.
            props_[phase_usage_.phase_pos[Aqua]].reset(new PvtConstCompr(0.5*Opm::prefix::centi*Opm::unit::Poise));
        }
    }

    // Oil PVT
    if (phase_usage_.phase_used[Liquid]) {
        const auto& pvdoTables = eclipseState->getPvdoTables();
        const auto& pvtoTables = eclipseState->getPvtoTables();
        if (!pvdoTables.empty()) {
            if (samples > 0) {
                std::shared_ptr<PvtDeadSpline> splinePvt(new PvtDeadSpline);
                splinePvt->initFromOil(pvdoTables, samples);
                props_[phase_usage_.phase_pos[Liquid]] = splinePvt;
            } else {
                std::shared_ptr<PvtDead> deadPvt(new PvtDead);
                deadPvt->initFromOil(pvdoTables);
                props_[phase_usage_.phase_pos[Liquid]] = deadPvt;
            }
        } else if (!pvtoTables.empty()) {
            props_[phase_usage_.phase_pos[Liquid]].reset(new PvtLiveOil(pvtoTables));
        } else if (deck->hasKeyword("PVCDO")) {
            std::shared_ptr<PvtConstCompr> pvcdo(new PvtConstCompr);
            pvcdo->initFromOil(deck->getKeyword("PVCDO"));

            props_[phase_usage_.phase_pos[Liquid]] = pvcdo;
        } else {
            OPM_THROW(std::runtime_error, "Input is missing PVDO, PVCDO or PVTO\n");
        }
    }
    // Gas PVT
    if (phase_usage_.phase_used[Vapour]) {
        const auto& pvdgTables = eclipseState->getPvdgTables();
        const auto& pvtgTables = eclipseState->getPvtgTables();
        if (!pvdgTables.empty()) {
            if (samples > 0) {
                std::shared_ptr<PvtDeadSpline> splinePvt(new PvtDeadSpline);
                splinePvt->initFromGas(pvdgTables, samples);
                props_[phase_usage_.phase_pos[Vapour]] = splinePvt;
            } else {
                std::shared_ptr<PvtDead> deadPvt(new PvtDead);
                deadPvt->initFromGas(pvdgTables);
                props_[phase_usage_.phase_pos[Vapour]] = deadPvt;
            }
        } else if (!pvtgTables.empty()) {
            props_[phase_usage_.phase_pos[Vapour]].reset(new PvtLiveGas(pvtgTables));
        } else {
            OPM_THROW(std::runtime_error, "Input is missing PVDG or PVTG\n");
        }
    }

    return props_;
}

void testmu(const double reltol, int n, int np, const std::vector<int> &pvtTableIdx, std::vector<double> p, std::vector<double> r,std::vector<double> z,
            std::vector<std::shared_ptr<PvtInterface> > props_, std::vector<Opm::PhasePresence> condition)
{
    std::vector<double> mu(n);
    std::vector<double> dmudp(n);
    std::vector<double> dmudr(n);
    std::vector<double> mu_new(n);
    double dmudp_diff;
    double dmudr_diff;
    double dmudp_diff_u;
    double dmudr_diff_u;

    // test mu
    for (int phase = 0; phase < np; ++phase) {
        props_[phase]->mu(n, &pvtTableIdx[0], &p[0], &r[0], &condition[0], &mu_new[0], &dmudp[0], &dmudr[0]);
        props_[phase]->mu(n, &pvtTableIdx[0], &p[0], &z[0], &mu[0]);
        dmudp_diff = (mu_new[1]-mu_new[0])/(p[1]-p[0]);
        dmudr_diff = (mu_new[2]-mu_new[0])/(r[2]-r[0]);
        dmudp_diff_u = (mu_new[4]-mu_new[3])/(p[4]-p[3]);
        dmudr_diff_u = (mu_new[5]-mu_new[3])/(r[5]-r[3]);

        for (int i = 0; i < n; ++i){
            BOOST_CHECK_CLOSE(mu_new[i],mu[i],reltol);
        }

        // saturated case
        BOOST_CHECK_CLOSE(dmudp_diff,dmudp[0],reltol);
        BOOST_CHECK_CLOSE(dmudr_diff,dmudr[0],reltol);

        // unsaturated case
        BOOST_CHECK_CLOSE(dmudp_diff_u,dmudp[3],reltol);
        BOOST_CHECK_CLOSE(dmudr_diff_u,dmudr[3],reltol);

    }
}

void testb(const double reltol, int n, int np, const std::vector<int> &pvtTableIdx, std::vector<double> p, std::vector<double> r,std::vector<double> z,
            std::vector<std::shared_ptr<PvtInterface> > props_, std::vector<Opm::PhasePresence> condition)
{
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

    for (int phase = 0; phase < np; ++phase) {
        props_[phase]->b(n, &pvtTableIdx[0], &p[0], &r[0], &condition[0], &b[0], &dbdp[0], &dbdr[0]);
        props_[phase]->dBdp(n, &pvtTableIdx[0], &p[0], &z[0], &B[0], &dBdp[0]);
        dbdp_diff = (b[1]-b[0])/(p[1]-p[0]);
        dbdr_diff = (b[2]-b[0])/(r[2]-r[0]);
        dbdp_diff_u = (b[4]-b[3])/(p[4]-p[3]);
        dbdr_diff_u = (b[5]-b[3])/(r[5]-r[3]);
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
}

void testrsSat(double reltol, int n, int np, const std::vector<int> &pvtTableIdx, std::vector<double> p, std::vector<std::shared_ptr<PvtInterface> > props_){
    // test bublepoint pressure
    std::vector<double> rs(n);
    std::vector<double> drsdp(n);
    double drsdp_diff;
    double drsdp_diff_u;

    for (int phase = 0; phase < np; ++phase) {
        props_[phase] ->rsSat(n, &pvtTableIdx[0], &p[0], &rs[0], &drsdp[0]);

        drsdp_diff = (rs[1]-rs[0])/(p[1]-p[0]);
        drsdp_diff_u = (rs[4]-rs[3])/(p[4]-p[3]);

        // saturated case
        BOOST_CHECK_CLOSE(drsdp_diff,drsdp[0], reltol);

        // unsaturad case
        BOOST_CHECK_CLOSE(drsdp_diff_u,drsdp[3], reltol);

    }
}

void testrvSat(double reltol, int n, int np, const std::vector<int> &pvtTableIdx, std::vector<double> p, std::vector<std::shared_ptr<PvtInterface> > props_){
    // test rv saturated
    std::vector<double> rv(n);
    std::vector<double> drvdp(n);
    double drvdp_diff;
    double drvdp_diff_u;

    for (int phase = 0; phase < np; ++phase) {
        props_[phase] ->rvSat(n, &pvtTableIdx[0], &p[0], &rv[0], &drvdp[0]);

        drvdp_diff = (rv[1]-rv[0])/(p[1]-p[0]);
        drvdp_diff_u = (rv[4]-rv[3])/(p[4]-p[3]);

        // saturated case
        BOOST_CHECK_CLOSE(drvdp_diff,drvdp[0], reltol);

        // unsaturad case
        BOOST_CHECK_CLOSE(drvdp_diff_u,drvdp[3], reltol);

    }
}

BOOST_AUTO_TEST_CASE(test_liveoil)
{


    // read eclipse deck
    const std::string filename = "liveoil.DATA";
    cout << "Reading deck: " << filename << endl;
    Opm::ParserPtr parser(new Opm::Parser());
    Opm::DeckConstPtr deck(parser->parseFile(filename));
    Opm::EclipseStateConstPtr eclipseState(new EclipseState(deck));

    // setup pvt interface
    PhaseUsage phase_usage_ = phaseUsageFromDeck(deck);
    std::vector<std::shared_ptr<PvtInterface> > props_ = getProps(deck, eclipseState, phase_usage_);


    // setup a test case. We will check 6 [p,r] pairs and compare them to both the [p,z] interface and a finite difference
    // approximation of the derivatives.
    const int n = 6;
    const int np = phase_usage_.num_phases;
    std::vector<int> pvtRegionIdx(n, 0);

    // the relative tolerance in percentage for acceptable difference in values
    const double reltolper = 1e-9;
    // the relative tolerance in percentage for acceptable difference in values for viscosity
    const double reltolpermu = 1e-1;

    std::vector<double> p(n);
    std::vector<double> r(n);
    std::vector<PhasePresence> condition(n);
    std::vector<double> z(n * np);


    // Used for forward difference calculations
    const double h_p = 1e4;
    const double h_r = 1;

    // saturated
    p[0] = 1e7;
    p[1] = p[0] + h_p;
    p[2] = p[0];

    r[0] = 200;
    r[1] = 200;
    r[2] = 200 + h_r;

    condition[0].setFreeGas();
    condition[1].setFreeGas();
    condition[2].setFreeGas();


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

    testmu(reltolpermu, n, np, pvtRegionIdx, p, r,z, props_, condition);

    testb(reltolper,n,np,pvtRegionIdx,p,r,z,props_,condition);

    testrsSat(reltolper,n,np,pvtRegionIdx,p,props_);

    testrvSat(reltolper,n,np,pvtRegionIdx,p,props_);


}

BOOST_AUTO_TEST_CASE(test_wetgas)
{


    // read eclipse deck

    const std::string filename = "wetgas.DATA";
    cout << "Reading deck: " << filename << endl;
    Opm::ParserPtr parser(new Opm::Parser());
    Opm::DeckConstPtr deck(parser->parseFile(filename));
    Opm::EclipseStateConstPtr eclipseState(new EclipseState(deck));

    // setup pvt interface
    PhaseUsage phase_usage_ = phaseUsageFromDeck(deck);
    std::vector<std::shared_ptr<PvtInterface> > props_ = getProps(deck,eclipseState,phase_usage_);


    // setup a test case. We will check 6 [p,r] pairs and compare them to both the [p,z] interface and a finite difference
    // approximation of the derivatives.
    const int n = 6;
    const int np = phase_usage_.num_phases;
    std::vector<int> pvtRegionIdx(n, 0);

    // the relative tolerance in percentage for acceptable difference in values
    const double reltolper = 1e-9;
    // the relative tolerance in percentage for acceptable difference in values for viscosity
    const double reltolpermu = 1e-1;

    std::vector<double> p(n);
    std::vector<double> r(n);
    std::vector<PhasePresence> condition(n);
    std::vector<double> z(n * np);


    // Used for forward difference calculations
    const double h_p = 1e4;
    const double h_r = 1e-7;

    // saturated
    p[0] = 1e7;
    p[1] = p[0] + h_p;
    p[2] = p[0];

    r[0] = 5e-5;
    r[1] = 5e-5;
    r[2] = 5e-5 + h_r;

    condition[0].setFreeOil();
    condition[1].setFreeOil();
    condition[2].setFreeOil();


    // undersaturated
    p[3] = p[0];
    p[4] = p[1];
    p[5] = p[2];

    r[3] = 1e-5;
    r[4] = 1e-5;
    r[5] = 1e-5 +h_r;

    // Corresponing z factors, used to compare with the [p,z] interface
    for (int i = 0; i < n; ++i) {
        z[0+i*np] = 0; z[1+i*np] = r[i];
        z[2+i*np] = 1;

    }

    testmu(reltolpermu, n, np, pvtRegionIdx, p, r,z, props_, condition);

    testb(reltolper,n,np,pvtRegionIdx,p,r,z,props_,condition);

    testrsSat(reltolper,n,np,pvtRegionIdx,p,props_);

    testrvSat(reltolper,n,np,pvtRegionIdx,p,props_);

}
