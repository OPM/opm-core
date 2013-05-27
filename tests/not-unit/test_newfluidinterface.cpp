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

#include <iostream>
#include <iterator>
#include <vector>
#include <string>
using namespace Opm;
using namespace std;


// The function object divides a Factor with an element
template <class Type>
class MultValue
{
   private:
      Type Factor;   // The value to multiply by
   public:
      // Constructor initializes the value to multiply by
      MultValue ( const Type& _Val ) : Factor ( _Val ) {
      }

      // The function call for the element to be multiplied
      int operator ( ) ( Type& elem ) const
      {
         return  Factor / elem;
      }
};
int main () {
    // read parameters from command-line
    const string filename = "../../opm-core/tests/not-unit/blackoil/SPE9small.DATA";
    cout << "Reading deck: " << filename << endl;
    const EclipseGridParser deck (filename);
    std::string mu_output = "mu_output";
    std::string b_output = "b_output";
    std::string rbub_output = "rbub_output";

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

    std::fstream bos(b_output.c_str(), std::fstream::out | std::fstream::trunc);
    bos << setiosflags(ios::scientific) << setprecision(12);
    if(!(bos.good())){
        std::cout << "Could not open"<< b_output << std::endl;
        exit(3);
    }

    std::fstream rbubos(rbub_output.c_str(), std::fstream::out | std::fstream::trunc);
    rbubos << setiosflags(ios::scientific) << setprecision(12);
    if(!(rbubos.good())){
        std::cout << "Could not open"<< rbub_output << std::endl;
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


        int n = 6;
        int np = 3; //phase_usage_.num_phases;
        double p[n];
        double r[n];
        double z[np*n];

        double mu[n];
        double dmudp[n];
        double dmudr[n];
        double mu_new[n];
        double dmudp_diff;
        double dmudr_diff;
        double dmudp_diff_u;
        double dmudr_diff_u;

        //double rf[n];

        double h = 1;
        double rh = 1;

        p[0] = 10000000;
        p[1] = p[0] + h;
        p[2] = 10000000;

        p[3] = 10000000;
        p[4] = p[0] + h;
        p[5] = 10000000;


        // saturated
        r[0] = 200;
        r[1] = 200;
        r[2] = 200 + rh;

        // undersaturated
        r[3] = 50;
        r[4] = 50;
        r[5] = 50 +rh;


        for (int i = 0; i < n; ++i) {
                z[0+i*np] = 0; z[1+i*np] = 1;
                z[2+i*np] = r[i];

       }



        // test mu
        for (int phase = 1; phase < 2; ++phase) {
            props_[phase]->mu(n, p, r, mu_new,dmudp,dmudr);
            props_[phase]->mu(n, p, z, mu);
            dmudp_diff = (mu_new[1]-mu_new[0])/h;
            dmudr_diff = (mu_new[2]-mu_new[0])/rh;
            dmudp_diff_u = (mu_new[4]-mu_new[3])/h;
            dmudr_diff_u = (mu_new[5]-mu_new[3])/rh;

            std::copy(mu,mu + n, std::ostream_iterator<double>(muos, " "));
            muos << "\n";
            std::copy(mu_new,mu_new + n, std::ostream_iterator<double>(muos, " "));
            muos << "\n";
            std::copy(dmudp,dmudp + n, std::ostream_iterator<double>(muos, " "));
            muos << "\n";
            muos << dmudp_diff << " " << dmudp_diff_u << "\n";
            std::copy(dmudr,dmudr + n, std::ostream_iterator<double>(muos, " "));
            muos << "\n";
            muos << dmudr_diff << " " << dmudr_diff_u << "\n";
        }

        // test b
        double b[n];
        double B[n];
        double invB[n];
        double dinvBdp[n];
        double dBdp[n];
        double dbdr[n];
        double dbdp[n];
        double dbdp_diff;
        double dbdr_diff;
        double dbdp_diff_u;
        double dbdr_diff_u;

        for (int phase = 1; phase < 2; ++phase) {
            props_[phase]->b(n, p, r, b,dbdp,dbdr);
            //props_[phase]->B(n, p, z, B);
            props_[phase]->dBdp(n, p, z, B, dBdp);
            dbdp_diff = (b[1]-b[0])/h;
            dbdr_diff = (b[2]-b[0])/rh;
            dbdp_diff_u = (b[4]-b[3])/h;
            dbdr_diff_u = (b[5]-b[3])/rh;
            for (int i = 0; i < n; ++i){
                invB[i] = 1/B[i];
                dinvBdp[i] = -1/pow(B[i],2) * dBdp[i];

            }
            std::copy(b,b + n, std::ostream_iterator<double>(bos, " "));
            bos << "\n";
            std::copy(invB,invB + n, std::ostream_iterator<double>(bos, " "));
            bos << "\n";
            std::copy(dinvBdp,dinvBdp + n, std::ostream_iterator<double>(bos, " "));
            bos << "\n";
            std::copy(dbdp,dbdp + n, std::ostream_iterator<double>(bos, " "));
            bos << "\n";
            bos << dbdp_diff << " " << dbdp_diff_u << "\n";
            std::copy(dbdr,dbdr + n, std::ostream_iterator<double>(bos, " "));
            bos << "\n";
            bos << dbdr_diff << " " << dbdr_diff_u << "\n";

        }

        // test rbub

        double rbub[n];
        double drbubdp[n];
        double drbubdp_diff;

        for (int phase = 1; phase < 2; ++phase) {
            props_[phase] ->rbub(n,p,rbub,drbubdp);

            drbubdp_diff = (rbub[1]-rbub[0])/h;
            std::copy(rbub,rbub + n, std::ostream_iterator<double>(rbubos, " "));
            rbubos << "\n";
            std::copy(drbubdp,drbubdp + n, std::ostream_iterator<double>(rbubos, " "));
            rbubos << drbubdp_diff;
            rbubos << "\n";


        }




}

