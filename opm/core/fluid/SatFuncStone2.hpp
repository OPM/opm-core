#ifndef SATFUNCSTONE2_HPP
#define SATFUNCSTONE2_HPP
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/fluid/blackoil/BlackoilPhases.hpp>
#include <vector>
namespace Opm
{
class SatFuncStone2: public BlackoilPhases
{
public:
    void init(const EclipseGridParser& deck, const int table_num, PhaseUsage phase_usg);
    void init(const EclipseGridParser& deck, const int table_num, PhaseUsage phase_usg,
              const int samples);
    void evalKr(const double* s, double* kr) const;
    void evalKrDeriv(const double* s, double* kr, double* dkrds) const;
    void evalPc(const double* s, double* pc) const;
    void evalPcDeriv(const double* s, double* pc, double* dpcds) const;
    double smin_[PhaseUsage::MaxNumPhases];
    double smax_[PhaseUsage::MaxNumPhases];
private:
    PhaseUsage phase_usage; // A copy of the outer class' phase_usage_.
    UniformTableLinear<double> krw_;
    UniformTableLinear<double> krow_;
    UniformTableLinear<double> pcow_;
    UniformTableLinear<double> krg_;
    UniformTableLinear<double> krog_;
    UniformTableLinear<double> pcog_;
    double krocw_; // = krow_(s_wc)
};
} // namespace Opm
#endif // SATFUNCSTONE2_HPP
