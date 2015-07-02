/*
  Copyright 2015 Andreas Lauser

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
/*!
 * \file
 * \copydoc Opm::SatFuncMultiplexer
 */
#ifndef OPM_SAT_FUNC_MULTIPLEXER_HPP
#define OPM_SAT_FUNC_MULTIPLEXER_HPP

#include "SatFuncGwseg.hpp"
#include "SatFuncSimple.hpp"
#include "SatFuncStone2.hpp"

#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <memory>

namespace Opm
{

/*!
 * \brief This is a multiplexer class which allows to saturation function to be used at
 *        runtime.
 */
class SatFuncMultiplexer
{
public:
    typedef NonuniformTableLinear<double> TableType;

    enum SatFuncType {
        Gwseg, // <- supposed to be the default
        Stone2,
        Simple
    };

// this is a helper macro which helps to save us from RSI in the following
#define OPM_MULTIPLEXER_CONST const
#define OPM_MULTIPLEXER_NON_CONST
#define OPM_MULTIPLEXER_SATFUNC_CALL(codeToCall, CONST)                 \
    switch (satFuncType_) {                                             \
    case Gwseg: {                                                       \
        typedef SatFuncGwseg<TableType> SatFunc;                        \
        __attribute__((unused)) CONST SatFunc& satFunc =                \
            *static_cast<SatFunc*>(satFunc_.get());                     \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case Stone2: {                                                      \
        typedef SatFuncStone2<TableType> SatFunc;                       \
        __attribute__((unused)) CONST SatFunc& satFunc =                \
            *static_cast<SatFunc*>(satFunc_.get());                     \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case Simple: {                                                      \
        typedef SatFuncSimple<TableType> SatFunc;                       \
        __attribute__((unused)) CONST SatFunc& satFunc =                \
            *static_cast<SatFunc*>(satFunc_.get());                     \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    };

    SatFuncMultiplexer()
    {}

    ~SatFuncMultiplexer()
    {}

    // this constructor not really a copy constructor, but it is
    // required to make std::unique_ptr happy
    SatFuncMultiplexer(const SatFuncMultiplexer&)
    {}

    // this operator does not do anything and is thus not a copy operator, but it is
    // required to make std::unique_ptr happy on old compilers
    SatFuncMultiplexer& operator=(const SatFuncMultiplexer& other)
    { return *this; }

    /*!
     * \brief Pick the correct saturation function type and initialize the object using
     *        an ECL deck.
     */
    void initFromDeck(Opm::EclipseStateConstPtr eclState,
                      size_t tableIdx,
                      SatFuncType satFuncType)
    {
        auto phaseUsage = phaseUsageFromDeck(eclState);

        setSatFuncType(satFuncType);

        OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.init(eclState,
                                                  tableIdx,
                                                  phaseUsage,
                                                  /*samples=*/0), OPM_MULTIPLEXER_NON_CONST);
    }

    /*!
     * \brief Returns the type of the actually selected saturation function.
     */
    void setSatFuncType(SatFuncType newSatFuncType)
    {
        // set the type of the saturation function
        satFuncType_ = newSatFuncType;

        // call the default constructor for the selected saturation function
        switch (satFuncType_) {
        case Gwseg: {
            typedef SatFuncGwseg<TableType> SatFunc;
            satFunc_.reset(new SatFunc);
            break;
        }
        case Stone2: {
            typedef SatFuncStone2<TableType> SatFunc;
            satFunc_.reset(new SatFunc);
            break;
        }
        case Simple: {
            typedef SatFuncSimple<TableType> SatFunc;
            satFunc_.reset(new SatFunc);
            break;
        }
        }
    }

    /*!
     * \brief Returns which saturation function is actually selected.
     */
    SatFuncType satFuncType() const
    { return satFuncType_; }

    /*!
     * \brief Returns a pointer to the base saturation function.
     *
     * All underlying saturation functions derive from this class. It stores the stuff
     * which is common for all saturation functions (e.g. the endpoint scaling
     * parameters).
     */
    const SatFuncBase<TableType>& getSatFuncBase() const
    { return *satFunc_; }

    SatFuncBase<TableType>& getSatFuncBase()
    { return *satFunc_; }

    /*!
     * \brief Return the raw "Gwseg" saturation function object.
     *
     * If the saturation function type is not "Gwseg", an assertation is triggered for
     * debug builds.
     */
    SatFuncGwseg<TableType>& getGwseg()
    {
        assert(satFuncType_ == Gwseg);
        return *static_cast<SatFuncGwseg<TableType>*>(satFunc_.get());
    }

    const SatFuncGwseg<TableType>& getGwseg() const
    {
        assert(satFuncType_ == Gwseg);
        return *static_cast<const SatFuncGwseg<TableType>*>(satFunc_.get());
    }

    /*!
     * \brief Return the raw "Simple" saturation function object.
     *
     * If the saturation function type is not "Simple", an assertation is triggered for
     * debug builds.
     */
    SatFuncSimple<TableType>& getSimple()
    {
        assert(satFuncType_ == Simple);
        return *static_cast<SatFuncSimple<TableType>*>(satFunc_.get());
    }

    const SatFuncSimple<TableType>& getSimple() const
    {
        assert(satFuncType_ == Simple);
        return *static_cast<const SatFuncSimple<TableType>*>(satFunc_.get());
    }

    /*!
     * \brief Return the raw "Stone2" saturation function object.
     *
     * If the saturation function type is not "Stone2", an assertation is triggered for
     * debug builds.
     */
    SatFuncStone2<TableType>& getStone2()
    {
        assert(satFuncType_ == Stone2);
        return *static_cast<SatFuncStone2<TableType>*>(satFunc_.get());
    }

    const SatFuncStone2<TableType>& getStone2() const
    {
        assert(satFuncType_ == Stone2);
        return *static_cast<const SatFuncStone2<TableType>*>(satFunc_.get());
    }

    template <class FluidState>
    void evalKr(const FluidState& fluidState, double* kr) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKr(fluidState, kr), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalKrDeriv(const FluidState& fluidState, double* kr, double* dkrds) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKrDeriv(fluidState, kr, dkrds), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalPc(const FluidState& fluidState, double* pc) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPc(fluidState, pc), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalPcDeriv(const FluidState& fluidState, double* pc, double* dpcds) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPcDeriv(fluidState, pc, dpcds), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalKr(const FluidState& fluidState, double* kr, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKr(fluidState, kr, epst), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalKr(const FluidState& fluidState, double* kr, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKr(fluidState, kr, epst, epst_hyst, sat_hyst), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalKrDeriv(const FluidState& fluidState, double* kr, double* dkrds, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKrDeriv(fluidState, kr, dkrds, epst), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalKrDeriv(const FluidState& fluidState, double* kr, double* dkrds, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKrDeriv(fluidState, kr, dkrds, epst, epst_hyst, sat_hyst), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalPc(const FluidState& fluidState, double* pc, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPc(fluidState, pc, epst), OPM_MULTIPLEXER_CONST); }

    template <class FluidState>
    void evalPcDeriv(const FluidState& fluidState, double* pc, double* dpcds, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPcDeriv(fluidState, pc, dpcds, epst), OPM_MULTIPLEXER_CONST); }

// undefine the helper macro here because we don't need it anymore
#undef OPM_MULTIPLEXER_SATFUNC_CALL
#undef OPM_MULTIPLEXER_NON_CONST
#undef OPM_MULTIPLEXER_CONST

private:
    SatFuncType satFuncType_;
    std::unique_ptr<SatFuncBase<TableType> > satFunc_;
};
} // namespace Opm

#endif // OPM_SAT_FUNC_MULTIPLEXER_HPP
