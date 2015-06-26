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
#define OPM_MULTIPLEXER_SATFUNC_CALL(codeToCall)    \
    switch (satFuncType_) {                         \
    case Gwseg: {                                   \
        typedef SatFuncGwseg<TableType> SatFunc;    \
        if (false) SatFunc dummyForWarning;         \
        auto& satFunc = getGwseg();                 \
        codeToCall;                                 \
        break;                                      \
    }                                               \
    case Stone2: {                                  \
        typedef SatFuncStone2<TableType> SatFunc;   \
        if (false) SatFunc dummyForWarning;         \
        auto& satFunc = getStone2();                \
        codeToCall;                                 \
        break;                                      \
    }                                               \
    case Simple: {                                  \
        typedef SatFuncSimple<TableType> SatFunc;   \
        if (false) SatFunc dummyForWarning;         \
        auto& satFunc = getSimple();                \
        codeToCall;                                 \
        break;                                      \
    }                                               \
    };

    SatFuncMultiplexer()
        : wasInitialized_(false)
    {}

    ~SatFuncMultiplexer()
    {
        if (!wasInitialized_)
            return;

        // call the destructor of the selected saturation function
        OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.~SatFunc());
    }

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
                                                  /*samples=*/0));
    }

    /*!
     * \brief Returns the type of the actually selected saturation function.
     */
    void setSatFuncType(SatFuncType newSatFuncType)
    {
        assert(!wasInitialized_);

        // set the type of the saturation function
        satFuncType_ = newSatFuncType;

        // call the default constructor for the selected saturation function
        OPM_MULTIPLEXER_SATFUNC_CALL(new(&satFunc) SatFunc());
        wasInitialized_ = true;
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
    {
        switch (satFuncType()) {
        case Gwseg:
            return getGwseg();
        case Stone2:
            return getStone2();
        case Simple:
            return getSimple();
        }
        OPM_THROW(std::logic_error, "Unhandled saturation function type");
    }

    SatFuncBase<TableType>& getSatFuncBase()
    {
        return const_cast<SatFuncBase<TableType>&>(static_cast<const SatFuncMultiplexer*>(this)->getSatFuncBase());
    }

// the strict aliasing compiler warnings need to be ignored here. This is safe because we
// deal with the alignment of the data ourselfs.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
    /*!
     * \brief Return the raw "Gwseg" saturation function object.
     *
     * If the saturation function type is not "Gwseg", an assertation is triggered for
     * debug builds.
     */
    SatFuncGwseg<TableType>& getGwseg()
    {
        assert(satFuncType_ == Gwseg);
        return *reinterpret_cast<SatFuncGwseg<TableType>*>(satFuncsData_.gwseg);
    }

    const SatFuncGwseg<TableType>& getGwseg() const
    {
        assert(satFuncType_ == Gwseg);
        return *reinterpret_cast<const SatFuncGwseg<TableType>*>(satFuncsData_.gwseg);
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
        return *reinterpret_cast<SatFuncSimple<TableType>*>(satFuncsData_.simple);
    }

    const SatFuncSimple<TableType>& getSimple() const
    {
        assert(satFuncType_ == Simple);
        return *reinterpret_cast<const SatFuncSimple<TableType>*>(satFuncsData_.simple);
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
        return *reinterpret_cast<SatFuncStone2<TableType>*>(satFuncsData_.stone2);
    }

    const SatFuncStone2<TableType>& getStone2() const
    {
        assert(satFuncType_ == Stone2);
        return *reinterpret_cast<const SatFuncStone2<TableType>*>(satFuncsData_.stone2);
    }
#pragma GCC diagnostic pop

    template <class FluidState>
    void evalKr(const FluidState& fluidState, double* kr) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKr(fluidState, kr)); }

    template <class FluidState>
    void evalKrDeriv(const FluidState& fluidState, double* kr, double* dkrds) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKrDeriv(fluidState, kr, dkrds)); }

    template <class FluidState>
    void evalPc(const FluidState& fluidState, double* pc) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPc(fluidState, pc)); }

    template <class FluidState>
    void evalPcDeriv(const FluidState& fluidState, double* pc, double* dpcds) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPcDeriv(fluidState, pc, dpcds)); }

    template <class FluidState>
    void evalKr(const FluidState& fluidState, double* kr, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKr(fluidState, kr, epst)); }

    template <class FluidState>
    void evalKr(const FluidState& fluidState, double* kr, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKr(fluidState, kr, epst, epst_hyst, sat_hyst)); }

    template <class FluidState>
    void evalKrDeriv(const FluidState& fluidState, double* kr, double* dkrds, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKrDeriv(fluidState, kr, dkrds, epst)); }

    template <class FluidState>
    void evalKrDeriv(const FluidState& fluidState, double* kr, double* dkrds, const EPSTransforms* epst, const EPSTransforms* epst_hyst, const SatHyst* sat_hyst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalKrDeriv(fluidState, kr, dkrds, epst, epst_hyst, sat_hyst)); }

    template <class FluidState>
    void evalPc(const FluidState& fluidState, double* pc, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPc(fluidState, pc, epst)); }

    template <class FluidState>
    void evalPcDeriv(const FluidState& fluidState, double* pc, double* dpcds, const EPSTransforms* epst) const
    { OPM_MULTIPLEXER_SATFUNC_CALL(satFunc.evalPcDeriv(fluidState, pc, dpcds, epst)); }

// undefine the helper macro here because we don't need it anymore
#undef OPM_MULTIPLEXER_SATFUNC_CALL

private:
    bool wasInitialized_;
    SatFuncType satFuncType_;
    union {
        char gwseg[sizeof(SatFuncGwseg<TableType>)];
        char simple[sizeof(SatFuncSimple<TableType>)];
        char stone2[sizeof(SatFuncStone2<TableType>)];
    } satFuncsData_  __attribute__((aligned(sizeof(double))));
};
} // namespace Opm

#endif // OPM_SAT_FUNC_MULTIPLEXER_HPP
