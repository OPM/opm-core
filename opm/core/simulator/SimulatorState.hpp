// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_SIMULATORSTATE_HEADER_INCLUDED
#define OPM_SIMULATORSTATE_HEADER_INCLUDED

#include <vector>

// forward declaration
struct UnstructuredGrid;

namespace Opm
{
    class SimulatorState
    {
    public:

        virtual void init(const UnstructuredGrid& g, int num_phases);

        virtual void init(int number_of_cells, int number_of_phases, int num_phases);
        
        enum ExtremalSat { MinSat, MaxSat };

    protected:
        /**
         * Initialize the first saturation to maximum value. This method
         * should be considered deprecated. Avoid to use it!
         *
         * \tparam Props Fluid and rock properties that pertain to this
         *               kind of simulation. Currently, only Blackoil-
         *               and IncompPropertiesInterface are supported.
         */
        template <typename Props>
        void setFirstSat(const std::vector<int>& cells,
                         const Props& props,
                         ExtremalSat es);
    public:
        int numPhases() const { return num_phases_; }

        std::vector<double>& pressure    () { return press_ ; }
        std::vector<double>& facepressure() { return fpress_; }
        std::vector<double>& faceflux    () { return flux_  ; }
        std::vector<double>& saturation  () { return sat_   ; }

        const std::vector<double>& pressure    () const { return press_ ; }
        const std::vector<double>& facepressure() const { return fpress_; }
        const std::vector<double>& faceflux    () const { return flux_  ; }
        const std::vector<double>& saturation  () const { return sat_   ; }

        /**
         * Compare this state with another, to see if they are different
         * only within a small margin.
         */
        virtual bool equals(const SimulatorState& other,
                            double epsilon = 1e-8) const;
    private:
        int num_phases_;
        std::vector<double> press_ ;
        std::vector<double> fpress_;
        std::vector<double> flux_  ;
        std::vector<double> sat_   ;

    protected:
        /**
         * Check if two vectors are equal within a margin.
         *
         * @param epsilon Relative difference that is tolerated for the
         *                vectors to still be considered equal.
         *
         * @return True if every element is within the margin, false if
         *         there is at least one that is not.
         */
        static bool vectorApproxEqual(const std::vector<double>& v1,
                                      const std::vector<double>& v2,
                                      double epsilon);
    };

} // namespace Opm

#endif // OPM_SIMULATORSTATE_HEADER_INCLUDED
