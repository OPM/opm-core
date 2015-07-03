// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_SIMULATORSTATE_HEADER_INCLUDED
#define OPM_SIMULATORSTATE_HEADER_INCLUDED

#include <vector>
#include <string>

// forward declaration
struct UnstructuredGrid;

namespace Opm
{
    class SimulatorState
    {
    public:

        virtual void init(const UnstructuredGrid& g, int num_phases);

        virtual void init(int number_of_cells, int number_of_faces, int num_phases);

        enum ExtremalSat { MinSat, MaxSat };

        /// \brief pressure per cell.
        static const int pressureId_ = 0;
        /// \brief temperature per cell.
        static const int temperatureId_ = 1;
        /// \brief The saturation of each phase per cell.
        static const int saturationId_ = 2;

        /// \brief pressure per face.
        static const int facePressureId_ = 0;
        /// \brief The fluxes at the faces.
        static const int faceFluxId_ = 1;

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

        std::vector<double>& pressure    () { return cellData_[ pressureId_ ]; }
        std::vector<double>& temperature () { return cellData_[ temperatureId_ ]; }
        std::vector<double>& facepressure() { return faceData_[ facePressureId_]; }
        std::vector<double>& faceflux    () { return faceData_[ faceFluxId_ ];        }
        std::vector<double>& saturation  () { return cellData_[ saturationId_ ];  }

        const std::vector<double>& pressure    () const { return cellData_[ pressureId_ ];    }
        const std::vector<double>& temperature () const { return cellData_[ temperatureId_ ]; }
        const std::vector<double>& facepressure() const { return faceData_[ facePressureId_]; }
        const std::vector<double>& faceflux    () const { return faceData_[ faceFluxId_ ];        }
        const std::vector<double>& saturation  () const { return cellData_[ saturationId_ ];  }

        /**
         * Compare this state with another, to see if they are different
         * only within a small margin.
         */
        virtual bool equals(const SimulatorState& other,
                            double epsilon = 1e-8) const;

        std::vector< std::vector<double> >& cellData() { return cellData_; }
        const std::vector< std::vector<double> >& cellData() const { return cellData_; }

        std::vector< std::vector<double> >& faceData() { return faceData_; }
        const std::vector< std::vector<double> >& faceData() const { return faceData_; }

    private:
        int num_phases_;

        /// \brief vector containing all registered cell data
        std::vector< std::vector< double > > cellData_;
        /// \brief vector containing all registered face data
        std::vector< std::vector< double > > faceData_;

        /// \brief names for the cell data
        std::vector< std::string > cellDataNames_;
        /// \brief names for the face data
        std::vector< std::string > faceDataNames_;

    protected:
        size_t registerCellData( const std::string& name, const int components, const double initialValue = 0.0 );
        size_t registerFaceData( const std::string& name, const int components, const double initialValue = 0.0 );

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
