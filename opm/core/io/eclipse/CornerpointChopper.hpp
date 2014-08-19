/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_CORNERPOINTCHOPPER_HEADER_INCLUDED
#define OPM_CORNERPOINTCHOPPER_HEADER_INCLUDED

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/Deck/DeckDoubleItem.hpp>
#include <opm/parser/eclipse/Deck/DeckIntItem.hpp>
#include <opm/parser/eclipse/Deck/DeckStringItem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

namespace Opm
{

    class CornerPointChopper
    {
    public:
        CornerPointChopper(const std::string& file)
        {
            Opm::ParserPtr parser(new Opm::Parser());
            deck_ = parser->parseFile(file);

            metricUnits_.reset(Opm::UnitSystem::newMETRIC());

            Opm::DeckRecordConstPtr specgridRecord = deck_->getKeyword("SPECGRID")->getRecord(0);
            dims_[0] = specgridRecord->getItem("NX")->getInt(0);
            dims_[1] = specgridRecord->getItem("NY")->getInt(0);
            dims_[2] = specgridRecord->getItem("NZ")->getInt(0);

            int layersz = 8*dims_[0]*dims_[1];
            const std::vector<double>& ZCORN = deck_->getKeyword("ZCORN")->getRawDoubleData();
            botmax_ = *std::max_element(ZCORN.begin(), ZCORN.begin() + layersz/2);
            topmin_ = *std::min_element(ZCORN.begin() + dims_[2]*layersz - layersz/2,
                                        ZCORN.begin() + dims_[2]*layersz);

            abszmax_ = *std::max_element(ZCORN.begin(), ZCORN.end());
            abszmin_ = *std::min_element(ZCORN.begin(), ZCORN.end());

            std::cout << "Parsed grdecl file with dimensions ("
                      << dims_[0] << ", " << dims_[1] << ", " << dims_[2] << ")" << std::endl;
        }




        const int* dimensions() const
        {
            return dims_;
        }




        const int* newDimensions() const
        {
            return new_dims_;
        }




        const std::pair<double, double> zLimits() const
        {
            return std::make_pair(botmax_, topmin_);
        }

        const std::pair<double, double> abszLimits() const
        {
            return std::make_pair(abszmin_, abszmax_);
        }


        void verifyInscribedShoebox(int imin, int ilen, int imax,
                                    int jmin, int jlen, int jmax,
                                    double zmin, double zlen, double zmax)
        {
            if (imin < 0) {
                std::cerr << "Error! imin < 0 (imin = " << imin << ")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (ilen > dims_[0]) {
                std::cerr << "Error! ilen larger than grid (ilen = " << ilen <<")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (imax > dims_[0]) {
                std::cerr << "Error! imax larger than input grid (imax = " << imax << ")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (jmin < 0) {
                std::cerr << "Error! jmin < 0 (jmin = " << jmin << ")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (jlen > dims_[1]) {
                std::cerr << "Error! jlen larger than grid (jlen = " << jlen <<")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (jmax > dims_[1]) {
                std::cerr << "Error! jmax larger than input grid (jmax = " << jmax << ")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (zmin < abszmin_) {
                std::cerr << "Error! zmin ("<< zmin << ") less than minimum ZCORN value ("<< abszmin_ << ")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (zmax > abszmax_) {
                std::cerr << "Error! zmax ("<< zmax << ") larger than maximal ZCORN value ("<< abszmax_ << ")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
            if (zlen > (abszmax_ - abszmin_)) {
                std::cerr << "Error! zlen ("<< zlen <<") larger than maximal ZCORN (" << abszmax_ << ") minus minimal ZCORN ("<< abszmin_ <<")\n";
                throw std::runtime_error("Inconsistent user input.");
            }
        }

        void chop(int imin, int imax, int jmin, int jmax, double zmin, double zmax, bool resettoorigin=true)
        {
            new_dims_[0] = imax - imin;
            new_dims_[1] = jmax - jmin;

            // Filter the coord field
            const std::vector<double>& COORD = deck_->getKeyword("COORD")->getRawDoubleData();
            int num_coord = COORD.size();
            if (num_coord != 6*(dims_[0] + 1)*(dims_[1] + 1)) {
                std::cerr << "Error! COORD size (" << COORD.size() << ") not consistent with SPECGRID\n";
                throw std::runtime_error("Inconsistent COORD and SPECGRID.");
            }
            int num_new_coord = 6*(new_dims_[0] + 1)*(new_dims_[1] + 1);
            double x_correction = COORD[6*((dims_[0] + 1)*jmin + imin)];
            double y_correction = COORD[6*((dims_[0] + 1)*jmin + imin) + 1];
            new_COORD_.resize(num_new_coord, 1e100);
            for (int j = jmin; j < jmax + 1; ++j) {
                for (int i = imin; i < imax + 1; ++i) {
                    int pos = (dims_[0] + 1)*j + i;
                    int new_pos = (new_dims_[0] + 1)*(j-jmin) + (i-imin);
                    // Copy all 6 coordinates for a pillar.
                    std::copy(COORD.begin() + 6*pos, COORD.begin() + 6*(pos + 1), new_COORD_.begin() + 6*new_pos);
                    if (resettoorigin) {
                        // Substract lowest x value from all X-coords, similarly for y, and truncate in z-direction
                      new_COORD_[6*new_pos]     -= x_correction;
                      new_COORD_[6*new_pos + 1] -= y_correction;
                      new_COORD_[6*new_pos + 2]  = 0;
                      new_COORD_[6*new_pos + 3] -= x_correction;
                      new_COORD_[6*new_pos + 4] -= y_correction;
                      new_COORD_[6*new_pos + 5]  = zmax-zmin;
                    }
                }
            }

            // Get the z limits, check if they must be changed to make a shoe-box.
            // This means that zmin must be greater than or equal to the highest
            // coordinate of the bottom surface, while zmax must be less than or
            // equal to the lowest coordinate of the top surface.
            int layersz = 8*dims_[0]*dims_[1];
            const std::vector<double>& ZCORN = deck_->getKeyword("ZCORN")->getRawDoubleData();
            int num_zcorn = ZCORN.size();
            if (num_zcorn != layersz*dims_[2]) {
                std::cerr << "Error! ZCORN size (" << ZCORN.size() << ") not consistent with SPECGRID\n";
                throw std::runtime_error("Inconsistent ZCORN and SPECGRID.");
            }

            zmin = std::max(zmin, botmax_);
            zmax = std::min(zmax, topmin_);
            if (zmin >= zmax) {
                std::cerr << "Error: zmin >= zmax (zmin = " << zmin << ", zmax = " << zmax << ")\n";
                throw std::runtime_error("zmin >= zmax");
            }
            std::cout << "Chopping subsample,  i: (" << imin << "--" << imax << ")  j: (" << jmin << "--" << jmax << ")   z: (" << zmin << "--" << zmax << ")" <<  std::endl;

            // We must find the maximum and minimum k value for the given z limits.
            // First, find the first layer with a z-coordinate strictly above zmin.
            int kmin = -1;
            for (int k = 0; k < dims_[2]; ++k) {
                double layer_max = *std::max_element(ZCORN.begin() + k*layersz, ZCORN.begin() + (k + 1)*layersz);
                if (layer_max > zmin) {
                    kmin = k;
                    break;
                }
            }
            // Then, find the last layer with a z-coordinate strictly below zmax.
            int kmax = -1;
            for (int k = dims_[2]; k > 0; --k) {
                double layer_min = *std::min_element(ZCORN.begin() + (k - 1)*layersz, ZCORN.begin() + k*layersz);
                if (layer_min < zmax) {
                    kmax = k;
                    break;
                }
            }
            new_dims_[2] = kmax - kmin;

            // Filter the ZCORN field, build mapping from new to old cells.
            double z_origin_correction = 0.0;
            if (resettoorigin) {
                z_origin_correction = zmin;
            }
            new_ZCORN_.resize(8*new_dims_[0]*new_dims_[1]*new_dims_[2], 1e100);
            new_to_old_cell_.resize(new_dims_[0]*new_dims_[1]*new_dims_[2], -1);
            int cellcount = 0;
            int delta[3] = { 1, 2*dims_[0], 4*dims_[0]*dims_[1] };
            int new_delta[3] = { 1, 2*new_dims_[0], 4*new_dims_[0]*new_dims_[1] };
            for (int k = kmin; k < kmax; ++k) {
                for (int j = jmin; j < jmax; ++j) {
                    for (int i = imin; i < imax; ++i) {
                        new_to_old_cell_[cellcount++] = dims_[0]*dims_[1]*k + dims_[0]*j + i;
                        int old_ix = 2*(i*delta[0] + j*delta[1] + k*delta[2]);
                        int new_ix = 2*((i-imin)*new_delta[0] + (j-jmin)*new_delta[1] + (k-kmin)*new_delta[2]);
                        int old_indices[8] = { old_ix, old_ix + delta[0],
                                               old_ix + delta[1], old_ix + delta[1] + delta[0],
                                               old_ix + delta[2], old_ix + delta[2] + delta[0],
                                               old_ix + delta[2] + delta[1], old_ix + delta[2] + delta[1] + delta[0] };
                        int new_indices[8] = { new_ix, new_ix + new_delta[0],
                                               new_ix + new_delta[1], new_ix + new_delta[1] + new_delta[0],
                                               new_ix + new_delta[2], new_ix + new_delta[2] + new_delta[0],
                                               new_ix + new_delta[2] + new_delta[1], new_ix + new_delta[2] + new_delta[1] + new_delta[0] };
                        for (int cc = 0; cc < 8; ++cc) {
                            new_ZCORN_[new_indices[cc]] = std::min(zmax, std::max(zmin, ZCORN[old_indices[cc]])) - z_origin_correction;
                        }
                    }
                }
            }

            filterIntegerField("ACTNUM", new_ACTNUM_);
            filterDoubleField("PORO", new_PORO_);
            filterDoubleField("NTG", new_NTG_);
            filterDoubleField("SWCR", new_SWCR_);
            filterDoubleField("SOWCR", new_SOWCR_);
            filterDoubleField("PERMX", new_PERMX_);
            filterDoubleField("PERMY", new_PERMY_);
            filterDoubleField("PERMZ", new_PERMZ_);
            filterIntegerField("SATNUM", new_SATNUM_);
        }

        /// Return a sub-deck with fields corresponding to the selected subset.
        Opm::DeckConstPtr subDeck()
        {
            Opm::DeckPtr subDeck(new Opm::Deck);

            Opm::DeckKeywordPtr specGridKw(new Opm::DeckKeyword("SPECGRID"));
            Opm::DeckRecordPtr specGridRecord(new Opm::DeckRecord());

            Opm::DeckIntItemPtr nxItem(new Opm::DeckIntItem("NX"));
            Opm::DeckIntItemPtr nyItem(new Opm::DeckIntItem("NY"));
            Opm::DeckIntItemPtr nzItem(new Opm::DeckIntItem("NZ"));
            Opm::DeckIntItemPtr numresItem(new Opm::DeckIntItem("NUMRES"));
            Opm::DeckStringItemPtr coordTypeItem(new Opm::DeckStringItem("COORD_TYPE"));

            nxItem->push_back(new_dims_[0]);
            nyItem->push_back(new_dims_[1]);
            nzItem->push_back(new_dims_[2]);
            numresItem->push_back(1);
            coordTypeItem->push_back("F");

            specGridRecord->addItem(nxItem);
            specGridRecord->addItem(nyItem);
            specGridRecord->addItem(nzItem);
            specGridRecord->addItem(numresItem);
            specGridRecord->addItem(coordTypeItem);

            specGridKw->addRecord(specGridRecord);

            subDeck->addKeyword(specGridKw);
            addDoubleKeyword_(subDeck, "COORD", /*dimension=*/"Length", new_COORD_);
            addDoubleKeyword_(subDeck, "ZCORN", /*dimension=*/"Length", new_ZCORN_);
            addIntKeyword_(subDeck, "ACTNUM", new_ACTNUM_);
            addDoubleKeyword_(subDeck, "PORO", /*dimension=*/"1", new_PORO_);
            addDoubleKeyword_(subDeck, "NTG", /*dimension=*/"1", new_NTG_);
            addDoubleKeyword_(subDeck, "SWCR", /*dimension=*/"1", new_SWCR_);
            addDoubleKeyword_(subDeck, "SOWCR", /*dimension=*/"1", new_SOWCR_);
            addDoubleKeyword_(subDeck, "PERMX", /*dimension=*/"Permeability", new_PERMX_);
            addDoubleKeyword_(subDeck, "PERMY", /*dimension=*/"Permeability", new_PERMY_);
            addDoubleKeyword_(subDeck, "PERMZ", /*dimension=*/"Permeability", new_PERMZ_);
            addIntKeyword_(subDeck, "SATNUM", new_SATNUM_);
            return subDeck;
        }
        void writeGrdecl(const std::string& filename)
        {
            // Output new versions of SPECGRID, COORD, ZCORN, ACTNUM, PERMX, PORO, SATNUM.
            std::ofstream out(filename.c_str());
            if (!out) {
                std::cerr << "Could not open file " << filename << "\n";
                throw std::runtime_error("Could not open output file.");
            }
            out << "SPECGRID\n" << new_dims_[0] << ' ' << new_dims_[1] << ' ' << new_dims_[2]
                << " 1 F\n/\n\n";

            outputField(out, new_COORD_, "COORD", /* nl = */ 3);
            outputField(out, new_ZCORN_, "ZCORN", /* nl = */ 4);
            outputField(out, new_ACTNUM_, "ACTNUM");
            outputField(out, new_PORO_, "PORO", 4);
            if (hasNTG()) {outputField(out, new_NTG_, "NTG", 4);}
            if (hasSWCR()) {outputField(out, new_SWCR_, "SWCR", 4);}
            if (hasSOWCR()) {outputField(out, new_SOWCR_, "SOWCR", 4);}
            outputField(out, new_PERMX_, "PERMX", 4);
            outputField(out, new_PERMY_, "PERMY", 4);
            outputField(out, new_PERMZ_, "PERMZ", 4);
            outputField(out, new_SATNUM_, "SATNUM");
        }
        bool hasNTG() const {return !new_NTG_.empty(); }
        bool hasSWCR() const {return !new_SWCR_.empty(); }
        bool hasSOWCR() const {return !new_SOWCR_.empty(); }

    private:
        Opm::DeckConstPtr deck_;
        std::shared_ptr<Opm::UnitSystem> metricUnits_;

        double botmax_;
        double topmin_;
        double abszmin_;
        double abszmax_;
        std::vector<double> new_COORD_;
        std::vector<double> new_ZCORN_;
        std::vector<int> new_ACTNUM_;
        std::vector<double> new_PORO_;
        std::vector<double> new_NTG_;
        std::vector<double> new_SWCR_;
        std::vector<double> new_SOWCR_;
        std::vector<double> new_PERMX_;
        std::vector<double> new_PERMY_;
        std::vector<double> new_PERMZ_;
        std::vector<int> new_SATNUM_;
        int dims_[3];
        int new_dims_[3];
        std::vector<int> new_to_old_cell_;

        void addDoubleKeyword_(Opm::DeckPtr subDeck,
                               const std::string& keywordName,
                               const std::string& dimensionString,
                               const std::vector<double>& data)
        {
            if (data.empty())
                return;

            Opm::DeckKeywordPtr dataKw(new Opm::DeckKeyword(keywordName));
            Opm::DeckRecordPtr dataRecord(new Opm::DeckRecord());
            Opm::DeckDoubleItemPtr dataItem(new Opm::DeckDoubleItem("DATA"));

            for (size_t i = 0; i < data.size(); ++i) {
                dataItem->push_back(data[i]);
            }

            std::shared_ptr<const Dimension> dimension = metricUnits_->parse(dimensionString);
            dataItem->push_backDimension(/*active=*/dimension, /*default=*/dimension);

            dataRecord->addItem(dataItem);
            dataKw->addRecord(dataRecord);
            subDeck->addKeyword(dataKw);
        }

        void addIntKeyword_(Opm::DeckPtr subDeck,
                               const std::string& keywordName,
                               const std::vector<int>& data)
        {
            if (data.empty())
                return;

            Opm::DeckKeywordPtr dataKw(new Opm::DeckKeyword(keywordName));
            Opm::DeckRecordPtr dataRecord(new Opm::DeckRecord());
            Opm::DeckIntItemPtr dataItem(new Opm::DeckIntItem("DATA"));

            for (size_t i = 0; i < data.size(); ++i) {
                dataItem->push_back(data[i]);
            }

            dataRecord->addItem(dataItem);
            dataKw->addRecord(dataRecord);
            subDeck->addKeyword(dataKw);
        }

        template <typename T>
        void outputField(std::ostream& os,
                         const std::vector<T>& field,
                         const std::string& keyword,
                         const typename std::vector<T>::size_type nl = 20)
        {
            if (field.empty()) return;

            os << keyword << '\n';

            typedef typename std::vector<T>::size_type sz_t;

            const sz_t n = field.size();
            for (sz_t i = 0; i < n; ++i) {
                os << field[i]
                   << (((i + 1) % nl == 0) ? '\n' : ' ');
            }
            if (n % nl != 0) {
                os << '\n';
            }
            os << "/\n\n";
        }



        template <typename T>
        void filterField(const std::vector<T>& field,
                                std::vector<T>& output_field)
        {
            int sz = new_to_old_cell_.size();
            output_field.resize(sz);
            for (int i = 0; i < sz; ++i) {
                output_field[i] = field[new_to_old_cell_[i]];
            }
        }

        void filterDoubleField(const std::string& keyword, std::vector<double>& output_field)
        {
            if (deck_->hasKeyword(keyword)) {
                const std::vector<double>& field = deck_->getKeyword(keyword)->getRawDoubleData();
                filterField(field, output_field);
            }
        }

        void filterIntegerField(const std::string& keyword, std::vector<int>& output_field)
        {
            if (deck_->hasKeyword(keyword)) {
                const std::vector<int>& field = deck_->getKeyword(keyword)->getIntData();
                filterField(field, output_field);
            }
        }

    };

}




#endif // OPM_CORNERPOINTCHOPPER_HEADER_INCLUDED
