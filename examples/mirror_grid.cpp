/*
  Copyright 2014 Statoil ASA.
  Copyright 2014 Andreas Lauser.

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

#include "config.h"

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <iostream>
#include <fstream>

/**
 * @file mirror_grid.cpp
 * @brief Mirror grid taken from grdecl file
 * 
 * The input grid is mirrored in either the x- or y-direction, resulting in a periodic grid.
 *
 */


/// Print init message in new grid filename
void printInitMessage(std::ofstream& out, const char* origfilename, std::string direction) {
    std::ifstream infile;
    infile.open(origfilename, std::ios::in);
    if (!infile) {
        std::cerr << "Can't open input file " << origfilename << std::endl;
        exit(1);
    }
    // Print init message and copy comments from original grid file
    out << "-- This grdecl file is generated from OPM application 'mirror_grid.cpp'." << std::endl
        << "-- The grid '" << origfilename << "' is mirrored around itself in the " << direction << " direction." << std::endl
        << "-- Thus, the resulting grid should be periodic in the " << direction << "-direction." << std::endl
        << "-- Original comments taken from '" << origfilename << "':" << std::endl;
    std::string nextInLine;
    while (getline(infile, nextInLine)) {
        if (nextInLine.substr(0,2) == "--") {
            out << nextInLine << std::endl;
        }
        else {
            break;
        }
    }
    out << std::endl;   
}

/// Write keyword values to file
template <class T>
void printKeywordValues(std::ofstream& out, std::string keyword, std::vector<T> values, int nCols) {
    out << keyword << std::endl;
    int col = 0;
    typename std::vector<T>::iterator iter;
    for (iter = values.begin(); iter != values.end(); ++iter) {
        out << *iter << " ";            
        ++col;
        // Break line for every nCols entry.
        if (col == nCols) {
            out << std::endl;
            col = 0;
        }
    }
    if (col != 0)
        out << std::endl;
    out << "/" << std::endl << std::endl;
}

// forward declaration
std::vector<double> getMapaxesValues(Opm::DeckConstPtr deck);

/// Mirror keyword MAPAXES in deck
void mirror_mapaxes(Opm::DeckConstPtr deck, std::string direction, std::ofstream& out) {
    // Assumes axis aligned with x/y-direction
    std::cout << "Warning: Keyword MAPAXES not fully understood. Result should be verified manually." << std::endl;
    if (deck->hasKeyword("MAPAXES")) {
        std::vector<double> mapaxes = getMapaxesValues(deck);
        std::vector<double> mapaxes_mirrored = mapaxes;
        // Double the length of the coordinate axis
        if (direction == "x") {
            mapaxes_mirrored[4] = (mapaxes[4]-mapaxes[2])*2 + mapaxes[2];
        }
        else if (direction == "y") {
            mapaxes_mirrored[1] = (mapaxes[1]-mapaxes[3])*2 + mapaxes[3];
        }
        printKeywordValues(out, "MAPAXES", mapaxes_mirrored, 2);
    }
}

/// Mirror keyword SPECGRID in deck
void mirror_specgrid(Opm::DeckConstPtr deck, std::string direction, std::ofstream& out) {
    // We only need to multiply the dimension by 2 in the correct direction.
    Opm::DeckRecordConstPtr specgridRecord = deck->getKeyword("SPECGRID")->getRecord(0);
    std::vector<int> dimensions(3);
    dimensions[0] = specgridRecord->getItem("NX")->getInt(0);
    dimensions[1] = specgridRecord->getItem("NY")->getInt(0);
    dimensions[2] = specgridRecord->getItem("NZ")->getInt(0);
    if (direction == "x")      {dimensions[0] *= 2;}
    else if (direction == "y") {dimensions[1] *= 2;}
    else                       {std::cerr << "Direction should be either x or y" << std::endl; exit(1);}
    out << "SPECGRID" << std::endl << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << " "
        << specgridRecord->getItem("NUMRES")->getInt(0) << " "
        << specgridRecord->getItem("COORD_TYPE")->getString(0) << " "
        << std::endl << "/" << std::endl << std::endl;
}

/// Mirror keyword COORD in deck
void mirror_coord(Opm::DeckConstPtr deck, std::string direction, std::ofstream& out) {
    // We assume uniform spacing in x and y directions and parallel top and bottom faces
    Opm::DeckRecordConstPtr specgridRecord = deck->getKeyword("SPECGRID")->getRecord(0);
    std::vector<int> dimensions(3);
    dimensions[0] = specgridRecord->getItem("NX")->getInt(0);
    dimensions[1] = specgridRecord->getItem("NY")->getInt(0);
    dimensions[2] = specgridRecord->getItem("NZ")->getInt(0);
    std::vector<double> coord = deck->getKeyword("COORD")->getRawDoubleData();
    const int entries_per_pillar = 6;
    std::vector<double> coord_mirrored;
    // Handle the two directions differently due to ordering of the pillars.
    if (direction == "x") {
        // Total entries in mirrored ZCORN. Number of pillars times 6
        const int entries = (2*dimensions[0] + 1) * (dimensions[1] + 1) * entries_per_pillar;
        // Entries per line in x-direction. Number of pillars in x-direction times 6
        const int entries_per_line = entries_per_pillar*(dimensions[0] + 1);
        coord_mirrored.assign(entries, 0.0);
        // Distance between pillars in x-directiion
        const double spacing = coord[entries_per_pillar]-coord[0];
        std::vector<double>::iterator it_new = coord_mirrored.begin();
        std::vector<double>::iterator it_orig;
        // Loop through each pillar line in the x-direction
        for (it_orig = coord.begin(); it_orig != coord.end(); it_orig += entries_per_line) {
            // Copy old pillars
            copy(it_orig, it_orig + entries_per_line, it_new);
            // Add new pillars in between
            it_new += entries_per_line;
            std::vector<double> next_vec(it_orig + entries_per_line - entries_per_pillar, it_orig + entries_per_line);
            for (int r=0; r < dimensions[0]; ++r) {
                next_vec[0] += spacing;
                next_vec[3] += spacing;
                copy(next_vec.begin(), next_vec.end(), it_new);
                it_new += entries_per_pillar;
            }
        }
    }
    else if (direction == "y") {
        // Total entries in mirrored ZCORN. Number of pillars times 6
        const int entries = (dimensions[0] + 1) * (2*dimensions[1] + 1) * entries_per_pillar;
        // Entries per line in y-direction. Number of pillars in y-direction times 6
        const int entries_per_line = entries_per_pillar*(dimensions[0] + 1);
        coord_mirrored.assign(entries, 0.0);
        // Distance between pillars in y-directiion
        const double spacing = coord[entries_per_line + 1]-coord[1];
        std::vector<double>::iterator it_new = coord_mirrored.begin();
        // Copy old pillars
        copy(coord.begin(), coord.end(), it_new);
        // Add new pillars at the end
        it_new += coord.size();
        std::vector<double> next_vec(coord.end() - entries_per_line, coord.end());
        for ( ; it_new != coord_mirrored.end(); it_new += entries_per_line) {
            for (int i = 1; i < entries_per_line; i += 3) {
                next_vec[i] += spacing;
            }
            copy(next_vec.begin(), next_vec.end(), it_new);
        }
    }
    else {
        std::cerr << "Direction should be either x or y" << std::endl;
        exit(1);
    }
    // Write new COORD values to output file
    printKeywordValues(out, "COORD", coord_mirrored, 6);
}

/// Mirror keyword ZCORN in deck
void mirror_zcorn(Opm::DeckConstPtr deck, std::string direction, std::ofstream& out) {
    Opm::DeckRecordConstPtr specgridRecord = deck->getKeyword("SPECGRID")->getRecord(0);
    std::vector<int> dimensions(3);
    dimensions[0] = specgridRecord->getItem("NX")->getInt(0);
    dimensions[1] = specgridRecord->getItem("NY")->getInt(0);
    dimensions[2] = specgridRecord->getItem("NZ")->getInt(0);
    std::vector<double> zcorn = deck->getKeyword("ZCORN")->getRawDoubleData();
    std::vector<double> zcorn_mirrored;
    // Handle the two directions differently due to ordering of the pillars.
    if (direction == "x") {
        // Total entries in mirrored ZCORN. Eight corners per cell.
        const int entries = dimensions[0]*2*dimensions[1]*dimensions[2]*8;
        zcorn_mirrored.assign(entries, 0.0);
        // Entries per line in x-direction. Two for each cell.
        const int entries_per_line = dimensions[0]*2;
        std::vector<double>::iterator it_new = zcorn_mirrored.begin();
        std::vector<double>::iterator it_orig = zcorn.begin();
        // Loop through each line and copy old corner-points and add new (which are the old reversed)
        for ( ; it_orig != zcorn.end(); it_orig += entries_per_line) {
            std::vector<double> next_vec(it_orig, it_orig + entries_per_line);
            std::vector<double> next_reversed = next_vec;
            reverse(next_reversed.begin(), next_reversed.end());
            // Copy old corner-points
            copy(it_orig, it_orig + entries_per_line, it_new);
            it_new += entries_per_line;
            // Add new corner-points
            copy(next_reversed.begin(), next_reversed.end(), it_new);
            it_new += entries_per_line;
        }
    }
    else if (direction == "y") {
        // Total entries in mirrored ZCORN. Eight corners per cell.
        const int entries = dimensions[0]*dimensions[1]*2*dimensions[2]*8;
        zcorn_mirrored.assign(entries, 0.0);
        // Entries per line in x-direction. Two for each cell.
        const int entries_per_line_x = dimensions[0]*2;
        // Entries per layer of corner-points. Four for each cell
        const int entries_per_layer = dimensions[0]*dimensions[1]*4;
        std::vector<double>::iterator it_new = zcorn_mirrored.begin();
        std::vector<double>::iterator it_orig = zcorn.begin();
        // Loop through each layer and copy old corner-points and add new (which are the old reordered) 
        for ( ; it_orig != zcorn.end(); it_orig += entries_per_layer) {
            // Copy old corner-points
            copy(it_orig, it_orig + entries_per_layer, it_new);
            it_new += entries_per_layer;
            // Add new corner-points
            std::vector<double> next_vec(it_orig, it_orig + entries_per_layer);
            std::vector<double> next_reordered(entries_per_layer, 0.0);
            std::vector<double>::iterator it_next = next_vec.end();
            std::vector<double>::iterator it_reordered = next_reordered.begin();
            // Reorder next entries
            for ( ; it_reordered != next_reordered.end(); it_reordered += entries_per_line_x) {
                copy(it_next - entries_per_line_x, it_next, it_reordered);
                it_next -= entries_per_line_x;
            }
            copy(next_reordered.begin(), next_reordered.end(), it_new);
            it_new += entries_per_layer;
        }
    }
    else {
        std::cerr << "Direction should be either x or y" << std::endl;
        exit(1);
    }
    // Write new ZCORN values to output file
    printKeywordValues(out, "ZCORN", zcorn_mirrored, 8);
}

std::vector<int> getKeywordValues(std::string keyword, Opm::DeckConstPtr deck, int /*dummy*/) {
    return deck->getKeyword(keyword)->getIntData();
}

std::vector<double> getKeywordValues(std::string keyword, Opm::DeckConstPtr deck, double /*dummy*/) {
    return deck->getKeyword(keyword)->getRawDoubleData();
}

std::vector<double> getMapaxesValues(Opm::DeckConstPtr deck)
{
    Opm::DeckRecordConstPtr mapaxesRecord = deck->getKeyword("MAPAXES")->getRecord(0);
    std::vector<double> result;
    for (size_t itemIdx = 0; itemIdx < mapaxesRecord->size(); ++itemIdx) {
        Opm::DeckItemConstPtr curItem = mapaxesRecord->getItem(itemIdx);

        for (size_t dataItemIdx = 0; dataItemIdx < curItem->size(); ++dataItemIdx) {
            result.push_back(curItem->getRawDouble(dataItemIdx));
        }
    }
    return result;
}

/// Mirror keywords that have one value for each cell
template <class T>
void mirror_celldata(std::string keyword, Opm::DeckConstPtr deck, std::string direction, std::ofstream& out) {
    if ( ! deck->hasKeyword(keyword)) {
        std::cout << "Ignoring keyword " << keyword << " as it was not found." << std::endl;
        return;
    }
    // Get data from eclipse deck
    Opm::DeckRecordConstPtr specgridRecord = deck->getKeyword("SPECGRID")->getRecord(0);
    std::vector<int> dimensions(3);
    dimensions[0] = specgridRecord->getItem("NX")->getInt(0);
    dimensions[1] = specgridRecord->getItem("NY")->getInt(0);
    dimensions[2] = specgridRecord->getItem("NZ")->getInt(0);
    std::vector<T> values = getKeywordValues(keyword, deck, T(0.0));
    std::vector<T> values_mirrored(2*dimensions[0]*dimensions[1]*dimensions[2], 0.0);
    // Handle the two directions differently due to ordering of the pillars.
    if (direction == "x") {
        typename std::vector<T>::iterator it_orig = values.begin();
        typename std::vector<T>::iterator it_new = values_mirrored.begin();
        // Loop through each line and copy old cell data and add new (which are the old reversed)
        for ( ; it_orig != values.end(); it_orig += dimensions[0]) {
            // Copy old cell data
            copy(it_orig, it_orig + dimensions[0], it_new);
            it_new += dimensions[0];
            // Add new cell data
            std::vector<double> next_vec(it_orig, it_orig + dimensions[0]);
            std::vector<double> next_reversed = next_vec;
            reverse(next_reversed.begin(), next_reversed.end());
            copy(next_reversed.begin(), next_reversed.end(), it_new);
            it_new += dimensions[0];
        }
    }
    else if (direction =="y") {
        typename std::vector<T>::iterator it_orig = values.begin();
        typename std::vector<T>::iterator it_new = values_mirrored.begin();
        // Entries per layer
        const int entries_per_layer = dimensions[0]*dimensions[1];
        // Loop through each layer and copy old cell data and add new (which are the old reordered) 
        for ( ; it_orig != values.end(); it_orig += entries_per_layer) {
            // Copy old cell data
            copy(it_orig, it_orig + entries_per_layer, it_new);
            it_new += entries_per_layer;
            // Add new cell data
            std::vector<T> next_vec(it_orig, it_orig + entries_per_layer);
            std::vector<T> next_reordered(entries_per_layer, 0.0);
            typename std::vector<T>::iterator it_next = next_vec.end();
            typename std::vector<T>::iterator it_reordered = next_reordered.begin();
            // Reorder next entries
            for ( ; it_reordered != next_reordered.end(); it_reordered += dimensions[0]) {
                copy(it_next - dimensions[0], it_next, it_reordered);
                it_next -= dimensions[0];
            }
            copy(next_reordered.begin(), next_reordered.end(), it_new);
            it_new += entries_per_layer;
        }
    }
    else {
        std::cerr << "Direction should be either x or y" << std::endl;
        exit(1);
    }
    // Write new keyword values to output file
    printKeywordValues(out, keyword, values_mirrored, 8);
}


int main(int argc, char** argv)
{
    // Set output precision
    int decimals = 16;

    // Process input parameters
    if (argc != 3) {
        std::cout << "Usage: mirror_grid filename.grdecl direction" << std::endl;
        std::cout << "(replace direction with either x or y)" << std::endl;
        exit(1);
    }
    const char* eclipsefilename = argv[1];
    std::string direction(argv[2]);
    if ( ! ((direction == "x") || (direction == "y")) ) {
        std::cerr << "Unrecognized input parameter for direction: '" << direction 
                  << "'. Should be either x or y (maybe also z later)." << std::endl;
        exit(1);
    }

    // Parse grdecl file
    std::cout << "Parsing grid file '" << eclipsefilename << "' ..." << std::endl;
    Opm::ParserPtr parser(new Opm::Parser);
    Opm::DeckConstPtr deck(parser->parseFile(eclipsefilename));
    if ( ! (deck->hasKeyword("SPECGRID") && deck->hasKeyword("COORD") && deck->hasKeyword("ZCORN")) ) {
        std::cerr << "Grid file " << eclipsefilename << "are missing keywords SPECGRID, COORD or ZCORN!" << std::endl;
        exit(1);
    }

    // Create new grid file
    std::string mirrored_eclipsefilename = std::string(eclipsefilename);
    std::string::size_type last_dot = mirrored_eclipsefilename.find_last_of('.');
    mirrored_eclipsefilename = mirrored_eclipsefilename.substr(0, last_dot) + "_mirrored-" + direction + ".grdecl";
    std::ofstream outfile;
    outfile.open(mirrored_eclipsefilename.c_str(), std::ios::out | std::ios::trunc);
    if (!outfile) {
        std::cerr << "Can't open output file " << mirrored_eclipsefilename << std::endl;
        exit(1);
    }
    outfile.precision(decimals);
    outfile.setf(std::ios::fixed);

    // Print init message
    printInitMessage(outfile, eclipsefilename, direction);
    
    // Mirror keywords
    mirror_mapaxes(deck, direction, outfile);
    mirror_specgrid(deck, direction, outfile);
    mirror_coord(deck, direction, outfile);
    mirror_zcorn(deck, direction, outfile);
    mirror_celldata<int>("ACTNUM", deck, direction, outfile);
    mirror_celldata<double>("PERMX", deck, direction, outfile);
    mirror_celldata<double>("PERMY", deck, direction, outfile);
    mirror_celldata<double>("PERMZ", deck, direction, outfile);
    mirror_celldata<double>("PORO", deck, direction, outfile);
    mirror_celldata<int>("SATNUM", deck, direction, outfile);
    mirror_celldata<double>("NTG", deck, direction, outfile);
    mirror_celldata<double>("SWCR", deck, direction, outfile);
    mirror_celldata<double>("SOWCR", deck, direction, outfile);
}
