/*
  Copyright 2015 Statoil ASA.

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

#ifndef OPM_PINCHPROCESSOR_HEADER_INCLUDED
#define OPM_PINCHPROCESSOR_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/PinchMode.hpp>
#include <opm/core/utility/Units.hpp>
#include <array>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <limits>

namespace Opm
{

    template <class Grid>
    class PinchProcessor
    {
    public:
        /// \brief Create a Pinch processor.
        /// \param[in]    minpvValue   value in MINPV keyword
        /// \param[in]    thickness    item 2 in PINCH keyword
        /// \param[in]    transMode    item 4 in PINCH keyword
        /// \param[in]    multzMode    item 5 in PINCH keyword
        PinchProcessor(const double minpvValue,
                       const double thickness,
                       const PinchMode::ModeEnum transMode,
                       const PinchMode::ModeEnum multzMode);
        /// Generate NNCs for cells which pv is less than MINPV.
        /// \param[in]    Grid    cpgrid or unstructured grid
        /// \param[in]    htrans  half cell transmissibility, size is number of cellfaces.
        /// \param[in]    multz   Z+ transmissibility multiplier for all active cells
        /// \param[in]    pv      pore volume for all the cartesian cells
        /// \param[in]    nnc     non-neighbor connection class
        /// Algorithm:
        /// 1. Mark all the cells which pv less than minpvValue.
        /// 2. Find out proper pinchouts column and associate top and bottom cells.
        /// 3. Compute transmissibility for nncs.
        /// 4. Apply multz due to different multz options.
        void process(const Grid& grid,
                     const std::vector<double>& htrans,
                     const std::vector<int>& actnum,
                     const std::vector<double>& multz,
                     const std::vector<double>& pv,
                     NNC& nnc);

    private:
        double minpvValue_;
        double thickness_;
        PinchMode::ModeEnum transMode_;
        PinchMode::ModeEnum multzMode_;
        
        /// Mark minpved cells.
        std::vector<int> getMinpvCells_(const std::vector<int>& actnum,
                                        const std::vector<double>& pv);

        /// Get the interface for two cells.
        int interface_(const Grid& grid,
                       const int cellIdx1,
                       const int cellIdx2);

        /// Get the proper face for one cell.
        int interface_(const Grid& grid,
                       const int cellIdx,
                       const Opm::FaceDir::DirEnum& faceDir);

        /// Get pinchouts column.
        std::vector<std::vector<int> >
        getPinchoutsColumn_(const Grid& grid,
                            const std::vector<int>& actnum,
                            const std::vector<double>& pv);
        
        /// Get global cell index.
        int getGlobalIndex_(const int i, const int j, const int k, const int* dims);

        /// Get cartesian index.
        std::array<int, 3> getCartIndex_(const int idx,
                                         const int* dims);

        /// Compute transmissibility for nnc.
        std::vector<double> transCompute_(const Grid& grid,
                                          const std::vector<double>& htrans,
                                          const std::vector<int>& pinCells,
                                          const std::vector<int>& pinFaces);

        /// Get map between half-trans index and the pair of face index and cell index.
        std::vector<int> getHfIdxMap_(const Grid& grid);
        
        /// Get active cell index.
        int getActiveCellIdx_(const Grid& grid,
                              const int globalIdx);

        /// Item 4 in PINCH keyword. 
        void transTopbot_(const Grid& grid,
                          const std::vector<double>& htrans,
                          const std::vector<int>& actnum,
                          const std::vector<double>& multz,
                          const std::vector<double>& pv,
                          NNC& nnc);
        
        /// Item 5 in PINCH keyword.
        std::unordered_multimap<int, double> multzOptions_(const Grid& grid,
                                                           const std::vector<int>& pinCells,
                                                           const std::vector<int>& pinFaces,
                                                           const std::vector<double>& multz,
                                                           const std::vector<std::vector<int> >& seg);

        /// Apply multz vector to face transmissibility.
        void applyMultz_(std::vector<double>& trans,
                         const std::unordered_multimap<int, double>& multzmap);

    };



    template <class Grid>
    inline PinchProcessor<Grid>::PinchProcessor(const double minpv,
                                                const double thickness,
                                                const PinchMode::ModeEnum transMode,
                                                const PinchMode::ModeEnum multzMode)
    {
        minpvValue_ = minpv;
        thickness_  = thickness;
        transMode_  = transMode;
        multzMode_  = multzMode;
    }

    

    template <class Grid>
    inline int PinchProcessor<Grid>::getGlobalIndex_(const int i, const int j, const int k, const int* dims)
    {
        return i + dims[0] * (j + dims[1] * k);
    }



    template <class Grid>
    inline std::array<int, 3> PinchProcessor<Grid>::getCartIndex_(const int idx,
                                                                  const int* dims)
    {
        std::array<int, 3> ijk;
        ijk[0] = (idx % dims[0]);
        ijk[1] = ((idx / dims[0]) % dims[1]);
        ijk[2] = ((idx / dims[0]) / dims[1]);

        return ijk;
    }



    template<class Grid>
    inline int PinchProcessor<Grid>::interface_(const Grid& grid,
                                                const int cellIdx1,
                                                const int cellIdx2)
    {
        const auto cell_faces = Opm::UgGridHelpers::cell2Faces(grid);
        int commonFace = -1;
        const int actCellIdx1 = getActiveCellIdx_(grid, cellIdx1);
        const int actCellIdx2 = getActiveCellIdx_(grid, cellIdx2);
        const auto cellFacesRange1 = cell_faces[actCellIdx1];
        const auto cellFacesRange2 = cell_faces[actCellIdx2];
        for (const auto& f1 : cellFacesRange1) {
            for (const auto& f2 : cellFacesRange2) {
                if (f1 == f2) {
                    commonFace = f1;
                    break;
                }
            }
        }

        if (commonFace == -1) {
            const auto dims = Opm::UgGridHelpers::cartDims(grid);
            const auto ijk1 = getCartIndex_(cellIdx1, dims);
            const auto ijk2 = getCartIndex_(cellIdx2, dims);

            OPM_THROW(std::logic_error, "Couldn't find the common face for cell " 
                      << cellIdx1<< "("<<ijk1[0]<<","<<ijk1[1]<<","<<ijk1[2]<<")"
                      << " and " << cellIdx2<<"("<<ijk2[0]<<","<<ijk2[1]<<","<<ijk2[2]<<")");
        }

        return commonFace;
    }



    template<class Grid>
    inline int PinchProcessor<Grid>::interface_(const Grid& grid,
                                                const int cellIdx,
                                                const Opm::FaceDir::DirEnum& faceDir)
    {
        const auto actCellIdx = getActiveCellIdx_(grid, cellIdx);
        const auto cell_faces = Opm::UgGridHelpers::cell2Faces(grid);
        const auto cellFacesRange = cell_faces[actCellIdx];
        int faceIdx = -1;
        for (auto cellFaceIter = cellFacesRange.begin(); cellFaceIter != cellFacesRange.end(); ++cellFaceIter) {
            int tag = Opm::UgGridHelpers::faceTag(grid, cellFaceIter);
            if ( (faceDir == Opm::FaceDir::ZMinus && tag == 4) || (faceDir == Opm::FaceDir::ZPlus && tag == 5) ) {
                faceIdx = *cellFaceIter;
            }
        }
        
        if (faceIdx == -1) {
            OPM_THROW(std::logic_error, "Couldn't find the face for cell ." << cellIdx);
        }
        
        return faceIdx;
    }



    template<class Grid>
    inline std::vector<int> PinchProcessor<Grid>::getMinpvCells_(const std::vector<int>& actnum,
                                                                 const std::vector<double>& pv)
    {
        std::vector<int> minpvCells(pv.size(), 0);
        for (int idx = 0; idx < static_cast<int>(pv.size()); ++idx) {
            if (actnum[idx]) {
                if (pv[idx] < minpvValue_) {
                    minpvCells[idx] = 1;
                }
            }
        }     
        return minpvCells;
    }



    template<class Grid>
    inline std::vector<int> PinchProcessor<Grid>::getHfIdxMap_(const Grid& grid)
    {
        std::vector<int> hf_ix(2*Opm::UgGridHelpers::numFaces(grid), -1);
        const auto& f2c = Opm::UgGridHelpers::faceCells(grid);
        const auto& cf = Opm::UgGridHelpers::cell2Faces(grid);
        
        for (int c = 0, i = 0; c < Opm::UgGridHelpers::numCells(grid); ++c) {
            for (const auto& f: cf[c]) {
                const auto off = 0 + (f2c(f, 0) != c);
                hf_ix[2*f + off] = i++;
            }
        }
        return hf_ix;
    }



    template<class Grid>
    inline int PinchProcessor<Grid>::getActiveCellIdx_(const Grid& grid,
                                                       const int globalIdx)
    {
        const int nc = Opm::UgGridHelpers::numCells(grid);
        const int* global_cell = Opm::UgGridHelpers::globalCell(grid);
        int idx = -1;
        for (int i = 0; i < nc; ++i) {
            if (global_cell[i] == globalIdx) {
                idx = i;
                break;
            }
        }
        return idx;
    }



    template<class Grid>
    inline std::vector<double> PinchProcessor<Grid>::transCompute_(const Grid& grid,
                                                                   const std::vector<double>& htrans,
                                                                   const std::vector<int>& pinCells,
                                                                   const std::vector<int>& pinFaces)
    {
        const int nc = Opm::UgGridHelpers::numCells(grid);
        const int nf = Opm::UgGridHelpers::numFaces(grid);
        std::vector<double> trans(nf, 0);
        int cellFaceIdx = 0;
        auto cell_faces = Opm::UgGridHelpers::cell2Faces(grid);
        const auto& hfmap = getHfIdxMap_(grid); 
        const auto& f2c = Opm::UgGridHelpers::faceCells(grid);
        for (int cellIdx = 0; cellIdx < nc; ++cellIdx) {
            auto cellFacesRange = cell_faces[cellIdx];
            for (auto cellFaceIter = cellFacesRange.begin(); cellFaceIter != cellFacesRange.end(); ++cellFaceIter, ++cellFaceIdx) {
                const int faceIdx = *cellFaceIter;
                const auto pos = std::find(pinFaces.begin(), pinFaces.end(), faceIdx);
                if (pos == pinFaces.end()) {
                    trans[faceIdx] += 1. / htrans[cellFaceIdx];
                } else {
                    const int idx1 = std::distance(std::begin(pinFaces), pos);
                    int idx2;
                    if (idx1 % 2 == 0) {
                        idx2 = idx1 + 1;
                    } else {
                        idx2 = idx1 - 1;
                    }
                    const int f1 = hfmap[2*pinFaces[idx1] + (f2c(pinFaces[idx1], 0) != getActiveCellIdx_(grid, pinCells[idx1]))];
                    const int f2 = hfmap[2*pinFaces[idx2] + (f2c(pinFaces[idx2], 0) != getActiveCellIdx_(grid, pinCells[idx2]))];
                    trans[faceIdx] = (1. / htrans[f1] + 1. / htrans[f2]);
                    trans[pinFaces[idx2]] = trans[faceIdx];
                }
            }
        }

        for (auto f = 0; f < nf; ++f) {
            trans[f] = 1. / trans[f];
        }

        return trans;
    }



    template<class Grid>
    inline std::vector<std::vector<int>> PinchProcessor<Grid>::getPinchoutsColumn_(const Grid& grid,
                                                                                   const std::vector<int>& actnum,
                                                                                   const std::vector<double>& pv)
    {
        const int* dims = Opm::UgGridHelpers::cartDims(grid);
        std::vector<int> minpvCells = getMinpvCells_(actnum, pv);
        std::vector<std::vector<int>> segment;
        for (int z = 0; z < dims[2]; ++z) {
            for (int y = 0; y < dims[1]; ++y) {
                for (int x = 0; x < dims[0]; ++x) {
                    const int c = getGlobalIndex_(x, y, z, dims);
                    std::vector<int> seg;
                    if (minpvCells[c]) {
                        seg.push_back(c);
                        minpvCells[c] = 0;
                        for (int zz = z+1; zz < dims[2]; ++zz) {
                            const int cc =  getGlobalIndex_(x, y, zz, dims);
                            if (minpvCells[cc]) {
                                seg.push_back(cc);
                                minpvCells[cc] = 0;
                            } else {
                                break;
                            }
                        }
                        segment.push_back(seg);
                    }
                }
            }
        }

        return segment;
    }



    template<class Grid>
    inline void PinchProcessor<Grid>::transTopbot_(const Grid& grid,
                                                   const std::vector<double>& htrans,
                                                   const std::vector<int>& actnum,
                                                   const std::vector<double>& multz,
                                                   const std::vector<double>& pv,
                                                   NNC& nnc)
    {
        const int* dims = Opm::UgGridHelpers::cartDims(grid);
        std::vector<int> pinFaces;
        std::vector<int> pinCells;
        std::vector<std::vector<int> > newSeg;
        auto minpvSeg = getPinchoutsColumn_(grid, actnum, pv);
        for (auto& seg : minpvSeg) {
            std::array<int, 3> ijk1 = getCartIndex_(seg.front(), dims);
            std::array<int, 3> ijk2 = getCartIndex_(seg.back(), dims);
            auto tmp = seg;
            if ((ijk1[2]-1) >= 0 && (ijk2[2]+1) < dims[2]) {
                int topCell = getGlobalIndex_(ijk1[0], ijk1[1], ijk1[2]-1, dims);
                int botCell = getGlobalIndex_(ijk2[0], ijk2[1], ijk2[2]+1, dims);
                /// for any segments, we need to find the active top and bottom cells.
                /// if the original segment's top and bottom is inactive, we need to lookup
                /// the column until they're found otherwise just ignore this segment.
                if (!actnum[topCell]) {
                    for (int topk = ijk1[2]-2; topk > 0; --topk) {
                        topCell = getGlobalIndex_(ijk1[0], ijk1[1], topk, dims);
                        if (actnum[topCell]) {
                            break;
                        }
                    }
                }
                pinFaces.push_back(interface_(grid, topCell, Opm::FaceDir::ZPlus));
                pinCells.push_back(topCell);

                tmp.insert(tmp.begin(), topCell);
                newSeg.push_back(tmp);
                if (!actnum[botCell]) {
                    for (int botk = ijk2[2]+2; botk < dims[2]; ++botk) {
                        botCell = getGlobalIndex_(ijk2[0], ijk2[1], botk, dims);
                        if (actnum[botCell]) {
                            break;

                        }
                    }
                }
                pinFaces.push_back(interface_(grid, botCell, Opm::FaceDir::ZMinus));
                pinCells.push_back(botCell);

            }
        }

        auto faceTrans = transCompute_(grid, htrans, pinCells, pinFaces);
        auto multzmap = multzOptions_(grid, pinCells, pinFaces, multz, newSeg);
        applyMultz_(faceTrans, multzmap);
        for (int i = 0; i < static_cast<int>(pinCells.size())/2; ++i) {
            nnc.addNNC(static_cast<int>(pinCells[2*i]), static_cast<int>(pinCells[2*i+1]), faceTrans[pinFaces[2*i]]);
        }
    }
        

 

    template<class Grid>
    inline std::unordered_multimap<int, double> PinchProcessor<Grid>::multzOptions_(const Grid& grid,
                                                                                    const std::vector<int>& pinCells,
                                                                                    const std::vector<int>& pinFaces,
                                                                                    const std::vector<double>& multz,
                                                                                    const std::vector<std::vector<int> >& segs)
    {
        std::unordered_multimap<int, double> multzmap;
        if (multzMode_ == PinchMode::ModeEnum::TOP) {
            for (int i = 0; i < static_cast<int>(pinFaces.size())/2; ++i) {
                multzmap.insert(std::make_pair(pinFaces[2*i], multz[getActiveCellIdx_(grid, pinCells[2*i])]));
                multzmap.insert(std::make_pair(pinFaces[2*i+1],multz[getActiveCellIdx_(grid, pinCells[2*i])]));
            }
        } else if (multzMode_ == PinchMode::ModeEnum::ALL) {
            for (auto& seg : segs) {
                //find the min multz in seg cells.
                auto multzValue = std::numeric_limits<double>::max();
                for (auto& cellIdx : seg) {
                    auto activeIdx = getActiveCellIdx_(grid, cellIdx);
                    if (activeIdx != -1) {
                        multzValue = std::min(multzValue, multz[activeIdx]);
                    }
                }
                //find the right face.
                auto index = std::distance(std::begin(pinCells), std::find(pinCells.begin(), pinCells.end(), seg.front()));
                multzmap.insert(std::make_pair(pinFaces[index], multzValue));
                multzmap.insert(std::make_pair(pinFaces[index+1], multzValue));
            }
        }

        return multzmap;
    }



    template<class Grid>
    inline void PinchProcessor<Grid>::applyMultz_(std::vector<double>& trans,
                                                  const std::unordered_multimap<int, double>& multzmap)
    {
        for (auto& x : multzmap) {
            trans[x.first] *= x.second;
        }
    }



    template<class Grid>
    inline void PinchProcessor<Grid>::process(const Grid& grid,
                                              const std::vector<double>& htrans,
                                              const std::vector<int>& actnum,
                                              const std::vector<double>& multz,
                                              const std::vector<double>& pv,
                                              NNC& nnc)
    {
        transTopbot_(grid, htrans, actnum, multz, pv, nnc);
    }

} // namespace Opm
#endif // OPM_PINCHPROCESSOR_HEADER_INCLUDED
