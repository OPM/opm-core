#include <memory>

namespace Opm
{
    template<class CentroidIterator>
    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           int number_of_cells,
                                                           const int* global_cell,
                                                           const int* cart_dims,
                                                           const CentroidIterator& begin_cell_centroids,
                                                           int dimension,
                                                           bool init_rock)
    {
        std::vector<int> compressedToCartesianIdx(number_of_cells);
        for (int cellIdx = 0; cellIdx < number_of_cells; ++cellIdx) {
            if (global_cell) {
                compressedToCartesianIdx[cellIdx] = global_cell[cellIdx];
            }
            else {
                compressedToCartesianIdx[cellIdx] = cellIdx;
            }
        }
        auto materialLawManager = std::make_shared<MaterialLawManager>();
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);

        init(deck, eclState, materialLawManager, number_of_cells, global_cell, cart_dims, begin_cell_centroids, dimension,
             init_rock);
    }

    template<class CentroidIterator>
    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           int number_of_cells,
                                                           const int* global_cell,
                                                           const int* cart_dims,
                                                           const CentroidIterator& begin_cell_centroids,
                                                           int dimension,
                                                           const parameter::ParameterGroup& param,
                                                           bool init_rock)
    {
        std::vector<int> compressedToCartesianIdx(number_of_cells);
        for (int cellIdx = 0; cellIdx < number_of_cells; ++cellIdx) {
            if (global_cell) {
                compressedToCartesianIdx[cellIdx] = global_cell[cellIdx];
            }
            else {
                compressedToCartesianIdx[cellIdx] = cellIdx;
            }
        }
        auto materialLawManager = std::make_shared<MaterialLawManager>();
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);

        init(deck,
             eclState,
             materialLawManager,
             number_of_cells,
             global_cell,
             cart_dims,
             begin_cell_centroids,
             dimension,
             param,
             init_rock);
    }

    template<class CentroidIterator>
    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           std::shared_ptr<MaterialLawManager> materialLawManager,
                                                           int number_of_cells,
                                                           const int* global_cell,
                                                           const int* cart_dims,
                                                           const CentroidIterator& begin_cell_centroids,
                                                           int dimension,
                                                           const parameter::ParameterGroup& param,
                                                           bool init_rock)
    {
        init(deck,
             eclState,
             materialLawManager,
             number_of_cells,
             global_cell,
             cart_dims,
             begin_cell_centroids,
             dimension,
             param,
             init_rock);
    }

    template<class CentroidIterator>
    inline void BlackoilPropertiesFromDeck::init(Opm::DeckConstPtr deck,
                                                 Opm::EclipseStateConstPtr eclState,
                                                 std::shared_ptr<MaterialLawManager> materialLawManager,
                                                 int number_of_cells,
                                                 const int* global_cell,
                                                 const int* cart_dims,
                                                 const CentroidIterator& begin_cell_centroids,
                                                 int dimension,
                                                 bool init_rock)
    {
        // retrieve the cell specific PVT table index from the deck
        // and using the grid...
        extractPvtTableIndex(cellPvtRegionIdx_, deck, number_of_cells, global_cell);

        if (init_rock){
           rock_.init(eclState, number_of_cells, global_cell, cart_dims);
        }
        pvt_.init(deck, eclState, /*numSamples=*/0);
        SaturationPropsFromDeck* ptr
            = new SaturationPropsFromDeck();
        ptr->init(phaseUsageFromDeck(deck), materialLawManager);
        satprops_.reset(ptr);

        if (pvt_.numPhases() != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }

    template<class CentroidIterator>
    inline void BlackoilPropertiesFromDeck::init(Opm::DeckConstPtr deck,
                                                 Opm::EclipseStateConstPtr eclState,
                                                 std::shared_ptr<MaterialLawManager> materialLawManager,
                                                 int number_of_cells,
                                                 const int* global_cell,
                                                 const int* cart_dims,
                                                 const CentroidIterator& begin_cell_centroids,
                                                 int dimension,
                                                 const parameter::ParameterGroup& param,
                                                 bool init_rock)
    {
        // retrieve the cell specific PVT table index from the deck
        // and using the grid...
        extractPvtTableIndex(cellPvtRegionIdx_, deck, number_of_cells, global_cell);

        if(init_rock){
            rock_.init(eclState, number_of_cells, global_cell, cart_dims);
        }

        const int pvt_samples = param.getDefault("pvt_tab_size", -1);
        pvt_.init(deck, eclState, pvt_samples);

        // Unfortunate lack of pointer smartness here...
        std::string threephase_model = param.getDefault<std::string>("threephase_model", "gwseg");
        if (deck->hasKeyword("ENDSCALE") && threephase_model != "gwseg") {
            OPM_THROW(std::runtime_error, "Sorry, end point scaling currently available for the 'gwseg' model only.");
        }

        SaturationPropsFromDeck* ptr
            = new SaturationPropsFromDeck();
        ptr->init(phaseUsageFromDeck(deck), materialLawManager);
        satprops_.reset(ptr);

        if (pvt_.numPhases() != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }
}

