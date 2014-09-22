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
        init(deck, eclState, number_of_cells, global_cell, cart_dims, begin_cell_centroids, dimension,
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
        init(deck,
             eclState,
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
        SaturationPropsFromDeck<SatFuncSimpleNonuniform>* ptr
            = new SaturationPropsFromDeck<SatFuncSimpleNonuniform>();
        satprops_.reset(ptr);
        ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids, dimension,
                  /*numSamples=*/0);

        if (pvt_.numPhases() != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }

    template<class CentroidIterator>
    inline void BlackoilPropertiesFromDeck::init(Opm::DeckConstPtr deck,
                                                 Opm::EclipseStateConstPtr eclState,
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
        const int sat_samples = param.getDefault("sat_tab_size", -1);
        std::string threephase_model = param.getDefault<std::string>("threephase_model", "gwseg");
        if (deck->hasKeyword("ENDSCALE") && threephase_model != "gwseg") {
            OPM_THROW(std::runtime_error, "Sorry, end point scaling currently available for the 'gwseg' model only.");
        }
        if (sat_samples > 1) {
            if (threephase_model == "stone2") {
                SaturationPropsFromDeck<SatFuncStone2Uniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncStone2Uniform>();
                satprops_.reset(ptr);
                ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "simple") {
                SaturationPropsFromDeck<SatFuncSimpleUniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncSimpleUniform>();
                satprops_.reset(ptr);
                ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "gwseg") {
                SaturationPropsFromDeck<SatFuncGwsegUniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncGwsegUniform>();
                satprops_.reset(ptr);
                ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else {
                OPM_THROW(std::runtime_error, "Unknown threephase_model: " << threephase_model);
            }
        } else {
            if (threephase_model == "stone2") {
                SaturationPropsFromDeck<SatFuncStone2Nonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncStone2Nonuniform>();
                satprops_.reset(ptr);
                ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "simple") {
                SaturationPropsFromDeck<SatFuncSimpleNonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncSimpleNonuniform>();
                satprops_.reset(ptr);
                ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "gwseg") {
                SaturationPropsFromDeck<SatFuncGwsegNonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncGwsegNonuniform>();
                satprops_.reset(ptr);
                ptr->init(deck, eclState, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else {
                OPM_THROW(std::runtime_error, "Unknown threephase_model: " << threephase_model);
            }
        }

        if (pvt_.numPhases() != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }
}

