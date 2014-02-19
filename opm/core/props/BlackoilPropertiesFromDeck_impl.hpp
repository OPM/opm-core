namespace Opm
{
    template<class T>
    void BlackoilPropertiesFromDeck::init(const EclipseGridParser& deck,
                                          int number_of_cells,
                                          const int* global_cell,
                                          const int* cart_dims,
                                          T begin_cell_centroids,
                                          int dimension,
                                          bool init_rock)
    {
        if (init_rock){
           rock_.init(deck, number_of_cells, global_cell, cart_dims);
        }
        pvt_.init(deck, 200);
        SaturationPropsFromDeck<SatFuncSimpleUniform>* ptr
            = new SaturationPropsFromDeck<SatFuncSimpleUniform>();
        satprops_.reset(ptr);
        ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimension, 200);

        if (pvt_.numPhases() != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }

    template<class T>
    void BlackoilPropertiesFromDeck::init(Opm::DeckConstPtr newParserDeck,
                                          int number_of_cells,
                                          const int* global_cell,
                                          const int* cart_dims,
                                          T begin_cell_centroids,
                                          int dimension,
                                          bool init_rock)
    {
        if (init_rock){
           rock_.init(newParserDeck, number_of_cells, global_cell, cart_dims);
        }
        pvt_.init(newParserDeck, /*numSamples=*/200);
        SaturationPropsFromDeck<SatFuncSimpleUniform>* ptr
            = new SaturationPropsFromDeck<SatFuncSimpleUniform>();
        satprops_.reset(ptr);
        ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids, dimension,
                  /*numSamples=*/200);

        if (pvt_.numPhases() != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }

    template<class T>
    void BlackoilPropertiesFromDeck::init(const EclipseGridParser& deck,
                                          int number_of_cells,
                                          const int* global_cell,
                                          const int* cart_dims,
                                          T begin_cell_centroids,
                                          int dimension,
                                          const parameter::ParameterGroup& param,
                                          bool init_rock)
    {
        if(init_rock){
            rock_.init(deck, number_of_cells, global_cell, cart_dims);
        }

        const int pvt_samples = param.getDefault("pvt_tab_size", 200);
        pvt_.init(deck, pvt_samples);

        // Unfortunate lack of pointer smartness here...
        const int sat_samples = param.getDefault("sat_tab_size", 200);
        std::string threephase_model = param.getDefault<std::string>("threephase_model", "simple");
        if (deck.hasField("ENDSCALE") && threephase_model != "simple") {
            OPM_THROW(std::runtime_error, "Sorry, end point scaling currently available for the 'simple' model only.");
        }
        if (sat_samples > 1) {
            if (threephase_model == "stone2") {
                SaturationPropsFromDeck<SatFuncStone2Uniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncStone2Uniform>();
                satprops_.reset(ptr);
                ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimension,
                          sat_samples);
            } else if (threephase_model == "simple") {
                SaturationPropsFromDeck<SatFuncSimpleUniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncSimpleUniform>();
                satprops_.reset(ptr);
                ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimension,
                          sat_samples);
            } else if (threephase_model == "gwseg") {
                SaturationPropsFromDeck<SatFuncGwsegUniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncGwsegUniform>();
                satprops_.reset(ptr);
                ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimension,
                          sat_samples);
            } else {
                OPM_THROW(std::runtime_error, "Unknown threephase_model: " << threephase_model);
            }
        } else {
            if (threephase_model == "stone2") {
                SaturationPropsFromDeck<SatFuncStone2Nonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncStone2Nonuniform>();
                satprops_.reset(ptr);
                ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimension,
                          sat_samples);
            } else if (threephase_model == "simple") {
                SaturationPropsFromDeck<SatFuncSimpleNonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncSimpleNonuniform>();
                satprops_.reset(ptr);
                ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimension,
                          sat_samples);
            } else if (threephase_model == "gwseg") {
                SaturationPropsFromDeck<SatFuncGwsegNonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncGwsegNonuniform>();
                satprops_.reset(ptr);
                ptr->init(deck, number_of_cells, global_cell, begin_cell_centroids, dimension,
                          sat_samples);
            } else {
                OPM_THROW(std::runtime_error, "Unknown threephase_model: " << threephase_model);
            }
        }

        if (pvt_.numPhases() != satprops_->numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_->numPhases() << ").");
        }
    }

    template<class T>
    void BlackoilPropertiesFromDeck::init(Opm::DeckConstPtr newParserDeck,
                                          int number_of_cells,
                                          const int* global_cell,
                                          const int* cart_dims,
                                          T begin_cell_centroids,
                                          int dimension,
                                          const parameter::ParameterGroup& param,
                                          bool init_rock)
    {

        if(init_rock){
            rock_.init(newParserDeck, number_of_cells, global_cell, cart_dims);
        }

        const int pvt_samples = param.getDefault("pvt_tab_size", 200);
        pvt_.init(newParserDeck, pvt_samples);

        // Unfortunate lack of pointer smartness here...
        const int sat_samples = param.getDefault("sat_tab_size", 200);
        std::string threephase_model = param.getDefault<std::string>("threephase_model", "simple");
        if (newParserDeck->hasKeyword("ENDSCALE") && threephase_model != "simple") {
            OPM_THROW(std::runtime_error, "Sorry, end point scaling currently available for the 'simple' model only.");
        }
        if (sat_samples > 1) {
            if (threephase_model == "stone2") {
                SaturationPropsFromDeck<SatFuncStone2Uniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncStone2Uniform>();
                satprops_.reset(ptr);
                ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "simple") {
                SaturationPropsFromDeck<SatFuncSimpleUniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncSimpleUniform>();
                satprops_.reset(ptr);
                ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "gwseg") {
                SaturationPropsFromDeck<SatFuncGwsegUniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncGwsegUniform>();
                satprops_.reset(ptr);
                ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else {
                OPM_THROW(std::runtime_error, "Unknown threephase_model: " << threephase_model);
            }
        } else {
            if (threephase_model == "stone2") {
                SaturationPropsFromDeck<SatFuncStone2Nonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncStone2Nonuniform>();
                satprops_.reset(ptr);
                ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "simple") {
                SaturationPropsFromDeck<SatFuncSimpleNonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncSimpleNonuniform>();
                satprops_.reset(ptr);
                ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids,
                          dimension, sat_samples);
            } else if (threephase_model == "gwseg") {
                SaturationPropsFromDeck<SatFuncGwsegNonuniform>* ptr
                    = new SaturationPropsFromDeck<SatFuncGwsegNonuniform>();
                satprops_.reset(ptr);
                ptr->init(newParserDeck, number_of_cells, global_cell, begin_cell_centroids,
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

