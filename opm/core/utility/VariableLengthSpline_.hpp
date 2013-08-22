// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2013 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::VariableLengthSpline_
 */
#ifndef OPM_VARIABLE_LENGTH_SPLINE_HH
#define OPM_VARIABLE_LENGTH_SPLINE_HH

#include "SplineCommon_.hpp"
#include "TridiagonalMatrix.hpp"

#include <array>
#include <vector>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Opm {

/*!
 * \brief The common code for all 3rd order polynomial splines with
 *        where the number of sampling points only known at run-time.
 */
template<class Scalar>
class VariableLengthSpline_
    : public SplineCommon_<Scalar,
                           VariableLengthSpline_<Scalar> >
{
    struct ComparatorX_
    {
        ComparatorX_(const std::vector<Scalar> &x)
            : x_(x)
        {};

        bool operator ()(int idxA, int idxB) const
        { return x_.at(idxA) < x_.at(idxB); }

        const std::vector<Scalar> &x_;
    };

    friend class SplineCommon_<Scalar, VariableLengthSpline_<Scalar> >;

    typedef std::vector<Scalar> Vector;
    typedef Opm::TridiagonalMatrix<Scalar> Matrix;

public:
    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples() const
    { return xPos_.size(); }


    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////
    // Full splines                      //
    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using C-style arrays.
     *
     * This method uses separate array-like objects for the values of
     * the X and Y coordinates. In this context 'array-like' means
     * that an access to the members is provided via the []
     * operator. (e.g. C arrays, std::vector, std::array, etc.)  Each
     * array must be of size 'nSamples' at least. Also, the number of
     * sampling points must be larger than 1.
     */
    template <class ScalarArrayX, class ScalarArrayY>
    void setXYArrays(int nSamples,
                     const ScalarArrayX &x,
                     const ScalarArrayY &y,
                     Scalar m0, Scalar m1,
                     bool sortInputs = false)
    {
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = x[i];
            yPos_[i] = y[i];
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using STL-compatible containers.
     *
     * This method uses separate STL-compatible containers for the
     * values of the sampling points' X and Y
     * coordinates. "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class ScalarContainerX, class ScalarContainerY>
    void setXYContainers(const ScalarContainerX &x,
                         const ScalarContainerY &y,
                         Scalar m0, Scalar m1,
                         bool sortInputs = false)
    {
        assert(x.size() == y.size());
        assert(x.size() > 1);

        setNumSamples_(x.size());

        std::copy(x.begin(), x.end(), xPos_.begin());
        std::copy(y.begin(), y.end(), yPos_.begin());

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using a C-style array.
     *
     * This method uses a single array of sampling points, which are
     * seen as an array-like object which provides access to the X and
     * Y coordinates.  In this context 'array-like' means that an
     * access to the members is provided via the [] operator. (e.g. C
     * arrays, std::vector, std::array, etc.)  The array containing
     * the sampling points must be of size 'nSamples' at least. Also,
     * the number of sampling points must be larger than 1.
     */
    template <class PointArray>
    void setArrayOfPoints(int nSamples,
                          const PointArray &points,
                          Scalar m0,
                          Scalar m1,
                          bool sortInputs = false)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = points[i][0];
            yPos_[i] = points[i][1];
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using a STL-compatible container of
     *        array-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be array-like objects storing the
     * X and Y coordinates.  "STL-compatible" means that the container
     * provides access to iterators using the begin(), end() methods
     * and also provides a size() method. Also, the number of entries
     * in the X and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfPoints(const XYContainer &points,
                              Scalar m0,
                              Scalar m1,
                              bool sortInputs = false)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(points.size() > 1);

        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++i, ++it) {
            xPos_[i] = (*it)[0];
            yPos_[i] = (*it)[1];
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        // make a full spline
        makeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of a
     *        full spline using a STL-compatible container of
     *        tuple-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be tuple-like objects storing the
     * X and Y coordinates.  "tuple-like" means that the objects
     * provide access to the x values via std::get<0>(obj) and to the
     * y value via std::get<1>(obj) (e.g. std::tuple or
     * std::pair). "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfTuples(const XYContainer &points,
                              Scalar m0,
                              Scalar m1,
                              bool sortInputs = false)
    {
        // resize internal arrays
        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++i, ++it) {
            xPos_[i] = std::get<0>(*it);
            yPos_[i] = std::get<1>(*it);
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        // make a full spline
        makeFullSpline_(m0, m1);
    }

    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////
    // Natural/Periodic splines          //
    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////
    /*!
     * \brief Set the sampling points natural spline using C-style arrays.
     *
     * This method uses separate array-like objects for the values of
     * the X and Y coordinates. In this context 'array-like' means
     * that an access to the members is provided via the []
     * operator. (e.g. C arrays, std::vector, std::array, etc.)  Each
     * array must be of size 'nSamples' at least. Also, the number of
     * sampling points must be larger than 1.
     */
    template <class ScalarArrayX, class ScalarArrayY>
    void setXYArrays(int nSamples,
                     const ScalarArrayX &x,
                     const ScalarArrayY &y,
                     SplineType splineType = NaturalSpline,
                     bool sortInputs = false)
    {
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = x[i];
            yPos_[i] = y[i];
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        if (splineType == PeriodicSpline)
            makePeriodicSpline_();
        else if (splineType == NaturalSpline)
            makeNaturalSpline_();
        else if (splineType == MonotonicSpline)
            this->makeMonotonicSpline_(slopeVec_);
        else
            OPM_THROW(std::runtime_error, "Spline type " << splineType << " not supported at this place");
    }

    /*!
     * \brief Set the sampling points of a natural spline using
     *        STL-compatible containers.
     *
     * This method uses separate STL-compatible containers for the
     * values of the sampling points' X and Y
     * coordinates. "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class ScalarContainerX, class ScalarContainerY>
    void setXYContainers(const ScalarContainerX &x,
                         const ScalarContainerY &y,
                         SplineType splineType = NaturalSpline,
                         bool sortInputs = false)
    {
        assert(x.size() == y.size());
        assert(x.size() > 1);

        setNumSamples_(x.size());
        std::copy(x.begin(), x.end(), xPos_.begin());
        std::copy(y.begin(), y.end(), yPos_.begin());

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        if (splineType == PeriodicSpline)
            makePeriodicSpline_();
        else if (splineType == NaturalSpline)
            makeNaturalSpline_();
        else if (splineType == MonotonicSpline)
            this->makeMonotonicSpline_(slopeVec_);
        else
            OPM_THROW(std::runtime_error, "Spline type " << splineType << " not supported at this place");
    }

    /*!
     * \brief Set the sampling points of a natural spline using a
     *        C-style array.
     *
     * This method uses a single array of sampling points, which are
     * seen as an array-like object which provides access to the X and
     * Y coordinates.  In this context 'array-like' means that an
     * access to the members is provided via the [] operator. (e.g. C
     * arrays, std::vector, std::array, etc.)  The array containing
     * the sampling points must be of size 'nSamples' at least. Also,
     * the number of sampling points must be larger than 1.
     */
    template <class PointArray>
    void setArrayOfPoints(int nSamples,
                          const PointArray &points,
                          SplineType splineType = NaturalSpline,
                          bool sortInputs = false)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(nSamples > 1);

        setNumSamples_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            xPos_[i] = points[i][0];
            yPos_[i] = points[i][1];
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        if (splineType == PeriodicSpline)
            makePeriodicSpline_();
        else if (splineType == NaturalSpline)
            makeNaturalSpline_();
        else if (splineType == MonotonicSpline)
            this->makeMonotonicSpline_(slopeVec_);
        else
            OPM_THROW(std::runtime_error, "Spline type " << splineType << " not supported at this place");
    }

    /*!
     * \brief Set the sampling points of a natural spline using a
     *        STL-compatible container of array-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be array-like objects storing the
     * X and Y coordinates.  "STL-compatible" means that the container
     * provides access to iterators using the begin(), end() methods
     * and also provides a size() method. Also, the number of entries
     * in the X and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfPoints(const XYContainer &points,
                              SplineType splineType = NaturalSpline,
                              bool sortInputs = false)
    {
        // a spline with no or just one sampling points? what an
        // incredible bad idea!
        assert(points.size() > 1);

        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++ i, ++it) {
            xPos_[i] = (*it)[0];
            yPos_[i] = (*it)[1];
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        if (splineType == PeriodicSpline)
            makePeriodicSpline_();
        else if (splineType == NaturalSpline)
            makeNaturalSpline_();
        else if (splineType == MonotonicSpline)
            this->makeMonotonicSpline_(slopeVec_);
        else
            OPM_THROW(std::runtime_error, "Spline type " << splineType << " not supported at this place");
    }

    /*!
     * \brief Set the sampling points of a natural spline using a
     *        STL-compatible container of tuple-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be tuple-like objects storing the
     * X and Y coordinates.  "tuple-like" means that the objects
     * provide access to the x values via std::get<0>(obj) and to the
     * y value via std::get<1>(obj) (e.g. std::tuple or
     * std::pair). "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfTuples(const XYContainer &points,
                              SplineType splineType = NaturalSpline,
                              bool sortInputs = false)
    {
        // resize internal arrays
        setNumSamples_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++i, ++it) {
            xPos_[i] = std::get<0>(*it);
            yPos_[i] = std::get<1>(*it);
        }

        if (sortInputs)
            sortInput_();
        else if (xPos_[0] > xPos_[numSamples() - 1])
            reverseSamplingPoints_();

        if (splineType == PeriodicSpline)
            makePeriodicSpline_();
        else if (splineType == NaturalSpline)
            makeNaturalSpline_();
        else if (splineType == MonotonicSpline)
            this->makeMonotonicSpline_(slopeVec_);
        else
            OPM_THROW(std::runtime_error, "Spline type " << splineType << " not supported at this place");
    }

protected:
    /*!
     * \brief Sort the sample points in ascending order of their x value.
     */
    void sortInput_()
    {
        int n = numSamples();

        // create a vector containing 0...n-1
        std::vector<int> idxVector(n);
        for (int i = 0; i < n; ++i)
            idxVector[i] = i;

        // sort the indices according to the x values of the sample
        // points
        ComparatorX_ cmp(xPos_);
        std::sort(idxVector.begin(), idxVector.end(), cmp);

        // reorder the sample points
        std::vector<Scalar> tmpX(n), tmpY(n);
        for (int i = 0; i < idxVector.size(); ++ i) {
            tmpX[i] = xPos_[idxVector[i]];
            tmpY[i] = yPos_[idxVector[i]];
        }
        xPos_ = tmpX;
        yPos_ = tmpY;
    }

    /*!
     * \brief Reverse order of the elements in the arrays which
     *        contain the sampling points.
     */
    void reverseSamplingPoints_()
    {
        // reverse the arrays
        int n = numSamples();
        for (int i = 0; i <= (n - 1)/2; ++i) {
            std::swap(xPos_[i], xPos_[n - i - 1]);
            std::swap(yPos_[i], yPos_[n - i - 1]);
        }
    }

    /*!
     * \brief Resizes the internal vectors to store the sample points.
     */
    void setNumSamples_(int nSamples)
    {
        xPos_.resize(nSamples);
        yPos_.resize(nSamples);
        slopeVec_.resize(nSamples);
    }

    /*!
     * \brief Create a natural spline from the already set sampling points.
     *
     * This creates a temporary matrix and right hand side vector.
     */
    void makeFullSpline_(Scalar m0, Scalar m1)
    {
        Matrix M(numSamples());
        std::vector<Scalar> d(numSamples());
        std::vector<Scalar> moments(numSamples());

        // create linear system of equations
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments (-> second derivatives)
        M.solve(moments, d);

        // convert the moments to slopes at the sample points
        this->setSlopesFromMoments_(slopeVec_, moments);
    }

    /*!
     * \brief Create a natural spline from the already set sampling points.
     *
     * This creates a temporary matrix and right hand side vector.
     */
    void makeNaturalSpline_()
    {
        Dune::DynamicMatrix<Scalar> M(numSamples(), numSamples());
        Dune::DynamicVector<Scalar> d(numSamples());
        Dune::DynamicVector<Scalar> moments(numSamples());

        // create linear system of equations
        this->makeNaturalSystem_(M, d);

        // solve for the moments (-> second derivatives)
        M.solve(moments, d);

        // convert the moments to slopes at the sample points
        this->setSlopesFromMoments_(slopeVec_, moments);
    }

    /*!
     * \brief Create a periodic spline from the already set sampling points.
     *
     * This creates a temporary matrix and right hand side vector.
     */
    void makePeriodicSpline_()
    {
        Matrix M(numSamples() - 1);
        Vector d(numSamples() - 1);
        Vector moments(numSamples() - 1);

        // create linear system of equations. This is a bit hacky,
        // because it assumes that std::vector internally stores its
        // data as a big C-style array, but it saves us from yet
        // another copy operation
        this->makePeriodicSystem_(M, d);

        // solve for the moments (-> second derivatives)
        M.solve(moments, d);

        moments.resize(numSamples());
        for (int i = numSamples() - 2; i >= 0; --i)
            moments[i+1] = moments[i];
        moments[0] = moments[numSamples() - 1];

        // convert the moments to slopes at the sample points
        this->setSlopesFromMoments_(slopeVec_, moments);
    }

    /*!
     * \brief Returns the x coordinate of the i-th sampling point.
     */
    Scalar x_(int i) const
    { return xPos_[i]; }

    /*!
     * \brief Returns the y coordinate of the i-th sampling point.
     */
    Scalar y_(int i) const
    { return yPos_[i]; }

    /*!
     * \brief Returns slope (i.e. first derivative) of the spline at
     *        the i-th sampling point.
     */
    Scalar slope_(int i) const
    { return slopeVec_[i]; }

    Vector xPos_;
    Vector yPos_;
    Vector slopeVec_;
};
}

#endif
