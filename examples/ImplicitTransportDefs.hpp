
#ifndef OPENRS_IMPLICITTRANSPORTDEFS_HEADER
#define OPENRS_IMPLICITTRANSPORTDEFS_HEADER

#include <dune/porsol/opmtransport/src/NormSupport.hpp>
#include <dune/porsol/opmtransport/src/ImplicitTransport.hpp>
template <class Ostream, class Collection>
Ostream&
operator<<(Ostream& os, const Collection& c)
{
    typedef typename Collection::value_type VT;

    os << "[ ";
    std::copy(c.begin(), c.end(), ::std::ostream_iterator<VT>(os, " "));
    os << "]";

    return os;
}
template <class Vector>
class MaxNorm {
	public:
	    static double
	    norm(const Vector& v) {
		return  Opm::ImplicitTransportDefault::AccumulationNorm <Vector,  Opm::ImplicitTransportDefault::MaxAbs>::norm(v);
	    }
};

template <int np = 2>
class ReservoirState { 
public:
    ReservoirState(const grid_t* g)
        : press_ (g->number_of_cells),
          fpress_(g->number_of_faces),
          flux_  (g->number_of_faces),
          sat_   (np * g->number_of_cells)
    {}

    ::std::vector<double>& pressure    () { return press_ ; }
    ::std::vector<double>& facepressure() { return fpress_; }
    ::std::vector<double>& faceflux    () { return flux_  ; }
    ::std::vector<double>& saturation  () { return sat_   ; }

    const ::std::vector<double>& faceflux    () const { return flux_; }
    const ::std::vector<double>& saturation  () const { return sat_ ; }
 
private:
    ::std::vector<double> press_ ;
    ::std::vector<double> fpress_;
    ::std::vector<double> flux_  ;
    ::std::vector<double> sat_   ;
};
class Rock {
public:
    Rock(::std::size_t nc, ::std::size_t dim)
        : dim_ (dim           ),
          perm_(nc * dim * dim),
          poro_(nc            ) {}

    const ::std::vector<double>& perm() const { return perm_; }
    const ::std::vector<double>& poro() const { return poro_; }

    void
    perm_homogeneous(double k) {
        setVector(0.0, perm_);

        const ::std::size_t d2 = dim_ * dim_;

        for (::std::size_t c = 0, nc = poro_.size(); c < nc; ++c) {
            for (::std::size_t i = 0; i < dim_; ++i) {
                perm_[c*d2 + i*(dim_ + 1)] = k;
            }
        }
    }

    void
    poro_homogeneous(double phi) {
        setVector(phi, poro_);
    }

private:
    void
    setVector(double x, ::std::vector<double>& v) {
        ::std::fill(v.begin(), v.end(), x);
    }

    ::std::size_t         dim_ ;
    ::std::vector<double> perm_;
    ::std::vector<double> poro_;
};

#endif // /OPENRS_IMPLICITTRANSPORTDEFS_HEADER
