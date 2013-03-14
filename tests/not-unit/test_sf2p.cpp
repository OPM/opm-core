#include <algorithm>
#include <tr1/array>
#include <iostream>
#include <iterator>

#include <opm/core/props/SimpleFluid2p.hpp>

template <class Ostream, typename T, size_t n>
Ostream&
operator<<(Ostream& os, const std::tr1::array<T, n>& a)
{
    os << "[ ";
    std::copy(a.begin(), a.end(), std::ostream_iterator<T>(os, " "));
    os << "]";

    return os;
}

template <int n>
void
test_simplefluid2p()
{
    using std::tr1::array;
    array<double, 2> mu  = { {1.0, 1.0} };
    array<double, 2> rho = { {0.0, 0.0} };

    Opm::SimpleFluid2p<n> fluid(mu, rho);

    std::cerr << "\\rho = [ " << fluid.density(0)
              << ", "         << fluid.density(1) << " ]\n";

    array<double, 2> sl = {{ 1.0, 0.0 }};

    array<double, 2>     mob;
    array<double, 2 * 2> dmob;

    fluid.mobility(0, sl, mob, dmob);
    std::cerr << "s = " << sl << ", m = " << mob << ", dm = " << dmob << '\n';

    array<double, 2> sm = {{ 0.5, 0.5 }};
    fluid.mobility(0, sm, mob, dmob);

    std::cerr << "s = " << sm << ", m = " << mob << ", dm = " << dmob << '\n';

    array<double, 2> sr = {{ 0.0, 1.0 }};
    fluid.mobility(0, sr, mob, dmob);

    std::cerr << "s = " << sr << ", m = " << mob << ", dm = " << dmob << '\n';
}

int main()
{
    std::cerr << "n = 1\n";
    test_simplefluid2p<1>();

    std::cerr << "n = 2\n";
    test_simplefluid2p<2> ();

    std::cerr << "n = 3\n";
    test_simplefluid2p<3>();

    return 0;
}
