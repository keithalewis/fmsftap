// main.cpp - test driver
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "fms_algebra.h"

/*
// value at time t on atom A for instrument I, and discount D
value(t,A,I,D)
{
    v = 0;
    [u, c] = I.next(t,A) // next cash flow time and cash flow for instrument given A

    for (B in atoms(A, u)) // atoms at time u contained in A
    {
        v += (c(B) + value(u, B, I, D))*D(B)
    }

    return v/D(A);
}
We want to write:
class<Algebra A, Instrument I>
I::cash_flow_type value(I::time_type t, A::atom_type a, ...)
{
    v = 0;
    auto [u, c] = ++i; // next time and cash flow of instrument
     
    for (b : atoms(a, u)) {
        v += (c(b) + value(u, b, ...)*D(b);
    }

    return v/D(a);
}
*/

constexpr size_t choose(size_t n, size_t k)
{
    if (2*k > n)
        return choose(n, n - k);
    if (n == 0 || k == 0)
        return 1;

    return (n*choose(n-1,k-1))/k;
}
template<size_t N> 
struct binomial_atom {
    size_t k; // level
    binomial_atom(size_t k = 0)
        : k(k)
    { }
    bool operator==(const binomial_atom& a) const
    {
        return k == a.k;
    }
    bool operator!=(const binomial_atom& a) const
    {
        return k != a.k;
    }
};

template<size_t N>
inline double prob(const binomial_atom<N>& atom)
{
    return choose(N,atom.k)*ldexp(1, -N);
}

template<size_t N>
class binomial_atom_iterator : public fms::atom_iterator<binomial_atom<N>> {
    binomial_atom<N> atom;
public:
    binomial_atom_iterator(size_t k = 0)
        : atom{k}
    { }
    bool operator==(const binomial_atom_iterator& a) const
    {
        return atom.k == a.atom.k;
    }
    bool operator!=(const binomial_atom_iterator& a) const
    {
        return !operator==(a);
    }
    /* !!!only return type can be covariant
    const bool op_equal(const binomial_atom_iterator& a) const override
    {
        return atom == a.atom;
    }
    */
    bool op_bool() const override
    {
        return atom.k <= N;
    }
    const binomial_atom<N>& op_star() const override
    {
        return atom;
    }
    binomial_atom_iterator& op_plus() override
    {
        ++atom.k;

        return *this;
    }
};

template<size_t N>
struct binomial_algebra : fms::algebra<binomial_atom<N>> {
    binomial_atom_iterator<N> b, e;
    binomial_algebra()
        : b(0), e(N+1)
    { }
    binomial_atom_iterator<N>& _begin() override
    {
        return b;
    }
    binomial_atom_iterator<N>& _end() override
    {
        return e;
    }
};

using namespace fms;

void test_binomial_atom_iterator()
{
    binomial_atom_iterator<3> a3;
    binomial_atom_iterator<3> a{a3};
    assert (a == a3);
    a3 = a;
    assert (a == a3);
    assert (a);
    assert (*a == *a3);
    assert (prob(*a) == 1./8);
    ++a; ++a3;
    assert (a);
    assert (*a == *a3);
    assert (prob(*a) == 3./8);
    ++a; assert(a);
    assert (prob(*a) == 3./8);
    ++a; assert(a);
    assert (prob(*a) == 1./8);
    ++a; assert(!a);
}
void test_binomial_algebra()
{
    binomial_algebra<3> a3;
    binomial_algebra<3> a{a3};
    a3 = a;
    auto& i = a.begin();
    assert (i != a.end());
    assert (prob(*i) == 1./8);
    ++i;
    assert (i != a.end());
    assert (prob(*i) == 3./8);
    ++i;
    assert (i != a.end());
    assert (prob(*i) == 3./8);
    ++i;
    assert (i != a.end());
    assert (prob(*i) == 1./8);
    ++i;
    assert (i == a.end());

    double one = 0;
    for (auto& i : a) {
        one += prob(*i);
    }
    assert (one == 1);
}
int main() 
{
    test_binomial_atom_iterator();
    test_binomial_algebra();

    return 0;
}