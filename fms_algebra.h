// fms_algebra.h
#pragma once
#include <iterator>

namespace fms {

    template<class A>
    struct atom_iterator {
        bool operator==(const atom_iterator& ai) const
        {
            return !op_bool(); // this knows end()
        }
        bool operator!=(const atom_iterator& ai) const
        {
            return !operator==(ai);
        }
        // can be derferenced?
        operator bool() const { return op_bool(); }
        const A& operator*() const { return op_star(); }
        atom_iterator& operator++() { return op_plus(); }
    private:
        // override in derived
        virtual bool op_bool() const = 0;
        virtual const A& op_star() const = 0;
        virtual atom_iterator& op_plus() = 0;
    };

    // An algebra as a partition into atoms.
    template<class A>
    struct algebra {
        // iterate over atoms
        atom_iterator<A>& begin()
        {
            return _begin();
        }
        atom_iterator<A>& end()
        {
            return _end();
        }
        // atoms of type B contained in atom
        template<class B>
        algebra<B> partition(const A& atom)
        {
            return algebra<B>{};
        }
    private:
        virtual atom_iterator<A>& _begin() = 0;
        virtual atom_iterator<A>& _end() = 0;
        // template cannot be virtual
    };

} // fms
