# Tinned

Tinned contains a set of nonnumerical routines for computational chemistry, in
which symbolic routines are built on top of
[SymEngine library](https://github.com/symengine/symengine). Theoretical
background can be found in reference [[1]](#1).

Currently, the following project uses Tinned library:

* [SymResponse](https://github.com/bingao/symresponse), a unified framework for
  response theory.

## License

Tinned is licensed under Mozilla Public License Version 2.0, see the
[LICENSE](LICENSE) for more information.

## Installation

Both SymEngine and Tinned require CMake and C++ compiler which supports C++11
standard.

Clone and build [forked SymEngine library](https://github.com/bingao/symengine)
first, which has implmented derivatives for different matrix expressions.

Then clone Tinned library and build it by setting `SymEngine_DIR` to the
SymEngine installation or build directory.

## Tinned APIs

Tinned currently provides C++ interface. Classes in Tinned that can be useful
for computational chemistry include:

* Class [`Perturbation`](include/Tinned/Perturbation.hpp), perturbations.
* Class [`PertTuple`](include/Tinned/PertTuple.hpp), `std::multiset` for
  perturbation tuples[[1]](#1).
* Class [`PertDependency`](include/Tinned/PertDependency.hpp), `std::set` for
  perturbations that an operator depends on and their maximum orders that can
  be differentiated.
* Class [`PerturbedParameter`](include/Tinned/PerturbedParameter.hpp),
  (perturbed) response parameter.
* Class [`OneElecDensity`](include/Tinned/OneElecDensity.hpp), one-electron
  spin-orbital density matrix derived from the abstract electronic state
  class [`ElectronicState`](include/Tinned/ElectronicState.hpp).
* Class [`OneElecOperator`](include/Tinned/OneElecOperator.hpp), one-electron
  like operators.
* Class [`TwoElecEnergy`](include/Tinned/TwoElecEnergy.hpp), two-electron like
  energies.
* Class [`TwoElecOperator`](include/Tinned/TwoElecOperator.hpp), two-electron
  like operators.
* Class [`ExchCorrEnergy`](include/Tinned/ExchCorrEnergy.hpp),
  exchange-correlation energy like functionals.
* Class [`ExchCorrPotential`](include/Tinned/ExchCorrPotential.hpp),
  exchange-correlation potential like operators.
* Class [`NonElecFunction`](include/Tinned/NonElecFunction.hpp), non-electron
  like functions.
* Class [`TemporumOperator`](include/Tinned/TemporumOperator.hpp), for
  operators $\langle\mathrm{i}\frac{\partial}{\partial t}\vert$ and
  $\vert\mathrm{i}\frac{\partial}{\partial t}\rangle$.
* Class [`TemporumOverlap`](include/Tinned/TemporumOverlap.hpp), for the operator
  $-\frac{1}{2}(\langle\mathrm{i}\frac{\partial}{\partial t}\chi_{\kappa}\vert\chi_{\lambda}\rangle
  +\langle\chi_{\kappa}\vert\mathrm{i}\frac{\partial}{\partial t}\chi_{\lambda}\rangle)$.
* Class [`AdjointMap`](include/Tinned/AdjointMap.hpp) represents an adjoint map
  $\prod_{j}\left(\text{ad}_{\tilde{\mathbf{X}}_{j}}\right)(\tilde{\mathbf{Y}})$
  where all $\tilde{\mathbf{X}}_{j}$'s and all their derivatives are commutative.
* Class [`ClusterConjHamiltonian`](include/Tinned/ClusterConjHamiltonian.hpp)
  is used for coupled-cluster theory and represents an exponential map in the
  form of $\text{e}^{\text{ad}_{-\hat{T}}}(\hat{H}^{b_{P}})$ or
  $\text{e}^{\text{ad}_{-\hat{T}}}\left(\prod_{j=1}^{j_{\max}}
    \left(\text{ad}_{\hat{T}^{b_{Q_{j}}}}\right)(\hat{H}^{b_{P}})\right)$ with
  $1\le j_{\max}\le3$.

**It is recommended that one uses helper functions in Tinned library to create
objects of any above classes.** These helper functions are defined and
implemented in the corresponding [header files](include/Tinned), and are named
as `make_...()`. Users can refer to comments in these header files for the use
of helper functions.

To facilitate the development of response theory, the following functions are
provided by Tinned:

* Function [`remove_zeros(x)`](include/Tinned/ZerosRemover.hpp) removes zero
  quantities from `x`.
* Function [`remove_if(x, symbols)`](include/Tinned/RemoveVisitor.hpp) removes
  given `symbols` from `x`.
* Function [`keep_if(x, symbols, remove_zero_quantities)`](include/Tinned/KeepVisitor.hpp)
  keeps given `symbols` in `x` while removing others, parameter
  `remove_zero_quantities` indicates if zero quantities will be removed from the
  output.
* Function [`replace(x, subs_dict, cache)`](include/Tinned/ReplaceVisitor.hpp)
  replaces classes defined in Tinned library in addition to those in
  `SymEngine::msubs()`;
* Function [`find_all(x, symbol)`](include/Tinned/FindAllVisitor.hpp)
  finds a given `symbol` and all its differentiated ones in `x`;
* Function [`eliminate(x, parameter, perturbations, min_order)`](include/Tinned/EliminationVisitor.hpp)
  eliminates a given response `parameter`'s derivatives from `x`. Maximum order
  of derivatives to be eliminated is the length of `perturbations`, and minimum
  order is specified by `min_order`. For wave function parameters, it should be
  greater than the floor function of the half length of `perturbations`, and for
  multipliers, it should be greater than or equal to the ceiling function of the
  half length of `perturbations` according to J. Chem. Phys. 129, 214103 (2008).
* Function [`clean_temporum(x)`](include/Tinned/TemporumCleaner.hpp) cleans
  `TemporumOperator` objects in `x`.
* Function [`latexify(x)`](include/Tinned/LaTeXifyVisitor.hpp) latexifies an
  expression `x`.
* Function [`stringify(x)`](include/Tinned/StringifyVisitor.hpp) stringifies an
  expression `x`.
* Function [`differentiate(expr, perturbations)`](include/Tinned/Utilities.hpp)
  can be used to do high-order differentiation, and to remove zero quantities.
* Function template [`replace_all<T>(x, subs_dict)`](include/Tinned/Utilities.hpp)
  replaces Tinned objects and their derivatives with SymEngine `Basic` symbols
  and corresponding derivatives. Template parameter `T` is the type of those
  Tinned objects.

**One should note that:**

* Functions `remove_zeros`, `remove_if`, `keep_if` and `eliminate` may
  return a null pointer `SymEngine::RCP<const SymEngine::Basic>()` when no
  symbols left after the action. Users can call `is_null()` to check the output
  of these functions.
* Functions `clean_temporum` and `differentiate` may return a zero output of
  type either `Tinned::ZeroOperator` or `SymEngine::zero`. Users can check if
  the output is a zero quantity by calling the function `Tinned::is_zero_quantity`.

## Examples

The following snippet shows how to use Tinned library to create (i) electrical
and geometrical perturbations with zero frequency, (ii) overlap $\mathbf{S}$
and density $\mathbf{D}$ matrices, (iii) the idempotency constraint
$\mathbf{Z}=\mathbf{DSD}-\mathbf{D}$, and finally take the second order
derivatives of $\mathbf{Z}$:

```
#include <iostream>

#include <symengine/constants.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>

#include "Tinned.hpp"

int main()
{
    // Create electrical and geometrical perturbations with zero frequency
    auto el = Tinned::make_perturbation(std::string("EL"));
    auto geo = Tinned::make_perturbation(std::string("GEO"));

    // Create overlap matrix which depends on geometrical perturbation
    auto S = Tinned::make_1el_operator(
        std::string("S"),
        Tinned::PertDependency({std::make_pair(geo, 99)})
    );

    // Create density matrix which depends on all perturbations
    auto D = Tinned::make_1el_density(std::string("D"));

    // Creat the idempotency constraint Z = DSD - D
    auto Z = SymEngine::matrix_add({
        SymEngine::matrix_mul({D, S, D}),
        SymEngine::matrix_mul({SymEngine::minus_one, D})
    });

    // Take the second order derivatives of the idempotency constraint
    auto K = Tinned::differentiate(Z, Tinned::PertTuple({el, geo}));
    std::cout << "Second order derivatives of the idempotency constraint: "
              << Tinned::stringify(K) << "\n";

    return 0;
}
```

More examples can be found in Tinned tests in the directory `tests`.

## TODO

* Add Boolean input for `differentiate` function to indicate if
  call `clean_temporum` after differentiation;
* Introduce `CoefficientMO` inherited from `ElectronicState`;
* Test newly added `TwoElecEnergy` and coupled-cluster classes, and visitors
  `EliminationVisitor`, `TemporumCleaner`, `ZerosRemover`.
* More tests on `ExchCorrPotential::get_potential_map`.
* Add tests for functions `ExchCorrEnergy::get_energy_terms` and
  `ExchCorrPotential::get_potential_terms`.

## References

<a id="1">[1]</a>
Bin Gao, "Tinned: A Symbolic Library for Response Theory and High-Order
Derivatives", J. Comput. Chem., DOI: 10.1002/jcc.27437.

**NOTE** the following changes have been made comparing to reference [[1]](#1):

* Function `find_dependency` is renamed `get_diff_order`;
* Class `TwoElecEnergy` stores a `TwoElecOperator` object `G_(inner_)` instead
  the density matrix `inner_`;
* Class `OneElecOperator` is used to represent basis functions on bra and ket
  in class `TemporumOverlap`;
* Class `TemporumOperator` does not support objects of class `NonElecFunction`
  as its target;
* Introduce class `PerturbedParameter` that replaces classes `LagMultiplier`
  and `StateVector`;
* Classes `AdjointMap` and `ExpAdjointHamiltonian` do not use `name` as their
  member variable, and class `ExpAdjointHamiltonian` is renamed
  `ClusterConjHamiltonian`;
* Add class `ConjugateTranspose`;
* Remove class `StateOperator`, such that the coupled-cluster operator is
  represented by the matrix multiplication of coupled-cluster amplitudes and
  the transpose of excitation operators.
* Type of perturbation frequency is changed to `SymEngine::RCP<const SymEngine::Basic>`.
