# Tinned

Tinned contains a set of nonnumerical routines for computational chemistry, in
which symbolic routines are built on top of
[SymEngine library](https://github.com/symengine/symengine).

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

## Try Tinned

Tinned currently provides C++ interface and the following classes and functions
for computational chemistry:

* [`class Perturbation`](include/Tinned/Perturbation.hpp), perturbations;
* [`class PertDependency`](include/Tinned/PertDependency.hpp), `std::set` for
  perturbations that an operator depends on and their maximum orders that can
  be differentiated;
* [`class OneElecDensity`](include/Tinned/OneElecDensity.hpp), one-electron
  spin-orbital density matrix derived from the abstract electronic state
  [`class ElectronicState`](include/Tinned/ElectronicState.hpp);
* [`class OneElecOperator`](include/Tinned/OneElecOperator.hpp), one-electron
  like operators;
* [`class TwoElecOperator`](include/Tinned/TwoElecOperator.hpp), two-electron
  like operators;
* [`class ExchCorrEnergy`](include/Tinned/ExchCorrEnergy.hpp),
  exchange-correlation energy like functionals;
* [`class ExchCorrPotential`](include/Tinned/ExchCorrPotential.hpp),
  exchange-correlation potential like operators;
* [`class NonElecFunction`](include/Tinned/NonElecFunction.hpp), non-electron
  like functions;
* [`class TemporumOperator`](include/Tinned/TemporumOperator.hpp), for
  operators $\langle\mathrm{i}\frac{\partial}{\partial t}\vert$ and
  $\vert\mathrm{i}\frac{\partial}{\partial t}\rangle$;
* [`class TemporumOverlap`](include/Tinned/TemporumOverlap.hpp), for the operator
  $-\frac{1}{2}(\langle\mathrm{i}\frac{\partial}{\partial t}\chi_{\kappa}\vert\chi_{\lambda}\rangle
  +\langle\chi_{\kappa}\vert\mathrm{i}\frac{\partial}{\partial t}\chi_{\lambda}\rangle)$;
* function [`keep_if(x, symbols)`](include/Tinned/KeepVisitor.hpp) keeps given
  `symbols` in `x` while removing others;
* function [`remove_if(x, symbols)`](include/Tinned/RemoveVisitor.hpp) removes
  given `symbols` from `x`;
* function [`replace(x, subs_dict, cache)`](include/Tinned/ReplaceVisitor.hpp)
  replaces classes defined in Tinned library in addition to those in
  `SymEngine::msubs()`;
* function [`stringify(x)`](include/Tinned/StringifyVisitor.hpp) stringifies
  symbols from SymEngine and additional ones defined in Tinned library.

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
    auto el = SymEngine::make_rcp<const Tinned::Perturbation>(
        std::string("EL")
    );
    auto geo = SymEngine::make_rcp<const Tinned::Perturbation>(
        std::string("GEO")
    );

    // Create overlap matrix which depends on geometrical perturbation
    auto S = SymEngine::make_rcp<const Tinned::OneElecOperator>(
        std::string("S"),
        Tinned::PertDependency({
            std::make_pair(geo, 99)
        })
    );

    // Create density matrix which depends on all perturbations
    auto D = SymEngine::make_rcp<const Tinned::OneElecDensity>(
        std::string("D")
    );

    // Creat the idempotency constraint Z = DSD - D
    auto Z = SymEngine::matrix_add({
        SymEngine::matrix_mul({D, S, D}),
        SymEngine::matrix_mul({SymEngine::minus_one, D})
    });

    // Take the second order derivatives of the idempotency constraint
    auto K = (Z->diff(geo))->diff(el);
    std::cout << "Second order derivatives of the idempotency constraint: "
              << Tinned::stringify(*K) << "\n";

    return 0;
}
```

More examples can be found in Tinned tests in the directory `tests`.
