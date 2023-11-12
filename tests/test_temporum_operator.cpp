#include <cstddef>
#include <map>
#include <string>
#include <utility>

#include <catch2/catch.hpp>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/real_double.h>
#include <symengine/symengine_rcp.h>

#include "Tinned.hpp"

using namespace Tinned;

TEST_CASE("Test TemporumOperator and make_dt_operator()", "[TemporumOperator]")
{
    auto D = make_1el_density(std::string("D"));
    auto Dbra = make_dt_operator(D, TemporumType::Bra);
    auto Dket = make_dt_operator(D, TemporumType::Ket);
    REQUIRE(SymEngine::eq(*D, *Dbra->get_args()[0]));
    REQUIRE(SymEngine::eq(*D, *Dket->get_args()[0]));
    REQUIRE(SymEngine::eq(*D, *Dbra->get_target()));
    REQUIRE(SymEngine::eq(*D, *Dket->get_target()));
    REQUIRE(Dbra->get_type() == TemporumType::Bra);
    REQUIRE(Dket->get_type() == TemporumType::Ket);
    REQUIRE(SymEngine::eq(*Dbra->get_frequency(), *SymEngine::zero));
    REQUIRE(SymEngine::eq(*Dket->get_frequency(), *SymEngine::zero));

    auto el = make_perturbation(std::string("EL"), SymEngine::real_double(0.5));
    auto geo = make_perturbation(std::string("GEO"), SymEngine::real_double(0.0));
    auto mag = make_perturbation(std::string("MAG"), SymEngine::real_double(1.5));
    auto Dp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        ((Dbra->diff(el))->diff(geo))->diff(mag)
    );
    auto sum_freq = SymEngine::addnum(
        SymEngine::addnum(el->get_frequency(), geo->get_frequency()),
        mag->get_frequency()
    );
    REQUIRE(SymEngine::eq(
        *Dp->get_frequency(),
        *SymEngine::subnum(SymEngine::real_double(0), sum_freq)
    ));
    REQUIRE(SymEngine::unified_eq(
        Dp->get_derivative(), SymEngine::multiset_basic({el, geo, mag})
    ));
    Dp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        ((Dket->diff(el))->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::eq(*Dp->get_frequency(), *sum_freq));
    REQUIRE(SymEngine::unified_eq(
        Dp->get_derivative(), SymEngine::multiset_basic({el, geo, mag})
    ));

    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto S = make_1el_operator(std::string("S"), dependencies);
    auto St = make_dt_operator(S);
    auto Sp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        (St->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::unified_eq(
        Sp->get_derivative(), SymEngine::multiset_basic({geo, mag})
    ));
    REQUIRE(SymEngine::eq(*Sp->diff(el), *make_zero_operator()));

    auto h = make_nonel_function(std::string("hnuc"), dependencies);
    auto ht = make_dt_operator(h);
    auto hp = SymEngine::rcp_dynamic_cast<const TemporumOperator>(
        (ht->diff(geo))->diff(mag)
    );
    REQUIRE(SymEngine::unified_eq(
        hp->get_derivative(), SymEngine::multiset_basic({geo, mag})
    ));
    REQUIRE(SymEngine::eq(*hp->diff(el), *SymEngine::zero));
}

TEST_CASE("Test TemporumOverlap and make_t_matrix()", "[TemporumOverlap]")
{
    auto el_freq = SymEngine::real_double(0.5);
    auto el = make_perturbation(std::string("EL"), el_freq);
    auto geo_freq = SymEngine::real_double(1.5);
    auto geo = make_perturbation(std::string("GEO"), geo_freq);
    auto mag_freq = SymEngine::real_double(2.5);
    auto mag = make_perturbation(std::string("MAG"), mag_freq);
    auto dependencies = PertDependency({
        std::make_pair(geo, 99), std::make_pair(mag, 99)
    });
    auto T = make_t_matrix(dependencies);
    REQUIRE(eq_dependency(dependencies, T->get_dependencies()));
    REQUIRE(T->size() == 1);
    REQUIRE(T->get_frequency(0)->is_zero());
    auto derivative = T->get_derivative(0);
    REQUIRE(SymEngine::unified_eq(derivative.first, SymEngine::multiset_basic({})));
    REQUIRE(SymEngine::unified_eq(derivative.second, SymEngine::multiset_basic({})));

    // = wg*S^{gg|0} + 0*2*S^{g|g} - wg*S^{0|gg}
    auto Tgg = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(
        (T->diff(geo))->diff(geo)
    );
    auto found = std::map<std::string, bool>({
        {std::string("gg|0"), false},
        {std::string("g|g"), false},
        {std::string("0|gg"), false}
    });
    REQUIRE(Tgg->size() == 3);
    for (std::size_t i = 0; i < Tgg->size(); ++i) {
        derivative = Tgg->get_derivative(i);
        if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({})
            )
        ) {
            REQUIRE(SymEngine::eq(*Tgg->get_frequency(i), *geo_freq));
            found[std::string("gg|0")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo})
            )
        ) {
            REQUIRE(Tgg->get_frequency(i)->is_zero());
            found[std::string("g|g")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, geo})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tgg->get_frequency(i),
                *SymEngine::subnum(SymEngine::real_double(0.0), geo_freq)
            ));
            found[std::string("0|gg")] = true;
        }
    }
    for (auto it = found.begin(); it != found.end(); ++it) REQUIRE(it->second);

    // = (wg+wb/2)*S^{ggb|0} + (wg-wb/2)*S^{gg|b} + wb*S^{gb|g}
    // - wb*S^{g|gb} + (wb/2-wg)*S^{b|gg} - (wg+wb/2)*S^{0|ggb}
    auto Tggb = SymEngine::rcp_dynamic_cast<const TemporumOverlap>(Tgg->diff(mag));
    found = std::map<std::string, bool>({
        {std::string("ggb|0"), false},
        {std::string("gg|b"), false},
        {std::string("gb|g"), false},
        {std::string("g|gb"), false},
        {std::string("b|gg"), false},
        {std::string("0|ggb"), false}
    });
    REQUIRE(Tggb->size() == 6);
    for (std::size_t i = 0; i < Tggb->size(); ++i) {
        derivative = Tggb->get_derivative(i);
        if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, geo, mag})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::addnum(
                    geo_freq,
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(0.5))
                )
            ));
            found[std::string("ggb|0")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({mag})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(
                    geo_freq,
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(0.5))
                )
            ));
            found[std::string("gg|b")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo, mag})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo})
            )
        ) {
            REQUIRE(SymEngine::eq(*Tggb->get_frequency(i), *mag_freq));
            found[std::string("gb|g")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({geo})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, mag})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(SymEngine::real_double(0.0), mag_freq)
            ));
            found[std::string("g|gb")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({mag})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, geo})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(0.5)),
                    geo_freq
                )
            ));
            found[std::string("b|gg")] = true;
        }
        else if (
            SymEngine::unified_eq(
                derivative.first,
                SymEngine::multiset_basic({})
            )
            &&
            SymEngine::unified_eq(
                derivative.second,
                SymEngine::multiset_basic({geo, geo, mag})
            )
        ) {
            REQUIRE(SymEngine::eq(
                *Tggb->get_frequency(i),
                *SymEngine::subnum(
                    SymEngine::mulnum(mag_freq, SymEngine::real_double(-0.5)),
                    geo_freq
                )
            ));
            found[std::string("0|ggb")] = true;
        }
    }
    for (auto it = found.begin(); it != found.end(); ++it) REQUIRE(it->second);

    REQUIRE(SymEngine::eq(*T->diff(el), *make_zero_operator()));
    REQUIRE(SymEngine::eq(*Tgg->diff(el), *make_zero_operator()));
    REQUIRE(SymEngine::eq(*Tggb->diff(el), *make_zero_operator()));
}
