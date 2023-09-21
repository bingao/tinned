#include <symengine/symengine_exception.h>

#include "TwoElectronOperator.h"

using namespace SymEngine;

TwoElectronOperator::TwoElectronOperator(const vec_basic &&arg) : MultiArgFunction(std::move(arg))
{
    SYMENGINE_ASSIGN_TYPEID()
    SYMENGINE_ASSERT(is_canonical(get_vec()))
}

bool TwoElectronOperator::is_canonical(const vec_basic &arg) const
{
    if (arg.size() < 2)
        return false;

    bool non_number_exists = false;

    for (const auto &p : arg) {
        if (is_a<Complex>(*p) or is_a<TwoElectronOperator>(*p))
            return false;
        if (not is_a_Number(*p))
            non_number_exists = true;
    }
    if (not std::is_sorted(arg.begin(), arg.end(), RCPBasicKeyLess()))
        return false;

    return non_number_exists; // all arguments cant be numbers
}

RCP<const Basic> TwoElectronOperator::create(const vec_basic &a) const
{
    return two_electron_operator(a);
}

RCP<const Basic> two_electron_operator(const vec_basic &arg)
{
    bool number_set = false;
    RCP<const Number> max_number, difference;
    set_basic new_args;

    for (const auto &p : arg) {
        if (is_a<Complex>(*p))
            throw SymEngineException("Complex can't be passed to max!");

        if (is_a_Number(*p)) {
            if (not number_set) {
                max_number = rcp_static_cast<const Number>(p);

            } else {
                if (eq(*p, *Inf)) {
                    return Inf;
                } else if (eq(*p, *NegInf)) {
                    continue;
                }
                difference = down_cast<const Number &>(*p).sub(*max_number);

                if (difference->is_zero() and not difference->is_exact()) {
                    if (max_number->is_exact())
                        max_number = rcp_static_cast<const Number>(p);
                } else if (difference->is_positive()) {
                    max_number = rcp_static_cast<const Number>(p);
                }
            }
            number_set = true;

        } else if (is_a<TwoElectronOperator>(*p)) {
            for (const auto &l : down_cast<const TwoElectronOperator &>(*p).get_args()) {
                if (is_a_Number(*l)) {
                    if (not number_set) {
                        max_number = rcp_static_cast<const Number>(l);

                    } else {
                        difference = rcp_static_cast<const Number>(l)->sub(
                            *max_number);

                        if (difference->is_zero()
                            and not difference->is_exact()) {
                            if (max_number->is_exact())
                                max_number = rcp_static_cast<const Number>(l);
                        } else if (difference->is_positive()) {
                            max_number = rcp_static_cast<const Number>(l);
                        }
                    }
                    number_set = true;
                } else {
                    new_args.insert(l);
                }
            }
        } else {
            new_args.insert(p);
        }
    }

    if (number_set)
        new_args.insert(max_number);

    vec_basic final_args(new_args.size());
    std::copy(new_args.begin(), new_args.end(), final_args.begin());

    if (final_args.size() > 1) {
        return make_rcp<const TwoElectronOperator>(std::move(final_args));
    } else if (final_args.size() == 1) {
        return final_args[0];
    } else {
        throw SymEngineException("Empty vec_basic passed to max!");
    }
}
