#ifndef SYMENGINE_FUNCTIONS_H
#define SYMENGINE_FUNCTIONS_H

#include <symengine/basic.h>
#include <symengine/symengine_casts.h>
#include <symengine/constants.h>
#include <symengine/functions.h>

// is a tensor contraction of ERI and density matrix
class TwoElectronOperator : public SymEngine::MultiArgFunction
{
private:
    SymEngine::RCP<const ElectronState> state_;
public:
    SymEngine::IMPLEMENT_TYPEID(SymEngine::SYMENGINE_MAX)
    //! TwoElectronOperator Constructor
    TwoElectronOperator(const SymEngine::vec_basic &&arg);
    //! \return `true` if canonical
    bool is_canonical(const SymEngine::vec_basic &arg) const;
    //! \return canonicalized TwoElectronOperator
    SymEngine::RCP<const SymEngine::Basic> create(const SymEngine::vec_basic &arg) const override;
};

//! Canonicalize TwoElectronOperator:
SymEngine::RCP<const SymEngine::Basic> two_electron_operator(const SymEngine::vec_basic &arg);

#endif
