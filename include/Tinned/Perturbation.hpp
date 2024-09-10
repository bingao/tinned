/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of perturbations.

   2024-09-05, Bin Gao:
   * change Perturbation's frequency to class
     SymEngine::RCP<const SymEngine::Basic> so that SymEngine::Number,
     SymEngine::Symbol, SymEngine::Add, etc. can be used as the frequency

   2023-10-09, Bin Gao:
   * change Perturbation's frequency to class
     SymEngine::RCP<const SymEngine::Number>. Users can use any derived class
     from SymEngine::Number as frequencies

   2023-09-16, Bin Gao:
   * first version
*/

#pragma once

#include <cstddef>
#include <set>
#include <string>

#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/constants.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

namespace Tinned
{
    class Perturbation: public SymEngine::Symbol
    {
        protected:
            // Frequency
            SymEngine::RCP<const SymEngine::Basic> frequency_;
            // Set of components
            std::set<std::size_t> components_;

        public:
            //! Constructor
            explicit Perturbation(
                const std::string& name,
                const SymEngine::RCP<const SymEngine::Basic>& frequency,
                const std::set<std::size_t>& components
            );

            SymEngine::hash_t __hash__() const override;
            bool __eq__(const SymEngine::Basic& o) const override;
            int compare(const SymEngine::Basic& o) const override;

            //! Get the frequency of the perturbation
            inline SymEngine::RCP<const SymEngine::Basic> get_frequency() const
            {
                return frequency_;
            }

            //! Get the set of components of the perturbation
            inline std::set<std::size_t> get_components() const
            {
                return components_;
            }
    };

    // Helper function to make a perturbation
    inline SymEngine::RCP<const Perturbation> make_perturbation(
        const std::string& name,
        const SymEngine::RCP<const SymEngine::Basic>& frequency = SymEngine::zero,
        const std::set<std::size_t>& components = std::set<std::size_t>()
    )
    {
        return SymEngine::make_rcp<const Perturbation>(name, frequency, components);
    }
}
