/**
 * \file RecoTHPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines an abstract base class for a plugin to recontruct a thq event.
 */

#pragma once

#include <Plugin.hpp>

#include <PhysicsObjects.hpp>

#include <string>
#include <utility>


/**
 * \class RecoTHPlugin
 * \brief Abstract base class for reconstruction of a thq event
 */
class RecoTHPlugin: public Plugin
{
public:
    /**
     * \brief Constructor
     * 
     * Accepts a name of the instance of the plugin.
     */
    RecoTHPlugin(std::string const &name_):
        Plugin(name_)
    {}
    
public:
    /// Returns the reconstructed recoil jet
    virtual Jet const &GetRecoilJet() const = 0;
    
    /// Returns the jet assigned to the b quark in the semileptonic decay of a top quark
    virtual Jet const &GetBTop() const = 0;
    
    /**
     * \brief Returns jets matched to the b quarks from decay of a Higgs boson
     * 
     * The first jet in the pair has larger transverse momentum.
     */
    virtual std::pair<Jet const &, Jet const &> GetBHiggs() const = 0;
};
