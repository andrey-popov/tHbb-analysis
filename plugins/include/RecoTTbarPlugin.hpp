/**
 * \file RecoTTbarPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines an abstract base class for a plugin to recontruct a semileptonic ttbar event.
 */

#pragma once

#include <Plugin.hpp>

#include <PhysicsObjects.hpp>

#include <string>
#include <utility>


/**
 * \class RecoTTbarPlugin
 * \brief Abstract base class for reconstruction of a semileptonic ttbar event
 */
class RecoTTbarPlugin: public Plugin
{
public:
    /**
     * \brief Constructor
     * 
     * Accepts a name of the instance of the plugin.
     */
    RecoTTbarPlugin(std::string const &name_):
        Plugin(name_)
    {}
    
public:
    /// Returns the jet assigned to the b quark in the semileptonic decay of a top quark
    virtual Jet const &GetBTopLep() const = 0;
    
    /// Returns the jet assigned to the b quark in the hadronic decay of a top quark
    virtual Jet const &GetBTopHad() const = 0;
    
    /**
     * \brief Returns jets matched to the light-flavour quarks from the hadronic decay of a top
     * 
     * The first jet in the pair has larger transverse momentum.
     */
    virtual std::pair<Jet const &, Jet const &> GetLightTopHad() const = 0;
};
