/**
 * \file RecoTPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines an abstract base class for a plugin to recontruct a t-channel single-top event.
 */

#pragma once

#include <Plugin.hpp>

#include <PhysicsObjects.hpp>

#include <string>
#include <utility>


/**
 * \class RecoTPlugin
 * \brief Abstract base class for reconstruction of a t-channel single-top event
 */
class RecoTPlugin: public Plugin
{
public:
    /**
     * \brief Constructor
     * 
     * Accepts a name of the instance of the plugin.
     */
    RecoTPlugin(std::string const &name_):
        Plugin(name_)
    {}
    
public:
    /// Returns the reconstructed recoil jet
    virtual Jet const &GetRecoilJet() const = 0;
    
    /// Returns the jet assigned to the b quark in the semileptonic decay of a top quark
    virtual Jet const &GetBTop() const = 0;
};
