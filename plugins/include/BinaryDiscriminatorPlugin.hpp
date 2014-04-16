/**
 * \file BinaryDiscriminatorPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines an abstract base class for plugins to perform binary discrimination of events.
 */

#pragma once

#include <Plugin.hpp>


/**
 * \class BinaryDiscriminatorPlugin
 * \brief An abstract base class for plugins to perform binary discrimination of events
 */
class BinaryDiscriminatorPlugin: public Plugin
{
public:
    /**
     * \brief Constructor
     * 
     * Forwards the given name of the instance to the base class
     */
    BinaryDiscriminatorPlugin(std::string const &name);
    
public:
    /**
     * \brief Returns discriminator's decision taken for the current event
     * 
     * Must be implemented by user. The method is expected to be fast and must not perform heavy
     * calculations.
     */
    virtual double GetResponse() const noexcept = 0;
};
