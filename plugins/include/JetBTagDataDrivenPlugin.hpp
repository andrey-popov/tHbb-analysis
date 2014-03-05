/**
 * \file JetBTagDataDrivenPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines an abstract base class for a plugin to perform a data-driven estimation based on jet
 * b-tags.
 */

#pragma once

#include <Plugin.hpp>

#include <BTaggerPlugin.hpp>
#include <EventWeightPlugin.hpp>
#include <PECReaderPlugin.hpp>
#include <PhysicsObjects.hpp>

#include <string>
#include <deque>


/**
 * \class JetBTagDataDrivenPlugin
 * \brief Abstract base class for a data-driven estimation based on b-tagging
 * 
 * The class provides an interface for a data-driven estimation that relies on b-tagging. A derived
 * class is expected to fill the mask jetTagBits so that it indicates which jets should be
 * considered as b-tagged. Several convenience functions to access b-tagging information are
 * provided.
 */
class JetBTagDataDrivenPlugin: public BTaggerPlugin, public EventWeightPlugin
{
public:
    /**
     * \brief Constructor
     * 
     * Forwards the given name to the constructor of Plugin.
     */
    JetBTagDataDrivenPlugin(std::string const &name) noexcept;
    
    /// Copy constructor
    JetBTagDataDrivenPlugin(JetBTagDataDrivenPlugin const &src) noexcept;
    
    /// Trivial destructor
    virtual ~JetBTagDataDrivenPlugin() noexcept;
     
public:
    /**
     * \brief Notifies *this  that a new dataset has been opened
     * 
     * In the default implementation the method only sets the pointer to the reader plugin. If the
     * method is overriden in a derived class, the pointer must be initialised properly or this
     * method should be called explicitly from its override.
     */
    virtual void BeginRun(Dataset const &dataset);
    
    /**
     * \brief Returns a mask that indicates which jets should be considered as b-tagged
     * 
     * The indices must correpond to the vector returned by PECReader::GetJets. Default
     * implementation simply returns a reference to jetTagBits.
     * 
     * The method might be called several times for each event. The second and subsequent calls
     * should be cheap.
     */
    virtual std::deque<bool> const &GetJetTagDecisions() const;
    
    /**
     * \brief Checks if the given jet should be considered as tagged
     * 
     * First, the method finds to which index in the collection PECReader::GetJets the given jet
     * corresponds. If not match is found, an exception is thrown. Otherwise, evaluation is
     * delegated to IsTagged(unsigned).
     */
    bool IsTagged(Jet const &jet) const;
    
    /**
     * \brief Checks if the given index corresponds to a tagged jet
     */
    bool IsTagged(unsigned index) const;
    
protected:
    /// Pointer to PECReaderPlugin
    PECReaderPlugin const *reader;
    
    /**
     * \brief A mask that specifies which jets should be considered as b-tagged
     * 
     * If PECReader::GetJets().at(i) should be b-tagged, jetTagBits.at(i) is true.
     */
    std::deque<bool> jetTagBits;
};
