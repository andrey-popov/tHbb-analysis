/**
 * \file JetTagDataDrivenPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines an abstract base class for a plugin to perform a data-driven estimation based on jet
 * tags of any kind.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <PhysicsObjects.hpp>

#include <string>
#include <vector>


/**
 * \class JetTagDataDrivenPlugin
 * \brief Abstract base class for a data-driven estimation
 * 
 * The class provides an interface for a data-driven estimation that relies on jet tagging. Specific
 * way of tagging is not important. A derived class can assign each event a weight and can specify
 * which jets should be considered tagged, this information is communicated to the outer world with
 * the dedicated methods GetWeight and GetTaggedJetIndices. It is not known in which order these
 * methods will be called and if they will be executed at all.
 */
class JetTagDataDrivenPlugin: public Plugin
{
public:
    /**
     * \brief Constructor
     * 
     * Forwards the given name to the constructor of Plugin.
     */
     JetTagDataDrivenPlugin(std::string const &name) noexcept;
    
    /// Trivial destructor
    virtual ~JetTagDataDrivenPlugin() noexcept;
     
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
     * \brief Returns event weight
     * 
     * The method might be called several times for each event. The second and subsequent calls
     * should be cheap.
     */
    virtual double GetWeight() const = 0;
    
    /**
     * \brief Returns a vector of indices of tagged jets
     * 
     * The indices must correpond to the vector returned by PECReader::GetJets. Default
     * implementation simply returns a reference to taggedJetIndices.
     * 
     * The returned vector must be sorted in incresing order. The method might be called several
     * times for each event. The second and subsequent calls should be cheap.
     */
    virtual std::vector<unsigned> const &GetTaggedJetIndices() const;
    
    /**
     * \brief Checks if the given jet should be considered as tagged
     * 
     * First the method finds to which index in the collection PECReader::GetJets the given jet
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
    
    /// Indices of jets that should be considered tagged
    std::vector<unsigned> taggedJetIndices;
};
