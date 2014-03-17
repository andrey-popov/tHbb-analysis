/**
 * \file ClassTTbarHFPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a plugin to classify a ttbar event depending on flavours of additional jets.
 */

#pragma once


#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <PhysicsObjects.hpp>

#include <list>
#include <vector>


/**
 * \class ClassTTbarHFPlugin
 * \brief Classifies a ttbar event depending on flavours of additional jets.
 * 
 * Classification is based on partons after parton shower and reconstructed jets. Supported
 * categories are listed in the enumeration Type. The plugin looks at partons with status 2. Partons
 * that are matched to b quarks from decays of top quarks are excluded from the list. In addition,
 * heavy-flavour quarks that are immediate daughters of beam particles are not considered.
 */
class ClassTTbarHFPlugin: public Plugin
{
public:
    /**
     * \brief Supported categories to classify an event
     * 
     * The classification is top-to-bottom: if an event matches some category, all the following
     * ones are not checked.
     */
    enum class Type
    {
        Unknown,
        TwoBottom,         ///< At least two b quarks matched to different reconstructed jets
        DoubleBottomJet,   ///< Same reconstructed jet matched to at least two b quarks
        SingleBottom,      ///< Only one b quark matched to a reconstructed jet
        TwoCharm,          ///< At least two c quarks matched to different reconstructed jets
        DoubleCharmJet,    ///< Same reconstructed jet matched to at least two c quarks
        SingleCharm,       ///< Only one c quark matched to a reconstructed jet
        UnmatchedBottom,   ///< There are b quarks, but none is matched to a reconstructed jet
        UnmatchedCharm,    ///< There are c quarks, but none is matched to a reconstructed jet
        NoHF               ///< There are no b or c quarks in the event
    };
    
private:
    /// Counts number of b and c quarks matched to a jet
    struct HFCounter
    {
        unsigned nb, nc;
    };
    
public:
    /// Constructor
    ClassTTbarHFPlugin(std::string const &name) noexcept;
    
public:
    /// Creates a newly-initialised copy of *this
    virtual Plugin *Clone() const;
    
    /// Saves a pointer to PECReaderPlugin for convenience
    virtual void BeginRun(Dataset const &dataset);
    
    /**
     * \brief Classifies the current event
     * 
     * Always returns true.
     */
    virtual bool ProcessEvent();
    
    /// Returns the classification decision
    Type GetDecision() const noexcept;
    
private:
    /// Pointer to the reader plugin
    PECReaderPlugin const *reader;
    
    /// Decision on classification of the current event
    Type type;
    
    /// A size of cone to match quarks to reconstructed jets
    double const matchConeSize;
    
    /**
     * \brief Collection of additional heavy-flavour quarks
     * 
     * Does not contain b quarks matched to the status-3 b quarks from decays of tops, nor
     * heavy-flavour quarks that are immediate daughters of beam particles.
     */
    std::list<ShowerParton const *> additionalPartons;
    
    /// Contains numbers of heavy-flavour partons matched to each reconstructed jet
    std::vector<HFCounter> partonCounts;
};
