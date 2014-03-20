/**
 * \file MassChi2RecoTPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a plugin to recontruct a t-channel single-top event.
 */

#pragma once

#include <RecoTPlugin.hpp>

#include <PECReaderPlugin.hpp>

#include <TMVA/Reader.h>

#include <string>
#include <vector>
#include <utility>


/**
 * \class MassChi2RecoTPlugin
 * \brief Reconstructs a t-channel single-top event with a simple mass chi2 minimisation
 * 
 * A central jet is matched to the b quark from decay of the top quark. The recoil jet is
 * reconstructed as the most forward of the remaining jets. Same jet cannot be matched to several
 * objects.
 */
class MassChi2RecoTPlugin: public RecoTPlugin
{
public:
    /**
     * \brief Constructor
     * 
     * The given name is forwarded to a constructor of Plugin.
     */
    MassChi2RecoTPlugin(std::string const &name);
    
public:
    /**
     * \brief Constructs a newly-initialised copy
     * 
     * Consult documentation of the overriden method for details.
     */
    virtual Plugin *Clone() const;
    
    /**
     * \brief Notifies this that a dataset has been opened
     * 
     * Only saves pointer to an instance of PECReaderPlugin for convenience.
     */
    virtual void BeginRun(Dataset const &);
    
    /**
     * \brief Processes the current event
     * 
     * Finds decay product of the top quark and the recoil jet. The method always returns
     * true. If some of the objects cannot be reconstructed, which would mean the event selection
     * is inappropriate, an exception of type std::logic_error is thrown.
     */
    virtual bool ProcessEvent();
    
    /// Returns the reconstructed recoil jet
    virtual Jet const &GetRecoilJet() const;
    
    /// Returns the jet assigned to the b quark in the semileptonic decay of a top quark
    virtual Jet const &GetBTop() const;
    
private:
    /// Pointer to PECReaderPlugin
    PECReaderPlugin const *reader;
    
    /// Index of the light-flavour recoil jet
    unsigned iRecoilJet;
    
    /// Index of the jet that corresponds to the b quark from semileptonic decay of a top quark
    unsigned iBTop;
};
