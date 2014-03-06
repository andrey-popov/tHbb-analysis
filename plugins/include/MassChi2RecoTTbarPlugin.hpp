/**
 * \file MassChi2RecoTTbarPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a plugin to reconstruct a semileptonic ttbar event.
 */

#pragma once

#include <RecoTTbarPlugin.hpp>

#include <PECReaderPlugin.hpp>

#include <TMVA/Reader.h>

#include <string>
#include <vector>
#include <utility>


/**
 * \class MassChi2RecoTTbarPlugin
 * \brief Reconstructs a ttbar event with a simple mass chi2 minimisation
 * 
 * Only central jets are matched b quark from decays of the top quarks. A jet cannot be matched to
 * several objects.
 */
class MassChi2RecoTTbarPlugin: public RecoTTbarPlugin
{
public:
    /**
     * \brief Constructor
     * 
     * The given name is forwarded to a constructor of Plugin.
     */
    MassChi2RecoTTbarPlugin(std::string const &name);
    
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
     * Identifies jets as decay products of top quark and Higgs boson. The method always returns
     * true. If some of the objects cannot be reconstructed, which would mean the event selection
     * is inappropriate, an exception of type std::logic_error is thrown.
     */
    virtual bool ProcessEvent();
    
    /// Returns the jet assigned to the b quark in the semileptonic decay of a top quark
    virtual Jet const &GetBTopLep() const;
    
    /// Returns the jet assigned to the b quark in the hadronic decay of a top quark
    virtual Jet const &GetBTopHad() const;
    
    /**
     * \brief Returns jets matched to the light-flavour quarks from the hadronic decay of a top
     * 
     * The first jet in the pair has larger transverse momentum.
     */
    virtual std::pair<Jet const &, Jet const &> GetLightTopHad() const;
    
private:
    /// Pointer to PECReaderPlugin
    PECReaderPlugin const *reader;
    
    /// Index of the jet that corresponds to the b quark from semileptonic decay of a top quark
    unsigned iBTopLep;
    
    /// Index of the jet that corresponds to the b quark from hadronic decay of a top quark
    unsigned iBTopHad;
    
    /// Indices of jets that correspond to light-flavour quarks from hadronic decay of a top
    unsigned iL1TopHad, iL2TopHad;
    
    /**
     * \brief Vector to keep indices of central jets
     * 
     * It is placed in the class description to avoid memory reallocation in each event.
     */
    std::vector<unsigned> centralJetIndices;
};
