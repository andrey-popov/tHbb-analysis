/**
 * \file MassChi2RecoTHPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a plugin to recontruct a tHq event.
 */

#pragma once

#include <RecoTHPlugin.hpp>

#include <PECReaderPlugin.hpp>

#include <TMVA/Reader.h>

#include <string>
#include <vector>
#include <utility>


/**
 * \class MassChi2RecoTHPlugin
 * \brief Reconstructs a thq event with a simple mass chi2 minimisation
 * 
 * Only central jets are matched to the top quark and the Higgs boson. The recoil jet is
 * reconstructed as the most forward of the remaining jets. A jet cannot be matched to several
 * objects.
 */
class MassChi2RecoTHPlugin: public RecoTHPlugin
{
public:
    /**
     * \brief Constructor
     * 
     * The given name is forwarded to a constructor of Plugin.
     */
    MassChi2RecoTHPlugin(std::string const &name);
    
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
    
    /// Returns the reconstructed recoil jet
    virtual Jet const &GetRecoilJet() const;
    
    /// Returns the jet assigned to the b quark in the semileptonic decay of a top quark
    virtual Jet const &GetBTop() const;
    
    /**
     * \brief Returns jets matched to the b quarks from decay of a Higgs boson
     * 
     * The first jet in the pair has larger transverse momentum.
     */
    virtual std::pair<Jet const &, Jet const &> GetBHiggs() const;
    
private:
    /// Pointer to PECReaderPlugin
    PECReaderPlugin const *reader;
    
    /// Index of the light-flavour recoil jet
    unsigned iRecoilJet;
    
    /// Index of the jet that corresponds to the b quark from semileptonic decay of a top quark
    unsigned iBTop;
    
    /// Indices of jets that correspond to b quarks from decay of a Higgs boson
    unsigned iB1Higgs, iB2Higgs;
    
    /**
     * \brief Vector to keep indices of central jets
     * 
     * It is placed in the class description to avoid memory reallocation in each event.
     */
    std::vector<unsigned> centralJetIndices;
};
