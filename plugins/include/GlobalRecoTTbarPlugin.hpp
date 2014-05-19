/**
 * \file GlobalRecoTTbarPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a plugin to recontruct a ttbar event.
 */

#pragma once

#include <RecoTTbarPlugin.hpp>

#include <PECReaderPlugin.hpp>
#include <BTaggerPlugin.hpp>

#include <TMVA/Reader.h>

#include <string>
#include <utility>


/**
 * \class GlobalRecoTTbarPlugin
 * \brief Reconstructs a ttbar event exploiting correlations between the objects
 * 
 * \warning It considers only jets from PECReader::GetJets(), but not the "additional" jets.
 */
class GlobalRecoTTbarPlugin: public RecoTTbarPlugin
{
public:
    /**
     * \brief Constructor
     * 
     * The arguments are: name for this plugin and the name of a plugin to perform b-tagging.
     */
    GlobalRecoTTbarPlugin(std::string const &name, std::string const &bTagPluginName);
    
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
     * Identifies jets as decay products of top quarks. The method always returns true. If some of
     * the objects cannot be reconstructed, which would mean the event selection is inappropriate,
     * an exception of type std::logic_error is thrown.
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
    /**
     * \brief Calculates input variables
     * 
     * 
     */
    void CalculateVariables(Jet const &bTopLep, Jet const &bTopHad, Jet const &q1TopHad,
     Jet const &q2TopHad);

private:
    /// Name of a plugin to perform b-tagging
    std::string bTagPluginName;
    
    /// Pointer to PECReaderPlugin
    PECReaderPlugin const *reader;
    
    /// Pointer to a b-tagging plugin
    BTaggerPlugin const *bTagger;
    
    /// MVA to peform event reconstruction
    TMVA::Reader mvaReco;
    
    /// Index of the jet that corresponds to the b quark from semileptonic decay of a top quark
    unsigned iBTopLep;
    
    /// Index of the jet that corresponds to the b quark from hadronic decay of a top quark
    unsigned iBTopHad;
    
    /// Indices of jets that correspond to light-flavour quarks from hadronic decay of a top
    unsigned iL1TopHad, iL2TopHad;
    
    
    // Buffers to store input variables for event reconstruction
    float LogDMass_TopHadWHad, LogMass_WHad, LogMass_BTopLepLep;
    float NumBTag_Light;
    float DeltaR_Light, DeltaR_BTopHadWHad, DeltaR_BTopLepWLep;
    float RelHt, LogPt_TopLep, LogPt_TopHad, AbsEta_TopHad;
    float DCharge_BTopHadBTopLep;  // (Q(b_thad) - Q(b_tlep)) * Q(l)
    float SumCharge_Light;
};
