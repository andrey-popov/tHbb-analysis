/**
 * \file GlobalRecoTHPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a plugin to recontruct a thq event.
 */

#pragma once

#include <RecoTHPlugin.hpp>

#include <PECReaderPlugin.hpp>
#include <BTaggerPlugin.hpp>

#include <TMVA/Reader.h>

#include <string>
#include <utility>


/**
 * \class GlobalRecoTHPlugin
 * \brief Reconstructs a tHq event exploiting correlations between different objects
 * 
 * \warning It considers only jets from PECReader::GetJets(), but not the "additional" jets.
 */
class GlobalRecoTHPlugin: public RecoTHPlugin
{
public:
    /**
     * \brief Constructor
     * 
     * The arguments are: name for this plugin and the name of a plugin to perform b-tagging.
     */
    GlobalRecoTHPlugin(std::string const &name, std::string const &bTagPluginName);
    
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
    
    /// Returns number of considered interpretations
    unsigned GetNumInterpretations() const;
    
    /// Returns MVA response for the chosen interpretation
    double GetMvaResponse() const;

private:
    /**
     * \brief Calculates input variables
     * 
     * 
     */
    void CalculateVariables(Jet const &bTop, Jet const &b1Higgs, Jet const &b2Higgs,
     Jet const &recoil);

private:
    /// Name of a plugin to perform b-tagging
    std::string bTagPluginName;
    
    /// Pointer to PECReaderPlugin
    PECReaderPlugin const *reader;
    
    /// Pointer to a b-tagging plugin
    BTaggerPlugin const *bTagger;
    
    /// MVA to reconstruct an event
    TMVA::Reader mvaReco;
    
    /// Index of the light-flavour recoil jet
    unsigned iRecoilJet;
    
    /// Index of the jet that corresponds to the b quark from semileptonic decay of a top quark
    unsigned iBTop;
    
    /// Indices of jets that correspond to b quarks from decay of a Higgs boson
    unsigned iB1Higgs, iB2Higgs;
    
    
    /// Technical information: number of considered interpretations
    unsigned numInterpretations;
    
    /// Technical information: highest MVA response
    double highestMvaResponse;
    
    
    // Buffers to store input variables for event reconstruction
    float LogMass_Higgs, LogMass_BTopLep;
    float AbsEta_Recoil;
    float PassBTag_BTop, NumBTag_Higgs;
    float DeltaR_BJetsHiggs, DeltaR_BTopW, DeltaR_TopHiggs;
    float RelHt;
    float MaxEta_BHiggs, LogMinPt_BHiggs;
    float Charge_BTop;
};
