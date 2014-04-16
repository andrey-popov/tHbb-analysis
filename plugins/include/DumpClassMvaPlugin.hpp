/**
 * \file DumpClassMvaPlugin.hpp
 * \author Andrey Popov
 * 
 * Defines a plugin to calculate response of the classification MVA. The plugin might be configured
 * to store the reponse and few auxiliary variables in a ROOT tree.
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <TopPtWeightPlugin.hpp>
#include <RecoTHPlugin.hpp>
#include <RecoTTbarPlugin.hpp>
#include <BTaggerPlugin.hpp>
#include <EventWeightPlugin.hpp>

#include <TMVA/Reader.h>

#include <string>
#include <memory>


/**
 * \class DumpClassMvaPlugin
 * \brief Calculates response of the classification MVA
 * 
 * Depending on the configuration, the plugin either stores the response and few additional
 * variables in a ROOT file or returns the response with the help of a dedicated method only.
 */
class DumpClassMvaPlugin: public Plugin
{
public:
    /**
     * \brief Constructor that configures the plugin to store the MVA response in a ROOT file
     * 
     * Accepts the name of the directory to host the produced files and a postfix to be applied
     * to the names of output files.
     */
    DumpClassMvaPlugin(std::string const &name, std::string const &bTagPluginName,
     std::string const &weightPluginName, std::string const &outDirectory,
     std::string const &filePostfix, bool saveSystWeights);
    
    /**
     * \brief Constructor that configures the plugin not to store the MVA response
     */
    DumpClassMvaPlugin(std::string const &name, std::string const &bTagPluginName);
    
public:
    /**
     * \brief Constructs a newly-initialised copy
     * 
     * Consult documentation of the overriden method for details.
     */
    Plugin *Clone() const;
    
    /**
     * \brief Notifies this that a dataset has been opened
     * 
     * Creates the output tree. Consult documentation of the overriden method for details.
     */
    void BeginRun(Dataset const &dataset);
    
    /**
     * \brief Notifies this that a dataset has been closed
     * 
     * Writes down the current output tree. Consult documentation of the overriden method for
     * details.
     */
    void EndRun();
    
    /**
     * \brief Processes the current event
     * 
     * Calculates and stores the input variables. The method always returns true, despite the fact
     * that the input variables are written for perfectly matched events only. Consult documentation
     * of the overriden method for details.
     */
    bool ProcessEvent();
    
    /// Returns response of the classification MVA for the current event
    float GetResponse() const noexcept;
    
private:
    /// Books an appropriate MVA
    void BookMVA();
    
private:
    /// Name of a b-tagging plugin
    std::string bTagPluginName;
    
    /**
     * \brief Name of a plugin to calculate event weight
     * 
     * Currently, used only for the data-driven ttbar.
     */
    std::string weightPluginName;
    
    /// Pointer to PECReaderPlugin
    PECReaderPlugin const *reader;
    
    /// Pointer to a plugin to perform top-pt reweighting (optional)
    TopPtWeightPlugin const *topPtReweighter;
    
    /// Pointer to a plugin that performs reconstruction of an event under thq hypothesis
    RecoTHPlugin const *builderTopHiggs;
    
    /// Pointer to a plugin that performs reconstruction of an event under ttbar hypothersis
    RecoTTbarPlugin const *builderTTbar;
    
    /// Pointer to a b-tagging plugin
    BTaggerPlugin const *bTagger;
    
    /// Pointer to a plugin for additional reweighting
    EventWeightPlugin const *reweighter;
    
    /// Indicates whether MVA output should be stored in a ROOT file
    bool storeResponse;
    
    /// Classification MVA
    TMVA::Reader mva;
    
    /// Directory to store output files
    std::string outDirectory;
    
    /// Postfix to be applied to the name of output file
    std::string filePostfix;
    
    /// Says if systematical weights should be saved
    bool saveSystWeights;
    
    /// Current output file
    std::unique_ptr<TFile> file;
    
    /// Current output tree
    std::unique_ptr<TTree> tree;
    
    
    // Output buffers
    ULong64_t bfEventNumber, bfRunNumber, bfLumiSection;
    
    Int_t bfNumWeights;
    Float_t bfWeights[16];
    //^ 0: central weight; 1, 2: pile-up up and down; 3, 4: tag rate up and down; 5, 6: mistag
    //rate up and down
    
    Int_t bfNumTag20;
    
    Float_t bfMvaResponse;
    
    
    // Input variables for MVA
    Float_t glb_Charge_Lep;
    Float_t glb_LogPt_J1;
    Float_t glb_LogSqrtSHat;  // mass of sum of all objects
    Float_t glb_Sphericity, glb_LogMET;
    
    
    Float_t thq_Cos_LepRecoil_TH;
    //^ Angle between three-momenta of the lepton and the light-flavour jet in the rest frame of
    //t+h system
    
    Float_t thq_DeltaR_BJetsHiggs;
    Float_t thq_AbsEta_Recoil;
    Float_t thq_NumBTag_Higgs;
    Float_t thq_LogPt_Higgs;
    Float_t thq_LogPt_Recoil;
    
    
    Float_t tt_Charge_BTopHad;  // multiplied by the charge of the lepton
    Float_t tt_DeltaR_Light;
    Float_t tt_LogDMass_TopHadWHad;
    Float_t tt_LogMass_TopHad;
    Float_t tt_MaxEta_Light;
    Float_t tt_LogMaxMass_BTopHadLight;
    Float_t tt_NumPassBTag_Light;
    Float_t tt_SumCharge_Light;  // multiplied by the charge of the lepton
};
