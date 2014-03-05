/**
 * \file TTbarDataDrivePlugin.hpp
 * \author Daniel Knowlton
 * 
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>
#include <PhysicsObjects.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TRandom3.h>

#include <string>
#include <vector>


/**
 * \class TTbarDataDrivePlugin
 */
class TTbarDataDrivePlugin: public Plugin
{
public:
  /**
   * \brief Constructor
   */
  TTbarDataDrivePlugin(std::string const &name, int region) noexcept;
    
  /// Trivial destructor
  virtual ~TTbarDataDrivePlugin() noexcept;
     
public:
  /**
   * \brief Creates a newly-initialized copy
   * 
   * Consult documentation of the overriden method for details.
   */
  Plugin *Clone() const;

  virtual void BeginRun(Dataset const &dataset);
    
  virtual double GetWeight() const = 0;
    
  virtual std::vector<int> const &GetTaggedJetIndices() const;

  /**
   * \brief Processes the current event
   * 
   * Consult documentation of the overriden method for details.
   */
  bool ProcessEvent();
    
protected:
  /// Pointer to PECReaderPlugin
  PECReaderPlugin const *reader;

  std::string const *name;
  int region;

  bool isMC;

  TFile *eff2D;
  TH1F *effbpt;
  TH1F *effcpt;
  TH1F *effudspt;
  TH1F *effgpt;
  TH1F *effbeta;
  TH1F *effceta;
  TH1F *effudseta;
  TH1F *effgeta;
  TH1F *effb2d;
  TH1F *effc2d;
  TH1F *effuds2d;
  TH1F *effg2d;
  TH1F *number_jets;

  // Flavor combinations and their corresponding probabilities
  // Can later be put into a nested map
  std::map < float, std::vector<int> > FlavFracs4, FlavFracs5, FlavFracs6, FlavFracs7, FlavFracs8, FlavFracs9, FlavFracs10;

  /// Vector off all jets in an event
  std::vector<Jet> allJets;

  TRandom3 r3;
  
  float ttweight3T = 0.;
  float ttweight4T = 0.;
  std::vector<int> taggedjets3t;
  std::vector<int> taggedjets4t;

};
