/**
 * \file TopPtWeightPlugin.hpp
 * \author Andrey Popov
 *
 * A draft of a plugin to perform top pt reweighting [1].
 * [1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
 */

#pragma once

#include <Plugin.hpp>

#include <PECReaderPlugin.hpp>


/**
 * \class TopPtWeightPlugin
 * \brief (Draft)
 * 
 * 
 */
class TopPtWeightPlugin: public Plugin
{
public:
    /**
     * \brief Constructor
     */
    TopPtWeightPlugin();
    
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
     * Only updates the pointer to the reader plugin.
     */
    void BeginRun(Dataset const &dataset);
    
    /// Calculates event weight
    bool ProcessEvent();
    
    /// Retrieves the calculated event weight
    double GetWeight() const;

private:
    /// Formula to calculate top-pt weight
    static double CalculateTopPtWeight(double pt1, double pt2, double a, double b);
    
private:
    /// Pointer to the reader plugin
    PECReaderPlugin const *reader;
    
    /**
     * \brief Shows if the current dataset is a ttbar one
     * 
     * If not, the plugin returns 1. as weight.
     */
    bool isTTbar;
    
    /// A buffer to hold event weight
    double weight;
};
