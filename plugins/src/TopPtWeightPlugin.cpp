#include <TopPtWeightPlugin.hpp>

#include <Processor.hpp>

#include <string>
#include <sstream>
#include <stdexcept>

#include <iostream>  // Delete me!


using namespace std;


TopPtWeightPlugin::TopPtWeightPlugin():
    Plugin("TopPtWeight")
{}


Plugin *TopPtWeightPlugin::Clone() const
{
    return new TopPtWeightPlugin();
}


void TopPtWeightPlugin::BeginRun(Dataset const &dataset)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    
    isTTbar = dataset.TestProcess(Dataset::Process::ttbar);
}


bool TopPtWeightPlugin::ProcessEvent()
{
    // First check if it is ttbar at all
    if (not isTTbar)
    {
        weight = 1.;
        return true;
    }
    
    
    unsigned nTopFound = 0;
    unsigned nLeptonFound = 0;  // to discriminate between semileptonic and dileptonic ttbar
    double topPt[2];
    
    // Loop over particles in the hard interaction
    for (auto const &p: (*reader)->GetHardGenParticles())
    {
        int const absPdgId = abs(p.GetPdgId());
        
        if (absPdgId == 6)
        {
            ++nTopFound;
            
            if (nTopFound == 3)
                throw runtime_error("TopPtWeightPlugin::ProcessEvent: Found more than two top "
                 "quarks in an event.");
            
            topPt[nTopFound - 1] = p.Pt();
        }
        
        if (absPdgId == 11 or absPdgId == 13 or absPdgId == 15)  // charged leptons
            ++nLeptonFound;
    }
    
    
    if (nTopFound < 2)
    {
        cout << "TopPtWeightPlugin: Found less than two top quarks. The weight will be set to 1.\n";
        //^ This is going to screw up the output!
        
        weight = 1.;
        return true;
    }
    
    if (nLeptonFound > 2)
        throw runtime_error("TopPtWeightPlugin::ProcessEvent: Found more than three leptons "
                 "in an event.");
    
    
    // Calculate the actual event weight
    if (nLeptonFound == 2)
        weight = CalculateTopPtWeight(topPt[0], topPt[1], 0.148, -0.00129);
    else if (nLeptonFound == 1)
        weight = CalculateTopPtWeight(topPt[0], topPt[1], 0.159, -0.00141);
    else
        weight = CalculateTopPtWeight(topPt[0], topPt[1], 0.156, -0.00137);
        //^ This is actually a combination of semileptonic and dileptonic tunning. But these is no
        //measurement for all-hadronic
    
    
    /// The plugin does not perform any filtering
    return true;
}


double TopPtWeightPlugin::GetWeight() const
{
    return weight;
}


double TopPtWeightPlugin::CalculateTopPtWeight(double pt1, double pt2, double a, double b)
{
    return sqrt(exp(a + b * pt1) * exp(a + b * pt2));
}

