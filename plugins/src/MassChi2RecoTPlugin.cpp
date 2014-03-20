#include <MassChi2RecoTPlugin.hpp>

#include <Processor.hpp>
#include <BTagSFInterface.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>


using namespace std;


MassChi2RecoTPlugin::MassChi2RecoTPlugin(string const &name_):
    RecoTPlugin(name_)
{}


Plugin *MassChi2RecoTPlugin::Clone() const
{
    return new MassChi2RecoTPlugin(name);
}


void MassChi2RecoTPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


bool MassChi2RecoTPlugin::ProcessEvent()
{
    auto const &jets = (*reader)->GetJets();
    
    
    // Make sure the event contains a sufficient number of jets
    if (jets.size() < 2)
        throw logic_error("MassChi2RecoTPlugin::ProcessEvent: Cannot reconstruct an event because "
         "it contains less than two jets.");
    
    
    // Reconstruct W boson
    Candidate const wBoson((*reader)->GetLeptons().front().P4() + (*reader)->GetNeutrino().P4());
    
    
    // Loop over all the possible ways to choose a central jet. Only a central jet can be assigned
    //to the top
    iBTop = -1;
    double minChi2 = numeric_limits<double>::infinity();
    
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        // Skip non-central jets
        if (fabs(jets.at(i).Eta()) > BTagSFInterface::GetMaxPseudorapidity())
            continue;
        
        double const chi2 = pow(((wBoson.P4() + jets.at(i).P4()).M() - 176.5) / 37.1, 2);
        //^ Mass mean and resoulution is cited according to
        //~aapopov/workspace/tHq/2012Bravo/2013.11.28_Mass-chi2-reconstruction/info.txt
        
        if (chi2 < minChi2)
        {
            iBTop = i;
            minChi2 = chi2;
        }
    }
    
    
    // A sanity check
    if (iBTop == unsigned(-1))
    //^ This could have happened if only there were no central jets in the event
        throw logic_error("MassChi2RecoTPlugin::ProcessEvent: Cannot reconstruct an event because "
         "it does not contain central jets.");
    
    
    // Find the recoil jet. It is chosen as the most forward jet excluding the jet already used to
    //reconstruct the top quark
    iRecoilJet = -1;
    double maxEta = -numeric_limits<double>::infinity();
    
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        if (i == iBTop)
            continue;
        
        auto const &jet = jets.at(i);
        
        if (fabs(jet.Eta()) > maxEta)
        {
            iRecoilJet = i;
            maxEta = fabs(jet.Eta());
        }
    }
    
    // Because an event is required to contain at least two jets, a recoil jet is always found
    
    
    return true;
}


Jet const &MassChi2RecoTPlugin::GetRecoilJet() const
{
    return (*reader)->GetJets().at(iRecoilJet);
}


Jet const &MassChi2RecoTPlugin::GetBTop() const
{
    return (*reader)->GetJets().at(iBTop);
}
