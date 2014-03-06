#include <MassChi2RecoTTbarPlugin.hpp>

#include <Processor.hpp>
#include <BTagSFInterface.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>


using namespace std;


MassChi2RecoTTbarPlugin::MassChi2RecoTTbarPlugin(string const &name_):
    RecoTTbarPlugin(name_)
{}


Plugin *MassChi2RecoTTbarPlugin::Clone() const
{
    return new MassChi2RecoTTbarPlugin(name);
}


void MassChi2RecoTTbarPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


bool MassChi2RecoTTbarPlugin::ProcessEvent()
{
    auto const &jets = (*reader)->GetJets();
    
    
    // Make sure the event contains a sufficient number of jets
    if (jets.size() < 4)
        throw logic_error("MassChi2RecoTTbarPlugin::ProcessEvent: Cannot reconstruct an event "
         "because it contains less than four jets.");
    
    
    // Reconstruct W boson that decays leptonically
    Candidate const wBoson((*reader)->GetLeptons().front().P4() + (*reader)->GetNeutrino().P4());
    
    
    // Only central jets are matched to the b quarks. Put their indices into a vector
    centralJetIndices.clear();
    
    for (unsigned i = 0; i < jets.size(); ++i)
        if (fabs(jets.at(i).Eta()) < BTagSFInterface::GetMaxPseudorapidity())
            centralJetIndices.push_back(i);
    
    if (centralJetIndices.size() < 2)
        throw logic_error("MassChi2RecoTTbarPlugin::ProcessEvent: Cannot reconstruct an event "
         "because it contains less than two central jets.");
    
    
    // Consider all the possible jet combinations
    iBTopLep = iBTopHad = iL1TopHad = iL2TopHad = -1;
    double minChi2 = numeric_limits<double>::infinity();
    
    // Loop over all the ways to choose the jet from semileptonic decay of a top quark
    for (unsigned icBLep = 0; icBLep < centralJetIndices.size(); ++icBLep)
    {
        // Precalculate the part of the chi2 that depends on the top quark that decays
        //semileptonically
        double const chi2TopLep =
         pow(((wBoson.P4() + jets.at(centralJetIndices.at(icBLep)).P4()).M() - 178.0) / 35.1, 2);
        //^ Mass mean and resoulution is cited according to
        //~aapopov/workspace/tHq/2012Bravo/2013.11.28_Mass-chi2-reconstruction/info.txt
        
        
        // Loop over all the ways to choose the b-jet from hadronically decaying top quark
        for (unsigned icBHad = 0; icBHad < centralJetIndices.size(); ++icBHad)
        {
            if (icBHad == icBLep)
                continue;
            
            
            // Now consider all the ways to choose the light-flavour quarks from decay of a W boson.
            //When a light-flavour jet happens to be central, just combination is symmetric against
            //a swap of this jet and the b-jet that corresponds to index icBHad. But such symmetry
            //is not expoited here, and the corresponding combination is evaluated several times.
            for (unsigned iq1 = 0; iq1 < jets.size() - 1; ++iq1)
            {
                if (iq1 == centralJetIndices.at(icBLep) or iq1 == centralJetIndices.at(icBHad))
                    continue;
                
                for (unsigned iq2 = iq1 + 1; iq2 < jets.size(); ++iq2)
                {
                    if (iq2 == centralJetIndices.at(icBLep) or iq2 == centralJetIndices.at(icBHad))
                        continue;
                    
                    
                    double const massTopHad = (jets.at(centralJetIndices.at(icBHad)).P4() +
                     jets.at(iq1).P4() + jets.at(iq2).P4()).M();
                    double const chi2 = chi2TopLep + pow((massTopHad - 175.9) / 30.3, 2);
                    //^ Mass mean and resoulution is cited according to
                    //~aapopov/workspace/tHq/2012Bravo/2013.11.28_Mass-chi2-reconstruction/info.txt
                    
                    if (chi2 < minChi2)
                    {
                        iBTopLep = centralJetIndices.at(icBLep);
                        iBTopHad = centralJetIndices.at(icBHad);
                        iL1TopHad = iq1;
                        iL2TopHad = iq2;
                        
                        minChi2 = chi2;
                    }
                }
            }
        }
    }
    
    return true;
}


Jet const &MassChi2RecoTTbarPlugin::GetBTopLep() const
{
    return (*reader)->GetJets().at(iBTopLep);
}


Jet const &MassChi2RecoTTbarPlugin::GetBTopHad() const
{
    return (*reader)->GetJets().at(iBTopHad);
}


pair<Jet const &, Jet const &> MassChi2RecoTTbarPlugin::GetLightTopHad() const
{
    auto const &jets = (*reader)->GetJets();
    
    return {jets.at(iL1TopHad), jets.at(iL2TopHad)};
}
