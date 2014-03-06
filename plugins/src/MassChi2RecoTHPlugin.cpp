#include <MassChi2RecoTHPlugin.hpp>

#include <Processor.hpp>
#include <BTagSFInterface.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>


using namespace std;


MassChi2RecoTHPlugin::MassChi2RecoTHPlugin(string const &name_):
    RecoTHPlugin(name_)
{}


Plugin *MassChi2RecoTHPlugin::Clone() const
{
    return new MassChi2RecoTHPlugin(name);
}


void MassChi2RecoTHPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


bool MassChi2RecoTHPlugin::ProcessEvent()
{
    auto const &jets = (*reader)->GetJets();
    
    
    // Make sure the event contains a sufficient number of jets
    if (jets.size() < 4)
        throw logic_error("MassChi2RecoTHPlugin::ProcessEvent: Cannot reconstruct an event because "
         "it contains less than four jets.");
    
    
    // Reconstruct W boson
    Candidate const wBoson((*reader)->GetLeptons().front().P4() + (*reader)->GetNeutrino().P4());
    
    
    // Find central jets. Only they are assigned to the top quark of the Higgs boson
    centralJetIndices.clear();
    
    for (unsigned i = 0; i < jets.size(); ++i)
        if (fabs(jets.at(i).Eta()) < BTagSFInterface::GetMaxPseudorapidity())
            centralJetIndices.push_back(i);
    
    
    if (centralJetIndices.size() < 3)
        throw logic_error("MassChi2RecoTHPlugin::ProcessEvent: Cannot reconstruct an event because "
         "it contains less than three central jets.");
    
    
    // Loop over all the possible ways to choose three central jets
    iBTop = iB1Higgs = iB2Higgs = -1;
    double minChi2 = numeric_limits<double>::infinity();
    
    // The outer loop checks all the possible ways to choose the jet from top decay
    for (unsigned it = 0; it < centralJetIndices.size(); ++it)
    {
        // Precalculate the part of the chi2 that depends on the top quark
        double const chi2Top =
         pow(((wBoson.P4() + jets.at(centralJetIndices.at(it)).P4()).M() - 176.5) / 37.1, 2);
        //^ Mass mean and resoulution is cited according to
        //~aapopov/workspace/tHq/2012Bravo/2013.11.28_Mass-chi2-reconstruction/info.txt
        
        
        // Loop over all the ways to choose the first b from higgs for each value of index it
        for (unsigned ih1 = 0; ih1 < centralJetIndices.size() - 1; ++ih1)
        {
            if (ih1 == it)
                continue;
            
            // Loop over all ways to choose the second b-jet from decay of higgs. Can assume that
            //this jet is softer than ih1
            for (unsigned ih2 = ih1 + 1; ih2 < centralJetIndices.size(); ++ih2)
            {
                if (ih2 == it)
                    continue;
                
                
                // Try the current combination of jets
                double const massHiggs = (jets.at(centralJetIndices.at(ih1)).P4() +
                 jets.at(centralJetIndices.at(ih2)).P4()).M();
                double const chi2 = chi2Top + pow((massHiggs - 113.5) / 20.6 , 2);
                //^ Mass mean and resoulution is cited according to
                //~aapopov/workspace/tHq/2012Bravo/2013.11.28_Mass-chi2-reconstruction/info.txt
                
                if (chi2 < minChi2)
                {
                    iBTop = centralJetIndices.at(it);
                    iB1Higgs = centralJetIndices.at(ih1);
                    iB2Higgs = centralJetIndices.at(ih2);
                    
                    minChi2 = chi2;
                }
            }
        }
    }
    
    
    // Find the recoil jet. It is chosen as the most forward jet excluding the jets already used to
    //reconstruct the top quark of the Higgs boson
    iRecoilJet = -1;
    double maxEta = -numeric_limits<double>::infinity();
    
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        if (i == iBTop or i == iB1Higgs or i == iB2Higgs)
            continue;
        
        auto const &jet = jets.at(i);
        
        if (fabs(jet.Eta()) > maxEta)
        {
            iRecoilJet = i;
            maxEta = fabs(jet.Eta());
        }
    }
    
    // Because an event is required to contain at least four jets, a recoil jet is always found
    
    
    return true;
}


Jet const &MassChi2RecoTHPlugin::GetRecoilJet() const
{
    return (*reader)->GetJets().at(iRecoilJet);
}


Jet const &MassChi2RecoTHPlugin::GetBTop() const
{
    return (*reader)->GetJets().at(iBTop);
}


pair<Jet const &, Jet const &> MassChi2RecoTHPlugin::GetBHiggs() const
{
    auto const &jets = (*reader)->GetJets();
    
    return {jets.at(iB1Higgs), jets.at(iB2Higgs)};
}
