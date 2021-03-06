#include <GlobalRecoTHPlugin.hpp>

#include <Processor.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <string>


using namespace std;


GlobalRecoTHPlugin::GlobalRecoTHPlugin(string const &name_, string const &bTagPluginName_):
    RecoTHPlugin(name_),
    bTagPluginName(bTagPluginName_), mvaReco("Silent")
{
    // Specify input variables and MVA for event reconstruction
    mvaReco.AddVariable("log(Mass_Higgs)", &LogMass_Higgs);
    mvaReco.AddVariable("log(Mass_BTopLep)", &LogMass_BTopLep);
    mvaReco.AddVariable("abs(Eta_Recoil)", &AbsEta_Recoil);
    mvaReco.AddVariable("PassBTag_BTop", &PassBTag_BTop);
    mvaReco.AddVariable("NumBTag_Higgs", &NumBTag_Higgs);
    mvaReco.AddVariable("DeltaR_BJetsHiggs", &DeltaR_BJetsHiggs);
    mvaReco.AddVariable("DeltaR_BTopW", &DeltaR_BTopW);
    mvaReco.AddVariable("DeltaR_TopHiggs", &DeltaR_TopHiggs);
    mvaReco.AddVariable("RelHt", &RelHt);
    mvaReco.AddVariable("MaxEta_BHiggs", &MaxEta_BHiggs);
    mvaReco.AddVariable("log(MinPt_BHiggs)", &LogMinPt_BHiggs);
    mvaReco.AddVariable("Charge_BTop", &Charge_BTop);
    
    
    // Read the install path from the environment
    string const path(getenv("TH_ANALYSIS_INSTALL"));
    
    mvaReco.BookMVA("Default", path + "/data/GlobalTHReco_BFGS_v3.weights.xml");
}


Plugin *GlobalRecoTHPlugin::Clone() const
{
    return new GlobalRecoTHPlugin(name, bTagPluginName);
}


void GlobalRecoTHPlugin::BeginRun(Dataset const &)
{
    // Save pointers to required plugins
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    bTagger = dynamic_cast<BTaggerPlugin const *>(processor->GetPluginBefore(bTagPluginName, name));
}


bool GlobalRecoTHPlugin::ProcessEvent()
{
    auto const &jets = (*reader)->GetJets();
    
    numInterpretations = 0;
    highestMvaResponse = -numeric_limits<double>::infinity();
    iBTop = iB1Higgs = iB2Higgs = iRecoilJet = -1;
    
    
    // Loop over all possible ways to choose four jets in the event
    for (unsigned i1 = 0; i1 < jets.size() - 3; ++i1)
        for (unsigned i2 = i1 + 1; i2 < jets.size() - 2; ++i2)
            for (unsigned i3 = i2 + 1; i3 < jets.size() - 1; ++i3)
                for (unsigned i4 = i3 + 1; i4 < jets.size(); ++i4)
                {
                    // Choose central jets out of these four
                    array<unsigned, 4> centralJetIndices;
                    unsigned nCentralJets = 0;
                    
                    for (unsigned const &index: {i1, i2, i3, i4})
                    {
                        if (fabs(jets.at(index).Eta()) < 2.4)
                        {
                            centralJetIndices.at(nCentralJets) = index;
                            ++nCentralJets;
                        }
                    }
                    
                    // Make sure there are at least three central jets
                    if (nCentralJets < 3)
                        continue;
                    
                    
                    // Loop over all the ways to choose the jet from decay of the top quark. The
                    //loop variable refers to an index in the centralJetIndices array
                    for (unsigned icBTop = 0; icBTop < nCentralJets; ++icBTop)
                    {
                        // Loop over the choice of the first jet from higgs decay
                        for (unsigned icB1Higgs = 0; icB1Higgs < nCentralJets - 1; ++icB1Higgs)
                        {
                            if (icB1Higgs == icBTop)
                                continue;
                            
                            // Loop over all the ways to choose the second jet from higgs decay.
                            //Note the jets are still ordered in pt, and this jet should be softer
                            //that icB1Higgs
                            for (unsigned icB2Higgs = icB1Higgs + 1; icB2Higgs < nCentralJets;
                             ++icB2Higgs)
                            {
                                if (icB2Higgs == icBTop)
                                    continue;
                                
                                
                                // Find index of the recoil jet
                                unsigned const iBTopProposed = centralJetIndices.at(icBTop);
                                unsigned const iB1HiggsProposed = centralJetIndices.at(icB1Higgs);
                                unsigned const iB2HiggsProposed = centralJetIndices.at(icB2Higgs);
                                unsigned iRecoilJetProposed = -1;
                                
                                for (unsigned const &index: {i1, i2, i3, i4})
                                {
                                    if (index == iBTopProposed or index == iB1HiggsProposed or
                                     index == iB2HiggsProposed)
                                        continue;
                                    
                                    iRecoilJetProposed = index;
                                    // Control will reach this line only once
                                }
                                
                                
                                // Make sure the recoil jet is not b-tagged
                                if (bTagger->IsTagged(jets.at(iRecoilJetProposed)))
                                    continue;
                                
                                
                                // Evaluate the MVA for the proposed match
                                ++numInterpretations;
                                CalculateVariables(jets.at(iBTopProposed),
                                 jets.at(iB1HiggsProposed), jets.at(iB2HiggsProposed),
                                 jets.at(iRecoilJetProposed));
                                double const response = mvaReco.EvaluateMVA("Default");
                                
                                if (response > highestMvaResponse)
                                {
                                    iBTop = iBTopProposed;
                                    iB1Higgs = iB1HiggsProposed;
                                    iB2Higgs = iB2HiggsProposed;
                                    iRecoilJet = iRecoilJetProposed;
                                    
                                    highestMvaResponse = response;
                                }
                            }
                        }
                    }
                }
    
    
    // Sanity check
    if (iBTop >= jets.size())
        throw logic_error("GlobalRecoTHPlugin::ProcessEvent: Failed to reconstruct current event.");
    
    return true;
}


Jet const &GlobalRecoTHPlugin::GetRecoilJet() const
{
    return (*reader)->GetJets().at(iRecoilJet);
}


Jet const &GlobalRecoTHPlugin::GetBTop() const
{
    return (*reader)->GetJets().at(iBTop);
}


pair<Jet const &, Jet const &> GlobalRecoTHPlugin::GetBHiggs() const
{
    auto const &jets = (*reader)->GetJets();
    
    return {jets.at(iB1Higgs), jets.at(iB2Higgs)};
}


unsigned GlobalRecoTHPlugin::GetNumInterpretations() const
{
    return numInterpretations;
}


double GlobalRecoTHPlugin::GetMvaResponse() const
{
    return highestMvaResponse;
}


void GlobalRecoTHPlugin::CalculateVariables(Jet const &bTop, Jet const &b1Higgs,
 Jet const &b2Higgs, Jet const &recoil)
{
    // Calculate variables related to the Higgs boson
    Candidate const higgs(b1Higgs.P4() + b2Higgs.P4());
    
    LogMass_Higgs = log(higgs.M());
    //bfPt_Higgs = higgs.Pt();
    //bfEta_Higgs = higgs.Eta();
    
    DeltaR_BJetsHiggs = b1Higgs.P4().DeltaR(b2Higgs.P4());
    //bfSumCharge_Higgs = fabs(b1Higgs.Charge() + b2Higgs.Charge());
    
    NumBTag_Higgs = 0 + bTagger->IsTagged(b1Higgs) + bTagger->IsTagged(b2Higgs);
    
    LogMinPt_BHiggs = log(min(b1Higgs.Pt(), b2Higgs.Pt()));
    MaxEta_BHiggs = max(fabs(b1Higgs.Eta()), fabs(b2Higgs.Eta()));
    
    
    // Calculate variables related to the top quark
    auto const &lepton = (*reader)->GetLeptons().front();
    Candidate const wBoson(lepton.P4() + (*reader)->GetNeutrino().P4());
    Candidate const top(wBoson.P4() + bTop.P4());
    
    //bfMass_Top = top.M();
    //bfPt_Top = top.Pt();
    //bfEta_Top = top.Eta();
    
    //bfPt_BTop = bTop.Pt();
    //bfEta_BTop = bTop.Eta();
    PassBTag_BTop = 0 + bTagger->IsTagged(bTop);
    Charge_BTop = bTop.Charge() * lepton.Charge();
    
    LogMass_BTopLep = log((bTop.P4() + lepton.P4()).M());
    DeltaR_BTopW = bTop.P4().DeltaR(wBoson.P4());
    //bfDEta_TopLep = fabs(top.Eta() - lepton.Eta());
    
    // Calculate the cosine
    // TVector3 b((wBoson.P4()).BoostVector());
    
    // TLorentzVector boostedLepton(lepton.P4());
    // boostedLepton.Boost(-b);
    // TVector3 p3Lepton(boostedLepton.Vect());
    
    // TLorentzVector boostedBJet(bTop.P4());
    // boostedBJet.Boost(-b);
    // TVector3 const p3BJet(boostedBJet.Vect());
    
    // bfCos_LepBTop_W = p3Lepton.Dot(p3BJet) / (p3Lepton.Mag() * p3BJet.Mag());
    
    
    // Variables related to the recoil quark
    //bfPt_Recoil = recoil.Pt();
    AbsEta_Recoil = fabs(recoil.Eta());
    //bfCharge_Recoil = recoil.Charge();
    
    
    // Calcualate variables with correlations between different objects
    double Ht = lepton.Pt() + (*reader)->GetMET().Pt();
    
    for (auto const &j: (*reader)->GetJets())
        Ht += j.Pt();
    
    RelHt = (top.Pt() + higgs.Pt() + recoil.Pt()) / Ht;
    
    DeltaR_TopHiggs = top.P4().DeltaR(higgs.P4());
}
