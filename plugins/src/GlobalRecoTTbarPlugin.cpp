#include <GlobalRecoTTbarPlugin.hpp>

#include <Processor.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <array>


using namespace std;


GlobalRecoTTbarPlugin::GlobalRecoTTbarPlugin(string const &name_ /*= string("")*/):
    RecoTTbarPlugin((name_.length() == 0) ? "GlobalRecoTTbar" : name_),
    mvaReco("Silent")
{
    // Specify input variables and MVA for event reconstruction
    mvaReco.AddVariable("log(Mass_TopHad - Mass_WHad)", &LogDMass_TopHadWHad);
    mvaReco.AddVariable("log(Mass_WHad)", &LogMass_WHad);
    mvaReco.AddVariable("log(Mass_BTopLepLep)", &LogMass_BTopLepLep);
    mvaReco.AddVariable("abs(Eta_TopHad)", &AbsEta_TopHad);
    mvaReco.AddVariable("PassBTag_BTopLep", &PassBTag_BTopLep);
    mvaReco.AddVariable("PassBTag_BTopHad", &PassBTag_BTopHad);
    mvaReco.AddVariable("NumBTag_Light", &NumBTag_Light);
    mvaReco.AddVariable("DeltaR_Light", &DeltaR_Light);
    mvaReco.AddVariable("DeltaR_BTopHadWHad", &DeltaR_BTopHadWHad);
    mvaReco.AddVariable("DeltaR_BTopLepWLep", &DeltaR_BTopLepWLep);
    mvaReco.AddVariable("RelHt", &RelHt);
    mvaReco.AddVariable("log(Pt_TopLep)", &LogPt_TopLep);
    mvaReco.AddVariable("log(Pt_TopHad)", &LogPt_TopHad);
    mvaReco.AddVariable("Charge_BTopHad - Charge_BTopLep", &DCharge_BTopHadBTopLep);
    mvaReco.AddVariable("SumCharge_Light", &SumCharge_Light);
    
    mvaReco.BookMVA("Default", "/afs/cern.ch/user/a/aapopov/workspace/tHq/2012Bravo/"
     "2014.02.24_Full-analysis/Step1_MVA/ttbar/Train/weights/"
     "GlobalTTbarReco_GlobalTTbarReco_BFGS.weights.xml");
}


Plugin *GlobalRecoTTbarPlugin::Clone() const
{
    return new GlobalRecoTTbarPlugin(name);
}


void GlobalRecoTTbarPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


bool GlobalRecoTTbarPlugin::ProcessEvent()
{
    auto const &jets = (*reader)->GetJets();
    
    
    // Loop over all possible ways to classify jets in the event and find the one with largest MVA
    //response
    iBTopLep = iBTopHad = iL1TopHad = iL2TopHad = -1;
    double highestRank = -numeric_limits<double>::infinity();
    
    for (unsigned i1 = 0; i1 < jets.size() - 3; ++i1)
        for (unsigned i2 = i1 + 1; i2 < jets.size() - 2; ++i2)
            for (unsigned i3 = i2 + 1; i3 < jets.size() - 1; ++i3)
                for (unsigned i4 = i3 + 1; i4 < jets.size(); ++i4)
                {
                    // Find central jets
                    array<unsigned, 4> centralJetIndices;
                    unsigned nCentralJets = 0;
                    
                    for (unsigned const &idx: {i1, i2, i3, i4})
                        if (fabs(jets.at(idx).Eta()) < 2.4)
                        {
                            centralJetIndices.at(nCentralJets) = idx;
                            ++nCentralJets;
                        }
                    
                    if (nCentralJets < 2)
                        continue;
                    
                    
                    // Consider all the ways to choose two b-jets out of the central ones
                    for (unsigned icBTopLep = 0; icBTopLep < nCentralJets; ++icBTopLep)
                        for (unsigned icBTopHad = 0; icBTopHad < nCentralJets; ++icBTopHad)
                        {
                            if (icBTopHad == icBTopLep)
                                continue;
                            
                            unsigned const iBTopLepProposed = centralJetIndices.at(icBTopLep);
                            unsigned const iBTopHadProposed = centralJetIndices.at(icBTopHad);
                            
                            
                            // The two chosen b-jets define the light-flavour jets in the chosen
                            //four uniquely
                            unsigned iL1TopHadProposed = -1, iL2TopHadProposed = -1;
                            
                            for (unsigned const &idx: {i1, i2, i3, i4})
                            {
                                if (idx == iBTopLepProposed or idx == iBTopHadProposed)
                                    continue;
                                
                                iL1TopHadProposed = idx;
                                break;
                            }
                            
                            for (unsigned const &idx: {i1, i2, i3, i4})
                            {
                                if (idx == iBTopLepProposed or idx == iBTopHadProposed or
                                 idx == iL1TopHadProposed)
                                    continue;
                                
                                iL2TopHadProposed = idx;
                                break;
                            }
                            
                            // Note that i1 < ... < i4; hence, the second light-flavour jet is
                            //softer than the first one by construction
                            
                            
                            CalculateVariables(jets.at(iBTopLepProposed), jets.at(iBTopHadProposed),
                             jets.at(iL1TopHadProposed), jets.at(iL2TopHadProposed));
                            double const rank = mvaReco.EvaluateMVA("Default");
                            
                            if (rank > highestRank)
                            {
                                iBTopLep = iBTopLepProposed;
                                iBTopHad = iBTopHadProposed;
                                iL1TopHad = iL1TopHadProposed;
                                iL2TopHad = iL2TopHadProposed;
                                
                                highestRank = rank;
                            }
                        }
                }
    
    
    // Sanity check
    if (iBTopLep >= jets.size())
        throw logic_error("GlobalRecoTTbarPlugin::ProcessEvent: Failed to reconstruct current "
         "event.");
    
    
    return true;
}


Jet const &GlobalRecoTTbarPlugin::GetBTopLep() const
{
    return (*reader)->GetJets().at(iBTopLep);
}


Jet const &GlobalRecoTTbarPlugin::GetBTopHad() const
{
    return (*reader)->GetJets().at(iBTopHad);
}


pair<Jet const &, Jet const &> GlobalRecoTTbarPlugin::GetLightTopHad() const
{
    auto const &jets = (*reader)->GetJets();
    
    return {jets.at(iL1TopHad), jets.at(iL2TopHad)};
}


void GlobalRecoTTbarPlugin::CalculateVariables(Jet const &bTopLep, Jet const &bTopHad,
 Jet const &q1TopHad, Jet const &q2TopHad)
{
    // Calculate variables related to the top quark decaying semileptonically
    auto const &lepton = (*reader)->GetLeptons().front();
    Candidate const wLep(lepton.P4() + (*reader)->GetNeutrino().P4());
    Candidate const topLep(wLep.P4() + bTopLep.P4());
    
    //bfMass_TopLep = topLep.M();
    LogPt_TopLep = log(topLep.Pt());
    //bfEta_TopLep = topLep.Eta();
    
    //bfPt_BTopLep = bTopLep.Pt();
    //bfEta_BTopLep = bTopLep.Eta();
    PassBTag_BTopLep = 0. + (bTopLep.CSV() > 0.898);
    //bfCharge_BTopLep = bTopLep.Charge() * lepton.Charge();
    
    LogMass_BTopLepLep = log((bTopLep.P4() + lepton.P4()).M());
    DeltaR_BTopLepWLep = bTopLep.P4().DeltaR(wLep.P4());
    //bfDEta_TopLepLep = fabs(topLep.Eta() - lepton.Eta());
    
    // Calculate the cosine
    // TVector3 b((wLep.P4()).BoostVector());
    
    // TLorentzVector boostedLepton(lepton.P4());
    // boostedLepton.Boost(-b);
    // TVector3 p3Lepton(boostedLepton.Vect());
    
    // TLorentzVector boostedBJet(bTopLep.P4());
    // boostedBJet.Boost(-b);
    // TVector3 const p3BJet(boostedBJet.Vect());
    
    // bfCos_LepBTopLep_WLep = p3Lepton.Dot(p3BJet) / (p3Lepton.Mag() * p3BJet.Mag());
    
    
    // Calculate variables related to the top quark decaying semileptonically
    Candidate const wHad(q1TopHad.P4() + q2TopHad.P4());
    Candidate const topHad(wHad.P4() + bTopHad.P4());
    
    //bfMass_TopHad = topHad.M();
    LogPt_TopHad = log(topHad.Pt());
    AbsEta_TopHad = fabs(topHad.Eta());
    
    LogMass_WHad = log(wHad.M());
    //bfPt_WHad = wHad.Pt();
    //bfEta_WHad= wHad.Eta();
    
    DeltaR_BTopHadWHad = bTopHad.P4().DeltaR(wHad.P4());
    
    //bfPt_BTopHad = bTopHad.Pt();
    //bfEta_BTopHad = bTopHad.Eta();
    PassBTag_BTopHad = 0. + (bTopHad.CSV() > 0.898);
    //bfCharge_BTopHad = bTopHad.Charge() * lepton.Charge();
    
    //bfMinPt_Light = min(fabs(q1TopHad.Pt()), fabs(q2TopHad.Pt()));
    //bfMaxEta_Light = max(fabs(q1TopHad.Eta()), fabs(q2TopHad.Eta()));
    SumCharge_Light = (q1TopHad.Charge() + q2TopHad.Charge()) * lepton.Charge();
    NumBTag_Light = 0. + (q1TopHad.CSV() > 0.898) + (q2TopHad.CSV() > 0.898);
    
    DeltaR_Light = q1TopHad.P4().DeltaR(q2TopHad.P4());
    // bfMaxMass_BTopHadLight =
    //  max((bTopHad.P4() + q1TopHad.P4()).M(), (bTopHad.P4() + q2TopHad.P4()).M());
    
    LogDMass_TopHadWHad = log(topHad.M() - wHad.M());
    
    
    // Calcualate variables with correlations between different objects
    double Ht = lepton.Pt() + (*reader)->GetMET().Pt();
    
    for (auto const &j: (*reader)->GetJets())
        Ht += j.Pt();
    
    RelHt = (topLep.Pt() + topHad.Pt()) / Ht;
    //bfMass_TopLepTopHad = (topLep.P4() + topHad.P4()).M();
    //bfPt_TopLepTopHad = (topLep.P4() + topHad.P4()).Pt();
    //bfEta_TopLepTopHad = (topLep.P4() + topHad.P4()).Eta();
    //bfDeltaR_TopLepTopHad = topLep.P4().DeltaR(topHad.P4());
    
    DCharge_BTopHadBTopLep = (bTopHad.Charge() - bTopLep.Charge()) * lepton.Charge();
}
