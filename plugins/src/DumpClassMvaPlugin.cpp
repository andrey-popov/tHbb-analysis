#include <DumpClassMvaPlugin.hpp>

#include <Processor.hpp>
#include <ROOTLock.hpp>
#include <Logger.hpp>

#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <boost/filesystem.hpp>


using namespace std;


DumpClassMvaPlugin::DumpClassMvaPlugin(string const &outDirectory_, string const &filePostfix_,
 bool saveSystWeights_):
    Plugin("DumpClassMva"),
    mva("Silent"),
    outDirectory(outDirectory_), filePostfix(filePostfix_),
    saveSystWeights(saveSystWeights_)
{
    // Make sure the directory path ends with a slash
    if (outDirectory.back() != '/')
        outDirectory += '/';
    
    // Create the output directory if it does not exist
    boost::filesystem::create_directories(outDirectory);
    
    
    // Specify input variables and MVA for event reconstruction
    mva.AddVariable("glb_Charge_Lep", &glb_Charge_Lep);
    mva.AddVariable("log(glb_SqrtSHat)", &glb_LogSqrtSHat);
    mva.AddVariable("log(glb_Pt_J1)", &glb_LogPt_J1);
    mva.AddVariable("log(glb_MET)", &glb_LogMET);
    mva.AddVariable("log(thq_Pt_Higgs)", &thq_LogPt_Higgs);
    mva.AddVariable("log(thq_Pt_Recoil)", &thq_LogPt_Recoil);
    mva.AddVariable("thq_Cos_LepRecoil_TH", &thq_Cos_LepRecoil_TH);
    mva.AddVariable("abs(thq_Eta_Recoil)", &thq_AbsEta_Recoil);
    mva.AddVariable("glb_Sphericity", &glb_Sphericity);
    mva.AddVariable("tt_MaxEta_Light", &tt_MaxEta_Light);
    mva.AddVariable("log(tt_Mass_TopHad)", &tt_LogMass_TopHad);
    mva.AddVariable("log(tt_DMass_TopHadWHad)", &tt_LogDMass_TopHadWHad);
    mva.AddVariable("log(tt_MaxMass_BTopHadLight)", &tt_LogMaxMass_BTopHadLight);
    mva.AddVariable("thq_NumBTag_Higgs", &thq_NumBTag_Higgs);
    mva.AddVariable("tt_NumPassBTag_Light", &tt_NumPassBTag_Light);
    mva.AddVariable("thq_DeltaR_BJetsHiggs", &thq_DeltaR_BJetsHiggs);
    mva.AddVariable("tt_DeltaR_Light", &tt_DeltaR_Light);
    mva.AddVariable("tt_SumCharge_Light", &tt_SumCharge_Light);
    mva.AddVariable("tt_Charge_BTopHad", &tt_Charge_BTopHad);
    
    mva.BookMVA("Default", "/afs/cern.ch/user/a/aapopov/workspace/tHq/2012Bravo/"
     "2014.02.24_Full-analysis/Step1_MVA/class/Train/weights/"
     "THvsBkg_Global_THvsBkg_Global_BFGS.weights.xml");
    
    
    // Set the number of weights in an event
    if (saveSystWeights)
        bfNumWeights = 7;
    else
        bfNumWeights = 1;
}


Plugin *DumpClassMvaPlugin::Clone() const
{
    return new DumpClassMvaPlugin(outDirectory, filePostfix, saveSystWeights);
}


void DumpClassMvaPlugin::BeginRun(Dataset const &dataset)
{
    // Save pointers to required plugins for convenience
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
    builderTopHiggs =
     dynamic_cast<RecoTHPlugin const *>(processor->GetPluginBefore("RecoTH", name));
    builderTTbar =
     dynamic_cast<RecoTTbarPlugin const *>(processor->GetPluginBefore("RecoTTbar", name));
    dataDrivenTTbar = dynamic_cast<JetTagDataDrivenPlugin const *>(
     processor->GetPluginBeforeQuiet("DataDrivenTTbar", name));
    
    // Save pointer to top-pt reweighter
    topPtReweighter = dynamic_cast<TopPtWeightPlugin const *>(
     processor->GetPluginBeforeQuiet("TopPtWeight", name));
    
    
    // Creation of ROOT objects is not thread-safe and must be protected
    ROOTLock::Lock();
    
    // Create the output file
    string fileName(dataset.GetFiles().front().GetBaseName());
    //fileName = fileName.substr(0, fileName.find_first_of('_'));
    string fullFileName = outDirectory + fileName;
    
    if (filePostfix.length() > 0)
        fullFileName += "_" + filePostfix;
    
    fullFileName += ".root";
    
    file.reset(new TFile(fullFileName.c_str(), "recreate"));
    
    // Create the tree
    tree.reset(new TTree("Vars", "Response of the classification MVA"));
    
    // End of critical block
    ROOTLock::Unlock();
    
    
    // Assign branch addresses
    tree->Branch("EventNumber", &bfEventNumber);
    tree->Branch("RunNumber", &bfRunNumber);
    tree->Branch("LumiSection", &bfLumiSection);
    
    tree->Branch("NumTag20", &bfNumTag20);
    
    tree->Branch("NumWeights", &bfNumWeights);
    tree->Branch("Weights", bfWeights, "Weights[NumWeights]/F");
    
    tree->Branch("MvaResponse", &bfMvaResponse);
}


void DumpClassMvaPlugin::EndRun()
{
    // Operations with ROOT objects performed here are not thread-safe and must be guarded
    ROOTLock::Lock();
    
    // Write the tree and close the file
    file->cd();
    tree->Write("", TObject::kOverwrite);
    
    // Delete the objects
    tree.reset();
    file.reset();
    
    ROOTLock::Unlock();
}


bool DumpClassMvaPlugin::ProcessEvent()
{
    // Save event ID
    auto const &eventID = (*reader)->GetEventID();
    
    bfRunNumber = eventID.Run();
    bfEventNumber = eventID.Event();
    bfLumiSection = eventID.LumiBlock();
    
    
    // Save event weight(s)
    bfWeights[0] = (*reader)->GetCentralWeight();
    
    if (saveSystWeights)
    {
        auto const &puWeights = (*reader)->GetSystWeight(SystTypeWeight::PileUp);
        bfWeights[1] = puWeights.front().up;
        bfWeights[2] = puWeights.front().down;
        
        auto const &tagWeights = (*reader)->GetSystWeight(SystTypeWeight::TagRate);
        bfWeights[3] = tagWeights.front().up;
        bfWeights[4] = tagWeights.front().down;
        
        auto const &mistagWeights = (*reader)->GetSystWeight(SystTypeWeight::MistagRate);
        bfWeights[5] = mistagWeights.front().up;
        bfWeights[6] = mistagWeights.front().down;
    }
    
    if (topPtReweighter)
    {
        double const topPtWeight = topPtReweighter->GetWeight();
        
        for (int i = 0; i < bfNumWeights; ++i)
            bfWeights[i] *= topPtWeight;
    }
    
    
    
    // Define some shortcuts and reconstruct W boson
    auto const &jets = (*reader)->GetJets();
    auto const &lepton = (*reader)->GetLeptons().front();
    
    Candidate const wLep(lepton.P4() + (*reader)->GetNeutrino().P4());
    
    
    // Count b-tagged jets
    bfNumTag20 = 0;
    
    for (auto const &j: jets)
        if (j.CSV() > 0.898)
            ++bfNumTag20;
    
    
    // Calculate observables that do not require event reconstruction
    glb_Charge_Lep = lepton.Charge();
    glb_LogPt_J1 = log(jets.at(0).Pt());
    //glb_Pt_J2 = jets.at(1).Pt();
    
    TLorentzVector p4AllJets;
    //glb_Ht = 0.;
    
    for (auto const &j: jets)
    {
        p4AllJets += j.P4();
        //glb_Ht += j.Pt();
    }
    
    glb_LogMET = log((*reader)->GetMET().Pt());
    //glb_Ht += lepton.Pt() + glb_MET;
    glb_LogSqrtSHat = log((p4AllJets + wLep.P4()).M());
    
    
    // Calculate sphericity
    TMatrixDSym sphericityTensor(3);
    double norm = 0.;
    
    for (auto const &p: {lepton.P4(), (*reader)->GetNeutrino().P4()})
    {
        TVector3 p3(p.Vect());
        norm += p3.Mag2();
        
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                sphericityTensor(i, j) += p3[i] * p3[j];
    }
    
    for (auto const &jet: jets)
    {
        TVector3 p3(jet.P4().Vect());
        norm += p3.Mag2();
        
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                sphericityTensor(i, j) += p3[i] * p3[j];
    }
    
    sphericityTensor *= 1. / norm;
    
    TMatrixDSymEigen eigenValCalc(sphericityTensor);
    TVectorD eigenVals(eigenValCalc.GetEigenValues());
    
    glb_Sphericity = 1.5 * (eigenVals[1] + eigenVals[2]);
    
    
    // Calculate variables reconstructed under thq hypothesis
    Jet const &bTop = builderTopHiggs->GetBTop();
    Jet const &b1Higgs = builderTopHiggs->GetBHiggs().first;
    Jet const &b2Higgs = builderTopHiggs->GetBHiggs().second;
    
    Candidate const top(bTop.P4() + wLep.P4());
    Candidate const higgs(b1Higgs.P4() + b2Higgs.P4());
    Jet const &recoilQuark = builderTopHiggs->GetRecoilJet();
    
    // thq_Mass_Top = top.M();
    // thq_Pt_Top = top.Pt();
    // thq_Eta_Top = top.Eta();
    
    // thq_Pt_BTop = bTop.Pt();
    // thq_Eta_BTop = bTop.Eta();
    // thq_PassBTag_BTop = 0. + (bTop.CSV() > 0.898);
    // thq_Charge_BTop = bTop.Charge() * lepton.Charge();
    
    // thq_DeltaR_BTopW = bTop.P4().DeltaR(wLep.P4());
    // thq_DEta_TopLep = fabs(top.Eta() - lepton.Eta());
    // thq_Mass_BTopLep = (bTop.P4() + lepton.P4()).M();
    
    // thq_Mass_Higgs = higgs.M();
    thq_LogPt_Higgs = log(higgs.Pt());
    // thq_Eta_Higgs = higgs.Eta();
    
    thq_NumBTag_Higgs = 0. + (b1Higgs.CSV() > 0.898) + (b2Higgs.CSV() > 0.898);
    // thq_AbsCharge_Higgs = fabs(b1Higgs.Charge() + b2Higgs.Charge());
    // thq_MinPt_BHiggs = min(b1Higgs.Pt(), b2Higgs.Pt());
    // thq_MaxEta_BHiggs = max(fabs(b1Higgs.Eta()), fabs(b2Higgs.Eta()));

    
    thq_LogPt_Recoil = log(recoilQuark.Pt());
    thq_AbsEta_Recoil = abs(recoilQuark.Eta());
    
    // thq_DeltaR_TopHiggs = higgs.P4().DeltaR(top.P4());
    thq_DeltaR_BJetsHiggs = b1Higgs.P4().DeltaR(b2Higgs.P4());
    
    // thq_Mass_TopHiggs = (higgs.P4() + top.P4()).M();
    // thq_RelHt = (higgs.Pt() + top.Pt() + recoilQuark.Pt()) / glb_Ht;
    
    
    // Calculate thq_Cos_LepRecoil_TH
    TVector3 b((higgs.P4() + top.P4()).BoostVector());
    
    TLorentzVector boostedLepton(lepton.P4());
    boostedLepton.Boost(-b);
    TVector3 const p3Lepton(boostedLepton.Vect());
    
    TLorentzVector boostedRecoil(recoilQuark.P4());
    boostedRecoil.Boost(-b);
    TVector3 const p3Recoil(boostedRecoil.Vect());
    
    thq_Cos_LepRecoil_TH = p3Lepton.Dot(p3Recoil) / (p3Lepton.Mag() * p3Recoil.Mag());
    
    
    // Calculate observables constructed under ttbar hypothesis
    Jet const &bTopLep = builderTTbar->GetBTopLep();
    Jet const &bTopHad = builderTTbar->GetBTopHad();
    Jet const &q1TopHad = builderTTbar->GetLightTopHad().first;
    Jet const &q2TopHad = builderTTbar->GetLightTopHad().second;
    
    Candidate const topLep(bTopLep.P4() + wLep.P4());
    Candidate const wHad(q1TopHad.P4() + q2TopHad.P4());
    Candidate const topHad(bTopHad.P4() + wHad.P4());
    
    // tt_Mass_TopLep = topLep.M();
    // tt_Pt_TopLep = topLep.Pt();
    // tt_Eta_TopLep = topLep.Eta();
    
    // tt_Pt_BTopLep = bTopLep.Pt();
    // tt_Eta_BTopLep = bTopLep.Eta();
    // tt_PassBTag_BTopLep = 0. + (bTopLep.CSV() > 0.898);
    // tt_Charge_BTopLep = bTopLep.Charge() * lepton.Charge();
    
    // tt_DeltaR_BTopLepWLep = bTopLep.P4().DeltaR(wLep.P4());
    // tt_DEta_TopLepLep = fabs(topLep.Eta() - lepton.Eta());
    // tt_Mass_BTopLepLep = (bTopLep.P4() + lepton.P4()).M();
    
    tt_LogMass_TopHad = log(topHad.M());
    // tt_Pt_TopHad = topHad.Pt();
    // tt_Eta_TopHad = topHad.Eta();
    
    // tt_Pt_BTopHad = bTopHad.Pt();
    // tt_Eta_BTopHad = bTopHad.Eta();
    // tt_PassBTag_BTopHad = 0. + (bTopHad.CSV() > 0.898);
    tt_Charge_BTopHad = bTopHad.Charge() * lepton.Charge();
    
    // tt_DeltaR_BTopHadWHad = bTopHad.P4().DeltaR(wHad.P4());
    
    // tt_Mass_WHad = wHad.M();
    // tt_Pt_WHad = wHad.Pt();
    // tt_Eta_WHad = wHad.Eta();
    
    tt_LogDMass_TopHadWHad = log(topHad.M() - wHad.M());
    // tt_RelHt = (topHad.Pt() + topLep.Pt()) / glb_Ht;
    
    tt_DeltaR_Light = q1TopHad.P4().DeltaR(q2TopHad.P4());
    tt_NumPassBTag_Light = 0. + (q1TopHad.CSV() > 0.898) + (q2TopHad.CSV() > 0.898);
    tt_MaxEta_Light = max(fabs(q1TopHad.Eta()), fabs(q2TopHad.Eta()));
    tt_SumCharge_Light = (q1TopHad.Charge() + q2TopHad.Charge()) * lepton.Charge();
    
    tt_LogMaxMass_BTopHadLight =
     log(max((bTopHad.P4() + q1TopHad.P4()).M(), (bTopHad.P4() + q2TopHad.P4()).M()));
    
    
    
    // Evaluate the MVA
    bfMvaResponse = mva.EvaluateMVA("Default");
    
    
    // Finally, fill the tree
    tree->Fill();
    
    return true;
}
