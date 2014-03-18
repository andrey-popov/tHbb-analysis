#include <ClassTTbarHFPlugin.hpp>

#include <Processor.hpp>

#include <array>
#include <limits>


using namespace std;


ClassTTbarHFPlugin::ClassTTbarHFPlugin(string const &name) noexcept:
    Plugin(name),
    matchConeSize(0.5)
{}


Plugin *ClassTTbarHFPlugin::Clone() const
{
    return new ClassTTbarHFPlugin(name);
}


void ClassTTbarHFPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


bool ClassTTbarHFPlugin::ProcessEvent()
{
    // Fill the collection of pointers to additional partons with appropriate partons
    additionalPartons.clear();
    
    for (ShowerParton const &p: (*reader)->GetShowerPartons())
    {
        int const absPdgId = abs(p.GetPdgId());
        
        if (absPdgId != 4 and absPdgId != 5)  // skip all but b and c quarks
            continue;
        
        if (p.GetOrigin() == ShowerParton::Origin::Proton)
            continue;
        //^ Skip immediate daughters of beam particles
        
        
        additionalPartons.push_back(&p);
    }
    
    
    // Find quarks whose descendants should be removed from the above list. They are quarks produced
    //in electroweak interactions, namely b quarks from decays of tops and c quarks from decays of
    //W bosons
    ewkPartons.clear();
    unsigned numBTopsFound = 0;
    
    for (GenParticle const &p: (*reader)->GetHardGenParticles())
    {
        if (abs(p.GetPdgId()) == 5 and abs(p.GetFirstMother()->GetPdgId()) == 6)
        //^ b quarks from decays of tops
        {
            ewkPartons.push_back(&p);
            ++numBTopsFound;
        }
        
        if (abs(p.GetPdgId()) == 4 and abs(p.GetFirstMother()->GetPdgId()) == 24)
        //^ Charms from decays of W bosons from top quarks
            ewkPartons.push_back(&p);
    }
    
    // A sanity check
    if (numBTopsFound != 2)
        throw runtime_error("ClassTTbarHFPlugin::ProcessEvent: Failed to find two b quarks from "
         "decays of the tops.");
    
    
    // Remove from the collection of additional partons those partons that are closest to the quarks
    //produced in electroweak interactions
    for (auto const &ewkParton: ewkPartons)
    {
        auto closestPartonIt = additionalPartons.end();
        double minDistance = numeric_limits<double>::infinity();
        
        for (auto pIt = additionalPartons.begin(); pIt != additionalPartons.end(); ++pIt)
        {
            if ((*pIt)->GetPdgId() != ewkParton->GetPdgId())
                continue;
            
            double const dR = (*pIt)->P4().DeltaR(ewkParton->P4());
            
            if (dR < minDistance)
            {
                closestPartonIt = pIt;
                minDistance = dR;
            }
        }
        
        // Sanity check
        if (closestPartonIt == additionalPartons.end())
            throw runtime_error("ClassTTbarHFPlugin::ProcessEvent: Cannot match a status-3 quark "
             "produced in electroweak interaction to any quark from the parton shower.");
        
        
        // Remove this parton
        additionalPartons.erase(closestPartonIt);
    }
    
    
    // Now match reconstructed jets to remaining partons
    auto const &jets = (*reader)->GetJets();
    partonCounts.clear();
    HFCounter unmatchedPartonCount{0, 0};
    
    for (unsigned i = 0; i < jets.size(); ++i)
        partonCounts.emplace_back(HFCounter{0, 0});
    
    for (auto const &p: additionalPartons)
    {
        // Find index of the closest jet
        unsigned iJetClosest = -1;
        double minDistance = numeric_limits<double>::infinity();
        
        for (unsigned iJet = 0; iJet < jets.size(); ++iJet)
        {
            double const dR = jets.at(iJet).P4().DeltaR(p->P4());
            
            if (dR < minDistance)
            {
                iJetClosest = iJet;
                minDistance = dR;
            }
        }
        
        
        // Check if the closest jet is close enough
        if (minDistance < matchConeSize)
        {
            if (abs(p->GetPdgId()) == 5)
                ++partonCounts.at(iJetClosest).nb;
            else
                ++partonCounts.at(iJetClosest).nc;
        }
        else
        {
            if (abs(p->GetPdgId()) == 5)
                ++unmatchedPartonCount.nb;
            else
                ++unmatchedPartonCount.nc;
        }
    }
    
    
    // Count different types of jets
    unsigned nJetB = 0, nJetDoubleB = 0, nJetC = 0, nJetDoubleC = 0;
    
    for (auto const &c: partonCounts)
    {
        if (c.nb > 0)
            ++nJetB;
        
        if (c.nb > 1)
            ++nJetDoubleB;
        
        if (c.nc > 0)
            ++nJetC;
        
        if (c.nc > 1)
            ++nJetDoubleC;
    }
    
    
    // Finally, classify the event
    if (nJetB >= 2)
        type = Type::TwoBottom;
    else if (nJetDoubleB > 0)
        type = Type::DoubleBottomJet;
    else if (nJetB > 0)
        type = Type::SingleBottom;
    else if (nJetC >= 2)
        type = Type::TwoCharm;
    else if (nJetDoubleC > 0)
        type = Type::DoubleCharmJet;
    else if (nJetC > 0)
        type = Type::SingleCharm;
    else if (unmatchedPartonCount.nb > 0)
        type = Type::UnmatchedBottom;
    else if (unmatchedPartonCount.nc > 0)
        type = Type::UnmatchedCharm;
    else
        type = Type::NoHF;
    
    
    return true;
}


ClassTTbarHFPlugin::Type ClassTTbarHFPlugin::GetDecision() const noexcept
{
    return type;
}
