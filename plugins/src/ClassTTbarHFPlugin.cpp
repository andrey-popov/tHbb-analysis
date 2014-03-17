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
    
    
    // Find the two b quarks from decays of tops
    array<GenParticle const *, 2> bTops;
    unsigned numBTopsFound = 0;
    
    for (GenParticle const &p: (*reader)->GetHardGenParticles())
    {
        if (abs(p.GetPdgId()) != 5)
            continue;
        
        if (abs(p.GetFirstMother()->GetPdgId()) == 6)
        {
            bTops.at(numBTopsFound) = &p;
            ++numBTopsFound;
        }
        
        if (numBTopsFound == 2)  // there is nothing else to check here
            break;
    }
    
    // A sanity check
    if (numBTopsFound != 2)
        throw runtime_error("ClassTTbarHFPlugin::ProcessEvent: Failed to find two b quarks from "
         "decays of the tops.");
    
    
    // Remove from the collection of additional partons two b quarks closest to the quarks from
    //decays of tops
    for (auto const &bTop: bTops)
    {
        auto closestPartonIt = additionalPartons.end();
        double minDistance = numeric_limits<double>::infinity();
        
        for (auto pIt = additionalPartons.begin(); pIt != additionalPartons.end(); ++pIt)
        {
            if (abs((*pIt)->GetPdgId()) != 5)
                continue;
            
            double const dR = (*pIt)->P4().DeltaR(bTop->P4());
            
            if (dR < minDistance)
            {
                closestPartonIt = pIt;
                minDistance = dR;
            }
        }
        
        // Sanity check
        if (closestPartonIt == additionalPartons.end())
            throw runtime_error("ClassTTbarHFPlugin::ProcessEvent: Cannot match a b quark from "
             "decay of a top to any b quark from the parton shower.");
        
        
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
