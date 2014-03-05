#include <JetBTagDataDrivenPlugin.hpp>

#include <Processor.hpp>

#include <algorithm>
#include <stdexcept>
#include <sstream>


using namespace std;


JetBTagDataDrivenPlugin::JetBTagDataDrivenPlugin(string const &name_) noexcept:
    Plugin(name_),
    //^ Must call constructor for the virtual base class explicitly
    BTaggerPlugin(name_), EventWeightPlugin(name_),
    reader(nullptr)
{}


JetBTagDataDrivenPlugin::JetBTagDataDrivenPlugin(JetBTagDataDrivenPlugin const &src) noexcept:
    Plugin(src),
    BTaggerPlugin(src), EventWeightPlugin(src),
    reader(src.reader), jetTagBits(src.jetTagBits)
{}


JetBTagDataDrivenPlugin::~JetBTagDataDrivenPlugin() noexcept
{}


void JetBTagDataDrivenPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


deque<bool> const &JetBTagDataDrivenPlugin::GetJetTagDecisions() const
{
    return jetTagBits;
}


bool JetBTagDataDrivenPlugin::IsTagged(Jet const &jet) const
{
    auto const &jets = (*reader)->GetJets();
    
    // Loop over all jets in the current event
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        if (&jet == &jets.at(i))
            return IsTagged(i);
    }
    
    
    throw logic_error("JetBTagDataDrivenPlugin::IsTagged: There is no such jet in the current "
     "event.");
    
    return false;
}


bool JetBTagDataDrivenPlugin::IsTagged(unsigned index) const
{
    if (index >= jetTagBits.size())
    {
        ostringstream ost;
        ost << "JetBTagDataDrivenPlugin::IsTagged: The given jet index " << index << " does not " <<
         "match the size of the mask of tagged jets, which is " << jetTagBits.size() << ". ";
        
        if (jetTagBits.size() != (*reader)->GetJets().size())
            ost << "The size of the mask does not agree with the size of the jet collection (" <<
             (*reader)->GetJets().size() << ") as well.";
        else
            ost << "The size of the mask agrees with the size of the jet collection.";
        
        throw logic_error(ost.str());
    }
    
    return jetTagBits.at(index);
}
