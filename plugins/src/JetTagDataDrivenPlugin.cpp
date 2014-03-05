#include <JetTagDataDrivenPlugin.hpp>

#include <algorithm>
#include <stdexcept>


using namespace std;


JetTagDataDrivenPlugin::JetTagDataDrivenPlugin(string const &name_) noexcept:
    Plugin(name_),
    reader(nullptr)
{}


JetTagDataDrivenPlugin::~JetTagDataDrivenPlugin() noexcept
{}


void JetTagDataDrivenPlugin::BeginRun(Dataset const &dataset)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


vector<unsigned> const &JetTagDataDrivenPlugin::GetTaggedJetIndices() const
{
    return taggedJetIndices;
}


bool JetTagDataDrivenPlugin::IsTagged(Jet const &jet) const
{
    auto const &jets = (*reader)->GetJets();
    
    // Loop over all jets in the current event
    for (unsigned i = 0; i < jets.size(); ++i)
    {
        if (&jet == &jets.at(i))
            return IsTagged(i);
    }
    
    
    throw logic_error("JetTagDataDrivenPlugin::IsTagged: There is no such jet in the current "
     "event.");
    
    return false;
}


bool JetTagDataDrivenPlugin::IsTagged(unsigned index) const
{
    auto const res = find(taggedJetIndices.begin(), taggedJetIndices.end(), index);
    
    return (res != taggedJetIndices.end());
}
