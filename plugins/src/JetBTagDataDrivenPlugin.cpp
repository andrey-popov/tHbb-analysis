#include <JetBTagDataDrivenPlugin.hpp>

#include <Processor.hpp>

#include <algorithm>
#include <stdexcept>


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
    reader(src.reader), taggedJetIndices(src.taggedJetIndices)
{}


JetBTagDataDrivenPlugin::~JetBTagDataDrivenPlugin() noexcept
{}


void JetBTagDataDrivenPlugin::BeginRun(Dataset const &)
{
    // Save pointer to the reader plugin
    reader = dynamic_cast<PECReaderPlugin const *>(processor->GetPluginBefore("Reader", name));
}


vector<unsigned> const &JetBTagDataDrivenPlugin::GetTaggedJetIndices() const
{
    return taggedJetIndices;
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
    auto const res = find(taggedJetIndices.begin(), taggedJetIndices.end(), index);
    
    return (res != taggedJetIndices.end());
}
