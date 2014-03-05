#include <TTbarDataDrivenPlugin.hpp>

#include <TFile.h>

#include <fstream>
#include <algorithm>
#include <stdexcept>


using namespace std;


TTbarDataDrivenPlugin::TTbarDataDrivenPlugin(string const &name_, Region region_) noexcept:
  Plugin(name_),
  JetBTagDataDrivenPlugin(name_),
  reader(nullptr),
  region(region_)
{}


Plugin *TTbarDataDrivenPlugin::Clone() const
{
  return new TTbarDataDrivenPlugin(name, region);
}


TTbarDataDrivenPlugin::~TTbarDataDrivenPlugin() noexcept
{

}


void TTbarDataDrivenPlugin::BeginRun(Dataset const &dataset)
{
    // Read in flavor combination fractions
    ifstream flavcombofile("flavpermfile_muon_tight_ttbar-mg_rev468_xYe.txt");
    if ( !flavcombofile.good() ) throw logic_error("flavcombofile not opened correctly");
    static float FlavCombos[11][5][5][5][5][5][5][5][5][5][5];
    for (int i0 = 4; i0 < 11; ++i0){  // N-jet
      for (int i1 = 1; i1 < 5; ++i1){  // Leading Jet
	for (int i2 = 1; i2 < 5; ++i2){  // Second Jet
	  for (int i3 = 1; i3 < 5; ++i3){  // ...
	    for (int i4 = 1; i4 < 5; ++i4){
	      for (int i5 = 0; i5 < 5; ++i5){
		for (int i6 = 0; i6 < 5; ++i6){
		  for (int i7 = 0; i7 < 5; ++i7){
		    for (int i8 = 0; i8 < 5; ++i8){
		      for (int i9 = 0; i9 < 5; ++i9){
			for (int i10 = 0; i10 < 5; ++i10){
			  flavcombofile >> FlavCombos[i0][i1][i2][i3][i4][i5][i6][i7][i8][i9][i10];
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    flavcombofile.close();

    // Fill maps for each N-jet category with non-zero fractions
    for (int i0 = 4; i0 < 11; ++i0){
      for (int i1 = 1; i1 < 5; ++i1){
	for (int i2 = 1; i2 < 5; ++i2){
	  for (int i3 = 1; i3 < 5; ++i3){
	    for (int i4 = 1; i4 < 5; ++i4){
	      for (int i5 = 0; i5 < 5; ++i5){
		for (int i6 = 0; i6 < 5; ++i6){
		  for (int i7 = 0; i7 < 5; ++i7){
		    for (int i8 = 0; i8 < 5; ++i8){
		      for (int i9 = 0; i9 < 5; ++i9){
			for (int i10 = 0; i10 < 5; ++i10){
			  float fracvalue = FlavCombos[i0][i1][i2][i3][i4][i5][i6][i7][i8][i9][i10];
			  if (i0 == 4 && fracvalue > 0.) FlavFracs4.insert(pair <float, vector<int> > (fracvalue,{i1,i2,i3,i4,i5,i6,i7,i8,i9,i10}));
			  if (i0 == 5 && fracvalue > 0.) FlavFracs5.insert(pair <float, vector<int> > (fracvalue,{i1,i2,i3,i4,i5,i6,i7,i8,i9,i10}));
			  if (i0 == 6 && fracvalue > 0.) FlavFracs6.insert(pair <float, vector<int> > (fracvalue,{i1,i2,i3,i4,i5,i6,i7,i8,i9,i10}));
			  if (i0 == 7 && fracvalue > 0.) FlavFracs7.insert(pair <float, vector<int> > (fracvalue,{i1,i2,i3,i4,i5,i6,i7,i8,i9,i10}));
			  if (i0 == 8 && fracvalue > 0.) FlavFracs8.insert(pair <float, vector<int> > (fracvalue,{i1,i2,i3,i4,i5,i6,i7,i8,i9,i10}));
			  if (i0 == 9 && fracvalue > 0.) FlavFracs9.insert(pair <float, vector<int> > (fracvalue,{i1,i2,i3,i4,i5,i6,i7,i8,i9,i10}));
			  if (i0 == 10 && fracvalue > 0.) FlavFracs10.insert(pair <float, vector<int> > (fracvalue,{i1,i2,i3,i4,i5,i6,i7,i8,i9,i10}));
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    // Read in Tagging Efficiencies
    eff2D = new TFile("effcalc_muon_tight_ttbar-mg_rev468_xYe.root");
    if ( !eff2D->IsOpen() ) throw logic_error("efficiency file not opened correctly");
    effbpt = (TH1F*)eff2D->Get("eff_bpt");
    effbeta = (TH1F*)eff2D->Get("eff_beta");
    effcpt = (TH1F*)eff2D->Get("eff_cpt");
    effceta = (TH1F*)eff2D->Get("eff_ceta");
    effudspt = (TH1F*)eff2D->Get("eff_udspt");
    effudseta = (TH1F*)eff2D->Get("eff_udseta");
    effgpt = (TH1F*)eff2D->Get("eff_gpt");
    effgeta = (TH1F*)eff2D->Get("eff_geta");
    effb2d = (TH1F*)eff2D->Get("eff_b2d");
    effc2d = (TH1F*)eff2D->Get("eff_c2d");
    effuds2d = (TH1F*)eff2D->Get("eff_uds2d");
    effg2d = (TH1F*)eff2D->Get("eff_g2d");

    if (dataset.IsMC()) isMC = true;

}


bool TTbarDataDrivenPlugin::ProcessEvent()
{
  // Put all the reconstructed jets into a single vector
  allJets.clear();
  allJets.insert(allJets.end(), (*reader)->GetJets().begin(), (*reader)->GetJets().end());

  // Count jets and tags
  int nJets = allJets.size();
  int NTagsTotal = 0;

  for (auto const &j: allJets){
    if (j.CSV() > 0.898 && fabs(j.Eta()) < 2.4)
      ++NTagsTotal;
  }

  if (NTagsTotal < 2 || NTagsTotal > 2) return false;

  // Define tagging efficiencies and btag SFs
  // This will later be done through the framework
  float TagEffAll[nJets][5][2];
  float TagEffAllE[nJets][5];
  float SFb = 0., SFc = 0., SFl = 0.;
    
  // Main loop over jets in the event
  for (int iJet = 0; iJet < nJets; ++iJet){

    float jetPT = allJets.at(iJet).Pt();
    float jetETA = fabs(allJets.at(iJet).Eta());
    if (jetPT >= 800.) jetPT = 799.; // just for SF uncertainty calculation

    int bin = 0;
    float ptmax[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
    while (bin < 15 and ptmax[bin] < jetPT)
      ++bin;
	    
    // Central
    SFb = ( 0.927563 + (1.55479e-05*jetPT) ) + ( -1.90666e-07*(jetPT*jetPT) );
    SFc = SFb;
    SFl = ((1.00462+(0.00325971*jetPT)) + (-7.79184e-06*(jetPT*jetPT))) + (5.22506e-09*(jetPT*(jetPT*jetPT)));

    // Tagging efficiencies 1-b's 2-c's 3-uds's 4-g's
    TagEffAll[iJet][1][1] = effb2d->GetBinContent(effbeta->FindBin(jetETA),effbpt->FindBin(jetPT));
    TagEffAll[iJet][2][1] = effc2d->GetBinContent(effceta->FindBin(jetETA),effcpt->FindBin(jetPT));
    TagEffAll[iJet][3][1] = effuds2d->GetBinContent(effudseta->FindBin(jetETA),effudspt->FindBin(jetPT));
    TagEffAll[iJet][4][1] = effg2d->GetBinContent(effgeta->FindBin(jetETA),effgpt->FindBin(jetPT));
      
    if (jetETA > 2.4){
      TagEffAll[iJet][1][1] = 0.;
      TagEffAll[iJet][2][1] = 0.;
      TagEffAll[iJet][3][1] = 0.;
      TagEffAll[iJet][4][1] = 0.;
    }

    if (!isMC){
      TagEffAll[iJet][1][1] = TagEffAll[iJet][1][1]*SFb;
      TagEffAll[iJet][2][1] = TagEffAll[iJet][2][1]*SFc;
      TagEffAll[iJet][3][1] = TagEffAll[iJet][3][1]*SFl;
      TagEffAll[iJet][4][1] = TagEffAll[iJet][4][1]*SFl;
    }

    // Tagging inefficiencies
    TagEffAll[iJet][1][0] = 1. - TagEffAll[iJet][1][1];
    TagEffAll[iJet][2][0] = 1. - TagEffAll[iJet][2][1];
    TagEffAll[iJet][3][0] = 1. - TagEffAll[iJet][3][1];
    TagEffAll[iJet][4][0] = 1. - TagEffAll[iJet][4][1];

    // Tagging efficency errors
    TagEffAllE[iJet][1] = effb2d->GetBinError(effbeta->FindBin(jetETA),effbpt->FindBin(jetPT));
    TagEffAllE[iJet][2] = effc2d->GetBinError(effceta->FindBin(jetETA),effcpt->FindBin(jetPT));
    TagEffAllE[iJet][3] = effuds2d->GetBinError(effudseta->FindBin(jetETA),effudspt->FindBin(jetPT));
    TagEffAllE[iJet][4] = effg2d->GetBinError(effgeta->FindBin(jetETA),effgpt->FindBin(jetPT));
    for (int checkinf = 0; checkinf < 5; ++checkinf){
      if (TagEffAllE[iJet][checkinf] != TagEffAllE[iJet][checkinf]) TagEffAllE[iJet][checkinf] = 1.;
    }

  }  // Done setting tagging efficiencies for jets

  float demCR = 0., numSR3T = 0., numSR4T = 0.;
  float effproduct = 0.;
  float num3TVar = 0.;
  float num4TVar = 0.;
  float demVar = 0.;
  ttweight3T = 0.; 
  float ttweight3TVar = 0.;
  ttweight4T = 0.; 
  float ttweight4TVar = 0.;

  // Used for choosing tagged jets
  int tagfirst = -1;
  int tagsecond = -1;
  int tagthird = -1;
  int tagfourth = -1;
  float tottagprob3t = 0;
  float tottagprob4t = 0;
  float tagcombo3t[10][10][10];
  float tagcombo4t[10][10][10][10];

  for (int tc1 = 0; tc1 < 10; ++tc1){
    for (int tc2 = 0; tc2 < 10; ++tc2){
      for (int tc3 = 0; tc3 < 10; ++tc3){
	tagcombo3t[tc1][tc2][tc3] = 0;
	for (int tc4 = 0; tc4 < 10; ++tc4){
	  tagcombo4t[tc1][tc2][tc3][tc4] = 0;
	}
      }
    }
  }

  // Choose a random flavor combination based of jet multiplicity
  map < float, vector<int> > njetflavcombo;
  if (nJets == 4) njetflavcombo.insert(FlavFracs4.begin(),FlavFracs4.end());
  if (nJets == 5) njetflavcombo.insert(FlavFracs5.begin(),FlavFracs5.end());
  if (nJets == 6) njetflavcombo.insert(FlavFracs6.begin(),FlavFracs6.end());
  if (nJets == 7) njetflavcombo.insert(FlavFracs7.begin(),FlavFracs7.end());
  if (nJets == 8) njetflavcombo.insert(FlavFracs8.begin(),FlavFracs8.end());
  if (nJets == 9) njetflavcombo.insert(FlavFracs9.begin(),FlavFracs9.end());
  if (nJets >= 10) njetflavcombo.insert(FlavFracs10.begin(),FlavFracs10.end());

  // Flavor Configurations
  vector<int> flavarray (10);  // Array that contains jet flavors ordered in pT
  float flavfraction = 0.;
  for(map <float, vector<int> >::iterator iter = njetflavcombo.begin(); iter != njetflavcombo.end(); ++iter){
    flavarray.clear();
    if (nJets < 10) flavarray.resize(nJets);
    else flavarray.resize(10);
    flavfraction = iter->first;
    copy ( (iter->second).begin(), (iter->second).end(), flavarray.begin() );
      
    // Tagged Configurations: not-tagged(0), tagged(1)
    for (int tagj1 = 0; tagj1 < 2; ++tagj1){
      for (int tagj2 = 0; tagj2 < 2; ++tagj2){
	for (int tagj3 = 0; tagj3 < 2; ++tagj3){
	  for (int tagj4 = 0; tagj4 < 2; ++tagj4){
	    for (int tagj5 = 0; tagj5 < 2; ++tagj5){
	      for (int tagj6 = 0; tagj6 < 2; ++tagj6){
		for (int tagj7 = 0; tagj7 < 2; ++tagj7){
		  for (int tagj8 = 0; tagj8 < 2; ++tagj8){
		    for (int tagj9 = 0; tagj9 < 2; ++tagj9){
		      for (int tagj10 = 0; tagj10 < 2; ++tagj10){

			int nJetmod = nJets;
			if (nJets > 10) nJetmod = 10;
			int jettags[10] = {tagj1,tagj2,tagj3,tagj4,tagj5,tagj6,tagj7,tagj8,tagj9,tagj10};
			int ctags = 0;
			for (int cjets = 0; cjets < nJetmod; ++cjets)
			  if (jettags[cjets] == 1) ++ctags;

			// Product of tagging (in)efficiencies
			effproduct = 1.;
			float effvarsum = 0.;
			for (int cjets = 0; cjets < nJetmod; ++cjets){
			  effproduct *= TagEffAll[cjets][flavarray[cjets]][jettags[cjets]];
			  if (TagEffAll[cjets][flavarray[cjets]][jettags[cjets]] > 0.)
			    effvarsum += pow(TagEffAllE[cjets][flavarray[cjets]]/TagEffAll[cjets][flavarray[cjets]][jettags[cjets]],2.);
			}
			float effprodVar = 0.;
			// Variance in this product
			if (effproduct > 0.) effprodVar = effproduct*effproduct*effvarsum;
				  
			// Two of the jets are tagged (2TCR)
			if (ctags == 2){
			  demVar += flavfraction*flavfraction*effprodVar;
			  demCR += flavfraction*effproduct;
			}
			  
			// Three of the jets are tagged (3TSR)
			if (ctags == 3){
			  num3TVar += flavfraction*flavfraction*effprodVar;
			  numSR3T += flavfraction*effproduct;
			  
			  // Finds all combinations with three tags
			  int counttags = 0;
			  for (int cjets = 0; cjets < nJetmod; ++cjets){
			    if (jettags[cjets] == 1){
			      if (counttags == 0) tagfirst = cjets;
			      if (counttags == 1) tagsecond = cjets;
			      if (counttags == 2) tagthird = cjets;
			      ++counttags;
			    }
			  }
			  tagcombo3t[tagfirst][tagsecond][tagthird] += flavfraction*effproduct;
			  tottagprob3t += flavfraction*effproduct;
			  
			}
			  
			// Four of the jets are tagged (4TSR)
			if (nJetmod > 4 && ctags == 4){
			  num4TVar += flavfraction*flavfraction*effprodVar;
			  numSR4T += flavfraction*effproduct;

			  // Finds all combinations with four tags
			  int counttags = 0;
			  for (int cjets = 0; cjets < nJetmod; ++cjets){
			    if (jettags[cjets] == 1){
			      if (counttags == 0) tagfirst = cjets;
			      if (counttags == 1) tagsecond = cjets;
			      if (counttags == 2) tagthird = cjets;
			      if (counttags == 3) tagfourth = cjets;
			      ++counttags;
			    }
			  }
			  tagcombo4t[tagfirst][tagsecond][tagthird][tagfourth] += flavfraction*effproduct;
			  tottagprob4t += flavfraction*effproduct;

			}

		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }  // End Tagging Combinations
  }  // End Flavor Combinations

  if (demCR > 0.){
    // Event weight = P(SR) / P(CR) for both 3T and 4T
    ttweight3T = numSR3T / demCR;
    ttweight4T = numSR4T / demCR;
    // Variance in the event weights (propagation of errors of numerator and denominator)
    if (numSR3T > 0. && demCR > 0.) ttweight3TVar = ttweight3T * ttweight3T * ( num3TVar / (numSR3T * numSR3T) + demVar / (demCR * demCR) );
    if (numSR4T > 0. && demCR > 0.) ttweight4TVar = ttweight4T * ttweight4T * ( num4TVar / (numSR4T * numSR4T) + demVar / (demCR * demCR) );
  }

  // Fills a map of the chosen tagging combinations for 3t and 4t
  // If P(3t) or P(4t) > 0
  map < float, vector <int> > choosetagconfig3t;
  map < float, vector <int> > choosetagconfig4t;
  for (int tc1 = 0; tc1 < 10; ++tc1){
    for (int tc2 = 0; tc2 < 10; ++tc2){
      for (int tc3 = 0; tc3 < 10; ++tc3){
	tagcombo3t[tc1][tc2][tc3] = tagcombo3t[tc1][tc2][tc3]/tottagprob3t;
	if (tagcombo3t[tc1][tc2][tc3] > 0.) choosetagconfig3t.insert(pair <float, vector<int> > (tagcombo3t[tc1][tc2][tc3],{tc1,tc2,tc3}));
	for (int tc4 = 0; tc4 < 10; ++tc4){
	  tagcombo4t[tc1][tc2][tc3][tc4] = tagcombo4t[tc1][tc2][tc3][tc4]/tottagprob4t;
	  if (tagcombo4t[tc1][tc2][tc3][tc4] > 0.) choosetagconfig4t.insert(pair <float, vector<int> > (tagcombo4t[tc1][tc2][tc3][tc4],{tc1,tc2,tc3,tc4}));
	}
      }
    }
  }
	
  if (region == Region::r3t)
  {
    // Randomly pick a tagging combination for the 3t model
    float picktagconfig = r3.Rndm();
    vector <int> tagarray3t (3);  // Contains the indices of tagged jets
    float sumprobs = 0.;
    for(map <float, vector<int> >::iterator iter = choosetagconfig3t.begin(); iter != choosetagconfig3t.end(); ++iter){
      sumprobs += iter->first;
      if (sumprobs > picktagconfig){
        copy ( (iter->second).begin(), (iter->second).end(), tagarray3t.begin() );
        break;
      }
    }
    sort(tagarray3t.begin(),tagarray3t.end());
    
    // Fill a vector that defines each jet being tagged (1) or not (0)
    for (int cjets = 0; cjets < nJets; ++cjets){
      bool foundtag = false;
      for (int ctags = 0; ctags < 3; ++ctags)
        if (tagarray3t[ctags] == cjets)
  	foundtag = true;
      if (foundtag) taggedjets3t.push_back(1);
      else taggedjets3t.push_back(0);
    }
  }
  else
  {
    // Randomly pick a tagging combination for the 4t model
    picktagconfig = r3.Rndm();
    vector <int> tagarray4t (4);
    sumprobs = 0.;
    for(map <float, vector<int> >::iterator iter = choosetagconfig4t.begin(); iter != choosetagconfig4t.end(); ++iter){
      sumprobs += iter->first;
      if (sumprobs > picktagconfig){
        copy ( (iter->second).begin(), (iter->second).end(), tagarray4t.begin() );
        break;
      }
    }
    sort(tagarray4t.begin(),tagarray4t.end());

    // Fill a vector that defines each jet being tagged (1) or not (0)
    for (int cjets = 0; cjets < nJets; ++cjets){
      bool foundtag = false;
      for (int ctags = 0; ctags < 4; ++ctags)
        if (tagarray4t[ctags] == cjets)
  	foundtag = true;
      if (foundtag) taggedjets4t.push_back(1);
      else taggedjets4t.push_back(0);
    }
  }
    
  return true;

}


vector<int> const &TTbarDataDrivenPlugin::GetTaggedJetIndices() const
{
  if (region == Region::r3t) return taggedjets3t;
  else return taggedjets4t;
}


double TTbarDataDrivenPlugin::GetWeight() const
{
  if (region == Region::r3t) return ttweight3T;
  else return ttweight4T;
}
