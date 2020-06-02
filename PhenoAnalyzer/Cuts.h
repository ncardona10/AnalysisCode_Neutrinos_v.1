/*
                               __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Cuts for single particles
*/

#include "./ROOTFunctions.h"
#include "./DelphesFunctions.h"
#include "./Physics.h"
#include "./Overlaps.h"
#include "./LeptonCounter.h"

bool noFilter(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr)
{
  return true;
}

bool met(ExRootTreeReader *treeReader,
         map<string, TClonesArray *> branchDict,
         int entry,
         vector<int> &cutsArr)
{
  // met>150

  if (cutsArr[entry])
  {
    treeReader->ReadEntry(entry);

    bool metBool = met(treeReader, branchDict, entry) > 150;

    cutsArr[entry] = metBool;
    return metBool;
  }

  return false;
}

bool bjets(ExRootTreeReader *treeReader,
           map<string, TClonesArray *> branchDict,
           int entry,
           vector<int> &cutsArr)
{
  // # bjets = 0

  treeReader->ReadEntry(entry);

  bool bJetsBool = true;

  if (cutsArr[entry])
  {

    for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->BTag == 1 && jet->PT > 30.0 && abs(jet->Eta) < 2.5)
      {
        bJetsBool = false;
        break;
      }
    }
  }
  else
  {
    bJetsBool = false;
  }
  cutsArr[entry] = bJetsBool;
  return bJetsBool;
}

// Verify if the leading jet satisfy the pt > 100 GeV and |eta| < 5.0
// Leading jet: Jet with the highest Pt
bool jetsCuts(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr)
{
  treeReader->ReadEntry(entry);
  if (cutsArr[entry])
  {
    bool ans = false;
    float max_pt = -1;

    // Create the leading jet
    Jet *leadingJet;

    // Find the leading jet
    for (int i = 0; i < branchDict["Jet"]->GetEntries(); i++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(i);
      if (jet->PT > max_pt)
      {
        max_pt = jet->PT;
        leadingJet = jet;
      }
    }

    // Check that at least one jet exists
    if (max_pt != -1)
    {
      // verify Pt and |eta| conditions
      ans = (leadingJet->PT > 100) && abs(leadingJet->Eta)<5.0;
    }

    cutsArr[entry] = ans;

    return ans;
  }
  return false;
}

//Verify if exist at least one electron and tau that satisfy the minimum pt and |eta| cuts
bool leptonsCuts(ExRootTreeReader *treeReader,
                 map<string, TClonesArray *> branchDict,
                 int entry,
                 vector<int> &cutsArr)
{
  treeReader->ReadEntry(entry);
  if (cutsArr[entry])
  {

    bool ansElec = false;
    bool ansTau = false;
    bool ans = false;

    int cont = 0;

    while (!ansElec && cont < branchDict["Electron"]->GetEntries())
    {
      // Get a lepton
      Electron *elec = (Electron *)branchDict["Electron"]->At(cont);

      // To create a muon use:
      //Muon *muon = (Muon *)branchDict["Muon"]->At(cont);

      if (elec->PT > 5.0 && abs(elec->Eta) < 2.4)
      {
        ansElec = true;
      }
      cont++;
    }

    //Reset counter
    cont = 0;

    while (!ansTau && cont < branchDict["Jet"]->GetEntries())
    {
      // Get a jet
      Jet *jet = (Jet *)branchDict["Jet"]->At(cont);

      // Verify if the jet is a tau
      if (jet->TauTag == 1)
      {
        if (jet->PT >= 20.0 && abs(jet->Eta) < 2.4)
        {
          ansTau = true;
        }
      }
      cont++;
    }

    ans = ansElec && ansTau;

    cutsArr[entry] = ans;
    return ans;
  }
  return false;
}