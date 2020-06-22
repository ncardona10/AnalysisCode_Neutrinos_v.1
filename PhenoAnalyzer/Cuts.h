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

//Baseline to create the lepton channels. Conside only an specific number for every lepton
bool nParticle(ExRootTreeReader *treeReader,
               map<string, TClonesArray *> branchDict,
               int entry,
               int n_electrons,
               int n_muon,
               int n_tau,
               vector<int> &cutsArr,
               bool checkElecPT,
               bool checkMuonPT,
               bool checkTauPT)
{

  /*
    must comply with VBF cuts  & cuts
    n_electrons electron with:
      pt>8 
      abs(eta) < 2.4
    n_tau taus with:
      pt>20
      abs(eta)<2.4
    n_muon muons with:
      pt>5
      abs(eta)<2.4
  */

  treeReader->ReadEntry(entry);

  if (cutsArr[entry])
  {
    // verify electron condition
    int nElectrons = 0;
    int i = 0;
    vector<TLorentzVector> particleCharacteristics;

    while (nElectrons <= n_electrons && i < branchDict["Electron"]->GetEntries())
    {
      Electron *elec = (Electron *)branchDict["Electron"]->At(i);
      if (elec->PT >= 8 && abs(elec->Eta) < 2.4)
      {
        if (checkElecPT)
        {
          if (elec->PT <= 40)
          {
            nElectrons++;
            particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
          }
        }

        else
        {
          nElectrons++;
          particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
        }
      }
      i++;
    }

    if (nElectrons == n_electrons)
    {
      //verify number of muons
      int nMuons = 0;
      int i = 0;

      particleCharacteristics.clear();

      while (nMuons <= n_muon && i < branchDict["Muon"]->GetEntries())
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(i);
        if (muon->PT >= 5 && abs(muon->Eta) < 2.4)
        {

          if (checkMuonPT)
          {
            if (muon->PT <= 40)
            {
              nMuons++;
              particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
            }
          }
          else
          {
            nMuons++;
            particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
          }
        }
        i++;
      }

      if (nMuons == n_muon)
      {
        //verify number of taus
        int nTaus = 0;
        int i = 0;
        particleCharacteristics.clear();

        while (nTaus <= n_tau && i < branchDict["Jet"]->GetEntries())
        {
          Jet *jet = (Jet *)branchDict["Jet"]->At(i);
          if (jet->TauTag == 1)
          {
            if (checkTauPT)
            {
              if (jet->PT >= 20.0 && abs(jet->Eta) < 2.4)
              {
                nTaus++;
                particleCharacteristics.push_back(createTLorentzVector(jet->PT, jet->Eta, jet->Mass, jet->Phi));
              }
            }
          }
          i++;
        }
        return nTaus == n_tau;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

/*
SINGLE PARTICLE CHANNEL
*/

bool mono_e(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry,
            vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 1, 0, 0, cutsArr, true, false, false);
}

bool mono_mu(ExRootTreeReader *treeReader,
             map<string, TClonesArray *> branchDict,
             int entry,
             vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 0, 1, 0, cutsArr, false, true, false);
}

bool mono_tau(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 0, 0, 1, cutsArr, false, false, true);
}

//Baseline to create the lepton channels. Consider at least an specific number for every lepton
bool nParticle_atLeast(ExRootTreeReader *treeReader,
               map<string, TClonesArray *> branchDict,
               int entry,
               int n_electrons,
               int n_muon,
               int n_tau,
               vector<int> &cutsArr,
               bool checkElecPT,
               bool checkMuonPT,
               bool checkTauPT)
{

  /*
    n_electrons electron with:
      pt>8 
      abs(eta) < 2.4
    n_tau taus with:
      pt>20
      abs(eta)<2.4
    n_muon muons with:
      pt>5
      abs(eta)<2.4
  */

  treeReader->ReadEntry(entry);

  if (cutsArr[entry])
  {
    // verify electron condition
    int nElectrons = 0;
    int i = 0;
    vector<TLorentzVector> particleCharacteristics;

    while (nElectrons <= n_electrons && i < branchDict["Electron"]->GetEntries())
    {
      Electron *elec = (Electron *)branchDict["Electron"]->At(i);
      if (elec->PT >= 8 && abs(elec->Eta) < 2.4)
      {
        if (checkElecPT)
        {
          if (elec->PT <= 40)
          {
            nElectrons++;
            particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
          }
        }

        else
        {
          nElectrons++;
          particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
        }
      }
      i++;
    }

    if (nElectrons >= n_electrons)
    {
      //verify number of muons
      int nMuons = 0;
      int i = 0;

      particleCharacteristics.clear();

      while (nMuons <= n_muon && i < branchDict["Muon"]->GetEntries())
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(i);
        if (muon->PT >= 5 && abs(muon->Eta) < 2.4)
        {

          if (checkMuonPT)
          {
            if (muon->PT <= 40)
            {
              nMuons++;
              particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
            }
          }
          else
          {
            nMuons++;
            particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
          }
        }
        i++;
      }

      if (nMuons >= n_muon)
      {
        //verify number of taus
        int nTaus = 0;
        int i = 0;
        particleCharacteristics.clear();

        while (nTaus <= n_tau && i < branchDict["Jet"]->GetEntries())
        {
          Jet *jet = (Jet *)branchDict["Jet"]->At(i);
          if (jet->TauTag == 1)
          {
            if (checkTauPT)
            {
              if (jet->PT >= 20.0 && abs(jet->Eta) < 2.4)
              {
                nTaus++;
                particleCharacteristics.push_back(createTLorentzVector(jet->PT, jet->Eta, jet->Mass, jet->Phi));
              }
            }
          }
          i++;
        }
        return nTaus >= n_tau;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

/*
MULTILEPTONS CHANNEL
*/

bool mu_mu_e_atLeast(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr)
{
  return nParticle(treeReader, branchDict, entry, 1, 2, 0, cutsArr, false, false, false);
}
