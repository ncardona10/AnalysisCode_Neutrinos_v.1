/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Counts the number of leptons in different PT ranges.                        
*/

#ifndef LEPTONCOUNTER_H
#define LEPTONCOUNTER_H

#include "./ROOTFunctions.h"
#include "./DelphesFunctions.h"
#include "../plots/MyHistograms.h"
#include "./Physics.h"
#include "./Overlaps.h"
#include <string>
#include <map>
#include <vector>
#include <set>

void fillHisto(TH1 *histo, float value)
{
  if (value > 0)
  {
    histo->Fill(value);
  }
}

bool inSet(int val, set<int> theSet)
{
  return theSet.count(val) > 0;
}

int ptEtaPhiMjjMt(ExRootTreeReader *treeReader,
                  map<string, TClonesArray *> branchDict,
                  vector<int> &cutsArr,
                  bool (*filter)(ExRootTreeReader *,
                                 map<string, TClonesArray *>,
                                 int,
                                 vector<int> &))
{

  int numEvents = 0;

  vector<string> variables = {"pt", "eta", "phi"};
  vector<string> particleTypes = {"electron",
                                  "muon",
                                  "tau",
                                  "jet"};

  // create histograms
  map<string, TH1 *> histos;
  for (int i = 0; (unsigned)i < variables.size(); i++)
  {
    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      int bins = 100;
      float x_min = 0.0;
      float x_max = 15.0;

      if (variables[i].compare("pt") == 0)
      {
        x_max = 500;
        bins = 500;
      }
      else if (variables[i].compare("phi") == 0)
      {
        x_max = TMath::Pi();
        x_min = -TMath::Pi();
      }
      else
      {
        x_min = -5;
        x_max = 5;
      }
      histos[variables[i] + particleTypes[j]] = blankHistogram(particleTypes[j] + " " + variables[i],
                                                               variables[i] + particleTypes[j],
                                                               bins, x_min, x_max); // check the histogram limits & bins
    }
  }

  histos["MET"] = blankHistogram("MET", "MET", 100, 0, 1000);

  Long64_t numberOfEntries = treeReader->GetEntries();

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);

    if (filter(treeReader, branchDict, entry, cutsArr))
    {
      numEvents += 1;

      // MET stuff
      MissingET *METPointer = (MissingET *)branchDict["MissingET"]->At(0);
      double MET = METPointer->MET;

      // taus
      for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
      {
        Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
        if (jet->TauTag == 1)
        {
          // its a tau!
          if ((abs(jet->Eta) < 2.4) && (jet->PT > 20.0))
          {
            histos["pttau"]->Fill(jet->PT);
            histos["etatau"]->Fill(jet->Eta);
            histos["phitau"]->Fill(normalizedDphi(jet->Phi));
          }
        }
      }

      cout << "electrons" << endl;

      // electrons
      for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
      {
        Electron *electron = (Electron *)branchDict["Electron"]->At(leaf);

        if ((abs(electron->Eta) < 2.4) && (electron->PT > 8.0))
        {
          histos["ptelectron"]->Fill(electron->PT);
          histos["etaelectron"]->Fill(electron->Eta);
          histos["phielectron"]->Fill(normalizedDphi(electron->Phi));
        }
      }

      cout << "muons" << endl;
      // muons
      for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(leaf);

        histos["ptmuon"]->Fill(muon->PT);
        histos["etamuon"]->Fill(muon->Eta);
        histos["phimuon"]->Fill(normalizedDphi(muon->Phi));
      }

      //jets
      for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
      {
        Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);

        if (!jet->TauTag)
        {
          histos["ptjet"]->Fill(jet->PT);
          histos["etajet"]->Fill(jet->Eta);
          histos["phijet"]->Fill(normalizedDphi(jet->Phi));
        }
      }

      //MET
      histos["MET"]->Fill(MET);
    }
  }

  for (int i = 0; (unsigned)i < variables.size(); i++)
  {
    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      histos[variables[i] + particleTypes[j]]->Write();
    }
  }

  histos["MET"]->Write();

  return numEvents;
}

#endif