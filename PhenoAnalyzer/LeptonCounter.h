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

  vector<string> variables = {"pt", "eta", "phi", "Mt", "leadingPT", "dphi_MET_"};
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

      if (variables[i].compare("pt") == 0 || variables[i].compare("leadingPT") == 0)
      {
        x_max = 500;
        bins = 500;
      }
      else if (variables[i].compare("phi") == 0)
      {
        x_max = TMath::Pi();
        x_min = -TMath::Pi();
      }
      else if (variables[i].compare("dphi_MET_") == 0)
      {
        x_max = TMath::Pi()*2;
      }
      else if (variables[i].compare("Mt") == 0)
      {
        x_max = 1000;
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
      float maxPt_tau = -1.0;
      Jet *leadingTau;
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

            double mtval = mt(jet->PT, MET, jet->Phi - METPointer->Phi);
            histos["Mttau"]->Fill(mtval);

            //leading pt
            if(jet->PT > maxPt_tau)
            {
              maxPt_tau = jet->PT;
              leadingTau = jet;
            }
          }
        }
      }
      if (maxPt_tau != -1.0)
      {
        histos["leadingPTtau"]->Fill(maxPt_tau);
        float dphi = abs(normalizedDphi(leadingTau->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_tau"]->Fill(dphi);
      }


      // electrons
      float maxPt_electrons = -1.0;
      Electron *leadingElectron;
      for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
      {
        Electron *electron = (Electron *)branchDict["Electron"]->At(leaf);

        if ((abs(electron->Eta) < 2.4) && (electron->PT > 8.0))
        {
          histos["ptelectron"]->Fill(electron->PT);
          histos["etaelectron"]->Fill(electron->Eta);
          histos["phielectron"]->Fill(normalizedDphi(electron->Phi));

          double mtval = mt(electron->PT, MET, electron->Phi - METPointer->Phi);
          histos["Mtelectron"]->Fill(mtval);

          //leading pt
          if(electron->PT > maxPt_electrons)
          {
            maxPt_electrons = electron->PT;
            leadingElectron = electron;
          }
        }
      }
      if (maxPt_electrons != -1.0)
      {
        histos["leadingPTelectron"]->Fill(maxPt_electrons);
        float dphi = abs(normalizedDphi(leadingElectron->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_electron"]->Fill(dphi);
      }

      // muons
      float maxPt_muons = -1.0;
      Muon *leadingMuon;
      for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(leaf);
        if ((abs(muon->Eta) < 2.4) && (muon->PT > 5.0))
        {
          histos["ptmuon"]->Fill(muon->PT);
          histos["etamuon"]->Fill(muon->Eta);
          histos["phimuon"]->Fill(normalizedDphi(muon->Phi));

          double mtval = mt(muon->PT, MET, muon->Phi - METPointer->Phi);
          histos["Mtmuon"]->Fill(mtval);

          //leading pt
          if(muon->PT > maxPt_muons)
          {
            maxPt_muons = muon->PT;
            leadingMuon = muon;
          }
        }
      }
      if (maxPt_muons != -1.0)
      {
        histos["leadingPTmuon"]->Fill(maxPt_muons);
        float dphi = abs(normalizedDphi(leadingMuon->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_muon"]->Fill(dphi);
      }

      //jets: jets minimum cuts are applied in jets cuts
      float maxPt_jet = -1.0;
      Jet *leadingJet;
      for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
      {
        Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
        if (!jet->TauTag)
        {
          histos["ptjet"]->Fill(jet->PT);
          histos["etajet"]->Fill(jet->Eta);
          histos["phijet"]->Fill(normalizedDphi(jet->Phi));

          //Doesnt make sense but needed to keep code symmetry
          double mtval = mt(jet->PT, MET, jet->Phi - METPointer->Phi);
          histos["Mtjet"]->Fill(mtval);

          //leading pt
          if(jet->PT > maxPt_jet)
          {
            maxPt_jet = jet->PT;
            leadingJet = jet;
          }
        }
      }
      if (maxPt_jet != -1.0)
      {
        histos["leadingPTjet"]->Fill(maxPt_jet);
        float dphi = abs(normalizedDphi(leadingJet->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_jet"]->Fill(dphi);        
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
