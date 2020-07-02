/*                             __
                           .--()°"."
 Author: Nathalia Cardona "|, . ,"
                           !_-(_\
*/

#include "PhenoAnalyzer.h"
#include "LeptonCounter.h"
#include "Cuts.h"
#include "Physics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>
#include <utility>

using namespace std;

void writeCsv(int count, string path, string cut)
{
  ofstream outfile;
  outfile.open("/home/n.cardonac/AnalysisCode_Neutrinos_v.1/PhenoAnalyzer/counts.csv", ios_base::app); // append instead of overwrite
  outfile << path << "," << cut << "," << count << "\n";
}

int main(int argc, char *argv[])
{
  cout << "Starting phenoanalyzer..." << endl;

  // Importing Delphes data
  TChain chain("Delphes");
  chain.Add(argv[1]);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  // output file manager
  TFile *HistoOutputFile = new TFile(argv[2], "RECREATE");

  // directory to store the histograms
  TDirectory *noCutsDirectory = HistoOutputFile->mkdir("noCuts");

  // directory to store the histograms after the mu_mu_e cuts
  TDirectory *mu_mu_eCutsDirectory = HistoOutputFile->mkdir("mu_mu_e");

  // directory to store the histograms after the e_e_e cuts
  TDirectory *e_e_eCutsDirectory = HistoOutputFile->mkdir("e_e_e");

  // directory to store the histograms after the mu_mu_mu cuts
  TDirectory *mu_mu_muCutsDirectory = HistoOutputFile->mkdir("mu_mu_mu");

  // directory to store the histograms after the e_e_mu cuts
  TDirectory *e_e_muCutsDirectory = HistoOutputFile->mkdir("e_e_mu");

  // get tree info
  vector<string> branches = {
      "Electron",
      "Muon",
      "Jet",
      "MissingET"};

  map<string, TClonesArray *> branchDict;

  // create a dictionary with the branches
  for (int i = 0; (unsigned)i < branches.size(); i++)
  {
    TClonesArray *branch = treeReader->UseBranch(branches[i].c_str());
    branchDict[branches[i]] = branch;
  }

  // boolean mask to avoid over computation
  vector<int> cutsArr;

  for (int i = 0; (unsigned)i < treeReader->GetEntries(); i++)
  {
    cutsArr.push_back(1);
  }

  int nEvents;

  // write number of events to csv
  nEvents = (int)treeReader->GetEntries();
  writeCsv(nEvents, string(argv[1]), "C0");

  // open output file
  HistoOutputFile->cd();

  // ---------------------------------No cuts--------------------------------------------

  noCutsDirectory->cd();
  cout << "No cuts" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, noFilter);
  cout << "No cuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "noCuts");

  // -----------------------------------MU MU E cuts---------------------------------------

  mu_mu_eCutsDirectory->cd();
  cout << "mu_mu_e cuts" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, mu_mu_e_atLeast);
  cout << "mu_mu_e cuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "mu_mu_e");

  // -----------------------------------E E E cuts---------------------------------------

  e_e_eCutsDirectory->cd();
  cout << "e_e_e cuts" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, e_e_e_atLeast);
  cout << "e_e_e cuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "e_e_e");

  // -----------------------------------MU MU MU cuts---------------------------------------

  mu_mu_muCutsDirectory->cd();
  cout << "mu_mu_mu cuts" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, mu_mu_mu_atLeast);
  cout << "mu_mu_mu cuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "mu_mu_mu");

  // -----------------------------------E E MU cuts---------------------------------------

  e_e_muCutsDirectory->cd();
  cout << "e_e_mu cuts" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr, e_e_mu_atLeast);
  cout << "e_e_mu cuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "e_e_mu");

  // ------------------------------------------------------------------------------------

  // close output file
  cout << "closing output file" << endl;
  HistoOutputFile->Close();

  cout << "DONE." << endl;
}
