// Raghav Kunnawalkam Elayavalli
// Feb 17th 2015
// Rutgers

// for questions or comments: raghav.k.e at CERN dot CH

//
// 
// Set up macro to read in different MC pthat files and generate a combined spectra 
// The cross section weights are taken from dijet samples and the values are given here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForestPA2013#Dijet_Cross_Sections_for_reweigh
// 
//


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

static const int nbins_pt = 39;
static const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 272, 300, 
  330, 362, 395, 430,
  468, 507, 548, 592,
  638, 686, 1000 
};

// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

static const int nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
//we have to multiply by 5, since centrality bin variable goes from 0-200. 

// finds the centrality bin of the given event. 
int findBin(int hiBin){
  int binNo = 0;

  for(int i = 0;i<nbins_cent;i++){
    if(hiBin>=5*boundaries_cent[i] && hiBin<5*boundaries_cent[i+1]) {
      binNo = i;
      break;
    }
  }

  return binNo;
}


// class JetData used to set branch address and variables for each pthat file (instead of loading them individually)
class JetData
{
public:
  JetData(char *fileName, char *jetTree) {

    cout <<"Open "<<fileName<<endl;

    tFile = new TFile(fileName,"read");
    tEvt = (TTree*)tFile->Get("hiEvtAnalyzer/HiTree");
    tSkim = (TTree*)tFile->Get("skimanalysis/HltTree");
    tHlt = (TTree*)tFile->Get("hltanalysis/HltTree");
    tJet = (TTree*)tFile->Get(jetTree);
    
    tJet->SetBranchAddress("jtpt" , jtpt );
    tJet->SetBranchAddress("jtpu", jtpu );
    tEvt->SetBranchAddress("hiNpix",&hiNpix);
    tEvt->SetBranchAddress("run",&run);
    tEvt->SetBranchAddress("evt",&evt);
    tEvt->SetBranchAddress("lumi",&lumi);
    tEvt->SetBranchAddress("hiNtracks",&hiNTracks);
    tEvt->SetBranchAddress("hiHF",&hiHF);

    tJet->SetBranchAddress("rawpt", rawpt);
    tJet->SetBranchAddress("trackMax" , trackMax );
    tJet->SetBranchAddress("chargedMax",chargedMax);
    tJet->SetBranchAddress("chargedSum",chargedSum);
    tJet->SetBranchAddress("neutralMax",neutralMax);
    tJet->SetBranchAddress("neutralSum",neutralSum);
    tJet->SetBranchAddress("photonSum",photonSum);
    tJet->SetBranchAddress("photonMax",photonMax);
    tJet->SetBranchAddress("eSum",eSum);
    tJet->SetBranchAddress("eMax",eMax);
    tJet->SetBranchAddress("muSum",muSum);
    tJet->SetBranchAddress("muMax",muMax);
    tJet->SetBranchAddress("refpt", refpt);
    tJet->SetBranchAddress("nref" ,&njets);
    tJet->SetBranchAddress("jteta", jteta);
    tJet->SetBranchAddress("jtphi", jtphi);
    tJet->SetBranchAddress("jtm", jtmass);
    tJet->SetBranchAddress("pthat", &pthat);
    tJet->SetBranchAddress("subid",&subid);
    tEvt->SetBranchAddress("hiBin",&bin);
    tEvt->SetBranchAddress("vz",&vz);
    tEvt->SetBranchAddress("vx",&vx);
    tEvt->SetBranchAddress("vy",&vy);
    tSkim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    tSkim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    tHlt->SetBranchAddress("HLT_HIJet80_v7",&jet80_1);
    tHlt->SetBranchAddress("HLT_HIJet65_v7",&jet65_1);
    tHlt->SetBranchAddress("HLT_HIJet55_v7",&jet55_1);
    tHlt->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_p_1);
    tHlt->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_p_1);
    tHlt->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_p_1);
    tJet->AddFriend(tEvt);
    tJet->AddFriend(tSkim);
    tJet->AddFriend(tHlt);
  };
  
  TFile *tFile;
  TTree *tJet;
  TTree *tEvt;
  TTree *tHlt;
  TTree *tSkim;
 
  //event varianbles
  int run;
  int evt;
  int lumi;
  int hiNTracks;
  Float_t hiHF;
  float vz;
  float vx;
  float vy;
  float pthat;
  int hiNpix;
  
  //jet variables 
  float jtpt[1000];
  float jtpu[1000];
  float rawpt[1000];
  float refpt[1000];
  float jteta[1000];
  float jtphi[1000];
  float jtmass[1000];
  float trackMax[1000];
  float chargedMax[1000];
  float neutralMax[1000];
  float chargedSum[1000];
  float neutralSum[1000];
  float photonSum[1000];
  float eSum[1000];
  float muSum[1000];
  float photonMax[1000];
  float eMax[1000];
  float muMax[1000];
  float subid[1000];

  int njets;
  int bin;     
  int pHBHENoiseFilter;
  int pcollisionEventSelection;
  int jet55_1;
  int jet65_1;
  int jet80_1;
  int jet55_p_1;
  int jet65_p_1;
  int jet80_p_1;

};



void MC_pthat_addition(char *algo = "Pu", int radius = 3, char *jet_type = "PF", int sub_id = 0){

  TDatime date;
  TStopwatch timer;
  timer.Start();

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);
  
  cout<<"Running for Algorithm "<<algo<<" "<<jet_type<<endl;
 
  bool printDebug = true;

  // get the pthat files: 
  const int nbins_pthat = 6;
  Double_t boundaries_pthat[nbins_pthat+1];
  char *fileName_pthat[nbins_pthat];
  Double_t xsection[nbins_pthat+1];

  boundaries_pthat[0] = 30;
  fileName_pthat[0] = "pthat30/HiForest_pthat30.root";
  xsection[0] = 1.075e-02;
  
  boundaries_pthat[1] = 50;
  fileName_pthat[1] = "pthat50/HiForest_pthat50.root";
  xsection[1] = 1.025e-03;
  
  boundaries_pthat[2] = 80;
  fileName_pthat[2] = "pthat80/HiForest_pthat80.root";
  xsection[2] = 9.865e-05;
  
  boundaries_pthat[3] = 120;
  fileName_pthat[3] = "pthat120/HiForest_pthat120.root";
  xsection[3] = 1.129e-05;
  
  boundaries_pthat[4] = 170;
  fileName_pthat[4] = "pthat170/HiForest_pthat170.root";
  xsection[4] = 1.465e-06;
  
  boundaries_pthat[5] = 220;
  fileName_pthat[5] = "pthat220/HiForest_pthat220.root";
  xsection[5] = 2.837e-07;

  boundaries_pthat[6] = 1000;
  xsection[6] = 0;

  
  // Vertex & centrality reweighting for PbPb
  TF1 *fVz;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);

  // get the centrality weight from the root file created in the plotting macro. 
  TFile *fcentin = TFile::Open("PbPb_DataMC_cent_ratio_20141117.root");
  TH1F *hCentWeight = (TH1F*)fcentin->Get("hCentRatio");

  // declare the histograms: sample jet pT here:
  TH1F * hPtHatRaw = new TH1F("hPtHatRaw","",nbins_pthat,boundaries_pthat);
  TH1F * hPtHat    = new TH1F("hPtHat","",nbins_pthat,boundaries_pthat);
  TH1F * hJetpT    = new TH1F("hJetPT","Jet pT spectra",1000,0,1000);
  
  
  JetData *PbPbMC[nbins_pthat];
  if(printDebug)cout<<"reading all the PbPb MC files"<<endl;
  for(int h = 0;h<nbins_pthat;h++){
    PbPbMC[h] = new JetData(fileName_pthat[h],Form("ak%s%d%sJetAnalyzer/t",algo,radius,jet_type));
    TH1F* hPtHatTmp = new TH1F("hPtHatTmp","",nbins_pthat,boundaries_pthat);
    PbPbMC[h]->tJet->Project("hPtHatTmp","pthat");
    hPtHatRaw->Add(hPtHatTmp);
    delete hPtHatTmp;
  }

  if(printDebug)cout<<"Filling PbPb MC"<<endl;

  for(int h = 0;h<nbins_pthat;h++){

    if(xsection[h]==0) continue;
    if(printDebug)cout <<"Loading pthat"<<boundaries_pthat[h]<<" sample, cross section = "<<xsection[h]<< Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[h],boundaries_pthat[h+1])<<endl;

    TEventList *el = new TEventList("el","el");

    double pthat_upper = boundaries_pthat[h+1];
    stringstream selection; selection<<"pthat<"<<pthat_upper;

    // get the number of events which are between the pthat ranges. 
    PbPbMC[h]->tJet->Draw(">>el",selection.str().c_str());
    double fentries = el->GetN();
    if(printDebug)cout<<"tree entries: "<<PbPbMC[h]->tJet->GetEntries()<<" elist: "<<fentries<<endl;
    delete el;
    int test_counter = 0; 
    Int_t nEntries = PbPbMC[h]->tJet->GetEntries();

    for(Long64_t nentry = 0;nentry < nEntries;nentry++){
 
      PbPbMC[h]->tEvt->GetEntry(nentry);
      PbPbMC[h]->tJet->GetEntry(nentry);
      PbPbMC[h]->tSkim->GetEntry(nentry);
      PbPbMC[h]->tHlt->GetEntry(nentry);
     
      int pthatBin = hPtHat->FindBin(PbPbMC[h]->pthat);
      double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/fentries;

      if(!PbPbMC[h]->pcollisionEventSelection) continue;
      int cBin = findBin(PbPbMC[h]->bin);
      //int cBin = nbins_cent-1;
      double weight_cent=1;
      double weight_vz=1;
      Float_t cent = cBin;
      
      weight_cent = hCentWeight->GetBinContent(hCentWeight->FindBin(PbPbMC[h]->bin));
      if(fabs(PbPbMC[h]->vz)>15) continue;
      weight_vz = fVz->Eval(PbPbMC[h]->vz);

      if(scale*weight_cent*weight_vz <=0 ) {
	cout<<"RED FLAG RED FLAG RED FLAG"<<endl;
	cout<<"pthat file = "<<boundaries_pthat[h]<<endl;
	cout<<"scale = "<<scale<<endl;
	cout<<"weight_cent = "<<weight_cent<<endl;
	cout<<"weight_vz = "<<weight_vz<<endl;
	continue;
      }

      for(int jentry = 0; jentry < PbPbMC[h]->njets; jentry++){

	// Make sure that we are getting the embedded events 
	if ( PbPbMC[h]->subid[jentry] != sub_id ) continue;
	if ( PbPbMC[h]->jtpt[jentry] > 2.*PbPbMC[h]->pthat) continue;

	// any other Jet ID/quality cuts should be placed here: 
	// eta selection here if youd like:
	if( fabs(PbPbMC[h]->jteta[jentry]) > 2 ) continue;

	hJetpT->Fill(PbPbMC[h]->jtpt[jentry],scale*weight_vz*weight_cent);

      }// jet loop

    }// entry loop
    
  }// pthat loop

  TFile fout(Form("jet_spectra_MC_ak%s%d%s_%d.root",algo,radius,jet_type,date.GetDate()),"RECREATE");
  fout.cd();
  hJetpT = (TH1F*)hJetpT->Rebin(nbins_pt,"hJetPT",boundaries_pt);
  divideBinWidth(hJetpT);
  hJetpT->Write();
  fout.Close();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
    

}// macro ends
