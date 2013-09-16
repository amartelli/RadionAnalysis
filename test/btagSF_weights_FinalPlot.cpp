
#include "PUreweightingUtils.h"
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "setTDRStyle.h"
#include "drawPlotsUtils.h"
#include "DiJetKinFitter.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THStack.h"
#include "TPad.h"
#include "TEfficiency.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TVector3.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

float nvtx;
float j1_eta;
float j2_eta;
float j3_eta;
float j4_eta;
float j1_betaStarClassic;
float j2_betaStarClassic;
float j3_betaStarClassic;
float j4_betaStarClassic;
float j1_dR2Mean;
float j2_dR2Mean;
float j3_dR2Mean;
float j4_dR2Mean;

bool passCutBasedJetId(int jet);
float weight_err(float jet1, float jet1_err,float jet2, float jet2_err,float weight, float rho);

int main(int argc, char** argv)
{

  // Set style options
  setTDRStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 
  //gStyle->SetOptStat(0000); 
  gStyle->SetOptFit(1); 

  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> Analyses::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");
  std::string inputWeightRoot = gConfigParser -> readStringOption("Input::inputWeightRoot");
  std::string inputTree = gConfigParser -> readStringOption("Input::inputTree");
  bool useMVA = gConfigParser -> readBoolOption("Input::useMVA");
 
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");
  
  std::map<int,TChain*> ntu;
  std::map<int, TFile*> Files;

  TChain* ntu_weight;
  ntu_weight= new TChain(inputTree.c_str());
  ntu_weight -> Add(inputWeightRoot.c_str());

  char trees[500];
  FILE *f_trees;
 
  f_trees = fopen(inputList.c_str(),"r");
  int pos = 0;
  int pos_total = 0;
  
  while(fscanf(f_trees,"%s \n", trees) !=EOF ){
    std::string TREES = std::string(trees);
 
    if(TREES.find("#") != std::string::npos) continue;
    std::cout << "\nReading input: " << trees << std::endl;
    ntu[pos] = new TChain(inputTree.c_str());
    ntu[pos] -> Add(trees);
   
    pos++;
    
  }

  pos_total = pos;
  
  for(int ii = 0; ii < pos_total; ii++){
      if(ntu[ii]->GetEntries() == 0 )
      {
         std::cout << "Error: input file" << ii << " is empty" << std::endl; 
         return -1;
      }
  }

  
  // vectors for jets
  float ptcorrjet[4], ecorrjet[4],ptuncorrjet[4];
  float etajet[4], phijet[4];              
  float btagjprobjet[4], btagcsvjet[4];
  float cefjet[4], chfjet[4];
  float vtxPtjet[4], vtx3dljet[4];
  float nconstjet[4];
  
  float fRegr_pt, fRegr_eta, fRegr_cef, fRegr_chf;
  float fRegr_nconst;
  float fRegr_vtxPt, fRegr_vtx3dl;
  float fRegr_met, fRegr_dPhiMet;
  
  TMVA::Reader *readerRegres = new TMVA::Reader( "!Color:!Silent" );
  readerRegres->AddVariable( "jet_pt",            &fRegr_pt);
  readerRegres->AddVariable( "jet_eta",           &fRegr_eta);
  readerRegres->AddVariable( "jet_emfrac",        &fRegr_cef);
  readerRegres->AddVariable( "jet_nConstituents", &fRegr_nconst);
  readerRegres->AddVariable( "jet_hadfrac",       &fRegr_chf);
  readerRegres->AddVariable( "jet_secVtxPt",      &fRegr_vtxPt);
  readerRegres->AddVariable( "jet_secVtx3dL",     &fRegr_vtx3dl);
  readerRegres->AddVariable( "ev_met_corr_pfmet",  &fRegr_met);
  readerRegres->AddVariable( "jet_dPhiMet",        &fRegr_dPhiMet);
  readerRegres->BookMVA("BDT","data/factoryJetRegGen2_globeinputs_BDT.weights.xml");
  
  float         event;
  float         lumis;
  float         evweight;
  float         weight;
  float         pu_weight;
  float         met_corr_pfmet;
  float         met_corr_phi_pfmet;
  float         met_corr_eta_pfmet;
  float         met_corr_e_pfmet;
  float         ph1_SCEta;
  float         ph1_pt;
  float         ph1_e;
  float         ph1_eta;
  float         ph1_phi;
  int           ph1_ciclevel;
  float         ph2_SCEta;
  float         ph2_pt;
  float         ph2_e;
  float         ph2_eta;
  float         ph2_phi;
  int           ph2_ciclevel;
  float         PhotonsMass;
  float         j1_e;
  float         j1_pt;
  float         j1_phi;
  float         j1_flavour;
  float         j1_weight;
  float         j1_weight_min;
  float         j1_weight_max;
  float         j1_weight_05min;
  float         j1_weight_05max;
  float         j1_weight_15min;
  float         j1_weight_15max;
  float         j1_weight_20min;
  float         j1_weight_20max;
  float         j1_csvBtag;
  float         j1_jetProbBtag; 
  float         j1_emfrac;
  float         j1_hadfrac;
  float         j1_secVtxPt;
  float         j1_secVtx3dL;
  int           j1_nNeutrals;
  int           j1_nCharged;
  float         j2_e;
  float         j2_pt;
  float         j2_phi;
  float         j2_flavour;
  float         j2_weight;
  float         j2_weight_min;
  float         j2_weight_max;
  float         j2_weight_05min;
  float         j2_weight_05max;
  float         j2_weight_15min;
  float         j2_weight_15max;
  float         j2_weight_20min;
  float         j2_weight_20max;
  float         j2_csvBtag;
  float         j2_jetProbBtag;
  float         j2_emfrac;
  float         j2_secVtxPt;
  float         j2_secVtx3dL;
  int           j2_nNeutrals;
  int           j2_nCharged;
  float         j2_hadfrac;
  float         j3_e;
  float         j3_pt;
  float         j3_phi;
  float         j3_flavour;
  float         j3_weight;
  float         j3_weight_min;
  float         j3_weight_max;
  float         j3_weight_05min;
  float         j3_weight_05max;
  float         j3_weight_15min;
  float         j3_weight_15max;
  float         j3_weight_20min;
  float         j3_weight_20max;
  float         j3_csvBtag;
  float         j3_jetProbBtag;
  float         j3_emfrac;
  float         j3_hadfrac;
  float         j3_secVtxPt;
  float         j3_secVtx3dL;
  int           j3_nNeutrals;
  int           j3_nCharged;
  float         j4_e;
  float         j4_pt;
  float         j4_phi;
  float         j4_flavour;
  float         j4_weight;
  float         j4_weight_min;
  float         j4_weight_max;
  float         j4_weight_05min;
  float         j4_weight_05max;
  float         j4_weight_15min;
  float         j4_weight_15max;
  float         j4_weight_20min;
  float         j4_weight_20max;
  float         j4_csvBtag;
  float         j4_jetProbBtag;
  float         j4_emfrac;
  float         j4_hadfrac;
  float         j4_secVtxPt;
  float         j4_secVtx3dL;
  int           j4_nNeutrals;
  int           j4_nCharged;
  float         dipho_cosThetaStar_CS;
  float         dipho_E;
  float         dipho_pt;
  float         dipho_eta;
  float         dipho_phi;
  int           njets_passing_kLooseID;

  
  ntu_weight->SetBranchStatus("*",0);
  ntu_weight->SetBranchStatus("event",1);
  ntu_weight->SetBranchStatus("lumis",1); 
  ntu_weight->SetBranchStatus("j1_phi",1);
  ntu_weight->SetBranchStatus("j1_eta",1);
  ntu_weight->SetBranchStatus("j1_pt",1);
  ntu_weight->SetBranchStatus("j1_e",1);
  ntu_weight->SetBranchStatus("j1_flavour",1);
  ntu_weight->SetBranchStatus("j1_weight",1);
  ntu_weight->SetBranchStatus("j1_weight_max",1);
  ntu_weight->SetBranchStatus("j1_weight_min",1);
  ntu_weight->SetBranchStatus("j1_weight_05max",1);
  ntu_weight->SetBranchStatus("j1_weight_05min",1);
  ntu_weight->SetBranchStatus("j1_weight_15max",1);
  ntu_weight->SetBranchStatus("j1_weight_15min",1);
  ntu_weight->SetBranchStatus("j1_weight_20max",1);
  ntu_weight->SetBranchStatus("j1_weight_20min",1);
  ntu_weight->SetBranchStatus("j2_phi",1);
  ntu_weight->SetBranchStatus("j2_eta",1);
  ntu_weight->SetBranchStatus("j2_pt",1);
  ntu_weight->SetBranchStatus("j2_e",1);
  ntu_weight->SetBranchStatus("j2_flavour",1);
  ntu_weight->SetBranchStatus("j2_weight",1);
  ntu_weight->SetBranchStatus("j2_weight_max",1);
  ntu_weight->SetBranchStatus("j2_weight_min",1);
  ntu_weight->SetBranchStatus("j2_weight_05max",1);
  ntu_weight->SetBranchStatus("j2_weight_05min",1);
  ntu_weight->SetBranchStatus("j2_weight_15max",1);
  ntu_weight->SetBranchStatus("j2_weight_15min",1);
  ntu_weight->SetBranchStatus("j2_weight_20max",1);
  ntu_weight->SetBranchStatus("j2_weight_20min",1);
  ntu_weight->SetBranchStatus("j3_phi",1);
  ntu_weight->SetBranchStatus("j3_eta",1);
  ntu_weight->SetBranchStatus("j3_pt",1);
  ntu_weight->SetBranchStatus("j3_e",1);
  ntu_weight->SetBranchStatus("j3_flavour",1);
  ntu_weight->SetBranchStatus("j3_weight",1);
  ntu_weight->SetBranchStatus("j3_weight_max",1);
  ntu_weight->SetBranchStatus("j3_weight_min",1);
  ntu_weight->SetBranchStatus("j3_weight_05max",1);
  ntu_weight->SetBranchStatus("j3_weight_05min",1);
  ntu_weight->SetBranchStatus("j3_weight_15max",1);
  ntu_weight->SetBranchStatus("j3_weight_15min",1);
  ntu_weight->SetBranchStatus("j3_weight_20max",1);
  ntu_weight->SetBranchStatus("j3_weight_20min",1);
  ntu_weight->SetBranchStatus("j4_phi",1);
  ntu_weight->SetBranchStatus("j4_eta",1);
  ntu_weight->SetBranchStatus("j4_pt",1);
  ntu_weight->SetBranchStatus("j4_e",1);
  ntu_weight->SetBranchStatus("j4_flavour",1);
  ntu_weight->SetBranchStatus("j4_weight",1);
  ntu_weight->SetBranchStatus("j4_weight_max",1);
  ntu_weight->SetBranchStatus("j4_weight_min",1);
  ntu_weight->SetBranchStatus("j4_weight_05max",1);
  ntu_weight->SetBranchStatus("j4_weight_05min",1);
  ntu_weight->SetBranchStatus("j4_weight_15max",1);
  ntu_weight->SetBranchStatus("j4_weight_15min",1);
  ntu_weight->SetBranchStatus("j4_weight_20max",1);
  ntu_weight->SetBranchStatus("j4_weight_20min",1);
  ntu_weight->SetBranchAddress("event",&event);
  ntu_weight->SetBranchAddress("lumis",&lumis);
  ntu_weight->SetBranchAddress("j1_phi",&j1_phi);
  ntu_weight->SetBranchAddress("j1_eta",&j1_eta);
  ntu_weight->SetBranchAddress("j1_pt",&j1_pt);
  ntu_weight->SetBranchAddress("j1_e",&j1_e);
  ntu_weight->SetBranchAddress("j1_flavour",&j1_flavour);
  ntu_weight->SetBranchAddress("j1_weight",&j1_weight);
  ntu_weight->SetBranchAddress("j1_weight_max",&j1_weight_max);
  ntu_weight->SetBranchAddress("j1_weight_min",&j1_weight_min);
  ntu_weight->SetBranchAddress("j1_weight_05max",&j1_weight_05max);
  ntu_weight->SetBranchAddress("j1_weight_05min",&j1_weight_05min);
  ntu_weight->SetBranchAddress("j1_weight_15max",&j1_weight_15max);
  ntu_weight->SetBranchAddress("j1_weight_15min",&j1_weight_15min);
  ntu_weight->SetBranchAddress("j1_weight_20max",&j1_weight_20max);
  ntu_weight->SetBranchAddress("j1_weight_20min",&j1_weight_20min);
  ntu_weight->SetBranchAddress("j2_phi",&j2_phi);
  ntu_weight->SetBranchAddress("j2_eta",&j2_eta);
  ntu_weight->SetBranchAddress("j2_pt",&j2_pt);
  ntu_weight->SetBranchAddress("j2_e",&j2_e);
  ntu_weight->SetBranchAddress("j2_flavour",&j2_flavour);
  ntu_weight->SetBranchAddress("j2_weight",&j2_weight);
  ntu_weight->SetBranchAddress("j2_weight_max",&j2_weight_max);
  ntu_weight->SetBranchAddress("j2_weight_min",&j2_weight_min);
  ntu_weight->SetBranchAddress("j2_weight_05max",&j2_weight_05max);
  ntu_weight->SetBranchAddress("j2_weight_05min",&j2_weight_05min);
  ntu_weight->SetBranchAddress("j2_weight_15max",&j2_weight_15max);
  ntu_weight->SetBranchAddress("j2_weight_15min",&j2_weight_15min);
  ntu_weight->SetBranchAddress("j2_weight_20max",&j2_weight_20max);
  ntu_weight->SetBranchAddress("j2_weight_20min",&j2_weight_20min);
  ntu_weight->SetBranchAddress("j3_phi",&j3_phi);
  ntu_weight->SetBranchAddress("j3_eta",&j3_eta);
  ntu_weight->SetBranchAddress("j3_pt",&j3_pt);
  ntu_weight->SetBranchAddress("j3_e",&j3_e);
  ntu_weight->SetBranchAddress("j3_flavour",&j3_flavour);
  ntu_weight->SetBranchAddress("j3_weight",&j3_weight);
  ntu_weight->SetBranchAddress("j3_weight_max",&j3_weight_max);
  ntu_weight->SetBranchAddress("j3_weight_min",&j3_weight_min);
  ntu_weight->SetBranchAddress("j3_weight_05max",&j3_weight_05max);
  ntu_weight->SetBranchAddress("j3_weight_05min",&j3_weight_05min);
  ntu_weight->SetBranchAddress("j3_weight_15max",&j3_weight_15max);
  ntu_weight->SetBranchAddress("j3_weight_15min",&j3_weight_15min);
  ntu_weight->SetBranchAddress("j3_weight_20max",&j3_weight_20max);
  ntu_weight->SetBranchAddress("j3_weight_20min",&j3_weight_20min);
  ntu_weight->SetBranchAddress("j4_phi",&j4_phi);
  ntu_weight->SetBranchAddress("j4_eta",&j4_eta);
  ntu_weight->SetBranchAddress("j4_pt",&j4_pt);
  ntu_weight->SetBranchAddress("j4_e",&j4_e);
  ntu_weight->SetBranchAddress("j4_flavour",&j4_flavour);
  ntu_weight->SetBranchAddress("j4_weight",&j4_weight);
  ntu_weight->SetBranchAddress("j4_weight_max",&j4_weight_max);
  ntu_weight->SetBranchAddress("j4_weight_min",&j4_weight_min);
  ntu_weight->SetBranchAddress("j4_weight_05max",&j4_weight_05max);
  ntu_weight->SetBranchAddress("j4_weight_05min",&j4_weight_05min);
  ntu_weight->SetBranchAddress("j4_weight_15max",&j4_weight_15max);
  ntu_weight->SetBranchAddress("j4_weight_15min",&j4_weight_15min);
  ntu_weight->SetBranchAddress("j4_weight_20max",&j4_weight_20max);
  ntu_weight->SetBranchAddress("j4_weight_20min",&j4_weight_20min);
  
  TLorentzVector p4;
  
  std::map<float,std::map<float,std::vector<TLorentzVector> > > PREjet_p4;
  std::map<float,std::map<float,std::vector<bool> > > PREbad_jet;
  std::map<float,std::map<float,std::vector<float> > > PREjet_flavour;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_max;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_min;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_05max;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_05min;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_15max;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_15min;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_20max;
  std::map<float,std::map<float,std::vector<float> > > PREjet_weight_20min;
  
  for(int ientry = 0; ientry < ntu_weight->GetEntries(); ientry++){
      ntu_weight->GetEntry(ientry);

      if(j1_pt < 0) PREbad_jet[lumis][event].push_back(true);
      if(j1_pt > 0) PREbad_jet[lumis][event].push_back(false);
      if(j2_pt < 0) PREbad_jet[lumis][event].push_back(true);
      if(j2_pt > 0) PREbad_jet[lumis][event].push_back(false);
      if(j3_pt < 0) PREbad_jet[lumis][event].push_back(true);
      if(j3_pt > 0) PREbad_jet[lumis][event].push_back(false);
      if(j4_pt < 0) PREbad_jet[lumis][event].push_back(true);
      if(j4_pt > 0) PREbad_jet[lumis][event].push_back(false);

      p4.SetPtEtaPhiE(j1_pt,j1_eta,j1_phi,j1_e);
      PREjet_p4[lumis][event].push_back(p4);
      p4.SetPtEtaPhiE(j2_pt,j2_eta,j2_phi,j2_e);
      PREjet_p4[lumis][event].push_back(p4);
      p4.SetPtEtaPhiE(j3_pt,j3_eta,j3_phi,j3_e);
      PREjet_p4[lumis][event].push_back(p4);
      p4.SetPtEtaPhiE(j4_pt,j4_eta,j4_phi,j4_e);
      PREjet_p4[lumis][event].push_back(p4);
      PREjet_flavour[lumis][event].push_back(j1_flavour);
      PREjet_flavour[lumis][event].push_back(j2_flavour);
      PREjet_flavour[lumis][event].push_back(j3_flavour);
      PREjet_flavour[lumis][event].push_back(j4_flavour);
      PREjet_weight[lumis][event].push_back(j1_weight);
      PREjet_weight[lumis][event].push_back(j2_weight);
      PREjet_weight[lumis][event].push_back(j3_weight);
      PREjet_weight[lumis][event].push_back(j4_weight);
      PREjet_weight_max[lumis][event].push_back(j1_weight_max);
      PREjet_weight_max[lumis][event].push_back(j2_weight_max);
      PREjet_weight_max[lumis][event].push_back(j3_weight_max);
      PREjet_weight_max[lumis][event].push_back(j4_weight_max);
      PREjet_weight_min[lumis][event].push_back(j1_weight_min);
      PREjet_weight_min[lumis][event].push_back(j2_weight_min);
      PREjet_weight_min[lumis][event].push_back(j3_weight_min);
      PREjet_weight_min[lumis][event].push_back(j4_weight_min);
      PREjet_weight_05max[lumis][event].push_back(j1_weight_05max);
      PREjet_weight_05max[lumis][event].push_back(j2_weight_05max);
      PREjet_weight_05max[lumis][event].push_back(j3_weight_05max);
      PREjet_weight_05max[lumis][event].push_back(j4_weight_05max);
      PREjet_weight_05min[lumis][event].push_back(j1_weight_05min);
      PREjet_weight_05min[lumis][event].push_back(j2_weight_05min);
      PREjet_weight_05min[lumis][event].push_back(j3_weight_05min);
      PREjet_weight_05min[lumis][event].push_back(j4_weight_05min);
      PREjet_weight_15max[lumis][event].push_back(j1_weight_15max);
      PREjet_weight_15max[lumis][event].push_back(j2_weight_15max);
      PREjet_weight_15max[lumis][event].push_back(j3_weight_15max);
      PREjet_weight_15max[lumis][event].push_back(j4_weight_15max);
      PREjet_weight_15min[lumis][event].push_back(j1_weight_15min);
      PREjet_weight_15min[lumis][event].push_back(j2_weight_15min);
      PREjet_weight_15min[lumis][event].push_back(j3_weight_15min);
      PREjet_weight_15min[lumis][event].push_back(j4_weight_15min);
      PREjet_weight_20max[lumis][event].push_back(j1_weight_20max);
      PREjet_weight_20max[lumis][event].push_back(j2_weight_20max);
      PREjet_weight_20max[lumis][event].push_back(j3_weight_20max);
      PREjet_weight_20max[lumis][event].push_back(j4_weight_20max);
      PREjet_weight_20min[lumis][event].push_back(j1_weight_20min);
      PREjet_weight_20min[lumis][event].push_back(j2_weight_20min);
      PREjet_weight_20min[lumis][event].push_back(j3_weight_20min);
      PREjet_weight_20min[lumis][event].push_back(j4_weight_20min);
      
     
  }
  
  for(int ii = 0; ii < pos_total; ii++){
      ntu[ii]->SetBranchStatus("*",0);
      ntu[ii]->SetBranchStatus("event",1);
      ntu[ii]->SetBranchStatus("lumis",1);
      ntu[ii]->SetBranchStatus("nvtx",1);
      ntu[ii]->SetBranchStatus("evweight",1);
      ntu[ii]->SetBranchStatus("weight",1);
      ntu[ii]->SetBranchStatus("pu_weight",1);
      ntu[ii]->SetBranchStatus("met_corr_pfmet",1);
      ntu[ii]->SetBranchStatus("met_corr_phi_pfmet",1);
      ntu[ii]->SetBranchStatus("met_corr_eta_pfmet",1);
      ntu[ii]->SetBranchStatus("met_corr_e_pfmet",1);
      ntu[ii]->SetBranchStatus("ph1_pt",1);
      ntu[ii]->SetBranchStatus("ph1_eta",1);
      ntu[ii]->SetBranchStatus("ph1_phi",1);
      ntu[ii]->SetBranchStatus("ph1_e",1);
      ntu[ii]->SetBranchStatus("ph1_SCEta",1);
      ntu[ii]->SetBranchStatus("ph1_ciclevel",1);
      ntu[ii]->SetBranchStatus("ph2_pt",1);
      ntu[ii]->SetBranchStatus("ph2_pt",1);
      ntu[ii]->SetBranchStatus("ph2_eta",1);
      ntu[ii]->SetBranchStatus("ph2_phi",1);
      ntu[ii]->SetBranchStatus("ph2_e",1);
      ntu[ii]->SetBranchStatus("ph2_SCEta",1);
      ntu[ii]->SetBranchStatus("ph2_ciclevel",1);
      ntu[ii]->SetBranchStatus("PhotonsMass",1);
      ntu[ii]->SetBranchStatus("PhotonsMass",1);
      ntu[ii]->SetBranchStatus("j1_e",1);
      ntu[ii]->SetBranchStatus("j1_pt",1);
      ntu[ii]->SetBranchStatus("j1_phi",1);
      ntu[ii]->SetBranchStatus("j1_eta",1);
      ntu[ii]->SetBranchStatus("j1_csvBtag",1);
      ntu[ii]->SetBranchStatus("j1_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j1_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j1_emfrac",1);
      ntu[ii]->SetBranchStatus("j1_hadfrac",1);
      ntu[ii]->SetBranchStatus("j1_secVtxPt",1);
      ntu[ii]->SetBranchStatus("j1_secVtx3dL",1);
      ntu[ii]->SetBranchStatus("j1_nNeutrals",1);
      ntu[ii]->SetBranchStatus("j1_nCharged",1);
      ntu[ii]->SetBranchStatus("j1_dR2Mean",1);
      ntu[ii]->SetBranchStatus("j2_e",1);
      ntu[ii]->SetBranchStatus("j2_pt",1);
      ntu[ii]->SetBranchStatus("j2_phi",1);
      ntu[ii]->SetBranchStatus("j2_eta",1);
      ntu[ii]->SetBranchStatus("j2_csvBtag",1);
      ntu[ii]->SetBranchStatus("j2_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j2_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j2_emfrac",1);
      ntu[ii]->SetBranchStatus("j2_hadfrac",1);
      ntu[ii]->SetBranchStatus("j2_secVtxPt",1);
      ntu[ii]->SetBranchStatus("j2_secVtx3dL",1);
      ntu[ii]->SetBranchStatus("j2_nNeutrals",1);
      ntu[ii]->SetBranchStatus("j2_nCharged",1);
      ntu[ii]->SetBranchStatus("j2_dR2Mean",1);
      ntu[ii]->SetBranchStatus("j3_e",1);
      ntu[ii]->SetBranchStatus("j3_pt",1);
      ntu[ii]->SetBranchStatus("j3_phi",1);
      ntu[ii]->SetBranchStatus("j3_eta",1);
      ntu[ii]->SetBranchStatus("j3_csvBtag",1);
      ntu[ii]->SetBranchStatus("j3_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j3_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j3_emfrac",1);
      ntu[ii]->SetBranchStatus("j3_hadfrac",1);
      ntu[ii]->SetBranchStatus("j3_secVtxPt",1);
      ntu[ii]->SetBranchStatus("j3_secVtx3dL",1);
      ntu[ii]->SetBranchStatus("j3_nNeutrals",1);
      ntu[ii]->SetBranchStatus("j3_nCharged",1);
      ntu[ii]->SetBranchStatus("j3_dR2Mean",1);
      ntu[ii]->SetBranchStatus("j4_e",1);
      ntu[ii]->SetBranchStatus("j4_pt",1);
      ntu[ii]->SetBranchStatus("j4_phi",1);
      ntu[ii]->SetBranchStatus("j4_eta",1);
      ntu[ii]->SetBranchStatus("j4_csvBtag",1);
      ntu[ii]->SetBranchStatus("j4_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j4_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j4_emfrac",1);
      ntu[ii]->SetBranchStatus("j4_hadfrac",1);
      ntu[ii]->SetBranchStatus("j4_secVtxPt",1);
      ntu[ii]->SetBranchStatus("j4_secVtx3dL",1);
      ntu[ii]->SetBranchStatus("j4_nNeutrals",1);
      ntu[ii]->SetBranchStatus("j4_nCharged",1);
      ntu[ii]->SetBranchStatus("j4_dR2Mean",1);
      ntu[ii]->SetBranchStatus("dipho_cosThetaStar_CS",1);
      ntu[ii]->SetBranchStatus("dipho_E",1);
      ntu[ii]->SetBranchStatus("dipho_pt",1);
      ntu[ii]->SetBranchStatus("dipho_eta",1);
      ntu[ii]->SetBranchStatus("dipho_phi",1);
      ntu[ii]->SetBranchStatus("njets_passing_kLooseID",1);
      ntu[ii]->SetBranchAddress("event",&event);
      ntu[ii]->SetBranchAddress("lumis",&lumis);
      ntu[ii]->SetBranchAddress("nvtx",&nvtx);
      ntu[ii]->SetBranchAddress("evweight",&evweight);
      ntu[ii]->SetBranchAddress("weight",&weight);
      ntu[ii]->SetBranchAddress("pu_weight",&pu_weight);
      ntu[ii]->SetBranchAddress("met_corr_pfmet",&met_corr_pfmet);
      ntu[ii]->SetBranchAddress("met_corr_phi_pfmet",&met_corr_phi_pfmet);
      ntu[ii]->SetBranchAddress("met_corr_eta_pfmet",&met_corr_eta_pfmet);
      ntu[ii]->SetBranchAddress("met_corr_e_pfmet",&met_corr_e_pfmet);
      ntu[ii]->SetBranchAddress("ph1_pt",&ph1_pt);
      ntu[ii]->SetBranchAddress("ph1_eta",&ph1_eta);
      ntu[ii]->SetBranchAddress("ph1_phi",&ph1_phi);
      ntu[ii]->SetBranchAddress("ph1_e",&ph1_e);
      ntu[ii]->SetBranchAddress("ph1_SCEta",&ph1_SCEta);
      ntu[ii]->SetBranchAddress("ph1_ciclevel",&ph1_ciclevel);
      ntu[ii]->SetBranchAddress("ph2_pt",&ph2_pt);
      ntu[ii]->SetBranchAddress("ph2_eta",&ph2_eta);
      ntu[ii]->SetBranchAddress("ph2_phi",&ph2_phi);
      ntu[ii]->SetBranchAddress("ph2_e",&ph2_e);
      ntu[ii]->SetBranchAddress("ph2_SCEta",&ph2_SCEta);
      ntu[ii]->SetBranchAddress("ph2_ciclevel",&ph2_ciclevel);
      ntu[ii]->SetBranchAddress("PhotonsMass",&PhotonsMass);
      ntu[ii]->SetBranchAddress("j1_e",&j1_e);
      ntu[ii]->SetBranchAddress("j1_pt",&j1_pt);
      ntu[ii]->SetBranchAddress("j1_phi",&j1_phi);
      ntu[ii]->SetBranchAddress("j1_eta",&j1_eta);
      ntu[ii]->SetBranchAddress("j1_csvBtag",&j1_csvBtag);
      ntu[ii]->SetBranchAddress("j1_betaStarClassic",&j1_betaStarClassic);
      ntu[ii]->SetBranchAddress("j1_jetProbBtag",&j1_jetProbBtag);
      ntu[ii]->SetBranchAddress("j1_emfrac",&j1_emfrac);
      ntu[ii]->SetBranchAddress("j1_hadfrac",&j1_hadfrac);
      ntu[ii]->SetBranchAddress("j1_secVtxPt",&j1_secVtxPt);
      ntu[ii]->SetBranchAddress("j1_secVtx3dL",&j1_secVtx3dL);
      ntu[ii]->SetBranchAddress("j1_nNeutrals",&j1_nNeutrals);
      ntu[ii]->SetBranchAddress("j1_nCharged",&j1_nCharged);
      ntu[ii]->SetBranchAddress("j1_dR2Mean",&j1_dR2Mean);
      ntu[ii]->SetBranchAddress("j2_e",&j2_e);
      ntu[ii]->SetBranchAddress("j2_pt",&j2_pt);
      ntu[ii]->SetBranchAddress("j2_phi",&j2_phi);
      ntu[ii]->SetBranchAddress("j2_eta",&j2_eta);
      ntu[ii]->SetBranchAddress("j2_csvBtag",&j2_csvBtag);
      ntu[ii]->SetBranchAddress("j2_betaStarClassic",&j2_betaStarClassic);
      ntu[ii]->SetBranchAddress("j2_jetProbBtag",&j2_jetProbBtag);
      ntu[ii]->SetBranchAddress("j2_emfrac",&j2_emfrac);
      ntu[ii]->SetBranchAddress("j2_hadfrac",&j2_hadfrac);
      ntu[ii]->SetBranchAddress("j2_secVtxPt",&j2_secVtxPt);
      ntu[ii]->SetBranchAddress("j2_secVtx3dL",&j2_secVtx3dL);
      ntu[ii]->SetBranchAddress("j2_nNeutrals",&j2_nNeutrals);
      ntu[ii]->SetBranchAddress("j2_nCharged",&j2_nCharged);
      ntu[ii]->SetBranchAddress("j2_dR2Mean",&j2_dR2Mean);
      ntu[ii]->SetBranchAddress("j3_e",&j3_e);
      ntu[ii]->SetBranchAddress("j3_pt",&j3_pt);
      ntu[ii]->SetBranchAddress("j3_phi",&j3_phi);
      ntu[ii]->SetBranchAddress("j3_eta",&j3_eta);
      ntu[ii]->SetBranchAddress("j3_csvBtag",&j3_csvBtag);
      ntu[ii]->SetBranchAddress("j3_betaStarClassic",&j3_betaStarClassic);
      ntu[ii]->SetBranchAddress("j3_jetProbBtag",&j3_jetProbBtag);
      ntu[ii]->SetBranchAddress("j3_emfrac",&j3_emfrac);
      ntu[ii]->SetBranchAddress("j3_hadfrac",&j3_hadfrac);
      ntu[ii]->SetBranchAddress("j3_secVtxPt",&j3_secVtxPt);
      ntu[ii]->SetBranchAddress("j3_secVtx3dL",&j3_secVtx3dL);
      ntu[ii]->SetBranchAddress("j3_nNeutrals",&j3_nNeutrals);
      ntu[ii]->SetBranchAddress("j3_nCharged",&j3_nCharged);
      ntu[ii]->SetBranchAddress("j3_dR2Mean",&j3_dR2Mean);
      ntu[ii]->SetBranchAddress("j4_e",&j4_e);
      ntu[ii]->SetBranchAddress("j4_pt",&j4_pt);
      ntu[ii]->SetBranchAddress("j4_phi",&j4_phi);
      ntu[ii]->SetBranchAddress("j4_eta",&j4_eta);
      ntu[ii]->SetBranchAddress("j4_csvBtag",&j4_csvBtag);
      ntu[ii]->SetBranchAddress("j4_betaStarClassic",&j4_betaStarClassic);
      ntu[ii]->SetBranchAddress("j4_jetProbBtag",&j4_jetProbBtag);
      ntu[ii]->SetBranchAddress("j4_emfrac",&j4_emfrac);
      ntu[ii]->SetBranchAddress("j4_hadfrac",&j4_hadfrac);
      ntu[ii]->SetBranchAddress("j4_secVtxPt",&j4_secVtxPt);
      ntu[ii]->SetBranchAddress("j4_secVtx3dL",&j4_secVtx3dL);
      ntu[ii]->SetBranchAddress("j4_nNeutrals",&j4_nNeutrals);
      ntu[ii]->SetBranchAddress("j4_nCharged",&j4_nCharged);
      ntu[ii]->SetBranchAddress("j4_dR2Mean",&j4_dR2Mean);
      ntu[ii]->SetBranchAddress("dipho_cosThetaStar_CS",&dipho_cosThetaStar_CS);
      ntu[ii]->SetBranchAddress("dipho_E",&dipho_E);
      ntu[ii]->SetBranchAddress("dipho_pt",&dipho_pt);
      ntu[ii]->SetBranchAddress("dipho_eta",&dipho_eta);
      ntu[ii]->SetBranchAddress("dipho_phi",&dipho_phi);
      ntu[ii]->SetBranchAddress("njets_passing_kLooseID",&njets_passing_kLooseID);    
  }


  TH1F* h_btagSF_weights = new TH1F("h_btagSF_weights","h_btagSF_weights",200,0,2);
  TH1F* h_btagSF_weights_up = new TH1F("h_btagSF_weights_up","h_btagSF_weights_up",200,0,2);
  TH1F* h_btagSF_weights_down = new TH1F("h_btagSF_weights_down","h_btagSF_weights_down",200,0,2);
  TH1F* h_btagSF_weights_1b = new TH1F("h_btagSF_weights_1b","h_btagSF_weights_1b",200,0,2);
  TH1F* h_btagSF_weights_up_1b = new TH1F("h_btagSF_weights_up_1b","h_btagSF_weights_up_1b",200,0,2);
  TH1F* h_btagSF_weights_down_1b = new TH1F("h_btagSF_weights_down_1b","h_btagSF_weights_down_1b",200,0,2);
  TH1F* h_btagSF_weights_2b = new TH1F("h_btagSF_weights_2b","h_btagSF_weights_2b",200,0,2);
  TH1F* h_btagSF_weights_up_2b = new TH1F("h_btagSF_weights_up_2b","h_btagSF_weights_up_2b",200,0,2);
  TH1F* h_btagSF_weights_down_2b = new TH1F("h_btagSF_weights_down_2b","h_btagSF_weights_down_2b",200,0,2);
  TH1F* h_diPhodiJet_InvMass_old = new TH1F("h_diPhodiJet_InvMass_old","h_diPhodiJet_InvMass_old",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass = new TH1F("h_diPhodiJet_InvMass","h_diPhodiJet_InvMass",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_up = new TH1F("h_diPhodiJet_InvMass_up","h_diPhodiJet_InvMass_up",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_down = new TH1F("h_diPhodiJet_InvMass_down","h_diPhodiJet_InvMass_down",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_old_1btag = new TH1F("h_diPhodiJet_InvMass_old_1btag","h_diPhodiJet_InvMass_old_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_1btag = new TH1F("h_diPhodiJet_InvMass_1btag","h_diPhodiJet_InvMass_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_up_1btag = new TH1F("h_diPhodiJet_InvMass_up_1btag","h_diPhodiJet_InvMass_up_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_05up_1btag = new TH1F("h_diPhodiJet_InvMass_05up_1btag","h_diPhodiJet_InvMass_05up_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_15up_1btag = new TH1F("h_diPhodiJet_InvMass_15up_1btag","h_diPhodiJet_InvMass_15up_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_20up_1btag = new TH1F("h_diPhodiJet_InvMass_20up_1btag","h_diPhodiJet_InvMass_20up_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_down_1btag = new TH1F("h_diPhodiJet_InvMass_down_1btag","h_diPhodiJet_InvMass_down_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_05down_1btag = new TH1F("h_diPhodiJet_InvMass_05down_1btag","h_diPhodiJet_InvMass_05down_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_15down_1btag = new TH1F("h_diPhodiJet_InvMass_15down_1btag","h_diPhodiJet_InvMass_15down_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_20down_1btag = new TH1F("h_diPhodiJet_InvMass_20down_1btag","h_diPhodiJet_InvMass_20down_1btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_old_2btag = new TH1F("h_diPhodiJet_InvMass_old_2btag","h_diPhodiJet_InvMass_old_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_2btag = new TH1F("h_diPhodiJet_InvMass_2btag","h_diPhodiJet_InvMass_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_up_2btag = new TH1F("h_diPhodiJet_InvMass_up_2btag","h_diPhodiJet_InvMass_up_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_down_2btag = new TH1F("h_diPhodiJet_InvMass_down_2btag","h_diPhodiJet_InvMass_down_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_05up_2btag = new TH1F("h_diPhodiJet_InvMass_05up_2btag","h_diPhodiJet_InvMass_05up_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_15up_2btag = new TH1F("h_diPhodiJet_InvMass_15up_2btag","h_diPhodiJet_InvMass_15up_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_20up_2btag = new TH1F("h_diPhodiJet_InvMass_20up_2btag","h_diPhodiJet_InvMass_20up_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_05down_2btag = new TH1F("h_diPhodiJet_InvMass_05down_2btag","h_diPhodiJet_InvMass_05down_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_15down_2btag = new TH1F("h_diPhodiJet_InvMass_15down_2btag","h_diPhodiJet_InvMass_15down_2btag",100,0.,1000.);
  TH1F* h_diPhodiJet_InvMass_20down_2btag = new TH1F("h_diPhodiJet_InvMass_20down_2btag","h_diPhodiJet_InvMass_20down_2btag",100,0.,1000.);

  TH2F* h2_pt1pt2 = new TH2F("h2_pt1pt2","h2_pt1pt2",150,0.,150.,150,0.,150.);
  TH2F* h2_eta1eta2 = new TH2F("h2_eta1eta2","h2_eta1eta2",100,-2.5,2.5,100,-2.5,2.5);
  TH2F* h2_flav1flav2 = new TH2F("h2_flav1flav2","h2_flav1flav2",21,0.5,21.5,21,0.5,21.5);
  TH2F* h2_jw1jw2_1btag = new TH2F("h2_jw1jw2_1btag","h2_jw1jw2_1btag",200,0,2,200,0,2);
  TH2F* h2_jw1jw2_2btag = new TH2F("h2_jw1jw2_2btag","h2_jw1jw2_2btag",200,0,2,200,0,2);
  TH2F* h2_jw1jw2_1btag_min = new TH2F("h2_jw1jw2_1btag_min","h2_jw1jw2_1btag_min",200,0,2,200,0,2);
  TH2F* h2_jw1jw2_2btag_min = new TH2F("h2_jw1jw2_2btag_min","h2_jw1jw2_2btag_min",200,0,2,200,0,2);
  TH2F* h2_jw1jw2_1btag_max = new TH2F("h2_jw1jw2_1btag_max","h2_jw1jw2_1btag_max",200,0,2,200,0,2);
  TH2F* h2_jw1jw2_2btag_max = new TH2F("h2_jw1jw2_2btag_max","h2_jw1jw2_2btag_max",200,0,2,200,0,2);
  TH2F* h2_jwpt_1btag_b = new TH2F("h2_jwpt_1btag_b","h2_jwpt_1btag_b",200,0,2,150,0.,150.);
  TH2F* h2_jwpt_1btag_nob = new TH2F("h2_jwpt_1btag_nob","h2_jwpt_1btag_nob",200,0,2,150,0.,150.);
  TH2F* h2_jwpt_2btag_b1 = new TH2F("h2_jwpt_2btag_b1","h2_jwpt_2btag_b1",200,0,2,150,0.,150.);
  TH2F* h2_jwpt_2btag_b2 = new TH2F("h2_jwpt_2btag_b2","h2_jwpt_2btag_b2",200,0,2,150,0.,150.);
  TH2F* h2_jweta_1btag_b = new TH2F("h2_jweta_1btag_b","h2_jweta_1btag_b",200,0,2,100,-2.5,2.5);
  TH2F* h2_jweta_1btag_nob = new TH2F("h2_jweta_1btag_nob","h2_jweta_1btag_nob",200,0,2,100,-2.5,2.5);
  TH2F* h2_jweta_2btag_b1 = new TH2F("h2_jweta_2btag_b1","h2_jweta_2btag_b1",200,0,2,100,-2.5,2.5);
  TH2F* h2_jweta_2btag_b2 = new TH2F("h2_jweta_2btag_b2","h2_jweta_2btag_b2",200,0,2,100,-2.5,2.5);
  TH2F* h2_jwflav_1btag_b = new TH2F("h2_jwflav_1btag_b","h2_jwflav_1btag_b",200,0,2,21,0.5,21.5);
  TH2F* h2_jwflav_1btag_nob = new TH2F("h2_jwflav_1btag_nob","h2_jwflav_1btag_nob",200,0,2,21,0.5,21.5);
  TH2F* h2_jwflav_1btag_b_min = new TH2F("h2_jwflav_1btag_b_min","h2_jwflav_1btag_b_min",200,0,2,21,0.5,21.5);
  TH2F* h2_jwflav_1btag_nob_min = new TH2F("h2_jwflav_1btag_nob_min","h2_jwflav_1btag_nob_min",200,0,2,21,0.5,21.5);
  TH2F* h2_jwflav_1btag_b_max = new TH2F("h2_jwflav_1btag_b_max","h2_jwflav_1btag_b_max",200,0,2,21,0.5,21.5);
  TH2F* h2_jwflav_1btag_nob_max = new TH2F("h2_jwflav_1btag_nob_max","h2_jwflav_1btag_nob_max",200,0,2,21,0.5,21.5);
  TH2F* h2_jwflav_2btag_b1 = new TH2F("h2_jwflav_2btag_b1","h2_jwflav_2btag_b1",200,0,2,21,0.5,21.5);
  TH2F* h2_jwflav_2btag_b2 = new TH2F("h2_jwflav_2btag_b2","h2_jwflav_2btag_b2",200,0,2,21,0.5,21.5);
  TH2F* h2_csvflav = new TH2F("h2_csvflav","h2_csvflav",200.,0.,1.,21,0.5,21.5);

  TGraphErrors* g_1btag2btag_correletion_up = new TGraphErrors();
  TGraphErrors* g_1btag2btag_correletion_down = new TGraphErrors();
  
  std::map<int,TLorentzVector*> jetP4;
  std::map<int,TLorentzVector*> phoP4;
  std::map<int,float> jet_csvBtag;
  std::map<int,bool> bad_jet;
  std::map<int,bool> bad_jet_pre;
  std::map<int,float> jet_weight;
  std::map<int,float> jet_weight_max;
  std::map<int,float> jet_weight_min;
  std::map<int,float> jet_weight_05max;
  std::map<int,float> jet_weight_05min;
  std::map<int,float> jet_weight_15max;
  std::map<int,float> jet_weight_15min;
  std::map<int,float> jet_weight_20max;
  std::map<int,float> jet_weight_20min;
  std::map<int,float> jet_flavour;

  int n_event = 0;
  int n_event_match = 0;
  int n_good_jets_total = 0;
  int n_flav_jets = 0;
  
  std::vector<float> dR;
  std::map<float,int> dR_map;

  TLorentzVector* sumP4;
  std::vector<float> jjPt;
  std::map<float,std::pair<int,int> > jjPt_map;

  int n_out = 0;

  TLorentzVector* Hgg;
  TLorentzVector* Hjj; 
  TLorentzVector* Radion;
  //TLorentzVector* Hgg_Rstar;

  int n_1btag_out = 0;
  int n_1btag_in = 0;
  int n_1btag_out_c = 0;
  int n_1btag_in_c = 0;
  int n_2btag_out = 0;
  int n_2btag_in = 0;

  float Hmass = 125.;
  DiJetKinFitter* fitter_jetsH = new DiJetKinFitter( "fitter_jetsH", "fitter_jets", Hmass );
  
  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%1000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);

          float puRe = evweight;

          if( ientry == 4886) continue;

          //if(pu_weight == 0 || weight == 0 || evweight == 0) continue;

          int n_good_jets_pre = 0;
          int n_b_jets_pre = 0;
          int n_good_jets = 0;
          int n_b_jets = 0;

          float btag_weight = 1.;
          float btag_weight_u = 1.;
          float btag_weight_d = 1.;
          float btag_weight_20u = 1.;
          float btag_weight_20d = 1.;
          float btag_weight_05u = 1.;
          float btag_weight_05d = 1.;
          float btag_weight_15u = 1.;
          float btag_weight_15d = 1.;


          // --------------------------------------------------------------------------------
          //preparing vectors with the infos used later on
          ecorrjet[0]     = j1_e;           
          ecorrjet[1]     = j2_e;           
          ecorrjet[2]     = j3_e;           
          ecorrjet[3]     = j4_e;
          ptcorrjet[0]    = j1_pt;          
          ptcorrjet[1]    = j2_pt;          
          ptcorrjet[2]    = j3_pt;          
          ptcorrjet[3]    = j4_pt;
          ptuncorrjet[0]  = j1_pt;          
          ptuncorrjet[1]  = j2_pt;          
          ptuncorrjet[2]  = j3_pt;          
          ptuncorrjet[3]  = j4_pt;
          etajet[0]       = j1_eta;         
          etajet[1]       = j2_eta;         
          etajet[2]       = j3_eta;        
          etajet[3]       = j4_eta;
          phijet[0]       = j1_phi;         
          phijet[1]       = j2_phi;         
          phijet[2]       = j3_phi;        
          phijet[3]       = j4_phi;
          btagjprobjet[0] = j1_jetProbBtag; 
          btagjprobjet[1] = j2_jetProbBtag; 
          btagjprobjet[2] = j3_jetProbBtag; 
          btagjprobjet[3] = j4_jetProbBtag;
          btagcsvjet[0]   = j1_csvBtag;     
          btagcsvjet[1]   = j2_csvBtag;     
          btagcsvjet[2]   = j3_csvBtag;    
          btagcsvjet[3]   = j4_csvBtag;
          cefjet[0]       = j1_emfrac;      
          cefjet[1]       = j2_emfrac;      
          cefjet[2]       = j3_emfrac;      
          cefjet[3]       = j4_emfrac;
          chfjet[0]       = j1_hadfrac;     
          chfjet[1]       = j2_hadfrac;     
          chfjet[2]       = j3_hadfrac;     
          chfjet[3]       = j4_hadfrac;
          vtxPtjet[0]     = j1_secVtxPt;    
          vtxPtjet[1]     = j2_secVtxPt;    
          vtxPtjet[2]     = j3_secVtxPt;    
          vtxPtjet[3]     = j4_secVtxPt;
          vtx3dljet[0]    = j1_secVtx3dL;   
          vtx3dljet[1]    = j2_secVtx3dL;   
          vtx3dljet[2]    = j3_secVtx3dL;   
          vtx3dljet[3]    = j4_secVtx3dL;
          nconstjet[0]    = (float)(j1_nNeutrals + j1_nCharged);   
          nconstjet[1]    = (float)(j2_nNeutrals + j2_nCharged);
          nconstjet[2]    = (float)(j3_nNeutrals + j3_nCharged);
          nconstjet[3]    = (float)(j4_nNeutrals + j4_nCharged);

          //----------------------------------------------------------------------------------
          
          // applying the jet regression
          if(useMVA){
           
           TVector3 tempT3jet, t3met;
           t3met.SetPtEtaPhi(met_corr_pfmet, 0, met_corr_phi_pfmet);               
           for (int ii=0; ii<4; ii++) {
	       
               if ( ptcorrjet[ii]<-1) continue;
	
	       fRegr_pt      = ptcorrjet[ii];
	       fRegr_eta     = etajet[ii];
	       fRegr_cef     = cefjet[ii];
	       fRegr_nconst  = nconstjet[ii];
	       fRegr_chf     = chfjet[ii];
	       fRegr_vtxPt   = vtxPtjet[ii];
	       fRegr_vtx3dl  = vtx3dljet[ii];
	       fRegr_met = met_corr_pfmet;
	       tempT3jet.SetPtEtaPhi(ptcorrjet[ii], etajet[ii], phijet[ii]);
	       fRegr_dPhiMet = tempT3jet.DeltaPhi(t3met);
	
	       float thePtCorr = 10000;
	       thePtCorr = (float)(readerRegres->EvaluateMVA("BDT"));
                  
               // corrected pT and energy 
	       float correctionFactor = thePtCorr/ptcorrjet[ii]; 
	       ptcorrjet[ii] = thePtCorr;  
	       ecorrjet[ii] *= correctionFactor;
            
           }
          }
         
          phoP4[0] = new TLorentzVector();
          phoP4[0]->SetPtEtaPhiE(ph1_pt,ph1_eta,ph1_phi,ph1_e);
          
          phoP4[1] = new TLorentzVector();
          phoP4[1]->SetPtEtaPhiE(ph2_pt,ph2_eta,ph2_phi,ph2_e);

          if((TMath::Abs(ph1_SCEta)>1.4442&&TMath::Abs(ph1_SCEta)<1.566)||(TMath::Abs(ph2_SCEta)>1.4442&&TMath::Abs(ph2_SCEta)<1.566)
       || TMath::Abs(ph1_SCEta)>2.5 || TMath::Abs(ph2_SCEta)>2.5) continue;  //
           
          if(ph1_pt > ph2_pt && (ph2_pt < 25 || ph1_pt < 40.*PhotonsMass/120.) ) continue;
	  if(ph2_pt > ph2_pt && (ph1_pt < 25 || ph2_pt < 40.*PhotonsMass/120.) ) continue;

          if(ph1_ciclevel < 4) continue;
          if(ph2_ciclevel < 4) continue;

          if(PhotonsMass < 100. || PhotonsMass > 180.) continue;

          for(int ii = 0; ii < 4; ii++){
              bad_jet[ii] = false;
              bad_jet_pre[ii] = false;
          }

          for(int ii = 0; ii < 4; ii++){
              jet_weight[ii] = 1.;
              jet_weight_max[ii] = 1.;
              jet_weight_min[ii] = 1.;
          }
              

          for (int ij=0; ij<4; ij++) { 
               if ( ptcorrjet[ij]<-1)               bad_jet_pre[ij] = true;
               if ( ptcorrjet[ij]<-1)               bad_jet[ij] = true;
               if ( btagcsvjet[ij]<=0)              bad_jet[ij] = true;
               if ( ptcorrjet[ij] < 25. )           bad_jet[ij] = true;
               if ( fabs(etajet[ij])>2.5 )          bad_jet[ij] = true;
               if ( !passCutBasedJetId(ij) )        bad_jet[ij] = true;
          }

          for(int ii = 0; ii < 4; ii++){
              if(bad_jet_pre[ii] == true) continue;
              n_good_jets_pre++;
          }

          for(int ii = 0; ii < 4; ii++){
              if(bad_jet[ii] == true) continue;
              jetP4[ii] = new TLorentzVector();
              jetP4[ii]->SetPtEtaPhiE(ptcorrjet[ii],etajet[ii],phijet[ii],ecorrjet[ii]);
              n_good_jets++;
          }
          
          jet_csvBtag[0] = j1_csvBtag;
          jet_csvBtag[1] = j2_csvBtag;
          jet_csvBtag[2] = j3_csvBtag;
          jet_csvBtag[3] = j4_csvBtag;
      
          for(int ii = 0; ii < 4; ii++)
              if(jet_csvBtag[ii] > 0.679 && bad_jet_pre[ii] == false) n_b_jets_pre++;

          for(int ii = 0; ii < 4; ii++)
              if(jet_csvBtag[ii] > 0.679 && bad_jet[ii] == false) n_b_jets++;
          
          if(n_good_jets_pre < 2) continue;
          if(n_b_jets_pre < 1) continue;
          if(n_good_jets < 2) continue; 
          
          if(n_b_jets < 1) continue;

          float invMass_Mjj = 0;
          float invMass_Mggjj = 0;
          float maxPt = 0;
          bool btagged[4] ={0,0,0,0};
          int j1,j2;
         
         for(int ii = 0; ii < 4; ii++){

              if(bad_jet[ii] == true) continue;
              
              if(jet_csvBtag[ii] > 0.679) btagged[ii] = true;
            
          }

          for(int ii = 0; ii < 4; ii++){
            
              if(bad_jet[ii] == true) continue;

              for(int jj = 0; jj < 4; jj++){

                  if(bad_jet[jj] == true) continue;

                  if(ii == jj) continue;

                  sumP4 = new TLorentzVector();

                  *sumP4 = *jetP4[ii]+*jetP4[jj];
                   
                  if(n_b_jets == 1)
                     if(btagged[ii] == true || btagged[jj] == true){
                        jjPt.push_back(sumP4->Pt());
                        jjPt_map.insert(std::pair<float,std::pair<int,int> >(sumP4->Pt(),std::pair<int,int>(ii,jj)));
                     }

                  if(n_b_jets >= 2)
                     if(btagged[ii] == true && btagged[jj] == true){
                        jjPt.push_back(sumP4->Pt());
                        jjPt_map.insert(std::pair<float,std::pair<int,int> >(sumP4->Pt(),std::pair<int,int>(ii,jj)));
                     }

                  delete sumP4;
              }
              
          }

          std::sort(jjPt.begin(),jjPt.end());

          j1 = jjPt_map[jjPt.at(jjPt.size()-1)].first;
          j2 = jjPt_map[jjPt.at(jjPt.size()-1)].second;

          jjPt.clear();
          jjPt_map.clear();

          Hgg = new TLorentzVector();
          *Hgg = *phoP4[0]+*phoP4[1];  
          Hjj = new TLorentzVector();
          *Hjj = *jetP4[j1]+*jetP4[j2];
          Radion = new TLorentzVector();
          *Radion = *Hgg + *Hjj;
          Hgg->Boost(-Radion->BoostVector());

          if(fabs(Hgg->CosTheta()) >= 0.9) continue;
         
          delete Hgg;
          delete Hjj;
          delete Radion;

          sumP4 = new TLorentzVector();

          *sumP4 = *jetP4[j1]+*jetP4[j2];
           
          invMass_Mjj = sumP4->M();

          delete sumP4;

          sumP4 = new TLorentzVector();

          *sumP4 = *jetP4[j1]+*jetP4[j2]+*phoP4[0]+*phoP4[1];
          invMass_Mggjj = sumP4->M();
           
          delete sumP4;

          if(n_b_jets == 1 && (invMass_Mjj < 90. || invMass_Mjj > 150.)) continue;
          if(n_b_jets > 1 && (invMass_Mjj < 95. || invMass_Mjj > 140.)) continue;


          TVector3 t3phot1, t3phot2;
          t3phot1.SetPtEtaPhi(phoP4[0]->Pt(),phoP4[0]->Eta(),phoP4[0]->Phi());
          t3phot2.SetPtEtaPhi(phoP4[1]->Pt(),phoP4[1]->Eta(),phoP4[1]->Phi());
          TLorentzVector t4phot1, t4phot2;
          t4phot1.SetPtEtaPhiE(phoP4[0]->Pt(),phoP4[0]->Eta(),phoP4[0]->Phi(),phoP4[0]->Energy());
          t4phot2.SetPtEtaPhiE(phoP4[1]->Pt(),phoP4[1]->Eta(),phoP4[1]->Phi(),phoP4[1]->Energy());
          
          TVector3 t3jet1, t3jet2;
          t3jet1.SetPtEtaPhi(jetP4[j1]->Pt(),jetP4[j1]->Eta(),jetP4[j1]->Phi());
          t3jet2.SetPtEtaPhi(jetP4[j2]->Pt(),jetP4[j2]->Eta(),jetP4[j2]->Phi());
          TLorentzVector t4jet1, t4jet2;
          t4jet1.SetPtEtaPhiE(jetP4[j1]->Pt(),jetP4[j1]->Eta(),jetP4[j1]->Phi(),jetP4[j1]->Energy());
          t4jet2.SetPtEtaPhiE(jetP4[j2]->Pt(),jetP4[j2]->Eta(),jetP4[j2]->Phi(),jetP4[j2]->Energy());

          std::pair<TLorentzVector,TLorentzVector> jets_kinfitH = fitter_jetsH->fit(t4jet1, t4jet2);
          TLorentzVector jet1_kinfit, jet2_kinfit;
          jet1_kinfit = jets_kinfitH.first;
          jet2_kinfit = jets_kinfitH.second;

          TLorentzVector t4diPhot = t4phot1+t4phot2;
          TLorentzVector dijet_kinfit = jet1_kinfit + jet2_kinfit;
          TLorentzVector Vstar_kinfit = dijet_kinfit + t4diPhot;
           
          float radMassKF = Vstar_kinfit.M();
          if(n_b_jets == 1 && (radMassKF<260. || radMassKF>335.)) continue;
          if(n_b_jets > 1 && (radMassKF<255. || radMassKF>320.)) continue;

          float deltaR_gg    = t4phot1.DeltaR(t4phot2);
          float deltaR_g1b1  = t4phot1.DeltaR(jet1_kinfit);
          float deltaR_g1b2  = t4phot1.DeltaR(jet2_kinfit);
          float deltaR_g2b1  = t4phot2.DeltaR(jet1_kinfit);
          float deltaR_g2b2  = t4phot2.DeltaR(jet2_kinfit);
          float minDeltaR_gb = deltaR_g1b1;
          if ( deltaR_g1b2 < minDeltaR_gb ) minDeltaR_gb = deltaR_g1b2;
          if ( deltaR_g2b1 < minDeltaR_gb ) minDeltaR_gb = deltaR_g2b1;
          if ( deltaR_g2b2 < minDeltaR_gb ) minDeltaR_gb = deltaR_g2b2;

          if (minDeltaR_gb<1) continue;

          if (njets_passing_kLooseID>=3) continue;
          
          n_good_jets_total = n_good_jets_total + n_good_jets;
          n_event++;  
        
            
          //if(n_b_jets <= 1) continue;
          //if(n_b_jets != 1) continue;
    
          for(int jj = 0; jj < 4; jj++){

             if(bad_jet[jj] == true) continue;
             
             for(unsigned int ii = 0; ii < PREjet_p4[lumis][event].size(); ii++){

                  if(PREbad_jet[lumis][event].at(ii) == true) continue;
                  
                  //if(PREjet_p4[lumis][event].at(ii).Pt()/ptuncorrjet[jj] < 0.4) continue;

                  dR.push_back(jetP4[jj]->DeltaR(PREjet_p4[lumis][event].at(ii)));
                  dR_map[jetP4[jj]->DeltaR(PREjet_p4[lumis][event].at(ii))] = ii;
              
              }
              
              std::sort(dR.begin(),dR.end());

              if(dR.size() == 0) continue;
            
              if(dR.at(0) < 0.3) jet_weight[jj] = PREjet_weight[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_max[jj] = PREjet_weight_max[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_min[jj] = PREjet_weight_min[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_05max[jj] = PREjet_weight_05max[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_05min[jj] = PREjet_weight_05min[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_15max[jj] = PREjet_weight_15max[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_15min[jj] = PREjet_weight_15min[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_20max[jj] = PREjet_weight_20max[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_weight_20min[jj] = PREjet_weight_20min[lumis][event].at(dR_map[dR.at(0)]);
              if(dR.at(0) < 0.3) jet_flavour[jj] = PREjet_flavour[lumis][event].at(dR_map[dR.at(0)]);
              
              dR.clear();
              dR_map.clear();
              
          }
        
          float rho = 0;
          if((jet_flavour[j1] == 5 && jet_flavour[j2] == 4) || (jet_flavour[j1] == 4 && jet_flavour[j2] == 5) || (jet_flavour[j1] == 5 && jet_flavour[j2] == 5) || (jet_flavour[j1] == 4 && jet_flavour[j2] == 4) || (jet_flavour[j1] != 5 && jet_flavour[j2] != 4) || (jet_flavour[j1] != 4 && jet_flavour[j2] != 5) || (jet_flavour[j1] != 5 && jet_flavour[j2] != 5) || (jet_flavour[j1] != 4 && jet_flavour[j2] != 4)) rho = 1.;
            
          btag_weight = jet_weight[j1]*jet_weight[j2];
          float err1_up = fabs(jet_weight_max[j1]-jet_weight[j1]);
          float err1_down = fabs(jet_weight_min[j1]-jet_weight[j1]);
          float err2_up = fabs(jet_weight_max[j2]-jet_weight[j2]);
          float err2_down = fabs(jet_weight_min[j2]-jet_weight[j2]);
          float weight_error_up = weight_err(jet_weight[j1],err1_up,jet_weight[j2],err2_up,btag_weight,rho);
          float weight_error_down = weight_err(jet_weight[j1],err1_down,jet_weight[j2],err2_down,btag_weight,rho);
             
          btag_weight_u = btag_weight+weight_error_up;
          btag_weight_d = btag_weight-weight_error_down;

          err1_up = fabs(jet_weight_05max[j1]-jet_weight[j1]);
          err1_down = fabs(jet_weight_05min[j1]-jet_weight[j1]);
          err2_up = fabs(jet_weight_05max[j2]-jet_weight[j2]);
          err2_down = fabs(jet_weight_05min[j2]-jet_weight[j2]);
          weight_error_up = weight_err(jet_weight[j1],err1_up,jet_weight[j2],err2_up,btag_weight,rho);
          weight_error_down = weight_err(jet_weight[j1],err1_down,jet_weight[j2],err2_down,btag_weight,rho);
  
          btag_weight_05u = btag_weight+weight_error_up;
          btag_weight_05d = btag_weight-weight_error_down;

          //if( ientry == 4886) std::cout << ientry << " - " << err1_up << " - " << err1_down << " - " << err2_up << " - " << err2_down << " - " << weight_error_up << " - " << weight_error_down << " - " << btag_weight << " - " << btag_weight_05u << " - " << btag_weight_05d << std::endl;
          //if( ientry == 4886) std::cout << ientry << " - " << jet_weight[j1] << " - " << jet_weight_max[j1] << " - " << jet_weight_min[j1] << " - " << jet_weight_05max[j1] << " - " << jet_weight_05min[j1] << " - " << jet_weight[j2] << " - " << jet_weight_max[j2] << " - " << jet_weight_min[j2] << " - " << jet_weight_05max[j2] << " - " << jet_weight_05min[j2] << " - " << weight_error_up << " - " << weight_error_down << " - " << btag_weight << " - " << btag_weight_05u << " - " << btag_weight_05d << std::endl;

          err1_up = fabs(jet_weight_15max[j1]-jet_weight[j1]);
          err1_down = fabs(jet_weight_15min[j1]-jet_weight[j1]);
          err2_up = fabs(jet_weight_15max[j2]-jet_weight[j2]);
          err2_down = fabs(jet_weight_15min[j2]-jet_weight[j2]);
          weight_error_up = weight_err(jet_weight[j1],err1_up,jet_weight[j2],err2_up,btag_weight,rho);
          weight_error_down = weight_err(jet_weight[j1],err1_down,jet_weight[j2],err2_down,btag_weight,rho);
  
          btag_weight_15u = btag_weight+weight_error_up;
          btag_weight_15d = btag_weight-weight_error_down;

          err1_up = fabs(jet_weight_20max[j1]-jet_weight[j1]);
          err1_down = fabs(jet_weight_20min[j1]-jet_weight[j1]);
          err2_up = fabs(jet_weight_20max[j2]-jet_weight[j2]);
          err2_down = fabs(jet_weight_20min[j2]-jet_weight[j2]);
          weight_error_up = weight_err(jet_weight[j1],err1_up,jet_weight[j2],err2_up,btag_weight,rho);
          weight_error_down = weight_err(jet_weight[j1],err1_down,jet_weight[j2],err2_down,btag_weight,rho);
  
          btag_weight_20u = btag_weight+weight_error_up;
          btag_weight_20d = btag_weight-weight_error_down;
          
          if(n_b_jets == 1 && btagged[j1] == true && fabs(jet_flavour[j1]) != 5) n_1btag_out++;
          if(n_b_jets == 1 && btagged[j2] == true && fabs(jet_flavour[j2]) != 5) n_1btag_out++;
          if(n_b_jets == 1 && btagged[j1] == true && fabs(jet_flavour[j1]) == 4) n_1btag_out_c++;
          if(n_b_jets == 1 && btagged[j2] == true && fabs(jet_flavour[j2]) == 4) n_1btag_out_c++;
          if(n_b_jets == 1 && btagged[j1] == false && fabs(jet_flavour[j1]) == 4) n_1btag_in_c++;
          if(n_b_jets == 1 && btagged[j2] == false && fabs(jet_flavour[j2]) == 4) n_1btag_in_c++;
          if(n_b_jets == 1 && btagged[j1] != true && fabs(jet_flavour[j1]) == 5) n_1btag_in++;
          if(n_b_jets == 1 && btagged[j2] != true && fabs(jet_flavour[j2]) == 5) n_1btag_in++;
          if(n_b_jets >= 2 && fabs(jet_flavour[j1]) != 5) n_2btag_out++;
          if(n_b_jets >= 2 && fabs(jet_flavour[j2]) != 5) n_2btag_out++;
          if(n_b_jets >= 2 && fabs(jet_flavour[j1]) == 5) n_2btag_in++;
          if(n_b_jets >= 2 && fabs(jet_flavour[j2]) == 5) n_2btag_in++;
             
          h_btagSF_weights->Fill(btag_weight);
          h_btagSF_weights_up->Fill(btag_weight_u);
          h_btagSF_weights_down->Fill(btag_weight_d);          

          if(n_b_jets == 1){
             h_btagSF_weights_1b->Fill(btag_weight);
             h_btagSF_weights_up_1b->Fill(btag_weight_u);
             h_btagSF_weights_down_1b->Fill(btag_weight_d);    
          }
          if(n_b_jets >= 2){
             h_btagSF_weights_2b->Fill(btag_weight);
             h_btagSF_weights_up_2b->Fill(btag_weight_u);
             h_btagSF_weights_down_2b->Fill(btag_weight_d);    
          }

          h_diPhodiJet_InvMass_old->Fill(invMass_Mggjj,puRe);
          h_diPhodiJet_InvMass->Fill(invMass_Mggjj,puRe*btag_weight);
          h_diPhodiJet_InvMass_up->Fill(invMass_Mggjj,puRe*btag_weight_u);
          h_diPhodiJet_InvMass_down->Fill(invMass_Mggjj,puRe*btag_weight_d);
          
          if(n_b_jets == 1){
             h_diPhodiJet_InvMass_old_1btag->Fill(invMass_Mggjj,puRe);
             h_diPhodiJet_InvMass_1btag->Fill(invMass_Mggjj,puRe*btag_weight);
             h_diPhodiJet_InvMass_up_1btag->Fill(invMass_Mggjj,puRe*btag_weight_u);
             h_diPhodiJet_InvMass_down_1btag->Fill(invMass_Mggjj,puRe*btag_weight_d);
             //std::cout << ientry << " - " << puRe*btag_weight_05u << " - " << puRe*btag_weight_05d << std::endl;
             h_diPhodiJet_InvMass_05up_1btag->Fill(invMass_Mggjj,puRe*btag_weight_05u);
             h_diPhodiJet_InvMass_05down_1btag->Fill(invMass_Mggjj,puRe*btag_weight_05d);
             h_diPhodiJet_InvMass_15up_1btag->Fill(invMass_Mggjj,puRe*btag_weight_15u);
             h_diPhodiJet_InvMass_15down_1btag->Fill(invMass_Mggjj,puRe*btag_weight_15d);
             h_diPhodiJet_InvMass_20up_1btag->Fill(invMass_Mggjj,puRe*btag_weight_20u);
             h_diPhodiJet_InvMass_20down_1btag->Fill(invMass_Mggjj,puRe*btag_weight_20d);
          }

          if(n_b_jets >= 2){
             h_diPhodiJet_InvMass_old_2btag->Fill(invMass_Mggjj,puRe);
             h_diPhodiJet_InvMass_2btag->Fill(invMass_Mggjj,puRe*btag_weight);
             h_diPhodiJet_InvMass_up_2btag->Fill(invMass_Mggjj,puRe*btag_weight_u);
             h_diPhodiJet_InvMass_down_2btag->Fill(invMass_Mggjj,puRe*btag_weight_d);
             h_diPhodiJet_InvMass_05up_2btag->Fill(invMass_Mggjj,puRe*btag_weight_05u);
             h_diPhodiJet_InvMass_05down_2btag->Fill(invMass_Mggjj,puRe*btag_weight_05d);
             h_diPhodiJet_InvMass_15up_2btag->Fill(invMass_Mggjj,puRe*btag_weight_15u);
             h_diPhodiJet_InvMass_15down_2btag->Fill(invMass_Mggjj,puRe*btag_weight_15d);
             h_diPhodiJet_InvMass_20up_2btag->Fill(invMass_Mggjj,puRe*btag_weight_20u);
             h_diPhodiJet_InvMass_20down_2btag->Fill(invMass_Mggjj,puRe*btag_weight_20d);
          }

          h2_pt1pt2->Fill(jetP4[j1]->Pt(),jetP4[j2]->Pt());
          h2_eta1eta2->Fill(jetP4[j1]->Eta(),jetP4[j2]->Eta());
          h2_flav1flav2->Fill(fabs(jet_flavour[j1]),fabs(jet_flavour[j2]));
          if(n_b_jets == 1) h2_jw1jw2_1btag->Fill(jet_weight[j1],jet_weight[j2]);
          if(n_b_jets >= 2) h2_jw1jw2_2btag->Fill(jet_weight[j1],jet_weight[j2]);
          if(n_b_jets == 1) h2_jw1jw2_1btag_max->Fill(jet_weight_max[j1],jet_weight_max[j2]);
          if(n_b_jets >= 2) h2_jw1jw2_2btag_max->Fill(jet_weight_max[j1],jet_weight_max[j2]);
          if(n_b_jets == 1) h2_jw1jw2_1btag_min->Fill(jet_weight_min[j1],jet_weight_min[j2]);
          if(n_b_jets >= 2) h2_jw1jw2_2btag_min->Fill(jet_weight_min[j1],jet_weight_min[j2]);

          if(btagged[j1] == true && n_b_jets == 1) h2_jwpt_1btag_b->Fill(jet_weight[j1],jetP4[j1]->Pt());
          if(btagged[j2] == true && n_b_jets == 1) h2_jwpt_1btag_b->Fill(jet_weight[j2],jetP4[j2]->Pt());
          if(btagged[j1] == false && n_b_jets == 1) h2_jwpt_1btag_nob->Fill(jet_weight[j1],jetP4[j1]->Pt());
          if(btagged[j2] == false && n_b_jets == 1) h2_jwpt_1btag_nob->Fill(jet_weight[j2],jetP4[j2]->Pt());
          if(btagged[j1] == true && n_b_jets >= 2) h2_jwpt_2btag_b1->Fill(jet_weight[j1],jetP4[j1]->Pt());
          if(btagged[j2] == true && n_b_jets >= 2) h2_jwpt_2btag_b2->Fill(jet_weight[j2],jetP4[j2]->Pt());
          if(btagged[j1] == true && n_b_jets == 1) h2_jweta_1btag_b->Fill(jet_weight[j1],jetP4[j1]->Eta());
          if(btagged[j2] == true && n_b_jets == 1) h2_jweta_1btag_b->Fill(jet_weight[j2],jetP4[j2]->Eta());
          if(btagged[j1] == false && n_b_jets == 1) h2_jweta_1btag_nob->Fill(jet_weight[j1],jetP4[j1]->Eta());
          if(btagged[j2] == false && n_b_jets == 1) h2_jweta_1btag_nob->Fill(jet_weight[j2],jetP4[j2]->Eta());
          if(btagged[j1] == true && n_b_jets >= 2) h2_jweta_2btag_b1->Fill(jet_weight[j1],jetP4[j1]->Eta());
          if(btagged[j2] == true && n_b_jets >= 2) h2_jweta_2btag_b2->Fill(jet_weight[j2],jetP4[j2]->Eta());
          if(btagged[j1] == true && n_b_jets == 1) h2_jwflav_1btag_b->Fill(jet_weight[j1],fabs(jet_flavour[j2]));
          if(btagged[j2] == true && n_b_jets == 1) h2_jwflav_1btag_b->Fill(jet_weight[j2],fabs(jet_flavour[j2]));
          if(btagged[j1] == false && n_b_jets == 1) h2_jwflav_1btag_nob->Fill(jet_weight[j1],fabs(jet_flavour[j2]));
          if(btagged[j2] == false && n_b_jets == 1) h2_jwflav_1btag_nob->Fill(jet_weight[j2],fabs(jet_flavour[j2]));
          if(btagged[j1] == true && n_b_jets == 1) h2_jwflav_1btag_b_min->Fill(jet_weight_min[j1],fabs(jet_flavour[j2]));
          if(btagged[j2] == true && n_b_jets == 1) h2_jwflav_1btag_b_min->Fill(jet_weight_min[j2],fabs(jet_flavour[j2]));
          if(btagged[j1] == false && n_b_jets == 1) h2_jwflav_1btag_nob_min->Fill(jet_weight_min[j1],fabs(jet_flavour[j2]));
          if(btagged[j2] == false && n_b_jets == 1) h2_jwflav_1btag_nob_min->Fill(jet_weight_min[j2],fabs(jet_flavour[j2]));
          if(btagged[j1] == true && n_b_jets == 1) h2_jwflav_1btag_b_max->Fill(jet_weight_max[j1],fabs(jet_flavour[j2]));
          if(btagged[j2] == true && n_b_jets == 1) h2_jwflav_1btag_b_max->Fill(jet_weight_max[j2],fabs(jet_flavour[j2]));
          if(btagged[j1] == false && n_b_jets == 1) h2_jwflav_1btag_nob_max->Fill(jet_weight_max[j1],fabs(jet_flavour[j2]));
          if(btagged[j2] == false && n_b_jets == 1) h2_jwflav_1btag_nob_max->Fill(jet_weight_max[j2],fabs(jet_flavour[j2]));
          if(btagged[j1] == true && n_b_jets >= 2) h2_jwflav_2btag_b1->Fill(jet_weight[j1],fabs(jet_flavour[j2]));
          if(btagged[j2] == true && n_b_jets >= 2) h2_jwflav_2btag_b2->Fill(jet_weight[j2],fabs(jet_flavour[j2]));
          
          h2_csvflav->Fill(jet_csvBtag[j1],fabs(jet_flavour[j1]));
          h2_csvflav->Fill(jet_csvBtag[j2],fabs(jet_flavour[j2]));
          
          for(int ii = 0; ii < 4; ii++)
              bad_jet[ii] = false;

          for(int ii = 0; ii < 4; ii++)
              bad_jet_pre[ii] = false;
           
          for(int ii = 0; ii < 4; ii++)
              jet_csvBtag[ii] = -1;

          for(int ii = 0; ii < 4; ii++)
              jet_flavour[ii] = 0.;
          
      } 
  }

  std::cout << "SELECTION: " << n_event << " - " << n_good_jets_total  << " - " << n_flav_jets << " - " << n_out << std::endl;
  //std::cout << "SELECTION: " << n_1btag_out << " - " << n_1btag_in << " - " << n_1btag_out_c << " - " << n_1btag_in_c <<  " - " << n_2btag_out << " - " << n_2btag_in << std::endl;
  std::cout << "Integral UP 1b: " << h_diPhodiJet_InvMass_up_1btag->Integral() << std::endl;
  std::cout << "Integral Medium 1b: " << h_diPhodiJet_InvMass_1btag->Integral() << std::endl;
  std::cout << "Integral DOWN 1b: " << h_diPhodiJet_InvMass_down_1btag->Integral() << std::endl;
  std::cout << "Integral UP 2b: " << h_diPhodiJet_InvMass_up_2btag->Integral() << std::endl;
  std::cout << "Integral Medium 2b: " << h_diPhodiJet_InvMass_2btag->Integral() << std::endl;
  std::cout << "Integral DOWN 2b: " << h_diPhodiJet_InvMass_down_2btag->Integral() << std::endl;
  std::cout << "Integral old 1b: " << h_diPhodiJet_InvMass_old_1btag->Integral() << std::endl;
  std::cout << "Integral old 2b: " << h_diPhodiJet_InvMass_old_2btag->Integral() << std::endl;

  std::cout << "Integrals UP-1btag: " << h_diPhodiJet_InvMass_05up_1btag->Integral() << " - " << h_diPhodiJet_InvMass_up_1btag->Integral() << " - " << h_diPhodiJet_InvMass_15up_1btag->Integral() << " - " << h_diPhodiJet_InvMass_20up_1btag->Integral() << std::endl;
  std::cout << "Integrals DOWN-1btag: " << h_diPhodiJet_InvMass_05down_1btag->Integral() << " - " << h_diPhodiJet_InvMass_down_1btag->Integral() << " - " << h_diPhodiJet_InvMass_15down_1btag->Integral() << " - " << h_diPhodiJet_InvMass_20down_1btag->Integral() << std::endl;
  std::cout << "Integrals UP-2btag: " << h_diPhodiJet_InvMass_05up_2btag->Integral() << " - " << h_diPhodiJet_InvMass_up_2btag->Integral() << " - " << h_diPhodiJet_InvMass_15up_2btag->Integral() << " - " << h_diPhodiJet_InvMass_20up_2btag->Integral() << std::endl;
  std::cout << "Integrals DOWN-2btag: " << h_diPhodiJet_InvMass_05down_2btag->Integral() << " - " << h_diPhodiJet_InvMass_down_2btag->Integral() << " - " << h_diPhodiJet_InvMass_15down_2btag->Integral() << " - " << h_diPhodiJet_InvMass_20down_2btag->Integral() << std::endl;
  
  g_1btag2btag_correletion_up->SetPoint(0,h_diPhodiJet_InvMass_05up_1btag->Integral(),h_diPhodiJet_InvMass_05up_2btag->Integral());
  g_1btag2btag_correletion_up->SetPoint(1,h_diPhodiJet_InvMass_up_1btag->Integral(),h_diPhodiJet_InvMass_up_2btag->Integral());
  g_1btag2btag_correletion_up->SetPoint(2,h_diPhodiJet_InvMass_15up_1btag->Integral(),h_diPhodiJet_InvMass_15up_2btag->Integral());
  g_1btag2btag_correletion_up->SetPoint(3,h_diPhodiJet_InvMass_20up_1btag->Integral(),h_diPhodiJet_InvMass_20up_2btag->Integral());

  g_1btag2btag_correletion_down->SetPoint(0,h_diPhodiJet_InvMass_05down_1btag->Integral(),h_diPhodiJet_InvMass_05down_2btag->Integral());
  g_1btag2btag_correletion_down->SetPoint(1,h_diPhodiJet_InvMass_down_1btag->Integral(),h_diPhodiJet_InvMass_down_2btag->Integral());
  g_1btag2btag_correletion_down->SetPoint(2,h_diPhodiJet_InvMass_15down_1btag->Integral(),h_diPhodiJet_InvMass_15down_2btag->Integral());
  g_1btag2btag_correletion_down->SetPoint(3,h_diPhodiJet_InvMass_20down_1btag->Integral(),h_diPhodiJet_InvMass_20down_2btag->Integral());

  
  g_1btag2btag_correletion_up->Fit("pol1");
  g_1btag2btag_correletion_down->Fit("pol1");

  TFile* output = new TFile(outputName.c_str(),"RECREATE");
  output->cd();
  h_btagSF_weights->Write();
  h_btagSF_weights_up->Write();
  h_btagSF_weights_down->Write();  
  h_btagSF_weights_1b->Write();
  h_btagSF_weights_up_1b->Write();
  h_btagSF_weights_down_1b->Write();
  h_btagSF_weights_2b->Write();
  h_btagSF_weights_up_2b->Write();
  h_btagSF_weights_down_2b->Write();
  h_diPhodiJet_InvMass_old->Write();
  h_diPhodiJet_InvMass->Write();
  h_diPhodiJet_InvMass_up->Write();
  h_diPhodiJet_InvMass_down->Write();
  h_diPhodiJet_InvMass_old_1btag->Write();
  h_diPhodiJet_InvMass_1btag->Write();
  h_diPhodiJet_InvMass_up_1btag->Write();
  h_diPhodiJet_InvMass_down_1btag->Write();
  h_diPhodiJet_InvMass_old_2btag->Write();
  h_diPhodiJet_InvMass_2btag->Write();
  h_diPhodiJet_InvMass_up_2btag->Write();
  h_diPhodiJet_InvMass_down_2btag->Write();
  h2_pt1pt2->Write();
  h2_eta1eta2->Write();
  h2_flav1flav2->Write();
  h2_jw1jw2_1btag->Write();
  h2_jw1jw2_2btag->Write();
  h2_jw1jw2_1btag_min->Write();
  h2_jw1jw2_1btag_max->Write();
  h2_jw1jw2_2btag_min->Write();
  h2_jw1jw2_2btag_max->Write();
  h2_jwpt_1btag_b->Write();
  h2_jwpt_1btag_nob->Write();
  h2_jwpt_2btag_b1->Write();
  h2_jwpt_2btag_b2->Write();
  h2_jweta_1btag_b->Write();
  h2_jweta_1btag_nob->Write();
  h2_jweta_2btag_b1->Write();
  h2_jweta_2btag_b2->Write();
  h2_jwflav_1btag_b->Write();
  h2_jwflav_1btag_nob->Write();
  h2_jwflav_1btag_b_min->Write();
  h2_jwflav_1btag_nob_min->Write();
  h2_jwflav_1btag_b_max->Write();
  h2_jwflav_1btag_nob_max->Write();
  h2_jwflav_2btag_b1->Write();
  h2_jwflav_2btag_b2->Write();
  h2_csvflav->Write();
  g_1btag2btag_correletion_up->Write("g_1btag2btag_correletion_up");
  g_1btag2btag_correletion_down->Write("g_1btag2btag_correletion_down");
  output->Close();
}

bool passCutBasedJetId(int jet) {

  bool isGood = true;

  float etajet[4],thebetastarjet[4], thermsjet[4];  
  etajet[0]         = j1_eta;  
  etajet[1]         = j2_eta; 
  etajet[2]         = j3_eta;   
  etajet[3]         = j4_eta;  
  thebetastarjet[0] = j1_betaStarClassic;
  thebetastarjet[1] = j2_betaStarClassic;
  thebetastarjet[2] = j3_betaStarClassic;
  thebetastarjet[3] = j4_betaStarClassic;
  thermsjet[0]      = j1_dR2Mean;
  thermsjet[1]      = j2_dR2Mean;
  thermsjet[2]      = j3_dR2Mean;
  thermsjet[3]      = j4_dR2Mean;

  if ( fabs(etajet[jet]) < 2.5 ) {
    if ( thebetastarjet[jet] > 0.2 * log( nvtx - 0.64) )  isGood = false;
    if (thermsjet[jet] > 0.06)                            isGood = false;
  } 
  else if (fabs(etajet[jet]) < 3.){
    if ( thermsjet[jet] > 0.05)  isGood =false;
  } 
  else {
    if ( thermsjet[jet] > 0.055) isGood =false;
  }

  return isGood;
}

float weight_err(float jet1, float jet1_err,float jet2, float jet2_err,float weight, float rho){
      
      float rel_err1 = jet1_err/jet1;
      float rel_err2 = jet2_err/jet2;
      return weight*sqrt(rel_err1*rel_err1+rel_err2*rel_err2+2*rel_err1*rel_err2*rho);
}
















