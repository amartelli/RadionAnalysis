
#include "PUreweightingUtils.h"
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "setTDRStyle.h"
#include "drawPlotsUtils.h"

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
float jet_weight_err(float eff, float eff_err, float sf, float sf_err);


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
  std::string inputFlav = gConfigParser -> readStringOption("Input::inputFlav");
  std::string inputTree = gConfigParser -> readStringOption("Input::inputTree");
  std::string inputSF = gConfigParser -> readStringOption("Input::inputSF");
  std::string inputEff = gConfigParser -> readStringOption("Input::inputEff");
  std::string inputEff_tot = gConfigParser -> readStringOption("Input::inputEffTot");
  bool useMVA = gConfigParser -> readBoolOption("Input::useMVA");
 
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");
  
  std::map<int,TChain*> ntu;
  std::map<int, TFile*> Files;

  char trees[500];
  FILE *f_trees;
  FILE *f_flav;
  
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

  std::map<int,std::map<int,std::vector<TLorentzVector> > > AODjet_p4;
  std::map<int,std::map<int,std::vector<int> > > AODjet_Flav;

  f_flav = fopen(inputFlav.c_str(),"r");
 
  int eventId,lumiId,flav,njets,njets_f;
  float pt,eta,phi,energy;
  
  TLorentzVector p4;

  std::cout << "Reading jet flavour" << std::endl;
  
  while(fscanf(f_flav,"%d %d %f %f %f %f %d %d %d \n", &lumiId, &eventId, &pt, &eta, &phi, &energy, &flav, &njets, &njets_f) !=EOF ){
            p4.SetPtEtaPhiE(pt,eta,phi,energy);

            AODjet_p4[lumiId][eventId].push_back(p4);  
            AODjet_Flav[lumiId][eventId].push_back(flav);          
      }
 
  
  for(int ii = 0; ii < pos_total; ii++){
      if(ntu[ii]->GetEntries() == 0 )
      {
         std::cout << "Error: input file" << ii << " is empty" << std::endl; 
         return -1;
      }
  }

  TFile* btagSF_File = new TFile(inputSF.c_str(),"READ");

  TF1* SFb_CSVM = (TF1*)btagSF_File->Get("SFb_CSVM");
  TH1F* h1_SFb_CSVM = (TH1F*)btagSF_File->Get("h1_SFb_CSVM");

  TF1* SFudsg_CSVM_00_08_max = (TF1*)btagSF_File->Get("SFudsg_CSVM_00_08_max");
  TF1* SFudsg_CSVM_00_08_mean = (TF1*)btagSF_File->Get("SFudsg_CSVM_00_08_mean");
  TF1* SFudsg_CSVM_00_08_min = (TF1*)btagSF_File->Get("SFudsg_CSVM_00_08_min");
   
  TF1* SFudsg_CSVM_08_16_max = (TF1*)btagSF_File->Get("SFudsg_CSVM_08_16_max");
  TF1* SFudsg_CSVM_08_16_mean = (TF1*)btagSF_File->Get("SFudsg_CSVM_08_16_mean");
  TF1* SFudsg_CSVM_08_16_min = (TF1*)btagSF_File->Get("SFudsg_CSVM_08_16_min");
   
  TF1* SFudsg_CSVM_16_24_max = (TF1*)btagSF_File->Get("SFudsg_CSVM_16_24_max");
  TF1* SFudsg_CSVM_16_24_mean = (TF1*)btagSF_File->Get("SFudsg_CSVM_16_24_mean");
  TF1* SFudsg_CSVM_16_24_min = (TF1*)btagSF_File->Get("SFudsg_CSVM_16_24_min");

  TFile* btagEff_File = new TFile(inputEff.c_str(),"READ");
  TFile* btagEff_File_tot = new TFile(inputEff_tot.c_str(),"READ");
  
  TH2F* h2_BTaggingEff_Denom_b_M = (TH2F*)btagEff_File->Get("h2_BTaggingEff_Denom_b_M");
  TH2F* h2_BTaggingEff_Num_b_M = (TH2F*)btagEff_File->Get("h2_BTaggingEff_Num_b_M");
  TH2F* h2_BTaggingEff_Denom_c_M = (TH2F*)btagEff_File->Get("h2_BTaggingEff_Denom_c_M");
  TH2F* h2_BTaggingEff_Num_c_M = (TH2F*)btagEff_File->Get("h2_BTaggingEff_Num_c_M");
  TH2F* h2_BTaggingEff_Denom_udsg_M = (TH2F*)btagEff_File->Get("h2_BTaggingEff_Denom_udsg_M");
  TH2F* h2_BTaggingEff_Num_udsg_M = (TH2F*)btagEff_File->Get("h2_BTaggingEff_Num_udsg_M");

  TH2F* h2_BTaggingEff_Denom_b_M_tot = (TH2F*)btagEff_File_tot->Get("h2_BTaggingEff_Denom_b_M");
  TH2F* h2_BTaggingEff_Num_b_M_tot = (TH2F*)btagEff_File_tot->Get("h2_BTaggingEff_Num_b_M");
  TH2F* h2_BTaggingEff_Denom_c_M_tot = (TH2F*)btagEff_File_tot->Get("h2_BTaggingEff_Denom_c_M");
  TH2F* h2_BTaggingEff_Num_c_M_tot = (TH2F*)btagEff_File_tot->Get("h2_BTaggingEff_Num_c_M");
  TH2F* h2_BTaggingEff_Denom_udsg_M_tot = (TH2F*)btagEff_File_tot->Get("h2_BTaggingEff_Denom_udsg_M");
  TH2F* h2_BTaggingEff_Num_udsg_M_tot = (TH2F*)btagEff_File_tot->Get("h2_BTaggingEff_Num_udsg_M");

  float ptmin_b[17] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  float etamin_b[4] = {0.,0.8,1.6,2.5};
  //float ptmin_b[14] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320,350.};
  //float etamin_b[6] = {0.,0.5,1.,1.5,2.,2.5};

  TH2F* h2_BTaggingEff_b_M = new TH2F("h2_BTaggingEff_b_M","h2_BTaggingEff_b_M",16, ptmin_b, 3,etamin_b);
  TH2F* h2_BTaggingEff_b_M_tot = new TH2F("h2_BTaggingEff_b_M_tot","h2_BTaggingEff_b_M_tot",16, ptmin_b, 3,etamin_b);
  //TH2F* h2_BTaggingEff_b_M = new TH2F("h2_BTaggingEff_b_M","h2_BTaggingEff_b_M",13, ptmin_b, 5,etamin_b);
  
  float ptmin_c[17] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  float etamin_c[4] = {0.,0.8,1.6,2.5};
  //float ptmin_c[14] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320,350.};
  //float etamin_c[6] = {0.,0.5,1.,1.5,2.,2.5};

  TH2F* h2_BTaggingEff_c_M = new TH2F("h2_BTaggingEff_c_M","h2_BTaggingEff_c_M",16, ptmin_c, 3,etamin_c);
  TH2F* h2_BTaggingEff_c_M_tot = new TH2F("h2_BTaggingEff_c_M_tot","h2_BTaggingEff_c_M_tot",16, ptmin_c, 3,etamin_c);
  //TH2F* h2_BTaggingEff_c_M = new TH2F("h2_BTaggingEff_c_M","h2_BTaggingEff_c_M",13, ptmin_c, 5,etamin_c);
  
  float ptmin_udsg_M[17] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  float etamin_udsg_M[4] = {0.,0.8,1.6,2.5};
  //float ptmin_udsg_M[14] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320,350.};
  //float etamin_udsg_M[6] = {0.,0.5,1.,1.5,2.,2.5};
  
  TH2F* h2_BTaggingEff_udsg_M = new TH2F("h2_BTaggingEff_udsg_M","h2_BTaggingEff_udsg_M",16, ptmin_udsg_M, 3,etamin_udsg_M);
  TH2F* h2_BTaggingEff_udsg_M_tot = new TH2F("h2_BTaggingEff_udsg_M_tot","h2_BTaggingEff_udsg_M_tot",16, ptmin_udsg_M, 3,etamin_udsg_M);
  //TH2F* h2_BTaggingEff_udsg_M = new TH2F("h2_BTaggingEff_udsg_M","h2_BTaggingEff_udsg_M",13, ptmin_udsg_M, 5,etamin_udsg_M);
  
  h2_BTaggingEff_b_M->Divide(h2_BTaggingEff_Num_b_M,h2_BTaggingEff_Denom_b_M,1,1,"B");
  h2_BTaggingEff_c_M->Divide(h2_BTaggingEff_Num_c_M,h2_BTaggingEff_Denom_c_M,1,1,"B");
  h2_BTaggingEff_udsg_M->Divide(h2_BTaggingEff_Num_udsg_M,h2_BTaggingEff_Denom_udsg_M,1,1,"B");
  h2_BTaggingEff_b_M_tot->Divide(h2_BTaggingEff_Num_b_M_tot,h2_BTaggingEff_Denom_b_M_tot,1,1,"B");
  h2_BTaggingEff_c_M_tot->Divide(h2_BTaggingEff_Num_c_M_tot,h2_BTaggingEff_Denom_c_M_tot,1,1,"B");
  h2_BTaggingEff_udsg_M_tot->Divide(h2_BTaggingEff_Num_udsg_M_tot,h2_BTaggingEff_Denom_udsg_M_tot,1,1,"B");
  
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
  float         event_outtree;
  float         lumis_outtree;
  float         evweight;
  float         weight;
  float         pu_weight;
  float         met_corr_pfmet;
  float         met_corr_phi_pfmet;
  float         met_corr_eta_pfmet;
  float         met_corr_e_pfmet;
  float         ph1_SCEta;
  float         ph1_pt;
  int           ph1_ciclevel;
  float         ph2_SCEta;
  float         ph2_pt;
  int           ph2_ciclevel;
  float         PhotonsMass;
  float         j1_e;
  float         j1_e_outtree;
  float         j1_pt;
  float         j1_pt_outtree;
  float         j1_phi;
  float         j1_phi_outtree;
  float         j1_eta_outtree;
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
  float         j2_e_outtree;
  float         j2_pt;
  float         j2_pt_outtree;
  float         j2_phi;
  float         j2_phi_outtree;
  float         j2_eta_outtree;
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
  float         j3_e_outtree;
  float         j3_pt;
  float         j3_pt_outtree;
  float         j3_phi;
  float         j3_phi_outtree;
  float         j3_eta_outtree;
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
  float         j4_e_outtree;
  float         j4_pt;
  float         j4_pt_outtree;
  float         j4_phi;
  float         j4_phi_outtree;
  float         j4_eta_outtree;
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

  TTree* outTree = new TTree(inputTree.c_str(),"Radion_m300_8TeV_nm");
  outTree->SetDirectory(0);

   outTree -> Branch("event",    &event_outtree,        "event/F");
   outTree -> Branch("lumis",    &lumis_outtree,        "lumis/F");
   outTree -> Branch("j1_pt",    &j1_pt_outtree,        "j1_pt/F");
   outTree -> Branch("j1_e",    &j1_e_outtree,          "j1_e/F");
   outTree -> Branch("j1_eta",    &j1_eta_outtree,      "j1_eta/F");
   outTree -> Branch("j1_phi",    &j1_phi_outtree,      "j1_phi/F");
   outTree -> Branch("j1_flavour",    &j1_flavour,      "j1_flavour/F");
   outTree -> Branch("j1_weight",    &j1_weight,        "j1_weight/F");
   outTree -> Branch("j1_weight_min",  &j1_weight_min,  "j1_weight_min/F");
   outTree -> Branch("j1_weight_max",  &j1_weight_max,  "j1_weight_max/F");
   outTree -> Branch("j1_weight_05min",  &j1_weight_05min,  "j1_weight_05min/F");
   outTree -> Branch("j1_weight_05max",  &j1_weight_05max,  "j1_weight_05max/F");
   outTree -> Branch("j1_weight_15min",  &j1_weight_15min,  "j1_weight_15min/F");
   outTree -> Branch("j1_weight_15max",  &j1_weight_15max,  "j1_weight_15max/F");
   outTree -> Branch("j1_weight_20min",  &j1_weight_20min,  "j1_weight_20min/F");
   outTree -> Branch("j1_weight_20max",  &j1_weight_20max,  "j1_weight_20max/F");
   outTree -> Branch("j2_pt",    &j2_pt_outtree,        "j2_pt/F");
   outTree -> Branch("j2_e",    &j2_e_outtree,          "j2_e/F");
   outTree -> Branch("j2_eta",    &j2_eta_outtree,      "j2_eta/F");
   outTree -> Branch("j2_phi",    &j2_phi_outtree,      "j2_phi/F");
   outTree -> Branch("j2_flavour",    &j2_flavour,      "j2_flavour/F");
   outTree -> Branch("j2_weight",    &j2_weight,        "j2_weight/F");
   outTree -> Branch("j2_weight_min",  &j2_weight_min,  "j2_weight_min/F");
   outTree -> Branch("j2_weight_max",  &j2_weight_max,  "j2_weight_max/F");
   outTree -> Branch("j2_weight_05min",  &j2_weight_05min,  "j2_weight_05min/F");
   outTree -> Branch("j2_weight_05max",  &j2_weight_05max,  "j2_weight_05max/F");
   outTree -> Branch("j2_weight_15min",  &j2_weight_15min,  "j2_weight_15min/F");
   outTree -> Branch("j2_weight_15max",  &j2_weight_15max,  "j2_weight_15max/F");
   outTree -> Branch("j2_weight_20min",  &j2_weight_20min,  "j2_weight_20min/F");
   outTree -> Branch("j2_weight_20max",  &j2_weight_20max,  "j2_weight_20max/F");
   outTree -> Branch("j3_pt",    &j3_pt_outtree,        "j3_pt/F");
   outTree -> Branch("j3_e",    &j3_e_outtree,          "j3_e/F");
   outTree -> Branch("j3_eta",    &j3_eta_outtree,      "j3_eta/F");
   outTree -> Branch("j3_phi",    &j3_phi_outtree,      "j3_phi/F");
   outTree -> Branch("j3_flavour",    &j3_flavour,      "j3_flavour/F");
   outTree -> Branch("j3_weight",    &j3_weight,        "j3_weight/F");
   outTree -> Branch("j3_weight_min",  &j3_weight_min,  "j3_weight_min/F");
   outTree -> Branch("j3_weight_max",  &j3_weight_max,  "j3_weight_max/F");
   outTree -> Branch("j3_weight_05min",  &j3_weight_05min,  "j3_weight_05min/F");
   outTree -> Branch("j3_weight_05max",  &j3_weight_05max,  "j3_weight_05max/F");
   outTree -> Branch("j3_weight_15min",  &j3_weight_15min,  "j3_weight_15min/F");
   outTree -> Branch("j3_weight_15max",  &j3_weight_15max,  "j3_weight_15max/F");
   outTree -> Branch("j3_weight_20min",  &j3_weight_20min,  "j3_weight_20min/F");
   outTree -> Branch("j3_weight_20max",  &j3_weight_20max,  "j3_weight_20max/F");
   outTree -> Branch("j4_pt",    &j4_pt_outtree,        "j4_pt/F");
   outTree -> Branch("j4_e",    &j4_e_outtree,          "j4_e/F");
   outTree -> Branch("j4_eta",    &j4_eta_outtree,      "j4_eta/F");
   outTree -> Branch("j4_phi",    &j4_phi_outtree,      "j4_phi/F");
   outTree -> Branch("j4_flavour",    &j4_flavour,      "j4_flavour/F");
   outTree -> Branch("j4_weight",    &j4_weight,        "j4_weight/F");
   outTree -> Branch("j4_weight_min",  &j4_weight_min,  "j4_weight_min/F");
   outTree -> Branch("j4_weight_max",  &j4_weight_max,  "j4_weight_max/F");
   outTree -> Branch("j4_weight_05min",  &j4_weight_05min,  "j4_weight_05min/F");
   outTree -> Branch("j4_weight_05max",  &j4_weight_05max,  "j4_weight_05max/F");
   outTree -> Branch("j4_weight_15min",  &j4_weight_15min,  "j4_weight_15min/F");
   outTree -> Branch("j4_weight_15max",  &j4_weight_15max,  "j4_weight_15max/F");
   outTree -> Branch("j4_weight_20min",  &j4_weight_20min,  "j4_weight_20min/F");
   outTree -> Branch("j4_weight_20max",  &j4_weight_20max,  "j4_weight_20max/F");
   
  
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
      ntu[ii]->SetBranchStatus("ph1_SCEta",1);
      ntu[ii]->SetBranchStatus("ph1_ciclevel",1);
      ntu[ii]->SetBranchStatus("ph2_pt",1);
      ntu[ii]->SetBranchStatus("ph2_SCEta",1);
      ntu[ii]->SetBranchStatus("ph2_ciclevel",1);
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
      ntu[ii]->SetBranchAddress("ph1_SCEta",&ph1_SCEta);
      ntu[ii]->SetBranchAddress("ph1_ciclevel",&ph1_ciclevel);
      ntu[ii]->SetBranchAddress("ph2_pt",&ph2_pt);
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
  }


  TH1F* h_btagSF_weights = new TH1F("h_btagSF_weights","h_btagSF_weights",200,0,2);
  TH1F* h_btagSF_weights_up = new TH1F("h_btagSF_weights_up","h_btagSF_weights_up",200,0,2);
  TH1F* h_btagSF_weights_down = new TH1F("h_btagSF_weights_down","h_btagSF_weights_down",200,0,2);
  
  float eff_error_b[h2_BTaggingEff_b_M->GetNbinsX()][h2_BTaggingEff_b_M->GetNbinsY()];
  for(int ii = 1; ii <= h2_BTaggingEff_b_M->GetNbinsX(); ii++)
      for(int jj = 1; jj <= h2_BTaggingEff_b_M->GetNbinsY(); jj++){
          float delta = fabs(h2_BTaggingEff_b_M->GetBinContent(ii,jj)-h2_BTaggingEff_b_M_tot->GetBinContent(ii,jj));
          float w1 = h2_BTaggingEff_b_M->GetBinError(ii,jj)/(h2_BTaggingEff_b_M_tot->GetBinError(ii,jj)-delta);
          float w2 = h2_BTaggingEff_b_M_tot->GetBinError(ii,jj)/(h2_BTaggingEff_b_M->GetBinError(ii,jj)-delta);
          if(h2_BTaggingEff_b_M->GetBinError(ii,jj) != 0 && h2_BTaggingEff_b_M_tot->GetBinError(ii,jj) != 0) eff_error_b[ii][jj] = (h2_BTaggingEff_b_M->GetBinError(ii,jj)*(w1*w1)+h2_BTaggingEff_b_M_tot->GetBinError(ii,jj)*(w2*w2))/(w1*w1+w2*w2);
          if(h2_BTaggingEff_b_M->GetBinError(ii,jj) == 0) eff_error_b[ii][jj] = h2_BTaggingEff_b_M_tot->GetBinError(ii,jj);
          if(h2_BTaggingEff_b_M_tot->GetBinError(ii,jj) == 0) eff_error_b[ii][jj] = h2_BTaggingEff_b_M->GetBinError(ii,jj);

          //std::cout << delta << " - " << h2_BTaggingEff_b_M_tot->GetBinError(ii,jj) << " - " << h2_BTaggingEff_b_M->GetBinError(ii,jj) << " - " << w1 << " - " << w2 << " - " << eff_error_b[ii][jj] << std::endl;
      }

  float eff_error_c[h2_BTaggingEff_c_M->GetNbinsX()][h2_BTaggingEff_c_M->GetNbinsY()];
  for(int ii = 1; ii <= h2_BTaggingEff_c_M->GetNbinsX(); ii++)
      for(int jj = 1; jj <= h2_BTaggingEff_c_M->GetNbinsY(); jj++){
          float delta = fabs(h2_BTaggingEff_c_M->GetBinContent(ii,jj)-h2_BTaggingEff_c_M_tot->GetBinContent(ii,jj));
          float w1 = h2_BTaggingEff_c_M->GetBinError(ii,jj)/(h2_BTaggingEff_c_M_tot->GetBinError(ii,jj)-delta);
          float w2 = h2_BTaggingEff_c_M_tot->GetBinError(ii,jj)/(h2_BTaggingEff_c_M->GetBinError(ii,jj)-delta);
          if(h2_BTaggingEff_c_M->GetBinError(ii,jj) != 0 && h2_BTaggingEff_c_M_tot->GetBinError(ii,jj) != 0) eff_error_c[ii][jj] = (h2_BTaggingEff_c_M->GetBinError(ii,jj)*(w1*w1)+h2_BTaggingEff_c_M_tot->GetBinError(ii,jj)*(w2*w2))/(w1*w1+w2*w2);
          if(h2_BTaggingEff_c_M->GetBinError(ii,jj) == 0) eff_error_c[ii][jj] = h2_BTaggingEff_c_M_tot->GetBinError(ii,jj);
          if(h2_BTaggingEff_c_M_tot->GetBinError(ii,jj) == 0) eff_error_c[ii][jj] = h2_BTaggingEff_c_M->GetBinError(ii,jj);
      }

  float eff_error_udsg[h2_BTaggingEff_udsg_M->GetNbinsX()][h2_BTaggingEff_udsg_M->GetNbinsY()];
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M->GetNbinsX(); ii++)
      for(int jj = 1; jj <= h2_BTaggingEff_udsg_M->GetNbinsY(); jj++){
          float delta = fabs(h2_BTaggingEff_udsg_M->GetBinContent(ii,jj)-h2_BTaggingEff_udsg_M_tot->GetBinContent(ii,jj));
          float w1 = h2_BTaggingEff_udsg_M->GetBinError(ii,jj)/(h2_BTaggingEff_udsg_M_tot->GetBinError(ii,jj)-delta);
          float w2 = h2_BTaggingEff_udsg_M_tot->GetBinError(ii,jj)/(h2_BTaggingEff_udsg_M->GetBinError(ii,jj)-delta);
          if(h2_BTaggingEff_udsg_M->GetBinError(ii,jj) != 0 && h2_BTaggingEff_udsg_M_tot->GetBinError(ii,jj) != 0) eff_error_udsg[ii][jj] = (h2_BTaggingEff_udsg_M->GetBinError(ii,jj)*(w1*w1)+h2_BTaggingEff_udsg_M_tot->GetBinError(ii,jj)*(w2*w2))/(w1*w1+w2*w2);
          if(h2_BTaggingEff_udsg_M->GetBinError(ii,jj) == 0) eff_error_udsg[ii][jj] = h2_BTaggingEff_udsg_M_tot->GetBinError(ii,jj);
          if(h2_BTaggingEff_udsg_M_tot->GetBinError(ii,jj) == 0) eff_error_udsg[ii][jj] = h2_BTaggingEff_udsg_M->GetBinError(ii,jj);
      }
  
   TH2F* h2_BTaggingEff_b_M_old = (TH2F*)h2_BTaggingEff_b_M->Clone("h2_BTaggingEff_b_M_old"); 
   /*h2_BTaggingEff_b_M = (TH2F*)h2_BTaggingEff_b_M_tot->Clone("h2_BTaggingEff_b_M");
   for(int ii = 1; ii <= h2_BTaggingEff_b_M->GetNbinsX(); ii++)
      for(int jj = 1; jj <= h2_BTaggingEff_b_M->GetNbinsY(); jj++)
          h2_BTaggingEff_b_M->SetBinError(ii,jj,eff_error_b[ii][jj]);

      
  h2_BTaggingEff_c_M = (TH2F*)h2_BTaggingEff_c_M_tot->Clone("h2_BTaggingEff_c_M");
  for(int ii = 1; ii <= h2_BTaggingEff_c_M->GetNbinsX(); ii++)
      for(int jj = 1; jj <= h2_BTaggingEff_c_M->GetNbinsY(); jj++)
          h2_BTaggingEff_c_M->SetBinError(ii,jj,eff_error_c[ii][jj]);
    
  h2_BTaggingEff_udsg_M = (TH2F*)h2_BTaggingEff_udsg_M_tot->Clone("h2_BTaggingEff_udsg_M");
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M->GetNbinsX(); ii++)
      for(int jj = 1; jj <= h2_BTaggingEff_udsg_M->GetNbinsY(); jj++)
          h2_BTaggingEff_udsg_M->SetBinError(ii,jj,eff_error_udsg[ii][jj]);*/
  
  TGraphErrors* g_BTaggingEff_b_M_vs_pt_00_08 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_b_M->GetNbinsX(); ii++){
      g_BTaggingEff_b_M_vs_pt_00_08->SetPoint(ii-1,h2_BTaggingEff_b_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_b_M->GetBinContent(ii,1));
      g_BTaggingEff_b_M_vs_pt_00_08->SetPointError(ii-1,0.,h2_BTaggingEff_b_M->GetBinError(ii,1));
  }

  TGraphErrors* g_BTaggingEff_b_M_vs_pt_00_08_old = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_b_M_old->GetNbinsX(); ii++){
      g_BTaggingEff_b_M_vs_pt_00_08_old->SetPoint(ii-1,h2_BTaggingEff_b_M_old->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_b_M_old->GetBinContent(ii,1));
      g_BTaggingEff_b_M_vs_pt_00_08_old->SetPointError(ii-1,0.,h2_BTaggingEff_b_M_old->GetBinError(ii,1));
  }
  
  TGraphErrors* g_BTaggingEff_b_M_vs_pt_08_16 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_b_M->GetNbinsX(); ii++){
      g_BTaggingEff_b_M_vs_pt_08_16->SetPoint(ii-1,h2_BTaggingEff_b_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_b_M->GetBinContent(ii,2));
      g_BTaggingEff_b_M_vs_pt_08_16->SetPointError(ii-1,0.,h2_BTaggingEff_b_M->GetBinError(ii,2));
  }
  
  TGraphErrors* g_BTaggingEff_b_M_vs_pt_16_25 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_b_M->GetNbinsX(); ii++){
      g_BTaggingEff_b_M_vs_pt_16_25->SetPoint(ii-1,h2_BTaggingEff_b_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_b_M->GetBinContent(ii,3));
      g_BTaggingEff_b_M_vs_pt_16_25->SetPointError(ii-1,0.,h2_BTaggingEff_b_M->GetBinError(ii,3));
  }

  TGraphErrors* g_BTaggingEff_c_M_vs_pt_00_08 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_c_M->GetNbinsX(); ii++){
      g_BTaggingEff_c_M_vs_pt_00_08->SetPoint(ii-1,h2_BTaggingEff_c_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_c_M->GetBinContent(ii,1));
      g_BTaggingEff_c_M_vs_pt_00_08->SetPointError(ii-1,0.,h2_BTaggingEff_c_M->GetBinError(ii,1));
  }


  TGraphErrors* g_BTaggingEff_c_M_vs_pt_08_16 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_c_M->GetNbinsX(); ii++){
      g_BTaggingEff_c_M_vs_pt_08_16->SetPoint(ii-1,h2_BTaggingEff_c_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_c_M->GetBinContent(ii,2));
      g_BTaggingEff_c_M_vs_pt_08_16->SetPointError(ii-1,0.,h2_BTaggingEff_c_M->GetBinError(ii,2));
  }

  TGraphErrors* g_BTaggingEff_c_M_vs_pt_16_25 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_c_M->GetNbinsX(); ii++){
      g_BTaggingEff_c_M_vs_pt_16_25->SetPoint(ii-1,h2_BTaggingEff_c_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_c_M->GetBinContent(ii,3));
      g_BTaggingEff_c_M_vs_pt_16_25->SetPointError(ii-1,0.,h2_BTaggingEff_c_M->GetBinError(ii,3));
  }
   
  TGraphErrors* g_BTaggingEff_udsg_M_vs_pt_00_08 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M->GetNbinsX(); ii++){
      g_BTaggingEff_udsg_M_vs_pt_00_08->SetPoint(ii-1,h2_BTaggingEff_udsg_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_udsg_M->GetBinContent(ii,1));
      g_BTaggingEff_udsg_M_vs_pt_00_08->SetPointError(ii-1,0.,h2_BTaggingEff_udsg_M->GetBinError(ii,1));
  }

  TGraphErrors* g_BTaggingEff_udsg_M_vs_pt_08_16 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M->GetNbinsX(); ii++){
      g_BTaggingEff_udsg_M_vs_pt_08_16->SetPoint(ii-1,h2_BTaggingEff_udsg_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_udsg_M->GetBinContent(ii,2));
      g_BTaggingEff_udsg_M_vs_pt_08_16->SetPointError(ii-1,0.,h2_BTaggingEff_udsg_M->GetBinError(ii,2));
  }

  TGraphErrors* g_BTaggingEff_udsg_M_vs_pt_16_25 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M->GetNbinsX(); ii++){
      g_BTaggingEff_udsg_M_vs_pt_16_25->SetPoint(ii-1,h2_BTaggingEff_udsg_M->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_udsg_M->GetBinContent(ii,3));
      g_BTaggingEff_udsg_M_vs_pt_16_25->SetPointError(ii-1,0.,h2_BTaggingEff_udsg_M->GetBinError(ii,3));
  }
  
  TGraphErrors* g_BTaggingEff_b_M_tot_vs_pt_00_08 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_b_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_b_M_tot_vs_pt_00_08->SetPoint(ii-1,h2_BTaggingEff_b_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_b_M_tot->GetBinContent(ii,1));
      g_BTaggingEff_b_M_tot_vs_pt_00_08->SetPointError(ii-1,0.,h2_BTaggingEff_b_M_tot->GetBinError(ii,1));
  }
  
  TGraphErrors* g_BTaggingEff_b_M_tot_vs_pt_08_16 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_b_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_b_M_tot_vs_pt_08_16->SetPoint(ii-1,h2_BTaggingEff_b_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_b_M_tot->GetBinContent(ii,2));
      g_BTaggingEff_b_M_tot_vs_pt_08_16->SetPointError(ii-1,0.,h2_BTaggingEff_b_M_tot->GetBinError(ii,2));
  }
  
  TGraphErrors* g_BTaggingEff_b_M_tot_vs_pt_16_25 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_b_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_b_M_tot_vs_pt_16_25->SetPoint(ii-1,h2_BTaggingEff_b_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_b_M_tot->GetBinContent(ii,3));
      g_BTaggingEff_b_M_tot_vs_pt_16_25->SetPointError(ii-1,0.,h2_BTaggingEff_b_M_tot->GetBinError(ii,3));
  }

  TGraphErrors* g_BTaggingEff_c_M_tot_vs_pt_00_08 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_c_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_c_M_tot_vs_pt_00_08->SetPoint(ii-1,h2_BTaggingEff_c_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_c_M_tot->GetBinContent(ii,1));
      g_BTaggingEff_c_M_tot_vs_pt_00_08->SetPointError(ii-1,0.,h2_BTaggingEff_c_M_tot->GetBinError(ii,1));
  }

  TGraphErrors* g_BTaggingEff_c_M_tot_vs_pt_08_16 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_c_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_c_M_tot_vs_pt_08_16->SetPoint(ii-1,h2_BTaggingEff_c_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_c_M_tot->GetBinContent(ii,2));
      g_BTaggingEff_c_M_tot_vs_pt_08_16->SetPointError(ii-1,0.,h2_BTaggingEff_c_M_tot->GetBinError(ii,2));
  }

  TGraphErrors* g_BTaggingEff_c_M_tot_vs_pt_16_25 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_c_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_c_M_tot_vs_pt_16_25->SetPoint(ii-1,h2_BTaggingEff_c_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_c_M_tot->GetBinContent(ii,3));
      g_BTaggingEff_c_M_tot_vs_pt_16_25->SetPointError(ii-1,0.,h2_BTaggingEff_c_M_tot->GetBinError(ii,3));
  }
   
  TGraphErrors* g_BTaggingEff_udsg_M_tot_vs_pt_00_08 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_udsg_M_tot_vs_pt_00_08->SetPoint(ii-1,h2_BTaggingEff_udsg_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_udsg_M_tot->GetBinContent(ii,1));
      g_BTaggingEff_udsg_M_tot_vs_pt_00_08->SetPointError(ii-1,0.,h2_BTaggingEff_udsg_M_tot->GetBinError(ii,1));
  }

  TGraphErrors* g_BTaggingEff_udsg_M_tot_vs_pt_08_16 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_udsg_M_tot_vs_pt_08_16->SetPoint(ii-1,h2_BTaggingEff_udsg_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_udsg_M_tot->GetBinContent(ii,2));
      g_BTaggingEff_udsg_M_tot_vs_pt_08_16->SetPointError(ii-1,0.,h2_BTaggingEff_udsg_M_tot->GetBinError(ii,2));
  }

  TGraphErrors* g_BTaggingEff_udsg_M_tot_vs_pt_16_25 = new TGraphErrors();
  for(int ii = 1; ii <= h2_BTaggingEff_udsg_M_tot->GetNbinsX(); ii++){
      g_BTaggingEff_udsg_M_tot_vs_pt_16_25->SetPoint(ii-1,h2_BTaggingEff_udsg_M_tot->GetXaxis()->GetBinCenter(ii),h2_BTaggingEff_udsg_M_tot->GetBinContent(ii,3));
      g_BTaggingEff_udsg_M_tot_vs_pt_16_25->SetPointError(ii-1,0.,h2_BTaggingEff_udsg_M_tot->GetBinError(ii,3));
  }
  
  std::map<int,TLorentzVector*> jetP4;
  std::map<int,float> jet_csvBtag;
  std::map<int,bool> bad_jet;
  std::map<int,int> jet_flavour;

  int n_event = 0;
  int n_event_match = 0;
  int n_good_jets_total = 0;
  int n_flav_jets = 0;
  
  std::vector<float> dR;
  std::map<float,int> dR_map;

  std::vector<float> SFw;

  int n_out = 0;

  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%1000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);

          event_outtree = -999;
          lumis_outtree = -999;
          j1_pt_outtree = -999.;
          j1_e_outtree = -999.;
          j1_eta_outtree = -999.;
          j1_phi_outtree = -999.;
          j1_flavour = -999.;
          j1_weight_min = 1.;
          j1_weight = 1.;
          j1_weight_max = 1.;
          j2_pt_outtree = -999.;
          j2_e_outtree = -999.;
          j2_eta_outtree = -999.;
          j2_phi_outtree = -999.;
          j2_flavour = -999.;
          j2_weight_min = 1.;
          j2_weight = 1.;
          j2_weight_max = 1.;
          j3_pt_outtree = -999.;
          j3_e_outtree = -999.;
          j3_eta_outtree = -999.;
          j3_phi_outtree = -999.;
          j3_flavour = -999.;
          j3_weight_min = 1.;
          j3_weight = 1.;
          j3_weight_max = 1.;
          j4_pt_outtree = -999.;
          j4_e_outtree = -999.;
          j4_eta_outtree = -999.;
          j4_phi_outtree = -999.;
          j4_flavour = -999.;
          j4_weight_min = 1.;
          j4_weight = 1.;
          j4_weight_max = 1.;
          
          float puRe = evweight/(weight/pu_weight);

          //if(pu_weight == 0 || weight == 0 || evweight == 0) continue;

          int n_good_jets = 0;
          int n_b_jets = 0;

          float btag_weight = 1.;
          float btag_weight_u = 1.;
          float btag_weight_d = 1.;


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


          /*if((TMath::Abs(ph1_SCEta)>1.4442&&TMath::Abs(ph1_SCEta)<1.566)||(TMath::Abs(ph2_SCEta)>1.4442&&TMath::Abs(ph2_SCEta)<1.566)
       || TMath::Abs(ph1_SCEta)>2.5 || TMath::Abs(ph2_SCEta)>2.5) continue;  //
           
          if(ph1_pt > ph2_pt && (ph2_pt < 25 || ph1_pt < 40.*PhotonsMass/120.) ) continue;
	  if(ph2_pt > ph2_pt && (ph1_pt < 25 || ph2_pt < 40.*PhotonsMass/120.) ) continue;

          if(ph1_ciclevel < 4) continue;
          if(ph2_ciclevel < 4) continue;

          if(PhotonsMass < 100. || PhotonsMass > 180.) continue;*/

          for(int ii = 0; ii < 4; ii++)
              bad_jet[ii] = false;

          for(int ii = 0; ii < 4; ii++)
              jet_flavour[ii] = 0;

          for (int ij=0; ij<4; ij++) { 
               if ( ptcorrjet[ij]<-1)               bad_jet[ij] = true;
               if ( btagcsvjet[ij]<=0)              bad_jet[ij] = true;
               //if ( ptcorrjet[ij] < 25. )           bad_jet[ij] = true;
               //if ( fabs(etajet[ij])>2.5 )          bad_jet[ij] = true;
               //if ( !passCutBasedJetId(ij) )        bad_jet[ij] = true;
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
              if(jet_csvBtag[ii] > 0.679 && bad_jet[ii] == false) n_b_jets++;
                 
          //if(n_good_jets < 2) continue; 
          //if(n_b_jets < 1) continue; 

          event_outtree = event;
          lumis_outtree = lumis;
           
          n_good_jets_total = n_good_jets_total + n_good_jets;
          n_event++;
          
               
         for(int jj = 0; jj < 4; jj++){

             if(bad_jet[jj] == true) continue;
             
             for(unsigned int ii = 0; ii < AODjet_p4[lumis][event].size(); ii++){
                  
                  if(AODjet_p4[lumis][event].at(ii).Pt()/ptuncorrjet[jj] < 0.4) continue;

                  dR.push_back(jetP4[jj]->DeltaR(AODjet_p4[lumis][event].at(ii)));
                  dR_map[jetP4[jj]->DeltaR(AODjet_p4[lumis][event].at(ii))] = ii;
              
              }
              
              std::sort(dR.begin(),dR.end());
            
              if(dR.at(0) < 0.3) jet_flavour[jj] = AODjet_Flav[lumis][event].at(dR_map[dR.at(0)]);
              else jet_flavour[jj] = 0;
               
              if(jet_flavour[jj] != 0) n_flav_jets++;
                
              dR.clear();
              dR_map.clear();
              
          }
        
          float sf[4]={1.,1.,1.,1.};
          float sf_err_down[4]={1.,1.,1.,1.}; 
          float sf_err_up[4]={1.,1.,1.,1.}; 
          float eff[4]={1.,1.,1.,1.}; 
          float eff_err[4]={1.,1.,1.,1.}; 
          float weight[4] = {1.,1.,1.,1.};
          float weight_err_up[4] = {1.,1.,1.,1.};
          float weight_err_down[4] = {1.,1.,1.,1.};
          float weight_err_05up[4] = {1.,1.,1.,1.};
          float weight_err_05down[4] = {1.,1.,1.,1.};
          float weight_err_15up[4] = {1.,1.,1.,1.};
          float weight_err_15down[4] = {1.,1.,1.,1.};
          float weight_err_20up[4] = {1.,1.,1.,1.};
          float weight_err_20down[4] = {1.,1.,1.,1.};
          
          for(int jj = 0; jj < 4; jj++){

              if(bad_jet[jj] == true) continue;
              
              if(abs(jet_flavour[jj]) == 5){
                 sf[jj] = SFb_CSVM->Eval(jetP4[jj]->Pt());
                 sf_err_up[jj] = h1_SFb_CSVM->GetBinError(h1_SFb_CSVM->FindBin(jetP4[jj]->Pt()));
                 sf_err_down[jj] = h1_SFb_CSVM->GetBinError(h1_SFb_CSVM->FindBin(jetP4[jj]->Pt()));
                 eff[jj] = h2_BTaggingEff_b_M->GetBinContent(h2_BTaggingEff_b_M->FindBin(jetP4[jj]->Pt(),fabs(jetP4[jj]->Eta())));
                 eff_err[jj] = h2_BTaggingEff_b_M->GetBinError(h2_BTaggingEff_b_M->FindBin(jetP4[jj]->Pt(),fabs(jetP4[jj]->Eta())));
              }
              if(abs(jet_flavour[jj]) == 4){
                 sf[jj] = SFb_CSVM->Eval(jetP4[jj]->Pt());
                 sf_err_up[jj] = 2*h1_SFb_CSVM->GetBinError(h1_SFb_CSVM->FindBin(jetP4[jj]->Pt()));
                 sf_err_down[jj] = 2*h1_SFb_CSVM->GetBinError(h1_SFb_CSVM->FindBin(jetP4[jj]->Pt()));
                 eff[jj] = h2_BTaggingEff_c_M->GetBinContent(h2_BTaggingEff_c_M->FindBin(jetP4[jj]->Pt(),fabs(jetP4[jj]->Eta())));
                 eff_err[jj] = h2_BTaggingEff_c_M->GetBinError(h2_BTaggingEff_c_M->FindBin(jetP4[jj]->Pt(),fabs(jetP4[jj]->Eta())));
              }
              if(abs(jet_flavour[jj]) == 1 || abs(jet_flavour[jj]) == 2 || abs(jet_flavour[jj]) == 3 || abs(jet_flavour[jj]) == 21){ 
                 eff[jj] = h2_BTaggingEff_udsg_M->GetBinContent(h2_BTaggingEff_udsg_M->FindBin(jetP4[jj]->Pt(),fabs(jetP4[jj]->Eta())));
                 eff_err[jj] = h2_BTaggingEff_udsg_M->GetBinError(h2_BTaggingEff_udsg_M->FindBin(jetP4[jj]->Pt(),fabs(jetP4[jj]->Eta())));   
                 if(fabs(jetP4[jj]->Eta()) > 0. && fabs(jetP4[jj]->Eta()) < 0.8){
                    sf[jj] = SFudsg_CSVM_00_08_mean->Eval(jetP4[jj]->Pt());
                    sf_err_up[jj] = SFudsg_CSVM_00_08_max->Eval(jetP4[jj]->Pt());
                    sf_err_down[jj] = SFudsg_CSVM_00_08_min->Eval(jetP4[jj]->Pt());
                 }
                 if(fabs(jetP4[jj]->Eta()) > 0.8 && fabs(jetP4[jj]->Eta()) < 1.6){
                    sf[jj] = SFudsg_CSVM_08_16_mean->Eval(jetP4[jj]->Pt());
                    sf_err_up[jj] = SFudsg_CSVM_08_16_max->Eval(jetP4[jj]->Pt());
                    sf_err_down[jj] = SFudsg_CSVM_08_16_min->Eval(jetP4[jj]->Pt());  
                 }
                 if(fabs(jetP4[jj]->Eta()) > 1.6 && fabs(jetP4[jj]->Eta()) < 2.4){
                    sf[jj] = SFudsg_CSVM_16_24_mean->Eval(jetP4[jj]->Pt());
                    sf_err_up[jj] = SFudsg_CSVM_16_24_max->Eval(jetP4[jj]->Pt());
                    sf_err_down[jj] = SFudsg_CSVM_16_24_min->Eval(jetP4[jj]->Pt());  
                 }
                 if(fabs(jetP4[jj]->Eta()) >= 2.4){
                    sf[jj] = SFudsg_CSVM_16_24_mean->Eval(jetP4[jj]->Pt());
                    sf_err_up[jj] = 2*SFudsg_CSVM_16_24_max->Eval(jetP4[jj]->Pt());
                    sf_err_down[jj] = 2*SFudsg_CSVM_16_24_min->Eval(jetP4[jj]->Pt());  
                 }     
              }
              
              if(jet_csvBtag[jj] > 0.679){
                 weight[jj] = sf[jj];
                 weight_err_up[jj] = sf_err_up[jj];
                 weight_err_down[jj] = sf_err_down[jj];
                 weight_err_05up[jj] = 0.5*sf_err_up[jj];
                 weight_err_05down[jj] = 0.5*sf_err_down[jj];
                 weight_err_15up[jj] = 1.5*sf_err_up[jj];
                 weight_err_15down[jj] = 1.5*sf_err_down[jj];
                 weight_err_20up[jj] = 2.*sf_err_up[jj];
                 weight_err_20down[jj] = 2.*sf_err_down[jj];
              }else{
                 weight[jj] = (1.-sf[jj]*eff[jj])/(1.-eff[jj]);
                 weight_err_up[jj] = jet_weight_err(eff[jj], eff_err[jj], sf[jj], sf_err_up[jj]);
                 weight_err_down[jj] = jet_weight_err(eff[jj], eff_err[jj], sf[jj], sf_err_down[jj]);
                 weight_err_05up[jj] = jet_weight_err(eff[jj], 0.5*eff_err[jj], sf[jj], 0.5*sf_err_up[jj]);
                 weight_err_05down[jj] = jet_weight_err(eff[jj], 0.5*eff_err[jj], sf[jj], 0.5*sf_err_down[jj]);
                 weight_err_15up[jj] = jet_weight_err(eff[jj], 1.5*eff_err[jj], sf[jj], 1.5*sf_err_up[jj]);
                 weight_err_15down[jj] = jet_weight_err(eff[jj], 1.5*eff_err[jj], sf[jj], 1.5*sf_err_down[jj]);
                 weight_err_20up[jj] = jet_weight_err(eff[jj], 2.*eff_err[jj], sf[jj], 2.*sf_err_up[jj]);
                 weight_err_20down[jj] = jet_weight_err(eff[jj], 2.*eff_err[jj], sf[jj], 2.*sf_err_down[jj]);

                 //if(sf_err_down[jj] != sf_err_up[jj])std::cout << ientry << " - " << jj << " - " << eff[jj] << " - " << eff_err[jj] << " - " << sf[jj] << " - " << sf_err_up[jj] << " - " << weight[jj] << " - " << weight_err_up[jj] << " - " << weight_err_down[jj] << " - " << sf_err_up[jj] << " - " << sf_err_down[jj] << std::endl;
              }
                
              if(fabs(jet_flavour[jj]) == 0){
                 weight[jj] = 1.;
                 weight_err_up[jj] = 1.;
                 weight_err_down[jj] = 1.;
              }
          }
            
          for(int jj = 0; jj < 4; jj++){
               
              if(bad_jet[jj] == true) continue;

              
              if(weight[jj] == 1 && weight_err_up[jj] == 1 && weight_err_down[jj] == 1) n_out++;
               
              if(jj == 0){
                 j1_pt_outtree = jetP4[jj]->Pt();
                 j1_e_outtree = jetP4[jj]->Energy();
                 j1_eta_outtree = jetP4[jj]->Eta();
                 j1_phi_outtree = jetP4[jj]->Phi();
                 j1_flavour = jet_flavour[jj];
                 j1_weight_max = weight[jj]+weight_err_up[jj];
                 j1_weight_05max = weight[jj]+weight_err_05up[jj];
                 j1_weight_15max = weight[jj]+weight_err_15up[jj];
                 j1_weight_20max = weight[jj]+weight_err_20up[jj];
                 j1_weight = weight[jj];
                 j1_weight_min = weight[jj]-weight_err_down[jj];
                 j1_weight_05min = weight[jj]-weight_err_05down[jj];
                 j1_weight_15min = weight[jj]-weight_err_15down[jj];
                 j1_weight_20min = weight[jj]-weight_err_20down[jj];
              }
              if(jj == 1){
                 j2_pt_outtree = jetP4[jj]->Pt();
                 j2_e_outtree = jetP4[jj]->Energy();
                 j2_eta_outtree = jetP4[jj]->Eta();
                 j2_phi_outtree = jetP4[jj]->Phi();
                 j2_flavour = jet_flavour[jj];
                 j2_weight_max = weight[jj]+weight_err_up[jj];
                 j2_weight_05max = weight[jj]+weight_err_05up[jj];
                 j2_weight_15max = weight[jj]+weight_err_15up[jj];
                 j2_weight_20max = weight[jj]+weight_err_20up[jj];
                 j2_weight = weight[jj];
                 j2_weight_min = weight[jj]-weight_err_down[jj];
                 j2_weight_05min = weight[jj]-weight_err_05down[jj];
                 j2_weight_15min = weight[jj]-weight_err_15down[jj];
                 j2_weight_20min = weight[jj]-weight_err_20down[jj];
              }
              if(jj == 2){
                 j3_pt_outtree = jetP4[jj]->Pt();
                 j3_e_outtree = jetP4[jj]->Energy();
                 j3_eta_outtree = jetP4[jj]->Eta();
                 j3_phi_outtree = jetP4[jj]->Phi();
                 j3_flavour = jet_flavour[jj];
                 j3_weight_max = weight[jj]+weight_err_up[jj];
                 j3_weight_05max = weight[jj]+weight_err_05up[jj];
                 j3_weight_15max = weight[jj]+weight_err_15up[jj];
                 j3_weight_20max = weight[jj]+weight_err_20up[jj];
                 j3_weight = weight[jj];
                 j3_weight_min = weight[jj]-weight_err_down[jj];
                 j3_weight_05min = weight[jj]-weight_err_05down[jj];
                 j3_weight_15min = weight[jj]-weight_err_15down[jj];
                 j3_weight_20min = weight[jj]-weight_err_20down[jj];
              }
              if(jj == 3){
                 j4_pt_outtree = jetP4[jj]->Pt();
                 j4_e_outtree = jetP4[jj]->Energy();
                 j4_eta_outtree = jetP4[jj]->Eta();
                 j4_phi_outtree = jetP4[jj]->Phi();
                 j4_flavour = jet_flavour[jj];
                 j4_weight_max = weight[jj]+weight_err_up[jj];
                 j4_weight_05max = weight[jj]+weight_err_05up[jj];
                 j4_weight_15max = weight[jj]+weight_err_15up[jj];
                 j4_weight_20max = weight[jj]+weight_err_20up[jj];
                 j4_weight = weight[jj];
                 j4_weight_min = weight[jj]-weight_err_down[jj];
                 j4_weight_05min = weight[jj]-weight_err_05down[jj];
                 j4_weight_15min = weight[jj]-weight_err_15down[jj];
                 j4_weight_20min = weight[jj]-weight_err_20down[jj];
              }
             
          }

          for(int ii = 0; ii < 4; ii++)
              bad_jet[ii] = false;
           
          for(int ii = 0; ii < 4; ii++)
              jet_csvBtag[ii] = -1;
          
          for(int ii = 0; ii < 4; ii++)
              jet_flavour[ii] = 0;

          outTree->Fill();
           
      } 
  }

  std::cout << "SELECTION: " << n_event << " - " << n_good_jets_total  << " - " << n_flav_jets << " - " << n_out << std::endl;
  
  TFile* output = new TFile(outputName.c_str(),"RECREATE");
  output->cd();
  outTree->Write();
  h2_BTaggingEff_b_M->Write();
  h2_BTaggingEff_c_M->Write();
  h2_BTaggingEff_udsg_M->Write();
  g_BTaggingEff_b_M_vs_pt_00_08_old->Write("g_BTaggingEff_b_M_vs_pt_00_08_old");
  g_BTaggingEff_b_M_vs_pt_00_08->Write("g_BTaggingEff_b_M_vs_pt_00_08");
  g_BTaggingEff_b_M_vs_pt_08_16->Write("g_BTaggingEff_b_M_vs_pt_08_16");
  g_BTaggingEff_b_M_vs_pt_16_25->Write("g_BTaggingEff_b_M_vs_pt_16_25");
  g_BTaggingEff_c_M_vs_pt_00_08->Write("g_BTaggingEff_c_M_vs_pt_00_08");
  g_BTaggingEff_c_M_vs_pt_08_16->Write("g_BTaggingEff_c_M_vs_pt_08_16");
  g_BTaggingEff_c_M_vs_pt_16_25->Write("g_BTaggingEff_c_M_vs_pt_16_25");
  g_BTaggingEff_udsg_M_vs_pt_00_08->Write("g_BTaggingEff_udsg_M_vs_pt_00_08");
  g_BTaggingEff_udsg_M_vs_pt_08_16->Write("g_BTaggingEff_udsg_M_vs_pt_08_16");
  g_BTaggingEff_udsg_M_vs_pt_16_25->Write("g_BTaggingEff_udsg_M_vs_pt_16_25");
  g_BTaggingEff_b_M_tot_vs_pt_00_08->Write("g_BTaggingEff_b_M_tot_vs_pt_00_08");
  g_BTaggingEff_b_M_tot_vs_pt_08_16->Write("g_BTaggingEff_b_M_tot_vs_pt_08_16");
  g_BTaggingEff_b_M_tot_vs_pt_16_25->Write("g_BTaggingEff_b_M_tot_vs_pt_16_25");
  g_BTaggingEff_c_M_tot_vs_pt_00_08->Write("g_BTaggingEff_c_M_tot_vs_pt_00_08");
  g_BTaggingEff_c_M_tot_vs_pt_08_16->Write("g_BTaggingEff_c_M_tot_vs_pt_08_16");
  g_BTaggingEff_c_M_tot_vs_pt_16_25->Write("g_BTaggingEff_c_M_tot_vs_pt_16_25");
  g_BTaggingEff_udsg_M_tot_vs_pt_00_08->Write("g_BTaggingEff_udsg_M_tot_vs_pt_00_08");
  g_BTaggingEff_udsg_M_tot_vs_pt_08_16->Write("g_BTaggingEff_udsg_M_tot_vs_pt_08_16");
  g_BTaggingEff_udsg_M_tot_vs_pt_16_25->Write("g_BTaggingEff_udsg_M_tot_vs_pt_16_25");
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

float jet_weight_err(float eff, float eff_err, float sf, float sf_err)
{
   float first = eff/(1-eff);
   float second = (1-sf)/((1-eff)*(1-eff));
   return sqrt(first*first*sf_err*sf_err+second*second*eff_err*eff_err);
}
















