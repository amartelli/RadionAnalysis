
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
#include "TGraphErrors.h"
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

  //  bool printOk = true;
  bool printOk = false;
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc <= 2)
  {
    std::cerr << ">>>>> Program usage: " << argv[0] << " configFileName btagWP" << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");
  std::string inputTree = gConfigParser -> readStringOption("Input::inputTree");
  
  bool useBJetRegression = gConfigParser -> readBoolOption("Input::useBJetRegression");
  std::string inputRegFile = gConfigParser -> readStringOption("Input::inputRegFile");
  
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");
  
  //  bool isReweighted = gConfigParser -> readBoolOption("Input::isReweighted");
  std::string inputMCPUHisto = gConfigParser -> readStringOption("Input::inputMCPUHisto");
  std::string inputDataPUHisto = gConfigParser -> readStringOption("Input::inputDataPUHisto");

  char* btagWP = argv[2];
   
  std::cout << "btagWP: " << btagWP << std::endl;

  if ( strcmp(btagWP,"No")!=0 && strcmp(btagWP,"Loose")!=0 && strcmp(btagWP,"Medium")!=0 && strcmp(btagWP,"Tight")!=0   )
  {
    std::cout << "CHK-STB Error: unknown  btag working point " << btagWP << std::endl;
    std::cout << "CHK-STB Select among No, Loose, Medium and Tight! " << std::endl;
    return -1;
  }

  float btag_wp = -1;
  
  if ( strcmp(btagWP,"No") == 0 ) btag_wp = 0.;
  if ( strcmp(btagWP,"Loose") == 0 ) btag_wp = 0.244;
  if ( strcmp(btagWP,"Medium") == 0 ) btag_wp = 0.679;
  if ( strcmp(btagWP,"Tight") == 0 ) btag_wp = 0.898;

  std::map<int,TChain*> ntu;
  std::map<int, TFile*> Files;

  std::map<int, std::map<int,float> > weight;
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
    if(printOk)    std::cout << " >>>> pos = " << pos << std::endl;
  }
  pos_total = pos;

  TTree *newtree;

  for(int ii = 0; ii < pos_total; ii++){
      if(ntu[ii]->GetEntries() == 0 )
      {
         std::cout << "Error: input file" << ii << " is empty" << std::endl; 
         return -1;
      }
      else newtree = ntu[ii]->CloneTree(0);
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
  if(useBJetRegression){
  readerRegres->AddVariable( "jet_pt",            &fRegr_pt);
  readerRegres->AddVariable( "jet_eta",           &fRegr_eta);
  readerRegres->AddVariable( "jet_emfrac",        &fRegr_cef);
  readerRegres->AddVariable( "jet_nConstituents", &fRegr_nconst);
  readerRegres->AddVariable( "jet_hadfrac",       &fRegr_chf);
  readerRegres->AddVariable( "jet_secVtxPt",      &fRegr_vtxPt);
  readerRegres->AddVariable( "jet_secVtx3dL",     &fRegr_vtx3dl);
  readerRegres->AddVariable( "ev_met_corr_pfmet",  &fRegr_met);
  readerRegres->AddVariable( "jet_dPhiMet",        &fRegr_dPhiMet);
  readerRegres->BookMVA("BDT",inputRegFile.c_str());
  }
  
  float         pu_n;
  float         weight_an;
  float         evweight;
  float         pu_weight;

  float         met_corr_pfmet;
  float         met_corr_phi_pfmet;
  float         met_corr_eta_pfmet;
  float         met_corr_e_pfmet;

  int           ph1_ciclevel;
  int           ph2_ciclevel;
  float         ph1_pfchargedisogood03;
  float         ph2_pfchargedisogood03;
  float         ph1_pfchargedisobad03;
  float         ph2_pfchargedisobad03;
  float         ph1_pfchargedisobad04;
  float         ph2_pfchargedisobad04;
  float         ph1_ecaliso;
  float         ph2_ecaliso;
  float         ph1_ecalisobad;
  float         ph2_ecalisobad;

  float         ph1_badvtx_Et;
  float         ph2_badvtx_Et;
                
  int           category;
  int           ph1_isEB;
  int           ph2_isEB;
  float         rho;

  float         ph1_e;
  float         ph2_e;
  float         ph1_pt;
  float         ph2_pt;
  float         ph1_phi;
  float         ph2_phi;
  float         ph1_eta;
  float         ph2_eta;
  float         ph1_SCEta;
  float         ph2_SCEta;
  float         ph1_r9;
  float         ph2_r9;
  float         ph1_hoe;
  float         ph2_hoe;
  float         ph1_sieie;
  float         ph2_sieie;
  float         ph1_sieip;
  float         ph2_sieip;
  float         ph1_sipip;
  float         ph2_sipip;
  //  float         ph1_isconv;
  //  float         ph2_isconv;
  float         PhotonsMass;
  float         dipho_E;
  float         dipho_pt;
  float         dipho_eta;
  float         dipho_phi;
  float         dipho_cosThetaStar_CS;
  float         dipho_tanhYStar;
  float         dipho_Y;
  int           vtx_ind;
  float         j1_e;
  float         j1_pt;
  float         j1_phi;
  float         j1_csvBtag;
  float         j1_jetProbBtag; 
  float         j1_csvMvaBtag;
  float         j1_tcheBtag;
  float         j1_emfrac;
  float         j1_hadfrac;
  float         j1_secVtxPt;
  float         j1_secVtx3dL;
  int           j1_nNeutrals;
  int           j1_nCharged;
  float         j2_e;
  float         j2_pt;
  float         j2_phi;
  float         j2_csvBtag;
  float         j2_jetProbBtag;
  float         j2_csvMvaBtag;
  float         j2_tcheBtag;
  float         j2_emfrac;
  float         j2_secVtxPt;
  float         j2_secVtx3dL;
  int           j2_nNeutrals;
  int           j2_nCharged;
  float         j2_hadfrac;
  float         j3_e;
  float         j3_pt;
  float         j3_phi;
  float         j3_csvBtag;
  float         j3_jetProbBtag;
  float         j3_csvMvaBtag;
  float         j3_tcheBtag;
  float         j3_emfrac;
  float         j3_hadfrac;
  float         j3_secVtxPt;
  float         j3_secVtx3dL;
  int           j3_nNeutrals;
  int           j3_nCharged;
  float         j4_e;
  float         j4_pt;
  float         j4_phi;
  float         j4_csvBtag;
  float         j4_jetProbBtag;
  float         j4_csvMvaBtag;
  float         j4_tcheBtag;
  float         j4_emfrac;
  float         j4_hadfrac;
  float         j4_secVtxPt;
  float         j4_secVtx3dL;
  int           j4_nNeutrals;
  int           j4_nCharged;
  float         JetsMass;
  float         dijet_E;
//   float         dijet_pt;
//   float         dijet_eta;
//   float         dijet_phi;
//  float         RadMass;
  float         radion_E;
//   float         radion_pt;
//   float         radion_eta;
//   float         radion_phi;
  int           njets_passing_kLooseID;
  
  for(int ii = 0; ii < pos_total; ii++){
      ntu[ii]->SetBranchStatus("*",0);
      ntu[ii]->SetBranchStatus("pu_n",1);
      ntu[ii]->SetBranchStatus("weight",1);
      ntu[ii]->SetBranchStatus("evweight",1);
      ntu[ii]->SetBranchStatus("pu_weight",1);
      ntu[ii]->SetBranchStatus("nvtx",1);

      ntu[ii]->SetBranchStatus("met_corr_pfmet",1);
      ntu[ii]->SetBranchStatus("met_corr_phi_pfmet",1);
      ntu[ii]->SetBranchStatus("met_corr_eta_pfmet",1);
      ntu[ii]->SetBranchStatus("met_corr_e_pfmet",1);

      ntu[ii]->SetBranchStatus("ph1_ciclevel",1);
      ntu[ii]->SetBranchStatus("ph2_ciclevel",1);
      ntu[ii]->SetBranchStatus("ph1_pfchargedisogood03",1);
      ntu[ii]->SetBranchStatus("ph2_pfchargedisogood03",1);
      ntu[ii]->SetBranchStatus("ph1_pfchargedisobad03",1);
      ntu[ii]->SetBranchStatus("ph2_pfchargedisobad03",1);
      ntu[ii]->SetBranchStatus("ph1_ecaliso",1);
      ntu[ii]->SetBranchStatus("ph2_ecaliso",1);
      ntu[ii]->SetBranchStatus("ph1_ecalisobad",1);
      ntu[ii]->SetBranchStatus("ph2_ecalisobad",1);

      ntu[ii]->SetBranchStatus("ph1_badvtx_Et",1);
      ntu[ii]->SetBranchStatus("ph2_badvtx_Et",1);
      ntu[ii]->SetBranchStatus("ph1_pfchargedisobad04",1);
      ntu[ii]->SetBranchStatus("ph2_pfchargedisobad04",1);
      ntu[ii]->SetBranchStatus("category",1);
      ntu[ii]->SetBranchStatus("ph1_isEB",1);
      ntu[ii]->SetBranchStatus("ph2_isEB",1);
      ntu[ii]->SetBranchStatus("rho",1);

      ntu[ii]->SetBranchStatus("ph1_e",1);
      ntu[ii]->SetBranchStatus("ph2_e",1);
      ntu[ii]->SetBranchStatus("ph1_pt",1);
      ntu[ii]->SetBranchStatus("ph2_pt",1);
      ntu[ii]->SetBranchStatus("ph1_phi",1);
      ntu[ii]->SetBranchStatus("ph2_phi",1);
      ntu[ii]->SetBranchStatus("ph1_eta",1);
      ntu[ii]->SetBranchStatus("ph2_eta",1);
      ntu[ii]->SetBranchStatus("ph1_SCEta",1);
      ntu[ii]->SetBranchStatus("ph2_SCEta",1);
      ntu[ii]->SetBranchStatus("ph1_r9",1);
      ntu[ii]->SetBranchStatus("ph2_r9",1);
      ntu[ii]->SetBranchStatus("ph1_hoe",1);
      ntu[ii]->SetBranchStatus("ph2_hoe",1);
      ntu[ii]->SetBranchStatus("ph1_sieie",1);
      ntu[ii]->SetBranchStatus("ph2_sieie",1);
      ntu[ii]->SetBranchStatus("ph1_sieip",1);
      ntu[ii]->SetBranchStatus("ph2_sieip",1);
      ntu[ii]->SetBranchStatus("ph1_sipip",1);
      ntu[ii]->SetBranchStatus("ph2_sipip",1);
      ntu[ii]->SetBranchStatus("PhotonsMass",1);
      ntu[ii]->SetBranchStatus("dipho_E",1);
      ntu[ii]->SetBranchStatus("dipho_pt",1);
      ntu[ii]->SetBranchStatus("dipho_eta",1);
      ntu[ii]->SetBranchStatus("dipho_phi",1);
      ntu[ii]->SetBranchStatus("dipho_cosThetaStar_CS",1);
      ntu[ii]->SetBranchStatus("dipho_tanhYStar",1);
      ntu[ii]->SetBranchStatus("dipho_Y",1);
      ntu[ii]->SetBranchStatus("vtx_ind",1);
      ntu[ii]->SetBranchStatus("j1_e",1);
      ntu[ii]->SetBranchStatus("j1_pt",1);
      ntu[ii]->SetBranchStatus("j1_phi",1);
      ntu[ii]->SetBranchStatus("j1_eta",1);
      ntu[ii]->SetBranchStatus("j1_csvBtag",1);
      ntu[ii]->SetBranchStatus("j1_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j1_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j1_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j1_tcheBtag",1);
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
      ntu[ii]->SetBranchStatus("j2_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j2_tcheBtag",1);
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
      ntu[ii]->SetBranchStatus("j3_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j3_tcheBtag",1);
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
      ntu[ii]->SetBranchStatus("j4_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j4_tcheBtag",1);
      ntu[ii]->SetBranchStatus("j4_emfrac",1);
      ntu[ii]->SetBranchStatus("j4_hadfrac",1);
      ntu[ii]->SetBranchStatus("j4_secVtxPt",1);
      ntu[ii]->SetBranchStatus("j4_secVtx3dL",1);
      ntu[ii]->SetBranchStatus("j4_nNeutrals",1);
      ntu[ii]->SetBranchStatus("j4_nCharged",1);
      ntu[ii]->SetBranchStatus("j4_dR2Mean",1);
      ntu[ii]->SetBranchStatus("JetsMass",1);
      ntu[ii]->SetBranchStatus("dijet_E",1);
//       ntu[ii]->SetBranchStatus("dijet_pt",1);
//       ntu[ii]->SetBranchStatus("dijet_phi",1);
//       ntu[ii]->SetBranchStatus("dijet_eta",1);
      ntu[ii]->SetBranchStatus("radion_E",1);
//       ntu[ii]->SetBranchStatus("radion_pt",1);
//       ntu[ii]->SetBranchStatus("radion_eta",1);
//       ntu[ii]->SetBranchStatus("radion_phi",1);
      ntu[ii]->SetBranchStatus("njets_passing_kLooseID",1);
      ntu[ii]->SetBranchAddress("pu_n",&pu_n);
      ntu[ii]->SetBranchAddress("weight",&weight_an);
      ntu[ii]->SetBranchAddress("evweight",&evweight);
      ntu[ii]->SetBranchAddress("pu_weight",&pu_weight);
      ntu[ii]->SetBranchAddress("nvtx",&nvtx);

      ntu[ii]->SetBranchAddress("met_corr_pfmet",&met_corr_pfmet);
      ntu[ii]->SetBranchAddress("met_corr_phi_pfmet",&met_corr_phi_pfmet);
      ntu[ii]->SetBranchAddress("met_corr_eta_pfmet",&met_corr_eta_pfmet);
      ntu[ii]->SetBranchAddress("met_corr_e_pfmet",&met_corr_e_pfmet);

      ntu[ii]->SetBranchAddress("ph1_ciclevel",&ph1_ciclevel);
      ntu[ii]->SetBranchAddress("ph2_ciclevel",&ph2_ciclevel);
      ntu[ii]->SetBranchAddress("ph1_pfchargedisogood03",&ph1_pfchargedisogood03);
      ntu[ii]->SetBranchAddress("ph2_pfchargedisogood03",&ph2_pfchargedisogood03);
      ntu[ii]->SetBranchAddress("ph1_pfchargedisobad03",&ph1_pfchargedisobad03);
      ntu[ii]->SetBranchAddress("ph2_pfchargedisobad03",&ph2_pfchargedisobad03);
      ntu[ii]->SetBranchAddress("ph1_ecaliso",&ph1_ecaliso);
      ntu[ii]->SetBranchAddress("ph2_ecaliso",&ph2_ecaliso);
      ntu[ii]->SetBranchAddress("ph1_ecalisobad",&ph1_ecalisobad);
      ntu[ii]->SetBranchAddress("ph2_ecalisobad",&ph2_ecalisobad);

      ntu[ii]->SetBranchAddress("ph1_badvtx_Et",&ph1_badvtx_Et);
      ntu[ii]->SetBranchAddress("ph2_badvtx_Et",&ph2_badvtx_Et);
      ntu[ii]->SetBranchAddress("ph1_pfchargedisobad04",&ph1_pfchargedisobad04);
      ntu[ii]->SetBranchAddress("ph2_pfchargedisobad04",&ph2_pfchargedisobad04);
      ntu[ii]->SetBranchAddress("category",&category);
      ntu[ii]->SetBranchAddress("ph1_isEB",&ph1_isEB);
      ntu[ii]->SetBranchAddress("ph2_isEB",&ph2_isEB);
      ntu[ii]->SetBranchAddress("rho",&rho);

      ntu[ii]->SetBranchAddress("ph1_e",&ph1_e);
      ntu[ii]->SetBranchAddress("ph2_e",&ph2_e);
      ntu[ii]->SetBranchAddress("ph1_pt",&ph1_pt);
      ntu[ii]->SetBranchAddress("ph2_pt",&ph2_pt);
      ntu[ii]->SetBranchAddress("ph1_phi",&ph1_phi);
      ntu[ii]->SetBranchAddress("ph2_phi",&ph2_phi);
      ntu[ii]->SetBranchAddress("ph1_eta",&ph1_eta);
      ntu[ii]->SetBranchAddress("ph2_eta",&ph2_eta);
      ntu[ii]->SetBranchAddress("ph1_SCEta",&ph1_SCEta);
      ntu[ii]->SetBranchAddress("ph2_SCEta",&ph2_SCEta);
      ntu[ii]->SetBranchAddress("ph1_r9",&ph1_r9);
      ntu[ii]->SetBranchAddress("ph2_r9",&ph2_r9);
      ntu[ii]->SetBranchAddress("ph1_hoe",&ph1_hoe);
      ntu[ii]->SetBranchAddress("ph2_hoe",&ph2_hoe);
      ntu[ii]->SetBranchAddress("ph1_sieie",&ph1_sieie);
      ntu[ii]->SetBranchAddress("ph2_sieie",&ph2_sieie);
      ntu[ii]->SetBranchAddress("ph1_sieip",&ph1_sieip);
      ntu[ii]->SetBranchAddress("ph2_sieip",&ph2_sieip);
      ntu[ii]->SetBranchAddress("ph1_sipip",&ph1_sipip);
      ntu[ii]->SetBranchAddress("ph2_sipip",&ph2_sipip);
      ntu[ii]->SetBranchAddress("PhotonsMass",&PhotonsMass);
      ntu[ii]->SetBranchAddress("dipho_E",&dipho_E);
      ntu[ii]->SetBranchAddress("dipho_pt",&dipho_pt);
      ntu[ii]->SetBranchAddress("dipho_eta",&dipho_eta);
      ntu[ii]->SetBranchAddress("dipho_phi",&dipho_phi);
      ntu[ii]->SetBranchAddress("dipho_cosThetaStar_CS",&dipho_cosThetaStar_CS);
      ntu[ii]->SetBranchAddress("dipho_tanhYStar",&dipho_tanhYStar);
      ntu[ii]->SetBranchAddress("dipho_Y",&dipho_Y);
      ntu[ii]->SetBranchAddress("vtx_ind",&vtx_ind);
      ntu[ii]->SetBranchAddress("j1_e",&j1_e);
      ntu[ii]->SetBranchAddress("j1_pt",&j1_pt);
      ntu[ii]->SetBranchAddress("j1_phi",&j1_phi);
      ntu[ii]->SetBranchAddress("j1_eta",&j1_eta);
      ntu[ii]->SetBranchAddress("j1_csvBtag",&j1_csvBtag);
      ntu[ii]->SetBranchAddress("j1_betaStarClassic",&j1_betaStarClassic);
      ntu[ii]->SetBranchAddress("j1_jetProbBtag",&j1_jetProbBtag);
      ntu[ii]->SetBranchAddress("j1_csvMvaBtag",&j1_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j1_tcheBtag",&j1_tcheBtag);
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
      ntu[ii]->SetBranchAddress("j2_csvMvaBtag",&j2_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j2_tcheBtag",&j2_tcheBtag);
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
      ntu[ii]->SetBranchAddress("j3_csvMvaBtag",&j3_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j3_tcheBtag",&j3_tcheBtag);
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
      ntu[ii]->SetBranchAddress("j4_csvMvaBtag",&j4_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j4_tcheBtag",&j4_tcheBtag);
      ntu[ii]->SetBranchAddress("j4_emfrac",&j4_emfrac);
      ntu[ii]->SetBranchAddress("j4_hadfrac",&j4_hadfrac);
      ntu[ii]->SetBranchAddress("j4_secVtxPt",&j4_secVtxPt);
      ntu[ii]->SetBranchAddress("j4_secVtx3dL",&j4_secVtx3dL);
      ntu[ii]->SetBranchAddress("j4_nNeutrals",&j4_nNeutrals);
      ntu[ii]->SetBranchAddress("j4_nCharged",&j4_nCharged);
      ntu[ii]->SetBranchAddress("j4_dR2Mean",&j4_dR2Mean);
      ntu[ii]->SetBranchAddress("JetsMass",&JetsMass);
      ntu[ii]->SetBranchAddress("dijet_E",&dijet_E);
//       ntu[ii]->SetBranchAddress("dijet_pt",&dijet_pt);
//       ntu[ii]->SetBranchAddress("dijet_phi",&dijet_phi);
//       ntu[ii]->SetBranchAddress("dijet_eta",&dijet_eta);
      ntu[ii]->SetBranchAddress("radion_E",&radion_E);
//       ntu[ii]->SetBranchAddress("radion_pt",&radion_pt);
//       ntu[ii]->SetBranchAddress("radion_eta",&radion_eta);
//       ntu[ii]->SetBranchAddress("radion_phi",&radion_phi);
      ntu[ii]->SetBranchAddress("njets_passing_kLooseID",&njets_passing_kLooseID);   
  }

  int nSTEP = 20;

  std::vector<TH1F*>  h_vtx_ind;
  std::vector<TH1F*>  h_PhoEcalIso;
  std::vector<TH1F*>  h_PhoChargedIso;
  std::vector<TH1F*>  h_PhoEcalIsoBad;
  std::vector<TH1F*>  h_PhoChargedIsoBad;
  std::vector<TH1F*>  h_ciclevel;

  std::vector<TH1F*>  h_PhoR9;
  std::vector<TH1F*>  h_PhoHoE;
  std::vector<TH1F*>  h_PhoSieie;
  std::vector<TH1F*>  h_PhoSipip;
  std::vector<TH1F*>  h_PhoSieip;
  std::vector<TH1F*>  h_PhoPt;
  std::vector<TH1F*>  h_PhoEta;
  std::vector<TH1F*>  h_PhoPhi;
  std::vector<TH1F*>  h_diPhoDeltaEta;
  std::vector<TH1F*>  h_diPhoDeltaPhi;
  std::vector<TH1F*>  h_diPhoDeltaR;
  std::vector<TH1F*>  h_diPhoInvMass;
  std::vector<TH1F*>  h_diPhoInvMass_CR1;
  std::vector<TH1F*>  h_diPhoInvMass_CR2;
  std::vector<TH1F*>  h_JetNum;
  std::vector<TH1F*>  h_JetNum_btagged;
  std::vector<TH1F*>  h_JetcsvBtag;
  std::vector<TH1F*>  h_JetcsvMvaBtag;
  std::vector<TH1F*>  h_JetProbBtag;
  std::vector<TH1F*>  h_JettcheBtag;
  std::vector<TH1F*>  h_JetPt;
  std::vector<TH1F*>  h_JetEta;
  std::vector<TH1F*>  h_JetPhi;
  std::vector<TH1F*>  h_diJetDeltaEta;
  std::vector<TH1F*>  h_diJetDeltaPhi;
  std::vector<TH1F*>  h_diJetDeltaR;
  std::vector<TH1F*>  h_diJetInvMass;
  std::vector<TH1F*>  h_diJetDeltaEta_maxpt;
  std::vector<TH1F*>  h_diJetDeltaPhi_maxpt;
  std::vector<TH1F*>  h_diJetDeltaR_maxpt;
  std::vector<TH1F*>  h_diJetInvMass_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaEta;
  std::vector<TH1F*>  h_PhoJetDeltaPhi;
  std::vector<TH1F*>  h_PhoJetDeltaR;
  std::vector<TH1F*>  h_PhoJetDeltaEta_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaPhi_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaR_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaEta_maxpt_btag;
  std::vector<TH1F*>  h_PhoJetDeltaPhi_maxpt_btag;
  std::vector<TH1F*>  h_PhoJetDeltaR_maxpt_btag;
  std::vector<TH1F*>  h_diPhodiJetDeltaEta;
  std::vector<TH1F*>  h_diPhodiJetDeltaPhi;
  std::vector<TH1F*>  h_diPhodiJetDeltaR;
  std::vector<TH1F*>  h_diPhodiJetInvMass;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands;
  std::vector<TH1F*>  h_diPhodiJetInvMass_maxpt;
  std::vector<TH1F*>  h_diPhodiJetInvMass_maxpt_CR1;
  std::vector<TH1F*>  h_diPhodiJetInvMass_maxpt_CR2;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core_maxpt;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands_maxpt;
  std::vector<TH1F*>  h_diPhodiJetInvMass_btag;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core_btag;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands_btag;
  std::vector<TH1F*>  h_diPhodiJetInvMass_bveto;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core_bveto;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands_bveto;
  std::vector<TH1F*>  h_diPhodiJetInvMass_kinfit;

  /////////////////////////////////////////////////////
  for(int nstep = 0; nstep<nSTEP; ++nstep){
    char histoName[50];
    sprintf(histoName, "_step%d", nstep);
    std::string stepName = std::string(histoName);
    h_vtx_ind.push_back(new TH1F(("h_vtx_ind"+stepName).c_str(), "", 40, 0., 40));
    h_PhoEcalIso.push_back(new TH1F(("h_PhoEcalIso"+stepName).c_str(),"",100,0.,20.));
    h_PhoChargedIso.push_back(new TH1F(("h_PhoChargedIso"+stepName).c_str(),"",100,0.,20.));
    h_PhoEcalIsoBad.push_back(new TH1F(("h_PhoEcalIsoBad"+stepName).c_str(),"",100,0.,20.));
    h_PhoChargedIsoBad.push_back(new TH1F(("h_PhoChargedIsoBad"+stepName).c_str(),"",100,0.,20.));
    h_ciclevel.push_back(new TH1F(("h_ciclevel"+stepName).c_str(),"",100,0.,20.));
    
    h_PhoR9.push_back(new TH1F(("h_PhoR9"+stepName).c_str(),"",50,0.2,1.));
    h_PhoHoE.push_back(new TH1F(("h_PhoHoE"+stepName).c_str(),"",50,0.,0.1));
    h_PhoSieie.push_back(new TH1F(("h_PhoSieie"+stepName).c_str(),"",50,0.,0.07));
    h_PhoSipip.push_back(new TH1F(("h_PhoSipip"+stepName).c_str(),"",50,0.,0.0045));
    h_PhoSieip.push_back(new TH1F(("h_PhoSieip"+stepName).c_str(),"",50,0.,0.0004));
    h_PhoPt.push_back(new TH1F(("h_PhoPt"+stepName).c_str(),"", 180,0.,180.));
    h_PhoEta.push_back(new TH1F(("h_PhoEta"+stepName).c_str(),"",700,-3.5,3.5));
    h_PhoPhi.push_back(new TH1F(("h_PhoPhi"+stepName).c_str(),"",700,-3.5,3.5));
    h_diPhoDeltaEta.push_back(new TH1F(("h_diPhoDeltaEta"+stepName).c_str(),"", 1000,-5.,5.));
    h_diPhoDeltaPhi.push_back(new TH1F(("h_diPhoDeltaPhi"+stepName).c_str(),"", 1000,-4.,4.));
    h_diPhoDeltaR.push_back(new TH1F(("h_diPhoDeltaR"+stepName).c_str(),"", 60,0.,6.));
    h_diPhoInvMass.push_back(new TH1F(("h_diPhoInvMass"+stepName).c_str(),"", 100,100.,300.));
    h_diPhoInvMass_CR1.push_back(new TH1F(("h_diPhoInvMass_CR1"+stepName).c_str(),"", 100,100.,300.));
    h_diPhoInvMass_CR2.push_back(new TH1F(("h_diPhoInvMass_CR2"+stepName).c_str(),"", 100,100.,300.));
    h_JetNum.push_back(new TH1F(("h_JetNum"+stepName).c_str(),"", 10,0.,10.));
    h_JetNum_btagged.push_back(new TH1F(("h_JetNum_btagged"+stepName).c_str(),"", 10,0.,10.));
    h_JetcsvBtag.push_back(new TH1F(("h_JetcsvBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JetcsvMvaBtag.push_back(new TH1F(("h_JetcsvMvaBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JetProbBtag.push_back(new TH1F(("h_JetProbBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JettcheBtag.push_back(new TH1F(("h_JettcheBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JetPt.push_back(new TH1F(("h_JetPt"+stepName).c_str(),"", 180,0.,180.));
    h_JetEta.push_back(new TH1F(("h_JetEta"+stepName).c_str(),"",800,-4,4));
    h_JetPhi.push_back(new TH1F(("h_JetPhi"+stepName).c_str(),"",700,-3.5,3.5));
    h_diJetDeltaEta.push_back(new TH1F(("h_diJetDeltaEta"+stepName).c_str(),"", 100,-7.,7.));
    h_diJetDeltaPhi.push_back(new TH1F(("h_diJetDeltaPhi"+stepName).c_str(),"", 100,-4.,4.));
    h_diJetDeltaR.push_back(new TH1F(("h_diJetDeltaR"+stepName).c_str(),"", 70,0.,7.));
    h_diJetInvMass.push_back(new TH1F(("h_diJetInvMass"+stepName).c_str(),"", 80,0.,800));
    h_diJetDeltaEta_maxpt.push_back(new TH1F(("h_diJetDeltaEta_maxpt"+stepName).c_str(),"", 100,-7.,7.));
    h_diJetDeltaPhi_maxpt.push_back(new TH1F(("h_diJetDeltaPhi_maxpt"+stepName).c_str(),"", 100,-4.,4.));
    h_diJetDeltaR_maxpt.push_back(new TH1F(("h_diJetDeltaR_maxpt"+stepName).c_str(),"",70,0.,7.));
    h_diJetInvMass_maxpt.push_back(new TH1F(("h_diJetInvMass_maxpt"+stepName).c_str(),"", 80,0.,800));
    h_PhoJetDeltaEta.push_back(new TH1F(("h_PhoJetDeltaEta"+stepName).c_str(),"",100,-5.,5.));
    h_PhoJetDeltaPhi.push_back(new TH1F(("h_PhoJetDeltaPhi"+stepName).c_str(),"",100,-4.,4.));
    h_PhoJetDeltaR.push_back(new TH1F(("h_PhoJetDeltaR"+stepName).c_str(),"",60,0.,6.));
    h_PhoJetDeltaEta_maxpt.push_back(new TH1F(("h_PhoJetDeltaEta_maxpt"+stepName).c_str(),"",100,-5.,5.));
    h_PhoJetDeltaPhi_maxpt.push_back(new TH1F(("h_PhoJetDeltaPhi_maxpt"+stepName).c_str(),"",100,-4.,4.));
    h_PhoJetDeltaR_maxpt.push_back(new TH1F(("h_PhoJetDeltaR_maxpt"+stepName).c_str(),"",60,0.,6.));
    h_PhoJetDeltaEta_maxpt_btag.push_back(new TH1F(("h_PhoJetDeltaEta_maxpt_btag"+stepName).c_str(),"",100,-5.,5.));
    h_PhoJetDeltaPhi_maxpt_btag.push_back(new TH1F(("h_PhoJetDeltaPhi_maxpt_btag"+stepName).c_str(),"",100,-4.,4.));
    h_PhoJetDeltaR_maxpt_btag.push_back(new TH1F(("h_PhoJetDeltaR_maxpt_btag"+stepName).c_str(),"",60,0.,6.));
    h_diPhodiJetDeltaEta.push_back(new TH1F(("h_diPhodiJetDeltaEta"+stepName).c_str(),"",120,-9.,9.));
    h_diPhodiJetDeltaPhi.push_back(new TH1F(("h_diPhodiJetDeltaPhi"+stepName).c_str(),"",100,-4.,4.));
    h_diPhodiJetDeltaR.push_back(new TH1F(("h_diPhodiJetDeltaR"+stepName).c_str(),"",90,0.,9.));
    h_diPhodiJetInvMass.push_back(new TH1F(("h_diPhodiJetInvMass"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core.push_back(new TH1F(("h_diPhodiJetInvMass_core"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_maxpt.push_back(new TH1F(("h_diPhodiJetInvMass_maxpt"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_maxpt_CR1.push_back(new TH1F(("h_diPhodiJetInvMass_maxpt_CR1"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_maxpt_CR2.push_back(new TH1F(("h_diPhodiJetInvMass_maxpt_CR2"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core_maxpt.push_back(new TH1F(("h_diPhodiJetInvMass_core_maxpt"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands_maxpt.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands_maxpt"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_btag.push_back(new TH1F(("h_diPhodiJetInvMass_btag"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core_btag.push_back(new TH1F(("h_diPhodiJetInvMass_core_btag"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands_btag.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands_btag"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_bveto.push_back(new TH1F(("h_diPhodiJetInvMass_bveto"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core_bveto.push_back(new TH1F(("h_diPhodiJetInvMass_core_bveto"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands_bveto.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands_bveto"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_kinfit.push_back(new TH1F(("h_diPhodiJetInvMass_kinfit"+stepName).c_str(),"",120,0.,1200));
  }
  /////////////////////////////////////////////////////

  //  TTree *newtree =  ntu[ii]->>CloneTree(0);
  
  std::map<int,TLorentzVector*> phoP4;
  std::map<int,TLorentzVector*> jetP4; 
  std::map<int,float> jet_csvBtag;
  std::map<int,float> jet_csvMvaBtag;
  std::map<int,float> jet_jetProbBtag;
  std::map<int,float> jet_tcheBtag;
  TLorentzVector* sumP4;
  TLorentzVector* PhoSumP4;
  std::map<int,TLorentzVector*> JetSumP4;

  std::map<int,bool> bad_jet_pre;
  std::map<int,bool> bad_jet;
  
  std::map<int,bool> isBTagCouple;
  std::map<int,bool> isBVetoCouple;

  std::vector<float>* pho_pt;
  std::vector<float>* jet_pt;

  std::map<float,int> phoPt_map;
  std::map<float,int> jetPt_map;

  int jet_couple_num; 
  int numPassingEvents[nSTEP];
  int numPassingEvents_CR1[nSTEP];
  int numPassingEvents_CR2[nSTEP];
  for(int nLocStep = 0; nLocStep<nSTEP; ++nLocStep) {
    numPassingEvents[nLocStep] = 0;
    numPassingEvents_CR1[nLocStep] = 0;
    numPassingEvents_CR2[nLocStep] = 0;
  }
  std::vector<int> event_selSTEP;
  std::vector<int> event_selSTEP_CR1;
  std::vector<int> event_selSTEP_CR2;

  std::vector<float> jjPt;
  std::map<float,std::pair<int,int> > jjPt_map;

  TLorentzVector* Hgg;
  TLorentzVector* Hjj; 
  TLorentzVector* Radion;

  float Hmass = 125.;
  DiJetKinFitter* fitter_jetsH = new DiJetKinFitter( "fitter_jetsH", "fitter_jets", Hmass );
 
  for(int nn = 0; nn < pos_total; nn++){

  std::cout << " ############################ " << std::endl;
  std::cout << " ################### Events = " << ntu[nn]->GetEntries() << std::endl;
  std::cout << " ############################ " << std::endl;

      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%1000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);

// 	  if(ientry < 31490) continue;
// 	  if(ientry > 31490) return 10;

          bool isfinal_event = true;

          int n_good_jets_pre = 0;
          int n_good_jets = 0;
          int n_b_jets_pre = 0;
          int n_b_jets = 0;

          bool isbad = false;
	  int Step_Counter = 0;
// 	  int Step_Counter_CR1 = 0;
// 	  int Step_Counter_CR2 = 0;
	  if(printOk) 	  std::cout << " Step_Counter = " << Step_Counter << std::endl;
// 	  if(printOk) 	  std::cout << " Step_Counter_CR1 = " << Step_Counter_CR1 << std::endl;
// 	  if(printOk) 	  std::cout << " Step_Counter_CR2 = " << Step_Counter_CR2 << std::endl;

	  if(printOk)	  std::cout << "--- Reading file_" << nn << " entry = " << ientry << std::endl;
	  //	  if(ientry > 10)return 10;

          float puRe = 1.;
          puRe = evweight;

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
          if(useBJetRegression){
           
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

	  //ALLEVTS (0)
	  ++numPassingEvents[Step_Counter];
          event_selSTEP.push_back(Step_Counter);  
	  ++numPassingEvents_CR1[Step_Counter];
          event_selSTEP_CR1.push_back(Step_Counter);  
	  ++numPassingEvents_CR2[Step_Counter];
          event_selSTEP_CR2.push_back(Step_Counter);  

          pho_pt = new std::vector<float>;
          jet_pt = new std::vector<float>;


	  //superTight
	  //EB 
	  bool diPhotonSuperTIGHT_passe = false;
	  if(ph1_ciclevel > 3 && ph2_ciclevel > 3) diPhotonSuperTIGHT_passe = true;

	  /*
	  bool Photon1SuperTIGHT_passe = false;
	  bool Photon2SuperTIGHT_passe = false;
	  //	  std::cout << " category = " << category << std::endl;
	  float val_isosumoet_1 = (ph1_pfchargedisogood03 + ph1_ecaliso + 2.5 - rho * 0.09) * 50. / ph1_pt;
	  float val_isosumoetbad_1 = (ph1_pfchargedisobad04 + ph1_ecalisobad + 2.5 - rho * 0.23) * 50. / ph1_badvtx_Et;
	  float val_trkisooet_1 = ph1_pfchargedisogood03 * 50. / ph1_pt;
	  float val_isosumoet_2 = (ph2_pfchargedisogood03 + ph2_ecaliso + 2.5 - rho * 0.09) * 50. / ph2_pt;
	  float val_isosumoetbad_2 = (ph2_pfchargedisobad04 + ph2_ecalisobad + 2.5 - rho * 0.23) * 50. / ph2_badvtx_Et;
	  float val_trkisooet_2 = ph2_pfchargedisogood03 * 50. / ph2_pt;
	  */

	  /*
	  if( category == 0 && // EBEB HR9 HR9
	      (ph1_isEB == 1 && ph1_r9 > 0.94 && 
	       ph1_sieie < 0.0125 && ph1_hoe < 0.141 && //ph1_ciclevel > 3 &&
	       val_isosumoet_1 < 6.3 && val_isosumoetbad_1 < 18.9 && val_trkisooet_1 < 4.5) &&
	      (ph2_isEB == 1 && ph2_r9 > 0.94 && 
	       ph2_sieie < 0.0125 && ph2_hoe < 0.141 && //ph2_ciclevel > 3 &&
	       val_isosumoet_2 < 6.3 && val_isosumoetbad_2 < 18.9 && val_trkisooet_2 < 4.5) ) diPhotonSuperTIGHT_passe = true;
	  else if( category == 1 && // EBEB !HH
		   (ph1_isEB == 1 && ph1_r9 > 0.33 && 
		    ph1_sieie < 0.0103 && ph1_hoe < 0.138 && //ph1_ciclevel > 3 &&
		    val_isosumoet_1 < 5.6 && val_isosumoetbad_1 < 8. && val_trkisooet_1 < 2.8 &&
		    ph2_isEB == 1 && ph2_r9 > 0.33 && 
		    ph2_sieie < 0.0103 && ph2_hoe < 0.138 && //ph2_ciclevel > 3 &&
		    val_isosumoet_2 < 5.6 && val_isosumoetbad_2 < 8. && val_trkisooet_2 < 2.8) ) diPhotonSuperTIGHT_passe = true;
	  else if( category == 2 && // !EBEB HH
		   (ph1_r9 > 0.94 && 
		    ph1_sieie < 0.029 && ph1_hoe < 0.12 && //ph1_ciclevel > 3 &&
		    val_isosumoet_1 < 5.8 && val_isosumoetbad_1 < 10. && val_trkisooet_1 < 4.  &&
		    ph2_r9 > 0.94 && 
		    ph2_sieie < 0.029 && ph2_hoe < 0.12 && //ph2_ciclevel > 3 &&
		    val_isosumoet_2 < 5.8 && val_isosumoetbad_2 < 10. && val_trkisooet_2 < 4.) ) diPhotonSuperTIGHT_passe = true;
	  else if( category == 3 && // !EBEB !HH
		   (ph1_r9 > 0.37  &&
		    ph1_sieie < 0.028 && ph1_hoe < 0.091 && ph1_ciclevel > 3 && 
		    val_isosumoet_1 < 5.1 &&  val_isosumoetbad_1 < 6.2 && val_trkisooet_1 < 1.62 && 
		    ph2_r9 > 0.37  && 
		    ph2_sieie < 0.028 && ph2_hoe < 0.091 && ph2_ciclevel > 3 && 
		    val_isosumoet_2 < 5.1 &&  val_isosumoetbad_2 < 6.2 && val_trkisooet_2 < 1.62) ) diPhotonSuperTIGHT_passe = true;
	  */
	  //	  if(Photon1SuperTIGHT_passe && Photon2SuperTIGHT_passe) diPhotonSuperTIGHT_passe = true;

         if((TMath::Abs(ph1_SCEta)>1.4442&&TMath::Abs(ph1_SCEta)<1.566)||(TMath::Abs(ph2_SCEta)>1.4442&&TMath::Abs(ph2_SCEta)<1.566)
       || TMath::Abs(ph1_SCEta)>2.5 || TMath::Abs(ph2_SCEta)>2.5) continue;
          
	  ///////////////////////////////////////////////////////////
	  TLorentzVector diPhoton1;
	  diPhoton1.SetPtEtaPhiE(ph1_pt,ph1_eta,ph1_phi,ph1_e);
	  TLorentzVector diPhoton2;
	  diPhoton2.SetPtEtaPhiE(ph2_pt,ph2_eta,ph2_phi,ph2_e);
	  TLorentzVector diPhoton = diPhoton1 + diPhoton2;	  
	  ///////////////////////////////////////////////////////////

	  if(printOk)     std::cout << "--- diPhoton.M() = " << diPhoton.M() << std::endl;

	  if(ph1_pt > ph2_pt && (ph2_pt < 25 || ph1_pt < diPhoton.M()/3.) ) isbad = true;
	  if(ph2_pt > ph1_pt && (ph1_pt < 25 || ph2_pt < diPhoton.M()/3.) ) isbad = true;
          
          //PhotonPt cut (1)
          ++Step_Counter;
          if(isbad == false){
	     ++numPassingEvents[Step_Counter];
             event_selSTEP.push_back(Step_Counter);
	     ++numPassingEvents_CR1[Step_Counter];
             event_selSTEP_CR1.push_back(Step_Counter);
	     ++numPassingEvents_CR2[Step_Counter];
             event_selSTEP_CR2.push_back(Step_Counter);
          }
         
          isbad = false;

          if(printOk)     std::cout << "--- prima di diPhotonSuperTIGHT_passe = " << diPhotonSuperTIGHT_passe << std::endl;
	  if(diPhotonSuperTIGHT_passe == false) isbad = true;

          //CIC Hgg + PhotonPt cut(2)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++)
//                  if(event_selSTEP.at(ii) == Step_Counter-1) isbad = false;
//              if(isbad == false){

	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-1) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
		++numPassingEvents_CR1[Step_Counter];
		event_selSTEP_CR1.push_back(Step_Counter);
		++numPassingEvents_CR2[Step_Counter];
		event_selSTEP_CR2.push_back(Step_Counter);
             }
          }
          
          isbad = false; 
	  float massLowCut = 120.;
	  if(diPhoton.M() < 119.9999 || diPhoton.M() > 130.) isbad = true;
	  if(printOk)     std::cout << "--- diPhoton.M() = " << diPhoton.M() << std::endl;
	  if(printOk)     std::cout << "--- PhotonsMass = " << PhotonsMass << std::endl;

          //CIC Hgg + DiMass cut 100-180 + PhotonPt cut (3)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-1) isbad = false;
//              }
//              if(isbad == false){

	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-1) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
		if(printOk)     std::cout << "--- Passed cut(3) SR " << std::endl;
             }
          }

	  // //CIC Hgg + DiMass cut < 100 CR1
	  if(diPhoton.M() < 119.9999 && event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-1){
		++numPassingEvents_CR1[Step_Counter];
		event_selSTEP_CR1.push_back(Step_Counter);
		if(printOk)     std::cout << "--- Passed cut(3) CR1 " << std::endl;
		if(printOk && diPhoton.M() < massLowCut)     std::cout << "--- Passed cut(3) CR1 && diPhoton.M() < 120." << std::endl;
	    }
	  // //CIC Hgg + DiMass cut > 180 CR2
	  if(isbad == true && diPhoton.M() > 130. && diPhoton.M() < 300. && event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-1){
	    ++numPassingEvents_CR2[Step_Counter];
	    event_selSTEP_CR2.push_back(Step_Counter);
		if(printOk)     std::cout << "--- Passed cut(3) CR2 " << std::endl;
	  }

          isbad = false;  

	  if(printOk) {
	  std::cout << " h_vtx_ind.size() = " << h_vtx_ind.size() << std::endl;
	  std::cout << " Step_Counter = " << Step_Counter << std::endl;
	  }
	  //	  if(printOk)     std::cout << "--- before new variables " << std::endl


	  //save photons
          phoP4[0] = new TLorentzVector();
	  phoP4[0]->SetPtEtaPhiE(ph1_pt,ph1_eta,ph1_phi,ph1_e);
          
          phoP4[1] = new TLorentzVector();
	  phoP4[1]->SetPtEtaPhiE(ph2_pt,ph2_eta,ph2_phi,ph2_e);
      
          for(int ii = 0; ii < 2; ii++){
	    pho_pt->push_back(phoP4[ii]->Pt());
	    phoPt_map[phoP4[ii]->Pt()] = ii;
          }
          
	  std::map<float,int>::reverse_iterator itPhoMap = phoPt_map.rbegin();
	  int itMaxPho,itNextToMaxPho;
	  itMaxPho = itPhoMap->second;
	  ++itPhoMap;
	  itNextToMaxPho = itPhoMap->second;

	  if(printOk){
	  std::cout << " >>> itMaxPho = " << itMaxPho << std::endl;
	  std::cout << " >>> itNextToMaxPho = " << itNextToMaxPho << std::endl;
	  
	  for(std::map<float,int>::iterator it=phoPt_map.begin(); it!=phoPt_map.end(); ++it)
	      std::cout << it->first << " => " << it->second << '\n';
	  }

         //------------------- Jets Selections -------------------//
        
          jet_csvBtag[0] = j1_csvBtag;
          jet_csvMvaBtag[0] = j1_csvMvaBtag;
          jet_jetProbBtag[0] = j1_jetProbBtag;
          jet_tcheBtag[0] = j1_tcheBtag;

          jet_csvBtag[1] = j2_csvBtag;
          jet_csvMvaBtag[1] = j2_csvMvaBtag;
          jet_jetProbBtag[1] = j2_jetProbBtag;
          jet_tcheBtag[1] = j2_tcheBtag;

          jet_csvBtag[2] = j3_csvBtag;
          jet_csvMvaBtag[2] = j3_csvMvaBtag;
          jet_jetProbBtag[2] = j3_jetProbBtag;
          jet_tcheBtag[2] = j3_tcheBtag;

          jet_csvBtag[3] = j4_csvBtag;
          jet_csvMvaBtag[3] = j4_csvMvaBtag;
          jet_jetProbBtag[3] = j4_jetProbBtag;
          jet_tcheBtag[3] = j4_tcheBtag;


          for(int ii = 0; ii < 4; ii++)
              bad_jet_pre[ii] = false;
         
          for (int ij=0; ij<4; ij++){ 
               if ( ptcorrjet[ij]< -1)               bad_jet_pre[ij] = true;
               //if ( btagcsvjet[ij]<=0)               bad_jet_pre[ij] = true;
          }
          for (int ij=0; ij<4; ij++){ 
               if(bad_jet_pre[ij] == true) continue;
               n_good_jets_pre++;
          }
          for(int ii = 0; ii < 4; ii++)
              bad_jet[ii] = false;

          for (int ij=0; ij<4; ij++) { 
               if ( ptcorrjet[ij]<-1)               bad_jet[ij] = true;
               if ( btagcsvjet[ij]<=0)              bad_jet[ij] = true;
               if ( ptcorrjet[ij] < 25. )           bad_jet[ij] = true;
               if ( fabs(etajet[ij])>2.5 )          bad_jet[ij] = true;
               if ( !passCutBasedJetId(ij) )        bad_jet[ij] = true;
          }

          for(int ii = 0; ii < 4; ii++){
              jetP4[ii] = new TLorentzVector();
              jetP4[ii]->SetPtEtaPhiE(ptcorrjet[ii],etajet[ii],phijet[ii],ecorrjet[ii]);
              if(bad_jet[ii] == true) continue;
              n_good_jets++;
          }
           
          for(int ii = 0; ii < 4; ii++)
              if(jet_csvBtag[ii] > btag_wp && bad_jet_pre[ii] == false) n_b_jets_pre++;

          for(int ii = 0; ii < 4; ii++)
              if(jet_csvBtag[ii] > btag_wp && bad_jet[ii] == false) n_b_jets++;

          if(n_good_jets_pre < 2) isbad = true; 
           
          //NJets >= 2 (4)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-1) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-1) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }

          isbad = false;

          if(n_b_jets_pre < 1 || n_good_jets_pre < 2) isbad = true; 

	  //NBJets >= 1 (5)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-1) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-1) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }

          isbad = false;

          if(n_good_jets < 2) isbad = true; 
           
          //NJets >= 2 (6)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-1) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-1) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }

          isbad = true;

          if(n_b_jets >= 1 && n_good_jets >= 2) isbad = false; 

	  //NBJets >= 1 (7)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-1) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-1) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }

          isbad = true;

          if(n_b_jets == 1 && n_good_jets  >= 2) isbad = false; 
          
	  //NBJets == 1 (8)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-1) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-1) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-1){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }

          isbad = true;
            
          if(n_b_jets >= 2 && n_good_jets  >= 2) isbad = false; 

          //NBJets >= 2 (9)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }

          isbad = true;

          float invMass_Mjj = 0;
          
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
         
          if(n_b_jets >= 1 && n_good_jets >= 2){
            j1 = jjPt_map[jjPt.at(jjPt.size()-1)].first;
            j2 = jjPt_map[jjPt.at(jjPt.size()-1)].second;
          }else{
           j1 = 0;
           j2 = 1;
           isfinal_event = false;
          }
         
          jjPt.clear();
          jjPt_map.clear();

          Hgg = new TLorentzVector();
          *Hgg = *phoP4[0]+*phoP4[1];  
          Hjj = new TLorentzVector();
          *Hjj = *jetP4[j1]+*jetP4[j2];
          Radion = new TLorentzVector();
          *Radion = *Hgg + *Hjj;
          Hgg->Boost(-Radion->BoostVector());
          
          isbad = true;

          if(fabs(Hgg->CosTheta()) < 0.9 && n_b_jets == 1 && n_good_jets  >= 2) isbad = false;
	  //remove cut ARA ===>>> will be removed anyway
	  isbad = false;

          //1BJet and |costheta*| < 0.9 (10)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }
           
          isbad = true;

          if(fabs(Hgg->CosTheta()) < 0.9 && n_b_jets >= 2 && n_good_jets  >= 2) isbad = false;
	  //remove cut ARA ===>>> will be removed anyway
	  isbad = false;

          //2BJet and |costheta*| < 0.9 (11)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }
          
          isbad = false;
         
          delete Hgg;
          delete Hjj;
          delete Radion;
             
          sumP4 = new TLorentzVector();

          *sumP4 = *jetP4[j1]+*jetP4[j2];
           
          invMass_Mjj = sumP4->M();

          delete sumP4;

	  bool cloneTree = false;

	  if(n_b_jets == 1 && (invMass_Mjj < 90. || invMass_Mjj > 170.)) isbad = true;
	  //  if(n_b_jets == 1 && (invMass_Mjj < 50. || invMass_Mjj > 250.)) isbad = true;
	  //	  else cloneTree = true;         

          //1BJet and mjj cut (12)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
		//		cloneTree = true;
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	       //	       cloneTree = true;
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	       //	       cloneTree = true;
	     }
          }
          
          isbad = false;

	  if(n_b_jets >= 2 && (invMass_Mjj < 90. || invMass_Mjj > 170.)) isbad = true;
	  //          if(n_b_jets >= 2 && (invMass_Mjj < 50. || invMass_Mjj > 250.)) isbad = true;
	  //	  else cloneTree = true;         

          //2BJet and mjj cut (13)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
		//		cloneTree = true;
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	       //	       cloneTree = true;
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	       //	       cloneTree = true;
	     }
          }
          

          isbad = false;

          TVector3 t3phot1, t3phot2;
          t3phot1.SetPtEtaPhi(phoP4[0]->Pt(),phoP4[0]->Eta(),phoP4[0]->Phi());
          t3phot2.SetPtEtaPhi(phoP4[1]->Pt(),phoP4[1]->Eta(),phoP4[1]->Phi());
          TLorentzVector t4phot1, t4phot2;
          t4phot1.SetPtEtaPhiE(phoP4[0]->Pt(),phoP4[0]->Eta(),phoP4[0]->Phi(),phoP4[0]->Energy());
          t4phot2.SetPtEtaPhiE(phoP4[1]->Pt(),phoP4[1]->Eta(),phoP4[1]->Phi(),phoP4[1]->Energy());
          
          TVector3 t3jet1, t3jet2;
          TLorentzVector t4jet1, t4jet2;
          if(n_b_jets >= 1 && n_good_jets >= 2){
             t3jet1.SetPtEtaPhi(jetP4[j1]->Pt(),jetP4[j1]->Eta(),jetP4[j1]->Phi());
             t3jet2.SetPtEtaPhi(jetP4[j2]->Pt(),jetP4[j2]->Eta(),jetP4[j2]->Phi());
             t4jet1.SetPtEtaPhiE(jetP4[j1]->Pt(),jetP4[j1]->Eta(),jetP4[j1]->Phi(),jetP4[j1]->Energy());
             t4jet2.SetPtEtaPhiE(jetP4[j2]->Pt(),jetP4[j2]->Eta(),jetP4[j2]->Phi(),jetP4[j2]->Energy());
          }
          
          std::pair<TLorentzVector,TLorentzVector> jets_kinfitH;
          if(n_b_jets >= 1 && n_good_jets >= 2) jets_kinfitH = fitter_jetsH->fit(t4jet1, t4jet2);
          TLorentzVector jet1_kinfit, jet2_kinfit;
          jet1_kinfit = jets_kinfitH.first;
          jet2_kinfit = jets_kinfitH.second;

          TLorentzVector t4diPhot = t4phot1+t4phot2;
          TLorentzVector dijet_kinfit = jet1_kinfit + jet2_kinfit;
          TLorentzVector Vstar_kinfit = dijet_kinfit + t4diPhot;
           
          float radMassKF = Vstar_kinfit.M();

          if(n_b_jets == 1 && (radMassKF<200. || radMassKF>400.)) isbad = true;

          //1BJet and mggjj cut (14)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
		//		cloneTree = true; 
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	       //	       cloneTree = true; 
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	       //	       cloneTree = true; 
	     }
          }
          
          isbad = false;

          if(n_b_jets >= 2 && (radMassKF<200. || radMassKF>400.)) isbad = true;

          //2BJet and mggjj cut (15)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
		cloneTree = true; 
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	       cloneTree = true; 
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	       cloneTree = true; 
	     }
          }
	  //if 15 
	  if(cloneTree) newtree->Fill();

          isbad = false;
          
          float deltaR_g1b1;
          float deltaR_g1b2;
          float deltaR_g2b1;
          float deltaR_g2b2;
          if(n_b_jets >= 1 && n_good_jets >= 2){
             deltaR_g1b1  = t4phot1.DeltaR(jet1_kinfit);
             deltaR_g1b2  = t4phot1.DeltaR(jet2_kinfit);
             deltaR_g2b1  = t4phot2.DeltaR(jet1_kinfit);
             deltaR_g2b2  = t4phot2.DeltaR(jet2_kinfit);
          }
          float minDeltaR_gb;
          if(n_b_jets >= 1 && n_good_jets >= 2){
             minDeltaR_gb = deltaR_g1b1;
             if ( deltaR_g1b2 < minDeltaR_gb ) minDeltaR_gb = deltaR_g1b2;
             if ( deltaR_g2b1 < minDeltaR_gb ) minDeltaR_gb = deltaR_g2b1;
             if ( deltaR_g2b2 < minDeltaR_gb ) minDeltaR_gb = deltaR_g2b2;
          }
          if (n_b_jets == 1 && minDeltaR_gb<1) isbad = true;

          //1BJet and deltaR(g,j) cut (16)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }
          
          isbad = false;
            
          if (n_b_jets >= 2 && minDeltaR_gb<1) isbad = true;

          //2BJet and deltaR(g,j) cut (17)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }
          
          isbad = false;

          if (n_b_jets == 1 && njets_passing_kLooseID>=3) isbad = true;

          //1BJet and njets cut (18)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }

          isbad = false;

          if (n_b_jets >= 2 && njets_passing_kLooseID>=3) isbad = true;

          //2BJet and njets cut (19)
          ++Step_Counter;
	  if(isbad == false){
//              isbad = true;
//              for(unsigned int ii = 0; ii < event_selSTEP.size(); ii++){
//                  if(event_selSTEP.at(ii) == Step_Counter-2) isbad = false;
//              }
//              if(isbad == false){
	    if(event_selSTEP.at(event_selSTEP.size()-1) == Step_Counter-2) {
	        ++numPassingEvents[Step_Counter];
                event_selSTEP.push_back(Step_Counter);
             }
	     if(event_selSTEP_CR1.at(event_selSTEP_CR1.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR1[Step_Counter];
	       event_selSTEP_CR1.push_back(Step_Counter);
	     }
	     if(event_selSTEP_CR2.at(event_selSTEP_CR2.size()-1) == Step_Counter-2){
	       ++numPassingEvents_CR2[Step_Counter];
	       event_selSTEP_CR2.push_back(Step_Counter);
	     }
          }
      
          isbad = false;


	  /////////////////////////////////////////////////////////////////////////////////
	  for(unsigned int nn = 0; nn < event_selSTEP_CR1.size(); nn++){
	    Step_Counter = event_selSTEP_CR1.at(nn);

	    if(Step_Counter > 9)
	      if(isfinal_event == false) continue; 
           
	    if(Step_Counter == 4 && PhotonsMass >= 120.) {
	      std::cout << " ientry = " << ientry << std::endl;
	      std::cout << " inputTree = " << inputTree << std::endl;
	    }
	    h_diPhoInvMass_CR1.at(Step_Counter)->Fill(PhotonsMass,puRe);


	    for(int ii = 0; ii < 4; ii++){
	      if(Step_Counter <= 3)
		if(bad_jet_pre[ii] == true) continue;
              if(Step_Counter > 3)
		if(bad_jet[ii] == true) continue;
// 	      if(printOk){
// 		std::cout << " >>>> jetP4[ii] Pt " << jetP4[ii]->Pt() << std::endl;
// 		std::cout << " >>>> jetP4[ii] Eta " << jetP4[ii]->Eta() << std::endl;
// 		std::cout << " >>>> jetP4[ii] Phi " << jetP4[ii]->Phi() << std::endl;
// 		std::cout << " >>>> jetP4[ii] E " << jetP4[ii]->E() << std::endl;
// 	      }
	      jetPt_map[jetP4[ii]->Pt()] = ii;
	    }

	    int itMaxJet_CR1,itNextToMaxJet_CR1;
	    if(jetPt_map.size() > 1.){
	      std::map<float,int>::reverse_iterator itJetMap = jetPt_map.rbegin();
	      itMaxJet_CR1 = itJetMap->second;
	      ++itJetMap;
	      itNextToMaxJet_CR1 = itJetMap->second;

	      sumP4 = new TLorentzVector();
	      //	      *sumP4 = *jetP4[itMaxJet_CR1] + *jetP4[itNextToMaxJet_CR1] + *PhoSumP4;
	      *sumP4 = *jetP4[itMaxJet_CR1] + *jetP4[itNextToMaxJet_CR1] + *phoP4[0]+ *phoP4[1];
	      h_diPhodiJetInvMass_maxpt_CR1.at(Step_Counter)->Fill(sumP4->M(),puRe);
	      delete sumP4;
	    }
	  }
	  event_selSTEP_CR1.clear();
	  jetPt_map.clear();

	  /* 
	  for(unsigned int nn = 0; nn < event_selSTEP_CR2.size(); nn++){
	    Step_Counter = event_selSTEP_CR2.at(nn);

	    if(Step_Counter > 9)
	      if(isfinal_event == false) continue; 
           
	    h_diPhoInvMass_CR2.at(Step_Counter)->Fill(PhotonsMass,puRe);

	    for(int ii = 0; ii < 4; ii++){
	      if(Step_Counter <= 3)
		if(bad_jet_pre[ii] == true) continue;
              if(Step_Counter > 3)
		if(bad_jet[ii] == true) continue;
// 	      if(printOk){
// 		std::cout << " >>>> jetP4[ii] Pt " << jetP4[ii]->Pt() << std::endl;
// 		std::cout << " >>>> jetP4[ii] Eta " << jetP4[ii]->Eta() << std::endl;
// 		std::cout << " >>>> jetP4[ii] Phi " << jetP4[ii]->Phi() << std::endl;
// 		std::cout << " >>>> jetP4[ii] E " << jetP4[ii]->E() << std::endl;
// 	      }
	      jetPt_map[jetP4[ii]->Pt()] = ii;
	    }

	    int itMaxJet_CR2,itNextToMaxJet_CR2;
	    if(jetPt_map.size() > 1.){
	      std::map<float,int>::reverse_iterator itJetMap = jetPt_map.rbegin();
	      itMaxJet_CR2 = itJetMap->second;
	      ++itJetMap;
	      itNextToMaxJet_CR2 = itJetMap->second;

	      sumP4 = new TLorentzVector();
	      *sumP4 = *jetP4[itMaxJet_CR2] + *jetP4[itNextToMaxJet_CR2] + *PhoSumP4;
	      h_diPhodiJetInvMass_maxpt_CR2.at(Step_Counter)->Fill(sumP4->M(),puRe);
	      delete sumP4;
	    }
	  }
	  jetPt_map.clear();
	  event_selSTEP_CR2.clear();
	  */
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          for(unsigned int nn = 0; nn < event_selSTEP_CR2.size(); nn++){
	    Step_Counter = event_selSTEP_CR2.at(nn);

          if(Step_Counter > 9)
             if(isfinal_event == false) continue; 
           
	  h_diPhoInvMass_CR2.at(Step_Counter)->Fill(PhotonsMass,puRe);

          for(int ii = 0; ii < 4; ii++){
	    if(Step_Counter <= 3)
	      if(bad_jet_pre[ii] == true) continue;
	    if(Step_Counter > 3)
	      if(bad_jet[ii] == true) continue;
	    jetPt_map[jetP4[ii]->Pt()] = ii;
          }

	  //	  if(printOk)	  std::cout << " >>>> CR2 STEP == " << Step_Counter << " jetPt_map.size() = " << jetPt_map.size() << std::endl;
	  //	  if(Step_Counter == 0)	  std::cout << " >>>> CR2 STEP == " << Step_Counter << " jetPt_map.size() = " << jetPt_map.size() << std::endl;

	  int itMaxJet,itNextToMaxJet;
	  if(jetPt_map.size() > 1.){
	  std::map<float,int>::reverse_iterator itJetMap = jetPt_map.rbegin();
	  itMaxJet = itJetMap->second;
	  ++itJetMap;
	  itNextToMaxJet = itJetMap->second;
 
	  sumP4 = new TLorentzVector();
	  //          *sumP4 = *jetP4[itMaxJet] + *jetP4[itNextToMaxJet] + *PhoSumP4;
	  *sumP4 = *jetP4[itMaxJet] + *jetP4[itNextToMaxJet] + *phoP4[0]+ *phoP4[1];
	  h_diPhodiJetInvMass_maxpt_CR2.at(Step_Counter)->Fill(sumP4->M(),puRe);
// 	  if(Step_Counter == 0)	  std::cout << " >>>> PhoSumP4->M() = " << (*phoP4[0]+ *phoP4[1]).M() << " >>>> sumP4->M() = " << sumP4->M() << std::endl;
// 	  if(Step_Counter == 0)	  std::cout << " >>>> itMaxJet = " << itMaxJet << " >>>> itNextToMaxJet = " << itNextToMaxJet << std::endl;
          delete sumP4;
          }
          }

          jetPt_map.clear();
          event_selSTEP_CR2.clear();

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          for(unsigned int nn = 0; nn < event_selSTEP.size(); nn++){
               
          Step_Counter = event_selSTEP.at(nn);

          if(Step_Counter > 9)
             if(isfinal_event == false) continue; 
           
          h_vtx_ind.at(Step_Counter)->Fill(vtx_ind, puRe);

          h_JetNum.at(Step_Counter)->Fill(n_good_jets, puRe);
          h_JetNum_btagged.at(Step_Counter)->Fill(n_b_jets, puRe);

	  h_PhoEcalIso.at(Step_Counter)->Fill(ph1_ecaliso/ph1_pt, puRe);
	  h_PhoEcalIso.at(Step_Counter)->Fill(ph2_ecaliso/ph2_pt, puRe);

	  h_PhoEcalIsoBad.at(Step_Counter)->Fill(ph1_pfchargedisogood03 * 50. / ph1_pt, puRe);
	  h_PhoEcalIsoBad.at(Step_Counter)->Fill(ph2_pfchargedisogood03 * 50. / ph2_pt, puRe);

	  h_PhoChargedIso.at(Step_Counter)->Fill( (ph1_pfchargedisogood03+ph1_ecaliso+2.5 - rho * 0.09)*50./ph1_pt, puRe);
	  h_PhoChargedIso.at(Step_Counter)->Fill( (ph2_pfchargedisogood03+ph2_ecaliso+2.5 - rho * 0.09)*50./ph2_pt, puRe);

	  h_PhoChargedIsoBad.at(Step_Counter)->Fill( (ph1_pfchargedisobad04+ph1_ecalisobad+2.5 - rho * 0.23) * 50. / ph1_badvtx_Et, puRe);
	  h_PhoChargedIsoBad.at(Step_Counter)->Fill( (ph2_pfchargedisobad04+ph2_ecalisobad+2.5 - rho * 0.23) * 50. / ph2_badvtx_Et, puRe);

	  h_ciclevel.at(Step_Counter)->Fill(ph1_ciclevel, puRe);
	  h_ciclevel.at(Step_Counter)->Fill(ph2_ciclevel, puRe);

	  h_PhoR9.at(Step_Counter)->Fill(ph1_r9,puRe);
          h_PhoR9.at(Step_Counter)->Fill(ph2_r9,puRe);

          h_PhoHoE.at(Step_Counter)->Fill(ph1_hoe,puRe);
          h_PhoHoE.at(Step_Counter)->Fill(ph2_hoe,puRe);
           
          h_PhoSieie.at(Step_Counter)->Fill(ph1_sieie,puRe);
          h_PhoSieie.at(Step_Counter)->Fill(ph2_sieie,puRe);

          h_PhoSipip.at(Step_Counter)->Fill(ph1_sipip,puRe);
          h_PhoSipip.at(Step_Counter)->Fill(ph2_sipip,puRe);

          h_PhoSieip.at(Step_Counter)->Fill(ph1_sieip,puRe);
          h_PhoSieip.at(Step_Counter)->Fill(ph1_sieip,puRe);
         
	  for(int ii = 0; ii < 2; ii++){
	    h_PhoPt.at(Step_Counter)->Fill(phoP4[ii]->Pt(),puRe); 
	    h_PhoEta.at(Step_Counter)->Fill(phoP4[ii]->Eta(),puRe);
	    h_PhoPhi.at(Step_Counter)->Fill(phoP4[ii]->Phi(),puRe);
	  }


          h_diPhoDeltaEta.at(Step_Counter)->Fill(ph1_eta-ph2_eta,puRe);
	  h_diPhoDeltaPhi.at(Step_Counter)->Fill(phoP4[0]->DeltaPhi(*phoP4[1]),puRe);
          h_diPhoDeltaR.at(Step_Counter)->Fill(phoP4[0]->DeltaR(*phoP4[1]),puRe); 
	  h_diPhoInvMass.at(Step_Counter)->Fill(PhotonsMass,puRe);

	  if(printOk)
	  std::cout << " >>>> prima di sumP4 " << std::endl;

          sumP4 = new TLorentzVector();

          *sumP4 = *phoP4[0]+*phoP4[1];

	  if(printOk)
	  std::cout << " >>>> dopo di sumP4 " << std::endl;

          PhoSumP4 = (TLorentzVector*)sumP4->Clone();
          delete sumP4;

	  if(printOk){
	  std::cout << " >>>> PhoSumP4 Pt " << PhoSumP4->Pt() << std::endl;
	  std::cout << " >>>> PhoSumP4 Eta " << PhoSumP4->Eta() << std::endl;
	  std::cout << " >>>> PhoSumP4 Phi " << PhoSumP4->Phi() << std::endl;
	  std::cout << " >>>> PhoSumP4 E " << PhoSumP4->E() << std::endl;
	  }

          for(int ii = 0; ii < 4; ii++){
	      if(Step_Counter <= 3)
                 if(bad_jet_pre[ii] == true) continue;
              if(Step_Counter > 3)
                 if(bad_jet[ii] == true) continue;
	    if(printOk){
	  std::cout << " >>>> jetP4[ii] Pt " << jetP4[ii]->Pt() << std::endl;
	  std::cout << " >>>> jetP4[ii] Eta " << jetP4[ii]->Eta() << std::endl;
	  std::cout << " >>>> jetP4[ii] Phi " << jetP4[ii]->Phi() << std::endl;
	  std::cout << " >>>> jetP4[ii] E " << jetP4[ii]->E() << std::endl;
	    }
           
              h_JetPt.at(Step_Counter)->Fill(jetP4[ii]->Pt(),puRe); 
              h_JetEta.at(Step_Counter)->Fill(jetP4[ii]->Eta(),puRe);
              h_JetPhi.at(Step_Counter)->Fill(jetP4[ii]->Phi(),puRe);
              
              if(jet_csvBtag[ii] > 0) h_JetcsvBtag.at(Step_Counter)->Fill(jet_csvBtag[ii],puRe);
              if(jet_csvMvaBtag[ii] > 0) h_JetcsvMvaBtag.at(Step_Counter)->Fill(jet_csvMvaBtag[ii],puRe);
              if(jet_jetProbBtag[ii] > 0) h_JetProbBtag.at(Step_Counter)->Fill(jet_jetProbBtag[ii],puRe);
              if(jet_tcheBtag[ii] > 0) h_JettcheBtag.at(Step_Counter)->Fill(jet_tcheBtag[ii],puRe);

              jet_pt->push_back(jetP4[ii]->Pt());
              jetPt_map[jetP4[ii]->Pt()] = ii;
          }

          jet_couple_num = 0;

          if(Step_Counter <= 9){
          for(int ii = 0; ii < 4; ii++){

              if(Step_Counter <= 3)
                 if(bad_jet_pre[ii] == true) continue;
              if(Step_Counter > 3)
                 if(bad_jet[ii] == true) continue;

              for(int jj = 0; jj < 4; jj++){

                  if(Step_Counter <= 3)
                     if(bad_jet_pre[jj] == true) continue;
                  if(Step_Counter > 3)
                     if(bad_jet[jj] == true) continue;    

                  if(ii == jj) continue;
                  
                  h_diJetDeltaEta.at(Step_Counter)->Fill(jetP4[ii]->Eta()-jetP4[jj]->Eta(), puRe);
                  h_diJetDeltaPhi.at(Step_Counter)->Fill(jetP4[ii]->DeltaPhi(*jetP4[jj]),puRe);
                  h_diJetDeltaR.at(Step_Counter)->Fill(jetP4[ii]->DeltaR(*jetP4[jj]),puRe); 
                  
                  sumP4 = new TLorentzVector();
                  *sumP4 = *jetP4[ii]+*jetP4[jj];
                  
                  JetSumP4[jet_couple_num] = (TLorentzVector*)sumP4->Clone();
                  h_diJetInvMass.at(Step_Counter)->Fill(sumP4->M(),puRe);

                  isBTagCouple[jet_couple_num] = false;
                  isBVetoCouple[jet_couple_num] = false;

                  if(jet_csvBtag[ii] > btag_wp) isBTagCouple[jet_couple_num] = true;
                  if(jet_csvBtag[jj] > btag_wp) isBTagCouple[jet_couple_num] = true;
                  if(jet_csvBtag[ii] <= 1-btag_wp && jet_csvBtag[jj] <= btag_wp) isBVetoCouple[jet_couple_num] = true;
                  if(jet_csvBtag[ii] <= btag_wp && jet_csvBtag[jj] <= 1-btag_wp) isBVetoCouple[jet_couple_num] = true;
                   
                  delete sumP4;
              
                 jet_couple_num++;
              }
          }
          }else{
              h_diJetDeltaEta.at(Step_Counter)->Fill(jetP4[j1]->Eta()-jetP4[j2]->Eta(), puRe);
              h_diJetDeltaPhi.at(Step_Counter)->Fill(jetP4[j1]->DeltaPhi(*jetP4[j2]),puRe);
              h_diJetDeltaR.at(Step_Counter)->Fill(jetP4[j1]->DeltaR(*jetP4[j2]),puRe); 
                  
              sumP4 = new TLorentzVector();
              *sumP4 = *jetP4[j1]+*jetP4[j2];
                 
              h_diJetInvMass.at(Step_Counter)->Fill(sumP4->M(),puRe);

              delete sumP4;
          }


	  if(printOk)	  std::cout << " >>>> jetPt_map.size() = " << jetPt_map.size() << std::endl;
	  //	  if(Step_Counter == 0)	  std::cout << " >>>> SR STEP == " << Step_Counter << " jetPt_map.size() = " << jetPt_map.size() << std::endl;

	  int itMaxJet,itNextToMaxJet;
	  if(jetPt_map.size() > 1.){
	  std::map<float,int>::reverse_iterator itJetMap = jetPt_map.rbegin();
	  itMaxJet = itJetMap->second;
	  ++itJetMap;
	  itNextToMaxJet = itJetMap->second;
 
	  if(printOk){
	  std::cout << " >>> itMaxJet = " << itMaxJet << std::endl;
	  std::cout << " >>> itNextToMaxJet = " << itNextToMaxJet << std::endl;

	  for(std::map<float,int>::iterator it=jetPt_map.begin(); it!=jetPt_map.end(); ++it)
	      std::cout << it->first << " => " << it->second << '\n';
	  }

          h_diJetDeltaEta_maxpt.at(Step_Counter)->Fill(jetP4[itMaxJet]->Eta()-jetP4[itNextToMaxJet]->Eta(), puRe);
          h_diJetDeltaPhi_maxpt.at(Step_Counter)->Fill(jetP4[itMaxJet]->DeltaPhi(*jetP4[itNextToMaxJet]), puRe);
          h_diJetDeltaR_maxpt.at(Step_Counter)->Fill(jetP4[itMaxJet]->DeltaR(*jetP4[itNextToMaxJet]), puRe); 

          sumP4 = new TLorentzVector();

          *sumP4 = *jetP4[itMaxJet] + *jetP4[itNextToMaxJet];
                  
          h_diJetInvMass_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);

          delete sumP4;

          for(int ii = 0; ii < 2; ii++)
              for(int jj = 0; jj < 4; jj++){
                   
                  if(Step_Counter <= 3)
                     if(bad_jet_pre[jj] == true) continue;
                  if(Step_Counter > 3)
                     if(bad_jet[jj] == true) continue;

                  h_PhoJetDeltaEta.at(Step_Counter)->Fill(phoP4[ii]->Eta()-jetP4[jj]->Eta(),puRe);
                  h_PhoJetDeltaPhi.at(Step_Counter)->Fill(phoP4[ii]->DeltaPhi(*jetP4[jj]),puRe);
                  h_PhoJetDeltaR.at(Step_Counter)->Fill(phoP4[ii]->DeltaR(*jetP4[jj]),puRe);  
              }
          
          h_PhoJetDeltaEta_maxpt.at(Step_Counter)->Fill(phoP4[itMaxPho]->Eta()-jetP4[itMaxJet]->Eta(),puRe);
          h_PhoJetDeltaPhi_maxpt.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaPhi(*jetP4[itMaxJet]),puRe);
          h_PhoJetDeltaR_maxpt.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaR(*jetP4[itMaxJet]),puRe);
           
          for(int ii = 0; ii < 4; ii++){
               
	      if(Step_Counter <= 3)
                 if(bad_jet_pre[ii] == true) continue;
              if(Step_Counter > 3)
                 if(bad_jet[ii] == true) continue;
	      if(jet_csvBtag[ii] < btag_wp) continue;

              h_PhoJetDeltaEta_maxpt_btag.at(Step_Counter)->Fill(phoP4[itMaxPho]->Eta() - jetP4[ii]->Eta(),puRe);
              h_PhoJetDeltaPhi_maxpt_btag.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaPhi(*jetP4[ii]),puRe);
              h_PhoJetDeltaR_maxpt_btag.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaR(*jetP4[ii]),puRe);

          }
 
          if(Step_Counter <= 9){
          for(unsigned int jj = 0; jj < JetSumP4.size(); jj++){
      
                  h_diPhodiJetDeltaEta.at(Step_Counter)->Fill(PhoSumP4->Eta()-JetSumP4[jj]->Eta(),puRe);
                  h_diPhodiJetDeltaPhi.at(Step_Counter)->Fill(PhoSumP4->DeltaPhi(*JetSumP4[jj]),puRe);
                  h_diPhodiJetDeltaR.at(Step_Counter)->Fill(PhoSumP4->DeltaR(*JetSumP4[jj]),puRe); 

                  sumP4 = new TLorentzVector();
                   
                  *sumP4 = *PhoSumP4+*JetSumP4[jj];

                  h_diPhodiJetInvMass.at(Step_Counter)->Fill(sumP4->M(),puRe);

                  if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core.at(Step_Counter)->Fill(sumP4->M(),puRe);
                  if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands.at(Step_Counter)->Fill(sumP4->M(),puRe);

                  if(isBTagCouple[jj] == true){
                     h_diPhodiJetInvMass_btag.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core_btag.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands_btag.at(Step_Counter)->Fill(sumP4->M(),puRe);
                  }
                  
                  if(isBVetoCouple[jj] == true){
                     h_diPhodiJetInvMass_bveto.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core_bveto.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands_bveto.at(Step_Counter)->Fill(sumP4->M(),puRe);
                  }
                  
                  delete sumP4;
          }
          }else{
              sumP4 = new TLorentzVector();
              *sumP4 = *jetP4[j1] + *jetP4[j2];
               h_diPhodiJetDeltaEta.at(Step_Counter)->Fill(PhoSumP4->Eta()-sumP4->Eta(),puRe);
               h_diPhodiJetDeltaPhi.at(Step_Counter)->Fill(PhoSumP4->DeltaPhi(*sumP4),puRe);
               h_diPhodiJetDeltaR.at(Step_Counter)->Fill(PhoSumP4->DeltaR(*sumP4),puRe); 
               delete sumP4;
               sumP4 = new TLorentzVector();
              *sumP4 = *jetP4[j1] + *jetP4[j2] + *PhoSumP4;
               h_diPhodiJetInvMass.at(Step_Counter)->Fill(sumP4->M(),puRe);
               if(Step_Counter >= 14) h_diPhodiJetInvMass_kinfit.at(Step_Counter)->Fill(radMassKF,puRe);
               delete sumP4;
          }

          sumP4 = new TLorentzVector();
           
          *sumP4 = *jetP4[itMaxJet] + *jetP4[itNextToMaxJet] + *PhoSumP4;

          h_diPhodiJetInvMass_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);
// 	  if(Step_Counter == 0)   std::cout << " >>>> PhoSumP4->M() = " << PhoSumP4->M() << " >>>> sumP4->M() = " << sumP4->M() << std::endl;
//           if(Step_Counter == 0)   std::cout << " >>>> itMaxJet = " << itMaxJet << " >>>> itNextToMaxJet = " << itNextToMaxJet << std::endl;

          if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);
          if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);
	  
          delete sumP4;
          }
           
          }

	  /////////////////////////////////////////////////////////////////////////////////

          for(int ii = 0; ii < 2; ii++)
              delete phoP4[ii];

          for(int ii = 0; ii < 4; ii++)
              delete jetP4[ii];

          for(int ii = 0; ii < 4; ii++)
              bad_jet[ii] = false;
          
          for(int ii = 0; ii < 4; ii++)
              bad_jet_pre[ii] = false;

          delete pho_pt; 
          delete jet_pt; 
    
	  phoPt_map.clear();
	  jetPt_map.clear();    
          event_selSTEP.clear(); 
      }
  }

  std::ofstream myfile(("EventCounts/Events_"+inputTree+".txt").c_str(), std::ios::out); 

  for(int nLocStep=0; nLocStep<nSTEP; ++ nLocStep){
    std::cout << " >>> numPassingEvents[" << nLocStep << "] = " << numPassingEvents[nLocStep] << std::endl;
    std::cout << " >>> numPassingEvents_CR1[" << nLocStep << "] = " << numPassingEvents_CR1[nLocStep] << std::endl;
    std::cout << " >>> numPassingEvents_CR2[" << nLocStep << "] = " << numPassingEvents_CR2[nLocStep] << std::endl;

    myfile << " >>> numPassingEvents[" << nLocStep << "] = " << numPassingEvents[nLocStep] << std::endl;
    myfile << " >>> numPassingEvents_CR1[" << nLocStep << "] = " << numPassingEvents_CR1[nLocStep] << std::endl;
    myfile << " >>> numPassingEvents_CR2[" << nLocStep << "] = " << numPassingEvents_CR2[nLocStep] << std::endl;
  }
  myfile.close();


  TFile* output = new TFile(outputName.c_str(),"RECREATE");
  output->cd();

  newtree->SetName("EeventTree");
  newtree->Print();
  newtree->AutoSave();
  newtree->Write();

  //  return 100;
  //  continue;
  for(int nLocStep=0; nLocStep<nSTEP; ++ nLocStep){

  h_vtx_ind.at(nLocStep)->Write();

  h_PhoEcalIso.at(nLocStep)->Write();
  h_PhoEcalIsoBad.at(nLocStep)->Write();
  h_PhoChargedIso.at(nLocStep)->Write();
  h_PhoChargedIsoBad.at(nLocStep)->Write();
  h_ciclevel.at(nLocStep)->Write();

  h_PhoR9.at(nLocStep)->Write();
  h_PhoHoE.at(nLocStep)->Write();
  h_PhoSieie.at(nLocStep)->Write();
  h_PhoSipip.at(nLocStep)->Write();
  h_PhoSieip.at(nLocStep)->Write();
  h_PhoPt.at(nLocStep)->Write();
  h_PhoEta.at(nLocStep)->Write();
  h_PhoPhi.at(nLocStep)->Write();
  h_diPhoDeltaEta.at(nLocStep)->Write();
  h_diPhoDeltaPhi.at(nLocStep)->Write();
  h_diPhoDeltaR.at(nLocStep)->Write();
  h_diPhoInvMass.at(nLocStep)->Write();

  h_diPhoInvMass_CR1.at(nLocStep)->Write();
  h_diPhoInvMass_CR2.at(nLocStep)->Write();

  h_JetNum.at(nLocStep)->Write();
  h_JetNum_btagged.at(nLocStep)->Write();
  h_JetcsvBtag.at(nLocStep)->Write();
  h_JetcsvMvaBtag.at(nLocStep)->Write();
  h_JetProbBtag.at(nLocStep)->Write();
  h_JettcheBtag.at(nLocStep)->Write();
  h_JetPt.at(nLocStep)->Write();
  h_JetEta.at(nLocStep)->Write();
  h_JetPhi.at(nLocStep)->Write();
  h_diJetDeltaEta.at(nLocStep)->Write();
  h_diJetDeltaPhi.at(nLocStep)->Write();
  h_diJetDeltaR.at(nLocStep)->Write();
  h_diJetInvMass.at(nLocStep)->Write();
  h_diJetDeltaEta_maxpt.at(nLocStep)->Write();
  h_diJetDeltaPhi_maxpt.at(nLocStep)->Write();
  h_diJetDeltaR_maxpt.at(nLocStep)->Write();
  h_diJetInvMass_maxpt.at(nLocStep)->Write();
  h_PhoJetDeltaEta.at(nLocStep)->Write();
  h_PhoJetDeltaPhi.at(nLocStep)->Write();
  h_PhoJetDeltaR.at(nLocStep)->Write();
  h_PhoJetDeltaEta_maxpt.at(nLocStep)->Write();
  h_PhoJetDeltaPhi_maxpt.at(nLocStep)->Write();
  h_PhoJetDeltaR_maxpt.at(nLocStep)->Write(); 
  h_PhoJetDeltaEta_maxpt_btag.at(nLocStep)->Write();
  h_PhoJetDeltaPhi_maxpt_btag.at(nLocStep)->Write();
  h_PhoJetDeltaR_maxpt_btag.at(nLocStep)->Write();
  h_diPhodiJetDeltaEta.at(nLocStep)->Write();
  h_diPhodiJetDeltaPhi.at(nLocStep)->Write();
  h_diPhodiJetDeltaR.at(nLocStep)->Write();
  h_diPhodiJetInvMass.at(nLocStep)->Write();
  h_diPhodiJetInvMass_core.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands.at(nLocStep)->Write();
  h_diPhodiJetInvMass_maxpt.at(nLocStep)->Write();
  h_diPhodiJetInvMass_maxpt_CR1.at(nLocStep)->Write();
  h_diPhodiJetInvMass_maxpt_CR2.at(nLocStep)->Write();

  h_diPhodiJetInvMass_core_maxpt.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands_maxpt.at(nLocStep)->Write();
  h_diPhodiJetInvMass_btag.at(nLocStep)->Write();
  h_diPhodiJetInvMass_core_btag.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands_btag.at(nLocStep)->Write();
  h_diPhodiJetInvMass_bveto.at(nLocStep)->Write();
  h_diPhodiJetInvMass_core_bveto.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands_bveto.at(nLocStep)->Write();
  h_diPhodiJetInvMass_kinfit.at(nLocStep)->Write();
  }
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


















