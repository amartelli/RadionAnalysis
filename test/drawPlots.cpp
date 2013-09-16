
#include "setTDRStyle.h"
#include "PUreweightingUtils.h"
#include "drawPlotsUtils.h"
#include "ConfigParser.h"
#include "ParserUtils.h"

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
#include "THStack.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

int main(int argc, char** argv)
{
  // Set style options
  setTDRStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  //gStyle->SetOptStat(1110); 
  gStyle->SetOptStat(0000); 
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
  
  std::string inputListData = gConfigParser -> readStringOption("Input::inputListData");
  std::string inputListMC = gConfigParser -> readStringOption("Input::inputListMC");
  
  std::string inputHistos = gConfigParser -> readStringOption("Input::inputHistos");

  float inputDataLumi = gConfigParser -> readFloatOption("Input::inputDataLumi");

  std::string inputLabelData = gConfigParser -> readStringOption("Input::inputLabelData");

  std::string inputCompType = gConfigParser -> readStringOption("Input::inputCompType");

  bool isRadionComp = gConfigParser -> readBoolOption("Input::isRadionComp");
  
  std::string outputDir = gConfigParser -> readStringOption("Output::outputDir");
  
  std::map<int, std::string> HistoName;
  std::map<int, TFile*> FileMC;
  std::map<int, TFile*> FileDA;
  std::map<int, std::map<int, TH1F*> > HistoMC;
  std::map<int, std::map<int, TH1F*> > HistoDA;
  std::map<int, TH1F*> HistoDA_total;
  std::map<int, std::map<int, TH1F*> > HistoMC_total;
  
  int histo_pos = 0;
  int fileMC_pos = 0;
  int fileDA_pos = 0;
  int histo_tot = 0;
  int fileMC_tot = 0;
  int fileDA_tot = 0;
  float scales;
  char label[500];
  int color; 
  int add;
  
  std::vector<float> scale;
  std::map<int, std::string> label_MC;
  std::map<int,int> color_MC; 
  std::map<int,int> addMC;
  std::map<int,std::vector<int> > antiAddMC;
 
  char histo[500];
  char file[500];
  FILE *f_histo;
  FILE *f_file_DA;
  FILE *f_file_MC;

  f_histo = fopen(inputHistos.c_str(),"r");
  f_file_DA = fopen(inputListData.c_str(),"r");
  f_file_MC = fopen(inputListMC.c_str(),"r");
  
  while(fscanf(f_histo,"%s \n", histo) !=EOF ){
    std::string HISTO = std::string(histo);

    if(HISTO.find("#") != std::string::npos) continue;
    std::cout << "\nReading Histo: " << histo << std::endl;
    HistoName[histo_pos] = histo;
    histo_pos++;
  }
  
  histo_tot = histo_pos;
  
  while(fscanf(f_file_DA,"%s \n", file) !=EOF ){
    std::string FILE = std::string(file);

    if(FILE.find("#") != std::string::npos) continue;
    std::cout << "\nReading input_DATA: " << file << std::endl;
    
    FileDA[fileDA_pos] = TFile::Open(file);
    
    fileDA_pos++;
  }
  
  fileDA_tot = fileDA_pos;

  while(fscanf(f_file_MC,"%s %f %s %d %d \n", file, &scales , label, &color, &add) !=EOF ){
    std::string FILE = std::string(file);

    if(FILE.find("#") != std::string::npos) continue;
    std::cout << "\nReading input_MC: " << file << " - " << scales << " - " << label << " - " << color << " - " << add <<std::endl;

    FileMC[fileMC_pos] = TFile::Open(file);
    scale.push_back(scales);
    label_MC[fileMC_pos] = std::string(label);
    color_MC[fileMC_pos] = color;
    addMC[fileMC_pos] = add;
    antiAddMC[add].push_back(fileMC_pos);
     
    fileMC_pos++;
  
  }
 
  fileMC_tot = fileMC_pos;
  
  for(int ii = 0; ii < fileMC_tot; ii++)
      for(int jj = 0; jj < histo_tot; jj++){
          HistoMC[ii][jj] = (TH1F*)FileMC[ii]->Get((HistoName[jj]).c_str());
          std::string histo_name = HistoMC[ii][jj]->GetName();
	  //          if(histo_name.find("Num") == std::string::npos && histo_name.find("Vtx") == std::string::npos && isRadionComp == false) HistoMC[ii][jj]->Rebin(5);
          HistoMC[ii][jj]->SetName((HistoName[jj]+"_total").c_str());
      }
 
  for(int ii = 0; ii < fileDA_tot; ii++)
      for(int jj = 0; jj < histo_tot; jj++){
          HistoDA[ii][jj] = (TH1F*)FileDA[ii]->Get((HistoName[jj]).c_str());
          std::string histo_name = HistoDA[ii][jj]->GetName();
	  //          if(histo_name.find("Num") == std::string::npos && histo_name.find("Vtx") == std::string::npos && isRadionComp == false) HistoDA[ii][jj]->Rebin(5);
          HistoDA[ii][jj]->SetName((HistoName[jj]+"_total").c_str());
      }
 
  for(int ii = 0; ii < histo_tot; ii++){
      HistoDA_total[ii] = (TH1F*)HistoDA[0][ii]->Clone();
      for(int jj = 1; jj < fileDA_tot; jj++)
          HistoDA_total[ii]->Add(HistoDA[jj][ii]);
  } 
 
  for(int ii = 0; ii < fileMC_tot; ii++)
      for(int jj = 0; jj < histo_tot; jj++)
	//          HistoMC[ii][jj]->Scale(5.537/19.62);
	HistoMC[ii][jj]->Scale(inputDataLumi*scale.at(ii));

  
  std::map<int,std::string> label_Final_MC;
  std::map<int,int> color_Final_MC;

  for(int jj = 0; jj < histo_tot; jj++)
      for(int ii = 0; ii < fileMC_tot; ii++)
          HistoMC_total[jj][ii] = (TH1F*)HistoMC[ii][jj]->Clone();
   
  for(int ii = 0; ii < histo_tot; ii++){
      HistoName[ii].erase(0,2);
  
      //      HistoDA_total[ii]->Scale(5.537/19.62);
   
      if(isRadionComp == false) compareHStacks(HistoDA_total[ii], HistoMC_total[ii] , HistoName[ii] ,inputLabelData, label_MC, color_MC, inputCompType , outputDir);
      if(isRadionComp == true) compareHStacks_vs_Radion(HistoDA_total[ii], HistoMC_total[ii] , HistoName[ii] ,inputLabelData, label_MC, color_MC, inputCompType , outputDir);
  }
  
}


















