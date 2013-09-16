#ifndef drawPlotsUtils_h
#define drawPlotsUtils_h

#include "setTDRStyle.h"

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
#include "TPaveLabel.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

void compareHistos(TH1F* , TH1F* , std::string , std::string, std::string , std::string);

void compareGraph(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, float x_min, float x_max, float y_min, float y_max , std::string xTitle, std::string yTitle, std::string label_1, std::string label_2, std::string outputName);

void compare4Graph(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TGraphAsymmErrors* g3, TGraphAsymmErrors* g4,float x_min, float x_max, float y_min, float y_max , std::string xTitle, std::string yTitle, std::string label_1, std::string label_2 , std::string label_3, std::string label_4, std::string outputName);

void compare6Graph(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TGraphAsymmErrors* g3, TGraphAsymmErrors* g4, TGraphAsymmErrors* g5, TGraphAsymmErrors* g6,float x_min, float x_max, float y_min, float y_max , std::string xTitle, std::string yTitle, std::string label_1, std::string label_2 , std::string label_3, std::string label_4, std::string label_5, std::string label_6,std::string outputName);

void compareHStacks(TH1F* histo_DA, std::map<int,TH1F*> histo_MC, std::string xTitle, std::string label_DA, std::map<int,std::string> label_MC, std::map<int,int> color_MC,  std::string comp_type , std::string outputDir);

void compareHStacks_vs_Radion(TH1F* histo_DA, std::map<int,TH1F*> histo_MC, std::string xTitle, std::string label_DA, std::map<int,std::string> label_MC, std::map<int,int> color_MC,  std::string comp_type , std::string outputDir);

void compareNHistos(std::string input_histo_variables, const float x_min, const float x_max ,std::string xTitle , std::string outputName);

void compareNHistosOther(std::string input_histo_variables, const float x_min, const float x_max , const float maximum, const int rebin, std::string xTitle , std::string outputName);

void compareProfile(TProfile*, TProfile*, std::string ,std::string, std::string, std::string , const float y_min, const float y_max, std::string , std::string);

void compareNProfiles(std::string, const float, const float, const float, const float, std::string, std::string, std::string);

void drawHisto(TH1F* , const float , const float , std::string , std::string , std::string, std::string  );

void drawProfile(TProfile*, const float , const float, const float , const float , std::string , std::string, std::string, std::string );

double myChiSquareTest(TH1F*);

double EfficiencyError(int n1, int n2);

TLorentzVector Get_pho_p4(const TVector3* caloPosition, const float vtx_x, const float vtx_y, const float vtx_z, const float energy);

#endif
