//g++ -Wall -o doControlBkg_Tree_Unbinned `root-config --cflags --glibs` FCNObject.cc doControlBkg_Tree_Unbinned.cpp

#include "FCNObject.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TChain.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>



FCNObject* myFCNObject;

#define CR1_Xlow 100.
#define CR1_Xmax 110.
#define CR2_Xlow 140.
#define CR2_Xmax 180.


Double_t fitNoSR(Double_t *x, Double_t *par)
{
  if(x[0] > CR1_Xmax && x[0] < CR2_Xlow){
    TF1::RejectPoint();
    return 0;
  }
  return (par[0]* exp(par[1] * x[0]));  
}

std::string labels2(int& bin)
{
  std::string Xaxis[25] = {
    "before #gamma ID",
    "#gamma pT cut",
    "#gamma CiC cut",
    /*4*/    "M#gamma#gamma cut",
    /*5*/    "njet >= 2",
    /*6*/    "nbjet >= 1",
    /*7*/    "njet >=2 + j selec.",
    /*8*/    "nbjet >=1 + j selec.",
    /*9*/    "1btag",
    /*10*/    "2btag",
    /*11*/    "1bt |cos #theta*|<.9",
    /*12*/    "2bt |cos #theta*|<.9",
    /*13*/    "1bt mjj cut", //" (90/170 for 1btag and 100/160 for 2btags)",    per ora tutto 90-170
    /*14*/    "2bt mjj cut", //" (90/170 for 1btag and 100/160 for 2btags)",
    /*15*/    "1bt kin fit M#gamma#gammajj cut", //" (260/335 and 255/320)",    per ora tutto 200-400
    /*16*/    "2bt kin fit M#gamma#gammajj cut", //" (260/335 and 255/320)",
    "1bt #DeltaR(#gamma,j) >= 1",
    "2bt #DeltaR(#gamma,j) >= 1",
    "1bt njets < 4",
    "2bt njets < 4"
  };
  return Xaxis[bin];
}


float ComputeChi2(TH1F* hist, TF1* func, TGraphErrors* grint, int& CR){

  int grintMax = (CR1_Xmax - CR1_Xlow) * 100.;
  float binMin = CR1_Xlow;
  float binMax = CR1_Xmax;
  if(CR == 2) {
    grintMax = (CR2_Xmax - CR2_Xlow) * 100.;
    binMin = CR2_Xlow;
    binMax = CR2_Xmax;
  }

  float chi2Value = 0.; 
  for(int xBin=1; xBin<hist->GetNbinsX(); ++xBin){
    if(hist->GetBinCenter(xBin) < binMin || hist->GetBinCenter(xBin) > binMax ||  hist->GetBinContent(xBin) == 0.) continue; 

    //    std::cout << " >>>>> hist->GetBinCenter(xBin) = " << hist->GetBinCenter(xBin) << std::endl;

    for(int iGr=0; iGr<grintMax; iGr++){
      double x,y;
      double Ex,Ey;
      grint->GetPoint(iGr, x, y);
      Ey = grint->GetErrorY(iGr);
      //      std::cout << " >>>>> x = " << x << " >>>> Ey = " << Ey << std::endl;
      if(hist->GetBinCenter(xBin) > x){
// 	std::cout << " >>>>> hist->GetBinCenter(xBin) = " << hist->GetBinCenter(xBin) 
// 		  << " >>>>> x = " << x 
// 		  << " >>>> Ey = " << Ey
// 		  << " >>>> hist->GetBinError(xBin) = " << hist->GetBinError(xBin) << std::endl;
	chi2Value += ( (pow(hist->GetBinContent(xBin) - func->Eval(hist->GetBinCenter(xBin)),2.)) /
		       (pow(hist->GetBinError(xBin),2.) + pow(Ey,2.)) );
	break;
      }
    }
  }
  return sqrt(chi2Value);  
}

void doFit(TTree* ntu, TH1F** histo, const std::string& name, const std::string& mode, double& p1, double& p1Err){
  std::cout << ">>> doFit <<<" << mode << std::endl;

  int nPar = 1;
  const int fitMinRangeMax = CR1_Xmax;
  const int fitMaxRangeMin = CR2_Xlow;
  const int fitMaxRangeMax = CR2_Xmax;
  myFCNObject = new FCNObject(1, ntu, nPar, fitMinRangeMax, fitMaxRangeMin, fitMaxRangeMax, mode);
  void (*myFCNPtr)(int&, double*, double&, double*, int) = &FCNObject::FCNWrapper;

  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,nPar);
  fitter -> SetFCN(myFCNPtr);
  //  fitter -> SetParameter(0, "p0", 12., 0.5,11.5,12.5);
  //  fitter -> SetParameter(0, "p0", 0.01, 0.005, 0.0, 0.02);
  //  fitter -> SetParameter(0, "p1", 0.01, 0.001, 0.005, 0.1);
  //    fitter -> SetParameter(0, "p1", p1, p1Err, p1-p1Err, p1+p1Err);
    if(p1-p1Err > 0.)  fitter -> SetParameter(0, "p1", p1, p1Err, p1-p1Err, p1+p1Err);
    else  fitter -> SetParameter(0, "p1", p1, p1Err, p1, p1+p1Err);
  //  fitter -> SetParameter(0, "p1", p1, p1Err, p1, p1+p1Err);

  // set print level                                                   
  double arglist[100];

  arglist[0] = 1;
  fitter -> ExecuteCommand("SET PRINT", arglist, 0);

  arglist[0] = 0.5;
  fitter -> ExecuteCommand("SET ERR", arglist, 1);

  // minimize                                     
  arglist[0] = 10000; // number of function calls 
  arglist[1] = 0.1*p1Err; // tolerance                 

  fitter -> ExecuteCommand("MIGRAD", arglist, 2);


  myFCNObject->passHisto(histo, name);
  std::cout << "histo->GetEntries()= " << (*histo)->GetEntries() << std::endl;
  delete myFCNObject;

  p1 = fitter->GetParameter(0);
  p1Err = fitter->GetParError(0);
  return;
}


int main(){
  //  gROOT->ProcessLine(".x ~/public/setTDRStyle.C");   
  gROOT->ProcessLine(".x ~/public/style.C");   
  std::cout << " ci sono " << std::endl;


  std::string nStep = "step12";
  std::string cutName = "1btag_Mjjcut";
  std::string cutNameLabel = "1btag Mjj in [90,170] cut";
  double p1_DA = 0.0173;          //-1.73320e-02;
  double p1Err_DA = 0.0024;       //2.40755e-030;
  double p1_MCDiPho = 0.0051;     //-5.09748e-03;
  double p1Err_MCDiPho = 0.0029;   //2.91910e-030;
  double p1_MCtt = 0.0105;         //1.05323e-02   
  double p1Err_MCtt = 0.0022;      //2.16894e-03 
  double p1_MC = 0.0193;          // -1.93267e-02 
  double p1Err_MC = 0.0016;       // 1.62999e-03
  //4.56572e+00   4.87092e-01           
  p1_DA = 4.566;
  p1Err_DA = 0.487;
  //3.25809e+00   2.21697e+00           
  p1_MC = 3.258;
  p1Err_MC = 2.217;
  //2.54628e+00   6.88845e-01           
  p1_MCDiPho = 2.546;
  p1Err_MCDiPho = 0.689;
  //                                    
  p1_MCtt =  1.;
  p1Err_MCtt = 4.;



  if(nStep == "step13"){
    cutName = "2btag_Mjjcut";
    cutNameLabel = "2btag Mjj in [90,170] cut";
    //2.00709e+00   9.31653e-01
    p1_DA = 2.007;
    p1Err_DA = 0.9317;

    //1.22847e-01   2.70784e+00
    p1_MC = 0.1;
    p1Err_MC = 1.;

    //6.82074e+00   3.99996e+00
    p1_MCDiPho = 6.821;
    p1Err_MCDiPho = 3.999;

    //-5.47964e+00   4.09947e+01
    p1_MCtt = 5.479;
    p1Err_MCtt = 2.;
    }
  if(nStep == "step14"){
    cutName = "1btag_Mjjggcut";
    cutNameLabel = "1btag Mjj in [90,170] cut M#gamma#gammajj in [200,400]";
    //6.98280e+00   1.54484e+00
    p1_DA = 6.983;
    p1Err_DA = 1.545;

    //5.09054e+00   3.18668e+00
    p1_MC = 5.091;
    p1Err_MC = 3.187;

    //6.82074e+00   3.99996e+00
    p1_MCDiPho = 6.821;
    p1Err_MCDiPho = 3.999;

    //-5.22793e-01   4.04227e+01
    p1_MCtt = 1.;
    p1Err_MCtt = 4;
  }
  if(nStep == "step15"){
    cutName = "2btag_Mjjggcut";
    cutNameLabel = "2btag Mjj in [90,170] cut M#gamma#gammajj in [200,400]";
    //2.87452e+00   4.51589e+00   
    p1_DA = 2.;
    p1Err_DA = 4.;

    //
    p1_MCDiPho =  0.001;
    p1Err_MCDiPho = 0.007;

    p1_MC = 0.0053;                      //-5.29939e-03
    p1Err_MC = 0.0084;                   //8.44790e-03

    //-6.68006e+00   6.82757e+01 
    p1_MCtt = 6.68;
    p1Err_MCtt = 6.8;
  }



  std::string outPlotDir = "UnbinnedFit_plot";

  TFile* fDA = new TFile(("output_"+nStep+"/Plots_tree_Data.root").c_str(), "read");
  //  TFile* fMC = new TFile(("output_"+nStep+"/Plots_tree_MCBackground_TTJets.root").c_str(), "read");
  TFile* fDiPho = new TFile(("output_"+nStep+"/Plots_tree_MCBackgrounds_DiPhotonJets-madgraph_diphojet_8TeV.root").c_str(), "read");

  TH1F* h_diPhoMass_DA;
  TH1F* h_diPhoMass_MC;
  TH1F* h_diPhoMass_MCDiPho;
  TH1F* h_diPhoMass_MCtt;

  TTree* ntu_DA = (TTree*)fDA->Get("EeventTree");
  TTree* ntu_MCDiPho = (TTree*)fDiPho->Get("EeventTree"); 
  TChain* ntu_MCtt = new TChain("EeventTree");
  ntu_MCtt->Add(("output_"+nStep+"/Plots_tree_MCBackgrounds_ttgg_Zgg_ttgg_dR02.root").c_str());
  TChain* ntu_MCgj = new TChain("EeventTree");
  ntu_MCgj->Add(("output_"+nStep+"/Plots_tree_MCBackground_TTGJets.root").c_str());
  TChain* ntu_MCH = new TChain("EeventTree");
  ntu_MCH->Add(("output_"+nStep+"/Plots_tree_MCBackground_ggh_m125_8TeV.root").c_str());
  ntu_MCH->Add(("output_"+nStep+"/Plots_tree_MCBackground_tth_m125_8TeV.root").c_str());
  ntu_MCH->Add(("output_"+nStep+"/Plots_tree_MCBackground_vbf_m125_8TeV.root").c_str());
  ntu_MCH->Add(("output_"+nStep+"/Plots_tree_MCBackground_wzh_m125_8TeV.root").c_str());
  TChain* ntu_MC = new TChain("EeventTree");
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_TTGJets.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_TTJets.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_ZZJetsTo2L2Q.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_dipho_Box_250_8TeV.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_dipho_Box_25_8TeV.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_ggh_m125_8TeV.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_gjets_20_8TeV_pf.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_gjets_40_8TeV_pf.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_qcd_30_8TeV_ff.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_qcd_30_8TeV_pf.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_qcd_40_8TeV_ff.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_qcd_40_8TeV_pf.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_tth_m125_8TeV.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_vbf_m125_8TeV.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackground_wzh_m125_8TeV.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackgrounds_DY_WG_DYJetsToLL.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackgrounds_DY_WG_WGToLNuG.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackgrounds_DiPhotonJets-madgraph_diphojet_8TeV.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackgrounds_ttgg_Zgg_Zgg_dR02.root").c_str());
  ntu_MC->Add(("output_"+nStep+"/Plots_tree_MCBackgrounds_ttgg_Zgg_ttgg_dR02.root").c_str());



  std::cout << " ntu_DA->GetEntries() = " << ntu_DA->GetEntries() << std::endl;
  std::cout << " ntu_MC->GetEntries() = " << ntu_MC->GetEntries() << std::endl;
  std::cout << " ntu_MCDiPho->GetEntries() = " << ntu_MCDiPho->GetEntries() << std::endl;
  std::cout << " ntu_MCtt->GetEntries() = " << ntu_MCtt->GetEntries() << std::endl;
  std::cout << " ntu_MCgj->GetEntries() = " << ntu_MCgj->GetEntries() << std::endl;
  std::cout << " ntu_MCH->GetEntries() = " << ntu_MCH->GetEntries() << std::endl;

  ntu_DA->SetDirectory(0);
  ntu_MC->SetDirectory(0);
  ntu_MCDiPho->SetDirectory(0);
  ntu_MCtt->SetDirectory(0);

  fDA->Close();
  //  fMC->Close();
  fDiPho->Close();
  std::cout << " >>>> Fine 1 " << std::endl;

  //  return 100; 

  //////// Fit method  

  float DAValueFit_CR1 = 0.;
  float DAValueFit_CR2 = 0.;
  float DAValueFit_SR = 0.;

  float MCDiPhoValueFit_CR1 = 0.;
  float MCDiPhoValueFit_CR2 = 0.;
  float MCDiPhoValueFit_SR = 0.;

  float MCttValueFit_CR1 = 0.;
  float MCttValueFit_CR2 = 0.;
  float MCttValueFit_SR = 0.;

  float MCValueFit_CR1 = 0.;
  float MCValueFit_CR2 = 0.;
  float MCValueFit_SR = 0.;

  float DAEstimFitDiPho_SR_fCR1 = 0.;
  float DAEstimFitDiPho_SR_fCR2 = 0.;

  float DAEstimFit_SR_fCR1 = 0.;
  float DAEstimFit_SR_fCR2 = 0.;

  float DAEstimFittt_SR_fCR1 = 0.;
  float DAEstimFittt_SR_fCR2 = 0.;

  /////////////////////////
  std::cout << " >>> normal binned fit " << std::endl;
  float DAValueB_CR1 = 0.;
  float DAValueB_CR2 = 0.;
  float DAValueB_SR = 0.;

  float MCDiPhoValueB_CR1 = 0.;
  float MCDiPhoValueB_CR2 = 0.;
  float MCDiPhoValueB_SR = 0.;

  float MCttValueB_CR1 = 0.;
  float MCttValueB_CR2 = 0.;
  float MCttValueB_SR = 0.;

  float MCValueB_CR1 = 0.;
  float MCValueB_CR2 = 0.;
  float MCValueB_SR = 0.;

  float DAEstimBDiPho_SR_fCR1 = 0.;
  float DAEstimBDiPho_SR_fCR2 = 0.;

  float DAEstimBtt_SR_fCR1 = 0.;
  float DAEstimBtt_SR_fCR2 = 0.;

  float DAEstimB_SR_fCR1 = 0.;
  float DAEstimB_SR_fCR2 = 0.;


  std::string histoName;

  TF1* DA_Fit;
  TF1* MC_Fit;
  TF1* MCDiPho_Fit;
  TF1* MCtt_Fit;

  char nomeFunc[100];
  int maxRangeFuncMax = int(CR2_Xmax);
  int maxRangeFuncMin = int(CR2_Xlow);
  int minRangeFuncMax = int(CR1_Xmax);
  //  sprintf(nomeFunc, "[0]*[1]*exp(-x*[1])/(exp(-100.*[1]) - exp(-[1]*%d))",maxRangeFunc);
  sprintf(nomeFunc, "[0]*[1]*exp(-x*[1])/(exp(-100.*[1]) + exp(-[1]*%d) - exp(-[1]*%d) - exp(-[1]*%d))",maxRangeFuncMin, minRangeFuncMax, maxRangeFuncMax);

  DA_Fit = new TF1("DA_Fit", nomeFunc, CR1_Xlow, CR2_Xmax);
  MC_Fit = new TF1("MC_Fit", nomeFunc, CR1_Xlow, CR2_Xmax);
  MCDiPho_Fit = new TF1("MCDiPho_Fit", nomeFunc, CR1_Xlow, CR2_Xmax);
  MCtt_Fit    = new TF1("MCtt_Fit", nomeFunc, CR1_Xlow, CR2_Xmax);

  //DA
  histoName = "h_diPhoMass_DA";
  doFit(ntu_DA, &(h_diPhoMass_DA), histoName, "likelihood3", p1_DA, p1Err_DA);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> DA p1 = " << p1_DA << " p1Err = " << p1Err_DA << std::endl;

  //MC All
  histoName = "h_diPhoMass_MC";
  doFit(ntu_MC, &(h_diPhoMass_MC), histoName, "likelihood3", p1_MC, p1Err_MC);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> MC p1 = " << p1_MC << " p1Err = " << p1Err_MC << std::endl;

  //MC DiPho
  histoName = "h_diPhoMass_MCDiPho";
  doFit(ntu_MCDiPho, &(h_diPhoMass_MCDiPho), histoName, "likelihood3", p1_MCDiPho, p1Err_MCDiPho);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> MCDiPho p1 = " << p1_MCDiPho << " p1Err = " << p1Err_MCDiPho << std::endl;

  //MC tt
  histoName = "h_diPhoMass_MCtt";
  doFit(ntu_MCtt, &(h_diPhoMass_MCtt), histoName, "likelihood3", p1_MCtt, p1Err_MCtt);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> MCtt p1 = " << p1_MCtt << " p1Err = " << p1Err_MCtt << std::endl;


//   h_diPhoMass_DA->Scale(1./h_diPhoMass_DA->Integral());
//   h_diPhoMass_MC->Scale(1./h_diPhoMass_MC->Integral());
//   h_diPhoMass_MCDiPho->Scale(1./h_diPhoMass_MCDiPho->Integral());
//   h_diPhoMass_MCtt->Scale(1./h_diPhoMass_MCtt->Integral());


  ////////// normal binned
  TF1* DA_FitB = new TF1("DA_FitB", "[0]*exp(-x*[1])", CR1_Xlow, CR2_Xmax);
  TF1* MC_FitB = new TF1("MC_FitB", "[0]*exp(-x*[1])", CR1_Xlow, CR2_Xmax);
  TF1* MCDiPho_FitB = new TF1("MCDiPho_FitB", "[0]*exp(-x*[1])", CR1_Xlow, CR2_Xmax);
  TF1* MCtt_FitB    = new TF1("MCtt_FitB", "[0]*exp(-x*[1])", CR1_Xlow, CR2_Xmax);

  DA_FitB->SetParameters(h_diPhoMass_DA->Integral(), p1_DA);
  MC_FitB->SetParameters(h_diPhoMass_MC->Integral(), p1_MC);
  MCDiPho_FitB->SetParameters(h_diPhoMass_MCDiPho->Integral(), p1_MCDiPho);
  MCtt_FitB->SetParameters(h_diPhoMass_MCtt->Integral(), p1_MCtt);


  h_diPhoMass_DA->Fit("DA_FitB", "L");
  h_diPhoMass_MC->Fit("MC_FitB", "L");
  h_diPhoMass_MCDiPho->Fit("MCDiPho_FitB", "L");
  h_diPhoMass_MCtt->Fit("MCtt_FitB", "L");

  //////// Setting Functions 1
  p1_DA = DA_FitB->GetParameter(1);
  p1Err_DA = DA_FitB->GetParError(1);
  p1_MC = MC_FitB->GetParameter(1);
  p1Err_MC = MC_FitB->GetParError(1);
  p1_MCDiPho = MCDiPho_FitB->GetParameter(1);
  p1Err_MCDiPho = MCDiPho_FitB->GetParError(1);
  p1_MCtt = MCtt_FitB->GetParameter(1);
  p1Err_MCtt = MCtt_FitB->GetParError(1);


  //refit unbinned
  //DA                                                                                                                               
  h_diPhoMass_DA->Reset();
  histoName = "h_diPhoMass_DA";
  doFit(ntu_DA, &(h_diPhoMass_DA), histoName, "likelihood3", p1_DA, p1Err_DA);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> DA p1 = " << p1_DA << " p1Err = " << p1Err_DA << std::endl;
  //MC All                                                                                                                           
  histoName = "h_diPhoMass_MC";
  doFit(ntu_MC, &(h_diPhoMass_MC), histoName, "likelihood3", p1_MC, p1Err_MC);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> MC p1 = " << p1_MC << " p1Err = " << p1Err_MC << std::endl;
  //MC DiPho                                                                                                                         
  histoName = "h_diPhoMass_MCDiPho";
  doFit(ntu_MCDiPho, &(h_diPhoMass_MCDiPho), histoName, "likelihood3", p1_MCDiPho, p1Err_MCDiPho);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> MCDiPho p1 = " << p1_MCDiPho << " p1Err = " << p1Err_MCDiPho << std::endl;
  //MC tt                                                                                                                            
  histoName = "h_diPhoMass_MCtt";
  doFit(ntu_MCtt, &(h_diPhoMass_MCtt), histoName, "likelihood3", p1_MCtt, p1Err_MCtt);
  std::cout << " %%%%%%%%%%%%  >>>>>>>>>>> MCtt p1 = " << p1_MCtt << " p1Err = " << p1Err_MCtt << std::endl;


  //////// Setting Functions 2
  //  DA_Fit->FixParameter(0, h_diPhoMass_DA->Integral()*5.);
  //DA_Fit->SetParameter(0, h_diPhoMass_DA->GetEntries());
  DA_Fit->SetParameter(0, h_diPhoMass_DA->Integral("width"));
  //  DA_Fit->FixParameter(0, h_diPhoMass_DA->GetBinContent(1));
  //DA_Fit->FixParameter(0, 1.);
  DA_Fit->SetParameter(1, p1_DA);

  std::cout << " ################ diPhoMass_DA->Integral width = " << h_diPhoMass_DA->Integral("width") 
 	    << " DA_FitB->GetParameter(0) = " << DA_FitB->GetParameter(0) << " \n "
 	    << " DA_Fit->GetParameter(0) = " << DA_Fit->GetParameter(0) << " \n "
 	    << " DA_FitB->GetParameter(1) = " << DA_FitB->GetParameter(1) << " \n "
 	    << " DA_Fit->GetParameter(1) = " << DA_Fit->GetParameter(1) << std::endl;
  
  //MC_Fit->FixParameter(0, h_diPhoMass_MC->GetEntries());
  // MC_Fit->FixParameter(0, h_diPhoMass_MC->Integral()*5.);
    MC_Fit->FixParameter(0, h_diPhoMass_MC->Integral("width"));
  //MC_Fit->FixParameter(0, h_diPhoMass_MC->GetBinContent(1));
  //  MC_Fit->FixParameter(0, 1.);
  MC_Fit->FixParameter(1, p1_MC);
  
  //MCDiPho_Fit->FixParameter(0, h_diPhoMass_MCDiPho->GetEntries());
  //MCDiPho_Fit->FixParameter(0, h_diPhoMass_MCDiPho->Integral()*5.);
  MCDiPho_Fit->FixParameter(0, h_diPhoMass_MCDiPho->Integral("width"));
  //  MCDiPho_Fit->FixParameter(0, h_diPhoMass_MCDiPho->GetBinContent(1));
  //MCDiPho_Fit->FixParameter(0, 1.);
  MCDiPho_Fit->FixParameter(1, p1_MCDiPho);

  //MCtt_Fit->FixParameter(0, h_diPhoMass_MCtt->GetEntries());
  //  MCtt_Fit->FixParameter(0, h_diPhoMass_MCtt->Integral()*5.);
  MCtt_Fit->FixParameter(0, h_diPhoMass_MCtt->Integral("width"));
  //  MCtt_Fit->FixParameter(0, h_diPhoMass_MCtt->GetBinContent(1));
  //MCtt_Fit->FixParameter(0, 1.);
  MCtt_Fit->FixParameter(1, p1_MCtt);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Estimate from fit
	
  DAValueFit_CR1 = DA_Fit->Integral(CR1_Xlow, CR1_Xmax);
  DAValueFit_CR2 = DA_Fit->Integral(CR2_Xlow, CR2_Xmax);
  DAValueFit_SR = DA_Fit->Integral(120., 130.);
	
  MCValueFit_CR1 = MC_Fit->Integral(CR1_Xlow, CR1_Xmax);
  MCValueFit_CR2 = MC_Fit->Integral(CR2_Xlow, CR2_Xmax);
  MCValueFit_SR  = MC_Fit->Integral(120., 130.);
	
  MCDiPhoValueFit_CR1 = MCDiPho_Fit->Integral(CR1_Xlow, CR1_Xmax);
  MCDiPhoValueFit_CR2 = MCDiPho_Fit->Integral(CR2_Xlow, CR2_Xmax);
  MCDiPhoValueFit_SR  = MCDiPho_Fit->Integral(120., 130.);

  MCttValueFit_CR1 = MCtt_Fit->Integral(CR1_Xlow, CR1_Xmax);
  MCttValueFit_CR2 = MCtt_Fit->Integral(CR2_Xlow, CR2_Xmax);
  MCttValueFit_SR  = MCtt_Fit->Integral(120., 130.);
	
  std::cout << " >>>>>>>> Measured by fit ****** " << nStep << std::endl;
  std::cout << " >>>> DA_FitB p0 = " << DA_FitB->GetParameter(0) << std::endl;
  std::cout << " >>>> h_diPhoMass_DA->Integral() = " << h_diPhoMass_DA->Integral() << std::endl;
  std::cout << " >>>> h_diPhoMass_DA->GetEntries() = " << h_diPhoMass_DA->GetEntries() << std::endl;
//   std::cout << " >>>> MC_Fit p0 = " << MC_Fit->GetParameter(0) << std::endl;
//   std::cout << " >>>> MC_Fit p1 = " << MC_Fit->GetParameter(1) << std::endl;
//   std::cout << " >>>> MCDiPho_Fit p0 = " << MCDiPho_Fit->GetParameter(0) << std::endl;
//   std::cout << " >>>> MCDiPho_Fit p1 = " << MCDiPho_Fit->GetParameter(1) << std::endl;
	
//   std::cout << " >>>>>>>> Measured by fit UNBINNED ****** " << nStep << std::endl;
//   std::cout << " MC_SR = " << MCValueFit_SR << std::endl;
//   std::cout << " MCDiPho_SR = " << MCDiPhoValueFit_SR << std::endl;
	
//   std::cout << " DA_CR1 = " << DAValueFit_CR1 << std::endl;
//   std::cout << " MC_CR1 = " << MCValueFit_CR1 << std::endl;
//   std::cout << " MCDiPho_CR1 = " << MCDiPhoValueFit_CR1 << std::endl;
	
//   std::cout << " DA_CR2 = " << DAValueFit_CR2 << std::endl;
//   std::cout << " MC_CR2 = " << MCValueFit_CR2 << std::endl;
//   std::cout << " MCDiPho_CR2 = " << MCDiPhoValueFit_CR2 << std::endl;
	
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  DAValueB_CR1 = DA_FitB->Integral(CR1_Xlow, CR1_Xmax);
  DAValueB_CR2 = DA_FitB->Integral(CR2_Xlow, CR2_Xmax);
  DAValueB_SR = DA_FitB->Integral(120., 130.);
  
  MCValueB_CR1 = MC_FitB->Integral(CR1_Xlow, CR1_Xmax);
  MCValueB_CR2 = MC_FitB->Integral(CR2_Xlow, CR2_Xmax);
  MCValueB_SR  = MC_FitB->Integral(120., 130.);
	
  MCDiPhoValueB_CR1 = MCDiPho_FitB->Integral(CR1_Xlow, CR1_Xmax);
  MCDiPhoValueB_CR2 = MCDiPho_FitB->Integral(CR2_Xlow, CR2_Xmax);
  MCDiPhoValueB_SR  = MCDiPho_FitB->Integral(120., 130.);

  MCttValueB_CR1 = MCtt_FitB->Integral(CR1_Xlow, CR1_Xmax);
  MCttValueB_CR2 = MCtt_FitB->Integral(CR2_Xlow, CR2_Xmax);
  MCttValueB_SR  = MCtt_FitB->Integral(120., 130.);
	
//   std::cout << " >>>>>>>> Measured by fitBinned ****** " << nStep << std::endl;
//   std::cout << " MC_SR = " << MCValueB_SR << std::endl;
//   std::cout << " MCDiPho_SR = " << MCDiPhoValueB_SR << std::endl;
//   std::cout << " DA_CR1 = " << DAValueB_CR1 << std::endl;
//   std::cout << " MC_CR1 = " << MCValueB_CR1 << std::endl;
//   std::cout << " MCDiPho_CR1 = " << MCDiPhoValueB_CR1 << std::endl;
//   std::cout << " DA_CR2 = " << DAValueB_CR2 << std::endl;
//   std::cout << " MC_CR2 = " << MCValueB_CR2 << std::endl;
//   std::cout << " MCDiPho_CR2 = " << MCDiPhoValueB_CR2 << std::endl;
  
	
  ////// ESTIMATE DA SR unbinned fit
  if(MCDiPhoValueFit_CR1 != 0.)  DAEstimFitDiPho_SR_fCR1 = MCDiPhoValueFit_SR * (DAValueFit_CR1/MCDiPhoValueFit_CR1);
  if(MCDiPhoValueFit_CR2 != 0.)  DAEstimFitDiPho_SR_fCR2 = MCDiPhoValueFit_SR * (DAValueFit_CR2/MCDiPhoValueFit_CR2);
  std::cout << " ESTIMATE DA in SR => DA func fitUNbinned ****** step " << nStep << std::endl;
  std::cout << " DA_CR2/DA_CR1  = " << DAValueFit_CR2/DAValueFit_CR1 << std::endl;
  std::cout << " DA_CR2 = " << DAValueFit_CR2 << std::endl;
  std::cout << " DA_CR1  = " << DAValueFit_CR1 << std::endl;
  std::cout << " DA_SR  = " << DAValueFit_SR << std::endl;

  std::cout << " ESTIMATE DA in SR => DiPho MC func fitUNbinned ****** step " << nStep << std::endl;
  std::cout << " MC_CR2/MC_CR1  = " << MCDiPhoValueFit_CR2/MCDiPhoValueFit_CR1 << std::endl;
  std::cout << " MC_CR2  = " << MCDiPhoValueFit_CR2 << std::endl;
  std::cout << " MC_CR1  = " << MCDiPhoValueFit_CR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR1  = " << DAEstimFitDiPho_SR_fCR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR2  = " << DAEstimFitDiPho_SR_fCR2 << std::endl;

  if(MCttValueFit_CR1 != 0.)  DAEstimFittt_SR_fCR1 = MCttValueFit_SR * (DAValueFit_CR1/MCttValueFit_CR1);
  if(MCttValueFit_CR2 != 0.)  DAEstimFittt_SR_fCR2 = MCttValueFit_SR * (DAValueFit_CR2/MCttValueFit_CR2);
  std::cout << " ESTIMATE DA in SR => tt MC func fitUNbinned ****** step " << nStep << std::endl;
  std::cout << " MC_CR2/MC_CR1  = " << MCttValueFit_CR2/MCttValueFit_CR1 << std::endl;
  std::cout << " MC_CR2  = " << MCttValueFit_CR2 << std::endl;
  std::cout << " MC_CR1  = " << MCttValueFit_CR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR1  = " << DAEstimFittt_SR_fCR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR2  = " << DAEstimFittt_SR_fCR2 << std::endl;
	
  if(MCValueFit_CR1 != 0.)  DAEstimFit_SR_fCR1 = MCValueFit_SR * (DAValueFit_CR1/MCValueFit_CR1);
  if(MCValueFit_CR2 != 0.)  DAEstimFit_SR_fCR2 = MCValueFit_SR * (DAValueFit_CR2/MCValueFit_CR2);
  std::cout << " ESTIMATE DA in SR => All MC func fitUNbinned ****** step " << nStep << std::endl;
  std::cout << " MC_CR2/MC_CR1  = " << MCValueFit_CR2/MCValueFit_CR1 << std::endl;
  std::cout << " MC_CR2  = " << MCValueFit_CR2 << std::endl;
  std::cout << " MC_CR1  = " << MCValueFit_CR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR1  = " << DAEstimFit_SR_fCR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR2  = " << DAEstimFit_SR_fCR2 << std::endl;
  

  ////// ESTIMATE DA SR   Binned	
  if(MCDiPhoValueB_CR1 != 0.)  DAEstimBDiPho_SR_fCR1 = MCDiPhoValueB_SR * (DAValueB_CR1/MCDiPhoValueB_CR1);
  if(MCDiPhoValueB_CR2 != 0.)  DAEstimBDiPho_SR_fCR2 = MCDiPhoValueB_SR * (DAValueB_CR2/MCDiPhoValueB_CR2);
  std::cout << " ESTIMATE DA in SR => DiPho MC func fitBinned ****** step " << nStep << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR1  = " << DAEstimBDiPho_SR_fCR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR2  = " << DAEstimBDiPho_SR_fCR2 << std::endl;

  if(MCttValueB_CR1 != 0.)  DAEstimBtt_SR_fCR1 = MCttValueB_SR * (DAValueB_CR1/MCttValueB_CR1);
  if(MCttValueB_CR2 != 0.)  DAEstimBtt_SR_fCR2 = MCttValueB_SR * (DAValueB_CR2/MCttValueB_CR2);
  std::cout << " ESTIMATE DA in SR => tt MC func fitBinned ****** step " << nStep << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR1  = " << DAEstimBtt_SR_fCR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR2  = " << DAEstimBtt_SR_fCR2 << std::endl;

  if(MCValueB_CR1 != 0.)  DAEstimB_SR_fCR1 = MCValueB_SR * (DAValueB_CR1/MCValueB_CR1);
  if(MCValueB_CR2 != 0.)  DAEstimB_SR_fCR2 = MCValueB_SR * (DAValueB_CR2/MCValueB_CR2);
  std::cout << " ESTIMATE DA in SR => All MC func fitBinned ****** step " << nStep << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR1  = " << DAEstimB_SR_fCR1 << std::endl;
  std::cout << " MC_SR * (DA/MC)_CR2  = " << DAEstimB_SR_fCR2 << std::endl;
	
	
  ////// DRAW 
  //   h_diPhoMass_DA->Rebin();
  //   h_diPhoMass_MC->Rebin();
  //   h_diPhoMass_MCDiPho->Rebin();
	
  DA_Fit->SetLineColor(kBlue);
  DA_Fit->SetLineWidth(1);
  MCDiPho_Fit->SetLineColor(kCyan);
  MCDiPho_Fit->SetLineWidth(1);
  MCtt_Fit->SetLineColor(kGreen);
  MCtt_Fit->SetLineWidth(1);
  MC_Fit->SetLineColor(kRed);
  MC_Fit->SetLineWidth(1);
  
  DA_FitB->SetLineColor(kBlue);
  DA_FitB->SetLineWidth(1);
  MCDiPho_FitB->SetLineColor(kCyan);
  MCDiPho_FitB->SetLineWidth(1);
  MCtt_FitB->SetLineColor(kGreen);
  MCtt_FitB->SetLineWidth(1);
  MC_FitB->SetLineColor(kRed);
  MC_FitB->SetLineWidth(1);
	
  h_diPhoMass_DA->SetMarkerColor(kBlue);
  h_diPhoMass_DA->SetLineColor(kBlue);
  h_diPhoMass_DA->SetMarkerStyle(20);
  h_diPhoMass_MCDiPho->SetMarkerColor(kCyan);
  h_diPhoMass_MCDiPho->SetLineColor(kCyan);
  h_diPhoMass_MCDiPho->SetMarkerStyle(21);
  h_diPhoMass_MCtt->SetMarkerColor(kGreen);
  h_diPhoMass_MCtt->SetLineColor(kGreen);
  h_diPhoMass_MCtt->SetMarkerStyle(21);
  h_diPhoMass_MC->SetMarkerColor(kRed);
  h_diPhoMass_MC->SetLineColor(kRed);
  h_diPhoMass_MC->SetMarkerStyle(21);
  
  TLegend *legF = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  legF->SetTextFont(42);
  legF->SetFillColor(kWhite);
  legF->SetLineColor(kWhite);
  legF->SetShadowColor(kWhite);
  
  legF->AddEntry(h_diPhoMass_DA, "DA","pl");
  legF->AddEntry(h_diPhoMass_MC, "MC all ","pl");
  legF->AddEntry(h_diPhoMass_MCDiPho, "MCDiPho","pl");
  legF->AddEntry(h_diPhoMass_MCtt, "MCtt","pl");
  
//   h_diPhoMass_DA->Scale(1./h_diPhoMass_DA->Integral());
//   h_diPhoMass_MC->Scale(1./h_diPhoMass_MC->Integral());
//   h_diPhoMass_MCDiPho->Scale(1./h_diPhoMass_MCDiPho->Integral());
//   h_diPhoMass_MCtt->Scale(1./h_diPhoMass_MCtt->Integral());

  TCanvas* c_cout = new TCanvas();
  c_cout->cd();
  h_diPhoMass_DA->GetXaxis()->SetTitle(("M_{#gamma#gamma} - "+cutNameLabel).c_str());
  //h_diPhoMass_DA->GetYaxis()->SetRangeUser(0., 0.1);
  h_diPhoMass_DA->Draw("e");
  h_diPhoMass_MCDiPho->Draw("samee");
  h_diPhoMass_MCtt->Draw("samee");
  h_diPhoMass_MC->Draw("samee");
  DA_Fit->Draw("same");
  MCDiPho_Fit->Draw("same");
  MCtt_Fit->Draw("same");
  MC_Fit->Draw("same");
  legF->Draw("same");
  c_cout->Print((outPlotDir+"/DiPho_"+cutName+".png").c_str(), "png");

  TFile outFile((outPlotDir+"/outFile_DiPho_"+cutName+".root").c_str(), "recreate");
  h_diPhoMass_DA->Write();
  h_diPhoMass_MCDiPho->Write();
  h_diPhoMass_MCtt->Write();
  h_diPhoMass_MC->Write();
  MCDiPho_Fit->Write();
  MCtt_Fit->Write();
  MC_Fit->Write();
  DA_Fit->Write();
  outFile.Close();
  
  ////// DRAW
  TCanvas* c_coutB = new TCanvas();
  c_coutB->cd();
  h_diPhoMass_DA->GetXaxis()->SetTitle(("M_{#gamma#gamma} - "+cutNameLabel).c_str());
  h_diPhoMass_DA->Draw("e");
  h_diPhoMass_MCDiPho->Draw("samee");
  h_diPhoMass_MCtt->Draw("samee");
  h_diPhoMass_MC->Draw("samee");
  DA_FitB->Draw("same");
  MCDiPho_FitB->Draw("same");
  MCtt_FitB->Draw("same");
  MC_FitB->Draw("same");
  legF->Draw("same");
  c_coutB->Print((outPlotDir+"/DiPho_"+cutName+"_FitBinned.png").c_str(), "png");

  TFile outFileB((outPlotDir+"/outFile_DiPho_"+cutName+"_FitBinned.root").c_str(), "recreate");
  h_diPhoMass_DA->Write();
  h_diPhoMass_MCDiPho->Write();
  h_diPhoMass_MCtt->Write();
  h_diPhoMass_MC->Write();
  MCDiPho_FitB->Write();
  MC_FitB->Write();
  MCtt_FitB->Write();
  DA_FitB->Write();
  outFileB.Close();


  return 0;
}
