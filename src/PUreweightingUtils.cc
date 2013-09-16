#include "PUreweightingUtils.h"

void makePUDistributionMC(std::string inputTree, std::string inputFile, std::string inputBranch, std::string outputName ){

  float PUit_TrueNumInteractions;
  
  TChain* ntu_MC = new TChain(inputTree.c_str());
  ntu_MC -> Add(inputFile.c_str());
  ntu_MC->SetBranchStatus("*",0);
  ntu_MC->SetBranchStatus(inputBranch.c_str(),1);
  ntu_MC->SetBranchAddress(inputBranch.c_str(),&PUit_TrueNumInteractions);

  TH1F* h_PUit_TrueNumInteractions = new TH1F("h_PUit_TrueNumInteractions","h_PUit_TrueNumInteractions", 100,-0.5,99.5);
  
  for(int entry = 0; entry < ntu_MC->GetEntries(); entry++)
    {
	ntu_MC->GetEntry(entry);
        h_PUit_TrueNumInteractions->Fill(PUit_TrueNumInteractions);
    }

  TFile* output = new TFile(outputName.c_str(),"RECREATE");

  output->cd();
  h_PUit_TrueNumInteractions->Write();
  output->Close();

}

void ComputePUweights(std::string inputFileDA, std::string inputFileMC, std::string inputHistoNameDA, std::string inputHistoNameMC, std::string outputName){

   // histograms
  TH1F* h_mc = new TH1F("h_mc","",100,-0.5,99.5);
  TH1F* h_da = new TH1F("h_da","",100,-0.5,99.5);
  TH1F* temp = NULL;
  
  // files
  TFile* f_mc = TFile::Open(inputFileMC.c_str());
  TFile* f_da = TFile::Open(inputFileDA.c_str());
  
  
  
  // fill mc histogram
  std::cout << "\n\n\n>>> MC pileup histogram" << std::endl;
  
  temp = (TH1F*)( f_mc->Get(inputHistoNameMC.c_str()) );
  for(int bin = 1; bin <= temp->GetNbinsX(); ++bin)
  {
    float binCenter  = temp -> GetBinCenter(bin);
    float binContent = temp -> GetBinContent(bin);
    std::cout << "bin: " << bin << "   binCenter: " << binCenter << "   binContent: " << binContent << std::endl;
    
    h_mc -> Fill(binCenter,binContent);
  }
  h_mc -> Scale(1./h_mc->Integral());
  std::cout << "Integral: " << h_mc -> Integral() << std::endl;
  
  
  // fill da histogram
  std::cout << "\n\n\n>>> DA pileup histogram" << std::endl;
  
  temp = (TH1F*)( f_da->Get(inputHistoNameDA.c_str()) );
  for(int bin = 1; bin <= temp->GetNbinsX(); ++bin)
  {
    float binCenter  = temp -> GetBinCenter(bin) - 0.5;
    float binContent = temp -> GetBinContent(bin);
    std::cout << "bin: " << bin << "   binCenter: " << binCenter << "   binContent: " << binContent << std::endl;
    h_da -> Fill(binCenter,binContent);
  }
  h_da -> Scale(1./h_da->Integral());
  std::cout << "Integral: " << h_da -> Integral() << std::endl;
  
  
  // save PU weights
  TFile* f = new TFile(outputName.c_str(),"RECREATE");
  f -> cd();
  
  TH1F* h_PUweights = (TH1F*)(h_da->Clone("h_PUweights"));
  h_PUweights -> Divide(h_mc);
  
  h_da -> Write();
  h_mc -> Write();
  h_PUweights -> Write();
  
  f -> Close();

}


// default constructor, reading the corrections from file
TPileupReweighting::TPileupReweighting(std::string fname, std::string hname) {
  TFile weightsFile(fname.c_str(),"READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get(hname.c_str());
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();
  return; 
}

TPileupReweighting::~TPileupReweighting() { return;}

double TPileupReweighting::GetWeight(int npu){
  return (w[npu]);
}

