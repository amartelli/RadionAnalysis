#ifndef FCNObject_h
#define FCNObject_h

#include <iostream>
#include <map>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"



class FCNObject
{
 public:
  
  //! ctor
  FCNObject(const float& scale,
            TTree* ntu,
            const int& nParams,
	    const int& fitMinRangeMax,
	    const int& fitMaxRangeMin,
	    const int& fitMaxRangeMax,
	    const std::string& fitMode = "likelihood");
  
  //! dtor
  ~FCNObject();
  
  
  //! the actual function
  void FCN(int& /*nPar*/,
           double* /*grad*/,
           double& fval,
           double* par,
           int /*iflag */);
  
  
  //! the function wrapper
  static void FCNWrapper(int& /*nPar*/,
                         double* /*grad*/,
                         double& fval,
                         double* par,
                         int /*iflag */);
  
  void passHisto(TH1F** histo, 
		 const std::string& name);

 private:
  
  TTree* m_ntu;
  
  Float_t m_var;
  Float_t m_evweight;
  float m_mu;
  float m_sigma;
  
  int m_nParams;
  int m_fitMinRangeMax;
  int m_fitMaxRangeMin;
  int m_fitMaxRangeMax;
  std::string m_fitMode;


  long double m_fval;
  int m_nCalls;
  
  //  TFile* m_outFile;

  TH1F* m_histoFit;  
  
  std::map<int, TGraph*> m_fvalMap;
};

#endif
