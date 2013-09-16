#include "FCNObject.h"
#include "TMath.h"


extern FCNObject* myFCNObject;
double sqrt2Pi = 2.506628275;



//! ctor
FCNObject::FCNObject(const float& scale,
                     TTree* ntu,
                     const int& nParams,
		     const int& fitMinRangeMax,                     
		     const int& fitMaxRangeMin,                     
		     const int& fitMaxRangeMax,                     
		     const std::string& fitMode):
  m_ntu(ntu),
  m_mu(90.),
  m_sigma(2.),
  m_nParams(nParams),
  m_fitMinRangeMax(fitMinRangeMax),
  m_fitMaxRangeMin(fitMaxRangeMin),
  m_fitMaxRangeMax(fitMaxRangeMax),
  m_fitMode(fitMode),
  m_fval(0.),
  m_nCalls(0)
{
  std::cout << ">>> FCNObject::ctor <<<" << std::endl;
  
  
  // define the tree variables
  m_ntu -> SetBranchAddress("PhotonsMass",&m_var);
  m_ntu -> SetBranchAddress("evweight",&m_evweight);

  // define the fval map        
  for(int parIt = 0; parIt < m_nParams; ++parIt)
    m_fvalMap[parIt] = new TGraph();
  
  m_histoFit = new TH1F("m_histoFit", "", 100, 100., 300.);
}
// ---------------------------------------------------






//! dtor
FCNObject::~FCNObject()
{
  std::cout << ">>> FCNObject::dtor <<<" << std::endl;
  
  std::cout << ">>> m_nParams = " << m_nParams << std::endl;
//   for(int parIt = 0; parIt < m_nParams; ++parIt)
//   {
//     char graphName[50];
//     sprintf(graphName, "g_par%d_fval", parIt);
//     
//     m_fvalMap[parIt] -> Write(graphName);
//   }
  

//   m_outFile -> cd();
//   histo->Write();
//   m_outFile -> Close();

//   TCanvas* cout = new TCanvas();
//   histo->Draw();
//   cout->Print("histo.png", "png");

  delete m_histoFit;

}

// ---------------------------------------------------






//! the actual function
void FCNObject::FCN(int& /*nPar*/,
                    double* /*grad*/,
                    double& fval,
                    double* par, 
                    int /*iflag */)
{
  if( (m_nCalls % 100) == 0 )
    std::cout << ">>> FCNObject::FCN::Call " << m_nCalls << std::endl;
    
  //-----------------
  // compute the fval
  m_fval = 0.;
  
  for(int entry = 0; entry < m_ntu->GetEntries(); ++entry)
  {
    m_ntu -> GetEntry(entry);
  
    if(m_fitMode == "likelihood")
    {
       m_fval += log(m_sigma/par[0]*sqrt2Pi) + 0.5*(m_var-m_mu/par[0])*(m_var-m_mu/par[0])/((m_sigma/par[0])*(m_sigma/par[0]));
    }
    
    if(m_fitMode == "likelihood2")
    {
       m_fval += log(m_sigma/sqrt2Pi) + 0.5*(par[0]*m_var-m_mu)*(par[0]*m_var-m_mu)/(m_sigma*m_sigma);
    }

    if(m_fitMode == "likelihood3"){
//       std::cout << " >>> m_fitMode = " << m_fitMode << std::endl;
//       std::cout << " >>> m_var = " << m_var << std::endl;
      if( m_var < m_fitMinRangeMax || (m_var > m_fitMaxRangeMin && m_var < m_fitMaxRangeMax )){
	//	std::cout << " >>> m_var OK = " << m_var << std::endl;
       	//      m_fval += log(TMath::Exp(par[0] * m_var) * par[0] );
m_fval += log( pow(TMath::Exp(-1.*par[0]*m_var)*par[0]/(TMath::Exp(-100.*par[0])+TMath::Exp(-1.*m_fitMaxRangeMin*par[0])
							-TMath::Exp(-1.*m_fitMinRangeMax*par[0])-TMath::Exp(-1.*m_fitMaxRangeMax*par[0])), m_evweight) );
	//      m_fval += log( par[0]*TMath::Exp(par[1] * m_var) );
	if(m_nCalls == 0)  m_histoFit->Fill(m_var, m_evweight); 
      } 
    }
  }
  
  fval = m_fval;
  
  //std::cout << "\ncall: " << m_nCalls << "   fval: " << fval << "   par = " << par[8] << std::endl;
  //for(int i = 5; i < 8; ++i)
  //  std::cout << "par[" << i << "] = " << par[i] << std::endl;
  
  for(int parIt = 0; parIt < m_nParams; ++parIt)
    m_fvalMap[parIt] -> SetPoint(m_nCalls, par[parIt], fval);
  
  ++m_nCalls;
}

// ---------------------------------------------------







//! the function wrapper
void FCNObject::FCNWrapper(int& nPar,
                           double* grad,
                           double& fval,
                           double* par,
                          int iflag)
{
  FCNObject* dummy = (FCNObject*)myFCNObject;
  dummy -> FCN(nPar,
               grad,
               fval,
               par,
               iflag);
}



void FCNObject::passHisto(TH1F** histo,
			  const std::string& name){


  (*histo) = new TH1F(*m_histoFit);
  (*histo)->SetName(name.c_str());
//   std::cout << " >>>>>   (*histo)->GetEntries() =  " << (*histo)->GetEntries() << std::endl;
//   std::cout << " >>>>>   (m_histoFit)->GetEntries() =  " << m_histoFit->GetEntries() << std::endl;
}
