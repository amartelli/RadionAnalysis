// ------------------------------------------------
// 
// GlobalFitter - kinematic fit for X->lvjj
//
// ------------------------------------------------
#ifndef GlobalFitter_h
#define GlobalFitter_h

#include "TKinFitter.h"
#include "TLorentzVector.h"


class MissingEnergy1{

 public:

  void SetSumEt(float SumEt);
  void SetNeutrino(TLorentzVector Neutrino);

  float SumEt_;
  TLorentzVector Neutrino_;

};



class GlobalFitter : public TKinFitter {

 public:

  GlobalFitter( const TString&  name, const TString& title, double massH, double massW = 80.399 );

  void set_Wmass( double massW ) { massW_ = massW; };
  void set_Hmass( double massH ) { massH_ = massH; };

  double ErrEt_jet(double Et, double Eta);
  double ErrEta_jet(double Et, double Eta);
  double ErrPhi_jet(double Et, double Eta);
  double ErrEt_lept(double SumEt, double Eta);
  double ErrEta_lept(double SumEt, double Eta);
  double ErrPhi_lept(double SumEt, double Eta);
  double ErrEt_neut(double SumEt, MissingEnergy1 neut);
  double ErrEta_neut(double SumEt, MissingEnergy1 neut);
  double ErrPhi_neut(double SumEt, MissingEnergy1 neut);

  std::vector<TLorentzVector> fit( TLorentzVector jet1, TLorentzVector jet2, TLorentzVector Lept, MissingEnergy1 neut );

 private:

  double massW_;
  double massH_;
  TString name_;
  TString title_;

};

#endif
