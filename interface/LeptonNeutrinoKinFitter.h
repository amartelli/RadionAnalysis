// ------------------------------------------------
// 
// LeptonNeutrinoKinFitter - kinematic fit for X->lv
//
// ------------------------------------------------

#ifndef LeptonNeutrinoKinFitter_h
#define LeptonNeutrinoKinFitter_h

#include "TKinFitter.h"
#include "TLorentzVector.h"


class MissingEnergy{

 public:

  void SetuncorrEnergyAK5Jet(float uncorrEnergyAK5Jet);
  void SetNeutrino(TLorentzVector Neutrino);

  float uncorrEnergyAK5Jet_;
  TLorentzVector Neutrino_;

};



class LeptonNeutrinoKinFitter : public TKinFitter {

 public:
  
  LeptonNeutrinoKinFitter( const TString&  name, const TString& title, double mass = 80.399 );

  void set_mass( double mass ) { mass_ = mass; };

  double ErrEt_lept(double SumEt, double Eta); 
  double ErrEta_lept(double SumEt, double Eta); 
  double ErrPhi_lept(double SumEt, double Eta); 
  double ErrEt_neut(double SumEt, double Eta); 
  double ErrEta_neut(double SumEt, double Eta); 
  double ErrPhi_neut(double SumEt, double Eta); 

  std::pair<TLorentzVector,TLorentzVector> fit( TLorentzVector jet1, MissingEnergy neut );

 private:

  double mass_;
  TString name_;
  TString title_;

};


#endif
