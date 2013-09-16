#include "LeptonNeutrinoKinFitter.h"

#include "TFitConstraintM.h"
#include "TFitParticleEtEtaPhi.h"
#include "TKinFitter.h"

#include <iostream>

void MissingEnergy::SetuncorrEnergyAK5Jet(float uncorrEnergyAK5Jet){
  
  uncorrEnergyAK5Jet_=uncorrEnergyAK5Jet;
    
    };
void MissingEnergy::SetNeutrino(TLorentzVector Neutrino){
  
  Neutrino_=Neutrino;
    
    };

LeptonNeutrinoKinFitter::LeptonNeutrinoKinFitter( const TString&  name, const TString& title, double mass ) : TKinFitter( name, title ) {
  
  name_ = name;
  title_ = title;
  mass_ = mass;

}



std::pair<TLorentzVector,TLorentzVector> LeptonNeutrinoKinFitter::fit(TLorentzVector lept, MissingEnergy neut ) {

  this->reset();

  TMatrixD m_lept(3,3);
  TMatrixD m_neut(3,3);

  m_lept(0,0) = this->ErrEt_lept(lept.Et(), lept.Eta()); // et
  m_lept(1,1) = this->ErrEta_lept(lept.Et(), lept.Eta()); // eta
  m_lept(2,2) = this->ErrPhi_lept(lept.Et(), lept.Eta()); // phi
  m_neut(0,0) = this->ErrEt_neut(neut.uncorrEnergyAK5Jet_, neut.Neutrino_.Eta()); // et
  m_neut(1,1) = this->ErrEta_neut(neut.uncorrEnergyAK5Jet_, neut.Neutrino_.Eta()); // eta
  m_neut(2,2) = this->ErrPhi_neut(neut.uncorrEnergyAK5Jet_, neut.Neutrino_.Eta()); // phi

  TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &lept, &m_lept );
  TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &neut.Neutrino_, &m_neut );
  
  TFitConstraintM *mCons_jets = new TFitConstraintM( "ZMassConstraint_jets", "ZMass-Constraint", 0, 0 , mass_);
  mCons_jets->addParticles1( fitJet1, fitJet2 );

  
  this->addMeasParticle( fitJet1 );
  this->addMeasParticle( fitJet2 );
  this->addConstraint( mCons_jets );
    

  //Set convergence criteria
  this->setMaxNbIter( 30 );
  this->setMaxDeltaS( 1e-2 );
  this->setMaxF( 1e-1 );
  this->setVerbosity(0);

  //Perform the fit
  TKinFitter::fit();

  TLorentzVector jet1_kinfit(*fitJet1->getCurr4Vec());
  TLorentzVector jet2_kinfit(*fitJet2->getCurr4Vec());

  std::pair<TLorentzVector,TLorentzVector> jetpair;
  jetpair.first  = jet1_kinfit;
  jetpair.second = jet2_kinfit;

  return jetpair;

}



// lepton resolutions. Setted low just for now
double LeptonNeutrinoKinFitter::ErrEt_lept( double Et, double Eta) {
  
  double InvPerr2 = 0.01*0.01;

  return InvPerr2;

}

// Lepton resolution. Setted low just for now
double LeptonNeutrinoKinFitter::ErrEta_lept( double Et, double Eta) {

  double InvPerr2 = 0.01*0.01;
 
  return InvPerr2;

}

// Lepton resolution. Setted low just for now
double LeptonNeutrinoKinFitter::ErrPhi_lept( double Et, double Eta) {

  double InvPerr2 = 0.01*0.01;
 
  return InvPerr2;

}

// Neutrino resolution.
double LeptonNeutrinoKinFitter::ErrEt_neut( double SumEt, double Eta) {
  
  double InvPerr2 = 0.5*0.5*SumEt*SumEt;

  return InvPerr2;

}

// Neutrino resolution
double LeptonNeutrinoKinFitter::ErrEta_neut( double SumEt, double Eta) {

  double InvPerr2=1.;

  return InvPerr2;

}

//Neutrino resolution
double LeptonNeutrinoKinFitter::ErrPhi_neut( double SumEt, double Eta) {
  
  double InvPerr2;
  
  double a,b,c,d,e;

    a = 926.978;
    b = 2.52747;
    c = 0.0304001;
    d = -926.224;
    e = -1.94117;
 
  InvPerr2 = 2*sqrt(a*a/(SumEt * SumEt) + b*b/SumEt + c*c) + d/SumEt + e/pow(SumEt,1.5);//la raddoppio
  InvPerr2 *= InvPerr2;

  return InvPerr2;

}
