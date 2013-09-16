#include "GlobalFitter.h"

#include "TFitConstraintM.h"
#include "TFitConstraintMGaus.h"
#include "TFitConstraintEp.h"
#include "TFitParticleEtEtaPhi.h"
#include "TKinFitter.h"

#include <iostream>

void MissingEnergy1::SetSumEt(float SumEt){
  SumEt_=SumEt;
    };

void MissingEnergy1::SetNeutrino(TLorentzVector Neutrino){
  Neutrino_=Neutrino;
    };

GlobalFitter::GlobalFitter( const TString&  name, const TString& title, double massH, double massW ) : TKinFitter( name, title ) {
  name_ = name;
  title_ = title;
  massH_ = massH;
  massW_ = massW;
}


std::vector<TLorentzVector> GlobalFitter::fit(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector lept, MissingEnergy1 neut ) {

  this->reset();

  TMatrixD m_jet1(3,3);
  TMatrixD m_jet2(3,3);
  TMatrixD m_lept(3,3);
  TMatrixD m_neut(3,3);
  float p=-0.0185205; //correaltion

  m_jet1(0,0) = this->ErrEt_jet(jet1.Et(), jet1.Eta()); // et
  m_jet1(1,1) = this->ErrEta_jet(jet1.Et(), jet1.Eta()); // eta
  m_jet1(2,2) = this->ErrPhi_jet(jet1.Et(), jet1.Eta()); // phi
  m_jet2(0,0) = this->ErrEt_jet(jet2.Et(), jet2.Eta()); // et
  m_jet2(1,1) = this->ErrEta_jet(jet2.Et(), jet2.Eta()); // eta
  m_jet2(2,2) = this->ErrPhi_jet(jet2.Et(), jet2.Eta()); // phi
  m_lept(0,0) = this->ErrEt_lept(lept.Et(), lept.Eta()); // et
  m_lept(1,1) = this->ErrEta_lept(lept.Et(), lept.Eta()); // eta
  m_lept(2,2) = this->ErrPhi_lept(lept.Et(), lept.Eta()); // phi
  m_neut(0,0) = this->ErrEt_neut(neut.SumEt_, neut ); // et
  m_neut(1,1) = this->ErrEta_neut(neut.SumEt_, neut ); // eta
  m_neut(2,2) = this->ErrPhi_neut(neut.SumEt_, neut ); // phi
//  m_neut(0,2) = p*this->ErrEt_neut(neut.SumEt_, neut ); // correlazione
//  m_neut(2,0) = m_neut(0,2); // correlazione


  TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &jet1, &m_jet1 );
  TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &jet2, &m_jet2 );

  TFitParticleEtEtaPhi *fitLept = new TFitParticleEtEtaPhi( "Lept", "Lept", &lept, &m_lept );
  TFitParticleEtEtaPhi *fitNeut = new TFitParticleEtEtaPhi( "Neut", "Neut", &neut.Neutrino_, &m_neut );

  TFitConstraintMGaus *mCons_jets = new TFitConstraintMGaus( "WMassConstraint_jets", "WMassConstraintJets", 0, 0 , massW_, 2.085);
  mCons_jets->addParticles1( fitJet1, fitJet2 );

  TFitConstraintMGaus *mCons_lept = new TFitConstraintMGaus( "WMassConstraint_lept", "WMassConstraintLept", 0, 0 , massW_, 2.085);
  mCons_lept->addParticles1( fitLept, fitNeut );

  TFitConstraintMGaus *mCons_Tot = new TFitConstraintMGaus( "HMassConstraint_totGaus" , "HMassConstraintTotGaus", 0, 0, massH_, 5.);
  mCons_Tot->addParticles1( fitJet1, fitJet2, fitLept, fitNeut );
  //TFitConstraintM *mCons_Tot = new TFitConstraintM( "HMassConstraint_tot", "HMassConstraintTot", 0, 0 , massH_);
  //mCons_Tot->addParticles1( fitJet1, fitJet2, fitLept, fitNeut );

//  TFitConstraintEp::component Pt;
//  TFitConstraintEp *mCons_Etot = new TFitConstraintEp( "EConstraint_all" , "EConstraintAll", 0, Pt, 0.);
//  mCons_Etot->addParticles( fitJet1, fitJet2, fitLept, fitNeut );

  this->addMeasParticle( fitJet1 );
  this->addMeasParticle( fitJet2 );
  this->addMeasParticle( fitLept );
  this->addMeasParticle( fitNeut );
  this->addConstraint( mCons_jets );
  this->addConstraint( mCons_lept );
  this->addConstraint( mCons_Tot );

  //Set convergence criteria
  this->setMaxNbIter( 30 );
  this->setMaxDeltaS( 1e-2 );
  this->setMaxF( 1e-1 );
  this->setVerbosity(0);

  //Perform the fit
  TKinFitter::fit();

  TLorentzVector Jet1_kinfit(*fitJet1->getCurr4Vec());
  TLorentzVector Jet2_kinfit(*fitJet2->getCurr4Vec());
  TLorentzVector Lept_kinfit(*fitLept->getCurr4Vec());
  TLorentzVector Neut_kinfit(*fitNeut->getCurr4Vec());

  std::vector<TLorentzVector> AllFitted;
  AllFitted.push_back( Jet1_kinfit );
  AllFitted.push_back( Jet2_kinfit );
  AllFitted.push_back( Lept_kinfit );
  AllFitted.push_back( Neut_kinfit );
 
 if(AllFitted.size() > 4){std::cout<<"Error in the size of TLorentzVector fitted"<<std::endl; exit(13);}

  return AllFitted;

}



// lepton resolutions. Setted low just for now
double GlobalFitter::ErrEt_lept( double Et, double Eta) {

  double InvPerr2 = 0.01*0.01;
//std::cout<<InvPerr2<<std::endl;
  return InvPerr2;
}

// Lepton resolution. Setted low just for now
double GlobalFitter::ErrEta_lept( double Et, double Eta) {

  double InvPerr2 = 0.01*0.01;
//std::cout<<InvPerr2<<std::endl;
  return InvPerr2;
}

// Lepton resolution. Setted low just for now
double GlobalFitter::ErrPhi_lept( double Et, double Eta) {

  double InvPerr2 = 0.01*0.01;
//std::cout<<InvPerr2<<std::endl;
  return InvPerr2;
}

// Neutrino resolution.
double GlobalFitter::ErrEt_neut( double SumEt, MissingEnergy1 neut ) {

  double InvPerr2 = pow( (0.5)*sqrt(SumEt)*(neut.Neutrino_.Px()+neut.Neutrino_.Py())/(2*neut.Neutrino_.Pt()), 2); // 0.5*sqrt(SumEt)*0.5*sqrt(SumEt);
//std::cout<<"Et_sqrt(InvPer): "<<sqrt(InvPerr2)<<std::endl;
  return InvPerr2;

}

// Neutrino resolution
double GlobalFitter::ErrEta_neut( double SumEt, MissingEnergy1 neut ) {

  double InvPerr2=1.;

//std::cout<<"Eta_sqrt(InvPer): "<<sqrt(InvPerr2)<<std::endl;
  return InvPerr2;

}

//Neutrino resolution
double GlobalFitter::ErrPhi_neut( double SumEt, MissingEnergy1 neut ) {

  double InvPerr2;

  InvPerr2 =  pow( ((0.6)*1-(neut.Neutrino_.Py()/neut.Neutrino_.Px()))*(sqrt(SumEt)*0.5)/(neut.Neutrino_.Px()*(1+(pow(neut.Neutrino_.Py(),2)/pow(neut.Neutrino_.Px(),2)))) ,2);

//std::cout<<"Phi_sqrt(InvPer): "<<sqrt(InvPerr2)<<std::endl;

 return InvPerr2;

}

// pfjet resolutions. taken from AN-2010-371
double GlobalFitter::ErrEt_jet( double Et, double Eta) {

  double InvPerr2;

  double N, S, C, m;
  if(fabs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( fabs(Eta) < 1. ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( fabs(Eta) < 1.5 ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( fabs(Eta) < 2. ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( fabs(Eta) < 2.5 ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
  } else if( fabs(Eta) < 3. ) {
    N = -3.33814;
    S = 0.73360;
    C = 0.;
    m = 0.08264;
  } else if( fabs(Eta) < 5. ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }
  // this is the absolute resolution (squared), not sigma(pt)/pt
  // so have to multiply by pt^2, thats why m+1 instead of m-1
  InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;

 //std::cout<<"Et Jet: "<<InvPerr2<<std::endl;
  return InvPerr2*1.2;

}

//pfjet position resolutions taken from AN 2010-104
double GlobalFitter::ErrEta_jet( double Et, double Eta) {

  double InvPerr2;

  double a,b,c,d,e;
  if( fabs(Eta) < 0.5 ) {
    a = 487.702;
    b = 1.72868;
    c = 0.0297405;
    d = -487.197;
    e = -1.16389;
  } else if( fabs(Eta) < 1.0 ) {
    a = 277.114;
    b = 1.31746;
    c = 0.0232343;
    d = -276.588;
    e = -1.07289;
  } else if( fabs(Eta) < 1.5 ) {
    a = 19.7603;
    b = 0.406775;
    c = 0.0056006;
    d = -19.1144;
    e = -1.24005;
  } else if( fabs(Eta) < 2.0 ) {
    a = 41.55;
    b = 0.556349;
    c = 0.0094941;
    d = -40.8018;
    e = -1.39179;
  } else if( fabs(Eta) < 2.5 ) {
    a = 0.833363;
    b = 0.0786743;
    c = 0.0036199;
    d = 0.0507317;
    e = -1.5492;
  } else if( fabs(Eta) < 3. ) {
    a = 3.4712;
    b = 0.383594;
    c = 2.6831e-7;
    d = -2.93791;
    e = -0.259687;
  } else {
    std::cout << "Not implemented Eta Err for eta > 3. Exiting." << std::endl;
    exit(1123);
  }

  InvPerr2 = sqrt(a*a/(Et * Et) + b*b/Et + c*c) + d/Et + e/pow(Et,1.5);
  InvPerr2 *= InvPerr2;
 //std::cout<<"Eta Jet: "<<InvPerr2<<std::endl;
  return InvPerr2*1.2;

}


double GlobalFitter::ErrPhi_jet( double Et, double Eta) {

  double InvPerr2;

  double a,b,c,d,e;
  if( fabs(Eta) < 0.5 ) {
    a = 926.978;
    b = 2.52747;
    c = 0.0304001;
    d = -926.224;
    e = -1.94117;
  } else if( fabs(Eta) < 1.0 ) {
    a = 3.3251e-6;
    b = 0.063941;
    c = 0.0038759;
    d = 0.301932;
    e = -0.825352;
  } else if( fabs(Eta) < 1.5 ) {
    a = 0.38469;
    b = 0.0755727;
    c = 0.0044353;
    d = 0.453887;
    e = -1.8947;
  } else if( fabs(Eta) < 2.0 ) {
    a = 2.9200e-7;
    b = 0.0718389;
    c = 0.0038557;
    d = 0.403668;
    e = -0.62698;
  } else if( fabs(Eta) < 2.5 ) {
    a = 0.0033663;
    b = 0.0880209;
    c = 0.0023084;
    d = 0.214304;
    e = -0.416353;
  } else if( fabs(Eta) < 3. ) {
    a = 11.1957;
    b = 0.643236;
    c = 0.0071142;
    d = -10.7613;
    e = 0.280927;
  } else {
    std::cout << "Not implemented Phi Err for eta > 3. Exiting." << std::endl;
    exit(1123);
  }

  InvPerr2 = sqrt(a*a/(Et * Et) + b*b/Et + c*c) + d/Et + e/pow(Et,1.5);
  InvPerr2 *= InvPerr2;
 //std::cout<<"Phi Jet: "<<InvPerr2<<std::endl;

  return InvPerr2*1.2;

}



