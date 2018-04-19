#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include <fstream>


TString from_env(TString env_var,TString default_val){
  if(gSystem->Getenv(env_var)){
    return gSystem->Getenv(env_var);
  } 
  return default_val;
}

Float_t from_env_float(TString env_var,TString default_val){
  TString val = default_val;
  if(gSystem->Getenv(env_var)){
    val = gSystem->Getenv(env_var);
  } 
  return val.Atof();
}

Int_t from_env_int(TString env_var,TString default_val){
  TString val = default_val;
  if(gSystem->Getenv(env_var)){
    val = gSystem->Getenv(env_var);
  } 
  return val.Atoi();
}



void gen_UV_charge(){

  

TString outfile_str = from_env("outfile","input_tracks.txt");  
  ofstream outfile;
  outfile.open (outfile_str);
  
// TH1F* muon_rate = new TH1F("muon_rate","muon_rate,theta (rad),intensity",360,0,
// Float_t mu_i=pow(cos(theta),2); // follows cos^2 law  

Int_t events = from_env_int("events","1");
Int_t charges = from_env_int("charges","1");
Float_t x = from_env_float("x","0");
Float_t y = from_env_float("y","0");
Float_t z = from_env_float("z","0");

Float_t n_photons = from_env_float("n_photons","2");

Float_t sigma_x = from_env_float("sigma_x","0.005"); // mm
Float_t sigma_y = from_env_float("sigma_y","0.005"); // mm
Float_t sigma_z = from_env_float("sigma_z","0.2045"); // mm // rayleigh length in mm
Float_t rayleigh = from_env_float("rayleigh","0.2045"); // mm // rayleigh length in mm


  TF1* charge_density_x = new TF1("charge_density_x", "gaus",-.5,.5);
  charge_density_x->SetParameter(0,1/(TMath::Sqrt(TMath::TwoPi())*sigma_x) );
  charge_density_x->SetParameter(1,0);
  charge_density_x->SetNpx(1000);
  charge_density_x->SetParameter(2,sigma_x);
  
  TF1* charge_density_y = new TF1("charge_density_y", "gaus",-.5,.5);
  charge_density_y->SetParameter(0,1/(TMath::Sqrt(TMath::TwoPi())*sigma_y) );
  charge_density_y->SetParameter(1,0);
  charge_density_y->SetNpx(1000);
  charge_density_y->SetParameter(2,sigma_y);
  
  
  // xingming distribution
  TF1* charge_density_z = new TF1("charge_density_z", "[0]*TMath::Power(1/(1+ (x/[1])*(x/[1]) ),[2]-1)",-3,3);
  charge_density_z->SetParameter(0,1); // just scaling constant
  charge_density_z->SetParameter(1,rayleigh); // rayleigh length in mm
  charge_density_z->SetParameter(2,n_photons); // n for n-photon absorbtion 
  charge_density_z->SetNpx(1000);
  
//   TF1* charge_density_z_gaus = new TF1("charge_density_z_gaus", "gaus",-3,3);
//   charge_density_z_gaus->SetParameter(0,1/(TMath::Sqrt(TMath::TwoPi())*sigma_z) );
//   charge_density_z_gaus->SetParameter(0,1);
//   charge_density_z_gaus->SetParameter(1,0);
//   charge_density_z_gaus->SetParameter(2,sigma_z);
//   charge_density_z_gaus->SetLineColor(1);
//   charge_density_z_gaus->SetNpx(1000);
  
//   TF1* charge_density_z = new TF1("charge_density_z", "gaus",-3,3);
//   charge_density_z->SetParameter(0,1/(TMath::Sqrt(TMath::TwoPi())*sigma_z) );
//   charge_density_z->SetParameter(1,0);
//   charge_density_z->SetParameter(2,sigma_z);


gRandom->SetSeed(0);

// TH1F* dummy = new TH1F("dummy","dummy",100,-5,5);

for (Int_t i = 0 ; i< events; i++){
  outfile << "Say \"new UV event\"" << endl;
  outfile << "Say \"laser position: "<< x << " " << y << " " << z << "\"" << endl;
  for (Int_t j = 0 ; j< charges; j++){
    
    Float_t charge_x = (charge_density_x->GetRandom()+x)/10.; // in cm  
    Float_t charge_y = (charge_density_y->GetRandom()+y)/10.;
    Float_t charge_z = (charge_density_z->GetRandom()+z)/10.;
  //   Float_t charge_x = 0;
  //   Float_t charge_y = 0;
  //   Float_t charge_z = 0;
//     dummy->Fill(charge_density_z->GetRandom());
    outfile << "area -0.40 -0.40 -0.40 0.40 0.40 0.40" << endl;
    outfile << Form("track %f %f %f %f %f %f lines 1", charge_x, charge_y, charge_z-0.0001, charge_x, charge_y, charge_z+0.0001)<< endl;
    outfile << "INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 1000" << endl;
    outfile << "DRIFT TRACK MONTE-CARLO-DRIFT LINE-PRINT NOLINE-PLOT" << endl;
    
  }
}

// new TCanvas();
// // dummy->Draw();
// charge_density_z_gaus->Draw();
// charge_density_z_gaus->GetXaxis()->SetTitle("z pos (mm)");
// charge_density_z_gaus->GetYaxis()->SetTitle("Ionization density (a.u.)");
// charge_density_z->DrawClone("same");
//   charge_density_z->SetParameter(2,3); // n for n-photon absorbtion 
// 
// charge_density_z->DrawClone("same");
//   charge_density_z->SetParameter(2,4); // n for n-photon absorbtion 
// 
// charge_density_z->DrawClone("same");


  
  
  
}