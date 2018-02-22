#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include <fstream>



void gen_cosmic_tracks(
  TString outfile_str     = "input_tracks.txt",
  TString number_str     = "5",
  TString z_length_str     = "0.6",
  TString displacement_x_str = "0",
  TString displacement_y_str = "0"
){

  ofstream outfile;
  outfile.open (outfile_str);
  
// TH1F* muon_rate = new TH1F("muon_rate","muon_rate,theta (rad),intensity",360,0,
// Float_t mu_i=pow(cos(theta),2); // follows cos^2 law  

Int_t number = number_str.Atoi();
Float_t displacement_x = displacement_x_str.Atof();
Float_t displacement_y = displacement_y_str.Atof();
Float_t z_length = z_length_str.Atof();

Float_t x_min = -0.35;
Float_t x_max = 0.35;
Float_t y_min = -0.35;
Float_t y_max = 0.35;
Float_t z_min = -0.35;
Float_t z_max = 0.35;


//http://pdg.lbl.gov/2011/reviews/rpp2011-rev-cosmic-rays.pdf

// TF1* muon_rate_theta = new TF1("muon_rate","cos(x)*cos(x)",0,TMath::Pi()/2);
// TF1* muon_rate_theta = new TF1("muon_rate","cos(x)*cos(x)",0,TMath::TwoPi()*60/360); // limit at 60 deg
TF2 *muon_rate_theta_energy = new TF2("muon_rate_theta_energy",
           "0.14*pow(y,-2.7)*( 1/(1+1.1*y*cos(x)/115.0) + 0.054/(1+1.1*y*cos(x)/850.0) )"
           ,0,TMath::TwoPi()/4,0,1000);
TF1 *muon_rate_energy = new TF1("muon_rate_theta_energy",
           "0.14*pow(x,-2.7)*( 1/(1+1.1*x*cos(0)/115.0) + 0.054/(1+1.1*x*cos(0)/850.0) )"
           ,0,100);

gRandom->SetSeed(0);

// muon_rate_theta_energy->Draw("surf");
muon_rate_energy->Draw();

Double_t theta,ekin;
muon_rate_theta_energy->GetRandom2(theta,ekin);




/*
// new TCanvas();
// muon_rate->Draw();

gRandom->SetSeed(0);


for (Int_t i = 0 ; i< number; i++){
  Float_t theta = muon_rate_theta->GetRandom();
  
//   theta = TMath::Pi()/2-0.1; 
//   theta = 0;
  
  Float_t phi   = gRandom->Uniform(0,TMath::TwoPi());
  
//   phi = TMath::Pi();
//   phi = 0;
//   cout << "theta : " << theta << " phi : " << phi <<endl;
  
  Float_t length = z_length/TMath::Cos(theta);
  
  Float_t dir_x = TMath::Cos(phi)*TMath::Sin(theta);
  Float_t dir_y = TMath::Sin(phi)*TMath::Sin(theta);
  Float_t dir_z = TMath::Cos(theta);
//   cout << " dirs " << dir_x << " " << dir_y << " "<< dir_z << endl;
  
  Float_t x_0 = gRandom->Uniform(-displacement_x,displacement_x);
  Float_t y_0 = gRandom->Uniform(-displacement_y,displacement_y);
  Float_t z_0 = 0;
  
  Float_t length_a = length/2.;
  Float_t length_b = length/2.;
  
//   cout << "length a, length b: " << length_a << " " << length_b << endl;
  
  Float_t x_a = x_0 - length_a*dir_x;
  Float_t y_a = y_0 - length_a*dir_y;
  Float_t z_a = z_0 - length_a*dir_z;
  Float_t x_b = x_0 + length_b*dir_x;
  Float_t y_b = y_0 + length_b*dir_y;
  Float_t z_b = z_0 + length_b*dir_z;
  
//   cout << "a pos: " << x_a << " " << y_a << " " << z_a << endl;
//   cout << "b pos: " << x_b << " " << y_b << " " << z_b << endl;
  
  if(x_a > x_max){
    length_a = -(x_max - x_0)/dir_x;
  }
  if(x_a < x_min){
    length_a = -(x_min - x_0)/dir_x;
  }
  if(x_b > x_max){
    length_b = (x_max - x_0)/dir_x;
  }
  if(x_b < x_min){
    length_b = (x_min - x_0)/dir_x;
  }
  
  if(y_a > y_max){
    length_a = -(y_max - y_0)/dir_y;
  }
  if(y_a < y_min){
    length_a = -(y_min - y_0)/dir_y;
  }
  if(y_b > y_max){
    length_b = (y_max - y_0)/dir_y;
  }
  if(y_b < y_min){
    length_b = (y_min - y_0)/dir_y;
  }
  
  
  // correct length, if track goes out of cell
  x_a = x_0 - length_a*dir_x;
  y_a = y_0 - length_a*dir_y;
  z_a = z_0 - length_a*dir_z;
  x_b = x_0 + length_b*dir_x;
  y_b = y_0 + length_b*dir_y;
  z_b = z_0 + length_b*dir_z;
  
  
//   outfile << "area -0.40 -0.40 -0.40 0.40 0.40 0.40 view x=0 3d" << endl;
  outfile << "area -0.40 -0.40 -0.40 0.40 0.40 0.40" << endl;
  outfile << "track "<< x_a << " " << y_a << " " << z_a 
  << " " << x_b << " " << y_b << " " << z_b 
  << " HEED electron energy  0.70 GeV" << endl;
  outfile << "INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 1000" << endl;
  outfile << "DRIFT TRACK MONTE-CARLO-DRIFT LINE-PRINT" << endl;
  
}



// TH1F* theta_test = new TH1F("theta_test","theta_test",40,0,TMath::Pi()/2);
// 
// for (Int_t i = 0 ; i< 1000000; i++){
//   theta_test->Fill(muon_rate->GetRandom());
//   
// }
// 
// new TCanvas();
// theta_test->Draw();
  */
  
  
}