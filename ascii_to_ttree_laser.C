
#include "TH1.h"
#include "TF1.h"
#include "TArray.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

// for compatibility for analysis scripts from 2016 ...
void hist_to_tarrayf(TH1* hist, TArrayF* xarr, TArrayF* yarr){
  
  Float_t x,y;
  
  for (Int_t i = 1 ; i <= hist->GetEntries(); i++){
  
    x = hist->GetXaxis()->GetBinCenter(i);
    y = hist->GetBinContent(i);
    xarr->AddAt(x,i-1);
    yarr->AddAt(y,i-1);
  
  }
  xarr->Set(hist->GetEntries());
  yarr->Set(hist->GetEntries());

  
}

TH1* fftconvolve(TH1D* h1, TH1D* h2){
  
  
  Int_t samples = h1->GetEntries();
  
  TVirtualFFT *fft_h1 = TVirtualFFT::FFT(1, &samples, "R2C M K");
  fft_h1->SetPoints(h1->GetArray());
  fft_h1->Transform();
  
  TVirtualFFT *fft_h2 = TVirtualFFT::FFT(1, &samples, "R2C M K");
  fft_h2->SetPoints(h2->GetArray());
  fft_h2->Transform();
 
  Double_t *re_h1 = new Double_t[samples];
  Double_t *im_h1 = new Double_t[samples];
  fft_h1->GetPointsComplex(re_h1,im_h1);
  
  Double_t *re_h2 = new Double_t[samples];
  Double_t *im_h2 = new Double_t[samples];
  fft_h2->GetPointsComplex(re_h2,im_h2);
 
  
  Double_t *re_conv = new Double_t[samples];
  Double_t *im_conv = new Double_t[samples];
 
  for (Int_t i = 0; i<samples; i++){
    // complex multiplication ... element wise
    re_conv[i] = re_h1[i]*re_h2[i] - im_h1[i]*im_h2[i];
    im_conv[i] = im_h1[i]*re_h2[i] + re_h1[i]*im_h2[i];
  }
  
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &samples, "C2R M K");
  fft_back->SetPointsComplex(re_conv,im_conv);
  fft_back->Transform();
  
  TH1 *hback = 0;
  hback = TH1::TransformHisto(fft_back,hback, "RE");
  hback->GetXaxis()->SetLimits(
    h1->GetXaxis()->GetXmin(),
    h1->GetXaxis()->GetXmax()
  );
  delete fft_h1;
  delete fft_h2;
  delete fft_back;
  
  return hback;
}






void ascii_to_ttree_laser(TString infile) {
  
  Bool_t draw_pulses = true;
  
  // for fish
  Int_t nth_electron = 0;
  Int_t elno_max = 20;
  
  TFile* f_out = new TFile(Form("%s%s",gSystem->DirName(infile),"/f_out.root"),"RECREATE");
  
  TFile* infile_root = new TFile(infile+".root");
  
  
  
  
  
  TTree* garfield_tree;
  
  if(infile_root->IsZombie()){ // does not exist yet
    
    garfield_tree = new TTree("garfield_tree", "garfield_tree");
    
    cout << "convert ascii data ..." << endl;
    garfield_tree->ReadFile(infile, "laser_x:laser_y:laser_z:n:x:y:z:t_drift:wire");
    cout << "done" << endl;
    
    TFile* tree_out = new TFile(infile+".root","RECREATE");
    garfield_tree->Write();
    tree_out->Close();
    
    
  } else {
    cout << "found already converted root tree in " << infile+".root" << endl;
  }
  
  infile_root = new TFile(infile+".root");
  garfield_tree = (TTree*) infile_root->Get("garfield_tree");
  f_out->cd();
  
//   garfield_tree->Draw("Energy_1:Energy_2”, "”, "");
//   garfield_tree->Draw("Energy_1”);
  new TCanvas();
  garfield_tree->Draw("t_drift >> tdrift_h()");
  new TCanvas();
//   garfield_tree->Draw("t_drift:x >> tdrift_vs_x()","","colz");
  garfield_tree->Draw("t_drift*1000:x*10 >> tdrift_vs_x(602,-3,3,256,0,102.3)","","colz");
  TH2F* tdrift_vs_x = (TH2F*) f_out->Get("tdrift_vs_x");
  tdrift_vs_x->GetXaxis()->SetTitle("x pos (mm)");
  tdrift_vs_x->GetYaxis()->SetTitle("drift time (ns)");
  tdrift_vs_x->SetTitle("drift time vs x start position");
  tdrift_vs_x->Draw("colz");
// nice binning:
//    garfield_tree->Draw("t_drift:x >> tdrift_vs_x(200,-0.3,0.3,256,0,0.1023)","","colz");
  
  new TCanvas();
  // nice bins!
  garfield_tree->Draw("t_drift*1000:z*10 >> tdrift_vs_z(602,-3,3,256,0,102.3)","","colz");
  TH2F* tdrift_vs_z = (TH2F*) f_out->Get("tdrift_vs_z");
  tdrift_vs_z->GetXaxis()->SetTitle("z pos (mm)");
  tdrift_vs_z->GetYaxis()->SetTitle("drift time (ns)");
  tdrift_vs_z->SetTitle("drift time vs z start position");
  tdrift_vs_z->Draw("colz");
  
  
//   garfield_tree->Draw("t_drift:z >> tdrift_vs_z()","","colz");
//   garfield_tree->Draw("t_drift:z >> tdrift_vs_z(200,-0.3,0.3,200,-0.02,0.1)","","colz");
// nice binning:
//  garfield_tree->Draw("t_drift:z >> tdrift_vs_z(602,-0.3,0.3,256,0,0.1023)","","colz");
  

//   new TBrowser();
  new TCanvas();
//   garfield_tree->Draw("t_drift:y >> tdrift_vs_y()","","colz");
  garfield_tree->Draw("t_drift*1000:y*10 >> tdrift_vs_y(602,-3,3,256,0,102.3)","","colz");
//   garfield_tree->Draw("t_drift:y >> tdrift_vs_y(200,-0.3,0.3,200,-0.02,0.1)","","colz");
//   nice binning
//   garfield_tree->Draw("t_drift:y >> tdrift_vs_y(602,-0.3,0.3,256,0,0.1023)","","colz");
  TH2F* tdrift_vs_y = (TH2F*) f_out->Get("tdrift_vs_y");
  tdrift_vs_y->GetXaxis()->SetTitle("y pos (mm)");
  tdrift_vs_y->GetYaxis()->SetTitle("drift time (ns)");
  tdrift_vs_y->SetTitle("drift time vs y start position");
  tdrift_vs_y->Draw("colz");
  
  
  /*
  
  /// load the kernel/chamber_IR function
  
  TGraph* tg_kern = new TGraph("chamber_IR.csv","%lg, %lg");
  
//   Float_t IR_y_scaler = 1./230./100.;
//   Float_t IR_y_scaler = 1./230; // divide through the number of charges of one Fe55 pulse
  Float_t IR_y_scaler = 1./100; // divide through the number of charges of one Fe55 pulse
//   Float_t IR_y_scaler = 1.;
  */
  Float_t sample_width = 1.6e-6/10;
  Int_t samples = 3200/10;
  
  /*
  
  TH1D* th_kern = new TH1D("th_kern","th_kern;t(ns)",samples,0,sample_width);
  
  for( Int_t i = 1; i <= samples; i++) {
    Float_t t = sample_width * (Float_t) i /(Float_t) samples;
    th_kern->SetBinContent(i,tg_kern->Eval(t-0.1e-6)*IR_y_scaler);
  }
  new TCanvas();
  tg_kern->Draw();
  
  new TCanvas();
  th_kern->Draw();
  
  // tree output compatible with earlier analysis scripts ... a bit bulky
  
  
  TH1D* th_fake_pmt = new TH1D("th_fake_pmt","th_fake_pmt;t(ns)",samples,0,sample_width);
  
  for( Int_t i = 1; i <= samples; i++) {
    Float_t t = sample_width * (Float_t) i /(Float_t) samples;
    th_fake_pmt->SetBinContent(i,tg_kern->Eval(t-0.1e-6)*600);
  }
  
  TTree* pulse_mem = new TTree("pulse_mem","Acquired single pulse samples");
 
  TArrayF* signal_xarr = new TArrayF(samples);
  TArrayF* signal_yarr = new TArrayF(samples);
  pulse_mem->Branch("signal_x",signal_xarr);
  pulse_mem->Branch("signal_y",signal_yarr);
  TArrayF* trigger_yarr = new TArrayF(samples);
  TArrayF* trigger_xarr = new TArrayF(samples);
  
  hist_to_tarrayf(th_fake_pmt,trigger_xarr,trigger_yarr);
  pulse_mem->Branch("trigger_x",trigger_xarr);
  pulse_mem->Branch("trigger_y",trigger_yarr);
  
  
  
  
  // process the garfield tracks
 */ 
  
  
  
  TH1D* th_esig = new TH1D("th_esig","th_esig;t(s)",samples,0,sample_width);
  TH1*  th_esig_cum = 0;
 
//   TH1* th_conv = 0;
  
//   Float_t weight = 1./((Float_t) samples ); // somehow I need that factor, so the convolution preserves the absolute Y values, because electron signals are no dirac peaks
  
  Float_t t_drift;
  Float_t t_drift_first;
  Float_t t_drift_second;
  Float_t t_drift_third;
  Float_t t_drift_fourth;
  Float_t n,x,y,z,last_x,last_y,last_z,wire;
  Float_t last_n = 1;
  garfield_tree->SetBranchAddress("t_drift",&t_drift);
  garfield_tree->SetBranchAddress("n",&n);
  garfield_tree->SetBranchAddress("laser_x",&x);
  garfield_tree->SetBranchAddress("laser_y",&y);
  garfield_tree->SetBranchAddress("laser_z",&z);
  garfield_tree->SetBranchAddress("wire",&wire);
 
//   new TCanvas();
//   th_kern->Draw();
  
  
  new TCanvas();
  Float_t t_drift_a = 1000;
  Float_t t_drift_b = 1000;
  Int_t   elno = 0;
  TTree* fish_tree = new TTree("fish_tree","fish_tree");
  fish_tree->Branch("t_drift_a",&t_drift_a);
  fish_tree->Branch("t_drift_b",&t_drift_b);
  fish_tree->Branch("x",&last_x);
  fish_tree->Branch("y",&last_y);
  fish_tree->Branch("z",&last_z);
  fish_tree->Branch("elno",&elno);
  
  std::vector<Float_t> t_drift_a_vec;
  std::vector<Float_t> t_drift_b_vec;
  
  
  
  // remember the y laser positions
  std::map<Float_t,Int_t> laser_x_positions;
  std::map<Float_t,Int_t> laser_ypositions;
  std::map<Float_t,Int_t> laser_z_positions;
  
  std::map<Float_t, std::map<Float_t,  std::map<Float_t,Int_t>>> laser_positions;
  
  
  
  Int_t primaries = garfield_tree->GetEntries();
//   primaries = 1;
  for (Int_t i = 0 ; i < primaries + 1; i++){
    
    if(i < primaries){
      garfield_tree->GetEntry(i);
    } else {
      n++; // to trigger last processing
    }
    
    if (wire == 1){
      t_drift_a_vec.push_back(t_drift*1000);
    } else if (wire == 2){
      t_drift_b_vec.push_back(t_drift*1000);
    }
    
    
//     cout << "n: " << n << endl;
    if (n > last_n){
      laser_x_positions[last_x]++;
      laser_ypositions[last_y]++;
      laser_z_positions[last_z]++;
      laser_positions[last_x][last_y][last_z]++;
//       cout << "new N! " << endl;
//       new TCanvas();
/*      
      delete th_conv;
      th_conv = fftconvolve(th_kern,th_esig);
      
      th_conv->GetXaxis()->SetRangeUser(-0.1e-6,0.5e-6);
      th_conv->GetYaxis()->SetRangeUser(-2e-3,0.5e-3); */
//       th_esig->DrawClone();

      t_drift_a = 1000;
      t_drift_b = 1000;
      std::sort(t_drift_a_vec.begin(),t_drift_a_vec.end());
      std::sort(t_drift_b_vec.begin(),t_drift_b_vec.end());
      
      for(elno = 0; elno < elno_max; ++elno){
        if(t_drift_a_vec.size() > elno){
          t_drift_a = t_drift_a_vec[elno];
        }
        if(t_drift_b_vec.size() > elno){
          t_drift_b = t_drift_b_vec[elno];
        }
        fish_tree->Fill();
      }
      
      
      t_drift_a_vec.clear();
      t_drift_b_vec.clear();
      
      if(draw_pulses && n < 100){
        if(last_n == 1){
//           th_conv->DrawClone();
          th_esig->DrawClone();
        } else {
//           th_conv->DrawClone("same");
//           th_esig->DrawClone("same");
        }
      }
//       gROOT->cd();
//       th_esig_cum = th_esig->GetCumulative();
//       th_conv->SetName(Form("%d pulse",i));
// //       th_conv->Write();
//       hist_to_tarrayf(th_conv,signal_xarr,signal_yarr);
//       pulse_mem->Fill();
      th_esig->Reset();
    }
    
    th_esig->Fill(t_drift*1e-6);
//     th_esig->Fill(0.01*1e-6,weight);
    
    last_n = n;
    last_x = x;
    last_y = y;
    last_z = z;
  }
//   pulse_mem->Write();

// new TCanvas();
// th2_first_e->Draw("colz");

// cout << "all laser x positions: " << endl;
// for(auto const& dummy: laser_x_positions){
//   cout << "  " << dummy.first << endl;
// }
// 
// cout << "all laser y positions: " << endl;
// for(auto const& dummy: laser_ypositions){
//   cout << "  " << dummy.first << endl;
// }
// 
// cout << "all laser z positions: " << endl;
// for(auto const& dummy: laser_z_positions){
//   cout << "  " << dummy.first << endl;
// }

/*
cout << "all laser positions: " << endl;
for(auto const& xmap: laser_positions){
//   cout << "  " << xmap.first << endl;
  for(auto const& ymap: xmap.second){
//     cout << "  " << ymap.first << endl;
    for(auto const& zmap: ymap.second){
      cout << "x " << xmap.first << " y " << ymap.first << " z " << zmap.first << endl;
      
    }
  }
}
*/


f_out->cd();

  TMultiGraph * mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  
  
Int_t graphno = 0;



  Double_t t1=0;   
  Double_t t1_mean=0;   
  Double_t counts=1;    
  Double_t ref_counts=1;
  Double_t efficiency=1;
  Double_t t1_std=0;    
  Double_t tot_mean=0;  
  Double_t tot_std=0;  
  Double_t xpos=0;
  Double_t ypos=0;
  Double_t zpos=0;
  Double_t polar_y=0;
  Double_t polar_z=0;
  Double_t radius=0;
  Double_t phi=0;
  Double_t pol_cent_x=0;
  Double_t pol_cent_y=0;
  Double_t pol_cent_z=0;
  
  TTree* scan_data_tree = new TTree("scan_data_tree","scan_data_tree");
  
  scan_data_tree->Branch("t1",&t1);
  scan_data_tree->Branch("t1_mean",&t1_mean);
  scan_data_tree->Branch("counts",&counts);
  scan_data_tree->Branch("ref_counts",&ref_counts);
  scan_data_tree->Branch("t1_std",&t1_std);
  scan_data_tree->Branch("tot_means",&tot_mean);
  scan_data_tree->Branch("tot_std",&tot_std);
  scan_data_tree->Branch("xpos",&xpos);
  scan_data_tree->Branch("ypos",&ypos);
  scan_data_tree->Branch("zpos",&zpos);
  scan_data_tree->Branch("efficiency",&efficiency);
  scan_data_tree->Branch("radius",&radius);
  scan_data_tree->Branch("phi",&phi);
  scan_data_tree->Branch("polar_y",&polar_y);
  scan_data_tree->Branch("polar_z",&polar_z);
  scan_data_tree->Branch("pol_cent_x",&pol_cent_x);
  scan_data_tree->Branch("pol_cent_y",&pol_cent_y);
  scan_data_tree->Branch("pol_cent_z",&pol_cent_z);



for(Int_t my_elno = 0; my_elno < 1; my_elno+=2){
  graphno++;

TGraphErrors* tge_drift_time = new TGraphErrors();
tge_drift_time->SetName(Form("tge_drift_time_elno_%d",my_elno));
tge_drift_time->SetTitle(Form("Simulated drift time, electron no %d",my_elno));
tge_drift_time->GetXaxis()->SetTitle("y pos (um)");
tge_drift_time->GetYaxis()->SetTitle("drift time (ns)");

tge_drift_time->SetLineColor(graphno);

TGraphErrors* tge_drift_time_uncert = new TGraphErrors();
tge_drift_time_uncert->SetName(Form("tge_drift_time_uncert_elno_%d",my_elno));
tge_drift_time_uncert->SetTitle(Form("Simulated drift time uncertainty, electron no %d",my_elno));
tge_drift_time_uncert->GetXaxis()->SetTitle("y pos (um)");
tge_drift_time_uncert->GetYaxis()->SetTitle("drift time uncertainty (ns)");

Int_t tgepoint = 0;
  
  
for(auto const& xmap: laser_positions){
//   cout << "  " << xmap.first << endl;
  for(auto const& ymap: xmap.second){
//     cout << "  " << ymap.first << endl;
    for(auto const& zmap: ymap.second){
      cout << "x " << xmap.first << " y " << ymap.first << " z " << zmap.first << endl;
      
  xpos = xmap.first * 1000.0;
  ypos = ymap.first * 1000.0;
  zpos = zmap.first * 1000.0;
  
//   cout << "  " << dummy.first << endl;
//   gROOT->cd();
  f_out->cd();
  fish_tree->Draw(Form("t_drift_a >> %05.0f_%05.0f_%05.0f_t1_hist",xpos,ypos,zpos),Form("abs(%f - x*1000)<1 && abs(%f - y*1000)<1 && abs(%f - z*1000)<1 && elno ==%d",xpos,ypos,zpos,my_elno),"");
//   TH1F* dummy_hist = (TH1F*) gROOT->Get(Form("%05.0f_%05.0f_%05.0f_t1_hist",xpos,ypos,zpos));
  TH1F* dummy_hist = (TH1F*) f_out->Get(Form("%05.0f_%05.0f_%05.0f_t1_hist",xpos,ypos,zpos));
  Float_t drift_time = dummy_hist->GetMean();
  Float_t drift_time_stdev = dummy_hist->GetStdDev();
  tge_drift_time->SetPoint(tgepoint,ypos,drift_time);
  tge_drift_time->SetPointError(tgepoint,0,drift_time_stdev);
  tge_drift_time_uncert->SetPoint(tgepoint,ypos,drift_time_stdev);
  tgepoint++;
//   cout << "y pos: " << ypos << " drift time: " << drift_time << endl;
  
  
  /// this is for exporting the data for comparison with rossendorf measurements
  
  
    polar_y = ypos - pol_cent_y;
    polar_z = zpos - pol_cent_z;

  
    Double_t pi = TMath::Pi();
    
    radius = TMath::Sqrt(TMath::Power(pol_cent_y-ypos,2) +  TMath::Power(pol_cent_z-zpos,2));
    
    if (polar_y >= 0 && polar_z >= 0) {
      phi = TMath::ATan( polar_z/polar_y);
    } else if ( polar_y < 0 && polar_z >= 0) {
      phi = TMath::ATan( (-polar_y)/polar_z) + pi / 2.0;
    } else if ( polar_y < 0 && polar_z < 0) {
      phi = TMath::ATan( polar_z/polar_y) + pi;
    } else if ( polar_y >= 0 && polar_z < 0) {
      phi = TMath::ATan( polar_y/(-polar_z)) + pi * 3.0/2.0;
    }
  
    scan_data_tree->Fill();
  
  
    }
  }
}
  f_out->cd();
  tge_drift_time->Write();
  tge_drift_time_uncert->Write();

  mg->Add(tge_drift_time);
  leg->AddEntry(tge_drift_time,Form("electron %d",my_elno),"l");
// new TCanvas();
// tge_drift_time->Draw("APL");
// new TCanvas();
// tge_drift_time_uncert->Draw("APL");
  
}
f_out->cd();

new TCanvas();
  mg->Draw("APL");
  mg->SetName("multigraph");
  mg->SetTitle("multigraph");
  mg->GetXaxis()->SetTitle("y pos (um)");
  mg->GetYaxis()->SetTitle("drift time (ns)");
  leg->Draw();

// TCanvas* fish_cavas = new TCanvas();
// fish_cavas->Divide(3,2);
// // draw a fish
// for (Int_t i = 0; i < 6; ++i){
//   fish_cavas->cd(i+1);
//   fish_tree->Draw(Form("(t_drift_b-t_drift_a):(t_drift_b+t_drift_a)>>fish%d(500,-200,300,200,-100,100)",i),
//                   Form("t_drift_b <1000 && elno == %d",i),"colz");
//   
// }




// th2_first_e_0->Draw(); // draw the first fit parameter (constant, in this case)
// new TCanvas();
// th2_first_e_1->Draw(); // draw the second fit parameter (mean, in this case)
// new TCanvas();
// th2_first_e_2->Draw(); // draw the third fit parameter (sigma, in this case)

f_out->Write();

  
}
