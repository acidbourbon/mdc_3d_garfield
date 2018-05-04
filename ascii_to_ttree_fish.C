
#include "TH1.h"
#include "TF1.h"
#include "TArray.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include <boost/math/distributions/skew_normal.hpp>

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


Int_t mcol(Int_t i){
  const Int_t wheel_size = 7;
  Int_t wheel[wheel_size] = {kBlack, kRed, kBlue, kGreen, kMagenta +1, kCyan +1, kYellow +1};
//   i+=wheel_size;
  if(i % wheel_size == 0){ // a blackish color
    if(i/wheel_size == 0)
      return kBlack;
    if(i/wheel_size == 1)
      return kGray+3;
    if(i/wheel_size == 2)
      return kGray+2;
    if(i/wheel_size == 3)
      return kGray+1;
    if(i/wheel_size == 4)
      return kGray;
  }
  return wheel[i%wheel_size] + i/wheel_size; 
}

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


Float_t skew_norm(Float_t x, Float_t location, Float_t scale, Float_t shape){
  
  if(scale <0.001)
    scale=0.001;
  
  boost::math::skew_normal_distribution<Float_t> dist(location,scale,shape);
  return boost::math::pdf(dist,x);
  
}

Float_t my_skew_norm(Float_t x, Float_t location, Float_t scale, Float_t shape){
  return TMath::Gaus( (x-location) - TMath::Exp(-(-1+shape*(x-location))),0,scale);
  
}


TObjArray* my_fit_slices_y(TH2F* hist, TF1* fit){
  
  TObjArray* return_objects = new TObjArray();
  TH1F* means = new TH1F(Form("%s%s",hist->GetName(),"_mean"),Form("%s%s",hist->GetName(),"_mean"),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  TH1F* stddevs = new TH1F(Form("%s%s",hist->GetName(),"_stddev"),Form("%s%s",hist->GetName(),"_stddev"),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  TH1F* fit_means = new TH1F(Form("%s%s",hist->GetName(),"_fit_mean"),Form("%s%s",hist->GetName(),"_fit_mean"),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  TH1F* fit_stddevs = new TH1F(Form("%s%s",hist->GetName(),"_fit_stddev"),Form("%s%s",hist->GetName(),"_fit_stddev"),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  
  
  
  for (Int_t i = 1; i<= hist->GetNbinsX(); ++i){
    TH1D *py = (TH1D*) hist->ProjectionY("py", i,i); // where firstXbin = 0 and lastXbin = 9
  //   fit->SetParameter(0,py->GetBinContent(py->GetMaximumBin()) );
  //   fit->SetParLimits(0,0,py->GetBinContent(py->GetMaximumBin())*2 );
    fit->SetParameter(1,py->GetMean());
    fit->SetParLimits(1,py->GetMean()-20,py->GetMean()+20);
    fit->SetParameter(2,py->GetStdDev());
    fit->SetParLimits(2,py->GetStdDev()/2,py->GetStdDev()*2);
    fit->FixParameter(3,2.2);
    py->Fit("fit","WW q");
    fit->ReleaseParameter(3);
    py->Fit("fit","WW M q");
    
    
    Float_t location = fit->GetParameter(1);
    Float_t scale = fit->GetParameter(2);
    Float_t shape = fit->GetParameter(3);
    
    Float_t delta = shape/TMath::Sqrt(1+shape*shape);
    Float_t fit_mean  = location + scale*delta*TMath::Sqrt(2/TMath::Pi());
    Float_t fit_stddev = scale*TMath::Sqrt(1-2*delta*delta/TMath::Pi());
    
    means->SetBinContent(i,py->GetMean());
//     means->SetBinError(i,fit->GetParError(1));
//     stddev->SetBinContent(i,fit->GetParameter(2));
//     stddev->SetBinError(i,fit->GetParError(2));
    stddevs->SetBinContent(i,py->GetStdDev());
//     stddevs->SetBinError(i,fit->GetParError(2));
    
    fit_means->SetBinContent(i,fit_mean);
    fit_means->SetBinError(i,0);
    fit_stddevs->SetBinContent(i,fit_stddev);
    fit_stddevs->SetBinError(i,0);
  }
  
  return_objects->AddLast(fit_means);
  return_objects->AddLast(fit_stddevs);
  return return_objects;
}


void ascii_to_ttree_fish(TString infile) {
  
  
  // this is the fit function also used in the beamtime analysis
    
    Float_t backgnd_scale_0     = 1.509e-01;
    Float_t backgnd_shift_1 = 1.0354;
    Float_t backgnd_scale_2 = 3.124;
    Float_t x_shift         = 5;
  
  
  if (1) { // is PASTTREC data
    backgnd_scale_0 = 0.2;
    backgnd_shift_1 = 1.05;
    backgnd_scale_2 = 2.6;
    x_shift         = 0;
  }    
  
//  TF1 *fit  = new TF1 ("fit",
//                                Form("[0]*TMath::Gaus(x,[1],[2])+%f*[0]*TMath::Gaus(x,[1]+%f*[2],%f*[2])",backgnd_scale_0, backgnd_shift_1, backgnd_scale_2)
// //                                ,-100,+200);
//         fit->SetParameter(1,20);
//         fit->SetParameter(2,3);
//         fit->SetParLimits(1,0,60);
//         fit->SetParLimits(2,0.8,8);
//         
//   TF1 *fit  = new TF1 ("fit","[0]*TMath::Gaus( (x-[1]) - TMath::Exp(-(-1+[3]*(x-[1]))),0,[2])"
//                                ,-100,+200);
//         fit->SetParLimits(0,0,1000);
//         fit->SetParLimits(1,0,100);
//         fit->SetParLimits(2,0.2,100);
//         fit->SetParLimits(3,0.001,1);
//         fit->SetParameter(0,100);
//         fit->SetParameter(1,10);
//         fit->SetParameter(2,1);
//         fit->SetParameter(3,0.01);

  TF1 *fit  = new TF1 ("fit","[0]*skew_norm(x,[1],[2],[3])"
                               ,-100,+200);
//         fit->SetParLimits(0,0,1000);
        fit->SetParLimits(1,0,100);
        fit->SetParLimits(2,0.2,100);
        fit->SetParLimits(3,0.001,10);
//         fit->SetParameter(0,100);
//         fit->SetParameter(1,10);
//         fit->SetParameter(2,1);
//         fit->SetParameter(3,0.01);
  
  // gaus fit
//   fit = new TF1("fit","[0]*TMath::Gaus(x,[1],[2])");
  
  
  
  
  Bool_t draw_pulses = true;
  
  // for fish
  Int_t nth_electron = 0;
  Int_t elno_max = 20;
  
  gStyle->SetOptFit();
  
  gRandom->SetSeed(0);
  
  Float_t t1_noise=from_env_float("t1_noise","0");
  
  TFile* f_out = new TFile("f_out.root","RECREATE");
  
  TFile* infile_root = new TFile(infile+".root");
  
  TTree* garfield_tree;
  
  if(infile_root->IsZombie()){ // does not exist yet
    
    garfield_tree = new TTree("garfield_tree", "garfield_tree");
    
    cout << "convert ascii data ..." << endl;
    garfield_tree->ReadFile(infile, "n:x:y:z:t_drift:wire");
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
  garfield_tree->Draw("t_drift*1000:y*10 >> tdrift_vs_y(602,-3,3,256,0,102.3)","wire == 1","colz");
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
  Int_t samples = 3200/10/2;
  
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
  
  Int_t y_pos_bins = 30;
  
  TH1D* th_esig = new TH1D("th_esig","th_esig;t(s)",samples,0,sample_width);
  TH1*  th_esig_cum = 0;
  
  Float_t t1_offset = 20;
  TH2F* th2_first_e = new TH2F("th2_first_e","th2_first_e;y pos (mm);tdrift first electron (ns)",y_pos_bins,-3,3,samples,0-t1_offset,sample_width*1e9-t1_offset);
  TH2F* th2_second_e = new TH2F("th2_second_e","th2_second_e;y pos (mm);tdrift second electron (ns)",y_pos_bins,-3,3,samples,0-t1_offset,sample_width*1e9-t1_offset);
  TH2F* th2_third_e = new TH2F("th2_third_e","th2_third_e;y pos (mm);tdrift third electron (ns)",y_pos_bins,-3,3,samples,0-t1_offset,sample_width*1e9-t1_offset);
  TH2F* th2_fourth_e = new TH2F("th2_fourth_e","th2_fourth_e;y pos (mm);tdrift fourth electron (ns)",y_pos_bins,-3,3,samples,0-t1_offset,sample_width*1e9-t1_offset);
  TH2F* th2_fifth_e = new TH2F("th2_fifth_e","th2_fifth_e;y pos (mm);tdrift fifth electron (ns)",y_pos_bins,-3,3,samples,0-t1_offset,sample_width*1e9-t1_offset);
  TH2F* th2_sixth_e = new TH2F("th2_sixth_e","th2_sixth_e;y pos (mm);tdrift sixth electron (ns)",y_pos_bins,-3,3,samples,0-t1_offset,sample_width*1e9-t1_offset);
 
//   TH1* th_conv = 0;
  
//   Float_t weight = 1./((Float_t) samples ); // somehow I need that factor, so the convolution preserves the absolute Y values, because electron signals are no dirac peaks
  
  Float_t t_drift;
  Float_t t_drift_first;
  Float_t t_drift_second;
  Float_t t_drift_third;
  Float_t t_drift_fourth;
  Float_t n,x,y,z,last_y,wire;
  Float_t last_n = 1;
  garfield_tree->SetBranchAddress("t_drift",&t_drift);
  garfield_tree->SetBranchAddress("n",&n);
  garfield_tree->SetBranchAddress("x",&x);
  garfield_tree->SetBranchAddress("y",&y);
  garfield_tree->SetBranchAddress("z",&z);
  garfield_tree->SetBranchAddress("wire",&wire);
 
//   new TCanvas();
//   th_kern->Draw();
  
  
  Float_t t_drift_a = 1000;
  Float_t t_drift_b = 1000;
  Int_t   elno = 0;
  TTree* fish_tree = new TTree("fish_tree","fish_tree");
  fish_tree->Branch("t_drift_a",&t_drift_a);
  fish_tree->Branch("t_drift_b",&t_drift_b);
  fish_tree->Branch("elno",&elno);
  
  std::vector<Float_t> t_drift_a_vec;
  std::vector<Float_t> t_drift_b_vec;
  
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
          t_drift_a = t_drift_a_vec[elno] + gRandom->Gaus(0,t1_noise);
        }
        if(t_drift_b_vec.size() > elno){
          t_drift_b = t_drift_b_vec[elno] + gRandom->Gaus(0,t1_noise);
        }
        fish_tree->Fill();
      }
      
      if(t_drift_a_vec.size() > 0){
        th2_first_e->Fill(last_y*10.,t_drift_a_vec[0] + gRandom->Gaus(0,t1_noise));
      }
      if(t_drift_a_vec.size() > 1){
        th2_second_e->Fill(last_y*10.,t_drift_a_vec[1] + gRandom->Gaus(0,t1_noise));
      }
      if(t_drift_a_vec.size() > 2){
        th2_third_e->Fill(last_y*10.,t_drift_a_vec[2] + gRandom->Gaus(0,t1_noise));
      }
      if(t_drift_a_vec.size() > 3){
        th2_fourth_e->Fill(last_y*10.,t_drift_a_vec[3] + gRandom->Gaus(0,t1_noise));
      }
      if(t_drift_a_vec.size() > 4){
        th2_fifth_e->Fill(last_y*10.,t_drift_a_vec[4] + gRandom->Gaus(0,t1_noise));
      }
      if(t_drift_a_vec.size() > 5){
        th2_sixth_e->Fill(last_y*10.,t_drift_a_vec[5] + gRandom->Gaus(0,t1_noise));
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
    last_y = y;
  }
//   pulse_mem->Write();

// new TCanvas();
// th2_first_e->Draw("colz");

f_out->cd();

Float_t max_y_sliced_means = 40.;


/*
TString fit_opt = "q";


th2_first_e->FitSlicesY(fit,1,-1,0,fit_opt);
// new TCanvas();
TH1* sliced_sigmas_first_e = ((TH1*) f_out->Get("th2_first_e_2"));
sliced_sigmas_first_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_first_e->DrawClone();
// new TCanvas();
TH1* sliced_means_first_e = ((TH1*) f_out->Get("th2_first_e_1"));
sliced_means_first_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_first_e->DrawClone();

th2_second_e->FitSlicesY(fit,1,-1,0,fit_opt);
// new TCanvas();
TH1* sliced_sigmas_second_e = ((TH1*) f_out->Get("th2_second_e_2"));
sliced_sigmas_second_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_second_e->DrawClone();
// new TCanvas();
TH1* sliced_means_second_e = ((TH1*) f_out->Get("th2_second_e_1"));
sliced_means_second_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_second_e->DrawClone();

th2_third_e->FitSlicesY(fit,1,-1,0,fit_opt);
// new TCanvas();
TH1* sliced_sigmas_third_e = ((TH1*) f_out->Get("th2_third_e_2"));
sliced_sigmas_third_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_third_e->DrawClone();
// new TCanvas();
TH1* sliced_means_third_e = ((TH1*) f_out->Get("th2_third_e_1"));
sliced_means_third_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_third_e->DrawClone();

th2_fourth_e->FitSlicesY(fit,1,-1,0,fit_opt);
// new TCanvas();
TH1* sliced_sigmas_fourth_e = ((TH1*) f_out->Get("th2_fourth_e_2"));
sliced_sigmas_fourth_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_fourth_e->DrawClone();
// new TCanvas();
TH1* sliced_means_fourth_e = ((TH1*) f_out->Get("th2_fourth_e_1"));
sliced_means_fourth_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_fourth_e->DrawClone();

th2_fifth_e->FitSlicesY(fit,1,-1,0,fit_opt);
// new TCanvas();
TH1* sliced_sigmas_fifth_e = ((TH1*) f_out->Get("th2_fifth_e_2"));
sliced_sigmas_fifth_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_fifth_e->DrawClone();
// new TCanvas();
TH1* sliced_means_fifth_e = ((TH1*) f_out->Get("th2_fifth_e_1"));
sliced_means_fifth_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_fifth_e->DrawClone();

th2_sixth_e->FitSlicesY(fit,1,-1,0,fit_opt);
// new TCanvas();
TH1* sliced_sigmas_sixth_e = ((TH1*) f_out->Get("th2_sixth_e_2"));
sliced_sigmas_sixth_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_sixth_e->DrawClone();
// new TCanvas();
TH1* sliced_means_sixth_e = ((TH1*) f_out->Get("th2_sixth_e_1"));
sliced_means_sixth_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_sixth_e->DrawClone();

*/

TObjArray* first_e_fit_results = my_fit_slices_y(th2_first_e, fit);
TH1F* sliced_means_first_e = (TH1F*) first_e_fit_results->At(0);
TH1F* sliced_sigmas_first_e = (TH1F*) first_e_fit_results->At(1);
sliced_sigmas_first_e->GetYaxis()->SetRangeUser(0,10);
sliced_means_first_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);

TObjArray* second_e_fit_results = my_fit_slices_y(th2_second_e, fit);
TH1F* sliced_means_second_e = (TH1F*) second_e_fit_results->At(0);
TH1F* sliced_sigmas_second_e = (TH1F*) second_e_fit_results->At(1);
sliced_sigmas_second_e->GetYaxis()->SetRangeUser(0,10);
sliced_means_second_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);

TObjArray* third_e_fit_results = my_fit_slices_y(th2_third_e, fit);
TH1F* sliced_means_third_e = (TH1F*) third_e_fit_results->At(0);
TH1F* sliced_sigmas_third_e = (TH1F*) third_e_fit_results->At(1);
sliced_sigmas_third_e->GetYaxis()->SetRangeUser(0,10);
sliced_means_third_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);

TObjArray* fourth_e_fit_results = my_fit_slices_y(th2_fourth_e, fit);
TH1F* sliced_means_fourth_e = (TH1F*) fourth_e_fit_results->At(0);
TH1F* sliced_sigmas_fourth_e = (TH1F*) fourth_e_fit_results->At(1);
sliced_sigmas_fourth_e->GetYaxis()->SetRangeUser(0,10);
sliced_means_fourth_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);

TObjArray* fifth_e_fit_results = my_fit_slices_y(th2_fifth_e, fit);
TH1F* sliced_means_fifth_e = (TH1F*) fifth_e_fit_results->At(0);
TH1F* sliced_sigmas_fifth_e = (TH1F*) fifth_e_fit_results->At(1);
sliced_sigmas_fifth_e->GetYaxis()->SetRangeUser(0,10);
sliced_means_fifth_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);

TObjArray* sixth_e_fit_results = my_fit_slices_y(th2_sixth_e, fit);
TH1F* sliced_means_sixth_e = (TH1F*) sixth_e_fit_results->At(0);
TH1F* sliced_sigmas_sixth_e = (TH1F*) sixth_e_fit_results->At(1);
sliced_sigmas_sixth_e->GetYaxis()->SetRangeUser(0,10);
sliced_means_sixth_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);

sliced_means_first_e->SetLineColor(mcol(1-1));
sliced_sigmas_first_e->SetLineColor(mcol(1-1));
sliced_means_second_e->SetLineColor(mcol(2-1));
sliced_sigmas_second_e->SetLineColor(mcol(2-1));
sliced_means_third_e->SetLineColor(mcol(3-1));
sliced_sigmas_third_e->SetLineColor(mcol(3-1));
sliced_means_fourth_e->SetLineColor(mcol(4-1));
sliced_sigmas_fourth_e->SetLineColor(mcol(4-1));
sliced_means_fifth_e->SetLineColor(mcol(5-1));
sliced_sigmas_fifth_e->SetLineColor(mcol(5-1));
sliced_means_sixth_e->SetLineColor(mcol(6-1));
sliced_sigmas_sixth_e->SetLineColor(mcol(6-1));

sliced_means_first_e->SetMarkerColor(mcol(1-1));
sliced_sigmas_first_e->SetMarkerColor(mcol(1-1));
sliced_means_second_e->SetMarkerColor(mcol(2-1));
sliced_sigmas_second_e->SetMarkerColor(mcol(2-1));
sliced_means_third_e->SetMarkerColor(mcol(3-1));
sliced_sigmas_third_e->SetMarkerColor(mcol(3-1));
sliced_means_fourth_e->SetMarkerColor(mcol(4-1));
sliced_sigmas_fourth_e->SetMarkerColor(mcol(4-1));
sliced_means_fifth_e->SetMarkerColor(mcol(5-1));
sliced_sigmas_fifth_e->SetMarkerColor(mcol(5-1));
sliced_means_sixth_e->SetMarkerColor(mcol(6-1));
sliced_sigmas_sixth_e->SetMarkerColor(mcol(6-1));

sliced_means_first_e->SetMarkerStyle(20+1);
sliced_sigmas_first_e->SetMarkerStyle(20+1);
sliced_means_second_e->SetMarkerStyle(20+2);
sliced_sigmas_second_e->SetMarkerStyle(20+2);
sliced_means_third_e->SetMarkerStyle(20+3);
sliced_sigmas_third_e->SetMarkerStyle(20+3);
sliced_means_fourth_e->SetMarkerStyle(20+1);
sliced_sigmas_fourth_e->SetMarkerStyle(20+1);
sliced_means_fifth_e->SetMarkerStyle(20+2);
sliced_sigmas_fifth_e->SetMarkerStyle(20+2);
sliced_means_sixth_e->SetMarkerStyle(20+3);
sliced_sigmas_sixth_e->SetMarkerStyle(20+3);

sliced_means_first_e->SetTitle("first electron");
sliced_sigmas_first_e->SetTitle("first electron");
sliced_means_second_e->SetTitle("second electron");
sliced_sigmas_second_e->SetTitle("second electron");
sliced_means_third_e->SetTitle("third electron");
sliced_sigmas_third_e->SetTitle("third electron");
sliced_means_fourth_e->SetTitle("fourth electron");
sliced_sigmas_fourth_e->SetTitle("fourth electron");
sliced_means_fifth_e->SetTitle("fifth electron");
sliced_sigmas_fifth_e->SetTitle("fifth electron");
sliced_means_sixth_e->SetTitle("sixth electron");
sliced_sigmas_sixth_e->SetTitle("sixth electron");

// sliced_means_first_e->SetLineColor(1);
// sliced_sigmas_first_e->SetLineColor(1);
// sliced_means_second_e->SetLineColor(2);
// sliced_sigmas_second_e->SetLineColor(2);
// sliced_means_third_e->SetLineColor(3);
// sliced_sigmas_third_e->SetLineColor(3);
// sliced_means_fourth_e->SetLineColor(4);
// sliced_sigmas_fourth_e->SetLineColor(4);

TFile* tf_t1_comp = new TFile("Link_to_Drift_Time.root");
TMultiGraph* mg_t1_comp = (TMultiGraph*) tf_t1_comp->Get("multigraph");
// new TCanvas();
// mg_t1_comp->DrawClone("APL");
f_out->cd();


TCanvas* c_means_first_to_fourth = new TCanvas();
mg_t1_comp->DrawClone("APL");
sliced_means_first_e->DrawClone("same P");
sliced_means_second_e->DrawClone("same P");
sliced_means_third_e->DrawClone("same P");
sliced_means_fourth_e->DrawClone("same P");
sliced_means_fifth_e->DrawClone("same P");
sliced_means_sixth_e->DrawClone("same P");

c_means_first_to_fourth->BuildLegend();

TFile* tf_sigma_comp = new TFile("Link_to_drift_time_uncertainty.root");
TMultiGraph* mg_sigma_comp = (TMultiGraph*) tf_sigma_comp->Get("multigraph");
// new TCanvas();
// mg_sigma_comp->DrawClone("APL");
f_out->cd();

TCanvas* c_sigmas_first_to_fourth = new TCanvas();
mg_sigma_comp->DrawClone("APL");
sliced_sigmas_first_e->DrawClone("same P");
sliced_sigmas_second_e->DrawClone("same P");
sliced_sigmas_third_e->DrawClone("same P");
sliced_sigmas_fourth_e->DrawClone("same P");
sliced_sigmas_fifth_e->DrawClone("same P");
sliced_sigmas_sixth_e->DrawClone("same P");

c_sigmas_first_to_fourth->BuildLegend();



TCanvas* c_first_to_fourth = new TCanvas();
c_first_to_fourth->Divide(3,2);
c_first_to_fourth->cd(1);
th2_first_e->Draw("colz");
c_first_to_fourth->cd(2);
th2_second_e->Draw("colz");
c_first_to_fourth->cd(3);
th2_third_e->Draw("colz");
c_first_to_fourth->cd(4);
th2_fourth_e->Draw("colz");
c_first_to_fourth->cd(5);
th2_fifth_e->Draw("colz");
c_first_to_fourth->cd(6);
th2_sixth_e->Draw("colz");

TCanvas* fish_canvas = new TCanvas();
fish_canvas->Divide(3,2);
// draw a fish
for (Int_t i = 0; i < 6; ++i){
  fish_canvas->cd(i+1);
  fish_tree->Draw(Form("(t_drift_b-t_drift_a):(t_drift_b+t_drift_a)>>fish%d(250,-200,300,200,-100,100)",i),
                  Form("t_drift_b <1000 && elno == %d",i),"colz");
  TH2F* this_fish = (TH2F*) f_out->Get(Form("fish%d",i));
  this_fish->SetMinimum(0);
  this_fish->SetMaximum(30);
  this_fish->GetXaxis()->SetRangeUser(-50,200);
  
}

TCanvas* fish_proj_canvas = new TCanvas();
fish_proj_canvas->Divide(3,2);
// draw a fish
for (Int_t i = 0; i < 6; ++i){
  fish_proj_canvas->cd(i+1);
  fish_tree->Draw(Form("(t_drift_b+t_drift_a)>>fish_proj%d(250,-200,300,200,-100,100)",i),
                  Form("abs(t_drift_b - t_drift_a) < 5 && t_drift_b <1000 && elno == %d",i),"colz");
  TH1F* this_fish = (TH1F*) f_out->Get(Form("fish_proj%d",i));
  this_fish->GetXaxis()->SetRangeUser(-50,200);
  this_fish->Fit("fit","WW q M");
  
}



TCanvas* vw_canv = new TCanvas("vw_canv","vw_canv",1600,700); 
vw_canv->Divide(2,1);
vw_canv->cd(1);

mg_t1_comp->DrawClone("APL");
sliced_means_third_e->DrawClone("same P");
sliced_means_fourth_e->DrawClone("same P");
sliced_means_fifth_e->DrawClone("same P");
// TH1F* fourth_scaled = (TH1F*) sliced_means_fourth_e->Clone();
// fourth_scaled->Scale(0.9);
// for (Int_t i = 1; i < fourth_scaled->GetNbinsX(); ++i){
//   fourth_scaled->SetBinError(i,0);
// }
// fourth_scaled->Draw("same");
vw_canv->cd(1)->BuildLegend();

vw_canv->cd(2);
mg_sigma_comp->DrawClone("AP");
sliced_sigmas_third_e->DrawClone("same P");
sliced_sigmas_fourth_e->DrawClone("same P");
sliced_sigmas_fifth_e->DrawClone("same P");
vw_canv->cd(2)->BuildLegend();

// th2_first_e_0->Draw(); // draw the first fit parameter (constant, in this case)
// new TCanvas();
// th2_first_e_1->Draw(); // draw the second fit parameter (mean, in this case)
// new TCanvas();
// th2_first_e_2->Draw(); // draw the third fit parameter (sigma, in this case)



TCanvas* slice_fit_canvas = new TCanvas();
slice_fit_canvas->Divide(4,5);
// draw a fish
Int_t slice_fit_canvas_pad=1;
for (Int_t i = 10; i < 30; ++i){
  slice_fit_canvas->cd(slice_fit_canvas_pad++);
  TH1D *py = (TH1D*) th2_first_e->ProjectionY("py", i,i)->Clone(); // where firstXbin = 0 and lastXbin = 9
  py->Draw();
//   fit->SetParameter(0,py->GetBinContent(py->GetMaximumBin()) );
//   fit->SetParLimits(0,0,py->GetBinContent(py->GetMaximumBin())*2 );
  fit->SetParameter(1,py->GetMean());
  fit->SetParLimits(1,py->GetMean()-20,py->GetMean()+20);
  fit->SetParameter(2,py->GetStdDev());
  fit->SetParLimits(2,py->GetStdDev()/2,py->GetStdDev()*2);
  fit->FixParameter(3,2.2);
  py->Fit("fit","WW q");
  fit->ReleaseParameter(3);
  py->Fit("fit","WW q M");
  
}

f_out->Write();


  
}
