
#include "TH1.h"
#include "TArray.h"

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






void ascii_to_ttree(TString infile) {
  
  Bool_t draw_pulses = true;
  
  TFile* f_out = new TFile("f_out.root","RECREATE");
  
  TTree* garfield_tree = new TTree("garfield_tree", "garfield_tree");
  garfield_tree->ReadFile(infile, "n:x:y:z:t_drift");
//   garfield_tree->Draw("Energy_1:Energy_2”, "”, "");
//   garfield_tree->Draw("Energy_1”);
  new TCanvas();
  garfield_tree->Draw("t_drift >> tdrift_h()");
  new TCanvas();
//   garfield_tree->Draw("t_drift:x >> tdrift_vs_x()","","colz");
  garfield_tree->Draw("t_drift:x >> tdrift_vs_x(100,-0.3,0.3,100,-0.02,0.1)","","colz");
  new TCanvas();
//   garfield_tree->Draw("t_drift:z >> tdrift_vs_z()","","colz");
  garfield_tree->Draw("t_drift:z >> tdrift_vs_z(100,-0.3,0.3,100,-0.02,0.1)","","colz");
//   new TBrowser();
  new TCanvas();
//   garfield_tree->Draw("t_drift:y >> tdrift_vs_y()","","colz");
  garfield_tree->Draw("t_drift:y >> tdrift_vs_y(100,-0.3,0.3,100,-0.02,0.1)","","colz");
  
  
  
  /// load the kernel/chamber_IR function
  
  TGraph* tg_kern = new TGraph("chamber_IR.csv","%lg, %lg");
  
//   Float_t IR_y_scaler = 1./230./100.;
//   Float_t IR_y_scaler = 1./230; // divide through the number of charges of one Fe55 pulse
  Float_t IR_y_scaler = 1./100; // divide through the number of charges of one Fe55 pulse
//   Float_t IR_y_scaler = 1.;
  
  Float_t sample_width = 1.6e-6;
  Int_t samples = 3200;
  
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
  
  TH1D* th_esig = new TH1D("th_esig","th_esig;t(s)",samples,0,sample_width);
  
  
  TH1* th_conv = 0;
  
  Float_t weight = 1./((Float_t) samples ); // somehow I need that factor, so the convolution preserves the absolute Y values, because electron signals are no dirac peaks
  
  Float_t t_drift;
  Float_t n;
  Float_t last_n = 1;
  garfield_tree->SetBranchAddress("t_drift",&t_drift);
  garfield_tree->SetBranchAddress("n",&n);
 
  new TCanvas();
//   th_kern->Draw();
  
  Int_t primaries = garfield_tree->GetEntries();
//   primaries = 1;
  for (Int_t i = 0 ; i < primaries + 1; i++){
    if(i < primaries){
      garfield_tree->GetEntry(i);
    } else {
      n++; // to trigger last processing
    }
    
//     cout << "n: " << n << endl;
    if (n > last_n){
//       cout << "new N! " << endl;
//       new TCanvas();
      delete th_conv;
      th_conv = fftconvolve(th_kern,th_esig);
      
      th_conv->GetXaxis()->SetRangeUser(-0.1e-6,0.5e-6);
      th_conv->GetYaxis()->SetRangeUser(-2e-3,0.5e-3);
//       th_esig->DrawClone();
      
      if(draw_pulses && n < 100){
        if(last_n == 1){
          th_conv->DrawClone();
        } else {
          th_conv->DrawClone("same");
        }
      }
      th_conv->SetName(Form("%d pulse",i));
//       th_conv->Write();
      hist_to_tarrayf(th_conv,signal_xarr,signal_yarr);
      pulse_mem->Fill();
      th_esig->Reset();
    }
    
    th_esig->Fill(t_drift*1e-6,weight);
//     th_esig->Fill(0.01*1e-6,weight);
    
    last_n = n;
  }
  pulse_mem->Write();
  
}
