#include "TH1.h"
#include "TTree.h"
#include "TVirtualFFT.h"


void fftconvolve(TH1D* h1, TH1D* h2, TH1D* h_out){
  
}


void convolution(void) {

  TFile* f_out = new TFile("f_out.root","RECREATE");
  
  TTree* MyTree = new TTree("MyTree", "MyTree");
  MyTree->ReadFile("track_drift_line_data.txt", "x:y:z:t_drift");
//   MyTree->Draw("Energy_1:Energy_2”, "”, "");
//   MyTree->Draw("Energy_1”);
  new TCanvas();
  MyTree->Draw("t_drift >> tdrift_h()");
  
//   TTree* KernTree = new TTree("KernTree", "KernTree");
//   KernTree->ReadFile("avg_waveform_meas08.csv", "t:u");
  
  new TCanvas();
  TGraph* tg_kern = new TGraph("avg_waveform_meas08.csv","%lg, %lg");
  tg_kern->Draw(); 
  
  
  Float_t sample_width = 1.6e-6;
  Int_t samples = 3200;
  
  TH1D* th_esig = new TH1D("th_esig","th_esig;t(ns)",samples,0,sample_width);
  
  Float_t t_drift;
  MyTree->SetBranchAddress("t_drift",&t_drift);
  
  for (Int_t i = 0 ; i < MyTree->GetEntries(); i++){
    MyTree->GetEntry(i);
    th_esig->Fill(t_drift*1e-6);
  }
  
  new TCanvas();
  th_esig->Draw();
  
  TH1D* th_kern = new TH1D("th_kern","th_kern;t(ns)",samples,0,sample_width);
  
  for( Int_t i = 1; i <= samples; i++) {
    Float_t t = sample_width * (Float_t) i /(Float_t) samples;
    th_kern->SetBinContent(i,tg_kern->Eval(t-0.1e-6));
  }
  
  new TCanvas();
  th_kern->Draw();
   
  TVirtualFFT *fft_kern = TVirtualFFT::FFT(1, &samples, "R2C M K");
  fft_kern->SetPoints(th_kern->GetArray());
  fft_kern->Transform();
  
  TVirtualFFT *fft_esig = TVirtualFFT::FFT(1, &samples, "R2C M K");
  fft_esig->SetPoints(th_esig->GetArray());
  fft_esig->Transform();

  
  new TCanvas();
  TH1 *hr = 0;
  hr = TH1::TransformHisto(fft_kern, hr, "RE");
  hr->SetTitle("Real part of the tranfsorm");
  hr->Draw();
  
  new TCanvas();
  TH1 *him = 0;
  him = TH1::TransformHisto(fft_kern, him, "IM");
  him->SetTitle("Im. part of the  transform");
  him->Draw();
  
  Double_t *re_kern = new Double_t[samples];
  Double_t *im_kern = new Double_t[samples];
  fft_kern->GetPointsComplex(re_kern,im_kern);
  
  Double_t *re_esig = new Double_t[samples];
  Double_t *im_esig = new Double_t[samples];
  fft_esig->GetPointsComplex(re_esig,im_esig);
 
  
  Double_t *re_conv = new Double_t[samples];
  Double_t *im_conv = new Double_t[samples];
 
  for (Int_t i = 0; i<samples; i++){
    // complex multiplication ... element wise
    re_conv[i] = re_kern[i]*re_esig[i] - im_kern[i]*im_esig[i];
    im_conv[i] = im_kern[i]*re_esig[i] + re_kern[i]*im_esig[i];
  }
  
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &samples, "C2R M K");
  fft_back->SetPointsComplex(re_conv,im_conv);
  fft_back->Transform();
  
  new TCanvas();
  TH1 *hback = 0;
  hback = TH1::TransformHisto(fft_back,hback, "RE");
  hback->GetXaxis()->SetLimits(0,sample_width);
  hback->Draw();
  
//   TH1 *hr = 0;
//   hr = TH1::TransformHisto(fft_kern, hr, "RE");
   // https://root-forum.cern.ch/t/convolution-of-two-th1f-histograms/21723/3
   
   
   // https://root.cern.ch/doc/master/FFT_8C.html
}  
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
