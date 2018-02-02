#include "TH1.h"
#include "TTree.h"
#include "TVirtualFFT.h"


TH1* fftconvolve(TH1D* h1, TH1D* h2){
  
  
  Int_t samples = h1->GetEntries();
  
  TVirtualFFT *fft_h1 = TVirtualFFT::FFT(1, &samples, "R2C M K");
  fft_h1->SetPoints(h1->GetArray());
  fft_h1->Transform();
  
  TVirtualFFT *fft_h2 = TVirtualFFT::FFT(1, &samples, "R2C M K");
  fft_h2->SetPoints(h2->GetArray());
  fft_h2->Transform();
  
  
 /* 
  new TCanvas();
  TH1 *hr = 0;
  hr = TH1::TransformHisto(fft_h1, hr, "RE");
  hr->SetTitle("Real part of the tranfsorm");
  hr->Draw();
  
  new TCanvas();
  TH1 *him = 0;
  him = TH1::TransformHisto(fft_h1, him, "IM");
  him->SetTitle("Im. part of the  transform");
  him->Draw();
  */
 
 
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
   
  
  TH1* th_conv = fftconvolve(th_kern,th_esig);
  
  new TCanvas();
  th_conv->Draw();
  
  
//   TH1 *hr = 0;
//   hr = TH1::TransformHisto(fft_kern, hr, "RE");
   // https://root-forum.cern.ch/t/convolution-of-two-th1f-histograms/21723/3
   
   
   // https://root.cern.ch/doc/master/FFT_8C.html
}  
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
