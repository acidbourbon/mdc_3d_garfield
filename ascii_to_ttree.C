
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
  
  Bool_t draw_pulses = false;
  
  TFile* f_out = new TFile("f_out.root","RECREATE");
  
  TTree* garfield_tree = new TTree("garfield_tree", "garfield_tree");
  garfield_tree->ReadFile(infile, "n:x:y:z:t_drift");
//   garfield_tree->Draw("Energy_1:Energy_2”, "”, "");
//   garfield_tree->Draw("Energy_1”);
  new TCanvas();
  garfield_tree->Draw("t_drift >> tdrift_h()");
  new TCanvas();
  garfield_tree->Draw("t_drift:x >> tdrift_vs_x()","","colz");
  new TCanvas();
  garfield_tree->Draw("t_drift:z >> tdrift_vs_z()","","colz");
//   new TBrowser();
  new TCanvas();
//   garfield_tree->Draw("t_drift:y >> tdrift_vs_y()","","colz");
  garfield_tree->Draw("t_drift:y >> tdrift_vs_y(100,-0.3,0.3,100,-0.02,0.1)","","colz");
  
  
  
  /// load the kernel/chamber_IR function
  
  TGraph* tg_kern = new TGraph("chamber_IR.csv","%lg, %lg");
  
  Float_t IR_y_scaler = 1./230.;
  
  Float_t sample_width = 1.6e-6;
  Int_t samples = 3200;
  
  TH1D* th_kern = new TH1D("th_kern","th_kern;t(ns)",samples,0,sample_width);
  
  for( Int_t i = 1; i <= samples; i++) {
    Float_t t = sample_width * (Float_t) i /(Float_t) samples;
    th_kern->SetBinContent(i,tg_kern->Eval(t-0.1e-6)*IR_y_scaler);
  }
  
  
  // process the garfield tracks
  
  TH1D* th_esig = new TH1D("th_esig","th_esig;t(s)",samples,0,sample_width);
  
  TH1* th_conv = 0;
  
  Float_t t_drift;
  Float_t n;
  Float_t last_n = 1;
  garfield_tree->SetBranchAddress("t_drift",&t_drift);
  garfield_tree->SetBranchAddress("n",&n);
 
  new TCanvas();
  
  Int_t tracks = garfield_tree->GetEntries();
  for (Int_t i = 0 ; i < tracks + 1; i++){
    if(i < tracks){
      garfield_tree->GetEntry(i);
    } else {
      n++; // to trigger last processing
    }
    
//     cout << "n: " << n << endl;
    if (n > last_n){
//       cout << "new N! " << endl;
//       new TCanvas();
      th_conv = fftconvolve(th_kern,th_esig);
      
      th_conv->GetXaxis()->SetRangeUser(-0.1e-6,0.5e-6);
      th_conv->GetYaxis()->SetRangeUser(-3.5,0.5);
//       th_esig->DrawClone();
      
      if(draw_pulses){
        if(last_n == 1){
          th_conv->DrawClone();
        } else {
          th_conv->DrawClone("same");
        }
      }
      th_conv->SetName(Form("%d pulse",i));
      th_conv->Write();
      th_esig->Reset();
    }
    
    th_esig->Fill(t_drift*1e-6);
    
    last_n = n;
  }
  
  
}
