
#include "TH1.h"
#include "TF1.h"
#include "TArray.h"
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






void ascii_to_ttree_fish(TString infile) {
  
  Bool_t draw_pulses = true;
  
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
  tdrift_vs_x->GetXaxis()->SetTitle("z pos (mm)");
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
  TH2D* th2_first_e = new TH2D("th2_first_e","th2_first_e;y pos (mm);tdrift first electron (ns)",60,-3,3,samples,0,sample_width*1e9);
  TH2D* th2_second_e = new TH2D("th2_second_e","th2_second_e;y pos (mm);tdrift second electron (ns)",60,-3,3,samples,0,sample_width*1e9);
  TH2D* th2_third_e = new TH2D("th2_third_e","th2_third_e;y pos (mm);tdrift third electron (ns)",60,-3,3,samples,0,sample_width*1e9);
  TH2D* th2_fourth_e = new TH2D("th2_fourth_e","th2_fourth_e;y pos (mm);tdrift fourth electron (ns)",60,-3,3,samples,0,sample_width*1e9);
 
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
  TTree* fish_tree = new TTree("fish_tree","fish_tree");
  fish_tree->Branch("t_drift_a",&t_drift_a);
  fish_tree->Branch("t_drift_b",&t_drift_b);
  
  Int_t primaries = garfield_tree->GetEntries();
//   primaries = 1;
  for (Int_t i = 0 ; i < primaries + 1; i++){
    
    if(i < primaries){
      garfield_tree->GetEntry(i);
    } else {
      n++; // to trigger last processing
    }
    
    if (wire == 1){
      if (t_drift*1000 < t_drift_a){
        t_drift_a = t_drift*1000;
      }
    } else if (wire == 2){
      if (t_drift*1000 < t_drift_b){
        t_drift_b = t_drift*1000;
      }
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

      fish_tree->Fill();
      t_drift_a = 1000;
      t_drift_b = 1000;
      
      if(draw_pulses && n < 100){
        if(last_n == 1){
//           th_conv->DrawClone();
          th_esig->DrawClone();
        } else {
//           th_conv->DrawClone("same");
//           th_esig->DrawClone("same");
        }
      }
      gROOT->cd();
      th_esig_cum = th_esig->GetCumulative();
      t_drift_first  = th_esig->GetXaxis()->GetBinCenter(  th_esig_cum->FindFirstBinAbove(0)  ) ;
      t_drift_second = th_esig->GetXaxis()->GetBinCenter(  th_esig_cum->FindFirstBinAbove(1)  ) ;
      t_drift_third  = th_esig->GetXaxis()->GetBinCenter(  th_esig_cum->FindFirstBinAbove(2)  ) ;
      t_drift_fourth = th_esig->GetXaxis()->GetBinCenter(  th_esig_cum->FindFirstBinAbove(3)  ) ;
      th2_first_e->Fill(last_y*10.,t_drift_first*1e9);
      th2_second_e->Fill(last_y*10.,t_drift_second*1e9);
      th2_third_e->Fill(last_y*10.,t_drift_third*1e9);
      th2_fourth_e->Fill(last_y*10.,t_drift_fourth*1e9);
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


th2_first_e->FitSlicesY(0,1,-1,0,"WW q");
// new TCanvas();
TH1* sliced_sigmas_first_e = ((TH1*) f_out->Get("th2_first_e_2"));
sliced_sigmas_first_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_first_e->DrawClone();
// new TCanvas();
TH1* sliced_means_first_e = ((TH1*) f_out->Get("th2_first_e_1"));
sliced_means_first_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_first_e->DrawClone();

th2_second_e->FitSlicesY(0,1,-1,0,"WW q");
// new TCanvas();
TH1* sliced_sigmas_second_e = ((TH1*) f_out->Get("th2_second_e_2"));
sliced_sigmas_second_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_second_e->DrawClone();
// new TCanvas();
TH1* sliced_means_second_e = ((TH1*) f_out->Get("th2_second_e_1"));
sliced_means_second_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_second_e->DrawClone();

th2_third_e->FitSlicesY(0,1,-1,0,"WW q");
// new TCanvas();
TH1* sliced_sigmas_third_e = ((TH1*) f_out->Get("th2_third_e_2"));
sliced_sigmas_third_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_third_e->DrawClone();
// new TCanvas();
TH1* sliced_means_third_e = ((TH1*) f_out->Get("th2_third_e_1"));
sliced_means_third_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_third_e->DrawClone();

th2_fourth_e->FitSlicesY(0,1,-1,0,"WW q");
// new TCanvas();
TH1* sliced_sigmas_fourth_e = ((TH1*) f_out->Get("th2_fourth_e_2"));
sliced_sigmas_fourth_e->GetYaxis()->SetRangeUser(0,10);
// sliced_sigmas_fourth_e->DrawClone();
// new TCanvas();
TH1* sliced_means_fourth_e = ((TH1*) f_out->Get("th2_fourth_e_1"));
sliced_means_fourth_e->GetYaxis()->SetRangeUser(0,max_y_sliced_means);
// sliced_means_fourth_e->DrawClone();


sliced_means_first_e->SetLineColor(1);
sliced_sigmas_first_e->SetLineColor(1);
sliced_means_second_e->SetLineColor(2);
sliced_sigmas_second_e->SetLineColor(2);
sliced_means_third_e->SetLineColor(3);
sliced_sigmas_third_e->SetLineColor(3);
sliced_means_fourth_e->SetLineColor(4);
sliced_sigmas_fourth_e->SetLineColor(4);

sliced_means_first_e->SetMarkerColor(1);
sliced_sigmas_first_e->SetMarkerColor(1);
sliced_means_second_e->SetMarkerColor(2);
sliced_sigmas_second_e->SetMarkerColor(2);
sliced_means_third_e->SetMarkerColor(3);
sliced_sigmas_third_e->SetMarkerColor(3);
sliced_means_fourth_e->SetMarkerColor(4);
sliced_sigmas_fourth_e->SetMarkerColor(4);

sliced_means_first_e->SetMarkerStyle(20+1);
sliced_sigmas_first_e->SetMarkerStyle(20+1);
sliced_means_second_e->SetMarkerStyle(20+2);
sliced_sigmas_second_e->SetMarkerStyle(20+2);
sliced_means_third_e->SetMarkerStyle(20+3);
sliced_sigmas_third_e->SetMarkerStyle(20+3);
sliced_means_fourth_e->SetMarkerStyle(20+4);
sliced_sigmas_fourth_e->SetMarkerStyle(20+4);

sliced_means_first_e->SetTitle("first electron");
sliced_sigmas_first_e->SetTitle("first electron");
sliced_means_second_e->SetTitle("second electron");
sliced_sigmas_second_e->SetTitle("second electron");
sliced_means_third_e->SetTitle("third electron");
sliced_sigmas_third_e->SetTitle("third electron");
sliced_means_fourth_e->SetTitle("fourth electron");
sliced_sigmas_fourth_e->SetTitle("fourth electron");

// sliced_means_first_e->SetLineColor(1);
// sliced_sigmas_first_e->SetLineColor(1);
// sliced_means_second_e->SetLineColor(2);
// sliced_sigmas_second_e->SetLineColor(2);
// sliced_means_third_e->SetLineColor(3);
// sliced_sigmas_third_e->SetLineColor(3);
// sliced_means_fourth_e->SetLineColor(4);
// sliced_sigmas_fourth_e->SetLineColor(4);


TCanvas* c_means_first_to_fourth = new TCanvas();
sliced_means_first_e->DrawClone();
sliced_means_second_e->DrawClone("same");
sliced_means_third_e->DrawClone("same");
sliced_means_fourth_e->DrawClone("same");

c_means_first_to_fourth->BuildLegend();

TCanvas* c_sigmas_first_to_fourth = new TCanvas();
sliced_sigmas_first_e->DrawClone();
sliced_sigmas_second_e->DrawClone("same");
sliced_sigmas_third_e->DrawClone("same");
sliced_sigmas_fourth_e->DrawClone("same");

c_sigmas_first_to_fourth->BuildLegend();



TCanvas* c_first_to_fourth = new TCanvas();
c_first_to_fourth->Divide(2,2);
c_first_to_fourth->cd(1);
th2_first_e->Draw("colz");
c_first_to_fourth->cd(2);
th2_second_e->Draw("colz");
c_first_to_fourth->cd(3);
th2_third_e->Draw("colz");
c_first_to_fourth->cd(4);
th2_fourth_e->Draw("colz");

new TCanvas();
// draw a fish
fish_tree->Draw("(t_drift_b-t_drift_a):(t_drift_b+t_drift_a)>>fish(100,0,100,100,-100,100)","t_drift_b <1000","colz");




// th2_first_e_0->Draw(); // draw the first fit parameter (constant, in this case)
// new TCanvas();
// th2_first_e_1->Draw(); // draw the second fit parameter (mean, in this case)
// new TCanvas();
// th2_first_e_2->Draw(); // draw the third fit parameter (sigma, in this case)

f_out->Write();

  
}
