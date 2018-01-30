void ascii_to_ttree(TString infile) {

  TTree* MyTree = new TTree("MyTree", "MyTree");
  MyTree->ReadFile(infile, "x:y:z:t_drift");
//   MyTree->Draw("Energy_1:Energy_2”, "”, "");
//   MyTree->Draw("Energy_1”);
  new TCanvas();
  MyTree->Draw("t_drift");
  new TCanvas();
  MyTree->Draw("t_drift:y","","colz");
  new TCanvas();
  MyTree->Draw("t_drift:z","","colz");
  

}
