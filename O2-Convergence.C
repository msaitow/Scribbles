{
  c = new TCanvas("c3","O2 PNO",200,10,700,500);
  c->SetGrid();
  c->SetLogx();  
      
  TGraph *normal = new TGraph("O2-Convergence.csv", "%lg %lg %*s %*s", ",");  
  normal->SetMarkerColor(kRed);
  normal->SetLineColor(kRed);  
  normal->SetMarkerStyle(21);
  normal->SetMarkerSize(2.0);  
  normal->SetLineWidth(2.0);

  TGraph *core = new TGraph("O2-Convergence.csv", "%lg %*s %lg %*s", ",");  
  core->SetMarkerColor(kBlue);
  core->SetLineColor(kBlue);  
  core->SetMarkerStyle(22);
  core->SetMarkerSize(2.0);  
  core->SetLineWidth(2.0);

  TGraph *coresomo = new TGraph("O2-Convergence.csv", "%lg %*s %*s %lg", ",");  
  coresomo->SetMarkerColor(kOrange);
  coresomo->SetLineColor(kOrange);  
  coresomo->SetMarkerStyle(23);
  coresomo->SetMarkerSize(2.0);  
  coresomo->SetLineWidth(2.0);
  
  //TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Errors in Isotropic HFCCs in MHz");
  TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Accuracy of Isotropic HFCC (%)");
  mg->Add(normal);
  mg->Add(core);
  mg->Add(coresomo);  
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(normal  , "Normal");
  leg->AddEntry(core    , "Tighter Core");
  leg->AddEntry(coresomo, "Tighter Core and SOMO");  
  leg->Draw();
}
