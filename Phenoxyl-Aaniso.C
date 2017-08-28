{
  
  c1 = new TCanvas("Phenoxyl-Fgrad","Error in EFG",200,10,700,500);
  c1->SetGrid();
  c1->SetLogx();  
  
  TGraph *x = new TGraph("Phenoxyl-Fgrad.csv", "%lg %lg %*s %*s", ",");
  x->SetMarkerColor(kBlue);
  x->SetMarkerStyle(20);
  x->SetMarkerSize(2.0);
  x->SetLineColor(kBlue);
  x->SetLineWidth(2.0);  

  TGraph *y = new TGraph("Phenoxyl-Fgrad.csv", "%lg %*s %lg %*s", ",");
  y->SetMarkerColor(kBlue);
  y->SetMarkerStyle(21);
  y->SetMarkerSize(2.0);
  y->SetLineColor(kBlue);
  y->SetLineWidth(2.0);  

  TGraph *z = new TGraph("Phenoxyl-Fgrad.csv", "%lg %*s %*s %lg", ",");
  z->SetMarkerColor(kBlue);
  z->SetMarkerStyle(21);
  z->SetMarkerSize(2.0);
  z->SetLineColor(kBlue);
  z->SetLineWidth(2.0);  

  TGaxis::SetMaxDigits(-2);
  
  TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Deviations from canonical values in MHz");
  mg->Add(x);
  mg->Add(y);
  mg->Add(z);     
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(x, "x-component");
  leg->AddEntry(y, "y-component");
  leg->AddEntry(z, "z-component");  
  leg->Draw();
}
