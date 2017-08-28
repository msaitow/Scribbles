{

  c1 = new TCanvas("Phenoxyl-Aiso","Error in Aiso",200,10,700,500);
  c1->SetGrid();
  c1->SetLogx();  
  
  TGraph *AAE = new TGraph("Phenoxyl-Aiso.csv", "%lg %lg %*s", ",");
  AAE->SetMarkerColor(kBlue);
  AAE->SetMarkerStyle(20);
  AAE->SetMarkerSize(2.0);
  AAE->SetLineColor(kBlue);
  AAE->SetLineWidth(2.0);  
  
  TGraph *Max = new TGraph("Phenoxyl-Aiso.csv", "%lg %*s %lg", ",");
  Max->SetMarkerColor(kBlue);
  Max->SetLineColor(kBlue);  
  Max->SetMarkerStyle(21);
  Max->SetMarkerSize(2.0);  
  Max->SetLineWidth(2.0);
    
  TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Deviations from canonical values in MHz");
  mg->Add(AAE);
  mg->Add(Max);   
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(AAE, "Average absolute deviations");
  leg->AddEntry(Max, "Maximum deviations");
  leg->Draw();
}
