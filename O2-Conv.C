{
  c3 = new TCanvas("c3","O2 PNO",200,10,700,500);
  c3->SetGrid();
  c3->SetLogx();  
      
  TGraph *epr2 = new TGraph("O2-Conv.csv", "%lg %lg %*s", ",");
  epr2->SetMarkerColor(kRed);
  epr2->SetLineColor(kRed);  
  epr2->SetMarkerStyle(21);
  epr2->SetMarkerSize(2.0);  
  epr2->SetLineWidth(2.0);
    
  TGraph *epr3 = new TGraph("O2-Conv.csv", "%lg %*s %lg", ",");
  epr3->SetMarkerColor(kBlue);
  epr3->SetLineColor(kBlue);    
  epr3->SetMarkerStyle(23);
  epr3->SetMarkerSize(2.0);  
  epr3->SetLineWidth(2.0);
    
  TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Errors in Isotropic HFCCs in MHz");
  mg->Add(epr2);
  mg->Add(epr3);
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(epr2, "With EPR2 (Canonical=55.77 MHz)");
  leg->AddEntry(epr3, "With EPR3 (Canonical=52.33 MHz)");
  leg->Draw();
}
