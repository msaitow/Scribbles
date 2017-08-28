{
  c3 = new TCanvas("c3","O2 PNO",200,10,700,500);
  c3->SetGrid();
  c3->SetLogx();  
      
  //TGraph *epr3 = new TGraph("O2-Tight-ScaleCore-1em3-ScaleSOMO-1em1.csv", "%lg %lg", ",");
  TGraph *epr3 = new TGraph("O2-Tight-ScaleCore-1em3-ScaleSOMO-1em1-Percentage.csv", "%lg %lg", ",");  
  epr3->SetMarkerColor(kRed);
  epr3->SetLineColor(kRed);  
  epr3->SetMarkerStyle(21);
  epr3->SetMarkerSize(2.0);  
  epr3->SetLineWidth(2.0);
        
  //TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Errors in Isotropic HFCCs in MHz");
  TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Accuracy of Isotropic HFCC (%)");
  mg->Add(epr3);
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(epr3, "With EPR3 (Canonical=52.3 MHz)");
  leg->Draw();
}
