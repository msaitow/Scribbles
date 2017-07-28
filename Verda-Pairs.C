{

  c1 = new TCanvas("Verda-Pairs","Pairs-Pairs",200,10,700,500);
  c1->SetGrid();
  c1->SetLogx();  
  
  TGraph *scale_100 = new TGraph("Pairs-Pairs.csv", "%lg %lg %*s %*s %*s", ",");
  scale_100->SetMarkerColor(kBlue);
  scale_100->SetMarkerStyle(20);
  scale_100->SetMarkerSize(2.0);
  scale_100->SetLineColor(kBlue);
  scale_100->SetLineWidth(2.0);  

///  scale_100->SetMaximum(100.05);
///  scale_100->SetMinimum(100.00);
///  scale_100->GetXaxis()->SetLimits(1,10);
  
  TGraph *scale_050 = new TGraph("Pairs-Pairs.csv", "%lg %*s %lg %*s %*s", ",");
  scale_050->SetMarkerColor(kGreen);
  scale_050->SetLineColor(kGreen);  
  scale_050->SetMarkerStyle(21);
  scale_050->SetMarkerSize(2.0);  
  scale_050->SetLineWidth(2.0);
  
///  scale_050->SetMaximum(100.05);
///  scale_050->SetMinimum(100.00);
///  scale_050->GetXaxis()->SetLimits(1,10);
  
  TGraph *scale_033 = new TGraph("Pairs-Pairs.csv", "%lg %*s %*s %lg %*s", ",");
  scale_033->SetMarkerColor(kRed);
  scale_033->SetLineColor(kRed);    
  scale_033->SetMarkerStyle(22);
  scale_033->SetMarkerSize(2.0);  
  scale_033->SetLineWidth(2.0);
  
///  scale_033->SetMaximum(100.05);
///  scale_033->SetMinimum(100.00);
///  scale_033->GetXaxis()->SetLimits(1,10);
  
  TGraph *scale_010 = new TGraph("Pairs-Pairs.csv", "%lg %*s %*s %*s %lg", ",");
  scale_010->SetMarkerColor(kOrange);
  scale_010->SetLineColor(kOrange);    
  scale_010->SetMarkerStyle(23);
  scale_010->SetMarkerSize(2.0);  
  scale_010->SetLineWidth(2.0);
  
///  scale_010->SetMaximum(100.05);
///  scale_010->SetMinimum(100.00);
///  scale_010->GetXaxis()->SetLimits(1,10);
  
  TMultiGraph *mg = new TMultiGraph("", ";TCutPairs; Isotropic HFCCs in MHz");
  mg->Add(scale_100);
  mg->Add(scale_050);
  mg->Add(scale_033);
  mg->Add(scale_010); 
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(scale_100, "N1");
  leg->AddEntry(scale_050, "N2");
  leg->AddEntry(scale_033, "N3");
  leg->AddEntry(scale_010, "N4");
  leg->Draw();
}
