{

  c2 = new TCanvas("c2","Spin Chains:BH",200,10,700,500);
  c2->SetGrid();
  
  TGraph *scale_100 = new TGraph("BH.latest.csv", "%lg %lg %*s %*s %*s", ",");
  scale_100->SetMarkerColor(kBlue);
  scale_100->SetMarkerStyle(20);
  scale_100->SetMarkerSize(2.0);
  scale_100->SetLineColor(kBlue);
  scale_100->SetLineWidth(2.0);  

//  scale_100->SetMaximum(100.0);
//  scale_100->SetMinimum( 99.8);
  scale_100->GetXaxis()->SetLimits(1,10);
  
  TGraph *scale_050 = new TGraph("BH.latest.csv", "%lg %*s %lg %*s %*s", ",");
  scale_050->SetMarkerColor(kGreen);
  scale_050->SetLineColor(kGreen);  
  scale_050->SetMarkerStyle(21);
  scale_050->SetMarkerSize(2.0);  
  scale_050->SetLineWidth(2.0);
  
//  scale_050->SetMaximum(100.0);
//  scale_050->SetMinimum( 99.8);
  scale_050->GetXaxis()->SetLimits(1,10);
  
  TGraph *scale_033 = new TGraph("BH.latest.csv", "%lg %*s %*s %lg %*s", ",");
  scale_033->SetMarkerColor(kRed);
  scale_033->SetLineColor(kRed);    
  scale_033->SetMarkerStyle(22);
  scale_033->SetMarkerSize(2.0);  
  scale_033->SetLineWidth(2.0);
  
//  scale_033->SetMaximum(100.0);
//  scale_033->SetMinimum( 99.8);
  scale_033->GetXaxis()->SetLimits(1,10);
  
  TGraph *scale_010 = new TGraph("BH.latest.csv", "%lg %*s %*s %*s %lg", ",");
  scale_010->SetMarkerColor(kOrange);
  scale_010->SetLineColor(kOrange);    
  scale_010->SetMarkerStyle(23);
  scale_010->SetMarkerSize(2.0);  
  scale_010->SetLineWidth(2.0);
  
//  scale_010->SetMaximum(100.0);
//  scale_010->SetMinimum( 99.8);
  scale_010->GetXaxis()->SetLimits(1,10);
  
  TMultiGraph *mg = new TMultiGraph("", ";Number of BH Monomers; Correlation Energy Recovered (%)");
  mg->Add(scale_100);
  mg->Add(scale_050);
  mg->Add(scale_033);
  mg->Add(scale_010); 
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(scale_100, "TscaleTCutPairs=1.00");
  leg->AddEntry(scale_050, "TscaleTCutPairs=0.50");
  leg->AddEntry(scale_033, "TscaleTCutPairs=0.33");
  leg->AddEntry(scale_010, "TscaleTCutPairs=0.10");
  leg->Draw();
}
