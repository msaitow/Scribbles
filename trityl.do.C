{

  trityl = new TCanvas("trityl.do","TCutDO.Trityl",200,10,700,500);
  trityl->Divide(2,2);
    
  TGraph *trityl_do = new TGraph("trityl.do.csv", "%lg %lg", ",");
  trityl_do->SetMarkerColor(kBlue);
  trityl_do->SetMarkerStyle(20);
  trityl_do->SetMarkerSize(2.0);
  trityl_do->SetLineColor(kBlue);
  trityl_do->SetLineWidth(2.0);  
  trityl_do->GetXaxis()->SetTitle("TCutDO");
  trityl_do->GetYaxis()->SetTitle("Correlation Energy Recovered (%)");
  
  trityl->cd(1);
  trityl->cd(1)->SetGrid();
  trityl->cd(1)->SetLogx();  
  trityl_do->Draw("APL");
  
  TGraph *trityl_mkn = new TGraph("trityl.mkn.csv", "%lg %lg", ",");
  trityl_mkn->SetMarkerColor(kBlue);
  trityl_mkn->SetMarkerStyle(20);
  trityl_mkn->SetMarkerSize(2.0);
  trityl_mkn->SetLineColor(kBlue);
  trityl_mkn->SetLineWidth(2.0);  
  trityl_mkn->GetXaxis()->SetTitle("TCutMKN");
  trityl_mkn->GetYaxis()->SetTitle("Correlation Energy Recovered (%)");

  trityl->cd(2);
  trityl->cd(2)->SetGrid();
  trityl->cd(2)->SetLogx();    
  trityl_mkn->Draw("APL");

  TGraph *trityl_strong = new TGraph("trityl.pairs.csv", "%lg %lg %*s", ",");
  trityl_strong->SetMarkerColor(kBlue);
  trityl_strong->SetMarkerStyle(20);
  trityl_strong->SetMarkerSize(2.0);
  trityl_strong->SetLineColor(kBlue);
  trityl_strong->SetLineWidth(2.0);  

  TGraph *trityl_total = new TGraph("trityl.pairs.csv", "%lg %*s %lg", ",");
  trityl_total->SetMarkerColor(kRed);
  trityl_total->SetMarkerStyle(20);
  trityl_total->SetMarkerSize(2.0);
  trityl_total->SetLineColor(kRed);
  trityl_total->SetLineWidth(2.0);  

  trityl->cd(3);
  trityl->cd(3)->SetGrid();
  trityl->cd(3)->SetLogx();    
  
  TMultiGraph *mg = new TMultiGraph("", ";TCutPairs; Correlation Energy Recovered (%)");
  mg->Add(trityl_strong);
  mg->Add(trityl_total);  
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(trityl_strong, "strong pairs");
  leg->AddEntry(trityl_total , "total pairs");  
  leg->Draw();
    
}
