{

  c_alkane_t0 = new TCanvas("alkane_t0","alkane triplet",200,10,700,500);
  c_alkane_t0->SetGrid();
  
  TGraph *zero_ext = new TGraph("alkane_triplet_integral_size.csv", "%lg %lg %*s %*s %*s %*s %*s", ",");
  zero_ext->SetMarkerColor(kBlue);
  zero_ext->SetMarkerStyle(20);
  zero_ext->SetMarkerSize(2.0);
  zero_ext->SetLineColor(kBlue);
  zero_ext->SetLineWidth(2.0); 
  
  TGraph *one_ext = new TGraph("alkane_triplet_integral_size.csv", "%lg %*s %lg %*s %*s %*s %*s", ",");
  one_ext->SetMarkerColor(kGreen);
  one_ext->SetLineColor(kGreen);  
  one_ext->SetMarkerStyle(21);
  one_ext->SetMarkerSize(2.0);  
  one_ext->SetLineWidth(2.0);

  TGraph *two_ext = new TGraph("alkane_triplet_integral_size.csv", "%lg %*s %*s %lg %*s %*s %*s", ",");
  two_ext->SetMarkerColor(kRed);
  two_ext->SetLineColor(kRed);    
  two_ext->SetMarkerStyle(22);
  two_ext->SetMarkerSize(2.0);  
  two_ext->SetLineWidth(2.0);

  TGraph *three_ext = new TGraph("alkane_triplet_integral_size.csv", "%lg %*s %*s %*s %lg %*s %*s", ",");
  three_ext->SetMarkerColor(kRed);
  three_ext->SetLineColor(kRed);    
  three_ext->SetMarkerStyle(23);
  three_ext->SetMarkerSize(2.0);  
  three_ext->SetLineWidth(2.0);

  TGraph *threeF_ext = new TGraph("alkane_triplet_integral_size.csv", "%lg %*s %*s %*s %*s %lg %*s", ",");
  threeF_ext->SetMarkerColor(kRed);
  threeF_ext->SetLineColor(kRed);    
  threeF_ext->SetMarkerStyle(24);
  threeF_ext->SetMarkerSize(2.0);  
  threeF_ext->SetLineWidth(2.0);

  TGraph *four_ext = new TGraph("alkane_triplet_integral_size.csv", "%lg %*s %*s %*s %*s %*s %lg", ",");
  four_ext->SetMarkerColor(kRed);
  four_ext->SetLineColor(kRed);    
  four_ext->SetMarkerStyle(25);
  four_ext->SetMarkerSize(2.0);  
  four_ext->SetLineWidth(2.0);
  
  TGaxis::SetMaxDigits(4);
  
  TMultiGraph *mg = new TMultiGraph("", ";Number of Carbon Atoms; Size of the PNO integral file (in MB)");
  mg->Add(zero_ext);
  mg->Add(two_ext);
  mg->Add(three_ext);
  mg->Add(threeF_ext);
  mg->Add(four_ext);
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(zero_ext  , "0-ext");
  leg->AddEntry(one_ext   , "1-ext");
  leg->AddEntry(two_ext   , "2-ext");
  leg->AddEntry(three_ext , "3-ext");
  leg->AddEntry(threeF_ext, "full 3-ext");
  leg->AddEntry(four_ext  , "4-ext");
  leg->Draw();
}
