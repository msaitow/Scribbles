{

  c_alkane_s0_ratio = new TCanvas("alkane_s0_ratio","alkane singlet",200,10,700,500);
  c_alkane_s0_ratio->SetGrid();
  
  TGraph *gen_ri = new TGraph("alkane_singlet_uhf_div_rhf.csv", "%lg %lg %*s %*s %*s", ",");
  gen_ri->SetMarkerColor(kBlue);
  gen_ri->SetMarkerStyle(20);
  gen_ri->SetMarkerSize(2.0);
  gen_ri->SetLineColor(kBlue);
  gen_ri->SetLineWidth(2.0);  
  
  TGraph *pno_trafo = new TGraph("alkane_singlet_uhf_div_rhf.csv", "%lg %*s %lg %*s %*s", ",");
  pno_trafo->SetMarkerColor(kGreen);
  pno_trafo->SetLineColor(kGreen);  
  pno_trafo->SetMarkerStyle(21);
  pno_trafo->SetMarkerSize(2.0);  
  pno_trafo->SetLineWidth(2.0);
    
  TGraph *guess = new TGraph("alkane_singlet_uhf_div_rhf.csv", "%lg %*s %*s %lg %*s", ",");
  guess->SetMarkerColor(kRed);
  guess->SetLineColor(kRed);    
  guess->SetMarkerStyle(22);
  guess->SetMarkerSize(2.0);  
  guess->SetLineWidth(2.0);

  TGraph *sigma = new TGraph("alkane_singlet_uhf_div_rhf.csv", "%lg %*s %*s %*s %lg", ",");
  sigma->SetMarkerColor(kMagenta);
  sigma->SetLineColor(kMagenta);    
  sigma->SetMarkerStyle(23);
  sigma->SetMarkerSize(2.0);  
  sigma->SetLineWidth(2.0);
  
  TGaxis::SetMaxDigits(4);
  
  TMultiGraph *mg = new TMultiGraph("", ";Number of Carbon Atoms; Ratio of UHF and RHF Computational Times");
  mg->Add(gen_ri);
  mg->Add(pno_trafo);
  mg->Add(guess);
  mg->Add(sigma);  
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(gen_ri   , "RI 3-Index Integral Generations");
  leg->AddEntry(pno_trafo, "RI-PNO Integral Transformations");
  leg->AddEntry(guess    , "PNO Generation");
  leg->AddEntry(sigma    , "CCSD Iterations");  
  leg->Draw();
}
