{

  c_alkane_t0_steps = new TCanvas("alkane_t0_steps","alkane triplet",200,10,700,500);
  c_alkane_t0_steps->Divide(2);
  
  TGraph *gen_ri = new TGraph("alkane_triplet_steps_per_atom.csv", "%lg %lg %*s %*s %*s %*s %*s %*s", ",");
  gen_ri->SetMarkerColor(kBlue);
  gen_ri->SetMarkerStyle(20);
  gen_ri->SetMarkerSize(2.0);
  gen_ri->SetLineColor(kBlue);
  gen_ri->SetLineWidth(2.0);  
  
  TGraph *pno_trafo = new TGraph("alkane_triplet_steps_per_atom.csv", "%lg %*s %lg %*s %*s %*s %*s %*s", ",");
  pno_trafo->SetMarkerColor(kGreen);
  pno_trafo->SetLineColor(kGreen);  
  pno_trafo->SetMarkerStyle(21);
  pno_trafo->SetMarkerSize(2.0);  
  pno_trafo->SetLineWidth(2.0);
    
  TGraph *guess = new TGraph("alkane_triplet_steps_per_atom.csv", "%lg %*s %*s %lg %*s %*s %*s %*s", ",");
  guess->SetMarkerColor(kRed);
  guess->SetLineColor(kRed);    
  guess->SetMarkerStyle(22);
  guess->SetMarkerSize(2.0);  
  guess->SetLineWidth(2.0);

  TGraph *sigma = new TGraph("alkane_triplet_steps_per_atom.csv", "%lg %*s %*s %*s %lg %*s %*s %*s", ",");
  sigma->SetMarkerColor(kMagenta);
  sigma->SetLineColor(kMagenta);    
  sigma->SetMarkerStyle(23);
  sigma->SetMarkerSize(2.0);  
  sigma->SetLineWidth(2.0);

  TGraph *loc = new TGraph("alkane_triplet_steps_per_atom.csv", "%lg %*s %*s %*s %*s %lg %*s %*s", ",");
  loc->SetMarkerColor(kMagenta);
  loc->SetLineColor(kMagenta);    
  loc->SetMarkerStyle(24);
  loc->SetMarkerSize(2.0);  
  loc->SetLineWidth(2.0);

  TGraph *fock = new TGraph("alkane_singlet_steps_per_atoms.csv", "%lg %*s %*s %*s %*s %*s %lg %*s", ",");
  fock->SetMarkerColor(kMagenta);
  fock->SetLineColor(kMagenta);    
  fock->SetMarkerStyle(25);
  fock->SetMarkerSize(2.0);  
  fock->SetLineWidth(2.0);

  TGraph *doi = new TGraph("alkane_triplet_steps_per_atom.csv", "%lg %*s %*s %*s %*s %*s %*s %lg", ",");
  doi->SetMarkerColor(kMagenta);
  doi->SetLineColor(kMagenta);    
  doi->SetMarkerStyle(26);
  doi->SetMarkerSize(2.0);  
  doi->SetLineWidth(2.0);

  //
  TGraph *crude_iaV = new TGraph("alkane_triplet_ri3index_steps_per_stom.csv", "%lg %lg %*s %*s %*s %*s", ",");
  crude_iaV->SetMarkerColor(kMagenta);
  crude_iaV->SetLineColor(kMagenta);    
  crude_iaV->SetMarkerStyle(27);
  crude_iaV->SetMarkerSize(2.0);  
  crude_iaV->SetLineWidth(2.0);
  //
  TGraph *fine_iaV = new TGraph("alkane_triplet_ri3index_steps_per_stom.csv", "%lg %*s %lg %*s %*s %*s", ",");
  fine_iaV->SetMarkerColor(kMagenta);
  fine_iaV->SetLineColor(kMagenta);    
  fine_iaV->SetMarkerStyle(28);
  fine_iaV->SetMarkerSize(2.0);  
  fine_iaV->SetLineWidth(2.0);
  //
  TGraph *abV = new TGraph("alkane_triplet_ri3index_steps_per_stom.csv", "%lg %*s %*s %lg %*s", ",");
  abV->SetMarkerColor(kMagenta);
  abV->SetLineColor(kMagenta);    
  abV->SetMarkerStyle(29);
  abV->SetMarkerSize(2.0);  
  abV->SetLineWidth(2.0);
  //
  TGraph *ijV = new TGraph("alkane_triplet_ri3index_steps_per_stom.csv", "%lg %*s %*s %*s%lg", ",");
  ijV->SetMarkerColor(kMagenta);
  ijV->SetLineColor(kMagenta);    
  ijV->SetMarkerStyle(30);
  ijV->SetMarkerSize(2.0);  
  ijV->SetLineWidth(2.0);
  
  //TGaxis::SetMaxDigits(4);
  
  TMultiGraph *mg = new TMultiGraph("", ";Number of Carbon Atoms; Computational Times Per Atom (sec.)");
  mg->Add(gen_ri);
  mg->Add(pno_trafo);
  mg->Add(guess);
  mg->Add(sigma);
  mg->Add(loc);
  mg->Add(fock);
  mg->Add(doi);

  c_alkane_t0_steps->cd(1);
  c_alkane_t0_steps->cd(1)->SetGrid();
  
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(gen_ri   , "RI 3-Index Integral Generations");
  leg->AddEntry(pno_trafo, "RI-PNO Integral Transformations");
  leg->AddEntry(guess    , "PNO Generation");
  leg->AddEntry(sigma    , "CCSD Iterations");
  leg->AddEntry(loc      , "Localization of Occupied Orbitals");
  leg->AddEntry(fock     , "Fock Matrix Generations");
  leg->AddEntry(doi      , "Differential Overlap Integral Generations");    
  leg->Draw();

  c_alkane_t0_steps->cd(2);  
  c_alkane_t0_steps->cd(2)->SetGrid();
  
  TMultiGraph *mg2 = new TMultiGraph("", ";Number of Carbon Atoms; Computational Times Per Atom (sec.)");
  mg2->Add(crude_iaV);
  mg2->Add(fine_iaV);
  mg2->Add(abV);  
  mg2->Add(ijV);
  mg2->Draw("APL");
  
  leg2 = new TLegend(0.1,0.7,0.48,0.9);
  leg2->AddEntry(crude_iaV, "1-external for crude guess");
  leg2->AddEntry(fine_iaV , "1-external for fine guess and CCSD");
  leg2->AddEntry(abV      , "2-external");
  leg2->AddEntry(ijV      , "0-external");
  leg2->Draw();
  
}
