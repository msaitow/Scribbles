{

  c_alkane_s0 = new TCanvas("alkane_s0","alkane singlet",200,10,700,500);
  c_alkane_s0->SetGrid();

  //TF1 *f1 = new TF1("f1","[0]*pow(x,[1])");
  TF1 *f1 = new TF1("f1","[0]*pow(3*x+2,[1])");  
  f1->SetParameter(0,5);
  f1->SetParameter(1,1);  
  
  TGraph *rijcosx_uhf = new TGraph("alkane_singlet_scaling.csv", "%lg %lg %*s %*s", ",");
  rijcosx_uhf->SetMarkerColor(kBlue);
  rijcosx_uhf->SetMarkerStyle(20);
  rijcosx_uhf->SetMarkerSize(2.0);
  rijcosx_uhf->SetLineColor(kBlue);
  rijcosx_uhf->SetLineWidth(2.0);
  
  rijcosx_uhf->Fit("f1");
  
  TGraph *uhf_dlpno = new TGraph("alkane_singlet_scaling.csv", "%lg %*s %lg %*s", ",");
  uhf_dlpno->SetMarkerColor(kGreen);
  uhf_dlpno->SetLineColor(kGreen);  
  uhf_dlpno->SetMarkerStyle(21);
  uhf_dlpno->SetMarkerSize(2.0);  
  uhf_dlpno->SetLineWidth(2.0);

  uhf_dlpno->Fit("f1");
  
  TGraph *rhf_dlpno = new TGraph("alkane_singlet_scaling.csv", "%lg %*s %*s %lg", ",");
  rhf_dlpno->SetMarkerColor(kRed);
  rhf_dlpno->SetLineColor(kRed);    
  rhf_dlpno->SetMarkerStyle(22);
  rhf_dlpno->SetMarkerSize(2.0);  
  rhf_dlpno->SetLineWidth(2.0);

  rhf_dlpno->Fit("f1");
  
  TGaxis::SetMaxDigits(4);
  
  TMultiGraph *mg = new TMultiGraph("", ";Number of Carbon Atoms; Computational Time (sec.)");
  mg->Add(rijcosx_uhf);
  mg->Add(uhf_dlpno);
  mg->Add(rhf_dlpno);
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(rijcosx_uhf, "RIJCOSX-UHF");
  leg->AddEntry(uhf_dlpno  , "UHF-DLPNO-CCSD");
  leg->AddEntry(rhf_dlpno  , "RHF-DLPNO-CCSD");
  leg->Draw();
}
