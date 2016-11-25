//vold reversegraphXaxis ()
//{
//   const Int_t n = 20;
//   Double_t x[n], y[n];
//   for (Int_t i=0;i<n;i++) {
//      x[i] = i*0.1+5;
//      y[i] = 10*sin(x[i]+0.2);
//   }
//   g = new TGraph(n,x,y);
//   g->SetMarkerStyle(21);
//   g->Draw("APL");
//
//   ReverseXAxis(g);
//   ReverseXGraph(g);
//}
//
//vold ReverseXAxis (TGraph *g)
//{
//   // Remove the current axis
//   g->GetXaxis()->SetLabelOffset(999);
//   g->GetXaxis()->SetTickLength(0);
//
//   // Redraw the new axis
//   gPad->Update();
//   TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
//                                gPad->GetUymin(),
//                                gPad->GetUxmin(),      
//                                gPad->GetUymin(),
//                                g->GetXaxis()->GetXmin(),
//                                g->GetXaxis()->GetXmax(),
//                                510,"-SDH");  
//   newaxis->SetLabelOffset(-0.03);    
//   newaxis->Draw();
//}
//
//vold ReverseXGraph (TGraph *g)
//{
//   // Create a new graph
//   Int_t n = g->GetN();
//   Double_t *x = g->GetX();
//   Double_t *y = g->GetY();
//   Double_t xr[100];
//   Double_t dx = g->GetXaxis()->GetXmin()+g->GetXaxis()->GetXmax();
//   for (Int_t i=0; i<n; i++) {
//      xr[i] = -x[i]+dx;
//   }
//
//   gr = new TGraph(n,xr,y);  
//   gr->SetMarkerStyle(20);
//   gr->SetLineColor(kRed);
//   gr->SetMarkerColor(kRed);
//   gr->Draw("PL");
//}

{
  c3 = new TCanvas("c3","semiquinone PNO",200,10,700,500);
  c3->SetGrid();
  c3->SetLogx();  
  
  TGraph *dlpno_strong = new TGraph("semiquinone.csv", "%lg %lg %*s %*s %*s", ",");
  dlpno_strong->SetMarkerColor(kOrange);
  dlpno_strong->SetMarkerStyle(20);
  dlpno_strong->SetMarkerSize(2.0);
  dlpno_strong->SetLineColor(kOrange);
  dlpno_strong->SetLineWidth(2.0);  

  //dlpno_strong->SetMaximum(100.05);
  //dlpno_strong->SetMinimum(100.00);
  //dlpno_strong->GetXaxis()->SetLimits(1,10);

  //ReverseXGraph(dlpno_strong)
    
  TGraph *dlpno_total = new TGraph("semiquinone.csv", "%lg %*s %lg %*s %*s", ",");
  dlpno_total->SetMarkerColor(kRed);
  dlpno_total->SetLineColor(kRed);  
  dlpno_total->SetMarkerStyle(21);
  dlpno_total->SetMarkerSize(2.0);  
  dlpno_total->SetLineWidth(2.0);
  
//  dlpno_total->SetMaximum(100.05);
//  dlpno_total->SetMinimum(100.00);
//  dlpno_total->GetXaxis()->SetLimits(1,10);
  
  TGraph *lpno_strong = new TGraph("semiquinone.csv", "%lg %*s %*s %lg %*s", ",");
  lpno_strong->SetMarkerColor(kGreen);
  lpno_strong->SetLineColor(kGreen);    
  lpno_strong->SetMarkerStyle(22);
  lpno_strong->SetMarkerSize(2.0);  
  lpno_strong->SetLineWidth(2.0);
  
//  lpno_strong->SetMaximum(100.05);
//  lpno_strong->SetMinimum(100.00);
//  lpno_strong->GetXaxis()->SetLimits(1,10);
  
  TGraph *lpno_total = new TGraph("semiquinone.csv", "%lg %*s %*s %*s %lg", ",");
  lpno_total->SetMarkerColor(kBlue);
  lpno_total->SetLineColor(kBlue);    
  lpno_total->SetMarkerStyle(23);
  lpno_total->SetMarkerSize(2.0);  
  lpno_total->SetLineWidth(2.0);
  
//  lpno_total->SetMaximum(100.05);
//  lpno_total->SetMinimum(100.00);
//  lpno_total->GetXaxis()->SetLimits(1,10);
  
  TMultiGraph *mg = new TMultiGraph("", ";TCutPNO; Correlation Energy Recovered (%)");
  mg->Add(dlpno_strong);
  mg->Add(dlpno_total );
  mg->Add( lpno_strong);
  mg->Add( lpno_total ); 
  mg->Draw("APL");

  leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(dlpno_strong, "DLPNO-CCSD (strong pairs)");
  leg->AddEntry(dlpno_total , "DLPNO-CCSD (total pairs)");
  leg->AddEntry(lpno_strong , "LPNO-CCSD (strong pairs)");
  leg->AddEntry(lpno_total  , "LPNO-CCSD (total pairs)");
  leg->Draw();
}
