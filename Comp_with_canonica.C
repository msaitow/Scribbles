{
   int i;
   
//   const Int_t nx = 8;
//   string os_X[nx]   = {"8","32","128","512","2048","8192","32768","131072"};
//   float d_35_0[nx] = {0.75, -3.30, -0.92, 0.10, 0.08, -1.69, -1.29, -2.37};
//   float d_35_1[nx] = {1.01, -3.02, -0.65, 0.37, 0.34, -1.42, -1.02, -2.10};


   const Int_t nx = 7;
   string func[nx]    = {"Default1","Default2","B3LYP","PBE","M06-2X","CAM-B3LYP","B2PLYP"};
   float Aiso[nx]     = {2.1634,1.0304,12.4139,24.5983,37.4135,12.7360,11.7639};
   float Aaniso_1[nx] = {0.2696,0.2287,3.6649 ,4.9592 ,6.2357 ,4.9461 ,3.0394 };
   float Aaniso_2[nx] = {0.2412,0.2289,3.3933 ,4.8864 ,5.8047 ,5.1971 ,2.6797 };
   float Aaniso_3[nx] = {0.3824,0.3179,6.2785 ,8.7344 ,7.9097 ,5.7116 ,5.1611 };
   
   TCanvas *cb = new TCanvas("cb","cb",600,400);
   cb->SetGrid();
   gStyle->SetHistMinimumZero();
   TH1F *h1b = new TH1F("h1b","Average Absolute Errors from canonical CCSD in MHz",nx,0,nx);
   h1b->SetFillColor(4);
   h1b->SetBarWidth(0.2);
   h1b->SetBarOffset(0.1);
   h1b->SetStats(0);
   for (i=1; i<=nx; i++) {
      h1b->SetBinContent(i, Aiso[i-1]);
      h1b->GetXaxis()->SetBinLabel(i,func[i-1].c_str());
   }
   //h1b->GetXaxis()->SetName("AAEs");
   h1b->Draw("b");
   
   TH1F *h2b = new TH1F("h2b","h2b",nx,0,nx);
   h2b->SetFillColor(38);
   h2b->SetBarWidth(0.2);
   h2b->SetBarOffset(0.3);
   h2b->SetStats(0);
   for (i=1;i<=nx;i++) h2b->SetBinContent(i, Aaniso_1[i-1]);
   h2b->Draw("b same");

   TH1F *h3b = new TH1F("h3b","h3b",nx,0,nx);
   h3b->SetFillColor(20);
   h3b->SetBarWidth(0.2);
   h3b->SetBarOffset(0.5);
   h3b->SetStats(0);
   for (i=1;i<=nx;i++) h3b->SetBinContent(i, Aaniso_2[i-1]);
   h3b->Draw("b same");

   TH1F *h4b = new TH1F("h3b","h3b",nx,0,nx);
   h4b->SetFillColor(50);
   h4b->SetBarWidth(0.2);
   h4b->SetBarOffset(0.7);
   h4b->SetStats(0);
   for (i=1;i<=nx;i++) h4b->SetBinContent(i, Aaniso_3[i-1]);
   h4b->Draw("b same");
   
   leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->AddEntry(h1b, "Isotropic HFCC");
   leg->AddEntry(h2b, "x-component of anisotropic HFC tensor");
   leg->AddEntry(h3b, "y-component of anisotropic HFC tensor");
   leg->AddEntry(h4b, "z-component of anisotropic HFC tensor");      
   leg->Draw();
   
   return cb;
}
