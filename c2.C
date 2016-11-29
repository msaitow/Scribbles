void c2()
{
//=========Macro generated from canvas: c2/Spin Chains:CH2
//=========  (Fri Nov 25 15:53:26 2016) by ROOT version6.08/00
   TCanvas *c2 = new TCanvas("c2", "Spin Chains:CH2",2471,189,1632,1144);
   c2->Range(-0.06000008,99.782,10.06,100.024);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetGridx();
   c2->SetGridy();
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle(";Number of BH Monomers; Correlation Energy Recovered (%)");
   
   Double_t Graph_fx1[9] = {
   1,
   2,
   3,
   4,
   5,
   5.999583,
   6.999166,
   8,
   9};
   Double_t Graph_fy1[9] = {
   99.991,
   99.987,
   99.982,
   99.976,
   99.953,
   99.89921,
   99.90116,
   99.855,
   99.815};
   TGraph *graph = new TGraph(9,Graph_fx1,Graph_fy1);
   graph->SetName("Graph");
   graph->SetTitle("BH.latest.csv");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(2.4);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","BH.latest.csv",100,1,10);
   Graph_Graph1->SetMinimum(99.8);
   Graph_Graph1->SetMaximum(100);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   multigraph->Add(graph,"");
   
   Double_t Graph_fx2[9] = {
   0.9954601,
   1.995043,
   2.994626,
   3.994209,
   5,
   5.999583,
   6.999166,
   7.998749,
   8.998331};
   Double_t Graph_fy2[9] = {
   99.99107,
   99.98717,
   99.98219,
   99.97612,
   99.95315,
   99.92412,
   99.90116,
   99.88318,
   99.89509};
   graph = new TGraph(9,Graph_fx2,Graph_fy2);
   graph->SetName("Graph");
   graph->SetTitle("BH.latest.csv");
   graph->SetFillColor(1);
   graph->SetLineColor(6);
   graph->SetLineStyle(2);
   graph->SetLineWidth(2);
   graph->SetMarkerColor(6);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(2);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","BH.latest.csv",100,1,10);
   Graph_Graph2->SetMinimum(99.8);
   Graph_Graph2->SetMaximum(100);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph2);
   
   multigraph->Add(graph,"");
   
   Double_t Graph_fx3[9] = {
   0.9954601,
   1.995043,
   2.994626,
   3.994209,
   5,
   5.999583,
   6.999166,
   7.998749,
   8.998331};
   Double_t Graph_fy3[9] = {
   99.99107,
   99.98717,
   99.98219,
   99.97612,
   99.95315,
   99.92412,
   99.93019,
   99.90917,
   99.89509};
   graph = new TGraph(9,Graph_fx3,Graph_fy3);
   graph->SetName("Graph");
   graph->SetTitle("BH.latest.csv");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);
   graph->SetLineStyle(5);
   graph->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(2.4);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","BH.latest.csv",100,1,10);
   Graph_Graph3->SetMinimum(99.8);
   Graph_Graph3->SetMaximum(100);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3->SetLineColor(ci);
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph3);
   
   multigraph->Add(graph,"");
   
   Double_t Graph_fx4[9] = {
   1,
   2,
   3,
   4,
   5,
   5.999583,
   7,
   8,
   9};
   Double_t Graph_fy4[9] = {
   99.991,
   99.987,
   99.982,
   99.977,
   99.9681,
   99.94406,
   99.934,
   99.903,
   99.888};
   graph = new TGraph(9,Graph_fx4,Graph_fy4);
   graph->SetName("Graph");
   graph->SetTitle("BH.latest.csv");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ffcc00");
   graph->SetLineColor(ci);
   graph->SetLineStyle(6);
   graph->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(2.4);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","BH.latest.csv",100,1,10);
   Graph_Graph4->SetMinimum(99.8);
   Graph_Graph4->SetMaximum(100);
   Graph_Graph4->SetDirectory(0);
   Graph_Graph4->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph4->SetLineColor(ci);
   Graph_Graph4->GetXaxis()->SetLabelFont(42);
   Graph_Graph4->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetXaxis()->SetTitleFont(42);
   Graph_Graph4->GetYaxis()->SetLabelFont(42);
   Graph_Graph4->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetYaxis()->SetTitleFont(42);
   Graph_Graph4->GetZaxis()->SetLabelFont(42);
   Graph_Graph4->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph4->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph4->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph4);
   
   multigraph->Add(graph,"");
   multigraph->Draw("APL");
   multigraph->GetXaxis()->SetTitle("Number of BH Monomers");
   multigraph->GetXaxis()->SetRange(5,96);
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetLabelSize(0.035);
   multigraph->GetXaxis()->SetTitleSize(0.05);
   multigraph->GetXaxis()->SetTickLength(0);
   multigraph->GetXaxis()->SetTitleOffset(0.88);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle(" Correlation Energy Recovered (%)");
   multigraph->GetYaxis()->CenterTitle(true);
   multigraph->GetYaxis()->SetDecimals();
   multigraph->GetYaxis()->SetLabelFont(42);
   multigraph->GetYaxis()->SetLabelSize(0.035);
   multigraph->GetYaxis()->SetTitleSize(0.05);
   multigraph->GetYaxis()->SetTitleOffset(1.15);
   multigraph->GetYaxis()->SetTitleFont(42);
   
   TLegend *leg = new TLegend(0.1484663,0.1960609,0.5282209,0.396598,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","TscaleTCutPairs=1.00","lpf");
   entry->SetFillColor(1);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(2.4);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","TscaleTCutPairs=0.50","lpf");
   entry->SetFillColor(1);
   entry->SetFillStyle(1001);
   entry->SetLineColor(6);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(6);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","TscaleTCutPairs=0.33","lpf");
   entry->SetFillColor(1);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(5);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(2.4);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","TscaleTCutPairs=0.10","lpf");
   entry->SetFillColor(1);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(6);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(2.4);
   entry->SetTextFont(42);
   leg->Draw();
   c2->Modified();
   c2->cd();
   c2->SetSelected(c2);
}
