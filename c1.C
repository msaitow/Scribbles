void c1()
{
//=========Macro generated from canvas: c1/Spin Chains:CH2
//=========  (Fri Nov 25 15:42:28 2016) by ROOT version6.08/00
   TCanvas *c1 = new TCanvas("c1", "Spin Chains:CH2",2281,138,1596,1191);
   c1->Range(-0.6755811,99.99195,10.7112,100.0524);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetLeftMargin(0.142409);
   c1->SetRightMargin(0.05771644);
   c1->SetTopMargin(0.06935909);
   c1->SetBottomMargin(0.1299385);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineWidth(3);
   c1->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle(";Number of CH2 Monomers; Correlation Energy Recovered (%)");
   
   Double_t Graph_fx1[10] = {
   0.9960048,
   1.996099,
   3.003337,
   3.996287,
   4.996381,
   5.996475,
   6.996569,
   7.996664,
   8.996758,
   9.996852};
   Double_t Graph_fy1[10] = {
   100.02,
   100.042,
   100.046,
   100.016,
   100.013,
   100.0111,
   100.012,
   100.0111,
   100.037,
   100.039};
   TGraph *graph = new TGraph(10,Graph_fx1,Graph_fy1);
   graph->SetName("Graph");
   graph->SetTitle("CH2.latest.csv");
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
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","CH2.latest.csv",100,1,10);
   Graph_Graph1->SetMinimum(100);
   Graph_Graph1->SetMaximum(100.05);
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
   
   Double_t Graph_fx2[10] = {
   1.003148,
   2.003242,
   2.989049,
   4.003431,
   5.003525,
   6.003619,
   7.003713,
   8.003807,
   9.003901,
   10.004};
   Double_t Graph_fy2[10] = {
   100.0201,
   100.0421,
   100.032,
   100.0121,
   100.0071,
   100.0041,
   100.0031,
   100.0021,
   100.0271,
   100.0291};
   graph = new TGraph(10,Graph_fx2,Graph_fy2);
   graph->SetName("Graph");
   graph->SetTitle("CH2.latest.csv");
   graph->SetFillColor(1);
   graph->SetLineColor(6);
   graph->SetLineStyle(2);
   graph->SetLineWidth(2);
   graph->SetMarkerColor(6);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(2.4);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","CH2.latest.csv",100,1,10);
   Graph_Graph2->SetMinimum(100);
   Graph_Graph2->SetMaximum(100.05);
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
   
   Double_t Graph_fx3[10] = {
   1,
   2,
   2.974762,
   4,
   5,
   6,
   7,
   8,
   9,
   10};
   Double_t Graph_fy3[10] = {
   100.02,
   100.042,
   100.0322,
   100.012,
   100.007,
   100.004,
   100.003,
   100.002,
   100.027,
   100.029};
   graph = new TGraph(10,Graph_fx3,Graph_fy3);
   graph->SetName("Graph");
   graph->SetTitle("CH2.latest.csv");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);
   graph->SetLineStyle(5);
   graph->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(2.4);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","CH2.latest.csv",100,1,10);
   Graph_Graph3->SetMinimum(100);
   Graph_Graph3->SetMaximum(100.05);
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
   
   Double_t Graph_fx4[10] = {
   1,
   2,
   2.996193,
   4,
   5,
   6,
   7,
   8,
   9,
   10};
   Double_t Graph_fy4[10] = {
   100.02,
   100.042,
   100.032,
   100.012,
   100.007,
   100.004,
   100.003,
   100.002,
   100.027,
   100.029};
   graph = new TGraph(10,Graph_fx4,Graph_fy4);
   graph->SetName("Graph");
   graph->SetTitle("CH2.latest.csv");
   graph->SetFillColor(1);
   graph->SetLineColor(50);
   graph->SetLineStyle(6);
   graph->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(2.4);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","CH2.latest.csv",100,1,10);
   Graph_Graph4->SetMinimum(100);
   Graph_Graph4->SetMaximum(100.05);
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
   multigraph->GetXaxis()->SetTitle("Number of CH2 Monomers");
   multigraph->GetXaxis()->SetRange(5,96);
   multigraph->GetXaxis()->CenterTitle(true);
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetLabelSize(0.035);
   multigraph->GetXaxis()->SetTitleSize(0.05);
   multigraph->GetXaxis()->SetTickLength(0);
   multigraph->GetXaxis()->SetTitleOffset(0.93);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle(" Correlation Energy Recovered (%)");
   multigraph->GetYaxis()->CenterTitle(true);
   multigraph->GetYaxis()->SetDecimals();
   multigraph->GetYaxis()->SetLabelFont(42);
   multigraph->GetYaxis()->SetLabelSize(0.03);
   multigraph->GetYaxis()->SetTitleSize(0.05);
   multigraph->GetYaxis()->SetTickLength(0.02);
   multigraph->GetYaxis()->SetTitleOffset(1.19);
   multigraph->GetYaxis()->SetTitleFont(42);
   
   TLegend *leg = new TLegend(0.412798,0.5856014,0.7885822,0.8551361,NULL,"brNDC");
   leg->SetBorderSize(1);

   ci = TColor::GetColor("#000000");
   leg->SetTextColor(ci);
   leg->SetLineColor(0);
   leg->SetLineStyle(4);
   leg->SetLineWidth(5);

   ci = TColor::GetColor("#ffffff");
   leg->SetFillColor(ci);
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
   entry->SetMarkerSize(2.4);
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
   entry->SetLineColor(50);
   entry->SetLineStyle(6);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(2.4);
   entry->SetTextFont(42);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->ToggleToolBar();
}
