///// triobject pt plot

void ctpt()
{
//=========Macro generated from canvas: c1_n8/c1_n8
//=========  (Tue Oct  3 15:49:42 2017) by ROOT version6.08/00
   TCanvas *c1_n8 = new TCanvas("c1_n8", "c1_n8",0,45,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n8->Range(0,0,1,1);
   c1_n8->SetFillColor(0);
   c1_n8->SetBorderMode(0);
   c1_n8->SetBorderSize(2);
   c1_n8->SetTickx(1);
   c1_n8->SetTicky(1);
   c1_n8->SetLeftMargin(0.16);
   c1_n8->SetRightMargin(0.04);
   c1_n8->SetTopMargin(0.08);
   c1_n8->SetBottomMargin(0.12);
   c1_n8->SetFrameFillStyle(0);
   c1_n8->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: tpt
   TPad *tpt = new TPad("tpt", "tpt",0,0.42,1,1);
   tpt->Draw();
   tpt->cd();
   tpt->Range(0,0,1,1);
   tpt->SetFillColor(0);
   tpt->SetBorderMode(0);
   tpt->SetBorderSize(2);
   tpt->SetLogy();
   tpt->SetTickx(1);
   tpt->SetTicky(1);
   tpt->SetLeftMargin(0.16);
   tpt->SetRightMargin(0.04);
   tpt->SetTopMargin(0.08);
   tpt->SetBottomMargin(0);
   tpt->SetFrameFillStyle(0);
   tpt->SetFrameBorderMode(0);
   
   TH1D *hTriobjectPt__36 = new TH1D("hTriobjectPt__36","",50,0,200);
   hTriobjectPt__36->SetBinContent(1,8834);
   hTriobjectPt__36->SetBinContent(2,13318);
   hTriobjectPt__36->SetBinContent(3,9407);
   hTriobjectPt__36->SetBinContent(4,6544);
   hTriobjectPt__36->SetBinContent(5,4472);
   hTriobjectPt__36->SetBinContent(6,3315);
   hTriobjectPt__36->SetBinContent(7,2505);
   hTriobjectPt__36->SetBinContent(8,1706);
   hTriobjectPt__36->SetBinContent(9,1338);
   hTriobjectPt__36->SetBinContent(10,1022);
   hTriobjectPt__36->SetBinContent(11,784);
   hTriobjectPt__36->SetBinContent(12,634);
   hTriobjectPt__36->SetBinContent(13,459);
   hTriobjectPt__36->SetBinContent(14,431);
   hTriobjectPt__36->SetBinContent(15,287);
   hTriobjectPt__36->SetBinContent(16,228);
   hTriobjectPt__36->SetBinContent(17,164);
   hTriobjectPt__36->SetBinContent(18,123);
   hTriobjectPt__36->SetBinContent(19,97);
   hTriobjectPt__36->SetBinContent(20,79);
   hTriobjectPt__36->SetBinContent(21,52);
   hTriobjectPt__36->SetBinContent(22,47);
   hTriobjectPt__36->SetBinContent(23,25);
   hTriobjectPt__36->SetBinContent(24,24);
   hTriobjectPt__36->SetBinContent(25,18);
   hTriobjectPt__36->SetBinContent(26,12);
   hTriobjectPt__36->SetBinContent(27,8);
   hTriobjectPt__36->SetBinContent(28,6);
   hTriobjectPt__36->SetBinContent(29,7);
   hTriobjectPt__36->SetBinContent(30,7);
   hTriobjectPt__36->SetBinContent(31,3);
   hTriobjectPt__36->SetBinContent(32,5);
   hTriobjectPt__36->SetBinContent(33,7);
   hTriobjectPt__36->SetBinContent(34,1);
   hTriobjectPt__36->SetBinContent(35,4);
   hTriobjectPt__36->SetBinContent(38,2);
   hTriobjectPt__36->SetBinContent(39,4);
   hTriobjectPt__36->SetBinContent(40,2);
   hTriobjectPt__36->SetBinContent(41,1);
   hTriobjectPt__36->SetBinContent(44,3);
   hTriobjectPt__36->SetBinContent(46,1);
   hTriobjectPt__36->SetBinContent(48,1);
   hTriobjectPt__36->SetBinContent(50,1);
   hTriobjectPt__36->SetBinContent(51,8);
   hTriobjectPt__36->SetEntries(55996);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cc0000");
   hTriobjectPt__36->SetFillColor(ci);
   hTriobjectPt__36->SetLineStyle(0);
   hTriobjectPt__36->SetMarkerStyle(20);
   hTriobjectPt__36->SetMarkerSize(0);
   hTriobjectPt__36->GetXaxis()->SetTitle("p_{T}_{ll} [GeV/c]");
   hTriobjectPt__36->GetYaxis()->SetTitle("Entries");
   hTriobjectPt__36->GetXaxis()->SetNdivisions(505);
   hTriobjectPt__36->GetXaxis()->SetLabelFont(42);
   hTriobjectPt__36->GetXaxis()->SetLabelOffset(0.007);
   hTriobjectPt__36->GetXaxis()->SetLabelSize(0.05);
   hTriobjectPt__36->GetXaxis()->SetTitleSize(0.06);
   hTriobjectPt__36->GetXaxis()->SetTitleOffset(0.9);
   hTriobjectPt__36->GetXaxis()->SetTitleFont(42);
   hTriobjectPt__36->GetYaxis()->SetLabelFont(42);
   hTriobjectPt__36->GetYaxis()->SetLabelOffset(0.007);
   hTriobjectPt__36->GetYaxis()->SetLabelSize(0.05);
   hTriobjectPt__36->GetYaxis()->SetTitleSize(0.085);
   hTriobjectPt__36->GetYaxis()->SetTitleOffset(0.65);
   hTriobjectPt__36->GetYaxis()->SetTitleFont(42);
   hTriobjectPt__36->GetZaxis()->SetNdivisions(505);
   hTriobjectPt__36->GetZaxis()->SetLabelFont(42);
   hTriobjectPt__36->GetZaxis()->SetLabelOffset(0.007);
   hTriobjectPt__36->GetZaxis()->SetLabelSize(0.05);
   hTriobjectPt__36->GetZaxis()->SetTitleSize(0.06);
   hTriobjectPt__36->GetZaxis()->SetTitleFont(42);
   hTriobjectPt__36->Draw("e2");
   
   TH1D *hTriobjectPt__37 = new TH1D("hTriobjectPt__37","",50,0,200);
   hTriobjectPt__37->SetBinContent(1,9146);
   hTriobjectPt__37->SetBinContent(2,13580);
   hTriobjectPt__37->SetBinContent(3,9527);
   hTriobjectPt__37->SetBinContent(4,6474);
   hTriobjectPt__37->SetBinContent(5,4498);
   hTriobjectPt__37->SetBinContent(6,3288);
   hTriobjectPt__37->SetBinContent(7,2481);
   hTriobjectPt__37->SetBinContent(8,1860);
   hTriobjectPt__37->SetBinContent(9,1344);
   hTriobjectPt__37->SetBinContent(10,983);
   hTriobjectPt__37->SetBinContent(11,763);
   hTriobjectPt__37->SetBinContent(12,624);
   hTriobjectPt__37->SetBinContent(13,481);
   hTriobjectPt__37->SetBinContent(14,400);
   hTriobjectPt__37->SetBinContent(15,288);
   hTriobjectPt__37->SetBinContent(16,230);
   hTriobjectPt__37->SetBinContent(17,168);
   hTriobjectPt__37->SetBinContent(18,141);
   hTriobjectPt__37->SetBinContent(19,118);
   hTriobjectPt__37->SetBinContent(20,72);
   hTriobjectPt__37->SetBinContent(21,61);
   hTriobjectPt__37->SetBinContent(22,32);
   hTriobjectPt__37->SetBinContent(23,38);
   hTriobjectPt__37->SetBinContent(24,17);
   hTriobjectPt__37->SetBinContent(25,13);
   hTriobjectPt__37->SetBinContent(26,11);
   hTriobjectPt__37->SetBinContent(27,11);
   hTriobjectPt__37->SetBinContent(28,10);
   hTriobjectPt__37->SetBinContent(29,6);
   hTriobjectPt__37->SetBinContent(30,8);
   hTriobjectPt__37->SetBinContent(31,7);
   hTriobjectPt__37->SetBinContent(32,6);
   hTriobjectPt__37->SetBinContent(33,1);
   hTriobjectPt__37->SetBinContent(36,6);
   hTriobjectPt__37->SetBinContent(37,2);
   hTriobjectPt__37->SetBinContent(38,4);
   hTriobjectPt__37->SetBinContent(39,1);
   hTriobjectPt__37->SetBinContent(40,2);
   hTriobjectPt__37->SetBinContent(41,3);
   hTriobjectPt__37->SetBinContent(42,1);
   hTriobjectPt__37->SetBinContent(43,3);
   hTriobjectPt__37->SetBinContent(44,1);
   hTriobjectPt__37->SetBinContent(48,1);
   hTriobjectPt__37->SetBinContent(50,1);
   hTriobjectPt__37->SetBinContent(51,13);
   hTriobjectPt__37->SetEntries(56725);

   ci = TColor::GetColor("#00cc00");
   hTriobjectPt__37->SetFillColor(ci);
   hTriobjectPt__37->SetLineStyle(0);
   hTriobjectPt__37->SetMarkerStyle(20);
   hTriobjectPt__37->SetMarkerSize(0);
   hTriobjectPt__37->GetXaxis()->SetTitle("p_{T}(ll#gamma) [GeV/c]");
   hTriobjectPt__37->GetXaxis()->SetNdivisions(505);
   hTriobjectPt__37->GetXaxis()->SetLabelFont(42);
   hTriobjectPt__37->GetXaxis()->SetLabelOffset(0.007);
   hTriobjectPt__37->GetXaxis()->SetLabelSize(0.05);
   hTriobjectPt__37->GetXaxis()->SetTitleSize(0.06);
   hTriobjectPt__37->GetXaxis()->SetTitleOffset(0.9);
   hTriobjectPt__37->GetXaxis()->SetTitleFont(42);
   hTriobjectPt__37->GetYaxis()->SetLabelFont(42);
   hTriobjectPt__37->GetYaxis()->SetLabelOffset(0.007);
   hTriobjectPt__37->GetYaxis()->SetLabelSize(0.05);
   hTriobjectPt__37->GetYaxis()->SetTitleSize(0.06);
   hTriobjectPt__37->GetYaxis()->SetTitleOffset(1.25);
   hTriobjectPt__37->GetYaxis()->SetTitleFont(42);
   hTriobjectPt__37->GetZaxis()->SetNdivisions(505);
   hTriobjectPt__37->GetZaxis()->SetLabelFont(42);
   hTriobjectPt__37->GetZaxis()->SetLabelOffset(0.007);
   hTriobjectPt__37->GetZaxis()->SetLabelSize(0.05);
   hTriobjectPt__37->GetZaxis()->SetTitleSize(0.06);
   hTriobjectPt__37->GetZaxis()->SetTitleFont(42);
   hTriobjectPt__37->Draw("e2same");
   
   TH1D *hTriobjectPt__38 = new TH1D("hTriobjectPt__38","",50,0,200);
   hTriobjectPt__38->SetBinContent(1,8532);
   hTriobjectPt__38->SetBinContent(2,12669);
   hTriobjectPt__38->SetBinContent(3,9178);
   hTriobjectPt__38->SetBinContent(4,6427);
   hTriobjectPt__38->SetBinContent(5,4533);
   hTriobjectPt__38->SetBinContent(6,3249);
   hTriobjectPt__38->SetBinContent(7,2543);
   hTriobjectPt__38->SetBinContent(8,1782);
   hTriobjectPt__38->SetBinContent(9,1315);
   hTriobjectPt__38->SetBinContent(10,987);
   hTriobjectPt__38->SetBinContent(11,784);
   hTriobjectPt__38->SetBinContent(12,652);
   hTriobjectPt__38->SetBinContent(13,455);
   hTriobjectPt__38->SetBinContent(14,365);
   hTriobjectPt__38->SetBinContent(15,300);
   hTriobjectPt__38->SetBinContent(16,207);
   hTriobjectPt__38->SetBinContent(17,172);
   hTriobjectPt__38->SetBinContent(18,131);
   hTriobjectPt__38->SetBinContent(19,83);
   hTriobjectPt__38->SetBinContent(20,82);
   hTriobjectPt__38->SetBinContent(21,63);
   hTriobjectPt__38->SetBinContent(22,33);
   hTriobjectPt__38->SetBinContent(23,40);
   hTriobjectPt__38->SetBinContent(24,22);
   hTriobjectPt__38->SetBinContent(25,23);
   hTriobjectPt__38->SetBinContent(26,11);
   hTriobjectPt__38->SetBinContent(27,7);
   hTriobjectPt__38->SetBinContent(28,8);
   hTriobjectPt__38->SetBinContent(29,11);
   hTriobjectPt__38->SetBinContent(30,4);
   hTriobjectPt__38->SetBinContent(31,5);
   hTriobjectPt__38->SetBinContent(32,4);
   hTriobjectPt__38->SetBinContent(33,3);
   hTriobjectPt__38->SetBinContent(34,4);
   hTriobjectPt__38->SetBinContent(35,1);
   hTriobjectPt__38->SetBinContent(36,4);
   hTriobjectPt__38->SetBinContent(38,1);
   hTriobjectPt__38->SetBinContent(39,2);
   hTriobjectPt__38->SetBinContent(43,1);
   hTriobjectPt__38->SetBinContent(44,1);
   hTriobjectPt__38->SetBinContent(45,1);
   hTriobjectPt__38->SetBinContent(46,1);
   hTriobjectPt__38->SetBinContent(47,1);
   hTriobjectPt__38->SetBinContent(48,2);
   hTriobjectPt__38->SetBinContent(50,3);
   hTriobjectPt__38->SetBinContent(51,8);
   hTriobjectPt__38->SetMinimum(10);
   hTriobjectPt__38->SetMaximum(5000);
   hTriobjectPt__38->SetEntries(54710);

   ci = TColor::GetColor("#3399ff");
   hTriobjectPt__38->SetFillColor(ci);
   hTriobjectPt__38->SetLineStyle(0);
   hTriobjectPt__38->SetMarkerStyle(20);
   hTriobjectPt__38->SetMarkerSize(0);
   hTriobjectPt__38->GetXaxis()->SetTitle("p_{T}(ll#gammma) [GeV/c]");
   hTriobjectPt__38->GetXaxis()->SetNdivisions(505);
   hTriobjectPt__38->GetXaxis()->SetLabelFont(42);
   hTriobjectPt__38->GetXaxis()->SetLabelOffset(0.007);
   hTriobjectPt__38->GetXaxis()->SetLabelSize(0.05);
   hTriobjectPt__38->GetXaxis()->SetTitleSize(0.06);
   hTriobjectPt__38->GetXaxis()->SetTitleOffset(0.9);
   hTriobjectPt__38->GetXaxis()->SetTitleFont(42);
   hTriobjectPt__38->GetYaxis()->SetLabelFont(42);
   hTriobjectPt__38->GetYaxis()->SetLabelOffset(0.007);
   hTriobjectPt__38->GetYaxis()->SetLabelSize(0.05);
   hTriobjectPt__38->GetYaxis()->SetTitleSize(0.06);
   hTriobjectPt__38->GetYaxis()->SetTitleOffset(1.25);
   hTriobjectPt__38->GetYaxis()->SetTitleFont(42);
   hTriobjectPt__38->GetZaxis()->SetNdivisions(505);
   hTriobjectPt__38->GetZaxis()->SetLabelFont(42);
   hTriobjectPt__38->GetZaxis()->SetLabelOffset(0.007);
   hTriobjectPt__38->GetZaxis()->SetLabelSize(0.05);
   hTriobjectPt__38->GetZaxis()->SetTitleSize(0.06);
   hTriobjectPt__38->GetZaxis()->SetTitleFont(42);
   hTriobjectPt__38->Draw("e2same");
   
   TLegend *leg = new TLegend(0.5,0.6,0.85,0.85,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.07);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hTriobjectPt__38","pp #rightarrow ee#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hTriobjectPt__37","pp #rightarrow #mu#mu#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hTriobjectPt__36","pp #rightarrow #tau#tau#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   tpt->Modified();
   c1_n8->cd();
  
// ------------>Primitives in pad: tptr
   TPad *tptr = new TPad("tptr", "tptr",0,0.0,1,0.42);
   tptr->Draw();
   tptr->cd();
   tptr->Range(0,0,1,1);
   tptr->SetFillColor(0);
   tptr->SetBorderMode(0);
   tptr->SetBorderSize(2);
   tptr->SetTickx(1);
   tptr->SetTicky(1);
   tptr->SetLeftMargin(0.16);
   tptr->SetRightMargin(0.04);
   tptr->SetTopMargin(0.08);
   tptr->SetBottomMargin(0.2);
   tptr->SetFrameFillStyle(0);
   tptr->SetFrameBorderMode(0);
   
   TH1D *htptr1__39 = new TH1D("htptr1__39","",50,0,200);
   htptr1__39->SetBinContent(1,1.035396);
   htptr1__39->SetBinContent(2,1.051227);
   htptr1__39->SetBinContent(3,1.024951);
   htptr1__39->SetBinContent(4,1.018204);
   htptr1__39->SetBinContent(5,0.9865431);
   htptr1__39->SetBinContent(6,1.020314);
   htptr1__39->SetBinContent(7,0.985057);
   htptr1__39->SetBinContent(8,0.9573513);
   htptr1__39->SetBinContent(9,1.01749);
   htptr1__39->SetBinContent(10,1.035461);
   htptr1__39->SetBinContent(11,1);
   htptr1__39->SetBinContent(12,0.9723926);
   htptr1__39->SetBinContent(13,1.008791);
   htptr1__39->SetBinContent(14,1.180822);
   htptr1__39->SetBinContent(15,0.9566667);
   htptr1__39->SetBinContent(16,1.101449);
   htptr1__39->SetBinContent(17,0.9534884);
   htptr1__39->SetBinContent(18,0.9389313);
   htptr1__39->SetBinContent(19,1.168675);
   htptr1__39->SetBinContent(20,0.9634146);
   htptr1__39->SetBinContent(21,0.8253968);
   htptr1__39->SetBinContent(22,1.424242);
   htptr1__39->SetBinContent(23,0.625);
   htptr1__39->SetBinContent(24,1.090909);
   htptr1__39->SetBinContent(25,0.7826087);
   htptr1__39->SetBinContent(26,1.090909);
   htptr1__39->SetBinContent(27,1.142857);
   htptr1__39->SetBinContent(28,0.75);
   htptr1__39->SetBinContent(29,0.6363636);
   htptr1__39->SetBinContent(30,1.75);
   htptr1__39->SetBinContent(31,0.6);
   htptr1__39->SetBinContent(32,1.25);
   htptr1__39->SetBinContent(33,2.333333);
   htptr1__39->SetBinContent(34,0.25);
   htptr1__39->SetBinContent(35,4);
   htptr1__39->SetBinContent(38,2);
   htptr1__39->SetBinContent(39,2);
   htptr1__39->SetBinContent(44,3);
   htptr1__39->SetBinContent(46,1);
   htptr1__39->SetBinContent(48,0.5);
   htptr1__39->SetBinContent(50,0.3333333);
   htptr1__39->SetBinContent(51,1);
   htptr1__39->SetBinError(1,0.01571637);
   htptr1__39->SetBinError(2,0.0130462);
   htptr1__39->SetBinError(3,0.01503782);
   htptr1__39->SetBinError(4,0.01788117);
   htptr1__39->SetBinError(5,0.02079286);
   htptr1__39->SetBinError(6,0.02518845);
   htptr1__39->SetBinError(7,0.02772963);
   htptr1__39->SetBinError(8,0.0324277);
   htptr1__39->SetBinError(9,0.03951009);
   htptr1__39->SetBinError(10,0.0462104);
   htptr1__39->SetBinError(11,0.05050763);
   htptr1__39->SetBinError(12,0.05423676);
   htptr1__39->SetBinError(13,0.06673633);
   htptr1__39->SetBinError(14,0.08399547);
   htptr1__39->SetBinError(15,0.07899109);
   htptr1__39->SetBinError(16,0.1057442);
   htptr1__39->SetBinError(17,0.1040636);
   htptr1__39->SetBinError(18,0.1178861);
   htptr1__39->SetBinError(19,0.1747451);
   htptr1__39->SetBinError(20,0.1518818);
   htptr1__39->SetBinError(21,0.1546464);
   htptr1__39->SetBinError(22,0.3234618);
   htptr1__39->SetBinError(23,0.1593444);
   htptr1__39->SetBinError(24,0.321996);
   htptr1__39->SetBinError(25,0.2462841);
   htptr1__39->SetBinError(26,0.4553712);
   htptr1__39->SetBinError(27,0.5914848);
   htptr1__39->SetBinError(28,0.4050463);
   htptr1__39->SetBinError(29,0.3076779);
   htptr1__39->SetBinError(30,1.096871);
   htptr1__39->SetBinError(31,0.438178);
   htptr1__39->SetBinError(32,0.8385255);
   htptr1__39->SetBinError(33,1.610153);
   htptr1__39->SetBinError(34,0.2795085);
   htptr1__39->SetBinError(35,4.472136);
   htptr1__39->SetBinError(38,2.44949);
   htptr1__39->SetBinError(39,1.732051);
   htptr1__39->SetBinError(44,3.464102);
   htptr1__39->SetBinError(46,1.414214);
   htptr1__39->SetBinError(48,0.6123724);
   htptr1__39->SetBinError(50,0.3849002);
   htptr1__39->SetBinError(51,0.5);
   htptr1__39->SetMinimum(0);
   htptr1__39->SetMaximum(2);
   htptr1__39->SetEntries(46.04084);

   ci = TColor::GetColor("#cc0000");
   htptr1__39->SetFillColor(ci);
   htptr1__39->SetLineStyle(0);
   htptr1__39->SetMarkerStyle(20);
   htptr1__39->SetMarkerSize(0);
   htptr1__39->GetXaxis()->SetTitle("p_{T}(ll#gamma) [GeV/c]");
   htptr1__39->GetYaxis()->SetTitle("Ratio");
   htptr1__39->GetXaxis()->CenterTitle(kTRUE);
   htptr1__39->GetXaxis()->SetNdivisions(505);
   htptr1__39->GetXaxis()->SetLabelFont(42);
   htptr1__39->GetXaxis()->SetLabelOffset(0.007);
   htptr1__39->GetXaxis()->SetLabelSize(0.05);
   htptr1__39->GetXaxis()->SetTitleSize(0.1);
   htptr1__39->GetXaxis()->SetTitleOffset(0.75);
   htptr1__39->GetXaxis()->SetTitleFont(42);
   htptr1__39->GetYaxis()->SetLabelFont(42);
   htptr1__39->GetYaxis()->SetLabelOffset(0.007);
   htptr1__39->GetYaxis()->SetLabelSize(0.05);
   htptr1__39->GetYaxis()->SetTitleSize(0.09);
   htptr1__39->GetYaxis()->SetTitleOffset(0.6);
   htptr1__39->GetYaxis()->SetTitleFont(42);
   htptr1__39->GetZaxis()->SetNdivisions(505);
   htptr1__39->GetZaxis()->SetLabelFont(42);
   htptr1__39->GetZaxis()->SetLabelOffset(0.007);
   htptr1__39->GetZaxis()->SetLabelSize(0.05);
   htptr1__39->GetZaxis()->SetTitleSize(0.06);
   htptr1__39->GetZaxis()->SetTitleFont(42);
   htptr1__39->Draw("e2");
   
   TH1D *htptr2__40 = new TH1D("htptr2__40","",50,0,200);
   htptr2__40->SetBinContent(1,1.071964);
   htptr2__40->SetBinContent(2,1.071908);
   htptr2__40->SetBinContent(3,1.038026);
   htptr2__40->SetBinContent(4,1.007313);
   htptr2__40->SetBinContent(5,0.9922788);
   htptr2__40->SetBinContent(6,1.012004);
   htptr2__40->SetBinContent(7,0.9756193);
   htptr2__40->SetBinContent(8,1.043771);
   htptr2__40->SetBinContent(9,1.022053);
   htptr2__40->SetBinContent(10,0.9959473);
   htptr2__40->SetBinContent(11,0.9732143);
   htptr2__40->SetBinContent(12,0.9570552);
   htptr2__40->SetBinContent(13,1.057143);
   htptr2__40->SetBinContent(14,1.09589);
   htptr2__40->SetBinContent(15,0.96);
   htptr2__40->SetBinContent(16,1.111111);
   htptr2__40->SetBinContent(17,0.9767442);
   htptr2__40->SetBinContent(18,1.076336);
   htptr2__40->SetBinContent(19,1.421687);
   htptr2__40->SetBinContent(20,0.8780488);
   htptr2__40->SetBinContent(21,0.968254);
   htptr2__40->SetBinContent(22,0.969697);
   htptr2__40->SetBinContent(23,0.95);
   htptr2__40->SetBinContent(24,0.7727273);
   htptr2__40->SetBinContent(25,0.5652174);
   htptr2__40->SetBinContent(26,1);
   htptr2__40->SetBinContent(27,1.571429);
   htptr2__40->SetBinContent(28,1.25);
   htptr2__40->SetBinContent(29,0.5454545);
   htptr2__40->SetBinContent(30,2);
   htptr2__40->SetBinContent(31,1.4);
   htptr2__40->SetBinContent(32,1.5);
   htptr2__40->SetBinContent(33,0.3333333);
   htptr2__40->SetBinContent(36,1.5);
   htptr2__40->SetBinContent(38,4);
   htptr2__40->SetBinContent(39,0.5);
   htptr2__40->SetBinContent(43,3);
   htptr2__40->SetBinContent(44,1);
   htptr2__40->SetBinContent(48,0.5);
   htptr2__40->SetBinContent(50,0.3333333);
   htptr2__40->SetBinContent(51,1.625);
   htptr2__40->SetBinError(1,0.01613451);
   htptr2__40->SetBinError(2,0.01324015);
   htptr2__40->SetBinError(3,0.01518221);
   htptr2__40->SetBinError(4,0.01773722);
   htptr2__40->SetBinError(5,0.0208833);
   htptr2__40->SetBinError(6,0.02503402);
   htptr2__40->SetBinError(7,0.02753079);
   htptr2__40->SetBinError(8,0.03459911);
   htptr2__40->SetBinError(9,0.03964333);
   htptr2__40->SetBinError(10,0.04487807);
   htptr2__40->SetBinError(11,0.04949181);
   htptr2__40->SetBinError(12,0.05359771);
   htptr2__40->SetBinError(13,0.06913427);
   htptr2__40->SetBinError(14,0.07932707);
   htptr2__40->SetBinError(15,0.07919596);
   htptr2__40->SetBinError(16,0.1064508);
   htptr2__40->SetBinError(17,0.1059501);
   htptr2__40->SetBinError(18,0.1306132);
   htptr2__40->SetBinError(19,0.2036676);
   htptr2__40->SetBinError(20,0.1418097);
   htptr2__40->SetBinError(21,0.1739262);
   htptr2__40->SetBinError(22,0.2405807);
   htptr2__40->SetBinError(23,0.2152034);
   htptr2__40->SetBinError(24,0.24953);
   htptr2__40->SetBinError(25,0.1961242);
   htptr2__40->SetBinError(26,0.4264014);
   htptr2__40->SetBinError(27,0.7597759);
   htptr2__40->SetBinError(28,0.5929271);
   htptr2__40->SetBinError(29,0.2768287);
   htptr2__40->SetBinError(30,1.224745);
   htptr2__40->SetBinError(31,0.8197561);
   htptr2__40->SetBinError(32,0.9682458);
   htptr2__40->SetBinError(33,0.3849002);
   htptr2__40->SetBinError(36,0.9682458);
   htptr2__40->SetBinError(38,4.472136);
   htptr2__40->SetBinError(39,0.6123724);
   htptr2__40->SetBinError(43,3.464102);
   htptr2__40->SetBinError(44,1.414214);
   htptr2__40->SetBinError(48,0.6123724);
   htptr2__40->SetBinError(50,0.3849002);
   htptr2__40->SetBinError(51,0.7302076);
   htptr2__40->SetMinimum(0);
   htptr2__40->SetMaximum(2);
   htptr2__40->SetEntries(50.70308);

   ci = TColor::GetColor("#00cc00");
   htptr2__40->SetFillColor(ci);
   htptr2__40->SetFillStyle(3001);
   htptr2__40->SetLineStyle(0);
   htptr2__40->SetMarkerStyle(20);
   htptr2__40->SetMarkerSize(0);
   htptr2__40->GetXaxis()->SetTitle("p_{T}(ll#gamma) [GeV/c]");
   htptr2__40->GetXaxis()->SetNdivisions(505);
   htptr2__40->GetXaxis()->SetLabelFont(42);
   htptr2__40->GetXaxis()->SetLabelOffset(0.007);
   htptr2__40->GetXaxis()->SetLabelSize(0.05);
   htptr2__40->GetXaxis()->SetTitleSize(0.06);
   htptr2__40->GetXaxis()->SetTitleOffset(0.9);
   htptr2__40->GetXaxis()->SetTitleFont(42);
   htptr2__40->GetYaxis()->SetLabelFont(42);
   htptr2__40->GetYaxis()->SetLabelOffset(0.007);
   htptr2__40->GetYaxis()->SetLabelSize(0.05);
   htptr2__40->GetYaxis()->SetTitleSize(0.06);
   htptr2__40->GetYaxis()->SetTitleOffset(1.25);
   htptr2__40->GetYaxis()->SetTitleFont(42);
   htptr2__40->GetZaxis()->SetNdivisions(505);
   htptr2__40->GetZaxis()->SetLabelFont(42);
   htptr2__40->GetZaxis()->SetLabelOffset(0.007);
   htptr2__40->GetZaxis()->SetLabelSize(0.05);
   htptr2__40->GetZaxis()->SetTitleSize(0.06);
   htptr2__40->GetZaxis()->SetTitleFont(42);
   htptr2__40->Draw("e2same");

   TLegend *leg2 = new TLegend(0.2,0.65,0.55,0.85,NULL,"brNDC");
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.07);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);
   TLegendEntry *entry2=leg2->AddEntry("htptr2__40","#mu#mu#gamma / ee#gamma","f");
   entry2->SetLineColor(1);
   entry2->SetLineStyle(1);
   entry2->SetLineWidth(1);
   entry2->SetMarkerColor(1);
   entry2->SetMarkerStyle(21);
   entry2->SetMarkerSize(1);
   entry2=leg2->AddEntry("htptr1__39","#tau#tau#gamma / ee#gamma","f");
   entry2->SetLineColor(1);
   entry2->SetLineStyle(1);
   entry2->SetLineWidth(1);
   entry2->SetMarkerColor(1);
   entry2->SetMarkerStyle(21);
   entry2->SetMarkerSize(1);
   leg2->Draw();
   TLine *line;
   line = new TLine(0,1,200,1);   
   line->Draw();
   
   
   tptr->Modified();
   c1_n8->cd();
   c1_n8->Modified();
   c1_n8->cd();
   c1_n8->SetSelected(c1_n8);
   c1_n8->SaveAs("plots/TriobjectPt_sherpa.png");
   c1_n8->SaveAs("plots/TriobjectPt_sherpa.pdf");
}
