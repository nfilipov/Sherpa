////// dilepton pT plot

void cdpt()
{
//=========Macro generated from canvas: c1_n7/c1_n7
//=========  (Tue Oct  3 15:49:42 2017) by ROOT version6.08/00
   TCanvas *c1_n7 = new TCanvas("c1_n7", "c1_n7",0,45,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n7->Range(0,0,1,1);
   c1_n7->SetFillColor(0);
   c1_n7->SetBorderMode(0);
   c1_n7->SetBorderSize(2);
   c1_n7->SetTickx(1);
   c1_n7->SetTicky(1);
   c1_n7->SetLeftMargin(0.16);
   c1_n7->SetRightMargin(0.04);
   c1_n7->SetTopMargin(0.08);
   c1_n7->SetBottomMargin(0.12);
   c1_n7->SetFrameFillStyle(0);
   c1_n7->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: dpt
   TPad *dpt = new TPad("dpt", "dpt",0,0.42,1,1);
   dpt->Draw();
   dpt->cd();
   dpt->Range(0,0,1,1);
   dpt->SetFillColor(0);
   dpt->SetBorderMode(0);
   dpt->SetBorderSize(2);
   dpt->SetLogy();
   dpt->SetTickx(1);
   dpt->SetTicky(1);
   dpt->SetLeftMargin(0.16);
   dpt->SetRightMargin(0.04);
   dpt->SetTopMargin(0.08);
   dpt->SetBottomMargin(0);
   dpt->SetFrameFillStyle(0);
   dpt->SetFrameBorderMode(0);
   
   TH1D *hDileptonPt__31 = new TH1D("hDileptonPt__31","",50,0,200);
   hDileptonPt__31->SetBinContent(1,1032);
   hDileptonPt__31->SetBinContent(2,3774);
   hDileptonPt__31->SetBinContent(3,7716);
   hDileptonPt__31->SetBinContent(4,10117);
   hDileptonPt__31->SetBinContent(5,9359);
   hDileptonPt__31->SetBinContent(6,7636);
   hDileptonPt__31->SetBinContent(7,5350);
   hDileptonPt__31->SetBinContent(8,3444);
   hDileptonPt__31->SetBinContent(9,2309);
   hDileptonPt__31->SetBinContent(10,1516);
   hDileptonPt__31->SetBinContent(11,995);
   hDileptonPt__31->SetBinContent(12,674);
   hDileptonPt__31->SetBinContent(13,526);
   hDileptonPt__31->SetBinContent(14,375);
   hDileptonPt__31->SetBinContent(15,280);
   hDileptonPt__31->SetBinContent(16,182);
   hDileptonPt__31->SetBinContent(17,141);
   hDileptonPt__31->SetBinContent(18,99);
   hDileptonPt__31->SetBinContent(19,83);
   hDileptonPt__31->SetBinContent(20,61);
   hDileptonPt__31->SetBinContent(21,46);
   hDileptonPt__31->SetBinContent(22,43);
   hDileptonPt__31->SetBinContent(23,35);
   hDileptonPt__31->SetBinContent(24,28);
   hDileptonPt__31->SetBinContent(25,23);
   hDileptonPt__31->SetBinContent(26,13);
   hDileptonPt__31->SetBinContent(27,16);
   hDileptonPt__31->SetBinContent(28,14);
   hDileptonPt__31->SetBinContent(29,14);
   hDileptonPt__31->SetBinContent(30,9);
   hDileptonPt__31->SetBinContent(31,12);
   hDileptonPt__31->SetBinContent(32,9);
   hDileptonPt__31->SetBinContent(33,7);
   hDileptonPt__31->SetBinContent(34,8);
   hDileptonPt__31->SetBinContent(35,5);
   hDileptonPt__31->SetBinContent(36,2);
   hDileptonPt__31->SetBinContent(37,3);
   hDileptonPt__31->SetBinContent(38,5);
   hDileptonPt__31->SetBinContent(39,1);
   hDileptonPt__31->SetBinContent(41,3);
   hDileptonPt__31->SetBinContent(42,1);
   hDileptonPt__31->SetBinContent(43,1);
   hDileptonPt__31->SetBinContent(44,3);
   hDileptonPt__31->SetBinContent(45,2);
   hDileptonPt__31->SetBinContent(46,1);
   hDileptonPt__31->SetBinContent(47,1);
   hDileptonPt__31->SetBinContent(48,1);
   hDileptonPt__31->SetBinContent(49,1);
   hDileptonPt__31->SetBinContent(51,20);
   hDileptonPt__31->SetEntries(55996);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cc0000");
   hDileptonPt__31->SetFillColor(ci);
   hDileptonPt__31->SetLineStyle(0);
   hDileptonPt__31->SetMarkerStyle(20);
   hDileptonPt__31->SetMarkerSize(0);
   hDileptonPt__31->GetXaxis()->SetTitle("p_{T}_{ll} [GeV/c]");
   hDileptonPt__31->GetXaxis()->SetNdivisions(505);
   hDileptonPt__31->GetXaxis()->SetLabelFont(42);
   hDileptonPt__31->GetXaxis()->SetLabelOffset(0.007);
   hDileptonPt__31->GetXaxis()->SetLabelSize(0.05);
   hDileptonPt__31->GetXaxis()->SetTitleSize(0.06);
   hDileptonPt__31->GetXaxis()->SetTitleOffset(0.9);
   hDileptonPt__31->GetXaxis()->SetTitleFont(42);
   hDileptonPt__31->GetYaxis()->SetLabelFont(42);
   hDileptonPt__31->GetYaxis()->SetLabelOffset(0.007);
   hDileptonPt__31->GetYaxis()->SetLabelSize(0.05);
   hDileptonPt__31->GetYaxis()->SetTitleSize(0.085);
   hDileptonPt__31->GetYaxis()->SetTitleOffset(0.65);
   hDileptonPt__31->GetYaxis()->SetTitleFont(42);
   hDileptonPt__31->GetZaxis()->SetNdivisions(505);
   hDileptonPt__31->GetZaxis()->SetLabelFont(42);
   hDileptonPt__31->GetZaxis()->SetLabelOffset(0.007);
   hDileptonPt__31->GetZaxis()->SetLabelSize(0.05);
   hDileptonPt__31->GetZaxis()->SetTitleSize(0.06);
   hDileptonPt__31->GetZaxis()->SetTitleFont(42);
   hDileptonPt__31->GetYaxis()->SetTitle("Entries");
   hDileptonPt__31->Draw("e2");
   
   TH1D *hDileptonPt__32 = new TH1D("hDileptonPt__32","",50,0,200);
   hDileptonPt__32->SetBinContent(1,1047);
   hDileptonPt__32->SetBinContent(2,3618);
   hDileptonPt__32->SetBinContent(3,7814);
   hDileptonPt__32->SetBinContent(4,10274);
   hDileptonPt__32->SetBinContent(5,9577);
   hDileptonPt__32->SetBinContent(6,7553);
   hDileptonPt__32->SetBinContent(7,5565);
   hDileptonPt__32->SetBinContent(8,3647);
   hDileptonPt__32->SetBinContent(9,2281);
   hDileptonPt__32->SetBinContent(10,1472);
   hDileptonPt__32->SetBinContent(11,1005);
   hDileptonPt__32->SetBinContent(12,735);
   hDileptonPt__32->SetBinContent(13,519);
   hDileptonPt__32->SetBinContent(14,363);
   hDileptonPt__32->SetBinContent(15,277);
   hDileptonPt__32->SetBinContent(16,207);
   hDileptonPt__32->SetBinContent(17,147);
   hDileptonPt__32->SetBinContent(18,130);
   hDileptonPt__32->SetBinContent(19,100);
   hDileptonPt__32->SetBinContent(20,62);
   hDileptonPt__32->SetBinContent(21,60);
   hDileptonPt__32->SetBinContent(22,51);
   hDileptonPt__32->SetBinContent(23,23);
   hDileptonPt__32->SetBinContent(24,23);
   hDileptonPt__32->SetBinContent(25,29);
   hDileptonPt__32->SetBinContent(26,17);
   hDileptonPt__32->SetBinContent(27,18);
   hDileptonPt__32->SetBinContent(28,13);
   hDileptonPt__32->SetBinContent(29,9);
   hDileptonPt__32->SetBinContent(30,8);
   hDileptonPt__32->SetBinContent(31,6);
   hDileptonPt__32->SetBinContent(32,2);
   hDileptonPt__32->SetBinContent(33,5);
   hDileptonPt__32->SetBinContent(34,9);
   hDileptonPt__32->SetBinContent(35,1);
   hDileptonPt__32->SetBinContent(36,3);
   hDileptonPt__32->SetBinContent(37,4);
   hDileptonPt__32->SetBinContent(38,3);
   hDileptonPt__32->SetBinContent(39,3);
   hDileptonPt__32->SetBinContent(40,10);
   hDileptonPt__32->SetBinContent(42,1);
   hDileptonPt__32->SetBinContent(43,2);
   hDileptonPt__32->SetBinContent(44,5);
   hDileptonPt__32->SetBinContent(45,2);
   hDileptonPt__32->SetBinContent(46,1);
   hDileptonPt__32->SetBinContent(47,1);
   hDileptonPt__32->SetBinContent(50,5);
   hDileptonPt__32->SetBinContent(51,18);
   hDileptonPt__32->SetEntries(56725);

   ci = TColor::GetColor("#00cc00");
   hDileptonPt__32->SetFillColor(ci);
   hDileptonPt__32->SetLineStyle(0);
   hDileptonPt__32->SetMarkerStyle(20);
   hDileptonPt__32->SetMarkerSize(0);
   hDileptonPt__32->GetXaxis()->SetTitle("p_{T}_{ll} [GeV/c]");
   hDileptonPt__32->GetXaxis()->SetNdivisions(505);
   hDileptonPt__32->GetXaxis()->SetLabelFont(42);
   hDileptonPt__32->GetXaxis()->SetLabelOffset(0.007);
   hDileptonPt__32->GetXaxis()->SetLabelSize(0.05);
   hDileptonPt__32->GetXaxis()->SetTitleSize(0.06);
   hDileptonPt__32->GetXaxis()->SetTitleOffset(0.9);
   hDileptonPt__32->GetXaxis()->SetTitleFont(42);
   hDileptonPt__32->GetYaxis()->SetLabelFont(42);
   hDileptonPt__32->GetYaxis()->SetLabelOffset(0.007);
   hDileptonPt__32->GetYaxis()->SetLabelSize(0.05);
   hDileptonPt__32->GetYaxis()->SetTitleSize(0.06);
   hDileptonPt__32->GetYaxis()->SetTitleOffset(1.25);
   hDileptonPt__32->GetYaxis()->SetTitleFont(42);
   hDileptonPt__32->GetZaxis()->SetNdivisions(505);
   hDileptonPt__32->GetZaxis()->SetLabelFont(42);
   hDileptonPt__32->GetZaxis()->SetLabelOffset(0.007);
   hDileptonPt__32->GetZaxis()->SetLabelSize(0.05);
   hDileptonPt__32->GetZaxis()->SetTitleSize(0.06);
   hDileptonPt__32->GetZaxis()->SetTitleFont(42);
   hDileptonPt__32->Draw("e2same");
   
   TH1D *hDileptonPt__33 = new TH1D("hDileptonPt__33","",50,0,200);
   hDileptonPt__33->SetBinContent(1,1058);
   hDileptonPt__33->SetBinContent(2,3733);
   hDileptonPt__33->SetBinContent(3,7502);
   hDileptonPt__33->SetBinContent(4,9756);
   hDileptonPt__33->SetBinContent(5,9270);
   hDileptonPt__33->SetBinContent(6,7427);
   hDileptonPt__33->SetBinContent(7,5276);
   hDileptonPt__33->SetBinContent(8,3489);
   hDileptonPt__33->SetBinContent(9,2157);
   hDileptonPt__33->SetBinContent(10,1389);
   hDileptonPt__33->SetBinContent(11,915);
   hDileptonPt__33->SetBinContent(12,674);
   hDileptonPt__33->SetBinContent(13,514);
   hDileptonPt__33->SetBinContent(14,369);
   hDileptonPt__33->SetBinContent(15,305);
   hDileptonPt__33->SetBinContent(16,172);
   hDileptonPt__33->SetBinContent(17,134);
   hDileptonPt__33->SetBinContent(18,107);
   hDileptonPt__33->SetBinContent(19,83);
   hDileptonPt__33->SetBinContent(20,65);
   hDileptonPt__33->SetBinContent(21,49);
   hDileptonPt__33->SetBinContent(22,27);
   hDileptonPt__33->SetBinContent(23,31);
   hDileptonPt__33->SetBinContent(24,30);
   hDileptonPt__33->SetBinContent(25,12);
   hDileptonPt__33->SetBinContent(26,19);
   hDileptonPt__33->SetBinContent(27,24);
   hDileptonPt__33->SetBinContent(28,22);
   hDileptonPt__33->SetBinContent(29,11);
   hDileptonPt__33->SetBinContent(30,5);
   hDileptonPt__33->SetBinContent(31,13);
   hDileptonPt__33->SetBinContent(32,5);
   hDileptonPt__33->SetBinContent(33,5);
   hDileptonPt__33->SetBinContent(34,6);
   hDileptonPt__33->SetBinContent(35,5);
   hDileptonPt__33->SetBinContent(36,5);
   hDileptonPt__33->SetBinContent(37,4);
   hDileptonPt__33->SetBinContent(39,2);
   hDileptonPt__33->SetBinContent(40,1);
   hDileptonPt__33->SetBinContent(41,2);
   hDileptonPt__33->SetBinContent(42,1);
   hDileptonPt__33->SetBinContent(43,4);
   hDileptonPt__33->SetBinContent(44,2);
   hDileptonPt__33->SetBinContent(46,5);
   hDileptonPt__33->SetBinContent(47,3);
   hDileptonPt__33->SetBinContent(49,2);
   hDileptonPt__33->SetBinContent(51,20);
   hDileptonPt__33->SetMinimum(10);
   hDileptonPt__33->SetMaximum(5000);
   hDileptonPt__33->SetEntries(54710);

   ci = TColor::GetColor("#3399ff");
   hDileptonPt__33->SetFillColor(ci);
   hDileptonPt__33->SetLineStyle(0);
   hDileptonPt__33->SetMarkerStyle(20);
   hDileptonPt__33->SetMarkerSize(0);
   hDileptonPt__33->GetXaxis()->SetTitle("p_{T}_{ll} [GeV/c]");
   hDileptonPt__33->GetXaxis()->SetNdivisions(505);
   hDileptonPt__33->GetXaxis()->SetLabelFont(42);
   hDileptonPt__33->GetXaxis()->SetLabelOffset(0.007);
   hDileptonPt__33->GetXaxis()->SetLabelSize(0.05);
   hDileptonPt__33->GetXaxis()->SetTitleSize(0.06);
   hDileptonPt__33->GetXaxis()->SetTitleOffset(0.9);
   hDileptonPt__33->GetXaxis()->SetTitleFont(42);
   hDileptonPt__33->GetYaxis()->SetLabelFont(42);
   hDileptonPt__33->GetYaxis()->SetLabelOffset(0.007);
   hDileptonPt__33->GetYaxis()->SetLabelSize(0.05);
   hDileptonPt__33->GetYaxis()->SetTitleSize(0.06);
   hDileptonPt__33->GetYaxis()->SetTitleOffset(1.25);
   hDileptonPt__33->GetYaxis()->SetTitleFont(42);
   hDileptonPt__33->GetZaxis()->SetNdivisions(505);
   hDileptonPt__33->GetZaxis()->SetLabelFont(42);
   hDileptonPt__33->GetZaxis()->SetLabelOffset(0.007);
   hDileptonPt__33->GetZaxis()->SetLabelSize(0.05);
   hDileptonPt__33->GetZaxis()->SetTitleSize(0.06);
   hDileptonPt__33->GetZaxis()->SetTitleFont(42);
   hDileptonPt__33->Draw("e2same");
   
   TLegend *leg = new TLegend(0.5,0.5,0.9,0.8,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.07);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hDileptonPt__33","pp #rightarrow ee#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hDileptonPt__32","pp #rightarrow #mu#mu#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hDileptonPt__31","pp #rightarrow #tau#tau#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   dpt->Modified();
   c1_n7->cd();
  
// ------------>Primitives in pad: dptr
   TPad *dptr = new TPad("dptr", "dptr",0,0.0,1,0.42);
   dptr->Draw();
   dptr->cd();
   dptr->Range(0,0,1,1);
   dptr->SetFillColor(0);
   dptr->SetBorderMode(0);
   dptr->SetBorderSize(2);
   dptr->SetTickx(1);
   dptr->SetTicky(1);
   dptr->SetLeftMargin(0.16);
   dptr->SetRightMargin(0.04);
   dptr->SetTopMargin(0.08);
   dptr->SetBottomMargin(0.25);
   dptr->SetFrameFillStyle(0);
   dptr->SetFrameBorderMode(0);
   
   TH1D *hdptr1__34 = new TH1D("hdptr1__34","",50,0,200);
   hdptr1__34->SetBinContent(1,0.9754253);
   hdptr1__34->SetBinContent(2,1.010983);
   hdptr1__34->SetBinContent(3,1.028526);
   hdptr1__34->SetBinContent(4,1.037003);
   hdptr1__34->SetBinContent(5,1.009601);
   hdptr1__34->SetBinContent(6,1.028141);
   hdptr1__34->SetBinContent(7,1.014026);
   hdptr1__34->SetBinContent(8,0.9871023);
   hdptr1__34->SetBinContent(9,1.070468);
   hdptr1__34->SetBinContent(10,1.091433);
   hdptr1__34->SetBinContent(11,1.087432);
   hdptr1__34->SetBinContent(12,1);
   hdptr1__34->SetBinContent(13,1.023346);
   hdptr1__34->SetBinContent(14,1.01626);
   hdptr1__34->SetBinContent(15,0.9180328);
   hdptr1__34->SetBinContent(16,1.05814);
   hdptr1__34->SetBinContent(17,1.052239);
   hdptr1__34->SetBinContent(18,0.9252336);
   hdptr1__34->SetBinContent(19,1);
   hdptr1__34->SetBinContent(20,0.9384615);
   hdptr1__34->SetBinContent(21,0.9387755);
   hdptr1__34->SetBinContent(22,1.592593);
   hdptr1__34->SetBinContent(23,1.129032);
   hdptr1__34->SetBinContent(24,0.9333333);
   hdptr1__34->SetBinContent(25,1.916667);
   hdptr1__34->SetBinContent(26,0.6842105);
   hdptr1__34->SetBinContent(27,0.6666667);
   hdptr1__34->SetBinContent(28,0.6363636);
   hdptr1__34->SetBinContent(29,1.272727);
   hdptr1__34->SetBinContent(30,1.8);
   hdptr1__34->SetBinContent(31,0.9230769);
   hdptr1__34->SetBinContent(32,1.8);
   hdptr1__34->SetBinContent(33,1.4);
   hdptr1__34->SetBinContent(34,1.333333);
   hdptr1__34->SetBinContent(35,1);
   hdptr1__34->SetBinContent(36,0.4);
   hdptr1__34->SetBinContent(37,0.75);
   hdptr1__34->SetBinContent(39,0.5);
   hdptr1__34->SetBinContent(41,1.5);
   hdptr1__34->SetBinContent(42,1);
   hdptr1__34->SetBinContent(43,0.25);
   hdptr1__34->SetBinContent(44,1.5);
   hdptr1__34->SetBinContent(46,0.2);
   hdptr1__34->SetBinContent(47,0.3333333);
   hdptr1__34->SetBinContent(49,0.5);
   hdptr1__34->SetBinContent(51,1);
   hdptr1__34->SetBinError(1,0.04267608);
   hdptr1__34->SetBinError(2,0.02333712);
   hdptr1__34->SetBinError(3,0.01667667);
   hdptr1__34->SetBinError(4,0.01471465);
   hdptr1__34->SetBinError(5,0.01479415);
   hdptr1__34->SetBinError(6,0.01675593);
   hdptr1__34->SetBinError(7,0.01967454);
   hdptr1__34->SetBinError(8,0.0237105);
   hdptr1__34->SetBinError(9,0.03205502);
   hdptr1__34->SetBinError(10,0.04053865);
   hdptr1__34->SetBinError(11,0.04980771);
   hdptr1__34->SetBinError(12,0.05447347);
   hdptr1__34->SetBinError(13,0.06346947);
   hdptr1__34->SetBinError(14,0.07451823);
   hdptr1__34->SetBinError(15,0.07598138);
   hdptr1__34->SetBinError(16,0.1125239);
   hdptr1__34->SetBinError(17,0.1269459);
   hdptr1__34->SetBinError(18,0.1290255);
   hdptr1__34->SetBinError(19,0.1552301);
   hdptr1__34->SetBinError(20,0.1672939);
   hdptr1__34->SetBinError(21,0.1927288);
   hdptr1__34->SetBinError(22,0.3910548);
   hdptr1__34->SetBinError(23,0.2784604);
   hdptr1__34->SetBinError(24,0.2452512);
   hdptr1__34->SetBinError(25,0.6825368);
   hdptr1__34->SetBinError(26,0.2462727);
   hdptr1__34->SetBinError(27,0.2151657);
   hdptr1__34->SetBinError(28,0.2175611);
   hdptr1__34->SetBinError(29,0.5127964);
   hdptr1__34->SetBinError(30,1.003992);
   hdptr1__34->SetBinError(31,0.3695265);
   hdptr1__34->SetBinError(32,1.003992);
   hdptr1__34->SetBinError(33,0.8197561);
   hdptr1__34->SetBinError(34,0.7200823);
   hdptr1__34->SetBinError(35,0.6324555);
   hdptr1__34->SetBinError(36,0.334664);
   hdptr1__34->SetBinError(37,0.572822);
   hdptr1__34->SetBinError(39,0.6123724);
   hdptr1__34->SetBinError(41,1.369306);
   hdptr1__34->SetBinError(42,1.414214);
   hdptr1__34->SetBinError(43,0.2795085);
   hdptr1__34->SetBinError(44,1.369306);
   hdptr1__34->SetBinError(46,0.219089);
   hdptr1__34->SetBinError(47,0.3849002);
   hdptr1__34->SetBinError(49,0.6123724);
   hdptr1__34->SetBinError(51,0.3162278);
   hdptr1__34->SetMinimum(0);
   hdptr1__34->SetMaximum(2);
   hdptr1__34->SetEntries(166.4223);

   ci = TColor::GetColor("#cc0000");
   hdptr1__34->SetFillColor(ci);
   hdptr1__34->SetLineStyle(0);
   hdptr1__34->SetMarkerStyle(20);
   hdptr1__34->SetMarkerSize(0);
   hdptr1__34->GetXaxis()->SetTitle("p_{T}(ll) [GeV/c]");
   hdptr1__34->GetYaxis()->SetTitle("Ratio");
   hdptr1__34->GetXaxis()->CenterTitle(kTRUE);
   hdptr1__34->GetXaxis()->SetNdivisions(505);
   hdptr1__34->GetXaxis()->SetLabelFont(42);
   hdptr1__34->GetXaxis()->SetLabelOffset(0.007);
   hdptr1__34->GetXaxis()->SetLabelSize(0.05);
   hdptr1__34->GetXaxis()->SetTitleSize(0.11);
   hdptr1__34->GetXaxis()->SetTitleOffset(0.75);
   hdptr1__34->GetXaxis()->SetTitleFont(42);
   hdptr1__34->GetYaxis()->SetLabelFont(42);
   hdptr1__34->GetYaxis()->SetLabelOffset(0.007);
   hdptr1__34->GetYaxis()->SetLabelSize(0.05);
   hdptr1__34->GetYaxis()->SetTitleSize(0.08);
   hdptr1__34->GetYaxis()->SetTitleOffset(0.6);
   hdptr1__34->GetYaxis()->SetTitleFont(42);
   hdptr1__34->GetZaxis()->SetNdivisions(505);
   hdptr1__34->GetZaxis()->SetLabelFont(42);
   hdptr1__34->GetZaxis()->SetLabelOffset(0.007);
   hdptr1__34->GetZaxis()->SetLabelSize(0.05);
   hdptr1__34->GetZaxis()->SetTitleSize(0.06);
   hdptr1__34->GetZaxis()->SetTitleFont(42);
   hdptr1__34->Draw("e2");
   
   TH1D *hdptr2__35 = new TH1D("hdptr2__35","",50,0,200);
   hdptr2__35->SetBinContent(1,0.989603);
   hdptr2__35->SetBinContent(2,0.9691937);
   hdptr2__35->SetBinContent(3,1.041589);
   hdptr2__35->SetBinContent(4,1.053096);
   hdptr2__35->SetBinContent(5,1.033118);
   hdptr2__35->SetBinContent(6,1.016965);
   hdptr2__35->SetBinContent(7,1.054776);
   hdptr2__35->SetBinContent(8,1.045285);
   hdptr2__35->SetBinContent(9,1.057487);
   hdptr2__35->SetBinContent(10,1.059755);
   hdptr2__35->SetBinContent(11,1.098361);
   hdptr2__35->SetBinContent(12,1.090504);
   hdptr2__35->SetBinContent(13,1.009728);
   hdptr2__35->SetBinContent(14,0.9837398);
   hdptr2__35->SetBinContent(15,0.9081967);
   hdptr2__35->SetBinContent(16,1.203488);
   hdptr2__35->SetBinContent(17,1.097015);
   hdptr2__35->SetBinContent(18,1.214953);
   hdptr2__35->SetBinContent(19,1.204819);
   hdptr2__35->SetBinContent(20,0.9538462);
   hdptr2__35->SetBinContent(21,1.22449);
   hdptr2__35->SetBinContent(22,1.888889);
   hdptr2__35->SetBinContent(23,0.7419355);
   hdptr2__35->SetBinContent(24,0.7666667);
   hdptr2__35->SetBinContent(25,2.416667);
   hdptr2__35->SetBinContent(26,0.8947368);
   hdptr2__35->SetBinContent(27,0.75);
   hdptr2__35->SetBinContent(28,0.5909091);
   hdptr2__35->SetBinContent(29,0.8181818);
   hdptr2__35->SetBinContent(30,1.6);
   hdptr2__35->SetBinContent(31,0.4615385);
   hdptr2__35->SetBinContent(32,0.4);
   hdptr2__35->SetBinContent(33,1);
   hdptr2__35->SetBinContent(34,1.5);
   hdptr2__35->SetBinContent(35,0.2);
   hdptr2__35->SetBinContent(36,0.6);
   hdptr2__35->SetBinContent(37,1);
   hdptr2__35->SetBinContent(39,1.5);
   hdptr2__35->SetBinContent(40,10);
   hdptr2__35->SetBinContent(42,1);
   hdptr2__35->SetBinContent(43,0.5);
   hdptr2__35->SetBinContent(44,2.5);
   hdptr2__35->SetBinContent(46,0.2);
   hdptr2__35->SetBinContent(47,0.3333333);
   hdptr2__35->SetBinContent(51,0.9);
   hdptr2__35->SetBinError(1,0.04313908);
   hdptr2__35->SetBinError(2,0.02261104);
   hdptr2__35->SetBinError(3,0.01683619);
   hdptr2__35->SetBinError(4,0.01488684);
   hdptr2__35->SetBinError(5,0.01505276);
   hdptr2__35->SetBinError(6,0.01661864);
   hdptr2__35->SetBinError(7,0.02026796);
   hdptr2__35->SetBinError(8,0.02475391);
   hdptr2__35->SetBinError(9,0.03176004);
   hdptr2__35->SetBinError(10,0.03964236);
   hdptr2__35->SetBinError(11,0.05018824);
   hdptr2__35->SetBinError(12,0.05815798);
   hdptr2__35->SetBinError(13,0.0628332);
   hdptr2__35->SetBinError(14,0.07272258);
   hdptr2__35->SetBinError(15,0.07537922);
   hdptr2__35->SetBinError(16,0.1241687);
   hdptr2__35->SetBinError(17,0.1310251);
   hdptr2__35->SetBinError(18,0.1585879);
   hdptr2__35->SetBinError(19,0.1788992);
   hdptr2__35->SetBinError(20,0.1693276);
   hdptr2__35->SetBinError(21,0.2357734);
   hdptr2__35->SetBinError(22,0.449559);
   hdptr2__35->SetBinError(23,0.2041824);
   hdptr2__35->SetBinError(24,0.2124809);
   hdptr2__35->SetBinError(25,0.8295051);
   hdptr2__35->SetBinError(26,0.2987069);
   hdptr2__35->SetBinError(27,0.2338536);
   hdptr2__35->SetBinError(28,0.2067149);
   hdptr2__35->SetBinError(29,0.3677454);
   hdptr2__35->SetBinError(30,0.9121403);
   hdptr2__35->SetBinError(31,0.2277914);
   hdptr2__35->SetBinError(32,0.334664);
   hdptr2__35->SetBinError(33,0.6324555);
   hdptr2__35->SetBinError(34,0.7905694);
   hdptr2__35->SetBinError(35,0.219089);
   hdptr2__35->SetBinError(36,0.438178);
   hdptr2__35->SetBinError(37,0.7071068);
   hdptr2__35->SetBinError(39,1.369306);
   hdptr2__35->SetBinError(40,10.48809);
   hdptr2__35->SetBinError(42,1.414214);
   hdptr2__35->SetBinError(43,0.4330127);
   hdptr2__35->SetBinError(44,2.09165);
   hdptr2__35->SetBinError(46,0.219089);
   hdptr2__35->SetBinError(47,0.3849002);
   hdptr2__35->SetBinError(51,0.2924038);
   hdptr2__35->SetMinimum(0);
   hdptr2__35->SetMaximum(2);
   hdptr2__35->SetEntries(23.70369);

   ci = TColor::GetColor("#00cc00");
   hdptr2__35->SetFillColor(ci);
   hdptr2__35->SetFillStyle(3001);
   hdptr2__35->SetLineStyle(0);
   hdptr2__35->SetMarkerStyle(20);
   hdptr2__35->SetMarkerSize(0);
   hdptr2__35->GetXaxis()->SetTitle("p_{T}_{ll} [GeV/c]");
   hdptr2__35->GetXaxis()->SetNdivisions(505);
   hdptr2__35->GetXaxis()->SetLabelFont(42);
   hdptr2__35->GetXaxis()->SetLabelOffset(0.007);
   hdptr2__35->GetXaxis()->SetLabelSize(0.05);
   hdptr2__35->GetXaxis()->SetTitleSize(0.06);
   hdptr2__35->GetXaxis()->SetTitleOffset(0.9);
   hdptr2__35->GetXaxis()->SetTitleFont(42);
   hdptr2__35->GetYaxis()->SetLabelFont(42);
   hdptr2__35->GetYaxis()->SetLabelOffset(0.007);
   hdptr2__35->GetYaxis()->SetLabelSize(0.05);
   hdptr2__35->GetYaxis()->SetTitleSize(0.06);
   hdptr2__35->GetYaxis()->SetTitleOffset(1.25);
   hdptr2__35->GetYaxis()->SetTitleFont(42);
   hdptr2__35->GetZaxis()->SetNdivisions(505);
   hdptr2__35->GetZaxis()->SetLabelFont(42);
   hdptr2__35->GetZaxis()->SetLabelOffset(0.007);
   hdptr2__35->GetZaxis()->SetLabelSize(0.05);
   hdptr2__35->GetZaxis()->SetTitleSize(0.06);
   hdptr2__35->GetZaxis()->SetTitleFont(42);
   hdptr2__35->Draw("e2same");

   TLegend *leg2 = new TLegend(0.2,0.67,0.6,0.88,NULL,"brNDC");
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.07);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);
   TLegendEntry *entry2=leg2->AddEntry("hdptr2__35","#mu#mu#gamma / ee#gamma","f");
   entry2->SetLineColor(1);
   entry2->SetLineStyle(1);
   entry2->SetLineWidth(1);
   entry2->SetMarkerColor(1);
   entry2->SetMarkerStyle(21);
   entry2->SetMarkerSize(1);
   entry2=leg2->AddEntry("hdptr1__34","#tau#tau#gamma / ee#gamma","f");
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
   
   dptr->Modified();
   c1_n7->cd();
   c1_n7->Modified();
   c1_n7->cd();
   c1_n7->SetSelected(c1_n7);
   c1_n7->SaveAs("plots/DileptonPt_sherpa.png");
   c1_n7->SaveAs("plots/DileptonPt_sherpa.pdf");
}
