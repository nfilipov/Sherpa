/// Photon Eta plot

void cpeta()
{
//=========Macro generated from canvas: c1_n3/c1_n3
//=========  (Tue Oct  3 15:49:41 2017) by ROOT version6.08/00
   TCanvas *c1_n3 = new TCanvas("c1_n3", "c1_n3",0,45,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n3->Range(0,0,1,1);
   c1_n3->SetFillColor(0);
   c1_n3->SetBorderMode(0);
   c1_n3->SetBorderSize(2);
   c1_n3->SetTickx(1);
   c1_n3->SetTicky(1);
   c1_n3->SetLeftMargin(0.16);
   c1_n3->SetRightMargin(0.04);
   c1_n3->SetTopMargin(0.08);
   c1_n3->SetBottomMargin(0.12);
   c1_n3->SetFrameFillStyle(0);
   c1_n3->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: peta
   TPad *peta = new TPad("peta", "peta",0,0.42,1,1);
   peta->Draw();
   peta->cd();
   peta->Range(0,0,1,1);
   peta->SetFillColor(0);
   peta->SetBorderMode(0);
   peta->SetBorderSize(2);
   peta->SetLogy();
   peta->SetTickx(1);
   peta->SetTicky(1);
   peta->SetLeftMargin(0.16);
   peta->SetRightMargin(0.04);
   peta->SetTopMargin(0.08);
   peta->SetBottomMargin(0);
   peta->SetFrameFillStyle(0);
   peta->SetFrameBorderMode(0);
   
   TH1D *hPhotonEta__11 = new TH1D("hPhotonEta__11","",50,-3,3);
   hPhotonEta__11->SetBinContent(0,129);
   hPhotonEta__11->SetBinContent(1,43);
   hPhotonEta__11->SetBinContent(2,69);
   hPhotonEta__11->SetBinContent(3,189);
   hPhotonEta__11->SetBinContent(4,602);
   hPhotonEta__11->SetBinContent(5,922);
   hPhotonEta__11->SetBinContent(6,1015);
   hPhotonEta__11->SetBinContent(7,1134);
   hPhotonEta__11->SetBinContent(8,1240);
   hPhotonEta__11->SetBinContent(9,1129);
   hPhotonEta__11->SetBinContent(10,1290);
   hPhotonEta__11->SetBinContent(11,1320);
   hPhotonEta__11->SetBinContent(12,1261);
   hPhotonEta__11->SetBinContent(13,1276);
   hPhotonEta__11->SetBinContent(14,1343);
   hPhotonEta__11->SetBinContent(15,1320);
   hPhotonEta__11->SetBinContent(16,1326);
   hPhotonEta__11->SetBinContent(17,1370);
   hPhotonEta__11->SetBinContent(18,1337);
   hPhotonEta__11->SetBinContent(19,1398);
   hPhotonEta__11->SetBinContent(20,1394);
   hPhotonEta__11->SetBinContent(21,1438);
   hPhotonEta__11->SetBinContent(22,1452);
   hPhotonEta__11->SetBinContent(23,1450);
   hPhotonEta__11->SetBinContent(24,1364);
   hPhotonEta__11->SetBinContent(25,1367);
   hPhotonEta__11->SetBinContent(26,1368);
   hPhotonEta__11->SetBinContent(27,1330);
   hPhotonEta__11->SetBinContent(28,1383);
   hPhotonEta__11->SetBinContent(29,1384);
   hPhotonEta__11->SetBinContent(30,1343);
   hPhotonEta__11->SetBinContent(31,1326);
   hPhotonEta__11->SetBinContent(32,1367);
   hPhotonEta__11->SetBinContent(33,1319);
   hPhotonEta__11->SetBinContent(34,1429);
   hPhotonEta__11->SetBinContent(35,1294);
   hPhotonEta__11->SetBinContent(36,1323);
   hPhotonEta__11->SetBinContent(37,1298);
   hPhotonEta__11->SetBinContent(38,1259);
   hPhotonEta__11->SetBinContent(39,1276);
   hPhotonEta__11->SetBinContent(40,1314);
   hPhotonEta__11->SetBinContent(41,1246);
   hPhotonEta__11->SetBinContent(42,1184);
   hPhotonEta__11->SetBinContent(43,1191);
   hPhotonEta__11->SetBinContent(44,1130);
   hPhotonEta__11->SetBinContent(45,1112);
   hPhotonEta__11->SetBinContent(46,965);
   hPhotonEta__11->SetBinContent(47,584);
   hPhotonEta__11->SetBinContent(48,140);
   hPhotonEta__11->SetBinContent(49,66);
   hPhotonEta__11->SetBinContent(50,40);
   hPhotonEta__11->SetBinContent(51,147);
   hPhotonEta__11->SetEntries(55996);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cc0000");
   hPhotonEta__11->SetFillColor(ci);
   hPhotonEta__11->SetLineStyle(0);
   hPhotonEta__11->SetMarkerStyle(20);
   hPhotonEta__11->SetMarkerSize(0);
   hPhotonEta__11->GetXaxis()->SetTitle("#eta_{#gamma}");
   hPhotonEta__11->GetXaxis()->SetNdivisions(505);
   hPhotonEta__11->GetXaxis()->SetLabelFont(42);
   hPhotonEta__11->GetXaxis()->SetLabelOffset(0.007);
   hPhotonEta__11->GetXaxis()->SetLabelSize(0.05);
   hPhotonEta__11->GetXaxis()->SetTitleSize(0.06);
   hPhotonEta__11->GetXaxis()->SetTitleOffset(0.9);
   hPhotonEta__11->GetXaxis()->SetTitleFont(42);
   hPhotonEta__11->GetYaxis()->SetLabelFont(42);
   hPhotonEta__11->GetYaxis()->SetLabelOffset(0.007);
   hPhotonEta__11->GetYaxis()->SetLabelSize(0.05);
   hPhotonEta__11->GetYaxis()->SetTitleSize(0.09);
   hPhotonEta__11->GetYaxis()->SetTitleOffset(0.65);
   hPhotonEta__11->GetYaxis()->SetTitle("Entries");
   hPhotonEta__11->GetYaxis()->SetTitleFont(42);
   hPhotonEta__11->GetZaxis()->SetNdivisions(505);
   hPhotonEta__11->GetZaxis()->SetLabelFont(42);
   hPhotonEta__11->GetZaxis()->SetLabelOffset(0.007);
   hPhotonEta__11->GetZaxis()->SetLabelSize(0.05);
   hPhotonEta__11->GetZaxis()->SetTitleSize(0.06);
   hPhotonEta__11->GetZaxis()->SetTitleFont(42);
   hPhotonEta__11->Draw("e2");
   
   TH1D *hPhotonEta__12 = new TH1D("hPhotonEta__12","",50,-3,3);
   hPhotonEta__12->SetBinContent(0,145);
   hPhotonEta__12->SetBinContent(1,48);
   hPhotonEta__12->SetBinContent(2,75);
   hPhotonEta__12->SetBinContent(3,145);
   hPhotonEta__12->SetBinContent(4,595);
   hPhotonEta__12->SetBinContent(5,1016);
   hPhotonEta__12->SetBinContent(6,1041);
   hPhotonEta__12->SetBinContent(7,1073);
   hPhotonEta__12->SetBinContent(8,1235);
   hPhotonEta__12->SetBinContent(9,1234);
   hPhotonEta__12->SetBinContent(10,1213);
   hPhotonEta__12->SetBinContent(11,1315);
   hPhotonEta__12->SetBinContent(12,1355);
   hPhotonEta__12->SetBinContent(13,1374);
   hPhotonEta__12->SetBinContent(14,1355);
   hPhotonEta__12->SetBinContent(15,1357);
   hPhotonEta__12->SetBinContent(16,1345);
   hPhotonEta__12->SetBinContent(17,1385);
   hPhotonEta__12->SetBinContent(18,1347);
   hPhotonEta__12->SetBinContent(19,1352);
   hPhotonEta__12->SetBinContent(20,1334);
   hPhotonEta__12->SetBinContent(21,1503);
   hPhotonEta__12->SetBinContent(22,1380);
   hPhotonEta__12->SetBinContent(23,1397);
   hPhotonEta__12->SetBinContent(24,1397);
   hPhotonEta__12->SetBinContent(25,1448);
   hPhotonEta__12->SetBinContent(26,1416);
   hPhotonEta__12->SetBinContent(27,1369);
   hPhotonEta__12->SetBinContent(28,1366);
   hPhotonEta__12->SetBinContent(29,1334);
   hPhotonEta__12->SetBinContent(30,1345);
   hPhotonEta__12->SetBinContent(31,1407);
   hPhotonEta__12->SetBinContent(32,1416);
   hPhotonEta__12->SetBinContent(33,1385);
   hPhotonEta__12->SetBinContent(34,1394);
   hPhotonEta__12->SetBinContent(35,1344);
   hPhotonEta__12->SetBinContent(36,1349);
   hPhotonEta__12->SetBinContent(37,1272);
   hPhotonEta__12->SetBinContent(38,1358);
   hPhotonEta__12->SetBinContent(39,1197);
   hPhotonEta__12->SetBinContent(40,1355);
   hPhotonEta__12->SetBinContent(41,1248);
   hPhotonEta__12->SetBinContent(42,1206);
   hPhotonEta__12->SetBinContent(43,1279);
   hPhotonEta__12->SetBinContent(44,1171);
   hPhotonEta__12->SetBinContent(45,1128);
   hPhotonEta__12->SetBinContent(46,921);
   hPhotonEta__12->SetBinContent(47,615);
   hPhotonEta__12->SetBinContent(48,156);
   hPhotonEta__12->SetBinContent(49,64);
   hPhotonEta__12->SetBinContent(50,48);
   hPhotonEta__12->SetBinContent(51,118);
   hPhotonEta__12->SetEntries(56725);

   ci = TColor::GetColor("#00cc00");
   hPhotonEta__12->SetFillColor(ci);
   hPhotonEta__12->SetLineStyle(0);
   hPhotonEta__12->SetMarkerStyle(20);
   hPhotonEta__12->SetMarkerSize(0);
   hPhotonEta__12->GetXaxis()->SetTitle("#eta_{#gamma}");
   hPhotonEta__12->GetXaxis()->SetNdivisions(505);
   hPhotonEta__12->GetXaxis()->SetLabelFont(42);
   hPhotonEta__12->GetXaxis()->SetLabelOffset(0.007);
   hPhotonEta__12->GetXaxis()->SetLabelSize(0.05);
   hPhotonEta__12->GetXaxis()->SetTitleSize(0.06);
   hPhotonEta__12->GetXaxis()->SetTitleOffset(0.9);
   hPhotonEta__12->GetXaxis()->SetTitleFont(42);
   hPhotonEta__12->GetYaxis()->SetLabelFont(42);
   hPhotonEta__12->GetYaxis()->SetLabelOffset(0.007);
   hPhotonEta__12->GetYaxis()->SetLabelSize(0.05);
   hPhotonEta__12->GetYaxis()->SetTitleSize(0.06);
   hPhotonEta__12->GetYaxis()->SetTitleOffset(1.25);
   hPhotonEta__12->GetYaxis()->SetTitleFont(42);
   hPhotonEta__12->GetZaxis()->SetNdivisions(505);
   hPhotonEta__12->GetZaxis()->SetLabelFont(42);
   hPhotonEta__12->GetZaxis()->SetLabelOffset(0.007);
   hPhotonEta__12->GetZaxis()->SetLabelSize(0.05);
   hPhotonEta__12->GetZaxis()->SetTitleSize(0.06);
   hPhotonEta__12->GetZaxis()->SetTitleFont(42);
   hPhotonEta__12->Draw("e2same");
   
   TH1D *hPhotonEta__13 = new TH1D("hPhotonEta__13","",50,-3,3);
   hPhotonEta__13->SetBinContent(0,160);
   hPhotonEta__13->SetBinContent(1,32);
   hPhotonEta__13->SetBinContent(2,74);
   hPhotonEta__13->SetBinContent(3,153);
   hPhotonEta__13->SetBinContent(4,585);
   hPhotonEta__13->SetBinContent(5,937);
   hPhotonEta__13->SetBinContent(6,1058);
   hPhotonEta__13->SetBinContent(7,1094);
   hPhotonEta__13->SetBinContent(8,1221);
   hPhotonEta__13->SetBinContent(9,1175);
   hPhotonEta__13->SetBinContent(10,1209);
   hPhotonEta__13->SetBinContent(11,1143);
   hPhotonEta__13->SetBinContent(12,1205);
   hPhotonEta__13->SetBinContent(13,1312);
   hPhotonEta__13->SetBinContent(14,1326);
   hPhotonEta__13->SetBinContent(15,1316);
   hPhotonEta__13->SetBinContent(16,1268);
   hPhotonEta__13->SetBinContent(17,1291);
   hPhotonEta__13->SetBinContent(18,1383);
   hPhotonEta__13->SetBinContent(19,1310);
   hPhotonEta__13->SetBinContent(20,1339);
   hPhotonEta__13->SetBinContent(21,1335);
   hPhotonEta__13->SetBinContent(22,1356);
   hPhotonEta__13->SetBinContent(23,1292);
   hPhotonEta__13->SetBinContent(24,1379);
   hPhotonEta__13->SetBinContent(25,1362);
   hPhotonEta__13->SetBinContent(26,1363);
   hPhotonEta__13->SetBinContent(27,1441);
   hPhotonEta__13->SetBinContent(28,1355);
   hPhotonEta__13->SetBinContent(29,1409);
   hPhotonEta__13->SetBinContent(30,1314);
   hPhotonEta__13->SetBinContent(31,1374);
   hPhotonEta__13->SetBinContent(32,1384);
   hPhotonEta__13->SetBinContent(33,1294);
   hPhotonEta__13->SetBinContent(34,1309);
   hPhotonEta__13->SetBinContent(35,1275);
   hPhotonEta__13->SetBinContent(36,1282);
   hPhotonEta__13->SetBinContent(37,1364);
   hPhotonEta__13->SetBinContent(38,1248);
   hPhotonEta__13->SetBinContent(39,1259);
   hPhotonEta__13->SetBinContent(40,1207);
   hPhotonEta__13->SetBinContent(41,1219);
   hPhotonEta__13->SetBinContent(42,1140);
   hPhotonEta__13->SetBinContent(43,1100);
   hPhotonEta__13->SetBinContent(44,1147);
   hPhotonEta__13->SetBinContent(45,1087);
   hPhotonEta__13->SetBinContent(46,894);
   hPhotonEta__13->SetBinContent(47,550);
   hPhotonEta__13->SetBinContent(48,165);
   hPhotonEta__13->SetBinContent(49,65);
   hPhotonEta__13->SetBinContent(50,39);
   hPhotonEta__13->SetBinContent(51,111);
   hPhotonEta__13->SetMinimum(10);
   hPhotonEta__13->SetMaximum(5000);
   hPhotonEta__13->SetEntries(54710);

   ci = TColor::GetColor("#3399ff");
   hPhotonEta__13->SetFillColor(ci);
   hPhotonEta__13->SetLineStyle(0);
   hPhotonEta__13->SetMarkerStyle(20);
   hPhotonEta__13->SetMarkerSize(0);
   hPhotonEta__13->GetXaxis()->SetTitle("#eta_{#gamma}");
   hPhotonEta__13->GetXaxis()->SetNdivisions(505);
   hPhotonEta__13->GetXaxis()->SetLabelFont(42);
   hPhotonEta__13->GetXaxis()->SetLabelOffset(0.007);
   hPhotonEta__13->GetXaxis()->SetLabelSize(0.05);
   hPhotonEta__13->GetXaxis()->SetTitleSize(0.06);
   hPhotonEta__13->GetXaxis()->SetTitleOffset(0.9);
   hPhotonEta__13->GetXaxis()->SetTitleFont(42);
   hPhotonEta__13->GetYaxis()->SetLabelFont(42);
   hPhotonEta__13->GetYaxis()->SetLabelOffset(0.007);
   hPhotonEta__13->GetYaxis()->SetLabelSize(0.05);
   hPhotonEta__13->GetYaxis()->SetTitleSize(0.06);
   hPhotonEta__13->GetYaxis()->SetTitleOffset(1.25);
   hPhotonEta__13->GetYaxis()->SetTitleFont(42);
   hPhotonEta__13->GetZaxis()->SetNdivisions(505);
   hPhotonEta__13->GetZaxis()->SetLabelFont(42);
   hPhotonEta__13->GetZaxis()->SetLabelOffset(0.007);
   hPhotonEta__13->GetZaxis()->SetLabelSize(0.05);
   hPhotonEta__13->GetZaxis()->SetTitleSize(0.06);
   hPhotonEta__13->GetZaxis()->SetTitleFont(42);
   hPhotonEta__13->Draw("e2same");
   
   TLegend *leg = new TLegend(0.4,0.2,0.7,0.5,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.08);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hPhotonEta__13","pp #rightarrow ee#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hPhotonEta__12","pp #rightarrow #mu#mu#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hPhotonEta__11","pp #rightarrow #tau#tau#gamma","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   peta->Modified();
   c1_n3->cd();
  
// ------------>Primitives in pad: petar
   TPad *petar = new TPad("petar", "petar",0,0.0,1,0.42);
   petar->Draw();
   petar->cd();
   petar->Range(0,0,1,1);
   petar->SetFillColor(0);
   petar->SetBorderMode(0);
   petar->SetBorderSize(2);
   petar->SetTickx(1);
   petar->SetTicky(1);
   petar->SetLeftMargin(0.16);
   petar->SetRightMargin(0.04);
   petar->SetTopMargin(0.08);
   petar->SetBottomMargin(0.25);
   petar->SetFrameFillStyle(0);
   petar->SetFrameBorderMode(0);
   
   TH1D *hpetar1__14 = new TH1D("hpetar1__14","",50,-3,3);
   hpetar1__14->SetBinContent(0,0.80625);
   hpetar1__14->SetBinContent(1,1.34375);
   hpetar1__14->SetBinContent(2,0.9324324);
   hpetar1__14->SetBinContent(3,1.235294);
   hpetar1__14->SetBinContent(4,1.02906);
   hpetar1__14->SetBinContent(5,0.9839915);
   hpetar1__14->SetBinContent(6,0.9593573);
   hpetar1__14->SetBinContent(7,1.036563);
   hpetar1__14->SetBinContent(8,1.015561);
   hpetar1__14->SetBinContent(9,0.9608511);
   hpetar1__14->SetBinContent(10,1.066998);
   hpetar1__14->SetBinContent(11,1.154856);
   hpetar1__14->SetBinContent(12,1.046473);
   hpetar1__14->SetBinContent(13,0.972561);
   hpetar1__14->SetBinContent(14,1.012821);
   hpetar1__14->SetBinContent(15,1.00304);
   hpetar1__14->SetBinContent(16,1.045741);
   hpetar1__14->SetBinContent(17,1.061193);
   hpetar1__14->SetBinContent(18,0.966739);
   hpetar1__14->SetBinContent(19,1.067176);
   hpetar1__14->SetBinContent(20,1.041075);
   hpetar1__14->SetBinContent(21,1.077154);
   hpetar1__14->SetBinContent(22,1.070796);
   hpetar1__14->SetBinContent(23,1.122291);
   hpetar1__14->SetBinContent(24,0.9891226);
   hpetar1__14->SetBinContent(25,1.003671);
   hpetar1__14->SetBinContent(26,1.003668);
   hpetar1__14->SetBinContent(27,0.9229702);
   hpetar1__14->SetBinContent(28,1.020664);
   hpetar1__14->SetBinContent(29,0.9822569);
   hpetar1__14->SetBinContent(30,1.02207);
   hpetar1__14->SetBinContent(31,0.9650655);
   hpetar1__14->SetBinContent(32,0.9877168);
   hpetar1__14->SetBinContent(33,1.01932);
   hpetar1__14->SetBinContent(34,1.091673);
   hpetar1__14->SetBinContent(35,1.014902);
   hpetar1__14->SetBinContent(36,1.031981);
   hpetar1__14->SetBinContent(37,0.9516129);
   hpetar1__14->SetBinContent(38,1.008814);
   hpetar1__14->SetBinContent(39,1.013503);
   hpetar1__14->SetBinContent(40,1.08865);
   hpetar1__14->SetBinContent(41,1.022149);
   hpetar1__14->SetBinContent(42,1.038596);
   hpetar1__14->SetBinContent(43,1.082727);
   hpetar1__14->SetBinContent(44,0.9851787);
   hpetar1__14->SetBinContent(45,1.022999);
   hpetar1__14->SetBinContent(46,1.079418);
   hpetar1__14->SetBinContent(47,1.061818);
   hpetar1__14->SetBinContent(48,0.8484848);
   hpetar1__14->SetBinContent(49,1.015385);
   hpetar1__14->SetBinContent(50,1.025641);
   hpetar1__14->SetBinContent(51,1.324324);
   hpetar1__14->SetBinError(0,0.09540339);
   hpetar1__14->SetBinError(1,0.3137183);
   hpetar1__14->SetBinError(2,0.1560432);
   hpetar1__14->SetBinError(3,0.1343405);
   hpetar1__14->SetBinError(4,0.05974338);
   hpetar1__14->SetBinError(5,0.04564527);
   hpetar1__14->SetBinError(6,0.04215064);
   hpetar1__14->SetBinError(7,0.04392767);
   hpetar1__14->SetBinError(8,0.04094427);
   hpetar1__14->SetBinError(9,0.04004344);
   hpetar1__14->SetBinError(10,0.04271088);
   hpetar1__14->SetBinError(11,0.04666055);
   hpetar1__14->SetBinError(12,0.04215737);
   hpetar1__14->SetBinError(13,0.03823902);
   hpetar1__14->SetBinError(14,0.03921001);
   hpetar1__14->SetBinError(15,0.03907292);
   hpetar1__14->SetBinError(16,0.04107503);
   hpetar1__14->SetBinError(17,0.04116171);
   hpetar1__14->SetBinError(18,0.03707805);
   hpetar1__14->SetBinError(19,0.04103656);
   hpetar1__14->SetBinError(20,0.03983642);
   hpetar1__14->SetBinError(21,0.04093855);
   hpetar1__14->SetBinError(22,0.04043825);
   hpetar1__14->SetBinError(23,0.04293623);
   hpetar1__14->SetBinError(24,0.03777232);
   hpetar1__14->SetBinError(25,0.03842557);
   hpetar1__14->SetBinError(26,0.0384114);
   hpetar1__14->SetBinError(27,0.03509523);
   hpetar1__14->SetBinError(28,0.03901386);
   hpetar1__14->SetBinError(29,0.03717381);
   hpetar1__14->SetBinError(30,0.03965891);
   hpetar1__14->SetBinError(31,0.03715126);
   hpetar1__14->SetBinError(32,0.03766391);
   hpetar1__14->SetBinError(33,0.03988326);
   hpetar1__14->SetBinError(34,0.04176605);
   hpetar1__14->SetBinError(35,0.04004828);
   hpetar1__14->SetBinError(36,0.04044378);
   hpetar1__14->SetBinError(37,0.03689944);
   hpetar1__14->SetBinError(38,0.04029658);
   hpetar1__14->SetBinError(39,0.04026018);
   hpetar1__14->SetBinError(40,0.04340334);
   hpetar1__14->SetBinError(41,0.04117768);
   hpetar1__14->SetBinError(42,0.04309597);
   hpetar1__14->SetBinError(43,0.04527718);
   hpetar1__14->SetBinError(44,0.04129292);
   hpetar1__14->SetBinError(45,0.04363356);
   hpetar1__14->SetBinError(46,0.05010685);
   hpetar1__14->SetBinError(47,0.06309122);
   hpetar1__14->SetBinError(48,0.09749627);
   hpetar1__14->SetBinError(49,0.1774343);
   hpetar1__14->SetBinError(50,0.2308058);
   hpetar1__14->SetBinError(51,0.1665268);
   hpetar1__14->SetMinimum(0);
   hpetar1__14->SetMaximum(2);
   hpetar1__14->SetEntries(8479.476);

   ci = TColor::GetColor("#cc0000");
   hpetar1__14->SetFillColor(ci);
   hpetar1__14->SetLineStyle(0);
   hpetar1__14->SetMarkerStyle(20);
   hpetar1__14->SetMarkerSize(0);
   hpetar1__14->GetXaxis()->SetTitle("#eta_{#gamma}");
   hpetar1__14->GetXaxis()->CenterTitle(kTRUE);
   hpetar1__14->GetXaxis()->SetNdivisions(505);
   hpetar1__14->GetXaxis()->SetLabelFont(42);
   hpetar1__14->GetXaxis()->SetLabelOffset(0.007);
   hpetar1__14->GetXaxis()->SetLabelSize(0.05);
   hpetar1__14->GetXaxis()->SetTitleSize(0.12);
   hpetar1__14->GetXaxis()->SetTitleOffset(0.65);
   hpetar1__14->GetXaxis()->SetTitleFont(42);
   hpetar1__14->GetYaxis()->SetTitle("Ratio");
   hpetar1__14->GetYaxis()->CenterTitle(kTRUE);
   hpetar1__14->GetXaxis()->SetNdivisions(505);
   hpetar1__14->GetYaxis()->SetLabelFont(42);   
   hpetar1__14->GetYaxis()->SetLabelOffset(0.007);
   hpetar1__14->GetYaxis()->SetLabelSize(0.05);
   hpetar1__14->GetYaxis()->SetTitleSize(0.1);
   hpetar1__14->GetYaxis()->SetTitleOffset(0.65);
   hpetar1__14->GetYaxis()->SetTitleFont(42);
   hpetar1__14->GetZaxis()->SetNdivisions(505);
   hpetar1__14->GetZaxis()->SetLabelFont(42);
   hpetar1__14->GetZaxis()->SetLabelOffset(0.007);
   hpetar1__14->GetZaxis()->SetLabelSize(0.05);
   hpetar1__14->GetZaxis()->SetTitleSize(0.06);
   hpetar1__14->GetZaxis()->SetTitleFont(42);
   hpetar1__14->Draw("e2");
   
   TH1D *hpetar2__15 = new TH1D("hpetar2__15","",50,-3,3);
   hpetar2__15->SetBinContent(0,0.90625);
   hpetar2__15->SetBinContent(1,1.5);
   hpetar2__15->SetBinContent(2,1.013514);
   hpetar2__15->SetBinContent(3,0.9477124);
   hpetar2__15->SetBinContent(4,1.017094);
   hpetar2__15->SetBinContent(5,1.084312);
   hpetar2__15->SetBinContent(6,0.9839319);
   hpetar2__15->SetBinContent(7,0.9808044);
   hpetar2__15->SetBinContent(8,1.011466);
   hpetar2__15->SetBinContent(9,1.050213);
   hpetar2__15->SetBinContent(10,1.003309);
   hpetar2__15->SetBinContent(11,1.150481);
   hpetar2__15->SetBinContent(12,1.124481);
   hpetar2__15->SetBinContent(13,1.047256);
   hpetar2__15->SetBinContent(14,1.02187);
   hpetar2__15->SetBinContent(15,1.031155);
   hpetar2__15->SetBinContent(16,1.060726);
   hpetar2__15->SetBinContent(17,1.072812);
   hpetar2__15->SetBinContent(18,0.9739696);
   hpetar2__15->SetBinContent(19,1.032061);
   hpetar2__15->SetBinContent(20,0.9962659);
   hpetar2__15->SetBinContent(21,1.125843);
   hpetar2__15->SetBinContent(22,1.017699);
   hpetar2__15->SetBinContent(23,1.081269);
   hpetar2__15->SetBinContent(24,1.013053);
   hpetar2__15->SetBinContent(25,1.063142);
   hpetar2__15->SetBinContent(26,1.038885);
   hpetar2__15->SetBinContent(27,0.9500347);
   hpetar2__15->SetBinContent(28,1.008118);
   hpetar2__15->SetBinContent(29,0.9467708);
   hpetar2__15->SetBinContent(30,1.023592);
   hpetar2__15->SetBinContent(31,1.024017);
   hpetar2__15->SetBinContent(32,1.023121);
   hpetar2__15->SetBinContent(33,1.070325);
   hpetar2__15->SetBinContent(34,1.064935);
   hpetar2__15->SetBinContent(35,1.054118);
   hpetar2__15->SetBinContent(36,1.052262);
   hpetar2__15->SetBinContent(37,0.9325513);
   hpetar2__15->SetBinContent(38,1.088141);
   hpetar2__15->SetBinContent(39,0.9507546);
   hpetar2__15->SetBinContent(40,1.122618);
   hpetar2__15->SetBinContent(41,1.02379);
   hpetar2__15->SetBinContent(42,1.057895);
   hpetar2__15->SetBinContent(43,1.162727);
   hpetar2__15->SetBinContent(44,1.020924);
   hpetar2__15->SetBinContent(45,1.037718);
   hpetar2__15->SetBinContent(46,1.030201);
   hpetar2__15->SetBinContent(47,1.118182);
   hpetar2__15->SetBinContent(48,0.9454545);
   hpetar2__15->SetBinContent(49,0.9846154);
   hpetar2__15->SetBinContent(50,1.230769);
   hpetar2__15->SetBinContent(51,1.063063);
   hpetar2__15->SetBinError(0,0.1039092);
   hpetar2__15->SetBinError(1,0.3423266);
   hpetar2__15->SetBinError(2,0.1660643);
   hpetar2__15->SetBinError(3,0.1098386);
   hpetar2__15->SetBinError(4,0.05921963);
   hpetar2__15->SetBinError(5,0.04911211);
   hpetar2__15->SetBinError(6,0.04295395);
   hpetar2__15->SetBinError(7,0.04214085);
   hpetar2__15->SetBinError(8,0.04082011);
   hpetar2__15->SetBinError(9,0.04280743);
   hpetar2__15->SetBinError(10,0.0407735);
   hpetar2__15->SetBinError(11,0.0465248);
   hpetar2__15->SetBinError(12,0.04452553);
   hpetar2__15->SetBinError(13,0.04042459);
   hpetar2__15->SetBinError(14,0.03947324);
   hpetar2__15->SetBinError(15,0.03989381);
   hpetar2__15->SetBinError(16,0.04151949);
   hpetar2__15->SetBinError(17,0.04150291);
   hpetar2__15->SetBinError(18,0.0372848);
   hpetar2__15->SetBinError(19,0.04001155);
   hpetar2__15->SetBinError(20,0.03853954);
   hpetar2__15->SetBinError(21,0.04234126);
   hpetar2__15->SetBinError(22,0.0389142);
   hpetar2__15->SetBinError(23,0.04173494);
   hpetar2__15->SetBinError(24,0.03845577);
   hpetar2__15->SetBinError(25,0.04013024);
   hpetar2__15->SetBinError(26,0.03942141);
   hpetar2__15->SetBinError(27,0.03585576);
   hpetar2__15->SetBinError(28,0.03865277);
   hpetar2__15->SetBinError(29,0.03616799);
   hpetar2__15->SetBinError(30,0.03970336);
   hpetar2__15->SetBinError(31,0.03883895);
   hpetar2__15->SetBinError(32,0.03867288);
   hpetar2__15->SetBinError(33,0.04138184);
   hpetar2__15->SetBinError(34,0.04098689);
   hpetar2__15->SetBinError(35,0.04120995);
   hpetar2__15->SetBinError(36,0.04104255);
   hpetar2__15->SetBinError(37,0.03634919);
   hpetar2__15->SetBinError(38,0.04266927);
   hpetar2__15->SetBinError(39,0.03838156);
   hpetar2__15->SetBinError(40,0.04443225);
   hpetar2__15->SetBinError(41,0.04122743);
   hpetar2__15->SetBinError(42,0.04369989);
   hpetar2__15->SetBinError(43,0.04781272);
   hpetar2__15->SetBinError(44,0.04241213);
   hpetar2__15->SetBinError(45,0.04410594);
   hpetar2__15->SetBinError(46,0.04836841);
   hpetar2__15->SetBinError(47,0.06562306);
   hpetar2__15->SetBinError(48,0.1055818);
   hpetar2__15->SetBinError(49,0.1733863);
   hpetar2__15->SetBinError(50,0.2653282);
   hpetar2__15->SetBinError(51,0.1405641);
   hpetar2__15->SetMinimum(0);
   hpetar2__15->SetMaximum(2);
   hpetar2__15->SetEntries(7880.083);

   ci = TColor::GetColor("#00cc00");
   hpetar2__15->SetFillColor(ci);
   hpetar2__15->SetFillStyle(3001);
   hpetar2__15->SetLineStyle(0);
   hpetar2__15->SetMarkerStyle(20);
   hpetar2__15->SetMarkerSize(0);
   hpetar2__15->GetXaxis()->SetTitle("#eta_{#gamma}");
   hpetar2__15->GetXaxis()->SetNdivisions(505);
   hpetar2__15->GetXaxis()->SetLabelFont(42);
   hpetar2__15->GetXaxis()->SetLabelOffset(0.007);
   hpetar2__15->GetXaxis()->SetLabelSize(0.05);
   hpetar2__15->GetXaxis()->SetTitleSize(0.06);
   hpetar2__15->GetXaxis()->SetTitleOffset(0.9);
   hpetar2__15->GetXaxis()->SetTitleFont(42);
   hpetar2__15->GetYaxis()->SetLabelFont(42);
   hpetar2__15->GetYaxis()->SetLabelOffset(0.007);
   hpetar2__15->GetYaxis()->SetLabelSize(0.05);
   hpetar2__15->GetYaxis()->SetTitleSize(0.06);
   hpetar2__15->GetYaxis()->SetTitleOffset(1.25);
   hpetar2__15->GetYaxis()->SetTitleFont(42);
   hpetar2__15->GetZaxis()->SetNdivisions(505);
   hpetar2__15->GetZaxis()->SetLabelFont(42);
   hpetar2__15->GetZaxis()->SetLabelOffset(0.007);
   hpetar2__15->GetZaxis()->SetLabelSize(0.05);
   hpetar2__15->GetZaxis()->SetTitleSize(0.06);
   hpetar2__15->GetZaxis()->SetTitleFont(42);
   hpetar2__15->Draw("e2same");

   TLegend *leg2 = new TLegend(0.4,0.66,0.7,0.9,NULL,"brNDC");
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.08);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);
   TLegendEntry *entry2=leg2->AddEntry("hpetar2__15","#mu#mu#gamma / ee#gamma","f");
   entry2->SetLineColor(1);
   entry2->SetLineStyle(1);
   entry2->SetLineWidth(1);
   entry2->SetMarkerColor(1);
   entry2->SetMarkerStyle(21);
   entry2->SetMarkerSize(1);
   entry2=leg2->AddEntry("hpetar1__14","#tau#tau#gamma / ee#gamma","f");
   entry2->SetLineColor(1);
   entry2->SetLineStyle(1);
   entry2->SetLineWidth(1);
   entry2->SetMarkerColor(1);
   entry2->SetMarkerStyle(21);
   entry2->SetMarkerSize(1);
   leg2->Draw();
   TLine *line;
   line = new TLine(-3,1,3,1);   
   line->Draw();
   
   petar->Modified();
   c1_n3->cd();
   c1_n3->Modified();
   c1_n3->cd();
   c1_n3->SetSelected(c1_n3);

   c1_n3->SaveAs("plots/PhotonEta_sherpa.png");
   c1_n3->SaveAs("plots/PhotonEta_sherpa.pdf");
}
