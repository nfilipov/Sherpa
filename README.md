# Sherpa

This package contains:

## SherpaKinematics.C
This code contains several void methods for building the physical objects, plotting their kinematics, overlaying them at will, overlaying between generators. The following serves as description of the various methods found in the code:
- SherpaKinematics() : This method takes Sherpa samples and builds objects (lepton, photon, lepton-lepton, lepton-photon, triobjects). Some cuts are applied. These objects are saved to blank histograms that you can later plot in a nicer way.
- overlays()
- draw_overlay(TH1D*, TH1D* , string, string , string)
- drawPlots();
- plotMG();

## Style macros
CMS_lumi.h, CMS_lumi.C, tdrstyle.C are used in conjunction to produce pretty plots complying with the CMS plotting guidelines (https://twiki.cern.ch/twiki/bin/view/CMS/Internal/FigGuidelines)

## "plots" directory
This is where the outputs of the methods of SherpaKinematics would go. Outputs can be either saved as .C, .pdf or .png formats, depending on needs. Looking for instances of the SaveAs("....") command in the SherpaKinematics macro should clarify this.
