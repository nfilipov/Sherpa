# Sherpa

Get the package:
```git clone https://github.com/nfilipov/Sherpa/ -d Sherpa```

## How to use SherpaKinematics.C?

```root -l```

```.L SherpaKinematics.C+```

```SherpaKinematics()``` 

This will save the selected events to raw histos. The other functions are for comparing MG and Sherpa, or just plotting things nicely.

```overlays()``` for example will plot madgraph and sherpa plots together and save them to /plots. Don't try this if you haven't got a /plots subdirectory in your work directory.

This package contains the following:

## SherpaKinematics.C
This code contains several void methods for building the physical objects, plotting their kinematics, overlaying them at will, overlaying between generators. The following serves as description of the various methods found in the code:
- SherpaKinematics() : This method takes Sherpa samples and builds objects (lepton, photon, lepton-lepton, lepton-photon, triobjects). Some cuts are applied. These objects are saved to blank histograms that you can later plot in a nicer way.
- plotMG() : Opens the Madgraph samples and builds objects. Same principle as SherpaKinematics() but for Madgraph. These objects are saved to blank histograms that you can later plot in a nicer way.
- overlays() : Opens the Madgraph output _histo_ file (where e,mu,tau kinematics are all together) and the Sherpa output _histo_ files to overlay kinematics together in nice plots.
- draw_overlay(TH1D*, TH1D* , string, string , string) : takes into arguments the raw histograms (for example, Sherpa and Madgraph's triobject pt in the mu channel and overlays the plots together. A ratio is computed automatically. The plots are saved as .C .pdf and .png into the /plots directory. Some basic style is applied; can be edited by editing the .C version of the plot. *No need to call this from within the ROOT session, the _overlays()_ macro will take care of it*
- drawPlots() : This is for plotting Sherpa histograms with style.


## Style macros
CMS_lumi.h, CMS_lumi.C, tdrstyle.C are used in conjunction to produce pretty plots complying with the CMS plotting guidelines (https://twiki.cern.ch/twiki/bin/view/CMS/Internal/FigGuidelines)

## "plots" directory
This is where the outputs of the methods of SherpaKinematics would go. Outputs can be either saved as .C, .pdf or .png formats, depending on needs. Looking for instances of the SaveAs("....") command in the SherpaKinematics macro should clarify this.

