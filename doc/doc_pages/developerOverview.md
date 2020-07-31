@page DEVDOCS Developer's documentation

The main program does very few things: intializes the paralllel environment, loads a subprogram (App) using an AppFactory, and launches it. Currently, supported apps are:
* PhononTransportApp
* ElectronWannierTransportApp
* ElectronEPATransportApp
* PhononDosApp
* ElectronWannierDosApp
* ElectronFourierDosApp
* PhononBandsApp
* ElectronWannierBandsApp
* ElectronFourierBandsApp
* ElectronPolarizationApp

To actually have a good overview of what are the main classes in phoebe and how these works, it's easier to start with a graphical sketch of the PhononTransportApp and its workflow.

![Fig.1: appScheme] (../images/TransportCode.png)

We refer to the class documentation for a detailed overview of all the classes.

We try to use camelCase notation for variable definitions. Classes should start with a capital letter. In general, the code should be formatted according to the clang-format style. Class and file documentation should follow Doxygen-compatible style, in order to render nicely in this page.