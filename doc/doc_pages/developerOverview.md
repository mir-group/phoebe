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

@section scheme Code schematics

To have a good overview of the main Phoebe classes and how these work, it's easier to start with a graphical sketch of the PhononTransportApp and its code workflow.

![Fig.1: appScheme] (../images/TransportCode.png)

We refer to the class documentation for a detailed overview of all the classes.

@section Format Code formatting

Use the camelCase notation for variable definitions.
Classes should start with a capital letter.
The code should be formatted according to the clang-format style.
Classes and files should be documented using a Doxygen-compatible style, in order to render nicely in this documentation.