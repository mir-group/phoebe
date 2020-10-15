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





@section Patch Quantum ESPRESSO Patch
Quantum espresso patch creation
We download Quantum ESPRESSO v6.6 (QE) and developed the modification in the Fortran code.

To create the patch, do the modifications of QE in a QE git repository.
Make sure you checkout the branch with the QE point release that we want to patch.
Once the code has been modified, compiles and works, create the patch by typing (from the root folder of the QE repository)

~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
git diff reference_branch..modified_branch > patchfile.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

and the patch file `"patchfile.txt"` will be created, containing the patch.

Alternatively, suppose we have a folder with an unmodified reference copy of QE `"./q-e-unmodified"`, and a modified copy of QE in `"./q-e-modified"`.
In the base folder that contains the two copies of QE, execute the command:

~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
diff -Naur q-e-unmodified q-e-modified > patchfile.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command will create a file `"patchfile.txt"` containing the patch.
Note however that diff will also create a difference between binary files: it is therefore recommended to do a make clean before applying the patch.

To apply the patch, go to a folder with a copy of quantum espresso to be patched. In that folder, copy `"patchfile.txt"` and execute the command

~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
patch -p1 < patchfile.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command executes successfully if the output looks like:

~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
patching file PHonon/PH/elphon.f90
patching file PHonon/PH/phq_readin.f90
Hunk #1 succeeded at 911 (offset 4 lines).
patching file PHonon/PH/run_nscf.f90
patching file PW/src/c_bands.f90
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The patch fails if the output looks like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
patching file PHonon/PH/elphon.f90
Hunk #1 succeeded at 1444 (offset -85 lines).
Hunk #2 succeeded at 1458 (offset -85 lines).
patching file PHonon/PH/phq_readin.f90
Hunk #1 FAILED at 907.
1 out of 1 hunk FAILED -- saving rejects to file PHonon/PH/phq_readin.f90.rej
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case, one should manually open the `"*.f90.rej"` file and manually apply the modification in the file `"*.f90"`.
Note that the patch consists in added code. In the file c_bands.f90, we simply add a new subroutine and its call. In the phonon code, we modify the behavior of the 'epa' calculation so that data are written in phoebe format.
Generally, we only added new lines to the code or added new variables, so, as a guideline to fix a failed patch, add the extra lines/variables that appear in the diff. Do so for every file where the patch has failed.

After the patch is successfully applied, compile QE and use as described in the tutorials.

The patch is tested to work in QE v6.5 and v6.6.
