@page DEVDOCS Developer's documentation

This section of the manual covers more technical details useful to a developer.

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

Name variables with sensible and clear names or the code can easily become a mess.

Use two spaces for indentation.
We are typically using CLion to format the source text, using its default code style.





@section Patch Quantum ESPRESSO Patch

**Code modifications**

First, download the repository phoebe-quantum-espresso from the MIR github page.
If necessary, you may need to pull from the official quantum-espresso repository into our fork.
Checkout the branch with the version of QE that we want to patch, for example, v6.6 is stored in the branch `qe-6.6`.
Create a new branch from `qe-6.6`, which we call `patched-qe-6.6`, and modify the source code with the lines of code needed by phoebe.
Alternatively, you may create the new branch `patched-qe-6.6` from an older patched branch, e.g. `patched-qe-6.5`, and pull the new commits of `qe-6.6` into `patched-qe-6.6`.

**Patch creation**

Assume now that the source code in branch `patched-qe-6.6` is ready to be distributed.
To create the patch, we simply save the `git diff` between our branch and the original branch released by the QE community.
Simply type
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
git diff qe-6.6..patched-qe-6.6 > patchfile.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This will create the file `"patchfile.txt"`, containing the patch.
You may inspect this file to make sure it didn't add unnecessary modifications.
Currently only 4 files are modified by the patch.


**Patch application**

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

After the patch is successfully applied, compile QE and use as described in the tutorials.
The patch is tested to work in QE v6.5 and v6.6.
The first patch has been developed with Quantum ESPRESSO v6.6.
The patch thus may or may not work with older versions, we don't provide any support with that.

Note also, the user could also directly download the source code of branch `patched-qe-6.6` and use that.



**Troubleshooting**

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





@section Debugging
Valgrind is a great programming tool for memory debugging, memory leak detection and profiling.
To launch it with phoebe, an example of syntax is:
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c}
valgrind --track-origins=yes --leak-check=full ./build/phoebe -in phoebe.in > phoebe.out &> valgrind.out
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This command launches the default valgrind tool, which is, "memcheck", which helps finding bugs in the code.

**Tip for running valgrind when the code is compiled with OpenMPI**:
Valgrind detects some false positive memory leaks, that originate from the MPI library and that you shouldn't worry about.
OpenMPI provides a file "openmpi-valgrind.supp", that can be used to suppress these wrong memory leaks messages.
The syntax is
~~~~~~~~~~~~~~~~~~~~~~~~~~~{.c }
valgrind --track-origins=yes --leak-check=full --suppressions=/usr/share/openmpi/openmpi-valgrind.supp ./build/phoebe -in phoebe.in > phoebe.out &> valgrind.out
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The "suppressions" flag is only relevant if you are compiling the code with MPI, and the path to the file "openmpi-valgrind.supp" should be modified to match the path on your computer.
