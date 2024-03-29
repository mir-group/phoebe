.. _mlPhononTransport:

Accelerated Lattice Thermal Conductivity Calculations using Machine Learning Force Fields
===========================================================================================

Synopsis
--------

In this tutorial, we will use Phoebe to compute the lattice thermal conductivity of silicon carbide (SiC) with a machine learning force field from `FLARE <https://github.com/mir-group/flare>`_.

This tutorial is based on the thermal conductivity calculation of SiC in `this paper <https://arxiv.org/abs/2203.03824>`_.

.. note::
  Support of phono3py is only available if Phoebe is built using `hdf5`. It's needed for this tutorial.

We will use

* Python (+ LAMMPS) to train a machine learning force field
* Phono3py to generate super cells with different displacements
* LAMMPS to compute forces for those super cells
* Phono3py to compute the force constants
* Phoebe to compute the lattice thermal conductivity.

Step 0: Train a Machine Learning Force Field
--------------------------------------------

Machine learning force field (MLFF) builds a model to predict energy/forces/stress from atomic coordinates, using machine learning algorithms. The model is usually trained from DFT data, and can reach first-principles accuracy while keep the linear scaling with system size in the computational cost. Therefore, the MLFF usually reaches 3-6 orders of magnitude speedup compared to DFT in production.

After a MLFF is trained, we can use it for large-scale molecular dynamics, and phonon transport calculations, where the latter is the subject of this tutorial.

Before any phonon transport calculations, we need to train a MLFF from the first-principles data. The MLFF model we use is `FLARE <https://github.com/mir-group/flare>`_, a Gaussian process regression model with many-body descriptors. The training can be efficiently done through Bayesian active learning, and we direct the users to the `FLARE tutorial <https://colab.research.google.com/drive/1qgGlfu1BlXQgSrnolS4c4AYeZ-2TaX5Y?usp=sharing>`_ to generate a force field before continuing.

For the phonon transport calculations, you need to prepare

* A LAMMPS executable for force calculation compiled with FLARE-LAMMPS code (in this tutorial the executable is named `lmp_mpi`)
* A coefficient file generated by FLARE for LAMMPS (in this tutorial the file is named `lmp.flare`).

Detailed explanation can be found in the tutorial linked above.


Step 1: Phono3py Installation
-----------------------------

See :ref:`phono3pyInstallation`.
Additionally, we recommend you use the unit cell that phonopy selects as the primitive cell. See a note on this in :ref:`forceConstantCalculation`.

Step 2: Relax the Atomic Structure
----------------------------------

As introduced in Step 0, suppose we have trained a good MLFF of SiC, compiled LAMMPS executable (`lmp_mpi`) with FLARE-LAMMPS code, and obtained a coefficient file (`lmp.flare`) for FLARE SiC force field.
Before computing force constants, we need to relax the unit cell with our MLFF. The relaxation is done by LAMMPS, with the Python interface provided by `ASE <https://wiki.fysik.dtu.dk/ase/ase/calculators/lammps.html>`_.

We first set up the ASE calculator of LAMMPS, where the LAMMPS input/output files will be saved in the `tmp` folder, and each LAMMPS calculation will run a relaxation of the cell and atomic positions. We run the relaxation `repeat=5` times to make sure the structure is relaxed sufficiently.

.. code:: python

    import numpy as np
    from ase.io import read, write
    import os, glob, shutil, sys
    from copy import deepcopy
    from ase.calculators.lammpsrun import LAMMPS


    def relax(struc, pot_file, lmp_command, specorder, repeat=5):
        """
        Relax cell and positions with LAMMPS
        """
        # create ASE calculator for LAMMPS
        calc = LAMMPS(
            label=f"tmp",
            keep_tmp_files=True,
            tmp_dir="tmp",
            files=[pot_file],
            specorder=specorder,
        )
        calc.set(
            command=lmp_command,                        # point to lammps executable
            pair_style="flare",                         # lammps command of using flare force field
            pair_coeff=[f"* * {pot_file}"],             # lammps command of using flare coefficients
            minimize="1.0e-4 1.0e-6 100 1000",          # lammps command of relaxing atomic positions
            fix=["1 all box/relax iso 0.0 vmax 0.001"], # lammps command of relaxing cell
        )

        # run relaxation multiple times
        for i in range(repeat):
            atoms = deepcopy(struc)
            atoms.calc = deepcopy(calc)
            atoms.get_forces()

            trj_files = glob.glob("tmp/trj*")
            assert len(trj_files) == 1
            struc = read(trj_files[0], specorder=specorder, format="lammps-dump-binary")
            os.remove(trj_files[0])

        return struc

We can start with a DFT relaxed unit cell structure or one downloaded from a site like the Materials Project. We then relax the structure with the ASE LAMMPS calculator and save the relaxed structure to file (`POSCAR-unitcell`).

.. code:: python

    # Read an atomic structure of unit cell from file (can be a DFT structure or downloaded from online)
    dft_struc = read(f"POSCAR", format="vasp")

    # Relax the structure using LAMMPS and FLARE force field
    struc = relax(
        dft_struc,
        pot_file="lmp.flare",
        lmp_command="./lmp_mpi",
        specorder=["Si", "C"],
        repeat=5,
    )

    # Write the relaxed unit cell into a file named "POSCAR-unitcell"
    write("POSCAR-unitcell", struc, format="vasp")


Step 3: Construct Force Constant Matrices
------------------------------------------

After obtaining a relaxed unit cell, we use the Python interface of phonopy and phono3py to compute force constants. As in the Step 2, we first define an ASE LAMMPS calculator for force calculation in later usage.

.. code:: python

    import numpy as np
    from ase import Atoms
    from ase.calculators.lammpsrun import LAMMPS

    def get_lmp_calc(pot_file, specorder, lmp_command):
        # create ASE calc for LAMMPS
        calc = LAMMPS(
            label=f"tmp",
            keep_tmp_files=True,
            tmp_dir="tmp",
            files=[pot_file],
            specorder=specorder,
        )
        calc.set(
            command=lmp_command,
            pair_style="flare",
            pair_coeff=[f"* * {pot_file}"],
        )
        return calc

Then we use `phono3py Python API <https://phonopy.github.io/phono3py/phono3py-api.html>`_ to make supercells and generate displacements.
Here we use different supercell sizes for 2nd order force constants (6x6x2) and 3rd order force constants (3x3x3). And for 3rd order force constants, we use a cutoff pair distance 2.5A.

.. code:: python

    from phonopy.interface.calculator import read_crystal_structure
    from phono3py import Phono3py
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
    from tqdm import tqdm

    # generate displacements
    unitcell, _ = read_crystal_structure("POSCAR-unitcell", interface_mode='vasp')
    ph3 = Phono3py(
        unitcell,
        supercell_matrix=[3, 3, 3],
        phonon_supercell_matrix=[6, 6, 2],
    )
    ph3.generate_displacements(cutoff_pair_distance=2.5)
    ph3.save("phono3py_disp.yaml")
    print("Generated displacements")

Next, we use ASE LAMMPS calculator to compute forces of the displaced supercells of 2nd order force constants. We put all forces into an array of shape `(n_displacements, n_atoms, 3)` and feed to phono3py.

.. code:: python

    pot_file = "lmp.flare"
    specorder = ["Si", "C"]
    lmp_command = "./lmp_mpi"

    # get forces for FC2
    print("Computing forces for FC2")
    forces = []
    for sc in tqdm(ph3.phonon_supercells_with_displacements):
        atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)
        atoms.calc = get_lmp_calc(pot_file, specorder, lmp_command)
        f = atoms.get_forces()
        forces.append(f)

    # compute 2nd order force constants
    ph3.phonon_forces = np.array(forces)

Then phono3py computes the 2nd order force constants and write to a file ``fc2.hdf5``.

.. code:: python

    print("Computing FC2")
    ph3.produce_fc2()
    write_fc2_to_hdf5(
        ph3.fc2,
        p2s_map=ph3.phonon_primitive.p2s_map,
        physical_unit="eV/angstrom^2",
    )

In the same way, the 3rd order force constants can be generated and written into ``fc3.hdf5``.

.. code:: python

    # get forces for FC3
    print("Computing forces for FC3")
    forces = []
    nat = len(ph3.supercells_with_displacements[0])
    for sc in tqdm(ph3.supercells_with_displacements):
        if sc is not None:
            atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)
            atoms.calc = get_lmp_calc(pot_file, specorder, lmp_command)
            f = atoms.get_forces()
        else:
            f = np.zeros((nat, 3))
        forces.append(f)

    # compute 3rd order force constants
    ph3.forces = np.array(forces)

    print("Computing FC3")
    ph3.produce_fc3()
    write_fc3_to_hdf5(
        ph3.fc3,
        p2s_map=ph3.primitive.p2s_map,
    )

The files ``fc2.hdf5`` and ``fc3.hdf5`` will be used by Phoebe.

If you want to check the phonon calculation, you can find instructions in either the :ref:`bands` or :ref:`harmonic_p3py` tutorials.


Step 4: Calculate Lattice Thermal Conductivity
------------------------------------------------

If this dispersion looks good, we are now ready to move on to phonon transport calculations using Phoebe.
See :ref:`thermalConductivityCalculation` of the :ref:`phononTransport` for instructions on how to use these files to generate lattice thermal conductivity.
