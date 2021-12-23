=======================
Rosetta docking plugin
=======================

This is a **Scipion** plugin that offers different **Rosetta molecular docking
tools**. These tools will make it possible to carry out molecular docking
experiments both **on the surface of proteins**, in order to avoid possible
protein-protein interactions, and **inside proteins**, in protein pockets
related to their biological activity.

Therefore, this plugin allows to use programs from the ROSETTA software suite
within the Scipion framework. **You need to download the Rosetta suite files
before installing the plugin, see section "Binary Files" for details**.

Rosetta is a software suite includes algorithms for computational modeling
and analysis of protein structures. In addition, it has implemented several
docking tools, that we used here
(see `Rosetta home page <https://www.rosettacommons.org/>`_ for details).

Current programs implemented:

    - score
    - make_ray_files
    - DARC


==========================
Install this plugin
==========================

You will need to use `Scipion3 <https://scipion-em.github.io/docs/docs/scipion
-modes/how-to-install.html>`_ to run these protocols.

1. **Binary files**

**Rosetta** binaries will **NOT** be downloaded automatically with the plugin.

The independent download of Rosetta software suite by the user is required
before running the programs.
The installation path is automatically found between all paths. In case of error, 
this path or any other of your preference has to be set in *ROSETTA_HOME* in
*scipion.conf*  file.

We recommend to install the last version of Rosetta, 3.12 .

The steps to download the Rosetta files are the following ones:

    - Go to  `Rosetta software page <https://www.rosettacommons.org/software>`_.
    - Rosetta is available to all non-commercial users for free and to commercial
      users for a fee. To download the software a license of Rosetta Software Suite
      is needed. You can get it
      `here <https://www.rosettacommons.org/software/license-and-download>`_.
    - With the license (an user and a key is sent to the email) you can download the
      file Rosetta <version> source + binaries for Linux - as one bundle.
    - Unzip the package and save it where you prefer.
    - Go to Rosetta/main and it contains the Rosetta source code, database, unit tests
      and integration tests. The source code is located in source/src can be compiled
      with SCons using the following commands (for more details see `getting started with Rosetta
      <https://www.rosettacommons.org/docs/latest/getting_started/Getting-Started#local
      -installation-and-use-of-rosetta>`_), with opencl as extra to allow to use GPU:

    .. code-block::

        cd path/to/Rosetta/main/source
        sudo apt install opencl-dev
        ./scons.py -j<NumOfJobs> mode=[debug/release*] extras=opencl [bin]


2. **AutoDock dependency**

The plugin scipion-chem-autodock is a requirement to run the DARC and make_ray_files programs.
Therefore, even though it is not a compulsory requirement, it is highly recommended:
(https://github.com/scipion-chem/scipion-chem-autodock)

3. **Install the plugin in Scipion**

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or

.. code-block::

    scipion3 installp -p scipion-chem-rosetta


- **Developer's version**

    1. Download repository:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-rosetta.git

    2. Install:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-rosetta --devel


