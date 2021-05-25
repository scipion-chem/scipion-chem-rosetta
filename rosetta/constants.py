# **************************************************************************
# *
# * Authors:
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


ROSETTA_HOME = 'ROSETTA_HOME' #Acceso by Plugin.getHome()

# Supported Versions
V3_12 = 'V3_12'

# Path to the binaries inside the ROSETTA folder
ROSETTA_BINARIES_PATH = "main/source/bin"
ROSETTA_DATABASE_PATH = "main/database"
ROSETTA_PARAMS_PATH = "main/source/scripts/python/public"


# Name of programs for linux
SCORE = 'score.static.linuxgccrelease'  # rescores PDBs and silent files, extracts, PDBs from silent files,
                                        # assembles PDBs into silent files.

PARAMS_FILE = 'molfile_to_params.py'
BATCH_PARAMS_FILE = 'batch_molfile_to_params.py'

MAKE_RAY_FILES = 'make_ray_files.static.linuxgccrelease'  # create a ray file to map the pocket or interface
MAKE_RAY_FILES_GPU = 'make_ray_files.opencl.linuxgccrelease'  # create a ray file to map the pocket or interface
                                                              # using GPU

DARC = 'DARC.static.linuxgccrelease'      # run DARC
DARC_GPU = 'DARC.opencl.linuxgccrelease'  # run DARC with GPU
