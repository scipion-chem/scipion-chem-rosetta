# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
This protocol uses a Rosetta suite program (named make_ray_files) to generate
a RAY file for the input protein. To generate this ray-file we need to input
the protein in PDB format and specify a target residue (or more than one)
at the interface.

The output will be a file named ray_<PDBname>_0001_<TargetResidue>.txt
"""

#IMPORT
from pyworkflow.protocol import params
from pyworkflow.protocol.params import (LEVEL_ADVANCED, GPU_LIST)
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pyworkflow.utils import createLink

import shutil
import os
import re

from rosetta import Plugin
from rosetta.constants import *
from rosetta.objects import RaysProtein, RaysStruct, GridAGD
from ..convert import adt2agdGrid



class Rosetta_make_rayFile(EMProtocol):
    """
    This protocol uses a Rosetta suite program (named make_ray_files) to generate
    a RAY file for the input protein. To generate this ray-file we need to input
    the protein in PDB format and specify a target residue (or more than one)
    at the interface.

    The output will be a file named ray_<PDBname>_0001_<TargetResidue>.txt
    """
    _label = 'Generate rays'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addHidden(GPU_LIST, params.StringParam,
                       label='Choose GPU ID',
                       default='1',
                       help="GPU may have several cores. Set it one if "
                            "you don't know what we are talking about but you have a GPU."
                            "For DARC, first core index is 1, second 2, and so on. Write 0 if you do not want"
                            "to use GPU")

        form.addParam("inputAtomStruct", params.PointerParam, pointerClass="AtomStruct",
                      label="Atomic structure",
                      important=True,
                      allowsNull=False,
                      help="Select the atomic structure of the prepared protein file")

        form.addParam("target_residue", params.IntParam,
                      label="Target residue",
                      important=True,
                      allowsNull=False,
                      help="Specify the number id of protein residue that you want to "
                           "use to generate the shape of pocket, around of that residue,"
                           " (preferably in the protein surface) using the rays method.")

        form.addParam("multiple_target", params.BooleanParam,
                      label='Multiple target residues',
                      default=False,
                      important=True,
                      help='Use to define the bound site for the small molecules')

        form.addParam('target_residues', params.StringParam, condition='multiple_target',
                      label='Additional target residues',
                      help='Write the additional number of protein residues used to throw the rays and create'
                           'the docking surface. The numbers have to be separated by semicolon, comma or space.\n\n'
                           'Examples: \n\n'
                           '1- 99,61 \n\n'
                           '2- 99;61 \n\n'
                           '3- 99 61  \n\n')

        form.addParam("change_origin", params.BooleanParam,
                      label='Change origin for ray generation',
                      default=False,
                      important=True,
                      help='If you want change the default () origin of rays generation, select "yes"')


        form.addParam("origin_residue", params.IntParam,
                      label="Residue origin",
                      important=True,
                      condition='change_origin',
                      allowsNull=True,
                      help="Specify the number id of protein residue that you want to"
                           " use as origin points for casting rays. This origin will be"
                           " the center of mass of that residue. ")

        form.addParam("electrostatics", params.BooleanParam,
                      label='Include electrostatics calculations',
                      default=False,
                      important=True,
                      help='')

        form.addParam("grid", params.PointerParam, pointerClass="GridADT",
                      condition="electrostatics",
                      label="Grid file",
                      important=True,
                      allowsNull=False,
                      help="Select the grid file")


        advanced = form.addGroup("Runs option",  expertLevel=LEVEL_ADVANCED)
        advanced.addParam("cseed", params.BooleanParam,
                       label='Use Constant Seed',
                       default=False,
                       help='Use this option to get reproducible results')

        advanced.addParam("seed", params.IntParam,
                       label='Set seed',
                       default=1111111,
                       condition="cseed",
                       help='Set a integer number as constant seed. The default one is 1111111 ')




    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('generate_ray')
        self._insertFunctionStep('createOutput')

    def convertInputStep(self):
        if self.electrostatics.get():
            adtGridName = self.grid.get().getFileName().split('/')[-1]
            self.agdGrid = adt2agdGrid(self.grid.get(), self._getExtraPath(adtGridName.replace('.e.map', '.agd')))

    def generate_ray(self):
        """Generate the txt and pdb file with the protein pocket mapping around a given residue
        """

        # Create the args of the program
        args = ""

        # Add protein file where the program will generate the rays (REQUIRED)
        pdb_file = self.inputAtomStruct.get().getFileName()
        pdb_file_extra = os.path.join(self._getExtraPath(), os.path.basename(pdb_file))
        createLink(pdb_file, pdb_file_extra)

        args += " -protein %s" % os.path.basename(pdb_file_extra)


        # Add Database path
        database_path = os.path.join(Plugin.getHome(), ROSETTA_DATABASE_PATH)
        args += " -database %s" % database_path


        # Add the specific residue that will be the center of ray generation (REQUIRED)
        # To use multiple origin points for casting rays or not
        if self.multiple_target.get():
            residues_string = self.target_residues.get()
            res = list(filter(None, re.split(",|;| ", residues_string.upper())))
            target_residues = ",".join(res)
            #res = list(map(int, res))
            args += " -central_relax_pdb_num %s,%s" %(self.target_residue.get(), target_residues)
            args += " -multiple_origin"
        else:
            args += " -central_relax_pdb_num %s" % self.target_residue.get()


        # By default the PocketGrid expands if pocket points are identified near the edge, this flag disables the autoexpanding feature.
        args += " -pocket_static_grid"


        # To use the mass center of a amino acid as origin of throwing tha rays
        if self.change_origin.get():

            args += " -set_origin "
            args += "-origin_res_num %i" % (self.origin_residue.get())  # number of residue in a chain

        # To include electrostatics calculations
        if self.electrostatics.get():
            # args += " -add_electrostatics"
            grid_file = self.agdGrid.getFileName()
            args += " -espGrid_file %s" % os.path.abspath(grid_file)

            # -add_electroestatics
            # -espGrid_file %.agd
        else:
            # Used to ignore electrostatic score and perform shape only calculation
            args += " -darc_shape_only"

        if self.cseed.get():
            args += " -run:constant_seed"
            args += " -run:jran %s" %self.seed.get()

        # Generate 2 file with different formats (pdb (rays are hetatm) and txt).
        # Run Make Ray Files w/wo GPU
        if GPU_LIST == 0:
            Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES), args, cwd=self._getExtraPath())
        else:
            args += " -gpu %s" % str(self.gpuList.get())
            Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES_GPU), args, cwd=self._getExtraPath())


    def createOutput(self):
        """ Create the greeting objects and give them the appropriate names
        """
        #Create filename

        pdb_ini = self.inputAtomStruct.get().getFileName()
        filename = os.path.splitext(os.path.basename(pdb_ini))[0]
        target_res = self.target_residue.get()


        if self.multiple_target.get():
            residues_string = self.target_residues.get()
            res = list(filter(None, re.split(",|;| ", residues_string.upper())))
            target_residues = ",".join(res)
            pdb_file_out_0 = self._getExtraPath('ray_%s_%s,%s.pdb' % (filename, target_res, target_residues))
            txt_file_out_0 = self._getExtraPath('ray_%s_%s,%s.txt' % (filename, target_res, target_residues))

            target_residues = "_".join(res)
            pdb_file_out = self._getExtraPath('ray_%s_%s_%s.pdb' % (filename, target_res, target_residues))
            txt_file_out = self._getExtraPath('ray_%s_%s_%s.txt' % (filename, target_res, target_residues))

            shutil.move(pdb_file_out_0, pdb_file_out)
            shutil.move(txt_file_out_0, txt_file_out)


        else:
            pdb_file_out = self._getExtraPath('ray_%s_%s.pdb' % (filename, target_res))
            txt_file_out = self._getExtraPath('ray_%s_%s.txt' % (filename, target_res))


        if os.path.exists(os.path.abspath(pdb_file_out)) and os.path.exists(os.path.abspath(txt_file_out)):
            target = RaysStruct(filename=pdb_file_out)
            rayfile = RaysProtein(filename=txt_file_out)
        else:
            raise Exception("Something wrong in output creation")


        if self.electrostatics.get():
            grid_file = self.agdGrid.getFileName()
            name = os.path.splitext(os.path.basename(grid_file))[0]
            agd_darcFile = self._getExtraPath("DARC_%s.agd" % name)
            if os.path.exists(os.path.abspath(agd_darcFile)):
                agd_darc = GridAGD(filename=agd_darcFile)
                outputdict = {'outputStructure': target, 'outputRay_TXT': rayfile, 'outputGRID_AGD': agd_darc}
        else:
            outputdict={'outputStructure': target, 'outputRay_TXT': rayfile}

        self._defineOutputs(**outputdict)
        self._defineSourceRelation(self.inputAtomStruct, target)




    # --------------------------- UTILS functions ------------------------------

    def length_residues_pdb(self, pdb):
        from Bio import PDB
        parser = PDB.PDBParser()

        name = os.path.splitext(os.path.basename(pdb))[0]
        structure = parser.get_structure(name, pdb)
        protein = structure[0]

        no_res = 0;  no_non_res = 0
        for chain in protein:
            for r in chain.get_residues():
                if r.id[0] == ' ':
                    no_res += 1
                else:
                    no_non_res += 1
        return no_res


    def _validate(self):
        """ Validate if the introduced residues is in the protein
        """
        errors = []

        no_res = self.length_residues_pdb(self.inputAtomStruct.get().getFileName())
        no_res_center = self.target_residue.get()
        no_res_origin = self.origin_residue.get()

        if (no_res_center > no_res) or (no_res_center < 1):
            errors.append('The given residue number is outside the expected '
                          'range of protein residues. The range of residues is'
                          ' *[1, %i]*.' % no_res)

        if self.change_origin.get():
            if (no_res_origin > no_res) or (no_res_origin < 1):
                errors.append('The given residue number is outside the expected '
                              'range of protein residues. The range of residues is'
                              ' *[1, %i]*.' % no_res)


        if self.multiple_target.get():
            residues_string = self.target_residues.get()
            res = list(filter(None, re.split(",|;| ", residues_string.upper())))

            for residue in res:
                if (int(residue) > no_res) or (int(residue) < 1):
                    errors.append('The given residue number is outside the expected '
                                  'range of protein residues. The range of residues is'
                                  ' *[1, %i]*.' % no_res)

        return errors