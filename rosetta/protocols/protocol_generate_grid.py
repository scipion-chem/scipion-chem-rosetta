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
This protocol uses a Rosetta suite program (named score) to build the missing atoms
in the protein, including the hydrogens. This is necessary to subsequently
generate the electrostatic potential grid.
"""

#IMPORT
import os
import math
import pandas as pd

from pyworkflow.protocol.params import  FloatParam, PointerParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED


from pwem.protocols import EMProtocol
from pwem.convert import AtomicStructHandler


from bioinformatics import Plugin as bioinformatics_plugin

from rosetta.objects import GridAGD


class Autodock_GridGeneration(EMProtocol):
    """
    Calls ADT to prepare a grid with an input that is the prepared target protein
    """

    _label = 'Grid generation with ADT'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputpdb', PointerParam, pointerClass="AtomStruct",
                      label='Input structure:',
                      allowsNull=False,
                      help="Select the prepared atomic structure of "
                           "the docking target protein")

        form.addParam('radius', FloatParam, label='Radius')

        form.addParam('spacing', FloatParam, default=0.375, label='Step size (A)',
                      expertLevel=LEVEL_ADVANCED,
                      help="Distance between each point in the electrostatic grid."
                           " This value is used to adjust the radius as number of "
                           "(x,y,z) points : radius/spacing = number of points along"
                           " 3 dimensions ")



    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('getpdbqt')
        self._insertFunctionStep('prepareGrid')
        self._insertFunctionStep('createOutput')



    def getpdbqt(self):
        """ Prepare the pdbqt used to generate the grid file e-map"""

        filename = self.inputpdb.get().getFileName()
        if filename.endswith('.cif'):
            fnIn = self._getTmpPath("atomStruct.pdb")
            name_protein = (os.path.basename(filename)).split(".")[0]
            aStruct1 = AtomicStructHandler(filename)
            aStruct1.write(fnIn)
        else:
            fnIn = filename
            name_protein = (os.path.basename(fnIn)).split(".")[0]

        fnOut = self._getExtraPath('%s.pdbqt' %name_protein)

        args = ' -v -r %s -o %s' % (fnIn, fnOut)

        program = "prepare_receptor4"
        self.runJob(bioinformatics_plugin.getMGLPath('bin/pythonsh'),
                    bioinformatics_plugin.getADTPath('Utilities24/%s.py' % program) + args)



    def prepareGrid(self):
        """
        """
        # Name of the grid
        inputpdb = self.inputpdb.get().getFileName()
        filename = os.path.abspath(inputpdb)
        name_protein = (os.path.basename(inputpdb)).split(".")[0]

        # pdbqt gasteiger path
        pdbqt = self._getExtraPath('%s.pdbqt' % name_protein)

        # Center of the molecule
        structure, x_center, y_center, z_center = self.calculate_centerMass(filename)

        # Create the GPF file, required by autogrid to build the grid and glg
        npts = (self.radius.get()+10)/self.spacing.get()  # x,y,z points of the grid

        gpf_file = self.generate_gpf(filename=pdbqt, spacing=self.spacing.get(),
                                     xc=x_center, yc=y_center, zc=z_center,
                                     npts=npts)

        glg_file = os.path.abspath(self._getExtraPath("library.glg"))
        open(glg_file, mode='a').close()

        args = "-p %s -l %s" % (gpf_file, glg_file)
        self.runJob(bioinformatics_plugin.getAutodockPath("autogrid4"), args, cwd=self._getExtraPath())

        # Build AGD format
        e_map = self._getExtraPath("%s.e.map" %name_protein)

        agdfile = self._getPath("%s.agd" %name_protein)


        with open(agdfile, "w") as agd:
            agd.write("Title: \n")
            agd.write("Mid:   %s     %s    %s\n" % (round(x_center,3), round(y_center,3), round(z_center,3)))
            agd.write("Dim:     %s     %s     %s\n" % (round(npts, None), round(npts, None), round(npts, None)))
            agd.write("Spacing:     %s\n" % (self.spacing.get()))
            agd.write("Values:\n")

            with open(e_map, "r") as emap:
                for line in emap.readlines():
                    if not line.startswith(("GRID", "MACROMOLECULE", "SPACING", "CENTER", "NELEMENTS")):
                        agd.write(line)


    def createOutput(self):
        fnIn = self.inputpdb.get().getFileName()
        name_protein = (os.path.basename(fnIn)).split(".")[0]
        fnOut = self._getPath("%s.agd" %name_protein)

        if os.path.exists(fnOut):
            grid = GridAGD(filename=fnOut)
            self._defineOutputs(outputGrid=grid)
            self._defineSourceRelation(self.inputpdb, grid)


    # --------------------------- Utils functions --------------------

    def calculate_centerMass(self, filename):
        """
        Returns the geometric center of mass of an Entity (anything with a get_atoms function in biopython).
        Geometric assumes all masses are equal
        """

        try:
            structureHandler = AtomicStructHandler()
            structureHandler.read(filename)
            center_coord = structureHandler.centerOfMass(geometric=True)
            structure = structureHandler.getStructure()

            return structure, center_coord[0], center_coord[1], center_coord[2]  # structure, x,y,z

        except Exception as e:
            print("ERROR: ", "A pdb file was not entered in the Atomic structure field. Please enter it.", e)
            return



    def calculate_OuterDistance(self, expand=10):
        """
        Returns the higher distance between the mass center and the different atoms of the protein
        """

        # Structure and this Center of mass.
        structure, x_center, y_center, z_center = self.calculate_centerMass()

        # First point to calculate euclidean distance
        dist = 0

        # iterate for all aminoacids
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        coord = atom.get_coord()
                        x = coord[0]
                        y = coord[1]
                        z = coord[2]

                        dist_new = (math.sqrt((x-0)**2 + (y-0)**2 + (z-0)**2)) /self.spacing.get()
                        #dist_new = (math.sqrt((x - x_center) ** 2 + (y - y_center) ** 2 + (z - z_center) ** 2)) * (1 / self.spacing.get())

                        if dist_new > dist:
                            dist=dist_new

        return (dist + expand)
        #dist_new = (math.sqrt((xmin - x_center) ** 2 + (ymin - y_center) ** 2 + (zmin - z_center) ** 2)) * (1 / self.spacing.get())





    def generate_gpf(self, filename, spacing, xc, yc, zc, npts):
        """

        :return:
        """

        filename = filename
        name_protein = (os.path.basename(filename)).split(".")[0]
        gpf_file = self._getExtraPath("%s.gpf" % name_protein)
        npts = round(npts, None)

        list_receptor = []
        index = 0
        # Read, parse and take receptor types atoms
        df = pd.read_csv(os.path.abspath(filename), sep='\s{1,}', header=None, engine='python')
        df.fillna(0)

        df_type_atoms = df[12].fillna("NA").tolist()
        for type in df_type_atoms:
            if not type in list_receptor:
                if not df[0][index] == "TER":
                    list_receptor.append(type)
            index += 1
        list_receptor = " ".join(list_receptor)



        with open(os.path.abspath(gpf_file), "w") as file:
            file.write("npts %s %s %s                        # num.grid points in xyz\n" % (npts, npts, npts))
            file.write("gridfld %s.maps.fld                # grid_data_file\n" % (name_protein))
            file.write("spacing %s                          # spacing(A)\n" % (spacing))
            file.write("receptor_types %s     # receptor atom types\n" % (list_receptor))
            file.write("ligand_types A C HD N NA OA SA       # ligand atom types\n")
            file.write("receptor %s                  # macromolecule\n" % (os.path.abspath(filename)))
            file.write("gridcenter %s %s %s           # xyz-coordinates or auto\n" % (xc, yc, zc))
            file.write("smooth 0.5                           # store minimum energy w/in rad(A)\n")
            file.write("map %s.A.map                       # atom-specific affinity map\n" % (name_protein))
            file.write("map %s.C.map                       # atom-specific affinity map\n" % (name_protein))
            file.write("map %s.HD.map                      # atom-specific affinity map\n" % (name_protein))
            file.write("map %s.N.map                       # atom-specific affinity map\n" % (name_protein))
            file.write("map %s.NA.map                      # atom-specific affinity map\n" % (name_protein))
            file.write("map %s.OA.map                      # atom-specific affinity map\n" % (name_protein))
            file.write("map %s.SA.map                      # atom-specific affinity map\n" % (name_protein))
            file.write("elecmap %s.e.map                   # electrostatic potential map\n" % (name_protein))
            file.write("dsolvmap %s.d.map                  # desolvation potential map\n" % (name_protein))
            file.write("dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant" )


        return os.path.abspath(gpf_file)