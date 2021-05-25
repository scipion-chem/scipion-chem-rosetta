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


import os

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, Viewer
from pwem.viewers.viewer_chimera import (Chimera)

from rosetta.protocols.protocol_generate_rays import Rosetta_make_rayFile



class RaysViewer(Viewer):
    """ Viewer for Rosetta program make ray file
    """
    _label = 'Rays Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [Rosetta_make_rayFile]


    def _visualize(self, obj, **args):
        # Get the rays pdb file
        rays = (self.protocol.outputStructure).getFileName()

        # Get the prepared protein pdb file
        base = os.path.basename(rays)
        name = (base).split("_")[1]
        print(name)
        pdb = self.protocol._getExtraPath("%s.pdb" % name)

        dim = 150.
        sampling = 1.
        bildFileName = self.protocol._getExtraPath("axis_output.bild")
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=bildFileName,
                                         sampling=sampling)


        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")
        f = open(fnCmd, 'w')
        # change to workingDir
        # If we do not use cd and the project name has an space the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())
        f.write("open %s\n" % bildFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates
        f.write("open %s\n" % rays)
        f.write("style stick\n")
        f.write("open %s\n" % pdb)
        f.write("show surfaces\n")

        # run in the background
        Chimera.runProgram(Chimera.getProgram(), os.path.abspath(fnCmd) + " &")
        return []