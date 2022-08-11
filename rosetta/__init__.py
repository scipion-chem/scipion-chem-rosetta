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

import pwem
import os
import fnmatch
import pyworkflow.utils as pwutils
from pwem import Config as emConfig

from pyworkflow.utils import Environ
from .constants import *



_version_ = "3.12"
_logo = "rosetta_icon.png"
_references = ['LeaverFay2011']


ROSETTA_DIC = {'name': 'rosetta', 'version': '3.12', 'home': 'ROSETTA_HOME'}

class Plugin(pwem.Plugin):
    _homeVar = ROSETTA_DIC['home']
    _pathVars = [ROSETTA_DIC['home']]
    _supportedVersions = [ROSETTA_DIC['version']]


    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. Set the Rosetta path on the computer
        """
        cls._defineVar(ROSETTA_DIC['home'], cls.getRosettaDir())


    @classmethod
    def defineBinaries(cls, env):
        pass


    @classmethod  #  Test that
    def getEnviron(cls):
        """ Setup the environment variables needed to launch programs from Rosetta """
        environ = pwutils.Environ(os.environ)
        environ.update({
            'ROSETTA': cls.getHome(''),
            'PATH': cls.getHome(''),
        }, position=Environ.BEGIN)
        return environ


    @classmethod
    def getProgram(cls, progName, path=ROSETTA_BINARIES_PATH):
        """ Return the program binary that will be used. """
        return os.path.join(cls.getHome(),
                            path,
                            progName)

    @classmethod
    def runRosettaProgram(cls, program, args=None, extraEnvDict=None, cwd=None):
        """ Internal shortcut function to launch a Rosetta program. """
        env = cls.getEnviron()
        if extraEnvDict is not None:
            env.update(extraEnvDict)
        pwutils.runJob(None, program, args, env=env, cwd=cwd)

    @classmethod
    def validateInstallation(cls):
        """ Check if the installation of this protocol is correct.
             Returning an empty list means that the installation is correct
             and there are not errors. If some errors are found, a list with
             the error messages will be returned.
             """
        missingPaths = []
        if not os.path.exists(os.path.expanduser(cls.getVar(ROSETTA_DIC['home']))):
            missingPaths.append("Path of Rosetta does not exist (%s) : %s " % (ROSETTA_DIC['home'],
                                                                               cls.getVar(ROSETTA_DIC['home'])))
        return missingPaths



    # ---------------------------------- Utils functions  -----------------------
    @staticmethod
    def find(path, pattern):  # This function is analogous to linux find (case of folders)
        paths = []
        for root, dirs, files in os.walk(path):
            for name in dirs:
                if fnmatch.fnmatch(name, pattern):
                    paths.append(os.path.join(root, name))
        return paths


    @classmethod
    def getRosettaDir(cls, fn=""):
        fileList = cls.find(emConfig.EM_ROOT, "rosetta_bin_linux*")
        if len(fileList) == 0:
            return None
        else:
            if fn == "":
                return fileList[0]
            else:
                return os.path.join(fileList[0], fn)
