# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
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

from pyworkflow.utils import Environ
from .constants import *



_version_ = "3.12" # Hacer como en phenix que toma cualquier versión, aunque tiene una por defecto. Recomendación COSS
_logo = "rosetta_icon.png"
_references = ['LeaverFay2011']



class Plugin(pwem.Plugin):
    _homeVar = ROSETTA_HOME
    _pathVars = [ROSETTA_HOME]
    _supportedVersions = [V3_12]


    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. Set the Rosetta path on the computer
        """
        cls._defineVar(ROSETTA_HOME, cls.getRosettaDir())
        cls._defineVar("OPENBABEL_ENV_ACTIVATION", 'conda activate openbabel-env')


    @classmethod
    def defineBinaries(cls, env):
        pass


    @classmethod  # Not used and can be deleted. Test
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
        return os.path.join(Plugin.getHome(),
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
    def getOpenbabelEnvActivation(cls):
        activation = cls.getVar("OPENBABEL_ENV_ACTIVATION")
        return activation

    @classmethod
    def addopenbabelPackage(cls, env, default=False):
        OPENBABEL_INSTALLED = 'openbabel_installed'

        # try to get CONDA activation command
        installationCmd = cls.getCondaActivationCmd()

        # Create the environment
        installationCmd += ' conda install openbabel -c conda-forge &&'

        # Flag installation finished
        installationCmd += ' touch %s' % OPENBABEL_INSTALLED

        openbabel_commands = [(installationCmd, OPENBABEL_INSTALLED)]

        envPath = os.environ.get('PATH', "")  # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage('openbabel',
                       tar='void.tgz',
                       commands=openbabel_commands,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def runOPENBABEL(cls, protocol, program="obabel", args=None, cwd=None):
        """ Run openbabel command from a given protocol. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getOpenbabelEnvActivation(), program)
        fullProgram = "obabel"
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd, numberOfThreads= 1)

    @classmethod
    def validateInstallation(cls):
        """ Check if the installation of this protocol is correct.
             Returning an empty list means that the installation is correct
             and there are not errors. If some errors are found, a list with
             the error messages will be returned.
             """
        missingPaths = []
        if not os.path.exists(os.path.expanduser(Plugin.getVar(ROSETTA_HOME))):
            missingPaths.append("Path of Rosetta does not exist (%s) : %s " % (ROSETTA_HOME, Plugin.getVar(ROSETTA_HOME)))
        return missingPaths



    # ---------------------------------- Utils functions to utils.py -----------------------
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
        fileList = Plugin.find("/home", "rosetta_bin_linux*")
        if len(fileList) == 0:
            return None
        else:
            if fn == "":
                return fileList[0]
            else:
                return os.path.join(fileList[0], fn)
