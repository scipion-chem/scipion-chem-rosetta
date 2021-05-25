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

import sys
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints


def gen_conformers(molecule, numConfs=200, maxAttempts=1000, pruneRmsThresh=0.5, useExpTorsionAnglePrefs=True,
                   useBasicKnowledge=True, enforceChirality=True, useRandomCoords=False,
                   numThreads=0, randomSeed=-1):
    ids = AllChem.EmbedMultipleConfs(molecule,
                                     numConfs=numConfs, # 200
                                     maxAttempts=maxAttempts, # 1000
                                     pruneRmsThresh=pruneRmsThresh, # 0.5
                                     useExpTorsionAnglePrefs=useExpTorsionAnglePrefs, # ET de ETKDG
                                     useBasicKnowledge=useBasicKnowledge, # K de ETKDG
                                     enforceChirality=enforceChirality, # True
                                     useRandomCoords=useRandomCoords, # False
                                     numThreads=numThreads, # 0 All available
                                     randomSeed=randomSeed, # -1 random
                                     )
    return list(ids)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print ("Usage: rdkit_conformer_generation.py <sdf input> <num conformers> <max attempts> "
               "<RMSD threshold>")
        exit()
    input_file = sys.argv[1]
    numConfs = int(sys.argv[2])
    maxAttempts = int(sys.argv[3])
    pruneRmsThresh = float(sys.argv[4])

    #https: // gist.github.com / tdudgeon / b061dc67f9d879905b50118408c30aac
