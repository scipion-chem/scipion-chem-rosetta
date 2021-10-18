# -*- coding: utf-8 -*-
# **************************************************************************
# Module to declare protocols
# Find documentation here: https://scipion-em.github.io/docs/docs/developer/creating-a-protocol
# **************************************************************************

#from .protocol_protein_preparation import Rosetta_protein_preparation

from .protocol_target_preparation import RosettaProteinPreparation
from .protocol_generate_rays import Rosetta_make_rayFile
from .protocol_darc import RosettaProtDARC

from .protocol_converter import ConvertStructures