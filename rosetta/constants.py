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

generateStructuresXML = '''<ROSETTASCRIPTS>
	<SCOREFXNS>
		<ScoreFunction name="cen" weights="score4_smooth_cart">
			<Reweight scoretype="elec_dens_fast" weight="15"/>
		</ScoreFunction>

		<ScoreFunction name="dens_soft" weights="soft_rep">
			<Reweight scoretype="cart_bonded" weight="0.5"/>
			<Reweight scoretype="pro_close" weight="0.0"/>
			<Reweight scoretype="fa_sol" weight="0.0"/> #REMOVE THIS LINE IF NOT A MEMBRANE PROTEIN
			<Reweight scoretype="elec_dens_fast" weight="25"/>
		</ScoreFunction>
		
		<ScoreFunction name="dens" weights="talaris2014_cart">
			<Reweight scoretype="elec_dens_fast" weight="25"/>
			<Reweight scoretype="fa_sol" weight="0.0"/> #REMOVE THIS LINE IF NOT A MEMBRANE PROTEIN
		<Set scale_sc_dens_byres="R:0.76,K:0.76,E:0.76,D:0.76,M:0.76, C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88, A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88"/> #These values were empirically determined by the ROSETTA group and SHOULD NOT BE CHANGED
		</ScoreFunction>
	</SCOREFXNS>
	<MOVERS>
		<SetupForDensityScoring name="setupdens"/>
		<LoadDensityMap name="loaddens" mapfile="../noMask0723.mrc"/>
		<SwitchResidueTypeSetMover name="tocen" set="centroid"/>
		<MinMover name="cenmin" scorefxn="cen" type="lbfgs_armijo_nonmonotone" max_iter="200" tolerance="0.00001" bb="1" chi="1" jump="ALL"/>
		#This section will use Z-scores to determing poorly fit areas for rebuilding with increasingly more strict cutoffs
		<CartesianSampler name="cen5_50" automode_scorecut="-0.5" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="2.0" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<CartesianSampler name="cen5_60" automode_scorecut="-0.3" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="1.5" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<CartesianSampler name="cen5_70" automode_scorecut="-0.1" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="1.5" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<CartesianSampler name="cen5_80" automode_scorecut="0.0" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="1.0" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<FastRelax name="relaxcart" scorefxn="dens" repeats="1" cartesian="1"/>
	</MOVERS>
	<PROTOCOLS>
		<Add mover="setupdens"/>
		<Add mover="loaddens"/>
		<Add mover="cenmin"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_50"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_60"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_70"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_80"/>
		<Add mover="relaxcart"/>
		<Add mover="relaxcart"/>
	</PROTOCOLS>
	<OUTPUT scorefxn="dens"/>
</ROSETTASCRIPTS>'''
