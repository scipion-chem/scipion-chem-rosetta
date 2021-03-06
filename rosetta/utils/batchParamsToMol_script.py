#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

'''
A wrapper around molfile_to_params.py for generating large numbers of params files and giving them unique names that don't conflict with anything in the standard Rosetta database.  The directory that this script produces can be directly parsed with rosetta using the -in:file:extra_res_batch path
Also: Various functions for manipulating params files.
Original code from RosettaCommons has just been adapted to Python3 and to solve a comflict in
get_disallowed_ligands function as stated in https://www.rosettacommons.org/node/11278

Author: Sam DeLuca
'''
import itertools
import os
import subprocess
import shutil
import fnmatch
from optparse import OptionParser
mol_to_params = "~/rosetta/rosetta_source/src/python/apps/public/molfile_to_params.py"

char_set = ['0','1','2','3','4','5',
            '6','7','8','9','A','B',
            'C','D','E','F','G','H',
            'I','J','K','L','M','N',
            'O','P','Q','R','S','T',
            'U','V','W','X','Y','Z']


def make_params(mol_path,ligand_name,base_name):
    '''Make the params file, append the conformers, and move the output into the proper directory'''
    logfile = "params/"+base_name+"/log.txt"
    log = open(logfile,'a')
    params_cmd = mol_to_params+" -n " +ligand_name + " " + mol_path
    params_child = subprocess.call(params_cmd,shell=True,stdout=log,stderr=log)
    log.close()
    if params_child != 0:
        return params_child

    pdb_path = "params/"+base_name+"/"+ligand_name+"_conformers.pdb"
    pdb_file = open(pdb_path, 'w')
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file,ligand_name+"_*.pdb"):
            conformer =open(file,'r')
            pdb_file.writelines(conformer.readlines())
            conformer.close()
            os.remove(file)

    pdb_file.close()

    if os.path.exists(ligand_name+".params"):
        paramsfile = open(ligand_name+".params",'a')
    else:
        return 1
    paramsfile.write("PDB_ROTAMERS "+ ligand_name+"_conformers.pdb\n")
    paramsfile.close()
    shutil.move(ligand_name+".params", "params/"+base_name+"/"+ligand_name+".params")
    return params_child

def get_name_from_params(path,database):
    '''Given the path to a params file, return the IO_STRING value.  if it doesn't exist, return None'''
    path=path.strip()
    #if path is commented, ignore it
    if len(path) == 0:
        return 0
    if path[0] == "#":
        return 0
    fullpath = database+path
    #print fullpath
    path = path.split("/")
    #is path a params file?
    type = path[-1].split(".").pop()
    if type == "params":
        if os.path.exists(fullpath):
            paramsfile = open(fullpath,'r')
            for line in paramsfile:
                line = line.split()
                if len(line) >0:
                    if line[0] == "IO_STRING":
                        return line[1]
        else:
            return None
    else:
        return None

def get_disallowed_ligands(database):
    '''Return a set of 3 letter names which are already assigned to ligands in the database'''
    residue_type_set_path = database+"chemical/residue_type_sets/"
    residue_set_list = os.listdir(residue_type_set_path)
    disallowed_ligands = set()
    for residue_set in residue_set_list:
        if residue_set[0] != "." and not residue_set == 'how_to_add_new_residue_types.txt':
            current_residue_set_path = residue_type_set_path+residue_set+"/"
            residue_types = open(current_residue_set_path+"residue_types.txt",'r')
            for line in residue_types:
                name = get_name_from_params(line, current_residue_set_path)
                #print name
                if name != None:
                    disallowed_ligands.add(name)
    return disallowed_ligands

def rename_param_file(param_path,new_name,new_conformer_path):
    '''Rename a param file residue and update the conformer file path'''
    param_file = open(param_path,'r')
    param_lines = [x.rstrip() for x in param_file]
    param_file.close()

    for index,line in enumerate(param_lines):
        fields = line.split()
        if fields[0] == "NAME":
            param_lines[index] = "NAME "+new_name
            continue
        if fields[0] == "IO_STRING":
            param_lines[index] = "IO_STRING "+new_name+" "+fields[2]
            continue
        if fields[0] == "PDB_ROTAMERS":
            param_lines[index] = "PDB_ROTAMERS "+new_conformer_path

    param_file = open(param_path,'w')
    for line in param_lines:
        param_file.write(line+"\n")


def rename_pdb_file(pdb_path,new_name):
    '''Renames all the HETATM resnames in the specified pdb to new_name'''
    pdb_file = open(pdb_path,'r')
    pdb_lines = [x.rstrip() for x in pdb_file]
    pdb_file.close()
    assert(len(new_name) == 3)
    for index,line in enumerate(pdb_lines):
        if len(line) < 6:
            continue
        if line[0:6] == "HETATM":
            pdb_lines[index] = line[:17]+new_name+line[20:]

    pdb_file = open(pdb_path,'w')
    for line in pdb_lines:
        pdb_file.write(line+"\n")

def getBatchMolToParamsPath():
    return __file__


if __name__ == "__main__":
    usage = "%prog -d path/to/rosetta/database --script_path=/path/to/molfile_to_params.py molfile_list.txt"
    parser=OptionParser(usage)
    parser.add_option("-d", dest ="database",help="path to minirosetta database",default=1)
    parser.add_option("--script_path",dest="script",help="location of the molfile_to_params script",default = "")
    parser.add_option("--exclusion_list",dest="excluded",help="list of ligand names to manually exclude",default="")
    (options, args) = parser.parse_args()

    if not os.path.exists("params"):
        os.mkdir("params")

    if options.script != "":
        mol_to_params = options.script


    if not os.path.exists(mol_to_params):
        parser.error("ERROR: make sure to specify the path to the molfile_to_params script with the option --params")

    if len(args) != 1:
        parser.error("ERROR: you must specify the path to a file containing a list of molfile paths")

    molfile_list_path = args[0]
    ligand_names = itertools.product(char_set,repeat=3)

    disallowed_ligands = get_disallowed_ligands(options.database+"/")

    if(options.excluded != ""):
        exclude_list = open(options.excluded,'r')
        for line in exclude_list:
            line = line.rstrip()
            disallowed_ligands.add(line)
        exclude_list.close()


    molfile_list = open(molfile_list_path, 'r')
    for molfile in molfile_list:
        molfile = molfile.strip()
        for ligand_name in ligand_names:
            ligand_name = "".join(ligand_name)
            if ligand_name not in disallowed_ligands:
                break
        mol_base_name = molfile.split("/").pop().split(".")[0]

        if os.path.exists("params/"+ mol_base_name +"/"+ ligand_name +".params"):
            print("ligand " + mol_base_name + " already processed, continuing")
            continue
        if not os.path.exists("params/" + mol_base_name):
            os.mkdir("params/" + mol_base_name)
        print(ligand_name, mol_base_name)
        params_status = make_params(molfile, ligand_name, mol_base_name)
        if params_status != 0:
            print("WARNING: failed to generate params for " + mol_base_name)

