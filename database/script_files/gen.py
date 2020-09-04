#!/usr/bin/env python3

import os, sys
import numpy as np
import math
from distutils.dir_util import copy_tree
import multiprocessing as mp
import distutils.spawn
from shutil import copyfile
import glob
import re
import copy
import ntpath
import g_var

def forcefield_selection():
    ##### forcefield selection
    if not g_var.info:
        if g_var.ff != None:
            if os.path.exists(g_var.ff):
                g_var.forcefield_location, g_var.forcefield = path_leaf(g_var.ff)
                folder_copy_and_check(g_var.ff, g_var.final_dir+forcefield)
                print('\nYou have chosen to use your own forcefield: '+g_var.forcefield+' in '+g_var.forcefield_location)
            elif path_leaf(g_var.ff)[1]+'.ff' in g_var.forcefield_available:
                forcefield_number = g_var.forcefield_available.index(path_leaf(g_var.ff)[1]+'.ff')
            elif path_leaf(g_var.ff)[1] in g_var.forcefield_available:
                forcefield_number = g_var.forcefield_available.index(path_leaf(g_var.ff)[1])
            else:
                print('Cannot find forcefield: '+g_var.ff+'  please select one from below\n')   
        if 'forcefield_number' not in locals() and 'forcefield' not in locals():
            forcefield_number = database_selection(g_var.forcefield_available, 'forcefields')    
        if 'forcefield' not in locals():
            print('\nYou have selected the forcefield: '+g_var.forcefield_available[forcefield_number].split('.')[0])
            folder_copy_and_check(g_var.database_dir+'/forcefields/'+g_var.forcefield_available[forcefield_number], g_var.final_dir+g_var.forcefield_available[forcefield_number])
            g_var.forcefield_location, g_var.forcefield=g_var.database_dir+'forcefields/', g_var.forcefield_available[forcefield_number]
            g_var.opt['ff'] = g_var.forcefield_available[forcefield_number]

def fragment_selection():
    ##### fragment selection
    frag_location, fragment_number, fragments_available_other = [],[],[]
    if g_var.fg != None:
        for frag_val, frag_path in enumerate(g_var.fg):
            if os.path.exists(frag_path):
                frag_loc, fragments = path_leaf(frag_path)
                frag_location.append(os.path.abspath(frag_loc)+'/')
                fragment_number.append(frag_val)
                fragments_available_other.append(fragments)
    if len(fragment_number) == 0:
        fragment_number = fetch_frag_number(g_var.fragments_available)
        frag_location = [g_var.database_dir+'fragments/']*len(g_var.fragments_available)
        if g_var.fg is None:
            g_var.opt['fg'] = ''
            for database in fragment_number:
                g_var.opt['fg'] += g_var.fragments_available[database]+' '
    else:
        g_var.fragments_available = fragments_available_other
    fetch_residues(frag_location, g_var.fragments_available, fragment_number)

def correct_number_cpus():
    if g_var.ncpus != None:
        if g_var.ncpus > mp.cpu_count():
            print('you have selected to use more CPU cores than are available: '+str(g_var.ncpus))
            print('defaulting to the maximum number of cores: '+str(mp.cpu_count()))
            g_var.ncpus = mp.cpu_count()
    else:
        if mp.cpu_count() >= 8:
            g_var.ncpus = 8
        else:
            g_var.ncpus = mp.cpu_count()
    g_var.opt['ncpus'] = g_var.ncpus

def path_leaf(path):
    head, tail = ntpath.split(path)
    if not tail:
        return path.replace(ntpath.basename(head)+'/', ''), ntpath.basename(head)
    else:
        return path.replace(tail, ''), tail

### finds gromacs installation
def find_gromacs():
    if g_var.gmx != None:
        g_var.gmx=distutils.spawn.find_executable(g_var.gmx)
    else:
        g_var.gmx=distutils.spawn.find_executable('gmx')
    if g_var.gmx is None or type(g_var.gmx) != str:
        if os.environ.get("GMXBIN") != None:
            for root, dirs, files in os.walk(os.environ.get("GMXBIN")):
                for file_name in files:
                    if file_name.startswith('gmx') and file_name.islower() and '.' not in file_name:
                        g_var.gmx=distutils.spawn.find_executable(file_name)
                        if type(g_var.gmx) == str and g_var.gmx != None :
                            break
                        else:
                            g_var.gmx=None
                break
        if g_var.gmx is None:
            sys.exit('Cannot find gromacs installation')
    g_var.opt['gromacs'] = g_var.gmx

def trunc_coord(xyz):
    xyz_new = []
    for coord in xyz:
        if len(str(coord)) > 8:
            if '.' in str(coord):
                xyz_new.append(np.round(coord, 7-len(str(int(coord)))))
            else:
                xyz_new.append(np.round(coord, 8-len(str(int(coord)))))
        else:
            xyz_new.append(coord)
    return xyz_new[0],xyz_new[1],xyz_new[2]

def calculate_distance(p1, p2):
    return np.sqrt(((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2)+((p1[2]-p2[2])**2))

def file_copy_and_check(file_in,file_out):
    if not os.path.exists(file_out) and os.path.exists(file_in):
        copyfile(file_in, file_out)

def folder_copy_and_check(folder_in,folder_out):
    if not os.path.exists(folder_out):
        copy_tree(folder_in, folder_out)

def flags_used():
    print('\nAll variables supplied have been saved in : \n'+g_var.input_directory+'script_inputs.dat')
    os.chdir(g_var.input_directory)
    with open('script_inputs.dat', 'w') as scr_input:
        scr_input.write('\n'+g_var.opt['input']+'\n')
        for var in g_var.opt:
            if var != 'input':
                scr_input.write('{0:9}{1:15}\n'.format(var,str(g_var.opt[var])))

def is_hydrogen(atom):
    if str.isdigit(atom[0]) and atom[1] != 'H':
        return False
    elif not str.isdigit(atom[0]) and not atom.startswith('H'):
        return False
    else:
        return True

def fetch_chain_groups():
    if g_var.group != None:
        if g_var.group[0] not in ['all','chain']:
            for group_val, group in enumerate(g_var.group):
                for chain in group.split(','):
                    g_var.group_chains[int(chain)]=group_val 
        else:
            g_var.group_chains =  g_var.group[0]   

def split_swap(swap):
    try:
        res_range = re.split(':', swap)[2].split(',')
        res_id = []
        for resid_section in res_range:
            if '-' in resid_section:
                spt = resid_section.split('-')
                for res in range(int(spt[0]), int(spt[1])+1):
                    res_id.append(res)
            else:
                res_id.append(int(resid_section))
        return res_range, res_id
    except BaseException:
        return 'ALL', 'ALL'

def sort_swap_group():
    if g_var.swap != None:
        for swap in g_var.swap:
            res_s = re.split(':', swap)[0].split(',')
            if re.split(':', swap)[1].split(',') is not type(int):
                res_e = re.split(':', swap)[1].split(',')
            else:
                sys.exit('swap layout is not correct')
            
            if len(res_s) == len(res_e):
                res_range, res_id = split_swap(swap)
                if res_s[0] not in g_var.swap_dict:
                    g_var.swap_dict[res_s[0]]={}
                if len(res_s) == 1:
                    g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]={'ALL':'ALL'}
                else:
                    g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]={}
                    for bead in range(1,len(res_s)):
                        g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]][res_s[bead]]=res_e[bead]
                g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]['resid']=res_id
                g_var.swap_dict[res_s[0]][res_s[0]+':'+res_e[0]]['range']=res_range
            else:
                sys.exit('The length of your swap groups do not match')
        print_swap_residues()


def create_ion_list(ion_pdb):
    with open(ion_pdb, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('['):
                header = strip_header(line)
                if header not in g_var.ions:
                    g_var.ions.append(header)


def print_swap_residues():
    print('\nYou have chosen to swap the following residues\n')
    print('{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}'.format('residue', 'bead', '     ', 'residue', 'bead', 'range'))
    print('{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}'.format('-------', '----', '     ', '-------', '----', '-----'))
    for residue in g_var.swap_dict:
        for swap in g_var.swap_dict[residue]:
            bead_s, bead_e='', ''
            for bead in g_var.swap_dict[residue][swap]:
                if bead not in ['resid', 'range']:
                    bead_s+=bead+' '
                    bead_e+=g_var.swap_dict[residue][swap][bead]+' '
                elif bead == 'range':
                    if g_var.swap_dict[residue][swap]['range'] != 'ALL':
                        ran=''
                        for resid_section in g_var.swap_dict[residue][swap]['range']:
                            ran += resid_section+', '
                        ran = ran[:-2]
                    else:
                        ran = g_var.swap_dict[residue][swap]['range']
            print('{0:^10}{1:^5}{2:^11}{3:^11}{4:^11}{5:^11}'.format(swap.split(':')[0], bead_s, ' --> ', swap.split(':')[1], bead_e, ran))
    

def new_box_vec(box_vec, box):
    print('box cutting only works for cubic boxes currently')
    box_vec_split = box_vec.split()[1:]
    box_vec_values, box_shift = [], []
    for xyz_val, xyz in enumerate(box):
        if xyz != 0:
            box_shift.append((float(box_vec_split[xyz_val])/2) - (float(xyz)/2))
            box_vec_values.append(np.round(float(xyz), 3))
        else:
            box_shift.append(0)
            box_vec_values.append(float(box_vec_split[xyz_val]))
    box_vec = g_var.box_line%(box_vec_values[0], box_vec_values[1], box_vec_values[2], box_vec_split[3],box_vec_split[4],box_vec_split[5])
    return box_vec, np.array(box_shift)

def strip_header(line):
    line = line.replace('[','')
    line = line.replace(']','')
    if len(line.split())>1:
        sys.exit('There is a issue in one of the fragment headers: \n', line)
    return line.strip()

def sep_fragments_topology(location):
    topology={}
    group = 1
    topology = copy.deepcopy(g_var.topology)
    if os.path.exists(location+'.top'):
        with open(location+'.top', 'r') as top_input:
            for line_nr, line in enumerate(top_input.readlines()):
                if len(line) > 0:
                    if not line.startswith('#'):
                        if line.startswith('['):
                            top = strip_header(line).upper()
                            if top in topology:
                                topology_function = top
                            else:
                                print('The topology header line is incorrect, therefore ignoring: \n'+location+'.top')
                                topology_function = ''                              
                        else:
                            line_sep = line.split()
                            if len(line_sep) > 0:
                                if topology_function == 'GROUPS':
                                    if not g_var.mod:
                                        for bead in line_sep:
                                            topology[topology_function][bead] = group
                                        group += 1
                                    topology[topology_function]['group_max'] = group
                                elif topology_function in ['STEER']:
                                    topology[topology_function] += line_sep
                                elif topology_function in ['N_TERMINAL', 'C_TERMINAL']:
                                    topology[topology_function] = ''.join(line_sep)
                                elif topology_function == 'CHIRAL':
                                    if len(line_sep) == 5:
                                        topology[topology_function]['atoms']+=line_sep
                                        topology[topology_function][line_sep[0]]={'m':line_sep[1], 'c1':line_sep[2], 'c2':line_sep[3], 'c3':line_sep[4]}
                                    else:
                                        print('The chiral group section is incorrect: \n'+location+'.top')
                                elif topology_function == 'CONNECT':
                                    if len(line_sep) == 4:
                                        topology[topology_function]['atoms'][line_sep[1]] = int(line_sep[3])
                                        if line_sep[0] in topology[topology_function]:
                                            topology[topology_function][line_sep[0]]['atom']+=[line_sep[1]]
                                            topology[topology_function][line_sep[0]]['Con_Bd']+=[line_sep[2]]
                                            topology[topology_function][line_sep[0]]['dir']+=[int(line_sep[3])]
                                        else:
                                            topology[topology_function][line_sep[0]]={'atom':[line_sep[1]], 'Con_Bd':[line_sep[2]], 'dir':[int(line_sep[3])]}
                                    else:
                                        print('The bead connection group section is incorrect: \n'+location+'.top')
    return topology

def get_fragment_topology(residue, location):
    topology = sep_fragments_topology(location[:-4])
    g_var.res_top[residue] = {'C_TERMINAL':topology['C_TERMINAL'], 'N_TERMINAL':topology['N_TERMINAL'], \
                            'STEER':topology['STEER'], 'CHIRAL':topology['CHIRAL'], 'GROUPS':{}, 'CONNECT':topology['CONNECT']}
    with open(location, 'r') as pdb_input:
        group=topology['GROUPS']['group_max']
        atom_list=[]
        grouped_atoms={}
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('['):
                bead = strip_header(line)
                if bead in topology['GROUPS']:
                    group_temp = topology['GROUPS'][bead]
                else:
                    group_temp = group
                    group+=1
                if group_temp not in grouped_atoms:
                    grouped_atoms[group_temp]={bead:[]}
                    g_var.res_top[residue]['GROUPS'][bead]=group_temp
                else:
                    grouped_atoms[group_temp][bead]=[]
                    g_var.res_top[residue]['GROUPS'][bead]=group_temp
            if line.startswith('ATOM'):
                line_sep = pdbatom(line)
                grouped_atoms[group_temp][bead].append(line_sep['atom_number'])
            ### return backbone info for each aminoacid residue
            if bead in g_var.res_top[residue]['CONNECT']:
                if line.startswith('ATOM'):
                    line_sep = pdbatom(line)
                    if not is_hydrogen(line_sep['atom_name']):
                        atom_list.append(line_sep['atom_name'])    ### list of backbone heavy atoms
                g_var.res_top[residue]['ATOMS']=atom_list
    return grouped_atoms

def switch_num_name(dictionary, input_val, num_to_letter):
    if num_to_letter:
        dictionary = {v: k for k, v in dictionary.items()}
    if input_val in dictionary:
        return dictionary[input_val]
    else:
        return input_val        

def sort_connectivity(atom_dict, heavy_bond):
    cut_group = {}
    if len(atom_dict) > 1:
        for group in atom_dict:
            cut_group[group]={}
            group_atoms = [atom for frag in atom_dict[group] for atom in atom_dict[group][frag] if atom in heavy_bond]
            non_self_group = [x for x in atom_dict.keys() if x != group]
            for atom in group_atoms:
                for bond in heavy_bond[atom]:
                    for group_2 in non_self_group:
                        for frag in atom_dict[group_2]:                          
                            if bond in atom_dict[group_2][frag]:
                                cut_group[group][atom] = [frag]
    return cut_group

def fetch_fragment():
#### fetches the Backbone heavy atoms and the connectivity with pre/proceeding residues 
    amino_acid_itp = fetch_amino_rtp_file_location(g_var.forcefield_location+g_var.forcefield) 
    g_var.at_mass = fetch_atom_masses(g_var.forcefield_location+g_var.forcefield)
    for residue_type in [g_var.p_directories, g_var.o_directories]:
        for directory in range(len(residue_type)):
            for residue in residue_type[directory][1:]:    
                if residue not in g_var.res_top:
                    location = fragment_location(residue)
                    g_var.hydrogen[residue], g_var.heavy_bond[residue], residue_list, at_mass, amide_h = fetch_bond_info(residue, amino_acid_itp, g_var.at_mass, location)
                    grouped_atoms = get_fragment_topology(residue, location)
                    g_var.sorted_connect[residue]  = sort_connectivity(grouped_atoms, g_var.heavy_bond[residue])
                    g_var.res_top[residue]['RESIDUE'] = residue_list
                    g_var.res_top[residue]['atom_masses'] = at_mass
                    if residue in g_var.p_residues:
                        g_var.res_top[residue]['amide_h'] = amide_h
    for directory in range(len(g_var.np_directories)):
        for residue in g_var.np_directories[directory][1:]:    
            if residue not in g_var.res_top:
                location = fragment_location(residue)
                if residue in ['SOL','ION']: 
                    at_mass = fetch_atoms_water(g_var.np_directories[directory][0]+residue+'/', g_var.at_mass)
                    if residue == 'ION':
                        create_ion_list(location[:-4]+'.pdb')
                    g_var.hydrogen[residue], g_var.heavy_bond[residue] = {},{}
                    residue_list = [residue]
                else:
                    g_var.hydrogen[residue], g_var.heavy_bond[residue], residue_list, at_mass, amide_h  = fetch_bond_info(residue, [location[:-4]+'.itp'], g_var.at_mass, location)

                grouped_atoms = get_fragment_topology(residue, location) 
                g_var.res_top[residue]['RESIDUE'] = residue_list
                g_var.res_top[residue]['atom_masses'] = at_mass
                if residue in ['SOL', 'ION']: 
                    g_var.sorted_connect[residue]={}
                else:
                    g_var.sorted_connect[residue]  = sort_connectivity(grouped_atoms, g_var.heavy_bond[residue])

def atom_bond_check(line_sep):
    if line_sep[1] == 'atoms':
        return True, False
    elif line_sep[1] == 'bonds':
        return False,True
    else:
        return False, False

def fetch_amino_rtp_file_location(forcefield_loc):
    rtp=[]
    for file in os.listdir(forcefield_loc):
        if file.endswith('.rtp'):#['aminoacids.rtp', 'merged.rtp']:
            rtp.append(forcefield_loc+'/'+file)
    return rtp

def fetch_atom_masses(forcefield_loc):
    at_mass = {}
    if os.path.exists(forcefield_loc+'/atomtypes.atp'):
        with open(forcefield_loc+'/atomtypes.atp', 'r') as itp_input:
            for line in itp_input.readlines():
                if len(line.split()) >= 2 and not line.startswith(';'):
                    line_sep = line.split()       
                    at_mass[line_sep[0]]=line_sep[1] 
    else:
        sys.exit('cannot find atomtypes.dat in the forcefield: '+forcefield_loc)
    return at_mass

def fetch_atoms_water(forcefield_loc, at_mass_p):
    at_mass = {}
    for file in os.listdir(forcefield_loc):
        if file.endswith('itp'):
            with open(forcefield_loc+file, 'r') as itp_input:
                for line in itp_input.readlines():
                    line_sep = line.split()
                    if len(line_sep) >= 3:
                        if line[0] not in [';', '#']:
                            line_sep = line.split()

                            if line_sep[0] == '[' and line_sep[1] != 'atoms':
                                strip_atoms = False
                            if strip_atoms:
                                at_mass[line_sep[4]] = float(at_mass_p[line_sep[1]])
                            if line_sep[0] == '[' and line_sep[1] == 'atoms':
                                strip_atoms = True
    return at_mass


def fetch_bond_info(residue, rtp, at_mass,location):
    bond_dict=[]
    heavy_dict, H_dict=[],[]
    residue_present = False
    atom_conversion = {}
    residue_list=[]
    res_at_mass = {}
    for rtp_file in rtp:
        with open(rtp_file, 'r') as itp_input:
            for line in itp_input.readlines():
                if len(line.split()) >= 2 and not line.startswith(';'):
                    line_sep = line.split()
                    if line_sep[1] == residue:
                        residue_present = True
                    elif line_sep[1] in ['HSE', 'HIE'] and residue == 'HIS': 
                        residue_present = True
                    elif residue_present or residue not in g_var.p_residues+g_var.o_residues:
                        if line_sep[0] == '[':
                            atoms, bonds = atom_bond_check(line_sep)
                        elif atoms:
                            if residue in g_var.p_residues or residue in g_var.o_residues:
                                if residue not in residue_list:
                                    residue_list.append(residue)
                                atom_conversion[line_sep[0]]=int(line_sep[3])+1
                                if is_hydrogen(line_sep[0]):
                                    H_dict.append(line_sep[0])
                                else:
                                    res_at_mass[line_sep[0]] = float(at_mass[line_sep[1]])
                                    heavy_dict.append(line_sep[0])
                            else:
                                if line_sep[3] not in residue_list:
                                    residue_list.append(line_sep[3])
                                if is_hydrogen(line_sep[4]):
                                    H_dict.append(int(line_sep[0]))
                                else:
                                    res_at_mass[line_sep[4]] = float(line_sep[7])
                                    heavy_dict.append(int(line_sep[0]))
                        elif bonds:
                            try:
                                bond_dict.append([int(line_sep[0]),int(line_sep[1])])
                            except BaseException:
                                bond_dict.append([line_sep[0],line_sep[1]])
                        elif not atoms and not bonds and residue in g_var.p_residues+g_var.o_residues:
                            break
        if len(heavy_dict) > 0:
            break
    bond_dict=np.array(bond_dict)
    hydrogen = {}
    heavy_bond = {}
    if residue not in g_var.mod_residues:
        atom_conversion = get_atomistic(location)
    for bond in bond_dict:
        hydrogen, amide_h = add_to_topology_list(bond[0], bond[1], hydrogen, heavy_dict, H_dict, atom_conversion, residue, g_var.p_residues+g_var.o_residues)
        if residue in g_var.p_residues and amide_h != None:
            amide_hydrogen = amide_h
        heavy_bond, amide_h= add_to_topology_list(bond[0], bond[1], heavy_bond, heavy_dict, heavy_dict, atom_conversion, residue, g_var.p_residues+g_var.o_residues)

    if 'amide_hydrogen' in locals():
        return hydrogen, heavy_bond, residue_list, res_at_mass, amide_hydrogen
    else:
        return hydrogen, heavy_bond, residue_list, res_at_mass, None

def get_atomistic(frag_location):
#### read in atomistic fragments into dictionary    
    residue = {} 
    with open(frag_location, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            if line.startswith('ATOM'):
                line_sep = pdbatom(line) ## splits up pdb line
                residue[line_sep['atom_name']] = line_sep['atom_number']
    return residue


def add_to_topology_list(bond_1, bond_2, top_list, dict1, dict2, conversion, residue, p_residues):
    amide_hydrogen = None

    for bond in [[bond_1, bond_2], [bond_2, bond_1]]:
        if bond[0] in dict1 and bond[1] in dict2:
            if residue in p_residues:
                if bond[0] == 'N' and is_hydrogen(bond[1]):
                    amide_hydrogen = bond[1]
                if bond[0] in conversion and bond[1] in conversion:
                    bond[0], bond[1] = conversion[bond[0]],conversion[bond[1]] 
            if bond[0] not in top_list:
                top_list[bond[0]]=[]
            top_list[bond[0]].append(bond[1])
    return top_list, amide_hydrogen

def fragment_location(residue):  
    # print(database_dir)
#### runs through dirctories looking for the atomistic fragments returns the correct location
    for res_type in [g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories]:
        # print(g_var.database_locations)
        for directory in range(len(res_type)):
            if os.path.exists(res_type[directory][0]+residue+'/'+residue+'.pdb'):
                return res_type[directory][0]+residue+'/'+residue+'.pdb'            
    sys.exit('cannot find fragment: '+residue+'/'+residue+'.pdb')


def read_database_directories():
#### Read in forcefields provided
    available_provided_database=[]
    for directory_type in ['forcefields', 'fragments']:
        if os.path.exists(g_var.database_dir+directory_type):
            for root, dirs, files in os.walk(g_var.database_dir+directory_type):
                available_provided = [x for x in sorted(dirs) if not x.startswith('_')]
                break
        else:
            sys.exit('no '+directory_type+' found')
            available_provided=[]
        available_provided_database.append(available_provided)
    g_var.forcefield_available, g_var.fragments_available = available_provided_database[0], available_provided_database[1]


def database_selection(provided, selection_type):
#### print out selection of forcefields
    print('\n\n{0:^45}\n'.format('Provided '+selection_type))
    print('{0:^20}{1:^30}'.format('Selection',selection_type))
    print('{0:^20}{1:^30}'.format('---------','----------'))
    for force_num_prov, line in enumerate(provided):
        print('{0:^20}{1:^30}'.format(force_num_prov,line.split('.')[0]))
    return ask_database(provided,  selection_type)

def ask_database(provided, selection_type):
#### ask which database to use
    while True:
        try:
            if len(provided)==1:
                print('\nOnly 1 '+selection_type[:-1]+' database is currently available, therefore you have no choice but to accept the following choice.')
                return 0
        #### if asking about fragments accept a list of libraries 
            if selection_type=='fragments': 
                number = np.array(input('\nplease select fragment libraries (in order of importance: eg. "1 0" then ENTER): ').split())
                number=number.astype(int)
                if len(number[np.where(number >= len(provided))]) == 0:
                    return number
        #### if forcefield only accept one selection
            else:
                number = int(input('\nplease select a forcefield: '))
                if number < len(provided):
                    return number
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except BaseException:
            print("Oops!  That was a invalid choice")

def fetch_frag_number(fragments_available):
    fragment_number = []
    if g_var.fg != None and len(g_var.fg) > 0 :
        for frag in g_var.fg:
            if frag in fragments_available:
                fragment_number.append(fragments_available.index(frag))
            elif not g_var.info:
                print('Cannot find fragment library: '+frag+' please select library from below\n')
                fragment_number += database_selection(fragments_available, 'fragments').tolist()
    else:
        fragment_number = database_selection(fragments_available, 'fragments')
    if len(fragment_number) > 0:
        return fragment_number 
    else:
        if g_var.info:
            return []
        sys.exit('no fragment databases selected')

def add_to_list(root, dirs, list_to_add, residues):
    list_to_add.append([])
    list_to_add[-1].append(root)
    list_to_add[-1] += dirs
    list1 = [x for x in list_to_add[-1] if not x.startswith('_')]
    list1.sort()
    list_to_add[-1] = list1
    residues += list_to_add[-1][1:]
    residues.sort()    

def fetch_residues(frag_dir, fragments_available_prov, fragment_number):
#### list of directories and water types  [[root, folders...],[root, folders...]]
#### run through selected fragments
    for database in fragment_number:
        if not g_var.info:
            if g_var.database_dir in frag_dir[database]:
                print('\nYou have selected the fragment library: '+fragments_available_prov[database])
            else:
                print('\nYou have chosen to use your own fragment library: '+fragments_available_prov[database]+' in '+frag_dir[database])
    #### separate selection between provided and user
        location = frag_dir[database]+ fragments_available_prov[database]
    #### runs through protein and non protein
        for directory_type in ['/non_protein/', '/protein/', '/other/']:
    #### adds non protein residues locations to np_directories
            if os.path.exists(location+directory_type):
                for root, dirs, files in os.walk(location+directory_type):
                    if directory_type =='/non_protein/':
                        add_to_list(root, dirs, g_var.np_directories, g_var.np_residues)
                    elif directory_type =='/other/':
                        add_to_list(root, dirs, g_var.o_directories, g_var.o_residues)
                    #### adds protein residues locations to p_directories
                    elif directory_type =='/protein/':
                        add_to_list(root, dirs, g_var.p_directories, g_var.p_residues)
                    #### adds modified residues to mod directories and removes MOD from p_directories
                        if os.path.exists(location+directory_type+'MOD/'):
                            g_var.p_directories[-1].remove('MOD')
                            g_var.p_residues.remove('MOD')
                            for root_mod, dirs_mod, files_mod in os.walk(location+directory_type+'MOD/'):
                                add_to_list(root_mod, dirs_mod, g_var.p_directories, g_var.p_residues)
                                add_to_list(root_mod, dirs_mod, g_var.mod_directories, g_var.mod_residues)
                                break
                    break
    

def print_water_selection(water, directory):
    if g_var.w != None:
        print('\nThe water type '+g_var.w+' doesn\'t exist')
    if len(water) == 0:
        sys.exit('\nCannot find any water models in: \n\n'+directory[0]+'SOL/'+'\n')
    print('\nPlease select a water molecule from below:\n')
    print('{0:^20}{1:^30}'.format('Selection','water_molecule'))
    print('{0:^20}{1:^30}'.format('---------','----------'))
    for selection, water_model in enumerate(water):
        print('{0:^20}{1:^30}'.format(selection,water_model))


def ask_for_water_model(directory, water):
    while True:
        try:
            number = int(input('\nplease select a water model: '))
            if number < len(water):
                return directory[0]+'SOL/', water[number]
        except KeyboardInterrupt:
            sys.exit('\nInterrupted')
        except BaseException:
            print("Oops!  That was a invalid choice")

def check_water_molecules():
    if any('SOL' in sublist for sublist in g_var.np_directories): 
        water=[]
        for directory in g_var.np_directories:
            if os.path.exists(directory[0]+'SOL/SOL.pdb'):
                g_var.water_info.append([directory[0]+'SOL'])
                with open(directory[0]+'SOL/SOL.pdb', 'r') as sol_input:
                    for line_nr, line in enumerate(sol_input.readlines()):
                        if line.startswith('['):
                            water.append(strip_header(line))
                            g_var.water_info[-1].append(strip_header(line))
        if not g_var.info:
            if g_var.w in water:
                print('\nYou have selected the water model: '+g_var.w)
                g_var.water_dir, g_var.water = directory[0]+'SOL/', g_var.w
            else:
                print_water_selection(water, directory)
                g_var.water_dir, g_var.water = ask_for_water_model(directory, water)
            if g_var.w is None:
                g_var.opt['w'] = water

############################################################################################## fragment rotation #################################################################################

def AnglesToRotMat(theta) :
#### rotation matrices for the rotation of fragments. theta is [x,y,z] in radians     
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
         
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
                                        
    R = np.dot(R_z, np.dot( R_y, R_x ))
    return R                          

def angle_clockwise(A, B):
#### find angle between vectors
    AB = np.linalg.norm(A)*np.linalg.norm(B)
    A_dot_B = A.dot(B)
    angle = np.degrees(np.arccos(A_dot_B/AB))
#### determinant of A, B
    determinant = np.linalg.det([A,B])

    if determinant < 0: 
        return angle
    else: 
        return 360-angle

############################################################################################## fragment rotation done #################################################################################

def pdbatom(line):
### get information from pdb file
### atom number, atom name, residue name,chain, resid,  x, y, z
    try:
        return dict([('atom_number',int(line[7:11].replace(" ", ""))),('atom_name',str(line[12:16]).replace(" ", "")),('residue_name',str(line[16:21]).replace(" ", "")),\
            ('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54]))])
    except BaseException:
        print(line[30:38],line[38:46],line[46:54])
        sys.exit('\npdb line is wrong:\t'+line) 

def create_pdb(file_name):
    pdb_output = open(file_name, 'w')
    pdb_output.write('TITLE     GENERATED BY CG2AT\nREMARK    Please don\'t explode\nREMARK    Good luck\n\
'+g_var.box_vec+'MODEL        1\n')
    return pdb_output

def mkdir_directory(directory):
#### checks if folder exists, if not makes folder
    if not os.path.exists(directory):
        os.mkdir(directory)

def clean():
#### cleans temp files from residue_types
    print()
    for residue_type in g_var.cg_residues:
        if residue_type not in ['SOL', 'ION']:
            print('Cleaning temp files from : '+residue_type)
            os.chdir(g_var.working_dir+residue_type)
            file_list = glob.glob('*temp*', recursive=True)
            file_list += glob.glob(residue_type+'*pdb', recursive=True)
            for file in file_list:
                if not file.endswith('.tpr') and not file.endswith('_merged.pdb') and 'gmx' not in file:
                    os.remove(file)
            os.chdir(g_var.working_dir+residue_type+'/MIN')
            file_list = glob.glob('*temp*', recursive=True)
            file_list += glob.glob('*trr', recursive=True)
            for file in file_list:
                if not file.endswith('.tpr'):
                    os.remove(file) 

def fix_time(t1, t2):
    minutes, seconds= divmod(t1-t2, 60)
    if minutes > 60:
        hours, minutes = divmod(minutes, 60)
    else:
        hours = 0
    return '{0:^3}{1:^6}{2:^3}{3:^4}{4:^3}{5:^4}'.format(int(np.round(hours)),'hours',int(np.round(minutes)),'min',int(np.round(seconds,0)),'sec') 

def print_script_timings():
    to_print=[]
    to_print.append('\n{:-<100}'.format(''))
    to_print.append('\n{0:^47}{1:^22}'.format('Job','Time'))
    to_print.append('{0:^47}{1:^22}'.format('---','----'))
    to_print.append('\n{0:47}{1}'.format('Initialisation: ', fix_time(g_var.tc['i_t_e'],g_var.tc['i_t'])))
    to_print.append('{0:47}{1}'.format('Read in CG system: ', fix_time(g_var.tc['r_i_t'],g_var.tc['i_t_e']))) 
    if 'PROTEIN' in g_var.system:
        to_print.append('{0:47}{1}'.format('Build protein systems: ',fix_time(g_var.tc['f_p_t'],g_var.tc['r_i_t'])))
    if 'OTHER' in g_var.system:
        to_print.append('{0:47}{1}'.format('Build other systems: ',fix_time(g_var.tc['f_o_t'],g_var.tc['f_p_t'])))        
    to_print.append('{0:47}{1}'.format('Build non protein system: ',fix_time(g_var.tc['n_p_t'],g_var.tc['f_o_t'])))
    to_print.append('{0:47}{1}'.format('merge and minimise de novo: ',fix_time(g_var.tc['m_t'],g_var.tc['n_p_t'])))
    to_print.append('{0:47}{1}'.format('NVT on de novo: ',fix_time(g_var.tc['eq_t'],g_var.tc['m_t'])))
    if g_var.o in ['all', 'align'] and g_var.a != None and g_var.user_at_input:
        to_print.append('{0:47}{1}'.format('Creating aligned system: ',fix_time(g_var.tc['a_e'],g_var.tc['a_s'])))
    to_print.append('{:-<69}'.format(''))
    to_print.append('{0:47}{1}'.format('Total run time: ',fix_time(g_var.tc['f_t'],g_var.tc['i_t'])))
    with open(g_var.final_dir+'script_timings.dat', 'w') as time_out:  
        for line in to_print:
            time_out.write(line+'\n')
            if g_var.v >= 1:
                print(line)
        print('\nAll script timings have been saved in: \n'+g_var.final_dir+'script_timings.dat\n')

def cg2at_header():
    print('{0:30}'.format('\nCG2AT2 is a fragment based conversion of coarsegrain to atomistic.\n'))
    print('{0:^90}\n'.format('CG2AT2 version: '+str(g_var.version)))
    print('{0:^90}'.format('CG2AT2 is written by Owen Vickery'))
    print('{0:^90}'.format('Project leader Phillip Stansfeld'))
    print('\n{0:^90}\n{1:^90}'.format('Contact email address:','owen.vickery@warwick.ox.ac.uk'))
    print('\n{0:^90}\n{1:^90}\n{2:^90}\n{3:-<90}'.format('Address:','School of Life Sciences, University of Warwick,','Gibbet Hill Road, Coventry, CV4 7AL, UK', ''))
    print('\n{0:^90}\n{1:^90}'.format('If you are using this script please acknowledge me (Dr Owen Vickery)','and cite the following DOI: 10.5281/zenodo.3890163'))    
    print('\n{0:^90}'.format('Executable: '+g_var.opt['input'].split()[0]))
    print('{0:^90}'.format('Database locations: '+g_var.database_dir))
    print('{0:^90}\n\n{1:-<90}'.format('Script locations: '+g_var.scripts_dir, ''))

def database_information():
    
    print('\n{0:^90}\n{1:-<90}\n'.format('The available forcefields within your database are (flag -ff):', ''))
    for forcefields in g_var.forcefield_available:
        print('{0:^90}'.format(forcefields))
    print('\n\n{0:^90}\n{1:-<90}\n'.format('The available fragment libraries within your database are (flag -fg):', ''))
    for fragments in g_var.fragments_available:
        print('{0:^90}'.format(fragments))    
    if g_var.fg != None :
        fragments_in_use()
    sys.exit('\n\"If all else fails, immortality can always be assured by spectacular error.\" (John Kenneth Galbraith)\n')

def fragments_in_use():
    protein_directories=[]
    if np.any([g_var.np_directories, protein_directories, g_var.mod_directories, g_var.o_directories, g_var.water_info]):
        for database_val, database in enumerate(sorted(g_var.fg)):
            if len(g_var.p_directories) > 0:
                protein_directories = [sublist for sublist in g_var.p_directories if path_leaf(sublist[0])[1] != 'MOD']
            print('\n\n{0:^90}\n{1:-<90}\n'.format('The following residues are available in the database: '+database,''))
            res_type_name = ['Non protein residues', 'Protein residues', 'Modified protein residues', 'Other linked residues', 'Water residues']
            for res_val, residue in enumerate([g_var.np_directories, protein_directories, g_var.mod_directories, g_var.o_directories, g_var.water_info]):
                try:
                    res_type = sorted(residue[database_val][1:])
                    print('\n{0:^90}\n{1:^90}'.format(res_type_name[res_val], '-'*len(res_type_name[res_val])))
                    if len(', '.join(map(str, res_type))) <= 80:
                        print('{0:^90}'.format(', '.join(map(str, res_type))))
                    else:
                        start, end = 0, 1                       
                        while end < len(res_type):
                            line = ', '.join(map(str, res_type[start:end]))
                            while len(line) <= 80:
                                if end < len(res_type):
                                    end+=1
                                    line = ', '.join(map(str, res_type[start:end]))
                                    if len(line) > 80:
                                        end-=1
                                        line = ', '.join(map(str, res_type[start:end]))
                                        break
                                else:
                                    break
                            print('{0:^90}'.format(line))
                            start = end
                except:
                    pass
        print('\n{0:-<90}\n'.format(''))

def write_system_components():
    print('\n{:-<100}'.format(''))
    print('{0:^100}'.format('Script has completed, time for a beer'))
    print('\n{0:^10}{1:^25}'.format('molecules','number'))
    print('{0:^10}{1:^25}'.format('---------','------'))
    for section in g_var.system:
        print('{0:^10}{1:^25}'.format(section, g_var.system[section]))
