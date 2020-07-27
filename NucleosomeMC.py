# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:43:39 2015

@author: John van Noort

Add functionalities for nucleoosmes to HelixMC
"""

import matplotlib as mpl

mpl.use(u'TkAgg')
mpl.interactive(False)

import numpy as np
from helixmc.pose import HelixPose
from helixmc.util import frames2params_3dna, frames2params
import matplotlib.pyplot as plt
import time
# ChromatinMC modules:
import FileIO as fileio

plt.interactive(False)
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})


def ofs2params(of1, of2, _3dna=False, flipx=[0,0]):
    o1 = of1[0]
    f1 = of1[1:, :] - o1
    o2 = of2[0]
    f2 = of2[1:, :] - o2

    if flipx[0]:
        f1[1] *= -1
        f1[2] *= -1
    if flipx[1]:
        f2[1] *= -1
        f2[2] *= -1

    if _3dna:
        params = frames2params_3dna(o1, o2, f1, f2)
    else:
        params = frames2params(o1, o2, f1, f2)
    return params


def find(lst, predicate):
    return (i for i, j in enumerate(lst) if predicate(j)).next()


def of2axis(of, length=[60, 90, 120]):
    """
    converts originframe to axis for plotting purposes
    """
    origin = of[0]
    frame = of[1:] - origin
    coord_out = []

    for ax_length, ax_direction in zip(length, frame):
        for j in np.linspace(0, ax_length, np.abs(ax_length)):
            coord_out.append(origin + j * ax_direction)
    return np.asarray(coord_out)


def apply_transformation(coord, trans_def):
    """
    Apply rigid transformation to coord using transformation
    parameters obtained with Kabsch algorithm
    
    Parameters
    ----------
    coord : ndarray of (N,3)
    trans_def : list containing center of rotation, rotation matrix and
                center after rotation
    
    Returns
    -------
    coord_out: ndarray of (N,3)    
    """
    N = coord.shape[0]
    coord_out = coord - np.tile(trans_def[0], (N, 1))
    coord_out = np.dot(coord_out, trans_def[1])
    coord_out = coord_out + np.tile(trans_def[2], (N, 1))

    return np.asarray(coord_out)


def get_transformation(start, target=None):
    """
    Align de coordinates defined by P such that
    they will overlap with Q, using Kabsch algorithm
    See http://en.wikipedia.org/wiki/Kabsch_algorithm
    """
    P = start
    Q = target
    if Q is None:
        Q = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    Pc = sum(P) / len(P)
    Qc = sum(Q) / len(Q)
    C = np.dot(np.transpose(P - Pc), Q - Qc)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return Pc, U, Qc


def join_o_f(origin, frame):
    of = []
    of.append(origin)
    for f in np.transpose(frame):
        of.append(origin + f)
    return np.asarray(of)


def get_nuc_of(coord, frames, dyad, nucl):
    """
    Calculate the center of mass and the reference frame (=of) of a nucleosome in dna
    """
    tf = get_transformation(nucl.dyad_of, target=get_of_2(coord, frames, dyad))
    n_of = apply_transformation(nucl.of, tf)
    return n_of


def get_of(dna, i):
    of = get_of_2(dna.coord, dna.frames, i)
    return of


def get_of_2(coord, frames, i):
    return np.concatenate(([coord[i]], frames[i].T+coord[i]), axis = 0)


# def get_wrap_params2(dna, dyad, fixed):
#     fixed_params = []
#     dyad_of = get_of(dna, dyad)
#     for i in fixed:
#         base_of = get_of(dna, dyad + i)
#         fixed_params.append((dyad_of, base_of))
#     return np.asarray(fixed_params).reshape((-1, 6))


def get_wrap_param(dna_coord, dna_frames, dyad, fixed):
    fixed_params = []
    for i in fixed:
        params = frames2params(dna_coord[dyad], dna_coord[dyad + i], dna_frames[dyad], dna_frames[dyad + i])
        fixed_params.append(params)
    return np.asarray(fixed_params)


def seq3_to_1(seq3):
    '''
    Turn a three letter protein into a one letter protein.
    Works also for DNA 'DA ' strings

    Parameters
    ----------
    seq3 : str(4) aminoacid or DNA code
    
    Returns
    -------
    seq1 : str

    '''
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER': '*',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'XAA': 'X',
         'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', '   ': ''}
    seq3 = seq3.strip()
    seq1 = ''
    for i in range(0, len(seq3), 4):
        seq1 += d[seq3[0 + i:3 + i].strip()]
    return seq1


def read_pdb(pdb_file):
    '''
    Get Chains, Chain type and coord of proteins from pdb file 

    Parameters
    ----------
    pdb_file : string
    
    Returns
    -------
    chain_dict: {chain, type, coord}

    '''
    #    print '###>'+ pdb_file
    f = open('PDBs\\'+pdb_file, 'r')
    pdb = f.readlines()
    f.close()

    cpd = [l for l in pdb if l[0:6] == "COMPND"]
    chains = []
    molecules = []
    keywords = ['DNA', 'H2A', 'H2B', 'H3', 'H4']
    for l in cpd:
        s = l[11:].split(':')
        if s[0] == 'MOLECULE':
            for k in keywords:
                if s[1].find(k) >= 0:
                    if k in chains:
                        chains.append(k + '*')
                    else:
                        chains.append(k)
        if s[0] == 'CHAIN':
            s = s[1].split(';')[0].split(',')
            i = 0
            for m in s:
                molecules.append(m.lstrip())
                if i > 0:
                    chains.append(chains[-1] + '*')
                i += 1

    chain_dict = dict([(c, ['', '', np.zeros((1, 3)), np.zeros(1)]) for c in chains])
    i = 0
    for c in chains:
        chain_dict[c][0] = molecules[i]
        i += 1
    chain_dict.pop('DNA*')

    # Use SEQRES for sequence 
    seq = [l for l in pdb if l[0:6] == "SEQRES"]
    for l in seq:
        for i in chain_dict:
            if chain_dict[i][0] == l[11]:
                chain_dict[i][1] += l[19:71]

    for i in chain_dict:
        chain_dict[i][1] = seq3_to_1(chain_dict[i][1])
        chain_dict[i][2] = np.zeros([len(chain_dict[i][1]), 3])
    # removed /0 from "chain_dict[i][2] = np.zeros([len(chain_dict[i][1]),3])/0"
    # -bert

    # Grab all CA from ATOM entries
    # Grab B -factor for Phosphates
    B_factor = np.zeros(len(chain_dict['DNA'][2]) - 1)
    P_I = np.zeros((len(chain_dict['DNA'][2]) - 1, 3))

    start_i = None
    for l in pdb:
        if l[0:6] == "ATOM  " and l[13:16] == "CA ":
            nc = np.fromstring(l[30:55], sep=' ')
            for i in chain_dict:
                if chain_dict[i][0] == l[21]:
                    chain_dict[i][2][int(l[22:26]) - 1, :] = nc
        if l[0:6] == "ATOM  " and l[13:16] == "P  ":
            if start_i == None:
                start_i = int(l[22:27])
            B_factor[int(l[22:27]) - start_i] += float(l[61:66])
            if l[21] == 'I':
                P_I[int(l[22:27]) - start_i] = np.fromstring(l[30:55], sep=' ')

    #    chain_dict['DNA'][3]=B_factor
    av_I = np.mean(P_I, axis=0)

    dI = np.sqrt(np.sum((P_I - av_I) ** 2, axis=1))
    chain_dict['DNA'][3] = dI
    #    plt.plot(np.arange(len(dI)), dI)
    #
    #    plt.plot(np.arange(len(B_factor)), B_factor)
    #    plt.show()
    return chain_dict


def crossproduct(u, v):
    w1 = u[1] * v[2] - u[2] * v[1]
    w2 = u[2] * v[0] - u[0] * v[2]
    w3 = u[0] * v[1] - u[1] * v[0]
    return np.array([w1, w2, w3])

class NucPose(object):
    '''
    NucPose object storing the state info of the nucleosome.

    Parameters
    ----------
    
    Attributes
    ----------
    nuc_type : string
        Name of the pdb file that contains teh coordinates
    step_list : list
        List of nucleotide-nucleotide steps (eg AC, AG etc)    
    dna : HelixMC pose
        nucleosomal DNA in HelixMC framework
    dyad : int
        index of dyad bp
    chains : dictionary
        dictionary of all chains {name, sequence, coord}
    fixed_i : list
        indices of fixed basepairs relative to dyad
    l_coord : ndarray of (4,3)
        The coordinates of Glu61 (H2A) and Asp24 (H4) that mediate nucleosome-
        nucleosome interactions through the tail of H4
    '''

    def __init__(self):
        return

    @classmethod
    def from_file(self, input_file):
        '''
        Load pose data from an input file.

        Parameters
        ----------
        input_file : str
        input file (.txt) obtained from w3DNA.

        '''
        self.nuc_type = input_file.split('.')[0]
        if input_file == None:
            input_file = "3LZ0.3DNA"
            print('Default file: ' + input_file)
        #   read protein coord from pdb file
        filename = input_file
        filename = fileio.change_extension(filename, 'pdb')
        chains = read_pdb(filename)

        #   read 3DNA file
        filename = fileio.change_extension(filename, '3DNA')
        with open('PDBs\\'+filename) as f:
            lines = f.read().splitlines()

        #   Basepair step parameters
        i = find(lines, lambda x:
        '    step       Shift     Slide      Rise      Tilt      Roll     Twist'
        in x) + 1

        j = find(lines[i:], lambda x: '~~~~~~~~~~~' in x)

        params = np.zeros((j, 6))
        for k in range(0, j):
            if not ('----' in lines[i + k]):
                tmp = np.asarray(map(float, lines[i + k].split()[2:]))
            params[k:] = np.append(tmp[0:3], np.radians(tmp[3:6]))

        #   get fixed frames
        B_factor = chains['DNA'][3]
        fixed = []
        for i in range(14):
            j = np.argmin(B_factor)
            fixed.append(j)
            B_factor[j - 6:j + 6] = np.max(B_factor)

        fixed = np.sort(fixed)
        self.dyad = int(round(np.mean(fixed)))
        self.fixed = fixed - self.dyad

        #   Centers of basepairs
        i = find(lines, lambda x:
        '      bp        Ox        Oy        Oz        Nx        Ny        Nz'
        in x) + 1

        coord = np.zeros((j, 3))
        frames = np.zeros((j, 3, 3))
        for k in range(0, j):
            tmp = np.asarray(map(float, lines[i + k].split()[2:]))
            coord[k, :] = tmp[0:3]
            frames[k][0] = tmp[3:6] / np.linalg.norm(tmp[3:6])
        chains['DNA'][2] = np.asarray(coord)

        # rotate all coord, such that frame0 = origin
        self.dna = HelixPose(params)
        tm = get_transformation(chains['DNA'][2][0:100], self.dna.coord[0:100])
        for chain in chains:
            chains[chain][2] = apply_transformation(chains[chain][2], tm)
        self.chains = chains

        # remove non-nucleosomal DNA
        if len(self.dna.coord) > 147:
            loose_ends = 147 - (self.fixed_i[-1] - self.fixed_i[0] + 1)
            if loose_ends % 2 == 1:
                l1 = loose_ends / 2
                l2 = l1 + 1
            else:
                l1 = loose_ends / 2
                l2 = l1
            start = self.fixed_i[0] - l1 + self.d_index
            end = self.fixed_i[-1] + l2 + self.d_index
            self.d_index = self.d_index - start
            params = params[start:end]
            self.dna = HelixPose(params)

        # get origin and frame of nucleosome
        cm = np.mean(self.dna.coord[self.fixed[4:-4]], axis=0)
        Nx = self.dna.coord[self.dyad] - cm
        Nx = Nx / np.linalg.norm(Nx)
        Nz = np.mean(self.dna.coord[self.fixed[:7], :], axis=0) - np.mean(self.dna.coord[self.fixed[7:], :], axis=0)
        Ny = np.cross(Nx, Nz)
        Ny = Ny / np.linalg.norm(Nz)
        Nz = np.cross(Nx, Ny)
        Nz = Nz / np.linalg.norm(Nz)
        origin = cm
        frame = np.array([Nx, Ny, Nz])
        self.of = join_o_f(origin, np.transpose(frame))
        self.dyad_of = get_of(self.dna, self.dyad)

        # get link coordinates Glu61 (H2A) and Asp24 (H4)
        # make dictionary with locations of Glu61 and Asp24
        int_dict = {'H2A': 60, 'H2A*': 60, 'H4': 23, 'H4*': 23}
        self.l_coord = []
        for locus in int_dict:
            self.l_coord.append(chains[locus][2][int_dict[locus]])
        # print('before: ', self.l_coord)
        # convert list into array
        # self.l_coord = np.asarray(self.l_coord)
        # print('after: ', self.l_coord)

def main():
    nuc = NucPose()
    nuc.from_file('1KX5.3DNA')

    coord = []
    for chain in nuc.chains:
        coord.append(nuc.chains[chain][2])

    # tf = get_transformation(nuc.of, target=np.asarray([[0,0,0],[0,0,-1],[0,1,0],[1,0,0]]))
    # tf = get_transformation(nuc.of, target=np.asarray([[0,0,0],[-1,0,0],[0,1,0],[0,0,-1]]))
    tf = get_transformation(nuc.of)

    n_coord = []
    for c in coord:
        n_coord.append(apply_transformation(c, tf))

    nuc_ax = apply_transformation(nuc.of, tf)
    n_coord.append(of2axis(nuc_ax))

    # add link coordinates to n_coords
    for l in nuc.l_coord:
        l.resize(1,3)
        n_coord.append(apply_transformation(l, tf))

    filename = fileio.get_filename(root='1nuc', ext='pov', incr=True)
    print(fileio.create_pov(filename, n_coord, range_A=[250, 350], offset_A=[0, 0, 150], show=True, width_pix=1500))

if __name__ == '__main__':
    main()
