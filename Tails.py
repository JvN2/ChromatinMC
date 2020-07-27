import numpy as np

# ChromatinMC modules:
import NucleosomeMC as nMC

def tf_dyad(dyads, dna, nucl):
    '''
    get transformation matrix for every dyad

    Parameters
    ----------
    dyads:  indices of dyads in fiber
    dna:    HelixPose
    nucl:   nucleosome pose

    Returns
    -------
    tf: transformation matrix, for each dyad

    '''
    of_d_nucl = nMC.get_of(nucl.dna, nucl.dyad)             # origin frame dyad in nucl pose
    of_d_fiber = []                                         # origin frame of dyad in fiber
    tf = []                                                 # transformation matrix
    for i, d in enumerate(dyads):
        # define origin frame of dyad in fiber
        of_d_fiber.append((nMC.get_of(dna, d)))
        # get transformation matrix of nucleosome dyad onto fiber dyad
        tf.append(nMC.get_transformation(of_d_nucl, of_d_fiber[i]))

    return tf

def get_histones(coord, dyads, dna, nucl):
    # type: (numpy.ndarray, numpy.ndarray, helixmc.pose.HelixPose, NucleosomeMC.NucPose) -> object
    '''
    projects histones coordinates onto every
    nucleosome in fiber pose, appends these
    coordinates to coord

    Parameters
    ----------
    coord:  coordinates of helixPose
    dyads:  indices of dyads in fiber
    dna:    HelixPose
    nucl:   nucleosome pose

    Returns
    -------
    coord:  [[DNA], n*[histones],[linker-coord]] n = number of nucleosomes
    radius: list
    color:  string

    '''
    # coordinates of histone proteins
    p_coord = []
    for chain in nucl.chains:
        if chain == 'DNA':
            pass
        else:
            p_coord.append(nucl.chains[chain][2])

    radius = [10]                                           # radius of DNA in POVray
    colors = 'o'                                            # color of DNA
    coord = [coord]
    tf = tf_dyad(dyads, dna, nucl)                          # transformation matrix for every dyad

    for tfm in tf:
        # apply transformation on coordinates of histone proteins
        for c in p_coord:
            coord.append(nMC.apply_transformation(c, tfm))
        # add link coordinates to coord
        for l in nucl.l_coord:
            coord.append(nMC.apply_transformation(l, tfm))
        # radius of histone proteins
        radius = np.append(radius, np.ones(8) * 4)
        # radius of linker-amino-acids
        radius = np.append(radius, np.ones(4) * 15)
        # colors of histone proteins and linker-amino-acids
        colors += 'bbggryrymzmz'           # color of DNA, 8 histone proteins + H2A, H2A*, H4, H4*

    return coord, radius, colors

def tail_dist(dyad_1, dyad_2, dyads, dna, nucl, orientation=None):
    """

    Parameters
    ----------
    dyad_1:         index of first dyad
    dyad_2:         index of second dyad
    orientation:    orientation of two nucleosomes
    dyads:          indices of all dyads
    dna:            helixpose
    nucl:           nucleosomepose

    Returns
    -------
    d_up:           distance from dyad 1 to dyad 2
    d_down:         distance from dyad 2 to dyad 1

    """
    if orientation == None:
        orientation = '*-'
    l_coord = np.asarray(nucl.l_coord)
    n_l_coord = []                                          # new coordinates of linker coordinates after transformation
    tf_matrix = tf_dyad(dyads, dna, nucl)

    for tf in tf_matrix:
        # apply transformation on coordinates of linker coordinates
        n_l_coord.append(nMC.apply_transformation(l_coord, tf))

    # orientation of tails and patches
    # * (star) and - (no star) refers to hist_int
    # '*-': top dyad 1 connects to bottom dyad 2
    # '-*': bottom dyad 1 connects to top dyad 2
    # '**': top dyad 1 connects to top dyad 2
    # '--': bottom dyad 1 connects to bottom dyad 2

    hist_int = {'H2A': 0, 'H2A*': 1, 'H4': 2, 'H4*': 3}

    tail_star_1 = n_l_coord[dyad_1][hist_int['H4*']]
    tail_stripe_1 = n_l_coord[dyad_1][hist_int['H4']]
    tail_star_2 = n_l_coord[dyad_2][hist_int['H4*']]
    tail_stripe_2 = n_l_coord[dyad_2][hist_int['H4']]

    patch_star_1 = n_l_coord[dyad_1][hist_int['H2A*']]
    patch_stripe_1 = n_l_coord[dyad_1][hist_int['H2A']]
    patch_star_2 = n_l_coord[dyad_2][hist_int['H2A*']]
    patch_stripe_2 = n_l_coord[dyad_2][hist_int['H2A']]

    if orientation == '*-':
        d_up = np.sqrt(np.sum((tail_star_1-patch_stripe_2)**2))
        d_down = np.sqrt(np.sum((tail_stripe_2-patch_star_1)**2))
    elif orientation == '-*':
        d_up = np.sqrt(np.sum((tail_stripe_1-patch_star_2)**2))
        d_down = np.sqrt(np.sum((tail_star_2-patch_stripe_1)**2))
    elif orientation == '**':
        d_up = np.sqrt(np.sum((tail_star_1-patch_star_2)**2))
        d_down = np.sqrt(np.sum((tail_star_2-patch_star_1)**2))
    elif orientation == '--':
        d_up = np.sqrt(np.sum((tail_stripe_1-patch_stripe_2)**2))
        d_down = np.sqrt(np.sum((tail_stripe_2-patch_stripe_1)**2))

    print('d up: ', d_up)
    print('d down: ', d_down)

    return d_up, d_down