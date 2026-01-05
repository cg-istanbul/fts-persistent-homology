############################################
# Homology computation (2-skeleton)
############################################

from sage.topology.simplicial_complex import SimplicialComplex

def homology_ranks_from_simplices(crosscut, maxdim=1):
    """
    simplices: list of tuples, e.g. [(0,), (1,), (0,1), (0,1,2)]
    Returns Betti numbers up to maxdim
    """
    if not crosscut:
        return {0: 0, 1: 0}

    K = SimplicialComplex(crosscut)
    betti = K.betti()

    return {
        0: betti.get(0, 0),
        1: betti.get(1, 0),
    }
