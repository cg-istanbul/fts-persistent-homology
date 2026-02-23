# scripts/topology/basis.sage

from sage.all import RR


def minimal_basis_of_topology_density_guided(
    n,
    A_sets,
    B_sets,
    drop_X_from_V=True
):
    n = int(n)
    X = frozenset(range(n))

    gens = []
    gens.extend([frozenset(A) for A in (A_sets or []) if A is not None])
    gens.extend([frozenset(B) for B in (B_sets or []) if B is not None])

    U_new = minimal_nonempty_intersections(gens, n)
    V_new = distinct_minimal_opens(U_new, drop_X=drop_X_from_V)

    
    return {"U": U_new, "V": V_new}

 
