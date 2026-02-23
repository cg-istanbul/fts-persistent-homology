# scripts/topology/anchors.sage

from sage.all import RR


def eligible_anchors_from_tau(rho, tau):
    t = RR(tau)
    anchors = [i for i, r in enumerate(rho) if RR(r) >= t]
    anchors.sort()
    return anchors


def anchor_neighborhood(D, a, r_k_a, lambda_scale=1.5):
    n = n_points(D)
    a = int(a)
    if r_k_a is None:
        return frozenset([a])
    rad = RR(lambda_scale) * RR(r_k_a)
    return frozenset([x for x in range(n) if dist(D, a, x) <= rad])


def build_A_V_at_tau(
    D, rho, r_k, tau, lambda_scale=1.5,
    drop_X_from_V=True
):
    n = n_points(D)
    X = frozenset(range(n))

    anchor_list = eligible_anchors_from_tau(rho, tau)

    Oa = {}
    A_sets = []
    for a in anchor_list:
        a = int(a)
        S = anchor_neighborhood(D, a, r_k[int(a)], lambda_scale=lambda_scale)
        Oa[a] = S
        A_sets.append(S)

    gens = list(A_sets)

    U_A = minimal_nonempty_intersections(gens, n)
    V_A = distinct_minimal_opens(U_A, drop_X=drop_X_from_V)

    covered_set = frozenset().union(*A_sets) if len(A_sets) > 0 else frozenset()
    uncovered_set = [x for x in range(n) if x not in covered_set]
    
    return dict(
        tau=RR(tau),
        anchor_list=list(anchor_list),
        Oa=Oa,
        A_sets=list(A_sets),
        U_A=U_A,
        V_A=V_A,
        covered_set=covered_set,
        uncovered_set=uncovered_set,
        diagnostics=dict(
            **{
                "anchors": int(len(anchor_list)),
                "covered": int(len(covered_set)),
                "uncovered": int(n - len(covered_set)),
            }
        ),
    )


def build_A_V_over_taus(
    D, rho, r_k, taus, lambda_scale=1.5,
    drop_X_from_V=True
):
    out = []
    for tau in taus:
        out.append(build_A_V_at_tau(
            D=D, rho=rho, r_k=r_k, tau=tau, lambda_scale=lambda_scale,
            drop_X_from_V=drop_X_from_V
        ))
    return out


