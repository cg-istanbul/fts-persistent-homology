# scripts/topology/uncovered.sage

def _fs(S):
    if S is None:
        return frozenset()
    return frozenset(int(x) for x in S)


def _min_dist_point_to_set(D, x, S):
    if not S:
        return float("inf")
    row = D[int(x)]
    best = float("inf")
    for y in S:
        d = float(row[int(y)])
        if d < best:
            best = d
            if best <= 0.0:
                break
    return best


def _inclusion_minimals(sets_list):
    uniq = list(dict.fromkeys([_fs(S) for S in sets_list if S]))
    uniq.sort(key=lambda s: (len(s), sorted(list(s))[:5]))
    mins = []
    for S in uniq:
        is_min = True
        for T in mins:
            if T.issubset(S):
                is_min = False
                break
        if is_min:
            mins.append(S)
    return mins


def _minimal_containing(curr_VA, subset):
    subset = _fs(subset)
    cands = [V for V in curr_VA if subset.issubset(V)]
    if not cands:
        return frozenset()
    cands.sort(key=lambda s: (len(s), sorted(list(s))[:5]))
    return _fs(cands[0])


def saturate_B_sets_against_curr_VA(B_sets, curr_VA, X_i):

    curr_VA = [_fs(S) for S in (curr_VA or []) if S]
    B_sets = [_fs(B) for B in (B_sets or []) if B]
    X_i = set(int(x) for x in (X_i or []))

    if not B_sets or not curr_VA:
        return B_sets

    pt_to_blocks = {}
    for V in curr_VA:
        for p in V:
            pt_to_blocks.setdefault(int(p), []).append(V)

    out = []
    for B in B_sets:
        covered = set(int(p) for p in B if int(p) in X_i)
        if not covered:
            out.append(B)
            continue

        union_blocks = set()
        for p in covered:
            for V in pt_to_blocks.get(p, []):
                union_blocks.update(V)

        uncovered_part = set(int(p) for p in B if int(p) not in X_i)
        B_sat = frozenset(uncovered_part | union_blocks)
        out.append(B_sat)

    return out


def build_B_sets(
    D,               
    prev_U,
    uncovered_points,
    prev_uncovered_points,
    prev_UA,
    curr_VA,
    curr_covered_set,
):
    uncovered_now = set(int(x) for x in (uncovered_points or []))
    uncovered_prev = set(int(x) for x in (prev_uncovered_points or []))
    newly_uncovered = sorted(list(uncovered_now - uncovered_prev))
    still_uncovered = sorted(list(uncovered_now & uncovered_prev))

    X_i = set(int(x) for x in (curr_covered_set or []))  

    prev_U  = [_fs(S) for S in prev_U]
    prev_UA = [_fs(S) for S in prev_UA]
    curr_VA = [_fs(S) for S in (curr_VA or []) if S]

    def Vprev_of_point(x):
        return prev_UA[int(x)]

    # ------------------------------------------------------------
    # Helper: "is x already contained by some B-set?"
    # ------------------------------------------------------------
    B_sets = []

    def _covered_by_existing_B(x):
        xx = int(x)
        for B in B_sets:
            if xx in B:
                return True
        return False

    def _add_if_needed(x, Cand):
        xx = int(x)
        if _covered_by_existing_B(xx):
            return
        Cand = _fs(Cand)
        if not Cand:
            Cand = frozenset([xx])
        B_sets.append(Cand)

    # ------------------------------------------------------------
    # Precompute regions for newly uncovered points and identify dead regions
    # ------------------------------------------------------------
    dead_regions = []
    region_of_new = {}
    for x in newly_uncovered:
        S = _fs(Vprev_of_point(x))
        if not S:
            continue
        region_of_new[int(x)] = S
        if len(S & X_i) == 0:
            dead_regions.append(S)

    dead_minimals = _inclusion_minimals(dead_regions)
    dead_minimals_sorted = sorted(dead_minimals, key=len)
    dead_minimals_set = set(dead_minimals)

    def B_dead_region(S):
        return _fs(S)

    # ------------------------------------------------------------
    # Define K, L, M
    # K = uncovered at step (i-1)
    # L = newly uncovered points whose dead region is MINIMAL dead
    # M = remaining uncovered points at step i excluding K and L
    # ------------------------------------------------------------
    K_points = set(still_uncovered)

    L_points = set()
    for x in newly_uncovered:
        S = region_of_new.get(int(x), frozenset())
        if S and (S in dead_minimals_set):
            L_points.add(int(x))

    M_points = sorted(list(uncovered_now - K_points - L_points))

    # ------------------------------------------------------------
    # Stage 1: process K (carry U^{i-1}(x)), add only if needed
    # ------------------------------------------------------------
    for x in sorted(K_points):
        Cand = prev_U[int(x)] if prev_U[int(x)] else frozenset([int(x)])
        _add_if_needed(x, Cand)

    # ------------------------------------------------------------
    # Stage 2: process L (minimal dead regions), add only if needed
    # ------------------------------------------------------------
    for x in sorted(L_points):
        S = region_of_new.get(int(x), _fs(Vprev_of_point(x)))
        Cand = B_dead_region(S)  # just S
        _add_if_needed(x, Cand)

    # ------------------------------------------------------------
    # Stage 3: process M (remaining uncovered), add only if needed
    # ------------------------------------------------------------
    B_for_region = {}

    for x in M_points:
        x = int(x)

        S = _fs(Vprev_of_point(x))
        if not S:
            Cand = prev_U[x] if prev_U[x] else frozenset([x])
            _add_if_needed(x, Cand)
            continue

        if S in B_for_region:
            _add_if_needed(x, B_for_region[S])
            continue

        if len(S & X_i) == 0:
            Cand = B_dead_region(S)  
            B_for_region[S] = Cand
            _add_if_needed(x, Cand)
            continue

        core = _fs(S & X_i)
        Vmin = _minimal_containing(curr_VA, core)  

        union_dead_B = frozenset()
        for T in dead_minimals_sorted:
            if T.issubset(S):
                union_dead_B |= B_dead_region(T)
                
        Cand = _fs(S | Vmin | union_dead_B)
        B_for_region[S] = Cand
        _add_if_needed(x, Cand)


    # ------------------------------------------------------------
    # NO-CARVING SATURATION (post-processing)
    # ------------------------------------------------------------
    B_sets = saturate_B_sets_against_curr_VA(B_sets, curr_VA, X_i)

    
    # ------------------------------------------------------------
    # Dedup B_sets (by equality) while preserving order
    # ------------------------------------------------------------
    seen = set()
    B_sets_dedup = []
    for S in B_sets:
        S = _fs(S)
        if S and S not in seen:
            seen.add(S)
            B_sets_dedup.append(S)
    B_sets = B_sets_dedup

    # ------------------------------------------------------------
    # Build B_by_point for compatibility (derived from B_sets)
    # ------------------------------------------------------------
    B_by_point = {}
    for x in sorted(uncovered_now):
        chosen = None
        for B in B_sets:
            if x in B:
                chosen = B
                break
        if chosen is None:
            chosen = prev_U[x] if prev_U[x] else frozenset([x])
        B_by_point[int(x)] = _fs(chosen)

    return B_by_point, B_sets
