# scripts/topology/schedule.sage

from sage.all import RR, floor

# ----------------------------
# Quantile utilities
# ----------------------------

def _quantile_value(sorted_vals, q):
    n = len(sorted_vals)
    if n == 0:
        return RR(0)
    if q <= 0:
        return RR(sorted_vals[0])
    if q >= 1:
        return RR(sorted_vals[-1])
    idx = int(floor(RR(q) * RR(n - 1)))
    idx = max(0, min(idx, n - 1))
    return RR(sorted_vals[idx])


def taus_from_rho_quantiles(rho, m, q_start=0.02, q_end=0.90, mode="linear_in_q", q_exponent=1.0):
    assert int(m) >= 2
    q_start = float(q_start)
    q_end = float(q_end)
    q_start = max(0.0, min(q_start, 1.0))
    q_end = max(0.0, min(q_end, 1.0))
    if q_end < q_start:
        q_start, q_end = q_end, q_start

    vals = sorted([RR(x) for x in rho])
    t_start = _quantile_value(vals, q_start)
    t_end = _quantile_value(vals, q_end)

    taus = []
    if str(mode) == "linear_in_tau":
        if t_end == t_start:
            return [t_start for _ in range(int(m))]
        for i in range(int(m)):
            a = RR(i) / RR(int(m) - 1)
            taus.append(t_start + a * (t_end - t_start))
    else:
        exp = float(q_exponent) if q_exponent is not None else 1.0
        if exp <= 0.0:
            exp = 1.0
        for i in range(int(m)):
            a = float(i) / float(int(m) - 1)
            q = q_start + (a ** exp) * (q_end - q_start)
            taus.append(_quantile_value(vals, q))

    for i in range(1, int(m)):
        if taus[i] < taus[i - 1]:
            taus[i] = taus[i - 1]

    uniq = []
    for t in taus:
        if len(uniq) == 0 or t != uniq[-1]:
            uniq.append(t)

    if len(uniq) >= int(m):
        return uniq[: int(m)]

    return taus_linear(t_start, t_end, int(m))


def taus_linear(tau_min, tau_max, m):
    assert int(m) >= 2
    t0 = RR(tau_min)
    t1 = RR(tau_max)
    if t1 < t0:
        t0, t1 = t1, t0
    if t1 == t0:
        return [t0 for _ in range(int(m))]
    taus = []
    for i in range(int(m)):
        a = RR(i) / RR(int(m) - 1)
        taus.append(t0 + a * (t1 - t0))
    return taus
