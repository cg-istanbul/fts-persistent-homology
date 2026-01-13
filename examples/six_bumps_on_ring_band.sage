import os
import random
from sage.all import RR, matrix, sqrt, cos, sin, pi, set_random_seed, Graphics, point

load('fts_pipeline.sage')

# ----------------------------
# Parameters (data)
# ----------------------------
SEED = 0

RING_RADIUS   = RR(3.0)     
RADIAL_STD    = RR(0.18)    
ANGULAR_JITTER_BG = RR(0.0) 

N_BACKGROUND  = 72          

NUM_BUMPS     = 6
PTS_PER_BUMP  = 20         

BUMP_ANG_STD  = RR(0.10)    
BUMP_RAD_STD  = RR(0.12)    

# ----------------------------
# Parameters (filtration)
# ----------------------------
m = 25
k_density = 6
lambda_scale = 2.0

schedule = "quantile"
tau_q_start = 0.10
tau_q_end   = 0.98
tau_min = None
tau_max = None

anchor_cap_topA = None
anchor_cap_until_step = None

include_X = False
tiny_k = 1

uncovered_k_every = 3
uncovered_k_max   = 8

do_saturation = True
saturation_passes = 1

coarsen_mode  = "none"     
coarsen_start = None
radius_start  = 0.0
radius_end    = None

merge_maximals = False
merge_every = 1
merge_k = 1
eps_q_start = 0.05
eps_q_end   = 0.90
keep_generators = True

check_continuity = True

# ----------------------------
# Geometry helpers
# ----------------------------
def euclidean_distance(p, q):
    return sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

def distance_matrix_from_points(points):
    n = len(points)
    D = matrix(RR, n, n)
    for i in range(n):
        D[i,i] = RR(0)
        for j in range(i+1, n):
            d = euclidean_distance(points[i], points[j])
            D[i,j] = d
            D[j,i] = d
    return D

def polar_to_cart(r, ang):
    return (r*cos(ang), r*sin(ang))

def ring_band_with_bumps(seed=0):
    set_random_seed(int(seed))
    random.seed(int(seed))

    pts = []

    for t in range(N_BACKGROUND):
        ang = RR(2*pi*t/N_BACKGROUND)
        if ANGULAR_JITTER_BG != 0:
            ang += RR(random.gauss(0.0, float(ANGULAR_JITTER_BG)))
        r = RR(random.gauss(float(RING_RADIUS), float(RADIAL_STD)))
        pts.append(polar_to_cart(r, ang))

    centers = [RR(2*pi*j/NUM_BUMPS) for j in range(NUM_BUMPS)]
    for a0 in centers:
        for _ in range(PTS_PER_BUMP):
            ang = RR(random.gauss(float(a0), float(BUMP_ANG_STD)))
            r   = RR(random.gauss(float(RING_RADIUS), float(BUMP_RAD_STD)))
            pts.append(polar_to_cart(r, ang))

    return pts

def save_scatter_plot(points, out_png):
    G = Graphics()
    G += point(points, size=18)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    G.save(out_png)
    print("[plot] saved to: {}".format(out_png))

# ----------------------------
# Main
# ----------------------------
points = ring_band_with_bumps(seed=SEED)
print("Total points n = {} (background = {}, bumps = {})".format(
    len(points), N_BACKGROUND, NUM_BUMPS*PTS_PER_BUMP
))

out_png = os.path.join(os.path.dirname(__file__), "figures", "six_bumps_on_ring_band.png")
save_scatter_plot(points, out_png)

D = distance_matrix_from_points(points)

run_density_filtration(
    D,
    m=m,
    k_density=k_density,
    lambda_scale=lambda_scale,
    schedule=schedule,
    tau_q_start=tau_q_start,
    tau_q_end=tau_q_end,
    tau_min=tau_min,
    tau_max=tau_max,
    anchor_cap_topA=anchor_cap_topA,
    anchor_cap_until_step=anchor_cap_until_step,
    include_X=include_X,
    tiny_k=tiny_k,
    uncovered_k_every=uncovered_k_every,
    uncovered_k_max=uncovered_k_max,
    do_saturation=do_saturation,
    saturation_passes=saturation_passes,
    coarsen_mode=coarsen_mode,
    coarsen_start=coarsen_start,
    radius_start=radius_start,
    radius_end=radius_end,
    merge_maximals=merge_maximals,
    merge_every=merge_every,
    merge_k=merge_k,
    eps_q_start=eps_q_start,
    eps_q_end=eps_q_end,
    keep_generators=keep_generators,
    check_continuity=check_continuity,
)
