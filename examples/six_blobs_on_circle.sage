import os
import random
from sage.all import RR, matrix, sqrt, cos, sin, pi, set_random_seed, Graphics, point

load('fts_pipeline.sage')

# ----------------------------
# Parameters (data)
# ----------------------------
SEED = 0
NUM_BLOBS = 6
PTS_PER_BLOB = 20          
CIRCLE_RADIUS = RR(3.0)    
BLOB_STD = RR(0.20)

# ----------------------------
# Parameters (filtration)
# ----------------------------
m = 25
k_density = 6
lambda_scale = 2.0

schedule = "quantile"
tau_q_start = 0.10
tau_q_end = 0.98
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

def six_blobs_on_circle(seed=0,
                        num_blobs=6,
                        pts_per_blob=20,
                        circle_radius=RR(3.0),
                        blob_std=RR(0.20)):
    set_random_seed(int(seed))
    random.seed(int(seed))

    centers = []
    for t in range(num_blobs):
        ang = RR(2*pi*t/num_blobs)
        centers.append((circle_radius*cos(ang), circle_radius*sin(ang)))

    pts = []
    for (cx, cy) in centers:
        for _ in range(pts_per_blob):
            x = RR(random.gauss(float(cx), float(blob_std)))
            y = RR(random.gauss(float(cy), float(blob_std)))
            pts.append((x, y))
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
points = six_blobs_on_circle(
    seed=SEED,
    num_blobs=NUM_BLOBS,
    pts_per_blob=PTS_PER_BLOB,
    circle_radius=CIRCLE_RADIUS,
    blob_std=BLOB_STD
)

out_png = os.path.join(os.path.dirname(__file__), "figures", "six_blobs_on_circle.png")
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
