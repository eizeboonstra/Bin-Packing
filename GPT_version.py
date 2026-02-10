# ae4446_2dbpp_gurobi.py
# Exact MILP for AE4446 vertical 2D bin packing with stability, fragility, incompatibility, LD3 cut
#
# Outputs (as required by assignment):
#   bins_used (list)
#   Items_in_Bin (dict: bin -> list(items))
#   I_info_solution (dict: item -> [x_i, z_i, w_i, h_i])
#
# Usage:
#   python ae4446_2dbpp_gurobi.py
#
# Requires:
#   pip install gurobipy matplotlib

import math
import pickle
from dataclasses import dataclass

import gurobipy as gp
from gurobipy import GRB

import matplotlib.pyplot as plt
import matplotlib.patches as patches


# -------------------------
# Data containers
# -------------------------
@dataclass
class Bin:
    ID: int
    btype: int
    W: float
    H: float
    cost: float
    a: float  # cut x-intercept (only if type 1, else -1)
    b: float  # cut z-intercept (only if type 1, else -1)

@dataclass
class Item:
    ID: int
    L: float
    H: float
    can_rotate: int
    fragile: int
    perishable: int
    radioactive: int


# -------------------------
# Helpers
# -------------------------
def load_inputs(b_path="B.pickle", i_path="I.pickle"):
    with open(b_path, "rb") as handle:
        B = pickle.load(handle)
    with open(i_path, "rb") as handle:
        I = pickle.load(handle)

    bins = []
    count= 0
    for k, v in B.items():
        if count < 2:
            count += 1
            btype = int(v[0])
            W, H, _count, cost, a, b = v[1]
            bins.append(Bin(ID=int(k), btype=btype, W=float(W), H=float(H), cost=float(cost), a=float(a), b=float(b)))

    items = []
    for k, v in I.items():
        if count < 20:
            count += 1
            L, H, can_rot, frag, per, rad = v
            items.append(Item(ID=int(k), L=float(L), H=float(H),
                            can_rotate=int(can_rot), fragile=int(frag),
                            perishable=int(per), radioactive=int(rad)))

    bins.sort(key=lambda bb: bb.ID)
    items.sort(key=lambda ii: ii.ID)
    return bins, items


def big_m(bins, items):
    # Safe-ish big-M: max bin dimension + max item dimension
    maxW = max(b.W for b in bins)
    maxH = max(b.H for b in bins)
    maxL = max(it.L for it in items)
    maxh = max(it.H for it in items)
    return maxW + maxL + maxH + maxh


# -------------------------
# Model build
# -------------------------
def build_model(bins, items, timelimit_sec=2 * 60 * 60, log_name="ae4446_exact.log"):
    m = gp.Model("AE4446_2DBPP_Exact")
    m.Params.TimeLimit = timelimit_sec
    m.Params.LogFile = log_name  # required to report incumbent/bound/gap/time :contentReference[oaicite:2]{index=2}
    m.Params.MIPFocus = 1
    m.Params.Presolve = 2

    B = [b.ID for b in bins]
    I = [it.ID for it in items]
    bin_by_id = {b.ID: b for b in bins}
    item_by_id = {it.ID: it for it in items}

    M = big_m(bins, items)
    eps = 1e-4  # small tolerance for strict-ish geometry

    # --------
    # Vars
    # --------
    # bin used
    u = m.addVars(B, vtype=GRB.BINARY, name="u")

    # assignment
    p = m.addVars(I, B, vtype=GRB.BINARY, name="p")

    # position (global coords, constrained by chosen bin via big-M)
    x = m.addVars(I, lb=0.0, vtype=GRB.CONTINUOUS, name="x")
    z = m.addVars(I, lb=0.0, vtype=GRB.CONTINUOUS, name="z")

    # rotation
    rho = m.addVars(I, vtype=GRB.BINARY, name="rho")

    # effective width/height
    w = m.addVars(I, lb=0.0, vtype=GRB.CONTINUOUS, name="w")
    h = m.addVars(I, lb=0.0, vtype=GRB.CONTINUOUS, name="h")

    # Non-overlap binaries per (i<j, b)
    # left[i,j,b]=1 means i left of j in bin b; below[i,j,b]=1 means i below j in bin b
    left = {}
    below = {}
    for ii in range(len(I)):
        for jj in range(ii + 1, len(I)):
            i = I[ii]
            j = I[jj]
            for b in B:
                left[i, j, b] = m.addVar(vtype=GRB.BINARY, name=f"L_{i}_{j}_{b}")
                left[j, i, b] = m.addVar(vtype=GRB.BINARY, name=f"L_{j}_{i}_{b}")
                below[i, j, b] = m.addVar(vtype=GRB.BINARY, name=f"B_{i}_{j}_{b}")
                below[j, i, b] = m.addVar(vtype=GRB.BINARY, name=f"B_{j}_{i}_{b}")

    # Stability vars:
    # g_i: on ground
    g = m.addVars(I, vtype=GRB.BINARY, name="g")

    # beta1_{i,j} supports left bottom corner of i by j
    # beta2_{i,j} supports right bottom corner of i by j
    beta1 = m.addVars(I, I, vtype=GRB.BINARY, name="beta1")
    beta2 = m.addVars(I, I, vtype=GRB.BINARY, name="beta2")

    # gamma_i: one corner supported by LD3 cut (only meaningful if item in a cut bin) 
    gamma = m.addVars(I, vtype=GRB.BINARY, name="gamma")

    # Incompatibility bin-level indicators
    Pbin = m.addVars(B, vtype=GRB.BINARY, name="Pbin")
    Rbin = m.addVars(B, vtype=GRB.BINARY, name="Rbin")

    # LD3 cut disjunction helper binaries per (i,b) for cut bins:
    # d1: x_i >= a, d2: z_i >= b, d3: x/a + z/b >= 1
    d1 = {}
    d2 = {}
    d3 = {}
    for b in B:
        bb = bin_by_id[b]
        if bb.a > 0 and bb.b > 0:  # cut present
            for i in I:
                d1[i, b] = m.addVar(vtype=GRB.BINARY, name=f"cut_d1_{i}_{b}")
                d2[i, b] = m.addVar(vtype=GRB.BINARY, name=f"cut_d2_{i}_{b}")
                d3[i, b] = m.addVar(vtype=GRB.BINARY, name=f"cut_d3_{i}_{b}")

    m.update()

    # --------
    # Objective: minimize cost of used ULDs :contentReference[oaicite:4]{index=4}
    # --------
    m.setObjective(gp.quicksum(bin_by_id[b].cost * u[b] for b in B), GRB.MINIMIZE)

    # --------
    # Core constraints
    # --------
    # Each item assigned to exactly one bin
    for i in I:
        m.addConstr(gp.quicksum(p[i, b] for b in B) == 1, name=f"assign_{i}")

    # Link bin usage
    for i in I:
        for b in B:
            m.addConstr(p[i, b] <= u[b], name=f"link_{i}_{b}")

    # Rotation feasibility + define w,h
    for i in I:
        it = item_by_id[i]
        if it.can_rotate == 0:
            m.addConstr(rho[i] == 0, name=f"no_rot_{i}")
        # w_i = L*(1-rho) + H*rho ; h_i = H*(1-rho) + L*rho
        m.addConstr(w[i] == it.L * (1 - rho[i]) + it.H * rho[i], name=f"def_w_{i}")
        m.addConstr(h[i] == it.H * (1 - rho[i]) + it.L * rho[i], name=f"def_h_{i}")

    # Fit inside chosen bin (big-M)
    for i in I:
        for b in B:
            bb = bin_by_id[b]
            m.addConstr(x[i] + w[i] <= bb.W + M * (1 - p[i, b]), name=f"inW_{i}_{b}")
            m.addConstr(z[i] + h[i] <= bb.H + M * (1 - p[i, b]), name=f"inH_{i}_{b}")

    # Non-overlap (only if both in same bin)
    for ii in range(len(I)):
        for jj in range(ii + 1, len(I)):
            i = I[ii]
            j = I[jj]
            for b in B:
                # If both in b, at least one relative placement must hold
                m.addConstr(left[i, j, b] + left[j, i, b] + below[i, j, b] + below[j, i, b]
                            >= p[i, b] + p[j, b] - 1,
                            name=f"disj_{i}_{j}_{b}")

                # Big-M implications
                m.addConstr(x[i] + w[i] <= x[j] + M * (1 - left[i, j, b]),
                            name=f"left_{i}_{j}_{b}")
                m.addConstr(x[j] + w[j] <= x[i] + M * (1 - left[j, i, b]),
                            name=f"left_{j}_{i}_{b}")
                m.addConstr(z[i] + h[i] <= z[j] + M * (1 - below[i, j, b]),
                            name=f"below_{i}_{j}_{b}")
                m.addConstr(z[j] + h[j] <= z[i] + M * (1 - below[j, i, b]),
                            name=f"below_{j}_{i}_{b}")

    # Incompatibility: perishable and radioactive cannot be in same ULD :contentReference[oaicite:5]{index=5}
    for b in B:
        # Link to assignments
        for i in I:
            it = item_by_id[i]
            if it.perishable == 1:
                m.addConstr(Pbin[b] >= p[i, b], name=f"linkP_{i}_{b}")
            if it.radioactive == 1:
                m.addConstr(Rbin[b] >= p[i, b], name=f"linkR_{i}_{b}")
        m.addConstr(Pbin[b] + Rbin[b] <= 1, name=f"incompat_{b}")

    # Fragility: no box stacked on top of fragile ones 
    # Enforce: fragile item j cannot support any other item's corners
    for j in I:
        if item_by_id[j].fragile == 1:
            for i in I:
                if i != j:
                    m.addConstr(beta1[i, j] == 0, name=f"frag_no_beta1_{i}_{j}")
                    m.addConstr(beta2[i, j] == 0, name=f"frag_no_beta2_{i}_{j}")

    # Stability constraints:
    # - item cannot float: both bottom corners supported by ground, other items, or (for LD3 cut) the cut. 
    # Given in assignment (with gamma_i for cut support): gamma_i + sum_j beta1_{i,j} + sum_j beta2_{i,j} + 2 g_i >= 2
    for i in I:
        m.addConstr(
            gamma[i] + gp.quicksum(beta1[i, j] for j in I if j != i) + gp.quicksum(beta2[i, j] for j in I if j != i) + 2 * g[i] >= 2,
            name=f"stability_{i}"
        )

        # On ground => z_i == 0 (if g_i=1)
        m.addConstr(z[i] <= M * (1 - g[i]), name=f"ground_z_ub_{i}")
        # (z >= 0 already)

        # A corner-support implies "i is above j": z_i = z_j + h_j
        for j in I:
            if j == i:
                m.addConstr(beta1[i, j] == 0, name=f"no_self_beta1_{i}")
                m.addConstr(beta2[i, j] == 0, name=f"no_self_beta2_{i}")
                continue

            # If beta is 1, enforce vertical contact
            m.addConstr(z[i] - (z[j] + h[j]) <= M * (1 - beta1[i, j]), name=f"b1_vert_ub_{i}_{j}")
            m.addConstr(z[i] - (z[j] + h[j]) >= -M * (1 - beta1[i, j]), name=f"b1_vert_lb_{i}_{j}")
            m.addConstr(z[i] - (z[j] + h[j]) <= M * (1 - beta2[i, j]), name=f"b2_vert_ub_{i}_{j}")
            m.addConstr(z[i] - (z[j] + h[j]) >= -M * (1 - beta2[i, j]), name=f"b2_vert_lb_{i}_{j}")
            for b in B:
                    # If beta is 1, enforce same-bin
                    m.addConstr(p[i, b] + p[j, b] >= 2 * beta1[i, j], name=f"b1_samebin_{i}_{j}_{b}")
                    m.addConstr(p[i, b] + p[j, b] >= 2 * beta2[i, j], name=f"b2_samebin_{i}_{j}_{b}")

            # Left corner x_i must lie within [x_j, x_j + w_j]
            m.addConstr(x[i] >= x[j] - M * (1 - beta1[i, j]), name=f"b1_xmin_{i}_{j}")
            m.addConstr(x[i] <= x[j] + w[j] + M * (1 - beta1[i, j]), name=f"b1_xmax_{i}_{j}")

            # Right corner x_i + w_i must lie within [x_j, x_j + w_j]
            m.addConstr(x[i] + w[i] >= x[j] - M * (1 - beta2[i, j]), name=f"b2_xmin_{i}_{j}")
            m.addConstr(x[i] + w[i] <= x[j] + w[j] + M * (1 - beta2[i, j]), name=f"b2_xmax_{i}_{j}")

            # Prevent impossible support if i and j are not in same bin:
            # beta <= sum_b p[i,b] AND beta <= sum_b p[j,b] AND beta <= 1 if same bin exists.
            # Stronger: beta <= sum_b y_{ijb} where y_{ijb} <= p[i,b], y_{ijb} <= p[j,b]
            # Keep it simpler but correct-enough:
            m.addConstr(beta1[i, j] <= gp.quicksum(p[i, b] for b in B), name=f"b1_assign_i_{i}_{j}")
            m.addConstr(beta1[i, j] <= gp.quicksum(p[j, b] for b in B), name=f"b1_assign_j_{i}_{j}")
            m.addConstr(beta2[i, j] <= gp.quicksum(p[i, b] for b in B), name=f"b2_assign_i_{i}_{j}")
            m.addConstr(beta2[i, j] <= gp.quicksum(p[j, b] for b in B), name=f"b2_assign_j_{i}_{j}")

    # LD3 cut constraints:
    # - For bins with cut (a,b): removed triangle with vertices (0,0), (a,0), (0,b).
    #   Feasible region is: x>=a OR z>=b OR x/a + z/b >= 1. Implemented as 3-way OR with binaries. 
    # - gamma_i can only be 1 if item is assigned to at least one cut bin, and then forces the lower-left corner to lie on the cut line.
    for i in I:
        # default: gamma must be 0 if no cut bins exist
        cut_bins = [b.ID for b in bins if b.a > 0 and b.b > 0]
        if not cut_bins:
            m.addConstr(gamma[i] == 0, name=f"no_cut_gamma_{i}")
        else:
            m.addConstr(gamma[i] <= gp.quicksum(p[i, b] for b in cut_bins), name=f"gamma_only_if_cutbin_{i}")

        # Also: if item is in a non-cut bin, gamma has no meaning; it doesn't help there.
        for b in B:
            bb = bin_by_id[b]
            if bb.a <= 0 or bb.b <= 0:
                # if assigned to non-cut bin, gamma must be 0 (otherwise it could “cheat” stability)
                m.addConstr(gamma[i] <= 1 - p[i, b], name=f"gamma_not_in_nocut_{i}_{b}")

    for b in B:
        bb = bin_by_id[b]
        if bb.a > 0 and bb.b > 0:
            a = bb.a
            c = bb.b

            for i in I:
                # activate only if item assigned to this bin
                m.addConstr(d1[i, b] + d2[i, b] + d3[i, b] >= p[i, b], name=f"cut_or_{i}_{b}")

                # d1 => x >= a
                m.addConstr(x[i] >= a - M * (1 - d1[i, b]), name=f"cut_d1_{i}_{b}")

                # d2 => z >= b
                m.addConstr(z[i] >= c - M * (1 - d2[i, b]), name=f"cut_d2_{i}_{b}")

                # d3 => x/a + z/b >= 1
                m.addConstr((x[i] / a) + (z[i] / c) >= 1 - M * (1 - d3[i, b]), name=f"cut_d3_{i}_{b}")

                # Make d's irrelevant unless assigned:
                m.addConstr(d1[i, b] <= p[i, b], name=f"cut_d1_link_{i}_{b}")
                m.addConstr(d2[i, b] <= p[i, b], name=f"cut_d2_link_{i}_{b}")
                m.addConstr(d3[i, b] <= p[i, b], name=f"cut_d3_link_{i}_{b}")

                # If gamma_i=1 and item in this cut bin -> force corner on cut line within [0,a]x[0,b]
                # x/a + z/b == 1
                m.addConstr(x[i] <= a + M * (1 - gamma[i]), name=f"gamma_xub_{i}_{b}")
                m.addConstr(z[i] <= c + M * (1 - gamma[i]), name=f"gamma_zub_{i}_{b}")
                m.addConstr((x[i] / a) + (z[i] / c) >= 1 - M * (1 - gamma[i]), name=f"gamma_line_lb_{i}_{b}")
                m.addConstr((x[i] / a) + (z[i] / c) <= 1 + M * (1 - gamma[i]), name=f"gamma_line_ub_{i}_{b}")

        else:
            # No cut => disallow d vars implicitly by not creating them; nothing to do.
            pass

    # Symmetry-breaking (lightweight): prefer using lower-index bins first (optional but helps)
    for b1, b2 in zip(B[:-1], B[1:]):
        m.addConstr(u[b1] >= u[b2], name=f"sym_u_{b1}_{b2}")

    m.update()
    return m, (u, p, x, z, w, h, rho, g, beta1, beta2, gamma)


# -------------------------
# Solve + extract outputs
# -------------------------
def extract_solution(bins, items, vars_pack):
    (u, p, x, z, w, h, rho, g, beta1, beta2, gamma) = vars_pack
    B = [b.ID for b in bins]
    I = [it.ID for it in items]

    bins_used = [b for b in B if u[b].X > 0.5]

    Items_in_Bin = {}
    for b in bins_used:
        Items_in_Bin[b] = [i for i in I if p[i, b].X > 0.5]

    I_info_solution = {}
    for i in I:
        I_info_solution[i] = [float(x[i].X), float(z[i].X), float(w[i].X), float(h[i].X)]

    return bins_used, Items_in_Bin, I_info_solution


def save_pickles(bins_used, Items_in_Bin, I_info_solution,
                out_bins="bins_used.pickle",
                out_items="Items_in_Bin.pickle",
                out_info="I_info_solution.pickle"):
    with open(out_bins, "wb") as f:
        pickle.dump(bins_used, f)
    with open(out_items, "wb") as f:
        pickle.dump(Items_in_Bin, f)
    with open(out_info, "wb") as f:
        pickle.dump(I_info_solution, f)


def plot_solution(bins, items, Items_in_Bin, I_info_solution, out_png="packing.png"):
    item_by_id = {it.ID: it for it in items}
    bin_by_id = {b.ID: b for b in bins}

    used_bins = sorted(list(Items_in_Bin.keys()))
    n = len(used_bins)
    if n == 0:
        print("No bins used? (unexpected)")
        return

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), squeeze=False)

    for idx, b in enumerate(used_bins):
        ax = axes[0][idx]
        bb = bin_by_id[b]

        ax.set_title(f"Bin {b} (type {bb.btype})")
        ax.set_xlim(0, bb.W)
        ax.set_ylim(0, bb.H)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x")
        ax.set_ylabel("z")

        # Draw bin boundary
        ax.add_patch(patches.Rectangle((0, 0), bb.W, bb.H, fill=False, linewidth=2))

        # Draw cut (if any)
        if bb.a > 0 and bb.b > 0:
            ax.plot([0, bb.a], [bb.b, 0], linewidth=2)
            # Optional: shade cut triangle
            cut_poly = patches.Polygon([(0, 0), (bb.a, 0), (0, bb.b)], closed=True, alpha=0.15)
            ax.add_patch(cut_poly)

        # Draw items
        for i in Items_in_Bin[b]:
            xi, zi, wi, hi = I_info_solution[i]
            fragile = item_by_id[i].fragile == 1

            rect = patches.Rectangle((xi, zi), wi, hi,
                                     fill=True, alpha=0.35 if fragile else 0.25,
                                     linewidth=1.5)
            ax.add_patch(rect)

            # label: id,(frag,per,rad) like example :contentReference[oaicite:9]{index=9}
            it = item_by_id[i]
            label = f"{i},({it.fragile},{it.perishable},{it.radioactive})"
            ax.text(xi + wi/2, zi + hi/2, label, ha="center", va="center", fontsize=8)

        ax.grid(True, linewidth=0.4, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    bins, items = load_inputs("B.pickle", "I.pickle")  # assignment input format 
    m, vars_pack = build_model(bins, items, timelimit_sec=2*60*60, log_name="ae4446_exact.log")  # 2 hours :contentReference[oaicite:11]{index=11}
    m.optimize()

    if m.SolCount == 0:
        print("No feasible solution found within time limit.")
        return

    bins_used, Items_in_Bin, I_info_solution = extract_solution(bins, items, vars_pack)
    print("bins_used:", bins_used)
    print("Items_in_Bin:", Items_in_Bin)

    save_pickles(bins_used, Items_in_Bin, I_info_solution,
                out_bins="bins_used.pickle",
                out_items="Items_in_Bin.pickle",
                out_info="I_info_solution.pickle")

    plot_solution(bins, items, Items_in_Bin, I_info_solution, out_png="packing.png")
    print("Wrote: ae4446_exact.log, bins_used.pickle, Items_in_Bin.pickle, I_info_solution.pickle, packing.png")


if __name__ == "__main__":
    main()
