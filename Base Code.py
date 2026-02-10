import gurobipy as gp
import pickle

import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Bin:
    def __init__(self, ID, length, height, cost, a, b):
        self.ID = ID
        self.length = length
        self.height = height
        self.cost = cost
        self.a = a
        self.b = b

class Item:
    def __init__(self, ID, length, height, rotate, fragile, perishable, radioactive):
        self.ID = ID
        self.length = length
        self.height = height
        self.rotate = rotate
        self.fragile = fragile
        self.perishable = perishable
        self.radioactive = radioactive

with open('B.pickle', 'rb') as handle:
    B = pickle.load(handle)
with open('I.pickle', 'rb') as handle:
    I = pickle.load(handle)

B_values = list(B.values())
Bins = []
for b in range(len(B_values)):
        Bins.append(Bin(b, B_values[b][1][0], B_values[b][1][1], B_values[b][1][3], B_values[b][1][4], B_values[b][1][5]))
Lmax = max([bin.length for bin in Bins])
Hmax = max([bin.height for bin in Bins])

I_values = list(I.values())
Items = []
for i in range(len(I_values)):
        Items.append(Item(i, I_values[i][0], I_values[i][1], I_values[i][2], I_values[i][3], I_values[i][4], I_values[i][5]))


model = gp.Model("Bin Packing Problem")


z = {}
for b in Bins:
    z[b.ID] = model.addVar(obj= b.cost, vtype=gp.GRB.BINARY, name="z_{}".format(b.ID))



p = {}
for i in Items:
    for bin in Bins:
        p[i.ID, bin.ID] = model.addVar(vtype=gp.GRB.BINARY, name="p_{}_{}".format(i.ID, bin.ID))

g = {}
rho = {}
gamma = {}
for i in Items:
    rho[i.ID] = model.addVar(vtype=gp.GRB.BINARY, name="rho_{}".format(i.ID))
    g[i.ID] = model.addVar(vtype=gp.GRB.BINARY, name="g_{}".format(i.ID))
    gamma[i.ID] = model.addVar(vtype=gp.GRB.BINARY, name="gamma_{}".format(i.ID))

radioactive = {}
for bin in Bins:
    radioactive[bin.ID] = model.addVar(vtype=gp.GRB.BINARY, name="radioactive_{}".format(bin.ID))

B1 = {}
B2 = {}
for i in Items:
    for j in Items:
        if i.ID != j.ID:
            B1[i.ID, j.ID] = model.addVar(vtype=gp.GRB.BINARY, name="B1_{}_{}".format(i.ID, j.ID))
            B2[i.ID, j.ID] = model.addVar(vtype=gp.GRB.BINARY, name="B2_{}_{}".format(i.ID, j.ID))

x = {}
y = {}
for i in Items:
    x[i.ID] = model.addVar(lb=0, vtype=gp.GRB.CONTINUOUS, name="x_{}".format(i.ID))
    y[i.ID] = model.addVar(lb=0, vtype=gp.GRB.CONTINUOUS, name="y_{}".format(i.ID))


def W(item_id):
    it = Items[item_id]
    return it.length * (1 - rho[item_id]) + it.height * rho[item_id]

def H(item_id):
    it = Items[item_id]
    return it.height * (1 - rho[item_id]) + it.length * rho[item_id]


# --- ONLY PAIRS i<j ---
pair_list = [(i.ID, j.ID) for idx_i, i in enumerate(Items) for j in Items[idx_i+1:]]

# w[i,j,k] = 1 iff i and j are both assigned to bin k
w = {}
left = {}
right = {}
below = {}
above = {}

for (i, j) in pair_list:
    for bin in Bins:
        k = bin.ID
        w[i, j, k] = model.addVar(vtype=gp.GRB.BINARY, name=f"w_{i}_{j}_{k}")
        left[i, j, k]  = model.addVar(vtype=gp.GRB.BINARY, name=f"L_{i}_{j}_{k}")
        right[i, j, k] = model.addVar(vtype=gp.GRB.BINARY, name=f"R_{i}_{j}_{k}")
        below[i, j, k] = model.addVar(vtype=gp.GRB.BINARY, name=f"B_{i}_{j}_{k}")
        above[i, j, k] = model.addVar(vtype=gp.GRB.BINARY, name=f"A_{i}_{j}_{k}")

        # linearize w = AND(p[i,k], p[j,k])
        model.addConstr(w[i, j, k] <= p[i, k], name=f"wub1_{i}_{j}_{k}")
        model.addConstr(w[i, j, k] <= p[j, k], name=f"wub2_{i}_{j}_{k}")
        model.addConstr(w[i, j, k] >= p[i, k] + p[j, k] - 1, name=f"wlb_{i}_{j}_{k}")

        # exactly one separation mode if in same bin; none if not
        model.addConstr(left[i, j, k] + right[i, j, k] + below[i, j, k] + above[i, j, k] == w[i, j, k],
                        name=f"disj_{i}_{j}_{k}")

        # INDICATORS (no big-M)
        # if left=1 then i is left of j: x_i + W_i <= x_j
        model.addGenConstrIndicator(left[i, j, k], 1, x[i] + W(i) <= x[j], name=f"indL_{i}_{j}_{k}")

        # if right=1 then j is left of i: x_j + W_j <= x_i
        model.addGenConstrIndicator(right[i, j, k], 1, x[j] + W(j) <= x[i], name=f"indR_{i}_{j}_{k}")

        # if below=1 then i is below j: y_i + H_i <= y_j
        model.addGenConstrIndicator(below[i, j, k], 1, y[i] + H(i) <= y[j], name=f"indB_{i}_{j}_{k}")

        # if above=1 then j is below i: y_j + H_j <= y_i
        model.addGenConstrIndicator(above[i, j, k], 1, y[j] + H(j) <= y[i], name=f"indA_{i}_{j}_{k}")
xbound = {}
ybound = {}
assigned = {}
for i in Items:
    xbound[i.ID] = model.addConstr(x[i.ID] + i.length * (1 - rho[i.ID]) + i.height * rho[i.ID] <= gp.quicksum(bin.length * p[i.ID, bin.ID] for bin in Bins), name="xbound_{}".format(i.ID))
    ybound[i.ID] = model.addConstr(y[i.ID] + i.height * (1 - rho[i.ID]) + i.length * rho[i.ID] <= gp.quicksum(bin.height * p[i.ID, bin.ID] for bin in Bins), name="ybound_{}".format(i.ID))
    assigned[i.ID] = model.addConstr(gp.quicksum(p[i.ID, bin.ID] for bin in Bins) == 1, name="assigned_{}".format(i.ID))

usage = {}
for bin in Bins:
    usage[bin.ID] = model.addConstr(gp.quicksum(p[i.ID, bin.ID] for i in Items) <= 1000 * z[bin.ID], name="bin_usage_{}".format(bin.ID))
model.addConstr(z[0]<=z[1])
model.addConstr(z[2]<=z[3])

grounding = {}
for i in Items:
    grounding[i.ID] = model.addConstr(y[i.ID] <= 10000 * (1-g[i.ID]), name="grounding_{}".format(i.ID))

cornerl = {}
corneru = {}
gamma_con = {}
for i in Items:
    for bin in Bins:
        if bin.b>0:
            cornerl[i.ID, bin.ID] = model.addConstr(y[i.ID]>=-bin.b/bin.a * x[i.ID] + bin.b - (1-p[i.ID, bin.ID])*(1000000), name="corner_{}_{}".format(i.ID, bin.ID))
            corneru[i.ID, bin.ID] = model.addConstr(y[i.ID]<=-bin.b/bin.a * x[i.ID] + bin.b + (2-p[i.ID, bin.ID] - gamma[i.ID])*(1000000), name="corner_{}_{}".format(i.ID, bin.ID))
        else:
            gamma_con[i.ID, bin.ID] = model.addConstr(gamma[i.ID] <= 1-p[i.ID, bin.ID], name="gamma_con_{}_{}".format(i.ID, bin.ID))


bcon = {}
rcon = {}
for i in Items:
    for j in Items:
        for bin in Bins:
            if i.ID != j.ID:
                bcon[i.ID, j.ID, bin.ID] = model.addConstr(y[i.ID] - (y[j.ID] + j.height * (1 - rho[j.ID]) + j.length * rho[j.ID]) <= 1000 * (3 - B1[i.ID, j.ID] - p[i.ID, bin.ID] - p[j.ID, bin.ID]), name="bcon_{}_{}_{}".format(i.ID, j.ID, bin.ID))
                rcon[i.ID, j.ID, bin.ID] = model.addConstr(y[i.ID] - (y[j.ID] + j.height * (1 - rho[j.ID]) + j.length * rho[j.ID]) <= 1000 * (3 - B2[i.ID, j.ID] - p[i.ID, bin.ID] - p[j.ID, bin.ID]), name="rcon_{}_{}_{}".format(i.ID, j.ID, bin.ID))

bcon2 = {}
rcon2 = {}
for i in Items:
    for j in Items:
        for bin in Bins:
            if i.ID != j.ID:
                bcon2[i.ID, j.ID, bin.ID] = model.addConstr(-y[i.ID] + (y[j.ID] + j.height * (1 - rho[j.ID]) + j.length * rho[j.ID]) <= 1000 * (3 - B1[i.ID, j.ID] - p[i.ID, bin.ID] - p[j.ID, bin.ID]), name="bcon2_{}_{}_{}".format(i.ID, j.ID, bin.ID))
                rcon2[i.ID, j.ID, bin.ID] = model.addConstr(-y[i.ID] + (y[j.ID] + j.height * (1 - rho[j.ID]) + j.length * rho[j.ID]) <= 1000 * (3 - B2[i.ID, j.ID] - p[i.ID, bin.ID] - p[j.ID, bin.ID]), name="rcon2_{}_{}_{}".format(i.ID, j.ID, bin.ID))

left = {}
left2 = {}
for i in Items:
    for j in Items:
        for bin in Bins:
            if i.ID != j.ID:
                left[i.ID, j.ID, bin.ID] = model.addConstr(x[j.ID]-x[i.ID] <= 1000000 * (3-p[i.ID, bin.ID] - p[j.ID, bin.ID] - B1[i.ID, j.ID]), name="left_{}_{}_{}".format(i.ID, j.ID, bin.ID))
                left2[i.ID, j.ID, bin.ID] = model.addConstr(-x[j.ID]-j.length*(1-rho[j.ID]) - j.height*rho[j.ID] + x[i.ID] <= 1000000 * (3-p[i.ID, bin.ID] - p[j.ID, bin.ID] - B1[i.ID, j.ID]), name="left2_{}_{}_{}".format(i.ID, j.ID, bin.ID))

right = {}
right2 = {}
for i in Items:
    for j in Items:
        for bin in Bins:
            if i.ID != j.ID:
                right[i.ID, j.ID, bin.ID] = model.addConstr(x[j.ID]-x[i.ID] - i.length*(1-rho[i.ID]) - i.height*rho[i.ID] <= 1000000 * (3-p[i.ID, bin.ID] - p[j.ID, bin.ID] - B2[i.ID, j.ID]), name="right_{}_{}_{}".format(i.ID, j.ID, bin.ID))
                right2[i.ID, j.ID, bin.ID] = model.addConstr(x[i.ID] + i.length*(1-rho[i.ID]) + i.height*rho[i.ID] - x[j.ID] - j.length*(1-rho[j.ID]) - j.height*rho[j.ID] <= 1000000 * (3-p[i.ID, bin.ID] - p[j.ID, bin.ID] - B2[i.ID, j.ID]), name="right2_{}_{}_{}".format(i.ID, j.ID, bin.ID))

stacking = {} #items can only stand on items in the same bin.
stacking2 = {}
for bin in Bins:
    for i in Items:
        for j in Items:
            if i.ID != j.ID:
                stacking[i.ID, j.ID, bin.ID] = model.addConstr(2 * B1[i.ID, j.ID] <= p[j.ID, bin.ID] + p[i.ID, bin.ID], name="stacking_{}_{}_{}".format(i.ID, j.ID, bin.ID))
                stacking2[i.ID, j.ID, bin.ID] = model.addConstr(2 * B2[i.ID, j.ID] <= p[j.ID, bin.ID] + p[i.ID, bin.ID], name="stacking2_{}_{}_{}".format(i.ID, j.ID, bin.ID))



rotation = {}
for i in Items:
    if i.rotate == 0:
        rotation[i.ID] = model.addConstr(rho[i.ID] == 0, name="rotation_{}".format(i.ID))

fragile = {}
for j in Items:
    if j.fragile == 1:
        fragile[j.ID] = model.addConstr(gp.quicksum(B1[i.ID, j.ID] + B2[i.ID, j.ID] for i in Items if i.ID != j.ID) == 0, name="fragile_{}".format(j.ID))

radiation_activation = {}
radiation_limiter = {}
for bin in Bins:
    radiation_activation[bin.ID] = model.addConstr(gp.quicksum(i.radioactive*p[i.ID, bin.ID] for i in Items) <= 1000 * radioactive[bin.ID], name="radiation_activation_{}".format(bin.ID))
    radiation_limiter[bin.ID] = model.addConstr(gp.quicksum(i.perishable*p[i.ID, bin.ID] for i in Items) <= 1000 * (1-radioactive[bin.ID]), name="radiation_limiter_{}".format(bin.ID))

feasibility = {}
for i in Items:
    feasibility[i.ID] = model.addConstr(gamma[i.ID] + gp.quicksum(B1[i.ID,j.ID] + B2[i.ID,j.ID] for j in Items if j.ID != i.ID) + 2 * g[i.ID] >= 2, name="feasibility_{}".format(i.ID))



model.update()
model.setParam("LogFile", "log_file")
model.update()
model.write("model.lp")
model.optimize()

model.write("model.sol")
print("Optimal value: {}".format(model.objVal))
print("Optimal solution:")
for i in Items:
    for bin in Bins:
        if p[i.ID, bin.ID].X > 0.5:
            print("Item {} is assigned to Bin {}".format(i.ID, bin.ID))
            print("Item {}: x = {}, y = {}, rho = {}".format(i.ID, x[i.ID].X, y[i.ID].X, rho[i.ID].X))
            for j in Items:
                if i.ID != j.ID:
                    if B1[i.ID, j.ID].X > 0.5:
                        print("Item {} is above Item {}".format(i.ID, j.ID))
                    if B2[i.ID, j.ID].X > 0.5:
                        print("Item {} is above Item {}".format(i.ID, j.ID))




for bin in Bins:
    if z[bin.ID].X > 0.5:
        fig, ax = plt.subplots()
        ax.set_xlim(0, bin.length)
        ax.set_ylim(0, bin.height)
        for i in Items:
            if p[i.ID, bin.ID].X > 0.5:
                rect = patches.Rectangle((x[i.ID].X, y[i.ID].X), i.length * (1 - rho[i.ID].X) + i.height * rho[i.ID].X, i.height * (1 - rho[i.ID].X) + i.length * rho[i.ID].X, edgecolor='black', facecolor='blue', alpha=0.5)
                ax.add_patch(rect)
                plt.text(x[i.ID].X + (i.length * (1 - rho[i.ID].X) + i.height * rho[i.ID].X)/2, y[i.ID].X + (i.height * (1 - rho[i.ID].X) + i.length * rho[i.ID].X)/2, str(i.ID), color='white', ha='center', va='center')
        plt.title("Bin {}".format(bin.ID))
        plt.xlabel("Length")
        plt.ylabel("Height")
        plt.grid()
plt.show()