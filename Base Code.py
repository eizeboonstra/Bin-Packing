import gurobipy as gp
import pickle

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

I_values = list(I.values())
Items = []
for i in range(len(I_values)):
    Items.append(Item(i, I_values[i][0], I_values[i][1], I_values[i][2], I_values[i][3], I_values[i][4], I_values[i][5]))


model = gp.Model("Bin Packing Problem")


z = {}
for b in Bins:
    z[b.ID] = model.addVar(obj= b.cost, vtype=gp.GRB.BINARY, name="z_{}".format(b.ID))

I = {}
for i in Items:
    for j in Items:
        if i.ID != j.ID:
            I[i.ID, j.ID] = model.addVar(vtype=gp.GRB.BINARY, name="I_{}_{}".format(i.ID, j.ID))

b = {}
for i in Items:
    for j in Items:
        if i.ID != j.ID:
            b[i.ID, j.ID] = model.addVar(vtype=gp.GRB.BINARY, name="b_{}_{}".format(i.ID, j.ID))

p = {}
for i in Items:
    for bin in Bins:
        p[i.ID, bin.ID] = model.addVar(vtype=gp.GRB.BINARY, name="p_{}_{}".format(i.ID, bin.ID))

g = {}
rho = {}
for i in Items:
    rho[i.ID] = model.addVar(vtype=gp.GRB.BINARY, name="rho_{}".format(i.ID))
    g[i.ID] = model.addVar(vtype=gp.GRB.BINARY, name="g_{}".format(i.ID))


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





overlapping = {}
for bin in Bins:
    for i in Items:
        for j in Items:
            if i.ID != j.ID:
                overlapping[bin.ID, i.ID, j.ID] = model.addConstr(I[i.ID, j.ID] + I[j.ID, i.ID] + b[i.ID, j.ID] + b[j.ID, i.ID] >= p[i.ID, bin.ID] + p[j.ID, bin.ID] - 1, name="overlapping_{}_{}_{}".format(bin.ID, i.ID, j.ID))

horizontal1 = {}
horizontal2 = {}

for i in Items:
    for j in Items:
        if i.ID != j.ID:
            horizontal1[i.ID, j.ID] = model.addConstr(x[j.ID] >= x[i.ID] + i.length * (1 - rho[i.ID]) + i.height * rho[i.ID] - Lmax * (1- I[i.ID, j.ID]), name="horizontal1_{}_{}".format(i.ID, j.ID))
            horizontal2[i.ID, j.ID] = model.addConstr(x[j.ID] + 0.00001 <= x[i.ID] + i.length * (1 - rho[i.ID]) + i.height * rho[i.ID] + Lmax *  I[j.ID, i.ID], name="horizontal2_{}_{}".format(i.ID, j.ID))

xbound = {}
ybound = {}
assigned = {}
for i in Items:
    xbound[i.ID] = model.addConstr(x[i.ID] + i.length * (1 - rho[i.ID]) + i.height * rho[i.ID] <= gp.quicksum(bin.length * p[i.ID, bin.ID] for bin in Bins), name="xbound_{}".format(i.ID))
    ybound[i.ID] = model.addConstr(y[i.ID] + i.height * (1 - rho[i.ID]) + i.length * rho[i.ID] <= gp.quicksum(bin.height * p[i.ID, bin.ID] for bin in Bins), name="ybound_{}".format(i.ID))
    assigned[i.ID] = model.addConstr(gp.quicksum(p[i.ID, bin.ID] for bin in Bins) == 1, name="assigned_{}".format(i.ID))

for bin in Bins:
    model.addConstr(gp.quicksum(p[i.ID, bin.ID] for i in Items) <= 1000 * z[bin.ID], name="bin_usage_{}".format(bin.ID))

grounding = {}
for i in Items:
    grounding[i.ID] = model.addConstr(y[i.ID] <= 1000 * (1-g[i.ID]), name="grounding_{}".format(i.ID))

bcon = {}
for i in Items:
    for j in Items:
        for bin in Bins:
            if i.ID != j.ID:
                bcon[i.ID, j.ID, bin.ID] = model.addConstr(y[i.ID] - (y[j.ID] + j.height * (1 - rho[j.ID]) + j.length * rho[j.ID]) <= 1000 * (3 - B1[i.ID, j.ID] - p[i.ID, bin.ID] - p[j.ID, bin.ID]), name="bcon_{}_{}_{}".format(i.ID, j.ID, bin.ID))







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

for i in Items:
    for j in Items:
        if i.ID != j.ID:
            if I[i.ID, j.ID].X > 0.5:
                print("Item {} is to the left of Item {}".format(i.ID, j.ID))
            if b[i.ID, j.ID].X > 0.5:
                print("Item {} is below Item {}".format(i.ID, j.ID))