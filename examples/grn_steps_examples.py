import jax.numpy as np

# Thresholds. If the gene is expressed above this number, it will switch on, otherwise it will remain off.
exp_thre = np.array([1.96876885, 0.0, 0.0, 0.0, 0.0])
# Exponential decay exponent lambda for each gene, the quantity of each gene at time tn will decay at this rate in tn+1
decay_lambdas = np.array([0.0, 1.42894643, 1.38985741, 1.65910196, 0.17701693])

# the grn, which is also the 'true_W'
grn = np.array(
    [
        [-1.86793303, 1.86367223, -0.37709589, 0.39320219, 0.0],
        [0.0, 0.0, 0.0, -0.35279367, 0.29892912],
        [-0.36162465, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [1.89152412, 0.0, -1.49580545, -1.76727241, 0.0],
    ]
)

# The development of the GRN through 15 steps, with the inital one being a standardized, user-defined seed.
dev = np.array(
    [
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [-0.0, 0.0, -0.0, 0.0, 1.89152412],
        [0.0, 0.47369919, 0.0, 0.0, 1.58465385],
        [0.0, 0.51032882, 0.0, 0.0, 1.32756849],
        [0.0, 0.45472132, 0.0, 0.0, 1.11219122],
        [0.0, 0.38746245, 0.0, 0.0, 0.93175555],
        [0.0, 0.32616287, 0.0, 0.0, 0.78059275],
        [0.0, 0.27362179, 0.0, 0.0, 0.65395376],
        [0.0, 0.22932046, 0.0, 0.0, 0.54785997],
        [0.0, 0.19213824, 0.0, 0.0, 0.45897824],
        [0.0, 0.16097195, 0.0, 0.0, 0.38451618],
        [0.0, 0.13485799, 0.0, 0.0, 0.32213443],
        [0.0, 0.11297968, 0.0, 0.0, 0.26987315],
        [0.0, 0.09465056, 0.0, 0.0, 0.22609044],
        [0.0, 0.079295, 0.0, 0.0, 0.18941079],
        [0.0, 0.06643063, 0.0, 0.0, 0.15868185],
    ]
)
# Equation to get the vector v1 from v0:
# take the starting vector as v0, this is identical always, and user defined (as I mentioned above).
v0 = dev[0]
# Decay the initial amount by one step
decayed_v0 = (lambda x, l: x * np.exp(-l))(v0, decay_lambdas)
# Calculate the regulatory effects of the decayed quantities
exp_change = np.matmul(grn, decayed_v0)
# add the decayed amounts to the regulatory effects
pre_thresholds = exp_change + decayed_v0
# a vector to rectify the resulting values to their thresholds
thresh_rectifier = (pre_thresholds > exp_thre).astype(int)
v1 = pre_thresholds * thresh_rectifier  # rectifying for the thresholds

X = dev[0:15]
y = dev[1:16]
true_W = grn
