import jax.numpy as np
from jax import grad, jit
from jax import random
import time

LR = 0.01


def loss(X, W, y):
    """Loss function for a linear regression. A forward pass of our model.

    Args:
        X: a matrix of input vectors.
        W: A matrix approximation.
        y: a matrix of target vectors.

    Returns:
        scalar: a cost of this solution.
    """

    y_hat = X.dot(W)  # Predict values.
    return ((y_hat - y) ** 2).mean()  # Return cost.


def gradient_descent(W_guess, inputs, targets):
    gradient = jit(grad(loss, argnums=1))
    print(f"Start loss {loss(inputs, W_guess, targets)}")
    for epoch in range(10000):
        W_guess -= LR * gradient(inputs, W_guess, targets)
        print(f"Iteration loss {loss(inputs, W_guess, targets)}")
    return W_guess


def find_linear_transform(inputs, targets):
    assert inputs.shape == targets.shape, "Inputs and targets have missmatching shapes"
    assert len(inputs.shape) == 2, "Your input vectors are wrong dimensions, Carlos!"
    W_shape = (inputs.shape[1], inputs.shape[1])
    key = random.PRNGKey(0)
    W_guess = random.normal(key, W_shape)  # <--- This is my initial incorrect guess about what the true W is.

    # This method will iteratively try to fix the initial guess to match the input/output pairs by minimising the
    # loss/mismatch.
    W_guess = gradient_descent(W_guess, inputs, targets)

    print(f"Final loss {loss(inputs, W_guess, targets)}")
    return W_guess


def main(X, true_W, y):
    W_guess = find_linear_transform(X, y)  # <--- this function should not take the correct W in any case, LoL
    print(f"Original W matrix {true_W}")
    print(f"final W matrix {W_guess}")
    return(W_guess)

exp_thre=np.array([1.96876885, 0.        , 0.        , 0.        , 0.        ]) # Thresholds. If the gene is expressed above this number, it will switch on, otherwise it will remain off.
decay_lambdas=np.array([0.        , 1.42894643, 1.38985741, 1.65910196, 0.17701693]) # Exponential decay exponent lambda for each gene, the quantity of each gene at time tn will decay at this rate in tn+1
grn=np.array([[-1.86793303,  1.86367223, -0.37709589,  0.39320219,  0.        ],
       [ 0.        ,  0.        ,  0.        , -0.35279367,  0.29892912],
       [-0.36162465,  0.        ,  0.        ,  0.        ,  0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
       [ 1.89152412,  0.        , -1.49580545, -1.76727241,  0.        ]]) # the grn, which is also the 'true_W'

dev=np.array([[ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ],
       [-0.        ,  0.        , -0.        ,  0.        ,  1.89152412],
       [ 0.        ,  0.47369919,  0.        ,  0.        ,  1.58465385],
       [ 0.        ,  0.51032882,  0.        ,  0.        ,  1.32756849],
       [ 0.        ,  0.45472132,  0.        ,  0.        ,  1.11219122],
       [ 0.        ,  0.38746245,  0.        ,  0.        ,  0.93175555],
       [ 0.        ,  0.32616287,  0.        ,  0.        ,  0.78059275],
       [ 0.        ,  0.27362179,  0.        ,  0.        ,  0.65395376],
       [ 0.        ,  0.22932046,  0.        ,  0.        ,  0.54785997],
       [ 0.        ,  0.19213824,  0.        ,  0.        ,  0.45897824],
       [ 0.        ,  0.16097195,  0.        ,  0.        ,  0.38451618],
       [ 0.        ,  0.13485799,  0.        ,  0.        ,  0.32213443],
       [ 0.        ,  0.11297968,  0.        ,  0.        ,  0.26987315],
       [ 0.        ,  0.09465056,  0.        ,  0.        ,  0.22609044],
       [ 0.        ,  0.079295  ,  0.        ,  0.        ,  0.18941079],
       [ 0.        ,  0.06643063,  0.        ,  0.        ,  0.15868185]]) # The development of the GRN through 15 steps, with the inital one being a standardized, user-defined seed.


# Equation to get the vector v1 from v0:
v0=dev[0] # take the starting vector as v0, this is identical always, and user defined (as I mentioned above).
decayed_v0=(lambda x, l: x*np.exp(-l))(vn,decay_lambdas) # Decay the initial amount by one step
exp_change = np.matmul(grn,decayed_v0) # Calculate the regulatory effects of the decayed quantities
pre_thresholds = exp_change + decayed_v0 # add the decayed amounts to the regulatory effects
thresh_rectifier = (pre_thresholds > exp_thre).astype(int) # a vector to rectify the resulting values to their thresholds
v1 = pre_thresholds * thresh_rectifier #rectifying for the thresholds

X=dev[0:15]
y=dev[1:16]

if __name__ == "__main__":
    X = np.array([  # These are my input vectors arranged in a matrix
        [4., 7., 3.],  # <--- ...
        [1., 8., 1.],  # <--- ...
        [-5., -6., -1],  # <--- ...
        [3., -1., -2.],
        [0., 9., 6.]
    ])

    true_W = np.array([  # This is a matrix I'm trying to learn, in this case I know the correct one since I defined
        [4., 7., 3.],  # it here. Normally, I don't have this matrix and I'm trying to infer it with this method
        [1., 8., -1.],
        [-5., -6., 11.],
    ])
    y = X.dot(true_W)  # These are the outputs of the transformation having the inputs from earlier. This is something
    #  I will have, and that I will use to infer the W matrix (the unknown transformation)
    main(X, true_W, y)
