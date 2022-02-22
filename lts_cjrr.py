import jax.numpy as np
from jax import grad, jit
from jax import random

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


def gradient_descent(W, inputs, targets):
    gradient = jit(grad(loss, argnums=1))
    print(f"Start loss {loss(inputs, W, targets)}")
    for epoch in range(10000):
        W -= LR * gradient(inputs, W, targets)
    return W


def find_linear_transform(inputs, targets, W):
    #W = np.array([
    #    [-1., 1., -1.],
    #    [1., 1., -1.],
    #    [-1., -1., 1.],
    #])
    W = gradient_descent(W, inputs, targets)

    print(f"Final loss {loss(inputs, W, targets)}")
    return W

X = np.array([
    [4., 7., 3.],
    [1., 8., 1.],
    [-5., -6., -1],
    [3., -1., -2.],
    [0., 9., 6.]
])

W = np.array([
    [4., 7., 3.],
    [1., 8., -1.],
    [-5., -6., 11.],
])

def main(X,W):
    y = X.dot(W)
    print(y)
    approximated_W = find_linear_transform(X, y, W)
    print(f"Original W matrix {W}")
    print(f"final W matrix {approximated_W}")


main()
