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
