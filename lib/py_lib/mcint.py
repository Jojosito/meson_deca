import numpy as np

def integral(func, bounds, N=10000):
    """
    Compute a definite Monte Carlo integral.

    Integrate func from bounds[i][0] to bounds[i][1] for i-th variable using 
    uniform random sampling.

    Parameters
    ----------
    func : function
        A Python function or method to integrate.
    bounds : list of 2-elemental lists
        List of elements of the type [a_i,b_i], a_i and b_i are floats, lower
        and upper limits of integration of i-th variable.
    N : int, optional
        2*N is the number of drawn samples.

    Returns
    -------
    res : float
        The integral of func from bounds[:][0] to bounds[:][1]
    abserr : float
        A crude estimate of the absolute error of the integration
    """
    n_vars = len(bounds)
    random_point = np.zeros(n_vars, dtype=float)
    total_volume = np.product([(bounds[i][1] - bounds[i][0]) for i in range(n_vars)])

    res1 =  sum(__eval_integrate(func, bounds, N, n_vars)) \
           / float(N) * total_volume
    res2 =  sum(__eval_integrate(func, bounds, N, n_vars)) \
           / float(N) * total_volume

    return [(res1+res2)/2., np.abs(res1-res2)/2]


def __eval_integrate(func, m, N_points, n_vars):
    # Generator yielding func evaluated at N_points random points.
    num = 0
    while num < N_points:
        yield func(*[np.random.uniform(m[i][0], m[i][1]) for i in range(n_vars)])
        num += 1


def integral_A(func, bounds, N=1000000):
    """
    Similar function, more oriented to our purposes.

    For the function func that takes a list of floats and
    returns a complex vector A, this integral_A returns the matrix
    I[i,j] = \int A*[i] A[j].
    """
    n_vars = len(bounds)
    random_point = np.zeros(n_vars, dtype=float)
    total_volume = np.product([(bounds[i][1] - bounds[i][0]) for i in range(n_vars)])

    res1 =  sum(__eval_integrate_A(func, bounds, N, n_vars)) \
           / float(N) * total_volume
    res2 =  sum(__eval_integrate_A(func, bounds, N, n_vars)) \
           / float(N) * total_volume

    return [(res1+res2)/2., np.abs(res1-res2)/2]

def __eval_integrate_A(func, m, N_points, n_vars):
    # Generator yielding I for func evaluated at N_points random points,
    # as described above.
    num = 0
    while num < N_points:
        A = func([np.random.uniform(m[i][0], m[i][1]) for i in range(n_vars)])
        n_res = A.shape[0]
        yield np.asarray([[np.conj(A[i1]) * A[i2] for i1 in range(n_res)] for i2 in range(n_res)])
        num += 1
