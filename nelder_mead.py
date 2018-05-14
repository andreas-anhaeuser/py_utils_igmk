#!/usr/bin/python
"""Minimization of scalar functions using Nelder-Mead algorithm.

    This solver is intended to be an improvement of the Nelder-Mead function
    in `scipy.optimize.minizme`.


    Improvements are:
    =================
    - The output now also includes the history of the function values on the
      vertices of each iteration step.
    - The codes has been completely re-written in order to make it more
      readable.
    - Documentation has been added.


    Terms
    =====
        Simplex
        -------
        A simplex is the higher-dimensional generalisation of a triangle on a
        plane or a tetrahedron in space. It exists in N-dimensional
        (hyper-)space and has (N+1) "corners".

        Vertex
        ------
        A "corner" of a simplex.

        Centroid
        --------
        The arithmetic mean of a simplex's vertices (i. e. its "center" point).


    Method
    ======
    The algorithm tries to minimize a scale-valued function in N-dimensional
    argument space. An N-dim simplex is initialized an the iteratively moved
    and reshaped depending on the function values on its vertices. If
    successful, this will eventually lead to a contraction of the simplex
    around the function minimum. For details see docstring of the main function.


    Author
    ======
    Andreas Anhaeuser (AA)
    University of Cologne
    <andreas.anhaeuser@posteo.net>
"""

from copy import deepcopy as copy
import numpy as np
from scipy.optimize import OptimizeResult

_alpha = 1.     # reflection factor
_gamma = 2.     # expansion factor
_rho = 0.5      # contraction factor
_sigma = 0.5    # shrinking factor

_fatol = 1e-4   # function difference tolerance
_xatol = 1e-4   # parameter difference tolerance
_maxiter = 200  # maximum number of iterations
_maxfev = 400   # maximum number of function evaluation

# simplex : ndarray, length N
# vertex : ndarray, shape (N+1) x N

###################################################
# MAIN                                            #
###################################################
def minimize_nelder_mead(
        fun,
        x0,
        args=(),
        callback=None,
        maxiter=None,
        maxfev=None,
        disp=False,
        return_all=False,
        initial_simplex=None,
        xatol=None,
        fatol=None,
        adaptive=False,
        x0_uncert=None,
    ):
    """Return a dict.

        Parameters
        ----------
        fun : callable
            Function that is to be minimized.
        x0 : array of shape (N,)
            Initial guess of the parameters to be optimized.
        args : tuble, optional
            Additional parameters (other than `x0`) passed to the objective
            function `fun`.
        callback : callable, optional
            Called after each iteration. If callback returns True the algorithm
            execution is terminated. The signature is:
            
                ``callback(xk)``
            
            where ``xk`` is the current parameter vector.
        maxiter : int, optional
            Maximum number of iterations.
        maxfev : int, optional
            Maximum number of function evaluations.
        disp : bool, optional
            Print convergence messages. Default: False
        return_all : bool, optional
            Return complete iteration record. Default: False
        initial_simplex : array of shape (N + 1, N), optional
            If given, overrides `x0` and `x0_uncert`. ``initial_simplex[j,:]``
            contains the coordinates of the j-th vertex of the ``N+1`` vertices
            in the simplex, where ``N`` is the dimension.
        xatol : array of length N, optional
            Maximum difference of the vertex coordinates upon termination
        fatol : float, optional
            Maximum difference of function value to the previous one upon
            termination
        adaptive : bool, optional
            Unused remnant of the original numpy implementation.
            Caution: If True, function will exit with ``NotImplementedError``.
            Default: False.
        x0_uncert : array of shape (N,), optional
            Initial uncertaintainties of the parameters in `x0`. If
            `initial_simplex` is given, `x0_uncert` remains unused. Otherwise,
            it is used to construct `initial_simplex`. It neither of the two is
            given, the function will still run (by just making up some
            `x0_uncert` itself), but this is discouraged as the function has no
            means to estimate `x0_uncert` is a sensible way.

        Returns
        -------
        scipy.optimize.OptimizeResult

        History
        -------
        2018-04-28 (AA): Created
        2018-05-14 (AA): Conformed to scipy.optimize._minimize_neldermead
    """
    ###################################################
    # DEFAULTS                                        #
    ###################################################
    if fatol is None:
        fatol = _fatol
    if xatol is None:
        xatol = _xatol 
    if maxiter is None:
        maxiter = _maxiter
    if maxfev is None:
        maxfev = _maxfev

    if x0_uncert is None:
        x0_uncert = make_up_uncertainties(x0)

    if initial_simplex is None:
        initial_simplex = get_initial_simplex(x0, x0_uncert)

    if adaptive:
        raise NotImplementedError('Adaptive parameters not implemented. ' + 
                'Set `adaptive` to False.')
    
    ###################################################
    # INITIALIZE                                      #
    ###################################################
    options = {
            'maxiter' : maxiter,
            'maxfev' : maxfev,
            'fatol' : fatol,
            'xatol' : xatol,
            }

    M, N = np.shape(initial_simplex)
    table = initialize_lookup_table(initial_simplex[0])
    record = OptimizeResult(
            {
            'all_simplices' : [],
            'best_vertices' : [],
            'worst_vertices' : [],
            'all_fun_values': [],
            'best_fun_values' : [],
            'worst_fun_values' : [],
            'lookup_table' : table,
            'action' : [],              # mode in which simplex is manipulated
            }
            )

    ###################################################
    # ABBREVIATIONS                                   #
    ###################################################
    f = lambda vertex: get_function_value(vertex, fun, args, table)
    F = lambda simplex: get_function_values(simplex, fun, args, table)

    # action flags
    _INITIALIZE = 0
    _SHRINK = 1
    _CONTRACT = 2
    _REFLECT = 3
    _EXPAND = 4

    ###################################################
    # ITERATE                                         #
    ################################################### 
    simplex = initial_simplex
    record['action'].append(_INITIALIZE)

    while True:
        ###################################################
        # EVALUATE                                        #
        ###################################################
        fun_values = F(simplex)

        # best
        m_best = np.argmin(fun_values)
        x_best = simplex[m_best]
        f_best = fun_values[m_best]                    

        # worst
        m_worst = np.argmax(fun_values)
        x_worst = simplex[m_worst]
        f_worst = fun_values[m_worst]

        ###################################################
        # RECORD                                          #
        ###################################################
        record['all_simplices'].append(copy(simplex))
        record['best_vertices'].append(copy(x_best))
        record['worst_vertices'].append(copy(x_worst))
        record['all_fun_values'].append(copy(fun_values))
        record['best_fun_values'].append(copy(f_best))
        record['worst_fun_values'].append(copy(f_worst))

        ###################################################
        # TERMINATE                                       #
        ###################################################
        if termination_reached(record, options):
            break

        ###################################################
        # FIND A BETTER SIMPLEX                           #
        ###################################################
        # centroid
        x_centroid = get_centroid_vertex(simplex, m_worst)

        # reflected
        x_reflected = get_reflected_vertex(x_worst, x_centroid)
        f_reflected = f(x_reflected)

        # check quality of reflected vertex
        if f_reflected < f_best:
            # x_reflected is best
            # -> try expansion
            x_expanded = get_expanded_vertex(x_reflected, x_centroid)
            f_expanded = f(x_expanded)
            if f_expanded < f_reflected:
                # expanded is better
                simplex[m_worst] = x_expanded
                record['action'].append(_EXPAND)
            else:
                # reflected is better
                simplex[m_worst] = x_reflected
                record['action'].append(_REFLECT)

        elif f_reflected > f_worst:
            # x_reflected is worst
            # -> try contraction
            x_contracted = get_contracted_vertex(x_worst, x_centroid)
            f_contracted = f(x_contracted)
            if f_contracted < f_worst:
                # contracted is better
                # -> use it
                simplex[m_worst] = x_contracted
                record['action'].append(_CONTRACT)
            else:
                # contracted is even worse (this can happen close to minima)
                # -> shrink
                simplex = get_shrunk_simplex(simplex, x_centroid, m_best)
                record['action'].append(_SHRINK)

        else:
            # x_reflected is neither worst nor best
            # -> use it
            simplex[m_worst] = x_reflected
            record['action'].append(_REFLECT)

        if callback is not None:
            terminate = callback(x_best)
            if terminate is True:
                break

    ###################################################
    # CREATE OptimizeResult OBJECT                    #
    ###################################################
    # list -> array
    for key in record.keys():
        if key in ['lookup_table']:
            continue
        record[key] = np.array(record[key])

    record['Nfeval'] = len(table['fun_values'])
    record['Niter'] = len(record['all_simplices'])

    record['x0'] = x0
    if x0_uncert is not None:
        record['x0_uncert'] = x0_uncert
    else:
        record['x0_uncert'] = np.ptp(initial_simplex, 0)
    record['f0'] = get_function_value(x0, fun, args, table)
    record['x'] = x_best
    record['f'] = f_best
    record['success'] = convergence_reached(record, options)
    record['message'] = get_termination_message(record, options)

    ###################################################
    # DISP                                            #
    ###################################################
    if disp:
        display_termination_message(record)

    return record

###################################################
# MISC                                            #
###################################################
def get_initial_simplex(x0, x0_uncert):
    """Return a (N+1 x N)-array."""
    # input check
    shape = np.shape(x0)
    assert len(shape) == 1
    assert np.shape(x0_uncert) == shape

    # initialize
    N = shape[0]
    simplex = np.nan * np.zeros((N+1, N))

    # ========== create vertices  ======================== #
    for n in range(-1, N):
        vertex = 1. * copy(x0)

        # manipulate n-th coordinate
        if n >= 0:
            vertex[n] += 0.5 * x0_uncert[n]

        # add vertex to simplex
        simplex[n] = vertex
    # ==================================================== #

    # output check
    assert np.sum(np.isnan(simplex)) == 0

    return simplex

def make_up_uncertainties(x0):
    unc = 0.1 * np.abs(x0) 
    unc[unc==0] = 0.1
    return unc

###################################################
# TERMINATION                                     #
###################################################
def display_termination_message(record, options=None):
    message = 'Succes                  : %s\n' % record['success'] + \
              'Terminatation reason    : %s\n' % record['message'] + \
              'Result                  : %s\n' % record['x'] + \
              '- uncertainty           : %s\n' % parameter_uncertainty(record,
                                                        options) + \
              'Function value (result) : %f\n' % record['f'] + \
              '- uncertainty           : %f\n' % function_uncertainty(record, 
                                                        options) + \
              'Ititial guess           : %s\n' % record['x0'] + \
              '- uncertainty           : %s\n' % record['x0_uncert'] + \
              'Ititial function value  : %f\n' % record['f0'] + \
              'Iterations              : %i\n' % record['Niter'] + \
              'Function evaluations    : %i\n' % record['Nfeval'] + \
              ''

    print(message)

def get_termination_message(record, options):
    """Return a str."""
    # convergence
    if convergence_reached(record, options):
        return 'Convergence reached.'
    
    # iteration limit
    itlim, message = iteration_limit_reachead(record, options)
    if itlim:
        return message

    # other
    return 'Termination reason unknown.'

def termination_reached(record, options):
    """Return a bool."""
    iteration_criterion = iteration_limit_reached(record, options)[0]
    convergence_criterion = convergence_reached(record, options)
    terminate = iteration_criterion or convergence_criterion
    return terminate

def iteration_limit_reached(record, options):
    """Return (bool, str)."""
    all_simplices = record['all_simplices']
    table = record['lookup_table']

    Niter = len(all_simplices)
    Nfeval = len(table['fun_values'])

    if Niter < 2:
        return False, ''

    # maxiter
    if Niter > options['maxiter']:
        return True, 'Maximum number of iterations reached.'

    # maxfev
    if Nfeval > options['maxfev']:
        return True, 'Maximum number of function evaluations reached.'

    return False, ''

def convergence_reached(record, options):
    """Return a bool."""
    # fatol
    span = function_uncertainty(record, options)
    fatol = options['fatol']
    if span < fatol:
        return True

    # xatol
    spans = parameter_uncertainty(record, options)
    xatol = options['xatol']
    if np.sum(spans >= xatol) == 0:
        return True

    return False

def function_uncertainty(record, options):
    """Return a float.

        Return maximum difference in function values of all points of the last
        two simplices. Return infinity if history is shorter than 2.
    """
    last_fun_values = record['all_fun_values'][-2:]

    # history too short:
    if len(last_fun_values) < 2:
        return np.inf

    return np.ptp(last_fun_values)

def parameter_uncertainty(record, options):
    """Return an array of shape (N,).

        Return maximum difference in parameter values of all points of the last
        two simplices. Return infinity if history is shorter than 2.
    """
    all_simplices = record['all_simplices']
    K, M, N = np.shape(all_simplices)

    # history too short:
    if K < 2:
        return np.ones(N) * np.inf

    # retrieve and reshape last two simplices
    last_simplices = np.array(all_simplices[-2:])
    last_simplices_reshaped = last_simplices.reshape((M*2, N))

    # return maximum differences
    return np.ptp(last_simplices_reshaped, 0)

###################################################
# MANIPULATE SIMPLEX                              #
###################################################
def get_centroid_vertex(simplex, exclude_index):
    """Return a vertex."""
    # ========== input check  ============================ #
    shape = np.shape(simplex)
    assert len(shape) == 2
    M, N = shape
    assert M == N + 1
    assert exclude_index < M
    # ==================================================== #

    # exclude
    m = exclude_index
    vertices = np.concatenate((simplex[:m], simplex[m+1:]), 0)

    # compute
    centroid = np.mean(vertices, 0)
    return centroid

def get_reflected_vertex(vertex, centroid):
    """Return a vertex."""
    assert _alpha > 0
    return centroid + _alpha * (centroid - vertex)

def get_expanded_vertex(vertex, centroid):
    """Return a vertex."""
    assert _gamma > 1
    return centroid + _gamma * (vertex - centroid)

def get_contracted_vertex(vertex, centroid):
    """Return a vertex."""
    assert 0 < _rho < 1
    return centroid + _rho * (vertex - centroid)

def get_shrunk_simplex(simplex, centroid, best_index):
    """Return a simplex."""
    # ========== input check  ============================ #
    shape = np.shape(simplex)
    assert len(shape) == 2
    M, N = shape
    assert M == N + 1
    assert 0 <= best_index < M

    assert 0 < _sigma < 1
    # ==================================================== #

    best_vertex = simplex[best_index]
    new_simplex = np.nan * np.ones_like(simplex)

    for m in range(M):
        vertex = simplex[m]

        if m != best_index:
            # regular case: shrink towards best_vertex
            new_vertex = best_vertex + _sigma * (vertex - best_vertex)
        else:
            # special case: best vertex
            new_vertex = copy(vertex)

        new_simplex[m] = new_vertex

    # output check
    assert np.sum(np.isnan(new_simplex)) == 0

    return new_simplex

###################################################
# LOOKUP TABLE                                    #
###################################################
def initialize_lookup_table(vertex):
    N = len(vertex)
    return {
            'vertices' : np.ones((0, N)),
            'fun_values' : np.ones(0),
            }

def get_function_value_from_lookup_table(lookup_table, vertex):
    """Return result if vertex is in lookup_table, else None."""
    index = find_index_in_lookup_table(lookup_table, vertex)
    if index is None:
        result = None
    else:
        result = lookup_table['fun_values'][index]

    return result

def find_index_in_lookup_table(lookup_table, vertex):
    """Return index if vertex is in lookup_table, else None."""
    vertices = lookup_table['vertices']
    K = len(vertices)
    index = None

    # try to find vertex in list of vertices
    for k in range(K):
        is_equal = np.sum(vertices[k] != vertex) == 0
        if is_equal:
            index = k
            break

    return index

def add_function_value_to_loopkup_table(lookup_table, vertex, result):
    index = find_index_in_lookup_table(lookup_table, vertex)
    if index is not None:
        return 

    table = lookup_table
    K = len(table['vertices'])
    table['vertices'] = np.insert(table['vertices'], K, vertex, 0)
    table['fun_values'] = np.insert(table['fun_values'], K, result, 0)

###################################################
# FORWARD FUNCTION                                #
###################################################
def get_function_values(simplex, fun, args, lookup_table):
    """Return an array."""
    # input check
    shape = np.shape(simplex)
    assert len(shape) == 2
    M, N = shape
    assert M == N + 1

    # initialize
    fvals = np.nan * np.ones(M)

    # fill array
    for m in range(M):
        vertex = simplex[m]
        fvals[m] = get_function_value(vertex, fun, args, lookup_table)

    return fvals

def get_function_value(vertex, fun, args, lookup_table):
    """Return a float."""
    # input check
    shape = np.shape(vertex)
    assert len(shape) == 1

    # check loopup_table
    result = get_function_value_from_lookup_table(lookup_table, vertex)

    # call function
    if result is None:
        result = fun(vertex, *args)
        add_function_value_to_loopkup_table(lookup_table, vertex, result)

    assert np.shape(result) == ()
    return result

###################################################
# TEST FUNCTIONS                                  #
###################################################
def himmelblau(a):
    x = a[0]
    y = a[1]
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def styblinski(a):
    N = len(a)
    return np.sum(a**4 - 16 * a**2 + 5 * a) + N * 2 * 39.16617

###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    x0 = - np.array([2.2, 1.5, 2.0])
    x0_uncert = np.array([5., 5., 5.])
    # f = himmelblau
    f = styblinski

    record = minimize_nelder_mead(f, x0, x0_uncert=x0_uncert, disp=True)
