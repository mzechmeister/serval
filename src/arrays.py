import numpy as np

class Arrays(np.ndarray):
    '''
    Broadcast operations on object type ndarrays.

    Subclassing see https://numpy.org/doc/stable/user/basics.subclassing.html.

    Examples
    --------
    >>> a = Arrays([[1,2,3,4], [7,8], [5,6,9,10]])
    >>> a
    Arrays([array([1, 2, 3, 4]), array([7, 8]), array([ 5,  6,  9, 10])],
           dtype=object)
    >>> a[0:2] + 1
    Arrays([array([2, 3, 4, 5]), array([8, 9])], dtype=object)
    >>> np.sum(a)
    55
    >>> np.sum(a, axis=1)
    array([10, 15, 30])
    >>> np.nansum(a, axis=1)
    array([10, 15, 30])
    >>> np.nansum(a)
    55
    >>> print(a[:,1:-1])
    [[2 3],
     [],
     [6 9]]
    >>> print(np.sin(a))
    [[ 0.84147098  0.90929743  0.14112001 -0.7568025 ],
     [0.6569866  0.98935825],
     [-0.95892427 -0.2794155   0.41211849 -0.54402111]]
    >>> print(a[0:2])
    [[1 2 3 4],
     [7 8]]
    >>> a[0:1]
    Arrays([array([1, 2, 3, 4])], dtype=object)
    >>> print(a[-1:])
    [[ 5  6  9 10]]
    >>> print(a + 2*a)
    [[ 3  6  9 12],
     [21 24],
     [15 18 27 30]]
    >>> print(a[a>6])
    [[],
     [7 8],
     [ 9 10]]

    '''
    def __new__(cls, input_array, dims=None):
        obj = np.array(list(map(np.array, input_array))).view(cls)
        return obj
    def ravel(self):
        return np.hstack(self)
    def astype(self, dtype):
        return Arrays(x.astype(dtype) for x in self)
    def __str__(self):
        return '[' + ",\n ".join(map(str, self)) + ']'
    def __setitem__(self, ij, val):
        if isinstance(ij, Arrays):
            # apply setitem to the subdimension
            [si.__setitem__(j, valj) for si,j,valj in zip(self, ij, val)]
        else:
            super(Arrays, self).__setitem__(ij, val)
    def __getitem__(self, ij):
        if isinstance(ij, tuple) and len(ij) > 1:
            i, j = ij
            # handle twodimensional slicing
            if isinstance(i, slice) or hasattr(i, '__iter__'):
                # [1:4,:] or [[1,2,3],[1,2]]
                return Arrays(arr[j] for arr in self[i])
            return self[i][j] # [1,:] np.array
        elif isinstance(ij, Arrays):
            # boolean indexing, e.g. a[a>0]
            # output is not flatten, dimensions are kept
            return Arrays(a[j] for a,j in zip(self,ij))
        return super(Arrays, self).__getitem__(ij)
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # this method is called whenever you use reduce (sum()), __call__ (sin), etc.
        '''this implementation of __array_ufunc__ makes sure that all custom attributes are maintained when a ufunc operation is performed on our class.'''
        axis = kwargs.pop('axis', None)
        # repeat scalar inputs before zipping
        pad_inputs = [arg if hasattr(arg, '__iter__') else [arg]*len(self) for arg in inputs]
        result = [np.ndarray.__array_ufunc__(self, ufunc, method, *x, **kwargs) for x in zip(*pad_inputs)]
        #print(len(inputs[0]), len(self))
        if method == 'reduce':
            # handle sum, min, max, etc.
            if axis == 1:
                return np.array(result)
            else:
                # repeat over remaining axis
                return np.ndarray.__array_ufunc__(self, ufunc, method, result, **kwargs)
        return Arrays(result)

# patch for nanfunction that cannot handle the object-ndarrays
def nanpatch(func):
    def wrapper(a, axis=None, **kwargs):
        if isinstance(a, Arrays):
            rowresults = [func(x, **kwargs) for x in a]
            if axis == 1:
                return np.array(rowresults)
            else:
                # repeat over remaining axis
                return func(rowresults)
        # otherwise keep the original version
        return func(a, axis=axis, **kwargs)
    return wrapper

np.nanmean = nanpatch(np.nanmean)
np.nansum = nanpatch(np.nansum)
np.nanmin = nanpatch(np.nanmin)
np.nanmax = nanpatch(np.nanmax)


