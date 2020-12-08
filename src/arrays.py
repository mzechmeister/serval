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

    '''
    def __new__(cls, input_array, dims=None):
        obj = np.array(list(map(np.array, input_array))).view(cls)
        return obj
    def __array_finalize__(self, obj):
        if obj is None:  # __new__ handles instantiation
            return
    def ravel(self):
        return np.hstack(self)
    def astype(self, typ):
        return self
    def __str__(self):
        return '['+",\n ".join(map(str, self))+']'
    def __getitem__(self, ij):
        if isinstance(ij, tuple) and len(ij) > 1:
            # handle twodimensional slicing
            if isinstance(ij[0],slice) or hasattr(ij[0], '__iter__'):
                # [1:4,:] or [[1,2,3],[1,2]]
                return Arrays([arr[ij[1]] for arr in self[ij[0]]])
            return self[ij[0]][ij[1]] # [1,:] np.array
        return super(Arrays, self).__getitem__(ij)
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):  # this method is called whenever you use a ufunc
        '''this implementation of __array_ufunc__ makes sure that all custom attributes are maintained when a ufunc operation is performed on our class.'''
        dimk = [hasattr(arg, '__iter__') and len(arg) or 1 for arg in inputs]
        dim = max(dimk)
        pad_inputs = [([i]*dim if (d<dim) else i) for d,i in zip(dimk, inputs)]
        result = [np.ndarray.__array_ufunc__(self, ufunc, method, *x, **kwargs) for x in   zip(*pad_inputs)]         
        if method == 'reduce':
           return np.array(result)
        return Arrays(result)

# patch for nanfunction that cannot handle the object-ndarrays along with second axis=-1
def nanpatch(func):
    def wrapper(a, axis=None, **kwargs):
        if isinstance(a, Arrays):
           if axis==1:
              #return func(a, **kwargs)
              return np.array([func(x, **kwargs) for x in a])
           else:
              #return func(func(a, **kwargs))
              return func(np.array([func(x, **kwargs) for x in a]))
        return func(a, axis=axis, **kwargs)
    return wrapper

np.sum = nanpatch(np.sum)
np.nanmean = nanpatch(np.nanmean)
np.nansum = nanpatch(np.nansum)
np.nanmin = nanpatch(np.nanmin)
np.nanmax = nanpatch(np.nanmax)



