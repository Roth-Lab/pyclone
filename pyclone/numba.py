try:
    from numba import *

except ImportError:
    def identity_decorator(*args, **kwargs):
        if len(args) == 1 and hasattr(args[0], '__call__'):
            return args[0]

        else:
            def _f(f):
                return f
            return _f

    jit = identity_decorator

    jitclass = identity_decorator

    int64 = None

    float64 = None
