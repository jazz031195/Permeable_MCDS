import numpy

def test_numpy_array_creation():
    arr = numpy.array([1, 2, 3])
    assert arr.tolist() == [1, 2, 3]
    assert arr.shape == (3,)
