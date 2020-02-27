import pytest
from Centroid import *
from mainF100 import *

def f_to_test_centroid():
    return (zcentroidtotal+crosssection.zc)/(-crosssection.zc) * 100

def test_centroid():
    assert f_to_test_centroid()<0.2

test_centroid()
