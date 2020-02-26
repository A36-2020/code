import pytest
from Moment_of_Inertia import *
from mainF100 import *

def f_to_test_Izz():
    return abs((crosssection.Izz-Izz_total)/crosssection.Izz*100)

def f_to_test_Iyy():
    return abs((crosssection.Iyy-Iyy_total)/crosssection.Iyy*100)

def test_Izz():
    assert f_to_test_Izz()<0.55

def test_Iyy():
    assert f_to_test_Iyy()<0.25

test_Izz()
test_Iyy()
