import numpy as np
import matplot.lib.pyplot as plt
import F100

def visualize(moments, shears, flows_l, flows_t):
    """ Takes in lists of internal forces, over the aileron and displays shit """

    moment = sum(moments)
    shear = sum(shears)
    q_l = sum(flows_l)
    q_t = sum(flows_t)

    n_sec = moment.shape[0]


