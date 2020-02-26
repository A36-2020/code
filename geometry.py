import F100
import numpy as np
from math import *
from Moment_of_Inertia import *

def Inertia_xx():
    ## Half circle in the front
    I_c = pi*(F100.h/2)**3
    ## spar in the middle
    I_s = F100.tsp * F100.h**3 / 12

    ## lentgh of the skin at the trailing edge
    l = sqrt((F100.h**2)/4+(F100.Ca-F100.h/2)**2)

    ## skin thickness assume all stringers are squished around the skin, like jello
    A_str = (F100.hst+F100.wst)*F100.tst
    A_total = F100.nst*A_str

    t = F100.tsk+A_total/(2*l+F100.h*pi)

    ## Moment of inertia of the plate at angle
    phi = asin(F100.h/2/l)
    I_sk = l**3*t*sin(phi)**2/12
    I_par = l*t*(F100.h/4)**2

    return Moment_of_inertia()[0]#I_c*t+I_s+2*(I_sk+I_par)

def Inertia_l():
    ## Half circle in the front
    I_c = pi*F100.tsk*(F100.h/2)**3
    ## spar in the middle
    I_s = F100.tsp * F100.h**3 / 12
    l = sqrt((F100.h**2)/4+(F100.Ca-F100.h/2)**2)
    ## skin thickness assume all stringers are squished around the skin, like jello
    A_str = (F100.hst+F100.wst)*F100.tst
    A_total = F100.nst*A_str

    t = F100.tsk+A_total/(2*l+F100.h*pi)
    return I_c*t+I_s


def Inertia_t():
    ## spar in the middle
    I_s = F100.tsp * F100.h**3 / 12

    ## lentgh of the skin at the trailing edge
    l = sqrt((F100.h**2)/4+(F100.Ca-F100.h/2)**2)

    ## skin thickness assume all stringers are squished around the skin, like jello
    A_str = (F100.hst+F100.wst)*F100.tst
    A_total = F100.nst*A_str

    t = F100.tsk+A_total/(2*l+F100.h*pi)

    ## Moment of inertia of the plate at angle
    phi = asin(F100.h/2/l)
    I_sk = l**3*t*sin(phi)**2/12
    I_par = l*t*(F100.h/4)**2

    return I_s+2*(I_sk+I_par)

def A_l():
    return 0.5*pi*(F100.h/2)**2

def A_t():
    return 0.5*F100.h*(F100.Ca-F100.h/2)

def S_al():
    return F100.h*F100.tsp + pi*F100.h*F100.tsk

def S_at():
    l = sqrt((F100.h**2)/4+(F100.Ca-F100.h/2)**2)
    return F100.h*F100.tsp+2*l*F100.tst

def S_l():
    return F100.h+ pi*F100.h/2

def S_t():
    l = sqrt((F100.h**2)/4+(F100.Ca-F100.h/2)**2)
    return F100.h+2*l

def l():
    return sqrt((F100.h**2)/4+(F100.Ca-F100.h/2)**2)

