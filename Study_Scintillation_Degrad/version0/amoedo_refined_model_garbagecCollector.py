import numpy as np
from amoedo_refined_model_simple import Pgamma_Ar3rd_refined,Pgamma_CF4_refined,Pgamma_CF3_refined
# ------------------------------------------------------------------------------

def Pgamma_vis_refined_GC(f_cf4, Pgamma_CF3_star_dir, P_Ar_dblstar,autoQ,C):
    a = Pgamma_CF3_refined(f_cf4, Pgamma_CF3_star_dir, P_Ar_dblstar,autoQ)
    return C*a

def Pgamma_UV_refined_GC(f_cf4, Pgamma_CF4_plus_star_dir, P_Ar_3rd, n,kcool,C):
    a = Pgamma_Ar3rd_refined(f_cf4, P_Ar_3rd, n)
    b = Pgamma_CF4_refined(f_cf4, Pgamma_CF4_plus_star_dir, P_Ar_3rd, n,kcool)
    return C*(a+b)