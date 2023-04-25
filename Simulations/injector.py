import math
import pint




def ox_injector_rad(Cd, p, dP, B, mDot, g):
    den = Cd * math.pi * math.sqrt(2 * g * p * dP) * (B + 1)
    num = mDot * B
    return math.sqrt(num / den)

def fu_injector_rad(mDot, Cd, g, p, dP, B, Cd_o, p_o, dP_o, w_oi):
    term1 = mDot / (Cd * math.sqrt(2 * g * p * dP) * (B + 1))
    term2 = (ox_injector_rad(Cd_o, p_o, dP_o, B, mDot, g) + w_oi) ** 2 * math.pi
    return math.sqrt((term1 + term2) / math.pi) 

def recess_len(w_oi, Cd, p, dP, B, mDot, g):
    return 4 * (w_oi + ox_injector_rad(Cd, p, dP, B, mDot, g))


def fu_massflow(B, mDot):
    return mDot / (B + 1)

def ox_massflow(B, mDot):
    return (B * mDot) / (B + 1)



