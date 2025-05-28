from numba import njit
import math

@njit
def heaviside(x):
    return 1.0 if x >= 0 else 0.0

@njit
def nan_to_num(x):
    if math.isnan(x):
        return 0
    if math.isinf(x):
        return 1e200
    return x

@njit
def safe_divide(a, b, epsilon=1e-200):
    """0'a bölmeyi önlemek için güvenli bölme."""
    if abs(b) < epsilon:
        return 1e200
    return a / b

@njit
def A(x, t, ti, m, a, b, U, D, K):
    A1 = heaviside(t - ti) * heaviside(x)
    a1 = safe_divide(U * x, 2 * D)  
    a2 = m * math.sqrt(safe_divide(U**2 * x**2, 4 * D**2) + safe_divide(K * x**2, D))
    A2 = math.exp(a1 + a2)
    A3 = a * safe_divide((t - ti) / 2 + (m * x / (2 * math.sqrt(U**2 + 4 * D * K))), 1) + (a * ti + b) / 2
    a11 = x + m * (t - ti) * math.sqrt(U**2 + 4 * D * K)
    a12 = 2 * math.sqrt(D * abs(t - ti))
    A4 = math.erfc(nan_to_num(safe_divide(a11, a12)))  
    A = A1 * A2 * A3 * A4
    return nan_to_num(A)

@njit
def B(x, t, m, U, D, K, Co):
    B1 = heaviside(x) * Co / 2
    B2 = math.exp(-K * t + (1 + m) * safe_divide(U * x, 2 * D))  
    B3 = math.erfc(nan_to_num(safe_divide(x + m * U * t, 2 * math.sqrt(D * t))))  
    B = B1 * B2 * B3
    return nan_to_num(B)

@njit
def main(t, x, x_source, source_t, a, b, U, D, K, Co, diameter):

    
    ti = source_t[0]
    ti2 = source_t[1]
    
    A_sum = (A(x, t, ti, -1, a, b, U, D, K) + A(x, t, ti, 1, a, b, U, D, K) - 
             A(x, t, ti2, -1, a, b, U, D, K) - A(x, t, ti2, 1, a, b, U, D, K))
    
    B_sum = (2 * B(0, t, 0, U, D, K, Co) - B(x, t, -1, U, D, K, Co) - B(x, t, 1, U, D, K, Co))
    
    C = A_sum + B_sum
    return C

@njit
def const_ab(source_t, source_ws):
    n = len(source_t) - 1
    a_list = []
    b_list = []
    
    for i in range(n):
        t1, t2 = source_t[i], source_t[i+1]
        w1, w2 = source_ws[i], source_ws[i+1]
        
        a = safe_divide((w2 - w1), (t2 - t1))  
        b = w1 - a * t1
        
        a_list.append(a)
        b_list.append(b)

    return a_list[0], b_list[0]

@njit
def solver(t, U, x, x_source, Co, source_t, source_ws, diameter, K):
    a, b = const_ab(source_t, source_ws)

    #Chlorine
    Reynold = (U * diameter/ 1.26E-6)
    D = (1E-4 * (Reynold)**0.66)

    #Florine
    #Reynold = (U * diameter/ 1.0533E-6)
    #D = (2E-5 * (Reynold)**0.83)
    
    result = main(t, x, x_source, source_t, a, b, U, D, K, Co, diameter)
    return result
