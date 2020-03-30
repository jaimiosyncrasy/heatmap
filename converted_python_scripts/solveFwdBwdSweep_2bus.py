import numpy as np
import cmath


def solveFwdBwdSweep_2bus(R12, X12, V1, P2, Q2):
    # Initialization
    print('~~~~~~~ Starting FBS Method for Solving PF')

    # Givens: z12, V1, P2, Q2
    S2 = complex(P2, Q2)  # per unit
    z12 = complex(R12, X12)  # per unit
    Vs = V1.copy()  # per unit

    # Init Cond
    V2 = []
    Vconv = []
    Vnom = Vs.copy()  # to check convergence
    tol = 0.0001
    k = 0
    V1[k] = 0
    V2.append(0)
    Vconv.append([0, 0])

    '''First Iteration'''
    k += 1

    # Fwd Sweep
    V1[k] = Vs
    V2.append(Vs)

    # Check convergence:
    Vconv.append([abs((abs(V1[k]) - abs(V1[k - 1]))) / Vnom, \
                  abs((abs(V2[k]) - abs(V2[k - 1]))) / Vnom])
    # Backward sweep
    I12 = np.conj(S2 / V2[k])

    '''Iterative Part'''
    while any(node >= tol for node in Vconv[k]):  # break when all nodes less than tol
        k += 1  # new iteration
        # Fwd sweep
        V1[k] = V1[k - 1]  # same as prev iter ZERO?
        V2.append(Vs - (z12 * I12))

        # Check convergence:
        Vconv.append([abs((abs(V1[k]) - abs(V1[k - 1]))) / Vnom, \
                      abs((abs(V2[k]) - abs(V2[k - 1]))) / Vnom])

        # print(Vconv) uncomment when debugging

        # Backward sweep
        I12 = np.conj(S2 / V2[k])

        if len(Vconv) > 30:
            print('Didnt converge')
            break  # break out of loop of iter too many times

    '''Output Results'''
    print('~~~~~~~ PF Results: ')
    Vsoln = [V1[-1], V2[-1]]  # /Vs to put into pu
    convergedIfZero = Vconv[-1]
    numIter = len(Vconv) - 1  # -1 because Vconv initialized at zero
    print('~~~~~~~ Finished FBS Method for SOlving PF')

    '''Polar to rect conversion for testing/probing'''
    mag = [abs(ele) for ele in Vsoln]
    ang = [np.degrees(cmath.phase(ele)) for ele in Vsoln]
    Vsoln_polarCoor = [mag, ang]  # Vpu, deg

    V2 = abs(Vsoln[1])
    del2 = np.degrees(cmath.phase(Vsoln[1]))

    return V2, del2
