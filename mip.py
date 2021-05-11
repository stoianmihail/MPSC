from ortools.linear_solver import pywraplp
import numpy as np
import math

def solveMIP():
    tmax = 3
    pmax = 3
    mmax = 3
    q = 3
    kmax = q - 1

    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')

    infinity = solver.infinity()
    B = {}
    for p in range(pmax):
        B[(p, 0)] = 250 * mmax / pmax
        for t in range(1, 1 + tmax):
            B[(p, t)] = solver.IntVar(0.0, infinity, 'B_%i_%i' % (p, t))
    I = {}
    for p in range(pmax):
        I[(p, 0)] = 500 * mmax / pmax
        for t in range(1, 1 + tmax):
            I[(p, t)] = solver.IntVar(0.0, infinity, 'I_%i_%i' % (p, t))
    sfo = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            sfo[(p, t)] = solver.IntVar(0.0, infinity, 'sfo_%i_%i' % (p, t))
    ssr = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            ssr[(p, t)] = solver.IntVar(0.0, infinity, 'ssr_%i_%i' % (p, t))
    u = {}
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                u[(p, m, t)] = solver.IntVar(0.0, 1.0, 'u_%i_%i_%i' % (p, m, t))
    x = {}
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                x[(p, m, t)] = solver.IntVar(0.0, infinity, 'x_%i_%i_%i' % (p, m, t))
    Cmin = {}
    Cmax = {}
    for m in range(mmax):
        for t in range(1, 1 + tmax):
            # 6720 hours
            Cmax[(m, t)] = 6720
            Cmin[(m, t)] = 0.1 * Cmax[(m, t)]
    cc = {}
    for p in range(pmax):
        for m in range(mmax):
            for k in range(1 + kmax):
                # 2 hours / wafer
                cc[(p, m, k)] = 2
    dfo = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            dfo[(p, t)] = 200 * mmax
    dsr = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            dsr[(p, t)] = 200 * mmax
    hc = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            hc[(p, t)] = 5
    lc = {}
    mc = {}
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                mc[(p, m, t)] = 10
                lc[(p, m, t)] = 375
    udc = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            udc[(p, t)] = 200
    rev = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            rev[(p, t)] = 80
    xi = {}
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                # Consider only in-house locations
                xi[(p, m, t)] = 400 / pmax
    delta = 0
    for m in range(mmax):
        for t in range(1, 1 + tmax):
            currMin = math.inf
            for p in range(pmax):
                sum = 0
                for k in range(min(kmax, tmax - t)):
                    sum += cc[(p, m, k)]
                if sum:
                    currMin = min(currMin, sum)
            if currMin != math.inf:
                delta = max(delta, Cmax[(m, t)] / currMin)
    # delta = max(delta, Cmax[(m, t)] / min([sum([cc[(p, m, k)] for k in range(min(kmax, tmax - t))]) for p in range(pmax)]))
    print(delta)

    for p in range(pmax):
        for t in range(1, 1 + tmax):
            sum1 = 0
            sum2 = 0
            for m in range(mmax):
                sum1 += x[(p, m ,t)]
                sum2 += xi[(p, m, t)]

                # (3.6)
                solver.Add(x[(p, m, t)] <= delta * u[(p, m, t)])
            # (3.2)
            solver.Add(I[(p, t)] + sfo[(p, t)] + ssr[(p, t)] == I[(p, t)] + sum1 + sum2)

            # (3.3)
            solver.Add(sfo[(p, t)] + B[(p, t)] == dfo[(p, t)]  + B[(p, t - 1)])

            # (3.4)
            solver.Add(ssr[(p, t)] <= dsr[(p, t)])

    # (3.5)
    if True:
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                sum = 0
                for p in range(pmax):
                    inner_sum = 0
                    for k in range(1 + min(kmax, tmax - t)):
                        inner_sum += cc[(p, m, k)] * (x[(p, m, t + k)] + xi[(p, m, t + k)])
                    sum += inner_sum
                solver.Add(Cmin[(m, t)] <= sum)
                solver.Add(sum <= Cmax[(m, t)])

    # Build the objective function
    expr = 0
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            expr += rev[(p, t)] * ssr[(p, t)]
            expr -= hc[(p, t)] * I[(p, t)]
            expr -= udc[(p, t)] * B[(p, t)]
            sum1 = 0
            sum2 = 0
            for m in range(mmax):
                sum1 += mc[(p, m, t)] * x[(p, m, t)]
                sum2 += lc[(p, m, t)] * u[(p, m, t)]
            expr -= sum1
            expr -= sum2
    solver.Maximize(expr)
    status = solver.Solve()

    def showSolution():
        if True:
            return
        for p in range(pmax):
            B[(p, 0)] = 250 * mmax / pmax
            for t in range(1, 1 + tmax):
                print('B_%i_%i' % (p, t) + " -> " + str(B[(p, t)].solution_value()))
        for p in range(pmax):
            for m in range(mmax):
                for t in range(1, 1 + tmax):
                    print('u_%i_%i_%i' % (p, m, t) + " -> " + str(u[(p, m, t)].solution_value()))
        for p in range(pmax):
            for m in range(mmax):
                for t in range(1, 1 + tmax):
                    print('x_%i_%i_%i' % (p, m, t) + " -> " + str(x[(p, m, t)].solution_value()))
        if False:
            for m in range(mmax):
                for t in range(1, 1 + tmax):
                    sum = 0
                    for p in range(pmax):
                        inner_sum = 0
                        for k in range(1 + min(kmax, tmax - t)):
                            inner_sum += cc[(p, m, k)] * (x[(p, m, t + k)].solution_value() + xi[(p, m, t + k)])
                        sum += inner_sum
                    print(str(Cmin[(m, t)]) + " <= " + str(sum) + " <= " + str(Cmax[(m, t)]))
                    print(sum)
    if status == pywraplp.Solver.OPTIMAL:
        print('Solution:')
        print('Objective value =', solver.Objective().Value())
        showSolution()
    else:
        if status == solver.FEASIBLE:
            print('A potentially suboptimal solution was found.')
        else:
            print('The solver could not solve the problem.')
    print('\nAdvanced usage:')
    print('Problem solved in %f milliseconds' % solver.wall_time())
    print('Problem solved in %d iterations' % solver.iterations())
    print('Problem solved in %d branch-and-bound nodes' % solver.nodes())
    pass

def main():
    solveMIP()
    pass

if __name__ == '__main__':
    main()
