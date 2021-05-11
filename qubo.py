import numpy as np
import dimod
import math
import os
import sys
import copy

from pyqubo import Binary
from dimod.generators.constraints import combinations
import dwave
import dimod
from dwave.system import LeapHybridSampler
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

class Wrapper(object):
    # TODO: Implement the best encoding from: https://arxiv.org/abs/1706.01945
    # in order to reduce the magnitude of the coefficients involved
    def __init__(self, name, bound = None):
    # Constructor
        self.name = name
        # Do we have an upper-bound for the value?
        # Use the log-trick from https://arxiv.org/abs/1302.5843
        if bound:
            M = math.floor(math.log2(bound))
            self.coefs = [2**n for n in range(M)]
            self.coefs.append(bound + 1 - 2**M)
            self.bits = [Binary(name + "-bit=" + str(index)) for index in range(M + 1)]

    def size(self):
    # Number of qubits in the wrapper
        return len(self.bits)

    def __mul__(self, other):
    # e.g. self * 2
        return other * self.expr()

    def __rmul__(self, other):
    # e.g. 2 * self
        return other * self.expr()

    def __rsub__(self, other):
    # This only works if the variable is a binary one
        assert len(self.bits) == 1
        return other - self.expr()

    def __str__(self):
    # Return the name
        return self.name

    def expr(self):
    # Return the bounded binary-encoding
        return sum(self.coefs[index] * self.bits[index] for index in range(len(self.coefs)))

    def sumOfBits(self):
    # Return the sum of bits: needed to test whether the represented number if zero
        return sum(self.bits[index] for index in range(len(self.coefs)))

    def valueOfBits(self, best):
        return sum(best[self.name + "-bit=" + str(index)] for index in range(len(self.coefs)))

    def value(self, best):
    # Return the value of the solution
        return sum(self.coefs[index] * best[self.name + "-bit=" + str(index)] for index in range(len(self.coefs)))

def solveQUBO(tmax, pmax, mmax, q, kmax, B, I, sfo, ssr, u, x, Cmin, Cmax, cc, dfo, dsr, hc, lc, mc, udc, rev, xi, delta):
    maxValue = delta
    numQubits = 0
    for p in range(pmax):
        # B[(p, 0)] has already been initialized
        for t in range(1, 1 + tmax):
            B[(p, t)] = Wrapper('B_%i_%i' % (p, t), bound=maxValue)
            numQubits += B[(p, t)].size()
    for p in range(pmax):
        # I[(p, 0)] has already been initialized
        for t in range(1, 1 + tmax):
            I[(p, t)] = Wrapper('I_%i_%i' % (p, t), bound=maxValue)
            numQubits += I[(p, t)].size()
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            ssr[(p, t)] = Wrapper('ssr_%i_%i' % (p, t), bound=dsr[(p, t)])
            numQubits += ssr[(p, t)].size()
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            sfo[(p, t)] = Wrapper('sfo_%i_%i' % (p, t), bound=maxValue)
            numQubits += sfo[(p, t)].size()
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                u[(p, m, t)] = Wrapper('u_%i_%i_%i' % (p, m, t), bound=1)
                numQubits += u[(p, m, t)].size()
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                x[(p, m, t)] = Wrapper('x_%i_%i_%i' % (p, m, t), bound=maxValue)
                numQubits += x[(p, m, t)].size()

    # Init the expressions of the constraints
    constraint2 = 0
    constraint3 = 0
    constraint4 = 0
    constraint5_left_side = 0
    constraint5_right_side = 0
    constraint6 = 0

    def buildSideOfConstraint(*elems):
        side = 0
        for elem in elems:
            side += elem * 1
        return side

    # Equality constraint modelling
    def modelEqConstraint(a, b):
        return (a - b)**2

    # And model the constraints
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            # Compute the sums used in (3.2) (see below)
            sum1 = 0
            sum2 = 0
            for m in range(mmax):
                sum1 += x[(p, m ,t)] * 1
                sum2 += xi[(p, m, t)]
            # (3.2)
            constraint2 += modelEqConstraint(buildSideOfConstraint(I[(p, t)], ssr[(p, t)], sfo[(p, t)]), buildSideOfConstraint(I[(p, t - 1)], sum1, sum2))

            # (3.3)
            constraint3 += modelEqConstraint(buildSideOfConstraint(sfo[(p, t)], B[(p, t)]), buildSideOfConstraint(dfo[(p, t)], B[(p, t - 1)]))

            # (3.4): already taken into consideration through the bounded binary-encoding
    # (3.5)
    if True:
        # Init slack variables
        l, r = {}, {}
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                l[(m, t)] =  Wrapper('l_%i_%i' % (m, t), bound=Cmax[(m, t)])
                r[(m, t)] =  Wrapper('r_%i_%i' % (m, t), bound=Cmax[(m, t)])
                numQubits += l[(m, t)].size()
                numQubits += r[(m, t)].size()

        # And build the constraint
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                variable_sum = 0
                const_sum = 0
                for p in range(pmax):
                    for k in range(1 + min(kmax, tmax - t)):
                        # Optimize the loop by splitting the variable part from the constant one
                        const_sum += cc[(p, m, k)] * xi[(p, m, t + k)]
                        variable_sum += cc[(p, m, k)] * x[(p, m, t + k)]
                constraint5_left_side += modelEqConstraint(Cmin[(m, t)] + l[(m, t)] * 1, const_sum + variable_sum * 1)
                constraint5_right_side += modelEqConstraint(const_sum + variable_sum * 1, Cmin[(m, t)] + r[(m, t)] * 1)

    # (3.6)
    if False:
        slack = {}
        for p in range(pmax):
            for m in range(mmax):
                for t in range(1, 1 + tmax):
                    slack[(p, m, t)] = Wrapper('slack_%i_%i' % (m, t), bound=maxValue)
                    constraint6 += modelEqConstraint(x[(p, m, t)] + slack[(p, m, t)] * 1, maxValue * u[(p, m, t)])

    # Build the objective function
    objFunc = 0
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            objFunc += rev[(p, t)] * ssr[(p, t)]
            objFunc -= hc[(p, t)] * I[(p, t)]
            objFunc -= udc[(p, t)] * B[(p, t)]
            sum1 = 0
            sum2 = 0
            for m in range(mmax):
                sum1 += mc[(p, m, t)] * x[(p, m, t)]
                sum2 += lc[(p, m, t)] * u[(p, m, t)]
            objFunc -= sum1
            objFunc -= sum2

    # Build the Hamiltonian
    H = 0

    # We want to maximize, so revert the sign
    H -= objFunc

    # Add the constraints multiplied by their strengths
    strengths = {constraint2 : 1, constraint3 : 1, constraint5_left_side : 1, constraint5_right_side : 1, constraint6 : 1}
    for c in strengths.keys():
        H += strengths[c] * c

    # And compile
    bqm = H.compile().to_bqm()

    def showSolution(best_solution):
        if False:
            return
        for p in range(pmax):
            for m in range(mmax):
                for t in range(1, 1 + tmax):
                    print('u_%i_%i_%i' % (p, m, t) + " -> " + str(u[(p, m, t)].value(best_solution)))
        for p in range(pmax):
            for m in range(mmax):
                for t in range(1, 1 + tmax):
                    print('x_%i_%i_%i' % (p, m, t) + " -> " + str(x[(p, m, t)].value(best_solution)))
                    print("sum of bits of: " + 'x_%i_%i_%i' % (p, m, t) + " -> " + str(x[(p, m, t)].valueOfBits(best_solution)))
        if True:
            for m in range(mmax):
                for t in range(1, 1 + tmax):
                    sum = 0
                    for p in range(pmax):
                        inner_sum = 0
                        for k in range(1 + min(kmax, tmax - t)):
                            inner_sum += cc[(p, m, k)] * (x[(p, m, t + k)].value(best_solution) + xi[(p, m, t + k)])
                        sum += inner_sum
                    print(str(Cmin[(m, t)]) + " <= " + str(sum) + " <= " + str(Cmax[(m, t)]))
                    print(sum)
    sampler = LeapHybridSampler()
    print("Sending problem...")
    sample_set = sampler.sample(bqm)
    print("Results from D-Wave:")
    if False:
        print(sample_set)
    best_solution = sample_set.lowest().first.sample
    if False:
        for item in best_solution.items():
            print(item)
            label, value = item
            print(label)
            print(value)

    Q, offset = bqm.to_qubo()
    print("#qubits=" + str(numQubits))
    print("Solution Energy: ", dimod.qubo_energy(best_solution, Q, offset))
    showSolution(best_solution)
    pass

def main():
    solveQUBO()
    pass

if __name__ == '__main__':
    main()
