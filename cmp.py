import sys
import mip
import qubo

def cmp(type=0):
    if type == 0:
        print("MIP vs QUBO..")
    elif type == 1:
        print("Run MIP only..")
    elif type == 2:
        print("Run QUBO only..")

    tmax = 2
    pmax = 2
    mmax = 2
    q = 2
    kmax = q - 1

    # Set the constants
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
    rev = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            rev[(p, t)] = 80
    dsr = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            dsr[(p, t)] = 200 * mmax
    hc = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            hc[(p, t)] = 5
    udc = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            udc[(p, t)] = 200
    lc = {}
    mc = {}
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                mc[(p, m, t)] = 10
                lc[(p, m, t)] = 375
    xi = {}
    for p in range(pmax):
        for m in range(mmax):
            for t in range(1, 1 + tmax):
                # Consider only in-house locations
                xi[(p, m, t)] = 400 / pmax
    dfo = {}
    for p in range(pmax):
        for t in range(1, 1 + tmax):
            dfo[(p, t)] = 200 * mmax

    # Init the variables
    B = {}
    for p in range(pmax):
        B[(p, 0)] = 250 * mmax / pmax
    I = {}
    for p in range(pmax):
        I[(p, 0)] = 500 * mmax / pmax
    ssr = {}
    sfo = {}
    u = {}
    x = {}

    # And compute `delta` := max(delta, Cmax[(m, t)] / min([sum([cc[(p, m, k)] for k in range(min(kmax, tmax - t))]) for p in range(pmax)]))
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
    print(delta)

    # And run
    if type == 0 or type == 1:
        mip.solveMIP(tmax=tmax, pmax=pmax, mmax=mmax, q=q, kmax=kmax, B=B, I=I, sfo=sfo, ssr=ssr, u=u, x=x, Cmin=Cmin, Cmax=Cmax, cc=cc, dfo=dfo, dsr=dsr, hc=hc, lc=lc, mc=mc, udc=udc, rev=rev, xi=xi, delta=delta)
    if type == 0 or type == 2:
        qubo.solveQUBO(tmax=tmax, pmax=pmax, mmax=mmax, q=q, kmax=kmax, B=B, I=I, sfo=sfo, ssr=ssr, u=u, x=x, Cmin=Cmin, Cmax=Cmax, cc=cc, dfo=dfo, dsr=dsr, hc=hc, lc=lc, mc=mc, udc=udc, rev=rev, xi=xi, delta=delta)
    pass

def main():
    # Run both the MIP and the QUBO
    if len(sys.argv) == 1:
        cmp()
    elif len(sys.argv) == 2:
        if sys.argv[1] == "mip":
            cmp(type=1)
        elif sys.argv[1] == "qubo":
            cmp(type=2)
        else:
            print("Option " + sys.argv[1] + " does not exist!")
            sys.exit(-1)
    else:
        sys.exit(-1)
    pass

if __name__ == '__main__':
    main()
