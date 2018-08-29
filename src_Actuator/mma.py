"""
Version September 2007 (and a small change August 2008)

Krister Svanberg <krille@math.kth.se>
Department of Mathematics KTH, SE-10044 Stockholm, Sweden.

Translated to python by A.J.J. Lagerweij TU Delft June 2018

"""
import numpy as np
from scipy.sparse import spdiags, csc_matrix

# MMA problem linearisation
def mma(m, n, itr, xval, xmin, xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c, d):
    '''
    This function mmasub performs one MMA-iteration, aimed at solving the
    nonlinear programming problem:

    Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
        subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m

                xmin_j <= x_j <= xmax_j,        j = 1,...,n

                z >= 0,   y_i >= 0,             i = 1,...,m
    Parameters
    _______
    m : int
        The number of general constraints.
    n : int
        The number of variables x_j
    itr : int
        Current iteration number ( =1 the first time mmasub is called).
    xval : 1-D array len(n)
        Vector with the current values of the variables x_j.
    xmin : 1-D array len(n)
        Vector with the lower bounds for the variables x_j.
    xmax : 1-D array len(n)
        Vector with the upper bounds for the variables x_j.
    xold1 : 1-D array len (n)
        xval, one iteration ago when iter>1, zero othewise.
    xold2 : 1-D array len (n)
        xval, two iteration ago when iter>2, zero othewise.
    f0val : float
        The value of the objective function f_0 at xval.
    df0dx : 1-D array len(n)
        Vector with the derivatives of the objective function f_0 with
        respect to the variables x_j, calculated at xval.
    fval : 1-D array len(m)
        Vector with the values of the constraint functions f_i,
        calculated at xval.
    dfdx : 2-D array size(m x n)
        (m x n)-matrix with the derivatives of the constraint functions f_i
        with respect to the variables x_j, calculated at xval
    low : 1-D array len(n)
        Vector with the lower asymptotes from the previous iteration
        (provided that iter>1).
    upp : 1-D array len(n)
        Vector with the upper asymptotes from the previous iteration
        (provided that iter>1).
    a0 : float
        The constants a_0 in the term a_0*z.
    a : 1-D array len(m)
        Vector with the constants a_i in the terms a_i*z.
    c : 1-D array len(m)
        Vector with the constants c_i in the terms c_i*y_i.
    d : 1-D array len(m)
        Vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

    Returns
    ------
    xmma : 1-D array len(n)
        Column vector with the optimal values of the variables x_j in the
        current MMA subproblem.
    low : 1-D array len(n)
        Column vector with the lower asymptotes, calculated and used in the
        current MMA subproblem
    upp : 1-D array len(n)
        Column vector with the upper asymptotes, calculated and used in the
        current MMA subproblem
    '''

    epsimin = 10**(-7)
    raa0 = 0.00001
    move = 1.0
    albefa = 0.1
    asyinit = 0.5
    asyincr = 1.2
    asydecr = 0.7
    eeen = np.ones((n))
    eeem = np.ones((m))
    zeron = np.zeros((n))

    # calculation of the upper and lower asymptotes
    if itr <= 2:
        low = xval - asyinit*(xmax-xmin)
        upp = xval + asyinit*(xmax-xmin)
    else:
        zzz = np.divide(xval-xold1, xold1-xold2, out=np.zeros_like(xval), where=xold1-xold2!=0)
        factor = eeen.copy()
        factor[np.where(zzz > 0)] = asyincr
        factor[np.where(zzz < 0)] = asydecr
        low = xval - factor*(xold1 - low)
        upp = xval + factor*(upp - xold1)
        lowmin = xval - 10*(xmax-xmin)
        lowmax = xval - 0.01*(xmax-xmin)
        uppmin = xval + 0.01*(xmax-xmin)
        uppmax = xval + 10*(xmax-xmin)
        low = np.minimum(lowmax, np.maximum(low, lowmin))
        upp = np.maximum(uppmin, np.minimum(upp, uppmax))

    # calculation of the bounds alfa and beta
    zzz1 = low + albefa*(xval - low)
    zzz2 = xval - move*(xmax - xmin)
    zzz = np.maximum(zzz1, zzz2)
    alfa = np.maximum(zzz, xmin)
    zzz1 = upp - albefa*(upp-xval)
    zzz2 = xval + move*(xmax - xmin)
    zzz = np.minimum(zzz1, zzz2)
    beta = np.minimum(zzz, xmax)

    # calculation of p0, q0
    xmami = xmax-xmin
    xmamieps = 0.00001*eeen
    xmami = np.maximum(xmami, xmamieps)
    xmamiinv = eeen/xmami

    ux1 = upp-xval
    ux2 = ux1*ux1
    uxinv = eeen/ux1
    xl1 = xval-low
    xl2 = xl1*xl1
    xlinv = eeen/xl1

    p0 = zeron
    q0 = zeron
    p0 = np.maximum(df0dx, 0)
    q0 = np.maximum(-df0dx, 0)
    pq0 = 0.001*(p0 + q0) + raa0*xmamiinv
    p0 = ((p0 + pq0)*ux2)
    q0 = ((q0 + pq0)*xl2)

    # calculation of P and Q and b
    P = np.maximum(dfdx, 0)
    Q = np.maximum(-dfdx, 0)
    PQ = 0.001*(P + Q) + raa0*(eeem*xmamiinv[:, np.newaxis]).T
    P = ((P + PQ)*spdiags(ux2, 0, n, n))
    Q = ((Q + PQ)*spdiags(xl2, 0, n, n))
    b = np.dot(P, uxinv) + np.dot(Q, xlinv) - fval

    # solving the simplified approximated problem
    xmma, ymma, zmma, lam, xsi, eta, mu, zet, s = solvemma(m, n, epsimin, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d)

    return xmma, low, upp


def solvemma(m, n, epsimin, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d):
    '''
    This function solves the MMA subproblem with a primal-dual Newton method:

    minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z + SUM[ ci*yi + 0.5*di*(yi)^2 ],

    subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
    alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
    '''
    epsi = 1
    een = np.ones((n))
    eem = np.ones((m))
    epsvecn = epsi*een
    epsvecm = epsi*eem
    x = 0.5*(alfa+beta)
    y = eem
    z = 1
    lam = eem
    xsi = een/(x-alfa)
    xsi = np.maximum(xsi, een)
    eta = een/(beta-x)
    eta = np.maximum(eta, een)
    mu = np.maximum(eem, 0.5*c)
    zet = 1
    s = eem
    itera = 0

    while epsi > epsimin:
        epsvecn = epsi*een
        epsvecm = epsi*eem
        ux1 = upp-x
        xl1 = x-low
        ux2 = ux1*ux1
        xl2 = xl1*xl1
        uxinv1 = een/ux1
        xlinv1 = een/xl1

        plam = p0 + np.dot(P.T, lam)
        qlam = q0 + np.dot(Q.T, lam)
        gvec = np.dot(P, uxinv1) + np.dot(Q, xlinv1)
        dpsidx = plam/ux2 - qlam/xl2
        rex = dpsidx - xsi + eta
        rey = c + d*y - mu - lam
        rez = a0 - zet - a*lam
        relam = gvec - a*z - y + s - b
        rexsi = xsi*(x-alfa) - epsvecn
        reeta = eta*(beta-x) - epsvecn
        remu = mu*y - epsvecm
        rezet = zet*z - epsi
        res = lam*s - epsvecm
        residu1 = np.hstack([rex.T, rey.T, rez])
        residu2 = np.hstack([relam.T, rexsi.T, reeta.T, remu.T, [rezet], res.T])
        residu = np.hstack([residu1, residu2])
        residunorm = np.sqrt(np.dot(residu, residu))
        residumax = np.max(np.abs(residu))
        ittt = 0
        while residumax > 0.9*epsi and ittt < 200:
            ittt = ittt + 1
            itera = itera + 1
            ux1 = upp-x
            xl1 = x-low
            ux2 = ux1*ux1
            xl2 = xl1*xl1
            ux3 = ux1*ux2
            xl3 = xl1*xl2
            uxinv1 = een/ux1
            xlinv1 = een/xl1
            uxinv2 = een/ux2
            xlinv2 = een/xl2
            plam = p0 + np.dot(P.T, lam)
            qlam = q0 + np.dot(Q.T, lam)
            gvec = np.dot(P, uxinv1) + np.dot(Q, xlinv1)
            GG = P*spdiags(uxinv2, 0, n, n) - Q*spdiags(xlinv2, 0, n, n)
            dpsidx = plam/ux2 - qlam/xl2
            delx = dpsidx - epsvecn/(x-alfa) + epsvecn/(beta-x)
            dely = c + d*y - lam - epsvecm/y
            delz = a0 - np.dot(a, lam) - epsi/z
            dellam = gvec - a*z - y - b + epsvecm/lam
            diagx = plam/ux3 + qlam/xl3
            diagx = 2*diagx + xsi/(x-alfa) + eta/(beta-x)
            diagxinv = een/diagx
            diagy = d + mu/y
            diagyinv = eem/diagy
            diaglam = s/lam
            diaglamyi = diaglam+diagyinv

            if m < n:
                blam = dellam + dely/diagy - np.dot(GG, delx/diagx)
                bb = np.hstack([blam, delz])
                Alam = spdiags(diaglamyi, 0, m, m) + np.dot(GG, spdiags(diagxinv, 0, n, n)*GG.T)
                AA = np.block([[Alam, a[np.newaxis].T], [a, (-zet/z)]])
                solut = np.linalg.solve(AA, bb)
                dlam = solut[0:m]
                dz = solut[m]
                dx = -delx/diagx - np.dot(GG.T, dlam)/diagx
            else:
                diaglamyiinv = eem/diaglamyi
                dellamyi = dellam + dely/diagy
                Axx = spdiags(diagx, 0, n, n) + np.dot(GG.T, spdiags(diaglamyiinv, 0, m, m)* GG)
                azz = zet/z + np.dot(a.T, (a/diaglamyi))
                axz = -np.dot(GG.T, (a/diaglamyi))
                bx = delx + np.dot(GG.T, (dellamyi/diaglamyi))
                bz = delz - np.dot(a.T, (dellamyi/diaglamyi))
                AA = np.block([[Axx, axz[np.newaxis].T], [axz, azz]])
                bb = np.hstack([-bx.T, -bz])
                solut = np.linalg.solve(AA, bb)
                dx = solut[0:n]
                dz = solut[n]
                dlam = np.dot(GG, dx)/diaglamyi - np.dot(dz, (a/diaglamyi)) + dellamyi/diaglamyi

            dy = -dely/diagy + dlam/diagy
            dxsi = -xsi + epsvecn/(x-alfa) - (xsi*dx)/(x-alfa)
            deta = -eta + epsvecn/(beta-x) + (eta*dx)/(beta-x)
            dmu = -mu + epsvecm/y - (mu*dy)/y
            dzet = -zet + epsi/z - zet*dz/z
            ds = -s + epsvecm/lam - (s*dlam)/lam
            xx = np.hstack([y.T, z, lam.T, xsi.T, eta.T, mu.T, zet, s.T])
            dxx = np.hstack([dy.T, dz, dlam.T, dxsi.T, deta.T, dmu.T, dzet, ds.T])

            stepxx = -1.01*dxx/xx  # check what is the correct formulation
            stmxx = np.max(stepxx)
            stepalfa = -1.01*dx/(x-alfa)
            stmalfa = np.max(stepalfa)
            stepbeta = 1.01*dx/(beta-x)
            stmbeta = np.max(stepbeta)
            stmalbe = np.maximum(stmalfa, stmbeta)
            stmalbexx = np.maximum(stmalbe, stmxx)
            stminv = np.maximum(stmalbexx, 1)
            steg = 1/stminv

            xold = x
            yold = y
            zold = z
            lamold = lam
            xsiold = xsi
            etaold = eta
            muold = mu
            zetold = zet
            sold = s

            itto = 0
            resinew = 2*residunorm
            while resinew > residunorm and itto < 50:
                itto = itto+1
                x = xold + steg*dx
                y = yold + steg*dy
                z = zold + steg*dz
                lam = lamold + steg*dlam
                xsi = xsiold + steg*dxsi
                eta = etaold + steg*deta
                mu = muold + steg*dmu
                zet = zetold + steg*dzet
                s = sold + steg*ds
                ux1 = upp - x
                xl1 = x - low
                ux2 = ux1*ux1
                xl2 = xl1*xl1
                uxinv1 = een/ux1
                xlinv1 = een/xl1
                plam = p0 + np.dot(P.T, lam)
                qlam = q0 + np.dot(Q.T, lam)
                gvec = np.dot(P, uxinv1) + np.dot(Q, xlinv1)
                dpsidx = plam/ux2 - qlam/xl2
                rex = dpsidx - xsi + eta
                rey = c + d*y - mu - lam
                rez = a0 - zet - np.dot(a.T, lam)
                relam = gvec - a*z - y + s - b
                rexsi = xsi*(x-alfa) - epsvecn
                reeta = eta*(beta-x) - epsvecn
                remu = mu*y - epsvecm
                rezet = np.dot(zet, z) - epsi
                res = lam*s - epsvecm
                residu1 = np.hstack([rex.T, rey.T, rez])
                residu2 = np.hstack([relam.T, rexsi.T, reeta.T, remu.T, rezet, res.T])
                residu = np.hstack([residu1, residu2])
                resinew = np.sqrt(np.dot(residu, residu))
                steg = steg/2

            residunorm = resinew
            residumax = np.max(np.abs(residu))
            steg = 2*steg

        if ittt >= 198:
            print('    MMA itteration runout')
            print('      ittt = ', ittt)
            print('      epsi = ', epsi)
        epsi = 0.1*epsi
    return x, y, z, lam, xsi, eta, mu, zet, s
