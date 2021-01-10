import numpy as np
from numpy.random import choice, binomial


def dict_quantile(cdf, q):
    if q > max(cdf.values()):
        q = max(cdf.values())
    key = 0
    val = 2
    for k in cdf.keys():
        if cdf[k] >= q:
            if cdf[k] <= val:
                key = k
                val = cdf[k]
    return key


def dict_quantiles(cdf, qs):
    quantiles = []
    for q in sorted(qs):
        quantiles.append(dict_quantile(cdf, q))
    return quantiles


def randomize_treatment_complete(Yl, Nlt):
    T = np.zeros(len(Yl))
    T[choice(range(len(Yl)), Nlt, replace=False)] = 1
    return np.array(Yl)[T == 1], np.array(Yl)[T == 0]


def randomize_treatment_bernoulli(Yl, el):
    T = np.array([binomial(1, el) for i in xrange(len(Yl))])
    return np.array(Yl)[T == 1], np.array(Yl)[T == 0]


def find_value(chi_min, chi_max, joint):
    vbest = 10000000
    kbest = (0, 0)
    for k in joint.keys():
        if abs(k[0] - chi_min) + abs(k[1] - chi_max) < vbest:
            vbest = abs(k[0] - chi_min) + abs(k[1] - chi_max)
            kbest = k
    return kbest, joint[kbest]


def gen_strata(X):
    strata = [() for i in xrange(0, len(X))]
    for i in xrange(len(X)):
        for j in xrange(len(X)):
            if np.array_equal(X[i], X[j]):
                strata[i] += tuple([j])
    return list(set(strata))


def coarsen(x, eps):
    ar = np.arange(min(x), max(x) + eps, eps)
    vals = [(ar[i], ar[i + 1]) for i in xrange(len(ar) - 1)]
    x_coarse = [() for i in xrange(len(x))]
    x_ind = [0 for i in xrange(len(x))]
    for i in xrange(len(x)):
        for j in xrange(len(vals)):
            if x[i] >= vals[j][0] and x[i] < vals[j][1]:
                x_coarse[i] = vals[j]
                x_ind[i] = j
        if x[i] >= vals[-1][1]:
            x_coarse[i] = (vals[-1][1], vals[-1][1] + eps)
            x_ind[i] = len(vals)
    return x_coarse, x_ind


def coarsen_matrix(X, eps):
    coarse_X = np.zeros(X.shape)
    for j in xrange(X.shape[1]):
        coarse_X[:, j] = coarsen(X[:, j], eps[j])[1]
    return coarse_X


def split_vectors(S, Y, T):
    Ylts = [[] for l in xrange(len(S))]
    Ylcs = [[] for l in xrange(len(S))]
    for l in xrange(len(S)):
        for i in S[l]:
            if T[i] == 0:
                Ylcs[l].append(Y[i])
            else:
                Ylts[l].append(Y[i])
    return Ylts, Ylcs


def compute_hsharp_stats(S, Ylts, Ylcs, remove_zeros=False):
    # quantities needed for null
    Nl1 = [0 for l in xrange(len(S))]
    Nl0 = [0 for l in xrange(len(S))]
    el = [0 for l in xrange(len(S))]

    for l in xrange(len(S)):
        Nl1[l] = int(sum(Ylts[l]) + sum(Ylcs[l]))
        Nl0[l] = int(len(Ylts[l]) - sum(Ylts[l]) + len(Ylcs[l]) - sum(Ylcs[l]))
        el[l] = float(len(Ylts[l])) / float(Nl1[l] + Nl0[l])

    Nl1_s = []
    Nl0_s = []
    el_s = []
    if remove_zeros:
        for l in xrange(len(S)):
            if(Nl1[l] >= 1 and Nl0[l] >= 1 and el[l] > 0 and el[l] < 1):
                Nl1_s.append(Nl1[l])
                Nl0_s.append(Nl0[l])
                el_s.append(el[l])
    else:
        Nl1_s = Nl1
        Nl0_s = Nl0
        el_s = el
    return Nl1_s, Nl0_s, el_s


def compute_hate_stats(S, Ylts, Ylcs, remove_zeros=True):
    Nlt = [0 for l in xrange(len(S))]
    Nlc = [0 for l in xrange(len(S))]
    pl = [0 for l in xrange(len(S))]

    for l in xrange(len(S)):
        Nlt[l] = len(Ylts[l])
        Nlc[l] = len(Ylcs[l])
        pl[l] = float(sum(Ylcs[l]) + sum(Ylts[l])) / \
            (Nlt[l] + Nlc[l]) if (Nlt[l] + Nlc[l]) > 0 else 0

    if remove_zeros:
        Nlt_s = []
        Nlc_s = []
        pl_s = []
        for l in xrange(len(S)):
            if(Nlt[l] >= 1 and Nlc[l] >= 1 and pl[l] > 0 and pl[l] < 1):
                Nlt_s.append(Nlt[l])
                Nlc_s.append(Nlc[l])
                pl_s.append(pl[l])

    return Nlt_s, Nlc_s, pl_s


def prune_excess(Ylt, Ylc, maximize=True):
    if maximize:
        if len(Ylc) > len(Ylt):
            return Ylt, np.sort(Ylc)[0:(len(Ylt))]
        elif len(Ylt) == len(Ylc):
            return Ylt, Ylc
        else:
            return np.sort(Ylt)[::-1][0:(len(Ylc))], Ylc
    else:
        v2, v1 = prune_excess(Ylc, Ylt, True)
        return v1, v2


def compute_max_mcnemar(Ylts, Ylcs):
    N = 0
    TE = 0.0
    SD = 0.0
    for l in xrange(len(Ylts)):
        Yltp, Ylcp = prune_excess(Ylts[l], Ylcs[l], True)
        N += len(Yltp)
        TE += (np.sum(Yltp) - np.sum(Ylcp))
    if TE > 0:
        for l in xrange(len(Ylts)):
            Yltp, Ylcp = prune_excess(Ylts[l], Ylcs[l], True)
            SD += abs(np.sum(Yltp) - np.sum(Ylcp))
    else:
        for l in xrange(len(Ylts)):
            Yltp, Ylcp = prune_excess(Ylts[l], Ylcs[l], True)
            SD += len(Ylcp) - abs(np.sum(Yltp) + np.sum(Ylcp) - len(Ylcp))

    return (TE - 1) / np.sqrt(SD + 1)


def compute_min_mcnemar(Ylts, Ylcs):
    N = 0
    TE = 0.0
    SD = 0.0
    for l in xrange(len(Ylts)):
        Yltp, Ylcp = prune_excess(Ylts[l], Ylcs[l], False)
        N += len(Yltp)
        TE += (np.sum(Yltp) - np.sum(Ylcp))
    if TE <= 0:
        for l in xrange(len(Ylts)):
            Yltp, Ylcp = prune_excess(Ylts[l], Ylcs[l], False)
            SD += abs(np.sum(Yltp) - np.sum(Ylcp))
    else:
        for l in xrange(len(Ylts)):
            Yltp, Ylcp = prune_excess(Ylts[l], Ylcs[l], False)
            SD += len(Ylcp) - abs(np.sum(Yltp) + np.sum(Ylcp) - len(Ylcp))
    if SD < TE: print TE, SD
    return (TE - 1) / np.sqrt(SD + 1)


def qq(pdf1, pdf2, quants=np.arange(0, 1.01, 0.01)):
    cdf1 = {k: sum([pdf1[v] for v in pdf1.keys() if v <= k])
            for k in pdf2.keys()}
    cdf2 = {k: sum([pdf2[v] for v in pdf2.keys() if v <= k])
            for k in pdf2.keys()}

    q1 = dict_quantiles(cdf1, quants)
    q2 = dict_quantiles(cdf2, quants)

    return q1, q2


def KL_div(pdf1, pdf2):
    KL = 0
    for k in pdf2.keys():
        if k in pdf1.keys():
            KL -= pdf2[k] * np.log(pdf1[k] / pdf2[k])
    return KL






