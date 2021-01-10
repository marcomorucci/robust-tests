# -*- coding: utf-8 -*-
# @Author: marco
# @Date:   2019-07-28 11:17:51
# @Last Modified by:   marco
# @Last Modified time: 2019-08-17 11:27:19

import numpy as np
from scipy.stats import norm
import pandas as pd
from utils import *

np.random.seed(42069)

sim_res = pd.DataFrame(columns=['dlt', "ATE", "L", "chi_max", "chi_min", "pval_avg"])
nsim = 100
N = 250
for dlt in [0.1, 0.5, 0.9]:
    print "dlt = ", dlt
    for sim in xrange(nsim):
        print "simulation", sim, 'of', nsim
        for te in [0, 1]:
            L = int(np.floor(N / 9))

            if te == 0:
                pl1 = np.random.beta(0.5, 0.5, L)
                pl0 = pl1
            else:
                pl1 = np.random.beta(0.5, 0.5, L)
                pl0 = np.random.beta(0.1, 0.5, L)
            
            el = np.random.uniform(0, 1, L)
            Sind = np.random.choice(range(L), N)
            S = [[] for l in xrange(L)]
            for i in xrange(N):
                S[Sind[i]].append(i)

            eps = [np.random.uniform(-min(pl0[Sind[i]], pl1[Sind[i]]),
                                     1 - max(pl0[Sind[i]], pl1[Sind[i]])) for i in xrange(N)]

            Y1 = [np.random.binomial(1, pl1[Sind[i]] + eps[i] * dlt) for i in xrange(N)]
            if (te != 0):
                Y0 = [np.random.binomial(1, pl0[Sind[i]] + eps[i] * dlt) for i in xrange(N)]
            else:
                Y0 = [np.random.binomial(1, pl1[Sind[i]] + eps[i] * dlt) for i in xrange(N)]

            T = [np.random.binomial(1, 0.2) for i in xrange(N)]

            Ylts = [[] for i in xrange(L)]
            Ylcs = [[] for i in xrange(L)]
            Yl = [[]for i in xrange(L)]

            for i in xrange(N):
                if T[i] == 1:
                    Ylts[Sind[i]].append(Y1[i])
                    Yl[Sind[i]].append(Y1[i])
                else:
                    Ylcs[Sind[i]].append(Y0[i])
                    Yl[Sind[i]].append(Y0[i])

            # exclude empty strata
            Ylts = [Ylts[i] for i in xrange(L) if len(Yl[i]) != 0]
            Ylcs = [Ylcs[i] for i in xrange(L) if len(Yl[i]) != 0]
            S = [S[l] for l in xrange(L) if len(Yl[l]) != 0]
            Yl = [Yl[i] for i in xrange(L) if len(Yl[i]) != 0]

            chi_max = compute_max_mcnemar(Ylts, Ylcs)
            chi_min = compute_min_mcnemar(Ylts, Ylcs)

            pval_avg = min(1, 2 * min(norm.cdf(chi_max, 0, 1), 1 - norm.cdf(chi_min, 0, 1)))
            
            # Match largest control first
            M = 0
            b = 0
            c = 0
            for l in xrange(len(S)):
                Ml = min(len(Ylts[l]), len(Ylcs[l]))
                if Ml == 0:
                    continue
                Yc = np.array(sorted(Ylcs[l], reverse=True)[0:Ml])
                Yt = np.random.choice(Ylts[l], Ml, replace=False)
                b += sum((Yt - Yc) == 1)
                c += sum((Yt - Yc) == -1)
                chi = (b - c - 1) / np.sqrt(b + c + 1)
                M += sum((Yt - Yc) == 1) + sum((Yt - Yc) == -1)
            max_p = 2 * min(norm.cdf(chi, 0, 1), 1 - norm.cdf(chi, 0, 1))
            # max_p = min(1, 2 * min(binom.cdf(b, M, .5), 1 - binom.cdf(b - 1, M, 0.5)))

            # Match smallest control first
            M = 0
            b = 0
            c = 0
            for l in xrange(len(S)):
                Ml = min(len(Ylts[l]), len(Ylcs[l]))
                if Ml == 0:
                    continue
                Yc = np.array(sorted(Ylcs[l], reverse=False)[0:Ml])
                Yt = np.random.choice(Ylts[l], Ml, replace=False)
                b += sum((Yt - Yc) == 1)
                M += sum((Yt - Yc) == 1) + sum((Yt - Yc) == -1)
                chi = (b - c - 1) / np.sqrt(b + c + 1)
            min_p = 2 * min(norm.cdf(chi, 0, 1), 1 - norm.cdf(chi, 0, 1))
            # min_p = min(1, 2 * min(binom.cdf(b, M, .5), 1 - binom.cdf(b - 1, M, 0.5)))

            # Match at random
            M = 0
            b = 0
            c = 0
            for l in xrange(len(S)):
                Ml = min(len(Ylts[l]), len(Ylcs[l]))
                if Ml == 0:
                    continue
                Yt = np.random.choice(Ylts[l], Ml, replace=False)
                Yc = np.random.choice(Ylcs[l], Ml, replace=False)
                b += sum((Yt - Yc) == 1)
                M += sum((Yt - Yc) == 1) + sum((Yt - Yc) == -1)
                chi = (b - c - 1) / np.sqrt(b + c + 1)
            rnd_p = 2 * min(norm.cdf(chi, 0, 1), 1 - norm.cdf(chi, 0, 1))
                                               
            sim_res = sim_res.append({"dlt": dlt, "ATE": np.mean(pl1 - pl0), "L": L,
                                      "chi_max": chi_max, "chi_min": chi_min,
                                      "pval_avg": pval_avg, 'min_p': min_p,
                                      "max_p": max_p, "rnd_p": rnd_p}, ignore_index=True)

sim_res.to_csv('confounded_method_comparison3.csv')
