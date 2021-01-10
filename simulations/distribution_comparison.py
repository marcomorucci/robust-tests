import numpy as np
from McNemar_conditional_distribution import compute_joint, compute_margins
from McNemar_sharp_distribution import compute_joint_sharp
from utils import *
import matplotlib.pyplot as plt
np.random.seed(42069)


N = 200
L = 25
pl1 = np.random.beta(5, 2, L)
pl0 = np.random.beta(0.2, 2, L)
S = np.random.choice(range(L), N)
Y1 = [np.random.binomial(1, pl1[S[i]]) for i in xrange(N)]
Y0 = [np.random.binomial(1, pl0[S[i]]) for i in xrange(N)]
T = [np.random.binomial(1, 0.2) for i in xrange(N)]
Nlt = [0 for i in xrange(L)]
Nlc = [0 for i in xrange(L)]
Ylts = [[] for i in xrange(L)]
Ylcs = [[] for i in xrange(L)]
strata = [[] for i in xrange(L)]
for i in xrange(N):
    Nlt[S[i]] += T[i]
    Nlc[S[i]] += 1 - T[i]
    strata[S[i]].append(i)
    if T[i] == 1:
        Ylts[S[i]].append(Y1[i])
    else:
        Ylcs[S[i]].append(Y0[i])

Nlt, Nlc, pl_hat = compute_hate_stats(strata, Ylts, Ylcs, True)
joint_avg_te = compute_joint(Nlt, Nlc, pl_hat)
chim_avg_te, chip_avg_te = compute_margins(joint_avg_te)

pl0 = pl1
Y1 = [np.random.binomial(1, pl1[S[i]]) for i in xrange(N)]
Y0 = [np.random.binomial(1, pl0[S[i]]) for i in xrange(N)]
Nlt = [0 for i in xrange(L)]
Nlc = [0 for i in xrange(L)]
Ylts = [[] for i in xrange(L)]
Ylcs = [[] for i in xrange(L)]
strata = [[] for i in xrange(L)]
for i in xrange(N):
    Nlt[S[i]] += T[i]
    Nlc[S[i]] += 1 - T[i]
    strata[S[i]].append(i)
    if T[i] == 1:
        Ylts[S[i]].append(Y1[i])
    else:
        Ylcs[S[i]].append(Y0[i])

Nlt, Nlc, pl_hat = compute_hate_stats(strata, Ylts, Ylcs, True)
joint_avg_note = compute_joint(Nlt, Nlc, pl_hat)
chim_avg_note, chip_avg_note = compute_margins(joint_avg_note)

fig = plt.figure(1)
fig.set_figwidth(15)
fig.set_figheight(15)
plt.subplot(2, 2, 1)
plt.hist(chim_avg_te.keys(), bins=len(chim_avg_te.keys()), weights=chim_avg_te.values(),
         align="mid", alpha=0.5, color="b")
plt.xlim(-5, 5)
plt.xlabel(r'$\chi^-$', fontsize=20)
plt.title(r"Marginal distributions of $\chi^-$ with SATE > 0.")
plt.subplot(2, 2, 2)
plt.hist(chim_avg_note.keys(), bins=len(chim_avg_note.keys()), weights=chim_avg_note.values(),
         align="mid", alpha=0.5, color="r")
plt.xlim(-5, 5)
plt.xlabel(r'$\chi^-$', fontsize=20)
plt.title(r"Marginal distributions of $\chi^-$ with SATE = 0.")
plt.subplot(2, 2, 3)
plt.hist(chip_avg_te.keys(), bins=len(chip_avg_te.keys()), weights=chip_avg_te.values(),
         align="mid", alpha=0.5, color="b")
plt.xlabel(r'$\chi^+$', fontsize=20)
plt.title(r"Marginal distributions of $\chi^+$ with SATE > 0.")
plt.xlim(-5, 5)
plt.subplot(2, 2, 4)
plt.hist(chip_avg_note.keys(), bins=len(chip_avg_note.keys()), weights=chip_avg_note.values(),
         align="mid", alpha=0.5, color="r")
plt.xlabel(r'$\chi^+$', fontsize=20)
plt.title(r"Marginal distributions of $\chi^+$ with SATE = 0.")
plt.xlim(-5, 5)
plt.savefig("HSATE_null.png", dpi=300, bbox_inches='tight')


Y1 = [np.random.binomial(1, pl1[S[i]]) for i in xrange(N)]
Y0 = Y1
Nlt = [0 for i in xrange(L)]
Nlc = [0 for i in xrange(L)]
Ylts = [[] for i in xrange(L)]
Ylcs = [[] for i in xrange(L)]
strata = [[] for i in xrange(L)]
for i in xrange(N):
    Nlt[S[i]] += T[i]
    Nlc[S[i]] += 1 - T[i]
    strata[S[i]].append(i)
    if T[i] == 1:
        Ylts[S[i]].append(Y1[i])
    else:
        Ylcs[S[i]].append(Y0[i])

Nl1, Nl0, el_hat = compute_hsharp_stats(strata, Ylts, Ylcs, True)
joint_shp = compute_joint_sharp(Nl1, Nl0, el_hat)
chim_shp, chip_shp = compute_margins(joint_shp)

fig = plt.figure(2)
fig.set_figwidth(15)
fig.set_figheight(7)
plt.subplot(1, 2, 1)
plt.hist(chim_shp.keys(), bins=len(chim_shp.keys()), weights=chim_shp.values(),
         align="mid", alpha=0.5, color="b")
plt.xlim(-5, 5)
plt.xlabel(r'$\chi^-$', fontsize=20)
plt.title(r"Marginal distribution of $\chi^-$ under $H_0^{Sharp}$.")
plt.subplot(1, 2, 2)
plt.hist(chip_shp.keys(), bins=len(chip_shp.keys()), weights=chip_shp.values(),
         align="mid", alpha=0.5, color="b")
plt.xlim(-5, 5)
plt.xlabel(r'$\chi^+$', fontsize=20)
plt.title(r"Marginal distribution of $\chi^+$ under $H_0^{Sharp}$.")
plt.savefig("Hsharp_null.png", dpi=300, bbox_inches='tight')






