import numpy as np
from scipy.stats import binom
from sparse_array import SparsePmf, SparseSingle
import sys


def u_plus_pmf(x, Nlt, Nlc, pl):
    if x > Nlt:
        return 0.0
    elif x == Nlc:
        return 1 - binom.cdf(Nlc - 1, Nlt, pl)
    elif x >= 0 and x < Nlc:
        return binom.pmf(x, Nlt, pl)
    else:
        return 0.0


def u_minus_pmf(x, Nlt, Nlc, pl):
    tdiff = max(Nlt - Nlc, 0.0)
    if x > Nlt:
        return 0.0
    elif x == 0:
        return binom.cdf(tdiff, Nlt, pl)
    elif x > 0 and x <= Nlt - tdiff:
        return binom.pmf(x + tdiff, Nlt, pl)
    else:
        return 0.0


def u_plus_minus_pmf(x, y, Nlt, Nlc, pl):
    tdiff = max(Nlt - Nlc, 0.0)
    if x == y + tdiff and x >= 0 and x < Nlc and y > 0 and y <= Nlt - tdiff:
        return binom.pmf(x, Nlt, pl)
    elif y == 0 and x >= 0 and x <= tdiff and x < Nlc:
        return binom.pmf(x, Nlt, pl)
    elif x == Nlc and y + tdiff >= Nlc and y > 0:
        return binom.pmf(y + tdiff, Nlt, pl)
    elif x == Nlc and y == 0 and tdiff >= Nlc:
        return binom.cdf(tdiff, Nlt, pl) - binom.cdf(Nlc - 1, Nlt, pl)
    else:
        return 0.0


def eta_plus_pmf(x, Nlt, Nlc, pl):
    return u_minus_pmf(x, Nlc, Nlt, pl)


def eta_minus_pmf(x, Nlt, Nlc, pl):
    return u_plus_pmf(x, Nlc, Nlt, pl)


def eta_plus_minus_pmf(x, y, Nlt, Nlc, pl):
    return u_plus_minus_pmf(y, x, Nlc, Nlt, pl)


def g1(a, b, c, d, Nlt, Nlc, plt, plc):
    if(abs(a) != b):
        return 0.0

    Ml = min(Nlt, Nlc)
    p = 0.0
    e1 = 0.5 * (2 * Ml - d - c)
    e2 = 0.5 * (d - c)
    for j in xrange(0, Ml + 1):
        # print "j = %d, e1 = %d, e2 = %d" % (j, e1, e2)
        add = u_plus_minus_pmf(e1 + c, a + j, Nlt, Nlc, plt) * \
            eta_plus_minus_pmf(e1, j, Nlt, Nlc, plc)
        # print "p1 = ", add
        if(e1 != e2):
            add += u_plus_minus_pmf(e2 + c, a + j, Nlt, Nlc, plt) * \
                eta_plus_minus_pmf(e2, j, Nlt, Nlc, plc)
        # print "p2 = ", add
        p += add

    # print "p =", p
    return p


def g2(a, b, c, d, Nlt, Nlc, plt, plc):
    if(abs(a) != b or abs(c) != d):
        return 0.0

    Ml = min(Nlt, Nlc)
    p = 0.0

    for k in xrange(0, Ml + 1):
        for j in xrange(0, Ml + 1):
            p += u_plus_minus_pmf(c + j, a + k, Nlt, Nlc, plt) * \
                eta_plus_minus_pmf(j, k, Nlt, Nlc, plc)
            # print "k = %d, j = %d, p = %f" % (k, j, p)
    # print p
    return p


def g3(a, b, c, d, Nlt, Nlc, plt, plc):
    if(abs(c) != d):
        return 0.0

    Ml = min(Nlt, Nlc)
    p = 0.0
    e1 = 0.5 * (2 * Ml - b - a)
    e2 = 0.5 * (b - a)

    for j in xrange(0, Ml + 1):
        add = u_plus_minus_pmf(c + j, e1 + a, Nlt, Nlc, plt) * \
            eta_plus_minus_pmf(j, e1, Nlt, Nlc, plc)

        if (e1 != e2):
            add += u_plus_minus_pmf(c + j, e2 + a, Nlt, Nlc, plt) *\
                eta_plus_minus_pmf(j, e2, Nlt, Nlc, plc)
        # print "j = %d, e1 = %d, e2 = %d, add = %f" % (j, e1, e2, add)
        p += add
    # print "p = ", p
    return p


def create_lookup(Nlt, Nlc, plt, plc, g):
    Ml = min(Nlt, Nlc)
    tab = np.array([[[[g(a, b, c, d, Nlt, Nlc, plt, plc) for d in xrange(0, Ml + 1)]
                      for c in xrange(-Ml, Ml + 1)] for b in xrange(0, Ml + 1)]
                    for a in xrange(-Ml, Ml + 1)])
    if np.sum(tab > 0) == 1:
        return SparseSingle()
    else:
        return SparsePmf(tab, [-Ml, 0, -Ml, 0], [Ml, Ml, Ml, Ml])


def convolve_many(pmfs):
    out_shp = np.ones(len(pmfs[0].shape))
    out_shp += sum([np.array(a.shape) - 1 for a in pmfs])
    out_shp = tuple(int(x) for x in out_shp)
    conv_fft = np.ones(out_shp)

    cnt = 0
    for a in pmfs:
        cnt += 1
        sys.stdout.write("%.1f %% done \n" % ((float(cnt) / len(pmfs)) * 100))
        sys.stdout.flush()
        fft = np.fft.fftn(a.to_np(), s=out_shp)
        conv_fft = np.prod([fft, conv_fft], axis=0)
        del fft

    conv = np.fft.ifftn(conv_fft)
    lo = (np.array(np.shape(conv))) / 2
    lo[1] = 0
    lo[3] = 0
    cv = {tuple(i - lo): conv.real[tuple(i)]
          for i in np.squeeze(np.dstack(np.where(conv.real > 1e-100)))}
    return cv


def create_joint_Z(f1, f2, f3):
    sys.stdout.write("Computing joint distribution... ")
    joint = {}
    tot_vals = len(f1.keys()) + len(f2.keys()) + len(f3.keys())
    cnt = 0
    for vs in f1.keys():
        if vs[0] <= 1 and vs[2] <= 1:
            z_m = (vs[0] - 1) / np.sqrt(vs[1] + 1)
            z_p = (vs[2] - 1) / np.sqrt(vs[3] + 1)
            if (z_m, z_p) in joint:
                joint[(z_m, z_p)] += f1[vs]
            else:
                joint[(z_m, z_p)] = f1[vs]
        cnt += 1

    for vs in f2.keys():
        if vs[0] <= 1 and vs[2] > 1:
            z_m = (vs[0] - 1) / np.sqrt(vs[1] + 1)
            z_p = (vs[2] - 1) / np.sqrt(vs[3] + 1)
            if (z_m, z_p) in joint:
                joint[(z_m, z_p)] += f2[vs]
            else:
                joint[(z_m, z_p)] = f2[vs]
        cnt += 1

    for vs in f3.keys():
        if vs[0] > 1 and vs[2] > 1:
            z_m = (vs[0] - 1) / np.sqrt(vs[1] + 1)
            z_p = (vs[2] - 1) / np.sqrt(vs[3] + 1)
            if (z_m, z_p) in joint:
                joint[(z_m, z_p)] += f3[vs]
            else:
                joint[(z_m, z_p)] = f3[vs]
        cnt += 1

    return joint


def compute_margins(joint):
    chip = {}
    chim = {}
    for m, p in joint.keys():
        if m in chim:
            chim[m] += joint[(m, p)]
        else:
            chim[m] = joint[(m, p)]
        if p in chip:
            chip[p] += joint[(m, p)]
        else:
            chip[p] = joint[(m, p)]
    return chim, chip


def create_range_Z(joint):
    rd = {}
    for m, p in joint.keys():
        if p - m in rd:
            rd[p - m] += joint[(m, p)]
        else:
            rd[p - m] = joint[(m, p)]
    return rd


def compute_joint(Nlt, Nlc, pl):
    sys.stdout.write("Computing ATE distribution\n")
    pmfs1 = [create_lookup(Nlt[l], Nlc[l], pl[l], pl[l], g1)
             for l in xrange(len(Nlt))]
    pmfs2 = [create_lookup(Nlt[l], Nlc[l], pl[l], pl[l], g2)
             for l in xrange(len(Nlt))]
    pmfs3 = [create_lookup(Nlt[l], Nlc[l], pl[l], pl[l], g3)
             for l in xrange(len(Nlt))]
    sys.stdout.write("Convolving f1... ")
    f1 = convolve_many(pmfs1)
    sys.stdout.write("Convolving f2... ")
    f2 = convolve_many(pmfs2)
    sys.stdout.write("Convolving f3... ")
    f3 = convolve_many(pmfs3)
    return create_joint_Z(f1, f2, f3)


def joint_pvalue(value, joint):
    pass


def range_pvalue(value, range_dist):
    pass

if __name__ == "__main__":
    pass
