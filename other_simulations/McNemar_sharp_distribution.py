import numpy as np
from scipy.stats import binom
from sparse_array import SparsePmf, SparseSingle
import sys

def inner_joint_lookup(Nl1, Nl0, el):
    Nl = Nl1 + Nl0
    pmf = {}
    for ul in xrange(Nl1 + 1):
        for vl in xrange(Nl0 + 1):
            etal = Nl1 - ul
            Nlt = ul + vl
            Nlc = Nl - Nlt
            Glp = max(Nlc - Nlt, 0)
            Glm = max(Nlt - Nlc, 0)
            ulp = ul - max(ul - Nlc, 0)
            ulm = max(ul - Glm, 0)
            etalp = max(etal - Glp, 0)
            etalm = etal - max(etal - Nlt, 0)
            Ml = min(Nlt, Nlc)
            # print "ul = %d, vl = %d, Nlt = %d, Nlc = %d, \
            #       ulp = %d, ulm = %d, elp = %d, elm = %d" % \
            #    (ul, vl, Nlt, Nlc, ulp, ulm, etalp, etalm)
            if((ulp, ulm, etalp, etalm, Ml) not in pmf):
                pmf[(ulp, ulm, etalp, etalm, Ml)] = (binom.pmf(ul, Nl1, el) *
                                                     binom.pmf(vl, Nl0, el))
            else:
                pmf[(ulp, ulm, etalp, etalm, Ml)] += (binom.pmf(ul, Nl1, el) *
                                                      binom.pmf(vl, Nl0, el))
    return pmf


def g1(a, b, c, d, Nl, pmfl):
    if(abs(a) != b):
        return 0.0
    p = 0
    for Ml in xrange(0, int(np.floor(Nl / 2)) + 1):
        e1 = 0.5 * (2 * Ml - d - c)
        e2 = 0.5 * (d - c)
        if np.any(np.abs((a, b, c, d)) > Ml):
            continue
        for j in xrange(0, Ml + 1):
            add = 0
            idx1 = (e1 + c, a + j, e1, j, Ml)
            idx2 = (e2 + c, a + j, e2, j, Ml)
            # print "e1 + c = %d, a +j = %d, e1 = %d, j= %d, Ml=%d" % idx1
            if idx1 in pmfl:
                add = pmfl[(e1 + c, a + j, e1, j, Ml)]
            # print "p1 = ", add
            # print "e2 + c = %d, a +j = %d, e2 = %d, j= %d, Ml=%d" % idx2
            if(e1 != e2) and idx2 in pmfl:
                add += pmfl[e2 + c, a + j, e2, j, Ml]
            # print "p2 = ", add
            p += add
    # print "p =", p
    return p


def g2(a, b, c, d, Nl, pmfl):
    if(abs(a) != b or abs(c) != d):
        return 0.0

    p = 0.0
    for Ml in xrange(0, int(np.floor(Nl / 2)) + 1):
        if np.any(np.abs((a, b, c, d)) > Ml):
            continue
        for k in xrange(0, Ml + 1):
            for j in xrange(0, Ml + 1):
                if (c + j, a + k, j, k, Ml) in pmfl:
                    p += pmfl[c + j, a + k, j, k, Ml]
            # print "k = %d, j = %d, p = %f" % (k, j, p)
    # print p
    return p


def g3(a, b, c, d, Nl, pmfl):
    if(abs(c) != d):
        return 0.0

    p = 0.0
    for Ml in xrange(0, int(np.floor(Nl / 2)) + 1):
        if np.any(np.abs((a, b, c, d)) > Ml):
            continue
        for j in xrange(0, Ml + 1):
            e1 = 0.5 * (2 * Ml - b - a)
            e2 = 0.5 * (b - a)
            add = 0
            idx1 = (c + j, e1 + a, j, e1, Ml)
            idx2 = (c + j, e2 + a, j, e2, Ml)
            if idx1 in pmfl:
                add = pmfl[idx1]
            if (e1 != e2) and idx2 in pmfl:
                add += pmfl[idx2]
            # print "j = %d, e1 = %d, e2 = %d, add = %f" % (j, e1, e2, add)
            p += add
        # print "p = ", p
    return p


def create_lookup(Nl1, Nl0, el, g):
    pmfl = inner_joint_lookup(Nl1, Nl0, el)
    Nl = Nl1 + Nl0
    Ml = min(Nl1, Nl0)  # int(np.floor(Nl / 2))
    if Ml == 0:
        print "Ml = 0! for ", Nl1, Nl0, el, g
    tab = np.array([[[[g(a, b, c, d, Nl, pmfl) for d in xrange(0, Ml + 1)]
                    for c in xrange(-Ml, Ml + 1)] for b in xrange(0, Ml + 1)]
                    for a in xrange(-Ml, Ml + 1)])

    if np.sum(tab > 0) == 1:
        return SparseSingle()
    else:
        return SparsePmf(tab, [-Ml, 0, -Ml, 0], [Ml, Ml, Ml, Ml])


def convolve_many(pmfs, lo):
    if len(pmfs) < 1:
        print pmfs
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
    cv = {tuple(i - lo): conv.real[tuple(i)]
          for i in np.squeeze(np.dstack(np.where(conv.real > 1e-100)))}
    return cv


def create_joint_Z(f1, f2, f3):
    joint = {}
    for vs in f1.keys():
        if vs[0] < 1 and vs[2] < 1:
            z_m = (vs[0] - 1) / np.sqrt(vs[1] + 1)
            z_p = (vs[2] - 1) / np.sqrt(vs[3] + 1)
            if (z_m, z_p) in joint:
                joint[(z_m, z_p)] += f1[vs]
            else:
                joint[(z_m, z_p)] = f1[vs]

    for vs in f2.keys():
        if vs[0] < 1 and vs[2] >= 1:
            z_m = (vs[0] - 1) / np.sqrt(vs[1] + 1)
            z_p = (vs[2] - 1) / np.sqrt(vs[3] + 1)
            if (z_m, z_p) in joint:
                joint[(z_m, z_p)] += f2[vs]
            else:
                joint[(z_m, z_p)] = f2[vs]

    for vs in f3.keys():
        if vs[0] >= 1 and vs[2] >= 1:
            z_m = (vs[0] - 1) / np.sqrt(vs[1] + 1)
            z_p = (vs[2] - 1) / np.sqrt(vs[3] + 1)
            if (z_m, z_p) in joint:
                joint[(z_m, z_p)] += f3[vs]
            else:
                joint[(z_m, z_p)] = f3[vs]

    return joint


def create_range_Z(joint):
    rd = {}
    for m, p in joint.keys():
        if p - m in rd:
            rd[p - m] += joint[(m, p)]
        else:
            rd[p - m] = joint[(m, p)]
    return rd


def compute_joint_sharp(Nl1, Nl0, el):
    sys.stdout.write("Computing sharp distribution\n")
    pmfs1 = [create_lookup(Nl1[l], Nl0[l], el[l], g1)
             for l in xrange(len(el))]
    pmfs2 = [create_lookup(Nl1[l], Nl0[l], el[l], g2)
             for l in xrange(len(el))]
    pmfs3 = [create_lookup(Nl1[l], Nl0[l], el[l], g3)
             for l in xrange(len(el))]
    Nl = sum([min(Nl1[l], Nl0[l]) for l in xrange(len(Nl1))])
    lo = np.array([Nl, 0, Nl, 0])
    sys.stdout.write("Convolving f1\n")
    f1 = convolve_many(pmfs1, lo)
    sys.stdout.write("Convolving f2\n")
    f2 = convolve_many(pmfs2, lo)
    sys.stdout.write("Convolving f3\n")
    f3 = convolve_many(pmfs3, lo)
    return create_joint_Z(f1, f2, f3)


def joint_pvalue(value, joint):
    pass


def range_pvalue(value, range_dist):
    pass

if __name__ == "__main__":
    pass
