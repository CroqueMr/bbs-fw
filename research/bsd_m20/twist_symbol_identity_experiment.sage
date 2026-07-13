from sage.all import *
from sage.version import version as sage_version
import itertools, json, os, random, time

BASE_AINVS = [1,-1,1,-1,1]  # 62a1
TWIST_Q = 5
P = 7
SAMPLES = 30000


def modp(x):
    x = QQ(x)
    den = int(x.denominator()) % P
    if den == 0:
        raise ZeroDivisionError(x)
    return int(x.numerator()) * inverse_mod(den, P) % P


def endpoints(c, d):
    g, x, y = xgcd(ZZ(c), ZZ(d))
    if g == 0:
        raise ValueError((c,d))
    alpha = oo if d // g == 0 else QQ(-x) / QQ(d // g)
    beta = oo if c // g == 0 else QQ(y) / QQ(c // g)
    return alpha, beta


def path_value(symbol, c, d):
    alpha, beta = endpoints(c, d)
    # symbol(r) is {infinity,r}; hence {alpha,beta}=symbol(beta)-symbol(alpha).
    va = QQ(0) if alpha is oo else symbol(alpha)
    vb = QQ(0) if beta is oo else symbol(beta)
    return modp(vb - va)


def legendre_mod_q(x):
    x %= TWIST_Q
    if x == 0:
        return 0
    return 1 if kronecker_symbol(x, TWIST_Q) == 1 else P-1


def main():
    E = EllipticCurve(BASE_AINVS)
    T = E.quadratic_twist(TWIST_Q).global_minimal_model()
    N = int(E.conductor())
    NT = int(T.conductor())
    base = E.modular_symbol(sign=+1, implementation='eclib')
    twist = T.modular_symbol(sign=+1, implementation='eclib')
    invq = inverse_mod(TWIST_Q, N)

    transforms = {
        'c_d': lambda c,d:(c,d),
        'qc_d': lambda c,d:(TWIST_Q*c,d),
        'iqc_d': lambda c,d:(invq*c,d),
        'c_qd': lambda c,d:(c,TWIST_Q*d),
        'c_iqd': lambda c,d:(c,invq*d),
        'd_c': lambda c,d:(d,c),
        'qd_c': lambda c,d:(TWIST_Q*d,c),
        'd_qc': lambda c,d:(d,TWIST_Q*c),
        'c_dpc': lambda c,d:(c,d+c),
        'c_dmc': lambda c,d:(c,d-c),
        'cpd_d': lambda c,d:(c+d,d),
        'cmd_d': lambda c,d:(c-d,d),
    }
    chars = {
        '1': lambda c,d:1,
        'chi_c': lambda c,d:legendre_mod_q(c),
        'chi_d': lambda c,d:legendre_mod_q(d),
        'chi_cpd': lambda c,d:legendre_mod_q(c+d),
        'chi_cmd': lambda c,d:legendre_mod_q(c-d),
        'chi_cd': lambda c,d:(legendre_mod_q(c)*legendre_mod_q(d))%P,
    }

    rng = random.Random(20260713)
    pairs=[]
    while len(pairs)<SAMPLES:
        c=rng.randrange(0,NT); d=rng.randrange(0,NT)
        if gcd(gcd(c,d),NT)==1:
            pairs.append((c,d))
    # Add structured representatives for all local q-types.
    for c in range(0,5*TWIST_Q):
        for d in range(0,5*TWIST_Q):
            if gcd(gcd(c,d),NT)==1:
                pairs.append((c,d))

    started=time.time()
    rows=[]
    for c,d in pairs:
        vt=path_value(twist,c,d)
        bv={}
        for name,f in transforms.items():
            cc,dd=f(c,d)
            bv[name]=path_value(base,int(cc),int(dd))
        rows.append((c,d,vt,bv))

    scores=[]
    for tname in transforms:
        for cname,ch in chars.items():
            inferred=[]
            for c,d,vt,bv in rows:
                x=(ch(c,d)*bv[tname])%P
                if x:
                    inferred.append(vt*inverse_mod(x,P)%P)
            candidates=range(P)
            for scalar in candidates:
                ok=0
                for c,d,vt,bv in rows:
                    pred=scalar*ch(c,d)*bv[tname]%P
                    ok += int(pred==vt)
                scores.append({'transform':tname,'character':cname,'scalar':int(scalar),'matches':ok,'total':len(rows),'rate':ok/len(rows)})
    scores.sort(key=lambda x:(x['matches'],-x['scalar']),reverse=True)

    local={}
    for rc in range(TWIST_Q):
        for rd in range(TWIST_Q):
            subset=[r for r in rows if r[0]%TWIST_Q==rc and r[1]%TWIST_Q==rd]
            if subset:
                local[f'{rc},{rd}']={'count':len(subset),'twist_zero_rate':sum(r[2]==0 for r in subset)/len(subset)}

    out={
        'base_ainvs':BASE_AINVS,'base_conductor':N,'twist_q':TWIST_Q,
        'twist_ainvs':[int(x) for x in T.ainvs()],'twist_conductor':NT,
        'p':P,'samples':len(rows),'top_candidates':scores[:50],
        'best_exact':scores[0]['matches']==len(rows),'local_q_patterns':local,
        'elapsed_s':time.time()-started,'sage_version':sage_version,
    }
    os.makedirs('research/bsd_m20/results',exist_ok=True)
    open('research/bsd_m20/results/twist_symbol_identity_experiment.json','w').write(json.dumps(out,indent=2,sort_keys=True))
    print(json.dumps(out,sort_keys=True))


if __name__=='__main__':
    main()
