from sage.all import *
from sage.version import version as sage_version
import hashlib, json, os, random

CASES = {
    "62a1_p7": {"ainvs": [1,-1,1,-1,1], "p": 7},
    "141c1_p7": {"ainvs": [1,0,0,-2,3], "p": 7},
    "105a1_p11": {"ainvs": [1,0,1,-3,1], "p": 11},
    "66b1_p31": {"ainvs": [1,1,1,-2,-1], "p": 31},
}


def mod_p_of_rational(x, p):
    x = QQ(x)
    den = int(x.denominator()) % p
    if den == 0:
        raise ZeroDivisionError((p, x))
    return (int(x.numerator()) % p) * inverse_mod(den, p) % p


def cusp_expr(num, den):
    return "oo" if den == 0 else str(QQ(num) / QQ(den))


def elementary_path(c, d):
    g, x, y = xgcd(ZZ(c), ZZ(d))
    if g == 0:
        raise ValueError("zero pair")
    return cusp_expr(-x, d // g), cusp_expr(y, c // g)


def eval_table(num, den, N, p, table, loverp):
    ans = table[1]
    c, d, a, b = 0, 1, int(num), int(den)
    while b:
        f = b
        b = a % b
        quotient = (a - b) // f
        a = -f
        e = d
        d = -c
        c = quotient * c + e
        value = table[(c % N) * N + (d % N)]
        if value is None:
            raise AssertionError((c, d))
        ans = (ans + value) % p
    return (ans + loverp) % p


def export_case(name, cfg, outdir):
    ainvs, p = cfg["ainvs"], cfg["p"]
    E = EllipticCurve(ainvs)
    N = int(E.conductor())
    z = pari.msfromell(pari.ellinit(ainvs), 1)
    M, xp = z[0], z[1]
    pari("MM=%s;XX=%s" % (M, xp))
    table = [None] * (N * N)
    valid = 0
    for c in range(N):
        for d in range(N):
            if gcd(gcd(c, d), N) != 1:
                continue
            alpha, beta = elementary_path(c, d)
            value = pari(f"mseval(MM,XX,[{alpha},{beta}])")
            table[c*N+d] = mod_p_of_rational(value, p)
            valid += 1
    loverp_q = pari("mseval(MM,XX,[oo,0])")
    loverp = mod_p_of_rational(loverp_q, p)
    rng = random.Random(20260713 + N + p)
    checks = [(0,1),(1,1),(-1,1),(1,2),(-3,7),(17,20329)]
    while len(checks) < 506:
        den = rng.randrange(1, 2_000_000)
        num = rng.randrange(-2_000_000, 2_000_001)
        g = gcd(num, den)
        checks.append((num//g, den//g))
    validation = []
    for num, den in checks:
        r = QQ(num)/QQ(den)
        direct = mod_p_of_rational(pari(f"mseval(MM,XX,[oo,{r}])"), p)
        rebuilt = eval_table(num, den, N, p, table, loverp)
        if direct != rebuilt:
            raise AssertionError((name,num,den,direct,rebuilt))
        validation.append([num,den,direct])
    payload = {
        "case": name, "ainvs": ainvs, "conductor": N, "p": p,
        "valid_projective_pairs": valid, "loverp_rational": str(loverp_q),
        "loverp_mod_p": int(loverp), "table_mod_p": table,
        "validation_count": len(validation),
        "validation_sha256": hashlib.sha256(json.dumps(validation,separators=(",",":")).encode()).hexdigest(),
        "pari_version": str(pari("version()")), "sage_version": sage_version,
    }
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name+"_manin_table.json")
    with open(path,"w",encoding="utf-8") as fh:
        json.dump(payload,fh,indent=2,sort_keys=True)
    print(json.dumps({k:payload[k] for k in ["case","conductor","p","valid_projective_pairs","validation_count","validation_sha256"]},sort_keys=True))


if __name__ == "__main__":
    for name,cfg in CASES.items():
        export_case(name,cfg,"research/bsd_m20/results")
