from sage.all import *
from sage.version import version as sage_version
import argparse, json, os, time

# Independent reproduction of M17.  This deliberately does not use the
# precomputed Manin table or the Numba evaluator used by M17.
AINVS = [1, 0, 1, -3, 1]  # 105a1
Q = 11117
P = 11
ELL1 = 463
ELL2 = 947
ETA1 = 3
ETA2 = 2
N = ELL1 * ELL2
NSHARDS = 128


def discrete_logs(ell, eta):
    out = [-1] * ell
    x = 1
    for k in range(ell - 1):
        out[x] = k
        x = (x * eta) % ell
    assert x == 1 and all(v >= 0 for v in out[1:])
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("shard", type=int)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    if not 0 <= args.shard < NSHARDS:
        raise ValueError("shard must lie in 0..127")

    l1 = discrete_logs(ELL1, ETA1)
    l2 = discrete_logs(ELL2, ETA2)
    lo = 1 + (N - 1) * args.shard // NSHARDS
    hi = (N - 1) * (args.shard + 1) // NSHARDS
    weights = []
    nonzero = 0
    for a in range(lo, hi + 1):
        r1, r2 = a % ELL1, a % ELL2
        w = 0 if (r1 == 0 or r2 == 0) else (l1[r1] * l2[r2]) % P
        weights.append(w)
        nonzero += int(w != 0)

    # Direct rational modular-symbol evaluation in PARI.  No table lookup,
    # no rational reconstruction and no BSD-predicted Sha value enter here.
    E = pari.ellinit(AINVS)
    Z = pari.msfromell(E, 1)
    weights_gp = "[" + ",".join(map(str, weights)) + "]"
    pari(
        f"M={Z[0]};S=mseval(M,{Z[1]});q={Q};p={P};n={N};"
        f"lo={lo};hi={hi};w={weights_gp};"
        "part()=sum(a=lo,hi,my(ww=w[a-lo+1]);"
        "if(ww,ww*sum(u=1,q-1,kronecker(q,u)*"
        "Mod(mseval(M,S,[oo,a/n+u/q]),p)),0))"
    )
    started = time.time()
    residue = pari("part()")
    elapsed = time.time() - started

    out = {
        "case": "105a1_q11117_p11",
        "engine": "direct PARI msfromell/mseval",
        "shard": args.shard,
        "nshards": NSHARDS,
        "lo": lo,
        "hi": hi,
        "p": P,
        "q": Q,
        "ell1": ELL1,
        "ell2": ELL2,
        "eta1": ETA1,
        "eta2": ETA2,
        "n": N,
        "nonzero_weights": nonzero,
        "residue": int(residue.lift()),
        "residue_pari": str(residue),
        "elapsed_s": elapsed,
        "pari_version": str(pari("version()")),
        "sage_version": sage_version,
    }
    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, f"p11_direct_shard_{args.shard:03d}.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()
