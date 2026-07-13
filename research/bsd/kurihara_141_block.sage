from sage.all import *
from sage.version import version as sage_version
import argparse, concurrent.futures, json, multiprocessing, os, time

CASE = {"case":"141c1_q3461_p7", "ainvs":[1,0,0,-2,3], "q":3461, "p":7, "ell1":337, "ell2":743}
NSHARDS=128
SHARDS_PER_BLOCK=8


def dlogs(ell):
    g=int(primitive_root(ell)); table={}; x=1
    for k in range(ell-1):
        table[x]=k; x=x*g%ell
    return g,table


def one(shard):
    c=CASE; P=pari; P.allocatemem(128_000_000)
    z=P.msfromell(P.ellinit(c['ainvs']),1)
    q,p,e1,e2=c['q'],c['p'],c['ell1'],c['ell2']; n=e1*e2
    g1,l1=dlogs(e1); g2,l2=dlogs(e2)
    lo=1+(n-1)*shard//NSHARDS; hi=(n-1)*(shard+1)//NSHARDS
    w=[]; nz=0
    for a in range(lo,hi+1):
        r1=a%e1; r2=a%e2
        v=0 if (not r1 or not r2) else l1[r1]*l2[r2]%p
        w.append(v); nz+=bool(v)
    ws='['+','.join(map(str,w))+']'
    P(f"M={z[0]};S=mseval(M,{z[1]});q={q};p={p};lo={lo};hi={hi};w={ws};part()=sum(a=lo,hi,my(ww=w[a-lo+1]);if(ww,ww*sum(u=1,q-1,kronecker(q,u)*Mod(mseval(M,S,[oo,a/{n}+u/q]),p)),0))")
    started=time.time(); residue=P('part()')
    return {'case':c['case'],'shard':shard,'nshards':NSHARDS,'lo':lo,'hi':hi,'lift':int(residue.lift()),'result':str(residue),'elapsed_s':time.time()-started,'nonzero_weights':nz,'eta1':g1,'eta2':g2,'pari_version':str(P('version()'))}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('block',type=int); ap.add_argument('--out',required=True); args=ap.parse_args()
    if not 0<=args.block<16: raise ValueError('block')
    ids=list(range(args.block*SHARDS_PER_BLOCK,(args.block+1)*SHARDS_PER_BLOCK)); started=time.time()
    ctx=multiprocessing.get_context('spawn')
    with concurrent.futures.ProcessPoolExecutor(max_workers=4,mp_context=ctx) as ex:
        rows=list(ex.map(one,ids))
    c=CASE
    out={'case':c['case'],'block':args.block,'blocks':16,'shards':rows,'partial_mod_p':sum(r['lift'] for r in rows)%c['p'],'p':c['p'],'q':c['q'],'ell1':c['ell1'],'ell2':c['ell2'],'n':c['ell1']*c['ell2'],'elapsed_wall_s':time.time()-started,'sage_version':sage_version}
    os.makedirs(args.out,exist_ok=True); path=os.path.join(args.out,f"{c['case']}_block_{args.block:02d}.json")
    with open(path,'w',encoding='utf-8') as fh: json.dump(out,fh,indent=2,sort_keys=True)
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__': main()
