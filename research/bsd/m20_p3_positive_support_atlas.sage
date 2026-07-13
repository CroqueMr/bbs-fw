from sage.all import *
from sage.version import version as sage_version
import hashlib, json, os, random, time
import numpy as np
from numba import njit, prange, set_num_threads, get_num_threads
AINVS=[1,-1,1,-1,1];N=62;P=3
QS=[197,509,613,677,821,941,1277,1381,1997,2069,2213,2693,2749,2917,2957,3037,3301,3413,3613,3797,3989,4357,4493,4549,4733,4909]
EXPECTED_SHA={q:(81 if q in (4357,4493) else 9) for q in QS}
pari.allocatemem(768_000_000);Epari=pari.ellinit(AINVS);z=pari.msfromell(Epari,1);M,xp=z[0],z[1]
pari('M20_M=%s; M20_X=%s' % (M,xp));pari('M20_T(q)=sum(a=0,q-1,kronecker(q,a)*mseval(M20_M,M20_X,[oo,a/q]))')
def modp(x,p):
 x=QQ(x);d=int(x.denominator())%p
 if d==0:raise ZeroDivisionError((x,p))
 return int(x.numerator())%p*int(inverse_mod(d,p))%p
def cusp(n,d):return 'oo' if d==0 else str(QQ(n)/QQ(d))
def path(c,d):
 g,x,y=xgcd(ZZ(c),ZZ(d));return cusp(-x,d//g),cusp(y,c//g)
def evpy(num,den,table,loverp):
 ans=int(table[1]);c,d,a,b=0,1,int(num),int(den)
 while b:
  f=b;b=a%b;q=(a-b)//f;a=-f;e=d;d=-c;c=q*c+e;v=int(table[(c%N)*N+(d%N)])
  if v<0:raise AssertionError((c,d))
  ans=(ans+v)%P
 return (ans+loverp)%P
def table_build():
 table=np.full(N*N,-1,dtype=np.int16);valid=0
 for c in range(N):
  for d in range(N):
   if gcd(gcd(c,d),N)!=1:continue
   a,b=path(c,d);table[c*N+d]=modp(pari(f'mseval(M20_M,M20_X,[{a},{b}])'),P);valid+=1
 lq=pari('mseval(M20_M,M20_X,[oo,0])');loverp=modp(lq,P);rng=random.Random(20260713+N+P);rs=[(0,1),(1,1),(-1,1),(1,2),(-3,7),(17,20329)]
 while len(rs)<506:
  d=rng.randrange(1,2000000);n=rng.randrange(-2000000,2000001);g=gcd(n,d);rs.append((n//g,d//g))
 checks=[]
 for n,d in rs:
  r=QQ(n)/QQ(d);direct=modp(pari(f'mseval(M20_M,M20_X,[oo,{r}])'),P);recon=evpy(n,d,table,loverp)
  if direct!=recon:raise AssertionError((n,d,direct,recon))
  neg=modp(pari(f'mseval(M20_M,M20_X,[oo,{-r}])'),P)
  if direct!=neg:raise AssertionError(('sym',n,d,direct,neg))
  checks.append([n,d,direct])
 return table,loverp,{'valid_projective_pairs':valid,'loverp_rational':str(lq),'loverp_mod_p':int(loverp),'validation_count':len(checks),'validation_sha256':hashlib.sha256(json.dumps(checks,separators=(',',':')).encode()).hexdigest()}
def firstkp(q):
 out=[]
 for ell in prime_range(2,100001):
  if ell in (P,q) or N%ell==0 or q%ell==0 or (ell-1)%P:continue
  ab=ZZ(pari.ellap(Epari,ell));at=ZZ(kronecker(q,ell))*ab
  if (at-ell-1)%P==0:
   o=ZZ(ell+1-at);out.append({'ell':int(ell),'a_base':int(ab),'a_twist':int(at),'group_order':int(o),'v3_group_order':int(valuation(o,3))})
   if len(out)==2:return out
 raise RuntimeError((q,out))
def logs(ell,eta):
 a=np.full(ell,-1,dtype=np.int64);x=1
 for k in range(ell-1):a[x]=k;x=x*eta%ell
 if x!=1 or np.any(a[1:]<0):raise AssertionError((ell,eta))
 return a
def chis(q):
 a=np.zeros(q,dtype=np.int8)
 for u in range(1,q):a[u]=1 if pow(u,(q-1)//2,q)==1 else -1
 return a
@njit(cache=True)
def ev(num,den,N,p,table,loverp):
 ans=int(table[1]);c=0;d=1;a=num;b=den
 while b!=0:
  f=b;b=a%b;q=(a-b)//f;a=-f;e=d;d=-c;c=q*c+e;v=int(table[(c%N)*N+(d%N)])
  if v<0:return -999999
  ans=(ans+v)%p
 return (ans+loverp)%p
@njit(parallel=True,cache=True)
def kd(q,p,e1,e2,l1,l2,chi,N,table,loverp):
 n=e1*e2;den=n*q;vals=np.zeros(n,dtype=np.int64)
 for a in prange(1,n):
  r1=a%e1;r2=a%e2
  if r1==0 or r2==0:continue
  w=(l1[r1]*l2[r2])%p
  if w==0:continue
  inn=0
  for u in range(1,q):
   s=ev(a*q+u*n,den,N,p,table,loverp)
   if s<0:return vals
   inn+=int(chi[u])*s
  vals[a]=(w*(inn%p))%p
 return vals
def direct(q,e1,e2,eeta1,eeta2):
 l1={};x=1
 for k in range(e1-1):l1[x]=k;x=x*eeta1%e1
 l2={};x=1
 for k in range(e2-1):l2[x]=k;x=x*eeta2%e2
 n=e1*e2;acc=0
 for a in range(1,n):
  if a%e1==0 or a%e2==0:continue
  w=(l1[a%e1]*l2[a%e2])%P
  if not w:continue
  inn=0
  for u in range(1,q):inn+=int(kronecker(q,u))*modp(pari(f'mseval(M20_M,M20_X,[oo,{QQ(a,n)+QQ(u,q)}])'),P)
  acc=(acc+w*inn)%P
 return acc
def main():
 set_num_threads(min(24,os.cpu_count() or 1));st0=time.time();table,loverp,tcert=table_build();tsha=hashlib.sha256(table.tobytes()).hexdigest();targets=[]
 for q in QS:
  raw=QQ(pari(f'M20_T({q})'))
  if raw.denominator()!=1 or int(raw)!=EXPECTED_SHA[q]:raise AssertionError(('quotient',q,raw,EXPECTED_SHA[q]))
  targets.append({'q':q,'L_over_Omega_plus':str(raw),'sha_an':int(raw),'v3_sha_an':int(valuation(ZZ(raw),3)),'kp':firstkp(q)})
 _=kd(5,3,7,13,logs(7,primitive_root(7)),logs(13,primitive_root(13)),chis(5),N,table,int(loverp));rows=[]
 for t in targets:
  q=t['q'];e1=t['kp'][0]['ell'];e2=t['kp'][1]['ell'];eta1=int(primitive_root(e1));eta2=int(primitive_root(e2));l1=logs(e1,eta1);l2=logs(e2,eta2);chi=chis(q);st=time.time();vals=kd(q,P,e1,e2,l1,l2,chi,N,table,int(loverp));delta=int(np.sum(vals)%P);nz=sum(1 for a in range(1,e1*e2) if a%e1 and a%e2 and (int(l1[a%e1])*int(l2[a%e2]))%P);r={**t,'ell1':e1,'ell2':e2,'eta1':eta1,'eta2':eta2,'n':e1*e2,'delta_n_mod_p':delta,'nonzero':bool(delta),'nonzero_outer_weights':nz,'evaluations':nz*(q-1),'elapsed_s':time.time()-st};rows.append(r);print('RESULT',json.dumps(r,sort_keys=True),flush=True)
 c=next(r for r in rows if r['q']==677);dc=direct(677,c['ell1'],c['ell2'],c['eta1'],c['eta2'])
 if dc!=c['delta_n_mod_p']:raise AssertionError(('direct',dc,c['delta_n_mod_p']))
 s={'date':'2026-07-13','curve':'62a1','p':3,'q_bound':5000,'positive_support_count':len(rows),'nonzero_count':sum(r['nonzero'] for r in rows),'rows':rows,'table_certificate':tcert,'table_sha256':tsha,'direct_pari_control_q677':dc,'direct_pari_control_matches':True,'numba_threads':get_num_threads(),'sage_version':sage_version,'pari_version':str(pari('version()')),'elapsed_total_s':time.time()-st0};tmp=json.dumps(s,indent=2,sort_keys=True);s['certificate_sha256']=hashlib.sha256(tmp.encode()).hexdigest();os.makedirs('research/bsd/results',exist_ok=True);open('research/bsd/results/m20_p3_positive_support_atlas.json','w').write(json.dumps(s,indent=2,sort_keys=True));open('research/bsd/results/m20_p3_positive_support_atlas.csv','w').write('q,sha_an,v3_sha_an,ell1,ell2,n,delta_n_mod_p,nonzero,evaluations,elapsed_s\n'+''.join(f"{r['q']},{r['sha_an']},{r['v3_sha_an']},{r['ell1']},{r['ell2']},{r['n']},{r['delta_n_mod_p']},{int(r['nonzero'])},{r['evaluations']},{r['elapsed_s']}\n" for r in rows));print(json.dumps(s,sort_keys=True),flush=True)
if __name__=='__main__':main()
