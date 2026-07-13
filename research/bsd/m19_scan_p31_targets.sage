from sage.all import *
from sage.version import version as sage_version
import csv,hashlib,json,os,random,time
import numpy as np
from numba import njit
AINVS=[1,1,1,-2,-1];N=66;P=31;QMAX=300000;MAXPOS=30;KBOUND=300000
pari.allocatemem(768_000_000);E=pari.ellinit(AINVS);Z=pari.msfromell(E,1);M,xp=Z[0],Z[1];pari(f'M19_M={M};M19_X={xp};M19_T(q)=sum(a=0,q-1,kronecker(q,a)*mseval(M19_M,M19_X,[oo,a/q]))')
def modp(x):
 x=QQ(x);d=int(x.denominator())%P
 if d==0:raise ZeroDivisionError((x,P))
 return int(x.numerator())%P*int(inverse_mod(d,P))%P
def cusp(n,d):return 'oo' if d==0 else str(QQ(n)/QQ(d))
def path(c,d):
 g,x,y=xgcd(ZZ(c),ZZ(d));return cusp(-x,d//g),cusp(y,c//g)
def evpy(num,den,t,l):
 ans=int(t[1]);c,d,a,b=0,1,int(num),int(den)
 while b:
  f=b;b=a%b;q=(a-b)//f;a=-f;e=d;d=-c;c=q*c+e;v=int(t[(c%N)*N+(d%N)])
  if v<0:raise AssertionError((c,d))
  ans=(ans+v)%P
 return (ans+l)%P
def build():
 t=np.full(N*N,-1,dtype=np.int16);valid=0
 for c in range(N):
  for d in range(N):
   if gcd(gcd(c,d),N)!=1:continue
   a,b=path(c,d);t[c*N+d]=modp(pari(f'mseval(M19_M,M19_X,[{a},{b}])'));valid+=1
 lq=pari('mseval(M19_M,M19_X,[oo,0])');l=modp(lq);rng=random.Random(20260713+N+P);rs=[(0,1),(1,1),(-1,1),(1,2),(-3,7),(17,20329)]
 while len(rs)<506:
  d=rng.randrange(1,2000000);n=rng.randrange(-2000000,2000001);g=gcd(n,d);rs.append((n//g,d//g))
 checks=[]
 for n,d in rs:
  r=QQ(n)/QQ(d);x=modp(pari(f'mseval(M19_M,M19_X,[oo,{r}])'));y=evpy(n,d,t,l)
  if x!=y:raise AssertionError((n,d,x,y))
  checks.append([n,d,x])
 return t,l,{'valid_projective_pairs':valid,'loverp_rational':str(lq),'loverp_mod_p':int(l),'validation_count':len(checks),'validation_sha256':hashlib.sha256(json.dumps(checks,separators=(',',':')).encode()).hexdigest(),'table_sha256':hashlib.sha256(t.tobytes()).hexdigest()}
@njit(cache=True)
def ev(num,den,N,p,t,l):
 ans=int(t[1]);c=0;d=1;a=num;b=den
 while b!=0:
  f=b;b=a%b;q=(a-b)//f;a=-f;e=d;d=-c;c=q*c+e;v=int(t[(c%N)*N+(d%N)])
  if v<0:return -999999
  ans=(ans+v)%p
 return (ans+l)%p
@njit(cache=True)
def powmod(a,e,m):
 r=1
 while e:
  if e&1:r=(r*a)%m
  a=(a*a)%m;e//=2
 return r
@njit(cache=True)
def twist_mod(q,N,p,t,l):
 s=0
 for a in range(1,q):
  chi=1 if powmod(a,(q-1)//2,q)==1 else -1
  v=ev(a,q,N,p,t,l)
  if v<0:return -1
  s+=chi*v
 return s%p
def eligible(q):return gcd(q,3*N)==1 and q%8==5 and kronecker(q,3)*kronecker(q,11)==-1
def vp(n,p):
 n=ZZ(n);e=0
 while n and n%p==0:n//=p;e+=1
 return e
def firstkp(q):
 out=[]
 for ell in prime_range(2,KBOUND+1):
  if ell in (P,q) or N%ell==0 or q%ell==0 or (ell-1)%P:continue
  ab=ZZ(pari.ellap(E,ell));at=ZZ(kronecker(q,ell))*ab
  if (at-ell-1)%P==0:
   out.append({'ell':int(ell),'a_base':int(ab),'a_twist':int(at),'group_order':int(ell+1-at)})
   if len(out)==2:return out
 raise RuntimeError((q,out))
def work(q,e1,e2):
 a=(e1-1)*(P-1)//P;b=(e2-1)*(P-1)//P;nz=a*b
 if nz%2:raise ArithmeticError(nz)
 return {'n':e1*e2,'nonzero_outer_weights':nz,'half_nonzero_outer_weights':nz//2,'evaluations':nz//2*(q-1)}
def main():
 st=time.time();table,loverp,tcert=build();_=twist_mod(5,N,P,table,int(loverp));rows=[];checked=0;modzero=0;exact_candidates=[]
 for q in prime_range(3,QMAX+1):
  if not eligible(q):continue
  checked+=1;m=int(twist_mod(int(q),N,P,table,int(loverp)))
  if m<0:raise ArithmeticError(('table',q))
  if m:continue
  modzero+=1;raw=QQ(pari(f'M19_T({int(q)})'));sha=raw/2
  if sha.denominator()!=1:raise ArithmeticError((q,raw,sha))
  si=ZZ(sha);e=vp(abs(si),P);exact_candidates.append({'q':int(q),'raw':str(raw),'sha':int(si),'valuation':e})
  if e<=0:continue
  k=firstkp(int(q));r={'q':int(q),'L_over_Omega_plus':str(raw),'sha_an':int(si),'v31_sha_an':e,'ell1':k[0]['ell'],'ell2':k[1]['ell'],'a_ell1_twist':k[0]['a_twist'],'a_ell2_twist':k[1]['a_twist'],**work(int(q),k[0]['ell'],k[1]['ell']),'checked_eligible_index':checked,'elapsed_s':time.time()-st};rows.append(r);print('FOUND',json.dumps(r,sort_keys=True),flush=True)
  if len(rows)>=MAXPOS:break
 rows.sort(key=lambda r:(r['evaluations'],r['n'],r['q']));s={'curve':'66b1','p':P,'q_min':3,'q_max':QMAX,'eligible_checked':checked,'mod31_zero_candidates':modzero,'exact_candidates':exact_candidates,'positive_found':len(rows),'best':rows[0] if rows else None,'rows':rows,'table_certificate':tcert,'sage_version':sage_version,'pari_version':str(pari('version()')),'elapsed_s':time.time()-st,'reference_q53861_evaluations':476984160000}
 if rows:s['best_speedup_vs_q53861']=476984160000/rows[0]['evaluations']
 os.makedirs('research/bsd/results',exist_ok=True);open('research/bsd/results/m19_p31_target_scan.json','w').write(json.dumps(s,indent=2,sort_keys=True));fields=['q','L_over_Omega_plus','sha_an','v31_sha_an','ell1','ell2','a_ell1_twist','a_ell2_twist','n','nonzero_outer_weights','half_nonzero_outer_weights','evaluations','checked_eligible_index','elapsed_s'];f=open('research/bsd/results/m19_p31_target_scan.csv','w',newline='');w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows);f.close();print(json.dumps(s,sort_keys=True),flush=True)
if __name__=='__main__':main()
