from sage.all import *
from sage.version import version as sage_version
import csv,json,os,time
AINVS=[1,1,1,-2,-1];N=66;P=31;QMAX=300000;MAXPOS=30;KBOUND=300000
pari.allocatemem(768_000_000);E=pari.ellinit(AINVS);Z=pari.msfromell(E,1);pari(f'M={Z[0]};X={Z[1]};T(q)=sum(a=0,q-1,kronecker(q,a)*mseval(M,X,[oo,a/q]))')
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
 nz1=(e1-1)*(P-1)//P;nz2=(e2-1)*(P-1)//P;nz=nz1*nz2
 if nz%2:raise ArithmeticError(nz)
 return {'n':e1*e2,'nonzero_outer_weights':nz,'half_nonzero_outer_weights':nz//2,'evaluations':nz//2*(q-1)}
st=time.time();rows=[];checked=0
for q in prime_range(3,QMAX+1):
 if not eligible(q):continue
 checked+=1;raw=QQ(pari(f'T({int(q)})'));sha=raw/2
 if sha.denominator()!=1:raise ArithmeticError((q,raw,sha))
 si=ZZ(sha);e=vp(abs(si),P)
 if e<=0:continue
 k=firstkp(int(q));r={'q':int(q),'L_over_Omega_plus':str(raw),'sha_an':int(si),'v31_sha_an':e,'ell1':k[0]['ell'],'ell2':k[1]['ell'],'a_ell1_twist':k[0]['a_twist'],'a_ell2_twist':k[1]['a_twist'],**work(int(q),k[0]['ell'],k[1]['ell']),'checked_eligible_index':checked,'elapsed_s':time.time()-st};rows.append(r);print('FOUND',json.dumps(r,sort_keys=True),flush=True)
 if len(rows)>=MAXPOS:break
rows.sort(key=lambda r:(r['evaluations'],r['n'],r['q']));s={'curve':'66b1','p':P,'q_min':3,'q_max':QMAX,'eligible_checked':checked,'positive_found':len(rows),'best':rows[0] if rows else None,'rows':rows,'sage_version':sage_version,'pari_version':str(pari('version()')),'elapsed_s':time.time()-st,'reference_q53861_evaluations':476984160000}
if rows:s['best_speedup_vs_q53861']=476984160000/rows[0]['evaluations']
os.makedirs('research/bsd/results',exist_ok=True);open('research/bsd/results/m19_p31_target_scan.json','w').write(json.dumps(s,indent=2,sort_keys=True));fields=['q','L_over_Omega_plus','sha_an','v31_sha_an','ell1','ell2','a_ell1_twist','a_ell2_twist','n','nonzero_outer_weights','half_nonzero_outer_weights','evaluations','checked_eligible_index','elapsed_s'];f=open('research/bsd/results/m19_p31_target_scan.csv','w',newline='');w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows);f.close();print(json.dumps(s,sort_keys=True),flush=True)
