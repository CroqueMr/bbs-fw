import glob,json,os,sys
root,outdir=sys.argv[1],sys.argv[2]; os.makedirs(outdir,exist_ok=True)
case='141c1_q3461_p7'
files=sorted(glob.glob(os.path.join(root,'**',f'{case}_block_*.json'),recursive=True))
if len(files)!=16: raise SystemExit(f'expected 16 blocks, found {len(files)}')
data=[json.load(open(f,encoding='utf-8')) for f in files]
assert [x['block'] for x in data]==list(range(16))
rows=[r for b in data for r in b['shards']]
assert [r['shard'] for r in rows]==list(range(128))
assert rows[0]['lo']==1 and rows[-1]['hi']==data[0]['n']-1
assert all(rows[i]['hi']+1==rows[i+1]['lo'] for i in range(127))
for key in ['case','p','q','ell1','ell2','n','sage_version']:
    assert len({str(b[key]) for b in data})==1
p=data[0]['p']; residue=sum(r['lift'] for r in rows)%p
out={k:data[0][k] for k in ['case','p','q','ell1','ell2','n','sage_version']}
out.update({'blocks':16,'nshards':128,'delta_n_mod_p':residue,'nonzero':residue!=0,'block_partial_residues':[b['partial_mod_p'] for b in data],'partial_sum_integer':sum(r['lift'] for r in rows),'total_nonzero_weights':sum(r['nonzero_weights'] for r in rows),'elapsed_cpu_sum_s':sum(r['elapsed_s'] for r in rows),'max_block_wall_s':max(b['elapsed_wall_s'] for b in data)})
with open(os.path.join(outdir,case+'_aggregate.json'),'w',encoding='utf-8') as fh: json.dump(out,fh,indent=2,sort_keys=True)
print(json.dumps(out,sort_keys=True))
