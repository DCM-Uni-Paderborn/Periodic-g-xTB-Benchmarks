#!/usr/bin/env python3
from pathlib import Path
import re
import math

CAMPAIGN = Path(__file__).resolve().parent.parent
R3 = CAMPAIGN / 'raw' / 'cache_runtime_r3' / 'runs'
B1 = CAMPAIGN / 'raw' / 'baseline_symmetry_v2_runs'
B2 = CAMPAIGN / 'raw' / 'baseline_partial_pbc_runs'

flt_re = re.compile(r'(?<![A-Za-z_])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')

def floats(s):
    return [float(x.replace('D','E').replace('d','e')) for x in flt_re.findall(s)]

def energy(lines):
    return [floats(x)[-1] for x in lines if 'ENERGY| Total FORCE_EVAL' in x]

def tagged(lines, tag):
    out=[]
    for x in lines:
        if tag in x:
            v=floats(x.split('|',1)[-1])
            if v: out.extend(v)
    return out

def debug_force(lines):
    # Numerical/analytical force pairs from CP2K DEBUG tables.
    out=[]
    active=False
    for x in lines:
        if 'DEBUG| Atom    E(' in x:
            active=True
            continue
        if 'DEBUG| Atom  Coordinate' in x:
            active=True
            continue
        if active:
            if not x.strip():
                active=False
                continue
            if 'DEBUG|' in x:
                v=floats(x.split('|',1)[-1])
                if len(v)>=2:
                    out.extend(v)
    return out

def pv(lines):
    out=[]
    active=False
    for x in lines:
        if any(k in x for k in ('DEBUG| Numerical pv_virial', 'DEBUG| Analytical pv_virial',
                                'DEBUG| Difference pv_virial')):
            active=True
            continue
        if active:
            if not x.strip():
                active=False
                continue
            if 'DEBUG|' in x:
                out.extend(floats(x.split('|',1)[-1]))
    return out

def pv_block(lines, header):
    out=[]; active=False
    for x in lines:
        if header in x:
            active=True
            continue
        if active:
            if not x.strip(): break
            if 'DEBUG|' in x:
                v=floats(x.split('|',1)[-1])
                if len(v)==3: out.extend(v)
    return out

def analytical_stress_tensor(lines):
    out=[]; active=False
    for x in lines:
        if 'STRESS| Analytical stress tensor' in x:
            active=True
            continue
        if active:
            if not x.strip(): break
            part=x.split('|',1)[-1]
            if re.match(r'\s*[xyz]\s+', part):
                v=floats(part)
                if len(v)==3: out.extend(v)
    return out

def force_summary(lines):
    out=[]; active=False
    row=re.compile(r'DEBUG\|\s+\d+\s+[xyz]\s+')
    for x in lines:
        if 'BEGIN OF SUMMARY' in x: active=True; continue
        if 'END OF SUMMARY' in x: active=False
        if active and row.search(x):
            v=floats(x.split('|',1)[-1])
            # atom, numerical, analytical, difference, optional percent error
            if len(v)>=4: out.append((v[1],v[2],v[3]))
    return out

def ksummary(lines):
    nirr=[]
    grid=[]
    for x in lines:
        if 'BRILLOUIN| List of Kpoints' in x:
            v=floats(x)
            if v: nirr.append(int(v[-1]))
        if 'BRILLOUIN| K-Point grid' in x:
            v=floats(x)
            if len(v)>=3: grid.append(tuple(int(z) for z in v[-3:]))
    return (grid[0] if grid else None, nirr[0] if nirr else None)

def maxdiff(a,b):
    if len(a)!=len(b): return math.inf, len(a), len(b)
    return max((abs(x-y) for x,y in zip(a,b)), default=0.0),len(a),len(b)

def load(p): return p.read_text(errors='replace').splitlines()

print('BASELINE_COMPARISON')
print('case\trc\tgrid\tnirr\tdE_Ha\tdForcePrint\tdForceDebug\tdStressPrint\tdVirial\tcounts')
for d in sorted(R3.iterdir()):
    if not d.is_dir(): continue
    name=d.name
    base=(B1/name if (B1/name).is_dir() else B2/name)
    r=load(d/'cp2k.out'); b=load(base/'cp2k.out')
    vals=[]
    labs=[]
    for lab,func in [('E',energy),('FP',lambda x: tagged(x,'FORCES|')),
                     ('FD',debug_force),('SP',lambda x: tagged(x,'STRESS|')),('PV',pv)]:
        md,nr,nb=maxdiff(func(r),func(b)); vals.append(md); labs.append(f'{lab}:{nr}/{nb}')
    grid,nirr=ksummary(r)
    rc=(d/'returncode.txt').read_text().strip()
    print(f'{name}\t{rc}\t{grid}\t{nirr}\t'+'\t'.join(f'{x:.3e}' if math.isfinite(x) else 'LEN' for x in vals)+'\t'+','.join(labs))

print('\nSYMMETRY_PAIRS')
print('reference\treduced\tgrid\tnk_ref/nk_red\tdE0_Ha\tdF_print\tdStress_print\tstress_debug_sum_red\tforce_debug_sum_red')
pairs=[
 ('ch4_full_debug','ch4_spglib_debug'),
 ('h2_full_force_stress','h2_time_reversal_force_stress'),
 ('gxtb_1d_x_k211_force_stress_full','gxtb_1d_x_k211_force_stress'),
 ('gxtb_1d_x_k211_force_stress_full','gxtb_1d_x_k211_force_stress_spglib'),
 ('gxtb_2d_xz_k212_force_stress_full','gxtb_2d_xz_k212_force_stress'),
 ('gxtb_2d_xz_k212_force_stress_full','gxtb_2d_xz_k212_force_stress_spglib'),
 ('si_shifted_full_energy','si_shifted_spglib_energy'),
]
def first_force(lines):
    out=[]; active=False
    for x in lines:
        if 'FORCES| Atomic forces' in x: active=True; continue
        if active:
            if not x.strip(): break
            if 'FORCES|' in x: out.extend(floats(x.split('|',1)[-1]))
    return out
def first_stress(lines):
    out=[]; active=False
    for x in lines:
        if 'STRESS| Analytical stress tensor' in x: active=True; continue
        if active:
            if not x.strip(): break
            if 'STRESS|' in x: out.extend(floats(x.split('|',1)[-1]))
    return out
def debug_sum(lines, periodic=False):
    key='Periodic-subspace sum of differences' if periodic else 'DEBUG| Sum of differences'
    vals=[]
    for x in lines:
        if key in x:
            v=floats(x)
            if v: vals.append(v[-1])
    return vals[0] if vals else float('nan')
for a,bn in pairs:
    la=load(R3/a/'cp2k.out'); lb=load(R3/bn/'cp2k.out')
    ga,na=ksummary(la); gb,nb=ksummary(lb)
    e0=abs(energy(la)[0]-energy(lb)[0])
    df=maxdiff(first_force(la),first_force(lb))[0]
    ds=maxdiff(first_stress(la),first_stress(lb))[0]
    print(f'{a}\t{bn}\t{ga}\t{na}/{nb}\t{e0:.3e}\t'+
          (f'{df:.3e}' if math.isfinite(df) else 'LEN')+'\t'+
          (f'{ds:.3e}' if math.isfinite(ds) else 'LEN')+f'\t{debug_sum(lb):.3e}\t{debug_sum(lb, periodic=True):.3e}')

print('\nDERIVATIVE_SUMMARIES')
print('case\tstress_sum_Ha\tperiodic_stress_sum_Ha\tforce_sum_HaBohr')
for d in sorted(R3.iterdir()):
    lines=load(d/'cp2k.out')
    ss=debug_sum(lines); ps=debug_sum(lines,periodic=True)
    # Force summary is the colon-bearing line, distinct from stress line.
    fs=float('nan')
    for x in lines:
        if 'DEBUG| Sum of differences:' in x:
            v=floats(x)
            if v: fs=v[0]
    if any(math.isfinite(v) for v in (ss,ps,fs)):
        print(f'{d.name}\t{ss:.12e}\t{ps:.12e}\t{fs:.12e}')

print('\nBASELINE_PHYSICAL_COMPONENTS')
print('case\tdE0_Ha\tdFcentral_HaBohr\tdFnum_HaBohr\tdFana_HaBohr\tdVirNum_Ha\tdVirAna_Ha\tdStressTensor_bar')
for d in sorted(R3.iterdir()):
    if not d.is_dir(): continue
    base=(B1/d.name if (B1/d.name).is_dir() else B2/d.name)
    r=load(d/'cp2k.out'); b=load(base/'cp2k.out')
    fr=force_summary(r); fb=force_summary(b)
    fnum=maxdiff([x[0] for x in fr],[x[0] for x in fb])[0]
    fana=maxdiff([x[1] for x in fr],[x[1] for x in fb])[0]
    vals=[abs(energy(r)[0]-energy(b)[0]),
          maxdiff(first_force(r),first_force(b))[0],fnum,fana,
          maxdiff(pv_block(r,'DEBUG| Numerical pv_virial'),pv_block(b,'DEBUG| Numerical pv_virial'))[0],
          maxdiff(pv_block(r,'DEBUG| Analytical pv_virial'),pv_block(b,'DEBUG| Analytical pv_virial'))[0],
          maxdiff(analytical_stress_tensor(r),analytical_stress_tensor(b))[0]]
    print(d.name+'\t'+'\t'.join(f'{v:.3e}' if math.isfinite(v) else 'LEN' for v in vals))

print('\nSYMMETRY_PHYSICAL_COMPONENTS')
print('reference\treduced\tdE0_Ha\tdFcentral_HaBohr\tdFana_HaBohr\tdVirAna_Ha\tdStressTensor_bar')
for a,bn in pairs:
    la=load(R3/a/'cp2k.out'); lb=load(R3/bn/'cp2k.out')
    fa=force_summary(la); fb=force_summary(lb)
    vals=[abs(energy(la)[0]-energy(lb)[0]),maxdiff(first_force(la),first_force(lb))[0],
          maxdiff([x[1] for x in fa],[x[1] for x in fb])[0],
          maxdiff(pv_block(la,'DEBUG| Analytical pv_virial'),pv_block(lb,'DEBUG| Analytical pv_virial'))[0],
          maxdiff(analytical_stress_tensor(la),analytical_stress_tensor(lb))[0]]
    print(f'{a}\t{bn}\t'+'\t'.join(f'{v:.3e}' if math.isfinite(v) else 'LEN' for v in vals))

print('\nBASELINE_SYMMETRY_METADATA')
print('case\tgrid_r3/base\tnirr_r3/base')
for d in sorted(R3.iterdir()):
    if not d.is_dir(): continue
    base=(B1/d.name if (B1/d.name).is_dir() else B2/d.name)
    kr=ksummary(load(d/'cp2k.out')); kb=ksummary(load(base/'cp2k.out'))
    print(f'{d.name}\t{kr[0]}/{kb[0]}\t{kr[1]}/{kb[1]}')
