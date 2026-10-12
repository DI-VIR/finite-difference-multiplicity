"""Matched initialization and budget ablation; no superiority is presumed."""
from collections import Counter
from pathlib import Path
import json
import math
import random
import statistics
import time

from constructions.modular_comparison import PAIRS, FAMILIES, select_seed
from extremal import delta

ROOT = Path(__file__).resolve().parents[1]
BUDGETS = (0,3000,30000)


def checkpoints(L,n,initial,seed,anneal=False):
    rng = random.Random(seed)
    current = list(initial)
    empty = sorted(set(range(L))-set(current))
    bits = sum(1<<t for t in current)
    value = best = delta(bits,L)
    witness = sorted(current)
    output = {0:(best,witness)}
    for j in range(1,max(BUDGETS)+1):
        i,k = rng.randrange(n),rng.randrange(L-n)
        x,y = current[i],empty[k]
        candidate = bits^(1<<x)^(1<<y)
        score = delta(candidate,L,stop_above=None if anneal else value)
        temperature = 1.0*(0.02**((j-1)/(max(BUDGETS)-1)))
        accept = score <= value or (anneal and rng.random() < math.exp((value-score)/temperature))
        if accept:
            bits,value = candidate,score
            current[i],empty[k] = y,x
            if value < best:
                best,witness = value,sorted(current)
        if not anneal:
            best,witness = value,sorted(current)
        if j in BUDGETS:
            assert max(Counter(b-a for i,a in enumerate(witness) for b in witness[i+1:]).values())==best
            output[j]=(best,witness.copy())
    return output


def run():
    start=time.perf_counter()
    records=[]
    for L,n in PAIRS:
        for family in (*FAMILIES,'Random-SA'):
            initial = None if family.startswith('Random') else select_seed(L,n,family)[0]
            for s in range(30):
                seed=2026100200+s
                taps=random.Random(seed).sample(range(L),n) if initial is None else initial
                for budget,(value,witness) in checkpoints(L,n,taps,seed,family=='Random-SA').items():
                    records.append(dict(L=L,n=n,family=family,seed=seed,budget=budget,
                        value=value,witness=witness,initial=sorted(taps)))
            print(L,n,family,'best30k',min(r['value'] for r in records if (r['L'],r['n'],r['family'],r['budget'])==(L,n,family,30000)),flush=True)
        (ROOT/'data/seed_budget_runs.json').write_text(json.dumps(records,indent=2)+'\n')
    summary=[]
    for L,n in PAIRS:
        for family in (*FAMILIES,'Random-SA'):
            for budget in BUDGETS:
                values=[r['value'] for r in records if (r['L'],r['n'],r['family'],r['budget'])==(L,n,family,budget)]
                summary.append(dict(L=L,n=n,family=family,budget=budget,runs=len(values),
                    mean=statistics.mean(values),SD=statistics.stdev(values),best=min(values)))
    # Conditional on the deterministic seeds and shared RNG identifiers, not independent algebraic seed populations.
    contrasts=[]
    for L,n in PAIRS:
        for family in FAMILIES[:-1]:
            for budget in (3000,30000):
                cells={f:{r['seed']:r['value'] for r in records if (r['L'],r['n'],r['family'],r['budget'])==(L,n,f,budget)} for f in (family,'Random')}
                diff=[cells[family][s]-cells['Random'][s] for s in sorted(cells[family])]
                rng=random.Random(2026100500+L+budget)
                boots=sorted(statistics.mean(rng.choices(diff,k=30)) for _ in range(10000))
                contrasts.append(dict(L=L,n=n,family=family,budget=budget,
                    algebraic_minus_random_mean=statistics.mean(diff),
                    paired_bootstrap_95_low=boots[249],paired_bootstrap_95_high=boots[9749],
                    bootstrap_replicates=10000,independent_blocks=30,
                    interpretation='Exploratory pointwise percentile interval; no multiplicity-adjusted superiority test'))
    for name,value in [('seed_budget_summary',summary),('seed_budget_contrasts',contrasts)]:
        (ROOT/'data'/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
    report=dict(status='PASS',streams=9*6*30,checkpoint_records=len(records),budgets=BUDGETS,
        seconds=round(time.perf_counter()-start,3),SA_schedule='T(j)=0.02^((j-1)/29999), best incumbent retained',
        protocol='All nine original pairs; thirty shared identifiers per arm; direct maximum multiplicity objective; initial and final full witnesses retained.')
    (ROOT/'data/seed_budget_protocol.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report,flush=True)


if __name__=='__main__':
    run()
