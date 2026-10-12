"""Build revised CSV tables directly from checked structured records."""
from pathlib import Path
import csv
import json
from regenerate_tables import read_tables, write_tables

HERE = Path(__file__).resolve().parents[1]


def load(name):
    return json.loads((HERE/'data'/(name+'.json')).read_text(encoding='utf-8'))


def run():
    catalog = read_tables()
    for row in catalog['S7']['rows']:
        if row[0] == 'This paper':
            row[-1] = ('Boundary certificate; slack identity; joint relaxation (Theorem 6); '
                       'binary chain (Theorem 7); hierarchy obstruction (Theorem 8); '
                       'fixed-cap candidate NP-completeness (Theorem 10)')
        elif 'Bogart' in str(row[0]) or 'Cuellar' in str(row[0]):
            row[-1] = row[-1].replace('Theorem 6', 'Theorem 8')
        elif row[0] == 'Singer [8]':
            row[-1] = 'Declared trace-zero seed of Section 8.3; not every Singer-family variant'
        elif row[0] == 'Martos et al. [7]':
            row[-1] = 'Published heuristic context; original algorithm not reproduced'

    def add(key, caption, headers, rows):
        assert all(len(r) == len(headers) for r in rows), key
        catalog[key] = dict(caption=caption, headers=headers, rows=rows)

    catalog['1']['rows'] = [r if 'B_{\\mathrm{pair}}' not in r[0] else
        [r'$B_{\mathrm{ph}}$', 'Pair-counting certificate', 'Universal lower bound'] for r in catalog['1']['rows']]
    add('2', 'Model-level comparison, not a proved LP equivalence or ordering. Classical length LPs retain pair labels and fractional assignments and optimize span. B_bd discards pair labels, minimizes a capacitated window second moment, and converts it into a multiplicity lower certificate.',
        ['Method','Primary variables','Fixed data','Objective','Retained information'], [
        ['Assignment length LP',r'$z_i,w_{ij,d}$','n, cap one','Minimum span','Pair labels; fractional distance assignments'],
        ['Flow-equivalent LPs [18]',r'$\varphi_{ij}$','n, cap one','Length lower bound','Equivalent classical LP formulations'],
        ['Pair counting',r'$r_d$','L,n','Multiplicity floor','Total number of pairs'],
        ['B_CS',r'$a_x\in\mathbb R$','L,n,u','Minimum square sum','Mass only'],
        ['B_int',r'$a_x\in\mathbb Z$','L,n,u','Minimum square sum','Integer mass'],
        ['B_bd',r'$a_x\in\mathbb Z$','L,n,u','Minimum square sum','Window second moments; tent capacities; pair labels discarded'],
        ['B_walk',r'$a_x\in\mathbb Z$','L,n,u','Minimum square sum','Capacities and adjacent increments'],
        ['B_real',r'$b_i\in\{0,1\}$','L,n,u','Minimum square sum','One binary window profile']])
    catalog['3']['rows'] = [r for r in catalog['3']['rows'] if r[0] not in ('B_walk','B_real')]
    for bound,relax,comp in [('B_walk','Capacities, endpoint conditions and adjacent increments',r'$B_{\mathrm{bd}}\leq B_{\mathrm{walk}}\leq B_{\mathrm{real}}$'),
                             ('B_real','Exact binary convolution minimum at each u',r'$B_{\mathrm{real}}\leq\Delta$')]:
        catalog['3']['rows'].append([bound,relax,comp])
    cases = load('seven_case_hierarchy')
    add('4', 'All seven boundary failures. Window triples are u_bd/u_walk/u_real (first maximizers); square-sum triples are Q_bd/Q_walk/Q_real, all at u_real. Closed means B_real equals Delta. Obstruction labels are defined in Section 5.3.',
        ['(L,n)',r'$\Delta$','B_bd','B_walk','B_real','Window triple','Square-sum triple','Closed by B_real?','Obstruction'],
        [[f"({r['L']},{r['n']})",r['Delta'],r['B_bd'],r['B_walk'],r['B_real'],
          f"{r['u_bd']}/{r['u_walk']}/{r['u_real']}",f"{r['Q_bd']}/{r['Q_walk']}/{r['Q_real']}",
          'Yes' if r['gap_closed'] else 'No','Binary coupling' if r['gap_closed'] else 'Scalar moment'] for r in cases])
    milp = load('milp_complete')
    from extremal import bounds
    add('6', 'Four archived SciPy/HiGHS feasibility solves, independently rechecked in this revision. '
        'A feasible witness at g plus the meeting lower certificate proves exact Delta. '
        'Times and nodes belong to the recorded solves, not new reruns; a timeout would not prove infeasibility.',
        ['(L,n)','Target g','Lower certificate','Decision','Seconds','Nodes',r'Exact $\Delta$'],
        [[f"({r['L']},{r['n']})",r['target'],f"B_bd = {bounds(r['L'],r['n'])['B_bd']}",
          r['decision'],f"{r['seconds']:.3f}",r['node_count'],r['exact_Delta']] for r in milp])
    add('7','Input models for complexity. Only the explicit fixed-cap selection theorem is claimed here.',
        ['Problem','Input','Cap','Complexity and contribution scope'],[
        ['Complete-interval extremum','L,n encoded succinctly','Variable','No hardness classification asserted here'],
        ['Classical ruler construction [4,5]','Model-specific marks, length or admissible coordinates','1','Prior algorithmic and hardness results; not a new case'],
        ['Golomb subruler [4,5]','Explicit candidate set, target size','1','Prior distinct-difference selection setting'],
        ['Bounded-overlap selection','Explicit candidate set, target size',r'$g\geq2$','NP-complete for each fixed g (Theorem 10)']])
    summary = load('modular_summary')
    constructors = load('modular_constructors')
    pairs = sorted({(r['L'],r['n']) for r in summary})
    rows = []
    families = ['Paley','Singer','Ruzsa-Spence','RDS','Random']
    for L,n in pairs:
        cells = {r['family']:r for r in summary if (r['L'],r['n']) == (L,n)}
        assert set(cells) == set(families)
        lower = cells['Paley']['B_comb']
        values = [cells[f]['best'] for f in families]
        rows.append([f'({L},{n})',lower,*values,min(values)-lower])
    add('9','Four declared modular seed variants and a random arm, with an equal budget of 30 streams of 3000 proposals per cell. Best is the observed best upper witness, not a solver optimum. Gap is the best listed witness minus B_comb.',
        ['(L,n)',r'$B_{\mathrm{comb}}$','Paley','Singer','Ruzsa-Spence','RDS','Random','Gap'],rows)
    add('S9','Executed verification domains. Archived timings and solver outcomes are not represented as fresh reruns.',
        ['Domain','Records','Verification'],[
        ['Nontrivial finite map','252','Complete re-enumeration of lengths 4 through 24'],
        ['Binary-window minima','4046','Every window minimum independently regenerated'],
        ['Full-interval trivial rows','21','Full-set witness and exact bound checked'],
        ['Seven exceptional walk cases','126 windows','Dense and sparse integer DP agree'],
        ['Original / strengthened residual vectors','5 / 5','Exact pair, moment, capacity and parity checks'],
        ['Literature grid','462','Source-cell inversion and 274 / 188 coverage split'],
        ['Principal medium witnesses','9','Full coordinate and difference validation'],
        ['Modular construction refinement','1350 runs','Pairwise histogram validation of all returned witnesses'],
        ['Marker-reduction fixtures','406','Finite construction checks; not a universal proof'],
        ['Cross-platform export','LF / CRLF / CR','Normalization regression PASS on Windows; Linux not run']])
    add('S10','Model mapping. A closed-form window certificate is not asserted equivalent to, or stronger than, every classical length LP.',
        ['Method','Domain and variables','Fixed data','Finite tent / binary profile','Relation to target'],[
        ['Classical assignment LP',r'$z_i,w_{ij,d}\geq0$','n, distance horizon H, cap one','Not explicit window tent / no','Length lower bound excludes cap one when above L-1'],
        ['Equivalent flow descriptions [18]',r'$\varphi_{ij}$ on mark-pair arcs','n, cap one','Not explicit window tent / no','Equivalent LP length bounds; not independent improvements'],
        ['Pair counting','Unlabelled nonnegative multiplicities','L,n','No / no','Mass divided among L-1 positive distances'],
        ['B_CS','Real window counts','L,n,u','No / no','Closed-form lower bound on maximum multiplicity'],
        ['B_int','Integer window counts','L,n,u','No / no','Integer tightening of real balancing'],
        ['B_bd','Capacitated integer window counts','L,n,u','Yes / no','Closed-form layer filling, slack and strictness'],
        ['B_walk','Increment-compatible integer counts','L,n,u','Yes / no','Dynamic-programming relaxation; shared bits discarded'],
        ['B_real','Counts from an n-element binary set','L,n,u','Yes / yes','Exact scalar binary minimum; joint distance structure not enforced']])
    add('S11a','All 45 construction cells. Each uses 30 seeds and 3000 proposals. After spread is sample SD; no inferential significance or optimality is implied by these descriptive values.',
        ['(L,n)','Family','Before mean','After mean +/- SD','Best',r'$B_{\mathrm{comb}}$','Gap'],
        [[f"({r['L']},{r['n']})",r['family'],f"{r['before_mean']:.3f}",f"{r['after_mean']:.3f} +/- {r['after_SD']:.3f}",r['best'],r['B_comb'],r['remaining_gap']] for r in summary])
    add('S11b','Every algebraic constructor. N is the original modulus, m the final folded modulus, g the modular cap and loss the actual collision deletion count. Eligibility is actual for explicitly counted Singer rows.',
        ['(L,n)','Family','q','N','m','g','Loss','Eligibility'],
        [[f"({r['L']},{r['n']})",r['family'],r['q'],r['original_modulus'],r['modulus'],r['modular_cap'],
          r.get('discarded_collisions',0),'Actual' if r.get('cardinality_check')=='explicit' else 'Guaranteed']
         for r in constructors if r['family'] != 'Random'])
    add('S12a','First maximizing windows for every layer of the complete seven-case hierarchy. All ties and per-window records are archived separately.',
        ['(L,n)',r'$\Delta$','B_bd','B_walk','B_real',r'$u_{\mathrm{bd}}$',r'$u_{\mathrm{walk}}$',r'$u_{\mathrm{real}}$'],
        [[f"({r['L']},{r['n']})",r['Delta'],r['B_bd'],r['B_walk'],r['B_real'],r['u_bd'],r['u_walk'],r['u_real']] for r in cases])
    add('S12b','Square sums evaluated at the common reference window u_real. A larger square sum need not change its ceiling bound.',
        ['(L,n)',r'$u_{\mathrm{real}}$',r'$Q_{\mathrm{bd}}$',r'$Q_{\mathrm{walk}}$',r'$Q_{\mathrm{real}}$','Obstruction','Gap closed'],
        [[f"({r['L']},{r['n']})",r['u_real'],r['Q_bd'],r['Q_walk'],r['Q_real'],
          'Binary coupling' if r['gap_closed'] else 'Scalar moment','Yes' if r['gap_closed'] else 'No'] for r in cases])
    residual = load('residual_diagnostics')
    add('S13a','Structural defects of the original vectors and odd-distance mass of the new vectors satisfying both parity and tail capacities.',
        ['(L,n)','Original odd mass','Original parity valid','Original tail violations d','New odd mass'],
        [[f"({r['L']},{r['n']})",r['old_odd_mass'],'Yes' if r['old_parity_feasible'] else 'No',
          ','.join(map(str,r['old_tail_violations'])) or 'None',sum(r['strengthened_witness'][::2])] for r in residual])
    add('S13b','Final witnesses used in main Theorem 8. These strengthened moment-feasible vectors satisfy tail capacity and parity; each has maximum B_bd below the completely enumerated optimum. Original vectors remain archived separately.',
        ['(L,n)',r'$r=(r_1,\ldots,r_{L-1})$'],
        [[f"({r['L']},{r['n']})",'$('+','.join(map(str,r['strengthened_witness']))+')$'] for r in residual])
    family = load('infinite_family_checks')
    add('S14','Exact arithmetic checks of Theorem 4: the first 20 consecutive k values and three larger checks. The proved lower gain may be nonpositive at small k; strict positivity is guaranteed for k>=14. These finite checks do not prove infinitude.',
        ['k','L','n','B_int','B_bd','Actual gain','Proved lower gain','Proved upper gain'],
        [[r['k'],r['L'],r['n'],r['B_int'],r['B_bd'],r['actual_gain'],r['proved_gain'],r['k']] for r in family if 3 <= r['k'] <= 22 or r['k'] in (50,80,100)])
    add('S15','Claim, evidence and scope. Computed domains do not restrict a proved universal theorem, and finite checks do not replace human proof review.',
        ['Claim','Evidence','Scope','Not claimed'],[
        ['B_bd lower certificate','Theorem 1',r'$2\leq n\leq L$','Global ordering relative to prior LP or diameter bounds'],
        [r'$\Theta(\sqrt L)$ gain','Theorem 4',r'$L=4k^2,n=2k^2$; positive for $k\geq14$',r'Second-order asymptotics of $\Delta$'],
        ['B_real benchmark','Theorem 7; enumeration',r'Universal validity; evaluated for $L\leq24$','Scalable general algorithm'],
        ['Five residual obstructions','Theorem 8; Table S13b','Five stated parameter pairs','Infinite obstruction family or impossibility for all window methods'],
        ['Medium exact values','Table 6; checked witnesses','Lower and upper certificates meet on listed rows','Optimality from an incumbent or timeout alone'],
        ['Fixed-cap NP-completeness','Theorem 10; Lemmas 11--12',r'Explicit candidates; each fixed $g\geq2$','Hardness for succinct complete intervals']])
    write_tables(catalog)
    # CSV profiles retain complete rows; figure and document tables are not the data archive.
    data_sets = {'seven_case_hierarchy':cases,'modular_summary':summary,'residual_diagnostics':residual,
                 'realizability_audit':load('realizability_audit'),
                 'infinite_family_checks':family}
    data_sets['seven_case_profiles'] = [dict(L=r['L'],n=r['n'],**p) for r in load('seven_case_profiles') for p in r['per_u']]
    for name, data in data_sets.items():
        fields = list(dict.fromkeys(k for r in data for k in r))
        with (HERE/'data'/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:
            writer = csv.DictWriter(f,fields,lineterminator='\n');writer.writeheader()
            writer.writerows({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in data)
    print('Built',len(catalog),'tables; 45 construction cells and 126 complete profile records.')


if __name__ == '__main__':
    run()
