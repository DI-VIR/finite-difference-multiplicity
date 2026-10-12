"""Executable finite fixtures for the complete candidate-set reduction.

The manuscript proves the theorem for all fixed g >= 2. These small fixtures
check implementation of the reduction, not the universal theorem itself.
"""
from itertools import combinations
from collections import Counter


def marker(s,i):return 2*s*i*i+i


def construct(N,edges,g,k):
    vertices=[marker(N,i) for i in range(N)]
    D=2*vertices[-1]+1;K=g*len(edges)
    marks=vertices.copy()
    for l,(u,w) in enumerate(edges):
        for t in range(1,g+1):
            offset=D*marker(K+1,g*l+t)
            marks.extend((offset+vertices[u],offset+vertices[w]))
    return sorted(marks),k+2*g*len(edges)


def multiplicities(marks):
    return Counter(b-a for i,a in enumerate(marks) for b in marks[i+1:])


def verify():
    cases=0
    for N in (3,4):
        possible=list(combinations(range(N),2))
        for mask in range(1,1<<len(possible)):
            edges=[e for j,e in enumerate(possible) if (mask>>j)&1]
            for g in (2,3):
                vertices=[marker(N,i) for i in range(N)]
                A,_=construct(N,edges,g,1)
                cc=multiplicities(A)
                edge_lengths={abs(vertices[u]-vertices[w]) for u,w in edges}
                assert all(cc[d]==g+1 for d in edge_lengths)
                assert all(r<=2 for d,r in cc.items() if d not in edge_lengths)
                for k in range(1,N):
                    independent=any(all(not ({u,w}<=set(W)) for u,w in edges)
                                    for W in combinations(range(N),k))
                    A,target=construct(N,edges,g,k)
                    # Enumerate complements; at most N-k marks can be omitted.
                    feasible=any(max(multiplicities([a for a in A if a not in omitted]).values(),default=0)<=g
                                 for omitted in combinations(A,len(A)-target))
                    assert feasible==independent,(N,edges,g,k)
                    cases+=1
    return cases


if __name__=='__main__':print('Candidate-set reduction fixtures:',verify())
