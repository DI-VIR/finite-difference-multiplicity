"""Instantiate the Paley modular g-Golomb construction at the main parameters."""

from __future__ import annotations

from pathlib import Path

import csv

import json

import random

import statistics

import sys

from extremal import bounds, delta

SOURCE = "https://arxiv.org/html/2607.07931#S5"

PAIRS = [(24,11),(32,15),(40,19),(48,22),(64,29),(80,37),
         (160,71),(256,118),(320,142)]

def prime(x):
    if x < 2:
        return False
    return all(x % q for q in range(2, int(x ** .5) + 1))

def paley_circle(p):
    residues = sorted({(x*x) % p for x in range(1,p)})
    assert len(residues) == (p-1)//2
    g = (p-1)//4
    modular = [0]*p
    rset=set(residues)
    for a in residues:
        for b in residues:
            if a != b:
                modular[(b-a)%p] += 1
    assert max(modular[1:]) <= g
    gaps=[(residues[(i+1)%len(residues)]-residues[i])%p
          for i in range(len(residues))]
    max_gap=max(gaps)
    cuts=[]
    for i,gap in enumerate(gaps):
        if gap == max_gap:
            start=residues[(i+1)%len(residues)]
            cuts.append(sorted((r-start)%p for r in residues))
    return min(cuts),g,max_gap

def subset(full, n):
    if len(full)==n:
        return full
    # A deterministic thinning; it is not asserted to be Paley-optimal.
    chosen=[full[(i*(len(full)-1))//(n-1)] for i in range(n)]
    assert len(chosen)==len(set(chosen))==n
    return chosen

def generate(length,n):
    for p in range(2*n+1, 4*length):
        if not prime(p):
            continue
        full,g,gap=paley_circle(p)
        if len(full)<n:
            continue
        taps=subset(full,n)
        if max(taps)<length:
            return p,g,gap,taps
    raise ValueError((length,n))
