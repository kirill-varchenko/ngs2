#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import copy
from itertools import tee
from Bio import SeqIO
import numpy as np
import argparse

def k_mers(s, k):
    """Генератор к-меров"""
    for i in range(len(s) - k + 1):
        yield s[i:i+k]

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def mergestrings(a, b) :
    """Объединяет строки по суффиксу из а и префиксу из b"""
    j = len(b)
    for i in range(len(a)) :
        if a[i:] == b[:j] :
            break
        else:
            j -= 1
    return a + b[j:]

class DeBruijnGraph:
    def __init__(self, read=None, k=None):
        self._outgoing = defaultdict(set)
        self._incoming = defaultdict(set)
        self._edge_data = defaultdict(list)
        self._node_coverage = defaultdict(int)

        if read is not None and k is not None:
            self.from_read(read, k+1)

    def from_read(self, read, k):
        for k_mer in k_mers(read, k):
            left = k_mer[:-1]
            right = k_mer[1:]
            self._node_coverage[left] += 1
            self.add_edge(left, right, label=k_mer)
        self._node_coverage[right] += 1

    def add_edge(self, u, v, **newdata):
        if 'label' not in newdata:
            newdata['label'] = ''
        if 'multiplicity' not in newdata:
            newdata['multiplicity'] = 1

        self._outgoing[u].add(v)
        self._incoming[v].add(u)
        for data in self._edge_data[(u, v)]:
            if data['label'] == newdata['label']:
                data['multiplicity'] += newdata['multiplicity']
                break
        else:
            self._edge_data[(u, v)].append(copy.copy(newdata))
            
    def remove_edge(self, u, v):
        if (u, v) not in self._edge_data:
            return
        
        del self._edge_data[(u, v)]
        self._outgoing[u].remove(v)
        self._incoming[v].remove(u)
        
    def remove_isolated_nodes(self):
        Q = {node for node in self._incoming if len(self._incoming[node]) == 0} & {node for node in self._outgoing if len(self._outgoing[node]) == 0}
        while Q:
            node = Q.pop()
            if node in self._incoming:
                del self._incoming[node]
            if node in self._outgoing:
                del self._outgoing[node]

    def condense(self):
        """Сжатие графа"""
        for u, v in self._edge_data:
            for data in self._edge_data[(u, v)]:
                if 'mean_coverage' not in data:
                    data['mean_coverage'] = (self._node_coverage[u] + self._node_coverage[v]) / 2
        
        # Очередь вершин с единичной входящей и исходящей степенью
        Q = {node for node in self._incoming if len(self._incoming[node]) == 1 and len(self._outgoing[node]) == 1}
        single_paths = []
        # Объединение вершин
        while Q:
            node = next(iter(Q))
            
            out_path = []
            u = node
            while u in Q:
                u = list(self._outgoing[u])[0]
                out_path.append(u)
            
            in_path = []
            u = node
            while u in Q:
                u = list(self._incoming[u])[0]
                in_path.append(u)
            
            single_path = list(reversed(in_path)) + [node] + out_path
            single_paths.append(single_path)
            Q -= set(single_path)
        
        for L in single_paths:
            # Объединение лейблов по всем рёбрам в пути
            label = ''
            for e in pairwise(L):
                if label == '':
                    label = self._edge_data[e][0]['label']
                else:
                    label = mergestrings(label, self._edge_data[e][0]['label'])
            
            # Среднее покрытие нового ребра
            mc = np.mean([self._edge_data[p][0]['mean_coverage'] for p in pairwise(L)])
            # Кратность нового ребра
            mm = min([self._edge_data[p][0]['multiplicity'] for p in pairwise(L)])          

            self.add_edge(L[0], L[-1], label=label, multiplicity=mm, mean_coverage=mc)

            for i in range(len(L) - 1):
                del self._edge_data[(L[i], L[i+1])]
                if i > 0:
                    del self._outgoing[L[i]]
                    del self._incoming[L[i]]

            self._outgoing[L[0]].remove(L[1])
            self._incoming[L[-1]].remove(L[-2])


    def write_dot(self, filename):
        nodes = list(self._node_coverage.keys())
        with open(filename, 'w') as fo:
            fo.write('digraph G {\nrankdir=LR\n')
            for u, v in self._edge_data:
                for data in self._edge_data[(u, v)]:
                    fo.write(f'"node_{nodes.index(u)}" -> "node_{nodes.index(v)}" [label="L = {len(data["label"])}, MC = {data.get("mean_coverage", -1):.2f} ({data["multiplicity"]})"];\n')
            fo.write('}\n')
        
    def drop_low_coveraged_edges(self, tiponly=False, q=0.6):
        """Удаление рёбер с небольшим покрытием
        q --- граница (квантиль) для покрытия 
        tiponly --- удалять только рёбра, одна из вершин которой является висящей"""
        
        # Если удалять не только висящие рёбра с небольшим покрытием, то можно
        # удалить не только ошибку, но и какой-нибудь редкий вариант
        
        mc = np.array([[x['mean_coverage'] for x in self._edge_data[e]] for e in self._edge_data]).flatten()
        cutoff_coverage = np.quantile(mc, q)
        
        edges2delete = []
        for e in self._edge_data:
            self._edge_data[e] = [x for x in self._edge_data[e] if x['mean_coverage'] >= cutoff_coverage]
            if len(self._edge_data[e]) == 0:
                edges2delete.append(e)
        
        for e in edges2delete:
            if tiponly:
                u_deg = len(self._incoming[e[0]]) + len(self._outgoing[e[0]])
                v_deg = len(self._incoming[e[1]]) + len(self._outgoing[e[1]])
                if u_deg > 1 and v_deg > 1:
                    continue
            self.remove_edge(*e)
        
        self.remove_isolated_nodes()
        
        
    def write_fasta(self, filename):
        with open(filename, 'w') as fo:
            for i, data in enumerate(self._edge_data.values()):
                fo.write(f'>Edge{i}\n{data[0]["label"]}\n')

    def __add__(self, other):
        res = copy.deepcopy(self)
        for u, v in other._edge_data:
            for data in other._edge_data[(u, v)]:
                res.add_edge(u, v, **data)
        for node, coverage in other._node_coverage.items():
            res._node_coverage[node] += coverage

        return res

parser = argparse.ArgumentParser()
parser.add_argument('input', help='input fasta/fastq file')
parser.add_argument('-k', '--k-mer', help='k-mer size', default=55, type=int)
parser.add_argument('-t', '--clean-tip-only', help='remove only low coveraged tips (default: all edges)', action='store_true')
parser.add_argument('-f', '--fasta', help='output fasta from graph', action='store_true')
args = parser.parse_args()

input_file_name = args.input
base_file_name, _, input_file_type = input_file_name.rpartition('.')
k = args.k_mer 

dbg = DeBruijnGraph()
for seq in SeqIO.parse(input_file_name, input_file_type):
    if len(seq) <= k:
        continue
    dbg += DeBruijnGraph(str(seq.seq), k)
    dbg += DeBruijnGraph(str(seq.seq.reverse_complement()), k)

dbg.condense()

dbg.write_dot(f'{base_file_name}_condensed.gv')

dbg.drop_low_coveraged_edges(tiponly=args.clean_tip_only)
dbg.condense()

dbg.write_dot(f'{base_file_name}_cleaned.gv')

if args.fasta:
    dbg.write_fasta(f'{base_file_name}_cleaned.fasta')
