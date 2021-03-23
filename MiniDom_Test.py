# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:02:11 2021

@author: theodore
"""

import xml.etree.cElementTree as ET
import networkx
import logging
import pylab

class KeggPathway(networkx.DiGraph):
    title = ''
    labels = {}
    reactions = {}
    def __repr__(self):
        return self.title + ' pathway' # TODO: __init__ method to make sure self.title exists
    def get_genes(self):
        """
        return a subgraph composed only by the genes
        >>> p = KeggPathway()
        >>> p.add_node('gene1', data={'type': 'gene'})
        >>> p.add_node('compound1', data={'type': 'compound'})
        >>> subgraph = p.get_genes()
        >>> print subgraph.nodes()
        ['gene1']
        """
#        subgraph = self.subgraph([node for node in self.nodes() if self.node[node]['type'] == 'gene'])
        genes = []
        labels = {}
        for n in self.nodes():
            try:
                if self.node[n]['type'] == 'gene':
                    genes.append(n)
                    labels[n] = self.node[n]
            except:
                pass
#            else:
#                self.labels.pop(node)
        subgraph = self.subgraph(genes)
        subgraph.title = self.title + ' (genes)'
        subgraph.labels = labels
        return subgraph

def KGML2Graph(xmlfile, filter_by = ()):
    pathway = KeggPathway()
    nodes = {}
    genes = []
    pathway.reactions = {}
    pathway.relations = {}
    pathway.labels = {}     # dictionary to keep node labels (gene name?)
    tree = ET.parse(xmlfile)
    organism = tree.getroot().get('org')
    if organism == 'ko':
        entriestype = ('ortholog', 'map', 'compound',)
    elif organism == 'ec':
        raise NotImplementedError('Didn\'t implement EC pathways yet')
    else:   # this is an organism-specific pathway
        entriestype = ('gene', 'compound', 'map')
    pathway.title = tree.getroot().get('title')
    pathway.name = tree.getroot().get('name')
    pathway.id = tree.getroot().get('id')
    #print(entriestype)
    #print(pathway.title)
    #print(pathway.name)
    #print(pathway.id)
    for entry in tree.getiterator('entry'):
        logging.debug(entry.get('type') + ' ' + entry.get('id'))
        node_type = entry.get('type')
        name = entry.get('name')
        node_id = entry.get('id')
#        if nodes.has_key(id):
#                raise TypeError('over writing a key')
        #print(node_type,name,node_id)  # can be ('gene', 'compound', 'map'..)
        graphics = entry.find('graphics')
        node_title = graphics.get('name')
        node_x = int(graphics.get('x'))  # Storing the original X and Y to recreate KEGG layout
        node_y = int(graphics.get('y'))
        logging.debug(node_title)
        nodes[node_id] = (name,node_title,node_type)
        pathway.labels[node_id] = node_title
        pathway.add_node(node_id,data={'label': node_title, 'type': node_type, 'xy': (node_x, node_y)})
        print((node_x,node_y))
       # pathway.add_node(node_id, data={'label': node_title, 'type': node_type, 'xy': (node_x, node_y)})
        #print(name,node_title,node_type)
    for rel in tree.getiterator('relation'):
        e1 = rel.get('entry1')
        e2 = rel.get('entry2')
        pathway.add_edge(e1, e2)
        pathway.relations[e1+'_'+e2] = rel
    for reaction in tree.getiterator('reaction'):
        Id = reaction.get('name')
        substrates = []
        products = []
        for sub in reaction.getiterator('substrate'):
            substrates.append(sub.get('name'))

        for prod in reaction.getiterator('product'):
            products.append(sub.get('name'))

        pathway.reactions[Id] = {'reaction': reaction, 'substrates': substrates, 'products': products}
    print(nodes)
    print(pathway)
    return tree, pathway, nodes, genes

def plot_starlike(pathway):
    pylab.figure()
    networkx.draw_circular(pathway, labels=pathway.labels)
    pylab.title(pathway.title)
    title = pathway.title.replace('/', '-') # TODO: which is the proper way to remove / in a filename?
    pylab.show()

def search_Route(query):
    if not isinstance(query,str):
        query = str(query)
    DFS_result = list(networkx.dfs_edges(pathway,source=query))
    print(DFS_result)
    Tabs = []
    List = dict()
    for idx,itm in enumerate(DFS_result):
        Tabs.append(itm)
        #print(Tab)
        if itm[0] == '46':
            print(itm,Tabs)
            List[idx] = Tabs[:-1]
            Tabs = [itm[1]]
    return List
tree, pathway, nodes, genes = KGML2Graph('ko00950.xml')
T = search_Route(46)
#plot_starlike(pathway)    
#plot_starlike(pathway.get_genes()) 