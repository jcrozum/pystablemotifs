import StableMotifs as sm
import networkx as nx

def format_reduction_label(s):
    """Helper function to make graph labels more readable. Removes single quotes
    and square brackets from the input string.

    Parameters
    ----------
    s : str
        Input label.

    Returns
    -------
    str
        Label with characters ', [, and ] removed.

    """
    return s.replace("'","").replace('[','').replace(']','')

def expanded_network(primes, single_parent_composites = False):
    """Produce the expanded network for given input update rules.

    Parameters
    ----------
    primes : PyBoolNet primes dictionary
        The update rules for which to construct the expanded network.
    single_parent_composites : bool
        Whether to insert composite nodes between virtual nodes when one is a prime
        implicant of the other. If False, the number of nodes is decreased; if
        True, then the expanded network is bipartite (the default is False).

    Returns
    -------
    networkx.DiGraph
        Digraph representing the expanded network. Nodes have a 'type' attribute
        that can be either 'virtual' or 'composite'.

    """
    G = nx.DiGraph()
    cnode_id = 0
    for p in primes:
        for v in [0,1]:
            name = '('+str(p)+','+str(v)+')'
            G.add_node(name)
            G.nodes[name]['label'] = name
            G.nodes[name]['type'] = 'virtual'

            for hedge in primes[p][v]:
                G.add_node(cnode_id)
                G.nodes[cnode_id]['type'] = 'composite'
                G.add_edge(cnode_id,name)

                for k in hedge:
                    parent = '(' + str(k) + ',' + str(hedge[k]) + ')'
                    G.add_edge(parent,cnode_id)
                cnode_id += 1

    # If we want to remove composite nodes of "size" one
    if not single_parent_composites:
        for i in range(cnode_id):
            if G.in_degree(i) == 1:
                pre = list(G.predecessors(i))[0]
                suc = G.successors(i)
                for j in suc:
                    G.add_edge(pre,j)
                G.remove_node(i)


    return G

def networkx_succession_diagram(ar,include_attractors_in_diagram=True,use_compressed_diagram=True):
    """Label the succesion diagram and (optionally) attractors of the input attractor
    repertoire according to the conventions of Rozum et al. (2021). This is an
    alias for the function Export.networkx_succession_diagram_reduced_network_based.

    Parameters
    ----------
    ar : AttractorRepertoire
        Attractor repertoire object for which to build the diagram.
    include_attractors_in_diagram : bool
        Whether attractors should be represented as nodes in the diagram (the
        default is True).
    use_compressed_diagram : bool
        Whether to use the (potentially compressed) succession diagram stored in
        ar.succession_digraph instead of the complete one ar.succession_diagram.digraph.
        These are equivelent unless ar.simplify_diagram is called. See
        AttractorRepertoire.py for additional details. The default is True.

    Returns
    -------
    networkx.DiGraph
        A labeled digraph that represents the succession diagram.

    """
    return networkx_succession_diagram_reduced_network_based(ar,
        include_attractors_in_diagram=include_attractors_in_diagram,
        use_compressed_diagram=use_compressed_diagram)

def networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=True,use_compressed_diagram=True):
    """Label the succesion diagram and (optionally) attractors of the input attractor
    repertoire according to the conventions of Rozum et al. (2021).

    Parameters
    ----------
    ar : AttractorRepertoire
        Attractor repertoire object for which to build the diagram.
    include_attractors_in_diagram : bool
        Whether attractors should be represented as nodes in the diagram (the
        default is True).
    use_compressed_diagram : bool
        Whether to use the (potentially compressed) succession diagram stored in
        ar.succession_digraph instead of the complete one ar.succession_diagram.digraph.
        These are equivelent unless ar.simplify_diagram is called. See
        AttractorRepertoire.py for additional details. The default is True.

    Returns
    -------
    networkx.DiGraph
        A labeled digraph that represents the succession diagram.

    """

    if use_compressed_diagram:
        G_reduced_network_based=ar.succession_digraph.copy()
        relevant_nodes = ar.relevant_nodes
    else:
        G_reduced_network_based=ar.succession_diagram.digraph.copy()
        relevant_nodes = ar.primes.keys()

    has_nodes = False
    for i in G_reduced_network_based.nodes():
        has_nodes = True
        states={}
        for k,v in ar.succession_diagram.motif_reduction_dict[i].logically_fixed_nodes.items():
            if k in relevant_nodes:
                states[k]=v
        G_reduced_network_based.nodes[i]['index'] = i
        G_reduced_network_based.nodes[i]['states']=states

        contraction_indices = [i]
        if use_compressed_diagram and 'contraction' in ar.succession_digraph.nodes[i].keys():
            contraction_indices += list(ar.succession_digraph.nodes[i]['contraction'].keys())

        histories = []
        for c in contraction_indices:
            histories.append(ar.succession_diagram.motif_reduction_dict[c].motif_history)
        G_reduced_network_based.nodes[i]['history']=histories

    if not has_nodes:
        G_reduced_network_based.add_node(0)
        G_reduced_network_based.nodes[0]['index'] = 0
        G_reduced_network_based.nodes[0]['states'] = {}
        G_reduced_network_based.nodes[0]['history'] = [[]]

    for u,v in G_reduced_network_based.edges():
        ufixed = G_reduced_network_based.nodes[u]['states']
        vfixed = G_reduced_network_based.nodes[v]['states']
        lock_in = {}
        for k,val in vfixed.items():
            if k not in ufixed.keys():
                lock_in[k] = val
        G_reduced_network_based.edges[u,v]['states'] = lock_in

    if include_attractors_in_diagram and not use_compressed_diagram:
        for a_index,a in enumerate(ar.attractors):
            astr = 'A'+str(a_index)
            G_reduced_network_based.add_node(astr)
            G_reduced_network_based.nodes[astr]['index']=astr
            G_reduced_network_based.nodes[astr]['states']=a.attractor_dict

            for r in a.reductions:
                r_ind = list(ar.succession_diagram.motif_reduction_dict.values()).index(r)
                r_key=list(ar.succession_diagram.motif_reduction_dict.keys())[r_ind]
                if r_key in G_reduced_network_based.nodes():
                    G_reduced_network_based.add_edge(r_key,astr,states='')
    if include_attractors_in_diagram and use_compressed_diagram:
        for a_index,a in enumerate(ar.attractor_equivalence_classes):
            astr = 'A'+str(a_index)
            G_reduced_network_based.add_node(astr)
            G_reduced_network_based.nodes[astr]['index']=astr
            G_reduced_network_based.nodes[astr]['states']=a['states']

            for r_key in a['reductions']:
                if r_key in G_reduced_network_based.nodes():
                    G_reduced_network_based.add_edge(r_key,astr,states='')

    for n in G_reduced_network_based.nodes():

        if not use_compressed_diagram:
            if str(n)[0]=='A':
                G_reduced_network_based.nodes[n]['label']=str(G_reduced_network_based.nodes[n]['states'])
            else:
                G_reduced_network_based.nodes[n]['label']=str(G_reduced_network_based.nodes[n]['history'])
        else:
            if str(n)[0]=='A':
                G_reduced_network_based.nodes[n]['label']=str({k:v for k,v in G_reduced_network_based.nodes[n]['states'].items() if k in ar.relevant_nodes})
            else:
                G_reduced_network_based.nodes[n]['label']=str([[{k:v for k,v in d.items() if k in ar.relevant_nodes} for d in h] for h in G_reduced_network_based.nodes[n]['history']])

    return G_reduced_network_based

def networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=True):
    """Label the succesion diagram and (optionally) attractors of the input attractor
    repertoire according to the conventions of Zanudo and Albert (2015). If
    attractors are not included, this is the line graph of the succession diagram
    defined in Rozum et al. (2021). Does not support compression.

    Parameters
    ----------
    ar : AttractorRepertoire
        Attractor repertoire object for which to build the diagram.
    include_attractors_in_diagram : bool
        Whether attractors should be represented as nodes in the diagram (the
        default is True).

    Returns
    -------
    networkx.DiGraph
        A labeled digraph that represents the succession diagram.

    """
    G_reduced_network_based = networkx_succession_diagram_reduced_network_based(ar,
        include_attractors_in_diagram=False,
        use_compressed_diagram=False)

    G_motif_based = nx.line_graph(G_reduced_network_based)

    for e in G_reduced_network_based.edges():
        G_motif_based.nodes[e]['states'] = G_reduced_network_based.edges[e]['states']
        G_motif_based.nodes[e]['label'] = str(G_reduced_network_based.edges[e]['states'])
    for e in G_motif_based.edges():
        G_motif_based.edges[e]['states']=''

    if include_attractors_in_diagram:
        for a_index,a in enumerate(ar.attractors):
            astr = 'A'+str(a_index)
            G_motif_based.add_node(astr)
            G_motif_based.nodes[astr]['states']=a.attractor_dict
            for r in a.reductions:
                r_ind = list(ar.succession_diagram.motif_reduction_dict.values()).index(r)
                r_key=list(ar.succession_diagram.motif_reduction_dict.keys())[r_ind]
                for n in G_motif_based.nodes():
                    if type(n)==tuple:
                        i,j=n
                        if r_key==j:
                            G_motif_based.add_edge((i,j),astr,states='')
    return G_motif_based

def plot_nx_succession_diagram(G, pos=None, fig_dimensions=(None,None), nx_node_kwargs=None, nx_edge_kwargs=None,
    draw_node_labels=True, labeling_convention='label', draw_edge_labels=False, nx_node_label_kwargs=None, nx_edge_label_kwargs=None):
    """Plot the input succession diagram. Requires matplotlib. For finer control
    over plot appearance, it is recommended to plot g directly.

    Parameters
    ----------
    g : networkx.DiGraph
        Labeled succession diagram, e.g., as is output from
        Export.networkx_succession_diagram_reduced_network_based().
    fig_dimensions : (int,int)
        Dimensions of the output figure. If (None,None), then the dimensions are
        calculated based on the number of nodes in g (the default is (None,None)).
    pos : str or graphviz_layout
        Layout for the nodes; A dictionary with nodes as keys and positions as
        values. Positions should be sequences of length 2. If none, we attempt to
        use pydot/graphviz to construct a layout, otherwise we fall back to the
        networkx planar_layout function (succession diagrams are always planar).
    draw_node_labels : bool
        Whether node labels should be drawn (True) or left as metadata (False)
        (the default is True).
    draw_edge_labels : bool
        Whether edge labels should be drawn (True) or left as metadata (False)
        (the default is True).
    labeling_convention : str
        Whether edge labels should be just the stable motifs ('label') or all stabilized states ('states')
        (the default is 'label').
    nx_node_kwargs : dictionary
        Keword arguments passed to nx.draw_networkx_nodes (in addition to G and pos).
        If None, we pass {'node_size':50*G.number_of_nodes()} by default.
    nx_edge_kwargs : dictionary
        Keword arguments passed to nx.draw_networkx_edges (in addition to G and pos).
        If None, we pass {'arrowstyle':'-|>','width':2,'arrowsize':30} by default.
    nx_node_label_kwargs : dictionary
        Keword arguments passed to nx.draw_networkx_labels (in addition to G and pos).
        If None, we pass {'font_size':16} by default.
    nx_edge_label_kwargs : dictionary
        Keword arguments passed to nx.draw_networkx_edge_labels (in addition to G and pos).
        If None, we pass {'font_size':16} by default.

    """
    import matplotlib.pyplot as plt

    if fig_dimensions == (None,None):
        fig_dimensions=(2*(G.number_of_nodes()+2),G.number_of_nodes()+2)

    if pos is None:
        try:
            from networkx.drawing.nx_agraph import graphviz_layout
            pos = graphviz_layout(G, prog='dot')
        except ImportError:
            pos = nx.planar_layout(G)

    plt.figure(figsize=fig_dimensions)

    if nx_node_kwargs is None:
        nx_node_kwargs = {'node_size':50*G.number_of_nodes()}
    if nx_edge_kwargs is None:
        nx_edge_kwargs = {'arrowstyle':'-|>','width':2,'arrowsize':30}
    if nx_node_label_kwargs is None:
        nx_node_label_kwargs = {'font_size':16}
    if nx_edge_label_kwargs is None:
        nx_edge_label_kwargs = {'font_size':16}

    nx.drawing.draw_networkx_nodes(G, pos,**nx_node_kwargs)
    nx.draw_networkx_edges(G, pos,**nx_edge_kwargs)
    if draw_node_labels:
        if labeling_convention=='label':
            nx.drawing.draw_networkx_labels(G,pos, labels=dict(G.nodes('label')),**nx_node_label_kwargs)
        else:
            nx.drawing.draw_networkx_labels(G,pos, labels=dict(G.nodes('states')),**nx_node_label_kwargs)
    if draw_edge_labels:
        nx.drawing.draw_networkx_edge_labels(G,pos,labels=dict(G.edges('states')),**nx_edge_label_kwargs)
    plt.axis('off')
    plt.show()

# The next two methods are commented out because they are experimental methods
# for compressing the succession diagram and are under development.

# def networkx_succession_diagram_motif_based_simplified(ar, GM=None, include_attractors_in_diagram=True):
#     """Produce a compressed version of the succession diagram using the conventions
#     of Zanudo and Albert (2015). In this compression, stable motifs (potentially
#     ones that are only stable after some number of reductions) are represented
#     only once. If a single motif appears more than once in the succession diagram,
#     those nodes are merged.
#
#     Parameters
#     ----------
#     ar : AttractorRepertoire
#         Attractor repertoire object for which to build the diagram.
#     GM : networkx.DiGraph
#         Labeled motif-based succession diagram to simplify. If None, the diagram
#         is generated from ar (the defaule is None).
#     include_attractors_in_diagram : bool
#         Whether attractors should be represented as nodes in the diagram (the
#         default is True).
#
#     Returns
#     -------
#     networkx.DiGraph
#         Simplified motif-based succession diagram.
#
#     """
#     if GM==None:
#         GM=networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
#     motifs_list=get_motif_set(ar)
#     motifs_dict = dict(zip(range(len(motifs_list)),motifs_list))
#     merged_dict=motifs_dict.copy()
#     attractors_dict=dict()
#     if include_attractors_in_diagram:
#         for a_index,a in enumerate(ar.attractors):
#             merged_dict['A'+str(a_index)]=a.attractor_dict
#             attractors_dict['A'+str(a_index)]=a.attractor_dict
#     motif_keys=list(merged_dict.keys())
#     motif_values=list(merged_dict.values())
#     GMM=nx.DiGraph()
#     GMM.add_nodes_from(merged_dict)
#     for n in GMM.nodes():
#         GMM.nodes[n]['virtual_nodes']=merged_dict[n]
#     for i,j in GM.edges():
#         print(i,j)
#         source=motif_keys[motif_values.index(GM.nodes[i]['virtual_nodes'])]
#         target=motif_keys[motif_values.index(GM.nodes[j]['virtual_nodes'])]
#         GMM.add_edge(source,target)
#     return GMM
#
# def networkx_motif_attractor_bipartite_graph(ar):
#     """Produce a motif-attractor bipartite compression of the succession diagram
#     for the given attractor repertoire. There is an edge from a (conditionally)
#     stable motif to an attractor if the motif is compatible with that attractor.
#
#     Parameters
#     ----------
#     ar : AttractorRepertoire
#         Attractor repertoire object for which to build the diagram.
#
#     Returns
#     -------
#     networkx.DiGraph
#         Motif-attractor bipartite condensed succession diagram.
#
#     """
#
#     GMM=networkx_succession_diagram_motif_based_simplified(ar,include_attractors_in_diagram=True)
#
#     motifs_list=get_motif_set(ar)
#     motifs_dict = dict(zip(range(len(motifs_list)),motifs_list))
#
#     attractors_dict=dict()
#     for a_index,a in enumerate(ar.attractors):
#         attractors_dict['A'+str(a_index)]=a.attractor_dict
#
#     GM_bp=nx.DiGraph()
#     GM_bp.add_nodes_from(motifs_dict)
#     GM_bp.add_nodes_from(attractors_dict)
#     for m in motifs_dict:
#         GM_bp.nodes[m]['virtual_nodes']=motifs_dict[m]
#         for a in attractors_dict:
#             GM_bp.nodes[a]['virtual_nodes']=attractors_dict[a]
#             if nx.has_path(GMM,m,a):
#                 GM_bp.add_edge(m,a)
#
#     for a in attractors_dict:
#             GM_bp.nodes[a]['virtual_nodes']=attractors_dict[a]
#
#     return GM_bp

def attractor_dataframe(ar):
    """Summarize the input attractor repertoire in a pandas DataFrame (requires
    pandas).

    Parameters
    ----------
    ar : AttractorRepertoire
        Attractor repertoire to summarize.

    Returns
    -------
    pandas.DataFrame
        Summary of the attractors.

    """
    import pandas as pd
    df=pd.DataFrame()
    for a in ar.attractors:
        df=df.append(a.attractor_dict,ignore_index=True).astype(int, errors='ignore').astype(str)

    return df

def save_to_graphml(G,model_name):
    """Export a labeled succesion diagram to graphml format.

    Parameters
    ----------
    G : networkx.DiGraph
        Labeled succession diagram to export.
    model_name : str
        Name of file to save to (.graphml extension will be appended).

    """

    #Graphml does not support complex attribues so we create a copy and stringify everything
    G_ex=G.copy()
    for n in G_ex.nodes():
        for k,v in G_ex.nodes[n].items():
            G_ex.nodes[n][k]=str(v)
    for e in G_ex.edges():
        for k,v in G_ex.edges[e].items():
            G_ex.edges[e][k]=str(v)
    nx.write_graphml(G_ex, "%s.graphml"%model_name)
