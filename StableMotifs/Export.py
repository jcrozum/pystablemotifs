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

def networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=True):
    """Label the succesion diagram and (optionally) attractors of the input attractor
    repertoire according to the conventions of Rozum et al. (2020). Useful for
    plotting.

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

    G_reduced_network_based=ar.succession_diagram.digraph.copy()
    has_nodes = False
    for i in G_reduced_network_based.nodes():
        has_nodes = True
        G_reduced_network_based.nodes[i]['label']=format_reduction_label(str(ar.succession_diagram.motif_reduction_dict[i].motif_history))
        G_reduced_network_based.nodes[i]['virtual_nodes']=ar.succession_diagram.motif_reduction_dict[i].motif_history

    if not has_nodes:
        G_reduced_network_based.add_node(0)
        G_reduced_network_based.nodes[0]['label'] = '[]'
        G_reduced_network_based.nodes[0]['virtual_nodes'] = []

    if include_attractors_in_diagram:
        for a_index,a in enumerate(ar.attractors):
            G_reduced_network_based.add_node('A'+str(a_index))
            G_reduced_network_based.nodes['A'+str(a_index)]['label']=format_reduction_label(str(a.attractor_dict))
            G_reduced_network_based.nodes['A'+str(a_index)]['virtual_nodes']=a.attractor_dict

            for r in a.reductions:
                r_key=list(ar.succession_diagram.motif_reduction_dict.keys())[list(ar.succession_diagram.motif_reduction_dict.values()).index(r)]
                G_reduced_network_based.add_edge(r_key,'A'+str(a_index))

    return G_reduced_network_based

def plot_nx_succession_diagram(g, fig_dimensions=[], pos='pydot', detailed_labels=True, node_size=[], node_color='grey',
                              font_size=12, font_color='black'):
    """Plot the input succession diagram. Requires matplotlib. For finer control
    over plot appearance, it is recommended to plot g directly.

    Parameters
    ----------
    g : networkx.DiGraph
        Labeled succession diagram, e.g., as is output from
        Export.networkx_succession_diagram_reduced_network_based().
    fig_dimensions : (int,int)
        Dimensions of the output figure. If [], then the dimensions are calculated
        based on the number of nodes in g (the default is []).
    pos : str or graphviz_layout
        Layout for the node labels; if 'pydot', these are automatically computed
        using pydot (requires pydot) (the default is 'pydot').
    detailed_labels : bool
        Whether node labels should be drawn (True) or left as metadata (False)
        (the default is True).
    node_size : int
        Size of the nodes. If [], the size is proportional to the number of nodes
        in g (the default is []).
    node_color : color or array of colors
        Node color. Can be a single color or a sequence of colors with the same
        length as nodelist. (the default is 'grey').
    font_size : int
        Font size for labels (the default is 12).
    font_color : str
        Color for label text (the default is 'black').

    """

    from networkx.drawing.nx_agraph import graphviz_layout
    import matplotlib.pyplot as plt

    if fig_dimensions==[]:
        fig_dimensions=(2*(g.number_of_nodes()+2),g.number_of_nodes()+2)
    if node_size==[]:
        node_size=50*g.number_of_nodes()

    if pos=='pydot':
        pos=graphviz_layout(g, prog='dot')

    plt.figure(figsize=fig_dimensions)
    nx.drawing.draw_networkx_nodes(g, pos,node_shape='s',node_color=node_color, node_size=node_size)
    nx.draw_networkx_edges(g, pos, arrowstyle='fancy',arrowsize=10)
    if detailed_labels:
        nx.drawing.draw_networkx_labels(g,pos, labels=dict(g.nodes('label')),font_size=font_size, font_color=font_color)
    else:
        nx.drawing.draw_networkx_labels(g,pos,font_size=font_size, font_color=font_color)
    plt.axis('off')
    plt.show()

def networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=True):
    """Label the succesion diagram and (optionally) attractors of the input attractor
    repertoire according to the conventions of Zanudo and Albert (2015). Useful
    for plotting.

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
    G_reduced_network_based=networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=False)
    G_motif_based = nx.line_graph(G_reduced_network_based)
    for i,j in G_motif_based.nodes():


        node_motif=set([frozenset(k.items()) for k in ar.succession_diagram.motif_reduction_dict[j].motif_history])-set([frozenset(k.items()) for k in ar.succession_diagram.motif_reduction_dict[i].motif_history])
        node_label=format_reduction_label(str(dict(list(node_motif)[0])))
        G_motif_based.nodes[(i,j)]['label']=node_label
        G_motif_based.nodes[(i,j)]['virtual_nodes']=dict(list(node_motif)[0])

    if include_attractors_in_diagram:
        for a_index,a in enumerate(ar.attractors):
            G_motif_based.add_node('A'+str(a_index))
            G_motif_based.nodes['A'+str(a_index)]['label']=format_reduction_label(str(a.attractor_dict))
            G_motif_based.nodes['A'+str(a_index)]['virtual_nodes']=a.attractor_dict
            for r in a.reductions:
                r_key=list(ar.succession_diagram.motif_reduction_dict.keys())[list(ar.succession_diagram.motif_reduction_dict.values()).index(r)]
                for n in G_motif_based.nodes():
                    if type(n)==tuple:
                        i,j=n
                        if r_key==j:
                            G_motif_based.add_edge((i,j),'A'+str(a_index))
    return G_motif_based

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

def get_motif_set(ar):
    """Extract the stable motifs of a system and its reduced networks from its
    attractor repertoire. Notably, these include both the system's primary stable
    motifs and conditionally stable motifs (see, e.g., Deritei et al. 2019).

    Parameters
    ----------
    ar : AttractorRepertoire
        The attractor repertoire of the system.

    Returns
    -------
    list of dictionaries
        Stable motifs that appear in the system or during reduction.

    """
    GM_no_attr=networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=False)
    SM_set=set([])
    for n in GM_no_attr.nodes(data=True):
        SM_set.add(frozenset([i for i in n[1]['virtual_nodes'].items()]))
    return [dict(sm) for sm in SM_set]

def save_to_graphml(G,model_name):
    """Export a labeled succesion diagram to graphml format.

    Parameters
    ----------
    G : networkx.DiGraph
        Labeled succession diagram to export.
    model_name : str
        Name of file to save to (.graphml extension will be appended).

    """

    #Graphml does not support complex attribues so we create a copy with the virtual_nodes attribute:
    G_ex=G.copy()
    for n in G_ex.nodes():
        G_ex.nodes[n]['virtual_nodes']=str(G_ex.nodes[n]['virtual_nodes'])
    nx.write_graphml(G_ex, "%s.graphml"%model_name)
