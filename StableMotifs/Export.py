import StableMotifs as sm
import networkx as nx

def format_reduction_label(s):
    return s.replace("'","").replace('[','').replace(']','')

def networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=True):

    '''
    tbd
    '''

    G_reduced_network_based=ar.succession_diagram.digraph.copy()
    has_nodes = False
    for i in G_reduced_network_based.nodes():
        has_nodes = True
        G_reduced_network_based.nodes[i]['label']=format_reduction_label(str(ar.succession_diagram.motif_reduction_dict[i].motif_history))
        G_reduced_network_based.nodes[i]['node_states']=ar.succession_diagram.motif_reduction_dict[i].motif_history

    if not has_nodes:
        G_reduced_network_based.add_node(0)
        G_reduced_network_based.nodes[0]['label'] = '[]'
        G_reduced_network_based.nodes[0]['node_states'] = []

    if include_attractors_in_diagram:
        for a_index,a in enumerate(ar.attractors):
            G_reduced_network_based.add_node('A'+str(a_index))
            G_reduced_network_based.nodes['A'+str(a_index)]['label']=format_reduction_label(str(a.attractor_dict))
            G_reduced_network_based.nodes['A'+str(a_index)]['node_states']=a.attractor_dict

            for r in a.reductions:
                r_key=list(ar.succession_diagram.motif_reduction_dict.keys())[list(ar.succession_diagram.motif_reduction_dict.values()).index(r)]
                G_reduced_network_based.add_edge(r_key,'A'+str(a_index))

    return G_reduced_network_based

def plot_nx_succession_diagram(g, fig_dimensions=[], pos='pydot', detailed_labels=True, node_size=[], node_color='grey',
                              font_size=12, font_color='black'):

    '''
    tbd
    '''
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
    '''
    tbd
    '''
    G_reduced_network_based=networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=False)
    G_motif_based = nx.line_graph(G_reduced_network_based)
    for i,j in G_motif_based.nodes():


        node_motif=set([frozenset(k.items()) for k in ar.succession_diagram.motif_reduction_dict[j].motif_history])-set([frozenset(k.items()) for k in ar.succession_diagram.motif_reduction_dict[i].motif_history])
        node_label=format_reduction_label(str(dict(list(node_motif)[0])))
        G_motif_based.nodes[(i,j)]['label']=node_label
        G_motif_based.nodes[(i,j)]['node_states']=dict(list(node_motif)[0])

    if include_attractors_in_diagram:
        for a_index,a in enumerate(ar.attractors):
            G_motif_based.add_node('A'+str(a_index))
            G_motif_based.nodes['A'+str(a_index)]['label']=format_reduction_label(str(a.attractor_dict))
            G_motif_based.nodes['A'+str(a_index)]['node_states']=a.attractor_dict
            for r in a.reductions:
                r_key=list(ar.succession_diagram.motif_reduction_dict.keys())[list(ar.succession_diagram.motif_reduction_dict.values()).index(r)]
                for i,j in G_motif_based.nodes():
                    if r_key==j:
                        G_motif_based.add_edge((i,j),'A'+str(a_index))
    return G_motif_based

def attractor_dataframe(ar):
    '''
    tbd
    '''
    import pandas as pd
    df=pd.DataFrame()
    for a in ar.attractors:
        df=df.append(a.attractor_dict,ignore_index=True).astype(int, errors='ignore').astype(str)

    return df

def get_motif_set(ar):
    '''
    tbd
    '''
    GM_no_attr=networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=False)
    SM_set=set([])
    for n in GM_no_attr.nodes(data=True):
        SM_set.add(frozenset([i for i in n[1]['node_states'].items()]))
    return [dict(sm) for sm in SM_set]

def save_to_graphml(G,model_name):
    '''
    tbd
    '''
    #Graphml does not support complex attribues so we create a copy with the node_states attribute:
    G_ex=G.copy()
    for n in G_ex.nodes():
        G_ex.nodes[n]['node_states']=str(G_ex.nodes[n]['node_states'])
    nx.write_graphml(G_ex, "%s.graphml"%model_name)
