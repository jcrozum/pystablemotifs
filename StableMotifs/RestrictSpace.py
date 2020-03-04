import PyBoolNet
import re

from StableMotifs.DomainOfInfluence import single_drivers, logical_domain_of_influence

# def exclusion_space(fspace):
#     """
#     NOTE: This function is commented out because it provides only extraneous information and is slow
#     It was used to negate conserved function to create a new conserved function, but this is superfluous.
#
#     Returns a list of PyBoolNet prime implicant dictionaries that is satisfied
#     if and only if all rules in fspace are violated.
#
#     Input:
#     fspace - a list of prime implicant dictionaries, to be interpeted as a disjunctive normal form function
#
#     Output:
#     p2 - the negation of fspace, provided in the same form
#     """
#     if len(fspace) == 0: return [{}]
#     # Build the rule; TRUE means we are in at least one of the spaces
#     # So we will want the implicants for the FALSE value
#     nrule = '__reserved,\t'
#     rule_list = []
#     for d in fspace:
#         clause_list = []
#         for k in d:
#             if d[k]: clause_list.append(k)
#             else: clause_list.append('!'+k)
#         rule_list.append('&'.join(clause_list))
#     nrule += '|'.join(rule_list)
#
#     # Build the implicants, grab the ones for the FALSE value
#     p2 = PyBoolNet.FileExchange.bnet2primes(nrule)['__reserved'][0]
#
#     return p2

def attractor_space_candidates(maxts,trmaxts):
    """
    Merge the maximum trap spaces maxts and time-reverse maximum trap spaces
    to obtain a list of attractor-conserved quantities.
    """
    #L=[exclusion_space(maxts+trmaxts)] # use exclusion space to include negation of conserved quantities
    L = []
    for t in maxts+trmaxts: L.append([t])
    return L

def rspace(maxts,primes):
    """
    In order for none of the trap spaces to "lock in", we would require that
    their single-node drivers are all sustained in a negated state. We can use
    this idea to hone the exclusion space. rspace will return the region that
    1) does not have any maxts active and
    2) has the negations of 1-node drivers of each maxts active and
    3) has the update rules of these 1-node drivers taking the appropriate value
    The return value is a list L of lists of prime implicants. Each element of L
    is to be interpreted as a list of OR-separated prime implicants; L is to be
    interpreted as AND-separated. e.g.,
    L=[[{'A':0,'B':1},{C:0}],[{'B':0,'D':1},{'A':1}]]
    should be read as L = (!A&B | !C) & (!B & D | A)
    """
    nds = {}
    for ts in maxts:
        drivers = single_drivers(ts,primes)
        for d in drivers:
            node = list(d.keys())[0]
            if node in nds:
                if nds[node] == d[node]:
                    return [[{'0':1}]] # In this case, at least one SM MUST stabilize
            nds[node] = int(not d[node])
    implied,contradicted = logical_domain_of_influence(nds,primes)
    if len(contradicted) != 0:
        return [[{'0':1}]] # In this case, at least one SM MUST stabilize

    for k in implied: nds[k]=implied[k]

    L = [[nds]]
    for k in nds:
        L.append(primes[k][nds[k]])
    return L

def fixed_rspace_nodes(L,primes):
    """
    Finds the nodes that must have a fixed value in order for the rspace
    constraint L to be satisfied in the system given by primes.
    Relies on the special structure of L.
    Returns {'0':1} if reduction uncovers that L is a contradiction.
    """
    fL = [r[0] for r in L if len(r)==1]

    fd = {}
    for d in fL:
        for k,v in d.items():
            if k in fd and fd[k] != v:
                return {'0':1}
            else:
                fd[k] = v

    imp,con = logical_domain_of_influence(fd,primes)

    if len(con) > 0:
        return {'0':1}
    fd.update(imp)
    return fd

def reduce_rspace_string(s,fd):
    """
    Replaces variables in the string s with the fixed values given by the dictionary fd.
    """
    s2 = s
    for k,v in fd.items():
        s2 = re.sub(rf"\b{k}\b",str(v),s2)
    s2 = PyBoolNet.BooleanLogic.minimize_espresso(s2)
    return s2
