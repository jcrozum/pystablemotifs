import PyBoolNet
import re

from StableMotifs.DomainOfInfluence import single_drivers, logical_domain_of_influence
import StableMotifs.DomainOfInfluence as sm_doi

def attractor_space_candidates(maxts,trmaxts):
    """
    Merge the maximum trap spaces maxts and time-reverse maximum trap spaces
    to obtain a list of attractor-conserved quantities.
    """
    #L=[exclusion_space(maxts+trmaxts)] # use exclusion space to include negation of conserved quantities
    L = []
    for t in maxts+trmaxts: L.append([t])
    return L

def state_in_rspace(state,L):
    for clause in L:
        sat = False
        for p in clause:
            if sm_doi.fixed_implies_implicant(state,p):
                sat = True
                break
        if not sat:
            return False
    return True

def partial_state_contradicts_rspace(st,L):
    for clause in L:
        csat = False
        for p in clause:
            if not sm_doi.fixed_excludes_implicant(st,p):
                csat = True
                break
        if not csat:
            return False
    return True

def reduce_rspace(L,primes):
    fixed = fixed_rspace_nodes(L,primes)
    L2 = [[fixed]]
    for clause in L:
        sat = False
        for p in clause:
            if sm_doi.fixed_implies_implicant(fixed,p):
                sat = True
                break
        if not sat:
            L2.append(clause)
    return L2

def rspace(maxts,trmaxts,primes):
    """
    In order for none of the trap spaces to "lock in", we would require that
    their single-node drivers are all sustained in a negated state. We can use
    this idea to hone the exclusion space. rspace will return the region that


    1) has the negations of 1-node drivers of each maxts active and . . .
    2) has the update rules of these 1-node drivers taking the appropriate value

    In addition, a time-reverse trap space (trmaxts) describes a region that,
    once exited, cannot be reentered. Thus, if the LDOI of the region contains
    any contradiciton, the region cannot contain any attractor. Therefore,
    we include a third criterion for the rspace:

    3) is not in a state belonging to an attractor-free time-reversed trap space

    The return value is a list L of lists of prime implicants. Each element of L
    is to be interpreted as a list of OR-separated prime implicants; L is to be
    interpreted as AND-separated. e.g.,
    L=[[{'A':0,'B':1},{C:0}],[{'B':0,'D':1},{'A':1}]]
    should be read as L = (!A&B | !C) & (!B & D | A)
    """
    L = []
    nds = {} # negated 1-node drivers of the stable motifs in maxts
    for ts in maxts:
        # Make sure we aren't in the stable motif ts
        # clause = []
        # for k,v in ts.items():
        #     clause.append({k:int(not v)})
        # L.append(clause)

        # Now consider drivers of ts
        drivers = single_drivers(ts,primes)
        for d in drivers:
            node = list(d.keys())[0] # d only has 1 key because its a 1-node driver set
            if node in nds:
                if nds[node] == d[node]:
                    return [[{'0':1}]] # In this case, at least one SM MUST stabilize
            nds[node] = int(not d[node])
    implied,contradicted = logical_domain_of_influence(nds,primes)
    if len(contradicted) != 0:
        return [[{'0':1}]] # In this case, at least one SM MUST stabilize

    for k in implied: nds[k]=implied[k]

    L.append([nds])

    # also include the conditions for the update rules of the negated drivers
    for k in nds:
        L.append(primes[k][nds[k]])


    # Now, investigate the trmaxts:
    for gs in trmaxts:
        implied,contradicted = logical_domain_of_influence(gs,primes)
        if len(contradicted) > 0:
            L.append([{n:int(not v)} for n,v in gs.items()])

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
