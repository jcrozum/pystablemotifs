import PyBoolNet
import re
# Convert rules from BooleanNet format to PyBoolNet format
def booleannet2bnet(rules):
    """
    Converts BooleanNet rules to BNet format.
    e.g., an input of
    "A*=B or C and not D"
    returns
    A,  B | C & !D

    Also replaces ~ with !
    """
    s = re.sub("\s*\*\s*=\s*",",\t",rules)
    s = re.sub("\s+not\s+"," !",s, flags=re.IGNORECASE)
    s = re.sub("\(\s*not\s+","(!",s, flags=re.IGNORECASE)
    s = re.sub("\s*~\s*"," !",s, flags=re.IGNORECASE)
    s = re.sub("\s+and\s+"," & ",s, flags=re.IGNORECASE)
    s = re.sub("\s+or\s+"," | ",s, flags=re.IGNORECASE)

    return s

def remove_comment_lines(stream):
    """
    Removes commented out lines, i.e., those starting with '#'
    """
    lines = list(stream)
    lines = filter(lambda x: not x.startswith("#"), lines)
    rules = "".join(lines)
    return rules

def import_primes(fname, format='BooleanNet', remove_constants=False):
    # TODO: add more formats
    rules = remove_comment_lines(open(fname))
    if format == 'BooleanNet':
        rules = booleannet2bnet(rules)
    elif format == 'BNet':
        rules = rules
    else:
        raise ValueError('Unrecognized format',format)

    primes = PyBoolNet.FileExchange.bnet2primes(rules)

    if remove_constants:
        PyBoolNet.PrimeImplicants._percolation(primes,True)

    return primes

def statestring2dict(statestring,names):
    """
    Converts a state string, which specifies a node in an STG, to the corresponding dictionary representation.
    Inputs:
    statestring - a binary string, e.g., '01101'
    names - an ordered list of variable names; (alphabetical order is PyBoolNet's default, e.g. sorted(primes))
    Outpus:
    sd - a dictionary with key set given by names and values given by the corresponding value in statestring
    """
    sd = {}
    for i,c in enumerate(statestring):
        sd[names[i]]=int(c)
    return sd

def statedict2str(statedict):
    """
    Converts a state dictionary to a statestring using default sorting (see statestring2dict)
    """
    return ''.join([str(statedict[x]) for x in sorted(statedict)])

def implicant2bnet(sd):
    """
    Converts a PyBoolNet-formatted prime implicant dictionary (sd) to a BNet string
    e.g., {'A':1,'B':0} returns 'A & !B'
    """
    return ' & '.join(["!"+k for k in sd if not sd[k]]+[k for k in sd if sd[k]])

def pretty_print_primes(primes):
    """
    Prints PyBoolNet a prime dictionary in a more readable format
    """
    for k,v in primes.items():
        for p in v[0]: print(p,"=> !"+k)
        for p in v[1]: print(p,"=>  "+k)

def pretty_print_prime_rules(primes):
    """
    Prints PyBoolNet a prime dictionary as Boolean rules
    The output format is of the form:
    A* = B & C | !D, for example.
    Assumes the PyBoolNet default of disjunctive normal form.
    """
    for k,v in primes.items():
        s = k + "* = "
        sl = []
        for c in v[1]:
            sll = []
            for kk,vv in c.items():
                if vv: sli = kk
                else: sli = '!'+kk
                sll.append(sli)
            if len(sll) > 0:
                sl.append(' & '.join(sll))
        if len(sl) > 0:
            s += ' | '.join(sl)
        print(s)

def pretty_print_rspace(L,simplify=True,silent=True):
    """
    Returns a string representation of the Boolean rule describing the input
    rspace L (see RestrictSpace.rspace). The rule is simplified if simplify=True,
    and the string is printed to screen if silent=False.
    """
    u = []
    for x in L:
        t = []
        for y in x:
            s = []
            for k,v in y.items():
                if v: s.append(k)
                else: s.append('!'+k)
            t.append(' & '.join(s))
        u.append('( '+' | '.join(t)+' )')
    s=' & '.join(u)
    if simplify:
        s = PyBoolNet.BooleanLogic.minimize_espresso(s)
    if not silent:
        print(s)

    return s
