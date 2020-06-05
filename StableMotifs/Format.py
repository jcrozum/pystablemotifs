import PyBoolNet
import re
import subprocess
import os
import ast
import datetime

BASE = os.path.join(os.path.dirname(PyBoolNet.__file__))
config = PyBoolNet.Utility.Misc.myconfigparser.SafeConfigParser()
config.read(os.path.join(BASE, "Dependencies", "settings.cfg"))

CMD_BNET2PRIMES = os.path.normpath(os.path.join(BASE, "Dependencies", config.get("Executables", "bnet2prime")))

def longbnet2primes(BNET, remove_constants = False):
    """
    A modified version of PyBoolNet's bnet2primes that does not do path-checking,
    as this can cause errors if the bnet rules are very long. Assumes BNET is a
    bnet string, not a file.
    """
    cmd = [CMD_BNET2PRIMES]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate( input=BNET.encode() )
    proc.stdin.close()
    PyBoolNet.FileExchange._bnet2primes_error(proc, out, err, cmd)
    out = out.decode()

    out = out.replace('\x08','') # remove backspaces
    out = out.replace(' ','') # remove whitespaces

    primes = ast.literal_eval(out)

    if remove_constants:
        PyBoolNet.PrimeImplicants._percolation(primes,True)

    return primes

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
    s = re.sub("\s*\*\s*=\s*",",\t",rules) # replace "=" with ",\t"
    s = re.sub("\s+not\s+"," !",s, flags=re.IGNORECASE) # not -> !
    s = re.sub("\(\s*not\s+","(!",s, flags=re.IGNORECASE) # not -> ! (with parens)
    s = re.sub("\s*~\s*"," !",s, flags=re.IGNORECASE) # ~ -> !
    s = re.sub("\s+and\s+"," & ",s, flags=re.IGNORECASE) # and -> &
    s = re.sub("\s+or\s+"," | ",s, flags=re.IGNORECASE) # or -> |

    return s

def bnetDNF2list(bnet):
    if bnet == "0":
        return []
    elif bnet == "1":
        return [{}]

    bnetList = []
    bnet_trim = re.sub("[\s\(\)]*","",bnet) # remove all whitespace and parens
    LL = [x.split('&') for x in bnet_trim.split('|')]

    for L in LL:
        Ldict = {}
        contradiction = False
        for literal in L:
            if literal[0]=="!":
                n = literal[1:]
                s = 0
            else:
                n = literal
                s = 1
            if n not in Ldict:
                Ldict[n] = s
            elif Ldict[n] != s:
                contradiction = True
                break
        if not contradiction:
            bnetList.append(Ldict)
    return bnetList

def build_rule_using_bnetDNFs(expr0,expr1):
    return [bnetDNF2list(expr0),bnetDNF2list(expr1)]

def bnet2sympy(rule):
    # print("bnet:",rule)
    crule = re.sub("!","~",rule)
    crule = re.sub("[\b\(]1[\b\)]","(x | ~x)",crule)
    crule = re.sub("[\b\(]0[\b\)]","(x & ~x)",crule)
    crule = re.sub("True","(x | ~x)",crule)
    crule = re.sub("False","(x & ~x)",crule)
    # print("sympy:",crule)
    return crule

def sympy2bnet(rule):
    crule = re.sub("~","!",rule)
    crule = re.sub("True","1",crule)
    crule = re.sub("False","0",crule)
    return crule

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

def statelist2dict(names,c):
    """
    Converts a collection of statestrings to a dictionary. If a node takes the
    same value in every state, the corresponding dictionary value matches its
    fixed value; otherwise, the dictionary value is 'X'.
    """
    d = {}
    for i,n in enumerate(names):
        for cs in c:
            if n not in d:
                d[n] = cs[i]
                continue
            if cs[i] != d[n]:
                d[n] = 'X'
                break
    return d

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

def rule2bnet(rule):
    """
    Converts a PyBoolNet prime rule into a BNet string.
    e.g., [{'A':1,'B':0},{'C':0}] returns 'A & !B | !C'

    There are two special cases:
    [] is identically 0
    [{}] is identically 1
    """
    if rule == []: return '0'
    elif rule == [{}]: return '1'
    else: return ' | '.join([implicant2bnet(imp) for imp in rule])

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

    if primes is None:
        return ""

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
    L_trim = [x for x in L if x != [{}]]
    u = []
    for x in L_trim:
        t = []
        for y in x:
            s = []
            for k,v in y.items():
                if v: s.append(k)
                else: s.append('!'+k)
            t.append(' & '.join(s))
        if len(t) > 0: u.append('( '+' | '.join(t)+' )')
    s=' & '.join(u)
    if simplify:
        s = PyBoolNet.BooleanLogic.minimize_espresso(s)
    if not silent:
        print(s)

    return s
