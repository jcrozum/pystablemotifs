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
    """A modified version of PyBoolNet's bnet2primes that does not do path-checking,
    as this can cause errors if the bnet rules are very long. Assumes BNET is a
    bnet string, not a file.

    Parameters
    ----------
    BNET : str
        BNET formatted rules to convert.
    remove_constants : bool
        Whether or not to remove and percolate constant input values (the default
        is False).

    Returns
    -------
    primes : PyBoolNet primes dictionary
        Update rules in PyBoolNet format.

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
    """Converts BooleanNet rules to BNet format.
    e.g., an input of
    "A*=B or C and not D"
    returns
    A,  B | C & !D

    Also replaces ~ with !

    Parameters
    ----------
    rules : str
        BooleanNet formatted rules.

    Returns
    -------
    str
        BNET formatted rules.

    """
    s = re.sub("\s*\*\s*=\s*",",\t",rules) # replace "=" with ",\t"
    s = re.sub("\s+not\s+"," !",s, flags=re.IGNORECASE) # not -> !
    s = re.sub("\(\s*not\s+","(!",s, flags=re.IGNORECASE) # not -> ! (with parens)
    s = re.sub("\s*~\s*"," !",s, flags=re.IGNORECASE) # ~ -> !
    s = re.sub("\s+and\s+"," & ",s, flags=re.IGNORECASE) # and -> &
    s = re.sub("\s+or\s+"," | ",s, flags=re.IGNORECASE) # or -> |
    s = re.sub("False","0",s, flags=re.IGNORECASE) # False -> 0 (ignore case)
    s = re.sub("True","1",s, flags=re.IGNORECASE) # True -> 1 (ignore case)
    return s

def bnetDNF2list(bnet):
    """Converts a BNet string expression to a list of prime implicant dictionaries.


    Parameters
    ----------
    bnet : str
        BNET formatted rules.

    Returns
    -------
    list of partial state dictionaries
        Variable states specified by each dictionary are to be thought of as "AND"-
        connected, and the dictionaries as "OR"-connected.

    """

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
    """Converts a BNet string expression (expr1) and its negation (expr0) to
    a PyBoolNet rule list. Note that this function does not test for consistency
    between expr0 and expr1.

    Parameters
    ----------
    expr0 : str
        Rule, in BNET format, for the "OFF" state of a variable.
    expr1 : str
        Rule, in BNET format, for the "ON" state of a variable.

    Returns
    -------
    PyBoolNet rule list
        The complementary expressions as they would appear in a PyBoolNet primes
        dictionary for a variable whose update rule is given by expr1.

    """
    return [bnetDNF2list(expr0),bnetDNF2list(expr1)]

def bnet2sympy(rule):
    """Converts a BNet string expression to a sympy string expression.

    Parameters
    ----------
    rule : str
        Boolean expression in BNET format.

    Returns
    -------
    str
        Expression in sympy format.

    """

    crule = re.sub("!","~",rule)
    crule = re.sub("[\b\(]1[\b\)]","(x | ~x)",crule)
    crule = re.sub("[\b\(]0[\b\)]","(x & ~x)",crule)
    crule = re.sub("True","(x | ~x)",crule)
    crule = re.sub("False","(x & ~x)",crule)

    return crule

def sympy2bnet(rule):
    """Converts a sympy string expression to a BNET string expression.

    Parameters
    ----------
    rule : str
        Boolean expression in sympy format.

    Returns
    -------
    str
        Expression in BNET format.

    """
    crule = re.sub("~","!",rule)
    crule = re.sub("True","1",crule)
    crule = re.sub("False","0",crule)
    return crule

def remove_comment_lines(stream, comment_char="#"):
    """Removes commented out lines from stream, e.g., those starting with '#'.

    Parameters
    ----------
    stream : iterable of str
        Lines from which comments should be excluded.
    comment_char : str
        Lines beginning with this character will be excluded.

    Returns
    -------
    list of str
        Lines that do not begin with comment_char.

    """
    lines = list(stream)
    lines = filter(lambda x: not x.startswith(comment_char), lines)
    rules = "".join(lines)
    return rules

def import_primes(fname, format='BooleanNet', remove_constants=False):
    """Import boolean rules from file and return PyBoolNet formatted primes list.

    Parameters
    ----------
    fname : str
        Path to (plaintext) file containing Boolean rules in format specified
        by the 'format' option.
    format : str
         Boolean rule format; options are 'BooleanNet' or 'BNet' (the default is
         'BooleanNet').
    remove_constants : bool
         If True, variables that are constant are removed and their influence is
         percolated. Otherwise, they remain and we consider initial conditions
         in opposition to their values (the default is False).

    Returns
    -------
    PyBoolNet primes dictionary
        Update rules in PyBoolNet format.

    """
    # TODO: add more formats
    rules = remove_comment_lines(open(fname))
    if format == 'BooleanNet':
        rules = booleannet2bnet(rules)
    elif format == 'BNet':
        rules = rules
    else:
        raise ValueError('Unrecognized format',format)

    primes = longbnet2primes(rules,remove_constants=remove_constants)

    return primes

def statelist2dict(names,statestrings):
    """Converts a collection of statestrings to a dictionary.

    Parameters
    ----------
    names : list of str
        An ordered list of variable names; (alphabetical order is PyBoolNet's
        default, e.g. sorted(primes)).
    c : iterable of str
        Each element should be a binary string, with each position corresponding
        to the variable name at the same position in names.

    Returns
    -------
    dictionary
        Dictionary summarizing c. If a node takes the same value in every state,
        the corresponding dictionary value matches its fixed value; otherwise,
        the dictionary value is 'X'.

    """
    d = {}
    for i,n in enumerate(names):
        for cs in statestrings:
            if n not in d:
                d[n] = cs[i]
                continue
            if cs[i] != d[n]:
                d[n] = 'X'
                break
    return d

def statestring2dict(statestring,names):
    """Converts a state string, which specifies a node in an STG, to the
    corresponding dictionary representation.

    Parameters
    ----------
    statestring : str
        A binary string, e.g., '01101'.
    names : list of str
        An ordered list of variable names; (alphabetical order is PyBoolNet's
        default, e.g. sorted(primes)).

    Returns
    -------
    dictionary
        The keys are the elements of names and the values are the corresponding
        value in statestring.

    """
    sd = {}
    for i,c in enumerate(statestring):
        sd[names[i]]=int(c)
    return sd

def statedict2str(statedict):
    """Converts a state dictionary to a statestring using alphabetical sorting.

    Parameters
    ----------
    statedict : partial state dictionary
        State to convert to a binary string representation.

    Returns
    -------
    str
        A binary string, with each position corresponding
        to the variable name at the same position in the alphabetized keys in
        statedict.

    """
    return ''.join([str(statedict[x]) for x in sorted(statedict)])

def implicant2bnet(partial_state):
    """Converts a partial state dictionary to a BNet string
    e.g., {'A':1,'B':0} returns 'A & !B'

    Parameters
    ----------
    partial_state : partial state dictionary
        Partial state to convert.

    Returns
    -------
    str
        BNET representation of the partial state.

    """
    return ' & '.join(["!"+k for k in partial_state if not partial_state[k]]+[k for k in partial_state if partial_state[k]])

def rule2bnet(rule):
    """Converts a PyBoolNet prime rule into a BNet string.
    e.g., [{'A':1,'B':0},{'C':0}] returns 'A & !B | !C'

    Parameters
    ----------
    rule : list of PyBoolNet partial states
        Update rule to convert.

    Returns
    -------
    str
        BNET representation of Boolean expression.

    """
    # There are two special cases:
    # [] is identically 0
    # [{}] is identically 1
    if rule == []: return '0'
    elif rule == [{}]: return '1'
    else: return ' | '.join([implicant2bnet(imp) for imp in rule])

def primes2bnet(primes):
    """A simpler version of PyBoolNet's FileExchange.primes2bnet function with
    fewer options and less organized output. Should handle prime rules with
    tautologies better than the PyBoolNet version though.

    Parameters
    ----------
    primes : PyBoolNet primes dictionary
        Update rules to convert.

    Returns
    -------
    str
        BNET representation of update rules.

    """
    lines = []
    width = max([len(x) for x in primes]) + 3

    for name in primes:
        if primes[name][0] == [] or primes[name][1] == [{}]:
            expression = '1'
        elif primes[name][1] == [] or primes[name][0] == [{}]:
            expression = '0'
        else:
            expression = ' | '.join(['&'.join([x if term[x]==1 else '!'+x for x in term]) for term in primes[name][1]  ])

        lines+= [(name+',').ljust(width)+expression]
    lines+=['']

    return "\n".join(lines)


def pretty_print_primes(primes):
    """Prints PyBoolNet a prime dictionary in a more readable format. Prints both
    state updates (1 and 0).

    Parameters
    ----------
    primes : PyBoolNet primes dictionary
        Update rules to print.

    """
    for k,v in primes.items():
        for p in v[0]: print(p,"=> !"+k)
        for p in v[1]: print(p,"=>  "+k)

def pretty_print_prime_rules(primes):
    """Prints PyBoolNet a prime dictionary as Boolean rules
    The output format is of the form:
    A* = B & C | !D, for example.

    Parameters
    ----------
    primes : PyBoolNet primes dictionary
        Update rules to print.

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
        if v[1]==[]:
            s = k + "* = 0"
        if v[1]==[{}]:
            s = k + "* = 1"
        print(s)

def pretty_print_rspace(L,simplify=True,silent=True):
    """Produces string representation of the Boolean rule describing the input
    rspace L (see RestrictSpace.rspace).

    Parameters
    ----------
    L : rspace list
        Restrict space list (see RestrictSpace.rspace for details).
    simplify : bool
        Whether to simplify the rule (the default is True).
    silent : bool
        Whether to supress output of the rule (the default is True).

    Returns
    -------
    str
        BNET expression that is true in and only in the rspace specified by L.

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
