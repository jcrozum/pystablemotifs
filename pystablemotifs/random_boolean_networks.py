import pandas as pd
import random as rd
import numpy as np
from datetime import datetime
from uuid import uuid4
import os

class RandomBooleanNetworks:
    """Generator of random Boolean networks (RBN) and ensembles of RBN.
    The RandomBooleanNetworks class object is a Boolean model and stores
    information of how the Boolean model was generated. It has functions that
    generate ensembles of RBN by generating multiple RandomBooleanNetworks objects.

    Attributes
    ----------
    node_names : list of str
        List of variable names.
    node_inputs_dictionary : dictionary
        Each value is a (fixed order) list of the names of the nodes whose values
        are inputs into the key variable's update function.
    node_rules_binary_dictionary : dictionary
        Each value is a list of outputs for the key variable's update function,
        stored as a list in ascending order of the numerical representation of
        the input row.
    node_rules_decimal_dictionary : dictionary
        Decimal conversion of node_rules_binary_dictionary.
    node_rules_string_dictionary : dictionary
        BooleanNet (str) conversion of node_rules_binary_dictionary values.
    node_rules_string : str
        BooleanNet representation of update rules.
    random_boolean_type : str
        Descrpition of generative process. Currently only "Kauffman NK" is
        implemented.
    N : int
        Number of nodes in the Boolean network.
    random_boolean_Network_parameters : list
        For Kauffman NK generation -
            [K,p], where K is the in-degree and p is the bias. K is a positive
            integer less than or equal to N, and p is a float between 0 and 1
            (inclusive).
    random_seed : int
        Seed for random functions.
    filename : str
        Path to file where network data are stored. If None, no files are written.

    """

    def __init__(self):
        self.node_names = []
        self.node_inputs_dictionary = {}
        self.node_rules_binary_dictionary = {}
        self.node_rules_decimal_dictionary = {}
        self.node_rules_string_dictionary = {}
        self.node_rules_string = ""
        self.random_boolean_type = ""
        self.N=0
        self.random_boolean_Network_parameters= []
        self.random_seed = None
        self.filename = None

    def random_boolean_network(self,random_boolean_type,N,rbn_parameters,seed=None,filename=None):
        """Construct network using specified generative process.

        Parameters
        ----------
        random_boolean_type : str
            Descrpition of generative process. Currently only "Kauffman NK" is
            implemented.
        N : int
            Number of nodes in the Boolean network.
        random_boolean_Network_parameters : list
            For Kauffman NK generation -
                [K,p], where K is the in-degree and p is the bias. K is a positive
                integer less than or equal to N, and p is a float between 0 and 1
                (inclusive).
        random_seed : int
            Seed for random functions.
        filename : str
            Path to file where network data are stored. If None, no files are written.

        """
        self.random_boolean_type = random_boolean_type
        self.random_boolean_Network_parameters= rbn_parameters
        self.N=N
        if seed is not None:
            rd.seed(seed)
        self.node_names=["n"+str(i) for i in range(N)]

        if(random_boolean_type=="Kauffman NK"):
            K,p=self.random_boolean_Network_parameters
            self.node_inputs_dictionary={node_name:rd.sample(self.node_names,K) for node_name in self.node_names}
            self.node_rules_binary_dictionary={node_name:
                                         (rd.choices([1,0],weights=[p,1-p],k=np.power(2,len(node_input_list))),
                                          node_input_list)
                                         for node_name,node_input_list in self.node_inputs_dictionary.items()}
            if filename is not None:
                self.filename = filename
                self.node_rules_decimal_dictionary={node_name:
                                      (int("".join([str(i) for i in node_rule_input_list[0]]),2),
                                       node_rule_input_list[1])
                                      for node_name,node_rule_input_list in self.node_rules_binary_dictionary.items()}
                write_boolean_network_decimal(self.node_rules_decimal_dictionary,self.filename)

    def random_boolean_network_Rules(self):
        """Generate various conversions of the node_rules_binary_dictionary
        attribute.

        """
        self.node_rules_string_dictionary=String_Rules_From_Binary(self.node_rules_binary_dictionary)
        self.node_rules_string=""
        for n,f in self.node_rules_string_dictionary.items():
             self.node_rules_string=self.node_rules_string+n+" *= "+f+"\n"
        return(self.node_rules_string)

def write_boolean_network_decimal(node_rules_decimal_dictionary,filename):
    """Write the decimal format of the Boolean rules to file.

    Parameters
    ----------
    node_rules_decimal_dictionary : dictionary
        Update rule truth table in decimal format.
    filename : str
        Path to file for csv output of the truth table.

    """
    df_boolean_model=pd.DataFrame.from_dict(node_rules_decimal_dictionary,orient='index')
    df_boolean_model.to_csv(filename, header=False)

def read_boolean_network_decimal(filename):
    """Imports rules from csv in decimal format.

    Parameters
    ----------
    filename : str
        Path to csv from which to import decimal-formatted rules.

    Returns
    -------
    str
        Rules in BooleanNet format.

    """
    df_dict=pd.read_csv(filename,header=None,index_col=0).to_dict('index')
    node_rules_decimal_dictionary={key:(element[1],[x[1:-1] for x in element[2].strip('][').split(', ')])  for key,element in df_dict.items()}
    node_rules_binary_dictionary=Binary_Rules_From_Decimal(node_rules_decimal_dictionary)
    node_rules_string_dictionary=String_Rules_From_Binary(node_rules_binary_dictionary)

    rules=""
    for n,f in node_rules_string_dictionary.items():
        rules=rules+n+" *= "+f+"\n"
    return(rules)

def Binary_Rules_From_Decimal(node_rules_decimal_dictionary):
    """Construct Binary format rules from decimal format rules.

    Parameters
    ----------
    node_rules_decimal_dictionary : dictionary
        Rules in decimal format to convert.

    Returns
    -------
    dictionary
        Binary rules dictionary.

    """
    node_rules_binary_dictionary={node_name:
                                 (Binary_Rule_From_Decimal(node_rule_input_list[0],node_rule_input_list[1]),
                                            node_rule_input_list[1])
                                 for node_name,node_rule_input_list in node_rules_decimal_dictionary.items()}
    return(node_rules_binary_dictionary)

def Binary_Rule_From_Decimal(node_rule_decimal,node_input_list):
    """Convert single decimal rule to its binary form.

    Parameters
    ----------
    node_rule_decimal : int
        Decimal form of a truth table's output column.
    node_input_list : list of str
        Variable names that correspond to each column of the truth table.

    Returns
    -------
    list of int
        Binary rule list corresponding to an output column of a truth table.

    """
    binary_rule=[int(binary) for binary in np.binary_repr(node_rule_decimal, width=np.power(2,len(node_input_list)))]
    return(binary_rule)

def String_Rule_From_Binary(node_rule_binary,node_input_list):
    """Convert binary rule to BooleanNet format.

    Parameters
    ----------
    node_rule_binary : list of int
        Binary rule list corresponding to an output column of a truth table.
    node_input_list : list of str
        Variable names that correspond to each column of the truth table.

    Returns
    -------
    str
        BooleanNet representation of rule.

    """
    notConstant=False
    bit_previous=-1
    for i,bit in enumerate(node_rule_binary):
        if(i==0):
            bit_previous=bit
            continue
        if(bit!=bit_previous):
            notConstant=True
            break;

    if(not notConstant):
        rule_string=str(bit_previous)
    else:
        implicants_rule=[[int(binary) for binary in np.binary_repr(i,width=len(node_input_list))] for i,bit in enumerate(node_rule_binary) if bit==1]
        rule_string=" or ".join(
            ["("+" and ".join(
            [node_input_list[k] if(state==1) else "not "+node_input_list[k] for k,state in enumerate(implicant)]
        )+")"
            for implicant in implicants_rule]
        )
    return(rule_string)

def String_Rules_From_Binary(node_rules_binary_dictionary):
    """Convert from binary dictionary rule format to BooleanNet format.

    Parameters
    ----------
    node_rules_binary_dictionary : dictionary
        Binary dictionary representation of rules.

    Returns
    -------
    str
        BooleanNet representation of rules.

    """
    node_rules_string_dictionary={node_name:String_Rule_From_Binary(node_rule_input_list[0],node_rule_input_list[1])
                                  for node_name,node_rule_input_list in node_rules_binary_dictionary.items()}
    return(node_rules_string_dictionary)

def random_boolean_network_ensemble_kauffman(N,K,p,N_ensemble,seed=1000,write_boolean_network=False):
    """Generate a sample from the Kauffman NK RBN ensemble.

    Parameters
    ----------
    N : int
        Number of nodes of RBN.
    K : int
        Number of inputs of each node in the RBN.
    p : float
        Probability that each entry in each truth table output column is equal to 1.
    N_ensemble : int
        Number of networks to generate.
    seed : int
        Random seed for generating the RBN ensemble (the default is 1000).
    write_boolean_network : bool
        Whether to write each network in the ensemble as a CSV file in a new
        directory (the default is False).

    Returns
    -------
    RBN_ensemble_rules : list of str
        Each string are the Boolea rules of an ensemble in booleannet format.
        Each element in RBN_ensemble_rules can be used as an input for the
        format.booleannet2bnet function.

    """
    rd.seed(seed)
    random_boolean_type = "Kauffman NK"
    RBN_ensemble_rules=[]
    if(write_boolean_network):
        directory = "RBN"+"_N-"+str(N)+"_K-"+str(K)+"_p-"+str(p)+"_"+ datetime.now().strftime('%Y%m-%d%H-%M%S-') + str(uuid4())
        try:
            os.makedirs(directory)
        except FileExistsError:
            # directory already exists
            print ('Error: Creating directory ' +  directory)
            pass
    for n in range(N_ensemble):
        if(write_boolean_network):
            filename="RBN"+"_N-"+str(N)+"_K-"+str(K)+"_p-"+str(p)+"_Nensemble-"+str(n)+".csv"
            rbn=RandomBooleanNetworks()
            rbn.random_boolean_network(random_boolean_type,N,rbn_parameters=[K,p],filename=os.path.join(directory, filename))
            rules=rbn.random_boolean_network_Rules()
            RBN_ensemble_rules.append(rules)
        else:
            rbn=RandomBooleanNetworks()
            rbn.random_boolean_network(random_boolean_type,N,rbn_parameters=[K,p])
            rules=rbn.random_boolean_network_Rules()
            RBN_ensemble_rules.append(rules)
    return(RBN_ensemble_rules)

def get_criticality_K_Kauffman(p):
    """The Kauffman RBN is at criticality when K = 2/(p(1-p)).

    Parameters
    ----------
    p : float
        Probability that each entry in each truth table output column is equal to 1.

    Returns
    -------
    K_criticality : int
        Number of inputs of each node in the RBN.

    """


    K_criticality=[2.0/(p*(1-p))]
    return(K_criticality)

def get_criticality_p_Kauffman(K):
    """The Kauffman RBN is at criticality when K = 2/(p(1-p)).

    Parameters
    ----------
    K : int
        Number of inputs of each node in the RBN.


    Returns
    -------
    p_criticality : float
        Probability that each entry in each truth table output column is equal to 1.

    """
    p_criticality=[]
    if(K<2):
        p_criticality=[]
    elif(K==2):
        p_criticality=[0.5]
    elif(K>2):
        p_criticality=[0.5*(1+np.sqrt(1.0-2.0/K)),0.5*(1-np.sqrt(1.0-2.0/K))]
    return(p_criticality)
