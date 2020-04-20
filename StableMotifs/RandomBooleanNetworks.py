import pandas as pd
import random as rd
import numpy as np
from datetime import datetime
from uuid import uuid4
import os

class RandomBooleanNetworks:
    """
    Generator of random Boolean networks (RBN) and ensembles of RBN.
    The RandomBooleanNetworks class object is a Boolean model and stores information of how the Boolean model was generated
    The RandomBooleanNetworks class has functions that generate ensembles of RBN by generating multiple RandomBooleanNetworks objects

    Variables:

    Functions:
    """

    def __init__(self):
        self.node_names = []
        self.node_inputs_dictionary = {}
        self.node_rules_binary_dictionary = {}
        self.node_rules_decimal_dictionary = {}
        self.node_rules_string_dictionary = {}
        self.node_rules_string = ""
        self.random_Boolean_type = ""
        self.N=0
        self.random_Boolean_Network_parameters= []
        self.random_seed = None
        self.filename = None

    def Random_Boolean_Network(self,random_Boolean_type,N,rbn_parameters,seed=None,filename=None):
        self.random_Boolean_type = random_Boolean_type
        self.random_Boolean_Network_parameters= rbn_parameters
        self.N=N
        if seed is not None:
            rd.seed(seed)
        self.node_names=["n"+str(i) for i in range(N)]

        if(random_Boolean_type=="Kauffman NK"):
            K,p=self.random_Boolean_Network_parameters
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
                write_Boolean_network_decimal(self.node_rules_decimal_dictionary,self.filename)

    def Random_Boolean_Network_Rules(self):
        self.node_rules_string_dictionary=String_Rules_From_Binary(self.node_rules_binary_dictionary)
        self.node_rules_string=""
        for n,f in self.node_rules_string_dictionary.items():
             self.node_rules_string=self.node_rules_string+n+" *= "+f+"\n"
        return(self.node_rules_string)

def write_Boolean_network_decimal(node_rules_decimal_dictionary,filename):
    df_boolean_model=pd.DataFrame.from_dict(node_rules_decimal_dictionary,orient='index')
    df_boolean_model.to_csv(filename, header=False)

def read_Boolean_network_decimal(filename):
    df_dict=pd.read_csv(filename,header=None,index_col=0).to_dict('index')
    node_rules_decimal_dictionary={key:(element[1],[x[1:-1] for x in element[2].strip('][').split(', ')])  for key,element in df_dict.items()}
    node_rules_binary_dictionary=Binary_Rules_From_Decimal(node_rules_decimal_dictionary)
    node_rules_string_dictionary=String_Rules_From_Binary(node_rules_binary_dictionary)

    rules=""
    for n,f in node_rules_string_dictionary.items():
        rules=rules+n+" *= "+f+"\n"
    return(rules)

def Binary_Rules_From_Decimal(node_rules_decimal_dictionary):
    node_rules_binary_dictionary={node_name:
                                 (Binary_Rule_From_Decimal(node_rule_input_list[0],node_rule_input_list[1]),
                                            node_rule_input_list[1])
                                 for node_name,node_rule_input_list in node_rules_decimal_dictionary.items()}
    return(node_rules_binary_dictionary)

def Binary_Rule_From_Decimal(node_rule_decimal,node_input_list):
    binary_rule=[int(binary) for binary in np.binary_repr(node_rule_decimal, width=np.power(2,len(node_input_list)))]
    return(binary_rule)

def String_Rule_From_Binary(node_rule_binary,node_input_list):
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
    node_rules_string_dictionary={node_name:String_Rule_From_Binary(node_rule_input_list[0],node_rule_input_list[1])
                                  for node_name,node_rule_input_list in node_rules_binary_dictionary.items()}
    return(node_rules_string_dictionary)

def Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=1000,write_Boolean_network=True):
    """
    Obtain the string version of the motif_history of a reduced network
    Given the motif_history of a reduction.from motif_reduction_list, obtain a text version of this motif_history

    Inputs:
    N - int, Number of nodes of RBN
    K - int, Number of inputs of each node in the RBN
    p - double, probability that each entry in the RBN is equal to 1
    N_ensemble - int, Number of elements in the RBN ensemble
    seed - int, random seed for generating the RBN ensemble
    write_Boolean_network - boolean, if True, will write each network in the ensemble as a CSV file in a new directory

    Outputs:
    RBN_ensemble_rules - List of strings, each string are the Boolea rules of an ensemble in booleannet format.
        Each element in RBN_ensemble_rules can be used as an input for the Format.booleannet2bnet function

    """
    rd.seed(seed)
    random_Boolean_type = "Kauffman NK"
    RBN_ensemble_rules=[]
    if(write_Boolean_network):
        directory = "RBN"+"_N-"+str(N)+"_K-"+str(K)+"_p-"+str(p)+"_"+ datetime.now().strftime('%Y%m-%d%H-%M%S-') + str(uuid4())
        try:
            os.makedirs(directory)
        except FileExistsError:
            # directory already exists
            print ('Error: Creating directory ' +  directory)
            pass
    for n in range(N_ensemble):
        if(write_Boolean_network):
            filename="RBN"+"_N-"+str(N)+"_K-"+str(K)+"_p-"+str(p)+"_Nensemble-"+str(n)+".csv"
            rbn=RandomBooleanNetworks()
            rbn.Random_Boolean_Network(random_Boolean_type,N,rbn_parameters=[K,p],filename=os.path.join(directory, filename))
            rules=rbn.Random_Boolean_Network_Rules()
            RBN_ensemble_rules.append(rules)
        else:
            rbn=RandomBooleanNetworks()
            rbn=RandomBooleanNetworks.Random_Boolean_Network(random_Boolean_type,N,rbn_parameters=[K,p])
            rules=rbn.Random_Boolean_Network_Rules()
            RBN_ensemble_rules.append(rules)
    return(RBN_ensemble_rules)

def get_criticality_K_Kauffman(p):
    """
    The Kauffman RBN is at criticality when S = 2p(1-p)K = 1. Given p, K = 2/(p(1-p))

    Inputs:
    p - double, probability that each entry in the RBN is equal to 1
    Outputs:
    K_criticality - list, Number of inputs (int) per node from which the Kauffman RBN is at criticality
    """

    K_criticality=[2.0/(p*(1-p))]
    return(K_criticality)

def get_criticality_p_Kauffman(K):
    """
    The Kauffman RBN is at criticality when S = 2p(1-p)K = 1. Given p, K = 2/(p(1-p))

    Inputs:
    K - int, Number of inputs per node for Kauffman RBN
    Outputs:
    p_criticality - list, probability that each entry in the RBN is equal to 1 for which the Kauffman RBN is at criticality
    """
    p_criticality=[]
    if(K<2):
        p_criticality=[]
    elif(K==2):
        p_criticality=[0.5]
    elif(K>2):
        p_criticality=[0.5*(1+np.sqrt(1.0-2.0/K)),0.5*(1-np.sqrt(1.0-2.0/K))]
    return(p_criticality)
