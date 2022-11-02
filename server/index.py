from flask import Flask, jsonify, send_file, Response, send_file
from rdkit import Chem
from rdkit.Chem import AllChem
import pybel
from io import *
import os
import zipfile
from django.http import HttpResponse

import random
import sys
import pandas as pd
import numpy as np
import progressbar
from IPython.display import display
import matplotlib as mpl
from matplotlib import pyplot as plt
import pylab as pl
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
import plotly.graph_objects as go
import dash
import dash_table
import dash_html_components as html
import chart_studio.plotly as py
import seaborn as sns
import plotly.express as px
import cufflinks as cf
import re
app =  Flask(__name__)

class PXA:
    def __init__(self, smile='',group='',formula='', numCl=0,numF=0,numBr=0,numC=0,numH=0,numHalo=0):
        self._smile=smile
        self._formula=formula
        self._group=group
        self._numCl=numCl
        self._numBr=numBr
        self._numF=numF
        self._numC=numC
        self._numH=numH
        self._numHalo=numHalo
        self.positions=[]

    def get_NumCl(self):
        return self._numCl
    def get_NumBr(self):
        return self._numBr
    def get_NumF(self):
        return self._numF
    def get_NumH(self):
        return self._numH
    def get_NumC(self):
        return self._numC
    def get_NumHalo(self):
        return self._numHalo
    def get_Formula(self):
        return self._formula
    
    def split(self,word):
        return [char for char in word]

    #~~~~~~ This section contains functions needed to generate the SMILES structures of all linear isomers of a compound ~~~~~~~~~

    #Generates all permutations of a string without repetitions
    def permutation(self,nCl=0,nF=0,nBr=0):
        elementsList=[]
        numHalo=nCl+nF+nBr
        for i in range(0,self._numH): elementsList.append('H')
        if nCl > 0:
            for x in range(0,nCl): elementsList.append('C')
        if nF > 0:
            for x in range(0,nF): elementsList.append('F')
        if nBr > 0:
            for x in range(0,nBr): elementsList.append('B')
        perms = [[]]
        for n in elementsList:
            new_perm = []
            for perm in perms:
                for i in range(len(perm) + 1):
                    new_perm.append(perm[:i] + [n] + perm[i:])
                    # handle duplication
                    if i < len(perm) and perm[i] == n: break
            perms = new_perm
        return perms
    #Searches for positions in a string which a specific character appears
    def duplicates(self,seq,item):
        start_at = -1
        locs = []
        while True:
            try:
                loc = seq.index(item,start_at+1)
            except ValueError:
                break
            else:
                locs.append(loc)
                start_at = loc
        return locs
    # Replaces a hydrogen from the C-H SMILES structure with a halogen, guided by a permutation
    def replaceSmileFromPermutation(self,positions_1=[],positions_2=[],positions_3=[],maxH=0, halo_1="",halo_2="", halo_3="",smile="",cont=0):
        cont=0
        smile=self.split(smile)
        for index in range(0,len(smile)):
                if smile[index] != 'H':
                    continue
                else:
                    cont+=1

                    if positions_1:
                        if (cont-1) in positions_1:
                            smile[index]=halo_1
                    if positions_2:
                        if (cont-1) in positions_2:
                            smile[index]=halo_2
                    if positions_3:
                        if (cont-1) in positions_3:
                            smile[index]=halo_3


        smile="".join(smile)
        return smile
    # Generates all linear isomers of a compound which group corresponds to a single halogen
    def generateGroupOfOne(self,group=""):
        numH=self._numH
        allPossibleSmiles=[]
        elementToLook=''
        halo=""

        if group=="Cl":
            perms=self.permutation(nCl=self._numCl)
            elementToLook='C'
            halo="Cl"
        if group=="F":
            perms=self.permutation(nF=self._numF)
            elementToLook='F'
            halo="F"
        if group=="Br":
            perms=self.permutation(nBr=self._numBr)
            elementToLook='B'
            halo="Br"

        for p in perms:
            smile=self._smile
            elementInPerms=''.join(p)
            positions=self.duplicates(elementInPerms,elementToLook)

            smile=self.replaceSmileFromPermutation(positions_1=positions,maxH=numH,halo_1=halo,smile=smile)
            allPossibleSmiles.append(smile)

        return allPossibleSmiles

    # Generates all linear isomers of a compound which group corresponds to two halogens
    def generateGroupOfTwo(self,group=''):
        numH=self._numH
        allPossibleSmiles=[]
        elementToLook_1=''
        elementToLook_2=''
        halo_1=""
        halo_2=""

        if group =="BrCl":
            perms=self.permutation(nCl=self._numCl,nBr=self._numF)
            elementToLook_1='C'
            elementToLook_2='B'
            halo_1="Cl" 
            halo_2="Br"
        if group =="ClF":
            perms=self.permutation(nCl=self._numCl,nF=self._numF)
            elementToLook_1='C'
            elementToLook_2='F'
            halo_1="Cl" 
            halo_2="F"
        if group == "BrF":
            perms=self.permutation(nBr=self._numBr,nF=self._numF)
            elementToLook_1='B'
            elementToLook_2='F'
            halo_1="Br" 
            halo_2="F"

        for p in perms:
            smile=self._smile
            elementInPerms=''.join(p)
            positions_1=self.duplicates(elementInPerms,elementToLook_1)
            positions_2=self.duplicates(elementInPerms,elementToLook_2)

        smile=self.replaceSmileFromPermutation(positions_1=positions_1,positions_2=positions_2,maxH=numH,halo_1=halo_1,halo_2=halo_2,smile=smile)
        allPossibleSmiles.append(smile)

        return allPossibleSmiles

    # Generates all linear isomers of a compound which group corresponds to three halogens     
    def generateGroupOfThree(self,group=""):
        numH=self._numH
        allPossibleSmiles=[]
        elementToLook_1=''
        elementToLook_2=''
        elementToLook_3=''
        halo_1=""
        halo_2=""
        halo_3=""

        smile=self._smile

        perms=self.permutation(nCl=self._numCl,nBr=self._numBr,nF=self._numF)
        elementToLook_1='C'
        elementToLook_2='B'
        elementToLook_3='F'
        halo_1="Cl" 
        halo_2="Br"
        halo_3="F"

        for p in perms:
            smile=self._smile
            elementInPerms=''.join(p)
            positions_1=self.duplicates(elementInPerms,elementToLook_1)
            positions_2=self.duplicates(elementInPerms,elementToLook_2)
            positions_3=self.duplicates(elementInPerms,elementToLook_3)

        smile=self.replaceSmileFromPermutation(positions_1=positions_1,positions_2=positions_2,positions_3=positions_3,maxH=numH,halo_1=halo_1,halo_2=halo_2,halo_3=halo_3,smile=smile)
        allPossibleSmiles.append(smile)

        return allPossibleSmiles
    def generateAllPossible(self):
        allPossibleSmiles=[]
        self._smile=self.generate_CH_Smile()

        if self._group =="Cl" :
            allPossibleSmiles=self.generateGroupOfOne(self._group)
        elif self._group =="Br":
            allPossibleSmiles=self.generateGroupOfOne(self._group)
        elif self._group =="F":
            allPossibleSmiles=self.generateGroupOfOne(self._group)
        elif self._group =="BrCl":
            allPossibleSmiles=self.generateGroupOfTwo(self._group)
        elif self._group =="ClF":
            allPossibleSmiles=self.generateGroupOfTwo(self._group)
        elif self._group == "BrF":
            allPossibleSmiles=self.generateGroupOfTwo(self._group)
        elif self._group == "BrClF":
            allPossibleSmiles=self.generateGroupOfThree(self._group)
        else:
            print("Group not recognized")
            return

        return allPossibleSmiles 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Replaces multiple strings at a time
    def replaceMultiple(self, mainString, toBeReplaced, newString):
        for elem in toBeReplaced :
            if elem in mainString :
                mainString = mainString.replace(elem, newString)

        return  mainString

#Generates a SMILES structure in a format containing just carbon and hydrogen by reading a chemical formula (CnH2*n+2).
#  Example: For the chemical formula: C10H22, the SMILES: 'C(C(C(C(C(C(C(C(C(C(H)(H)H)(H)H)(H)H)(H)H)(H)H)(H)H)(H)H)(H)H)(H)H)(H)(H)H' is generated.
# A second option is added, where the SMILES can be generated by inputing the number of carbons and hydrogens instead of the chemical formula.
    def generate_CH_Smile(self):
        if self._numH > 0 and self._numC > 0:
            total_H=self._numH+self._numCl+self._numBr+self._numF
            mult_H=int((total_H - 2)/2)
            str_carbon= "C("*self._numC
            str_H="H)(H)"*mult_H
            newSmile=str_carbon+str_H+"(H)H"

        else:
            arr_formula=re.sub( r"([A-Z])", r" \1", self._formula).split()
            numC=int(arr_formula[0].replace('C',''))
            numH=int(arr_formula[1].replace('H',''))

            total_H=numH
            mult_H=int((total_H - 2)/2)
            str_carbon= "C("*numC
            str_H="H)(H)"*mult_H
            newSmile=str_carbon+str_H+"(H)H"

        return newSmile

    #Generates a halogenated SMILES structure by replacing a hydrogen ('H') with a halogen ('Cl','Br', or 'F')
    #  When a carboxylic acid is added to the original structure, it replaces the string 'X' with halogen instead of 'H', this way the hydrogen corresponding to the acid is not edited
    def generate(self,smile,num,halogen,toBeReplaced ='H'):
        smileStr=list(smile)

        if num > len(smileStr):
            smileStr=[halogen for i in smileStr]
        else:
            i=0
            while i in range(num):
                position=random.randrange(0,len(smileStr))
                if smileStr[position] == toBeReplaced:
                    while smileStr[position]==halogen:
                           position=random.randrange(0,len(smileStr))
                    smileStr[position]=halogen
                    self.positions.append(position)
                    i=i+1

        smileStr="".join(smileStr)
        return smileStr

    #Creates a halogenated structure containing just one halogen, i.e. Groups= "Cl","Br",or "F"
    # Two options: 1) Manipulates a SMILES structure containing a carboxylic acid (Change  variable 'carboxylic=True')
    #              2) Creates a new congener (Default, i.e 'carboxylic=False' )
    def add_Group_Of_One(self,smile,x,numC,numH,group,halo,carboxylic=False):
        non=''
        if carboxylic:
          # Returns new information of new structure with added Carboxylic acid
            newSmile=self.generate(smile,(x-1) if numH == 0 else x,halo,"X")
            num_halo= (x-1) if numH == 0 else x
            numH= (numH+1) if numH == 0 else numH
            numHalo=x
            formula=f"C{self._numC}H{self._numH}{halo}{x if x > 1 else non}O2"  
            return newSmile, num_halo, numH, numHalo, formula
        else:
          # Creates an individual row in dataframe for each new congener with columns ["SMILE","Group","Formula","Cl","F","Br","C","H","Halo"]
            individualList=[]
            individualList.extend([self.generate(smile,x,halo),group,f"C{numC}{'H'if numH-x > 0 else ''}{numH-x if numH-x > 1 else non }{halo}{x if x > 1 else non}", x if group =='Cl' else 0, x if group=='F' else 0, x if group=='Br' else 0,numC,numH-x,x]) 
            return individualList

    #Creates a halogenated structure containing two halogens, i.e. Groups= "ClF","BrCl",or "BrF"
    # Two options: 1) Manipulates a SMILES structure containing a carboxylic acid (Change  variable 'carboxylic=True')
    #              2) Creates a new congener (Default, i.e 'carboxylic=False' )
    def add_Group_Of_Two(self,smile,h1,h2,numC,numH,group,halo1,halo2,carboxylic=False):
        non=''
        if carboxylic:
          # Returns new information of new structure with added Carboxylic acid
            newSmile=self.generate(smile,(h1-1) if numH == 0 else h1,halo1,"X")
            newSmile=self.generate(newSmile,h2,halo2,"X")
            num_halo1=(h1-1) if self._numH == 0 else h1
            numH= (numH+1) if numH == 0 else numH
            numHalo=h1+h2
            formula=f"C{numC}H{numH}{halo1}{h1 if h1 > 1 else non}{halo2}{h2 if h2 > 1 else non}O2"
            return newSmile, num_halo1, numH, numHalo, formula
        else:
            # Creates and individual row in dataframe for each congener with columns ["SMILE","Group","Formula","Cl","F","Br","C","H","Halo"]
            individualList=[]
            smile_new=""
            smile_new=self.generate(smile,h1,halo1)
            smile_new=self.generate(smile_new,h2,halo2)
            individualList.extend([smile_new,group,f"C{numC}{'H'if numH-(h1+h2) > 0 else ''}{numH-(h1+h2) if numH-(h1+h2) > 1 else non}{halo1}{h1 if h1 > 1 else non}{halo2}{h2 if h2 > 1 else non}", h1 if halo1=='Cl' else 0, h1 if halo1=='F' else 0 or h2 if halo2=='F' else 0, h2 if halo2=='Br' else 0,numC,numH-(h1+h2),h1+h2])
            return individualList

    # Generates all halogenated congeners' SMILES for a corresponding chemical formula read from a dataframe
    def generate_Congeners(self,smile_CH,numH,numC):
        all_Possible=[]

        #Generates all halogenated SMILES containing one halogen, i.e. Groups= "Cl","Br",and "F"
        for x in range(1,numH+1):
            individualList_Cl=self.add_Group_Of_One(smile_CH,x,numC,numH,'Cl','Cl')
            individualList_F=self.add_Group_Of_One(smile_CH,x,numC,numH,'F','F')
            individualList_Br=self.add_Group_Of_One(smile_CH,x,numC,numH,'Br','Br')

            all_Possible.extend([individualList_Cl,individualList_F,individualList_Br])

        #Generates all halogenated SMILES containing two halogens, i.e. Groups= "ClF","BrCl",and "BrF"
        for h1 in range(1,numH+2):

            for h2 in range(1,(numH+1)-h1):
                individualList_ClF=self.add_Group_Of_Two(smile_CH,h1,h2,numC,numH,'ClF',"Cl","F")
                individualList_BrCl=self.add_Group_Of_Two(smile_CH,h1,h2,numC,numH,'BrCl',"Cl","Br")
                individualList_BrF=self.add_Group_Of_Two(smile_CH,h1,h2,numC,numH,'BrF',"F","Br")
                all_Possible.extend([individualList_ClF,individualList_BrCl,individualList_BrF])

        smile=""
        individualList=[]
        non=""

        #Generates all halogenated SMILES containing three halogen, i.e. Group= "BrClF"
        for h1 in range(1,numH+2):

            for h2 in range(1,(numH+1)-h1):

                for h3 in range(1,(numH+1)-(h1+h2)):
                    individualList=[]
                    smile=self.generate(smile_CH,h1,"Cl")
                    smile=self.generate(smile,h2,"F")
                    smile=self.generate(smile,h3,"Br")
                    individualList.extend([smile,"BrClF",f"C{numC}{'H'if numH-(h1+h2+h3) > 0 else ''}{numH-(h1+h2+h3) if numH-(h1+h2+h3) > 1 else non}Br{h3 if h3 > 1 else non}Cl{h1 if h1 > 1 else non}F{h2 if h2 > 1 else non}",h1,h2,h3,numC,numH-(h1+h2+h3),h1+h2+h3])
                    all_Possible.append(individualList)

        return all_Possible

    # Generates a new random SMILES structure for existent congener, i.e. generates a new SMILES with the same amount of carbons, hydrogens, chlorines, flourines, and bromines.
    def generateRandomSMILE(self):
        smile_CH=self.generate_CH_Smile()

        if self._numCl > 0:
            smile_CH=self.generate(smile_CH,self._numCl,"Cl")
        if self._numBr > 0:
            smile_CH=self.generate(smile_CH,self._numBr,"Br")
        if self._numF > 0:
            smile_CH=self.generate(smile_CH,self._numF,"F")
        return smile_CH

    # Generates a new SMILES structure with an added carboxylic acid for each congener, updating all the information of each compound, i.e. ["SMILE","Group","Formula","Cl","F","Br","C","H","Halo"]    
    def add_CarboxylicAcid(self):
        newSmile=None
        str_carbon= "C("*(self._numC+1)
        str_carboxylic='=O)OH)'
        str_X="(X)X)"*(self._numC-1)
        newSmile=str_carbon+str_carboxylic+str_X+"(X)(X)X"
        non=''
        self._numC+=1

        if self._group == 'Cl':
            newSmile, self._numCl, self._numH, self._numHalo, self._formula =self.add_Group_Of_One(newSmile,self._numCl,self._numC,self._numH,self._group,"Cl",carboxylic=True)
        if self._group == 'F':
            newSmile, self._numF, self._numH, self._numHalo, self._formula =self.add_Group_Of_One(newSmile,self._numF,self._numC,self._numH,self._group,"F",carboxylic=True)

        if self._group == 'Br':
            newSmile, self._numBr, self._numH, self._numHalo, self._formula =self.add_Group_Of_One(newSmile,self._numBr,self._numC,self._numH,self._group,"Br",carboxylic=True)

        if self._group == 'BrCl':
            newSmile, self._numBr, self._numH, self._numHalo, self._formula=self.add_Group_Of_Two(newSmile,self._numBr,self._numCl,self._numC,self._numH,self._group,"Br","Cl",carboxylic=True)

        if self._group == 'ClF':
            newSmile, self._numCl, self._numH, self._numHalo, self._formula=self.add_Group_Of_Two(newSmile,self._numCl,self._numF,self._numC,self._numH,self._group,"Cl","F",carboxylic=True)

        if self._group == 'BrF':
            newSmile, self._numBr, self._numH, self._numHalo, self._formula=self.add_Group_Of_Two(newSmile,self._numBr,self._numF,self._numC,self._numH,self._group,"Br","F",carboxylic=True)

        if self._group == 'BrClF':
            newSmile=self.generate(newSmile,(self._numBr-1) if self._numH == 0 else self._numBr,"Br","X")
            newSmile=self.generate(newSmile,self._numF,"F","X")
            newSmile=self.generate(newSmile,self._numCl,"Cl","X")
            self._numBr=self._numBr-1 if self._numH == 0 else self._numBr
            self._numH=self._numH+1 if self._numH == 0 else self._numH
            self._numHalo=self._numBr +self._numCl+self._numF
            self._formula=f"C{self._numC}H{self._numH}{'Br' if self._numBr > 0 else non }{self._numBr if self._numBr > 1 else non}{'Cl' if self._numCl > 0 else non }{self._numCl if self._numCl > 1 else non}{'F' if self._numF > 0 else non }{self._numF if self._numF > 1 else non}O2"
        if self._numH > 0:
            newSmile=self.generate(newSmile,(self._numH-1),"H","X")

        return newSmile



#~~~~~~~~~~~~~~~~~~~Functions used in main function for generating SMILES structures ~~~~~~~~~~~~~~~~~~~~#

# Generates all linear isomers of a single compound 
# Must enter the compound's chemical formula as a parameter
def getIsomers(compound_formula):
    num_C,num_H,num_Cl,num_F,num_Br,compound_group=getCompound_Information(compound_formula)
    pxa=PXA(group=compound_group,numC=num_C,numH=num_H,numCl=num_Cl,numF=num_F,numBr=num_Br)

    all_isomers=pxa.generateAllPossible()
    isomers_df = pd.DataFrame (all_isomers,columns=['SMILE'])

    return isomers_df

def getGroup(numCl,numF,numBr):
        
    if numCl > 0 and numF == 0 and numBr==0:
        group="Cl"

    if numF > 0 and numCl==0 and numBr==0:
        group="F"

    if numBr > 0 and numCl==0 and numF==0:
        group="Br"

    if numBr > 0 and numCl>0 and numF==0:
        group="ClBr"

    if numBr == 0 and numCl>0 and numF>0:
        group="ClF"

    if numBr > 0 and numCl==0 and numF>0:
        group="FBr"

    if numBr > 0 and numCl>0 and numF>0:
        group="BrClF"
    return group

def getCompound_Information(compound_formula):
    formula_Arr= re.sub( r"([A-Z])", r" \1", compound_formula).split()
    numCl= numF=numBr=numC=numH=0
    arr_Cl=[i for i in formula_Arr if "Cl" in i]
    arr_F=[i for i in formula_Arr if "F" in i]
    arr_Br=[i for i in formula_Arr if "Br" in i]
    arr_H=[i for i in formula_Arr if "H" in i]

    if arr_Cl:
        numCl=int(arr_Cl[0].replace('Cl',''))
        formula_Arr.remove(arr_Cl[0])
    if arr_F:
        numF=int(arr_F[0].replace('F',''))
        formula_Arr.remove(arr_F[0])
    if arr_Br:
        numBr=int(arr_Br[0].replace('Br',''))
        formula_Arr.remove(arr_Br[0])
    if arr_H:
        numH=int(arr_H[0].replace('H',''))
        formula_Arr.remove(arr_H[0])

    arr_C=[i for i in formula_Arr if "C" in i]
    if arr_C:
        numC=int(arr_C[0].replace('C',''))
        formula_Arr.remove(arr_C[0])

    group=getGroup(numCl,numF,numBr)

    return numC,numH,numCl,numF,numBr,group


# Generates all congeners corresponding to the 15 alkanes of carbon chains 10-to-25
# The dataframe(df) must contain the chemical formula ['Formula'] of the saturated hydrocarbons (alkanes), in the format: CnH2n+2
def generate_Congeners_From_Dataframe(df):
    total_rows = len(df.index)
    pd.set_option('display.max_colwidth', -1)
    bar = progressbar.ProgressBar(max_value=total_rows)
#     all_Possible=[["SMILE","Group","Formula","Cl","F","Br","C","H","Halo"]]
    all_Possible=[]
    for row in range(0,total_rows):
            df_formula=df.loc[row,"Formula"]
            pxa=PXA(formula=df_formula)
            
            smile=pxa.generate_CH_Smile()
            numH=smile.count("H")
            numC=smile.count("C")
            all_Possible+=pxa.generate_Congeners(smile,numH,numC)
            bar.update(row)
            
    data=pd.DataFrame(all_Possible)
            
    return data

#While reading the dataframe containing all congeners corresponding to the 15 alkanes of carbon chains 10-to-25, generates a new smile containing a carboxylic acid
#The dataframe (df) must contain columns ["SMILE","Group","Formula","Cl","F","Br","C","H","Halo"]
def all_generateWithCarboxylic(df):
    total_rows = len(df.index)
    pd.set_option('display.max_colwidth', -1)
    bar = progressbar.ProgressBar(max_value=total_rows)

    for row in range(0,total_rows):
        smile=df.loc[row,"SMILE"]
        group=df.loc[row,"Group"]  
        formula=df.loc[row,"Formula"] 
        numCl=df.loc[row,"Cl"]
        numBr=df.loc[row,"Br"]
        numF=df.loc[row,"F"]
        numC=df.loc[row,"C"]
        numH=df.loc[row,"H"]
        numHalo=df.loc[row,"Halo"]
        
        pxa=PXA(smile,group,formula,numCl,numF,numBr,numC,numH,numHalo)

        smileCX=pxa.add_CarboxylicAcid()

        df.loc[row,"SMILE GENERATED"] = smileCX
        df.loc[row,"Formula"]=pxa.get_Formula()
        df.loc[row,"Cl"] =pxa.get_NumCl()
        df.loc[row,"Br"] =pxa.get_NumBr()
        df.loc[row,"F"] =pxa.get_NumF()
        df.loc[row,"H"] =pxa.get_NumH()
        df.loc[row,"C"] =pxa.get_NumC()
        df.loc[row,"Halo"] =pxa.get_NumHalo()

        bar.update(row)
    return df       

#~~~~ This section generates a new random SMILES structure from existent congener ~~~~#
# The dataframe (df) must contain columns ["SMILE","Group","Formula","Cl","F","Br","C","H","Halo"]

#Generates a new random SMILES for all congeners in the dataframe at once 
def generateRandomSMILE_fromExistent(df):
    total_rows = len(df.index)
    bar = progressbar.ProgressBar(max_value=total_rows-1)

    for row in range(0,total_rows-1):
        numCl=df.loc[row,"Cl"]
        numBr=df.loc[row,"Br"]
        numF=df.loc[row,"F"]
        numC=df.loc[row,"C"]
        numH=df.loc[row,"H"]
        pxa=PXA(numCl,numBr,numF,numC,numH)

        smileCX=pxa.generateRandomSMILE()
        df.loc[row,"SMILE GENERATED"] = smileCX
        bar.update(row)
    return df       

# Generates a new random SMILES for one structure read as an input
def individual_Input():
    smile=input("\nEnter the SMILE structure: ")
    numCl=int(input("\nEnter the number of Chlorines: "))
    numBr=int(input("Enter the number of Bromines: "))
    numF=int(input("Enter the number of Florines: "))
    pxa=PXA(smile,numCl,numBr,numF)
    print("\n"+pxa.generateRandomSMILE())       

# Generates a new random SMILES for group of structure read as inputs
def group_Input():
    num=int(input("Enter the amount of structures:"))

    for x in range(0,num):
        smile=input("\nEnter the SMILE structure: ")
        numCl=int(input("\nEnter the number of Chlorines: "))
        numBr=int(input("Enter the number of Bromines: "))
        numF=int(input("Enter the number of Florines: "))
        pxa=PXA(smile,numCl,numBr,numF)
        print("\n"+pxa.generateRandomSMILE())      


#~~~~ This section generates the  figures : chemical partitioning space, Venn diagram and histogram ~~~~#    
# Visualized using a Jupyter Notebook
def define_Chain(arr, numC):

    if numC>=10 and numC<= 13:
        arr[0]+=1
    elif numC>=14 and numC<= 17:
        arr[1]+=1
    elif numC>=18 and numC<= 25:
        arr[2]+=1

def generateVennDiagram(df):
    total_rows = len(df.index)
    arr_Cl,arr_F,arr_Br,arr_ClF,arr_ClBr,arr_FBr,arr_BrClF=[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]

    for row in range(0,total_rows):
            smile=df.iloc[row,1] 
            numCl=df.loc[row,"Cl"]
            numBr=df.loc[row,"Br"]
            numF=df.loc[row,"F"]
            numC=df.loc[row,"C"]
            
            if pd.isna(smile):
                continue
        
            if numCl>0 and numF==0 and numBr==0:
                arr_Cl[3]+=1
                define_Chain(arr_Cl,numC)
                    
            elif numCl==0 and numF>0 and numBr==0:
                arr_F[3]+=1
                define_Chain(arr_F,numC)

            elif numCl==0 and numF==0 and numBr>0:
                arr_Br[3]+=1
                define_Chain(arr_Br,numC)

            elif numCl>0 and numF>0 and numBr==0:
                arr_ClF[3]+=1
                define_Chain(arr_ClF,numC)

            elif numCl>0 and numF==0 and numBr>0:
                arr_ClBr[3]+=1
                define_Chain(arr_ClBr,numC)

            elif numCl==0 and numF>0 and numBr>0:
                arr_FBr[3]+=1
                define_Chain(arr_FBr,numC)

            elif numCl>0 and numF>0 and numBr>0:
                arr_BrClF[3]+=1
                define_Chain(arr_BrClF,numC)
                
                

    df2 = pd.DataFrame([
        [40.73,31.29,'Cl (3-11)'],
        [3.59,2.32,'Br (2-6)'],
        [11.21,8.80,'Cl (1-9) \n Br (1-6)'],
        [56.36,38.02,'Cl (4-11) \n F (1-7)'],
        [4.91,2.77,'F (1-11) \n Br (2-7)'],
        [16.99,15.88,'Br (1-6) \n Cl (1-10) \n F (1-10)']
    ], columns=['Mean','Standard Deviation','Halogens Range'], index=['Cl','Br','ClBr','ClF','FBr','BrClF'])
    
    fig = plt.figure(figsize=(15,15))
    
    ax1 = fig.add_subplot(121)
    
    v = venn3_unweighted(subsets = {'100': arr_Cl[3], '010':arr_F[3],'110':arr_ClF[3],'001': arr_Br[3], '101':arr_ClBr[3],'011': arr_FBr[3], '111': arr_BrClF[3]}, set_labels = ('Cl', 'F', 'Br'))
    v.get_patch_by_id('100').set_color('green')
    v.get_patch_by_id('101').set_color('#d17600')
    v.get_patch_by_id('001').set_color('red')
    v.get_patch_by_id('111').set_color('#9e53e0')
    v.get_patch_by_id('110').set_color('#e8da10')
    v.get_patch_by_id('011').set_color('#00d1b9')
    v.get_patch_by_id('010').set_color('#0700d1')

    for t in v.set_labels:
        t.set_fontsize(22)
        t.set_color('blue')
    for t in v.subset_labels: t.set_fontsize(20)

    ax1.set_title("")
    ax1.annotate('Cl compounds', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'),fontsize=18)
    
    ax1.annotate('ClBr compounds', xy=v.get_label_by_id('101').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'),fontsize=18)
    
    ax1.annotate('F compounds', xy=v.get_label_by_id('010').get_position() - np.array([0, -0.05]), xytext=(70,70),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'),fontsize=18)
    
    ax1.annotate('ClF compounds', xy=v.get_label_by_id('110').get_position() - np.array([0.05, 0]), xytext=(-100,50),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'),fontsize=18)
    
    ax1.annotate('FBr compounds', xy=v.get_label_by_id('011').get_position() - np.array([-0.08, -0.01]), xytext=(110,40),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'),fontsize=18)
    
    ax1.annotate('BrClF compounds', xy=v.get_label_by_id('111').get_position() - np.array([-0.1, 0.10]), xytext=(110,-50),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='gray'),fontsize=18)
    
    ax1.annotate('Br compounds', xy=v.get_label_by_id('001').get_position() - np.array([0, 0.05]), xytext=(80,-20),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0',color='gray'),fontsize=18)
     
    plt.show()

def generateHistogram(df):
    total_rows = len(df.index)
    arr_Cl,arr_F,arr_Br,arr_ClF,arr_ClBr,arr_FBr,arr_BrClF=[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]

    for row in range(0,total_rows):
            smile=df.iloc[row,1] 
            numCl=df.loc[row,"Cl"]
            numBr=df.loc[row,"Br"]
            numF=df.loc[row,"F"]
            numC=df.loc[row,"C"]
            
            if pd.isna(smile):
                continue
        
            if numCl>0 and numF==0 and numBr==0:
                arr_Cl[3]+=1
                define_Chain(arr_Cl,numC)
                    
            elif numCl==0 and numF>0 and numBr==0:
                arr_F[3]+=1
                define_Chain(arr_F,numC)

            elif numCl==0 and numF==0 and numBr>0:
                arr_Br[3]+=1
                define_Chain(arr_Br,numC)

            elif numCl>0 and numF>0 and numBr==0:
                arr_ClF[3]+=1
                define_Chain(arr_ClF,numC)

            elif numCl>0 and numF==0 and numBr>0:
                arr_ClBr[3]+=1
                define_Chain(arr_ClBr,numC)

            elif numCl==0 and numF>0 and numBr>0:
                arr_FBr[3]+=1
                define_Chain(arr_FBr,numC)

            elif numCl>0 and numF>0 and numBr>0:
                arr_BrClF[3]+=1
                define_Chain(arr_BrClF,numC)
    labels = [' '*9 +'Short Chain'+'\n'+' '*9 +'(C=10-13)',' '*9 +'Medium Chain'+'\n'+' '*9 +'(C=14-17)',' '*12 +'Long Chain'+'\n'+' '*12 +'(C=18-25)',' '*18 +'Total']

    x = np.arange(len(labels))  # the label locations
    width = 0.10  # the width of the bars
    my_dpi=96
    fig, ax = plt.subplots()
    rects1 = ax.bar(x + (width+0.02), arr_Cl, width, label='Cl')
    rects2 = ax.bar(x + (width*2 +0.02), arr_F, width, label='F')
    rects3 = ax.bar(x + (width*3 +0.04), arr_Br, width, label='Br')
    rects4 = ax.bar(x + (width*4 +0.06), arr_ClF, width, label='ClF')
    rects5 = ax.bar(x + (width*5 +0.08), arr_ClBr,  width, label='ClBr')
    rects6 = ax.bar(x + (width*6 +0.10), arr_FBr, width, label='FBr')
    rects7 = ax.bar(x + (width*7 +0.12), arr_BrClF, width, label='BrClF')

    ax.set_xticks(x)
    ax.set_xticklabels(labels,fontsize = 15, ha= 'left' )
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     axes.set_xlim([xmin,xmax])
    ax.set_ylim([0,520])
    ax.set_yticks([0, 100, 300, 500])
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(22)
    ax.legend(fontsize=18)

    
    fig.set_size_inches(18.5, 9.5)
    
    rects_T=[rects1,rects2,rects3,rects4,rects5,rects6,rects7]
    
    for rects in rects_T:
        for rect in rects:
                height = rect.get_height()
                ax.annotate('{}'.format(height),
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 6),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom',fontsize=22 )
    
    fig.savefig("AC-BAP_Histogram_Carboxylic.png", dpi=600)

    fig.tight_layout()

def generate_Partitioning_Space(df):
    fig=px.scatter(df,x="logKoa", y="logKaw", color='Halo', hover_data=['C','Cl','Br','F'],range_x=[-5,90], range_y=[-45,22], category_orders={"Group": ['Cl','BrCl','Br','BrClF','ClF','BrF']},color_discrete_map={
                    "Cl": "green",
                    "BrCl": "#d17600",
                    "Br": "red",
                    "BrClF": "#9e53e0",
                    "ClF": "#e8da10",
                    "BrF": "#186aed",
                    "F": "#60d104"},opacity=0.7)
    fig.add_annotation(
            x=8.5,
            y=-2,
            xref="x",
            yref="y",
            text="<b>AC-BAP Compositions= 733</b>",
            showarrow=True,
            font=dict(
                family="Courier New, monospace",
                size=12,
                color="#ffffff"
                ),
            align="center",
            arrowhead=2,
            ay=49,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor="black",
            bordercolor="#6130db",
            borderwidth=2,
            borderpad=4,
            bgcolor="#865fe8",
            opacity=0.9
            )
    fig.add_annotation(
            x=50,
            y=10,
            xref="x",
            yref="y",
            text="<b>Total Compositions= 184,584</b>",
            showarrow=False,
            font=dict(
                family="Courier New, monospace",
                size=12,
                color="#ffffff"
                ),
            align="center",
            bordercolor="#bd4800",
            borderwidth=2,
            borderpad=4,
            bgcolor="#f26d00",
            opacity=0.9
            )

    fig.add_vline(x=6, line_width=3, line_color="#1746d1")
    fig.add_annotation(
            x=5,
            y=-38,
            xref="x",
            yref="y",
            text="<b> Log Koa = 6</b>",
            showarrow=False, textangle=-90,
            font=dict(
                family="Arial Black",
                size=12,
                color="black"
                )
            )
    fig.add_vline(x=12, line_width=3, line_color="#1746d1")
    fig.add_annotation(
            x=11,
            y=-37.5,
            xref="x",
            yref="y",
            text="<b> Log Koa = 12</b>",
            showarrow=False,textangle=-90,
            font=dict(
                family="Arial Black",
                size=12,
                color="black"
                )
            )
    fig.add_hline(y=0.5, line_width=3, line_color="#00ab23")
    fig.add_annotation(
            x=83,
            y=2,
            xref="x",
            yref="y",
            text="<b>Log Kaw = 0.5</b>",
            showarrow=False,
            font=dict(
                family="Arial Black",
                size=12,
                color="black"
                )
            )
    fig.add_hline(y=-5, line_width=3, line_color="#00ab23")
    fig.add_annotation(
            x=83,
            y=-3,
            xref="x",
            yref="y",
            text="<b>Log Kaw = -5</b>",
            showarrow=False,
            font=dict(
                family="Arial Black",
                size=12,
                color="black"
                )
            )
    fig.add_shape(type="line",
        x0=53.5, y0=-45, x1=-5, y1=13.8,
        line=dict(
            color="#de2c2c",
            width=3,
        )
                 )
    fig.add_annotation(
            x=40,
            y=-35,
            xref="x",
            yref="y",
            text="<b>Log Kow=3.5 </b>",
            showarrow=False,textangle=35,
            font=dict(
                family="Arial Black",
                size=12,
                color="black"
                )
            )
    fig.add_shape(type="line",
        x0=48, y0=-45, x1=-5, y1=8.546,
        line=dict(
            color="#de2c2c",width=3
        )
                 )
    fig.add_annotation(
            x=48.5,
            y=-38,
            xref="x",
            yref="y",
            text="<b>Log Kow=8.5 </b>",
            showarrow=False,textangle=35,
            font=dict(
                family="Arial Black",
                size=12,
                color="black"
                )
            )
    
    fig.show()

def smileToMol(smile):
  '''Convert SMILES to rdkit.Mol with 3D coordinates'''
  mol = Chem.MolFromSmiles(smile)
  if mol is not None:
      mol = Chem.AddHs(mol)
      AllChem.EmbedMolecule(mol)
      AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
      return mol
  else:
      return None 
    
def getFormat(smile,format):
  mol = pybel.readstring("smi", f"{smile}")
  mol.make3D()
  text = mol.write(f"{format}")
  return text

@app.route("/zip/<smiles_file_data>/<format>",methods = ['GET','POST'])
def getZIP(smiles_file_data,format):
  lines = smiles_file_data.split(",")
  lines.pop()
  output = BytesIO()
  with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as zipf:
      cont=0
      for line in lines:
        cont+=1
        smile=line
        format_text = getFormat(smile,format) 
        zipf.writestr(f'smile_{cont}.{format}', format_text)
  output.seek(0)
  return send_file(output,
                   mimetype="application/zip",
                   attachment_filename="test.zip",
                   as_attachment=True)


@app.route("/format/<smile>_<format>")
def smileToFormat(smile,format):
  text = getFormat(smile,format)
  response = jsonify(data=text)
  response.headers.add('Access-Control-Allow-Origin', '*')
  return response

@app.route("/download/<smile>/<format>")
def downloadFormat(smile,format):
  text = getFormat(smile,format)
  return Response(
    text, 
    mimetype=f'text/{format}',
     headers={
       "Content-disposition": f"attachment; filename={smile}.{format}",
       'Access-Control-Allow-Origin': '*'
       }
  )
  
@app.route("/3d/<smile>")
def get_3d_object(smile):
  mol = smileToMol(smile)
  mol_str = Chem.MolToMolBlock(mol)
  return Response(
    mol_str, 
    mimetype='text/sdf',
     headers={
       "Content-disposition": f"attachment; filename={smile}.sdf",
       'Access-Control-Allow-Origin': '*'
       }
  )
@app.route("/congeners/<smile>")
def getAlkanesCongeners(smile):
    numC = smile.count("C")
    numH=2*numC+2
    formula = f"C{numC}H{numH}"
    pxa=PXA(formula=formula)
    smile=pxa.generate_CH_Smile()
    all_Possible = [["SMILE","Group","Formula","Cl","F","Br","C","H","Halo"]]
    all_Possible += pxa.generate_Congeners(smile,numH,numC)
    df = pd.DataFrame(all_Possible)
    text = df.to_string()
    format="txt"
    print(text)
    return Response(
        text, 
        mimetype=f'text/{format}',
        headers={
        "Content-disposition": f"attachment; filename={smile}.{format}",
        'Access-Control-Allow-Origin': '*'
        }
    )
@app.route("/downloadCV")
def downloadCv():
    pass
if __name__ == '__main__':
  app.run(debug=True, port=8000)