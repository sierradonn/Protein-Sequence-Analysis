#!/usr/bin/env python3
# Name: Sierra Donn (sdonn)
''' The following program will calculate the pysical-chemical properties of a protein sequence

The program should look like the following when running properly:

inputs =  VLSPADKTNVKAAW
output:
Number of Amino Acids: 14
Molecular Weight: 1499.7
molar Extinction coefficient: 5500.00
mass Extinction coefficient: 3.67
Theoretical pI: 9.88
Amino acid composition:
A = 21.43%
C = 0.00%
D = 7.14%
E = 0.00%
F = 0.00%
G = 0.00%
H = 0.00%
I = 0.00%
K = 14.29%
L = 7.14%
M = 0.00%
N = 7.14%
P = 7.14%
Q = 0.00%
R = 0.00%
S = 7.14%
T = 7.14%
V = 14.29%
W = 7.14%
Y = 0.00%

'''
class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein): #added the protein attribute to be used for all the other functions
        self.protein = protein
        
        aaGroup = ''.join(protein).split() #concatenates characters in protein and splits them this removes whitespace
        self.proteinInput = ''.join(aaGroup).upper() #concatenates aa and converts to uppercase 
        self.aaComp = {} #this creates a dictionary that will store the counts of each amino acid entered from the input sequence 
        
        for aa in self.aa2mw.keys(): #loops through all the aa in the dictionary aa2mw
            self.aaComp[aa] = float(self.proteinInput.count(aa)) #will count each aa and store it in aaComp dictionary 
           
    '''the count function will count the number of amino acids in the input string 
    if you input VLSPADKTNVKAAW
    then the print result will be 14 '''

    def aaCount (self): #added the len method to count number of amino acids
        self.aaNum = 0 #initilizese the variable aaNum 
        for aa in self.protein:  #for the amino acids in the protein, loop through each one
            if aa in self.aa2mw: #if an amino acid is in the dictionary aa2wm do something to it 
                self.aaNum += 1 #what will be done is one will be added for each charcter in the dictionary 
        return self.aaNum #stores the count value in the variable aaNum to be used later 
    
    '''the composition function will call the function written in _init_ 
    this function will output percentages of aa found in the sequence ''' 
    
    def aaComposition (self) :
        return self.aaComp #returns dictionary set up in __init__ function 
   
    '''the pI and charge functions will be used to find the  theoretical isolelectric point 
    this will cahnge depending on ipnut and aa that are in the key 
    but an example is the sequence VLSPADKTNVKAAW will output Theoretical pI: 9.88 '''

    
    def pI (self):
        arbitraryCharge =  1000000 #allows for loop to run through infinite positive and update goodPH when any lower charge is found 
        goodPH = 0 
        thepH = 0

        while thepH < 14.01: #loop through all the cases where pH less than 14
            myCharge = self._charge_(thepH)#calculates the absolute value of the charge at currect pH 
            if abs(myCharge) < arbitraryCharge: #if the charge is less than the arbitrary charge 
                #which is set to a high number so the loop can iterate over a lot of values 
                #is this is true this means the current pH will result in a charge closer to zero 
                arbitraryCharge = abs(myCharge) #updates the arbitraryCharge with the current charge value 
                goodPH = thepH 
            thepH += 0.01 #makes sure the best pH value is down to 2 decimals 
        return goodPH
        
    def _charge_ (self, pH):
        posCharge = 0 #initilizes the variables 
        negCharge = 0
        nTerminous = (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)
        cTerminous = (10 ** pH) / (10 ** self.aaCterm + 10 ** pH)
        for aa in self.protein: #loops through the aa that are in the input protein
            if aa in self.aa2chargePos: #if aa are in the dictionary with the positive charges the keys will be searched 
                #countAminoAcids = self.aaComp[aa]
                posCharge += (10 ** self.aa2chargePos[aa]) / (10 ** self.aa2chargePos[aa] + 10 ** pH)
                #((self.aaComp[aa]) * (10 ** self.aa2chargePos[aa])/((10 ** self.aa2chargePos[aa]) + (10 ** pH))) 
                #nTerm = (10 ** self.aaNterm)/((10 ** self.aaNterm) + (10 ** pH))
                
            elif aa in self.aa2chargeNeg: #if aa are in the dictionary with the negative charges the keys will be searched 
                #countAminoAcids = self.aaComp[aa] 
                negCharge += (10 ** pH) / (10 ** self.aa2chargeNeg[aa] + 10 ** pH)
                #((self.aaComp[aa]) * (10 ** pH)/((10 ** self.aa2chargeNeg[aa]) + (10 ** pH)))
                #cTerm = (10 ** pH)/((10 ** self.aaCterm) + (10 ** pH)) 
        negCharge += cTerminous #adds charge to c-terminus
        posCharge += nTerminous #adds charge to n-terminus 

        netCharge = posCharge - negCharge
        return netCharge
    
    '''the molar Extinction function is used to find the molarExtinction and also the massExtinction 
    from the sequence VLSPADKTNVKAAW the out put for both are as follows molar Extinction coefficient: 5500.00
    mass Extinction coefficient: 3.67'''

    def molarExtinction (self): 
        if "W" in self.protein or "Y" in self.protein or "C" in self.protein: # check if the letters W,Y, or C are in the protein
            y = self.protein.count("Y") #initilize the variables to use in calculations 
            w = self.protein.count("W") #these are the count value for the letters y,w,c respectively
            c = self.protein.count("C")
            aaList = ["Y", "W", "C"] #make a list of the three amino acids that have absorbance at 280 nm 
            
            calculate = 0 #initilize the variable 
            for aa, count in zip(aaList, [y, w, c]): #using the zip function i combined elements from aaList and [y, w, c]
                #this creates a tuple to be able to connect the values in the dictionary to the count of letter int the sequence 
                calculate += count * self.aa2abs280.get(aa, 0) #for each aa found in absorbance dictionary
                #the product of letter count and absorbance for that specific letter are added to the variable calculate  
            return calculate
        else: #if the aa from the sequence can not be found zero will be returned 
            return 0
       
    
    def massExtinction (self):
        pass
        myMW =  self.molecularWeight() 
        return self.molarExtinction() / myMW if myMW else 0.0


    '''using the dictionary of aa weights this function will sum the averages of the aa in the 
    given sequence and also acount for water loss during bonding
    from VLSPADKTNVKAAW the Molecular Weight will be  1499.7''' 

    def molecularWeight (self):
        if self.aaNum == 0: #if there is no entry of sequence then the weight will be zero  
            print("0")
        else: 
            self.weight = 0 #initilize the variables 
            loss = self.aaNum - 1 #amount of amino acids minus one will give us the water lost in bonding 
            for aa in self.protein: #for the aa in the sequence 
                self.weight += self.aa2mw.get(aa, 0) #updates and adds the weight if the aa is found in the dictionary 
                #from the dictionary the molecular weights will be added 
                
            return self.weight - (loss * self.mwH2O) #this will store the massExtinction
            #loss times the molecular weight of water gives the molecular weight of water loss 
            #must subtract the molecular weight of water loss from the weight of the sequence to get total weight 
    

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?').upper()
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()
    

