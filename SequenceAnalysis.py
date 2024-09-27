import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
class NucParams:
    '''Take in the .fa file given to us and calculate the sequence length, GC content, codon frequence, and number of codons'''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''Instaniate every dictionary you need for the program, nucComp, codonComp, and aaComp'''
        self.inString = inString
        self.nucComp = {'A': 0, 'C': 0, 'T':0, 'G': 0, 'U': 0, 'N': 0}
        self.codonComp = {key:0 for key, value in self.rnaCodonTable.items()}
        self.sorted_dict = dict(sorted(self.rnaCodonTable.items()))
        self.noOrderaaComp = {self.rnaCodonTable[key]:0 for key, value in self.codonComp.items()}
        self.aaComp = dict(sorted(self.noOrderaaComp.items()))

    def addSequence (self, inSeq):
        #take out everything that is not ACGTU and replace it with N
        addSeq = ''.join(['N' if nuc not in 'ACGTU' else nuc for nuc in inSeq])

        #replace every T with U to make a RNA sequence
        nucleotideSeq = addSeq.replace('T','U')

        #for every nucleotide in the sequence, add to the nucComp dictionary
        for i in addSeq.upper(): #uppercase every input to make sure it is uppercase
            if i in self.nucComp:
                self.nucComp[i] += 1      
        
        #for all the nucleotides in the sequence, increment in groups of 3 for the whole sequence
        codonRNA = [nucleotideSeq[i:i+3] for i in range(0, len(nucleotideSeq), 3)]

        for i in codonRNA: #for every increment of three, codons, count in the dictionary
            if i in self.codonComp:
                self.codonComp[i] += 1
        
    def aaComposition(self):
        '''Return the amino acid dictionary under aaComposition'''
        return self.aaComp 
    def nucComposition(self):
        '''Return the nucleotide dictionary under nucComp'''
        return self.nucComp
    def codonComposition(self):
        '''Return the codon dictionary under codonComposition'''
        return self.codonComp
    def nucCount(self):
        '''Return the total count of nucleotides'''
        return sum(self.nucComp.values()) # sum all the values in the nucleotide dictionary
    
    def GCcount(self):
        '''Count all the total G and C nucleotides in sequence and then find the GC content'''
        wantGC = {'G':0, 'C':0} #instantiate the GC dictionary
        totalGC = 0 #instantiate the GC count
        for key in wantGC: # for all the GC in the nucleotide dictionary
            if key in self.nucComp:
                totalGC += self.nucComp[key] #add to total GC
        GCcontent = (totalGC/self.nucCount()) # GC count divided by the total nucleotide count
        return GCcontent
    
    def codonFreq(self):
        '''Caluculate the codon frequency, the amount of codons per amino acid'''
        codonFrequency = {} #instantiate the codonFrequency dictionary
        for key, value in self.rnaCodonTable.items(): 
            if self.aaComp[value] != 0: #if there are no aminoacids in the dictionary don't run
                codonFrequency[key] = (self.codonComp[key]/self.aaComp[value])*100 # number of codons divded by number of aminoacids
        return codonFrequency

    def printer(self):
        '''
        Print the wanted out string for every codon, aminoacid pair
        For aminoacid composition, if codonComp key has value from rnaCodonTable, increment count for that value in the aaComp
        '''
        for key, value in self.codonComp.items():
            if key in self.rnaCodonTable:
                self.aaComp[self.rnaCodonTable[key]] += value # add the value from codonComp to the aaComp value
       
        outStr = ""
        for i in self.aaComp: # for all aminoacids in aaComp
            for key, value in self.sorted_dict.items(): #for codons and aminoacids in rnaCodonTable
                # this will connect the codons from rnaCodonTable and aminoacids from aaComp
                if value == i: # if the aminoacids of rnaCodonTable = aminoacids of aaComp
                    if self.codonComp[key] != 0: # the number of codons must >= 1 or won't print
                        # this will equal "new line, codon, :, aminoacid, space, rounded CodonFrequency value, (, count of codons, )
                        outStr += "\n" + str(key) + " : " + str(value) + " " + str(round(self.codonFreq()[key], 1)) + " ( " + str(self.codonComp[key]) + " ) "
        # result will now give our sequence length, GC content, and the outStr 
        result = ("sequence length = {:.3} Mb\n\nGC content = {:.1%}\n{}".format(self.nucCount()/1000000, self.GCcount(), outStr))
        return result
    
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

    
    weight = 0.0
    totalWeight = 0.0
    mwH2O = 18.015
    
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
    
  
    def __init__ (self, protein):
        ''' Assign each parameter to be in __init__.'''
        self.protein = protein #self.protein can now be referred in every method
        
        self.aaProfile = {key:0 for key,value in self.aa2mw.items()} # I remake a new dictionary with values 0 that
        # counts the number of amino acids but keep the same keys

        # This will count each amino acid in protein and will mark +1 in the dictionary
        for char in self.protein.upper(): # turns all the inputed protein string to upper case
            if char in self.aaProfile:
                self.aaProfile[char] += 1
        
    def aaCount (self):
        
        '''I want to return the summed values of the new dictionary I made to keep count of the amino acids.'''
        return sum(self.aaProfile.values())
        
    def pI (self):
        
        '''
        Iterate through each value of pH and compare to the pH before it. Then find the next pH and
        and continuously compare to the one before it.
        '''
        
        pH = 0.0 # we want pH to start at 0.0
        pI_Return = 0.0 # we want our isoelectric point to be 0.0
        while pH <= 14.0: # we want a while loop to iterate through every value before pH is less than or equal to 14
            if pH > 0.0: # if the pH is greater than 0.0 but less than 14.
                # For value pH, input into our __charge__ method and run through that function to find
                # the best pH value for our isoelectric point(total charge = 0)
                # We compare our charge with pH to be less than our charge with pH -0.01
                # so that we find the lowest pH that give our charge to be closer to 0
                if abs(self._charge_(pH)) < abs(self._charge_(pH-.01)): 
                    pI_Return = pH # Set the pI to the pH that was valued
            pH += .01 # in the while loop add 0.01 to our pH and that will be tested until pH = 14
        return pI_Return

    def aaComposition (self) :
        '''Give us the new dictionary that will count the amount of aminoacids in the protein.'''
        return self.aaProfile 

    def _charge_ (self, pH):
        '''
        Create two dictionaries to split up postive charge count and negative charge count, iterate through our
        dictionary with the count of amino acids and if the amino acid in our positive charge dictionary, 
        do the math for each amino acid and then multiply that with the number of amino acids and then sum up
        all our those values.
        Same thing for our negative charge dictionary and then subract the positive from the negative to get
        the net charge.
        '''
        self.pH = pH
        
        posCharge = {} # create our postive dictionary
        for i in self.aaProfile.keys(): 
            # if the protein is in our profile and in our postive charge dictionary
            if i in self.aa2chargePos.keys():
                # multiply the amino acid in our postive charge dictionary with the number of total amino acids
                # as well as do the desired math for the charge
                posCharge[i] = self.aaProfile[i] * (10**self.aa2chargePos[i]/(10**self.aa2chargePos[i]+(10**self.pH)))
            # sum all the values in our new dictionary
            posSum = sum(posCharge.values())
            # add our Nterm value
            totalPos = posSum + (10**self.aaNterm/((10**self.aaNterm)+(10**self.pH)))

        negCharge = {}
        for i in self.aaProfile.keys():
            # if the protein is in our profile and in our negative charge dictionary
            if i in self.aa2chargeNeg.keys():
                # multiply the amino acid in our postive charge dictionary with the number of total amino acids
                # as well as do the desired math for the charge
                negCharge[i] = self.aaProfile[i] * (10**self.pH/((10**self.aa2chargeNeg[i])+(10**self.pH)))
            # sum all the values in our new dictionary
            negSum = sum(negCharge.values())
            # add out Cterm value
            totalNeg = negSum + (10**self.pH/((10**self.aaCterm)+(10**self.pH)))
        # now subtract each other to get total charge
        totalCharge = totalPos - totalNeg
        return totalCharge
    
    def molarExtinction (self):
        '''
        iterate through dictionary and then multiply the number of amino acids in protein with
        the value of extinction coeffiecient of each amino acid. Then sum our dictionary value.
        '''
        molarExt = {} # create an empty dictionary
        for key, val in self.aaProfile.items():
            if key in self.aa2abs280: # if the key is in both dictionaries
                # multiply the values of each dictionary
                molarExt[key] = val * self.aa2abs280[key]
            # sum the values of the new dictionary
            totalE = sum(molarExt.values())
        return totalE

    def massExtinction (self):
        '''divide our molarExtinction with molecularWeight '''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''
        iterate through the number of amino acids in our profile dictionary to multiply with the 
        values of molecular weight in our mw dictionary.
        '''
        # define weight and totalWeight in our molecularWeight
        weight = 0.0
        totalWeight = 0.0
       
        for i in self.aaProfile.keys(): # for all the amino acids in the profile dictionary
            # take the values of molecular weight associated to the amino acid, subtract 
            # to the water weight and multiply to the number of amino acids.
            weight += ((self.aa2mw[i] - self.mwH2O)* self.aaProfile[i]) 
        totalWeight += weight + self.mwH2O # our total weight is our weight plus the weight of water
        return totalWeight

class OrfFinder:
    def __init__(self, sequence, start = ['ATG','TTG','GTG'], stop = ['TAG','TAA','TGA']):
        self.sequence = sequence
        self.start = start
        self.stop = stop
        self.orfList = []
    
    def findOrfs(self):

        addSeq = ''.join(['N' if nuc not in 'ACGTU' else nuc for nuc in self.sequence])
        mainSeq = addSeq.replace('U','T')
        
        for frame in [0,1,2]:
            posStart = []

            for index in range(frame, len(mainSeq), 3):
                codon = [mainSeq[index:index+3]]

                if codon in self.start:
                        posStart.append(index)
                
                if codon in self.stop:
                        if codon in self.start:

                            starts = posStart[0] + 1 - frame
                            stops = index + 3
                            lengths = self.start - self.stop
                            self.saveOrf((frame % 3) + 1, starts, stops, lengths)
                            posStart = []

                if codon not in self.start:
                    if codon in self.stop:
                        start = 1
                        stop = index + 3
                        length = stop - start + 1
                        self.saveOrf((frame % 3) + 1, start, stop, length)
                        posStart = []
                    

            if codon not in self.stop:
                if codon in self.start:  # If no stop codon was found but start codon was found. - dangling start
            
                    #create an ORF with a dangling start codon
                    start = posStart[0] + 1
                    stop = len(mainSeq)
                    length = stop - start + 1
                    self.saveOrf((frame % 3) + 1, start, stop, length)
        
        return self.orfList
    
    def complimentarySequence(self):
        self.translateNuc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([self.translateNuc[base] for base in self.sequence[::-1]])

    def findRevOrfs(self):

        compSeq = self.complimentarySequence()

        for frame in [0,1,2]:
            posStart = []

            for index in range(frame, len(compSeq), 3):
                codon = compSeq[index:index+3]

                if codon in self.start:
                    posStart.append(index)

                if codon in self.stop:
                    if codon in self.start:
                        stop = len(compSeq) - startPos[0] #distance from start of reverse to start codon
                        start = len(compSeq) - (index + 2) #distance from end of reverse to position of last nuc, subtract from total of reverse
                        
                        #adjust stop position depending on the frame
                        if frame == 1: 
                            stop += 1
                        
                        elif frame == 2: 
                            stop += 2
                        
                        length = stop - start + 1 #calculate length after adjusting
                        
                        #save ORF and reset variables
                        #computes the reading frame of the ORF relative to the 3'-5' strand of the DNA sequence, but with a negative sign because reverse complement
                        self.saveOrf(-1 * ((frame%3) + 1), start, stop, length)
                        
                        startPos = []
                if codon not in self.start:
                    if codon in self.stop:
                        start = len(compSeq) - index - 2 #calculate start position
                        stop = len(compSeq)
                        length = stop - start + 1
                        
                        self.saveOrf(-1 * ((frame % 3) + 2), start, stop, length)
                        startPos = []

            if codon in self.start:
                start =  startPos[0] + 1
                stop = 1
                length = stop - start + 1 
                self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length)

        return self.orfList
    
    def saveOrf(self, frame, start, stop, length):
        """ Saves ORF info in following format:
        
        Adds list to the Orfs attribute, each list contains:
        - Frame: the reading frame of the ORF (0, 1, or 2)
        - Start: the index of the start codon for the ORF
        - Stop: the index of the stop codon for the ORF
        - Length: the length of the ORF in nucleotides
        """
        self.orfList.append([frame, start, stop, length]) #puts all together

import math
class FindPrimers:
    
    def __init__(self, sequence, readRanges, readLengths, optimalGC, optimalTemp, findPairs):
        self.sequence = sequence
        self.readRanges = readRanges
        self.readLengths = readLengths
        self.optimalGC = optimalGC
        self.optimalTemp = optimalTemp
        self.findPairs = findPairs
        
        self.rangeList = []
        self.sectionedLengths = []
        self.gcDict = {}
        self.tempDict = {}
        self.posDict = {}
        
        self.forRange = []
        self.forLength = []
        self.forGC = {}
        self.forTemp = {}
        self.forDict = {}
        
        self.finalList = []
        self.finalForList = []
        self.finalPairs = []
        
        #Instantiating all lists and dictionaries that are going to be manipulated throughout the code
        
        self.forSequence = self.ForwardStrand()
        nucPair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        
        #Establish the forward sequence
        
    def ForwardStrand(self):
        '''Reverse and find Compliment Strand of the File'''
        nucPair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([nucPair[base] for base in self.sequence[::-1]]) 

        #Populate the forward strand sequence
        
    def FindRange(self):
        '''Parses the input value into defined ranges that we will later manipulate'''
        readRanges = self.readRanges.split(', ')
        for index in readRanges:
            splitRanges = index.split('-')
            if len(splitRanges) == 2 and splitRanges[0] and splitRanges[1]:  
                start, end = int(splitRanges[0]), int(splitRanges[1])
                self.rangeList.append(self.sequence[start:end + 1])
                # Extract the sequence within the specified range and append to rangeList
            else:
                print("Invalid range:", index)
        #Take input ranges and iterate through the different ranges and isolate bases
                
        return self.rangeList
    
    def FindForRange(self):
        '''Finds the coordinate at which the primer latches'''
        seqLength = len(self.forSequence)
        readRanges = self.readRanges.split(', ')
        for index in readRanges:
            # Iterate through each range
            start, end = map(int, index.split('-'))
            # Calculate the positions on the forward strand
            forStart = seqLength - end - 1
            forEnd = seqLength - start
            forSequence = self.forSequence[forStart:forEnd]
            self.forRange.append(forSequence)
        return self.forRange


    def FindLength(self):
        '''the range will be iterated through 19-24 bases'''
        
        for index in self.readRanges.split(', '):
            start, end = map(int, index.split('-'))
            # Iterate through possible positions within the range
            for j in range(start, end - self.readLengths + 2):
                partition = self.sequence[j - 1:j + self.readLengths - 1]
                # Calculate start and end positions for the partition
                startPos = j
                endPos = j + self.readLengths - 1
                self.sectionedLengths.append(partition)
                # Store the start and end positions as a tuple with the partition as the key
                self.posDict[partition] = (startPos, endPos)
        return self.sectionedLengths, self.posDict
        
    def FindForLength(self):
        '''Finds the length of forward primers '''
        readRanges = self.readRanges.split(', ')
        forwardRanges = self.FindForRange()
        # Iterate through each range and corresponding sequence
        for ranges, sequence in zip(readRanges, forwardRanges):
            start, end = map(int, ranges.split('-'))
            # Adjust the start and end positions for the forward sequence
            forStart = len(self.sequence) - end
            forEnd = len(self.sequence) - start
            for j in range(forStart, forEnd - self.readLengths + 2):
                partition = sequence[j - forStart:j - forStart + self.readLengths]
                 # Calculate start and end positions for the partition
                startPos = j
                endPos = j + self.readLengths - 1
                self.forLength.append(partition)
                 # Store partition as key and its start and end positions as value in forDict
                self.forDict[partition] = (startPos, endPos)
        return self.forLength, self.forDict

    def FindGC(self):
        '''calculate the GC content of the range with length'''
        for sequence in self.sectionedLengths:
            nucDict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            # Initialize a dictionary to count nucleotides
            for base in sequence:
                if base in nucDict:
                    nucDict[base] += 1
                     # Count occurrences of each nucleotide in the sequence
            totalGC = nucDict['G'] + nucDict['C']
            totalBases = sum(nucDict.values())
            if totalBases != 0:
                GCcontent = (totalGC / totalBases) * 100
            if self.optimalGC - 4 <= GCcontent <= self.optimalGC + 5:
                # Check if GC content falls within the optimal range
                self.gcDict[sequence] = GCcontent
                
        for forward in self.forLength:
            forDict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            for base in forward:
                if base in forDict:
                    forDict[base] += 1
            totalGC = forDict['G'] + forDict['C']
            totalBases = sum(forDict.values())
            if totalBases != 0:
                GCcontent = (totalGC / totalBases) * 100
            if self.optimalGC - 4 <= GCcontent <= self.optimalGC + 5:
                self.forGC[forward] = GCcontent
        return self.gcDict, self.forGC
    
    def FindTemp(self):
        '''Find temperatures around the optimal temperature using Tm calculation'''
        for index in range(-50, 51):
            temp = round(self.optimalTemp + index * 0.1, 2)
            for sequence in self.sectionedLengths:
                nucDict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
                for base in sequence:
                    if base in nucDict:
                        nucDict[base] += 1
                # Calculate Tm using the provided formula
                NAconcentration = 0.2  
                # Na+ concentration in M
                GCcontent = (nucDict['G'] + nucDict['C']) / len(sequence) * 100
                meltTemp = 81.5 + 16.6 * math.log10(NAconcentration) + 0.41 * GCcontent - 600 / len(sequence)
                # Check if the calculated Tm is within the desired range
                if abs(meltTemp - temp) <= 1.0:
                    if sequence not in self.tempDict:
                        self.tempDict[sequence] = temp
                    else:
                        # If a temperature is already associated with the sequence, only update if the new temperature is closer to the optimal temperature
                        if abs(temp - self.optimalTemp) < abs(self.tempDict[sequence] - self.optimalTemp):
                            self.tempDict[sequence] = temp
                            
        for index in range(-50, 51):
            temp = round(self.optimalTemp + index * 0.1, 2)
            for forward in self.forLength:
                forDict = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
                for base in forward:
                    if base in forDict:
                        forDict[base] += 1
                # Calculate Tm using the provided formula
                NAconcentration = 0.2  
                # Na+ concentration in M
                GCcontent = (forDict['G'] + forDict['C']) / len(forward) * 100
                meltTemp = 81.5 + 16.6 * math.log10(NAconcentration) + 0.41 * GCcontent - 600 / len(forward)
                # Check if the calculated Tm is within the desired range
                if abs(meltTemp - temp) <= 1.0:
                    if forward not in self.forTemp:
                        self.forTemp[forward] = temp
                    else:
                        # If a temperature is already associated with the sequence, only update if the new temperature is closer to the optimal temperature
                        if abs(temp - self.optimalTemp) < abs(self.forTemp[forward] - self.optimalTemp):
                            self.forTemp[forward] = temp
                            
        return self.tempDict, self.forTemp
    
    def FindPairs(self):
        # Iterate through each pair of sequences in gcDict and tempDict
        for gcKey, gcValue in self.gcDict.items():
            for tempKey, tempValue in self.tempDict.items():
                if gcKey in tempKey:  
                    commonSequence = gcKey
                    self.finalList.append((commonSequence, "reverse", self.posDict[commonSequence], len(commonSequence), gcValue, tempValue))
        # Iterate through each pair of sequences in forGC and forTemp
        for gcKey, gcValue in self.forGC.items():
            for tempKey, tempValue in self.forTemp.items():
                if gcKey in tempKey:  
                    forSequence = gcKey
                    self.finalForList.append((forSequence, "forward", self.forDict[forSequence], len(forSequence), gcValue, tempValue))
         # Iterate through each pair of sequences in finalList and finalForList to find pairs with proper distance
        for i in self.finalList:
            for j in self.finalForList:
                 # Check if the distance between end of sequence from finalList and start of sequence from finalForList is greater than or equal to findPairs
                if abs(j[2][1] - i[2][0]) >= self.findPairs:
                    self.finalPairs.append((i, j))
                    # Append the pair to finalPairs list
                    break  # Exit the inner loop once a suitable pair is found

        return self.finalPairs


def main(inFile=None):
    myReader = FastAreader(inFile)
    sequence = ""
    for head, seq in myReader.readFasta():
        sequence += seq
    
    inputRanges = input("Enter comma-separated ranges you would like to iterate through (e.g., 1-5, 10-15): ")
    inputLengths = int(input("What is your ideal length (Ideal length is between 18-24bp)? (e.g, 20):"))
    optimalGC = float(input("What is your ideal gc content range (Ideal conditions have GC% between 40-60%)? (e.g, 45):"))
    optimalTemp = float(input("What is your ideal melting temp (Ideal conditions lie between 40-70 degrees)? (e.g, 60.0):"))
    inputPair = int(input("What is the minimum amount of base pairs you would like to have the primers read? (e.g, 100):"))
    
    
    primers = FindPrimers(sequence, inputRanges, inputLengths, optimalGC, optimalTemp, inputPair)

    ranges = primers.FindRange()
    forRange = primers.FindForRange()
    length, posDict = primers.FindLength()
    forLength, forDict = primers.FindForLength()
    gc, temp = primers.FindGC(), primers.FindTemp()

    pairs = primers.FindPairs()
    
    for reverse, forward in pairs:
        reverseStart = reverse[2][0]
        reverseEnd = reverse[2][1]
        forwardStart = forward[2][0]
        forwardEnd = forward[2][1]
        productSize = abs(forwardEnd - reverseStart)
        print(f"Pair: {reverse[0]} (Reverse) and {forward[0]} (Forward)")
        print(f"Reverse Sequence Position: {reverse[2]}, Length: {reverse[3]}, GC: {reverse[4]}%, Temp: {reverse[5]}°C")
        print(f"Forward Sequence Position: {forward[2]}, Length: {forward[3]}, GC: {forward[4]}%, Temp: {forward[5]}°C")
        print(f"Product Size: {productSize}")
        print()


if __name__ == "__main__":
    main("lambda_attB_primer_regions.fa")