from collections import defaultdict
from collections import Counter
import numpy as np
import re

class functions:

    def nucleotide_Count(DNA_String):
        counts = Counter(DNA_String)
        print(str(counts['A']) + " " + str(counts['C']) + " " + str(counts['G']) + " " + str(counts['T']))

    def DNAtoRNA(DNA_String):
        return DNA_String.replace('T','U')

    def DNA_ReverseComplement(DNA_String):
        reverse = DNA_String[::-1]
        complement = ""
        for char in reverse:
            if(char == 'A'):
                complement = complement +"T"
            if(char == 'T'):
                complement = complement +"A"
            if(char == 'C'):
                complement = complement +"G"
            if(char == 'G'):
                complement = complement +"C"

        return complement

    def Rabbit_Reccurence(n, k):
        if(n == 1):
            return 1
        if(n == 2):
            return k
        if (n <= 4):
            return (Rabbit_Reccurence(n-1,k) + Rabbit_Reccurence(n-2,k))

        return Rabbit_Reccurence(n-1,k) + k*Rabbit_Reccurence(n-2,k)

    def Mortal_Rabbits(n, m):
        bunnies = [1,1]
        months = 2

        while(months<n):
            if(months<m):
                bunnies.append(bunnies[-2] + bunnies[-1])
            elif months == m or months == m+1:
                bunnies.append(bunnies[-2] + bunnies[-1] -1)
            else:
                bunnies.append(bunnies[-2] + bunnies[-1] - bunnies[-(m+1)])
            months += 1

        print(bunnies[-1])

    def GC_Content(FASTA_String):
        id = FASTA_String.splitlines()[0][1:]
        sequence = "".join(FASTA_String.splitlines()[1:])
        content = (Counter(sequence)['C'] + Counter(sequence)['G'])/(len(sequence)) * 100
        print(id)
        print(content)

    @staticmethod
    def DNA_Motif_Find(inputstring, sub):
        start = 0
        while True:
            start = inputstring.find(sub,start)
            if start == -1: return
            yield start
            start += 1

    def DNA_Motif(inputstring, sub):
        idxlist = list(functions.DNA_Motif_Find(inputstring,sub))
        result = ""
        for ans in idxlist:
            result = result + str(ans+1) + " "
        print(result)

    def HammingDistance(DNA1, DNA2):
        distance = 0
        for idx in range(0, len(DNA1)):
            if(DNA1[idx] != DNA2[idx]):
                distance += 1
        return distance

    def DNA_Profile(file):
        """ Solution for Consensus and Profile """
        fastaDict = functions.readFasta(file)
        DNA_Strings = fastaDict.values()
        for strings in DNA_Strings:
            DNA_Array.append(list(strings))
        DNA_Array = np.array(DNA_Array)
        output_Dict = defaultdict(list)
        nucleotides = ['A', 'C', 'G', 'T']
        consensus_String = ""
        columnCounts = []

        for nucleotide in nucleotides:
            output_Dict[nucleotide] = [0]*(len(DNA_Array) + 1)

        for idx in range(0, len(DNA_Array) + 1):
            column = DNA_Array[:, idx]
            columnCounts.append(Counter(column))

        for idx in range(0,len(columnCounts)):
            consensus_String = consensus_String + columnCounts[idx].most_common(1)[0][0]
            counter = columnCounts[idx]
            for ele in counter:
                output_Dict[ele][idx] = counter[ele]

        print(consensus_String)
        print(dict(output_Dict))

    @staticmethod
    def readFasta(file):
        raw = file.readlines()
        fastaDict = {}
        DNA = ""
        for idx in range(0,len(raw)):
            if raw[idx][0] == '>':
                id = raw[idx][1:]
            else
                DNA = DNA + raw[idx]
            if raw[idx+1][0] == '>':
                fastaDict[id] = DNA
                DNA = ""
        return fastaDict
