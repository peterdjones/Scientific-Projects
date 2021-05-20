#script to carry out regular expression to identify a particular sequence motif (Heparin binding motif) in a protein sequence

import re
from colorama import Fore
import numpy as np

#set up our protein sequence
proteinSeq = "MATTAAGLYRRGRGYSGSGRRARALLABBBAAAAARRDRRD"
#regular expression to find Heparin binding motif XBBXBX, where B is a basic amino acid(R,L,H) and X is any amino acid
matches = re.finditer(r".[RLH][RLH].[RLH]", proteinSeq)

#define our array for storing matched motifs and locations
matchArray = []

#set printPos at 0 so printing starts at beginning of protein sequence
printPos = 0 
for match in matches:
    matchMotif = match.group() #matched sequence
    startPos = match.start() #start point of match
    endPos = match.end() #end point of match
    print(proteinSeq[printPos:startPos] + Fore.RED + proteinSeq[startPos:endPos] + Fore.RESET, end='') #print from current point in black, print matched sequence in red, then restore to black
    printPos = endPos #reset current position through protein sequence
    
    matchArray.append(np.array([matchMotif, startPos, endPos])) #append matched motif and locations to matchArray
    
print(proteinSeq[printPos:(len(proteinSeq))]) #print rest of the sequence after the final matched motif 

matchArray = np.asarray(matchArray) #convert matched data to a numpy array for printing

#loop through array, printing motif, and start and end positions of match
for x in range(matchArray.shape[0]):
    print('Motif: ' + matchArray[x,0] + ', Start Position: ' + matchArray[x,1] + ', End position:' + matchArray[x,2])
    