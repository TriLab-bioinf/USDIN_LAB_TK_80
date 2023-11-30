filename = input("Enter the name of the alignment file: ")
data = open(filename, 'r')
alignDict = {}
nameDict = {}
for line in data:
    # collect the alignments
    if line[:5] == "Query":
        field = line.strip().split()
        if field[0] in alignDict:
            alignDict[field[0]] = alignDict[field[0]] + field[2]
        else:
            alignDict[field[0]] = field[2]
    # map the queryID to local ID
    elif line[:7] == "Subject":
        field = line.strip().split()
        name = field[1].split(":")
        qId = field[3].split("|")
        nameDict[qId[1]] = name[1]
    else: pass
data.close()

# find the position of all the Cs in CpGs
alignto = alignDict["Query"].upper()
del alignDict["Query"]
isC = []
lenQuery = len(alignto)
for i in range(len(alignto)):
    if alignto[i] == "C":
        isC.append(i)
    else: pass


isPerfect = {}
notPerfect = {}
notPerfectlarge = {}
largedel = {}
for key, seq in alignDict.items():
    if len(seq) < 410:
        largedel[key] = seq
    elif seq[402] == '-' or seq[410] == '-':
        notPerfectlarge[key] = seq
    elif seq[403:409] == '------':
        isPerfect[key] = seq
    else:
        notPerfect[key] = seq

smalldelsize = []
for key, seq in notPerfectlarge.items():
    delcount = 0
    for base in seq:
        if seq[i] == '-':
            delcount += 1
    smalldelsize.append(delcount)
    
delsize = []
for key, seq in largedel.items():
    delsize.append(len(alignto)-len(seq))
    
# set up lists to later calculate % methylation at each position
pCpcnt = [0] * len(isC)
ptpcnt = [0] * len(isC)
ipCpcnt = [0] * len(isC)
iptpcnt = [0] * len(isC)
ldCpcnt = [0] * len(isC)
ldtpcnt = [0] * len(isC)

# count the number of alleles with a given % of methylation in bins of 10
pallPcnt = [0] * 10
ipallPcnt = [0] * 10
ldallPcnt = [0] * 10

# count the number of unique alleles
pallType = {} 
ipallType = {} 
ldallType = {} 

# then see if those positions are C or T in the test sequences
for id in sorted(isPerfect):
    checkC = isPerfect[id]
    # equalise the lengths of the test sequences
    if len(checkC) != lenQuery:
        checkC = checkC + ((lenQuery - len(checkC)) * "X")
    else: pass
    # collect the C or t value and update % methylation at each position
    pcCount = ""
    for i, posn in enumerate(isC):
        if checkC[posn] == ".":
            pcCount += "C"
            pCpcnt[i] +=1
        elif checkC[posn] == "T":
            pcCount += "t"
            ptpcnt[i] +=1
        else:
            pcCount += "-"
    # get the pcent methylation
    pnumC = pcCount.count("C")
    pnumt = pcCount.count("t")
    ppcent = pnumC/(pnumC + pnumt)
    ppcentIndex = int(10 * ppcent)
    if ppcentIndex == 10:
        ppcentIndex = 9
    else: pass
    pallPcnt[ppcentIndex] +=1
    # check for allele uniqueness
    if pcCount in pallType:
        pallType[pcCount] +=1
    else:
        pallType[pcCount] =1

# calculate % methylation at each position and prepare output
cpgNum = ""
pcpgVal = ""
for i in range(len(pCpcnt)):
    cpgNum += str(i + 1) + "\t"
    if pCpcnt[i] + ptpcnt[i] == 0:
        pcpgVal += str('none') + '\t'
    else:
        pcpgVal += str(round(pCpcnt[i]/(pCpcnt[i] + ptpcnt[i]), 2)) + "\t"

for id in sorted(notPerfect):
    checkC = notPerfect[id]
    # equalise the lengths of the test sequences
    if len(checkC) != lenQuery:
        checkC = checkC + ((lenQuery - len(checkC)) * "X")
    else: pass
    # collect the C or t value and update % methylation at each position
    ipcCount = ""
    for i, posn in enumerate(isC):
        if checkC[posn] == ".":
            ipcCount += "C"
            ipCpcnt[i] +=1
        elif checkC[posn] == "T":
            ipcCount += "t"
            iptpcnt[i] +=1
        else:
            ipcCount += "-"
    # get the pcent methylation
    ipnumC = ipcCount.count("C")
    ipnumt = ipcCount.count("t")
    ippcent = ipnumC/(ipnumC + ipnumt)
    ippcentIndex = int(10 * ippcent)
    if ippcentIndex == 10:
        ippcentIndex = 9
    else: pass
    ipallPcnt[ippcentIndex] +=1
    # check for allele uniqueness
    if ipcCount in ipallType:
        ipallType[ipcCount] +=1
    else:
        ipallType[ipcCount] =1

# calculate % methylation at each position and prepare output
ipcpgNum = ""
ipcpgVal = ""
for i in range(len(ipCpcnt)):
    ipcpgNum += str(i + 1) + "\t"
    if ipCpcnt[i] + iptpcnt[i] == 0:
        ipcpgVal += str('none') + '\t'
    else:
        ipcpgVal += str(round(ipCpcnt[i]/(ipCpcnt[i] + iptpcnt[i]), 2)) + "\t"

for id in sorted(largedel):
    checkC = largedel[id]
    # equalise the lengths of the test sequences
    if len(checkC) != lenQuery:
        checkC = checkC + ((lenQuery - len(checkC)) * "X")
    else: pass
    # collect the C or t value and update % methylation at each position
    ldcCount = ""
    for i, posn in enumerate(isC):
        if checkC[posn] == ".":
            ldcCount += "C"
            ldCpcnt[i] +=1
        elif checkC[posn] == "T":
            ldcCount += "t"
            ldtpcnt[i] +=1
        else:
            ldcCount += "-"
    # get the pcent methylation
    ldnumC = ldcCount.count("C")
    ldnumt = ldcCount.count("t")
    ldpcent = ldnumC/(ldnumC + ldnumt)
    ldpcentIndex = int(10 * ldpcent)
    if ldpcentIndex == 10:
        ldpcentIndex = 9
    else: pass
    ldallPcnt[ldpcentIndex] +=1
    # check for allele uniqueness
    if ldcCount in ldallType:
        ldallType[ldcCount] +=1
    else:
        ldallType[ldcCount] =1

# calculate % methylation at each position and prepare output
ldcpgNum = ""
ldcpgVal = ""
for i in range(len(ldCpcnt)):
    ldcpgNum += str(i + 1) + "\t"
    if ldCpcnt[i] + ldtpcnt[i] == 0:
        ldcpgVal += str('none') + '\t'
    else:
        ldcpgVal += str(round(ldCpcnt[i]/(ldCpcnt[i] + ldtpcnt[i]), 2)) + "\t"

outfile = open("Methylation_statistics.txt", 'w')
print("Number of alleles: ", len(alignDict), file=outfile)
print("Number of perfect deletions: ", len(isPerfect), file=outfile)
print("Number of imperfect deletions: ", len(notPerfect), file=outfile)
print('\t', 'Size of each deletion:', smalldelsize, file=outfile)
print("", file=outfile)
print("Number of large deletions: ", len(largedel), file=outfile)
print('\t', 'Size of each deletion:', delsize, file=outfile)
print("", file=outfile)
print("Percent methylation of each CpG", file=outfile)
print(cpgNum, file=outfile)
print(pcpgVal, file=outfile)
print(ipcpgVal, file=outfile)
print(ldcpgVal, file=outfile)
print("", file=outfile)
print("Number of alleles with methylation percent", file=outfile)
print("<10\t<20\t<30\t<40\t<50\t<60\t<70\t<80\t<90\t90+", file=outfile)
methpcent = ""
for value in pallPcnt:
    methpcent += str(value) + "\t"
print(methpcent, file=outfile)
methpcent = ""
for value in ipallPcnt:
    methpcent += str(value) + "\t"
print(methpcent, file=outfile)
methpcent = ""
for value in ldallPcnt:
    methpcent += str(value) + "\t"
print(methpcent, file=outfile)
print("", file=outfile)
print("Unique alleles found", file=outfile)
print("Frequency\tNumt\tNumdash\tNumC\tSequence", file=outfile)
print('Perfect deletions:', file=outfile)
sortAllele = sorted(pallType.items(), key=lambda x:x[1], reverse=True)
for allele in sortAllele:
    outLine = str(allele[1]) + "\t" + str(allele[0].count("t")) + "\t" + str(allele[0].count("-")) + "\t" + str(allele[0].count("C"))
    for char in allele[0]:
        outLine += "\t" + char
    print(outLine, file=outfile)
print("", file=outfile)
print('Imperfect deletions:', file=outfile)
sortAllele = sorted(ipallType.items(), key=lambda x:x[1], reverse=True)
for allele in sortAllele:
    outLine = str(allele[1]) + "\t" + str(allele[0].count("t")) + "\t" + str(allele[0].count("-")) + "\t" + str(allele[0].count("C"))
    for char in allele[0]:
        outLine += "\t" + char
    print(outLine, file=outfile)
print("", file=outfile)
print('Large deletions:', file=outfile)
sortAllele = sorted(ldallType.items(), key=lambda x:x[1], reverse=True)
for allele in sortAllele:
    outLine = str(allele[1]) + "\t" + str(allele[0].count("t")) + "\t" + str(allele[0].count("-")) + "\t" + str(allele[0].count("C"))
    for char in allele[0]:
        outLine += "\t" + char
    print(outLine, file=outfile)
outfile.close()
