import re

filename = "NR2BHRNB114-Alignment.txt"
data = open(filename, 'r')
alignDict = {}
nameDict = {}

# Query regex
Q = re.compile(r'^Query\s+\d+\s+([ACGT-]+)\s{2}\d+$')
query_len = 1

# alignment chunk counter
chunk = {} # chunk number using query_id as key
chunk_len = [] # records chunk lengths

for line in data:
    # parse line

    ##############################################################################
    # NOTE: Keep track of alignment chunks where subject sequence contains all gaps and 
    # therefore it does not appear in that chunk. Hence, it is necessary to add a
    # blank line to the missing subject sequence with the same lenght as the query
    ##############################################################################

    # Query line regex (it will depend upon Query sequence en in the alignment, which is variable)
    if Q.match(line):
        query_len = len(Q.match(line).group(1))

    # Alignment line regex
    q = re.compile(r'^(Query_?\d*)\s+(\d+)\s+(.{' + f'{query_len}' + '})\s{2}(\d+)$')
    # Subject (header) line regex
    s = re.compile(r'Subject\s#\d+:(\d+)\sID:\slcl\|(Query_\d+)')
    g = re.compile(r'\s', )

    # collect the alignments 
    if q.match(line):
        query_id = q.match(line).group(1)
        q_start = q.match(line).group(2)
        sequence = q.match(line).group(3)

        # Replace blank spaces with "_"
        sequence = g.sub('_', sequence)

        # Kepp track of chunk number
        chunk[query_id] = 1 + chunk.get(query_id, 0)

        if query_id == "Query":
            chunk_len.append(query_len) 

        # Check if the new subject chuck number is out of sync with the Query chunk number
        # If so, this means that subject sequence was missing in previous chunk
        if chunk[query_id] < chunk["Query"]:
            gap = ''
            for i in range(chunk[query_id]-1, chunk["Query"]-1):
                chunk[query_id] += 1
                gap += '_' * chunk_len[i]
            
            sequence = gap + sequence

        if query_id in alignDict:
            alignDict[query_id] += sequence
        else:
            alignDict[query_id] = sequence
    
    # map the queryID to local ID
    elif s.match(line):
        name = s.match(line).group(1)
        qId  = s.match(line).group(2)
        nameDict[qId] = name
        chunk[qId] = 0

    else: pass

data.close()

# The code below is just for printing the sequences captured from the alignment and can be removed.
for q_id in nameDict.keys():
    print(f">{q_id} length={len(alignDict[q_id])}")
    print(alignDict[q_id])

