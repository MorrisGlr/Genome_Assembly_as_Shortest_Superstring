f = open('rosalind_genome_assembly_as_shortest_superstring.txt', 'r')
rawdata = f.readlines()
data = []
for i in rawdata:
    data.append(i.strip('\n'))
fastanames = []
fastanamesindexes = []
dataindexes = []
for i in data:
    dataindexes.append(data.index(i))
    if '>' in i:
        fastanames.append(i)
for i in fastanames:
    index = data.index(i)
    fastanamesindexes.append(index)
dict ={}
for i in fastanamesindexes:
    y = data[i]
    dict[y]={'seq':'', 'connections':[]}
    dataindexes.pop(0)
    countforpop = 0
    for i in dataindexes:
        if i not in fastanamesindexes:
            dict[y]['seq'] = (dict[y]['seq']) + data[i]
            countforpop +=1
            continue
        if i in fastanamesindexes:
            for i in range(0,countforpop):
                dataindexes.pop(0)
            break
#The code above loads the fasta file sequences into a dicitonary in preparation for
#further manipulation


for i in dict:
    tail_name = i
    tail_seq = dict[i]['seq']
    for s in dict:
        if s != tail_name:
            tail_seq2 = tail_seq
            head_name = s
            head_seq = dict[s]['seq']
            while len(head_seq) > len(head_seq)/2:
                head_seq = head_seq[:-1]    #deletes last character in string/sequence
            while len(tail_seq2) > len(tail_seq2)/2:
                tail_seq2 = tail_seq2[1:]     #deletes first character in string/sequence
            if tail_seq2 == head_seq and len(tail_seq2)==3 and len(head_seq)==3: #graph of k=3!!!
                if head_name not in dict[tail_name]['connections']: #stops duplicate additions
                    dict[tail_name]['connections'].append(head_name)
