f = open('rosalind_genome_assembly_as_shortest_superstring.txt', 'r')
answercheck = 'ATTAGACCTGCCGGAATAC'
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
    dict[y]={'seq':''}
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


recursive_seq1 = []
for i in dict:
    recursive_seq1.append(dict[i]['seq'])


#now append the joinedread to the list and make the recursive assembler
print(recursive_seq1)


global should_restart
global end_of_tail_loop
global end_of_head_loop
recursive_seq2 =[]
final_assembly = []
whole_loop_counter = 0
joincount = 0



should_restart = True
while should_restart:
    recursive_seq_1_active = False
    recursive_seq_2_active = False
    should_restart = False
    removeread_list = []    #the list takes into account the reads that have joined to other reads so that they are not duplicated downstream
    active_data = []    #This list will be populated by the recursive_seqx that contains data, i.e. len > 0. Using a list list recursively in addition to ading and removing elements is problematic
    if len(recursive_seq1) > 0:
        active_data = recursive_seq1
        recursive_seq_1_active = True
        print('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/R1 active'+' | join count: '+str(joincount)+'\n'+'R1 = '+str(recursive_seq1))
        lenofactive_data = len(recursive_seq1)
    elif len(recursive_seq2) > 0:
        active_data = recursive_seq2
        recursive_seq_2_active = True
        print('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/R2 active'+' | join count: '+str(joincount)+'\n'+'R1 = '+str(recursive_seq2))
        lenofactive_data = len(recursive_seq2)
    for i in active_data:
        tail_i = active_data.index(i)
        tail_seq = i #sequence for manipulation
        tail_seq_ref = i #full sequence string for downstream reference
         #length of the list before joined reads are appended to this list. This value iff used to clear the

        while len(tail_seq) > (len(tail_seq_ref)/2):
            tail_seq = tail_seq[1:]
            print('chopping tail seq('+str(tail_i)+')    : '+tail_seq)

            for k in active_data:
                head_i = active_data.index(k)
                head_seq_ref = k
                head_seq = k
                if head_i != tail_i:
                    while len(head_seq) > (len(head_seq_ref)/2):
                        head_seq = head_seq[:-1]
                        print('chopping head seq    ('+str(head_i)+'): '+head_seq)
                        #print('line 107 | boolean tail: '+str(end_of_tail_loop)+' | '+'boolean head: '+str(end_of_head_loop)+' | '+ 'join count: '+str(joincount)+' | '+'Whole cycle count: '+str(whole_loop_counter))
                        if tail_seq != head_seq:
                            if recursive_seq_1_active == True:          #from this line to 5 lines below, adds the unjoined reads to the inactive
                                if tail_seq_ref not in recursive_seq2 or head_seq_ref not in recursive_seq2:
                                    if head_seq_ref not in recursive_seq2:
                                        recursive_seq2.append(head_seq_ref)
                                        print('sequences transfering over to R1: '+head_seq_ref)
                                        #if tail_i == len(active_data)-2: #checks if the last list item of the first FOR loop has been iterated
                                            #
                                            #print('-----Reached TRUE in end of TAIL loop')
                                    if tail_seq_ref not in recursive_seq2:
                                        recursive_seq2.append(tail_seq_ref)
                                        print('sequences transfering over to R1: '+tail_seq_ref)
                                    if head_i == len(active_data)-2 and tail_i ==len(active_data)-1 and len(head_seq)<=((len(head_seq_ref)/2)+1): #checks if the last list item of the first FOR loop has been iterated
                                        end_of_head_loop = True
                                        print('-----Reached True in end of Recursive List')
                                        if joincount >0:
                                            end_of_tail_loop = False
                                            end_of_head_loop = False
                                            joincount = 0
                                            whole_loop_counter +=1
                                            print('TT>0 | active status R1: '+str(recursive_seq_1_active)+' | active status R2: '+str(recursive_seq_2_active))
                                            print('active data: '+str(active_data)+' | R1 data: '+str(recursive_seq1)+'\n'+' whole cycle count: '+str(whole_loop_counter)+' | join count '+ str(joincount))
                                            if recursive_seq_1_active == True: #This checks what data list is active and then I delete its entire contents
                                                del recursive_seq1[:] #the joined reads will be added to the inactive list which is empty
                                                print('/\/\/\/\/DELETION: OF R1 | there should be nothing between these brackets '+'['+str(recursive_seq1)+']')
                                                joincount = 0 #resets the join count number
                                                for i in removeread_list:
                                                    if i in recursive_seq2:
                                                        recursive_seq2.remove(i)
                                                final_assembly = recursive_seq2
                                                should_restart = True
                                            elif recursive_seq_2_active == True:
                                                del recursive_seq2[:]
                                                joincount = 0 #resets the join count number
                                                print('/\/\/\/\/DELETION: OF R2 | there should be nothing between these brackets '+'['+str(recursive_seq2)+']')
                                                for i in removeread_list:
                                                    if i in recursive_seq1:
                                                        recursive_seq1.remove(i)
                                                final_assembly = recursive_seq1
                                                should_restart = True
                                else:
                                    if head_i == len(active_data)-2 and tail_i ==len(active_data)-1 and len(head_seq)<=((len(head_seq_ref)/2)+1): #checks if the last list item of the first FOR loop has been iterated
                                        end_of_head_loop = True
                                        print('-----Reached True in end of Recursive List')
                                        if joincount >0:
                                            end_of_tail_loop = False
                                            end_of_head_loop = False
                                            joincount = 0
                                            whole_loop_counter +=1
                                            print('TT>0 | active status R1: '+str(recursive_seq_1_active)+' | active status R2: '+str(recursive_seq_2_active))
                                            print('active data: '+str(active_data)+' | R1 data: '+str(recursive_seq1)+'\n'+' whole cycle count: '+str(whole_loop_counter)+' | join count '+ str(joincount))
                                            if recursive_seq_1_active == True: #This checks what data list is active and then I delete its entire contents
                                                del recursive_seq1[:] #the joined reads will be added to the inactive list which is empty
                                                print('/\/\/\/\/DELETION: OF R1 | there should be nothing between these brackets '+'['+str(recursive_seq1)+']')
                                                joincount = 0 #resets the join count number
                                                for i in removeread_list:
                                                    if i in recursive_seq2:
                                                        recursive_seq2.remove(i)
                                                final_assembly = recursive_seq2
                                                should_restart = True
                                            elif recursive_seq_2_active == True:
                                                del recursive_seq2[:]
                                                joincount = 0 #resets the join count number
                                                print('/\/\/\/\/DELETION: OF R2 | there should be nothing between these brackets '+'['+str(recursive_seq2)+']')
                                                for i in removeread_list:
                                                    if i in recursive_seq1:
                                                        recursive_seq1.remove(i)
                                                final_assembly = recursive_seq1
                                                should_restart = True
                            elif recursive_seq_2_active == True:
                                if tail_seq_ref not in recursive_seq1 or head_seq_ref not in recursive_seq1:
                                    if head_seq_ref not in recursive_seq1:
                                        recursive_seq1.append(head_seq_ref)
                                        print('sequences transfering over to R1: '+head_seq_ref)
                                    if tail_seq_ref not in recursive_seq1:
                                        recursive_seq1.append(tail_seq_ref)
                                        print('sequences transfering over to R1: '+tail_seq_ref)
                                    elif head_i == len(active_data)-2 and tail_i ==len(active_data)-1 and len(head_seq)<=((len(head_seq_ref)/2)+1): #checks if the last list item of the first FOR loop has been iterated
                                        end_of_head_loop = True
                                        print('-----Reached True in end of Recursive List')
                                        if joincount >0:
                                            end_of_tail_loop = False
                                            end_of_head_loop = False
                                            joincount = 0
                                            whole_loop_counter +=1
                                            print('TT>0 | active status R1: '+str(recursive_seq_1_active)+' | active status R2: '+str(recursive_seq_2_active))
                                            print('active data: '+str(active_data)+' | R1 data: '+str(recursive_seq1)+'\n'+' whole cycle count: '+str(whole_loop_counter)+' | join count '+ str(joincount))
                                            if recursive_seq_1_active == True: #This checks what data list is active and then I delete its entire contents
                                                del recursive_seq1[:] #the joined reads will be added to the inactive list which is empty
                                                print('/\/\/\/\/DELETION: OF R1 | there should be nothing between these brackets '+'['+str(recursive_seq1)+']')
                                                joincount = 0 #resets the join count number
                                                for i in removeread_list:
                                                    if i in recursive_seq2:
                                                        recursive_seq2.remove(i)
                                                final_assembly = recursive_seq2
                                                should_restart = True
                                            elif recursive_seq_2_active == True:
                                                del recursive_seq2[:]
                                                joincount = 0 #resets the join count number
                                                print('/\/\/\/\/DELETION: OF R2 | there should be nothing between these brackets '+'['+str(recursive_seq2)+']')
                                                for i in removeread_list:
                                                    if i in recursive_seq1:
                                                        recursive_seq1.remove(i)
                                                final_assembly = recursive_seq1
                                                should_restart = True
                                else:
                                    if head_i == len(active_data)-2 and tail_i ==len(active_data)-1 and len(head_seq)<=((len(head_seq_ref)/2)+1): #checks if the last list item of the first FOR loop has been iterated
                                        end_of_head_loop = True
                                        print('-----Reached True in end of Recursive List')
                                        if joincount >0:
                                            end_of_tail_loop = False
                                            end_of_head_loop = False
                                            joincount = 0
                                            whole_loop_counter +=1
                                            print('TT>0 | active status R1: '+str(recursive_seq_1_active)+' | active status R2: '+str(recursive_seq_2_active))
                                            print('active data: '+str(active_data)+' | R1 data: '+str(recursive_seq1)+'\n'+' whole cycle count: '+str(whole_loop_counter)+' | join count '+ str(joincount))
                                            if recursive_seq_1_active == True: #This checks what data list is active and then I delete its entire contents
                                                del recursive_seq1[:] #the joined reads will be added to the inactive list which is empty
                                                print('/\/\/\/\/DELETION: OF R1 | there should be nothing between these brackets '+'['+str(recursive_seq1)+']')
                                                joincount = 0 #resets the join count number
                                                for i in removeread_list:
                                                    if i in recursive_seq2:
                                                        recursive_seq2.remove(i)
                                                final_assembly = recursive_seq2
                                                should_restart = True
                                            elif recursive_seq_2_active == True:
                                                del recursive_seq2[:]
                                                joincount = 0 #resets the join count number
                                                print('/\/\/\/\/DELETION: OF R2 | there should be nothing between these brackets '+'['+str(recursive_seq2)+']')
                                                for i in removeread_list:
                                                    if i in recursive_seq1:
                                                        recursive_seq1.remove(i)
                                                final_assembly = recursive_seq1
                                                should_restart = True
                        elif tail_seq == head_seq:
                            joincount +=1 #this tracks how many joins were made. as longer as the counter in line 73 is != 0 then the program will proceed
                            head_len = len(head_seq)
                            tail_len = len(tail_seq)
                            head_join = head_seq_ref[head_len:]
                            tail_join = tail_seq_ref[:-tail_len]
                            joinedread = tail_join +head_seq+ head_join
                            removeread_list.append(head_seq_ref)
                            removeread_list.append(tail_seq_ref)
                            if recursive_seq_1_active == True:
                                if joinedread not in recursive_seq2:
                                    recursive_seq2.append(joinedread)
                                    print('             RECURSIVE LIST JOIN(R1)!----->           '+joinedread)
                                    #if tail_i == len(active_data)-2: #checks if the last list item of the first FOR loop has been iterated
                                        #
                                        #print('-----Reached TRUE in end of TAIL loop')
                                    if head_i == len(active_data)-2 and tail_i ==len(active_data)-1 and len(head_seq)<=((len(head_seq_ref)/2)+1): #checks if the last list item of the first FOR loop has been iterated
                                        end_of_head_loop = True
                                        print('-----Reached True in end of Recursive List')
                                        if joincount >0:
                                            end_of_tail_loop = False
                                            end_of_head_loop = False
                                            joincount = 0
                                            whole_loop_counter +=1
                                            print('TT>0 | active status R1: '+str(recursive_seq_1_active)+' | active status R2: '+str(recursive_seq_2_active))
                                            print('active data: '+str(active_data)+' | R1 data: '+str(recursive_seq1)+'\n'+' whole cycle count: '+str(whole_loop_counter)+' | join count '+ str(joincount))
                                            if recursive_seq_1_active == True: #This checks what data list is active and then I delete its entire contents
                                                del recursive_seq1[:] #the joined reads will be added to the inactive list which is empty
                                                print('/\/\/\/\/DELETION: OF R1 | there should be nothing between these brackets '+'['+str(recursive_seq1)+']')
                                                joincount = 0 #resets the join count number
                                                for i in removeread_list:
                                                    if i in recursive_seq2:
                                                        recursive_seq2.remove(i)
                                                final_assembly = recursive_seq2
                                                should_restart = True
                                            elif recursive_seq_2_active == True:
                                                del recursive_seq2[:]
                                                joincount = 0 #resets the join count number
                                                print('/\/\/\/\/DELETION: OF R2 | there should be nothing between these brackets '+'['+str(recursive_seq2)+']')
                                                for i in removeread_list:
                                                    if i in recursive_seq1:
                                                        recursive_seq1.remove(i)
                                                final_assembly = recursive_seq1
                                                should_restart = True

                            elif recursive_seq_2_active == True:
                                if joinedread not in recursive_seq1:
                                    recursive_seq1.append(joinedread)
                                    print('             RECURSIVE LIST JOIN(R2)!----->           '+joinedread)
                                    #if tail_i == len(active_data)-2: #checks if the last list item of the first FOR loop has been iterated
                                        #
                                        #print('-----Reached TRUE in end of TAIL loop')
                                    if head_i == len(active_data)-2 and tail_i ==len(active_data)-1 and len(head_seq)<=((len(head_seq_ref)/2)+1): #checks if the last list item of the first FOR loop has been iterated
                                        end_of_head_loop = True
                                        print('-----Reached True in end of Recursive List')
                                        if joincount >0:
                                            end_of_tail_loop = False
                                            end_of_head_loop = False
                                            joincount = 0
                                            whole_loop_counter +=1
                                            print('TT>0 | active status R1: '+str(recursive_seq_1_active)+' | active status R2: '+str(recursive_seq_2_active))
                                            print('active data: '+str(active_data)+' | R1 data: '+str(recursive_seq1)+'\n'+' whole cycle count: '+str(whole_loop_counter)+' | join count '+ str(joincount))
                                            if recursive_seq_1_active == True: #This checks what data list is active and then I delete its entire contents
                                                del recursive_seq1[:] #the joined reads will be added to the inactive list which is empty
                                                print('/\/\/\/\/DELETION: OF R1 | there should be nothing between these brackets '+'['+str(recursive_seq1)+']')
                                                joincount = 0 #resets the join count number
                                                final_assembly = recursive_seq2
                                                should_restart = True
                                            elif recursive_seq_2_active == True:
                                                del recursive_seq2[:]
                                                joincount = 0 #resets the join count number
                                                print('/\/\/\/\/DELETION: OF R2 | there should be nothing between these brackets '+'['+str(recursive_seq2)+']')
                                                final_assembly = recursive_seq1
                                                should_restart = True

#look at the infinite loops and see if it is making any joins at all after each whole cycle. maybe there is a loop hole

if len(recursive_seq1) > 0:
    print('R1: '+str(recursive_seq1))
if len(recursive_seq2) > 0:
    print('R2: '+str(recursive_seq2))
print('Final assembly ('+str(whole_loop_counter)+') : '+str(final_assembly))
print(answercheck)

#Note to self: fix the cycling issue. last tiem i figured out that the last item in the list was not being iterated because
