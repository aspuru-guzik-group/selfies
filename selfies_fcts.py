import math

def make_brackets_around_atoms(smiles):
    ii=0
    smiles=smiles.replace(' ','')
    current_new_smiles=''
    while ii<len(smiles):
        if smiles[ii]=='[':
            small_smiles=smiles[ii:]
            closing_bracket=small_smiles.find(']')
            in_bracket_content=small_smiles[0:closing_bracket+1]
            current_new_smiles+=in_bracket_content
            ii+=len(in_bracket_content)
            
        elif (smiles[ii]>='A' and smiles[ii]<='Z') or (smiles[ii]>='a' and smiles[ii]<='z') or smiles[ii]=='*':
            smiles=smiles+' '
            if smiles[ii:ii+2]=='Br' or smiles[ii:ii+2]=='Cl':
                current_new_smiles+='['+smiles[ii:ii+2]+']'
                ii+=2
            else:
                current_new_smiles+='['+smiles[ii:ii+1]+']'   
                ii+=1
            smiles=smiles.replace(' ','')
        elif smiles[ii]=='(' or smiles[ii]==')' or smiles[ii]=='.' or smiles[ii]=='%' or (smiles[ii]>='0' and smiles[ii]<='9'):
            current_new_smiles+=smiles[ii]
            ii+=1
        elif smiles[ii]=='=' or smiles[ii]=='#' or smiles[ii]=='\\' or smiles[ii]=='/' or smiles[ii]=='-':
            pre_symbol=smiles[ii]
            ii+=1
            if smiles[ii]=='[':
                small_smiles=smiles[ii:]
                closing_bracket=small_smiles.find(']')
                in_bracket_content=small_smiles[1:closing_bracket]
                current_new_smiles+='['+pre_symbol+in_bracket_content+']'
                ii+=len(in_bracket_content)+2
            elif (smiles[ii]>='A' and smiles[ii]<='Z') or (smiles[ii]>='a' and smiles[ii]<='z')  or smiles[ii]=='*':
                smiles=smiles+' '
                if smiles[ii:ii+2]=='Br' or smiles[ii:ii+2]=='Cl':
                    current_new_smiles+='['+pre_symbol+smiles[ii:ii+2]+']'
                    ii+=2
                else:
                    current_new_smiles+='['+pre_symbol+smiles[ii:ii+1]+']'   
                    ii+=1
                smiles=smiles.replace(' ','')
            elif smiles[ii]=='%' or (smiles[ii]>='0' and smiles[ii]<='9'): # explicit bond-number with ring                
                current_new_smiles+=pre_symbol+smiles[ii]
                ii+=1
            else:
                raise ValueError('make_brackets_around_atoms: Unknown symbol in the string.')

        else:
            raise ValueError('make_brackets_around_atoms: Unknown symbol in the string.')
    return current_new_smiles




def reconfigure_smiles_numbers1(smiles):
    available_nums=['@$aa','@$ab','@$ac','@$ad','@$ae','@$af','@$ag','@$ah','@$ai','@$aj','@$ak','@$al','@$am','@$an','@$ao','@$ap','@$aq','@$ar','@$as','@$at','@$au','@$av','@$aw','@$ax','@$ay','@$az','@$ba','@$bb','@$bc','@$bd','@$be','@$bf','@$bg','@$bh','@$bi','@$bj','@$bk','@$bl','@$bm','@$bn','@$bo','@$bp','@$bq','@$br','@$bs','@$bt','@$bu','@$bv','@$bw','@$bx','@$by','@$bz','@$ca','@$cb','@$cc','@$cd','@$ce','@$cf','@$cg','@$ch','@$ci','@$cj','@$ck','@$cl','@$cm','@$cn','@$co','@$cp','@$cq','@$cr','@$cs','@$ct','@$cu','@$cv','@$cw','@$cx','@$cy','@$cz','@$da','@$db','@$dc','@$dd','@$de','@$df','@$dg','@$dh','@$di','@$dj','@$dk','@$dl','@$dm','@$dn','@$do','@$dp','@$dq','@$dr','@$ds','@$dt','@$du','@$dv']
    all_smiles=smiles.split('.')
    all_smiles_new=''

    for current_smiles in all_smiles:
        current_smiles_empty=current_smiles
        in_bracket=0
        jj=0
        while jj<len(current_smiles_empty):
            if current_smiles_empty[jj]=='[':
                in_bracket=1
            elif current_smiles_empty[jj]==']':
                in_bracket=0
            else:
                if in_bracket==1:
                    current_smiles_empty=current_smiles_empty[0:jj]+' '+current_smiles_empty[jj+1:]
                    
            if current_smiles_empty[jj]>='0' and current_smiles_empty[jj]<='9':
                current_smiles_empty=current_smiles_empty[0:jj]+'%0'+current_smiles_empty[jj:]
                current_smiles=current_smiles[0:jj]+'%0'+current_smiles[jj:]
                jj+=2
                
            if current_smiles_empty[jj]=='%':
                jj+=2
            jj+=1
    
        an_idx=0
        cc=0
        loop_count=0
        while cc<len(current_smiles_empty):
            if cc<0:
                raise ValueError('reconfigure_smiles_numbers1: Malformed ring.')
            letter=current_smiles_empty[cc]            
            if letter=='%':
                current_num=current_smiles_empty[cc:cc+3]
                for jj in range(2):
                    pos_of_num=current_smiles_empty.find(current_num)
                    current_smiles=current_smiles[0:pos_of_num]+available_nums[an_idx]+current_smiles[pos_of_num+3:]
                    current_smiles_empty=current_smiles_empty[0:pos_of_num]+available_nums[an_idx]+current_smiles_empty[pos_of_num+3:]
                    cc-=2
                an_idx+=1
            cc+=1
            loop_count+=1
            if loop_count>10000:
                raise ValueError('reconfigure_smiles_numbers1: Malformed ring.')
        all_smiles_new=all_smiles_new+'.'+current_smiles
        
    return all_smiles_new[1:]



def reconfigure_smiles_numbers2(smiles):
    available_nums=['@$aa','@$ab','@$ac','@$ad','@$ae','@$af','@$ag','@$ah','@$ai','@$aj','@$ak','@$al','@$am','@$an','@$ao','@$ap','@$aq','@$ar','@$as','@$at','@$au','@$av','@$aw','@$ax','@$ay','@$az','@$ba','@$bb','@$bc','@$bd','@$be','@$bf','@$bg','@$bh','@$bi','@$bj','@$bk','@$bl','@$bm','@$bn','@$bo','@$bp','@$bq','@$br','@$bs','@$bt','@$bu','@$bv','@$bw','@$bx','@$by','@$bz','@$ca','@$cb','@$cc','@$cd','@$ce','@$cf','@$cg','@$ch','@$ci','@$cj','@$ck','@$cl','@$cm','@$cn','@$co','@$cp','@$cq','@$cr','@$cs','@$ct','@$cu','@$cv','@$cw','@$cx','@$cy','@$cz','@$da','@$db','@$dc','@$dd','@$de','@$df','@$dg','@$dh','@$di','@$dj','@$dk','@$dl','@$dm','@$dn','@$do','@$dp','@$dq','@$dr','@$ds','@$dt','@$du','@$dv']
    all_smiles=smiles.split('.')
    all_smiles_new=''

    for current_smiles in all_smiles:
        
        for jj in range(len(available_nums)):
            find_num_1=current_smiles.find(available_nums[jj])            
            
            if find_num_1>=0:
                tmp_current_smiles=current_smiles[0:find_num_1]+'X'+current_smiles[find_num_1+1:]
                find_num_2=tmp_current_smiles.find(available_nums[jj])  
                sub_smiles=current_smiles[find_num_1+4:find_num_2]
                ring_size=len(sub_smiles)-len(sub_smiles.replace('[',''))
                str_ring_size=str(ring_size)
                
                if len(str_ring_size)==1:
                    ring_sizeSymbol='%00'+str_ring_size
                elif len(str_ring_size)==2:
                    ring_sizeSymbol='%0'+str_ring_size
                elif len(str_ring_size)==3:
                    ring_sizeSymbol='%'+str_ring_size
                else:
                    raise ValueError('reconfigure_smiles_numbers2: Very long ring is not implemented.')
                    
                current_smiles=current_smiles[0:find_num_1]+current_smiles[find_num_1+4:find_num_2]+ring_sizeSymbol+current_smiles[find_num_2+4:]
            else:
                break
            
        all_smiles_new=all_smiles_new+'.'+current_smiles
        
    return all_smiles_new[1:]
  
    
def cumsum(int_list,cum_offset=0): # without numpy
    cum_list=[cum_offset]
    for x in int_list:
        cum_list.append(cum_list[-1]+x)
    return cum_list[1:]
    

def smiles_to_selfies(smiles):
    start_alphabet=['[epsilon]','[Ring1]','[Ring2]','[Branch1_1]','[Branch1_2]','[Branch1_3]','[Branch2_1]','[Branch2_2]','[Branch2_3]','[F]','[O]','[=O]','[N]','[=N]','[#N]','[C]','[=C]','[#C]','[S]','[=S]'];
    all_smiles=smiles.split('.')
    all_smiles_new=''

    for current_smiles in all_smiles:
        transitions=''
        tmp_smiles=current_smiles+'      '
        while True:
            current_symbol=tmp_smiles[0]
            tmp_smiles=tmp_smiles[1:]

            if current_symbol==' ':
                break

            if current_symbol=='[':
                pos_of_close=tmp_smiles.find(']')
                current_symbol=tmp_smiles[0:pos_of_close] # without the brackets
                tmp_smiles=tmp_smiles[pos_of_close+1:]
                transitions=transitions+'['+current_symbol+']'

            elif current_symbol=='%':
                next_symbol=tmp_smiles[0:3]
                tmp_smiles=tmp_smiles[3:]

                current_num=int(next_symbol)
                if current_num>=2 and current_num<=400:
                    if current_num<=20:
                        ring_symbol=start_alphabet[current_num-2]
                        transitions=transitions+'[Ring1]'+ring_symbol 
                    else:
                        ring_num_1=int(math.floor(current_num/20.))
                        ring_num_2=current_num%20
                        ring_symbol1=start_alphabet[ring_num_1-1]
                        ring_symbol2=start_alphabet[ring_num_2]

                        transitions=transitions+'[Ring2]'+ring_symbol1+ring_symbol2

                else:
                    if current_num<2:
                        raise ValueError('smiles_to_selfies: Malformed Ring.')                    
                    else:
                        raise ValueError('smiles_to_selfies: Very large ring is not implemented.')

            elif (current_symbol=='=' or current_symbol=='#' or current_symbol=='\\' or current_symbol=='/' or current_symbol=='-') and tmp_smiles[0]=='%':
                # explicit bond number for ring                
                pre_symbol=current_symbol
                next_symbol=tmp_smiles[1:4]
                tmp_smiles=tmp_smiles[4:]
                current_num=int(next_symbol)
                if current_num>=2 and current_num<=400:
                    if current_num<=20:
                        ring_symbol=start_alphabet[current_num-2]
                        transitions=transitions+'[Expl'+pre_symbol+'Ring1]'+ring_symbol 
                    else:
                        ring_num_1=int(math.floor(current_num/20.))
                        ring_num_2=current_num%20
                        ring_symbol1=start_alphabet[ring_num_1-1]
                        ring_symbol2=start_alphabet[ring_num_2]

                        transitions=transitions+'[Expl'+pre_symbol+'Ring2]'+ring_symbol1+ring_symbol2

                else:
                    if current_num<2:
                        raise ValueError('smiles_to_selfies: Malformed Ring.')                    
                    else:
                        raise ValueError('smiles_to_selfies: Very large ring is not implemented.')

            
            elif current_symbol=='(':
                open_vec=[0]*len(tmp_smiles)
                close_vec=[0]*len(tmp_smiles)
                tmp_smiles2=tmp_smiles
                while True:
                    open_bracket=tmp_smiles2.find('(')
                    if open_bracket>=0:
                        tmp_smiles2=tmp_smiles2[0:open_bracket]+' '+tmp_smiles2[open_bracket+1:]
                        open_vec[open_bracket]=1
                    else:
                        break
                
                while True:
                    close_bracket=tmp_smiles2.find(')')
                    if close_bracket>=0:
                        tmp_smiles2=tmp_smiles2[0:close_bracket]+' '+tmp_smiles2[close_bracket+1:]
                        close_vec[close_bracket]=1
                    else:
                        break

                difference_of_list=[x - y for x, y in zip(open_vec, close_vec)]                
                layer_of_brackets=cumsum(difference_of_list,1)
                
                if 0 in layer_of_brackets:
                    sub_smiles=tmp_smiles[0:layer_of_brackets.index(0)]
                    tmp_smiles=tmp_smiles[layer_of_brackets.index(0)+1:]
                    
                else:
                    sub_smiles=tmp_smiles
                    tmp_smiles=[]
                    
                sub_ds=smiles_to_selfies(sub_smiles)

                current_num=len(sub_ds)-len(sub_ds.replace('[',''))

                if current_num<=20:
                    ring_symbol=start_alphabet[current_num-1]
                    if tmp_smiles[1]=='=' or tmp_smiles[1]=='#' or tmp_smiles[0]=='(':
                        if sub_smiles[1]=='=' or sub_smiles[1]=='#':
                            transitions=transitions+'[Branch1_1]'+ring_symbol+sub_ds
                        else:
                            transitions=transitions+'[Branch1_2]'+ring_symbol+sub_ds
                    else:
                        transitions=transitions+'[Branch1_3]'+ring_symbol+sub_ds
                
                elif current_num<=400:
                    ring_num_1=int(math.floor(current_num/20.))
                    ring_num_2=current_num%20
                    ring_symbol1=start_alphabet[ring_num_1-1]
                    ring_symbol2=start_alphabet[ring_num_2]
                    
                    if len(tmp_smiles)<=2:
                        transitions=transitions+'[Branch2_3]'+ring_symbol1+ring_symbol2+sub_ds
                    else:
                        if tmp_smiles[1]=='=' or tmp_smiles[1]=='#' or tmp_smiles[0]=='(':
                            if sub_smiles[1]=='=' or sub_smiles[1]=='#':
                                transitions=transitions+'[Branch2_1]'+ring_symbol1+ring_symbol2+sub_ds
                            else:
                                transitions=transitions+'[Branch2_2]'+ring_symbol1+ring_symbol2+sub_ds
                        else:
                            transitions=transitions+'[Branch2_3]'+ring_symbol1+ring_symbol2+sub_ds
                            
                elif current_num<=8000:
                    ring_num_1=int(math.floor(current_num/400.))
                    current_num1=current_num-ring_num_1*400
                    ring_num_2=int(math.floor(current_num1/20.))
                    ring_num_3=current_num1-ring_num_2*20
                    ring_symbol1=start_alphabet[ring_num_1-1]
                    ring_symbol2=start_alphabet[ring_num_2]
                    ring_symbol3=start_alphabet[ring_num_3]

                    if len(tmp_smiles)<=2:
                        transitions=transitions+'[Branch3_3]'+ring_symbol1+ring_symbol2+ring_symbol3+sub_ds
                    else:
                        if tmp_smiles[1]=='=' or tmp_smiles[1]=='#' or tmp_smiles[0]=='(':
                            if sub_smiles[1]=='=' or sub_smiles[1]=='#':
                                transitions=transitions+'[Branch3_1]'+ring_symbol1+ring_symbol2+ring_symbol3+sub_ds
                            else:
                                transitions=transitions+'[Branch3_2]'+ring_symbol1+ring_symbol2+ring_symbol3+sub_ds
                        else:
                            transitions=transitions+'[Branch3_3]'+ring_symbol1+ring_symbol2+ring_symbol3+sub_ds                
                else:
                    raise ValueError('smiles_to_selfies: Very large branch is not implemented (current_num='+str(current_num)+').')
            
            else:
                raise ValueError('smiles_to_selfies: Unknown Symbol')

        all_smiles_new=all_smiles_new+'.'+transitions
        
    return all_smiles_new[1:]
      


def get_next_selfies_symbol(tmp_ds):
    next_symbol=''
    tmp_ds_new=tmp_ds
    if len(tmp_ds)<=2:
        return [next_symbol, tmp_ds_new]
    
    if tmp_ds[0]!='[':
        raise ValueError('get_next_selfies_symbol: Decoding Problem 1: '+tmp_ds)
    
    end_of_symbol=tmp_ds.find(']')
    
    if end_of_symbol==-1:
        raise ValueError('get_next_selfies_symbol: Decoding Problem 2: '+tmp_ds)
    else:
        next_symbol=tmp_ds[0:end_of_symbol+1]
        tmp_ds_new=tmp_ds_new[end_of_symbol+1:]
        
    return [next_symbol, tmp_ds_new]
    
    


def selfies_to_smiles_derive(selfies,smiles):    
    start_alphabet=['[epsilon]','[Ring1]','[Ring2]','[Branch1_1]','[Branch1_2]','[Branch1_3]','[Branch2_1]','[Branch2_2]','[Branch2_3]','[F]','[O]','[=O]','[N]','[=N]','[#N]','[C]','[=C]','[#C]','[S]','[=S]'];
            
    tmp_ds=selfies.replace('X','Z!') # X will be used as states of the derivation
    
    next_X=smiles.find('X');
    while next_X>=0:
        state=int(smiles[next_X+1])
        before_smiles=smiles[0:next_X]
        after_smiles=smiles[next_X+2:]
                
        [current_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
        
        if state==0:
            if current_symbol=='[epsilon]':
                new_smiles_symbol='X0'
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0:    
                new_smiles_symbol='X0'
                [_,tmp_ds]=get_next_selfies_symbol(tmp_ds)  # ignore next symbol                 
            elif current_symbol=='[F]':
                new_smiles_symbol='[F]X1'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]X1'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]X1'
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X2'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[O]X2'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X4'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[C]X4'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[C]X4'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'
            smiles=before_smiles+new_smiles_symbol+after_smiles

        if state==1:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''
            elif current_symbol.find('Branch')>=0:   
                new_smiles_symbol='X0'
                [_,tmp_ds]=get_next_selfies_symbol(tmp_ds)  # ignore next symbol               
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='2'
                
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')
            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='2'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[F]':
                new_smiles_symbol='[F]'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]'                
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X1'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[O]'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'                
            smiles=before_smiles+new_smiles_symbol+after_smiles

        if state==2:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''              
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                     
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='2'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X1'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X1'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X1'
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                      
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='2'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X1'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X1'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X1'
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]' or current_symbol=='[Branch1_2]' or current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'

            elif current_symbol=='[Branch2_1]' or current_symbol=='[Branch2_2]' or current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                

            elif current_symbol=='[Branch3_1]' or current_symbol=='[Branch3_2]' or current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol

                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'              
                
            elif current_symbol=='[F]':
                new_smiles_symbol='[F]'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]'                
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X1'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[=O]'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[=N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[=N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'
            smiles=before_smiles+new_smiles_symbol+after_smiles


        if state==3:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''             
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                  
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='2'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X2'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X2'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X2'
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                    
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='2'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X2'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X2'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X2'
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]' or current_symbol=='[Branch1_2]':
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
            elif current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X6')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'             
                
            elif current_symbol=='[Branch2_1]' or current_symbol=='[Branch2_2]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
                
            elif current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X6')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                
                

            elif current_symbol=='[Branch3_1]' or current_symbol=='[Branch3_2]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
                
            elif current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X6')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'                 
                
                
            elif current_symbol=='[F]':
                new_smiles_symbol='[F]'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]'                
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X1'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[=O]'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[=N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[#N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'

            smiles=before_smiles+new_smiles_symbol+after_smiles


        if state==4:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''            
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='2'                

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X4'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X4'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X4'
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                new_smiles_symbol='X4'
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='2'

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X4'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X4'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X4'
                else:
                    raise ValueError('selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]':
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X6')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'
                
            elif current_symbol=='[Branch1_2]':
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'              
                
            elif current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X7')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                
            elif current_symbol=='[Branch2_1]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X6')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'
                

            elif current_symbol=='[Branch2_2]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'            
                
                
            elif current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)   
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X7')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                
                
                

            elif current_symbol=='[Branch3_1]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X6')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'
                

            elif current_symbol=='[Branch3_2]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X5')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'         
                
                
            elif current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=selfies_to_smiles_derive(branch_str,'X7')
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'         
                
                
                
            elif current_symbol=='[F]':
                new_smiles_symbol='[F]'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]'                  
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X1'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[=O]'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[=N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[#N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'

            smiles=before_smiles+new_smiles_symbol+after_smiles

        if state==5:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0: 
                new_smiles_symbol='X5'
            elif current_symbol=='[F]':
                new_smiles_symbol='[F]'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]'                  
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X1'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[O]'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'
            smiles=before_smiles+new_smiles_symbol+after_smiles



        if state==6:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''          
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0: 
                new_smiles_symbol='X6' 
            elif current_symbol=='[F]':
                new_smiles_symbol='[F]'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]'                  
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X1'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[=O]'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[=N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[=N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'
            smiles=before_smiles+new_smiles_symbol+after_smiles




        if state==7:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''        
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0: 
                new_smiles_symbol='X7'     
            elif current_symbol=='[F]':
                new_smiles_symbol='[F]'
            elif current_symbol=='[Cl]':
                new_smiles_symbol='[Cl]'
            elif current_symbol=='[Br]':
                new_smiles_symbol='[Br]'                  
            elif current_symbol=='[O]':
                new_smiles_symbol='[O]X1'
            elif current_symbol=='[=O]':
                new_smiles_symbol='[=O]'
            elif current_symbol=='[N]':
                new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                new_smiles_symbol='[=N]X4'
            elif current_symbol=='[#N]':
                new_smiles_symbol='[#N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X4'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X4'

            smiles=before_smiles+new_smiles_symbol+after_smiles
    
        if len(tmp_ds)<=2:
            while True:
                non_terminal=smiles.find('X')
                if non_terminal>=0:
                    smiles=smiles[0:non_terminal]+smiles[non_terminal+2:]
                else:
                    break;

        next_X=smiles.find('X')
        

    return smiles.replace('Z!','X')



def selfies_to_smiles(selfies):
    all_selfies=selfies.split('.')
    all_selfies_new=''

    for current_smiles in all_selfies:
        all_selfies_new=all_selfies_new+'.'+selfies_to_smiles_derive(current_smiles,'X0')
            
    return all_selfies_new[1:] 




def insert_rings_to_smiles(smiles):
    # For simplicity, rings are inserted after the full derivation. However, this does not change anything in the concept of the Grammar.
    
    available_nums=['@$aa','@$ab','@$ac','@$ad','@$ae','@$af','@$ag','@$ah','@$ai','@$aj','@$ak','@$al','@$am','@$an','@$ao','@$ap','@$aq','@$ar','@$as','@$at','@$au','@$av','@$aw','@$ax','@$ay','@$az','@$ba','@$bb','@$bc','@$bd','@$be','@$bf','@$bg','@$bh','@$bi','@$bj','@$bk','@$bl','@$bm','@$bn','@$bo','@$bp','@$bq','@$br','@$bs','@$bt','@$bu','@$bv','@$bw','@$bx','@$by','@$bz','@$ca','@$cb','@$cc','@$cd','@$ce','@$cf','@$cg','@$ch','@$ci','@$cj','@$ck','@$cl','@$cm','@$cn','@$co','@$cp','@$cq','@$cr','@$cs','@$ct','@$cu','@$cv','@$cw','@$cx','@$cy','@$cz','@$da','@$db','@$dc','@$dd','@$de','@$df','@$dg','@$dh','@$di','@$dj','@$dk','@$dl','@$dm','@$dn','@$do','@$dp','@$dq','@$dr','@$ds','@$dt','@$du','@$dv']
    unique_nums=[]
    for ii in range(200):
        num_str=str(ii+1)
        if len(num_str)==1:
            num_str='%00'+num_str
        elif len(num_str)==2:
            num_str='%0'+num_str   
        elif len(num_str)==3:
            num_str='%'+num_str 
        unique_nums.append(num_str)
        
    smiles=smiles.replace('%','_') # to avoid conflicts when adding numbers
    
    ring_counter=1
    while True:        
        pos_of_ring_symbol=smiles.find('_')
        if pos_of_ring_symbol>=0:
            if (len(smiles[0:pos_of_ring_symbol])-len(smiles[0:pos_of_ring_symbol].replace('[','')))>2:
                    
                size_of_ring=int(smiles[pos_of_ring_symbol+1:pos_of_ring_symbol+4])+1
                
                num_of_elements=len(smiles[0:pos_of_ring_symbol])-len(smiles[0:pos_of_ring_symbol].replace(']',''))            
                r_count=pos_of_ring_symbol
                s_count=min(num_of_elements,size_of_ring)
                            
                while True:
                    if s_count==0: # found the target atom
                        break
                    
                    r_count-=1
                    if smiles[r_count]==']':
                        s_count-=1
    
                ring_symbol=available_nums[ring_counter]
    
                if smiles[r_count+1]=='@': # There is already an identifier, i'm the subsequent one
                    r_count+=4
                    if smiles[r_count+1]=='@':
                        r_count+=4
                
                smiles=smiles[0:r_count+1]+ring_symbol+smiles[r_count+1:pos_of_ring_symbol]+ring_symbol+smiles[pos_of_ring_symbol+4:]
                ring_counter+=1                
            
            else:    # this can be seen as two additional states of the derivation.  X0 -> X0A -> X0B, and in X0* rings are not allowed. 
                smiles=smiles[0:pos_of_ring_symbol]+smiles[pos_of_ring_symbol+4:]
            
        else:
            break
        
    cc=1
    while True:
        if cc>(len(smiles)-2):
            break
        
        if smiles[cc:cc+2]=='@$': # replace by unique identifier
            ring_id=smiles[cc:cc+4]            
            smiles=smiles.replace(ring_id,unique_nums[0])            
            unique_nums=unique_nums[1:]
        
        elif smiles[cc:cc+2]=='%0': # this identifier is free again for further use
            uniquering_id=smiles[cc:cc+4]
            unique_nums=[uniquering_id]+unique_nums
            unique_nums.sort()   # sort the list such that smallest identifier will be used next
            
        cc+=1
        
    smiles=smiles.replace('[=','=[')
    smiles=smiles.replace('[#','#[')
    smiles=smiles.replace('[\\','\\[')
    smiles=smiles.replace('[/','/[')
    smiles=smiles.replace('[-','-[')
            
    smiles=smiles.replace('[B]','B')
    smiles=smiles.replace('[I]','I')
    smiles=smiles.replace('[P]','P')
    smiles=smiles.replace('[C]','C')
    smiles=smiles.replace('[N]','N')
    smiles=smiles.replace('[O]','O')  
    smiles=smiles.replace('[F]','F')     
    smiles=smiles.replace('[S]','S')  
    smiles=smiles.replace('[Cl]','Cl')
    smiles=smiles.replace('[Br]','Br')
    smiles=smiles.replace('[c]','c')  
    smiles=smiles.replace('[n]','n')     
    smiles=smiles.replace('[o]','o')    
    smiles=smiles.replace('[s]','s')  
    smiles=smiles.replace('[p]','p')      
    
    
    smiles=smiles.replace('%00','')
    smiles=smiles.replace('%0','%')

    return(smiles)