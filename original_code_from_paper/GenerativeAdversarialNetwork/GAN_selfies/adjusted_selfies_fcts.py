import math



def _make_brackets_around_atoms(smiles): # first function in the encoder: All atoms are itemized via brackets. 
                                         # for example: C1=CO1O -> [C]1[=C][O][O]. 
                                         # Meaning that bond information is treated as part of the atom, to ensure semantical validity
                                         # Brackets define an element of the alphabet, and for each of them there is a rule vector in the grammar.
                                        
                                        
    ii=0
    smiles=smiles.replace(' ','')
    current_new_smiles=''
    while ii<len(smiles):
        if smiles[ii]=='[':
            small_smiles=smiles[ii:]
            closing_bracket=small_smiles.find(']')
            in_bracket_content=small_smiles[0:closing_bracket+1]
            current_new_smiles+=in_bracket_content[0:-1]+'expl]'
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
        elif smiles[ii]=='(' or smiles[ii]==')' or smiles[ii]=='.' or smiles[ii]=='%' or (smiles[ii]>='0' and smiles[ii]<='9'): # these symbols, for rings, brackets and divisions between molecules, will be outside the brackets
            current_new_smiles+=smiles[ii]
            ii+=1
        elif smiles[ii]=='=' or smiles[ii]=='#' or smiles[ii]=='\\' or smiles[ii]=='/' or smiles[ii]=='-': # these symbols will be part of the brackets. The sterochemical information could potentially be outside the bracket, by introducing a new rule vector for them. What is more advantageos need to be investigated later.
            pre_symbol=smiles[ii]
            ii+=1
            if smiles[ii]=='[':
                small_smiles=smiles[ii:]
                closing_bracket=small_smiles.find(']')
                in_bracket_content=small_smiles[1:closing_bracket]
                current_new_smiles+='['+pre_symbol+in_bracket_content+'expl]'
                ii+=len(in_bracket_content)+2
            elif (smiles[ii]>='A' and smiles[ii]<='Z') or (smiles[ii]>='a' and smiles[ii]<='z')  or smiles[ii]=='*': # upper and lower case letters, as well as the wildcard symbol is part of the bracket.
                smiles=smiles+' '
                if smiles[ii:ii+2]=='Br' or smiles[ii:ii+2]=='Cl':      # Chlor and Brom are usually written without brackets in SMILES strings, so need to be treated differently
                    current_new_smiles+='['+pre_symbol+smiles[ii:ii+2]+']'
                    ii+=2
                else:
                    current_new_smiles+='['+pre_symbol+smiles[ii:ii+1]+']'   # Add brackets around element
                    ii+=1
                smiles=smiles.replace(' ','')
            elif smiles[ii]=='%' or (smiles[ii]>='0' and smiles[ii]<='9'): # explicit bond-number with ring                
                current_new_smiles+=pre_symbol+smiles[ii]
                ii+=1
            else:
                raise ValueError('_make_brackets_around_atoms: Unknown symbol in the string.')

        else:
            raise ValueError('_make_brackets_around_atoms: Unknown symbol in the string.')
    return current_new_smiles




def _reconfigure_smiles_numbers1(smiles): # All rings get a unique identifiers
    # list of unique identifiers, for purpose of speed, this list is hardcoded
    available_nums=['@$aa','@$ab','@$ac','@$ad','@$ae','@$af','@$ag','@$ah','@$ai','@$aj','@$ak','@$al','@$am','@$an','@$ao','@$ap','@$aq','@$ar','@$as','@$at','@$au','@$av','@$aw','@$ax','@$ay','@$az','@$ba','@$bb','@$bc','@$bd','@$be','@$bf','@$bg','@$bh','@$bi','@$bj','@$bk','@$bl','@$bm','@$bn','@$bo','@$bp','@$bq','@$br','@$bs','@$bt','@$bu','@$bv','@$bw','@$bx','@$by','@$bz','@$ca','@$cb','@$cc','@$cd','@$ce','@$cf','@$cg','@$ch','@$ci','@$cj','@$ck','@$cl','@$cm','@$cn','@$co','@$cp','@$cq','@$cr','@$cs','@$ct','@$cu','@$cv','@$cw','@$cx','@$cy','@$cz','@$da','@$db','@$dc','@$dd','@$de','@$df','@$dg','@$dh','@$di','@$dj','@$dk','@$dl','@$dm','@$dn','@$do','@$dp','@$dq','@$dr','@$ds','@$dt','@$du','@$dv']
    all_smiles=smiles.split('.')
    all_smiles_new=''

    for current_smiles in all_smiles: # make all rings of the form '%NN'
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
        while cc<len(current_smiles_empty): # replace all rings with a unique identifier
            if cc<0:
                raise ValueError('_reconfigure_smiles_numbers1: Malformed ring.')
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
                raise ValueError('_reconfigure_smiles_numbers1: Malformed ring.')

        # bugfix for non canonical smiles - 26.08.2019
        for symbol in available_nums:
            find_identifier1=current_smiles.find(symbol)
            if find_identifier1==-1:
                break
            else:
                find_identifier1_single=current_smiles.find('-'+symbol)
                if find_identifier1_single==find_identifier1-1:
                    current_smiles_tmp=current_smiles[0:find_identifier1]+'x'+current_smiles[find_identifier1+1:]
                    
                    find_identifier2=current_smiles_tmp.find(symbol)
                    if current_smiles_tmp[find_identifier2-1]!='-':
                        all_smiles_new=current_smiles[0:find_identifier2]+'-'+current_smiles[find_identifier2:] # add explicit single bond from second identifier
                    current_smiles=current_smiles[0:find_identifier1_single]+current_smiles[find_identifier1:] # remove explicit single bond from first identifier
                    

                find_identifier1_single=current_smiles.find('='+symbol)
                if find_identifier1_single==find_identifier1-1:
                    current_smiles_tmp=current_smiles[0:find_identifier1]+'x'+current_smiles[find_identifier1+1:]
                    
                    find_identifier2=current_smiles_tmp.find(symbol)
                    if current_smiles_tmp[find_identifier2-1]!='=':
                        all_smiles_new=current_smiles[0:find_identifier2]+'-'+current_smiles[find_identifier2:] # add explicit single bond from second identifier
                    current_smiles=current_smiles[0:find_identifier1_single]+current_smiles[find_identifier1:] # remove explicit single bond from first identifier                    
        # ##
        
        
        
        all_smiles_new=all_smiles_new+'.'+current_smiles

    return all_smiles_new[1:]



def _reconfigure_smiles_numbers2(smiles): # pairs of unique identifiers will be transformed into symbols which referes to the relative distance between them. all of them are of the form '%NNN'
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
                    raise ValueError('_reconfigure_smiles_numbers2: Very long ring is not implemented.') # Rings larger than 999 Elements cannot be translated (can easily be extended if necessary)
                    
                current_smiles=current_smiles[0:find_num_1]+current_smiles[find_num_1+4:find_num_2]+ring_sizeSymbol+current_smiles[find_num_2+4:]
            else:
                break
            
        all_smiles_new=all_smiles_new+'.'+current_smiles
        
    return all_smiles_new[1:]
  
    
def _cumsum(int_list,cum_offset=0): # cumulative sum, without numpy (such that we dont need to include numpy at all)
    cum_list=[cum_offset]
    for x in int_list:
        cum_list.append(cum_list[-1]+x)
    return cum_list[1:]
    


def _get_next_selfies_symbol(tmp_ds): # get the next selfies symbol
    next_symbol=''
    tmp_ds_new=tmp_ds
    if len(tmp_ds)<=2:
        return [next_symbol, tmp_ds_new]
    
    if tmp_ds[0]!='[':
        raise ValueError('_get_next_selfies_symbol: Decoding Problem 1: '+tmp_ds)
    
    end_of_symbol=tmp_ds.find(']')
    
    if end_of_symbol==-1:
        raise ValueError('_get_next_selfies_symbol: Decoding Problem 2: '+tmp_ds)
    else:
        next_symbol=tmp_ds[0:end_of_symbol+1]
        tmp_ds_new=tmp_ds_new[end_of_symbol+1:]
        
    return [next_symbol, tmp_ds_new]
    
    


def __selfies_to_smiles_derive(selfies,smiles,N_restrict=True):    
    # Elements of start_alphabet, again, stand for integers (see comments in _smiles_to_selfies function for more details)
    start_alphabet=['[epsilon]','[F]','[=O]','[#N]','[O]','[N]','[=N]','[C]','[=C]','[#C]','[Branch1_1]','[Branch1_2]','[Branch1_3]','[Ring1]'];
    tmp_ds=selfies.replace('X','Z!') # X will be used as states of the derivation
    
    next_X=smiles.find('X');
    while next_X>=0:
        before_smiles=smiles[0:next_X] # smiles before the non-terminal
        if smiles[next_X+1]=='9':
            state=int(smiles[next_X+1:next_X+5])  # states after branches are called X999...
            after_smiles=smiles[next_X+6:] # smiles after the non-terminal            
        else:            
            state=int(smiles[next_X+1]) # the state is given by the nonterminal symbol X_n, where n=state        
            after_smiles=smiles[next_X+2:] # smiles after the non-terminal
            
        [current_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds) # the current selfies symbol gives the rule vector, and the current state indentifies the one specific, current rule.
        
        # The semantic informations of this set of rules could be significantly extended and more details could be added. Here, we have semantic rules for the most important molecules in organic chemistry, Carbon, Oxygen, Nitrogen, Flour.
        # Other elements get a generic (very weak) restrictions

        if state==0:
            if current_symbol=='[epsilon]':
                new_smiles_symbol='X0'
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0:
                new_smiles_symbol='X0'
                [_,tmp_ds]=_get_next_selfies_symbol(tmp_ds)  # ignore next symbol  
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
                if N_restrict:
                    new_smiles_symbol='[N]X3'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[N]X3'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[N]X3'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X4'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[C]X4'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[C]X4'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X6'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[S]X6'
            else:
                new_smiles_symbol=current_symbol+'X6'
            smiles=before_smiles+new_smiles_symbol+after_smiles

        if state==1:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''
            elif current_symbol.find('Branch')>=0:   
                new_smiles_symbol='X1'
                [_,tmp_ds]=_get_next_selfies_symbol(tmp_ds)  # ignore next symbol               
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='5'
                
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')
            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='5'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[S]X5'
            else:
                new_smiles_symbol=current_symbol+'X6'                
            smiles=before_smiles+new_smiles_symbol+after_smiles

        if state==2:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''              
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                     
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='5'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X1'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X1'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X1'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                      
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='5'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X1'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X1'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X1'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]' or current_symbol=='[Branch1_2]' or current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'

            elif current_symbol=='[Branch2_1]' or current_symbol=='[Branch2_2]' or current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                

            elif current_symbol=='[Branch3_1]' or current_symbol=='[Branch3_2]' or current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol

                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X6'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X6'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X6'
            smiles=before_smiles+new_smiles_symbol+after_smiles


        if state==3:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''             
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                  
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='5'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X2'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X2'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X2'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]                    
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='5'
                    
                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X2'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X2'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X2'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]' or current_symbol=='[Branch1_2]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
            elif current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'             
                
            elif current_symbol=='[Branch2_1]' or current_symbol=='[Branch2_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
                
            elif current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                
                

            elif current_symbol=='[Branch3_1]' or current_symbol=='[Branch3_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
                
            elif current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X6'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[#N]'
                else:
                    new_smiles_symbol='[#N]X6'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X6'

            smiles=before_smiles+new_smiles_symbol+after_smiles


        if state==4:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''            
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='5'                

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X3'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X3'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X3'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                new_smiles_symbol='X4'
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='5'

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X3'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X3'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X3'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
            elif current_symbol=='[Branch1_2]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'              
                
            elif current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                
            elif current_symbol=='[Branch2_1]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                

            elif current_symbol=='[Branch2_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'            
                
                
            elif current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)   
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X1'
                
            elif current_symbol=='[Branch3_1]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                

            elif current_symbol=='[Branch3_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'         
                
                
            elif current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X6'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[#N]'
                else:
                    new_smiles_symbol='[#N]X6'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X6'
            smiles=before_smiles+new_smiles_symbol+after_smiles



        if state==5:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''            
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='5'                

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X4'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X4'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X4'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                new_smiles_symbol='X5'
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='5'

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X4'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X4'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X4'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'
                
            elif current_symbol=='[Branch1_2]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1

                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'              

            elif current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1

                branch_str=''
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                
            elif current_symbol=='[Branch2_1]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'
                

            elif current_symbol=='[Branch2_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'            
                
                
            elif current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)   
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'
                                
            elif current_symbol=='[Branch3_1]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'
                
            elif current_symbol=='[Branch3_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'         
                
            elif current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X2'     

                
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X6'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[#N]'
                else:
                    new_smiles_symbol='[#N]X6'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X6'
            smiles=before_smiles+new_smiles_symbol+after_smiles




        if state==6:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''            
            elif current_symbol.find('Ring1]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    ring_num=str(start_alphabet.index(next_symbol)+2)
                else:
                    ring_num='5'                

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X5'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X5'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X5'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol.find('Ring2]')>=0:
                pre_symbol=''
                if current_symbol[1:5]=='Expl': # Explicit Bond Information
                    pre_symbol=current_symbol[5]
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                new_smiles_symbol='X6'
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    ring_num_1=(start_alphabet.index(next_symbol1)+1)*20
                    ring_num_2=start_alphabet.index(next_symbol2)                
                    ring_num=str(ring_num_1+ring_num_2)
                else:
                    ring_num='5'

                if len(ring_num)==1:
                    new_smiles_symbol=pre_symbol+'%00'+ring_num+'X5'
                elif len(ring_num)==2:
                    new_smiles_symbol=pre_symbol+'%0'+ring_num+'X5'
                elif len(ring_num)==3:
                    new_smiles_symbol=pre_symbol+'%'+ring_num+'X5'
                else:
                    raise ValueError('__selfies_to_smiles_derive: Problem with deriving very long ring.')

            elif current_symbol=='[Branch1_1]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'
                
            elif current_symbol=='[Branch1_2]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1

                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X5'              

            elif current_symbol=='[Branch1_3]':
                [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if next_symbol in start_alphabet:
                    branch_num=start_alphabet.index(next_symbol)+1
                else:
                    branch_num=1

                branch_str=''
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'
                
            elif current_symbol=='[Branch2_1]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'
                

            elif current_symbol=='[Branch2_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X5'            
                
                
            elif current_symbol=='[Branch2_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)   
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*20
                    branch_num2=start_alphabet.index(next_symbol2)
                    branch_num=branch_num1+branch_num2
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'
                                
            elif current_symbol=='[Branch3_1]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9992',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X4'
                
            elif current_symbol=='[Branch3_2]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9991',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X5'         
                
            elif current_symbol=='[Branch3_3]':
                [next_symbol1,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol2,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                [next_symbol3,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                if (next_symbol1 in start_alphabet) and (next_symbol2 in start_alphabet) and (next_symbol3 in start_alphabet):
                    branch_num1=(start_alphabet.index(next_symbol1)+1)*400
                    branch_num2=(start_alphabet.index(next_symbol2))*20
                    branch_num3=start_alphabet.index(next_symbol3)
                    branch_num=branch_num1+branch_num2+branch_num3
                else:
                    branch_num=1
                
                branch_str=''                
                for bii in range(branch_num):
                    [next_symbol,tmp_ds]=_get_next_selfies_symbol(tmp_ds)
                    branch_str+=next_symbol
                    
                branch_smiles=__selfies_to_smiles_derive(branch_str,'X9993',N_restrict)
                new_smiles_symbol=''
                if len(branch_smiles)>0:
                    new_smiles_symbol='('+branch_smiles+')X3'     

                
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X6'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X6'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[#N]'
                else:
                    new_smiles_symbol='[#N]X6'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X6'
            smiles=before_smiles+new_smiles_symbol+after_smiles



        if state==9991: # states 5-7 occure after branches are derived, because a branch or a ring directly after a branch is syntactically illegal.
                     # state  5 corresponds to state 1, state 6 corresponds to state 2, and state 7 corresponds to state 3, without branches & rings
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0: 
                new_smiles_symbol='X9991'
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X4'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[S]X5'
            else:
                new_smiles_symbol=current_symbol+'X6'
            smiles=before_smiles+new_smiles_symbol+after_smiles



        if state==9992:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''          
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0: 
                new_smiles_symbol='X9992' 
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X4'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X6'
            smiles=before_smiles+new_smiles_symbol+after_smiles




        if state==9993:
            if current_symbol=='[epsilon]':
                new_smiles_symbol=''        
            elif current_symbol.find('Ring')>=0 or current_symbol.find('Branch')>=0: 
                new_smiles_symbol='X9993'     
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
                if N_restrict:
                    new_smiles_symbol='[N]X2'
                else:
                    new_smiles_symbol='[N]X4'
            elif current_symbol=='[=N]':
                if N_restrict:
                    new_smiles_symbol='[=N]X1'
                else:
                    new_smiles_symbol='[=N]X4'
            elif current_symbol=='[#N]':
                if N_restrict:
                    new_smiles_symbol='[#N]'
                else:
                    new_smiles_symbol='[#N]X4'
            elif current_symbol=='[C]':
                new_smiles_symbol='[C]X3'
            elif current_symbol=='[=C]':
                new_smiles_symbol='[=C]X2'
            elif current_symbol=='[#C]':
                new_smiles_symbol='[#C]X1'
            elif current_symbol=='[S]':
                new_smiles_symbol='[S]X5'
            elif current_symbol=='[=S]':
                new_smiles_symbol='[=S]X4'
            else:
                new_smiles_symbol=current_symbol+'X6'

            smiles=before_smiles+new_smiles_symbol+after_smiles
    
        if len(tmp_ds)<=2: # if all selfies symbols are derived, the final non-terminals are removed
            while True:
                non_terminal=smiles.find('X')
                if non_terminal>=0:
                    if smiles[non_terminal+1]=='9':
                        smiles=smiles[0:non_terminal]+smiles[non_terminal+5:]
                    else:
                        smiles=smiles[0:non_terminal]+smiles[non_terminal+2:]
                else:
                    break;

        next_X=smiles.find('X')
        
    return smiles.replace('Z!','X')



def _selfies_to_smiles(selfies,N_restrict=True): # here we derive molecule by molecule
    all_selfies=selfies.split('.') # the dot symbol characterizes the start of a new, independent molecule
    all_selfies_new=''

    for current_smiles in all_selfies:
        all_selfies_new=all_selfies_new+'.'+__selfies_to_smiles_derive(current_smiles,'X0',N_restrict) # derive the current selfies string in __selfies_to_smiles_derive (this is a recursive function, and X0 is the starting state of the derivation)
            
    return all_selfies_new[1:] 




def _insert_rings_to_smiles(smiles,N_restrict=True,bilocal_ring_function=True):
    # For simplicity, rings are inserted after the full derivation. However, this does not
    # change anything in the concept of the Grammar - rings could likewise be included
    # directly at the stage of the grammar derivation. It is more compact afterwards.
    
    # No posteriori correction algorithm is included, this function consists only of
    # ring-inclusion with the additional property of verifying bond-accessability at the
    # destination of the ring. This is the only place where SELFIES could have been invalid.
    
    # This so-called bi-local ring-function (bond information at origin is guaranteed
    # by grammar, at destination by the extended ring-function) leads to extremly high
    # validities beyond 99.99% even for random strings.
    
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

    # one additional rule, which can easily encoded into the grammar itself (similar as states X0)
    # for the initial state, or X9991-X9993 for the states after a branch, one could have three initial state classes)
    

    atom1=smiles[1:].find(']')
    if len(smiles)>atom1+2:
        if smiles[atom1+2]=='_':    
            smiles=smiles[0:atom1+2]+smiles[atom1+6:] # after first atom
            atom1=smiles[1:].find(']')

    if atom1>=0:
        atom2=smiles[atom1+2:].find(']')
        if atom2>=0 and len(smiles)>atom1+atom2+3:
            if smiles[atom1+atom2+3]=='_':
                smiles=smiles[0:atom1+atom2+3]+smiles[atom1+atom2+7:] # after second atom
            elif smiles[atom1+atom2+3]=='(':
                pos_end=smiles[atom1+atom2+4:].find(')')
                if len(smiles)>atom1+atom2+5+pos_end:
                    if smiles[atom1+atom2+5+pos_end]=='_':
                        smiles=smiles[0:atom1+atom2+5+pos_end]+smiles[atom1+atom2+9+pos_end:]

    if len(smiles)>atom1+2:
        if smiles[atom1+2]=='(':
            b_cnt=1
            cc=atom1+3
            while b_cnt>0:
                if smiles[cc]==')':            
                    b_cnt-=1
                if smiles[cc]=='(':            
                    b_cnt+=1
                cc+=1
            bracket1=atom1+2
            bracket2=cc
            bracket_cont=smiles[bracket1:bracket2]
            smiles=smiles[0:atom1+2]+smiles[cc:]

            atom1=smiles[1:].find(']')
            if len(smiles)>atom1+2:
                if smiles[atom1+2]=='_':
                    smiles=smiles[0:atom1+2]+smiles[atom1+6:] # after first atom (branch1 removed)
                    atom1=smiles[1:].find(']')
                    
                    if len(smiles)>atom1+2:
                        if smiles[atom1+2]=='_':
                            smiles=smiles[0:atom1+2]+smiles[atom1+6:] # double ring here
                            atom1=smiles[1:].find(']')

            if atom1>=0:
                atom2=smiles[atom1+2:].find(']')
                if atom2>=0 and len(smiles)>atom1+atom2+3:
                    if smiles[atom1+atom2+3]=='_':
                        smiles=smiles[0:atom1+atom2+3]+smiles[atom1+atom2+7:] # after second atom (branch1 removed)
            smiles=smiles[0:bracket1]+bracket_cont+smiles[bracket1:]
        else:
            atom2=smiles[atom1+2:].find(']')
            if atom2>=0 and len(smiles)>atom1+atom2+3:
                
                if smiles[atom1+atom2+3]=='(':
                    
                    b_cnt=1
                    cc=atom1+atom2+4
                    while b_cnt>0:
                        if smiles[cc]==')':            
                            b_cnt-=1
                        if smiles[cc]=='(':            
                            b_cnt+=1
                        cc+=1                    
                    
                    pos_bracket=cc
                    if len(smiles)>pos_bracket:
                        if smiles[pos_bracket]=='_':
                            smiles=smiles[0:pos_bracket]+smiles[pos_bracket+4:]

    pos=smiles.find(')_')
    while pos>=0:
        b_cnt=1
        cc=smiles.find(')_')-1
        cc_init=cc
        while b_cnt>0 and cc>0:
            if smiles[cc]==')':
                b_cnt+=1
            if smiles[cc]=='(':
                b_cnt-=1
            cc-=1
        smiles=smiles[0:cc+1]+smiles[cc_init+2:cc_init+6]+smiles[cc+1:cc_init+2]+smiles[cc_init+6:]
        pos=smiles.find(')_')
         


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

                # this is bi-local ring-function extension. it checks the bond structure at
                # target, and inserts ring only if number of bonds is not exhausted
                allowed_bond_at_target=42   # some large number
                if bilocal_ring_function:
                    target_atom_type=smiles[r_count-1]
                    if target_atom_type=='F':
                        allowed_bond_at_target=1
                    elif target_atom_type=='O':
                        allowed_bond_at_target=2
                    elif target_atom_type=='N':
                        if N_restrict:
                            allowed_bond_at_target=3
                    elif target_atom_type=='C':
                        allowed_bond_at_target=4    
                    elif target_atom_type=='S':
                        allowed_bond_at_target=6                          
    
                    if smiles[1:r_count].find('[')>=0: # bond with previous 
                        if smiles[1:r_count].find('.')>=0:
                            tmp_smiles=smiles[smiles.find('.')+1:r_count]
                            if tmp_smiles[1:].find('[')>=0:
                                allowed_bond_at_target-=1
                        else:
                            allowed_bond_at_target-=1
                    if smiles[r_count-2]=='=':
                        allowed_bond_at_target-=1 # double bond with previous
                    if smiles[r_count-2]=='#':
                        allowed_bond_at_target-=2 # triple bond with previous
                    while smiles[r_count+1]=='@' and allowed_bond_at_target>0:                        
                        allowed_bond_at_target-=1 # already ring here
                        r_count+=4                    
                    if smiles[r_count+1]=='(':
                        allowed_bond_at_target-=1 # already branch here
                        
                        if smiles[r_count+3]=='=':
                            allowed_bond_at_target-=1 # double bond branch
                        if smiles[r_count+3]=='#':
                            allowed_bond_at_target-=2 # triple bond branch           
                        
                        tmp_str=smiles[r_count+1:]
                        b_cnt=1
                        gg=1
                        while b_cnt>0:
                            if tmp_str[gg]==')':            
                                b_cnt-=1
                            if tmp_str[gg]=='(':            
                                b_cnt+=1
                            gg+=1                                        
                        tmp_str2=tmp_str[gg:]
                        if len(tmp_str2)>0:
                            if tmp_str2[0]=='[':
                                allowed_bond_at_target-=1 # bond with atom afterwards 
                            if tmp_str2[0]=='_':
                                allowed_bond_at_target-=1 # ring after bond end
                            if tmp_str2[0]=='(':
                                allowed_bond_at_target-=1 # branch after bond end
                                if tmp_str2[2]=='=':
                                    allowed_bond_at_target-=1 # double branch after bond end
                                if tmp_str2[2]=='#':
                                    allowed_bond_at_target-=2 # double branch after bond end                                    
                                
                                b_cnt=1
                                gg=1
                                while b_cnt>0:
                                    if tmp_str2[gg]==')':            
                                        b_cnt-=1
                                    if tmp_str2[gg]=='(':            
                                        b_cnt+=1
                                    gg+=1
                                tmp_str3=tmp_str2[gg:]
                                if len(tmp_str3)>0:
                                    allowed_bond_at_target-=1 # branch, branch, something
                                    if len(tmp_str3)>1:
                                        if tmp_str3[1]=='=':
                                            allowed_bond_at_target-=1 # double branch after bond end
                                        if tmp_str3[1]=='#':
                                            allowed_bond_at_target-=2 # double branch after bond end  
                                    if tmp_str3[0]=='(':
                                        b_cnt=1
                                        gg=1
                                        while b_cnt>0:
                                            if tmp_str3[gg]==')':            
                                                b_cnt-=1
                                            if tmp_str3[gg]=='(':            
                                                b_cnt+=1
                                            gg+=1
                                        tmp_str4=tmp_str3[gg:]
                                        if len(tmp_str4)>0:
                                            allowed_bond_at_target-=1 # branch, branch, something
                                            if len(tmp_str4)>1:
                                                if tmp_str4[1]=='=':
                                                    allowed_bond_at_target-=1 # double branch after bond end
                                                if tmp_str4[1]=='#':
                                                    allowed_bond_at_target-=2 # double branch after bond end                                              
                                
                        if len(tmp_str2)>1:
                            if tmp_str2[1]=='=':
                                allowed_bond_at_target-=1 # double bond with atom afterwards                                         
                            if tmp_str2[1]=='#':
                                allowed_bond_at_target-=2 # triple bond with atom afterwards 

                    if smiles[r_count+1]=='[':
                        allowed_bond_at_target-=1 # bond with atom afterwards 
                    if smiles[r_count+2]=='=':
                        allowed_bond_at_target-=1 # double bond with atom afterwards  
                    if smiles[r_count+2]=='#':
                        allowed_bond_at_target-=2 # triple bond with atom afterwards    

                if allowed_bond_at_target>0:
                    ring_symbol=available_nums[ring_counter]
            
                    if smiles[r_count+1]=='@': # There is already an identifier, i'm the subsequent one
                        r_count+=4
                        if smiles[r_count+1]=='@':
                            r_count+=4                    
                    smiles=smiles[0:r_count+1]+ring_symbol+smiles[r_count+1:pos_of_ring_symbol]+ring_symbol+smiles[pos_of_ring_symbol+4:]
                    ring_counter+=1           
                else:
                    smiles=smiles[0:pos_of_ring_symbol]+smiles[pos_of_ring_symbol+4:]   # bi-local analysation of property at ring insertion, leads to non-inclusion of ring
            else:    # this can be seen as two additional states of the derivation.  X0 -> X0A -> X0B, and in X0* rings are not allowed. 
                smiles=smiles[0:pos_of_ring_symbol]+smiles[pos_of_ring_symbol+4:]

        else:
            break

    # extended ring function, which has access to the graph structure, i.e. branch information      
          
    # X12YYYYYYZ12 stands for a double ring between X and Z.
    # rewrite to match definition of double ring in RDKit as X1YYYYYYZ=1
    dd=0
    
    while True:
        if dd>len(smiles)-9:
            break

        if smiles[dd:dd+2]=='@$' and smiles[dd+4:dd+6]=='@$': # double ring input
            double_ring_idf=smiles[dd:dd+8]
            tmp_smiles=smiles[dd+5:]
            if tmp_smiles.find(double_ring_idf)>=0: # double ring between same two vertices
                smiles=smiles[0:dd+4]+smiles[dd+8:]
                pos2=smiles.find(double_ring_idf)
                smiles=smiles[0:pos2]+'='+smiles[pos2:pos2+4]+smiles[pos2+8:]
                dd=0
        dd+=1
                  
    tmp_smiles=smiles.replace('(','').replace(')','')
    for unique_id in available_nums:
        pos1=tmp_smiles.find(unique_id)
        if pos1>=0:
            tmp_smiles2=tmp_smiles[pos1+4:]
            pos2=tmp_smiles2.find(unique_id)                
            indices = [i for i, x in enumerate(tmp_smiles2[0:pos2]) if x == '[']
            if len(indices)<2:
                for twice in range(2):
                    pos1=smiles.find(unique_id)
                    if smiles[pos1-1]=='=':
                        smiles=smiles[0:pos1-1]+smiles[pos1+4:]
                    else:
                        smiles=smiles[0:pos1]+smiles[pos1+4:]
                        
    gg=0
    while gg<len(smiles):
        tmp_smiles0=smiles[gg:]
        if tmp_smiles0.find('(')==-1:
            gg=len(smiles)
        else:
            b_cnt=1
            cc=tmp_smiles0.find('(')+1
            cc_init=cc-1
            while b_cnt>0 and cc<len(tmp_smiles0):
                if tmp_smiles0[cc]==')':
                    b_cnt-=1
                if tmp_smiles0[cc]=='(':
                    b_cnt+=1
                cc+=1
            
            tmp_smiles=tmp_smiles0[0:cc_init]+tmp_smiles0[cc:]
    
            for unique_id in available_nums:
                pos1=tmp_smiles.find(unique_id)
                indices = [i for i, x in enumerate(tmp_smiles) if x =='@']
                if pos1>=0 and len(indices)>=2:
                    tmp_smiles2=tmp_smiles[pos1+4:]
                    pos2=tmp_smiles2.find(unique_id)
                    if pos2==-1:
                        break
                    else:
                        indices = [i for i, x in enumerate(tmp_smiles2[0:pos2]) if x == '[']                    
                        if len(indices)<2:
                            for twice in range(2):
                                pos1=smiles.find(unique_id)
                                if smiles[pos1-1]=='=':
                                    smiles=smiles[0:pos1-1]+smiles[pos1+4:]
                                else:
                                    smiles=smiles[0:pos1]+smiles[pos1+4:]
                                        
            indices = [i for i, x in enumerate(tmp_smiles0) if x == '@']
            if len(indices)<2:
                gg=len(smiles)
            gg+=1

    tmp_smiles=smiles
    while tmp_smiles.find('(')>=0:
        b_cnt=1
        cc=tmp_smiles.find('(')+1
        cc_init=cc-1
        while b_cnt>0 and cc<len(tmp_smiles):
            if tmp_smiles[cc]==')':
                b_cnt-=1
            if tmp_smiles[cc]=='(':
                b_cnt+=1
            cc+=1
        tmp_smiles=tmp_smiles[0:cc_init]+tmp_smiles[cc:]

    for unique_id in available_nums:
        pos1=tmp_smiles.find(unique_id)
        indices = [i for i, x in enumerate(tmp_smiles) if x =='@']
        if pos1>=0 and len(indices)>=2:
            tmp_smiles2=tmp_smiles[pos1+4:]
            pos2=tmp_smiles2.find(unique_id)
            if pos2>=0:
                indices = [i for i, x in enumerate(tmp_smiles2[0:pos2]) if x == '[']                    
                if len(indices)<2:
                    for twice in range(2):
                        pos1=smiles.find(unique_id)
                        if smiles[pos1-1]=='=':
                            smiles=smiles[0:pos1-1]+smiles[pos1+4:]
                        else:
                            smiles=smiles[0:pos1]+smiles[pos1+4:]

    # #########
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
        
    smiles=smiles.replace('[=','=[') # all special symbols for bond or sterochemical information are now outside of the bracket again
    smiles=smiles.replace('[#','#[')
    smiles=smiles.replace('[\\','\\[')
    smiles=smiles.replace('[/','/[')
    smiles=smiles.replace('[-','-[')
    
    smiles=smiles.replace('==','#') # RDKit syntax
            
    smiles=smiles.replace('[B]','B') # brackets are removed for these molecules (according to SMILES standard)
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
    
    smiles=smiles.replace('expl]',']')
    
    
    smiles=smiles.replace('%00','') # finally, make rings in a standard form.
    smiles=smiles.replace('%0','%')

    return(smiles)
    
    
    
    
    

def encoder(smiles,PrintErrorMessage=True): # encodes SMILES -> SELFIES
    """
    SELFIES: a robust representation of semantically constrained graphs with an example application in chemistry
                  v0.2.0, 02. September 2019
    by Mario Krenn, Florian Haese, AkshatKuman Nigam, Pascal Friederich, Alan Aspuru-Guzik

    SELFIES (SELF-referencIng Embedded Strings) is a general-purpose, sequence-based,
    robust representation of semantically constrained graphs. It is based on a Chomsky
    type-2 grammar, augmented with two self-referencing functions. A main objective is
    to use SELFIES as direct input into machine learning models, in particular
    in generative models, for the generation of outputs with high validity.

    The code presented here is a concrete application of SELIFES in chemistry, for
    the robust representation of molecule. We show the encoding and decoding of three
    molecules from various databases, and the generation of a new, random molecule
    with high semantical and syntactical validity.
    
    Example:
        test_molecule1='CN1C(=O)C2=C(c3cc4c(s3)-c3sc(-c5ncc(C#N)s5)cc3C43OCCO3)N(C)C(=O)C2=C1c1cc2c(s1)-c1sc(-c3ncc(C#N)s3)cc1C21OCCO1' # non-fullerene acceptors for organic solar cells
        selfies1=encoder(test_molecule1)
        smiles1=decoder(selfies1)
        
        selfies1 -> '[C][N][C][Branch1_3][epsilon][=O][C][=C][Branch2_3][epsilon][S][c][c][c][c][Branch1_3][Ring2][s]
                     [Ring1][Ring2][-c][s][c][Branch1_3][O][-c][n][c][c][Branch1_3][Ring1][C][#N][s][Ring1][Branch1_2]
                     [c][c][Ring1][F][C][Ring1][=N][O][C][C][O][Ring1][Ring2][N][Branch1_3][epsilon][C][C][Branch1_3]
                     [epsilon][=O][C][Ring2][epsilon][Branch2_3][=C][Ring2][epsilon][N][c][c][c][c][Branch1_3][Ring2]
                     [s][Ring1][Ring2][-c][s][c][Branch1_3][O][-c][n][c][c][Branch1_3][Ring1][C][#N][s][Ring1][Branch1_2]
                     [c][c][Ring1][F][C][Ring1][=N][O][C][C][O][Ring1][Ring2]'
        smiles1 -> 'CN1C(=O)C2=C(c3cc4c(s3)-c3sc(-c5ncc(C#N)s5)cc3C43OCCO3)N(C)C(=O)C2=C1c1cc2c(s1)-c1sc(-c3ncc(C#N)s3)cc1C21OCCO1'
        
    With the following alphabet
        my_alphabet=['[Branch1_1]','[Branch1_2]','[Branch1_3]','[Ring1]','[O]','[=O]','[N]','[=N]','[C]','[=C]','[#C]','[S]','[=S]','[P]','[F]'];
    random SELFIES lead to validity of >99.99% (for string lengths of 50 symbols)
    
    For comments, bug reports or feature ideas, please send an email to
    mario.krenn@utoronto.ca and alan@aspuru.com
    """
    try:
        preselfies1=_make_brackets_around_atoms(smiles)        # Itemize Atoms
        preselfies2=_reconfigure_smiles_numbers1(preselfies1)  # Unique Ringsymbols
        preselfies3=_reconfigure_smiles_numbers2(preselfies2)  # Rings as relative distance
        selfies=_smiles_to_selfies(preselfies3)                # Create selfies
    except ValueError as err:
        if PrintErrorMessage:
            print(err)
            print('Could not encode smiles string. Please contact authors.')
        return -1
        
    return selfies
    



def decoder(selfies,N_restrict=True,bilocal_ring_function=True,PrintErrorMessage=True): # decodes SELFIES -> SMILES
    smiles=-1
    if selfies!=-1:
        try:
            presmiles1=_selfies_to_smiles(selfies,N_restrict)     # Runs Grammar Rules
            smiles=_insert_rings_to_smiles(presmiles1,N_restrict,bilocal_ring_function)       # Inserts Rings
        except  ValueError as err:
            if PrintErrorMessage:
                print(err)            
                print('Could not decode selfies string. Please contact authors.')
            return -1
    
    return smiles
    
