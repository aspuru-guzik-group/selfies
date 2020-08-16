# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 14:17:55 2019

@author: MK
"""

def GrammarPlusToSMILES(grammarString,SMILES):
    # for initial input, SMILES has to be 'X0'
    # State can be 0-6, where i=0...4=Xi, and i=5: NUM1, and i=6: NUM2
    
    #print('---------------------------------')
    #print('---------------------------------')
    #print('   INPUT:')
    #print('       grammarString: ', grammarString)
    #print('       SMILES: ', SMILES) 
    #print(' ')

#    grammarString='HC'
#    SMILES='X0'

    nextX=SMILES.find('X')
    #print('nextX: ', nextX)
    #print('SMILES: ', SMILES)
    while nextX>=0 and len(grammarString)>0:
        curSMILES=''
        currentState=int(SMILES[nextX+1])
        nextSymbol=grammarString[0]
        #print(' ')
        #print('grammarString: ',grammarString)
        #print('SMILES: ', SMILES)
        #print('currentState: ',currentState)
        grammarString=grammarString[1:]
        if currentState==0:
            if nextSymbol=='A':
                # do nothing
                curSMILES='X0'
            elif nextSymbol=='B':
                curSMILES='FX1'
            elif nextSymbol=='C':
                curSMILES='OX2'
            elif nextSymbol=='D':
                curSMILES='NX3'
            elif nextSymbol=='E':
                curSMILES='OX2'
            elif nextSymbol=='F':
                curSMILES='NX3'
            elif nextSymbol=='G':
                curSMILES='NX3'
            elif nextSymbol=='H':
                curSMILES='CX4'
            elif nextSymbol=='I':
                curSMILES='CX4'
            elif nextSymbol=='J':
                curSMILES='CX4'
            elif nextSymbol=='K':
                curSMILES='X0'
                grammarString=grammarString[1:] # 'Y', delete next symbol
            elif nextSymbol=='L':
                curSMILES='X0'
                grammarString=grammarString[1:] # 'Y', delete next symbol
            elif nextSymbol=='M':
                curSMILES='X0'
                grammarString=grammarString[1:] # 'Y', delete next symbol
            elif nextSymbol=='N':
                curSMILES='X0'
                grammarString=grammarString[1:] # 'Y', delete next symbol
                
                
        elif currentState==1:
            if nextSymbol=='A':
                curSMILES=''
            elif nextSymbol=='B':
                curSMILES='F'
            elif nextSymbol=='C':
                curSMILES='O'
            elif nextSymbol=='D':
                curSMILES='N'
            elif nextSymbol=='E':
                curSMILES='OX1'
            elif nextSymbol=='F':
                curSMILES='NX2'
            elif nextSymbol=='G':
                curSMILES='NX2'
            elif nextSymbol=='H':
                curSMILES='CX3'
            elif nextSymbol=='I':
                curSMILES='CX3'
            elif nextSymbol=='J':
                curSMILES='CX3'
            elif nextSymbol=='K':
                curSMILES='X1'
                grammarString=grammarString[1:] # 'Y', delete next symbol
            elif nextSymbol=='L':
                curSMILES='X1'
                grammarString=grammarString[1:] # 'Y', delete next symbol
            elif nextSymbol=='M':
                curSMILES='X1'
                grammarString=grammarString[1:] # 'Y', delete next symbol                
            elif nextSymbol=='N':                
                if len(grammarString)>0:
                    nextSymbol=grammarString[0]                
                    curSMILES=chr(ord(nextSymbol.lower()))
                    grammarString=grammarString[1:]
                
        elif currentState==2:
            if nextSymbol=='A':
                curSMILES=''
            elif nextSymbol=='B':
                curSMILES='F'
            elif nextSymbol=='C':
                curSMILES='=O'
            elif nextSymbol=='D':
                curSMILES='=N'
            elif nextSymbol=='E':
                curSMILES='OX1'
            elif nextSymbol=='F':
                curSMILES='NX2'
            elif nextSymbol=='G':
                curSMILES='=NX1'
            elif nextSymbol=='H':
                curSMILES='CX3'
            elif nextSymbol=='I':
                curSMILES='=CX2'
            elif nextSymbol=='J':
                curSMILES='=CX2'
            elif nextSymbol=='K':            
                if len(grammarString)>0:
                    nextSymbol=grammarString[0]
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols+1]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X5')
                    grammarString=grammarString[NumOfSymbols+1:]               
                    curSMILES='('+BranchSMILES+')X1'                
            elif nextSymbol=='L':
                if len(grammarString)>0:                
                    nextSymbol=grammarString[0]
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols+1]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X5')
                    grammarString=grammarString[NumOfSymbols+1:]
                    curSMILES='('+BranchSMILES+')X1'    
            elif nextSymbol=='M':
                if len(grammarString)>0:                
                    nextSymbol=grammarString[0]
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols+1]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X5')
                    grammarString=grammarString[NumOfSymbols+1:]
                    curSMILES='('+BranchSMILES+')X1'                     
            elif nextSymbol=='N':
                if len(grammarString)>0:                
                    nextSymbol=grammarString[0]                 
                    curSMILES=chr(ord(nextSymbol.lower()))+'X1'
                    grammarString=grammarString[1:] 
                
        elif currentState==3:
            if nextSymbol=='A':
                curSMILES=''
            elif nextSymbol=='B':
                curSMILES='F'
            elif nextSymbol=='C':
                curSMILES='=O'
            elif nextSymbol=='D':
                curSMILES='#N'
            elif nextSymbol=='E':
                curSMILES='OX1'
            elif nextSymbol=='F':
                curSMILES='NX2'
            elif nextSymbol=='G':
                curSMILES='=NX1'
            elif nextSymbol=='H':
                curSMILES='CX3'
            elif nextSymbol=='I':
                curSMILES='=CX2'
            elif nextSymbol=='J':
                curSMILES='#CX1'
            elif nextSymbol=='K':            
                if len(grammarString)>0:
                    nextSymbol=grammarString[0] 
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols+1]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X5')
                    grammarString=grammarString[NumOfSymbols+1:]  
                    curSMILES='('+BranchSMILES+')X2'                
            elif nextSymbol=='L':
                if len(grammarString)>0:
                    nextSymbol=grammarString[0]
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols+1]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X6')
                    grammarString=grammarString[NumOfSymbols+1:]
                    curSMILES='('+BranchSMILES+')X1'
            elif nextSymbol=='M':            
                if len(grammarString)>0:
                    nextSymbol=grammarString[0] 
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols+1]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X5')
                    grammarString=grammarString[NumOfSymbols+1:]  
                    curSMILES='('+BranchSMILES+')X2'                       
            elif nextSymbol=='N':
                if len(grammarString)>0:
                    nextSymbol=grammarString[0]                 
                    curSMILES=chr(ord(nextSymbol.lower()))+'X2'
                    grammarString=grammarString[1:]      

        elif currentState==4:
            if nextSymbol=='A':
                curSMILES=''
            elif nextSymbol=='B':
                curSMILES='F'
            elif nextSymbol=='C':
                curSMILES='=O'
            elif nextSymbol=='D':
                curSMILES='#N'
            elif nextSymbol=='E':
                curSMILES='OX1'
            elif nextSymbol=='F':
                curSMILES='NX2'
            elif nextSymbol=='G':
                curSMILES='=NX1'
            elif nextSymbol=='H':
                curSMILES='CX3'
            elif nextSymbol=='I':
                curSMILES='=CX2'
            elif nextSymbol=='J':
                curSMILES='#CX1'
            elif nextSymbol=='K':            
                if len(grammarString)>0:                
                    nextSymbol=grammarString[0] 
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X5')
                    grammarString=grammarString[NumOfSymbols+1:]
                    curSMILES='('+BranchSMILES+')X3'                
            elif nextSymbol=='L':
                if len(grammarString)>0:                
                    nextSymbol=grammarString[0]
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X7')
                    grammarString=grammarString[NumOfSymbols+1:]
                    curSMILES='('+BranchSMILES+')X1' 
            elif nextSymbol=='M':
                if len(grammarString)>0:                
                    nextSymbol=grammarString[0]
                    NumOfSymbols=ord(nextSymbol.lower())-97+1
                    BranchGrammarString=grammarString[1:NumOfSymbols]
                    BranchSMILES=GrammarPlusToSMILES(BranchGrammarString,'X6')
                    grammarString=grammarString[NumOfSymbols+1:]
                    curSMILES='('+BranchSMILES+')X2'                     
            elif nextSymbol=='N':
                if len(grammarString)>0:
                    nextSymbol=grammarString[0]                 
                    curSMILES=chr(ord(nextSymbol.lower()))+'X3'
                    grammarString=grammarString[1:]




        elif currentState==5:
            if nextSymbol=='A':
                curSMILES='C'
            elif nextSymbol=='B':
                curSMILES='F'
            elif nextSymbol=='C':
                curSMILES='O'
            elif nextSymbol=='D':
                curSMILES='N'
            elif nextSymbol=='E':
                curSMILES='OX1'
            elif nextSymbol=='F':
                curSMILES='NX2'
            elif nextSymbol=='G':
                curSMILES='NX2'
            elif nextSymbol=='H':
                curSMILES='CX3'
            elif nextSymbol=='I':
                curSMILES='CX3'
            elif nextSymbol=='J':
                curSMILES='CX3'
            elif nextSymbol=='K':
                curSMILES='C'
            elif nextSymbol=='L':
                curSMILES='C'
            elif nextSymbol=='M':
                curSMILES='C'
            elif nextSymbol=='N':
                curSMILES='C'  
                
        elif currentState==6:
            if nextSymbol=='A':
                curSMILES='C'
            elif nextSymbol=='B':
                curSMILES='F'
            elif nextSymbol=='C':
                curSMILES='=O'
            elif nextSymbol=='D':
                curSMILES='=N'
            elif nextSymbol=='E':
                curSMILES='OX1'
            elif nextSymbol=='F':
                curSMILES='NX2'
            elif nextSymbol=='G':
                curSMILES='=NX1'
            elif nextSymbol=='H':
                curSMILES='CX3'
            elif nextSymbol=='I':
                curSMILES='=CX2'
            elif nextSymbol=='J':
                curSMILES='=CX2'
            elif nextSymbol=='K':
                curSMILES='C'
            elif nextSymbol=='L':
                curSMILES='C'
            elif nextSymbol=='M':
                curSMILES='C'
            elif nextSymbol=='N':
                curSMILES='C' 
                
        elif currentState==7:
            if nextSymbol=='A':
                curSMILES=''
            elif nextSymbol=='B':
                curSMILES='F'
            elif nextSymbol=='C':
                curSMILES='=O'
            elif nextSymbol=='D':
                curSMILES='#N'
            elif nextSymbol=='E':
                curSMILES='OX1'
            elif nextSymbol=='F':
                curSMILES='NX2'
            elif nextSymbol=='G':
                curSMILES='=NX1'
            elif nextSymbol=='H':
                curSMILES='CX3'
            elif nextSymbol=='I':
                curSMILES='=CX2'
            elif nextSymbol=='J':
                curSMILES='#CX1'
            elif nextSymbol=='K':
                curSMILES='C'
            elif nextSymbol=='L':
                curSMILES='C'
            elif nextSymbol=='M':
                curSMILES='C'
            elif nextSymbol=='N':
                curSMILES='C'      



        #print('end curSMILES: ', curSMILES)
        #print('end SMILES: ', SMILES)
        SMILES=SMILES[0:nextX]+curSMILES+SMILES[nextX+2:]
        nextX=SMILES.find('X')
        #print('SMILES: ', SMILES)
        #print('nextX: ', nextX)



    # remove remaining terminal symbols
    findXinSMILES=SMILES.find('X')
    while findXinSMILES>=0:
        SMILES=SMILES[0:findXinSMILES]+SMILES[findXinSMILES+2:]
        findXinSMILES=SMILES.find('X')
        
    
    #print('   OUTPUT:')
    #print('       grammarString: ', grammarString)
    #print('       SMILES: ', SMILES) 
    #print('---------------------------------')
    #print('---------------------------------')
    return SMILES
        
def IsLowerCase(onechar):
    rv=0    
    if ord(onechar)>=97 and ord(onechar)<=122:
        rv=1
        
    return rv


def IncludeRingsForSMILES(SMILES):
#    for ii in range(len(SMILES)-1):
#        if IsLowerCase(SMILES[ii]) and IsLowerCase(SMILES[ii+1]):
#            if ord(SMILES[ii])>ord(SMILES[ii+1]):
#                tmp=SMILES[ii]
#                SMILES=list(SMILES)
#                SMILES[ii]=SMILES[ii+1]
#                SMILES[ii+1]=tmp
#                SMILES=''.join(SMILES)
                
        
    #print(' ')
    #print(SMILES)
    RingNumber=1
    for ii in range(14):
        #print(chr(ii+97))
        findLetter=SMILES.find(chr(ii+97))
        while findLetter>=0:
            AtomIdxList=[]
            for jj in range(len(SMILES)):
                if SMILES[jj]=='C' or SMILES[jj]=='N' or SMILES[jj]=='O' or SMILES[jj]=='F':
                    AtomIdxList.append(jj)

            SMILES = SMILES[:findLetter] + str(RingNumber) + SMILES[(findLetter+1):]
            #print('AtomIdxList: ',AtomIdxList)
            
            for jj in range(findLetter)[::-1]:
                if SMILES[jj]=='C' or SMILES[jj]=='N' or SMILES[jj]=='O' or SMILES[jj]=='F':
                    findBeforeAtomPos=jj
                    break
            
            CurrentAtomPosInSMILES=AtomIdxList.index(findBeforeAtomPos)
            #print('CurrentAtomPosInSMILES: ',CurrentAtomPosInSMILES)           
            PositionForRingStart=AtomIdxList[max(0,CurrentAtomPosInSMILES-ii-2)]+1
            #print('PositionForRingStart: ',PositionForRingStart)
            
            SMILES = SMILES[:PositionForRingStart] + str(RingNumber) + SMILES[PositionForRingStart:]
            RingNumber+=1
            
            #print('   SMILES: ',SMILES)
            findLetter=SMILES.find(chr(ii+97))
            
            
    return SMILES


    
#res=IncludeRingsForSMILES(GrammarPlusToSMILES('FIFFGHKDHIMEG','X0'))
#res=IncludeRingsForSMILES('CC=CCCC(O)CCfc')
#print(' ')
#print('Final Result: ', res)