import copy 
import pickle 
import random 
from collections import defaultdict
from collections import Counter 
import threading 
import os





#CLASSES ================================================================================================================================

#Class stores the board and info like which player's turn it is.
class GMPO:
    def __init__(self,brd,plr,cst_rgh,ENPT,HMC,hist = {}):
        self.brd = brd 
        self.plr = plr 
        self.castling = cst_rgh 
        self.EnP = ENPT 
        self.HMC = HMC 
        self.hist = hist   
    def getbrd(self):
        return self.brd
    def setbrd(self,brd):
        self.brd = brd
    def getplr(self):
        return self.plr
    def setplr(self,plr):
        self.plr = plr
    def getCSRG(self):
        return self.castling
    def setCSRG(self,cst_rgh):
        self.castling = cst_rgh
    def getEnP(self):
        return self.EnP
    def setEnP(self, ENPT):
        self.EnP = ENPT
    def getHMC(self):
        return self.HMC
    def setHMC(self,HMC):
        self.HMC = HMC
    def chRE(self):
        return any(value>=3 for value in self.hist.itervalues())
    def gH(self):
        return self.hist
    def clone(self):
        clone = GMPO(copy.deepcopy(self.brd), self.plr,copy.deepcopy(self.castling),self.EnP,self.HMC)
        return clone

#Class stores information of a single chess piece 
class Unit:
    def __init__(self,pINF,chess_crd):
        Unit = pINF[0]
        cl = pINF[1]
        self.pINF = pINF
        self.chess_crd = chess_crd
        self.pos = (-1,-1)
    def getInfo(self):
        return [self.chess_crd, self.subsection,self.pos]
    def setpos(self,pos):
        self.pos = pos
    def getpos(self):
        return self.pos
    def setcrd(self,crd):
        self.chess_crd = crd
    def __repr__(self):
        return self.pINF









#CHESS FUNCTIONS ==========================================================================================================================

#Tells if given co-ordinate on board is not empty
def isOCC(brd,x,y):
    x = int(x)
    y = int(y)
    if brd[y][x] == 0:
        return False
    return True

#Gives true if given co-ordinate on board is occupied by passed color
def isOCCby(brd,x,y,cl):
    x = int(x)
    y = int(y)
    if brd[y][x]==0:
        return False
    if brd[y][x][1] == cl[0]:
        return True
    return False

#Removes co-ordinates with pieces of passed color from list that is passed into function
def FBC(brd,LOT,cl):
    flt_lst = []
    for pos in LOT:
        x = pos[0]
        y = pos[1]
        if x>=0 and x<=7 and y>=0 and y<=7 and not isOCCby(brd,x,y,cl):
            flt_lst.append(pos)
    return flt_lst

#Returns list of co-ordinates occupied by piece that is passed
def lf(brd,Unit):
    LOL = []
    for row in range(8):
        for col in range(8):
            if brd[row][col] == Unit:
                x = col
                y = row
                LOL.append((x,y))
    return LOL

#Checks whether passed square is being attacked by pieces of passed color
def IAB(pst,tx,ty,cl):
    brd = pst.getbrd()
    cl = cl[0]
    LOAS = []
    for x in range(8):
        for y in range(8):
            if brd[y][x]!=0 and brd[y][x][1]==cl:
                LOAS.extend(
                    FPS(pst,x,y,True))
    return (tx,ty) in LOAS 

#Returns list of possible squares that given piece can move to                
def FPS(pst,x,y,AtSe=False):
    x = int(x)
    y = int(y)
    brd = pst.getbrd()
    plr = pst.getplr()
    cst_rgh = pst.getCSRG()
    ENPT = pst.getEnP()
    if len(brd[y][x])!=2: 
        return [] 
    Unit = brd[y][x][0] 
    cl = brd[y][x][1] 
    enm_cl = opp(cl)
    LOT = []

    #pawn ----------------------------------------------------------------------------------------
    if Unit == 'P': 
        if cl=='w':
            if not isOCC(brd,x,y-1) and not AtSe:
                LOT.append((x,y-1))
                if y == 6 and not isOCC(brd,x,y-2):
                    LOT.append((x,y-2))
            if x!=0 and isOCCby(brd,x-1,y-1,'black'):
                LOT.append((x-1,y-1))
            if x!=7 and isOCCby(brd,x+1,y-1,'black'):
                LOT.append((x+1,y-1))
            if ENPT!=-1:
                if ENPT == (x-1,y-1) or ENPT == (x+1,y-1):
                    LOT.append(ENPT)    
        elif cl=='b': 
            if not isOCC(brd,x,y+1) and not AtSe:
                LOT.append((x,y+1))
                if y == 1 and not isOCC(brd,x,y+2):
                    LOT.append((x,y+2))
            if x!=0 and isOCCby(brd,x-1,y+1,'white'):
                LOT.append((x-1,y+1))
            if x!=7 and isOCCby(brd,x+1,y+1,'white'):
                LOT.append((x+1,y+1))
            if ENPT == (x-1,y+1) or ENPT == (x+1,y+1):
                LOT.append(ENPT)

    #Rook ----------------------------------------------------------------------------------------
    elif Unit == 'R': 
        for i in [-1,1]:
            kx = x 
            while True: 
                kx = kx + i 
                if kx<=7 and kx>=0:   
                    if not isOCC(brd,kx,y):
                        LOT.append((kx,y))
                    else:
                        if isOCCby(brd,kx,y,enm_cl):
                            LOT.append((kx,y))
                        break 
                else: 
                    break
        for i in [-1,1]:
            ky = y
            while True:
                ky = ky + i 
                if ky<=7 and ky>=0: 
                    if not isOCC(brd,x,ky):
                        LOT.append((x,ky))
                    else:
                        if isOCCby(brd,x,ky,enm_cl):
                            LOT.append((x,ky))
                        break
                else:
                    break
        
    #knight ----------------------------------------------------------------------------------------
    elif Unit == 'N':
        for dx in [-2,-1,1,2]:
            if abs(dx)==1:
                sy = 2
            else:
                sy = 1
            for dy in [-sy,+sy]:
                LOT.append((x+dx,y+dy))
        LOT = FBC(brd,LOT,cl)
    elif Unit == 'B': 
        for dx in [-1,1]: 
            for dy in [-1,1]: 
                kx = x 
                ky = y
                while True: 
                    kx = kx + dx 
                    ky = ky + dy 
                    if kx<=7 and kx>=0 and ky<=7 and ky>=0:
                        if not isOCC(brd,kx,ky):
                            LOT.append((kx,ky))
                        else:
                            if isOCCby(brd,kx,ky,enm_cl):
                                LOT.append((kx,ky))
                            break    
                    else:
                        break

    #Queen ----------------------------------------------------------------------------------------
    elif Unit == 'Q': 
        brd[y][x] = 'R' + cl
        l_rook = FPS(pst,x,y,True)
        brd[y][x] = 'B' + cl
        l_bishop = FPS(pst,x,y,True)
        LOT = l_rook + l_bishop
        brd[y][x] = 'Q' + cl
    elif Unit == 'K': 
        for dx in [-1,0,1]:
            for dy in [-1,0,1]:
                LOT.append((x+dx,y+dy))
        LOT = FBC(brd,LOT,cl)
        if not AtSe:
            right = cst_rgh[plr]
            if (right[0] and 
                brd[y][7]!=0 and 
                brd[y][7][0]=='R' and 
                not isOCC(brd,x+1,y) and 
                not isOCC(brd,x+2,y) and 
                not IAB(pst,x,y,enm_cl) and 
                not IAB(pst,x+1,y,enm_cl) and 
                not IAB(pst,x+2,y,enm_cl)       ):
                LOT.append((x+2,y))
            if (right[1] and 
                brd[y][0]!=0 and 
                brd[y][0][0]=='R' and 
                not isOCC(brd,x-1,y)and 
                not isOCC(brd,x-2,y)and 
                not isOCC(brd,x-3,y) and 
                not IAB(pst,x,y,enm_cl) and 
                not IAB(pst,x-1,y,enm_cl) and 
                not IAB(pst,x-2,y,enm_cl)       ):
                LOT.append((x-2,y)) 

    if not AtSe:
        n_list = []
        for tpq in LOT:
            x2 = tpq[0]
            y2 = tpq[1]
            temp_pos = pst.clone()
            makemove(temp_pos,x,y,x2,y2)
            if not isCHK(temp_pos,cl):
                n_list.append(tpq)
        LOT = n_list
    return LOT

#Returns opposite color
def opp(cl):
    cl = cl[0]
    if cl == 'w':
        oppcl = 'b'
    else:
        oppcl = 'w'
    return oppcl

#Returns true if king of this color is under check
def isCHK(pst,cl):
    brd = pst.getbrd()
    cl = cl[0]
    enm = opp(cl)
    Unit = 'K' + cl
    x,y = lf(brd,Unit)[0]
    return IAB(pst,x,y,enm)

#Returns true if king of this color is under checkmate
def isCHKmate(pst,cl=-1):
    
    if cl==-1:
        return isCHKmate(pst,'white') or isCHKmate(pst,'b')
    cl = cl[0]
    if isCHK(pst,cl) and allM(pst,cl)==[]:
            return True
    return False

#Returns true if current situation is checkmate
def isSTL(pst):
    plr = pst.getplr()
    if plr==0:
        cl = 'w'
    else:
        cl = 'b'
    if not isCHK(pst,cl) and allM(pst,cl)==[]:
        return True
    return False

#Returns list of all pieces of given color    
def getAU(pst,cl):
    brd = pst.getbrd()
    LOP = []
    for j in range(8):
        for i in range(8):
            if isOCCby(brd,i,j,cl):
                LOP.append((i,j))
    return LOP

#Returns list of all possible moves for pieces of given color
def allM(pst, cl):
    if cl==1:
        cl = 'white'
    elif cl ==-1:
        cl = 'black'
    cl = cl[0]
    listofUnits = getAU(pst,cl)
    moves = []
    for pos in listofUnits:

        targets = FPS(pst,pos[0],pos[1])
        for target in targets:
    
             moves.append([pos,target])
    return moves

#Creates piece instances for board pieces
def crtUNITS(brd):

    LOWU = []
    LOBU= []
  
    for i in range(8):
        for k in range(8):
            if brd[i][k]!=0:
               
                p = Unit(brd[i][k],(k,i))
                
                if brd[i][k][1]=='w':
                    LOWU.append(p)
                else:
                    LOBU.append(p)
   
    return [LOWU,LOBU]









#MAKE MOVES ===============================================================================================================================

#Moves piece from current co-ordinates to given co-ordinates
def makemove(pst,x,y,x2,y2):
    x = int(x)
    y = int(y)
    x2 = int(x2)
    y2 = int(y2)
    brd = pst.getbrd()
    Unit = brd[y][x][0]
    cl = brd[y][x][1]

    plr = pst.getplr()
    cst_rgh = pst.getCSRG()
    ENPT = pst.getEnP()
    HalfMoveClock = pst.getHMC()
    
    if isOCC(brd,x2,y2) or Unit=='P':
        HalfMoveClock = 0
    else:
        HalfMoveClock += 1

    brd[y2][x2] = brd[y][x]
    brd[y][x] = 0
   
    if Unit == 'K':
        cst_rgh[plr] = [False,False]
        if abs(x2-x) == 2:
            if cl=='w':
                l = 7
            else:
                l = 0
            if x2>x:
                brd[l][5] = 'R'+cl
                brd[l][7] = 0
            else:
                brd[l][3] = 'R'+cl
                brd[l][0] = 0

    if Unit=='R':
        if x==0 and y==0:       
            cst_rgh[1][1] = False
        elif x==7 and y==0:     
            cst_rgh[1][0] = False
        elif x==0 and y==7:      
            cst_rgh[0][1] = False
        elif x==7 and y==7:      
            cst_rgh[0][0] = False
    
    if Unit == 'P':  
        if ENPT == (x2,y2):
            if cl=='w':
                brd[y2+1][x2] = 0
            else:
                brd[y2-1][x2] = 0   
        if abs(y2-y)==2:
            ENPT = (x,(y+y2)/2)
        else:
            ENPT = -1   
        if y2==0:
            brd[y2][x2] = 'Qw'
        elif y2 == 7:
            brd[y2][x2] = 'Qb'
    else:
        ENPT = -1

    plr = 1 - plr    
        
    pst.setplr(plr)
    pst.setCSRG(cst_rgh)
    pst.setEnP(ENPT)
    pst.setHMC(HalfMoveClock)










#MINIMAX ALGORITHM =========================================================================================================================

def miniMax(pst,depth,alpha,beta,maximizingplr,bestMoveReturn,root):

    if depth==0:
        return maximizingplr*evaluate(pst)

    moves = allM(pst, maximizingplr) 
    if moves==[]:
        return maximizingplr*evaluate(pst)    
    random.shuffle(moves)
  
    #Player is white , so find advantage for black(By minimizing)
    if maximizingplr == -1:
        bValue = -100000
        maxEval = -1000000
        storage = []
        for item in moves:
            newpos = pst.clone()
            makemove(newpos,item[0][0],item[0][1],item[1][0],item[1][1])

            eval = miniMax(newpos,depth-1,alpha,beta,1,bestMoveReturn,False)
            if(root == True):
                storage.append(eval)
            maxEval = max(eval,maxEval)
            alpha = max(alpha,eval)
            if beta <= alpha:
                break
            if(root == True):
                #Find Min move
                counter = 0
                maxIndex = 0
                maxValue = storage[0]
                for x in storage:
                    if x < maxValue:
                        maxValue = x
                        maxIndex = counter
                    counter = counter + 1
                bestMove = moves[maxIndex]
                bestMoveReturn[:] = bestMove
                return
        return maxEval

    #Player is black so find advantage for white(By Maximizing)
    else:
        bValue = 100000
        minEval = 1000000
        storage = []
        for item in moves:
            newpos = pst.clone()
            makemove(newpos,item[0][0],item[0][1],item[1][0],item[1][1])

            eval = miniMax(newpos,depth-1,alpha,beta,-1,bestMoveReturn,False)

            if(root == True):
                storage.append(eval)
            minEval = min(eval,minEval)
            beta = min(beta,eval)
            if beta <= alpha:
                break
            if(root == True):
                #Find Max move
                counter = 0
                maxIndex = 0
                maxValue = storage[0]
                for x in storage:
                    if x > maxValue:
                        maxValue = x
                        maxIndex = counter
                    counter = counter + 1
                bestMove = moves[maxIndex]
                bestMoveReturn[:] = bestMove
                return
        return minEval        

    







#SCORING FUNCTIONS ========================================================================================================================

#Functions for scoring positions on multiple criteria (Score tables , possibility of checkmate ,etc)
def evaluate(pst):
    if isCHKmate(pst,'white'):
        return -20000
    if isCHKmate(pst,'black'):
        return 20000
 
    brd = pst.getbrd()
  
    flatbrd = [x for row in brd for x in row]
   
    c = Counter(flatbrd)
    Qw = c['Qw']
    Qb = c['Qb']
    Rw = c['Rw']
    Rb = c['Rb']
    Bw = c['Bw']
    Bb = c['Bb']
    Nw = c['Nw']
    Nb = c['Nb']
    Pw = c['Pw']
    Pb = c['Pb']
   
    wM = 9*Qw + 5*Rw + 3*Nw + 3*Bw + 1*Pw
    bM = 9*Qb + 5*Rb + 3*Nb + 3*Bb + 1*Pb
    movesAvailable = len(pst.gH())
    PHASE = 'opening'
    if movesAvailable>40 or (wM<14 and bM<14):
        PHASE = 'ending'
  
    Dw = dblPawns(brd,'white')
    Db = dblPawns(brd,'black')
    Sw = blkPawns(brd,'white')
    Sb = blkPawns(brd,'black')
    Iw = islPawns(brd,'white')
    Ib = islPawns(brd,'black')
  
    evaluation1 = 900*(Qw - Qb) + 500*(Rw - Rb) +330*(Bw-Bb )+320*(Nw - Nb) +100*(Pw - Pb) +-30*(Dw-Db + Sw-Sb + Iw- Ib )
    evaluation2 = UnitST(flatbrd,PHASE) 
    evaluation = evaluation1 + evaluation2
    return evaluation
    
def UnitST(flatbrd,PHASE):
    
    score = 0
    for i in range(64):
        if flatbrd[i]==0:       
            continue
       
        Unit = flatbrd[i][0]
        cl = flatbrd[i][1]
        sign = +1
    
        if cl=='b':
            i = (7-i/8)*8 + i%8
            sign = -1
        i = int(i)
        if Unit=='P':
            score += sign*pawn_values[i]
        elif Unit=='N':
            score+= sign*knight_values[i]
        elif Unit=='B':
            score+=sign*bishop_values[i]
        elif Unit=='R':
            score+=sign*rook_values[i]
        elif Unit=='Q':
            score+=sign*queen_values[i]
        elif Unit=='K':
            if PHASE == 'opening':
                score+=sign*king_values[i]
            else:
                score+=sign*king_endValues[i]
    return score  
    
def dblPawns(brd,cl):
    cl = cl[0]
    LPW = lf(brd,'P'+cl)
    repeats = 0
    temp = []
    for pawnpos in LPW:
        if pawnpos[0] in temp:
            repeats = repeats + 1
        else:
            temp.append(pawnpos[0])
    return repeats
    
def blkPawns(brd,cl):
    cl = cl[0]
    LPW = lf(brd,'P'+cl)
    blocked = 0 
    for pawnpos in LPW:
        if ((cl=='w' and isOCCby(brd,pawnpos[0],pawnpos[1]-1, 'black'))
            or (cl=='b' and isOCCby(brd,pawnpos[0],pawnpos[1]+1, 'white'))):
            blocked = blocked + 1
    return blocked

def islPawns(brd,cl):
    cl = cl[0]
    LPW = lf(brd,'P'+cl)
    
    xlist = [x for (x,y) in LPW]
    isolated = 0
    for x in xlist:
        if x!=0 and x!=7:
            if x-1 not in xlist and x+1 not in xlist:
                isolated+=1
        elif x==0 and 1 not in xlist:
            isolated+=1
        elif x==7 and 6 not in xlist:
            isolated+=1
    return isolated















#CHESS BOARD =================================================================================================================================

#Game board:
brd = [   ['Rb', 'Nb', 'Bb', 'Qb', 'Kb', 'Bb', 'Nb', 'Rb'],  #8
          ['Pb', 'Pb', 'Pb', 'Pb', 'Pb', 'Pb', 'Pb', 'Pb'],  #7
          [  0,    0,    0,    0,    0,    0,    0,    0],   #6
          [  0,    0,    0,    0,    0,    0,    0,    0],   #5
          [  0,    0,    0,    0,    0,    0,    0,    0],   #4
          [  0,    0,    0,    0,    0,    0,    0,    0],   #3
          ['Pw', 'Pw', 'Pw',  'Pw', 'Pw', 'Pw', 'Pw', 'Pw'], #2
          ['Rw', 'Nw', 'Bw',  'Qw', 'Kw', 'Bw', 'Nw', 'Rw'] ]#1
          # a      b     c     d     e     f     g     h

plr = 0 
cst_rgh = [[True, True],[True, True]]
E_P_T = -1 
HalfMoveClock = 0 
pst = GMPO(brd,plr,cst_rgh,E_P_T ,HalfMoveClock)

#Unit Square Tables
pawn_values = [  0,  0,  0,  0,  0,  0,  0,  0,
                50, 50, 50, 50, 50, 50, 50, 50,
                10, 10, 20, 30, 30, 20, 10, 10,
                 5,  5, 10, 25, 25, 10,  5,  5,
                 0,  0,  0, 20, 20,  0,  0,  0,
                 5, -5,-10,  0,  0,-10, -5,  5,
                 5, 10, 10,-20,-20, 10, 10,  5,
                 0,  0,  0,  0,  0,  0,  0,  0 ]
knight_values =[-50,-40,-30,-30,-30,-30,-40,-50,
                -40,-20,  0,  0,  0,  0,-20,-40,
                -30,  0, 10, 15, 15, 10,  0,-30,
                -30,  5, 15, 20, 20, 15,  5,-30,
                -30,  0, 15, 20, 20, 15,  0,-30,
                -30,  5, 10, 15, 15, 10,  5,-30,
                -40,-20,  0,  5,  5,  0,-20,-40,
                -50,-90,-30,-30,-30,-30,-90,-50 ]
bishop_values =[-20,-10,-10,-10,-10,-10,-10,-20,
                -10,  0,  0,  0,  0,  0,  0,-10,
                -10,  0,  5, 10, 10,  5,  0,-10,
                -10,  5,  5, 10, 10,  5,  5,-10,
                -10,  0, 10, 10, 10, 10,  0,-10,
                -10, 10, 10, 10, 10, 10, 10,-10,
                -10,  5,  0,  0,  0,  0,  5,-10,
                -20,-10,-90,-10,-10,-90,-10,-20 ]
rook_values = [  0,  0,  0,  0,  0,  0,  0,  0,
                 5, 10, 10, 10, 10, 10, 10,  5,
                -5,  0,  0,  0,  0,  0,  0, -5,
                -5,  0,  0,  0,  0,  0,  0, -5,
                -5,  0,  0,  0,  0,  0,  0, -5,
                -5,  0,  0,  0,  0,  0,  0, -5,
                -5,  0,  0,  0,  0,  0,  0, -5,
                 0,  0,  0,  5,  5,  0,  0,  0 ]
queen_values = [-20,-10,-10, -5, -5,-10,-10,-20,
                -10,  0,  0,  0,  0,  0,  0,-10,
                -10,  0,  5,  5,  5,  5,  0,-10,
                 -5,  0,  5,  5,  5,  5,  0, -5,
                  0,  0,  5,  5,  5,  5,  0, -5,
                -10,  5,  5,  5,  5,  5,  0,-10,
                -10,  0,  5,  0,  0,  0,  0,-10,
                -20,-10,-10, 70, -5,-10,-10,-20 ]
king_values =  [-30,-40,-40,-50,-50,-40,-40,-30,
                -30,-40,-40,-50,-50,-40,-40,-30,
                -30,-40,-40,-50,-50,-40,-40,-30,
                -30,-40,-40,-50,-50,-40,-40,-30,
                -20,-30,-30,-40,-40,-30,-30,-20,
                -10,-20,-20,-20,-20,-20,-20,-10,
                 20, 20,  0,  0,  0,  0, 20, 20,
                 20, 30, 10,  0,  0, 10, 30, 20 ]
king_endValues=[-50,-40,-30,-20,-20,-30,-40,-50,
                -30,-20,-10,  0,  0,-10,-20,-30,
                -30,-10, 20, 30, 30, 20,-10,-30,
                -30,-10, 30, 40, 40, 30,-10,-30,
                -30,-10, 30, 40, 40, 30,-10,-30,
                -30,-10, 20, 30, 30, 20,-10,-30,
                -30,-30,  0,  0,  0,  0,-30,-30,
                -50,-30,-30,-30,-30,-30,-30,-50 ]

LOWU,LOBU= crtUNITS(brd)
listofSds = []

openings = defaultdict(list)

isDraw = False 
chessEnded = False 
prevMove = [-1,-1,-1,-1] 
ax,ay=0,0
numm = 0











#PRINT ========================================================================================================================================

print("\n\t\t\tCHESS GAME\n")
print("* Player is White. ")
print("* Player chess pieces")
print("  + Pw : Pawn white")
print("  + Rw : Rook white")
print("  + Nw : Knight white")
print("  + Bw : Bishop white")
print("  + Qw : Queen white")
print("  + Kw : King white")
print("* Enter moves treating board as matrix with rows and columns from 0-7")
print("\n")
print("CHESS BOARD Cordinates:")
print("[ 00 ,  01 ,  02 ,  03 ,  04 ,  05 ,  06 ,  07 ]")
print("[ 10 ,  11 ,  12 ,  13 ,  14 ,  15 ,  16 ,  17 ]")
print("[ 20 ,  21 ,  22 ,  23 ,  24 ,  25 ,  26 ,  27 ]")
print("[ 30 ,  31 ,  32 ,  33 ,  34 ,  35 ,  36 ,  37 ]")
print("[ 40 ,  41 ,  42 ,  43 ,  44 ,  45 ,  46 ,  47 ]")
print("[ 50 ,  51 ,  52 ,  53 ,  54 ,  55 ,  56 ,  57 ]")
print("[ 60 ,  61 ,  62 ,  63 ,  64 ,  65 ,  66 ,  67 ]")
print("[ 70 ,  71 ,  72 ,  73 ,  74 ,  75 ,  76 ,  77 ]")












#GAME RUNNING ================================================================================================================================

while chessEnded == False:   
    counter = 0 
    r1check=0
    colll = 0
    print("\n\n\n\t\tCHESS BOARD")
    for x in brd:
        for r1 in range(0,8):
            if x[r1]==0:
                r1check=1
        if r1check==0:
            print(str(colll) + " " + str(x))
            colll = colll+1
        elif r1check==1:
            print( str(colll) + " [" , end="")
            colll+1
            for r2 in range(0,7):
                if x[r2]==0:
                    print("'  '," , end=" ")
                else:
                    print("'" + str(x[r2]) + "'," , end=" ")
            if x[7] == 0:
                print("'  '" , end="")
            else:
                print("'" + str(x[7]) + "'" , end="")
            print("]")
            colll = colll+1
            r1check=0
    print("    0  ,  1  ,  2  ,  3  ,  4  ,  5  ,  6  ,  7   ")

    #player move
    print("\nMake your MOVE ! ")
    row = input("Row Initial    :")
    column = input("Column Initial :")
    rowFinal = input("Row Final      :")
    columnFinal = input("Column Final   :")
    y = int(row)
    x = int(column)
    y2 = int(rowFinal)
    x2 = int(columnFinal)
    
    #Make sure player doesn't move the wrong colour
    while isOCCby(pst.brd,x,y,'white') != True:
        print("\nPlease pick a square with a white Unit")
        row = input("Row Initial    :")
        column = input("Column Initial :")
        y = int(row)
        x = int(column)
        rowFinal = input("Row Final      :")
        columnFinal = input("Column Final   :")
        y2 = int(rowFinal)
        x2 = int(columnFinal)
    
    #Make sure player doesn't move to a square with an existing Unit
    while isOCCby(pst.brd,x2,y2,'white') != False:
        print("\nPlease pick a valid square to move to")
        rowFinal = input("Row Final      :")
        columnFinal = input("Column Final   :")
        y2 = int(rowFinal)
        x2 = int(columnFinal)
    
    #Make sure player makes a valid move with the given piece
    possibilities = allM(pst,1)
    found = False

    while found != True:
        moveCount = 0
        for m in possibilities:
            if (m[0][0] == x) and (m[0][1]==y):
                moveCount = moveCount + 1
            if (m[0][0] == x) and (m[0][1]==y)and(m[1][0]==x2)and(m[1][1]==y2):
                found = True
        if (moveCount == 0):
            print("\nThis piece has no available moves. Please change your pick ! ")
            row = input("Row Initial    :")
            column = input("Column Initial :")
            y = int(row)
            x = int(column)
            rowFinal = input("Row Final      :")
            columnFinal = input("Column Final   :")
            y2 = int(rowFinal)
            x2 = int(columnFinal)
        elif(found == False):
            print("\nPlease pick a valid move for this piece !")
            rowFinal = input("RowFinal      :")
            columnFinal = input("ColumnFinal   :")
            y2 = int(rowFinal)
            x2 = int(columnFinal)
    
    makemove(pst,x,y,x2,y2)


    clsign=-1
    bestMoveReturn = []
    miniMax(pst,3,-1000000,1000000,clsign,bestMoveReturn,True)
    
    x = int(bestMoveReturn[0][0])
    y = int(bestMoveReturn[0][1])
    x2 = int(bestMoveReturn[1][0])
    y2 = int(bestMoveReturn[1][1])
    makemove(pst,x,y,x2,y2)    

    prevMove = [x,y,x2,y2]
    plr = pst.getplr()
    HMC = pst.getHMC()
    #pst.aTh(pst)
    if HMC>=100 or isSTL(pst):
        isDraw = True
        chessEnded = True
    if isCHKmate(pst,'white'):
        winner = 'b'
        chessEnded = True
    if isCHKmate(pst,'black'):
        winner = 'w'
        chessEnded = True
    










#RESULT ======================================================================================================================================

if isDraw == True:
    print("ITS A DRAW !")
elif winner == 'b':
    print("CHECKMATE ! BLACK WINS !\n\n")
elif winner == 'w':
    print("CHECKMATE! WHITE WINS !\n\n")
