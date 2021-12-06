import numpy as np
import math
import re

#Parser to Read mps files
def read_mps(path):
  mode = ""
  name = None
  cost = None
  restrictions_names = []
  types = []
  var_names = []
  var_types = []
  A = np.matrix([[]])
  c = np.array([])
  rhs_names = []
  rhs = {}
  bnd_names = []
  bnd = {}
  k = False

  with open(path, "r") as reader:
   for l in reader:
      l = re.split(" |\t", l)
      l = [x.strip() for x in l]
      l = list(filter(None, l))

      if l[0] == "ENDATA":
        break
      if l[0] == "*":
        continue
      if l[0] == "NAME":
        name = l[1]
      elif l[0] in ["ROWS", "COLUMNS"]:
        mode = l[0]
      elif l[0] == "RHS" and len(l) <= 2:
        if len(l) > 1:
          rhs_names.append(l[1])
          rhs[l[1]] = np.zeros(len(restrictions_names))
          mode = "RHS_NAME"
        else:
          print('RHS_NO_NAME')
          mode = "RHS_NO_NAME"
      elif l[0] == "BOUNDS" and len(l) <= 2:
        if len(l) > 1:
          bnd_names.append(l[1])
          bnd[l[1]] = {"LO": np.zeros(len(var_names)), "UP": np.repeat(math.inf, len(var_names))}
          mode = "BOUNDS_NAME"
        else:
           mode =  "BOUNDS_NO_NAME"
      elif mode == "ROWS":
        if l[0] == "N":
           cost = l[1]
        else:
           types.append(l[0])
           restrictions_names.append(l[1])
      elif mode == "COLUMNS":
        if len(l) > 1 and l[1] == "'MARKER'":
           if l[2] == "'INTORG'":
              k = True
           elif l[2] == "'INTEND'":
              k = False
           continue
        try:
           i = var_names.index(l[0])
        except:
           if A.shape[1] == 0:
               A = np.zeros((len(restrictions_names), 1))
           else:
               A = np.concatenate((A, np.zeros((len(restrictions_names), 1))), axis = 1)
           var_names.append(l[0])
           var_types.append(k * 'integral' + (not k) * 'continuous')
           c = np.append(c, 0)
           i = -1
        j = 1
        while j < len(l) - 1:
           if l[j] == cost:
             c[i] = float(l[j + 1])
           else:
             A[restrictions_names.index(l[j]), i] = float(l[j + 1])
           j = j + 2
      elif mode == "RHS_NO_NAME":
       print('RHS_N0_NAME2')
       try:
             i = rhs_names.index(l[0])
       except:
             rhs_names.append(l[0])
             rhs[l[0]] = np.zeros(len(restrictions_names))
             i = -1
       rhs[l[0]][restrictions_names.index(l[1])] = float(l[2])
      elif mode == "BOUNDS_NAME":
       if l[1] != bnd_names[-1]:
            raise Exception("Other BOUNDS name was given even though name was set after BOUNDS tag.")
       if l[0] in ["LO", "UP"]:
            bnd[l[1]][l[0]][var_names.index(l[2])] = float(l[3])
       elif l[0] == "FX":
            bnd[l[1]]["LO"][var_names.index(l[2])] = float(l[3])
            bnd[l[1]]["UP"][var_names.index(l[2])] = float(l[3])
       elif l[0] == "FR":
            bnd[l[1]]["LO"][var_names.index(l[2])] = -math.inf
      elif mode ==  "BOUNDS_NO_NAME":
       try:
            i = bnd_names.index(l[1])
       except:
            bnd_names.append(l[1])
            bnd[l[1]] = {"LO": np.zeros(len(var_names)), "UP": np.repeat(math.inf, len(var_names))}
            i = -1
       if l[0] in ["LO", "UP"]:
            bnd[l[1]][l[0]][var_names.index(l[2])] = float(l[3])
       elif l[0] == "FX":
            bnd[l[1]]["LO"][var_names.index(l[2])] = float(l[3])
            bnd[l[1]]["UP"][var_names.index(l[2])] = float(l[3])
       elif l[0] == "FR":
            bnd[l[1]]["LO"][var_names.index(l[2])] = -math.inf
  dim_A=np.shape(A)
  number_of_restrictions=dim_A[0]
  Eqin = [None] * (number_of_restrictions)
  for t in range(0,number_of_restrictions):
        if types[t] == 'L':
            Eqin[t]=-1
        elif types[t] == 'G':
            Eqin[t]=1
        else:
            Eqin[t]=0
  return name, cost, restrictions_names, var_names, var_types, types, c, A, rhs_names, rhs, bnd_names, bnd, Eqin

# Add your path!
ap=read_mps('*YOUR PATH HERE*/sdata1_100x100.mps')


A=ap[7]#the matrix which contains the coefficients of each variable in the set of constraints
print("Τhe matrix which contains the coefficients of each variable in the set of constraints,is the following:")
print(A)
import numpy 
dim_A=np.shape(A)#dimensions of A matrix
print('The number of constraints which exist in the problem is:',dim_A[0],"\n",'The number of variables in the problem is:',dim_A[1])
met=ap[3]#the variables of the problem
pl_met=dim_A[1]#number of variables
ty_per=ap[5]#types of constraints(LO,UP,E)
pl_per=len(ty_per)#number of constraints
names_res=ap[2]#names of constraints
name_rhs=list(ap[9].keys())[0]
b=ap[9][name_rhs]#nonzero right-hand side values of the constraints
print("nonzero right-hand side values of the constraints are the following:",b)
c=ap[6]#the coefficients of cost function
print("The coefficients of the cost function are:",c)
per=ap[2]
if not not ap[10]:#in case there are bounds for some variables
  name_bounds=list(ap[11].keys())[0]
  type(list(ap[11].values())[0])
  val_lo_bnd=list(ap[11].values())[0]['LO']
  val_up_bnd=list(ap[11].values())[0]['UP']
eqin=ap[12]
print('the kinds of the constraints are the following(-1 is for <=, 1 is for>=, 0 is for = ',eqin)

###################################################################################################################################################
#Exterior point algorithm

print(A.shape) 
w, h = len(A), len(A);
bd = np.zeros((w, h))

B = np.array([])
N = np.array([])
P = np.array([])
Q = np.array([])
L = np.array([])
S0 = np.array([])
dB = np.array([])
Sp = np.array([])
eis1 = np.array([])
eis2 = np.array([])
Sq = np.array([])
hj = np.array([])
Wt = np.array([])
SnP = np.array([])


bcols=0

for i in range(len(bd)):                #check inequalities, where eqin -1 replace with +1 to the bd matrix. Where 1 replace with -1. 
    for j in range(0, i+1):              
      if eqin[i] == -1:
          bd[i,bcols] = 1
      if eqin[i] == 1:
          bd[i,bcols] = -1
    bcols +=1 
    
Aa = np.concatenate((A,bd), axis=1)     # New A ( A + bd)  

for i in range(0, len(A)):              # Split Ct to 2 matrices B and N. -To B put the indicator Ct to which the values of vector C do not belong
    N= np.append(N,i)                   
for i in range(100,len(bd)*2):          # -To N put the indicator Ct where belong the values of matrix C. The values of colmun's of matrix A 
    B= np.append(B,i)                   # correspont to indecators of matrix A and the rows of matrix B to the values of matrix B

for i in range(0, len(bd)):             # Wt= loose variables to the matrix Ct
 if eqin[i] == -1 or eqin[i] == 1:
      Wt= np.append(Wt,0)

Ct = np.concatenate((c,Wt), axis=0)     # Vector Ct 

Xb = np.dot(bd,b)                       # Matrix Xb ( bd * b)
tempsn= np.dot(A,Wt)
Sn = np.subtract(c, tempsn)             # Matrix Sn (c - (A * Wt))

for i in range(0, len(Sn)):             # Create vectors P and Q, Where I find value < 0 into vector Sn, I save the indicator to the vector P
    if Sn[i] < 0:                       # Where I find value >= 0 into vector Sn, I save the indicator to the vector Q
        P = np.append(P,i)  
    else:
        Q = np.append(Q,i)

for i in range(0, len(P)):              # Create vector L
    L = np.append(L,1)


for i,x in enumerate(Sn):               
    if x < 0 :                          # i = indeces , x= values
        SnP = np.append(SnP,x)

S0 = -1
while S0 != 0:
    for i in range(0, len(SnP)):            # Multiply matrix Sn with vectro L to calculate S0
      S0 += L[i]*SnP[i]                     # (If S0=0 then optimal solution)
    print(S0)
    
    list1 = P.tolist()                      # Take every column of A where refered "P"
    list1 = list(map(int, list1))           
    PcolsA= Aa[:, list1]
    
    
    for i in range(0, len(P)):              # Calculate vector dB. Multiply Vector L with h and sum it to find vector dB
       h += np.dot(bd,PcolsA[:,i])          
    dB = np.dot(-(L)[i],h)                  
                                            
            
                
        
    # step 2.1  -Select which variable will export from vectro B
    
    
    negDBvalues = 0                         
    for i in range(0, len(dB)):
       if dB[i] < 0:
        negDBvalues += 1
  
    #while S0 !=0 and not P:                 # Start the loop. Eliminated when find the optimal solution. ( S0 = 0 and P= empty list) (I have it as commend because)
                                             # As mentioned the code is not correct!
           
    exer = np.array([])
    min_exer = [0]                                  # Find out the minimum value and his indicator  
    for i in range(0, negDBvalues):                 
       temp2 = Xb[i] / -(dB[i])                     
       exer = np.append(exer,temp2 )                     
        
    min_exer = np.amin(exer)                        
    r = np.where(exer == np.amin(exer))             
    r = int(r[0])
     
    for i,x in enumerate(B):                     
        if i== r :                              
            K = int(x)                          
       
    #vima 2.2                                # Select incoming variable
    Br =np.array([])
    for i,x in enumerate(bd):               # Take the values from Bd from the rows that r is refered to and save the values to matrix Br
        if i == r :                         
            Br = np.append(Br,x)
        
    for i in range(0, len(P)):              # Calculate to Hrp. 
        HrP = np.dot(Br,PcolsA)             
           
    list3 = Q.tolist()                      
    list3 = list(map(int, list3))           
    QcolsA= Aa[:, list3]
    
    for i in range(0, len(Aa)):              # Same as HrP
        HrQ = np.dot(Br,QcolsA)              #QcolsA, Br
            
    Sp= SnP                                  # Create Sp vector. Contains all negative values from Sn.
    
    
    eis1 = np.array([])
    min_eis1 = [0]                                    
    for i in range(0, len(Sp)):                       
       if HrP[i] != 0:                               
            temp3 = -(Sp[i]) / HrP[i]                 
            eis1 = np.append(eis1, temp3 )                  
       else:
            eis1 = np.append(eis1, 1000)     # When operations were performed with 0 the result was inf. I put a large number so as not to affect it 
                                                                                   
    min_eis1 = np.amin(eis1)                 
    t1 = np.where(eis1 == np.amin(eis1))     
    t1 = int(t1[0])
    
    SnQ = np.array([]) 
    for i,x in enumerate(Sn):                
        if x >= 0 :                         
            SnQ = np.append(SnQ,x)          
        
    eis2 = np.array([])
    min_eis2 = [0]                                    
    for i in range(0, len(Q)):
        if HrQ[i] != 0:                               
            temp4 = -(SnQ[i]) / HrQ[i]                
            eis2 = np.append(eis2, temp4 )                
        else:
            eis2 = np.append(eis2, 1000)              
                                                      
    min_eis2 = np.amin(eis2)                          
    t2 = np.where(eis2 == np.amin(eis2))              
    t2 = int(t2[0])    
     
           
    l1=0
    for i in range(0, len(Sp)):             # Indicator of potential incoming variable from matrix Sp
        if t1 == i:
              l1 = t1
    
    
    p = np.array([]) 
    for i,x in enumerate(P):                 
        if l1 == i :                        
            p = np.append(p,x)              
    
     
    l2=0
    for i in range(0, len(SnQ)):            # Indicator of potential incoming variable from matrix Sp
        if t2 == i:
              l2 = t2
    
    
    q = np.array([]) 
    for i,x in enumerate(Q):                
        if l2 == i :                        
            q = np.append(q,x) 
    pp= p
    qq= q  
    if min_eis1 <= min_eis2:                # If Thita1 <= Thita2 ( eis1 <= eis2) then the incoming value is the value from indicator p (l=p). Else q(l=q) 
      l = pp[0]                             
    else:                                   
      l = qq[0]
      
      
    
    # Step 2.3 pivoting
        
    temp = K 
    
    q = np.array([]) 
    for i,x in enumerate(B):                # Replace the value of indicator K and Vector B with value l
        if i == r :                         
            B[i] = l                                
    #print(B[28]) 
    #print("diktis t1:",t1,"    timi pp: ",pp[0])                                
    #print("diktis t2:",t2,"    timi qq: ",qq[0])
    if l == pp[0]:
     for i,x in enumerate(P):
       if i == t1:                              
        Q = np.append(Q, temp)                  
        P = np.delete(P,t1)   
    elif l==qq[0]:
      for i,x in enumerate(Q):
          if i==t2:
              P = np.append(P, temp)
              Q = np.delete(Q,t2)
              
              
    if l == pp[0]:                          
      for i,x in enumerate(B):              # If eqin == -1 and -1 for eqin == 1.          
        if i == t1 :                         
         if eqin[t1] == -1 :
            b[t1] = 1
         else:
            b[t1] = -1
    elif l==qq[0]:
     for i,x in enumerate(B):                
        if i == t2 :                         
         if eqin[t2] == -1 :
                  b[t2] = 1
         else:
                  b[t2] = -1
    
    
    antistrofosB = np.linalg.solve(A, b)         # New inverse of B
    
    list1 = P.tolist()                      
    list1 = list(map(int, list1))           
    PcolsA= Aa[:, list1]                    
    
    
    for i in range(0, len(P)):                  # Same commands as previous
       hj += np.dot(antistrofosB, PcolsA[:,i])  
    
    L=[1]   
    for i in range(0, len(P)-1):               
        L = np.append(L,1)     
        
    for i in range(0, len(L)):              
        dB = np.dot(-(L)[i],h)                