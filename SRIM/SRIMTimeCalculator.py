import matplotlib.pyplot as plt
import numpy as np
import os, pickle
import sys
import pandas as pd
from statistics import stdev
pd.options.display.max_seq_items = None
#pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)

file_name = 'EXYZArgon0.125keV.txt'
file_name_edit = file_name[0:-4] + '_edit.txt'
skip = 0
Data = 0
avo = 6.0221408e+23
amu = 0 
KE = 0
IonData = ''
file_lines = []

with open(file_name, 'r') as file: # Searches for the row in which data starts and
    done = 0                       # looks for ion data in the file (we want the atomic mass)
    for num, line in enumerate(file, 1):
        if '0000001' in line and done == 0:
            #print('found at line:', num)
            skip = num - 1
            done = 1
        if 'Ion Data' in line:
            Data = num
            
            
with open(file_name, 'r') as file:  # Adds a space after first column to allow
    n = 0                           # pandas to use double-space as separator    
    for line in file:
        if n >= skip:
            file_lines.append(''.join([line.strip()[:7], ' ', line.strip()[7:], '\n'])) 
        else:
            #file_lines.append(''.join([line.strip(), '\n']))
            file_lines.append(''.join([line]))
        n += 1
        
with open(file_name, 'r') as file: # Adds the ion data string to a variable
    n = 0
    for line in file:
        if n == Data:
            IonData = line.rstrip()
            IonData = (' '.join(IonData.split())).split()
        n += 1

amu = float(IonData[1]) # Set the atomic mass variable
mass = amu/avo * 0.001 # calculate mass of 1 atom in kg
#masskev = amu/avo * 0.001 * 5.6095886e32 # mass in keV
print(amu)
print(mass)
    
#with open(file_name_edit, 'w') as file: # Write the editied file for pandas to open
 #   file.writelines(file_lines)
    

df = pd.read_csv(file_name, sep='  | ', skiprows=skip, header=None, engine=('python'))
#print(df.head(5))
df.rename( columns={0 :'Ion Number'}, inplace=True )
df.rename( columns={1 :'Energy (keV)'}, inplace=True )
df.rename( columns={2 :'x (A)'}, inplace=True )
df.rename( columns={3 :'y (A)'}, inplace=True )
df.rename( columns={4 :'z (A)'}, inplace=True )
df.rename( columns={5 :'Electronic Stop (eV/A)'}, inplace=True )
df.rename( columns={6 :'Energy Lost to Last Recoil (eV)'}, inplace=True )

Posx = df['x (A)'].tolist()
Posy = df['y (A)'].tolist()
Posz = df['z (A)'].tolist()
Energy = df['Energy (keV)'].tolist()
#print(Energy)
Energy = [1.60218e-16*x for x in Energy] # Writing all energy in Joules
ElStop = df['Electronic Stop (eV/A)'].tolist()
IonNum = df['Ion Number'].tolist()
Events = max(IonNum)
KEInitial = (df.iloc[0][1]) # Set kinetic energy variable (in keV)
#print(KE)

displacement = [] # Want to record the displacement in each step (in meters)

for i in range(Events): # Calculate displacement and record it in a 2D array which is
    temp = []           # indexed by the event number
    for j in range(len(Posx)):
        if (IonNum[j] == i + 1):
            #if Posx[j] == 0:
            #    temp.append(0)
            if Posx[j] !=0:
                temp.append((1e-10)*(((Posx[j]-Posx[j-1])**2+(Posy[j]-Posy[j-1])**2+(Posz[j]-Posz[j-1])**2)**0.5))
    displacement.append(temp)
     
#print(displacement[slice(3)])
#print(len(displacement[1]))


velocity_temp = [(2*x/mass)**0.5 for x in Energy] # Want velocity from mass and energy for each step 
#print(velocity_temp)
velocity = []
finalpos = [] # Want to record the final position

for i in range(Events): # Write velocity in a nicer 2D list indexed by the event number
    temp = []           
    for j in range(len(Posx)):
        if (IonNum[j] == i + 1):
            if ElStop[j] != 0: 
                temp.append(velocity_temp[j])
            if ElStop[j] == 0:
                finalpos.append([Posx[j], Posy[j], Posz[j]])
    velocity.append(temp)

del displacement[-1] #Last event sometimes is not complete so remove it
del velocity[-1]

finaldisplacement = [(1e-10)*(((finalpos[i][0])**2+(finalpos[i][1])**2+(finalpos[i][2])**2)**0.5) for i in range(len(finalpos))]

#for i in range(len(displacement)):
#    print(len(displacement[i]), len(velocity[i]))
print(len(velocity_temp))
print(len(velocity))
print(len(displacement))

time = [[displacement[j][i]/velocity[j][i] for i in range(len(velocity[j]))] for j in range(len(velocity))] # Want to calculate total time for each event 
#print(time[0])

timetot = [sum(time[i]) for i in range(len(time))] # Sums the times at each step for each event
#print(timetot)

timeaverage = 0 
for i in range(len(timetot)): # Calculate average time over all events
    timeaverage += timetot[i]
timeaverage = timeaverage/Events

finaldisplacementaverage = 0
for i in range(len(finaldisplacement)): # Calculate average final displacement over all events
    finaldisplacementaverage += finaldisplacement[i]
finaldisplacementaverage = finaldisplacementaverage/Events

totaldisplacement = [sum(displacement[i]) for i in range(len(displacement))]
#print(totaldisplacement)

totaldisplacementaverage = 0
for i in range(len(totaldisplacement)): # Calculate average total displacement over all events
    totaldisplacementaverage += totaldisplacement[i]
totaldisplacementaverage = totaldisplacementaverage/Events

print('The average time is:', timeaverage, 's', 'with a stdev of:', stdev(timetot), 's')
print('The average final distance is:', finaldisplacementaverage, 'm', 'with a stdev of:', stdev(finaldisplacement), 'm')
print('The average total displacement is:', totaldisplacementaverage, 'm', 'with a stdev of:', stdev(totaldisplacement), 'm')

fig, ax = plt.subplots(1,2)

def scale(list1, x):
    temp = [x*y for y in list1]
    return temp
    
ax[0].hist(scale(timetot, 1e+15), 25, histtype = "step") 
ax[1].hist(scale(finaldisplacement, 1e+9), 25, histtype = "step") 
ax[0].set_title("Total Time Dist.")
fig.text(0.45, 0.97, "E= "+str(KEInitial)+"keV")
fig.text(0.30, 0.04, 'Total Time (fs)', ha='center', va='center')
fig.text(0.72, 0.04, 'Final Position (nm)', ha='center', va='center')
fig.text(0.06, 0.5, 'frequency', ha='center', va='center', rotation='vertical')
ax[1].set_title("Final Pos. Dist.")
plt.show()
print(max(timetot), np.argmax(timetot), max(finaldisplacement), finaldisplacement[np.argmax(finaldisplacement)])
#print(sorted(finaldisplacement))
print((1e-10)*(((finalpos[400][0])**2+(finalpos[i][1])**2+(finalpos[400][2])**2)**0.5))
print(finalpos[400][0], finalpos[400][1], finalpos[400][2])






