
from pylab import plot,axis,xlabel,ylabel,title,legend,savefig,show


input = 'vk3.8'

# reading file
Z=[];F=[];G=[];H=[];C=[];
file = open(input,"r")
lastLine=len(file.readlines())-3 # numero de linhas valido para plot
print lastLine
file = open(input,"r")
file.seek(13) # movendo ponteiro para 2 linha do arquivo de dados
for i in range(0,lastLine):
 aux = map(float,file.readline().split())
 Z.append(aux[0])
 F.append(aux[1])
 G.append(aux[2])
 H.append(aux[3])
 C.append(aux[4])

Za=[];Fa=[];Ga=[];Ha=[];Ca=[];
inputAnalytic = '../../../db/baseState/nuC/Sc2000/analiticoNuC.dat'
fileAnalytic = open(inputAnalytic,"r")
lastLine=len(fileAnalytic.readlines()) # numero de linhas valido para plot
for i in range(0,lastLine):
 aux = map(float,fileAnalytic.readline().split())
 Za.append(aux[0])
 Fa.append(aux[1])
 Ga.append(aux[2])
 Ha.append(aux[3])
 Ca.append(aux[4])

# plotting file
plot(Za,Fa,Za,Ga,Za,Ha,Za,Ca)
axis([0.0,10.0,-0.1,1.01])
xlabel('Z')
ylabel('r')
title('Velocities and Concentration profiles of rotating disk')
legend(('F','G','H','C'))
show()

