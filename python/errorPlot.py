
from pylab import plot,axis,xlabel,ylabel,title,legend,savefig,show


it = 2042

F1=[];G1=[];H1=[];time=[];
dt=0
for i in range(0,it):
 inputFile = 'vk12.' + str(i)
 # reading file
 file = open('../sim/'+inputFile,"r")
 lastLine=len(file.readlines())-3 # numero de linhas valido para plot
 file = open('../sim/'+inputFile,"r")
 file.readline()
 for i in range(0,lastLine):
  aux = map(float,file.readline().split())
  if i == lastLine-1:
   F1.append(aux[1])
   G1.append(aux[2])
   H1.append(aux[3]-0.866)
   dt+=5.1250850251e-03
   time.append(dt)

H2=[]
F2=[];G2=[];
for i in range(0,it):
 inputFile = 'vk13.' + str(i)
 # reading file
 file = open('../sim/'+inputFile,"r")
 lastLine=len(file.readlines())-3 # numero de linhas valido para plot
 file = open('../sim/'+inputFile,"r")
 file.readline()
 for i in range(0,lastLine):
  aux = map(float,file.readline().split())
  if i == lastLine-1:
   F2.append(aux[1])
   G2.append(aux[2])
   H2.append(aux[3]-0.866)

H3=[]
F3=[];G3=[];
for i in range(0,it):
 inputFile = 'vk8.' + str(i)
 # reading file
 file = open('../sim/'+inputFile,"r")
 lastLine=len(file.readlines())-3 # numero de linhas valido para plot
 file = open('../sim/'+inputFile,"r")
 file.readline()
 for i in range(0,lastLine):
  aux = map(float,file.readline().split())
  if i == lastLine-1:
   F3.append(aux[1])
   G3.append(aux[2])
   H3.append(aux[3]-0.866)

# plotting file
plot(time,F1,time,G1,time,H1)
#axis([0.0,10.0,-0.1,1.11])
xlabel('time')
ylabel('H')
legend(('R=80','R=86'))
show()

