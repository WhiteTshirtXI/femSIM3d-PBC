
from pylab import plot,axis,xlabel,ylabel,title,legend,savefig,yscale,show

it=2500

F1=[];G1=[];
H1=[];time=[];
dt=0
for i in range(0,it):
 inputFile = 'vkError.' + str(i)
 # reading file
 file = open('../sim/'+inputFile,"r")
 lastLine=len(file.readlines())-1 # numero de linhas valido para plot
 file = open('../sim/'+inputFile,"r")
 file.readline()
 for i in range(0,lastLine):
  aux = map(float,file.readline().split())
  if i == lastLine-14:
   F1.append(aux[1])
   G1.append(aux[2])
   H1.append(aux[3])
   dt+=5.1250850251e-03
   time.append(dt)

H2=[];F2=[];G2=[]
for i in range(0,it):
 inputFile = 'vkError.' + str(i)
 # reading file
 file = open('../sim/'+inputFile,"r")
 lastLine=len(file.readlines())-1 # numero de linhas valido para plot
 file = open('../sim/'+inputFile,"r")
 file.readline()
 for i in range(0,lastLine):
  aux = map(float,file.readline().split())
  if i == lastLine-13:
   F2.append(aux[1])
   G2.append(aux[2])
   H2.append(aux[3])

# plotting file
plot(time,F1,time,G1,time,H1)
plot(time,F2,time,G2,time,H2)
#axis([0.0,10.0,-0.1,1.11])
yscale('log')
xlabel('time')
ylabel('H')
legend(('F1','G1','H1','F2','G2','H2'))
show()

