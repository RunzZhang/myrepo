import pickle
import matplotlib.pyplot, numpy
el32 = 2**(-15)
with open('C:\\Users\\24230\\Downloads\\MFC1008_scaling.pkl', 'rb') as f:
    data = pickle.load(f)
print(data, "\n",el32)
x = [706,1521,4797,8065,11333,14612,16631]
x_dig = []
for i in x:
    x_dig.append(i*el32)
y = [0.5,1,3,5,7,9,11]
m,b = numpy.polyfit(x_dig, y, 1)
print(m,b)
