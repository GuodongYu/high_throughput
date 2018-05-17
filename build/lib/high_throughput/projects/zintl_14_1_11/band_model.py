import matplotlib
matplotlib.use('Tkagg')
from matplotlib import pyplot as plt
vbmx = 0; cbmx = 0
vbmy = 0; cbmy = 0.5
vbm2x = 2; vbm2y = -0.3
def cb(x):
    return 0.1*(x-cbmx)**2 + cbmy

def vb(x):
    return -0.3*(x-vbmx)**2 + vbmy

def vb2(x):
    return -0.02*(x-vbm2x)**2 + vbm2y

xs = [i*0.1 for i in range(-50,50)]
plt.plot(xs, [cb(i) for i in xs],linewidth=5,color='black')
plt.plot(xs, [vb(i) for i in xs],linewidth=5,color='black')
plt.plot(xs, [vb2(i) for i in xs],linewidth=5,color='black')
plt.plot([-1,3],[vbmy,vbmy],linestyle='dashed',linewidth=3,color='red')
plt.plot([-1,1],[cbmy,cbmy],linestyle='dashed',linewidth=3,color='red')
plt.plot([1,3],[vbm2y,vbm2y],linestyle='dashed',linewidth=3,color='red')
plt.text(0.0,-0.1,'VBM',fontsize=15,horizontalalignment='center',verticalalignment='center')
plt.text(0.0,0.55,'CBM',fontsize=15,horizontalalignment='center',verticalalignment='center')
plt.text(2,-0.35,'VBM2',fontsize=15,horizontalalignment='center',verticalalignment='center')
plt.text(0.0,0.25,'Gap',fontsize=15)
plt.text(2,-0.15,'D12',fontsize=15)
plt.annotate('',xy=(vbmx,vbmy),xytext=(cbmx,cbmy),arrowprops=dict(arrowstyle="<->",color='red',linestyle='solid',linewidth=2))
plt.annotate('',xy=(vbm2x,vbm2y),xytext=(vbm2x,vbmy),arrowprops=dict(arrowstyle="<->",color='red',linestyle='solid',linewidth=2))
#plt.annotate('',xy=(cbmx,cbmy),xytext=(cbmx,cbmy/2),arrowprops=dict(arrowstyle="->"))
plt.xlim(-3,5)
plt.ylim(-0.4,0.8)
plt.axis('off')
plt.savefig('band_model.svg',fmt='svg')
