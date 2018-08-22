import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

def velocitymaps(data, sys):
    gas = data[data['PTG'] == 1]
    stars = data[data['PTG'] == 2]
    plt.figure(figsize = (10,4))
    plt.subplot(121)
    plt.scatter(gas['X'], gas['Y'], c=gas['V']-sys, cmap = 'bwr')
    plt.title('Gas')
    plt.colorbar()
    plt.subplot(122)
    plt.scatter(stars['X'], stars['Y'], c=stars['V']-sys, cmap = 'bwr')
    plt.title('Stars')
    plt.colorbar()
    plt.tight_layout()    
    plt.show()

def rotationcurve(data, ulim = 300, color = 'b', label = ''):
    plt.plot(data['R'], (data['V']-sys)/np.cos(data['PHI']*np.pi/180), color+',')
    uniquer = np.unique(data['R'],return_index = True)[1]
    uniquei = np.zeros(len(data), dtype = bool)
    uniquei[uniquer] = 1
    plt.plot(data['R'][uniquei], (data['VMOD'][uniquei]-sys)/np.cos(data['PHI'][uniquei]*np.pi/180), 
            color+'.', label = label)
    plt.ylim(0,ulim)
    
def readheader(file):
    header = open(file).readlines()[11:24]
    readrange = (header[0].find('PARAMETER') + 9, header[0].find(' F'))
    header = [i[readrange[0]:readrange[1]].split(' ') for i in header]
    return np.asarray([[i for i in header[j] if i] for j in range(len(header))][1:], dtype = float)

def visualize(file):
    names = ('PTG', 'FIB', 'X', 'Y', 'R', 'PHI', 'V', 'Ve', 'Ve_ADJ', 'VMOD', 'CHI', 'WGT', 'FLG', 'O')
    data = np.genfromtxt(file, comments = '#', names = names)
    gas = data[data['PTG'] == 1]
    stars = data[data['PTG'] == 2]
    header = readheader('test.out')
    inc = header[0,0]
    sys = header[2,0]
    velocitymaps(data, sys)
    rotationcurve(stars, label = 'Stars')
    rotationcurve(gas, color = 'r', label = 'Gas')
    plt.legend()

if __name__=='__main__':
    visualize('/home/brian/Desktop/Code/AD/test.out')
    visualize('/home/brian/Downloads/kbw_release/test_c.out')
    visualize('/home/brian/Downloads/kbw_release/testic.out')
    '''


    # In[28]:


    #18: untie systemic velocities, n, 1,2
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[29]:


    #6: delete based on error limit, 300
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[30]:


    #7: delete based on chisq limit, 1000
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[31]:


    #14: tie pointings, 1,2 (why tie back together?)
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[32]:


    #8: undelete based on chisq limit, 100 (what does this mean? why?)
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[33]:


    #8: undelete based on chisq limit, 100 (why again?)
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[34]:


    #8: undelete based on chisq limit, 100 (why?)
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[35]:


    #7: delete based on chisq limit, 25 (isn't this just undoing previous steps?)
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[36]:


    #26,n,1,2: untie stochastic errors
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[37]:


    #10: adjust chisq based off error dist, 1
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[38]:


    #10: adjust chisq based off error dist, 1 (why again?)
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[39]:


    #8: undelete on chisq, 36
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[40]:


    #10: adjust chisq based on error dist, 1
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[41]:


    #3: delete points based on gaussian growth curve, 1 (what does this mean? Why does it take so long?)
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[42]:


    #3,1 again
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[43]:


    #3,1 again
    visualize('/home/brian/Downloads/kbw_release/example/test.out')


    # In[45]:



    visualize('/home/brian/Downloads/kbw_release/example/test.out')
    '''
