from activeEL import ActiveEL
import matplotlib.pyplot as plt
for i in range(1,5):
    #error_array_1 = ActiveEL(i,8,O_eta = 1)
    #plt.plot(error_array_1)
    #plt.show()
    error_array_2 = ActiveEL(i,nStep = 256, nex = 16,T = 1, O_eta = -1)
    error_array_1 = ActiveEL(i,nStep = 256, nex = 16,T = 1, O_eta = 1)

    

    
