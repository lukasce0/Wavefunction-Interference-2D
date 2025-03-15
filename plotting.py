
import numpy as np
import matplotlib.pyplot as plt



"""
Plotting the total probability (normalization)
"""



def readFilePTotal(filename):
    
    probs = []

    with open(filename) as file:
        for line in file:
            values = line.split()
            prob_line = 0
            
            for number in values:
                prob_line += float(number)
            
            probs.append(prob_line)
        
    return probs





P_tot_nowall = readFilePTotal("results_probability_no_wall.txt")
P_tot_wall = readFilePTotal("results_probability_wall.txt")
t = np.linspace(0, 0.08, len(P_tot_nowall))


plt.rcParams.update({'font.size': 18})

plt.figure(figsize=(8,6))
plt.xlabel(r"Time")
plt.ylabel(r"Sum of probability")
plt.plot(t, P_tot_nowall, label="Without wall")
plt.plot(t, P_tot_wall, label="With wall")
plt.legend()
plt.savefig("normlization_error.pdf")







"""
Plottig the probability of a double slit in room and real/imaginary components of wavefunction:
"""


dim = 199 #M-2




def readFilePLoc(filename):
    prob_t = []
    prob_all = []

    with open(filename) as file:
    
        for line in file:
            values = line.split()
        
            count = 0
            
            for i in range(dim):
                prob_t.append([])
            
                for j in range(dim):
                    prob_t[-1].append(float(values[count]))
                    count += 1

        
            prob_all.append(prob_t)
            prob_t = []
    
    return prob_all



def plotContour(metrix, time_index, fig_name):
    plt.figure(figsize=(8,6))
    plt.imshow(metrix[time_index])
    plt.xlabel("position x")
    plt.ylabel("position y")
    plt.xticks(dim*np.array([0,0.2,0.4,0.6,0.8,1]),[0,0.2,0.4,0.6,0.8,1])
    plt.yticks(dim*np.array([0,0.2,0.4,0.6,0.8]),[1,0.8,0.6,0.4,0.2])
    plt.colorbar()
    plt.savefig(fig_name)
        



prob_double = readFilePLoc("results_probability_contur.txt")
N_times = len(prob_double)
plotContour(prob_double, 0, "contur_P_2_t0.pdf")
plotContour(prob_double, N_times//2, "contur_P_2_12T.pdf") #time 0.001
plotContour(prob_double, -1, "contur_P_2_T.pdf") #time 0.002


prob_double_real = readFilePLoc("results_probability_contur_real.txt")
plotContour(prob_double_real, 0, "contur_rel_2_t0.pdf")
plotContour(prob_double_real, N_times//2, "contur_rel_2_12T.pdf") #time 0.001
plotContour(prob_double_real, -1, "contur_rel_2_T.pdf") #time 0.002


prob_double_imag = readFilePLoc("results_probability_contur_imag.txt")
plotContour(prob_double_imag, 0, "contur_imag_2_t0.pdf")
plotContour(prob_double_imag, N_times//2, "contur_imag_2_12T.pdf") #time 0.001
plotContour(prob_double_imag, -1, "contur_imag_2_T.pdf") #time 0.002





"""
Plottig the probability of detection along a line for double, single and triple slit:
"""
    
    

def plot1D(y, col, fig_name):
    plt.figure(figsize=(8,6))
    plt.xlabel("y")
    plt.ylabel("Probability mass function")
    plt.plot(x,y, color=col)
    plt.savefig(fig_name)
    
    

room_x = int(0.8*dim) #index of
time_index = -1
x = np.linspace(0,1,dim)

y_double = np.array(prob_double[time_index])[:,room_x]
y_double /= sum(y_double) #normalize
plot1D(y_double, "b", "1D_2.pdf")


y_single = np.array(readFilePLoc("results_probability_contur_single.txt")[time_index])[:,room_x]
y_single /= sum(y_single) #normalize
plot1D(y_single, "r", "1D_1.pdf")

y_triple = np.array(readFilePLoc("results_probability_contur_triple.txt")[time_index])[:,room_x]
y_triple /= sum(y_triple)
plot1D(y_triple, "g", "1D_3.pdf")


