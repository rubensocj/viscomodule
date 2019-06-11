# This script aims to simulate the interpolation function
# adopted to tensions for characterization of viscoelastic materials
# by CCMRT test

# Contains:
#
# tension
# prony
# plot
# getRelaxationTimes
# getRelaxationTimesEinf
# isAllPositiveArray
# readCSV
# writeCSV
# writeTexTable
# toFloatArray
# main

# import modules
import numpy as np
import scipy.linalg as sc
import math as mt
import matplotlib.pyplot as plt
import system as st
import csv
import collections as col

# define here auxiliary functions

"""
 Francisco adapted interpolation function
 
 time - time array
 kz - deformation constant rate
 eInf - single spring constant
 eArray - spring constant array (relaxation modulus)
 pArray - relaxation time array
 num - number of terms in Prony series
 
  Returns - array with values of tension
"""
def tension(time, kz, eInf, eArray, pArray, num):
        i = 0
        tension = []
        
        while i < len(time):
                j = 0
                
                x1 = 0
                while j < num:
                        # auxiliary variable
                        x1 = x1 + eArray[j]*pArray[j]*(1 - mt.exp(-time[i]/pArray[j]))
                        
                        # appends the result to sumArray
                        #sumArray.append(x1)
                        j = j+1
                
                # auxiliary variable, means the tension
                x2 = kz*(eInf*time[i] + x1)
                
                # appends the result to tension array
                tension.append(x2)
                i = i+1
                        
        return tension

"""
 The Prony series function
 
 time - time array
 eInf - single spring constant
 eArray - spring constant array (relaxation modulus)
 pArray - relaxation time array
 num - number of terms in Prony series
 
  Returns - array with values of tension
"""
def prony(time, eInf, eArray, pArray, num):
        # initializes the prony constants array
        pronyArray = np.empty(shape=(len(time)), dtype=float)
                              
        for i in np.arange(0, len(time), 1):

                # the summatory term of prony series equation
                ss = 0
                for j in np.arange(0, num, 1):
                        ss = ss + eArray[j]*mt.exp(-time[i]/pArray[j])

                # sums the summatory on the i-th time on the E_inf
                ee = eInf + ss

                # appends the i-th prony constant in function of time
                pronyArray[i] = ee

        return pronyArray

        
"""
 Plot time x tension curve
 
 time - time array
 tension - tension array
"""
def plot(time, tension):
        
        # Get the coefficients from a least squares polynomial fit with degree 2
        coef = np.polyfit(time, tension, 2)
        print("coeff: ", coef)

        # Build a a*x^2 + b*x + c like function where a, b, c are the coefficients from coef array
        poly = np.poly1d(coef)
        print("adjusted function: \n", np.poly1d(poly))
        
        # plots time x tension curve in blue
        # plots adjusted quadratic curve in red dashed
        lines = plt.plot(time, tension, time, poly(time))
        plt.setp(lines[0], color='b', linewidth=2.0, linestyle='-')
        plt.setp(lines[1], color='r', linewidth=2.0, linestyle=':')
        
        # custom plot info
        plt.xlabel('time, t')
        plt.ylabel(r'tension, $\sigma$')
        plt.title('tension x time output')

        # shows the plot
        plt.show()

"""
 Calculates the best relaxation times array and the best input deformation rate
 for the given Prony series number of terms and time and tension values array.
 Iterates through the relaxation times array in the defined exponent range
 and the deformation rate values from 1.0 to 4.0.

 time - time values array
 tension - tension values array
 kz - deformation ratio
 num - number of terms in Prony series
 
 Returns - a collection
"""
def getRelaxationTimes(time, tension, kz, num, step):

        # initializes relaxation times array
        # value range: 10e-6 to 10e10
        # set the geomspace numpy array
        # geomspace initial exponents
        ini = -6    # lowest possible exponent
        fin = 10    # highest possible exponent

        # the 0-th exponent when the last exponent = 10
        lastini = 11 - num

        # set kz values
        # value range: 1.0 to 4.0 by 1.0
##        kz = np.linspace(1.0, 3.0, 3)

        # set pp mantissa values
        # value range: 1.0 to 9.0 by 1
        # e.g. 1E3, 2E3, 3E4...
##        ppmant = np.linspace(1.0, 9.0, 9.0)
        ppmant = np.arange(1.0, 10.0, step=step)

        # initializes ee array to keep Prony series constants values
        ee = np.empty(shape=(num), dtype=float)

        eeopt = np.zeros(num)   # array to keep the optimum array    
        ecopt = 1               # variation coefficient as the minimum possible: 0
        rsopt = 0               # r-squared optimum value

        # initializes matrix to keep the ee values
        mtxe = np.empty(shape=((lastini-ini+1),len(ppmant),len(ee)+1))

        # loop through pp
        for i in np.arange(0, (lastini-ini+1), 1):
                
                # i-th start and stop exponents from i-th geomspace array
                start = ini + i
                stop = (num - 1) + start

                for j in np.arange(0, len(ppmant), 1):

                        startValue = ppmant[j]*np.float(10)**(start)
                        stopValue = ppmant[j]*np.float(10)**(stop)

                        # i-th geomspace array
                        # logaritmic decades
                        pp = np.geomspace(startValue, stopValue, num)
                        #print(pp)

                        # set vector B
                        vec = st.vectorB(tension, time, pp, num)

                        # loop through kz
##                        for j in np.arange(0, len(kz), 1):
                        
                        # set matrix
                        mtx = st.matrixA(kz, time, pp, num)

                        # solve
                        ee = st.solve(mtx, vec)

                        # store ee j-th array in 3-dimensional matrix
                        # each input in this matrix is a Prony series constant array
                        # next step, obtain the best values by the stop criteria
                        mtxe[i][j] = ee

                        if not np.array_equal(ee, np.ones((len(ee),), dtype=int)):

                                # check the variance and variance coefficient of the exponents
                                # get the exponent array
                                expvec = np.floor(np.log10(np.abs(ee[1:]))).astype(int)
                                expstd = np.std(expvec)
                                expmen = np.mean(expvec)
                                expcov = expstd/expmen
                                expinf = np.floor(np.log10(np.abs(ee[0]))).astype(int)
##                                print('expinf: {} | expmen: {}'.format(expinf, expmen))

                                # fits ee values to 3th degree polynomy
                                print('pp:', pp)
                                print('ee:',ee)
                                reg = polynomialRegression(pp, abs(ee[1:]), 3)
        ##                        print(isWellAdjusted(ee, reg))
##                                print('ecopt:', ecopt)
        ##                        if(isWellAdjusted(ee, reg) and reg.r_squared > rsopt):
        ##                        if(reg.r_squared > rsopt and reg.r_squared < 0.8):
                                if(reg.r_squared > rsopt and reg.r_squared < 0.8 and
                                   expcov < ecopt and expcov > 0.10 and
                                   0.7*expmen <= expinf <= 1.3*expmen):  
                                        rsopt = reg.r_squared
                                        eeopt = ee
        ##                                        kzopt = kz[j]
                                        ppopt = pp
                                        maopt = mtx
                                        vbopt = vec
                                        ecopt = expcov                        

        Relaxation = col.namedtuple('Relaxation',
                                    ['relaxation_times',
                                     'modules',
                                     'equilibrium_module',
                                     'deformation_rate',
                                     'matrixA',
                                     'vectorB',
                                     'all_modules'])
        rr = Relaxation(ppopt, eeopt[1:len(eeopt)], eeopt[0], kz, maopt, vbopt, mtxe)

        return rr

def getRelaxationTimesEinf(time, tension, eInf, kz, num, step):

        # initializes relaxation times array
        # value range: 10e-6 to 10e10
        # set the geomspace numpy array
        # geomspace initial exponents
        ini = -6    # lowest possible exponent
        fin = 10    # highest possible exponent

        # the 0-th exponent when the last exponent = 10
        lastini = 11 - num

        # set pp mantissa values
        # value range: 1.0 to 9.0 by 1
        # e.g. 1E3, 2E3, 3E4...
##        ppmant = np.linspace(1.0, 9.0, 9.0)
        ppmant = np.arange(1.0, 10.0, step=step)

        # initializes ee array to keep Prony series constants values
        ee = np.empty(shape=(num), dtype=float)

        eeopt = np.zeros(num)   # array to keep the optimum array    
        ecopt = 1               # variation coefficient as the minimum possible: 0
        rsopt = 0               # r-squared optimum value

        # initializes matrix to keep the ee values
        mtxe = np.empty(shape=((lastini-ini+1),len(ppmant),len(ee)))

        # loop through pp
        for i in np.arange(0, (lastini-ini+1), 1):
                
                # i-th start and stop exponents from i-th geomspace array
                start = ini + i
                stop = (num - 1) + start

                for j in np.arange(0, len(ppmant), 1):

                        startValue = ppmant[j]*np.float(10)**(start)
                        stopValue = ppmant[j]*np.float(10)**(stop)

                        # i-th geomspace array
                        # logaritmic decades
                        pp = np.geomspace(startValue, stopValue, num)

                        # set vector B
                        vec = st.vectorBred(tension,
                                            time,
                                            pp,
                                            eInf,
                                            kz,
                                            num)
            
                        # set matrix
                        mtx = st.matrixAred(kz,
                                            time,
                                            pp,
                                            num)

                        # solve
                        ee = st.solve(mtx, vec)

                        # store ee j-th array in 3-dimensional matrix
                        # each input in this matrix is a Prony series constant array
                        # next step, obtain the best values by the stop criteria
                        mtxe[i][j] = ee

                        if not np.array_equal(ee, np.ones((len(ee),), dtype=int)):

                                # check the variance and variance coefficient of the exponents
                                # get the exponent array
                                expvec = np.floor(np.log10(np.abs(ee))).astype(int)
                                expstd = np.std(expvec)
                                expmen = np.mean(expvec)
                                expcov = expstd/expmen
                                expinf = np.floor(np.log10(np.abs(eInf))).astype(int)
##                                print('expinf: {} | expmen: {}'.format(expinf, expmen))

                                # fits ee values to 3th degree polynomy
                                print('pp:', pp)
                                print('ee:',ee)
                                reg = polynomialRegression(pp, abs(ee), 3)
##                                print('r² optimum:', rsopt)
##                                if(reg.r_squared > rsopt and reg.r_squared < 0.8 and
                                if(0.7*expmen <= expinf <= 1.5*expmen):   
                                        rsopt = reg.r_squared
                                        eeopt = ee
                                        ppopt = pp
                                        maopt = mtx
                                        vbopt = vec
                                        ecopt = expcov

        Relaxation = col.namedtuple('Relaxation',
                                    ['relaxation_times',
                                     'modules',
                                     'equilibrium_module',
                                     'deformation_rate',
                                     'matrixA',
                                     'vectorB',
                                     'all_modules'])
        rr = Relaxation(ppopt, eeopt, eInf, kz, maopt, vbopt, mtxe)

        return rr

"""
 Fit the prony series modules to a model by 3th degre polynomial regression

 xx - x values (log(p) - log relaxation times prony series array)
 yy - y values (E - prony modules)
 degree - polynomial degree. Studies shows usualy is a 3th degre

 return - a collection {model, R-Squared, 1 degree derivative, 2 degree derivative]
"""
def polynomialRegression(xx, yy, degree):
        xk = np.log(xx)
        yk = yy
        wk = np.polyfit(xk, yk, degree)
        mk = np.poly1d(wk)

        rsqd = rsquared(yk, mk(xk))

        modelStudy = col.namedtuple('ModelStudy', ['fitted_model', 'r_squared'])
        ms = modelStudy(mk, rsqd)
        return ms

"""
 Computes the R-Squared to a model

 yy - observed values (prony series module array)
 ff - modeled values, a model obtained by polyfit1d

 return - R-Square value
"""
def rsquared(yy, ff):
    ssres=0
    sstot=0
    for i in range(len(yy)):
        ssres = ssres + (yy[i] - ff[i])**2
        sstot = sstot + (yy[i] - np.mean(yy))**2

    rsq = 1-(ssres/sstot)
    return rsq

"""
 Check if the model is well adjusted
 The second derivative must be negative, and then positive, as seen in tests
"""
def isWellAdjusted(array, polyreg):
        if (polyreg.derivative2(max(range(len(array)))*0.35) < 0 and polyreg.derivative2(max(range(len(array)))*0.9) > 0):
                return True
        else: return False

"""
 Tests if an input array is all positive. Returns True if all elements in
 the array is positive, or False otherwise.

 array - input array
 
 Returns - a boolean
"""
def isAllPositiveArray(array):
        test = True
        for i in array:
                if i < 0:
                        test = False
                        break

        return test

"""
 Reads a csv file with two columns. It must be p and E respectively and
 must not have the columns header in the file

 Filename - file name within the .csv extension
 
 Returns - a 2d array with p and E columns respectively
"""
def readCSV(fileName):
        with open(fileName) as csvFile:
                reader = csv.reader(csvFile, delimiter = ',')
                cstE = []
                cstP = []
                results = []
                for i in reader:
                        results.append(i)
                        # ignores column name
                        cstP.append(i[0])
                        cstE.append(i[1])
                        
        csvFile.close()
        
        return [cstP, cstE]

"""
 @deprecated
 
 Writes a csv file.
 
 fileName - name of the output file within the .csv extension
 tension - an array of tension values, where each line must be an array
 
"""
def writeCSV(fileName, tension):
        with open(fileName, 'w') as csvFile:
                writer = csv.writer(csvFile, delimiter = ',')
                writer.writerows(tension)
                
        csvFile.close()

""" 
 Writes a txt file with a table in LaTeX.
 
 fileName - name of the output file within the .csv extension
 table - an 2d-array with the values E(t) and p(t)
 e_inf - equilibrium module
 k_z - deformation rate
 
"""
def writeTexTable(fileName, table, e_inf, k_z):
        print('OOOOOOOOOOOOOOOOOOOOOOOOOO')
        with open(fileName, 'w') as textFile:
                textFile.write('\\begin{table}[htb] \n')
                textFile.write('\t \\centering \n')
                textFile.write('\t \\setlength{\\'+'tabcolsep}{2pt} \n')
                textFile.write('\t \\'+'renewcommand{\\arraystretch}{1} \n')
                textFile.write('\t \\captionsetup{width=11cm} \n')
                textFile.write('\t \\Caption{\\label{tab} Constantes da série de Prony do módulo de relaxação}{ \n')
                textFile.write('\t\t \\begin{tabular}{ccc} \n')
                textFile.write('\t\t\t \\toprule \n')
                textFile.write('\t\t\t i & $\\'+'rho_i$ (sec) & $E_i$ (MPa) \\\\ \n')
                textFile.write('\t\t\t \\midrule \\midrule \n')
                textFile.write('\t\t\t \\hline \n')
                
                for x in range(len(table[0])):
                        textFile.write('\t\t\t {} & {} & {} \\\\ \n'.format(x+1,
                                                                            str("{:.2E}".format(table[1][x])),
                                                                            str("{:.2E}".format(table[0][x]))))

                textFile.write('\t\t\t & $E_\\infty = {} $ & \\\\ \n'.format(str("{:.2E}".format(e_inf))))
                textFile.write('\t\t\t & $k = {} $ & \\\\ \n'.format(str("{:.2E}".format(k_z))))
                textFile.write('\t\t\t \\bottomrule \n')
                textFile.write('\t\t \\end{tabular} \n')
                textFile.write('\t }{ \n')
                textFile.write('\t \\Fonte{} \n')
                textFile.write('\t \\Nota{} \n')
                textFile.write('\t } \n')
                textFile.write('\\end{table} \n')


        textFile.close()
        
"""
 Converts an String array in to a float array
 
 a - a string array
 
 Returns - a float array
"""
def toFloatArray(a):
        f = []  
        i = 0
        while i < len(a):
                f.append(float(a[i]))
                i = i+1
        
        return f


"""
 The main method
"""
def main():
        # Kim, pp. 145

        # imports the csv file with p and E values
        input = readCSV('jun2019-KimGeneratedSerie_rate5e-6csv.csv')
        p = toFloatArray(input[0])
        e = toFloatArray(input[1])
                
        # time array
        # linspace(init, end, number os samples)
        t = np.linspace(0, 3600, 10000)
        # deformation rate
        k = 5e-6
        # single spring constant (relaxation modulus) in Pa
        eInf = 34.5 # kim
##        #eInf = 2.24e6 # park
        # number of terms in Prony series
        n = 11
##                
##        # tension
##        outputTension = tension(t, k, eInf, e, p, n)
##
##        # organizes the time x tension output matrix
##        o1 = np.asarray([t, outputTension])
##        output = o1.transpose()
##
##        print('o1: ', o1)
##        print('...')
##        print('output: ', output)
##
##        # creates the output.csv file
##        np.savetxt("outputkim.csv", output, delimiter = ",")

        # plot
        #plot(t, outputTension)
        print(e)
        print(p)

        serie = prony(t,eInf,np.abs(e),np.abs(p),len(e))
        print(serie)

        plt.style.use('seaborn-paper')
        plt.rcParams['font.sans-serif'] = "DejaVu Sans"
        plt.rcParams['font.family'] = "sans-serif"
        ax = plt.subplot(111) # create an axis
        ax.clear() # discards the old graph
        color = 'orangered'

        ax.plot(t, serie, color=color)

        ax.set_xlabel('Tempo - t (s)', fontsize=11)
        ax.set_ylabel('Módulo de relaxação - E(t) (MPa)', fontsize=11)
        ax.tick_params(axis='both', which='major', labelsize=10)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,4))

        plt.show()

##        k=1
##        einf=1
##        e=[10,20,30]
##        p=[0.1,1,10]
##        t=[1,2]
##
##        sigma = tension(t,k,einf,e,p,len(e))
##        print(sigma)


        
        
# the main module
if __name__ == "__main__":
  main()
        
