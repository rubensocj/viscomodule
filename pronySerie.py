import prony as pr
import matplotlib.pyplot as plt

class PronySerie:

    def setTestOutput(self, output):
        self.time = pr.toFloatArray(output[0])
        self.tension = pr.toFloatArray(output[1])

    # test time
    def setTime(self, time):
        self.time = time

    # test tension
    def setTension(self, tension):
        self.tension = tension

    # prony serie terms
    def setTerms(self, num):
        self.num = num

    # test deformation rate
    def setRate(self, kk):
        self.kk = kk

    # prony serie tyep: with a given E_inf or without
    def setGivenEinfSerieType(self, bo):
        self.einfType = bo

    # prony given E_inf value
    def setGivenEinf(self, inEinf):
        self.givenEinf = inEinf

    # relaxation method step
    def setStep(self, step):
        self.step = step

    # set prony serie characterization to test simulation
    def setSimulationInput(self, simulation):
        self.relaxatioTimes = pr.toFloatArray(simulation[0])
        self.modules = pr.toFloatArray(simulation[1])

    # get optimum relaxation times
    def runRelaxation(self):
        if self.einfType == True:
            self.results = pr.getRelaxationTimesEinf(self.time,
                                                     self.tension,
                                                     self.givenEinf,
                                                     self.kk,
                                                     self.num,
                                                     self.step)

        else:
            self.results = pr.getRelaxationTimes(self.time,
                                                 self.tension,                                                 
                                                 self.kk,
                                                 self.num,
                                                 self.step)
    # get prony serie characterization
    def runPronySerie(self):
        self.prony = pr.prony(self.time,
                              abs(self.results.equilibrium_module),
                              abs(self.results.modules),
                              abs(self.results.relaxation_times),
                              self.num)

    # get tension from test simulation
    def runSimulation(self):
        tension = pr.tension(self.time,
                             self.kk,
                             self.givenEinf,
                             self.modules,
                             self.relaxatioTimes,
                             self.num)
        self.setTension(tension)