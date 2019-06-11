import sys
 
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QSizePolicy, QAction
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout, QHBoxLayout, QGridLayout
from PyQt5.QtWidgets import QMessageBox, QPushButton, QLabel, QComboBox, QCheckBox, QSlider
from PyQt5.QtWidgets import QWidget, QDesktopWidget, QFileDialog, QLineEdit, QDialog, QProgressBar
from PyQt5.QtWidgets import QTableWidget,QTableWidgetItem, QTableView, QHeaderView

from PyQt5 import QtGui
from PyQt5.QtGui import QIntValidator, QDoubleValidator, QDesktopServices, QIcon
from PyQt5.QtCore import pyqtSlot, QSize, Qt, QAbstractTableModel, QUrl
 
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import prony as pr
import numpy as np
import pronySerie as ps

import threading
 
class App(QMainWindow):
 
    def __init__(self):
        super().__init__()

        # ceterlizes the window
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())

        # defines the title and dimensions
        self.title = 'Viscomodule 1.0.2'
        self.width = 800
        self.height = 500

        # initializes the UI
        self.initUI()

    '''
     Initializes the UI setting the menubar and the widgets in the main window
    '''
    def initUI(self):

        #initializes the horizontal layout
        self.verLayout = QVBoxLayout()

        # initializes the main layout as a VBoxLayout
        self.horLayout = QHBoxLayout()

        # initializes the central widget as a blank widget
        self.centralWidget = QWidget()

        # the input from csv file - the 2d array time x tension from test
        self.inputProny = [[],[]]
        self.inputTension = [[],[]]
        self.outputTension = []

        # the pronySerie object
        self.pronySerie = ps.PronySerie()
        self.simulationProny = ps.PronySerie()

        # the table
        self.initializeTable()

        # initializes and configures the left and right panels
        self.widLeft = self.configureLeftPanel()
        self.widRight = self.configureRightPanel()
        self.configureStyleSheets()
        
        # set the title
        # set the size
        self.setWindowTitle(self.title)
        self.resize(self.width, self.height)

        # creates the action from menubar
        menuHelp = QAction('Ajuda', self)
        menuHelp.triggered.connect(lambda: self.openUrl('https://rubensocj.github.io/viscoprony/ajuda.html'))
        menuAbout = QAction('Sobre', self)
        menuAbout.triggered.connect(lambda: self.openUrl('https://rubensocj.github.io/viscoprony/sobre.html'))

        # creates and configures the menubar and items
        menuMain = self.menuBar()
        menuMain.addAction(menuHelp)
        menuMain.addAction(menuAbout)

        # add the widgets
        self.horLayout.addWidget(self.widLeft, 3)
        self.horLayout.addWidget(self.widRight, 1)
        
        # defines central widget layout
        self.centralWidget.setLayout(self.horLayout)

        # sets the central widget
        self.setCentralWidget(self.centralWidget)

        # sets icon
        self.setWindowIcon(QIcon('icon5.png'))  

        # shows the application 
        self.show()

    '''
     Configures the widgets in the left side on the main window
    '''
    def configureLeftPanel(self):

        # labels
        self.lblProny = QLabel('Ensaio de creep estático')
        self.lblInput = QLabel('Caracterização pelo módulo de relaxação:')
        self.lblTerms = QLabel('Número de termos:')
        self.lblSpace = QLabel('Espaçamento entre tempos de relaxação:')
        self.lblSpace.setToolTip('Configura a mantissa dos valores de tempo de relaxação testados. \n' +
                                 'Quanto menor este valor, maiores serão o esforço computacional\n' +
                                 'e a precisão da busca pelo melhor resultado')
        self.lblFile = QLabel('Arquivo:')
        self.lblFileName = QLabel('...')
        self.lblFileName.setWordWrap(True)
        self.lblRate = QLabel('Taxa de deformação:')
        
        self.lblResults = QLabel('Resultados:')
        self.lblEinf = QLabel('Módulo de Equilíbrio (E_\u221E):')
        self.lblEinfValue = QLabel('')
        self.lblK = QLabel('Taxa de deformação com o tempo:')
        self.lblKValue = QLabel('')
        
        self.lblVoid1 = QLabel('')
        self.lblVoid2 = QLabel('')
        self.lblVoid3 = QLabel('')

        self.lblExport = QLabel('Exportar resultados da caracterização:')
        
        # combobox
        self.cbxTerms = QComboBox(self)
        self.cbxTerms.addItem('7')
        self.cbxTerms.addItem('8')
        self.cbxTerms.addItem('9')
        self.cbxTerms.addItem('10')
        self.cbxTerms.addItem('11')

        self.cbxSpace = QComboBox(self)
        self.cbxSpace.addItem('10')
        self.cbxSpace.addItem('1')
        self.cbxSpace.addItem('0.5')
        self.cbxSpace.addItem('0.25')
        self.cbxSpace.addItem('0.1')

        # textfields
        self.doubleValidator = QDoubleValidator()
        self.tfdRate = QLineEdit(self)
        self.tfdRate.setValidator(self.doubleValidator)
        self.tfdEinf = QLineEdit(self)
        self.tfdEinf.setValidator(self.doubleValidator)
        self.tfdEinf.setEnabled(False)

        # checkbox
        self.ckbEinf = QCheckBox('Fornecer Módulo de Equilíbrio (E_\u221E):')
        self.ckbEinf.setChecked(False)
        self.ckbEinf.stateChanged.connect(self.changeEinf)
        self.setEinf = False # boolean tells wich system method call

        # button
        self.btnCurve = QPushButton('Ver curva', self)
        self.btnCurve.clicked.connect(self.plotCreepCurve) # creep test curve 
        self.btnProny = QPushButton('Calcular', self)
        self.btnProny.clicked.connect(self.actionProny) # prony    
        self.btnFile = QPushButton('Importar')
        self.btnFile.clicked.connect(lambda: self.openFileNameDialog('prony')) # import     
        self.btnClear = QPushButton('Limpar tudo')
        self.btnClear.clicked.connect(self.actionClearProny) # clear results    
        self.btnImage = QPushButton('Gráfico')
        self.btnImage.clicked.connect(self.plotPronyCurve) # image      
        self.btnTable = QPushButton('Arquivo csv')
        self.btnTable.clicked.connect(lambda: self.actionTableExportData('csv')) # csv
        self.btnTexTable = QPushButton('Tabela tex')
        self.btnTexTable.clicked.connect(lambda: self.actionTableExportData('tex')) # tex   
        self.setEnabledLeftButtons(False)
        
        # left panel
        self.widLeft1 = QWidget()

        # final layout
        self.leftLayout = QHBoxLayout()

        # auxiliary panels
        self.leftLayout1 = QHBoxLayout()
        self.leftLayout2 = QHBoxLayout()
        self.leftLayout3 = QHBoxLayout()
        self.leftLayout4 = QHBoxLayout()        
        self.leftLayout5 = QHBoxLayout()        
        self.leftLayout6 = QHBoxLayout()        
        self.leftLayout7 = QHBoxLayout()        
        self.leftLayout8 = QHBoxLayout()        
        self.leftLayout9 = QHBoxLayout()        
        self.leftLayout10 = QHBoxLayout()

        # build the panels
        self.leftLayout1.addWidget(self.lblFile, 1)
        self.leftLayout1.addWidget(self.lblFileName, 3)        
        self.leftLayout1.addWidget(self.btnFile, 1)
        self.leftLayout2.addWidget(self.lblVoid1, 3)
        self.leftLayout2.addWidget(self.btnCurve, 1)
        self.leftLayout3.addWidget(self.lblTerms, 5)
        self.leftLayout3.addWidget(self.cbxTerms, 1)
        self.leftLayout4.addWidget(self.lblRate, 5)
        self.leftLayout4.addWidget(self.tfdRate, 1)
        self.leftLayout5.addWidget(self.ckbEinf, 4)
        self.leftLayout5.addWidget(self.tfdEinf, 1)
        self.leftLayout6.addWidget(self.lblSpace, 5)
        self.leftLayout6.addWidget(self.cbxSpace, 1)
        self.leftLayout7.addWidget(self.lblVoid2, 3)
        self.leftLayout7.addWidget(self.btnProny, 1)
        
        self.leftLayout8.addWidget(self.btnImage, 1)
        self.leftLayout8.addWidget(self.btnTable, 1)
        self.leftLayout8.addWidget(self.btnTexTable, 1)
        
        self.leftLayout9.addWidget(self.lblEinf, 3)
        self.leftLayout9.addWidget(self.lblEinfValue, 1)
        
        self.leftLayout10.addWidget(self.lblVoid3, 3)
        self.leftLayout10.addWidget(self.btnClear, 1)

        # Fit the panels
        self.midLeft4 = QVBoxLayout()
        self.midLeft4.addWidget(self.lblProny)
        self.midLeft4.addLayout(self.leftLayout1)
        self.midLeft4.addLayout(self.leftLayout2)
        self.midWid4 = QWidget()
        
        self.midLeft1 = QVBoxLayout()
        self.midLeft1.addWidget(self.lblInput)
        self.midLeft1.addLayout(self.leftLayout3)
        self.midLeft1.addLayout(self.leftLayout4)
        self.midLeft1.addLayout(self.leftLayout5)
        self.midLeft1.addLayout(self.leftLayout6)
        self.midLeft1.addLayout(self.leftLayout7)
        self.midWid1 = QWidget() # widget input
        self.midWid1.setLayout(self.midLeft1)
        self.midLeft4.addWidget(self.midWid1) # left layout: add input
        
        self.midLeft2 = QVBoxLayout()
        self.midLeft2.addWidget(self.lblExport)
        self.midLeft2.addLayout(self.leftLayout8)
        self.midWid2 = QWidget() # widget export
        self.midWid2.setLayout(self.midLeft2)
        self.midLeft4.addWidget(self.midWid2) # left layout: add export
        self.midLeft4.addStretch(0)

        self.midWid4.setLayout(self.midLeft4) # widget left: input + export
        
        self.midLeft3 = QVBoxLayout()
        self.midLeft3.addWidget(self.lblResults)
        self.midLeft3.addWidget(self.table)
        self.midLeft3.addLayout(self.leftLayout9)
        self.midLeft3.addLayout(self.leftLayout10)
        self.midWid3 = QWidget() # widget results
        self.midWid3.setLayout(self.midLeft3)

        self.leftLayout.addWidget(self.midWid4, 1)
        self.leftLayout.addWidget(self.midWid3, 4)

        # set the layout to left layout
        self.widLeft1.setLayout(self.leftLayout)

        return self.widLeft1

    '''
     Configures the widgets in the left side on the main window
    '''
    def configureRightPanel(self):
        # left panel
        self.widRight1 = QWidget()

        # final layout
        self.rightLayout = QVBoxLayout()

        # auxiliary panels
        self.rightLayout1 = QHBoxLayout()
        self.rightLayout2 = QHBoxLayout()
        self.rightLayout3 = QHBoxLayout()
        self.rightLayout4 = QHBoxLayout()        
        self.rightLayout5 = QHBoxLayout()        
        self.rightLayout6 = QHBoxLayout()        
        self.rightLayout7 = QHBoxLayout()        
        self.rightLayout8 = QHBoxLayout()

        # labels
        self.lblTension = QLabel('Simular ensaio de creep estático')
        self.lblTime = QLabel('Duração do ensaio (s):')
        self.lblTimeRate = QLabel('Número de leituras:')
        self.lblTensionEinf = QLabel('Módulo de Equilíbrio (E_\u221E):')
        self.lblTensionK = QLabel('Taxa de deformação:')
        self.lblTensionTerms = QLabel('Número de termos:')  
        self.lblTensionFile = QLabel('Arquivo:')
        self.lblTensionFileName = QLabel('...')
        self.lblTensionFileName.setWordWrap(True)      
        self.lblTensionVoid1 = QLabel('')      
        self.lblTensionVoid2 = QLabel('')      
        self.lblTensionVoid3 = QLabel('')

        # textfield
        self.validator = QIntValidator()
        self.tfdTensionTime = QLineEdit(self)
        self.tfdTensionTime.setValidator(self.validator)
        self.tfdTensionRate = QLineEdit(self)
        self.tfdTensionRate.setValidator(self.validator)
        self.tfdTensionEinf = QLineEdit(self)
        self.tfdTensionEinf.setValidator(self.doubleValidator)
        self.tfdTensionK = QLineEdit(self)
        self.tfdTensionK.setValidator(self.doubleValidator)

        # button
        self.btnTensionFile = QPushButton('Importar')
        self.btnTensionFile.clicked.connect(lambda: self.openFileNameDialog('tension')) # import
        self.btnTension = QPushButton('Simular', self)
        self.btnTension.clicked.connect(self.actionTension) # simulate      
        self.btnTensionImage = QPushButton('Ver curva')
        self.btnTensionImage.clicked.connect(self.plotSimulationCurve) # plot    
        self.btnTensionTable = QPushButton('Arquivo csv')
        self.btnTensionTable.clicked.connect(self.actionSimulationExportData) # export    
        self.btnTensionClear = QPushButton('Limpar tudo')
        self.btnTensionClear.clicked.connect(self.actionClearTension) # clear
        self.setEnabledRightButtons(False)

        # build the panels
        self.rightLayout1.addWidget(self.lblTime, 3)
        self.rightLayout1.addWidget(self.tfdTensionTime, 1)
        self.rightLayout2.addWidget(self.lblTimeRate, 3)
        self.rightLayout2.addWidget(self.tfdTensionRate, 1)
        self.rightLayout3.addWidget(self.lblTensionEinf, 3)
        self.rightLayout3.addWidget(self.tfdTensionEinf, 1)
        self.rightLayout4.addWidget(self.lblTensionK, 3)
        self.rightLayout4.addWidget(self.tfdTensionK, 1)
        self.rightLayout5.addWidget(self.lblTensionFile, 1)
        self.rightLayout5.addWidget(self.lblTensionFileName, 5)
        self.rightLayout5.addWidget(self.btnTensionFile, 2)
        self.rightLayout6.addWidget(self.lblTensionVoid1, 3)
        self.rightLayout6.addWidget(self.btnTension, 1)
        self.rightLayout7.addWidget(self.lblTensionVoid2, 1)
        self.rightLayout7.addWidget(self.btnTensionTable, 2)
        self.rightLayout7.addWidget(self.btnTensionImage, 2)
        self.rightLayout8.addWidget(self.lblTensionVoid3, 3)
        self.rightLayout8.addWidget(self.btnTensionClear, 1)

        self.rightLayout.addWidget(self.lblTension)
        self.rightLayout.addLayout(self.rightLayout1)
        self.rightLayout.addLayout(self.rightLayout2)
        self.rightLayout.addLayout(self.rightLayout3)
        self.rightLayout.addLayout(self.rightLayout4)
        self.rightLayout.addLayout(self.rightLayout5)
        self.rightLayout.addLayout(self.rightLayout6)
        self.rightLayout.addLayout(self.rightLayout7)
        self.rightLayout.addLayout(self.rightLayout8)

        self.rightLayout.addStretch(0)

        # set the layout to left layout
        self.widRight1.setLayout(self.rightLayout)

        return self.widRight1     

    '''
     Checks if the input is correct and then do the prony algoritm to get the
     Prony series constants and the relaxation times
    '''
    def actionProny(self):
        
        if self.lblFileName.text() == '...':
            QMessageBox.about(self, 'Aviso', 'Nenhum arquivo selecionado')
        elif len(self.tfdRate.text()) == 0 or self.tfdRate.text() == '0':
            QMessageBox.about(self, 'Aviso', 'Valor da taxa de deformação não informado')
        elif self.setEinf == True and len(self.tfdEinf.text()) == 0:
            QMessageBox.about(self, 'Aviso', 'Valor do módulo de equilíbrio não informado')
        else:

            self.pronySerie.setTerms(int(self.cbxTerms.currentText()))
            self.pronySerie.setRate(float(self.tfdRate.text()))
            self.pronySerie.setStep(float(self.cbxSpace.currentText()))

            if self.setEinf == True:
                self.pronySerie.setGivenEinfSerieType(True)
                self.pronySerie.setGivenEinf(float(self.tfdEinf.text()))
            else:
                self.pronySerie.setGivenEinfSerieType(False)

            try:
                QMessageBox.about(self, 'Aviso', 'Esta operação poderá demorar alguns minutos.\n' +
                                  'Pressione ''OK'' para prosseguir')
                self.pronySerie.runRelaxation()
            except UnboundLocalError:
                QMessageBox.about(self, 'Erro', 'Não foi possível encontrar um resultado com parâmetros informados')
            except:
                QMessageBox.about(self, 'Erro', 'Um erro ocorreu: pronySerie.runRelaxation')            
            else:
                try:
                    self.pronySerie.runPronySerie()
                    self.updateTable()

                    self.lblEinfValue.setText(str("{:.2E}".format(self.pronySerie.results.equilibrium_module)))
                    self.setEnabledLeftButtons(True)
                except:
                    QMessageBox.about(self, 'Erro', 'Um erro ocorreu: pronySerie.runPronySerie')

    '''
     Checks if the input is correct and then do the tension algoritm to get the
     tension values from simulation of CCRMT
    '''
    def actionTension(self):
        
        if self.lblTensionFileName.text() == '...':
            QMessageBox.about(self, "Aviso", "Nenhum arquivo selecionado")
        elif len(self.tfdTensionK.text()) == 0 or self.tfdTensionK.text() == '0':
            QMessageBox.about(self, 'Aviso', 'Valor da taxa de deformação não informado') 
        elif self.tfdTensionTime.text() == '' or self.tfdTensionRate.text() == '' or self.tfdTensionEinf.text() == '':
            QMessageBox.about(self, "Aviso", "Preencha todos os campos")        
        else:

            self.simulationProny.setTime(np.linspace(0, int(self.tfdTensionTime.text()), int(self.tfdTensionRate.text())))
            self.simulationProny.setRate(float(self.tfdTensionK.text()))
            self.simulationProny.setGivenEinf(float(self.tfdTensionEinf.text()))
            self.simulationProny.setTerms(len(self.simulationProny.modules))

            self.simulationProny.runSimulation()

            # enable results option buttons
            self.setEnabledRightButtons(True)

    '''
     Exports the result data in the table to a csv  or text table
    '''
    def actionTableExportData(self, fxt):

        # organizes the time x tension output matrix
        o1 = np.asarray([self.pronySerie.results.modules, self.pronySerie.results.relaxation_times])
        print('o1:', o1)
        output = o1.transpose()
        print('output:', output)

        if fxt == 'csv':        
            self.saveCSVFileDialog(output)
        elif fxt == 'tex':
            print('tex')
            self.saveTEXTableFileDialog(o1, self.pronySerie.results.equilibrium_module, self.pronySerie.kk)

    def actionSimulationExportData(self):
        # organizes the time x tension output matrix
        o1 = np.asarray([self.simulationProny.time, self.simulationProny.tension])
        output = o1.transpose()
        
        self.saveCSVFileDialog(output)        

    '''
     Shows file selector and deals with the selected file
    '''
    def openFileNameDialog(self, typeInput):    
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filePath, _ = QFileDialog.getOpenFileName(self,"Importar dados de arquivo CSV", "","csv(*.csv)", options=options)
        if filePath:
            print(typeInput, filePath)

            fileName = filePath.split("/")

            if typeInput == 'prony':
                # imports the csv and puts into a prony list
                self.pronySerie.setTestOutput(pr.readCSV(filePath))
                
                # updates the path to file label
                self.lblFileName.setText(fileName[len(fileName)-1])
                self.lblFileName.setToolTip(filePath)

            elif typeInput == 'tension':
                # imports the csv and puts into a prony list
                self.simulationProny.setSimulationInput(pr.readCSV(filePath))
                
                # updates the path to file label
                self.lblTensionFileName.setText(fileName[len(fileName)-1])
                self.lblTensionFileName.setToolTip(filePath)

    '''
     Shows dialog "Save as csv file" to save the result from simulation or from table
    '''
    def saveCSVFileDialog(self, output):    
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,'Salvar resultado como arquivo CSV','','csv(*.csv)', options=options)
        if fileName:
            print(fileName)
            np.savetxt(fileName + '.csv', output, delimiter = ',')

    '''
     Shows dialog "Save as tex file" to save the result from simulation or from table in a tex table
    '''
    def saveTEXTableFileDialog(self, output, eInf, k_z):    
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,'Salvar resultado como tabela TEX','','txt(*.txt)', options=options)
        if fileName:
            print(fileName)
            pr.writeTexTable(fileName + '.txt', output, eInf, k_z)

    '''
     Salva a imagem
    '''
    def saveImageFileDialog(self, typeInput):    
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,'Salvar imagem','','JPEG Image (*.jpg);;PNG Image (*.png)', options=options)
        if fileName:
            print(typeInput, fileName)

            #fileName = filePath.split("/")

            self.plot.savePlot(fileName, typeInput)

    '''
     Plots the creep test outpu in to a curve
    '''
    def plotCreepCurve(self):        
        if self.lblFileName.text() == '...':
            QMessageBox.about(self, 'Aviso', 'Nenhum arquivo selecionado')
        else:
            self.plotCreep = PlotWindow(self.pronySerie.time, self.pronySerie.tension, 'Tempo - t (s)', 'Tensão - $\sigma$(t) (MPa)')
            self.plotCreep.show()

    def plotPronyCurve(self):
        self.plot = PlotWindow(self.pronySerie.time, self.pronySerie.prony, 'Tempo - t (s)', 'Módulo de relaxação - E(t) (MPa)')
        self.plot.show()

    def plotSimulationCurve(self):
        self.plotSimulation = PlotWindow(self.simulationProny.time, self.simulationProny.tension, 'Tempo - t (s)', 'Tensão - $\sigma$(t) (MPa)')
        self.plotSimulation.show()

    '''
     Clear all the left panel data and the plot
    '''
    def actionClearProny(self):
        if self.table.rowCount() != 0:
            self.lblFileName.setText('...')
            self.lblEinfValue.setText('')
            self.tfdRate.setText('')
            self.tfdEinf.setText('')
            self.table.setRowCount(0)
            self.pronySerie = ps.PronySerie()
            self.setEnabledLeftButtons(False)

    def actionClearTension(self):
        self.lblTensionFileName.setText('...')
        self.tfdTensionTime.setText('')
        self.tfdTensionRate.setText('')
        self.tfdTensionK.setText('')
        self.tfdTensionEinf.setText('')
        self.simulationProny = ps.PronySerie()
        self.setEnabledRightButtons(False)

    '''
     Enables the buttons in the left and right panel
    '''
    def setEnabledLeftButtons(self, bo):        
        self.btnImage.setEnabled(bo)
        self.btnTable.setEnabled(bo)
        self.btnTexTable.setEnabled(bo)
        self.btnClear.setEnabled(bo)
        
    def setEnabledRightButtons(self, bo):        
        self.btnTensionImage.setEnabled(bo)
        self.btnTensionTable.setEnabled(bo)
        self.btnTensionClear.setEnabled(bo)

    def changeEinf(self, state):
        if state == Qt.Checked:
            self.tfdEinf.setEnabled(True)
            self.setEinf = True
            print('checked')
        else:
            self.tfdEinf.setEnabled(False)
            self.setEinf = False
            print('not checked')

    '''
     Initializes the empty table
    '''
    def initializeTable(self):
        self.table = QTableWidget()
        
        # set column count
        self.table.setColumnCount(2)
        self.table.setRowCount(0)

        # header label text
        self.table.setHorizontalHeaderLabels(["Tempos de \nrelaxação ("u"\u03C1)", "Constantes da \nsérie de Prony (E)"])

        # header width to fit in the widget
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.Stretch)

        # header label text alignment
        self.table.horizontalHeaderItem(0).setTextAlignment(Qt.AlignHCenter)
        self.table.horizontalHeaderItem(1).setTextAlignment(Qt.AlignHCenter)
        
        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()

    '''
     Set the result data in to the table
    '''
    def updateTable(self):
        self.table.setRowCount(self.pronySerie.num)
        i=0
        while i < self.pronySerie.num:
            
            itemTime = QTableWidgetItem(str("{:.2E}".format(abs(self.pronySerie.results.relaxation_times[i]))))
            itemTime.setTextAlignment(Qt.AlignCenter)
            itemTime.setFlags(Qt.ItemIsEnabled)
            self.table.setItem(i,0,itemTime)

            itemModule = QTableWidgetItem(str("{:.2E}".format(self.pronySerie.results.modules[i])))
            itemModule.setTextAlignment(Qt.AlignCenter)
            itemModule.setFlags(Qt.ItemIsEnabled)
            self.table.setItem(i,1,itemModule)
            
            i = i+1

    '''
     Configures the CSS atributtes of QWidgets
    '''
    def configureStyleSheets(self):
        
        self.lblTension.setObjectName('lblTension')
        self.lblTension.setStyleSheet('QWidget#lblTension {font-weight: bold}')
        
        self.lblProny.setObjectName('lblProny')
        self.lblProny.setStyleSheet('QWidget#lblProny {font-weight: bold}')

        self.lblInput.setObjectName('lblInput')
        self.lblInput.setStyleSheet('QWidget#lblInput {font-weight: bold}')
        self.lblResults.setObjectName('lblResults')
        self.lblResults.setStyleSheet('QWidget#lblResults {font-weight: bold}')
        self.lblExport.setObjectName('lblExport')
        self.lblExport.setStyleSheet('QWidget#lblExport {font-weight: bold}')
        
        self.midWid1.setObjectName('midWid1')
        self.midWid1.setStyleSheet('QWidget#midWid1 {border: 1px solid LightGray}')

        self.midWid2.setObjectName('midWid2')
        self.midWid2.setStyleSheet('QWidget#midWid2 {border: 1px solid LightGray}')

        self.midWid3.setObjectName('midWid3')
        self.midWid3.setStyleSheet('QWidget#midWid3 {border: 1px solid LightGray}')
        
        self.widLeft1.setObjectName("widLeft")
        self.widLeft1.setStyleSheet("QWidget#widLeft {border: 1px solid LightGray}")
        
        self.widRight1.setObjectName('widRight1')
        self.widRight1.setStyleSheet('QWidget#widRight1 {border: 1px solid LightGray}')

    '''
     Open a url
    '''
    def openUrl(self, url):
        QDesktopServices.openUrl(QUrl(url))

'''
 This class show a Dialog with a plot
'''
class PlotWindow(QDialog):
    def __init__(self, xx, yy, xlabel, ylabel, parent=None):
        super(PlotWindow, self).__init__(parent)
        
        plt.style.use('seaborn-paper')
        plt.rcParams['font.sans-serif'] = "DejaVu Sans"
        plt.rcParams['font.family'] = "sans-serif"

        self.setWindowTitle('Viscomodule')

        # a figure instance to plot on
        self.figure = Figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        ax = self.figure.add_subplot(111) # create an axis
        ax.clear() # discards the old graph
        color = 'orangered'
        ax.plot(xx, yy, color=color)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.tick_params(axis='both', which='major', labelsize=10)
##        ax.text(0.85, 0.9, '$E_\infty$ = ' + str(einf), ha='center', va='center', transform=ax.transAxes, fontsize=10,
##                bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
        
        self.canvas.draw() # refresh canvas

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,4))

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
                
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
