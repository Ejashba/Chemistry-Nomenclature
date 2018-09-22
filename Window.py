__author__ = 'Emerson Shoichet-Bartus'
__high_school__ = 'Upper Canada College'
__house__ = "Martland's"
__graduation_year__ = '2016'
__version__ = 1.1

"""
This file defines the GUI window for Nomenclature Boss

All buttons and GUI objects with which the user interacts are defined here

This took quite a long time to create
"""

import tkinter as tk
from Chem import Nomenclature # import all the nomenclature methods from the chemistry file

class GUI(tk.Tk):

    def __init__(self):
        tk.Tk.__init__(self)
        self.frame = tk.Frame(self)
        #self.geometry('{}x{}'.format(500, 100))
        self.frame.master.title('Nomenclature Boss')
        self.resizable(height=tk.FALSE, width=tk.FALSE)  # makes window NOT resizable
        self.bind("<Return>", self.displayAns)
        self.makeNamingMenu()
        self.makeFNRadios()
        self.makeEnterLabel()
        self.makeAnsLabel()
        self.makeAnsLabelText()
        self.makeNomSystemRadios()
        self.makeUserInField()
        self.makeSubmitButton()

    def makeNamingMenu(self):
        self.nomTypeVar = tk.StringVar(self)
        self.nomTypeVar.set('Select Nomenclature Type')
        self.namingMenu = tk.OptionMenu(self, self.nomTypeVar, 'Alk|anes/enes/ynes', 'Polyatomic Oxyanions', 'Acids')
        self.namingMenu.grid(row=2, column=3, padx=2, pady=2)

    # formula to name & name to formula radio buttons
    def makeFNRadios(self):
        self.FNOption = tk.IntVar()
        self.FNOption.set(-1)
        self.FtoNRadio = tk.Radiobutton(self, text='Formula to Name', value=0, variable=self.FNOption)
        self.FtoNRadio.grid(row=0, column=3)
        self.NtoFRadio = tk.Radiobutton(self, text='Name to Formula', value=1, variable=self.FNOption)
        self.NtoFRadio.grid(row=1, column=3)
        self.FtoNRadio.bind(sequence="<Button-1>", func=self.updateEnterText)  # Button-1 is left mouse button
        self.NtoFRadio.bind(sequence="<Button-1>", func=self.updateEnterText)

    def makeEnterLabel(self):
        self.enterText = tk.StringVar()
        self.enterText.set('Enter:')
        self.enterLabel = tk.Label(self, textvariable=self.enterText)
        self.enterLabel.grid(row=0, column=0)

    def makeAnsLabel(self):
        self.ansLabel = tk.Label(self, text='Answer:')
        self.ansLabel.grid(row=1, column=0)

    def makeAnsLabelText(self):
        self.ansLabelText = tk.StringVar()
        self.ans = tk.Label(self, textvariable=self.ansLabelText)
        self.ans.grid(column=1, row=1)

    # nomenclature system (common or IUPAC) radio buttons
    def makeNomSystemRadios(self):
        self.nomSystemInt = tk.IntVar()
        self.nomSystemInt.set(-1)
        self.commonRadio = tk.Radiobutton(self, text='Common', value=0, variable=self.nomSystemInt)
        self.commonRadio.grid(row=3, column=3)
        self.IUPACRadio = tk.Radiobutton(self, text='IUPAC\t', value=1, variable=self.nomSystemInt)
        self.IUPACRadio.grid(row=4, column=3)

    # field that user types in
    def makeUserInField(self):
        self.inVar = tk.StringVar()
        self.userIn = tk.Entry(self, textvariable=self.inVar)
        self.userIn.grid(row=0, column=1, pady=2, padx=2)

    def makeSubmitButton(self):
        self.submitButton = tk.Button(self, text='Submit')
        self.submitButton.grid(row=0, column=2)
        self.submitButton.bind(sequence="<Button-1>", func=self.displayAns)

    def displayAns(self, *args):
        if self.nomTypeVar.get() == 'Select Nomenclature Type':
            self.ansLabelText.set('Please select a Nomenclature Type')
        elif self.FNOption.get() == -1:
            self.ansLabelText.set('Please choose Formula to Name or Name to Formula')
        elif self.nomSystemInt.get() == -1:
            self.ansLabelText.set('Please choose either Common or IUPAC Nomenclature')
        else:
            textInput = self.inVar.get()
            self.funcMap =  {  # formula to name is 0, name to formula is 1
                            (0, 'Alk|anes/enes/ynes'): Nomenclature.alkneFtoN(Nomenclature, textInput),
                            (1, 'Alk|anes/enes/ynes'): Nomenclature.alkneNtoF(Nomenclature, textInput),
                            (0, 'Polyatomic Oxyanions'): Nomenclature.polyatomicOxyanionFtoN(Nomenclature, textInput),
                            (1, 'Polyatomic Oxyanions'): Nomenclature.polyatomicOxyanionNtoF(Nomenclature, textInput),
                            (0, 'Acids'): Nomenclature.acidFtoN(Nomenclature, textInput),
                            (1, 'Acids'): Nomenclature.acidNtoF(Nomenclature, textInput)
                            }
            if type(Nomenclature.polyatomicOxyanionFtoN(Nomenclature, textInput)) == list:
                self.funcMap.update({(0, 'Polyatomic Oxyanions'):
                                    Nomenclature.polyatomicOxyanionFtoN(Nomenclature, textInput)[self.nomSystemInt.get()]})
            self.ansLabelText.set(self.funcMap[(self.FNOption.get(), self.nomTypeVar.get())])
        self.update_idletasks()

    def updateEnterText(self, *args):
        if self.FNOption.get() == 0:
            self.enterText.set('Name:')
        elif self.FNOption.get() == 1:
            self.enterText.set('Formula:')
        self.update_idletasks()

gui = GUI()
gui.mainloop()