__author__ = 'Emerson Shoichet-Bartus'
__high_school__ = 'Upper Canada College'
__house__ = "Martland's"
__graduation_year__ = '2016'
__version__ = 1.0

import re

"""
This file contains all the chemistry nomenclature algorithms

'FtoN' signifies a method that converts a formula to a name
'NtoF' signifies a method that converts a name to a formula

The compounds that a method names is given explicitly in the function name

The commented out print functions in methods that print numbers are used for debugging

Planned expansions to this program are works in progress, existing incomplete functions can be seen at the bottom of the file
"""

class Nomenclature:


    errMsgF, errMsgN, convOne = 'Invalid formula', 'Invalid name', lambda num: '' if num == 1 else str(num)
    checkRegex = lambda regex, str: re.fullmatch(regex, str)

    dataFile = open('data.txt', 'r')
    elemDataset = {}
    dataFile.readline()  # get rid of legend line
    for line in dataFile.readlines():
        dataPieces = line.strip().split(', ')
        if '/' in dataPieces[0]:
            elemDataset.update({dataPieces[1].upper(): (tuple([spelling for spelling in dataPieces[0].split('/')]), int(dataPieces[2]), dataPieces[3])})
        elif '/' in dataPieces[2]:
            elemDataset.update({dataPieces[1].upper(): (dataPieces[0], tuple([int(n) for n in dataPieces[2].split('/')]), dataPieces[3])})
        else:
            elemDataset.update({dataPieces[1].upper(): (dataPieces[0], int(dataPieces[2]), dataPieces[3])})

    def alkneFtoN(self, formula):
        try:
            if type(formula) != str:
                raise TypeError('formula must be a character sequence')
            apology, formula, prefix = 'cannot name alk|anes/enes/ynes which exceed 20 carbon atoms', formula.upper(), ''
            prefixList1, prefixList2 = ['meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec'], \
                                       ['un', 'do', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona']
            if formula[0] != 'C':
                #print(0)
                return self.errMsgF
            if len(formula) == 6:
                if formula[3] != 'H':
                    #print(1)
                    return self.errMsgF
                Cnum, HnumRange = int(formula[1:3]), 4
            elif len(formula) == 3:
                if formula[1] != 'H':
                    #print(2)
                    return self.errMsgF
                Cnum, HnumRange, prefix = 1, 2, prefixList1[0]
            elif len(formula) < 6:
                HnumRange, Cnum = 3, int(formula[1])
            else:
                #print(3)
                return self.errMsgF
            Hnum = int(formula[HnumRange:])
            if Cnum > 20:
                return apology
            for i in range(len(prefixList1)):
                prefix = prefixList1[i] if Cnum == i+1 and Cnum <= 10 \
                    else prefixList2[i] + prefixList1[9] if Cnum == i+11 and 10 < Cnum < 20 else None
                if prefix is not None:
                    break
            if Cnum == 20:
                prefix = 'eicos'
            if prefix == '':
                #print(4)
                return self.errMsgF
            if Hnum not in {2*Cnum + 2: 'a', 2*Cnum: 'e', 2*Cnum - 2: 'y'}:
                #print(5)
                return self.errMsgF
            if Cnum == 1 and Hnum != 4:
                #print(6)
                return self.errMsgF
            return prefix + {2*Cnum + 2: 'a', 2*Cnum: 'e', 2*Cnum - 2: 'y'}[Hnum] + 'ne'
            # suffix is suffixDict[Hnum]
        except Exception:
            #print('error')
            return self.errMsgF

    # TODO fully implement regex
    def alkneNtoF(self, name):
        try:
            if type(name) != str:
                raise TypeError('hydrocarbon name must be a word')

            apology, name, Cnum = 'alk|ane/ene/yne exceeds 20 carbon atoms', name.lower(), 0
            prefix, suffixLetter = name[:len(name) - 3], name[len(name) - 3].lower()
            prefixDict1, prefixDict2 = {'meth': 1, 'eth': 2, 'prop': 3, 'but': 4, 'pent': 5, 'hex': 6, 'hept': 7,
                                        'oct': 8, 'non': 9, 'dec': 10}, \
                                       {'un': 1, 'do': 2, 'tri': 3, 'tetra': 4, 'penta': 5, 'hexa': 6, 'hepta': 7,
                                        'octa': 8, 'nona': 9}
            regex = re.compile('(((un|do|tri|tetra|penta|hexa|hepta|octa|nona)?dec)'
                               '|(meth|eth|prop|but|pent|hex|hept|oct|non|dec|eicos))([aey]ne)', re.VERBOSE)
            if not self.checkRegex(regex, name):
                return self.errMsgN
            if name[len(name)-2:] != 'ne':
                #print(0)
                return self.errMsgN
            if prefix in prefixDict1:
                Cnum = prefixDict1[prefix]
            else:
                for p in prefixDict2:
                    Cnum = prefixDict2[p] + 10 if prefix == p + 'dec' else 20 if prefix == 'eicos' else None
                    if Cnum is not None:
                        break
            if Cnum == 0:
                #print(1)
                return self.errMsgN + ' or ' + apology
            if suffixLetter not in 'aey':
                #print(2)
                return self.errMsgN
            Hnum = {'a': 2*Cnum + 2, 'e': 2*Cnum, 'y': 2*Cnum - 2}[suffixLetter]
            if Cnum == 1 and Hnum != 4:
                #print(3)
                return self.errMsgN
            else:
                return 'C' + self.convOne(Cnum) + 'H' + str(Hnum)
        except Exception:
            #print('error')
            return self.errMsgN

    # TODO regex testing
    def polyatomicOxyanionFtoN(self, formula):
        try:
            if type(formula) != str:
                raise TypeError('formula must be a character sequence')
            formula, numOsIndex = formula.upper().replace(' ', ''), 1
            dataset1 = {'P': ('phosph', 4, -3), 'S': ('sulf', 4, -2), 'N': ('nitr', 3, -1), 'C': ('carbon', 3, -2),
                        'I': ('iod', 3, -1), 'F': ('fluor', 3, -1), 'B': ('bor', 3, -3)}
            dataset2 = {'BR': ('brom', 3, -1), 'CL': ('chlor', 3, -1), 'AS': ('arsen', 4, -3), 'SE': ('selen', 4, -2),
                        'TC': ('technet', 3, -1), 'RE': ('rhen', 3, -2), 'CR': ('chrom', 4, -1)}
            romanNumerals = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII']
            regex = re.compile('([PS]|As|Se|Cr)O[2-5]|([FINCB]|Br|Cl|Tc|Re)O[1-4]?')
            if not self.checkRegex(regex):
                print('regex error')
                return self.errMsgF
            '''keyFlag = False
            for k in dataset2.keys():
                if k in formula:
                    keyFlag = True
                    name = dataset2[k][0]
            if not keyFlag:
                for k in dataset1.keys():
                    if k in formula:
                        name = dataset1[k][0]'''
            stan = -5
            if len(formula) == 4:
                key = formula[0:2]
                name = dataset2[key][0]
                stan = dataset2[key][1]
                mainElCharge = dataset2[key][2] + 2 * stan
                numOs = int(formula[3])
            elif len(formula) == 3:
                if formula[2].isdigit():
                    key = formula[0]
                    name = dataset1[key][0]
                    stan = dataset1[key][1]
                    mainElCharge = dataset1[key][2] + 2*stan
                    numOs = int(formula[2])
                elif not formula[2].isdigit():
                    key = formula[0:2]
                    name = dataset2[key][0]
                    stan = dataset2[key][1]
                    mainElCharge = dataset2[key][2] + 2*stan
                    numOs = 1
            prefSufDict = {stan+1: ('per', 'ate'), stan: ('', 'ate'), stan-1: ('', 'ite'), stan-2: ('hypo', 'ite')}
            '''if formula[0] in dataset1 and formula[0:2] not in dataset2: #if One (key in dataset)
                key, One = formula[0], True
                name, stan = dataset1[key][0], dataset1[key][1]
            elif formula[0:2] in dataset2: #if Two (key in dataset)
                key, One = formula[0:2], False
                numOsIndex, name, stan = 2, dataset2[key][0], dataset2[key][1]
            numOs = 1
            if formula[numOsIndex] != 'O':
                #print(0)
                return errMsgF
            if len(formula) == numOsIndex+2 and formula[numOsIndex+1].isdigit():
                    numOs = int(formula[numOsIndex+1])

            if numOs not in prefSufDict:
                #print(1)
                return errMsgF
            #mainElCharge = dataset1[key][2] + 2*numOs if One else dataset2[key][2] + 2*numOs if not One else None  # math to compute charge on central atom
            if mainElCharge > 8:
                #print(2)
                return errMsgF'''
            return [prefSufDict[numOs][0] + name + prefSufDict[numOs][1], name + 'ate (' + romanNumerals[mainElCharge-1] + ')']
        except Exception:
            #print('error')
            return self.errMsgF


    #INCOMPLETE regex missing roman numeral multiples
    # (hypo)?((phosph)|(sulf)|(sulph)|(nitr)|(carbon)|(iod)|(fluor)|(brom)|(chlor)|(bor)|(arsen)|(selen)|(technet)|(rhen)|(chrom))ite|((phosph)|(sulf)|(sulph)|(nitr)|(carbon)|(iod)|(fluor)|(brom)|(chlor)|(bor)|(arsen)|(selen)|(technet)|(rhen)|(chrom))(ate)( [(]((i{1,3})|(v(i{1,2})?)|(i?v))[)])?
    def polyatomicOxyanionNtoF(self, name):
        try:
            if type(name) != str:
                raise TypeError('name of polyatomic oxyanion must be a word')
            name = name.lower().replace(' ', '')
            num, suffix, convCharge = None, name[len(name)-3:] if len(name.split()) == 1 else 'ate' if len(name.split()) == 2 else None, \
                                      lambda s: -1 if s == '-' else -int(s[0])
            dataset = {'phosph': ('P', 4, '3-'), 'sulf': ('S', 4, '2-'), 'sulph': ('S', 4, '2-'), 'nitr': ('N', 3, '-'),
                       'carbon': ('C', 3, '2-'), 'iod': ('I', 3, '-'), 'fluor': ('F', 3, '-'), 'brom': ('Br', 3, '-'),
                       'chlor': ('Cl', 3, '-'), 'bor': ('B', 3, '3-'), 'arsen': ('As', 4, '3-'), 'selen': ('Se', 4, '2-'),
                       'technet': ('Tc', 3, '-'), 'rhen': ('Re', 3, '-'), 'chrom': ('Cr', 4, '2-')}
            prefix, start = 'hypo' if name[:4] == 'hypo' else 'per' if name[:3] == 'per' else '', 4 if name[:4] == 'hypo' \
                else 3 if name[:3] == 'per' else 0
            for el in dataset:
                polyName = el if name[start:len(name)-3] == el else None
                if polyName is None:
                    try:
                        polyName = el if name[start:name.index('(')-3] == el else None
                    except ValueError:
                        pass
                if polyName is not None:
                    break
            if polyName is None:
                #print(0)
                return self.errMsgN
            if '(' in name and ')' in name: # and name.split()[1][len(name.split()[1]-1]
                if name[:4] == 'hypo' or name[:3] == 'per' or name[len(polyName):name.index('(')] == 'ite':
                    #print(1)
                    return self.errMsgN
                romanNumConvChart = {'i': 1, 'ii': 2, 'iii': 3, 'iv': 4, 'v': 5, 'vi': 6, 'vii': 7, 'viii': 8}
                try:
                    mainElStan = .5*(romanNumConvChart[name[len(polyName)+4:len(name)-1]] - convCharge(dataset[polyName][2])) #last term is mainElCharge
                except KeyError:
                    #print(2)
                    return self.errMsgN
                if not mainElStan.is_integer():
                    #print(3)
                    return self.errMsgN
                return dataset[polyName][0] + 'O' + self.convOne(int(mainElStan)) + ' ' + dataset[polyName][2]
            if suffix != 'ate' and suffix != 'ite':
                #print(4)
                return self.errMsgN
            if (prefix, suffix) in [('hypo', 'ate'), ('per', 'ite')]:  # all BAD prefix & suffix combinations
                #print(5)
                return self.errMsgN
            stan = dataset[polyName][1]
            mainElStan = {('', 'ate'): stan, ('per', 'ate'): stan+1, ('', 'ite'): stan-1, ('hypo', 'ite'): stan-2}[(prefix, suffix)]
            if mainElStan is None:
                #print(6)
                return self.errMsgN
            return dataset[polyName][0] + 'O' + self.convOne(mainElStan) + ' ' + dataset[polyName][2]
            # return (letter symbol for element) + 'O' + (number of Os) + ' ' + (charge on oxyanion, which is in elementDict)
        except Exception as e:
            #print('error')
            e.with_traceback(0)
            return self.errMsgN

    # regex = H(([INF]|Br|Cl|Tc|Re)(O[1-4]?)?|2((S|Se|Cr)(O[2-5])|CO[1-4])|3((P|As)O[2-5]|BO[1-4]))
    def acidFtoN(self, formula):
        try:
            global start, end, Obool
            if type(formula) != str:
                raise TypeError('formula must be a character sequence')
            formula, Two, mainElSymbol, numOs, hydro = formula.upper(), False, None, 0, False
            dataset = {'P': ('phosphor', 4, -3), 'S': ('sulfur', 4, -2), 'N': ('nitr', 3, -1), 'C': ('carbon', 3, -2),
                        'I': ('iod', 3, -1), 'F': ('fluor', 3, -1), 'BR': ('brom', 3, -1), 'CL': ('chlor', 3, -1),
                        'B': ('bor', 3, -3), 'AS': ('arsen', 4, -3), 'SE': ('selen', 4, -2), 'TC': ('technet', 3, -1),
                        'RE': ('rhen', 3, -1), 'CR': ('chrom', 4, -2)}
            if formula[0] != 'H':
                #print(0)
                return self.errMsgF
            HnumBool = False if not formula[1].isdigit() else True
            start, end = 2 if HnumBool else 1, 4 if HnumBool else 3
            if formula[start:end] in dataset:
                mainElSymbol, Obool = formula[start:end], len(formula) >= end+1 and formula[end] == 'O'
            elif formula[start:end] not in dataset:
                if formula[1] not in dataset and not HnumBool:
                    #print(1)
                    return self.errMsgF
                mainElSymbol, Obool = formula[start], formula[start+1] == 'O'
            if mainElSymbol is None:
                #print(2)
                return self.errMsgF
            if HnumBool and int(formula[1]) != abs(dataset[mainElSymbol][2]):
                    #print(3)
                    return self.errMsgF
            if not HnumBool and dataset[mainElSymbol][2] != -1:
                    #print(4)
                    return self.errMsgF
            if mainElSymbol not in dataset:
                #print(5)
                return self.errMsgF
            elif dataset[mainElSymbol][2] == -1 and formula[1] == '1': #TODO prohibits C subscript from being '1' in formula
                #print(6)
                return self.errMsgF
            if not Obool:
                if formula[len(formula)-1].isdigit():
                    #print(7)
                    return self.errMsgF
                return 'hydro' + dataset[mainElSymbol][0] + 'ic acid'
            elif Obool:
                if not formula[len(formula)-1].isdigit() and dataset[mainElSymbol][1] == 3: #TODO Added second part of and clause
                    print(8)
                    return self.errMsgF
                numOs = int(formula[len(formula)-1])
                if dataset[mainElSymbol][1] == 3:
                    return 'hypo' + dataset[mainElSymbol][0] + 'ous acid'
                stan = int(dataset[mainElSymbol][1])
                prefSufs = {stan: ('', 'ic'), stan+1: ('per', 'ic'), stan-1: ('', 'ous'), stan-2: ('hypo', 'ous')}
                if numOs not in prefSufs.keys():
                    print(9)
                    return self.errMsgF
                elif formula[len(formula)-2] != 'O':
                    print(10)
                    return self.errMsgF
                return prefSufs[numOs][0] + dataset[mainElSymbol][0] + prefSufs[numOs][1] + ' acid'
                # return prefix + (main element name) + suffix + ' acid'
        except Exception:
            print('error')
            return self.errMsgF

    def acidNtoF(self, name):
        try:
            if type(name) != str:
                raise TypeError('acid name must be a word')
            name = name.replace(' ', '').lower()
            dataset = {'phosphor': ('P', 4, -3), 'sulfur': ('S', 4, -2), 'sulphur': ('S', 4, -2), 'nitr': ('N', 3, -1),
                       'carbon': ('C', 3, -2), 'iod': ('I', 3, -1), 'fluor': ('F', 3, -1), 'brom': ('Br', 3, -1),
                        'chlor': ('Cl', 3, -1), 'bor': ('B', 3, -3), 'arsen': ('As', 4, -3), 'selen': ('Se', 4, -2),
                       'technet': ('Tc', 3, -1), 'rhen': ('Re', 3, -1), 'chrom': ('Cr', 4, -2)}
            if name[len(name)-4:] != 'acid':
                return 'Must have "acid" at end of name'
            name, prefSufs = name[:len(name)-4], {'hydro': True, 'hypo': True, 'per': True, 'ic': False, 'ous': False} #prefixes are True, suffixes are False
            for pOrS, boolean in zip(prefSufs.keys(), prefSufs.values()): #pOrS is prefixOrSuffix
                prefix = pOrS if boolean and pOrS == name[:len(pOrS)] else ''
                start = len(pOrS) if boolean and pOrS == name[:len(pOrS)] else 0
                if start is not 0:
                    break
            for pOrS, boolean in zip(prefSufs.keys(), prefSufs.values()): #pOrS is prefixOrSuffix
                suffix = pOrS if not boolean and pOrS == name[len(name)-len(pOrS):] else ''
                end = len(pOrS) if not boolean and pOrS == name[len(name)-len(pOrS):] else None
                if end is not None:
                    break
            if end is None:
                #print(0)
                return self.errMsgN
            if (prefix, suffix) in [('hydro', 'ous'), ('hypo', 'ic'), ('per', 'ous')]: #all BAD prefix & suffix combinations
                #print(1)
                return self.errMsgN
            mainEl = name[start:len(name)-end]
            if mainEl not in dataset.keys():
                #print(2)
                return self.errMsgN
            returnStan = abs(dataset[mainEl][2])
            if prefix == 'hydro':
                return 'H' + self.convOne(returnStan) + dataset[mainEl][0]
            stan = dataset[mainEl][1]
            returnNumOs = {('', 'ic'): stan, ('per', 'ic'): stan+1, ('', 'ous'): stan-1, ('hypo', 'ous'): stan-2}[(prefix, suffix)]
            return 'H' + self.convOne(returnStan) + dataset[mainEl][0] + 'O' + self.convOne(returnNumOs)
        except Exception:
            print('error')
            return self.errMsgN


    #TODO work in progress
    def saltsAndCovalentsFtoN(self, formula):
        try:
            covalentPrefixMap = {1: 'mono', 2: 'di', 3: 'tri', 4: 'tetra', 5: 'penta', 6: 'hexa', 7: 'hepta', 8: 'octa',
                         9: 'nona', 10: 'deca'}
            romanNumerals = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII']
            name1, name2, formula, name2start = None, None, formula.upper(), 0
            name1end = 1 if formula[0] in self.elemDataset and formula[0:2] not in self.elemDataset else 2 if formula[0:2] in self.elemDataset else 0
            #print(formula[0] in elemDataset and formula[0:2] not in elemDataset)
            if name1end == 0:
                #print(0)
                return self.errMsgF
            name1key = formula[0:name1end]
            name1 = self.elemDataset[name1key][0]
            if name1 is None:
                #print(1)
                return self.errMsgF
            assert name1end != 0
            if formula[name1end:name1end+2].isdigit():
                name2start = name1end + 2
                name2digit = int(formula[name1end:name1end+2])
            elif formula[name1end].isdigit():
                name2start = name1end + 1
                name1digit = int(formula[name1end])
            elif formula[name1end].isalpha():
                name2start = name1end
                name1digit = 1
            else:
                #print(2)
                return self.errMsgF
            if formula[name2start] in self.elemDataset and formula[name2start:name2start + 2] not in self.elemDataset:
                name2end = name2start + 1
            elif formula[name2start:name2start + 2] in self.elemDataset:
                name2end = name2start + 2
            else:
                #print(3)
                return self.errMsgF
            name2key = formula[name2start:name2end]
            name2 = self.elemDataset[name2key][0]
            if (len(formula) == name2end or len(formula) == name2end + 1) and formula[name2end:].isdigit():
                name2digit = int(formula[name2end:])
            if name2 is None:
                #print(3)
                return self.errMsgF
            if self.elemDataset[name1key][3] == self.elemDataset[name2key][3] == 'nonmetal':
                name1prefix = ''
                if name2digit > 1:
                    name1prefix = covalentPrefixMap[name1digit]
                return name1prefix + name1 + ' ' + covalentPrefixMap[name2digit] + name2
            elif (self.elemDataset[name1key][3] == 'nonmetal' and self.elemDataset[name2key][3] == 'metal') or (self.elemDataset[name1key][3] == 'metal' and self.elemDataset[name2key][3] == 'nonmetal'):
                return name1 #TODO
            #if polyatomicOxyanionFtoN(formula[name2start:]) != 'Invalid formula':


        #if type(dataset[key][1]) == int:
        #    if not formula[end].isdigit() and abs(dataset[key][1]) > 1:
        #        print(4)
        #        return errMsgF
            return name1 + ' ' + name2
        except Exception:
            print('error')
            return self.errMsgF


    # TODO work in progress
    def covalentsFtoN(self, formula):
        try:
            covalentSet, formula = {}, formula.upper()
            prefixMap = {1: 'mono', 2: 'di', 3: 'tri', 4: 'tetra', 5: 'penta', 6: 'hexa', 7: 'hepta', 8: 'octa', 9: 'nona', 10: 'deca'}
            for elem in self.elemDataset:
                if self.elemDataset[elem][2] == 'nonmetal':
                    if type(self.elemDataset[elem][0]) == tuple:
                        covalentSet.update({elem: (self.elemDataset[elem][0][0], self.elemDataset[elem][1])})
                    else:
                        covalentSet.update({elem: (self.elemDataset[elem][0], self.elemDataset[elem][1])})
            mainEl1 = formula[0] if formula[0] in covalentSet and formula[0:2] not in covalentSet else formula[0:2] if formula[0:2] in covalentSet else None
            if mainEl1 is None:
                #print(0)
                return self.errMsgF
            mainEl1Prefix = prefixMap[int(formula[len(mainEl1)])] if formula[len(mainEl1)].isdigit() else ''
            mainEl2startIndex = len(mainEl1) + 1
            #mainEl2 = formula[] #TODO######################################################################
            if mainEl1Prefix == '':
                mainEl2 = formula[len(mainEl1)] if formula[len(mainEl1)] in covalentSet and formula[len(mainEl1):len(mainEl1)+2] not in covalentSet \
                else formula[len(mainEl1):len(mainEl1)+2] if formula[len(mainEl1):len(mainEl1)+2] else None
                mainEl2startIndex = len(mainEl1)
            if mainEl2 is None:
                #print(1)
                return self.errMsgF
            mainEl2Prefix = prefixMap[int(formula[len(mainEl2) + mainEl2startIndex].isdigit())] if formula[len(mainEl2) + mainEl2startIndex].isdigit() else 'mono'
            return mainEl1Prefix + mainEl1 + ' ' + mainEl2Prefix + mainEl2
        except Exception as e:
            e.with_traceback(None)
            #print('error')
            return self.errMsgF


#print(Nomenclature.covalentsFtoN('Br4O2')) # testing the new covalent nomenclature function


#Used for testing
#while 1 + 1 == 2: # infinite loop
    #print(Nomenclature.alkneFtoN(input()))
    #print(Nomenclature.alkneNtoF(input()))
    #print(Nomenclature.polyatomicOxyanionFtoN(input()))
    #print(Nomenclature.polyatomicOxyanionNtoF(input()))
    #print(Nomenclature.acidFtoN(input()))
    #print(Nomenclature.acidNtoF(input()))
#print(Nomenclature.saltsAndCovalentsFtoN(input()))
    #print(Nomenclature.saltsAndCovalentsNtoF(input())) #TODO future
    #print(Nomenclature.covalentsFtoN(input()))