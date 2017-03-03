from operator import itemgetter

inputfile = open("nist_v41_atom_data.txt")

atoms = []
new_atom = {}
for line in inputfile:
    fields = [x.strip() for x in line.split('=')]
        
    if fields[0] == 'Isotopic Composition' or fields[0] == 'Relative Atomic Mass':
        my_str = fields[1].split('(')[0]
            
        if my_str != '':
            fields[1] = float(my_str)
        else:
            fields[1] = "NAN"
        
        new_atom[fields[0]] = fields[1]
        
    elif fields[0] == "Notes":
        atoms.append(new_atom)
        new_atom = {}
        
    elif fields[0] == "Atomic Number" or fields[0] == "Mass Number":
        new_atom[fields[0]] = int(fields[1])
    
    elif fields[0] == "Atomic Symbol":
        new_atom[fields[0]] = fields[1]
        

atoms.sort(key=itemgetter("Atomic Number"))

for atom in atoms:
    if(atom['Atomic Number'] != 1):
        atom['name'] = str(atom['Mass Number']) + atom['Atomic Symbol']
    else:
        atom['name'] = atom['Atomic Symbol']


for i in range(0, len(atoms)):
    atom = atoms[i]
    if(i != len(atoms) - 1):
        print("{" + str(atom['Atomic Number']) + ", {" + "\"" + atom['name'] + "\"" + ", " + str(atom['Mass Number']) +
              ", " + str(atom['Relative Atomic Mass']) + ", " + str(atom['Isotopic Composition']) + "}},")
    else: 
        print("{" + str(atom['Atomic Number']) + ", {" + "\"" + atom['name']+ "\"" + ", " + str(atom['Mass Number']) +
              ", " + str(atom['Relative Atomic Mass']) + ", " + str(atom['Isotopic Composition']) + "}}")
         
