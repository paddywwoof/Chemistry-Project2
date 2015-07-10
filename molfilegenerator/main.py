import numpy as np
import random
from interactionmanager import InteractionManager

singlebond = [0,100]
doublebond = [100,160]
oxygen = [160,222]

def main():
    interaction_manager = InteractionManager(1000, 'interactions/cucumber.np')


    im = interaction_manager.interaction_matrix
    atoms = interaction_manager.atom_types
    shift_data = interaction_manager.shift_data

    bonds = getbonds(im,shift_data)

    #Add Oxygens
    for shift_index,shift_value in enumerate(shift_data):
        if between(shift_value, oxygen):
            atoms.append("O")
            bonds.append([shift_index+1,len(atoms),2])

    header = """Molecule Name
 Additional Information

 %s %s  0  0  0  0  0  0  0  0999 V2000
""" % (len(atoms), len(bonds))
    footer = "M  END"
    line = "    %s    %s    %s %s "+"  0"*12
    file = ""
    file += header
    for atom in atoms:
        i1 = str(random.uniform(0,1))[:6]
        i2 = str(random.uniform(0,1))[:6]
        i3 = str(random.uniform(0,1))[:6]
        file += line%(i1, i2, i3,atom)+"\n"
    bondline = " %s %s  %s  0  0  0  0"

    for bond in bonds:
        pass
        a1 = " "*(2-len(str(bond[0])))+str(bond[0])
        a2 = " "*(2-len(str(bond[1])))+str(bond[1])
        a3 = bond[2]
        b = bondline%(a1,a2,a3)
        file += (b+"\n")

    file += footer

    outfile = open('testing.mol','w+')
    outfile.write(file)
    outfile.close()
    print(file)

def between(x, interval):
    if interval[0] <= x and x < interval[1]:
        return True
    else:
        return False

def getbonds(im,shift_data):
    bonds = []
    for i in range(len(im)):
        for j in range(len(im)):
            if i < j:
                if im[i][j] in [2,4]:
                    if between(shift_data[i], singlebond) and between(shift_data[j], singlebond):
                        btype=1
                    if between(shift_data[i], doublebond) and between(shift_data[j], doublebond):
                        btype=2
                    else:
                        btype=1

                    bonds.append([i+1, j+1,btype])
    return bonds

main()