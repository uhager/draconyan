### Ulrike Hager, Apr 2010 ###

uInKeV = 931494.013


def read_nubtab(nuclName):
    """Read nuclide from nubtab.asc.

    Name format example: 17O; calculates mass in u, MeV and kg."""
    nubtab = open('nubtab03.asc', 'r')
    nucleus = {'name':nuclName, 'A':0, 'Z':0, 'ME':0, 'massU':0, 'massMeV':0}
    for line in nubtab:
        dataList = line.split()
        if dataList[2]==nuclName:
            nucleus['Z'] = int(round(float(dataList[1])/10.))
            nucleus['A'] = int(dataList[0])
            nucleus['ME'] = float(dataList[3])
            nucleus['massU'] = nucleus['A'] + nucleus['ME']/uInKeV
            nucleus['massMeV'] = (nucleus['A']*uInKeV + nucleus['ME'])/1000
            nucleus['massKg'] = nucleus['massU'] * 1.660538782e-27
            break
    else: return 0
    return nucleus

            
        
            
        
