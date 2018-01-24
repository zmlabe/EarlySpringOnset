"""
*Function converts given temperature units.*
"""

def tempunits(inT,inunits,outunits):
    """
    Parameters
    ----------
    inT : Temperature
    inunits : Celsius='C', Kelvin= 'K', Fahrenheit= 'F'
    outunits : Celsius='C', Kelvin= 'K', Fahrenheit= 'F'
    
    Returns
    ----------
    outT : Temperature with converted Units
    """
    if inunits =='C' and outunits =='F':
        outT = (9./5)*inT+32
    elif inunits =='F' and outunits =='C':
        outT = (inT-32)*(5./9)
    elif inunits =='K' and outunits =='F':
        outT = (9./5)*inT-459.67
    elif inunits =='F' and outunits =='K':
        outT = (inT+459.67)/(9./5)
    elif inunits == 'C' and outunits == 'K':
        outT = inT + 273.15
    elif inunits == 'K' and outunits == 'C':
        outT = inT - 273.15
    elif inunits == outunits:
        outT = inT
    return outT
