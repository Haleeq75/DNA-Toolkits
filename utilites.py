def colored(seq):
    bcolor= {
        'A': '\0033[92m ',
        'C': '\0033[94m ',
        'T': '\0033[91m ',
        'U': '\0033[91m ',
        'G': '\0033[93m ',
        'reset': '\003[O;Om'
    }
    tempStr= ""

    for nuc in seq:
        if nuc in bcolor:
            tempStr+=bcolor[nuc]+nuc
        else:
            tempStr+=bcolor['reset']+nuc
    return tempStr + '\003[O;Om'