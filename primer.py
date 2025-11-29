# Compute GC content and Tm(melting temperature) of DNA sequence 

import math

# GC CONTENT 
def gc_content(seq):
    seq = seq.upper()
    g = seq.count('G')
    c = seq.count('C')
    return (g+c)/len(seq) *100

# Tm calculation 
def tm_value(seq):
    seq=seq.upper()
    a = seq.count('A')
    g = seq.count('G')
    c = seq.count('C')
    t = seq.count('T')
    return 2*(a+t) + 4*(g+c)

# Primer design Function 

def design_primers(dna_seq, primer_len=20, tm_range=(50,65), gc_range=(40,60)):
    dna_seq = dna_seq.upper()
    primers = []
    
    for i in range(len(dna_seq) - primer_len+1):
        primer= dna_seq[i:i+primer_len]
        tm= tm_value(primer)
        gc= gc_content(primer)
        
        if tm_range[0] <= tm <=tm_range[1] and gc_range[0]<= gc<=gc_range[1]:
            primers.append({
                'Primer': primer, 
                'Start_Pos' : i+1, 
                'Tm' : tm, 
                'GC%': round(gc,2)     
            })
    return primers 

if __name__ == "__main__":
    dna_sequence = input("Enter the DNA sequence: ")

    primers = design_primers(
        dna_sequence,
        primer_len=20,
        gc_range=(45, 65),
        tm_range=(55, 65)
    )

    if primers:
        print(f"\nFound {len(primers)} valid primers:\n")
        for p in primers:
            print(f"Pos {p['Start_Pos']} | Tm: {p['Tm']}°C | GC: {p['GC%']}% | Primer: {p['Primer']}")
    else:
        print("No primers matched the criteria.")

    

