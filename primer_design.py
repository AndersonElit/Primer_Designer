from BioinfomaticsTools import SEARCH, PRIMER_DESIGN

keyword = "ORF7b coronavirus homo sapiens"

# searching
webpage = SEARCH.searcher(keyword)
page_content = SEARCH.scraper1(webpage)
links = SEARCH.links_results(page_content)
length = len(links)

if length == 0:

    print("multiple results, specify more its search")

else:

    fasta_seq = PRIMER_DESIGN.get_fasta(page_content)
    dictionary = PRIMER_DESIGN.PRIMER3(fasta_seq)
    Primers_Tm_GC = PRIMER_DESIGN.primers_tm_gc(dictionary)
    Hairpins_calc = PRIMER_DESIGN.Hairpins(Primers_Tm_GC)
    Homodimers_calc = PRIMER_DESIGN.Homodimers(Hairpins_calc)
    Heterodimers_calc = PRIMER_DESIGN.Heterodimers(Homodimers_calc)
    
    
    
    

