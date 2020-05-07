from Primer_Designer_Functions import SEARCH, PRIMER_DESIGN, IN_SILICO_PCR

keyword = "ORF7b coronavirus"

# searching
webpage = SEARCH.searcher(keyword)
page_content = SEARCH.scraper1(webpage)
results = SEARCH.links_results(page_content)
length = len(results)

if length > 1:

    SEARCH.list_creator(results)

else:

    # Design primers
    fasta_seq = PRIMER_DESIGN.get_fasta(page_content)
    dictionary = PRIMER_DESIGN.PRIMER3(fasta_seq)
    Primers_Tm_GC = PRIMER_DESIGN.primers_tm_gc(dictionary)
    Hairpins_calc = PRIMER_DESIGN.Hairpins(Primers_Tm_GC)
    Homodimers_calc = PRIMER_DESIGN.Homodimers(Hairpins_calc)
    Heterodimers_calc = PRIMER_DESIGN.Heterodimers(Homodimers_calc)

    # In_silico PCR
    primer_list = IN_SILICO_PCR.primers(Heterodimers_calc)
    pcr_products = IN_SILICO_PCR.in_silico_pcr(primer_list, fasta_seq)
    
    
    
    
    

