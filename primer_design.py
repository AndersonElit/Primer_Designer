from Primer_Designer_Functions import SEARCH, PRIMER_DESIGN, IN_SILICO_PCR, PRIMER_ANALYSIS

keyword = "orf1ab AND coronavirus AND human AND HKU1"

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

    # In_silico PCR
    pcr_products = IN_SILICO_PCR.in_silico_pcr(Primers_Tm_GC, fasta_seq)

    # oligoanalyzer
    thermo_analysis = PRIMER_ANALYSIS.oligoanalyzer(pcr_products)

    
    
    
    
    

