from Primer_Designer_Functions import SEARCH, PRIMER_DESIGN, IN_SILICO_PCR, PRIMER_ANALYSIS, COTIZER

def all_analysis():

    keyword = "orf1ab AND coronavirus AND human AND HKU1"

    # searching
    webpage = SEARCH.searcher(keyword)
    page_content = SEARCH.scraper1(webpage)
    # Design primers
    fasta_seq = PRIMER_DESIGN.get_fasta(page_content)
    dictionary = PRIMER_DESIGN.PRIMER3(fasta_seq)
    Primers_Tm_GC = PRIMER_DESIGN.primers_tm_gc(dictionary)

    # In_silico PCR
    pcr_products = IN_SILICO_PCR.in_silico_pcr(Primers_Tm_GC, fasta_seq)

    # oligoanalyzer
    thermo_analysis = PRIMER_ANALYSIS.oligoanalyzer(pcr_products)

    # cotizer
    IDT_prices = COTIZER.IDT(pcr_products)

    all_data = []
    count = 0

    while count <= 4:
        
        sets = [pcr_products[count], thermo_analysis[count], IDT_prices[count]]
        all_data.append(sets)
        count += 1

    return all_data
