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
    primers = PRIMER_DESIGN.PRIMER3()
    
    

