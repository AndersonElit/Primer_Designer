from BioinfomaticsTools import SEARCH

# search gene

keyword = "ORF7b coronavirus"
webpage = SEARCH.searcher(keyword)
page_content = SEARCH.scraper1(webpage)

# extract results links
links = SEARCH.links_results(page_content)

# get gene names, links genbank, fasta and protein
data_content_links = SEARCH.links_genbank_fasta_protein(links)



