import urllib.request
from bs4 import BeautifulSoup
import primer3
from bash import bash

class SEARCH:
    
    def searcher(keyword):

        webpage = "https://www.ncbi.nlm.nih.gov/gene/"
        words = keyword.split()
        length = len(words) -1
        search = ""
        cont = 0

        while cont <= length:

            if cont == length:
                conc = words[cont]
                search += conc
            else:
                conc = words[cont] + "+AND+"
                search += conc

            cont += 1

        tag = "?term="
        webpage += tag
        webpage += search

        return webpage

    def scraper1(webpage):

        content = urllib.request.urlopen(webpage)
        soup = BeautifulSoup(content, 'html5lib')
        
        return soup

    def links_results(page_content):

        table_results = page_content.find_all('tr', class_ = 'rprt')
        # get links
        links = []

        for item in table_results:

            link = 'https://www.ncbi.nlm.nih.gov'
            tds = item.find_all('td')
            href = tds[0].a['href']
            link += href
            links.append(link)

        return links

    def links_genbank_fasta_protein(links):

        def scraper(webpage):
            
            content = urllib.request.urlopen(webpage)
            soup = BeautifulSoup(content, 'html5lib')

            return soup

        data_set = []
        length = len(links)
        cont = 1

        for link in links:

            inc_link = 'https://www.ncbi.nlm.nih.gov'
            page_content = scraper(link)
            print(str(cont)+' de '+str(length))
            cont += 1
            
            # get gene name
            title = page_content.find_all('h1', class_ = 'title')
            ID = title[0].span.text
            name = title[0].em.text
            gene_name = ID+' - '+name

            # get other info
            content_link = page_content.find_all('ol')

            # get link genbank and fasta
            olgene = content_link[0]
            gene_info = olgene.find_all('a')
            href_genbank = gene_info[0]['href']
            href_fasta = gene_info[1]['href']
            link_genbank = inc_link + href_genbank
            link_fasta = inc_link + href_fasta

            # get protein link
            olprotein = content_link[1]
            protein_info = olprotein.find_all('a')
            href_protein = protein_info[0]['href']
            link_protein = inc_link + href_protein

            data = [gene_name, link_genbank, link_fasta, link_protein]
            data_set.append(data)

        return data_set

class PRIMER_DESIGN:

    def get_fasta(page_content):

        '''
        link_fasta = 'https://www.ncbi.nlm.nih.gov'
        '''

        # get gene name
        title = page_content.find_all('h1', class_ = 'title')
        ID = title[0].span.text
        name = title[0].em.text
        gene_name = ID+' - '+name

        # get sequence id and region
        geneid = page_content.find_all('dl', class_ = 'dl-chr-info')
        id_region = geneid[0].dd.text
        id_region_list = id_region.split(' ')
        ID = id_region_list[0]
        region = id_region_list[1]
        regionc1 = region.replace('(', '')
        regionc2 = regionc1.replace(')', '')
        region_list = regionc2.split('..')

        # extract sequence form nucleotide database
        command = 'esearch -db nucleotide -query "' + ID + '" | efetch -format fasta'
        seq = str(bash(command))
        title_seq = seq.split('\n')
        length = len(title_seq) - 1

        cont = 1
        seqt = ''

        while cont <= length:

            fragment = title_seq[cont]
            seqt += fragment
            cont += 1

        flank1 = int(region_list[0]) - 1
        flank2 = int(region_list[1])

        gene = seqt[flank1:flank2]

        return gene

        '''
        # get fasta link
        content_link = page_content.find_all('ol')
        olgene = content_link[0]
        gene_info = olgene.find_all('a')
        href_fasta = gene_info[1]['href']
        link_fasta += href_fasta
        link_fasta_content = scraper(link_fasta)
        seq = link_fasta_content.find_all('pre')
        seq_fasta = seq[0].text
        seq_fastac = seq_fasta.replace('\n\n', '')

        return seq_fastac
        '''

    def PRIMER3():

        primers = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': 'NC_045512.2',
                'SEQUENCE_TEMPLATE': 'ATGATTGAACTTTCATTAATTGACTTCTATTTGTGCTTTTTAGCCTTTCTGCTATTCCTTGTTTTAATTATGCTTATTATCTTTTGGTTCTCACTTGAACTGCAAGATCATAATGAAACTTGTCACGCCTAA'
            },
            {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_OPT_TM': 60.0,
                'PRIMER_MIN_TM': 57.0,
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_MIN_GC': 20.0,
                'PRIMER_MAX_GC': 80.0,
                'PRIMER_MAX_POLY_X': 100,
                'PRIMER_INTERNAL_MAX_POLY_X': 100,
                'PRIMER_SALT_MONOVALENT': 50.0,
                'PRIMER_DNA_CONC': 50.0,
                'PRIMER_MAX_NS_ACCEPTED': 0,
                'PRIMER_MAX_SELF_ANY': 12,
                'PRIMER_MAX_SELF_END': 8,
                'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                'PRIMER_PAIR_MAX_COMPL_END': 8,
                'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]]
            })

        return primers


        
        

        
            
            

            

        
        

        
    

    

    
