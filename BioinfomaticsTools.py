import urllib.request
from bs4 import BeautifulSoup
from selenium import webdriver

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

    def scraper2(website):

        # this require download geckodrive and move it to usr/local/bin
        driver = webdriver.Firefox()
        driver.get(website)
        content = driver.page_source
        soup = BeautifulSoup(content, 'html5lib')
        driver.close()

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

        def scraper(webpage):

            # this require download geckodrive and move it to usr/local/bin
            driver = webdriver.Firefox()
            driver.get(website)
            content = driver.page_source
            soup = BeautifulSoup(content, 'html5lib')
            driver.close()

            return soup

        link_fasta = 'https://www.ncbi.nlm.nih.gov'

        # get gene name
        title = page_content.find_all('h1', class_ = 'title')
        ID = title[0].span.text
        name = title[0].em.text
        gene_name = ID+' - '+name

        # get fasta link
        content_link = page_content.find_all('ol')
        olgene = content_link[0]
        gene_info = olgene.find_all('a')
        href_fasta = gene_info[1]['href']
        link_fasta += href_fasta
        link_fasta_content = scraper(link_fasta)
        seq = link_fasta_content.find_all('pre')
        seq_fasta = seq[0].text

        return seq_fasta
        

        
            
            

            

        
        

        
    

    

    
