import urllib.request
from bs4 import BeautifulSoup
import primer3
from bash import bash
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
import time

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
        data = []

        for item in table_results:

            # get link
            link = 'https://www.ncbi.nlm.nih.gov'
            tds = item.find_all('td')
            href = tds[0].a['href']
            link += href

            # get title_gene
            title_gene = tds[1].text

            # get gene_id
            gene_id = href[6:]

            # put data together
            data_set = [title_gene, gene_id, link]
            data.append(data_set)

        return data

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
'''
    def list_creator(results):

        # open home.html
        inf = open('templates/reference.html', 'r')
        soup = BeautifulSoup(inf)
        cont = 1

        for item in results:

            tr = BeautifulSoup("<tr></tr>")
            tr_tag = tr.tr

            # item num
            th = BeautifulSoup("<th scope='row'></th>")
            th_tag = th.th
            th_tag.string = str(cont)
            tr_tag.append(th_tag)

            # gene name
            gene_name = item[0]
            td_gene = BeautifulSoup("<td></td>")
            td_gene_tag = td_gene.td
            td_gene_tag.string = gene_name
            tr_tag.append(td_gene_tag)

            # gene id
            ID = item[1]
            link = item[2]
            td_id = BeautifulSoup("<td></td>")
            td_id_tag = td_id.td
            a = "<a href='" + link + "' target='_blank'>" + ID + "</a>"
            a_id = BeautifulSoup(a)
            a_id_tag = a_id.a
            td_id_tag.append(a_id_tag)
            tr_tag.append(td_id_tag)

            # button
            td_button = BeautifulSoup("<td></td>")
            td_button_tag = td_button.td
            button = BeautifulSoup("<button type='submit' class='btn btn-dark mr-2 font-weight-bolder wordstyle1 btn-sm'>Generar primers</button>")
            button_tag = button.button
            td_button_tag.append(button_tag)
            tr_tag.append(td_button_tag)

            # insert main html
            soup.tbody.append(tr_tag)

            cont += 1

        # create results.html e insertar html modificado
        new_html = open('templates/results.html', 'w+')
        new_html.write(str(soup))
'''      

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

        '''
        # extract sequence form nucleotide database, this is since linux terminal
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
        '''

        #extract sequence with biopython
        seq = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text")
        record = SeqIO.read(seq, "fasta")
        seq.close()
        sequence = str(record.seq)
        
        flank1 = int(region_list[0]) - 1
        flank2 = int(region_list[1])

        gene = sequence[flank1:flank2]

        # build fasta file
        all_record = str(record)
        title_seq = all_record.split('\n')
        description = title_seq[2]
        split = description.split(': ')
        seq_t = split[1]
        fasta_region = '>' + seq_t + '\n' + gene

        # create .fasta file with sequences

        f= open("seq.fasta","w+")
        f.write(fasta_region)
        f.close()

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

    def PRIMER3(sequence):

        primers = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': 'NC_045512.2',
                'SEQUENCE_TEMPLATE': sequence
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

    def primers_tm_gc(dictionary):

        primers_num = dictionary.get('PRIMER_LEFT_NUM_RETURNED') - 1
        cont = 0
        data_set = []

        while cont <= primers_num:

            # promer sequences
            left = 'PRIMER_LEFT_' + str(cont) + '_SEQUENCE'
            right = 'PRIMER_RIGHT_' + str(cont) + '_SEQUENCE'
            left_primer = dictionary.get(left)
            right_primer = dictionary.get(right)
            primer_pair = [left_primer, right_primer]
            
            # Tm
            Tml = 'PRIMER_LEFT_' + str(cont) + '_TM'
            Tmr = 'PRIMER_RIGHT_' + str(cont) + '_TM' 
            left_primer_tm = round(dictionary.get(Tml), 2)
            right_primer_tm = round(dictionary.get(Tmr), 2)
            Tm_pair = [left_primer_tm, right_primer_tm]

            # GC percent
            GCl = 'PRIMER_LEFT_' + str(cont) + '_GC_PERCENT'
            GCr = 'PRIMER_RIGHT_' + str(cont) + '_GC_PERCENT'
            left_primer_GC = dictionary.get(GCl)
            right_primer_GC = dictionary.get(GCr)
            GC_pair = [left_primer_GC, right_primer_GC]

            # data compile
            data = [primer_pair, Tm_pair, GC_pair] 
            data_set.append(data)
            cont += 1

        return data_set

class IN_SILICO_PCR:

    def in_silico_pcr(Primers_Tm_GC, fasta_seq):

        product_list = []

        for data in Primers_Tm_GC:

            left = data[0][0]
            right = data[0][1]
            start = fasta_seq.find(left)
            reverse_right = ''.join(reversed(right))
            seq_right = Seq(reverse_right)
            complement_right = str(seq_right.complement())
            end = fasta_seq.find(complement_right) + len(right)
            distance = end - start
            product_leght = str(distance) + ' bp'
            product = fasta_seq[start:end]
            seq_product = Seq(product)
            complement_product = str(seq_product.complement())
            
            cont = 0
            lines = []
            
            while cont < distance:
                line = '|'
                lines.append(line)
                cont += 1
            
            product_pair = [data, product_leght, product, complement_product, lines]
            product_list.append(product_pair)

        return product_list

################################oligoanalyzer##################################

class PRIMER_ANALYSIS:

    def oligoanalyzer(pcr_products):
        
        driver = webdriver.Firefox()
        driver.get('https://www.idtdna.com/site/account/login?returnurl=%2Fcalc%2Fanalyzer%2F')
        print('cargo pagina')
        wait = WebDriverWait(driver, 30)
        time.sleep(20)
        #login
        '''
        username = driver.find_element_by_id('UserName')
        password = driver.find_element_by_id('Password')
        button = driver.find_element_by_id('login-button')
        '''
        cookie = driver.find_elements_by_css_selector("a.cc-btn.cc-dismiss")
        lencookie = len(cookie)
        print(lencookie)
                
        if lencookie > 0:
            cookiebtn = driver.find_element_by_css_selector("a.cc-btn.cc-dismiss")
            time.sleep(10)
            cookiebtn.click()

        time.sleep(20)
        
        username = wait.until(ec.presence_of_element_located((By.ID, 'UserName')))
        username.send_keys("AndersonElit")
        password = wait.until(ec.presence_of_element_located((By.ID, 'Password')))
        password.send_keys("Anderlit89")
        button = wait.until(ec.presence_of_element_located((By.ID, 'login-button')))
        button.click()
        time.sleep(20)
        print('inicio de sesion')

        '''        
        cookie = driver.find_elements_by_css_selector("a.cc-btn.cc-dismiss")
        lencookie = len(cookie)
        print(lencookie)
                
        if lencookie > 0:
            cookiebtn = driver.find_element_by_css_selector("a.cc-btn.cc-dismiss")
            time.sleep(10)
            cookiebtn.click()
        '''
            
        driver.refresh()
        time.sleep(10)
        all_data = []
        
        for data_set in pcr_products:
            primers_data = []
            primers = data_set[0][0]
            for primer in primers:

                #textboxprimer = wait.until(ec.presence_of_element_located((By.ID, 'textarea-sequence')))
                driver.refresh()
                time.sleep(30)
                textboxprimerp = driver.find_elements_by_id('textarea-sequence')               
                found = len(textboxprimerp)
                count = 1

                while found == 0:

                    driver.refresh()
                    time.sleep(10)
                    textboxprimer2 = driver.find_elements_by_id('textarea-sequence')
                    found2 = len(textboxprimer2)
                    if found2 > 0:
                        found += 1
                    print('contador hairpin homodimer: ' + str(count))
                    count += 1

                print('textarea-sequence hairpin homodimers encontrado!!!!!!!!')               
                textboxprimer = driver.find_element_by_id('textarea-sequence')
                textboxprimer.send_keys(primer)
                # calc and extract hairpin results
                #hairpin_btn = driver.find_element_by_xpath('//button[text()="Hairpin"]')
                hairpin_btn = wait.until(ec.presence_of_element_located((By.XPATH, '//button[text()="Hairpin"]')))
                hairpin_btn.click()
                time.sleep(30)
                pagehairpin = driver.page_source
                soup = BeautifulSoup(pagehairpin, 'html5lib')
                tables_hairpin = soup.find_all('table', class_ = 'table')
                img_hairpin = tables_hairpin[3].find_all('img', class_ = 'imageThumb')
                imgs_hairpins_src = [img_hairpin[0]['src'], img_hairpin[3]['src']]
                trs = tables_hairpin[3].find_all('tr')
                tds1 = trs[1].find_all('td')
                tds2 = trs[2].find_all('td')
                dGs_hairpin = [tds1[2].text, tds2[2].text]
                hairpins_data = [imgs_hairpins_src, dGs_hairpin]
                                    
                # calc and extract homodimer results
                #homo_btn = driver.find_element_by_xpath('//button[text()="Self-Dimer"]')
                homo_btn = wait.until(ec.presence_of_element_located((By.XPATH, '//button[text()="Self-Dimer"]')))
                homo_btn.click()
                time.sleep(30)
                pagehomo =  driver.page_source
                souphomo = BeautifulSoup(pagehomo, 'html5lib')
                results_homo = souphomo.find_all('div', class_ = 'well')
                spans1 = results_homo[6].find_all('span')
                spans2 = results_homo[7].find_all('span')
                spans3 = results_homo[8].find_all('span')
                dGs_homo = [spans1[0].text, spans2[0].text, spans3[0].text]
                homodimers = [spans1[2].text, spans2[2].text, spans3[2].text]
                homodimers_data = [dGs_homo, homodimers]

                # primer data hairpin homodimer
                primer_data1 = [hairpins_data, homodimers_data]
                primers_data.append(primer_data1)
                textboxprimer.clear()
                
            # calc and extract heterodimer results
            driver.refresh()
            time.sleep(30)
            textboxprimerf = driver.find_element_by_id('textarea-sequence')
            textboxprimerf.send_keys(primers[0])
            #hetero_btn = driver.find_element_by_xpath('//button[text()="Hetero-Dimer"]')
            hetero_btn = wait.until(ec.presence_of_element_located((By.XPATH, '//button[text()="Hetero-Dimer"]')))
            hetero_btn.click()
            time.sleep(20)
            listextarea = driver.find_elements_by_tag_name('textarea')
            length = len(listextarea)
            print('componentes: ' + str(listextarea))
            print('longitud textarea: ' + str(length))
            
            count = 1

            while length == 2:

                driver.refresh()
                time.sleep(20)
                hetero_btn2 = wait.until(ec.presence_of_element_located((By.XPATH, '//button[text()="Hetero-Dimer"]')))
                hetero_btn2.click()
                time.sleep(20)
                listextarea2 = driver.find_elements_by_tag_name('textarea')
                length2 = len(listextarea)
                if length2 > 2:
                    length = length2
                print('contador: ' + str(count))
                count += 1

            print('tagname encontrado!!!!!')
            
            primerr = driver.find_elements_by_tag_name('textarea')[2]
            #print(primerr)
            primerr.send_keys(primers[1])        
            calc_hetero_btn = wait.until(ec.presence_of_element_located((By.XPATH, '//button[text()="Calculate"]')))
            #calc_hetero_btn = driver.find_element_by_xpath('//button[text()="Calculate"]')
            #driver.execute_script("arguments[0].click();", calc_hetero_btn)
            calc_hetero_btn.click()
            time.sleep(30)
            pagehetero = driver.page_source
            souphetero = BeautifulSoup(pagehetero, 'html5lib')
            results_hetero = souphetero.find_all('div', class_ = 'well')
            spanshetero1 = results_hetero[6].find_all('span')
            spanshetero2 = results_hetero[7].find_all('span')
            spanshetero3 = results_hetero[8].find_all('span')
            dGs_hetero = [spanshetero1[0].text, spanshetero2[0].text, spanshetero3[0].text]
            heterodimers = [spanshetero1[2].text, spanshetero2[2].text, spanshetero3[2].text]
            heterodimers_data = [dGs_hetero, heterodimers]

            # primer thermo data
            setdata = [primers_data, heterodimers_data]
            primers_data.append(setdata)
            all_data.append(primers_data)
            textboxprimerf.clear()
            print('\nresultados: ' + str(all_data))
                        
        time.sleep(30)
        driver.close()
        #return all_data
        print('analisis completado')
        return all_data
        

class COTIZER:

    def IDT(pcr_products):

        return pcr_products

        

    

          
    

    

    
