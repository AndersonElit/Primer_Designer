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
        fasta_region = title_seq[0] + '\n' + gene

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

    def Hairpins(Primers_Tm_GC):

        data_set = []

        for data in Primers_Tm_GC:

            primers = data[0]
            primer_left = primers[0]
            primer_right = primers[1]
            Hairpin_left = primer3.calcHairpin(primer_left)
            Hairpin_right = primer3.calcHairpin(primer_right)
            thermo_left = [round(Hairpin_left.tm, 2), round(Hairpin_left.dg, 2), round(Hairpin_left.dh, 2), round(Hairpin_left.ds, 2)]
            thermo_right = [round(Hairpin_right.tm, 2), round(Hairpin_right.dg, 2), round(Hairpin_right.dh, 2), round(Hairpin_right.ds, 2)]
            hairpin_data = [thermo_left, thermo_right]
            data.append(hairpin_data)
            data_set.append(data)

        return data_set

    def Homodimers(Hairpins_calc):

        data_set = []

        for data in Hairpins_calc:

            primers = data[0]
            primer_left = primers[0]
            primer_right = primers[1]
            Homodimer_left = primer3.calcHomodimer(primer_left)
            Homodimer_right = primer3.calcHomodimer(primer_right)
            thermo_left = [round(Homodimer_left.tm, 2), round(Homodimer_left.dg, 2), round(Homodimer_left.dh, 2), round(Homodimer_left.ds, 2)]
            thermo_right = [round(Homodimer_right.tm, 2), round(Homodimer_right.dg, 2), round(Homodimer_right.dh, 2), round(Homodimer_right.ds, 2)]
            Homodimer_data = [thermo_left, thermo_right]
            data.append(Homodimer_data)
            data_set.append(data)

        return data_set

    def Heterodimers(Homodimers_calc):

        data_set = []

        for data in Homodimers_calc:

            primers = data[0]
            primer_left = primers[0]
            primer_right = primers[1]
            Heterodimer = primer3.calcHeterodimer(primer_left, primer_right)
            thermo = [round(Heterodimer.tm, 2), round(Heterodimer.dg, 2), round(Heterodimer.dh, 2), round(Heterodimer.ds, 2)]
            data.append(thermo)
            data_set.append(data)

        return data_set

class IN_SILICO_PCR:

    # get primers

    def primers(Heterodimers_calc):

        primers_set = []

        for data in Heterodimers_calc:

            primers = data[0]
            primer_left = primers[0]
            primer_right = primers[1]
            primer_pair = [primer_left, primer_right]
            primers_set.append(primer_pair)

        return primers_set

    def in_silico_pcr(primer_list):

        count1 = 0
        product_list = []

        for primer_pair in primer_list:

            forward = primer_pair[0]
            reverse = primer_pair[1]
            ipcress_name = 'primers' + str(count) + '.ipcress'
            ipcress_content = 'experiment' + str(count) + ' ' + forward + ' ' + reverse + ' ' + '75 225'
            file= open(ipcress_name,"w+")
            file.write(ipcress_content)
            file.close()
            experiment = 'ipcress -i ' + ipcress_name + ' -s seq.fasta -P'
            insilicopcr = str(bash(experiment))
            pcr_list = insilicopcr.split('\n')

            # extract product length(bp)
            product_length = pcr_list[6]
            product_lengthc = product_length.replace('    ', '')
            product_lengthc_list = product_lengthc.split(' ')
            product_length_bp = product_lengthc_list[1] + ' ' + product_lengthc_list[2]

            #extract fasta product
            length = len(pcr_list) - 1
            fasta_list = pcr_list[16:length]
            fasta_sequence = fasta_list[0] + '\n'

            count2 = 1
            len_fasta = len(fasta_list) - 1

            while cont <= len_fasta:

                fragment = fasta_list[cont]
                fasta_sequence += fragment
                count2 += 1

            product_list.append(fasta_sequence)
            count1 += 1

        return experiment
            
        
            


        
        

        
            
            

            

        
        

        
    

    

    