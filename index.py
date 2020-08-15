from flask import Flask, render_template, request, redirect, url_for, flash
from Primer_Designer_Functions import SEARCH, PRIMER_DESIGN, IN_SILICO_PCR, PRIMER_ANALYSIS, COTIZER
from werkzeug.utils import secure_filename
import os

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('fastaseq.html')


@app.route('/keyword')
def keyword():
    return render_template('keyword.html')

@app.route('/geneid')
def geneid():
    return render_template('geneid.html')

@app.route('/fastaseq')
def fastaseq():
    return render_template('fastaseq.html')

# def function
def analysis(fasta_seq):
    # Design primers
    dictionary = PRIMER_DESIGN.PRIMER3(fasta_seq)
    Primers_Tm_GC = PRIMER_DESIGN.primers_tm_gc(dictionary)
    pcr_products = IN_SILICO_PCR.in_silico_pcr(Primers_Tm_GC, fasta_seq)
    thermo_analysis = PRIMER_ANALYSIS.oligoanalyzer(pcr_products)
    IDT_prices = COTIZER.IDT(pcr_products)
    all_data = []
    count = 0

    while count <= 4:
        
        sets = [pcr_products[count], thermo_analysis[count], IDT_prices[count]]
        all_data.append(sets)
        count += 1

    return all_data

@app.route('/results_fasta', methods=['POST'])
def results_fasta():
    
    all_results = [[[[['GCCTTGCGAATGAATGTGCT', 'CACTGCTAGTACCACCAGGC'], [59.83, 60.11], [50.0, 60.0]], '87 bp', 'GCCTTGCGAATGAATGTGCTCAAGTTTTGAGTGAAATAGTTATGTGTGGCGGTTGCTATTATGTTAAGCCTGGTGGTACTAGCAGTG', 'CGGAACGCTTACTTACACGAGTTCAAAACTCACTTTATCAATACACACCGCCAACGATAATACAATTCGGACCACCATGATCGTCAC', ['|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|']], [[[['https://www.idtdna.com//unafold/home/extras?img=1cda4a42-5062-4b9a-b5d8-f4aab9a85205/1cda4a42-5062-4b9a-b5d8-f4aab9a85205_1_thm.jpg&cdnblocker=E637330686923416362', 'https://www.idtdna.com//unafold/home/extras?img=1cda4a42-5062-4b9a-b5d8-f4aab9a85205/1cda4a42-5062-4b9a-b5d8-f4aab9a85205_2_thm.jpg&cdnblocker=E637330686923416362'], ['-0.29', '0.02']], [['-3.61 kcal/mole', '-3.14 kcal/mole', '-3.14 kcal/mole'], ["5'       GCCTTGCGAATGAATGTGCT\n           : : || : :        \n3' TCGTGTAAGTAAGCGTTCCG", "5'                   GCCTTGCGAATGAATGTGCT\n                     ||                  \n3' TCGTGTAAGTAAGCGTTCCG", "5'              GCCTTGCGAATGAATGTGCT\n                ||   ::             \n3' TCGTGTAAGTAAGCGTTCCG"]]], [[['https://www.idtdna.com//unafold/home/extras?img=ca81575c-e17b-4324-8531-398ab57cb492/ca81575c-e17b-4324-8531-398ab57cb492_1_thm.jpg&cdnblocker=E637330687898791323', 'https://www.idtdna.com//unafold/home/extras?img=ca81575c-e17b-4324-8531-398ab57cb492/ca81575c-e17b-4324-8531-398ab57cb492_2_thm.jpg&cdnblocker=E637330687898791323'], ['-0.55', '-0.52']], [['-4.16 kcal/mole', '-3.65 kcal/mole', '-3.14 kcal/mole'], ["5'       CACTGCTAGTACCACCAGGC\n            : |||| :         \n3' CGGACCACCATGATCGTCAC", "5' CACTGCTAGTACCACCAGGC\n     ::: : |||| : :::  \n3' CGGACCACCATGATCGTCAC", "5'           CACTGCTAGTACCACCAGGC\n                 ||              \n3' CGGACCACCATGATCGTCAC"]]], [[...], [['-7.81 kcal/mole', '-3.3 kcal/mole', '-3.3 kcal/mole'], ["5' GCCTTGCGAATGAATGTGCT\n   |||| : :     :   :  \n3' CGGACCACCATGATCGTCAC", "5'   GCCTTGCGAATGAATGTGCT\n         ::  : :    |||  \n3' CGGACCACCATGATCGTCAC", "5' GCCTTGCGAATGAATGTGCT\n                  ||| :          \n3'           CGGACCACCATGATCGTCAC"]]]], '$14.80 USD'], [[[['GCCTTGCGAATGAATGTGCT', 'ACTGCTAGTACCACCAGGCT'], [59.83, 60.25], [50.0, 55.0]], '86 bp', 'GCCTTGCGAATGAATGTGCTCAAGTTTTGAGTGAAATAGTTATGTGTGGCGGTTGCTATTATGTTAAGCCTGGTGGTACTAGCAGT', 'CGGAACGCTTACTTACACGAGTTCAAAACTCACTTTATCAATACACACCGCCAACGATAATACAATTCGGACCACCATGATCGTCA', ['|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|']], [[[['https://www.idtdna.com//unafold/home/extras?img=5940321a-eb26-4cfe-bbb4-60d90f539b8c/5940321a-eb26-4cfe-bbb4-60d90f539b8c_1_thm.jpg&cdnblocker=E637330689721866202', 'https://www.idtdna.com//unafold/home/extras?img=5940321a-eb26-4cfe-bbb4-60d90f539b8c/5940321a-eb26-4cfe-bbb4-60d90f539b8c_2_thm.jpg&cdnblocker=E637330689721866202'], ['-0.29', '0.02']], [['-3.61 kcal/mole', '-3.14 kcal/mole', '-3.14 kcal/mole'], ["5'       GCCTTGCGAATGAATGTGCT\n           : : || : :        \n3' TCGTGTAAGTAAGCGTTCCG", "5'                   GCCTTGCGAATGAATGTGCT\n                     ||                  \n3' TCGTGTAAGTAAGCGTTCCG", "5'              GCCTTGCGAATGAATGTGCT\n                ||   ::             \n3' TCGTGTAAGTAAGCGTTCCG"]]], [[['https://www.idtdna.com//unafold/home/extras?img=6a71f5d7-b304-4671-aa16-faa03c8ab08b/6a71f5d7-b304-4671-aa16-faa03c8ab08b_1_thm.jpg&cdnblocker=E637330690685822463', 'https://www.idtdna.com//unafold/home/extras?img=6a71f5d7-b304-4671-aa16-faa03c8ab08b/6a71f5d7-b304-4671-aa16-faa03c8ab08b_2_thm.jpg&cdnblocker=E637330690685822463'], ['-0.55', '-0.52']], [['-4.16 kcal/mole', '-3.65 kcal/mole', '-3.14 kcal/mole'], ["5'         ACTGCTAGTACCACCAGGCT\n             : |||| :          \n3' TCGGACCACCATGATCGTCA", "5'   ACTGCTAGTACCACCAGGCT\n      ::: : |||| : :::   \n3' TCGGACCACCATGATCGTCA", "5'             ACTGCTAGTACCACCAGGCT\n                  ||               \n3' TCGGACCACCATGATCGTCA"]]], [[...], [['-7.81 kcal/mole', '-3.3 kcal/mole', '-3.14 kcal/mole'], ["5'  GCCTTGCGAATGAATGTGCT\n    |||| : :     :   :  \n3' TCGGACCACCATGATCGTCA", "5' GCCTTGCGAATGAATGTGCT\n            :     ||| :         \n3'          TCGGACCACCATGATCGTCA", "5'                GCCTTGCGAATGAATGTGCT\n                  ||  :               \n3' TCGGACCACCATGATCGTCA"]]]], '$14.80 USD'], [[[['AGCAGTCTGTATGGCTTGCC', 'TGCAGGCGCATAACAAATCG'], [60.39, 59.9], [55.0, 50.0]], '88 bp', 'AGCAGTCTGTATGGCTTGCCAATTGTGACTTTGATATTGTAGTGGCTTGGCATGTAGTTCGTGATTCACGATTTGTTATGCGCCTGCA', 'TCGTCAGACATACCGAACGGTTAACACTGAAACTATAACATCACCGAACCGTACATCAAGCACTAAGTGCTAAACAATACGCGGACGT', ['|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|']], [[[['https://www.idtdna.com//unafold/home/extras?img=3746e073-bfe8-46e3-b27a-1e1ae20ceee2/3746e073-bfe8-46e3-b27a-1e1ae20ceee2_1_thm.jpg&cdnblocker=E637330692511568487', 'https://www.idtdna.com//unafold/home/extras?img=3746e073-bfe8-46e3-b27a-1e1ae20ceee2/3746e073-bfe8-46e3-b27a-1e1ae20ceee2_2_thm.jpg&cdnblocker=E637330692511568487'], ['-1', '-0.64']], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-4.74 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n               |||  :::            \n3'             CCGTTCGGTATGTCTGACGA", "5' AGCAGTCTGTATGGCTTGCC\n    |||  :  ::  :  ::: \n3' CCGTTCGGTATGTCTGACGA", "5'     AGCAGTCTGTATGGCTTGCC\n       |||  :    :  :::    \n3' CCGTTCGGTATGTCTGACGA"]]], [[['https://www.idtdna.com//unafold/home/extras?img=77448a9f-ad46-4bee-8871-6851f376841e/77448a9f-ad46-4bee-8871-6851f376841e_1_thm.jpg&cdnblocker=E637330693477242512', 'https://www.idtdna.com//unafold/home/extras?img=77448a9f-ad46-4bee-8871-6851f376841e/77448a9f-ad46-4bee-8871-6851f376841e_2_thm.jpg&cdnblocker=E637330693477242512'], ['-1.15', '-0.95']], [['-9.89 kcal/mole', '-7.05 kcal/mole', '-5.09 kcal/mole'], ["5'       TGCAGGCGCATAACAAATCG\n            : |||| :         \n3' GCTAAACAATACGCGGACGT", "5'                 TGCAGGCGCATAACAAATCG\n                   ||||                \n3' GCTAAACAATACGCGGACGT", "5'           TGCAGGCGCATAACAAATCG\n             |||    :::          \n3' GCTAAACAATACGCGGACGT"]]], [[...], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-5.09 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n          : : :     |||    \n3'     GCTAAACAATACGCGGACGT", "5'                 AGCAGTCTGTATGGCTTGCC\n                    |||                \n3' GCTAAACAATACGCGGACGT", "5' AGCAGTCTGTATGGCTTGCC\n    :   : :     :: ||| \n3' GCTAAACAATACGCGGACGT"]]]], '$14.80 USD'], [[[['AGCAGTCTGTATGGCTTGCC', 'GCAGGCGCATAACAAATCGT'], [60.39, 59.9], [55.0, 50.0]], '87 bp', 'AGCAGTCTGTATGGCTTGCCAATTGTGACTTTGATATTGTAGTGGCTTGGCATGTAGTTCGTGATTCACGATTTGTTATGCGCCTGC', 'TCGTCAGACATACCGAACGGTTAACACTGAAACTATAACATCACCGAACCGTACATCAAGCACTAAGTGCTAAACAATACGCGGACG', ['|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|']], [[[['https://www.idtdna.com//unafold/home/extras?img=797fc9b1-3f61-4fb4-b94a-46a2539553d7/797fc9b1-3f61-4fb4-b94a-46a2539553d7_1_thm.jpg&cdnblocker=E637330695295789916', 'https://www.idtdna.com//unafold/home/extras?img=797fc9b1-3f61-4fb4-b94a-46a2539553d7/797fc9b1-3f61-4fb4-b94a-46a2539553d7_2_thm.jpg&cdnblocker=E637330695295789916'], ['-1', '-0.64']], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-4.74 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n               |||  :::            \n3'             CCGTTCGGTATGTCTGACGA", "5' AGCAGTCTGTATGGCTTGCC\n    |||  :  ::  :  ::: \n3' CCGTTCGGTATGTCTGACGA", "5'     AGCAGTCTGTATGGCTTGCC\n       |||  :    :  :::    \n3' CCGTTCGGTATGTCTGACGA"]]], [[['https://www.idtdna.com//unafold/home/extras?img=3a45ab66-c412-4773-821b-a57dffb08a63/3a45ab66-c412-4773-821b-a57dffb08a63_1_thm.jpg&cdnblocker=E637330696273352576', 'https://www.idtdna.com//unafold/home/extras?img=3a45ab66-c412-4773-821b-a57dffb08a63/3a45ab66-c412-4773-821b-a57dffb08a63_2_thm.jpg&cdnblocker=E637330696273352576'], ['-0.54', '0.1']], [['-9.89 kcal/mole', '-3.61 kcal/mole', '-3.61 kcal/mole'], ["5'         GCAGGCGCATAACAAATCGT\n             : |||| :          \n3' TGCTAAACAATACGCGGACG", "5' GCAGGCGCATAACAAATCGT\n        ||  :    :  ::     \n3'     TGCTAAACAATACGCGGACG", "5' GCAGGCGCATAACAAATCGT\n                    ||                 \n3'                 TGCTAAACAATACGCGGACG"]]], [[...], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-5.09 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n      :   : : :     |||   \n3'    TGCTAAACAATACGCGGACG", "5'  AGCAGTCTGTATGGCTTGCC\n     :   : :     :: ||| \n3' TGCTAAACAATACGCGGACG", "5' AGCAGTCTGTATGGCTTGCC\n         :  : ::   |||      \n3'      TGCTAAACAATACGCGGACG"]]]], '$14.80 USD'], [[[['TCTTCAGCCGAACAGACAGC', 'CCGCGCAATTACCACCTAGA'], [60.32, 60.18], [55.0, 55.0]], '75 bp', 'TCTTCAGCCGAACAGACAGCCTATGCTCAAGCCATTACTGGAGCCCACAAGGTAATCTAGGTGGTAATTGCGCGG', 'AGAAGTCGGCTTGTCTGTCGGATACGAGTTCGGTAATGACCTCGGGTGTTCCATTAGATCCACCATTAACGCGCC', ['|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|']], [[[['https://www.idtdna.com//unafold/home/extras?img=1ea9529a-5620-4b42-aa0a-f13e6dcea7b9/1ea9529a-5620-4b42-aa0a-f13e6dcea7b9_1_thm.jpg&cdnblocker=E637330698110802233', 'https://www.idtdna.com//unafold/home/extras?img=1ea9529a-5620-4b42-aa0a-f13e6dcea7b9/1ea9529a-5620-4b42-aa0a-f13e6dcea7b9_2_thm.jpg&cdnblocker=E637330698110802233'], ['0.46', '0.53']], [['-3.61 kcal/mole', '-3.52 kcal/mole', '-3.17 kcal/mole'], ["5'   TCTTCAGCCGAACAGACAGC\n     : :     ||     : :  \n3' CGACAGACAAGCCGACTTCT", "5'       TCTTCAGCCGAACAGACAGC\n         : ||| :: ::: :      \n3' CGACAGACAAGCCGACTTCT", "5'     TCTTCAGCCGAACAGACAGC\n       |||          :::    \n3' CGACAGACAAGCCGACTTCT"]]], [[['https://www.idtdna.com//unafold/home/extras?img=c2512018-bea2-4ff4-9479-768367cac868/c2512018-bea2-4ff4-9479-768367cac868_1_thm.jpg&cdnblocker=E637330699070855097', 'https://www.idtdna.com//unafold/home/extras?img=c2512018-bea2-4ff4-9479-768367cac868/c2512018-bea2-4ff4-9479-768367cac868_2_thm.jpg&cdnblocker=E637330699070855097'], ['1.98', '2.58']], [['-10.36 kcal/mole', '-9.89 kcal/mole', '-5.36 kcal/mole'], ["5'               CCGCGCAATTACCACCTAGA\n                  ||||               \n3' AGATCCACCATTAACGCGCC", "5'             CCGCGCAATTACCACCTAGA\n                 ||||              \n3' AGATCCACCATTAACGCGCC", "5'     CCGCGCAATTACCACCTAGA\n           : |||| :        \n3' AGATCCACCATTAACGCGCC"]]], [[...], [['-3.61 kcal/mole', '-3.61 kcal/mole', '-3.17 kcal/mole'], ["5'          TCTTCAGCCGAACAGACAGC\n            :  :    ||          \n3' AGATCCACCATTAACGCGCC", "5'        TCTTCAGCCGAACAGACAGC\n            :     ||          \n3' AGATCCACCATTAACGCGCC", "5' TCTTCAGCCGAACAGACAGC\n   |||       ::  :   : \n3' AGATCCACCATTAACGCGCC"]]]], '$14.80 USD']]

    return render_template('analysis_fasta.html', all_results=all_results)
    
    '''
    seq = request.form.get('text')
    f = request.files['file']
    list_file = str(f)
    list_file2 = list_file.split(' ')

    def seq_assemble(fasta):
        list = fasta.split('\n')
        seq_title = list[0]
        fragments = list[1:]
        sequence = ''
        for fragment in fragments:
            found = fragment.find('\r')
            if found == -1:
                sequence += fragment
            else:
                new_fragment = fragment.replace('\r','')
                sequence += new_fragment
        return sequence

    if list_file2[1] == "''":
        full_seq = seq_assemble(seq)
        all_results = analysis(full_seq)
        return render_template('analysis_fasta.html', all_results=all_results)
    else:
        f.save(secure_filename(f.filename))
        lenf = len(list_file2[1]) - 1
        filen = list_file2[1][1:lenf]
        fileopen = open(filen, 'r')
        content = str(fileopen.read())
        fileopen.close()
        os.remove(filen)
        full_seq = seq_assemble(content)
        all_results = analysis(full_seq)
        return render_template('analysis_fasta.html', all_results=all_results)
    '''

if __name__ == '__main__':
    app.run(debug=True)


