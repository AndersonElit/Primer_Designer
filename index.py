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
    
    all_results = [[[[['GCCTTGCGAATGAATGTGCT', 'CACTGCTAGTACCACCAGGC'], [59.83, 60.11], [50.0, 60.0]], '87 bp', 'GCCTTGCGAATGAATGTGCTCAAGTTTTGAGTGAAATAGTTATGTGTGGCGGTTGCTATTATGTTAAGCCTGGTGGTACTAGCAGTG\n|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nCGGAACGCTTACTTACACGAGTTCAAAACTCACTTTATCAATACACACCGCCAACGATAATACAATTCGGACCACCATGATCGTCAC'], [[[['https://www.idtdna.com//unafold/home/extras?img=2ddfa896-84a1-490d-9771-b9797c92a94b/2ddfa896-84a1-490d-9771-b9797c92a94b_1_thm.jpg&cdnblocker=E637332229659134407', 'https://www.idtdna.com//unafold/home/extras?img=2ddfa896-84a1-490d-9771-b9797c92a94b/2ddfa896-84a1-490d-9771-b9797c92a94b_2_thm.jpg&cdnblocker=E637332229659134407'], ['-0.29', '0.02']], [['-3.61 kcal/mole', '-3.14 kcal/mole', '-3.14 kcal/mole'], ["5'       GCCTTGCGAATGAATGTGCT\n           : : || : :        \n3' TCGTGTAAGTAAGCGTTCCG", "5'                   GCCTTGCGAATGAATGTGCT\n                     ||                  \n3' TCGTGTAAGTAAGCGTTCCG", "5'              GCCTTGCGAATGAATGTGCT\n                ||   ::             \n3' TCGTGTAAGTAAGCGTTCCG"]]], [[['https://www.idtdna.com//unafold/home/extras?img=d46a3ccb-fd4e-4ee3-b25e-72e04403cad9/d46a3ccb-fd4e-4ee3-b25e-72e04403cad9_1_thm.jpg&cdnblocker=E637332230633702140', 'https://www.idtdna.com//unafold/home/extras?img=d46a3ccb-fd4e-4ee3-b25e-72e04403cad9/d46a3ccb-fd4e-4ee3-b25e-72e04403cad9_2_thm.jpg&cdnblocker=E637332230633702140'], ['-0.55', '-0.52']], [['-4.16 kcal/mole', '-3.65 kcal/mole', '-3.14 kcal/mole'], ["5'       CACTGCTAGTACCACCAGGC\n            : |||| :         \n3' CGGACCACCATGATCGTCAC", "5' CACTGCTAGTACCACCAGGC\n     ::: : |||| : :::  \n3' CGGACCACCATGATCGTCAC", "5'           CACTGCTAGTACCACCAGGC\n                 ||              \n3' CGGACCACCATGATCGTCAC"]]], [[...], [['-7.81 kcal/mole', '-3.3 kcal/mole', '-3.3 kcal/mole'], ["5' GCCTTGCGAATGAATGTGCT\n   |||| : :     :   :  \n3' CGGACCACCATGATCGTCAC", "5'   GCCTTGCGAATGAATGTGCT\n         ::  : :    |||  \n3' CGGACCACCATGATCGTCAC", "5' GCCTTGCGAATGAATGTGCT\n                  ||| :          \n3'           CGGACCACCATGATCGTCAC"]]]], '$14.80 USD'], [[[['GCCTTGCGAATGAATGTGCT', 'ACTGCTAGTACCACCAGGCT'], [59.83, 60.25], [50.0, 55.0]], '86 bp', 'GCCTTGCGAATGAATGTGCTCAAGTTTTGAGTGAAATAGTTATGTGTGGCGGTTGCTATTATGTTAAGCCTGGTGGTACTAGCAGT\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nCGGAACGCTTACTTACACGAGTTCAAAACTCACTTTATCAATACACACCGCCAACGATAATACAATTCGGACCACCATGATCGTCA'], [[[['https://www.idtdna.com//unafold/home/extras?img=7713d3c2-1473-40c1-ae75-4406b5db5237/7713d3c2-1473-40c1-ae75-4406b5db5237_1_thm.jpg&cdnblocker=E637332232483602955', 'https://www.idtdna.com//unafold/home/extras?img=7713d3c2-1473-40c1-ae75-4406b5db5237/7713d3c2-1473-40c1-ae75-4406b5db5237_2_thm.jpg&cdnblocker=E637332232483602955'], ['-0.29', '0.02']], [['-3.61 kcal/mole', '-3.14 kcal/mole', '-3.14 kcal/mole'], ["5'       GCCTTGCGAATGAATGTGCT\n           : : || : :        \n3' TCGTGTAAGTAAGCGTTCCG", "5'                   GCCTTGCGAATGAATGTGCT\n                     ||                  \n3' TCGTGTAAGTAAGCGTTCCG", "5'              GCCTTGCGAATGAATGTGCT\n                ||   ::             \n3' TCGTGTAAGTAAGCGTTCCG"]]], [[['https://www.idtdna.com//unafold/home/extras?img=fe97cdd0-93c5-4aad-89ce-e183042e66ab/fe97cdd0-93c5-4aad-89ce-e183042e66ab_1_thm.jpg&cdnblocker=E637332233477217199', 'https://www.idtdna.com//unafold/home/extras?img=fe97cdd0-93c5-4aad-89ce-e183042e66ab/fe97cdd0-93c5-4aad-89ce-e183042e66ab_2_thm.jpg&cdnblocker=E637332233477217199'], ['-0.55', '-0.52']], [['-4.16 kcal/mole', '-3.65 kcal/mole', '-3.14 kcal/mole'], ["5'         ACTGCTAGTACCACCAGGCT\n             : |||| :          \n3' TCGGACCACCATGATCGTCA", "5'   ACTGCTAGTACCACCAGGCT\n      ::: : |||| : :::   \n3' TCGGACCACCATGATCGTCA", "5'             ACTGCTAGTACCACCAGGCT\n                  ||               \n3' TCGGACCACCATGATCGTCA"]]], [[...], [['-7.81 kcal/mole', '-3.3 kcal/mole', '-3.14 kcal/mole'], ["5'  GCCTTGCGAATGAATGTGCT\n    |||| : :     :   :  \n3' TCGGACCACCATGATCGTCA", "5' GCCTTGCGAATGAATGTGCT\n            :     ||| :         \n3'          TCGGACCACCATGATCGTCA", "5'                GCCTTGCGAATGAATGTGCT\n                  ||  :               \n3' TCGGACCACCATGATCGTCA"]]]], '$14.80 USD'], [[[['AGCAGTCTGTATGGCTTGCC', 'TGCAGGCGCATAACAAATCG'], [60.39, 59.9], [55.0, 50.0]], '88 bp', 'AGCAGTCTGTATGGCTTGCCAATTGTGACTTTGATATTGTAGTGGCTTGGCATGTAGTTCGTGATTCACGATTTGTTATGCGCCTGCA\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nTCGTCAGACATACCGAACGGTTAACACTGAAACTATAACATCACCGAACCGTACATCAAGCACTAAGTGCTAAACAATACGCGGACGT'], [[[['https://www.idtdna.com//unafold/home/extras?img=895103e8-e441-4811-a81d-d2173f58d73b/895103e8-e441-4811-a81d-d2173f58d73b_1_thm.jpg&cdnblocker=E637332235397897316', 'https://www.idtdna.com//unafold/home/extras?img=895103e8-e441-4811-a81d-d2173f58d73b/895103e8-e441-4811-a81d-d2173f58d73b_2_thm.jpg&cdnblocker=E637332235397897316'], ['-1', '-0.64']], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-4.74 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n               |||  :::            \n3'             CCGTTCGGTATGTCTGACGA", "5' AGCAGTCTGTATGGCTTGCC\n    |||  :  ::  :  ::: \n3' CCGTTCGGTATGTCTGACGA", "5'     AGCAGTCTGTATGGCTTGCC\n       |||  :    :  :::    \n3' CCGTTCGGTATGTCTGACGA"]]], [[['https://www.idtdna.com//unafold/home/extras?img=d3fd96c1-3d78-49f4-81a5-11d22a2d0e18/d3fd96c1-3d78-49f4-81a5-11d22a2d0e18_1_thm.jpg&cdnblocker=E637332236391525377', 'https://www.idtdna.com//unafold/home/extras?img=d3fd96c1-3d78-49f4-81a5-11d22a2d0e18/d3fd96c1-3d78-49f4-81a5-11d22a2d0e18_2_thm.jpg&cdnblocker=E637332236391525377'], ['-1.15', '-0.95']], [['-9.89 kcal/mole', '-7.05 kcal/mole', '-5.09 kcal/mole'], ["5'       TGCAGGCGCATAACAAATCG\n            : |||| :         \n3' GCTAAACAATACGCGGACGT", "5'                 TGCAGGCGCATAACAAATCG\n                   ||||                \n3' GCTAAACAATACGCGGACGT", "5'           TGCAGGCGCATAACAAATCG\n             |||    :::          \n3' GCTAAACAATACGCGGACGT"]]], [[...], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-5.09 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n          : : :     |||    \n3'     GCTAAACAATACGCGGACGT", "5'                 AGCAGTCTGTATGGCTTGCC\n                    |||                \n3' GCTAAACAATACGCGGACGT", "5' AGCAGTCTGTATGGCTTGCC\n    :   : :     :: ||| \n3' GCTAAACAATACGCGGACGT"]]]], '$14.80 USD'], [[[['AGCAGTCTGTATGGCTTGCC', 'GCAGGCGCATAACAAATCGT'], [60.39, 59.9], [55.0, 50.0]], '87 bp', 'AGCAGTCTGTATGGCTTGCCAATTGTGACTTTGATATTGTAGTGGCTTGGCATGTAGTTCGTGATTCACGATTTGTTATGCGCCTGC\n|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nTCGTCAGACATACCGAACGGTTAACACTGAAACTATAACATCACCGAACCGTACATCAAGCACTAAGTGCTAAACAATACGCGGACG'], [[[['https://www.idtdna.com//unafold/home/extras?img=c2405a4d-5a0c-4a6d-b817-2e1322587cd1/c2405a4d-5a0c-4a6d-b817-2e1322587cd1_1_thm.jpg&cdnblocker=E637332238263299637', 'https://www.idtdna.com//unafold/home/extras?img=c2405a4d-5a0c-4a6d-b817-2e1322587cd1/c2405a4d-5a0c-4a6d-b817-2e1322587cd1_2_thm.jpg&cdnblocker=E637332238263299637'], ['-1', '-0.64']], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-4.74 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n               |||  :::            \n3'             CCGTTCGGTATGTCTGACGA", "5' AGCAGTCTGTATGGCTTGCC\n    |||  :  ::  :  ::: \n3' CCGTTCGGTATGTCTGACGA", "5'     AGCAGTCTGTATGGCTTGCC\n       |||  :    :  :::    \n3' CCGTTCGGTATGTCTGACGA"]]], [[['https://www.idtdna.com//unafold/home/extras?img=13bfb323-0f4f-4fa9-8340-dbbd4c2ddc98/13bfb323-0f4f-4fa9-8340-dbbd4c2ddc98_1_thm.jpg&cdnblocker=E637332239300206859', 'https://www.idtdna.com//unafold/home/extras?img=13bfb323-0f4f-4fa9-8340-dbbd4c2ddc98/13bfb323-0f4f-4fa9-8340-dbbd4c2ddc98_2_thm.jpg&cdnblocker=E637332239300206859'], ['-0.54', '0.1']], [['-9.89 kcal/mole', '-3.61 kcal/mole', '-3.61 kcal/mole'], ["5'         GCAGGCGCATAACAAATCGT\n             : |||| :          \n3' TGCTAAACAATACGCGGACG", "5' GCAGGCGCATAACAAATCGT\n        ||  :    :  ::     \n3'     TGCTAAACAATACGCGGACG", "5' GCAGGCGCATAACAAATCGT\n                    ||                 \n3'                 TGCTAAACAATACGCGGACG"]]], [[...], [['-6.21 kcal/mole', '-5.09 kcal/mole', '-5.09 kcal/mole'], ["5' AGCAGTCTGTATGGCTTGCC\n      :   : : :     |||   \n3'    TGCTAAACAATACGCGGACG", "5'  AGCAGTCTGTATGGCTTGCC\n     :   : :     :: ||| \n3' TGCTAAACAATACGCGGACG", "5' AGCAGTCTGTATGGCTTGCC\n         :  : ::   |||      \n3'      TGCTAAACAATACGCGGACG"]]]], '$14.80 USD'], [[[['TCTTCAGCCGAACAGACAGC', 'CCGCGCAATTACCACCTAGA'], [60.32, 60.18], [55.0, 55.0]], '75 bp', 'TCTTCAGCCGAACAGACAGCCTATGCTCAAGCCATTACTGGAGCCCACAAGGTAATCTAGGTGGTAATTGCGCGG\n|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nAGAAGTCGGCTTGTCTGTCGGATACGAGTTCGGTAATGACCTCGGGTGTTCCATTAGATCCACCATTAACGCGCC'], [[[['https://www.idtdna.com//unafold/home/extras?img=6d75dec3-f20c-4e18-a8eb-824d63aabcc0/6d75dec3-f20c-4e18-a8eb-824d63aabcc0_1_thm.jpg&cdnblocker=E637332241299171371', 'https://www.idtdna.com//unafold/home/extras?img=6d75dec3-f20c-4e18-a8eb-824d63aabcc0/6d75dec3-f20c-4e18-a8eb-824d63aabcc0_2_thm.jpg&cdnblocker=E637332241299171371'], ['0.46', '0.53']], [['-3.61 kcal/mole', '-3.52 kcal/mole', '-3.17 kcal/mole'], ["5'   TCTTCAGCCGAACAGACAGC\n     : :     ||     : :  \n3' CGACAGACAAGCCGACTTCT", "5'       TCTTCAGCCGAACAGACAGC\n         : ||| :: ::: :      \n3' CGACAGACAAGCCGACTTCT", "5'     TCTTCAGCCGAACAGACAGC\n       |||          :::    \n3' CGACAGACAAGCCGACTTCT"]]], [[['https://www.idtdna.com//unafold/home/extras?img=12ee9a9d-c94f-4c1f-850c-dd9cb6623fba/12ee9a9d-c94f-4c1f-850c-dd9cb6623fba_1_thm.jpg&cdnblocker=E637332242424675765', 'https://www.idtdna.com//unafold/home/extras?img=12ee9a9d-c94f-4c1f-850c-dd9cb6623fba/12ee9a9d-c94f-4c1f-850c-dd9cb6623fba_2_thm.jpg&cdnblocker=E637332242424675765'], ['1.98', '2.58']], [['-10.36 kcal/mole', '-9.89 kcal/mole', '-5.36 kcal/mole'], ["5'               CCGCGCAATTACCACCTAGA\n                  ||||               \n3' AGATCCACCATTAACGCGCC", "5'             CCGCGCAATTACCACCTAGA\n                 ||||              \n3' AGATCCACCATTAACGCGCC", "5'     CCGCGCAATTACCACCTAGA\n           : |||| :        \n3' AGATCCACCATTAACGCGCC"]]], [[...], [['-3.61 kcal/mole', '-3.61 kcal/mole', '-3.17 kcal/mole'], ["5'          TCTTCAGCCGAACAGACAGC\n            :  :    ||          \n3' AGATCCACCATTAACGCGCC", "5'        TCTTCAGCCGAACAGACAGC\n            :     ||          \n3' AGATCCACCATTAACGCGCC", "5' TCTTCAGCCGAACAGACAGC\n   |||       ::  :   : \n3' AGATCCACCATTAACGCGCC"]]]], '$14.80 USD']]

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


