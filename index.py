from flask import Flask, render_template, request, redirect, url_for
from Primer_Designer_Functions import SEARCH, PRIMER_DESIGN, IN_SILICO_PCR

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('keyword.html')

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
def analysis(page_content):
    # Design primers
    fasta_seq = PRIMER_DESIGN.get_fasta(page_content)
    dictionary = PRIMER_DESIGN.PRIMER3(fasta_seq)
    Primers_Tm_GC = PRIMER_DESIGN.primers_tm_gc(dictionary)
    Hairpins_calc = PRIMER_DESIGN.Hairpins(Primers_Tm_GC)
    Homodimers_calc = PRIMER_DESIGN.Homodimers(Hairpins_calc)
    Heterodimers_calc = PRIMER_DESIGN.Heterodimers(Homodimers_calc)
    # In_silico PCR
    primer_list = IN_SILICO_PCR.primers(Heterodimers_calc)
    pcr_products = IN_SILICO_PCR.in_silico_pcr(primer_list, fasta_seq)
    # primers
    pair0 = primer_list[0]
    left0 = "5' - " + pair0[0] + " - 3'"
    right0 = "3' - " + pair0[1] + " - 5'"
    pair1 = primer_list[1]
    left1 = "5' - " + pair1[0] + " - 3'"
    right1 = "3' - " + pair1[1] + " - 5'"
    pair2 = primer_list[2]
    left2 = "5' - " + pair2[0] + " - 3'"
    right2 = "3' - " + pair2[1] + " - 5'"
    pair3 = primer_list[3]
    left3 = "5' - " + pair3[0] + " - 3'"
    right3 = "3' - " + pair3[1] + " - 5'"
    pair4 = primer_list[4]
    left4 = "5' - " + pair4[0] + " - 3'"
    right4 = "3' - " + pair4[1] + " - 5'"
    # products
    products0 = pcr_products[0]
    forward0 = "5' - " + products0[0] + " - 3'"
    reverse0 = "3' - " + products0[1] + " - 5'"
    lines0 = products0[2]
    products1 = pcr_products[1]
    forward1 = "5' - " + products1[0] + " - 3'"
    reverse1 = "3' - " + products1[1] + " - 5'"
    lines1 = products1[2]
    products2 = pcr_products[2]
    forward2 = "5' - " + products2[0] + " - 3'"
    reverse2 = "3' - " + products2[1] + " - 5'"
    lines2 = products2[2]
    products3 = pcr_products[3]
    forward3 = "5' - " + products3[0] + " - 3'"
    reverse3 = "3' - " + products3[1] + " - 5'"
    lines3 = products3[2]
    products4 = pcr_products[4]
    forward4 = "5' - " + products4[0] + " - 3'"
    reverse4 = "3' - " + products4[1] + " - 5'"
    lines4 = products4[2]

    return left0, right0, left1, right1, left2, right2, left3, right3, left4, right4, forward0, reverse0, lines0, forward1, reverse1, lines1, forward2, reverse2, lines2, forward3, reverse3, lines3, forward4, reverse4, lines4  

@app.route('/results', methods=['POST'])
def results():

    keyword = request.form.get('text')
    webpage = SEARCH.searcher(keyword)
    page_content = SEARCH.scraper1(webpage)
    results = SEARCH.links_results(page_content)
    length = len(results)

    if length == 0:
        return render_template('error.html')
    elif length == 1:
        left0, right0, left1, right1, left2, right2, left3, right3, left4, right4, forward0, reverse0, lines0, forward1, reverse1, lines1, forward2, reverse2, lines2, forward3, reverse3, lines3, forward4, reverse4, lines4 = analysis(page_content)
        return render_template('analysis.html', left0=left0, right0=right0, left1=left1, right1=right1, left2=left2, right2=right2, left3=left3, right3=right3, left4=left4, right4=right4, forward0=forward0, reverse0=reverse0, lines0=lines0, forward1=forward1, reverse1=reverse1, lines1=lines1, forward2=forward2, reverse2=reverse2, lines2=lines2, forward3=forward3, reverse3=reverse3, lines3=lines3, forward4=forward4, reverse4=reverse4, lines4=lines4)
    else:
        SEARCH.list_creator(results)
        return render_template('reference.html')

if __name__ == '__main__':
    app.run(debug=True)