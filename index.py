from flask import Flask, render_template, request, redirect, url_for, flash
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

    return pcr_products
    
@app.route('/results', methods=['POST'])
def results():
    keyword = request.form.get('text')
    webpage = SEARCH.searcher(keyword)
    page_content = SEARCH.scraper1(webpage)
    results = SEARCH.links_results(page_content)
    length = len(results)

    if length == 0:
        flash("No items found!")
        return redirect(url_for('keyword'))
    elif length == 1:
        all_results = analysis(page_content)
        return render_template('analysis.html', all_results=all_results)
    else:
        if length > 10:
            results2 = results[0:10]
            length_results2 = len(results2)
            return render_template('table_results.html', results2=results2, length_results2=length_results2)
        else:
            results2 = results
            length_results2 = len(results2)
            return render_template('table_results.html', results2=results2, length_results2=length_results2)

if __name__ == '__main__':
    app.run(debug=True)


