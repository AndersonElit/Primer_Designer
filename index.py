from flask import Flask, render_template, request, redirect, url_for, flash
from Primer_Designer_Functions import SEARCH, PRIMER_DESIGN, IN_SILICO_PCR
from werkzeug.utils import secure_filename
import os

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('fastaseq.html')

'''
@app.route('/keyword')
def keyword():
    return render_template('keyword.html')

@app.route('/geneid')
def geneid():
    return render_template('geneid.html')

@app.route('/fastaseq')
def fastaseq():
    return render_template('fastaseq.html')
'''

# def function
def analysis(fasta_seq):
    # Design primers
    dictionary = PRIMER_DESIGN.PRIMER3(fasta_seq)
    Primers_Tm_GC = PRIMER_DESIGN.primers_tm_gc(dictionary)
    pcr_products = IN_SILICO_PCR.in_silico_pcr(Primers_Tm_GC, fasta_seq)

    return pcr_products

'''    
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
        fasta_seq = PRIMER_DESIGN.get_fasta(page_content)
        all_results = analysis(fasta_seq)
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
'''

@app.route('/results_fasta', methods=['POST'])
def results_fasta():
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

if __name__ == '__main__':
    app.run(debug=True)


