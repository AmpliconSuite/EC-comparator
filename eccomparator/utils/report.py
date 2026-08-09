"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 5:40 PM 08/29/23

Create report with all the results
"""
import os
import datetime
from jinja2 import Environment, FileSystemLoader
from xhtml2pdf import pisa
from io import BytesIO

from eccomparator.utils.utils import OUTFILES as o

def generate_report(path_s1,
                    path_s2,
                    total_cost,
                    outdir=None):

    if outdir:
        template_dir = os.path.dirname(os.path.abspath(__file__))
        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template("/report/results.html")


        html_content = template.render(page_title_text='EC-comparator report',
                               description_time=datetime.datetime.now(),
                               description_s1=path_s1,
                               description_s2=path_s2,
                               description_output=outdir,
                               total_cost=total_cost,
                               path_total_cost=o.TOTAL_COST_PNG,
                               path_total_table=o.TOTAL_COST_TABLE,
                               path_plot_cn_sv=o.COVERAGE_BREAKPOINTS_PROFILE)

        # 5. Write the template to an HTML file
        outdirabs = os.path.abspath(outdir)
        outfile = os.path.join(outdirabs, 'report.html')
        with open(outfile, 'w') as f:
            f.write(html_content)

        # 6. Convert report to pdf
        outfile = os.path.join(outdirabs, 'report.pdf')
        os.chdir(outdirabs)

        out_pdf_file_handle = BytesIO()
        pisa.CreatePDF(
            src=BytesIO(html_content.encode('utf-8')),  # HTML to convert
            dest=out_pdf_file_handle)  # File handle to receive result

        # Save to a file
        with open(outfile, "wb") as f:
            f.write(out_pdf_file_handle.getvalue())

if __name__ == "__main__":
    generate_report("/Users/madag/Projects/PhD/github/ecdna-compare/examples/COLODM320/COLO320DM_Hung2021_amplicon3_cycles.bed",
                    "/Users/madag/Projects/PhD/github/ecdna-compare/examples/COLODM320/COLO320DM_Wu2019_amplicon4_cycles.bed",
                    total_cost=0.5,
                    outdir="/Users/madag/Projects/PhD/github/ecdna-compare/examples/COLODM320/output")