from PyPDF2 import PdfMerger
import os

files = os.listdir()

pdfs = []

for f in files:
    if f.endswith(".pdf"):
        pdfs.append(f)


merger = PdfMerger()

for pdf in pdfs:
    merger.append(pdf)

merger.write("result.pdf")
merger.close()


