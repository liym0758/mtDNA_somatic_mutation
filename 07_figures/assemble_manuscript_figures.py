"""Assemble original vector panels, preserving their visual styles and source PDFs."""
from pathlib import Path
from io import BytesIO
import os
from pypdf import PdfReader,PdfWriter,PageObject,Transformation
from reportlab.pdfgen import canvas

ROOT=Path(os.environ['MITOEVO_FIGURE_ASSEMBLY_ROOT'])
PRIOR=Path(os.environ['MITOEVO_SOURCE_CODE_ROOT'])/'manuscript_revision_20261104_v1'
FIG=ROOT/'Figures'; SUP=ROOT/'Supplementary_candidates'
SP=['human','mouse','rat','chicken','cobra','sturgeon','stingray']

def put(dest,src,index,x,y,w,h):
 p=PdfReader(src).pages[index];sw=float(p.mediabox.width);sh=float(p.mediabox.height);s=min(w/sw,h/sh)
 dest.merge_transformed_page(p,Transformation().scale(s).translate(x+(w-sw*s)/2,y+(h-sh*s)/2))

def label(dest,items,w,h):
 b=BytesIO();c=canvas.Canvas(b,pagesize=(w,h))
 for x,y,t,size in items:c.setFont('Helvetica-Bold',size);c.drawString(x,y,t)
 c.save();b.seek(0);dest.merge_page(PdfReader(b).pages[0])

def save(page,path):
 w=PdfWriter();w.add_page(page)
 with path.open('wb') as f:w.write(f)

# Original Figure 1 slide, with only caption wording and reviewer-note removal in a copy.
src=SUP/'figure1_source_slide_caption_corrected.pdf';w=PdfWriter();w.add_page(PdfReader(src).pages[1])
with (FIG/'Figure1.pdf').open('wb') as f:w.write(f)

# Figure 2a-d: exact already-reviewed vector panel assembly; e: seven unchanged original cell types.
p=PageObject.create_blank_page(width=1420,height=1310)
put(p,PRIOR/'figures/Figure2_revised.pdf',0,15,80,695,1210)
hist=PRIOR/'baseline/AllSpecies_Reorder_celltype_VAF_hist_260128.pdf'
for i,sp in enumerate(SP):
 y=1110-i*165
 put(p,hist,i,750,y,650,160)
label(p,[(740,1281,'e  Cell-type VAF distributions (original all-positive analysis)',13)],1420,1310)
save(p,FIG/'Figure2.pdf')

# Figure 3a-g: unchanged lollipop, including original D-loop/intergenic track.
# Figure 3h-n: original R style, same all-positive counts, corrected y unit/label.
p=PageObject.create_blank_page(width=1200,height=2040)
lol=PRIOR/'baseline/7species_lolipop_SameYscale_DotScale_final_filter_260128.pdf'
for i,sp in enumerate(SP):
 y=1740-i*275
 put(p,lol,i,35,y,690,260)
 put(p,SUP/'figure3_vaf_original_style'/f'{sp}_VAF_all_positive_axis_corrected.pdf',0,730,y,445,260)
 label(p,[(20,y+246,chr(ord('a')+i),17),(730,y+246,chr(ord('h')+i),17)],1200,2040)
save(p,FIG/'Figure3.pdf')

# Figure 4a-g: original ggplot aesthetics from original 56-row table, with missing ns labels.
# Figure 4h: unchanged original overall-bar PDF and its existing significance logic.
p=PageObject.create_blank_page(width=1300,height=1240)
positions=[(0,3),(0,2),(0,1),(0,0),(1,3),(1,2),(1,1)]
for i,(col,row) in enumerate(positions):
 sp=SP[i];x=25+col*640;y=55+row*290
 put(p,SUP/'figure4_original_style_species'/f'{sp}_T_B_NK.pdf',0,x,y,285,260)
 put(p,SUP/'figure4_original_style_species'/f'{sp}_Mono.pdf',0,x+290,y,315,260)
 label(p,[(x,y+259,chr(ord('a')+i)+'  '+sp.title(),15)],1300,1240)
x=665;y=55
put(p,PRIOR/'baseline/AllSpecies_sitelevel_dNdS_noBin_cluster_barplot_withP.pdf',0,x,y,605,260)
label(p,[(x,y+259,'h  Overall dN/dS (original)',15)],1300,1240)
save(p,FIG/'Figure4.pdf')

for f in FIG.glob('*.pdf'):
 r=PdfReader(f);assert len(r.pages)==1 and r.pages[0].extract_text().strip();print(f.name,f.stat().st_size)
