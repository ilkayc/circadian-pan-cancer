import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Set high-resolution publication styling
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

# Super compact dimensions: 14.8 x 6.6 inches (occupies minimal vertical space in Word, text appears big and readable)
fig, ax = plt.subplots(figsize=(14.8, 6.6), dpi=300, facecolor='#fafbfe')
ax.set_xlim(0, 14.8)
ax.set_ylim(0, 6.6)
ax.axis('off')

# Pastel Colors
c_stage1 = '#4575b4' # Soft Ocean / Periwinkle
c_stage2 = '#d97736' # Warm Terracotta / Apricot
c_stage3 = '#389163' # Soft Sage Green
c_stage4 = '#c04d69' # Soft Dusty Rose

c_bg1 = '#f4f8fd'
c_bg2 = '#fef7f2'
c_bg3 = '#f2faf5'
c_bg4 = '#fdf3f5'

c_border1 = '#cbdff8'
c_border2 = '#fcd9c4'
c_border3 = '#c2ebd5'
c_border4 = '#f8cbd5'

# Compact Card Helper with Big Readable Text
def draw_card(ax, x, y, w, h, title, items, header_color, border_color):
    # Subtle soft drop shadow
    shadow = patches.FancyBboxPatch((x+0.025, y-0.025), w, h, boxstyle='round,pad=0.10,rounding_size=0.10',
                                   facecolor='#334155', alpha=0.04, edgecolor='none', zorder=1)
    ax.add_patch(shadow)
    
    # White card body with pastel border
    card = patches.FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.10,rounding_size=0.10',
                                 facecolor='#ffffff', edgecolor=border_color, linewidth=1.2, zorder=2)
    ax.add_patch(card)
    
    # Header badge
    hb_h = 0.48
    hb = patches.FancyBboxPatch((x, y + h - hb_h), w, hb_h, boxstyle='round,pad=0.06,rounding_size=0.08',
                               facecolor=header_color, edgecolor='none', zorder=3)
    ax.add_patch(hb)
    # Header text: 10.2pt Bold
    ax.text(x + w/2, y + h - hb_h/2, title, ha='center', va='center', 
            fontsize=10.0, fontweight='bold', color='#ffffff', zorder=4)
    
    # Items: 9.3pt Bold/Clear Slate
    cur_y = y + h - hb_h - 0.25
    for item in items:
        # Bullet dot
        dot = patches.Circle((x + 0.18, cur_y), 0.052, facecolor=header_color, edgecolor='none', zorder=4)
        ax.add_patch(dot)
        ax.text(x + 0.30, cur_y, item, ha='left', va='center', 
                fontsize=9.2, color='#090d16', fontweight='normal', zorder=4)
        cur_y -= 0.29

# Banner helper
def draw_stage_banner(ax, x, y, w, h, step_num, title, color):
    banner = patches.FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.08,rounding_size=0.10',
                                   facecolor=color, edgecolor='none', zorder=2)
    ax.add_patch(banner)
    ax.text(x + w/2, y + h/2, f'{step_num}   |   {title}', ha='center', va='center', 
            fontsize=10.5, fontweight='bold', color='#ffffff', zorder=3)

# Column setup
start_x = 0.32
col_w = 3.36
gap = 0.24
y_top = 6.05
y_bot = 0.20
col_h = y_top - y_bot

stages = [
    ('PHASE 1', 'DATA ACQUISITION', c_stage1, c_bg1, c_border1, start_x),
    ('PHASE 2', 'VARIANT FILTERING', c_stage2, c_bg2, c_border2, start_x + (col_w + gap)),
    ('PHASE 3', 'FUNCTIONAL DISCOVERY', c_stage3, c_bg3, c_border3, start_x + 2 * (col_w + gap)),
    ('PHASE 4', 'CLINICAL PROGNOSTICS', c_stage4, c_bg4, c_border4, start_x + 3 * (col_w + gap))
]

for step_code, step_title, color, bg, border, x_pos in stages:
    col_bg = patches.FancyBboxPatch((x_pos, y_bot), col_w, col_h, boxstyle='round,pad=0.12,rounding_size=0.14',
                                   facecolor=bg, edgecolor=border, linewidth=1.0, alpha=0.7, zorder=0)
    ax.add_patch(col_bg)
    draw_stage_banner(ax, x_pos + 0.10, y_top - 0.54, col_w - 0.20, 0.46, step_code, step_title, color)

card_w = col_w - 0.24

# Phase 1
x1 = start_x + 0.12
draw_card(ax, x1, 3.70, card_w, 1.66, 'TCGA Pan-Cancer Atlas',
          ['10,239 Primary Tumour Profiles', '33 Malignant Histologies', 'Curated MC3 Consensus VCFs'],
          c_stage1, c_border1)
draw_card(ax, x1, 1.94, card_w, 1.64, '16 Core Clock Genes',
          ['Core TTFL (PER1-3, CRY1-2, CLOCK, ARNTL)', 'Auxiliary Loop (NR1D1-2, RORA-C)', 'Kinases (CSNK1D, CSNK1E, TIMELESS)'],
          c_stage1, c_border1)
draw_card(ax, x1, 0.30, card_w, 1.52, 'Clinical Resource (CDR)',
          ['Standardized TCGA-CDR Database', 'Validated Overall Survival (OS)', 'Progression-Free Interval (PFI)'],
          c_stage1, c_border1)

# Phase 2
x2 = start_x + (col_w + gap) + 0.12
draw_card(ax, x2, 3.70, card_w, 1.66, 'Somatic Mutation Calls',
          ['2,503 Total Somatic Alterations', '595 Synonymous Variants (23.8%)', '1,908 Non-Silent Mutations (Table 1)'],
          c_stage2, c_border2)
draw_card(ax, x2, 1.94, card_w, 1.64, 'Coding Missense Spectrum',
          ['1,532 Missense Substitutions (61.2%)', '136 Nonsense & 115 Frameshift', '39 Splice Site & 101 Other/Reg.'],
          c_stage2, c_border2)
draw_card(ax, x2, 0.30, card_w, 1.52, 'Dual In Silico Pathogenicity',
          ['SIFT (PROVEAN Suite, < 0.05)', 'PolyPhen-2 (HumVar, Damaging)', '706 Consensus High-Conf. (46.1%)'],
          c_stage2, c_border2)

# Phase 3
x3 = start_x + 2 * (col_w + gap) + 0.12
draw_card(ax, x3, 3.70, card_w, 1.66, 'Hotspots & Topography',
          ['Lollipop Mutation Distribution', 'PER3 p.R316C (Melanoma UV Sig)', 'CSNK1E/D p.R127W Hotspot'],
          c_stage3, c_border3)
draw_card(ax, x3, 1.94, card_w, 1.64, 'Domain Disruption',
          ['bHLH-PAS Dimerization Targets', 'Photolyase / FAD-binding Targets', 'Kinase Phosphorylation Loops'],
          c_stage3, c_border3)
draw_card(ax, x3, 0.30, card_w, 1.52, 'Landscape & Exclusivity',
          ['Oncoplot Waterfall Profiling (Fig 3-4)', 'Pairwise Analysis (120 Gene Pairs)', 'Multiple Testing FDR (q >= 0.81)'],
          c_stage3, c_border3)

# Phase 4
x4 = start_x + 3 * (col_w + gap) + 0.12
draw_card(ax, x4, 3.70, card_w, 1.66, 'Kaplan-Meier Survival',
          ['Significant Protection in UCEC (HR=0.43)', 'BLCA Survival Benefit (HR=0.66)', 'Adverse Prognosis in KIRC (HR=1.98)'],
          c_stage4, c_border4)
draw_card(ax, x4, 1.94, card_w, 1.64, 'Multi-Cohort Hazard Ratios',
          ['Forest Plot Across 10 Tumors', 'Systematic Cox Proportional Hazards', 'All 33 TCGA Cohorts (Supp. S4)'],
          c_stage4, c_border4)
draw_card(ax, x4, 0.30, card_w, 1.52, 'TMB Confounder Modeling',
          ['Firth Penalized Regression', 'Adjusted for log1p Somatic Burden', 'Robust Signal (All FDR < 0.05)'],
          c_stage4, c_border4)

# Arrows
def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='-|>', color='#94a3b8', lw=1.4, mutation_scale=10,
                                shrinkA=2, shrinkB=2), zorder=5)

arrow_ys = [4.53, 2.76, 1.06]
for idx in range(3):
    x_from = start_x + idx * (col_w + gap) + col_w
    x_to = start_x + (idx + 1) * (col_w + gap)
    for y in arrow_ys:
        draw_arrow(ax, x_from + 0.02, y, x_to - 0.02, y)

# Top Title Banner
main_title_bg = patches.FancyBboxPatch((start_x, 6.18), 14.16, 0.34, boxstyle='round,pad=0.05,rounding_size=0.06',
                                      facecolor='#0f172a', edgecolor='none', zorder=2)
ax.add_patch(main_title_bg)
ax.text(7.40, 6.35, 'Pan-Cancer Genomic, Structural, and Translational Workflow of the Human Circadian Clock Machinery',
        ha='center', va='center', fontsize=10.8, fontweight='bold', color='#ffffff', zorder=4)

plt.tight_layout()
out_png = 'figures/Figure_1_Pipeline_Workflow.png'
plt.savefig(out_png, dpi=300, bbox_inches='tight')
print('Ultra-Compact Figure 1 successfully generated at:', out_png)
