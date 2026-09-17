import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

# 2x2 Grid Layout: Much narrower width, fits perfectly into Word without wide shrinking!
# Dimensions: 10.8 x 9.6 inches (Balanced 1.1:1 ratio)
fig, ax = plt.subplots(figsize=(10.8, 9.6), dpi=300, facecolor='#fafbfe')
ax.set_xlim(0, 10.8)
ax.set_ylim(0, 9.6)
ax.axis('off')

# Refined Pastel Colors
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

# Title banner at top
main_title_bg = patches.FancyBboxPatch((0.4, 9.0), 10.0, 0.42, boxstyle='round,pad=0.06,rounding_size=0.08',
                                      facecolor='#0f172a', edgecolor='none', zorder=2)
ax.add_patch(main_title_bg)
ax.text(5.4, 9.21, 'Pan-Cancer Genomic, Structural, and Translational Pipeline',
        ha='center', va='center', fontsize=11.2, fontweight='bold', color='#ffffff', zorder=4)

# Function to draw a Phase Box containing 3 sub-cards
def draw_phase_block(ax, x, y, w, h, phase_num, phase_title, cards_data, header_col, bg_col, border_col):
    # Phase container backdrop
    container = patches.FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.10,rounding_size=0.14',
                                      facecolor=bg_col, edgecolor=border_col, linewidth=1.2, zorder=1)
    ax.add_patch(container)
    
    # Phase top header
    ph_h = 0.44
    ph_bg = patches.FancyBboxPatch((x + 0.12, y + h - ph_h - 0.08), w - 0.24, ph_h, boxstyle='round,pad=0.06,rounding_size=0.08',
                                  facecolor=header_col, edgecolor='none', zorder=2)
    ax.add_patch(ph_bg)
    ax.text(x + w/2, y + h - ph_h/2 - 0.08, f'{phase_num}   |   {phase_title}', ha='center', va='center',
            fontsize=10.2, fontweight='bold', color='#ffffff', zorder=3)
    
    # 3 Sub-cards stacked inside
    sub_y = y + h - ph_h - 0.20
    card_h = (h - ph_h - 0.45) / 3.0
    
    for title, items in cards_data:
        sub_y -= card_h
        # Sub-card background
        c_box = patches.FancyBboxPatch((x + 0.16, sub_y), w - 0.32, card_h - 0.08, boxstyle='round,pad=0.06,rounding_size=0.08',
                                      facecolor='#ffffff', edgecolor=border_col, linewidth=1.0, zorder=3)
        ax.add_patch(c_box)
        
        # Sub-card header badge
        ch_h = 0.32
        ch_bg = patches.FancyBboxPatch((x + 0.16, sub_y + card_h - 0.08 - ch_h), w - 0.32, ch_h,
                                      boxstyle='round,pad=0.04,rounding_size=0.06',
                                      facecolor=header_col, edgecolor='none', zorder=4)
        ax.add_patch(ch_bg)
        ax.text(x + w/2, sub_y + card_h - 0.08 - ch_h/2, title, ha='center', va='center',
                fontsize=8.8, fontweight='bold', color='#ffffff', zorder=5)
        
        # Sub-card bullet items (LARGE and BOLD)
        cur_text_y = sub_y + card_h - 0.08 - ch_h - 0.16
        for item in items:
            dot = patches.Circle((x + 0.32, cur_text_y), 0.042, facecolor=header_col, edgecolor='none', zorder=5)
            ax.add_patch(dot)
            ax.text(x + 0.44, cur_text_y, item, ha='left', va='center',
                    fontsize=8.2, color='#0f172a', fontweight='normal', zorder=5)
            cur_text_y -= 0.21
        sub_y -= 0.04

# Layout coordinates: 2 columns, 2 rows
pw = 4.72
ph = 3.95
col1_x = 0.46
col2_x = 5.62
row1_y = 4.75 # Top row
row2_y = 0.50 # Bottom row

# Data
p1_cards = [
    ('TCGA Pan-Cancer Atlas', ['10,239 Primary Tumour Profiles', '33 Malignant Histologies (Curated MC3)']),
    ('16 Core Clock Genes', ['Core TTFL (PER1-3, CRY1-2, CLOCK, ARNTL)', 'Auxiliary Loop & Kinases (CSNK1D/E, etc.)']),
    ('Clinical Resource (CDR)', ['Standardized Overall Survival (OS) Records', 'Progression-Free Interval (PFI) Endpoints'])
]

p2_cards = [
    ('Somatic Mutation Calls', ['2,503 Total Somatic Alterations', '1,908 Non-Silent Mutations (Table 1)']),
    ('Coding Missense Spectrum', ['1,532 Missense Alterations (61.2%)', '136 Nonsense & 115 Frameshift Del/Ins']),
    ('Dual In Silico Pathogenicity', ['SIFT (< 0.05) & PolyPhen-2 (Damaging)', '706 Consensus High-Confidence Alleles'])
]

p3_cards = [
    ('Hotspots & Topography', ['Lollipop Mutation Distribution (16 Genes)', 'PER3 p.R316C (Melanoma) & CSNK1E p.R127W']),
    ('Conserved Domain Disruption', ['bHLH-PAS Dimerization Interface Targets', 'Photolyase/FAD & Kinase Catalytic Loops']),
    ('Landscape & Exclusivity', ['Oncoplot Waterfall Profiling (Fig 3-4)', 'Pairwise Analysis & FDR Correction (q >= 0.81)'])
]

p4_cards = [
    ('Kaplan-Meier Survival Analysis', ['Significant Protection in UCEC (HR=0.43)', 'BLCA Benefit (HR=0.66) vs. KIRC Risk (1.98)']),
    ('Multi-Cohort Hazard Ratios', ['Forest Plot Across 10 Solid Tumors', 'Systematic Cox Proportional Hazards Models']),
    ('TMB Confounder Modeling', ['Firth Bias-Reduced Penalized Regression', 'Adjusted for log1p Somatic Burden (FDR < 0.05)'])
]

# Draw blocks
draw_phase_block(ax, col1_x, row1_y, pw, ph, 'PHASE 1', 'DATA ACQUISITION', p1_cards, c_stage1, c_bg1, c_border1)
draw_phase_block(ax, col2_x, row1_y, pw, ph, 'PHASE 2', 'VARIANT FILTERING', p2_cards, c_stage2, c_bg2, c_border2)
draw_phase_block(ax, col1_x, row2_y, pw, ph, 'PHASE 3', 'FUNCTIONAL DISCOVERY', p3_cards, c_stage3, c_bg3, c_border3)
draw_phase_block(ax, col2_x, row2_y, pw, ph, 'PHASE 4', 'CLINICAL PROGNOSTICS', p4_cards, c_stage4, c_bg4, c_border4)

# Connecting Arrows
# Arrow 1: Phase 1 -> Phase 2 (Top horizontal)
ax.annotate('', xy=(col2_x - 0.04, row1_y + ph/2), xytext=(col1_x + pw + 0.04, row1_y + ph/2),
            arrowprops=dict(arrowstyle='-|>', color='#64748b', lw=2.0, mutation_scale=14))

# Arrow 2: Phase 2 -> Phase 4 (Vertical down) OR Phase 2 -> Phase 3:
# Let's show flow: Phase 1 -> Phase 2 -> Phase 3 -> Phase 4
# Downward from Phase 1 to Phase 3:
ax.annotate('', xy=(col1_x + pw/2, row2_y + ph + 0.04), xytext=(col1_x + pw/2, row1_y - 0.04),
            arrowprops=dict(arrowstyle='-|>', color='#64748b', lw=2.0, mutation_scale=14))

# Horizontal from Phase 3 -> Phase 4:
ax.annotate('', xy=(col2_x - 0.04, row2_y + ph/2), xytext=(col1_x + pw + 0.04, row2_y + ph/2),
            arrowprops=dict(arrowstyle='-|>', color='#64748b', lw=2.0, mutation_scale=14))

# Downward from Phase 2 -> Phase 4:
ax.annotate('', xy=(col2_x + pw/2, row2_y + ph + 0.04), xytext=(col2_x + pw/2, row1_y - 0.04),
            arrowprops=dict(arrowstyle='-|>', color='#64748b', lw=2.0, mutation_scale=14))

plt.tight_layout()
out_png = 'figures/Figure_1_Pipeline_Workflow_2X2.png'
plt.savefig(out_png, dpi=300, bbox_inches='tight')
print('2x2 Figure 1 successfully generated at:', out_png)
