"""Render the single-source Russian riboflavin DBTL report and its figures."""
from html import escape
from pathlib import Path
import json
import re

import fitz
import matplotlib.pyplot as plt
from PIL import Image as PILImage, ImageDraw, ImageFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate, Frame, HRFlowable, Image, KeepTogether, PageBreak,
    PageTemplate, Paragraph, Spacer, Table, TableStyle,
)

# Reuse established project fonts; the legacy report is not built on import.
from build_fmn_report import MPL_FONT, ROOT_DIR, crop_structure_png

ROOT = ROOT_DIR
SOURCE = ROOT / 'docs/riboflavin-dbtl-report.md'
OUT = ROOT / 'riboflavin-dbtl-report.pdf'
ASSETS = ROOT / 'assets'
SCRATCH = ROOT / '.context/riboflavin-report'
INK = '#1E282D'
MUTED = '#5F6C72'
TEAL = '#1F634E'
RULE = '#C9D8D0'
WIDTH = A4[0] - 37 * mm
MARKET_DATA = ROOT / 'docs/riboflavin-market-data.json'


def market_figures():
    data = json.loads(MARKET_DATA.read_text(encoding='utf-8'))
    plt.rcParams.update({'font.size': 12, 'axes.labelcolor': INK,
                         'text.color': INK, 'axes.spines.top': False,
                         'axes.spines.right': False})
    rows = data['fefac_regional_shares']['rows']
    fig, ax = plt.subplots(figsize=(10, 3.6), dpi=180)
    positions = list(range(len(rows)))
    for offset, key, label, color in [(-.18, 'capacity_pct', 'Мощности', TEAL),
                                      (.18, 'production_pct', 'Выпуск', '#88A6BE')]:
        bars = ax.bar([p + offset for p in positions], [r[key] for r in rows],
                      width=.33, label=label, color=color)
        ax.bar_label(bars, labels=[f'{r[key]}%' for r in rows], padding=4)
    ax.set_xticks(positions, [r['region'] for r in rows])
    ax.set_ylim(0, 79)
    ax.set_ylabel('Доля мирового объёма, %')
    ax.legend(frameon=False, loc='upper right')
    fig.tight_layout()
    fig.savefig(ASSETS / 'riboflavin-market-capacity.png')
    plt.close(fig)

    rows = data['world_exports_2024']['rows']
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), dpi=180, sharey=True)
    for ax, key, divisor, label, color in [
        (axes[0], 'usd_thousand', 1000, 'Экспорт, млн долларов', TEAL),
        (axes[1], 'kg', 1000, 'Экспорт, тонн', '#88A6BE')]:
        values = [r[key] / divisor for r in rows]
        bars = ax.barh(positions := list(range(len(rows))), values, color=color, height=.58)
        ax.bar_label(bars, labels=[f'{v:,.1f}'.replace(',', ' ') for v in values], padding=4, fontsize=11)
        ax.set_xlim(0, max(values) * 1.28)
        ax.set_xlabel(label)
        ax.set_yticks(positions, [r['country'] for r in rows])
    axes[0].invert_yaxis()
    fig.tight_layout(w_pad=2)
    fig.savefig(ASSETS / 'riboflavin-market-exporters.png')
    plt.close(fig)

    rows = data['china_exports_to_russia']['rows']
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), dpi=180)
    for ax, key, divisor, label, color in [
        (axes[0], 'kg', 1000, 'Экспорт Китай → РФ, тонн', TEAL),
        (axes[1], 'usd_thousand', 1000, 'Стоимость, млн долларов', '#88A6BE')]:
        values = [r[key] / divisor for r in rows]
        bars = ax.bar([str(r['year']) for r in rows], values, width=.5, color=color)
        ax.bar_label(bars, labels=[f'{v:.2f}' for v in values], padding=4)
        ax.set_ylim(0, max(values) * 1.23)
        ax.set_ylabel(label)
    fig.tight_layout(w_pad=2)
    fig.savefig(ASSETS / 'riboflavin-market-china-russia.png')
    plt.close(fig)

    rows = data['vega_price_observations']['rows']
    fig, ax = plt.subplots(figsize=(10, 3.3), dpi=180)
    labels = []
    for y, row in enumerate(rows):
        fraction = row['assay_pct'] / 100
        low, high = row['low'] / fraction, row['high'] / fraction
        color = TEAL if row['assay_pct'] == 80 else '#60849F'
        ax.plot([low, high], [y, y], color=color, lw=7, solid_capstyle='round')
        ax.plot([low, high], [y, y], '|', color=color, markersize=16)
        value = f'{low:.1f}' if low == high else f'{low:.1f}–{high:.1f}'
        ax.text(high + 1.2, y, value, va='center', fontsize=11)
        labels.append(f"{row['period']} · {row['assay_pct']}%")
    ax.set_yticks(range(len(rows)), labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 65)
    ax.set_ylim(len(rows) - .5, -.5)
    ax.set_xlabel('Долларов за кг действующего B2 (пересчёт по концентрации)')
    ax.grid(axis='x', alpha=.15)
    fig.tight_layout()
    fig.savefig(ASSETS / 'riboflavin-market-prices.png')
    plt.close(fig)


def sampling_figure():
    fig, ax = plt.subplots(figsize=(10, 4.1), dpi=180)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    def box(x, y, text, color='#E6EFEA'):
        ax.text(x, y, text, ha='center', va='center', fontsize=11.5, color=INK,
                bbox=dict(boxstyle='round,pad=.62', facecolor=color, edgecolor='#BFD3C7', lw=1.0))

    def arrow(start, end):
        ax.annotate('', xy=end, xytext=start,
                    arrowprops=dict(arrowstyle='->', color=TEAL, lw=1.8))

    box(.5, .86, 'Репрезентативная проба культуры\nФиксированный исходный объём · защита от света')
    box(.24, .54, 'Аликвота А\nЦентрифугирование')
    box(.77, .54, 'Аликвота Б\nПолное извлечение из целой культуры')
    box(.13, .20, 'Супернатант\nРастворённые флавины')
    box(.42, .20, 'Осадок\nКлетки + кристаллы B2')
    box(.77, .20, 'Суммарный B2\nКонтроль материального баланса')
    arrow((.39, .76), (.25, .65))
    arrow((.61, .76), (.76, .65))
    arrow((.19, .43), (.13, .31))
    arrow((.29, .43), (.41, .31))
    arrow((.77, .43), (.77, .31))
    fig.tight_layout()
    fig.savefig(ASSETS / 'riboflavin-sampling.png')
    plt.close(fig)


def denovo_pathway_figure():
    fig, ax = plt.subplots(figsize=(10.5, 4.4), dpi=180, facecolor='white')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    def node(x, y, text, bg='#F2F6F4', edge='#BFD3C7', bold=False, fs=10.2):
        ax.text(x, y, text, ha='center', va='center', fontsize=fs,
                fontweight='bold' if bold else 'normal', color=INK,
                bbox=dict(boxstyle='round,pad=.48', facecolor=bg, edgecolor=edge, lw=1.15))

    def gene(x, y, text, key=False):
        bg = '#D9ECE3' if key else '#F7F9F8'
        edge = TEAL if key else '#C5D2CC'
        ax.text(x, y, text, ha='center', va='center', fontsize=8.8,
                fontweight='bold' if key else 'normal', color=TEAL if key else '#3D4C53',
                bbox=dict(boxstyle='round,pad=.28', facecolor=bg, edgecolor=edge, lw=1.0))

    def arrow(start, end, color=TEAL, style='->'):
        ax.annotate('', xy=end, xytext=start,
                    arrowprops=dict(arrowstyle=style, color=color, lw=1.6))

    # Top title banner / stoichiometry
    ax.text(.01, .95, 'Пуриновая ветвь (1 × GTP)', fontsize=10.5, fontweight='bold', color=TEAL, ha='left')
    ax.text(.01, .39, 'Пентозофосфатная ветвь (2 × Ру-5-Ф)', fontsize=10.5, fontweight='bold', color='#4C7290', ha='left')
    ax.text(.99, .95, 'Активатор SEF1 (C. famata) — координированное усиление пути',
            fontsize=9.2, color=TEAL, ha='right',
            bbox=dict(boxstyle='round,pad=.32', facecolor='#EBF4EF', edgecolor=TEAL, lw=0.9))

    # Purine branch (top row)
    node(.09, .77, 'GTP', bg='#EBF3F7', edge='#9BB8CD', bold=True)
    node(.33, .77, 'DARPP')
    node(.58, .77, 'ArPP → ArP')
    arrow((.145, .77), (.27, .77))
    gene(.208, .855, 'RIB1 (ключ.)', key=True)
    arrow((.39, .77), (.505, .77))
    gene(.448, .855, 'Cf: RIB2, RIB3\nKp: RIB7, RIB2')

    # Pentose phosphate branch (bottom row)
    node(.14, .21, 'Рибулозо-5-фосфат', bg='#EBF3F7', edge='#9BB8CD', bold=True)
    node(.54, .21, 'DHBP', bold=True)
    arrow((.265, .21), (.485, .21), color='#4C7290')
    gene(.375, .295, 'DHBP-синтаза\nCf: RIB6 · Kp: RIB3', key=True)

    # Convergence to DMRL and Riboflavin
    node(.74, .49, 'DMRL\n(люмазин)', bold=True)
    node(.925, .49, 'Рибофлавин\n(B2)', bg='#D5EAE0', edge=TEAL, bold=True, fs=10.8)

    arrow((.635, .72), (.715, .565))
    arrow((.595, .25), (.70, .425), color='#4C7290')
    gene(.61, .50, 'Люмазинсинтаза\nCf: RIB5 · Kp: RIB4')

    arrow((.80, .49), (.86, .49))
    gene(.83, .635, 'B2-синтаза\nCf: RIB7 · Kp: RIB5', key=True)

    # Stoichiometry footer
    ax.text(.5, .04, 'Стехиометрия: 1 GTP + 2 рибулозо-5-фосфата → 1 рибофлавин (B2)   |   Зелёным выделены приоритетные мишени усиления в DBTL',
            ha='center', va='center', fontsize=9.2, color=MUTED)
    fig.tight_layout()
    fig.savefig(ASSETS / 'riboflavin-denovo-pathway.png')
    plt.close(fig)


def figures():
    market_figures()
    sampling_figure()
    denovo_pathway_figure()
    with PILImage.open(ASSETS / 'candida-kinetics-2022.png') as original:
        original.crop((0, 3600, original.width, original.height)).save(
            ASSETS / 'candida-kinetics-2022-ef.png')
    fig = plt.figure(figsize=(11, 3.3), dpi=180, facecolor='white')
    for x, name, label in [(0.02, 'riboflavin', 'Рибофлавин · B2'),
                           (0.34, 'fmn', 'FMN · 5′-фосфат'),
                           (0.67, 'fad', 'FAD')]:
        ax = fig.add_axes([x, .23, .29, .65])
        ax.imshow(crop_structure_png(ASSETS / 'structures' / f'{name}.png'))
        ax.axis('off')
        fig.text(x + .145, .94, label, ha='center', fontsize=13, color=INK)
    fig.text(.326, .49, '→', ha='center', fontsize=25, color=TEAL)
    fig.text(.658, .49, '→', ha='center', fontsize=25, color=TEAL)
    fig.text(.32, .12, 'АТФ → АДФ\nРибофлавинкиназа', ha='center', fontsize=10)
    fig.text(.66, .12, 'АТФ → PPi\nFAD-синтетаза', ha='center', fontsize=10)
    fig.savefig(ASSETS / 'riboflavin-pathway.png', bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 2.5), dpi=180)
    volumes = [1250, 125, 62.5]
    bars = ax.barh(['1 г/л', '10 г/л', '20 г/л'], volumes, color=TEAL, height=.52)
    ax.invert_yaxis()
    for bar, value in zip(bars, volumes):
        ax.text(value + 18, bar.get_y() + bar.get_height() / 2,
                f'{value:g} л', va='center', fontsize=12)
    ax.set_xlim(0, 1450)
    ax.set_xlabel('Объём культуры на 1 кг выделенного B2 (извлечение 80%)')
    ax.set_ylabel('Титр')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(ASSETS / 'riboflavin-volume.png')
    plt.close(fig)


def _inline_markup(fragment):
    escaped = escape(fragment)
    escaped = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', escaped)
    escaped = re.sub(r'\*(.+?)\*', r'<i>\1</i>', escaped)
    return escaped


def markup(text):
    parts = re.split(r'(https?://[^\s)]+)', text)
    rendered = []
    for part in parts:
        if part.startswith(('https://', 'http://')):
            url = part.rstrip(';.,')
            rendered.append(f'<link href="{escape(url, quote=True)}" color="{TEAL}">{escape(url)}</link>')
            rendered.append(_inline_markup(part[len(url):]))
        else:
            rendered.append(_inline_markup(part))
    return ''.join(rendered)


BODY = ParagraphStyle('body', fontName='Sans', fontSize=9.3, leading=12.9,
                      textColor=colors.HexColor(INK), spaceAfter=4.4)
SUBTITLE = ParagraphStyle('subtitle', parent=BODY, fontSize=9.3, leading=12.5,
                          textColor=colors.HexColor(MUTED), spaceAfter=8.5)
H1 = ParagraphStyle('title', parent=BODY, fontName='SansBold', fontSize=20,
                    leading=25, spaceAfter=6)
H2 = ParagraphStyle('section', parent=BODY, fontName='SansBold', fontSize=14.0,
                    leading=18.0, textColor=colors.HexColor(TEAL),
                    spaceBefore=13, spaceAfter=5.5, keepWithNext=True)
H2_TOP = ParagraphStyle('section_top', parent=H2, spaceBefore=0)
H3 = ParagraphStyle('subhead', parent=BODY, fontName='SansBold', fontSize=10.6,
                    leading=14.0, textColor=colors.HexColor(INK),
                    spaceBefore=7.5, spaceAfter=3.5, keepWithNext=True)
TOC_TITLE = ParagraphStyle('toc_title', parent=BODY, fontName='SansBold', fontSize=11.5,
                           leading=15, textColor=colors.HexColor(INK),
                           spaceBefore=6, spaceAfter=5, keepWithNext=True)
TOC_ITEM = ParagraphStyle('toc_item', parent=BODY, fontSize=9.3, leading=12.8, spaceAfter=0)
TOC_PAGE = ParagraphStyle('toc_page', parent=TOC_ITEM, alignment=TA_RIGHT,
                          textColor=colors.HexColor(MUTED))
CELL = ParagraphStyle('cell', parent=BODY, fontSize=7.9, leading=10.2, spaceAfter=0)
CELL_HEAD = ParagraphStyle('cell_head', parent=CELL, fontName='SansBold',
                           textColor=colors.HexColor(INK))
CAPTION = ParagraphStyle('caption', parent=BODY, fontSize=8.2, leading=11.1,
                         textColor=colors.HexColor(MUTED), spaceAfter=5.5)
SMALL = ParagraphStyle('small', parent=BODY, fontSize=8.1, leading=11.2, spaceAfter=3.5)
BULLET = ParagraphStyle('bullet', parent=BODY, leftIndent=11, firstLineIndent=-8, spaceAfter=2.8)


class Report(BaseDocTemplate):
    def __init__(self, filename):
        super().__init__(str(filename), pagesize=A4, leftMargin=18.5 * mm,
                         rightMargin=18.5 * mm, topMargin=15.5 * mm, bottomMargin=15.5 * mm,
                         title='Рибофлавин в дрожжах: от воспроизведения базового процесса к оптимизации DBTL',
                         author='DBTL project')
        self.pages = {}
        frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height,
                      leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        self.addPageTemplates(PageTemplate(id='page', frames=[frame], onPage=self.footer))

    @staticmethod
    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont('Sans', 8)
        canvas.setFillColor(colors.HexColor('#77817F'))
        canvas.drawCentredString(A4[0] / 2, 9.5 * mm, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if hasattr(flowable, 'section_key'):
            self.canv.bookmarkPage(flowable.section_key)
            self.canv.addOutlineEntry(flowable.section_title, flowable.section_key, level=0)
            self.pages[flowable.section_key] = self.page


def _build_toc(headings, pages):
    items = [
        HRFlowable(width='100%', thickness=0.7, color=colors.HexColor(RULE), spaceBefore=6, spaceAfter=6),
        Paragraph('Содержание', TOC_TITLE),
    ]
    rows = []
    for n, heading in enumerate(headings, 1):
        key = f'section-{n}'
        rows.append([
            Paragraph(f'<link href="#{key}" color="{TEAL}">{escape(heading)}</link>', TOC_ITEM),
            Paragraph(str(pages.get(key, '—')), TOC_PAGE),
        ])
    toc_table = Table(rows, colWidths=[WIDTH - 14 * mm, 14 * mm], hAlign='LEFT')
    toc_table.setStyle(TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('LEFTPADDING', (0, 0), (-1, -1), 0),
        ('RIGHTPADDING', (0, 0), (-1, -1), 0),
        ('TOPPADDING', (0, 0), (-1, -1), 2.8),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 3.6),
        ('LINEBELOW', (0, 0), (-1, -2), 0.35, colors.HexColor('#E4ECE7')),
    ]))
    items.append(toc_table)
    return items


def story(text, pages):
    lines = text.splitlines()
    headings = [line[3:] for line in lines if line.startswith('## ')]
    result = []
    index = 0
    first_section = True
    seen_title = False
    seen_subtitle = False
    section = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if line.startswith('# '):
            result += [
                HRFlowable(width='100%', thickness=2.0, color=colors.HexColor(TEAL), spaceBefore=0, spaceAfter=9),
                Paragraph(markup(line[2:]), H1),
            ]
            seen_title = True
        elif line.startswith('## '):
            if first_section:
                result += _build_toc(headings, pages)
                first_section = False
            result.append(PageBreak())
            section += 1
            heading = Paragraph(markup(line[3:]), H2_TOP)
            heading.section_key = f'section-{section}'
            heading.section_title = line[3:]
            result.append(heading)
        elif line == '<!-- pagebreak -->':
            result.append(PageBreak())
        elif line.startswith('### '):
            result.append(Paragraph(markup(line[4:]), H3))
        elif line.startswith('|'):
            rows = []
            while index < len(lines) and lines[index].startswith('|'):
                values = [v.strip() for v in lines[index].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[:\- ]+', v) for v in values):
                    cell_style = CELL_HEAD if not rows else CELL
                    rows.append([Paragraph(markup(v), cell_style) for v in values])
                index += 1
            count = len(rows[0])
            widths = [WIDTH / count] * count
            if count == 2:
                widths = [WIDTH * .30, WIDTH * .70]
            elif count == 3:
                widths = [WIDTH * .25, WIDTH * .375, WIDTH * .375]
            elif count == 4:
                widths = [WIDTH * .34, WIDTH * .22, WIDTH * .22, WIDTH * .22]
            table = Table(rows, colWidths=widths, repeatRows=1, hAlign='LEFT')
            table.setStyle(TableStyle([
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E6EFEA')),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F6F9F7')]),
                ('LEFTPADDING', (0, 0), (-1, -1), 6.0),
                ('RIGHTPADDING', (0, 0), (-1, -1), 6.0),
                ('TOPPADDING', (0, 0), (-1, -1), 3.4),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 3.4),
                ('LINEBELOW', (0, 0), (-1, 0), .6, colors.HexColor('#B4C9BD')),
                ('LINEBELOW', (0, -1), (-1, -1), .4, colors.HexColor('#DCE5DF')),
            ]))
            block = [table, Spacer(1, 5.5)]
            result.append(KeepTogether(block) if len(rows) <= 6 else table)
            if len(rows) > 6:
                result.append(Spacer(1, 5.5))
            continue
        elif line.startswith('!['):
            match = re.fullmatch(r'!\[(.*?)\]\((.*?)\)', line)
            path = (SOURCE.parent / match[2]).resolve()
            with PILImage.open(path) as image:
                ratio = image.height / image.width
            image_width = min(WIDTH * 0.94, 140 * mm / ratio)
            img_flow = Image(str(path), width=image_width, height=image_width * ratio)
            # Keep figure together with its immediately following caption paragraph.
            lookahead = index + 1
            while lookahead < len(lines) and not lines[lookahead].strip():
                lookahead += 1
            if lookahead < len(lines) and lines[lookahead].strip().startswith('Рисунок '):
                chunks = [lines[lookahead].strip()]
                lookahead += 1
                while lookahead < len(lines) and lines[lookahead].strip() and not lines[lookahead].startswith(('#', '|', '- ', '![', '<!--')):
                    chunks.append(lines[lookahead].strip())
                    lookahead += 1
                cap = Paragraph(markup(' '.join(chunks)), CAPTION)
                group = [Spacer(1, 3), img_flow, Spacer(1, 3), cap]
                if result and isinstance(result[-1], Paragraph) and result[-1].style is H3:
                    group.insert(0, result.pop())
                result.append(KeepTogether(group))
                index = lookahead
                continue
            result += [img_flow, Spacer(1, 6)]
        elif line.startswith('- ') or re.match(r'^\d+\.\s', line):
            if line.startswith('- '):
                result.append(Paragraph('• ' + markup(line[2:]), BULLET))
            else:
                result.append(Paragraph(markup(line), BULLET))
        else:
            chunks = [line]
            while index + 1 < len(lines) and lines[index + 1].strip() and not lines[index + 1].startswith(('#', '|', '- ', '![', '<!--')) and not re.match(r'^\d+\.\s', lines[index + 1].strip()):
                index += 1
                chunks.append(lines[index].strip())
            content = ' '.join(chunks)
            if seen_title and not seen_subtitle and first_section:
                result.append(Paragraph(markup(content), SUBTITLE))
                seen_subtitle = True
            else:
                is_reference = re.match(r'^\[\d+\]', content)
                paragraph = Paragraph(markup(content), SMALL if is_reference else BODY)
                result.append(KeepTogether([paragraph]) if is_reference else paragraph)
        index += 1
    return result


def render_review():
    for old in SCRATCH.glob('page-*.png'):
        old.unlink()
    for old in SCRATCH.glob('contact-*.png'):
        old.unlink()
    doc = fitz.open(OUT)
    thumbnails = []
    for i, page in enumerate(doc):
        pix = page.get_pixmap(matrix=fitz.Matrix(1.2, 1.2), alpha=False)
        path = SCRATCH / f'page-{i + 1:02}.png'
        pix.save(path)
        with PILImage.open(path) as image:
            thumb = image.copy()
        thumb.thumbnail((230, 325))
        thumbnails.append(thumb)
    for start in range(0, len(thumbnails), 12):
        batch = thumbnails[start:start + 12]
        sheet = PILImage.new('RGB', (1000, ((len(batch) + 3) // 4) * 355), '#dce3e0')
        draw = ImageDraw.Draw(sheet)
        for n, thumb in enumerate(batch):
            x, y = (n % 4) * 250 + 10, (n // 4) * 355 + 20
            sheet.paste(thumb, (x, y))
            draw.text((x, y - 16), str(start + n + 1), fill='black')
        sheet.save(SCRATCH / f'contact-{start // 12 + 1}.png')
    cover = PILImage.new('RGB', (1600, 900), '#edf3ef')
    draw = ImageDraw.Draw(cover)
    font = 'C:/Windows/Fonts/segoeui.ttf'
    bold = 'C:/Windows/Fonts/segoeuib.ttf'
    draw.rectangle([(0, 0), (1600, 18)], fill=TEAL)
    draw.text((95, 110), 'DBTL · БИОСИНТЕЗ ВИТАМИНА B2', font=ImageFont.truetype(font, 32), fill=TEAL)
    draw.multiline_text((90, 230), 'Рибофлавин в дрожжах:\nвоспроизведение и DBTL',
                        font=ImageFont.truetype(bold, 78), spacing=22, fill=INK)
    draw.multiline_text((95, 620), 'Выбор продуцента · сравнение 16 альтернатив · аналитика и рынок\nПлан работ и ресурсная оценка · Ноябрь 2026',
                        font=ImageFont.truetype(font, 34), spacing=15, fill=MUTED)
    cover.save(SCRATCH / 'cover.webp')
    links = sum(len(page.get_links()) for page in doc)
    print(f'Built {OUT.name}: {len(doc)} pages, {links} links')
    print('Review renders:', SCRATCH)


def main():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    figures()
    text = SOURCE.read_text(encoding='utf-8')
    first = Report(SCRATCH / 'first-pass.pdf')
    first.build(story(text, {}))
    final = Report(OUT)
    final.build(story(text, first.pages))
    if final.pages != first.pages:
        raise RuntimeError('TOC page numbers changed between passes')
    validate_report(text)
    render_review()


def validate_report(text):
    """Runnable checks: source coverage, market arithmetic and final PDF geometry."""
    data = json.loads(MARKET_DATA.read_text(encoding='utf-8'))
    for key in ('capacity_pct', 'production_pct'):
        assert sum(row[key] for row in data['fefac_regional_shares']['rows']) == 100
    selected = data['selected_exports_to_russia_2024']['rows']
    assert sum(row['kg'] for row in selected) == 283667
    assert abs(sum(row['usd_thousand'] for row in selected) - 4440.27) < .001
    assert abs(1000 / (20 * .8) - 62.5) < 1e-9
    assert round(20 * .8 * 8000 / (200 + 24)) == 571
    body_part, _, refs_part = text.partition('## 13. Список источников')
    assert refs_part, 'Missing section 13 references'
    citations = set()
    for group in re.findall(r'\[([\d, –-]+)\]', body_part):
        for item in group.split(','):
            bounds = re.split('[–-]', item.strip())
            citations.update(range(int(bounds[0]), int(bounds[-1]) + 1))
    references = set(map(int, re.findall(r'^\[(\d+)\]', refs_part, re.M)))
    assert citations == references == set(range(1, 43)), (citations - references, references - citations)
    assert 'Предложение научному руководителю' not in text
    assert 'Редакция 2' not in text
    with fitz.open(OUT) as doc:
        toc = doc.get_toc()
        assert len(toc) == 13, toc
        toc_links = doc[0].get_links()
        assert len(toc_links) == 13, len(toc_links)
        destinations = [int(link['page']) if link['kind'] == 4 else link['page'] + 1
                        for link in toc_links]
        assert destinations == [item[2] for item in toc], (destinations, toc)
        for page in doc:
            extracted = page.get_text()
            assert '<!--' not in extracted
            assert '\ufffd' not in extracted
            assert '**' not in extracted
            for block in page.get_text('blocks'):
                assert block[0] >= 0 and block[1] >= 0, (page.number, block)
                assert block[2] <= page.rect.width + .5 and block[3] <= page.rect.height + .5, (page.number, block)
        print(f'PASS: {len(references)} cited references, 13 TOC links on page 1, market calculations, {len(doc)} pages verified')


if __name__ == '__main__':
    main()
