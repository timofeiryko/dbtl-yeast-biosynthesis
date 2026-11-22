"""Render the current, single-source Russian report and its evidence-based figures."""
from html import escape
from pathlib import Path
import re

import fitz
import matplotlib.pyplot as plt
from PIL import Image as PILImage, ImageDraw, ImageFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, PageBreak, PageTemplate, Paragraph,
    Spacer, Table, TableStyle,
)

# Reuse established project fonts; the legacy report is not built on import.
from build_fmn_report import MPL_FONT, ROOT_DIR, crop_structure_png

ROOT = ROOT_DIR
SOURCE = ROOT / 'docs/riboflavin-dbtl-report.md'
OUT = ROOT / 'riboflavin-dbtl-report.pdf'
ASSETS = ROOT / 'assets'
SCRATCH = ROOT / '.context/riboflavin-report'
INK = '#202B30'
TEAL = '#236852'
WIDTH = A4[0] - 40 * mm


def figures():
    fig = plt.figure(figsize=(11, 3.4), dpi=180, facecolor='white')
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

    fig, ax = plt.subplots(figsize=(10, 2.6), dpi=180)
    titers = [1, 10, 20]
    volumes = [1250, 125, 62.5]
    bars = ax.barh(['1 г/л', '10 г/л', '20 г/л'], volumes, color=TEAL, height=.52)
    ax.invert_yaxis()
    for bar, value in zip(bars, volumes):
        ax.text(value + 18, bar.get_y() + bar.get_height()/2,
                f'{value:g} л', va='center', fontsize=12)
    ax.set_xlim(0, 1450)
    ax.set_xlabel('Объём культуры на 1 кг выделенного B2; извлечение 80%')
    ax.set_ylabel('Титр')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(ASSETS / 'riboflavin-volume.png')
    plt.close(fig)


def markup(text):
    parts = re.split(r'(https?://[^\s]+)', text)
    rendered = []
    for part in parts:
        if part.startswith(('https://', 'http://')):
            url = part.rstrip(';.,')
            rendered.append(f'<link href="{escape(url, quote=True)}" color="{TEAL}">{escape(url)}</link>')
            rendered.append(escape(part[len(url):]))
        else:
            rendered.append(escape(part))
    return ''.join(rendered)


BODY = ParagraphStyle('body', fontName='Sans', fontSize=10.1, leading=14.8,
                      textColor=colors.HexColor(INK), spaceAfter=7)
H1 = ParagraphStyle('title', parent=BODY, fontName='SansBold', fontSize=31,
                    leading=39, spaceAfter=20)
H2 = ParagraphStyle('section', parent=BODY, fontName='SansBold', fontSize=21,
                    leading=27, spaceAfter=16, keepWithNext=True)
H3 = ParagraphStyle('subhead', parent=BODY, fontName='SansBold', fontSize=12.2,
                    leading=17, spaceBefore=12, spaceAfter=7, keepWithNext=True)
CELL = ParagraphStyle('cell', parent=BODY, fontSize=8.2, leading=11.3, spaceAfter=0)
SMALL = ParagraphStyle('small', parent=BODY, fontSize=8.5, leading=12.5)
BULLET = ParagraphStyle('bullet', parent=BODY, leftIndent=11, firstLineIndent=-9)


class Report(BaseDocTemplate):
    def __init__(self, filename):
        super().__init__(str(filename), pagesize=A4, leftMargin=20*mm,
                         rightMargin=20*mm, topMargin=19*mm, bottomMargin=19*mm,
                         title='Рибофлавин: от воспроизведения к DBTL',
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
        canvas.drawCentredString(A4[0]/2, 10*mm, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if hasattr(flowable, 'section_key'):
            self.canv.bookmarkPage(flowable.section_key)
            self.canv.addOutlineEntry(flowable.section_title, flowable.section_key, level=0)
            self.pages[flowable.section_key] = self.page


def story(text, pages):
    lines = text.splitlines()
    headings = [line[3:] for line in lines if line.startswith('## ')]
    result = []
    index = 0
    first_section = True
    section = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if line.startswith('# '):
            result += [Spacer(1, 22*mm), Paragraph(markup(line[2:]), H1)]
        elif line.startswith('## '):
            if first_section:
                result += [PageBreak(), Paragraph('Содержание', H2), Spacer(1, 5*mm)]
                for n, heading in enumerate(headings, 1):
                    key = f'section-{n}'
                    row = Table([[Paragraph(f'<link href="#{key}" color="{TEAL}">{escape(heading)}</link>', BODY),
                                  Paragraph(str(pages.get(key, '—')), BODY)]],
                                colWidths=[WIDTH-12*mm, 12*mm])
                    row.setStyle(TableStyle([('VALIGN', (0, 0), (-1, -1), 'TOP'),
                                              ('BOTTOMPADDING', (0, 0), (-1, -1), 9)]))
                    result.append(row)
                first_section = False
            result.append(PageBreak())
            section += 1
            heading = Paragraph(markup(line[3:]), H2)
            heading.section_key = f'section-{section}'
            heading.section_title = line[3:]
            result.append(heading)
        elif line.startswith('### '):
            result.append(Paragraph(markup(line[4:]), H3))
        elif line.startswith('|'):
            rows = []
            while index < len(lines) and lines[index].startswith('|'):
                values = [v.strip() for v in lines[index].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[:\- ]+', v) for v in values):
                    rows.append([Paragraph(markup(v), CELL) for v in values])
                index += 1
            count = len(rows[0])
            widths = [WIDTH/count] * count
            if count == 3:
                widths = [WIDTH*.25, WIDTH*.375, WIDTH*.375]
            table = Table(rows, colWidths=widths, repeatRows=1, hAlign='LEFT')
            table.setStyle(TableStyle([
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E5EEE8')),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F6F8F7')]),
                ('LEFTPADDING', (0, 0), (-1, -1), 7),
                ('RIGHTPADDING', (0, 0), (-1, -1), 7),
                ('TOPPADDING', (0, 0), (-1, -1), 6.5),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 6.5),
                ('LINEBELOW', (0, 0), (-1, 0), .5, colors.HexColor('#B9CDC0')),
            ]))
            result += [table, Spacer(1, 10)]
            continue
        elif line.startswith('!['):
            match = re.fullmatch(r'!\[(.*?)\]\((.*?)\)', line)
            path = (SOURCE.parent / match[2]).resolve()
            with PILImage.open(path) as image:
                ratio = image.height/image.width
            result += [Image(str(path), width=WIDTH, height=WIDTH*ratio), Spacer(1, 8)]
        elif line.startswith('- '):
            result.append(Paragraph('• ' + markup(line[2:]), BULLET))
        else:
            chunks = [line]
            while index + 1 < len(lines) and lines[index+1].strip() and not lines[index+1].startswith(('#', '|', '- ', '![')):
                index += 1
                chunks.append(lines[index].strip())
            content = ' '.join(chunks)
            result.append(Paragraph(markup(content), SMALL if re.match(r'^\[\d+\]', content) else BODY))
        index += 1
    return result


def render_review():
    doc = fitz.open(OUT)
    thumbnails = []
    for i, page in enumerate(doc):
        pix = page.get_pixmap(matrix=fitz.Matrix(1.2, 1.2), alpha=False)
        path = SCRATCH / f'page-{i+1:02}.png'
        pix.save(path)
        with PILImage.open(path) as image:
            thumb = image.copy()
        thumb.thumbnail((230, 325))
        thumbnails.append(thumb)
    for start in range(0, len(thumbnails), 12):
        batch = thumbnails[start:start+12]
        sheet = PILImage.new('RGB', (1000, ((len(batch)+3)//4)*355), '#dce3e0')
        draw = ImageDraw.Draw(sheet)
        for n, thumb in enumerate(batch):
            x, y = (n % 4)*250+10, (n//4)*355+20
            sheet.paste(thumb, (x, y))
            draw.text((x, y-16), str(start+n+1), fill='black')
        sheet.save(SCRATCH / f'contact-{start//12+1}.png')
    # Feed-sized cover also remains useful for sharing outside the repository.
    cover = PILImage.new('RGB', (1600, 900), '#edf3ef')
    draw = ImageDraw.Draw(cover)
    font = 'C:/Windows/Fonts/segoeui.ttf'
    bold = 'C:/Windows/Fonts/segoeuib.ttf'
    draw.text((95, 100), 'DBTL / РИБОФЛАВИН', font=ImageFont.truetype(font, 32), fill=TEAL)
    draw.multiline_text((90, 235), 'Сначала воспроизвести.\nЗатем улучшать.',
                        font=ImageFont.truetype(bold, 83), spacing=20, fill=INK)
    draw.multiline_text((95, 620), 'Выбор дрожжей • протоколы • аналитика\nТитр + выход + скорость • ресурсы для запуска',
                        font=ImageFont.truetype(font, 36), spacing=15, fill=TEAL)
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
    render_review()


if __name__ == '__main__':
    main()