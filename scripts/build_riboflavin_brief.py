# -*- coding: utf-8 -*-
"""
Build concise publication-grade A4 analytical PDF brief on Variant B:
Vitamin B2 (Riboflavin, Item 32 of Order No. 392) -> FMN & FAD coenzymes for the DBTL platform.

Generates:
  - assets/fig_b2_pathway.png (composited with RDKit 2D chemical structures)
  - assets/fig_b2_economics.png
  - riboflavin-fmn-fad-dbtl-brief.pdf
"""

from pathlib import Path

import fitz
import matplotlib
matplotlib.use("Agg")
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as PILImage, ImageOps

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

ROOT_DIR = Path(__file__).resolve().parent.parent
ASSETS_DIR = ROOT_DIR / "assets"
STRUCT_DIR = ASSETS_DIR / "structures"
CTX_DIR = ROOT_DIR / ".context" / "riboflavin-brief"
ASSETS_DIR.mkdir(parents=True, exist_ok=True)
CTX_DIR.mkdir(parents=True, exist_ok=True)

OUT_PDF = ROOT_DIR / "riboflavin-fmn-fad-dbtl-brief.pdf"


def _register_fonts() -> str:
    candidates = [
        (
            r"C:\Windows\Fonts\segoeui.ttf",
            r"C:\Windows\Fonts\segoeuib.ttf",
            r"C:\Windows\Fonts\segoeuii.ttf",
            r"C:\Windows\Fonts\segoeuiz.ttf",
            "Segoe UI",
        ),
        (
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf",
            "DejaVu Sans",
        ),
    ]
    for reg, bold, ital, bital, mpl_name in candidates:
        if Path(reg).exists():
            pdfmetrics.registerFont(TTFont("Sans", reg))
            pdfmetrics.registerFont(TTFont("SansBold", bold))
            pdfmetrics.registerFont(TTFont("SansItalic", ital))
            pdfmetrics.registerFont(TTFont("SansBoldItalic", bital))
            pdfmetrics.registerFontFamily(
                "Sans",
                normal="Sans",
                bold="SansBold",
                italic="SansItalic",
                boldItalic="SansBoldItalic",
            )
            return mpl_name
    raise RuntimeError("No suitable Cyrillic TrueType font found.")


MPL_FONT = _register_fonts()
plt.rcParams["font.family"] = MPL_FONT
plt.rcParams["axes.edgecolor"] = "#C8D0D6"
plt.rcParams["axes.linewidth"] = 0.8

INK = colors.HexColor("#181D20")
SUBINK = colors.HexColor("#2E373D")
MUTED = colors.HexColor("#5E6A73")
LINE = colors.HexColor("#D5DCE0")
SOFT_BG = colors.HexColor("#F5F7F8")
ACCENT_TEAL = colors.HexColor("#1B5E4A")
ACCENT_SOFT = colors.HexColor("#EAF2EF")
SLATE_BLUE = colors.HexColor("#2F5573")

SUB_MAP = {
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
    "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9", "ₐ": "a",
}
SUP_MAP = {"⁺": "+", "⁻": "-", "²": "2", "³": "3", "⁴": "4", "¹": "1"}


def fmt_chem(text: str) -> str:
    out = []
    for ch in text:
        if ch in SUB_MAP:
            out.append(f'<sub size="6.5">{SUB_MAP[ch]}</sub>')
        elif ch in SUP_MAP:
            out.append(f'<super size="6.5">{SUP_MAP[ch]}</super>')
        else:
            out.append(ch)
    return "".join(out)


def P(text: str, style) -> Paragraph:
    return Paragraph(fmt_chem(text), style)


class BookmarkHeading(Paragraph):
    def __init__(self, text: str, style: ParagraphStyle, key: str, level: int, toc_title: str):
        anchored = f'<a name="{key}"/>' + fmt_chem(text)
        super().__init__(anchored, style)
        self.toc_key = key
        self.toc_level = level
        self.toc_title = toc_title


def crop_structure_png(png_path: Path, pad: int = 14) -> np.ndarray:
    img = PILImage.open(png_path).convert("RGB")
    inv = ImageOps.invert(img)
    bbox = inv.getbbox()
    if bbox:
        x0 = max(0, bbox[0] - pad)
        y0 = max(0, bbox[1] - pad)
        x1 = min(img.width, bbox[2] + pad)
        y1 = min(img.height, bbox[3] + pad)
        img = img.crop((x0, y0, x1, y1))
    return np.asarray(img)


def make_pathway_figure(out_path: Path) -> None:
    """Figure 1: Direct biosynthetic cascade Glucose -> Riboflavin (Step 1) -> FMN -> FAD (Step 2)."""
    fig = plt.figure(figsize=(11.0, 5.0), dpi=300)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 45.5)
    ax.axis("off")

    outer = FancyBboxPatch(
        (0.6, 0.6), 98.8, 44.3,
        boxstyle="round,pad=0.3,rounding_size=1.1",
        facecolor="#FAFCFD", edgecolor="#CFD8DE", linewidth=1.0, zorder=0
    )
    ax.add_patch(outer)

    ax.text(
        2.0, 42.8,
        "ПРЯМОЙ БИОХИМИЧЕСКИЙ КАСКАД ВАРИАНТА Б: Шаг 1 (Рибофлавин, п. 32 МСХ № 392) → Шаг 2А (ФМН) → Шаг 2Б (ФАД)",
        fontsize=9.2, fontweight="bold", color="#181D20", va="center"
    )

    # Feeder precursor box (PPP + Purine IMP fork)
    feeder = FancyBboxPatch(
        (2.0, 32.8), 96.0, 8.2,
        boxstyle="round,pad=0.2,rounding_size=0.8",
        facecolor="#EEF4F8", edgecolor="#9EB4C4", linewidth=0.9, zorder=1
    )
    ax.add_patch(feeder)
    ax.text(
        3.2, 39.2,
        "Центральный узел предшественников из глюкозы (задача балансировки потоков в контуре GEM + ML-DBTL):",
        fontsize=8.0, fontweight="bold", color="#2F5573", va="center"
    )
    ax.text(
        3.2, 36.8,
        "[Ветвь 1: Пентозофосфатный путь ZWF1 / GND1 / RPE] → Рибулозо-5-фосфат (2 Ru5P) + PRPP (предшественник пуринового кольца)",
        fontsize=7.25, color="#2E373D", va="center"
    )
    ax.text(
        3.2, 34.5,
        "[Ветвь 2: Развилка пуринов в точке ИМФ (IMP)] → (а) Гуаниловая ветвь guaB/guaA → 1 GTP (в каскад RIB1–RIB7 на сборку B₂)   |   "
        "(б) Адениловая ветвь purA/ADE12 → 2 ATP (на Шаг 2)",
        fontsize=7.15, color="#2E373D", va="center"
    )

    def place_card(x, y, w, h, img_key, title, formula, sub, border, bg, title_col):
        card = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.2,rounding_size=0.9",
            facecolor=bg, edgecolor=border, linewidth=1.35, zorder=2
        )
        ax.add_patch(card)
        ax.text(x + w / 2, y + h - 1.8, title, fontsize=8.2, fontweight="bold", color=title_col, ha="center", va="center", zorder=4)
        ax.text(x + w / 2, y + h - 3.7, formula, fontsize=6.9, color="#5E6A73", ha="center", va="center", zorder=4)

        arr = crop_structure_png(STRUCT_DIR / f"{img_key}.png", pad=16)
        ih, iw = arr.shape[:2]
        bot_pad = 4.4
        box_w = w - 2.0
        box_h = h - 4.9 - bot_pad
        aspect = iw / ih
        if box_w / box_h > aspect:
            dh = box_h
            dw = dh * aspect
        else:
            dw = box_w
            dh = dw / aspect
        cx = x + w / 2
        cy = y + bot_pad + box_h / 2
        ax.imshow(arr, extent=[cx - dw / 2, cx + dw / 2, cy - dh / 2, cy + dh / 2], zorder=3, interpolation="lanczos")
        ax.text(x + w / 2, y + 2.25, sub, fontsize=6.65, color="#2E373D", ha="center", va="center", linespacing=1.22, zorder=4)

    yc, hc = 2.4, 26.8
    place_card(
        2.0, yc, 25.2, hc, "riboflavin",
        "Шаг 1: Рибофлавин (Вит. B2)",
        "C₁₇H₂₀N₄O₆ · 376,4 г/моль",
        "Импорт РФ: 300–500 т/год (п. 32 МСХ № 392)\nТест за 10 с: супернатант 445 нм без реактивов",
        "#1B5E4A", "#EAF2EF", "#1B5E4A"
    )
    place_card(
        36.6, yc, 25.6, hc, "fmn",
        "Шаг 2А: 5'-ФМН (мононуклеотид)",
        "C₁₇H₂₁N₄O₉P · 456,3 г/моль",
        "Строго 5'-фосфат (>96% vs ~75% хим. POCl₃)\nАФС для «Цитофлавина», «Декстралинка», инъекций",
        "#2F5573", "#F4F8FB", "#2F5573"
    )
    place_card(
        71.6, yc, 26.4, hc, "fad",
        "Шаг 2Б: Кофермент ФАД",
        "C₂₇H₃₃N₉O₁₅P₂ · 785,6 г/моль",
        "Химсинтез невозможен ($1 200–3 000+/кг)\nОфтальмология и тест-полоски глюкометров FAD-GDH",
        "#2F5573", "#FFFFFF", "#2F5573"
    )

    def rxn_arrow(x0, x1, y, top_t, bot_t):
        arr = FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>,head_length=5,head_width=3", color="#1B5E4A", linewidth=1.6, zorder=5)
        ax.add_patch(arr)
        xm = (x0 + x1) / 2
        ax.text(xm, y + 1.4, top_t, fontsize=7.3, fontweight="bold", color="#1B5E4A", ha="center", va="bottom")
        ax.text(xm, y - 1.3, bot_t, fontsize=6.25, color="#2E373D", ha="center", va="top", linespacing=1.18)

    rxn_arrow(27.6, 36.2, 15.8, "+ ген FMN1", "+1 АТФ → АДФ\n(рибофлавин-\nкиназа)")
    rxn_arrow(62.6, 71.2, 15.8, "+ ген FAD1", "+1 АТФ → PPᵢ\n(ФАД-синтаза)")

    # Feeder vertical arrows
    fa1 = FancyArrowPatch((14.6, 32.7), (14.6, 29.5), arrowstyle="-|>,head_length=4.2,head_width=2.5", color="#5E6A73", linewidth=1.0, linestyle="--", zorder=4)
    fa2 = FancyArrowPatch((31.9, 32.7), (31.9, 21.0), arrowstyle="-|>,head_length=4.2,head_width=2.5", color="#5E6A73", linewidth=1.0, linestyle="--", zorder=4)
    fa3 = FancyArrowPatch((66.9, 32.7), (66.9, 21.0), arrowstyle="-|>,head_length=4.2,head_width=2.5", color="#5E6A73", linewidth=1.0, linestyle="--", zorder=4)
    for f_a in (fa1, fa2, fa3):
        ax.add_patch(f_a)
    ax.text(15.8, 31.0, "RIB1–RIB7 (1 GTP + 2 Ru5P)", fontsize=6.6, color="#2F5573", fontweight="bold", va="center")

    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_economics_figure(out_path: Path) -> None:
    """Figure 2: Value escalation per kg along the cascade & isomer purity advantage of bio-FMN."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 3.65), dpi=300)
    fig.patch.set_facecolor("white")

    stages = [
        "Кормовой B₂ (80%)\nШаг 1 · п. 32 МСХ",
        "Пищ. E101 / Фарм B₂\n(98%, ГФ XV)",
        "Био-5'-ФМН (>96%)\nШаг 2А · фарм АФС",
        "Кофермент ФАД\nШаг 2Б · фарм / диагн.",
    ]
    prices = [17.0, 48.0, 420.0, 1950.0]
    lbls = ["$12–22/кг\n(сотни т в РФ)", "$35–60/кг\n(десятки т в РФ)", "$250–600/кг\n(чистый 5'-изомер)", "$1 200–3 000+/кг\n(только биосинтез)"]
    cols = ["#7A8B96", "#2F5573", "#1B5E4A", "#1B5E4A"]

    y_pos = np.arange(len(stages))
    ax1.barh(y_pos, prices, height=0.56, color=cols, alpha=0.92, zorder=3)
    ax1.set_xscale("log")
    ax1.set_xlim(5, 16000)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(stages, fontsize=7.8, color="#2E373D")
    for idx, (val, note) in enumerate(zip(prices, lbls)):
        ax1.text(val * 1.25, idx, note, va="center", fontsize=7.2, fontweight="bold" if idx >= 2 else "normal", color="#1B5E4A" if idx >= 2 else "#2E373D")

    ax1.grid(True, axis="x", which="major", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax1.set_xlabel("Оптовая цена B2B в мире ($/кг, лог. шкала)", fontsize=8.2, color="#2E373D")
    ax1.set_title("А. Рост удельной стоимости вдоль прямого каскада (в ~100 раз)", fontsize=9.1, fontweight="bold", color="#181D20", pad=8)
    ax1.tick_params(axis="x", labelsize=8, colors="#4A555E")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Stacked bar comparing chemical POCl3 FMN vs enzymatic/microbial FMN1
    routes = ["Химсинтез ФМН\n(рибофлавин + POCl₃)", "Биосинтез Шага 2А\n(рибофлавинкиназа FMN1)"]
    main_5fmn = np.array([75.0, 96.5])
    side_isomers = np.array([25.0, 3.5])

    x = np.arange(len(routes))
    ax2.bar(x, main_5fmn, width=0.52, color="#1B5E4A", label="Целевой активный 5'-ФМН (%)", zorder=3)
    ax2.bar(x, side_isomers, width=0.52, bottom=main_5fmn, color="#C86454", label="Побочные 3'/4'-фосфаты, дифосфаты и B₂ (%)", zorder=3)

    ax2.text(0, 38, "~73–77%\n5'-ФМН", ha="center", va="center", fontsize=8.0, fontweight="bold", color="white", zorder=4)
    ax2.text(0, 87.5, "~23–27%\nпримесей", ha="center", va="center", fontsize=7.4, fontweight="bold", color="white", zorder=4)
    ax2.text(1, 48, ">95–98%\nчистый\n5'-ФМН", ha="center", va="center", fontsize=8.0, fontweight="bold", color="white", zorder=4)

    ax2.set_ylim(0, 118)
    ax2.set_xticks(x)
    ax2.set_xticklabels(routes, fontsize=8.0, color="#2E373D")
    ax2.set_ylabel("Доля компонентов в субстанции (%)", fontsize=8.2, color="#2E373D")
    ax2.set_title("Б. Качественное преимущество био-ФМН над химией POCl₃", fontsize=9.1, fontweight="bold", color="#181D20", pad=8)
    ax2.legend(frameon=False, fontsize=7.2, loc="upper center")
    ax2.grid(True, axis="y", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax2.tick_params(axis="y", labelsize=8, colors="#4A555E")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout(pad=1.2)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


class BriefDocTemplate(BaseDocTemplate):
    def __init__(self, filename, toc_map: dict, **kw):
        super().__init__(filename, **kw)
        self.toc_map = toc_map

    def afterFlowable(self, flowable):
        if isinstance(flowable, BookmarkHeading):
            self.canv.bookmarkPage(flowable.toc_key)
            self.canv.addOutlineEntry(flowable.toc_title, flowable.toc_key, level=flowable.toc_level - 1, closed=False)
            self.toc_map[flowable.toc_key] = self.page


def build_pdf() -> None:
    fig_p = ASSETS_DIR / "fig_b2_pathway.png"
    fig_e = ASSETS_DIR / "fig_b2_economics.png"
    if not fig_p.exists():
        make_pathway_figure(fig_p)
    if not fig_e.exists():
        make_economics_figure(fig_e)

    page_w, page_h = A4
    left_m = right_m = 18 * mm
    top_m = 18.5 * mm
    bot_m = 16 * mm
    usable_w = page_w - left_m - right_m

    st_kicker = ParagraphStyle("Kicker", fontName="SansBold", fontSize=8.0, leading=10.2, textColor=ACCENT_TEAL, spaceAfter=1.8 * mm)
    st_title = ParagraphStyle("Title", fontName="SansBold", fontSize=16.5, leading=20.8, textColor=INK, spaceAfter=2.2 * mm)
    st_sub = ParagraphStyle("Sub", fontName="Sans", fontSize=9.3, leading=13.0, textColor=MUTED, spaceAfter=3.8 * mm)

    st_h1 = ParagraphStyle("H1", fontName="SansBold", fontSize=12.2, leading=15.5, textColor=INK, spaceBefore=3.0 * mm, spaceAfter=1.8 * mm, keepWithNext=True)
    st_h2 = ParagraphStyle("H2", fontName="SansBold", fontSize=9.9, leading=13.0, textColor=ACCENT_TEAL, spaceBefore=2.2 * mm, spaceAfter=1.1 * mm, keepWithNext=True)

    st_body = ParagraphStyle("Body", fontName="Sans", fontSize=9.05, leading=13.0, textColor=INK, alignment=TA_JUSTIFY, spaceAfter=1.6 * mm)
    st_bullet = ParagraphStyle("Bullet", fontName="Sans", fontSize=8.95, leading=12.7, textColor=INK, leftIndent=4.5 * mm, firstLineIndent=-4.5 * mm, alignment=TA_JUSTIFY, spaceAfter=1.2 * mm)
    st_box = ParagraphStyle("Box", fontName="Sans", fontSize=8.95, leading=12.8, textColor=INK, alignment=TA_LEFT)
    st_caption = ParagraphStyle("Caption", fontName="SansItalic", fontSize=8.0, leading=11.0, textColor=MUTED, spaceBefore=1.4 * mm, spaceAfter=3.0 * mm)

    st_toc_h = ParagraphStyle("TOCH", fontName="SansBold", fontSize=10.2, leading=13.0, textColor=INK, spaceBefore=2.2 * mm, spaceAfter=1.5 * mm)
    st_toc_1 = ParagraphStyle("TOC1", fontName="SansBold", fontSize=8.85, leading=11.6, textColor=INK)
    st_toc_pg = ParagraphStyle("TOCPg", fontName="SansBold", fontSize=8.85, leading=11.6, textColor=ACCENT_TEAL, alignment=2)

    st_th = ParagraphStyle("TH", fontName="SansBold", fontSize=7.7, leading=9.9, textColor=INK)
    st_td = ParagraphStyle("TD", fontName="Sans", fontSize=7.6, leading=10.2, textColor=SUBINK)
    st_td_b = ParagraphStyle("TDB", fontName="SansBold", fontSize=7.6, leading=10.2, textColor=ACCENT_TEAL)
    st_ref = ParagraphStyle("Ref", fontName="Sans", fontSize=7.7, leading=10.4, textColor=SUBINK, leftIndent=6.0 * mm, firstLineIndent=-6.0 * mm, spaceAfter=0.8 * mm)

    def draw_hf(canvas, doc):
        canvas.saveState()
        if doc.page > 1:
            canvas.setStrokeColor(LINE)
            canvas.setLineWidth(0.5)
            canvas.line(left_m, page_h - 12.8 * mm, page_w - right_m, page_h - 12.8 * mm)
            canvas.setFont("Sans", 7.5)
            canvas.setFillColor(MUTED)
            canvas.drawString(left_m, page_h - 10.8 * mm, "Вариант Б: каскад Витамин B₂ (рибофлавин) → ФМН → ФАД для платформы DBTL")
            canvas.drawRightString(page_w - right_m, page_h - 10.8 * mm, "Ноябрь 2026")
        canvas.setStrokeColor(LINE)
        canvas.setLineWidth(0.45)
        canvas.line(left_m, 11.5 * mm, page_w - right_m, 11.5 * mm)
        canvas.setFont("Sans", 8.2)
        canvas.setFillColor(MUTED)
        canvas.drawCentredString(page_w / 2, 7.6 * mm, f"{doc.page}")
        canvas.restoreState()

    toc_spec = [
        ("s1", "1. Проверка по пяти критериям и сопоставление с платформой L-эрготионеина"),
        ("s2", "2. Готовые штаммы-референсы, отличие от базовых дрожжей и аудит работ по DBTL/ML"),
        ("s3", "3. Прямой биохимический каскад (B2 → ФМН → ФАД) и развилка пуринов для моделирования"),
        ("s4", "4. Быстрая дифференциальная фотометрия без масс-спектрометра (Test) и план циклов"),
        ("s5", "5. Рынок РФ (п. 32 МСХ № 392, реестр ГРЛС), преимущество над химией и экономика"),
        ("s6", "6. Список литературы и нормативных источников"),
    ]

    def callout(html_text: str) -> Table:
        t = Table([[P(html_text, st_box)]], colWidths=[usable_w])
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), ACCENT_SOFT),
                    ("LINEBEFORE", (0, 0), (0, -1), 2.6, ACCENT_TEAL),
                    ("LEFTPADDING", (0, 0), (-1, -1), 9),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 9),
                    ("TOPPADDING", (0, 0), (-1, -1), 6),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                ]
            )
        )
        return t

    def build_story(toc_pages: dict):
        story = []

        def H1(key: str, text: str):
            return BookmarkHeading(text, st_h1, key, 1, text)

        # ==================== PAGE 1 ====================
        story.append(P("АНАЛИТИЧЕСКАЯ ЗАПИСКА · АЛЬТЕРНАТИВНАЯ ПЛАТФОРМА DBTL (ВАРИАНТ Б) · НОЯБРЬ 2026", st_kicker))
        story.append(P("Витамин B2 (рибофлавин) → коферменты ФМН и ФАД: платформа DBTL с тоннажным первым шагом из Приказа Минсельхоза № 392", st_title))
        story.append(P("Краткое обоснование выбора, прямая химическая преемственность Шага 1 и Шага 2, архитектура моделирования и безреактивной фотометрии, анализ ГРЛС и экономические границы масштабирования в РФ", st_sub))

        hr = Table([[""]], colWidths=[usable_w], rowHeights=[0.5])
        hr.setStyle(TableStyle([("LINEABOVE", (0, 0), (-1, -1), 0.9, INK), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2 * mm)]))
        story.append(hr)

        story.append(
            callout(
                "<b>Итоговый вывод.</b> Связка <b>Витамин B₂ (рибофлавин, Шаг 1) → рибофлавин-5'-мононуклеотид (ФМН, Шаг 2А) → флавинадениндинуклеотид (ФАД, Шаг 2Б)</b> закрывает главные узкие места маршрута через эрготионеин: (1) уже на Шаге 1 импорт в РФ составляет <b>300–500 тонн/год</b> (п. 32 Приказа Минсельхоза № 392, субсидия до 20% CAPEX с 2027 г., свой выпуск в РФ — 0%); (2) <b>в базовых дрожжах уже есть все 6 генов пути B₂ (<i>RIB1–RIB7</i>) и оба гена Шага 2 (<i>FMN1, FAD1</i>)</b> — чужеродные гены не нужны: замена промотора только перед собственным <i>RIB1</i> (+<i>RIB3</i>) сразу окрашивает среду в жёлтый цвет (4–20 мг/л в колбах, 175 мг/л в реакторе), а в БРЦ ВКПМ (Курчатовский ин-т) доступны готовые штаммы <i>B. subtilis</i> (<i>B-2694, B-6797, B-8387</i>, до 12–21 г/л); (3) B₂ является <b>прямым химическим интермедиатом</b> высокомаржинальных <b>5'-ФМН ($250–600/кг)</b> и <b>ФАД ($1 200–3 000+/кг)</b> при замере супернатанта за <b>10 с на 445 нм без реактивов</b> и <b>нуле опубликованных работ по ML-DBTL</b> в дрожжах и на каскаде ФМН/ФАД."
            )
        )

        story.append(P("Содержание (разделы кликабельны)", st_toc_h))
        toc_rows = []
        for k, t_str in toc_spec:
            pg = toc_pages.get(k, "—")
            toc_rows.append([
                Paragraph(f'<a href="#{k}" color="#181D20">{fmt_chem(t_str)}</a>', st_toc_1),
                Paragraph(f'<a href="#{k}" color="#1B5E4A">{pg}</a>', st_toc_pg),
            ])
        toc_tbl = Table(toc_rows, colWidths=[usable_w - 12 * mm, 12 * mm])
        toc_tbl.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 1),
            ("RIGHTPADDING", (0, 0), (-1, -1), 1),
            ("TOPPADDING", (0, 0), (-1, -1), 1.8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.8),
            ("LINEBELOW", (0, 0), (-1, -1), 0.3, colors.HexColor("#E6ECEF")),
        ]))
        story.append(toc_tbl)

        # ==================== SECTION 1 ====================
        story.append(H1("s1", "1. Проверка по пяти критериям и сопоставление с платформой L-эрготионеина"))
        story.append(
            P(
                "Флавиновый каскад полностью укладывается в согласованную систему из пяти критериев отбора целевой молекулы для цикла DBTL, предлагая альтернативный баланс между масштабом рынка первого шага и ценовой конкуренцией с крупнотоннажной ферментацией КНР (Таблица 1):",
                st_body,
            )
        )

        comp_rows = [
            [
                P("Критерий отбора", st_th),
                P("<b>Базовый вариант:</b><br/>L-Эрготионеин → (S,S)-SAMe", st_th),
                P("<b>Вариант Б (флавиновый каскад):</b><br/>Витамин B₂ (рибофлавин) → 5'-ФМН → ФАД", st_th),
            ],
            [
                P("<b>1. Понятный Шаг 1 и референсы</b>", st_td_b),
                P("Вставка 2 чужеродных генов (<i>NcEgt1, CpEgt2</i>) в <i>S. cerevisiae</i>, сразу 46–154 мг/л.", st_td),
                P("<b>0 чужеродных генов в дрожжах:</b> свои <i>RIB1–RIB7</i> под <i>pGAP</i> сразу дают 4–175 мг/л [17]; в БРЦ ВКПМ есть готовые <i>B. subtilis</i> (до 21 г/л) [14].", st_td),
            ],
            [
                P("<b>2. Быстрый оптический тест без МС</b>", st_td_b),
                P("УФ 311 нм (нужна наработка фермента <i>TdETL</i>) или 343 нм (2-PDS при pH 1,0).", st_td),
                P("<b>Абсолютный минимум операций:</b> прямое поглощение жёлтого супернатанта при <b>445 нм за 10 с без реактивов</b>; для ФАД/ФМН — тест тушения флуоресценции 525 нм [11].", st_td),
            ],
            [
                P("<b>3. Рынок РФ на Шаге 1 и биосинтез</b>", st_td_b),
                P("Микрониша <b>единицы–десятки кг/год</b> в РФ ($400–580/кг); монетизация Шага 1 только через R&amp;D-пилот с косм. брендом.", st_td),
                P("Импорт в РФ <b>300–500 тонн/год</b> (100% микробный синтез в мире, химсинтез вытеснен); но цена кормового B₂ низкая ($12–22/кг) — требует крупного масштаба [12].", st_td),
            ],
            [
                P("<b>4. Связь с Шагом 2 на том же узле</b>", st_td_b),
                P("Непрямая: переносится конвейер и <i>SAM2/АТФ</i>, но ветви <i>HIS1</i> и <i>CYS4</i> нужно разворачивать.", st_td),
                P("<b>Прямая химическая:</b> B₂ + 1 ген (<i>FMN1</i>) → <b>5'-ФМН</b> ($250–600/кг) + 1 ген (<i>FAD1</i>) → <b>ФАД</b> ($1 200–3 000+/кг) без переделки базового потока [7, 8].", st_td),
            ],
            [
                P("<b>5. Госперечни РФ и масштаб Шага 2</b>", st_td_b),
                P("Шаг 2 (SAMe): 14–16 млрд ₽/год в РФ (ЖНВЛП «Фарма-2030», п. 24 МСХ № 392).", st_td),
                P("<b>Оба шага в списках:</b> п. 32 Приказа МСХ № 392 («Витамин B₂, его соединения и производные») [12] + АФС для «Цитофлавина», «Декстралинка» и сенсоров глюкозы [15, 16].", st_td),
            ],
        ]
        c_tbl = Table(comp_rows, colWidths=[30 * mm, 71 * mm, 73 * mm], repeatRows=1)
        c_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
            ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
            ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3.8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3.8),
            ("TOPPADDING", (0, 0), (-1, -1), 3.4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3.4),
        ]))
        story.append(c_tbl)
        story.append(P("Таблица 1. Прямое сравнение двух продуктовых платформ DBTL по пяти базовым критериям проекта.", st_caption))

        # ==================== SECTION 2: EXISTING STRAINS, GENES VS WT YEAST & DBTL AUDIT ====================
        story.append(PageBreak())
        story.append(H1("s2", "2. Готовые штаммы-референсы, отличие от базовых дрожжей и аудит работ по DBTL/ML"))

        story.append(Paragraph(fmt_chem("2.1. Какие гены вставлены по сравнению с базовыми дрожжами и есть ли готовый референс Шага 1"), st_h2))
        story.append(
            P(
                "Ключевое практическое преимущество рибофлавина для быстрого запуска «Цикла 0» состоит в том, что <b>любые базовые лабораторные и промышленные дрожжи (<i>S. cerevisiae, Komagataella phaffii / Pichia pastoris, Yarrowia lipolytica, Candida famata</i>) изначально содержат в собственном геноме полный путь синтеза и витамина B₂, и коферментов ФМН/ФАД</b> — синтезировать и вставлять чужеродные гены из других царств не требуется [1, 17]:",
                st_body,
            )
        )
        wt_pts = [
            "<b>Отличие продуцента от диких дрожжей — снятие промоторного тормоза с собственных генов <i>RIB</i>:</b> в хромосомах базовых дрожжей уже лежат все 6 генов сборки рибофлавина (<i>RIB1</i> — ГТФ-циклогидролаза II, первый и главный лимитирующий фермент; <i>RIB7</i> — редуктаза; <i>RIB2</i> — дезаминаза; <i>RIB3</i> [у <i>C. famata</i> — <i>RIB6</i>] — синтаза ДГБФ, вводящая пентозофосфатный поток Ru5P; <i>RIB4</i> — люмазинсинтаза; <i>RIB5</i> — рибофлавинсинтаза) плюс оба гена Шага 2 (<i>FMN1</i> и <i>FAD1</i>). В диком типе на глюкозе <i>RIB1</i> жёстко зажат на уровне &lt;0,5 мг/л. Как показано в классическом протоколе <i>Marx, Mattanovich &amp; Sauer (2008)</i> на <i>K. phaffii (P. pastoris X-33 / GS115)</i> [17], простая постановка <b>одного собственного гена <i>PpRIB1</i> под сильный конститутивный промотор <i>pGAP</i></b> (амплифицируется обычным ПЦР с геномной ДНК самих дрожжей за 1 день) сразу окрашивает колонии и среду в ярко-жёлтый цвет, давая <b>3,9–6,9 мг/л в колбах и 50 мг/л в супернатанте реактора</b>; пара <i>pGAP–RIB1 + RIB3</i> даёт <b>10,1 мг/л</b>, а дерегуляция всех 6 нативных генов <i>RIB1–RIB7</i> — <b>20,2 мг/л в колбах и 175 мг/л в супернатанте биореактора</b> ещё до оптимизации предшественников.",
            "<b>Готовые штаммы в российской коллекции БРЦ ВКПМ и мировые дрожжевые референсы (Таблица 2):</b> помимо одношаговой сборки дрожжевого референса, в фонде БРЦ ВКПМ НИЦ «Курчатовский институт» (Москва) физически депонирована отечественная линейка штаммов <i>Bacillus subtilis</i> с истёкшими сроками патентной защиты РФ (включая открытый каталожный штамм <b>ВКПМ B-2694</b> по цене 3 000 руб. и высокопродуктивные штаммы <b>ВКПМ B-6797, B-8386, B-8387</b> с выходом 12,4–21 г/л [14]).",
        ]
        for b in wt_pts:
            story.append(P("• " + b, st_bullet))

        str_rows = [
            [
                P("Шасси и конкретные штаммы", st_th),
                P("Генетические отличия от базового дикого штамма", st_th),
                P("Стартовый и рекордный выход B₂", st_th),
                P("Доступность для Цикла 0 и патенты РФ", st_th),
            ],
            [
                P("<b>Модельные дрожжи:</b><br/><i>K. phaffii</i> X-33 / GS115 [17, 18]; <i>Y. lipolytica</i> [4, 19]; <i>S. cerevisiae</i> CEN.PK", st_td_b),
                P("<b>0 чужеродных генов:</b> ПЦР-вставка своих генов <i>RIB1</i> (1 ген), <i>RIB1+RIB3</i> (2 гена) или всех 6 <i>RIB1–7</i> под промоторами <i>pGAP / pTEF1</i> (либо <i>dCas9-VP64</i> к <i>RIB1</i>).", st_td),
                P("<b>4–10 мг/л</b> (1–2 гена); <b>20,2 мг/л</b> (колбы) и <b>175 мг/л</b> (реактор, 6 генов <i>K.p.</i>) — в 30–300× выше порога 445 нм.", st_td),
                P("Штаммы X-33, GS115, CEN.PK и <i>Y.l.</i> стандартны; сборка с геномной ДНК за 1–2 нед.; патентов РФ нет.", st_td),
            ],
            [
                P("<b>Флавиногенные дрожжи:</b><br/><i>C. famata</i> AF-4 → BRP → BRPI [1–3, 9]; <i>M. guilliermondii</i> DSM 11947", st_td_b),
                P("Доп. копии своего транскрипционного активатора <i>SEF1</i> + <i>RIB1, RIB6, RIB7</i> + снятие ретроингибирования пуринов (<i>PRS3, ADE4, IMH3</i>), <i>GND1</i> и экспортёр <i>RFE1</i>.", st_td),
                P("Не-ГМО <i>M.g.</i>: <b>34–318 мг/л</b>; инж. <i>C. famata</i> BRPI: <b>653–822 мг/л</b> (колбы) и <b>16,4 г/л</b> (7-л реактор).", st_td),
                P("В БРЦ ВКПМ и ВКМ — сотни изолятов <i>Pichia / Candida / Ogataea</i> (вкл. ВКПМ Y-201 с высокой активностью <i>FMN1</i>).", st_td),
            ],
            [
                P("<b>Бактерии БРЦ ВКПМ:</b><br/><i>B. subtilis</i> <b>ВКПМ B-2694</b>, <b>B-6797</b>, <b>B-8386</b>, <b>B-8387</b> [14]", st_td_b),
                P("Мутации <i>ribC</i> и оператора <i>ribO</i> (снятие репрессии ФМН-рибопереключателя) + плазмида <i>pMX45</i> с <i>rib</i>-опероном (+ оперон <i>B. amyloliquefaciens</i> в хромосоме).", st_td),
                P("B-2694: <b>3,5–4,5 г/л</b>; B-6797: <b>12,4 г/л</b> (42 ч); B-8386/87: <b>8,8–9,6 г/л</b> (колбы) и <b>18–21 г/л</b> (70 ч, реактор).", st_td),
                P("<b>Физически в Москве</b> (каталог ВКПМ: B-2694 — 3 000 ₽); патенты RU 2081906 и RU 2261273 <b>истекли</b>.", st_td),
            ],
        ]
        st_tbl = Table(str_rows, colWidths=[38 * mm, 54 * mm, 43 * mm, 39 * mm], repeatRows=1)
        st_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
            ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
            ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3.5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3.5),
            ("TOPPADDING", (0, 0), (-1, -1), 2.5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
        ]))
        story.append(st_tbl)
        story.append(P("Таблица 2. Существующие штаммы-продуценты рибофлавина и их точные генетические отличия от базовых штаммов дикого типа.", st_caption))

        story.append(Paragraph(fmt_chem("2.2. Использовали ли уже подход DBTL и ML для синтеза рибофлавина? (Аудит 2014–2026 гг.)"), st_h2))
        dbtl_audit = [
            "<b>Что уже делали на бактериях и в случайном мутагенезе (без замкнутого ML-цикла Learn на титр):</b> (1) в <i>B. subtilis</i> группа Long Liu построила общую омиксную модель клетки <i>EM_iBsu1209-ME</i> [20], а в июле 2026 г. (<i>Li et al.</i> [20]) применила стехиометрический расчёт <i>OptKnock/MOMA</i> (без ML) и слепой 96-луночный скрининг случайной библиотеки <i>CRISPRi</i> (<i>↓guaC + ↓ilvD</i>, 20,2 г/л) — при этом авторы экспериментально зафиксировали <b>сильный отрицательный (неаддитивный) эпистаз</b>: механическое объединение лучших одиночных нокаутов снижало титр B₂; (2) в <i>E. coli</i> и <i>C. glutamicum</i> применяли однократный биофизический расчёт силы RBS (<i>RBS Calculator</i> в один проход без дообучения [10]) и однораундовый скрининг 91 мишени <i>sRNA</i> (20,5 г/л) [10]; (3) в дрожжах <i>Y. lipolytica</i> (<i>Wagner et al., 2018</i> [4]) сортировали в микрокаплях (FACS) случайные УФ/EMS-мутанты без генотипирования и без стадии Learn.",
            "<b>Что в мире НЕ сделано (наша чистая научная новизна Q1):</b> (1) <b>в дрожжах (<i>S. cerevisiae, K. phaffii, Y. lipolytica, C. famata</i>) нет ни одной статьи с замкнутым циклом ML-DBTL</b> — все работы вплоть до 2026 г. шли ручным поочерёдным перебором по одному гену [1–3, 17–19]; (2) <b>ни в одном организме мира не применялся ML-DBTL к сквозному каскаду «B₂ → 5'-ФМН → ФАД»</b>, где возникает идеальный для активного обучения метаболический конфликт на развилке пуринов (раздел 3).",
        ]
        for b in dbtl_audit:
            story.append(P("• " + b, st_bullet))

        # ==================== SECTION 3 ====================
        story.append(H1("s3", "3. Прямой биохимический каскад (B2 → ФМН → ФАД) и развилка пуринов для моделирования"))
        story.append(
            P(
                "Синтез одной молекулы рибофлавина (C₁₇H₂₀N₄O₆) объединяет в клетке две центральные метаболические магистрали (Рисунок 1): одну молекулу <b>гуанозин-5'-трифосфата (GTP)</b> из пути синтеза пуринов <i>de novo</i> и две молекулы <b>D-рибулозо-5-фосфата (Ru5P)</b> из окислительной ветви пентозофосфатного пути (каскад генов <i>RIB1, RIB7, RIB2, RIB3, RIB4, RIB5</i> у дрожжей либо оперон <i>ribDEAHT</i> у бацилл) [1, 5]. На втором шаге полученный рибофлавин не сбрасывается в боковую ветвь, а непосредственно достраивается двумя последовательными АТФ-зависимыми реакциями: монофункциональная рибофлавинкиназа <b>FMN1</b> селективно фосфорилирует 5'-гидроксил рибитильного хвоста с образованием <b>ФМН</b>, после чего ФАД-синтаза (ФМН-аденилилтрансфераза <b>FAD1</b>) переносит аденилильный остаток со второй молекулы АТФ, замыкая динуклеотид <b>ФАД</b> [7, 9].",
                st_body,
            )
        )

        img_p = Image(str(fig_p), width=usable_w, height=usable_w * (5.0 / 11.0))
        story.append(
            KeepTogether([
                img_p,
                P("Рисунок 1. Структурные химические формулы прямого каскада «Рибофлавин (Шаг 1) → 5'-ФМН (Шаг 2А) → ФАД (Шаг 2Б)» и развилка предшественников пентозофосфатного и пуринового путей, оптимизируемая моделью DBTL.", st_caption),
            ])
        )

        story.append(Paragraph(fmt_chem("Преимущество дрожжевого шасси на Шаге 2 и парадокс узла ИМФ для модели Design–Learn"), st_h2))
        sec2_pts = [
            "<b>Разделение ферментов <i>FMN1</i> и <i>FAD1</i> у дрожжей против слитого <i>RibC/RibF</i> у бактерий:</b> у бактерий (<i>B. subtilis, E. coli</i>) обе реакции Шага 2 ведёт один бифункциональный белок <i>RibC / RibF</i>, который трудно остановить на стадии чистого ФМН и который катализирует обратный распад ФАД [7, 14]. У дрожжей (<i>K. phaffii, C. famata, Y. lipolytica, S. cerevisiae</i>) <b>FMN1</b> и <b>FAD1</b> — два независимых монофункциональных белка: усиление только <i>FMN1</i> даёт чистый 5'-ФМН, а добавление <i>FAD1</i> переключает поток на ФАД без эндотоксинов и угрозы бактериофагов [7, 9, 19].",
            "<b>Нелинейный парадокс развилки ИМФ (<i>IMP → GTP vs ATP</i>) как ядро статьи ML-DBTL:</b> и гуаниновый предшественник рибофлавина (<b>GTP</b>), и энергетический донор второй ступени (<b>2×ATP</b>) образуются из общего узла <b>инозинмонофосфата (IMP)</b>. На Шаге 1 (чистый B₂) выгодно блокировать адениловую ветвь <i>IMP → AMP</i> (делеция <i>ADE12</i> дала рост B₂ в 5,4 раза [18]), но при переходе к Шагу 2 (ФМН и особенно ФАД) клетке одновременно нужны <b>1 моль GTP + 2 моль ATP</b> (поэтому в работе <i>Hou et al., февр. 2026</i> для получения 2,31 г/л ФАД пришлось, наоборот, сверхэкспрессировать ортолог <i>purA/ADE12</i> [7]!). Нахождение точного баланса силы промоторов на узле <i>[ZWF1/GND1/RPE + ADE4 + IMD vs ADE12 + RIB1/RIB3 + FMN1/FAD1]</i> через связку GEM (<i>ecYeast8 / iKp1260</i>) и байесовского активного обучения решает задачу, недоступную ручному OFAT-перебору.",
        ]
        for b in sec2_pts:
            story.append(P("• " + b, st_bullet))

        # ==================== SECTION 4 ====================
        story.append(H1("s4", "4. Быстрая дифференциальная фотометрия без масс-спектрометра (Test) и план циклов"))
        test_pts = [
            "<b>Шаг 1 (рибофлавин) — измерение за 10 секунд без единого реактива:</b> рибофлавин активно секретируется клетками в прозрачную минеральную среду и имеет интенсивную полосу поглощения в видимой синей области при <b>λ = 445 нм</b> (молярный коэффициент ε₄₄₅ = 12 500 М⁻¹·см⁻¹; ни аминокислоты, ни сахара, ни пурины при 445 нм не поглощают). Центрифугированный супернатант просто разбавляют водой в 10–50 раз в обычной кювете или 96-луночном планшете и сразу вычисляют концентрацию: <b>C (мг/л) = 30,1 × A₄₄₅ × K_разв</b> (при микроколичествах используется собственная ярко-зелёная флуоресценция Ex 450 нм / Em 525 нм).",
            "<b>Шаг 2 (раздельное определение B₂, ФМН и ФАД без ВЭЖХ-МС — метод Бесси–Лоури [11]):</b> все три флавина поглощают при 445 нм одинаково, но резко различаются по квантовому выходу флуоресценции и гидрофобности. Из-за внутримолекулярного стэкинга аденинового и изоаллоксазинового колец <b>молярная флуоресценция ФАД при pH 7,0 потушена почти в 10 раз (~9–10% от сигнала ФМН)</b>. Это позволяет за 20 минут в обычном планшетном флуориметре (Ex 450 / Em 525 нм) разделить продукты без хроматографа: (1) замер свечения нативной пробы при pH 7,0 (F₀); (2) прогрев аликвоты в 0,1 М HCl (100 °C, 15 мин — количественный гидролиз пирофосфатной связи <i>ФАД → ФМН</i>) и повторный замер при pH 7,0 (F₁) — скачок ΔF = F₁ − F₀ строго пропорционален концентрации <b>ФАД</b>; (3) однократная экстракция нейтрального рибофлавина бензиловым спиртом отделяет его от заряженных в воде фосфатов <b>ФМН + ФАД</b>. Лучшие 10–15% клонов верифицируются быстрым 5-минутным изократическим ВЭЖХ-УФ (колонка C18, метанол / ацетатный буфер, детектор 445 нм).",
        ]
        for b in test_pts:
            story.append(P("• " + b, st_bullet))

        rm_rows = [
            [
                P("Цикл конвейера", st_th),
                P("Срок", st_th),
                P("Экспериментальные работы (Build + Test)", st_th),
                P("Целевой выход и обучение модели (Learn)", st_th),
            ],
            [
                P("<b>Цикл 0–1 (Шаг 1)</b><br/>Узел B₂ (<i>GTP + Ru5P</i>)", st_td_b),
                P("Нед. 1–8", st_td),
                P("Сборка библиотеки промоторов перед узлами ПФП (<i>ZWF1, GND1, RPE</i>), входа в пурины (<i>ADE4/purF, IMD</i>), каскада <i>RIB1–RIB7</i> и экспортёра <i>RFE1</i>; прямой скрининг при 445 нм.", st_td),
                P("Массив >350 точек «вектор экспрессии → мг/л·ч B₂»; обучение GP / XGBoost модели баланса <i>Ru5P : GTP</i> без просадки роста.", st_td),
            ],
            [
                P("<b>Цикл 2–3 (Шаг 2)</b><br/>Переключение на ФМН и ФАД", st_td_b),
                P("Нед. 9–18", st_td),
                P("Интеграция <i>ScFMN1</i> (линия <b>2А — чистый 5'-ФМН</b>) и связки <i>ScFMN1 + CgFAD1</i> с усилением <i>ADE12/purA</i> и делецией 5'-нуклеотидаз (линия <b>2Б — ФАД</b>); дифференциальный флуоресцентный скрининг 525 нм.", st_td),
                P("Выход на >2,5–4,0 г/л региоселективного 5'-ФМН (>96%) и >1,5–2,5 г/л ФАД <i>de novo</i>; заявка на патент РФ, статья Q1 и образцы фармсубстанции.", st_td),
            ],
        ]
        rm_tbl = Table(rm_rows, colWidths=[31 * mm, 17 * mm, 66 * mm, 60 * mm], repeatRows=1)
        rm_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
            ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
            ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
            ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 3.8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3.8),
            ("TOPPADDING", (0, 0), (-1, -1), 3.6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3.6),
        ]))
        story.append(Spacer(1, 1.2 * mm))
        story.append(rm_tbl)
        story.append(P("Таблица 3. Компактный план отработки цикла DBTL с прямым переходом от рибофлавина к коферментам ФМН и ФАД.", st_caption))

        # ==================== SECTION 5 ====================
        story.append(H1("s5", "5. Рынок РФ (п. 32 МСХ № 392, реестр ГРЛС), преимущество над химией и экономика"))
        econ_pts = [
            "<b>Шаг 1 в РФ — гарантированный спрос в сотни тонн и прямая строка господдержки (с честной оговоркой по масштабу):</b> по данным аналитики <i>FEEDLOT</i> (сентябрь 2026 г.) [13], выпуск комбикормов в РФ достиг 45,1 млн т/год (премиксов — 586,6 тыс. т/год), при этом внутреннее производство витамина B₂ в РФ равно <b>0%</b>, а годовой импорт кормового (80%), пищевого (краситель E101) и фармацевтического рибофлавина (ГФ РФ XV <b>ФС.2.1.0559</b> прямо определяет субстанцию как продукт ферментации; в ГРЛС зарегистрированы субстанции <i>Long Sheng</i> КНР, <i>BASF</i> и <i>DSM</i> [15]) составляет <b>300–500 тонн/год</b>. Позиция <b>№ 32 Приказа Минсельхоза России от 06.06.2025 № 392</b> (<i>«Витамин B₂ (рибофлавин), его соединения и производные»</i>) открывает доступ к возмещению 20% капитальных затрат с 2027 г. [12]. Однако важно трезво учитывать цену: кормовой B₂ 80% в мире стоит <b>$12–22/кг</b> (пищевой/фарм 98% — <b>$35–60/кг</b>; Рисунок 2А), поскольку китайские заводы (<i>Hubei Guangji</i>) ведут процесс в ферментёрах >100 м³ при титрах 26–34 г/л [6]. Поэтому на опытном реакторе 0,5–2 м³ продавать кормовой B₂ убыточно — экономика малотоннажной линии сходится именно на продуктах Шага 2.",
            "<b>Шаг 2А в РФ — региоселективный био-5'-ФМН против грязного фосфорилирования POCl₃:</b> водорастворимая фармсубстанция <i>рибофлавина фосфат натрия (ФМН)</i> остро востребована в России в трёх сегментах: (1) в составе крупнейшего отечественного рецептурного нейрометаболического препарата <b>«Цитофлавин»</b> (ООО «НТФФ ПОЛИСАН», РУ ЛП-№(000973)-(РГ-RU)) [15]; (2) в инъекционном препарате <b>«Рибофлавин-мононуклеотид»</b> (ОАО «Фармстандарт-УфаВИТА», ЛП-№(008441)-(РГ-RU)) [15]; (3) в офтальмологических фотосенсибилизаторах для УФ-кросслинкинга роговицы при кератоконусе (российские изделия <b>«Декстралинк», «Риболинк», «Гиполинк»</b> Уфимского НИИ глазных болезней [16]). В ГРЛС по субстанции ФМН доминируют поставщики из Индии (<i>Harman Finochem, Danashmand</i>) и Швейцарии (<i>DSM</i>) при малотоннажном химическом синтезе у <i>ООО «Полисинтез»</i> [15]. При этом традиционное химическое фосфорилирование рибофлавина хлорокисью фосфора (POCl₃) неселективно и даёт <b>20–28% изомерных примесей</b> (3'- и 4'-фосфаты, дифосфаты; Рисунок 2Б), которые ингибируют клеточные флавопротеины и требуют дорогой очистки, тогда как ферментативное фосфорилирование киназой <i>FMN1</i> в воде даёт <b>>95–98% чистого природного 5'-изомера</b> по цене АФС <b>$250–600/кг</b>.",
            "<b>Шаг 2Б в РФ и мире — высокомаржинальный кофермент ФАД ($1 200–3 000+/кг):</b> флавинадениндинуклеотид невозможно синтезировать химически в промышленном масштабе. В медицине он применяется в парентеральных и глазных формах при дистрофиях сетчатки/роговицы и гепатопатиях, а в биотехнологической индустрии является <b>ключевым расходным кофактором современных тест-полосок глюкометров</b> на базе ФАД-зависимой глюкозодегидрогеназы (<i>FAD-GDH</i>, не чувствительной к кислороду крови и мальтозе), что формирует устойчивый высокомаржинальный спрос для компактного ферментационного участка.",
        ]
        for b in econ_pts:
            story.append(P("• " + b, st_bullet))

        img_e = Image(str(fig_e), width=usable_w, height=usable_w * (3.65 / 10.6))
        story.append(
            KeepTogether([
                Spacer(1, 1.2 * mm),
                img_e,
                P("Рисунок 2. Экономика и качество каскада Варианта Б: (А) кратный рост оптовой цены B2B при переходе от кормового B₂ к региоселективному 5'-ФМН и коферменту ФАД; (Б) сопоставление изомерной чистоты химического (POCl₃) и биотехнологического (FMN1) синтеза ФМН.", st_caption),
            ])
        )

        # ==================== SECTION 6: REFERENCES ====================
        story.append(H1("s6", "6. Список литературы и нормативных источников"))
        refs = [
            "[1] Dmytruk K. V., Ruchała J., Fedorovych D. V., Ostapiv R. D., Sibirny A. A. Modulation of the purine pathway for riboflavin production in flavinogenic recombinant strain of the yeast <i>Candida famata</i> (до 16,4 г/л). <i>Biotechnology Journal</i>. 2020;15(7):e1900468. doi:10.1002/biot.201900468.",
            "[2] Ruchała J. et al. Cheese whey supports high riboflavin synthesis by the engineered strains of the flavinogenic yeast <i>Candida famata</i>. <i>Microbial Cell Factories</i>. 2022;21:161. doi:10.1186/s12934-022-01888-0; Dzanaeva L. et al. <i>FEMS Yeast Research</i>. 2024;24:foae020.",
            "[3] Коэкспрессия генов <i>RFE1</i> (экспортёр рибофлавина), <i>GND1</i>, <i>RIB1</i> и <i>RIB6</i> в дрожжах <i>Candida famata</i>. <i>J. Fungi</i>. 2022;8(4):141 &amp; <i>Microbiology</i>. 2025; не-ГМО скрининг <i>Meyerozyma guilliermondii</i> DSM 11947 (317,6 мг/л). <i>Organic Agriculture</i>. 2023;13:399–410.",
            "[4] Wagner J. M., Liu L., Yuan S.-F., Venkataraman M. V., Abate A. R., Alper H. S. A comparative analysis of single cell and droplet-based FACS for improving production phenotypes: riboflavin overproduction in <i>Yarrowia lipolytica</i>. <i>Metabolic Engineering</i>. 2018;47:346–356. doi:10.1016/j.ymben.2018.04.015.",
            "[5] Liu C., Xia M., Fang H., Xu F., Wang S., Zhang D. De novo engineering riboflavin production in <i>Bacillus subtilis</i> by overexpressing the downstream genes in the purine biosynthesis pathway (<i>guaB, guaA, gmk, ndk, ribA</i>). <i>Microbial Cell Factories</i>. 2024;23:157. doi:10.1186/s12934-024-02426-w.",
            "[6] Wang S., Zhu Q., Liu C., Dong H., Xia M., Jin Z., Zhang D. Engineering riboflavin-overproducing <i>Bacillus subtilis</i> via pathway gene overexpression (30–34 г/л и анализ ростовых дефектов / нестабильности). <i>Synthetic and Systems Biotechnology</i>. 2025;11:91–101. doi:10.1016/j.synbio.2025.09.005.",
            "[7] Hou J., Hu W., Ding M., Liu C., Zou S. Metabolic engineering of <i>Escherichia coli</i> for efficient de novo production of redox cofactors flavin mononucleotide (3.46 g/L FMN) and flavin adenine dinucleotide (2.31 g/L FAD via <i>ScFMN1 + CgFMNAT</i> + <i>↑purA</i>). <i>J. Agric. Food Chem.</i> 2026 (Feb 2026). doi:10.1021/acs.jafc.5c07311.",
            "[8] Wu Y., Sun Y., Miao H., Xia X., Gao L. Engineering of riboflavin transporter (<i>RibM</i>) and enzymes in <i>Bacillus subtilis</i> facilitates high-efficiency FAD synthesis from riboflavin (конверсия 92%). <i>ACS Synthetic Biology</i>. 2026 (May 2026). doi:10.1021/acssynbio.6c00153.",
            "[9] Fedorovych D. V., Dmytruk K. V., Sibirny A. A. Recent advances in construction of the efficient producers of riboflavin and flavin nucleotides (FMN, FAD) in the yeast <i>Candida famata</i>. <i>Methods in Molecular Biology</i>. 2021;2280:19–33; Yatsyshyn V. Y. et al. <i>Metabolic Engineering</i>. 2009;11:163–167.",
            "[10] Lin Z. et al. <i>Microb. Cell Fact.</i> 2014;13:104 (RBS-модуляция <i>ribF</i> в <i>E. coli</i>, 2,7 г/л); Pérez-García F. et al. <i>Biofabrication</i>. 2024 (8,3 г/л в <i>C. glutamicum</i>); sRNA-аттенюация штамма <i>E. coli</i> SN49M (20,5 г/л, 2026); Liu S., Sun D. <i>World J. Microbiol. Biotechnol.</i> 2025;41.",
            "[11] Bessey O. A., Lowry O. H., Love R. H. The fluorometric measurement of the nucleotides of riboflavin and their concentration in tissues (эффект 10-кратного тушения флуоресценции ФАД относительно ФМН). <i>J. Biol. Chem.</i> 1949;180(2):755–769.",
            "[12] Приказ Минсельхоза России от 06.06.2025 № 392 (рег. Минюста РФ № 83535 от 12.09.2025), позиция 32: «Витамин B₂ (рибофлавин), его соединения и производные»; Постановление Правительства РФ о возмещении до 20% CAPEX с 01.01.2027.",
            "[13] Аналитический обзор агентства FEEDLOT «Инвестиции в производство кормовых добавок в России» (03.09.2026: выпуск 45,1 млн т комбикормов и 586,6 тыс. т премиксов при нулевой самообеспеченности витаминами группы B); ФС.2.1.0559 «Рибофлавин» ГФ РФ XV изд. (Приказ Минздрава РФ № 377).",
            "[14] Отечественные штаммы БРЦ ВКПМ НИЦ «Курчатовский институт» (ГосНИИгенетика): <i>B. subtilis</i> ВКПМ B-2694 (304/pMX45, открытый каталог ВКПМ), ВКПМ B-6797 (RU 2081906, 12,4 г/л, Козлов Ю. И., Перумов Д. А. и др.), ВКПМ B-8386 и B-8387 (RU 2261273 C2, до 21 г/л, Миронов А. С., Перумов Д. А., Дебабов В. Г. и др.); скрининг активности <i>FMN1</i> у <i>O. polymorpha</i> ВКПМ Y-201 (<i>Studia Biologica</i>, 2009).",
            "[15] ГРЛС Минздрава России (выгрузка ноябрь 2026 г.): субстанции рибофлавина и рибофлавин-5'-фосфата натрия (ФС-000117 Long Sheng, П N013263/01 DSM, ЛСР-007577/08 Harman Finochem, ЛС-002018 Danashmand, ЛСР-002295/07 ООО «Полисинтез»); готовые ЛП «Цитофлавин» (НТФФ «ПОЛИСАН», ЛП-№(000973)-(РГ-RU)) и «Рибофлавин-мононуклеотид» (ОАО «Фармстандарт-УфаВИТА», ЛП-№(008441)-(РГ-RU)).",
            "[16] Бикбов М. М., Халимов А. Р. и др. УФ-А-кросслинкинг роговицы с 0,1% растворами рибофлавина мононуклеотида («Декстралинк» РУ ФСР 2010/09071, «Риболинк», «Гиполинк», Уфимский НИИ глазных болезней). <i>Офтальмология</i>. 2024;21(4):802–808.",
            "[17] Marx H., Mattanovich D., Sauer M. Overexpression of the riboflavin biosynthetic pathway in <i>Pichia pastoris</i> (<i>Komagataella phaffii</i>: от одного гена <i>RIB1</i> до всех 6 нативных генов <i>RIB1–RIB7</i> под <i>pGAP</i>, 20,2–175 мг/л в супернатанте). <i>Microbial Cell Factories</i>. 2008;7:23. doi:10.1186/1475-2859-7-23; Santos M. A., Revuelta J. L. <i>J. Biol. Chem.</i> 1995;270:437–444 (гены <i>RIB1–7</i> в <i>S. cerevisiae</i>).",
            "[18] Baumschabl M., Prielhofer R., Mattanovich D., Steiger M. G. Fine-tuning of transcription in <i>Pichia pastoris</i> using dCas9 and RNA scaffolds (активация <i>RIB1</i>). <i>ACS Synthetic Biology</i>. 2020;9(12):3202–3209; Ruchała J., Sibirny A. et al. Regulation of riboflavin biosynthesis in microorganisms (вкл. эффект <i>ΔADE12</i>). <i>Int. J. Mol. Sci.</i> 2025;26(13):6243; Jiménez-Nava R. A. et al. <i>Pharmaceuticals</i>. 2026;19(3):389.",
            "[19] Elsharawy H., He Q., Zeng W., Zhou J. Systematic enzyme and cofactor engineering (<i>YlRIB1, RIB2, RIB3, FMN1, FAD1</i>) in <i>Yarrowia lipolytica</i>. 2026 (Jul 2026, PMC13400399); Wu Y. et al. (<i>MCH5, RIB1, FLX1</i> в <i>Y. lipolytica</i>). 2025 (PMC12767179).",
            "[20] Bi X., Liu L. et al. A multi-omics, machine learning-aware, genome-wide metabolic model of <i>Bacillus subtilis</i> (<i>EM_iBsu1209-ME</i>). <i>Advanced Science</i>. 2024;11(42):2408705; Li Z., Bi X., Liu L., Chen J. et al. Systems metabolic engineering of <i>B. subtilis</i> for high-level riboflavin production via dynamic regulation and CRISPRi screening (неаддитивность нокаутов и <i>↓guaC + ↓ilvD</i>, 20,2 г/л). <i>ENGINEERING Foods</i>. 2026;1(1):29–45. doi:10.2738/ENGF.2026.0003; Xu X. et al. <i>Nucleic Acids Res.</i> 2025 (<i>dhfCas13X</i>).",
        ]
        for r in refs:
            story.append(P(r, st_ref))

        return story

    def render_pass(out_file: Path, toc_in: dict) -> tuple[int, dict]:
        rec = {}
        doc = BriefDocTemplate(
            str(out_file),
            rec,
            pagesize=A4,
            leftMargin=left_m,
            rightMargin=right_m,
            topMargin=top_m,
            bottomMargin=bot_m,
            title="Вариант Б: Витамин B2 (рибофлавин) -> ФМН и ФАД для платформы DBTL",
            author="Проект DBTL",
        )
        frame = Frame(left_m, bot_m, usable_w, page_h - top_m - bot_m, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=draw_hf)])
        doc.build(build_story(toc_in))
        return doc.page, rec

    p1, exact_map = render_pass(OUT_PDF, {})
    p2, verify_map = render_pass(OUT_PDF, exact_map)
    assert p1 == p2 and exact_map == verify_map, f"TOC pass mismatch: {p1} vs {p2}"

    pdf_doc = fitz.open(str(OUT_PDF))
    bad_pages = [p.number + 1 for p in pdf_doc if chr(65533) in p.get_text()]
    assert not bad_pages, f"Missing glyphs on pages: {bad_pages}"
    links_p1 = len(pdf_doc[0].get_links())
    for old_png in CTX_DIR.glob("page_*.png"):
        old_png.unlink()
    for idx, page in enumerate(pdf_doc, start=1):
        pix = page.get_pixmap(dpi=160, alpha=False)
        pix.save(str(CTX_DIR / f"page_{idx:02d}.png"))
    pdf_doc.close()

    print(f"SUCCESS: {OUT_PDF} ({p2} pages, {links_p1} clickable TOC links on page 1, {OUT_PDF.stat().st_size} bytes)")


if __name__ == "__main__":
    build_pdf()
