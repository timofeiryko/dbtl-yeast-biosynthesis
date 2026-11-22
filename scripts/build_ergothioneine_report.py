# -*- coding: utf-8 -*-
"""
Build publication-grade A4 analytical report PDF on selecting L-Ergothioneine (Step 1)
and the SAM metabolic hub -> (S,S)-Ademetionine / SAMe (Step 2) for the yeast DBTL platform.

Generates:
  - assets/fig1_candidates.png
  - assets/fig2_pathway_structures.png (composited with RDKit 2D chemical structures)
  - assets/fig3_dbtl_loop.png
  - assets/fig4_markets.png
  - ergothioneine-dbtl-report.pdf
"""

from pathlib import Path

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
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
CTX_DIR = ROOT_DIR / ".context" / "ergothioneine-report"
ASSETS_DIR.mkdir(parents=True, exist_ok=True)
CTX_DIR.mkdir(parents=True, exist_ok=True)

OUT_PDF = ROOT_DIR / "ergothioneine-dbtl-report.pdf"


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

# Restrained academic editorial palette
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
    """Heading paragraph that records its exact page number for the clickable TOC."""

    def __init__(self, text: str, style: ParagraphStyle, key: str, level: int, toc_title: str):
        anchored = f'<a name="{key}"/>' + fmt_chem(text)
        super().__init__(anchored, style)
        self.toc_key = key
        self.toc_level = level
        self.toc_title = toc_title


def crop_structure_png(png_path: Path, pad: int = 12) -> np.ndarray:
    """Auto-crop whitespace around an RDKit 2D structure PNG for crisp high-fill placement."""
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


def make_chart_alternatives(out_path: Path) -> None:
    """Figure 1: Quantitative comparison of candidate molecules."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.35), dpi=300)
    fig.patch.set_facecolor("white")

    items = [
        ("Таурин", 2.3, 5.5, "#8A969E", (6, 16)),
        ("Пуллулан (E1204)", 23.0, 5.2, "#8A969E", (-68, -16)),
        ("Аннато (E160b)", 32.0, 4.5, "#8A969E", (-46, 13)),
        ("2'-FL (ОГМ)", 35.0, 14.5, "#2F5573", (8, 4)),
        ("Нарингенин", 75.0, 7.0, "#8A969E", (-52, 16)),
        ("Ресвератрол (98%)", 105.0, 9.2, "#8A969E", (-26, 15)),
        ("п-Кумаровая к-та", 160.0, 6.2, "#8A969E", (-56, -17)),
        ("Ликопин (10%)", 220.0, 4.4, "#B84A39", (-44, -22)),
        ("Биотин B7", 350.0, 4.8, "#8A969E", (12, -13)),
        ("Шафран (10% кроцин)", 380.0, 8.4, "#2F5573", (10, -2)),
        ("L-Эрготионеин 99% (КНР)", 490.0, 18.0, "#1B5E4A", (-122, 4)),
        ("(S,S)-SAMe (фарм АФС)", 680.0, 12.2, "#1B5E4A", (-106, 11)),
    ]

    for label, price, cagr, col, offset in items:
        is_target = "Эрготионеин" in label or "SAMe" in label
        sz = 82 if is_target else 44
        ax1.scatter(price, cagr, s=sz, color=col, alpha=0.92, zorder=4, edgecolors="white", linewidth=0.8)
        ax1.annotate(
            label,
            (price, cagr),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.3 if not is_target else 7.9,
            fontweight="bold" if is_target else "normal",
            color="#1B5E4A" if is_target else "#2E373D",
            arrowprops=dict(arrowstyle="-", color="#B8C2C8", lw=0.55, shrinkA=0, shrinkB=3),
        )

    ax1.set_xscale("log")
    ax1.set_xlim(1.3, 3500)
    ax1.set_ylim(0.5, 22.5)
    ax1.axhspan(11.0, 22.5, color="#EAF2EF", alpha=0.45, zorder=1)
    ax1.grid(True, which="major", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax1.set_xlabel("Оптовая цена B2B субстанции ($/кг, лог. шкала)", fontsize=8.4, color="#2E373D")
    ax1.set_ylabel("Реалистичный темп роста рынка, CAGR (%)", fontsize=8.4, color="#2E373D")
    ax1.set_title("А. Оптовая цена субстанции и реалистичный CAGR", fontsize=9.4, fontweight="bold", color="#181D20", pad=8)
    ax1.tick_params(labelsize=8, colors="#4A555E")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    names_b = [
        "Биксин / аннато (3Б)",
        "Кальция фолинат",
        "Таурин (S.c.)",
        "Кроцетин (S.c., 3А)",
        "Биотин B7 (дрожжи)",
        "L-Эрготионеин (S.c.)",
        "L-Эрготионеин (Y.l.)",
        "(S,S)-SAMe (дрожжи*)",
        "Ликопин (Y.l.)",
        "Ресвератрол (Y.l.)",
        "2'-Фукозиллактоза (S.c.)",
    ]
    titers = [0.0005, 0.0004, 0.038, 0.107, 0.002, 3.05, 8.03, 17.1, 17.6, 22.5, 26.0]
    bar_cols = [
        "#B84A39", "#B84A39", "#B84A39", "#2F5573", "#B84A39",
        "#1B5E4A", "#1B5E4A", "#1B5E4A", "#8A969E", "#8A969E", "#8A969E",
    ]
    notes_b = [
        "0 г/л (путь не собран)",
        "0,0004 г/л",
        "0,038 г/л",
        "0,11 г/л (0 ML-работ)",
        "~0,002 г/л",
        "3,05 г/л (Ma 2026, 0 ML)",
        "8,03 г/л (0 ML-работ)",
        "17,1 г/л (*с подачей Met)",
        "17,6 г/л (десятки ML-работ)",
        "22,5 г/л (ML опубликован)",
        "26,0 г/л (патенты DSM)",
    ]

    y_pos = np.arange(len(names_b))
    ax2.barh(y_pos, titers, height=0.62, color=bar_cols, alpha=0.9, zorder=3)
    ax2.set_xscale("log")
    ax2.set_xlim(0.0002, 520)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(names_b, fontsize=7.8, color="#2E373D")
    ax2.axvline(1.0, color="#5E6A73", linestyle=":", linewidth=0.9, zorder=2)
    ax2.text(1.15, 0.45, "Ориентир пилота\n~1–3 г/л", fontsize=7.0, color="#5E6A73", va="bottom")

    for idx, (val, note) in enumerate(zip(titers, notes_b)):
        is_egt = "Эрготионеин" in names_b[idx] or "SAMe" in names_b[idx]
        ax2.text(
            val * 1.22,
            idx,
            note,
            va="center",
            fontsize=7.0 if not is_egt else 7.4,
            fontweight="bold" if is_egt else "normal",
            color="#1B5E4A" if is_egt else "#2E373D",
        )

    ax2.grid(True, axis="x", which="major", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax2.set_xlabel("Максимальный опубликованный титр в дрожжах (г/л, лог.)", fontsize=8.4, color="#2E373D")
    ax2.set_title("Б. Достигнутые титры в дрожжах и статус ML-DBTL", fontsize=9.4, fontweight="bold", color="#181D20", pad=8)
    ax2.tick_params(axis="x", labelsize=8, colors="#4A555E")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout(pad=1.2)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_pathway_structures_figure(out_path: Path) -> None:
    """
    Figure 2: Comprehensive biochemical & 2D chemical-structure diagram showing:
      - Native yeast synthesis of precursors from glucose & sulfate
      - Step 1 enzymatic reactions (NcEgt1, CpEgt2) with RDKit 2D structures & formulas
      - Explicit transition from Step 1 engineered strain to Step 2 products (SAMe, Spermidine, Glutathione)
    """
    fig = plt.figure(figsize=(11.2, 7.6), dpi=300)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 68)
    ax.axis("off")

    # Outer subtle frame
    outer = FancyBboxPatch(
        (0.6, 0.6), 98.8, 66.8,
        boxstyle="round,pad=0.3,rounding_size=1.2",
        facecolor="#FAFCFD", edgecolor="#CFD8DE", linewidth=1.0, zorder=0
    )
    ax.add_patch(outer)

    # Top banner: Native yeast precursor synthesis + Step 1 heterologous pathway
    ax.text(
        2.0, 65.0,
        "ШАГ 1 ПЛАТФОРМЫ: биосинтез L-эрготионеина в дрожжах de novo из глюкозы и сульфата (2 вставляемых гена: NcEgt1, CpEgt2)",
        fontsize=9.6, fontweight="bold", color="#181D20", va="center"
    )

    # Native Yeast Precursor Box (top feeder bar)
    feeder_bg = FancyBboxPatch(
        (2.0, 55.6), 96.0, 7.6,
        boxstyle="round,pad=0.2,rounding_size=0.8",
        facecolor="#EEF4F8", edgecolor="#9EB4C4", linewidth=0.9, zorder=1
    )
    ax.add_patch(feeder_bg)
    ax.text(
        3.2, 61.5,
        "Эндогенная база дрожжей (S. cerevisiae и Y. lipolytica сами синтезируют все 3 предшественника на минеральной среде):",
        fontsize=8.1, fontweight="bold", color="#2F5573", va="center"
    )
    ax.text(
        3.2, 59.1,
        "[1] Глюкоза → пентозофосфатный путь → PRPP → L-Гистидин (гены HIS1–HIS7)        "
        "[2] Глюкоза + NH₄⁺ → L-Метионин + АТФ → кофактор SAM (SAM2; рецикл SAH1, MET6)",
        fontsize=7.25, color="#2E373D", va="center"
    )
    ax.text(
        3.2, 56.9,
        "[3] Неорганический сульфат SO₄²⁻ + 2 АТФ + 4 NADPH → восстановление серы (MET3/14/16/10) и транссульфурирование (CYS3/CYS4) → L-Цистеин",
        fontsize=7.25, color="#2E373D", va="center"
    )

    def place_mol_card(x, y, w, h, img_key, title, formula_str, subtitle, border="#B8C4CC", bg="white", title_col="#181D20", lw=1.0, sub_lines=1):
        card = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.2,rounding_size=0.9",
            facecolor=bg, edgecolor=border, linewidth=lw, zorder=2
        )
        ax.add_patch(card)
        ax.text(x + w / 2, y + h - 1.6, title, fontsize=8.1, fontweight="bold", color=title_col, ha="center", va="center", zorder=4)
        ax.text(x + w / 2, y + h - 3.4, formula_str, fontsize=6.85, color="#5E6A73", ha="center", va="center", zorder=4)

        arr = crop_structure_png(STRUCT_DIR / f"{img_key}.png", pad=18)
        ih, iw = arr.shape[:2]
        bot_pad = 3.9 if sub_lines == 2 else 2.6
        box_w = w - 1.8
        box_h = h - 4.5 - bot_pad
        aspect_img = iw / ih
        if box_w / box_h > aspect_img:
            draw_h = box_h
            draw_w = draw_h * aspect_img
        else:
            draw_w = box_w
            draw_h = draw_w / aspect_img
        cx = x + w / 2
        cy = y + bot_pad + box_h / 2
        ax.imshow(arr, extent=[cx - draw_w / 2, cx + draw_w / 2, cy - draw_h / 2, cy + draw_h / 2], zorder=3, interpolation="lanczos")
        ax.text(
            x + w / 2,
            y + (2.05 if sub_lines == 2 else 1.35),
            subtitle,
            fontsize=6.55 if sub_lines == 2 else 6.75,
            color="#2E373D",
            ha="center",
            va="center",
            linespacing=1.22,
            zorder=4,
        )

    # Row 1: Four Step-1 chemical structure cards with wider reaction gaps (8.2 units)
    y1, h1 = 32.0, 20.6
    place_mol_card(
        2.0, y1, 17.2, h1,
        "histidine", "L-Гистидин", "C₆H₉N₃O₂ · 155,2 г/моль",
        "Нативный предшественник", border="#9EB4C4"
    )
    place_mol_card(
        27.4, y1, 17.6, h1,
        "hercynine", "Герцинин", "C₉H₁₅N₃O₂ · 197,2 г/моль",
        "Триметил-L-гистидин (бетаин)", border="#1B5E4A", bg="#FBFEFD"
    )
    place_mol_card(
        53.2, y1, 20.4, h1,
        "cysher", "Сульфоксид Cys-HER", "C₁₂H₂₀N₄O₅S · 332,4 г/моль",
        "C2-тиоэфирный сульфоксид", border="#1B5E4A", bg="#FBFEFD"
    )
    place_mol_card(
        81.4, y1, 16.6, h1,
        "ergothioneine", "L-Эрготионеин", "C₉H₁₅N₃O₂S · 229,3 г/моль",
        "Продукт Шага 1 (тион C=S)", border="#1B5E4A", bg="#EAF2EF", title_col="#1B5E4A", lw=1.6
    )

    def draw_rxn_arrow(x0, x1, y, top_lbl, bot_lbl, col="#1B5E4A"):
        arr = FancyArrowPatch(
            (x0, y), (x1, y),
            arrowstyle="-|>,head_length=5,head_width=3",
            color=col, linewidth=1.5, zorder=5
        )
        ax.add_patch(arr)
        xm = (x0 + x1) / 2
        ax.text(xm, y + 1.4, top_lbl, fontsize=7.25, fontweight="bold", color=col, ha="center", va="bottom")
        ax.text(xm, y - 1.3, bot_lbl, fontsize=6.45, color="#2E373D", ha="center", va="top", linespacing=1.18)

    # Reaction arrows between the 4 cards (gaps: 19.2..27.4, 45.0..53.2, 73.6..81.4)
    draw_rxn_arrow(19.4, 27.2, 41.5, "NcEgt1 (MT)", "+3 SAM\n→ 3 SAH")
    draw_rxn_arrow(45.2, 53.0, 41.5, "NcEgt1 (SO)", "+L-Цистеин\n+O₂ (Fe²⁺)")
    draw_rxn_arrow(73.8, 81.2, 41.5, "CpEgt2", "PLP (вит. B₆)\n−пируват, −NH₃")

    # Clean vertical feeder arrows stopping well above reaction text
    fa1 = FancyArrowPatch((10.6, 55.5), (10.6, 52.8), arrowstyle="-|>,head_length=4,head_width=2.4", color="#5E6A73", linewidth=0.95, linestyle="--", zorder=4)
    fa2 = FancyArrowPatch((23.3, 55.5), (23.3, 46.8), arrowstyle="-|>,head_length=4,head_width=2.4", color="#5E6A73", linewidth=0.95, linestyle="--", zorder=4)
    fa3 = FancyArrowPatch((49.1, 55.5), (49.1, 46.8), arrowstyle="-|>,head_length=4,head_width=2.4", color="#5E6A73", linewidth=0.95, linestyle="--", zorder=4)
    for fa in (fa1, fa2, fa3):
        ax.add_patch(fa)

    # Divider & Transition Banner: How Step 1 relates to Step 2 (with honest biochemical reconfiguration)
    trans_bg = FancyBboxPatch(
        (2.0, 24.4), 96.0, 5.4,
        boxstyle="round,pad=0.2,rounding_size=0.8",
        facecolor="#EAF2EF", edgecolor="#1B5E4A", linewidth=1.0, zorder=2
    )
    ax.add_patch(trans_bg)
    ax.text(
        50.0, 28.0,
        "ПЕРЕХОД К ШАГУ 2: переносится конвейер DBTL и энергетический блок [↑SAM2 + регенерация АТФ], но потоки HIS1/CYS4 реконфигурируются",
        fontsize=7.65, fontweight="bold", color="#1B5E4A", ha="center", va="center"
    )
    ax.text(
        50.0, 25.8,
        "Для SAMe (2А): без Egt1/2, подавление слива ↓CYS4, подпитка дешёвым L-метионином и депо в вакуоли (альтернатива — прямой старт с 2А):",
        fontsize=7.15, color="#181D20", ha="center", va="center"
    )

    a_mid = FancyArrowPatch((89.7, 31.8), (89.7, 29.9), arrowstyle="-|>,head_length=4.5,head_width=2.8", color="#1B5E4A", linewidth=1.4, zorder=5)
    ax.add_patch(a_mid)

    # Row 2: Three Step-2 chemical structure cards with clean 2-line captions inside borders
    y2, h2 = 1.8, 19.6
    step2_specs = [
        (
            2.0, 31.2, "sam",
            "Шаг 2А: (S,S)-Адеметионин (SAMe)",
            "C₁₅H₂₃N₆O₅S⁺ · 399,4 г/моль",
            "Подпитка Met + депо в вакуоли (10–17 г/л)\nРынок РФ: 14–16 млрд ₽/год (ЖНВЛП; п. 24 МСХ)",
            "#2F5573", "#F4F8FB",
            "↓CYS4 + Met + вакуоль",
        ),
        (
            35.2, 29.6, "spermidine",
            "Шаг 2Б: Спермидин",
            "C₇H₁₉N₃ · 145,2 г/моль",
            "Декарбоксилированный SAM + путресцин\nСегмент БАД здорового долголетия и кормов",
            "#2F5573", "#FFFFFF",
            "+ гены SPE2, SPE3",
        ),
        (
            66.8, 31.2, "glutathione",
            "Шаг 2В: L-Глутатион (GSH)",
            "C₁₀H₁₇N₃O₆S · 307,3 г/моль",
            "Усиленный пул L-цистеина + L-глутамат + АТФ\nФармацевтика, косметика и пищевая индустрия",
            "#2F5573", "#FFFFFF",
            "↑CYS4 + GSH1, GSH2",
        ),
    ]

    for x, w, key, title, form, sub, brd, bg, branch_lbl in step2_specs:
        place_mol_card(x, y2, w, h2, key, title, form, sub, border=brd, bg=bg, title_col="#2F5573", lw=1.15, sub_lines=2)
        cx = x + w / 2
        ba = FancyArrowPatch((cx, 24.3), (cx, 21.5), arrowstyle="-|>,head_length=4.5,head_width=2.8", color="#2F5573", linewidth=1.3, zorder=5)
        ax.add_patch(ba)
        ax.text(cx + 0.8, 22.85, branch_lbl, fontsize=6.9, fontweight="bold", color="#2F5573", ha="left", va="center")

    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_dbtl_loop_figure(out_path: Path) -> None:
    """Figure 3: Clean visual schematic of the closed-loop DBTL workflow for yeast ergothioneine."""
    fig = plt.figure(figsize=(10.6, 3.25), dpi=300)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 31)
    ax.axis("off")

    stages = [
        (
            1.5, "1. DESIGN (Проектирование)", "#2F5573", "#EEF4F8",
            "• Стехиометрическая модель ecYeast8\n  (отсев тупиковых и токсичных мутаций)\n"
            "• Выбор комбинации промоторов для\n  6–8 узлов (NcEgt1, CpEgt2, HIS1,\n  MET14, SAM2, SAH1, Δspe2, AQR1)\n"
            "• Варьирование Fe²⁺, SO₄²⁻ и вит. B6",
        ),
        (
            26.2, "2. BUILD (Сборка штаммов)", "#1B5E4A", "#EAF2EF",
            "• Синтез кодон-оптимизированных\n  генов NcEgt1 и CpEgt2\n"
            "• Модульная сборка кассет\n  (Golden Gate / USER)\n"
            "• Безмаркерная CRISPR/Cas9-вставка\n  в хромосомы S. cerevisiae CEN.PK",
        ),
        (
            50.9, "3. TEST (Простой замер)", "#1B5E4A", "#EAF2EF",
            "• Выращивание 72 ч на прозрачной\n  минеральной среде Делфта\n"
            "• Нагрев осадка в воде (94 °C, 10 мин)\n"
            "• Спектрофотометр (311 или 343 нм)\n  в кювете / планшете за 5–10 мин\n"
            "• Контроль топ-15% клонов на ВЭЖХ",
        ),
        (
            75.6, "4. LEARN (Обучение модели)", "#2F5573", "#EEF4F8",
            "• Расчёт мг/л, г/(л·сут) и отношения\n  интермедиата [Герцинин] / [EGT]\n"
            "• Обучение ансамбля GP / XGBoost\n  (включая нулевые/слабые варианты)\n"
            "• Байесовская генерация набора\n  штаммов следующего цикла",
        ),
    ]

    bw, bh, y0 = 22.8, 22.0, 6.8
    for x, title, stroke, fill, body in stages:
        box = FancyBboxPatch(
            (x, y0), bw, bh,
            boxstyle="round,pad=0.25,rounding_size=0.9",
            facecolor=fill, edgecolor=stroke, linewidth=1.15, zorder=2
        )
        ax.add_patch(box)
        ax.text(x + bw / 2, y0 + bh - 2.3, title, fontsize=8.1, fontweight="bold", color=stroke, ha="center", va="center")
        ax.plot([x + 1.0, x + bw - 1.0], [y0 + bh - 4.1, y0 + bh - 4.1], color=stroke, linewidth=0.6, alpha=0.5)
        ax.text(x + 1.1, y0 + bh - 5.3, body, fontsize=6.95, color="#181D20", ha="left", va="top", linespacing=1.32)

    # Forward arrows between stages 1 -> 2 -> 3 -> 4
    for x_end in (24.3, 49.0, 73.7):
        arr = FancyArrowPatch(
            (x_end + 0.1, y0 + bh / 2), (x_end + 1.8, y0 + bh / 2),
            arrowstyle="-|>,head_length=4.5,head_width=2.6",
            color="#181D20", linewidth=1.3, zorder=4
        )
        ax.add_patch(arr)

    # Feedback return arrow from LEARN (stage 4) back to DESIGN (stage 1) along bottom
    ax.plot([87.0, 87.0, 12.9, 12.9], [y0 - 0.1, 2.3, 2.3, y0 - 1.6], color="#2F5573", linewidth=1.2, linestyle="--")
    fb_arr = FancyArrowPatch((12.9, 4.2), (12.9, y0 - 0.1), arrowstyle="-|>,head_length=4.5,head_width=2.6", color="#2F5573", linewidth=1.2, zorder=4)
    ax.add_patch(fb_arr)
    ax.text(
        50.0, 3.3,
        "Итеративная обратная связь (Циклы 1 → 3): уточнение узкого места и рекомендации новых комбинаций промоторов",
        fontsize=7.3, fontweight="bold", color="#2F5573", ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none")
    )

    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_chart_markets(out_path: Path) -> None:
    """Figure 4: Global EGT growth vs Lycopene + Russian hepatoprotector market (SAMe)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 3.95), dpi=300)
    fig.patch.set_facecolor("white")

    years = np.array([2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031])
    egt_real = np.array([63, 75, 88, 104, 123, 145, 171, 202])
    egt_promo = np.array([58, 78, 108, 146, 195, 255, 325, 432])
    lyco_mkt = np.array([152, 159, 166, 174, 182, 190, 199, 207])

    ax1.plot(years, egt_real, marker="o", color="#1B5E4A", linewidth=2.3, markersize=5, label="EGT: базовая оценка (~68 т, CAGR 18%)")
    ax1.fill_between(years, egt_real, color="#1B5E4A", alpha=0.10)
    ax1.plot(years, egt_promo, color="#6A9E8E", linewidth=1.3, linestyle=":", label="EGT: промо-прогнозы (CAGR 33%)")
    ax1.plot(years, lyco_mkt, marker="s", color="#B84A39", linewidth=1.7, linestyle="--", markersize=4.2, label="Ликопин (CAGR ~4,5%)")

    ax1.annotate("$75M (~68 т)", (2025, 75), xytext=(6, -15), textcoords="offset points", ha="left", fontsize=7.2, fontweight="bold", color="#1B5E4A")
    ax1.annotate("$202M", (2031, 202), xytext=(-14, -15), textcoords="offset points", ha="center", fontsize=7.3, fontweight="bold", color="#1B5E4A")
    ax1.annotate("$432M", (2031, 432), xytext=(-8, 6), textcoords="offset points", ha="center", fontsize=7.0, color="#6A9E8E")
    ax1.annotate("$159M", (2025, 159), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7.1, color="#B84A39")

    ax1.set_ylim(0, 495)
    ax1.grid(True, linestyle="--", linewidth=0.5, color="#E1E6EA")
    ax1.set_ylabel("Мировой рынок субстанции ($ млн)", fontsize=8.4, color="#2E373D")
    ax1.set_xlabel("Год (базовая оценка Nexchem / Emergen vs агрегаторы)", fontsize=8.1, color="#2E373D")
    ax1.set_title("А. Мировой рынок EGT (узкий тоннаж) и ликопина", fontsize=9.3, fontweight="bold", color="#181D20", pad=8)
    ax1.legend(frameon=False, fontsize=7.1, loc="upper left")
    ax1.tick_params(labelsize=8, colors="#4A555E")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    cats = [
        "Адеметионин (SAMe)\nГептрал, Гептор,\nСамеликс, Гепцифол",
        "Эссенциальные\nфосфолипиды\n(Эссенциале и др.)",
        "Глицирризиновая к-та\n+ фосфолипиды\n(Фосфоглив)",
        "Прочие средства A05B\n(Карсил, Гепа-Мерц,\nРемаксол, Лив-52)",
    ]
    rub_bn = [13.0, 5.0, 3.4, 6.2]
    shares = ["13,0 млрд ₽\n(47%)", "5,0 млрд ₽\n(18%)", "3,4 млрд ₽\n(12%)", "6,2 млрд ₽\n(22%)"]
    cols2 = ["#1B5E4A", "#7A8B96", "#9AA8B2", "#B8C2C8"]

    x_pos = np.arange(len(cats))
    bars = ax2.bar(x_pos, rub_bn, width=0.54, color=cols2, zorder=3)
    for bar, lbl, idx in zip(bars, shares, range(len(cats))):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.35,
            lbl,
            ha="center",
            va="bottom",
            fontsize=7.5,
            fontweight="bold" if idx == 0 else "normal",
            color="#1B5E4A" if idx == 0 else "#2E373D",
        )

    ax2.set_ylim(0, 16.2)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(cats, fontsize=7.0, color="#2E373D")
    ax2.grid(True, axis="y", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax2.set_ylabel("Аптечные продажи в РФ за 11 мес. 2025 г. (млрд ₽)", fontsize=8.2, color="#2E373D")
    ax2.set_title("Б. Рынок Шага 2 в РФ: доминирование SAMe (RNC Pharma)", fontsize=9.3, fontweight="bold", color="#181D20", pad=8)
    ax2.tick_params(axis="y", labelsize=8, colors="#4A555E")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout(pad=1.2)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


class ReportDocTemplate(BaseDocTemplate):
    """Custom DocTemplate that records exact bookmark page numbers for the clickable TOC."""

    def __init__(self, filename, toc_map: dict, **kw):
        super().__init__(filename, **kw)
        self.toc_map = toc_map

    def afterFlowable(self, flowable):
        if isinstance(flowable, BookmarkHeading):
            self.canv.bookmarkPage(flowable.toc_key)
            self.canv.addOutlineEntry(flowable.toc_title, flowable.toc_key, level=flowable.toc_level - 1, closed=False)
            self.toc_map[flowable.toc_key] = self.page


def build_pdf() -> None:
    fig1_path = ASSETS_DIR / "fig1_candidates.png"
    fig2_path = ASSETS_DIR / "fig2_pathway_structures.png"
    fig3_path = ASSETS_DIR / "fig3_dbtl_loop.png"
    fig4_path = ASSETS_DIR / "fig4_markets.png"

    make_chart_alternatives(fig1_path)
    make_pathway_structures_figure(fig2_path)
    make_dbtl_loop_figure(fig3_path)
    make_chart_markets(fig4_path)

    page_w, page_h = A4
    left_m = right_m = 18 * mm
    top_m = 19 * mm
    bot_m = 16.5 * mm
    usable_w = page_w - left_m - right_m

    # Typography styles — generous line height and paragraph spacing for high readability
    st_kicker = ParagraphStyle("Kicker", fontName="SansBold", fontSize=8.2, leading=10.5, textColor=ACCENT_TEAL, spaceAfter=2 * mm)
    st_title = ParagraphStyle("DocTitle", fontName="SansBold", fontSize=18.0, leading=22.5, textColor=INK, spaceAfter=2.5 * mm)
    st_subtitle = ParagraphStyle("DocSubtitle", fontName="Sans", fontSize=9.6, leading=13.6, textColor=MUTED, spaceAfter=4.5 * mm)

    st_h1 = ParagraphStyle("H1", fontName="SansBold", fontSize=13.2, leading=17.0, textColor=INK, spaceBefore=0, spaceAfter=3.2 * mm, keepWithNext=True)
    st_h2 = ParagraphStyle("H2", fontName="SansBold", fontSize=10.5, leading=14.0, textColor=ACCENT_TEAL, spaceBefore=4.2 * mm, spaceAfter=1.8 * mm, keepWithNext=True)
    st_h3 = ParagraphStyle("H3", fontName="SansBold", fontSize=9.5, leading=13.0, textColor=SUBINK, spaceBefore=2.8 * mm, spaceAfter=1.2 * mm, keepWithNext=True)

    st_body = ParagraphStyle("Body", fontName="Sans", fontSize=9.3, leading=13.8, textColor=INK, alignment=TA_JUSTIFY, spaceAfter=2.3 * mm)
    st_bullet = ParagraphStyle("Bullet", fontName="Sans", fontSize=9.2, leading=13.6, textColor=INK, leftIndent=4.8 * mm, firstLineIndent=-4.8 * mm, alignment=TA_JUSTIFY, spaceAfter=1.7 * mm)
    st_subbullet = ParagraphStyle("SubBullet", fontName="Sans", fontSize=8.95, leading=13.1, textColor=SUBINK, leftIndent=9.5 * mm, firstLineIndent=-4.2 * mm, alignment=TA_JUSTIFY, spaceAfter=1.3 * mm)

    st_box = ParagraphStyle("CalloutBox", fontName="Sans", fontSize=8.95, leading=12.8, textColor=INK, alignment=TA_LEFT)
    st_caption = ParagraphStyle("Caption", fontName="SansItalic", fontSize=8.1, leading=11.4, textColor=MUTED, spaceBefore=1.6 * mm, spaceAfter=3.5 * mm)

    st_toc_h = ParagraphStyle("TOCHead", fontName="SansBold", fontSize=10.6, leading=13.5, textColor=INK, spaceBefore=2.0 * mm, spaceAfter=1.6 * mm)
    st_toc_1 = ParagraphStyle("TOC1", fontName="SansBold", fontSize=8.95, leading=11.8, textColor=INK)
    st_toc_2 = ParagraphStyle("TOC2", fontName="Sans", fontSize=8.55, leading=11.3, textColor=SUBINK, leftIndent=4.5 * mm)
    st_toc_pg = ParagraphStyle("TOCPg", fontName="SansBold", fontSize=8.8, leading=11.3, textColor=ACCENT_TEAL, alignment=2)  # TA_RIGHT=2

    st_th = ParagraphStyle("TH", fontName="SansBold", fontSize=7.7, leading=9.9, textColor=INK)
    st_td = ParagraphStyle("TD", fontName="Sans", fontSize=7.6, leading=10.2, textColor=SUBINK)
    st_td_bold = ParagraphStyle("TDBold", fontName="SansBold", fontSize=7.6, leading=10.2, textColor=ACCENT_TEAL)
    st_ref = ParagraphStyle("Ref", fontName="Sans", fontSize=7.8, leading=10.6, textColor=SUBINK, leftIndent=6.2 * mm, firstLineIndent=-6.2 * mm, spaceAfter=0.85 * mm)

    def draw_header_footer(canvas, doc):
        canvas.saveState()
        if doc.page > 1:
            canvas.setStrokeColor(LINE)
            canvas.setLineWidth(0.5)
            canvas.line(left_m, page_h - 13 * mm, page_w - right_m, page_h - 13 * mm)
            canvas.setFont("Sans", 7.6)
            canvas.setFillColor(MUTED)
            canvas.drawString(left_m, page_h - 11 * mm, "Обоснование выбора L-эрготионеина и метаболического узла SAM для цикла DBTL")
            canvas.drawRightString(page_w - right_m, page_h - 11 * mm, "Ноябрь 2026")

        # Minimalist footer: ONLY page number, no extra text
        canvas.setStrokeColor(LINE)
        canvas.setLineWidth(0.45)
        canvas.line(left_m, 11.8 * mm, page_w - right_m, 11.8 * mm)
        canvas.setFont("Sans", 8.2)
        canvas.setFillColor(MUTED)
        canvas.drawCentredString(page_w / 2, 7.8 * mm, f"{doc.page}")
        canvas.restoreState()

    toc_entries_spec = [
        ("sec1", 1, "1. Итоговое решение и критерии отбора целевого продукта"),
        ("sec1_1", 2, "1.1. Выбор класса задачи и целевой метрики для цикла DBTL"),
        ("sec1_2", 2, "1.2. Пять согласованных критериев отбора целевого вещества"),
        ("sec1_3", 2, "1.3. Сводное обоснование связки L-эрготионеин → адеметионин (SAMe)"),
        ("sec2", 1, "2. Сравнительный анализ альтернативных вариантов"),
        ("sec2_1", 2, "2.1. Каротиноиды (ликопин, бета-каротин, астаксантин) и апокаротиноиды (шафран, аннато)"),
        ("sec2_2", 2, "2.2. Фенилпропаноиды (п-кумаровая кислота, ресвератрол, нарингенин)"),
        ("sec2_3", 2, "2.3. Промышленные кандидаты и позиции государственных перечней РФ"),
        ("sec2_4", 2, "2.4. Сводная таблица и количественное сопоставление кандидатов"),
        ("sec3", 1, "3. Биохимический путь L-эрготионеина, синтез предшественников и модель DBTL"),
        ("sec3_1", 2, "3.1. Нативный синтез предшественников в дрожжах и воспроизводимость первого шага"),
        ("sec3_2", 2, "3.2. Химическая схема пути (Шаг 1) и биохимические условия перехода к Шагу 2"),
        ("sec3_3", 2, "3.3. Топология метаболического узла и архитектура моделирования Design–Learn"),
        ("sec4", 1, "4. Пошаговый лабораторный протокол и план экспериментов"),
        ("sec4_1", 2, "4.1. Приготовление прозрачной минеральной среды Делфта (на 1 литр)"),
        ("sec4_2", 2, "4.2. Сборка генетической кассеты и трансформация дрожжей (Build)"),
        ("sec4_3", 2, "4.3. Выращивание культур и водная термоэкстракция при 94 °C"),
        ("sec4_4", 2, "4.4. Быстрое измерение концентрации на спектрофотометре (311 и 343 нм) и ВЭЖХ"),
        ("sec4_5", 2, "4.5. Календарный план итераций Циклов 0–3 и перехода к Шагу 2"),
        ("sec5", 1, "5. Экономика, регуляторика и критический стресс-тест рынков мира и РФ"),
        ("sec5_1", 2, "5.1. Мировой рынок L-эрготионеина, профицит мощностей в КНР и ценовая реальность"),
        ("sec5_2", 2, "5.2. Шаг 1 в РФ (микрониша и механика пилота с брендом) и масштаб Шага 2 (SAMe)"),
        ("sec5_3", 2, "5.3. Критический стресс-тест концепции и развилка траекторий (Сценарий А vs Сценарий Б)"),
        ("sec6", 1, "6. Список литературы и источников"),
    ]

    def callout_box(html_text: str, bg=ACCENT_SOFT, border=ACCENT_TEAL) -> Table:
        t = Table([[P(html_text, st_box)]], colWidths=[usable_w])
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), bg),
                    ("LINEBEFORE", (0, 0), (0, -1), 2.6, border),
                    ("LEFTPADDING", (0, 0), (-1, -1), 9),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 9),
                    ("TOPPADDING", (0, 0), (-1, -1), 6.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6.5),
                ]
            )
        )
        return t

    def build_story(toc_pages: dict):
        story = []

        def H1(key: str, text: str):
            return BookmarkHeading(text, st_h1, key, 1, text)

        def H2(key: str, text: str):
            return BookmarkHeading(text, st_h2, key, 2, text)

        # ==================== PAGE 1: COVER SUMMARY & CLICKABLE TOC ====================
        story.append(P("АНАЛИТИЧЕСКИЙ ОТЧЁТ · ПРОЕКТ DBTL · НОЯБРЬ 2026", st_kicker))
        story.append(
            P(
                "Обоснование выбора L-эрготионеина и метаболического узла SAM для отработки цикла Design–Build–Test–Learn",
                st_title,
            )
        )
        story.append(
            P(
                "Критерии выбора целевой задачи, сравнительный аудит 16 биотехнологических кандидатов, химическая схема пути, пошаговый лабораторный протокол на дрожжах и анализ рынков мира и РФ",
                st_subtitle,
            )
        )

        hr = Table([[""]], colWidths=[usable_w], rowHeights=[0.5])
        hr.setStyle(TableStyle([("LINEABOVE", (0, 0), (-1, -1), 0.9, INK), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5 * mm)]))
        story.append(hr)

        story.append(
            callout_box(
                "<b>Краткое резюме и главный вывод фактчека.</b> Сопоставлены 16 биотехнологических кандидатов для запуска цикла Design–Build–Test–Learn (DBTL). Поскольку молекулы с простым кюветным тестом без ВЭЖХ имеют в РФ малый физический тоннаж, базовый маршрут разделяет задачи двух этапов: <b>Шаг 1 — L-эрготионеин как модельный полигон конвейера DBTL, база статьи Q1 и R&amp;D-пилота с косметическим брендом РФ</b> (2 гена <i>NcEgt1, CpEgt2</i> дают 46–154 мг/л <i>de novo</i>, тест 311/343 нм, 0 ML-статей; рынок в РФ — единицы–десятки кг/год на фоне профицита КНР по $400–580/кг, но низкий регуляторный порог ТР ТС 009/2011 и 4-стадийная мембранно-угольная очистка позволяют наработать в 10-л реакторе 10–20 кг 1%-го концентрата на серию 15 000 сывороток и привлечь индустриального партнёра уровня «Гельтек» / «Ангиофарм» уже на Шаге 1). <b>Шаг 2 — (S,S)-адеметионин (SAMe)</b> с рынком в РФ <b>14–16 млрд руб./год</b> (ЖНВЛП «Фарма-2030», п. 24 МСХ № 392), куда переносится конвейер DBTL и блок <i>SAM2/АТФ</i> с разворотом <i>↓CYS4</i> под подачу метионина (либо прямой старт с SAMe по <b>Сценарию Б</b>)."
            )
        )
        story.append(Spacer(1, 3 * mm))

        story.append(P("Оглавление (разделы кликабельны)", st_toc_h))
        toc_rows = []
        for key, lvl, title in toc_entries_spec:
            pg_num = toc_pages.get(key, "—")
            st_item = st_toc_1 if lvl == 1 else st_toc_2
            link_html = f'<a href="#{key}" color="#181D20">{fmt_chem(title)}</a>'
            pg_html = f'<a href="#{key}" color="#1B5E4A">{pg_num}</a>'
            toc_rows.append([Paragraph(link_html, st_item), Paragraph(pg_html, st_toc_pg)])

        toc_tbl = Table(toc_rows, colWidths=[usable_w - 14 * mm, 14 * mm])
        toc_styles = [
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 1),
            ("RIGHTPADDING", (0, 0), (-1, -1), 1),
            ("TOPPADDING", (0, 0), (-1, -1), 1.45),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.45),
            ("LINEBELOW", (0, 0), (-1, -1), 0.3, colors.HexColor("#E6ECEF")),
        ]
        for idx, (_, lvl, _) in enumerate(toc_entries_spec):
            if lvl == 1 and idx > 0:
                toc_styles.append(("TOPPADDING", (0, idx), (-1, idx), 3.2))
        toc_tbl.setStyle(TableStyle(toc_styles))
        story.append(toc_tbl)

        # ==================== SECTION 1 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec1", "1. Итоговое решение и критерии отбора целевого продукта"))

        story.append(H2("sec1_1", "1.1. Выбор класса задачи и целевой метрики для цикла DBTL"))
        story.append(
            P(
                "Методология Design–Build–Test–Learn (проектирование — сборка — тестирование — обучение, DBTL) используется в инженерной биологии для разных задач: от мутагенеза одиночных ферментов и подбора питательных сред до перестройки клеточного метаболизма [1–3]. При запуске первого итеративного конвейера критично выбрать такую постановку, в которой результат каждого раунда выражается строгим числом и напрямую подаётся на вход алгоритма машинного обучения.",
                st_body,
            )
        )
        story.append(
            P(
                "По этой причине из возможных вариантов выбрана задача <b>биосинтеза конкретного низкомолекулярного вещества с максимизацией концентрации и объёмной скорости его накопления (г/л·сут)</b>. В отличие от сложных физиологических фенотипов, выбора самого вида микроорганизма или реологических характеристик полимеров, концентрация целевой молекулы задаёт <b>однозначную численную метрику качества, на которую мы опираемся в цикле DBTL</b>.",
                st_body,
            )
        )

        story.append(H2("sec1_2", "1.2. Пять согласованных критериев отбора целевого вещества"))
        story.append(
            P(
                "Чтобы проект был воспроизводим в стандартной молекулярно-биологической лаборатории, давал содержательную задачу для моделирования и одновременно имел коммерческую ценность, целевой продукт должен удовлетворять пяти критериям:",
                st_body,
            )
        )
        criteria = [
            "<b>1. Понятный первый шаг с опубликованным референсом и действующим промышленным производством.</b> Молекулу первого шага уже умеют синтезировать в биологических (клеточных) продуцентах — бактериях, дрожжах или мицелиальных грибах: точно известны вставляемые гены, описаны рабочие лабораторные методики и стартовые выходы. При этом в мире уже существует промышленное производство данного вещества, показатели которого (скорость синтеза, отказ от дорогих предшественников в среде, выход вещества наружу из клеток) ставится задача превзойти за счёт циклов DBTL.",
            "<b>2. Быстрое и простое измерение результата цикла в наших лабораторных условиях.</b> Концентрация продукта после культивирования должна быстро и надёжно измеряться простым оптическим методом без масс-спектрометрии (ВЭЖХ-МС / ГХ-МС) и длительной хроматографии каждого образца — например, на обычном лабораторном спектрофотометре в кюветах (или, при увеличении числа проверяемых штаммов, в 96-луночных планшетах). Хроматография применяется только для контрольной проверки лучших отобранных клонов.",
            "<b>3. Самостоятельный рынок в мире и в РФ и преимущество биосинтеза уже на первом шаге.</b> Первое же вещество должно иметь растущий платёжеспособный спрос в мире и подтверждённый импортозависимый рынок с понятными покупателями в России, причём биологический синтез клетками-продуцентами должен быть экономически и технологически выгоднее химического синтеза и выделения из природного сырья.",
            "<b>4. Перспективные и уникальные вещества следующих этапов («Шаг 2»).</b> Созданный и оптимизированный на первом шаге штамм с усиленным синтезом предшественников и кофакторов должен служить готовой платформой для перехода (добавлением 1–3 генов или переключением потока без смены продуцента) к более сложным высокомаржинальным молекулам с крупным рынком.",
            "<b>5. Вхождение целевых продуктов в государственные приоритетные перечни РФ.</b> Продукты первого и/или последующих шагов должны входить в отраслевые перечни критического импортозамещения — Приказ Минсельхоза России № 392 от 06.06.2025, стратегию Минпромторга «Фарма-2030» и продуктовые цепочки нацпроекта «Технологическое обеспечение биоэкономики» (Распоряжение Правительства РФ № 2777-р от 05.10.2026) [36].",
        ]
        for c_item in criteria:
            story.append(P("• " + c_item, st_bullet))

        story.append(H2("sec1_3", "1.3. Сводное обоснование связки L-эрготионеин → адеметионин (SAMe)"))
        story.append(
            P(
                "По результатам проверки шестнадцати кандидатов всем пяти критериям одновременно отвечает двухэтапная платформа на дрожжевом продуценте: <b>Шаг 1 — L-эрготионеин</b>, <b>Шаг 2 — (S,S)-адеметионин (SAMe), спермидин и глутатион</b>:",
                st_body,
            )
        )
        summary_pts = [
            "<b>Простой и гарантированный старт (Критерий 1):</b> дрожжи исходно синтезируют все три предшественника (гистидин, SAM и цистеин) из обычной глюкозы и сульфата; вставка всего двух грибных генов (<i>NcEgt1</i> и <i>CpEgt2</i>) сразу даёт 46–154 мг/л эрготионеина без нокаутов хозяйских генов [27, 29, 32], а на рынке КНР уже работает крупное ферментационное производство [40], задающее жёсткий технологический ориентир.",
            "<b>Быстрое измерение без масс-спектрометра (Критерий 2):</b> эрготионеин во внеклеточном супернатанте измеряется за 5–10 минут на обычном спектрофотометре при 311 нм (после разовой наработки рекомбинантной эрготионазы <i>TdETL</i>) или 343 нм (реактив 2-PDS при pH 1,0, где протонированные клеточные тиолы подавлены, с обязательной калибровкой по ВЭЖХ-УФ на внутриклеточных пробах) [33, 34].",
            "<b>Трезвая оценка рынка Шага 1 и преимущество биосинтеза (Критерий 3):</b> в плодовых телах грибов эрготионеина мало (0,1–2 мг/г), а химсинтез сложен и даёт примесь D-изомера, поэтому индустрия перешла на ферментацию (~$75 млн / ~68 т в мире в 2025 г., реалистичный рост ~18%/год) [40]. В РФ 100% сырья импортируется и нет патентов на рекомбинантные дрожжи с <i>Egt1/Egt2</i>, но из-за микродозировок в косметике (0,01–0,05%) физический спрос в РФ составляет единицы–десятки кг/год на фоне профицита мощностей КНР ($400–580/кг). Поэтому Шаг 1 позиционируется именно как <b>калибровочный полигон DBTL, база для статьи Q1 и мелкосерийных B2B-партий</b>, а не отдельный завод.",
            "<b>Продукты Шага 2, биохимическая перестройка и списки РФ (Критерии 4 и 5):</b> главный коммерческий целевой продукт платформы — <b>(S,S)-адеметионин (SAMe)</b> с рынком в РФ 14–16 млрд руб./год (ЖНВЛП «Фарма-2030», кормовые соединения метионина п. 24 Приказа МСХ № 392) [36, 37, 43]. При переходе от EGT к SAMe переносится весь конвейер DBTL и энергетический блок <i>SAM2/АТФ</i>, тогда как ветви <i>HIS1</i> и <i>CYS4</i> целенаправленно разворачиваются под подачу L-метионина (раздел 3.2; прямое сравнение со стартом сразу на SAMe дано в разделе 5.3).",
        ]
        for pt in summary_pts:
            story.append(P("• " + pt, st_bullet))

        # ==================== SECTION 2 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec2", "2. Сравнительный анализ альтернативных вариантов"))
        story.append(
            P(
                "Ниже приведены результаты анализа альтернативных молекул с точки зрения лабораторной сложности старта, удобства измерения, научной новизны применения ML-DBTL и реальной экономики на российском рынке.",
                st_body,
            )
        )

        story.append(H2("sec2_1", "2.1. Каротиноиды (ликопин, бета-каротин, астаксантин) и апокаротиноиды (шафран, аннато)"))
        story.append(
            P(
                "Каротиноиды удобны на первом шаге за счёт видимой окраски клеток (поглощение при 450–475 нм) и наличия готовых штаммов-референсов, однако проигрывают по научной новизне и экономике:",
                st_body,
            )
        )
        carot_bullets = [
            "<b>Исчерпанная научная новизна в DBTL:</b> именно на каротиноидах в <i>S. cerevisiae</i> и <i>Y. lipolytica</i> были отработаны почти все классические алгоритмы машинного обучения для метаболической инженерии — байесовская платформа <i>BioAutomata</i> на ликопине (2019) [4], нейросетевой конвейер <i>MiYA</i> (2018) [6], студенческие проекты iGEM <i>NEFU-China</i> (2022) и <i>NACIS Shanghai</i> (2025) [8], стекинг 24 моделей ML (<i>Shimazaki et al., 2024</i>) [10], девятираундовая ML-оптимизация бета-каротина в MIT (<i>Tu &amp; Stephanopoulos, 2026</i>) [11] и модель с ферментными ограничениями <i>eciYali5-GEM</i> (2026) [12]. Титры в дрожжах уже доведены до 17,6 г/л по ликопину, 39,5 г/л по бета-каротину и 2,82 г/л по астаксантину [9, 15].",
            "<b>Слабая экономика ликопина и астаксантина:</b> мировой рынок ликопина невелик (~$159 млн, рост ~4,5% в год) и занят дешёвым химическим синтезом (<i>BASF, dsm-firmenich</i>), томатным экстрактом и промышленной ферментацией мицелиального гриба <i>Blakeslea trispora</i>. По предварительному техноэкономическому расчёту для линии 1 000 л в РФ при титре 5 г/л себестоимость 10%-го ликопина (12,4 тыс. руб./кг) сравнима с ценой готового импорта.",
            "<b>Барьеры разворота из ликопина в апокаротиноиды — кроцин шафрана (3А) и аннато E160b (3Б):</b>",
        ]
        for b in carot_bullets:
            story.append(P("• " + b, st_bullet))

        story.append(
            P(
                "– <i>Аннато E160b (биксин / норбиксин):</i> семена тропического кустарника <i>Bixa orellana</i> содержат 2–5% пигмента при цене сырья $2,5–3/кг, поэтому готовый краситель E160b продаётся в России по 1 700–5 100 руб./кг. При этом гены <i>BoLCD</i> из ранней статьи 2003 г. не подтвердились при расшифровке генома растения [14], и работающего пути биксина в дрожжах не существует (0 мг/л).",
                st_subbullet,
            )
        )
        story.append(
            P(
                "– <i>Кроцетин и кроцины шафрана:</i> экстракты дороги ($300–3 500/кг), но в дрожжах возникает пространственное разобщение субстрата и фермента: промежуточный жирорастворимый зеаксантин накапливается внутри гидрофобных липидных капель клетки, куда не проникает водорастворимый фермент цитозоля диоксигеназа <i>CsCCD2</i>. В работе <i>Zhou et al. (2025)</i> из 1 575 мг/л зеаксантина в <i>Y. lipolytica</i> расщепилось менее 2% (30,2 мг/л кроцетина за 168 ч при вынужденном охлаждении реактора до 20 °C) [13], а собственные глюкозидазы дрожжей (<i>EXG1/EXG2</i>) разрушают целевой кроцин обратно. По этим причинам проект дрожжевого шафрана компании <i>Evolva</i> (WO 2014/036130) был закрыт.",
                st_subbullet,
            )
        )

        story.append(H2("sec2_2", "2.2. Фенилпропаноиды (п-кумаровая кислота, ресвератрол, нарингенин)"))
        phenyl_bullets = [
            "<b>Отсутствие экономики на первом шаге:</b> <i>п</i>-кумаровая кислота легко получается в дрожжах вставкой одного гена тирозинаммиаклиазы <i>TAL</i> (до 12,5 г/л, поглощение в среде при 310 нм) [17], но в промышленности её синтезируют в одну химическую стадию конденсацией Кнёвенагеля из дешёвого 4-гидроксибензальдегида по $90–250/кг.",
            "<b>Закрытая ниша на втором шаге:</b> добавление 2–3 генов (<i>4CL, STS / CHS</i>) даёт бесцветные антиоксиданты для БАД — ресвератрол (до 22,5 г/л [18]) и нарингенин (до 8,65 г/л). В 2024 г. группы Вагенингена и <i>dsm-firmenich</i> уже опубликовали полный цикл ML-DBTL на 18 генах этого пути [16]. Коммерчески дешёвый растительный экстракт горца японского из КНР ($60–120/кг за 98% ресвератрола) привёл пионера дрожжевой ферментации <i>Evolva AG</i> к многолетним убыткам (~CHF 14,6 млн/год) и продаже бизнеса в декабре 2023 г. за CHF 20 млн, а в РФ рынок антиоксидантных БАД занят отечественным дигидрокверцетином из сибирской лиственницы (АО «Аметис», >33 т/год по 40–90 тыс. руб./кг). Удлинение же пути до редких изофлаванов (глабридин, бакучиол) требует +6…+12 генов с падением титров в дрожжах до 1,2–9,3 мг/л [18].",
        ]
        for b in phenyl_bullets:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec2_3", "2.3. Промышленные кандидаты и позиции государственных перечней РФ"))
        story.append(
            P(
                "Отдельно проверены пять специализированных промышленных и фармацевтических молекул, четыре направления предварительного патентного отсева (эктоин, хондроитинсульфат, глюкозамин и бензоназоподобные нуклеазы, по которым в РФ уже выпускается аналог <i>Dreamnase</i> либо опубликованы дрожжевые работы 2025–2026 гг.), а также позиции Приказа Минсельхоза РФ № 392 [36]:",
                st_body,
            )
        )
        gov_bullets = [
            "<b>Пуллулан (пищевая добавка E1204, фармацевтические капсулы):</b> внеклеточный полисахарид продуцируется не модельными дрожжами, а полиморфным «чёрным грибом» <i>Aureobasidium pullulans</i>, параллельно синтезирующим тёмный пигмент меланин (титры до 140,2 г/л [19], заявка WO2026050525A1). Для отработки DBTL он не подходит по трём причинам: (1) культуральная жидкость превращается в густой гель вязкостью до 6 270 мПа·с, что затрудняет перемешивание и отбор проб; (2) ключевое коммерческое свойство полимера — распределение длины цепей (молекулярной массы), которое нельзя измерить на спектрофотометре (нужна гель-проникающая хроматография каждого образца) и нельзя предсказать стехиометрической моделью метаболизма; (3) оптовая цена пуллулана в КНР составляет $15–30/кг при небольшом мировом рынке ~$80 млн [20].",
            "<b>Витамин B7 (D-биотин, поз. 38 Приказа МСХ № 392):</b> полностью импортируется для кормовых премиксов, но биологический синтез ограничен природой последнего фермента (биотинсинтаза <i>Bio2</i> у дрожжей / <i>BioB</i> у бактерий). Это радикальный SAM-фермент одноразового («суицидного») действия: он расходует собственный железосерный кластер [2Fe–2S] каждые 1–3 оборота реакции [21]. Даже с недавно открытым обходным ферментом <i>CfBio1</i> (<i>PNAS, 2026</i> [22]) титры в дрожжах составляют единицы мг/л (в бактериях — до 209 мг/л), молекула не поглощает УФ/видимый свет и не выдерживает ценовой конкуренции с химсинтезом КНР.",
            "<b>Рекомбинантный химозин (поз. 12 Приказа МСХ № 392):</b> молокосвёртывающий фермент для сыроделия (рынок РФ ~730 млн руб./год) легко тестируется по скорости коагуляции молока, но представляет собой задачу белковой секреции одного гена в <i>Komagataella phaffii</i>, а не оптимизации метаболического пути; не имеет продуктов «Шага 2» и уже плотно закрыт в РФ патентом ФИЦ Биотехнологии РАН RU 2769175 [23] и линейкой химозинов консорциума АлтГУ, НГУ, ГНЦ «Вектор» и СибНИИС.",
            "<b>Таурин и фолиновая кислота (кальция фолинат):</b> таурин производится нефтехимическим синтезом из этиленоксида по $1,8–3,0/кг, тогда как в инженерных <i>S. cerevisiae</i> в 2025–2026 гг. получено лишь 2–38 мг/л [24, 25]. Фолиновая кислота (5-формилтетрагидрофолат) накапливается в дрожжах в виде нестабильной смеси полиглутаматов на уровне 0,36 мг/л [38] при малом объёме онкологического рынка и жёстких требованиях фарм-GMP к разделению изомеров.",
            "<b>Перспективные резервные мишени из госперечней — 7-дегидрохолестерин → 25-OH-D3 (поз. 29 МСХ № 392) и 2'-фукозиллактоза (нацпроект «Биоэкономика»):</b> 100% импорт в РФ (покупатель 2'-FL — завод детских смесей «Инфаприм»), но для первого цикла DBTL они неудобны: спектр поглощения 7-дегидрохолестерина (282 нм) полностью перекрывается собственным эргостеролом дрожжей (без ВЭЖХ штаммы неразличимы), а по 2'-фукозиллактозе нет прямого УФ-сигнала и действуют широкие патентные блоки <i>dsm-firmenich</i> и <i>Chr. Hansen</i>.",
        ]
        for b in gov_bullets:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec2_4", "2.4. Сводная таблица и количественное сопоставление кандидатов"))

        tbl_data = [
            [
                P("Кандидат / платформа", st_th),
                P("Сложность старта и гены", st_th),
                P("Быстрое измерение без МС / ВЭЖХ", st_th),
                P("Статус ML-DBTL в дрожжах", st_th),
                P("Экономика Шага 1 и рынок РФ", st_th),
                P("Шаг 2 и госперечни РФ", st_th),
            ],
            [
                P("<b>L-Эрготионеин → SAMe</b><br/>(Сценарии А / Б)", st_td_bold),
                P("<b>Низкая:</b> 2 гена (<i>NcEgt1, CpEgt2</i>), сразу 46–154 мг/л", st_td),
                P("<b>Да (для EGT):</b> 311 или 343 нм; для SAMe — быстрый ВЭЖХ 260 нм", st_td),
                P("<b>0 ML-статей по EGT</b>; нет патентов РФ на рекомб. дрожжи", st_td),
                P("<b>Нишевая на Шаге 1:</b> ~68 т в мире, единицы–десятки кг в РФ, опт КНР $400–580/кг", st_td),
                P("<b>Сильный Шаг 2:</b> SAMe — 14–16 млрд ₽ в РФ (ЖНВЛП Фарма-2030, п. 24 МСХ); реконфигурация <i>↓CYS4</i>", st_td),
            ],
            [
                P("Ликопин → Астаксантин", st_td),
                P("Низкая: 3–5 генов (есть штаммы ВКПМ)", st_td),
                P("Да: поглощение 450–475 нм (после экстракции клеток)", st_td),
                P("Десятки статей (MIT 2026, MiYA, BioAutomata)", st_td),
                P("Слабая: рынок $159M, занят химсинтезом BASF и <i>B. trispora</i>", st_td),
                P("Астаксантин занят синтетикой и водорослью <i>H. pluvialis</i>", st_td),
            ],
            [
                P("Ликопин → Кроцин (3А) / Аннато (3Б)", st_td),
                P("Высокая: +4–5 генов, охлаждение реактора до 20 °C", st_td),
                P("Частично: полоса 440 нм в водной фазе", st_td),
                P("0 статей, но конверсия из капель <2%", st_td),
                P("Аннато в РФ 1,7–5 тыс. ₽/кг; проект шафрана Evolva закрыт", st_td),
                P("Нет платформенного расширения", st_td),
            ],
            [
                P("<i>п</i>-Кумаровая к-та → Ресвератрол", st_td),
                P("Низкая: 1–3 гена (<i>TAL, 4CL, STS</i>)", st_td),
                P("Да: УФ 310 / 390 нм на спектрофотометре", st_td),
                P("Опубликовано (Wageningen / DSM 2024)", st_td),
                P("Убыточно против синтеза Кнёвенагеля и горца КНР; продажа Evolva", st_td),
                P("Глабридин / бакучиол дают 1,2–9,3 мг/л (+6…12 генов)", st_td),
            ],
            [
                P("Пуллулан (E1204)", st_td),
                P("Сложная: гриб <i>A. pullulans</i>, побочный меланин", st_td),
                P("Нет: вязкость 6 270 мПа·с, нужна хроматография ГПХ", st_td),
                P("Модель метаболизма не предсказывает длину цепи", st_td),
                P("Узкая маржа: опт КНР $15–30/кг, весь рынок ~$80M", st_td),
                P("Одиночный гомополимер без Шага 2", st_td),
            ],
            [
                P("Биотин B7 / Таурин / Фолинат", st_td),
                P("Средняя–высокая (одноразовый фермент <i>Bio2</i> у B7)", st_td),
                P("Нет удобного оптического сигнала без ВЭЖХ-МС", st_td),
                P("Титры в дрожжах лишь 0,4–38 мг/л", st_td),
                P("Себестоимость на 2–3 порядка хуже химсинтеза КНР", st_td),
                P("Биотин — п. 38 МСХ № 392, но экономика не сходится", st_td),
            ],
            [
                P("Рекомбинантный химозин", st_td),
                P("Низкая: секреция 1 белка в <i>K. phaffii</i>", st_td),
                P("Да: тест по свёртыванию молока", st_td),
                P("Секреция белка, а не метаболический путь", st_td),
                P("Рынок РФ ~730 млн ₽ занят патентами ФИЦ Биотех и АлтГУ", st_td),
                P("П. 12 МСХ № 392; нет перехода к Шагу 2", st_td),
            ],
        ]

        col_w = [26 * mm, 28 * mm, 27 * mm, 25 * mm, 34 * mm, 34 * mm]
        comp_tbl = Table(tbl_data, colWidths=col_w, repeatRows=1)
        comp_tbl.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
                    ("BACKGROUND", (0, 1), (-1, 1), ACCENT_SOFT),
                    ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
                    ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 3.5),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 3.5),
                    ("TOPPADDING", (0, 0), (-1, -1), 3.6),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 3.6),
                ]
            )
        )
        story.append(comp_tbl)
        story.append(P("Таблица 1. Сравнение рассмотренных кандидатов по пяти критериям пригодности для цикла DBTL и рынка РФ.", st_caption))

        img1 = Image(str(fig1_path), width=usable_w, height=usable_w * (4.35 / 10.6))
        story.append(
            KeepTogether(
                [
                    img1,
                    P(
                        "Рисунок 1. Количественное сопоставление кандидатов: (А) оптовая цена B2B и реалистичный темп роста рынка (CAGR); (Б) максимальные опубликованные титры в дрожжах и наличие работ по ML-DBTL (*для промышленного SAMe ферментация ведётся с подпиткой L-метионином).",
                        st_caption,
                    ),
                ]
            )
        )

        # ==================== SECTION 3 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec3", "3. Биохимический путь L-эрготионеина, синтез предшественников и модель DBTL"))

        story.append(H2("sec3_1", "3.1. Нативный синтез предшественников в дрожжах и воспроизводимость первого шага"))
        story.append(
            P(
                "Важнейший практический вопрос перед стартом первого цикла DBTL — <b>насколько легко получить измеримое количество продукта на базовом («нулевом») шаге без сложной перестройки генома</b>. Путь L-эрготионеина в дрожжах полностью опирается на готовый клеточный метаболизм хозяина:",
                st_body,
            )
        )
        prec_bullets = [
            "<b>Все три исходных предшественника дрожжи синтезируют сами по умолчанию:</b> при росте на простой минеральной среде (глюкоза + сульфат аммония) дикие штаммы <i>Saccharomyces cerevisiae</i> и <i>Yarrowia lipolytica</i> самостоятельно производят <b>L-гистидин</b> (ветвь пентозофосфатного пути, гены <i>HIS1–HIS7</i>), кофактор <b>S-аденозил-L-метионин (SAM)</b> (из метионина и АТФ синтазами <i>SAM1/SAM2</i>) и <b>L-цистеин</b> (восстановление сульфата SO₄²⁻ через <i>MET3 → MET14 → MET16 → MET10</i> и путь транссульфурирования <i>CYS3/CYS4</i>). При этом пекарские дрожжи от природы обладают самым высоким среди промышленных микроорганизмов внутриклеточным пулом SAM и глутатиона.",
            "<b>Короткий эукариотический модуль из 2 генов вместо 5 бактериальных:</b> бактериям (<i>Mycobacterium smegmatis</i>) для сборки эрготионеина необходимы пять отдельных ферментов (<i>EgtA–EgtE</i>) с промежуточным расходом глутамата. У нитчатых грибов ту же работу выполняют всего <b>два фермента</b> (<i>NcEgt1</i> и <i>CpEgt2</i>), которые напрямую соединяют свободный клеточный L-гистидин, SAM и L-цистеин [27].",
            "<b>Подтверждённый выход с первой же трансформации (включая студенческие команды):</b> как показано в Таблице 2, простая вставка базовой кассеты <i>NcEgt1 + CpEgt2</i> без единого нокаута хозяйских генов стабильно даёт в <i>S. cerevisiae</i> <b>46–84 мг/л</b> (а в <i>Y. lipolytica</i> — <b>154 мг/л</b>) уже в малом лабораторном объёме за 72 ч, что в 90–300 раз превышает порог чувствительности спектрофотометра (~0,5 мг/л).",
        ]
        for b in prec_bullets:
            story.append(P("• " + b, st_bullet))

        ref_tbl_data = [
            [
                P("Публикация / группа", st_th),
                P("Штамм дрожжей", st_th),
                P("Генетическое изменение на базовом шаге", st_th),
                P("Базовый выход (колбы / планшеты, 72–84 ч)", st_th),
                P("Результат дальнейшей ручной инженерии (без ML)", st_th),
            ],
            [
                P("<b>van der Hoek et al., 2019</b> (DTU Biosustain) [27]", st_td),
                P("<i>S. cerevisiae</i> CEN.PK113-7D", st_td),
                P("Вставка только 2 генов (<i>NcEgt1 + CpEgt2</i>), 0 нокаутов хозяина", st_td),
                P("<b>46–84 мг/л</b> на минеральной среде (~360 мкМ)", st_td_bold),
                P("598 мг/л в 1-л реакторе (84 ч, с подпиткой аминокислотами)", st_td),
            ],
            [
                P("<b>Команда iGEM HBUT-China, 2022</b> [32]", st_td),
                P("<i>S. cerevisiae</i> CEN.PK113-5D", st_td),
                P("Вставка 1 и 2 копий <i>NcEgt1 + CpEgt2</i> в локус <i>PDC1–STU2</i>", st_td),
                P("<b>Воспроизведено с 1-й сборки</b> (+40,6% от 2-й копии)", st_td_bold),
                P("Подтверждена реализуемость силами студенческой группы за 4 нед.", st_td),
            ],
            [
                P("<b>van der Hoek et al., 2022</b> [28]; <b>Ma et al., 2026</b> [41]", st_td),
                P("<i>S. cerevisiae</i> CEN.PK / BY4741", st_td),
                P("2 копии <i>Egt1 + Egt2</i> + модульная инженерия трёх ветвей (His / SAM / S)", st_td),
                P("<b>100–170 мг/л</b> в глубоких планшетах без аминокислот", st_td_bold),
                P("<b>2,39 г/л</b> (DTU, 2022) → <b>3,05 г/л</b> <i>de novo</i> без аминокислот (BUCT, апр. 2026; без ML)", st_td),
            ],
            [
                P("<b>van der Hoek et al., 2022</b> [29]; <b>Chen et al., 2026</b> [31]", st_td),
                P("<i>Yarrowia lipolytica</i> (и справка по <i>E. coli</i> [42])", st_td),
                P("Вставка только 2 генов <i>NcEgt1 + CpEgt2</i>, 0 нокаутов хозяина", st_td),
                P("<b>154 мг/л</b> в малом объёме на минеральной среде", st_td_bold),
                P("<b>8,03 г/л</b> <i>de novo</i> в <i>Y.l.</i> (2026); в <i>E. coli</i> — 7,2–12,16 г/л [42] (промышленный штамм GeneIII [40])", st_td),
            ],
        ]
        ref_tbl = Table(ref_tbl_data, colWidths=[35 * mm, 29 * mm, 38 * mm, 35 * mm, 37 * mm], repeatRows=1)
        ref_tbl.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
                    ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
                    ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 3.5),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 3.5),
                    ("TOPPADDING", (0, 0), (-1, -1), 3.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
                ]
            )
        )
        story.append(Spacer(1, 1.5 * mm))
        story.append(ref_tbl)
        story.append(P("Таблица 2. Опубликованные экспериментальные референсы первого шага в дрожжах, подтверждающие гарантированный стартовый выход без сложного редактирования генома.", st_caption))

        story.append(H2("sec3_2", "3.2. Химическая схема пути (Шаг 1) и биохимические условия перехода к Шагу 2"))
        story.append(
            P(
                "Вставляемый двухгенный каскад проводит три последовательные реакции (Рисунок 2): (1) N-концевой метилтрансферазный домен фермента <b>NcEgt1</b> (из <i>Neurospora crassa</i>, GenBank XP_956324.3) переносит три метильные группы от трёх молекул SAM на α-аминогруппу L-гистидина с образованием бетаина <b>герцинина</b> (C₉H₁₅N₃O₂) и трёх молекул S-аденозил-L-гомоцистеина (SAH); (2) C-концевой Fe²⁺-зависимый сульфоксидсинтазный домен того же <b>NcEgt1</b> при участии растворённого кислорода O₂ пришивает атом серы L-цистеина ко 2-му углероду имидазольного кольца с образованием <b>сульфоксида Cys-HER</b> (C₁₂H₂₀N₄O₅S); (3) пиридоксальфосфат-(PLP)-зависимая лиаза <b>CpEgt2</b> (из <i>Claviceps purpurea</i>) отщепляет пируват и аммиак, формируя стабильный тионный таутомер <b>L-эрготионеина</b> (C₉H₁₅N₃O₂S).",
                st_body,
            )
        )
        story.append(
            P(
                "<b>Критическое биохимическое уточнение по переходу «Шаг 1 → Шаг 2 (SAMe)»:</b> L-эрготионеин не является химическим интермедиатом адеметионина, а при переходе к промышленному накоплению SAMe переносятся <b>сам программно-лабораторный конвейер DBTL и энергетический модуль регенерации АТФ / синтазы <i>SAM2</i></b>, тогда как две периферийные ветви подлежат целенаправленной инверсии. Во-первых, на Шаге 1 усиливается транссульфурирование (<i>↑CYS3/CYS4</i>) для снабжения <i>NcEgt1</i> цистеином, тогда как в продуцентах SAMe ген цистатионин-β-синтазы <i>CYS4</i>, напротив, ослабляют (<i>↓CYS4</i>), чтобы предотвратить утечку гомоцистеина в цистеин и глутатион. Во-вторых, разгон <i>↑HIS1</i> на Шаге 1 расходует пул PRPP, который на Шаге 2 необходим для синтеза аденозинового кольца самого SAMe. Наконец, промышленный синтез SAMe ведут не <i>de novo</i> из сульфата (сжигающего 2 АТФ и 4 NADPH на моль серы), а с подпиткой дешёвым L/DL-метионином ($2,0–2,8/кг).",
                st_body,
            )
        )

        img2 = Image(str(fig2_path), width=usable_w, height=usable_w * (7.6 / 11.2))
        story.append(
            KeepTogether(
                [
                    img2,
                    P(
                        "Рисунок 2. Структурные химические формулы интермедиатов биосинтеза L-эрготионеина в дрожжах (Шаг 1) и схема перехода к продуктам Шага 2 с указанием переносимого блока [SAM2 + регенерация АТФ] и реконфигурируемых узлов (↓CYS4 и подпитка метионином для SAMe).",
                        st_caption,
                    ),
                ]
            )
        )

        story.append(H2("sec3_3", "3.3. Топология метаболического узла и архитектура моделирования Design–Learn"))
        story.append(
            P(
                "В отличие от линейных путей (таких как каротиноидный), сборка одной молекулы L-эрготионеина требует одновременного сведения в строгой пропорции <b>трёх конкурирующих за ресурсы потоков клетки</b>: 1 моль L-гистидина + 1 моль L-цистеина + 3 моль SAM (на регенерацию которых тратится ~9–10 молекул АТФ). Ручная и модульная инженерия без ML (включая выход на 2,39 г/л в DTU [28] и 3,05 г/л в BUCT в апреле 2026 г. [41]) неизбежно упирается в нелинейные перекрёстные ограничения (эпистаз), делающие этот узел показательным полигоном для активного обучения:",
                st_body,
            )
        )
        model_bullets = [
            "<b>Конкуренция за АТФ и PRPP:</b> снятие обратной связи по гистидину (мутантная форма фермента <i>HIS1</i>, не тормозящаяся избытком гистидина) перетягивает общий предшественник PRPP от синтеза адениновых нуклеотидов, снижая внутриклеточный уровень АТФ и тормозя образование SAM.",
            "<b>Ингибирование побочным продуктом SAH:</b> ускорение метилирования гистидина без одновременного усиления фермента гидролиза S-аденозилгомоцистеина (<i>SAH1</i>) и регенерации метионина (<i>MET6</i>) приводит к накоплению SAH, который блокирует активный центр <i>NcEgt1</i>.",
            "<b>Сульфидная токсичность:</b> чрезмерное усиление восстановления сульфата (<i>MET3/MET14/MET16</i>) при отставании <i>NcEgt1</i> вызывает сброс избытка серы в сероводород (H₂S), угнетающий дыхание дрожжей [28]. Именно поэтому рекорд 9,3 г/л в <i>Y. lipolytica</i> [30] требовал дорогой подпитки 20 г/л готового L-метионина, а <i>de novo</i> выход в <i>S. cerevisiae</i> остановился на 2,39 г/л [28].",
        ]
        for b in model_bullets:
            story.append(P("• " + b, st_bullet))

        story.append(
            P(
                "Для нахождения сбалансированного режима вычислительный контур DBTL объединяет два уровня (Рисунок 3): (1) геномную стехиометрическую модель дрожжей с ферментными ограничениями (<i>ecYeast8 / Yeast9</i>) с добавленными реакциями <i>NcEgt1/CpEgt2</i> для расчёта допустимых потоков; (2) вероятностную суррогатную модель машинного обучения (Gaussian Process / XGBoost в контуре байесовского активного обучения по схеме ART [5]), предсказывающую объёмную производительность (мг/л·ч) по вектору силы промоторов 6–8 генов (<i>NcEgt1, CpEgt2, HIS1, MET14, SAM2, SAH1, Δspe2</i> и экспортёров <i>AQR1/TPO1</i>) и концентрациям кофакторов среды (Fe²⁺, SO₄²⁻, пиридоксин).",
                st_body,
            )
        )

        img3 = Image(str(fig3_path), width=usable_w, height=usable_w * (3.25 / 10.6))
        story.append(
            KeepTogether(
                [
                    img3,
                    P("Рисунок 3. Схема работы замкнутого контура Design–Build–Test–Learn на задаче оптимизации биосинтеза L-эрготионеина.", st_caption),
                ]
            )
        )

        # ==================== SECTION 4 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec4", "4. Пошаговый лабораторный протокол и план экспериментов"))
        story.append(
            P(
                "Ниже приведён пошаговый рабочий регламент нулевого и первого циклов на пекарских дрожжах <i>Saccharomyces cerevisiae</i>, опирающийся на проверенные методики группы И. Бородиной [27, 28], протокол <i>HBUT-China</i> [32] и селективные спектрофотометрические тесты [33, 34]. Все процедуры выполняются на базовом лабораторном оборудовании (термостат-шейкер, ПЦР-амплификатор, настольная центрифуга, обычный кюветный или планшетный спектрофотометр).",
                st_body,
            )
        )

        story.append(H2("sec4_1", "4.1. Приготовление прозрачной минеральной среды Делфта (на 1 литр)"))
        story.append(
            P(
                "Вместо обычной богатой среды YPD с дрожжевым экстрактом и пептоном используется прозрачная минеральная среда Делфта по двум практическим причинам: (1) дрожжевой экстракт уже содержит готовые аминокислоты (гистидин, метионин, цистеин), что мешает оценить способность штамма синтезировать продукт из сахара; (2) жёлто-коричневые компоненты пептона и дрожжевого экстракта сильно поглощают ультрафиолет при 310–350 нм и забивают сигнал продукта на спектрофотометре. Минеральная среда полностью бесцветна и прозрачна в УФ.",
                st_body,
            )
        )
        media_steps = [
            "<b>Шаг 1 (солевая основа):</b> растворить в 850 мл деионизированной воды <b>7,5 г (NH₄)₂SO₄</b> (источник азота и сульфатной серы), <b>14,4 г KH₂PO₄</b> и <b>0,5 г MgSO₄·7H₂O</b>. Довести pH до <b>6,0</b> раствором 2 М KOH, довести объём водой до 900 мл и проавтоклавировать (121 °C, 20 мин).",
            "<b>Шаг 2 (раствор глюкозы):</b> отдельно растворить 200 г D-глюкозы в 500 мл воды (40%-й сток) и простерилизовать автоклавированием или фильтрацией через мембрану 0,22 мкм.",
            "<b>Шаг 3 (сборка готовой среды после остывания до &lt;40 °C):</b> в стерильных условиях прилить к 900 мл солевой основы <b>50 мл</b> 40%-й глюкозы (конечная концентрация <b>20 г/л</b>), <b>2,0 мл</b> стерильного раствора микроэлементов 500× (обязательно содержащего <b>FeSO₄·7H₂O 3,0 г/л в стоке</b> для работы активного центра Fe²⁺ в <i>NcEgt1</i>), <b>1,0 мл</b> фильтрованного раствора витаминов 1000× (с увеличенной до <b>10 г/л в стоке</b> навеской пиридоксина HCl — витамина B₆, из которого клетка образует кофактор PLP для <i>CpEgt2</i>) и довести стерильной водой до 1 000 мл [28].",
        ]
        for s in media_steps:
            story.append(P("• " + s, st_bullet))

        story.append(H2("sec4_2", "4.2. Сборка генетической кассеты и трансформация дрожжей (Build)"))
        build_steps = [
            "<b>Шаг 1 (выбор реципиентного штамма):</b> используется стандартный лабораторный штамм <i>S. cerevisiae</i> <b>CEN.PK113-7D</b> (прототроф, растущий на минимальной среде без добавок) либо его производное <b>CEN.PK113-5D</b> (<i>ura3-52</i>, нуждающееся в урациле до вставки плазмиды с маркером <i>URA3</i>) [27, 32].",
            "<b>Шаг 2 (синтез и сборка ДНК):</b> заказываются оптимизированные под кодоны дрожжей последовательности генов <b>NcEgt1</b> (2 631 п.н., GenBank XP_956324.3) и <b>CpEgt2</b> (~1 400 п.н.). Методом Golden Gate, Гибсона или USER собирается двунаправленная кассета с сильными конститутивными промоторами дрожжей <i>pTEF1–NcEgt1–tADH1 + pPGK1–CpEgt2–tCYC1</i> внутри интеграционной плазмиды набора <i>EasyClone-MarkerFree</i> с плечами гомологии по 500 п.н. к проверенным межгенным участкам хромосом (сайты <i>XI-3</i> и <i>X-4</i> либо локус <i>PDC1–STU2</i>). Перед трансформацией 2 мкг плазмиды разрезают рестриктазой <i>NotI</i> (37 °C, 1 ч) для получения линейного фрагмента вставки.",
            "<b>Шаг 3 (трансформация ацетатом лития):</b> ночную культуру дрожжей подращивают в 25 мл YPD до оптической плотности OD₆₀₀ = 0,8–1,0 (30 °C, 3,5 ч), осаждают клетки (3 000 × g, 5 мин) и промывают водой. К 50 мкл осадка клеток последовательно добавляют: <b>240 мкл</b> 50%-го ПЭГ-3350, <b>36 мкл</b> 1,0 М ацетата лития (LiAc), <b>10 мкл</b> прокипячённой одноцепочечной ДНК-носителя (10 мг/мл), <b>500 нг</b> плазмиды Cas9/gRNA и <b>1 000 нг</b> линейного фрагмента кассеты; доводят водой до 360 мкл, перемешивают и инкубируют на водяной бане при <b>42 °C в течение 35 мин</b>.",
            "<b>Шаг 4 (отбор и проверка клонов):</b> клетки осаждают, выдерживают 2 ч в 1 мл свежей среды при 30 °C и высевают на селективные чашки (YPD с антибиотиками G418 200 мг/л и нурсеотрицином 100 мг/л либо агар без урацила SC-Ura) на 48 ч при 30 °C. Вставку подтверждают ПЦР с колоний по краям локуса, после чего сбрасывают временную плазмиду Cas9 пересевом без селекции (или на чашках с 1 г/л 5-фтороротовой кислоты). Для Цикла 1 аналогичным образом собирают панель из 24–48 штаммов с заменой промоторов (<i>pRPL18B, pRAD27, pACT1, pTEF1, pTDH3</i>) перед генами <i>HIS1, MET14, SAM2, SAH1</i>, делецией <i>Δspe2</i> и экспортёрами <i>AQR1/TPO1</i>.",
        ]
        for s in build_steps:
            story.append(P("• " + s, st_bullet))

        story.append(H2("sec4_3", "4.3. Выращивание культур и водная термоэкстракция при 94 °C"))
        cult_steps = [
            "<b>Шаг 1 (посев и ферментация 72 ч):</b> одиночную колонию проверенного штамма (в 3 повторах + пустой родительский штамм CEN.PK как холостой контроль) засевают в 5 мл минеральной среды Делфта на 24 ч (30 °C, 250 об/мин). Затем измеряют мутность культуры при 600 нм (OD₆₀₀) и пересевают в свежую минеральную среду Делфта до стартовой плотности <b>OD₆₀₀ = 0,1</b> (по 5–25 мл в пробирках/колбах либо по 500 мкл в глубоком 96-луночном планшете). Выращивают ровно <b>72 ч при 30 °C и 250 об/мин</b>.",
            "<b>Шаг 2 (отделение культуральной жидкости от клеток):</b> по окончании 72 ч записывают финальную плотность OD₆₀₀, отбирают <b>1,0 мл</b> суспензии в микропробирку и центрифугируют при <b>3 000 × g в течение 5 мин</b>. Верхнюю прозрачную жидкость над осадком (супернатант) аккуратно переносят пипеткой в чистую пробирку — это <i>внеклеточная проба</i> (то, что дрожжи выделили наружу).",
            "<b>Шаг 3 (извлечение продукта изнутри клеток кипячением в воде):</b> клеточный осадок на дне однократно промывают 1 мл воды, затем приливают к нему ровно <b>1,0 мл чистой дистиллированной воды</b>, взбалтывают и ставят закрытую пробирку в термоблок при <b>94 °C на 10 мин</b> (при нагреве мембраны открываются, белки денатурируют, а термостойкий эрготионеин полностью переходит в воду [27, 28]). Пробирку встряхивают ещё 20 мин, центрифугируют при <b>10 000 × g (5 мин)</b> и отбирают прозрачный водный экстракт — это <i>внутриклеточная проба</i>.",
            "<b>Важное правило безопасности методики:</b> для осаждения белков категорически нельзя добавлять трихлоруксусную кислоту (ТХУ) — она за несколько минут окисляет и уничтожает более 90% эрготионеина [33].",
        ]
        for s in cult_steps:
            story.append(P("• " + s, st_bullet))

        story.append(H2("sec4_4", "4.4. Быстрое измерение концентрации на спектрофотометре (311 и 343 нм) и ВЭЖХ"))
        story.append(
            P(
                "Оба приведённых ниже метода выполняются при комнатной температуре за 5–10 минут в стандартной УФ-прозрачной кювете объёмом 1 мл (длина оптического пути 1 см) на любом обычном спектрофотометре; при работе в 96-луночном УФ-планшете все объёмы просто уменьшаются в 5 раз (до суммы 200 мкл):",
                st_body,
            )
        )
        assay_steps = [
            "<b>Метод А (основной, абсолютная избирательность, λ = 311 нм; требует предварительной наработки фермента):</b> фермент эрготионаза <i>TdETL</i> не продаётся готовым реагентом, поэтому на подготовительном этапе (2–3 нед. до Цикла 0) его ген клонируют в <i>pET28a</i>, экспрессируют в <i>E. coli</i> BL21(DE3) и очищают на колонке Ni-NTA до запаса на тысячи анализов [34]. Для замера в кювету наливают <b>900 мкл</b> 50 мМ калий-фосфатного буфера (pH 8,0), вносят <b>100 мкл</b> пробы, считывают фон A₀ при 311 нм, добавляют <b>5 мкл</b> <i>TdETL</i> (до 400 нМ) и через 10 мин измеряют пик тиоурокановой кислоты A₁ (ε₃₁₁ = 22 500 М⁻¹·см⁻¹; цистеин, глутатион, гистидин и герцинин сигнала не дают): <b>C (мг/л) = 101,9 × ΔA₃₁₁</b>.",
            "<b>Метод Б (без наработки фермента, реактив Карлссона 2-PDS, λ = 343 нм — оптимален для внеклеточного супернатанта):</b> смешивают <b>500 мкл</b> прозрачной культуральной жидкости и <b>500 мкл</b> 0,5 мМ раствора 2,2'-дипиридилдисульфида (2-PDS) в 0,25 М HCl (итоговый <b>pH ≈ 1,0</b>) и через 5 мин измеряют пик 2-тиопиридона при <b>343 нм</b> [33]. Важное практическое ограничение: во внеклеточной среде Делфта тест работает чисто, но внутри клеток <i>S. cerevisiae</i> содержится 5–10 мМ (1,5–3 г/л) глутатиона (20–50-кратный молярный избыток над EGT на ранних циклах), поэтому внутриклеточные экстракты требуют как пред-окисления 0,1 мМ CuSO₄ (pH 9,5, 15 мин), так и обязательной калибровочной сверки с ВЭЖХ-УФ.",
            "<b>Контрольная проверка топ-15% штаммов и калибровка методом ВЭЖХ-УФ:</b> лучшие штаммы каждого раунда и калибровочную подвыборку проверяют на обычном хроматографе с УФ-детектором на колонке C18 (0,1% муравьиная кислота / ацетонитрил, 9 мин) на двух длинах волн — <b>254 нм</b> (пик L-эрготионеина) и <b>220 нм</b> (пик герцинина) [27, 28]. Отношение [Герцинин]/[EGT] показывает модели, какой из двух ферментов (<i>NcEgt1</i> или <i>CpEgt2</i>) лимитирует поток.",
        ]
        for s in assay_steps:
            story.append(P("• " + s, st_bullet))

        roadmap_data = [
            [
                P("Этап / Цикл", st_th),
                P("Срок", st_th),
                P("Лабораторные работы (Build + Test)", st_th),
                P("Целевой результат и выход для модели Learn", st_th),
            ],
            [
                P("<b>Цикл 0</b><br/>Базовый референс", st_td_bold),
                P("Нед. 1–4", st_td),
                P("Вставка 1 и 2 копий <i>NcEgt1 + CpEgt2</i> в <i>S. cerevisiae</i> CEN.PK; проверка экстракции при 94 °C и калибровка спектрофотометра (311 / 343 нм) против ВЭЖХ.", st_td),
                P("Воспроизводимый выход 50–150 мг/л на минеральной среде; линейность калибровки спектрофотометра R² &gt; 0,95.", st_td),
            ],
            [
                P("<b>Цикл 1</b><br/>Обучающая выборка", st_td_bold),
                P("Нед. 5–9", st_td),
                P("Сборка 36–48 штаммов по дробному факторному плану (разная сила промоторов перед <i>HIS1, MET14, SAM2, SAH1</i> + <i>Δspe2</i> + экспортёры <i>AQR1/TPO1</i>) × 3 варианта среды.", st_td),
                P("Размеченный массив &gt;350–400 измерений (включая слабые и нулевые комбинации); первичное обучение суррогатной модели.", st_td),
            ],
            [
                P("<b>Циклы 2–3</b><br/>ML-оптимизация", st_td_bold),
                P("Нед. 10–18", st_td),
                P("Два раунда активного обучения: балансировка регенерации 3×SAM и секреции наружу без добавки аминокислот; проверка лидеров в 1–5 л биореакторе.", st_td),
                P("Выход на &gt;2,0–3,0+ г/л <i>de novo</i> на глюкозе и сульфате при доле секреции в среду &gt;80%; заявка на патент РФ и методологическая статья Q1.", st_td),
            ],
            [
                P("<b>Переход в Шаг 2</b><br/>Косм. пилот + SAMe", st_td_bold),
                P("Мес. 6+", st_td),
                P("Наработка 10–20 кг 1%-го водного концентрата EGT (ультрафильтрация 3 кДа + уголь) под R&amp;D-коллаборацию с брендом РФ; реконфигурация штамма (<i>↓CYS4</i> + подача Met) под <b>(S,S)-SAMe</b>.", st_td),
                P("Индустриальный партнёр, хоздоговор и акт внедрения на Шаге 1 (серия ~15 тыс. флаконов сыворотки) + регламент глубокой ферментации SAMe.", st_td),
            ],
        ]
        rm_tbl = Table(roadmap_data, colWidths=[28 * mm, 18 * mm, 66 * mm, 62 * mm], repeatRows=1)
        rm_tbl.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
                    ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
                    ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ("TOPPADDING", (0, 0), (-1, -1), 4),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ]
            )
        )
        story.append(
            KeepTogether(
                [
                    H2("sec4_5", "4.5. Календарный план итераций Циклов 0–3 и перехода к Шагу 2"),
                    rm_tbl,
                    P("Таблица 3. Поэтапный план проведения циклов DBTL и масштабирования метаболического узла.", st_caption),
                ]
            )
        )

        # ==================== SECTION 5 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec5", "5. Экономика, регуляторика и критический стресс-тест рынков мира и РФ"))

        story.append(H2("sec5_1", "5.1. Мировой рынок L-эрготионеина, профицит мощностей в КНР и ценовая реальность"))
        mkt_world = [
            "<b>Физический объём мирового рынка против промо-прогнозов:</b> по физическому балансу мировое потребление субстанции L-эрготионеина в 2025 г. составило всего <b>~68 метрических тонн (~$75 млн)</b> с реалистичным темпом роста <b>~18% в год</b> (против рекламных оценок CAGR 28–38% у маркетинговых агрегаторов; Рисунок 4А). Интерес косметических концернов (<i>Estée Lauder, L'Oréal, Shiseido, Proya</i>) обусловлен наличием у человека специфического переносчика <b>OCTN1</b> (<i>SLC22A4</i>) и устойчивостью тионной формы (C=S) к автоокислению при pH 7,4, а низкое содержание в грибах (0,1–2 мг/г) и рацемизация при химсинтезе закрепили переход отрасли на ферментацию.",
            "<b>Навес производственных мощностей в КНР и обвал оптовых цен:</b> за 2024–2026 гг. китайские биотех-компании ввели мощности, кратно перекрывающие текущий мировой спрос. Лидер рынка <i>Jiangsu GeneIII Biotech</i> запустил ферментёры по 30 м³ с выпуском <b>3–5 тонн/мес (36–60 т/год)</b> порошка чистотой 99,99% и в январе 2026 г. получил статус <b>FDA GRAS (GRN 001270)</b> на штамме <i>E. coli</i> K-12 MG1655 [40]. Параллельно введены и строятся линии <i>Bloomage Biotech</i> (1 т/год), <i>Yili Chuanning</i> (0,5 т/год; компания отчиталась о чистом убытке своего синбио-дивизиона в 92,8 млн юаней за 2025 г. из-за падения цен), <i>Baikaisheng</i> (10 т) и <i>Zhuhai United Lab</i> (20 т). В результате оптовая цена 98–99% порошка FOB КНР к осени 2026 г. снизилась до <b>$400–580/кг</b>.",
            "<b>Ограничение гипотезы «высокомаржинального фильтрата»:</b> жидкий косметический актив <i>Bloomage Bioyouth™-EGT</i> (соферментат мицелия грибов мацутакэ и ежовика) содержит лишь <b>0,013–0,05% (130–500 мг/кг)</b> чистого эрготионеина и продаётся малыми канистрами как брендированный грибной экстракт (вся косметическая линейка EGT принесла <i>Bloomage</i> лишь чуть более 1 млн юаней выручки). При доступности на рынке бесцветного 99% кристаллического порошка без запаха по $400–580/кг рассчитывать на массовую продажу небрендированного бульона пекарских дрожжей по пересчётной цене $5 000–10 000/кг актива нереалистично.",
        ]
        for b in mkt_world:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec5_2", "5.2. Шаг 1 в РФ (микрониша и механика пилота с брендом) и масштаб Шага 2 (SAMe)"))
        mkt_ru = [
            "<b>Шаг 1 в России — патентная чистота дрожжей, но малый физический тоннаж для отдельного завода:</b> в реестре ФИПС нет патентов РФ/ЕАПО на рекомбинантные дрожжи с генами <i>Egt1/Egt2</i> — действующие иностранные патенты закрывают только природный мицелий <i>Hericium erinaceus</i> HT-3 (<b>RU 2788417 C1</b>, <i>Bloomage</i>, 2023) и бромный химсинтез (<b>RU 2550202 C2</b>, <i>Tetrahedron</i>) [39]. Однако при рабочей дозировке чистого EGT в косметике <b>0,01–0,05%</b> (1 кг порошка хватает на 33 000–66 000 флаконов сыворотки по 30 мл) физический объём рынка РФ невелик: сырьё возят четыре нишевых поставщика — <b>ООО «БИОХИМ»</b> (СПб, декларация ЕАЭС июля 2025 г. на <i>ErgoX™</i> от <i>Effepharm</i>), <b>ООО «Лотос Трейдинг»</b> (<i>Bioyouth™-EGT</i>), <b>Witec</b> (складская фасовка по <b>100 г</b>) и <b>GEON</b> (от 5 кг), а среди конечных линеек подтверждены <b>«Ангиофарм»</b> (<i>Antioxidant Complex Serum</i>, IV кв. 2026 г.), <b>Alga Spa</b>, <b>A.M.&amp;P.M.</b> и грибные капсулы <b>FUNGILINE</b>. Весь спрос РФ составляет <b>единицы–десятки кг/год</b>, поэтому ценовая конкуренция с китайским оптом ($400–580/кг из реакторов 30 м³) не окупит капитальный цех.",
            "<b>Как превратить Шаг 1 в раннего индустриального партнёра (R&amp;D-коллаборация с косметическим брендом РФ):</b> вместо биржевой конкуренции с Китаем по цене сухого порошка малый физический расход EGT становится преимуществом для запуска <b>совместного спецпроекта («первый российский биотехнологический актив, созданный на ИИ-платформе DBTL»)</b> с наукоёмкими отечественными марками (потенциальные партнёры — <i>«Гельтек», «Ангиофарм», Openface, «Краснополянская косметика», Art&amp;Fact</i>). Для промышленной серии в <b>15 000 флаконов сыворотки по 30 мл</b> (при вводе 0,03% EGT) заводу требуется всего <b>135 г чистого эрготионеина</b> (= 13,5 кг 1%-го водного концентрата, нарабатываемого за 2–3 цикла в обычном 10-л лабораторном биореакторе). Порог входа минимален и по очистке, и по документам: (1) <b>настольная 4-стадийная очистка без хроматографа</b> — микрофильтрация 0,22 мкм (100% удаление ГМ-клеток и ДНК) → ультрафильтрация через мембрану <b>3 кДа</b> (отсечение всех белков и маннанов дрожжей) → осветление 1–2% фарм-углём и смешанной смолой (удаление запаха и солей среды Делфта) → стабилизация 5% косметического <i>Pentylene Glycol</i> и стерилизующая фильтрация 0,22 мкм (чистый состав INCI: <i>Aqua, Pentylene Glycol, Ergothioneine</i>); (2) <b>быстрый регуляторный допуск</b> — само B2B-сырьё не требует СГР или декларации по ТР ТС 009/2011 (декларируется уже готовая сыворотка на стороне бренда за 2–3 нед.), а полный пакет поставщика для входного контроля завода по ГОСТ ISO 22716 (внутренние ТУ + Паспорт безопасности ПБХП по ГОСТ 30333 + сертификат анализа партии и протокол аккредитованной лаборатории на микробиологию, тяжёлые металлы и отсутствие раздражающего действия <i>in vitro</i>) укладывается в <b>~70–130 тыс. руб. и 1–1,5 мес.</b> В формате эксклюзивного R&amp;D-договора (0,3–1,5 млн руб.) это даёт проекту реального индустриального заказчика, акт внедрения и серийный продукт на полке уже по итогам Шага 1.",
            "<b>Шаг 2 в России — реальный фармацевтический и ветеринарный рынок SAMe на 14–16 млрд руб./год:</b> в отличие от эрготионеина, <b>(S,S)-адеметионин (SAMe)</b> формирует более 47% всего рынка гепатопротекторов РФ (>13 млрд руб. в аптеках за 11 мес. 2025 г. по данным <i>RNC Pharma</i> + госпитальный сегмент ЖНВЛП и ветеринарные гепатопротекторы; Рисунок 4Б) [37]. Сверка карточек <b>ГРЛС Минздрава РФ (ноябрь 2026 г.)</b> [43] подтверждает доминирование зарубежных дрожжевых субстанций: у «Гептрала» (Abbott) — <i>Gnosis Bioresearch</i> (Швейцария/Италия), у российских дженериков («Адеметионин Велфарм», «Самеликс», «Гептиэль») — китайские заводы <i>Zhejiang Hisun</i>, <i>Wuxi Fortune</i>, <i>Shandong Jincheng</i> и <i>Shandong Ruiying</i>. В 2026 г. в ГРЛС появились первые российские записи на стадию АФС (<b>ООО «БЕРАХИМ»</b>, г. Обнинск, ФС-002473 от 01.07.2026, и площадка <b>ООО «Велфарм»</b> в Долгопрудном), что делает <b>собственный высокопродуктивный штамм глубокой ферментации на отечественном сырье</b> ключевым технологическим преимуществом для правила «второй лишний» в торгах ЖНВЛП.",
            "<b>Точная привязка к государственным перечням РФ:</b> в <b>Приказе Минсельхоза РФ № 392</b> [36] п. 21 (<i>«Бетаин, его соли и соединения»</i>) ориентирован на крупнотоннажный кормовой триметилглицин и не применим для субсидирования косметического эрготионеина, а п. 24 (<i>«Метионин, его соли и соединения»</i>) релевантен для кормовых и ветеринарных форм SAMe; основной же объём господдержки и спроса по SAMe относится к стратегии Минпромторга <b>«Фарма-2030» (перечень ЖНВЛП)</b> и нацпроекту «Технологическое обеспечение биоэкономики».",
        ]
        for b in mkt_ru:
            story.append(P("• " + b, st_bullet))

        img4 = Image(str(fig4_path), width=usable_w, height=usable_w * (3.95 / 10.6))
        story.append(
            KeepTogether(
                [
                    Spacer(1, 1.5 * mm),
                    img4,
                    P(
                        "Рисунок 4. Сравнение масштаба рынков двух шагов: (А) базовая физическая оценка мирового рынка L-эрготионеина (~68 т, $75 млн, CAGR 18%) на фоне промо-прогнозов и ликопина; (Б) аптечные продажи гепатопротекторов в РФ за 11 мес. 2025 г. (RNC Pharma, доля SAMe >47%).",
                        st_caption,
                    ),
                ]
            )
        )

        story.append(H2("sec5_3", "5.3. Критический стресс-тест концепции и развилка траекторий (Сценарий А vs Сценарий Б)"))
        story.append(
            P(
                "Сквозная проверка всех 16 кандидатов выявила фундаментальную развилку: ни одна молекула не сочетает одновременно фоново-чистый кюветный тест без хроматографии и незанятый многомиллиардный рынок Шага 1 в РФ. В зависимости от приборного оснащения лаборатории и приоритета первых 6 месяцев проект может стартовать по одному из двух честных сценариев (Таблица 4):",
                st_body,
            )
        )

        scen_data = [
            [
                P("Параметр выбора", st_th),
                P("<b>Сценарий А (базовый двухшаговый):</b><br/>Модельный старт на <b>L-эрготионеине</b> → реконфигурация под <b>SAMe</b>", st_th),
                P("<b>Сценарий Б (прямой индустриальный):</b><br/>Старт цикла DBTL сразу на дрожжевом <b>(S,S)-адеметионине (SAMe)</b>", st_th),
            ],
            [
                P("<b>Условие выбора</b>", st_td_bold),
                P("Приоритет первых 4–6 мес. — отладить с нуля конвейер DBTL силами малой/учебной группы без очередей на ВЭЖХ и выпустить статью Q1 по ML.", st_td),
                P("В лаборатории доступен обычный ВЭЖХ-УФ с автосамплером; приоритет — с первого дня работать на крупный рынок АФС и ветпрепаратов РФ.", st_td),
            ],
            [
                P("<b>Старт (Цикл 0) и стадия Test</b>", st_td_bold),
                P("Вставка 2 генов (<i>NcEgt1, CpEgt2</i>), сразу 46–84 мг/л <i>de novo</i>. Спектрофотометр 311 / 343 нм в супернатанте + ВЭЖХ для топ-15%.", st_td),
                P("<b>0 чужеродных генов на старте</b> (нативный <i>SAM2</i> в <i>S. cerevisiae / K. phaffii</i> + подача дешёвого L-метионина даёт 0,3–0,8 г/л). Быстрый изократический ВЭЖХ-УФ 254–260 нм (3–4 мин/проба).", st_td),
            ],
            [
                P("<b>Задача для ML (Design–Learn)</b>", st_td_bold),
                P("Балансировка 3 конкурирующих ветвей <i>de novo</i> (<i>HIS + CYS + 3 SAM</i>); <b>0 статей по ML-DBTL в мире</b>.", st_td),
                P("Многофакторный баланс регенерации АТФ, экспрессии <i>SAM2</i>, блока слива <i>↓CYS4</i>, полифосфатного депо вакуоли (<i>VTC4</i>) и подпитки Met.", st_td),
            ],
            [
                P("<b>Экономика в РФ и главный риск</b>", st_td_bold),
                P("Завод на Шаге 1 не окупится, но <b>реален R&amp;D-пилот с брендом РФ</b> (13,5 кг 1% р-ра на 15 тыс. сывороток; допуск ~100 тыс. ₽); для SAMe нужен разворот <i>↓CYS4</i>.", st_td),
                P("Прямая работа на рынок <b>14–16 млрд ₽/год</b> (ЖНВЛП «Фарма-2030», п. 24 МСХ № 392) без потери 6 мес.; нагрузка 3–6 ч ВЭЖХ на раунд.", st_td),
            ],
        ]
        scen_tbl = Table(scen_data, colWidths=[30 * mm, 72 * mm, 72 * mm], repeatRows=1)
        scen_tbl.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
                    ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
                    ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ("TOPPADDING", (0, 0), (-1, -1), 3.8),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 3.8),
                ]
            )
        )
        story.append(scen_tbl)
        story.append(P("Таблица 4. Сравнение двух рабочих сценариев запуска проекта DBTL по результатам критического фактчека.", st_caption))

        # ==================== SECTION 6: REFERENCES (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec6", "6. Список литературы и источников"))
        refs = [
            "[1] Nielsen J., Keasling J. D. Engineering cellular metabolism. <i>Cell</i>. 2016;164(6):1185–1197. doi:10.1016/j.cell.2016.02.004.",
            "[2] Lawson C. E., Martí J. M., Radivojevic T. et al. Machine learning for metabolic engineering: a review. <i>Metabolic Engineering</i>. 2021;63:34–60. doi:10.1016/j.ymben.2020.10.005.",
            "[3] Carbonell P., Jervis A. J., Robinson C. J. et al. An automated Design-Build-Test-Learn pipeline for enhanced microbial production of fine chemicals. <i>Communications Biology</i>. 2018;1:66. doi:10.1038/s42003-018-0076-9.",
            "[4] HamediRad M., Chao R., Weisberg S. et al. Towards a fully automated algorithm driven platform for biosystems design (BioAutomata, lycopene). <i>Nature Communications</i>. 2019;10:5150. doi:10.1038/s41467-019-13189-z.",
            "[5] Radivojević T., Costello Z., Workman K., Garcia Martin H. A machine learning Automated Recommendation Tool for synthetic biology. <i>Nature Communications</i>. 2020;11:4879. doi:10.1038/s41467-020-18008-4.",
            "[6] Zhou Y., Li G., Dong J., Xing X., Dai J., Zhang C. MiYA, an efficient machine-learning workflow in conjunction with the YeastFab assembly strategy for combinatorial optimization of heterologous metabolic pathways. <i>Metabolic Engineering</i>. 2018;47:294–302. doi:10.1016/j.ymben.2018.03.020.",
            "[7] Watcharawipas A., Runguphan W. et al. Integrating yeast biodiversity and machine learning for predictive metabolic engineering. <i>FEMS Yeast Research</i>. 2025;25:foaf072. doi:10.1093/femsyr/foaf072.",
            "[8] Команды iGEM NEFU-China (2022, ликопин в <i>S. cerevisiae</i>, 2022.igem.wiki/nefu-china) и NACIS Shanghai (2025, бета-каротин, 2025.igem.wiki/nacis-shanghai).",
            "[9] Ma Y., Liu N., Greisen P. et al. Removal of lycopene substrate inhibition enables high carotenoid productivity in <i>Yarrowia lipolytica</i> (17.6 g/L lycopene, 39.5 g/L β-carotene). <i>Nature Communications</i>. 2022;13:572. doi:10.1038/s41467-022-28277-w.",
            "[10] Shimazaki S. et al. Stacking machine-learning framework for combinatorial promoter optimization of carotenoid biosynthesis in <i>Saccharomyces cerevisiae</i>. 2024 (цит. по <i>FEMS Yeast Res.</i> 2025;25:foaf072).",
            "[11] Tu W., Xu J., Ma Y., Katsimpouras C., Stephanopoulos G. Machine learning-driven optimization of metabolic balance for β-carotene production in <i>Yarrowia lipolytica</i>. <i>Metabolic Engineering</i>. 2026;94:192–201. doi:10.1016/j.ymben.2025.12.002.",
            "[12] De Biaggi J. S., Park Y.-K., Kerkhoven E. J., Ledesma-Amaro R., Lahtvee P.-J. Enzyme-constrained genome-scale model of <i>Yarrowia lipolytica</i> (eciYali5-GEM) predicts growth-phase specific metabolic engineering targets. <i>Appl. Microbiol. Biotechnol.</i> 2026;110. doi:10.1007/s00253-026-13791-4.",
            "[13] Zhou L. et al. (Ledesma-Amaro lab). Engineering <i>Yarrowia lipolytica</i> for de novo crocetin production and extracellular secretion (30.17 mg/L crocetin vs 1575 mg/L trapped zeaxanthin). <i>Biotechnol. Biofuels Bioprod.</i> 2025;18.",
            "[14] Liang N. et al. Tunnel engineering of carotenoid cleavage dioxygenase CsCCD2 for crocetin production in <i>S. cerevisiae</i> (107 mg/L). <i>J. Agric. Food Chem.</i> 2021;69(39):11626–11636; Cárdenas-Conejo Y. et al. <i>BMC Genomics</i>. 2015;16:877 (ревизия генов биксина <i>Bixa orellana</i>); Lee J. et al. <i>Microb. Cell Fact.</i> 2024;23 (кроцины 15,8 мг/л в <i>E. coli</i>).",
            "[15] Abdullah C. N., Liu M., Zhou J. et al. Efficient production of astaxanthin in <i>Yarrowia lipolytica</i> (2.82 g/L). <i>Synthetic and Systems Biotechnology</i>. 2025;10(3):737–750. doi:10.1016/j.synbio.2025.02.014.",
            "[16] Moreno-Paz S. et al. Machine learning-guided metabolic engineering of <i>Saccharomyces cerevisiae</i> for <i>p</i>-coumaric acid production. <i>ACS Synthetic Biology</i>. 2024.",
            "[17] Rodriguez A., Borodina I., Nielsen J. et al. <i>Metabolic Engineering</i>. 2015;31:181–188; Liu Q., Yu T., Li X., Chen Y., Nielsen J., Chen Y. Rewiring carbon metabolism in yeast for high level production of aromatic chemicals (12.5 g/L <i>p</i>-coumaric acid). <i>Nature Communications</i>. 2019;10:4976. doi:10.1038/s41467-019-12961-5.",
            "[18] Sáez-Sáez J., Wang G., Marella E. R., Sudarsan S., Cernuda Pastor M., Borodina I. Engineering the oleaginous yeast <i>Yarrowia lipolytica</i> for high-level resveratrol production (12.4 g/L). <i>Metabolic Engineering</i>. 2020;62:51–61; Liu M. et al. <i>Bioresour. Technol.</i> 2022 (22.5 g/L); Deng Y. et al. <i>Nature Communications</i>. 2026 (глабридин 1,2 мг/л).",
            "[19] Chen S., Zheng H., Gao J., Song H., Bai W. High-level production of pullulan (140.2 g/L) and its biosynthesis regulation in <i>Aureobasidium pullulans</i> BL06. <i>Front. Bioeng. Biotechnol.</i> 2023;11:1131875. doi:10.3389/fbioe.2023.1131875.",
            "[20] Syrigou M., Tsouko E. Pullulan-based gels with food-related orientation (рынок ~$80 млн, опт $25–30/кг). <i>Gels</i>. 2026;12(7):631. doi:10.3390/gels12070631; West T. P. <i>J. Fungi</i>. 2026;12(2):84.",
            "[21] Perli T., Wronska A. K., Pronk J. T., Daran J.-M. Vitamin requirements and biosynthesis in <i>Saccharomyces cerevisiae</i>. <i>Yeast</i>. 2020;37(4):283–304; Ma D. et al. Advances and prospects in microbial production of biotin (механизм инактивации BioB). <i>Microbial Cell Factories</i>. 2024;23:135. doi:10.1186/s12934-024-02413-1.",
            "[22] Ihsan R. M., Cronan J. E. Bypass of the canonical biotin synthesis pathway in yeast (<i>CfBio1</i>). <i>PNAS</i>. 2026;123. doi:10.1073/pnas.2621848123.",
            "[23] Патент РФ RU 2769175 C1 (ФИЦ Биотехнологии РАН: Эльдаров М. А., Филькин С. Ю., Фёдоров А. Н. и др.). Способ микробиологической продукции химозина быка с использованием рекомбинантного штамма <i>Pichia pastoris</i> с коэкспрессией фактора HAC1; Филькин С. Ю., Липкин А. В., Фёдоров А. Н. <i>Food Systems</i>. 2024;7(2):46.",
            "[24] Moon S., Smanski M. J. Biosynthesis of tauro-ursodeoxycholic acid and taurine (2–3 mg/L) in <i>Saccharomyces cerevisiae</i>. <i>Microbial Cell Factories</i>. 2025;24:218. doi:10.1186/s12934-025-02829-3.",
            "[25] Cui H., Liang Y., Dai Z., Zhang L. et al. Direct synthesis of taurine from CO₂ via a tandem electro-biosystem in <i>S. cerevisiae</i> (38.17 mg/L). <i>Chemical Engineering Journal</i>. 2026. doi:10.1016/j.cej.2026.160821.",
            "[26] Li L., Li C., Chen Y., Sun Z., Qu G. Engineering a hybrid chemical-biological system for de novo taurine production in <i>E. coli</i> (4.08 g/L). 2026. PMID:41865819.",
            "[27] van der Hoek S. A., Darbani B., Zugaj K. E., Kell D. B., Borodina I. Engineering the yeast <i>Saccharomyces cerevisiae</i> for the production of L-(+)-ergothioneine. <i>Front. Bioeng. Biotechnol.</i> 2019;7:262. doi:10.3389/fbioe.2019.00262.",
            "[28] van der Hoek S. A., Rusnák M., Wang G., Jessop-Fabre M. M., Kell D. B., Borodina I. Engineering precursor supply for the high-level production of ergothioneine in <i>Saccharomyces cerevisiae</i> (2.39 g/L de novo). <i>Metabolic Engineering</i>. 2022;70:129–142. doi:10.1016/j.ymben.2022.01.012.",
            "[29] van der Hoek S. A., Rusnák M., Jacobsen I. H., Martínez J. L., Kell D. B., Borodina I. Engineering ergothioneine production in <i>Yarrowia lipolytica</i> (1.63 g/L). <i>FEBS Letters</i>. 2022;596(10):1356–1364. doi:10.1002/1873-3468.14239.",
            "[30] Hu L., Liu M., Chen Q., Zeng W., Zhou J., Zhang C., Xu S. Enhanced production of ergothioneine in <i>Yarrowia lipolytica</i> via enzyme engineering and precursor perturbation. <i>ACS Agric. Sci. Technol.</i> 2025;5. doi:10.1021/acsagscitech.4c00730.",
            "[31] Chen Q. et al. (Guo lab). De novo biosynthesis of ergothioneine in <i>Yarrowia lipolytica</i> via TrEgt1–TrEgt2 fusion (8.03 g/L). <i>Biotechnol. Bioeng.</i> 2026; Sato S. et al. <i>Appl. Microbiol. Biotechnol.</i> 2025 (PMC11993508).",
            "[32] Команда iGEM HBUT-China (2022). CRISPR/Cas9 integration of <i>NcEgt1</i> and <i>CpEgt2</i> in <i>S. cerevisiae</i> CEN.PK113-5D (протокол экспериментов, 2022.igem.wiki/hbut-china/experiments).",
            "[33] Carlsson J., Kierstan M. P. J., Brocklehurst K. A convenient spectrophotometric assay for the determination of L-ergothioneine (2,2'-dipyridyl disulphide at pH 1.0, 343 nm). <i>Biochemical Journal</i>. 1974;139(1):237–242. doi:10.1042/bj1390237; Nguyen T. H. et al. <i>Food Chemistry</i>. 2012;133:585–591.",
            "[34] Maurer A. The mechanism of enzyme-catalyzed ergothioneine degradation (эрготионаза TdETL, ε₃₁₁ = 22 500 M⁻¹cm⁻¹). PhD Diss., Univ. of Basel, 2020. doi:10.5451/unibas-ep80474; Beliaeva M. A., Leisinger F., Seebeck F. P. <i>ACS Chem. Biol.</i> 2021;16(2):397–403. doi:10.1021/acschembio.0c00968.",
            "[35] Kayrouz C. M., Huang J., Hauser N., Seyedsayamdost M. R. Biosynthesis of selenium-containing small molecules in diverse microorganisms (<i>senABC</i>). <i>Nature</i>. 2022;610:199–204. doi:10.1038/s41586-022-05174-2.",
            "[36] Приказ Минсельхоза России от 06.06.2025 № 392 (рег. Минюста РФ № 83535 от 12.09.2025) «Об утверждении перечня критически важных ферментных препаратов, пищевых и кормовых добавок, технологических вспомогательных средств»; Распоряжение Правительства РФ от 05.10.2026 № 2777-р «Стратегия развития биоэкономики в РФ до 2036 г.».",
            "[37] Аналитические сводки RNC Pharma («Аудит розничных продаж ЛП в России, total sell-out» от 14.01.2026 и итоги I кв. 2026 г.: группа гепатопротекторов A05B и доля препаратов МНН адеметионин); данные ТР ТС 009/2011 и ТР ТС 021/2011 ЕАЭС.",
            "[38] Revuelta J. L., Serrano-Amatriain C., Ledesma-Amaro R., Jiménez A. Formation of folates by microorganisms: towards the biotechnological production of this vitamin. <i>Appl. Microbiol. Biotechnol.</i> 2018;102:8613–8620. doi:10.1007/s00253-018-9266-0.",
            "[39] Патенты РФ на иные способы получения эрготионеина: RU 2788417 C1 (Bloomage Biotechnology Corp. Ltd., 2023: штамм базидиомицета <i>Hericium erinaceus</i> HT-3 CCTCC M 2018567) и RU 2550202 C2 (Tetrahedron SAS, Эрдельмейер И.: химический синтез производных бетаина с бромом); международная заявка DTU на рекомбинантные дрожжи WO 2020/221795 A1 (без национальной фазы в РФ).",
            "[40] FDA GRAS Notice GRN 001270 (Jiangsu GeneIII Biotech Co., Ltd., Jan 2026: L-ergothioneine ≥99.9% produced by engineered <i>E. coli</i> K-12 MG1655, 30 m³ scale); аналитические обзоры тоннажа и мощностей КНР: Emergen Research (2026, ~68 t global volume), Nexchem Intelligence (2026, ~$75M, CAGR ~18%), Onerchem Market Price Assessment (Oct 2026, $400–580/kg); годовой отчёт Yili Chuanning Biotechnology Co., Ltd. за 2025 г.",
            "[41] Ma Y. et al. (Beijing University of Chemical Technology, BUCT). Modular metabolic engineering of <i>Saccharomyces cerevisiae</i> for high-level de novo L-ergothioneine production without amino acid supplementation (3.05 g/L). 2026 (SSRN preprint, Apr 2026); Khan M. A. et al. <i>Acta Microbiologica Sinica</i>. 2026;66(8) (1.31 g/L).",
            "[42] Yan F. et al. Amino-acid-free high-level biosynthesis of ergothioneine in <i>Escherichia coli</i> using betaine and inorganic sulfate (7.2 g/L). <i>Trends in Biotechnology</i>. 2025 (Oct 2025); Zhang H. et al. Metabolic engineering of <i>E. coli</i> for L-ergothioneine production (12.16 g/L in 5-L bioreactor). <i>Journal of Biotechnology</i>. 2026 (Jun 2026); заявка Цзяннаньского ун-та CN 120099041 A (2025, 10.3 g/L в <i>Y. lipolytica</i>).",
            "[43] Государственный реестр лекарственных средств (ГРЛС) Минздрава России, выгрузка ноябрь 2026 г. по МНН «Адеметионин»: регистрационные удостоверения «Гептрал» (П N011968/02, ЛП-№(004983)-(РГ-RU) — АФС Gnosis Bioresearch S.A.), «Адеметионин Велфарм» / «Гептиэль» (ЛП-007320, ЛП-№(014846)-(РГ-RU) — АФС Zhejiang Hisun, Wuxi Fortune, Shandong Jincheng, Shandong Ruiying, ООО «Велфарм», ООО «БЕРАХИМ»), «Самеликс» (ЛП-№(006234)-(РГ-RU)) и реестровая запись фармсубстанции ФС-002473 от 01.07.2026 (ООО «Лиоф Плюс» / ООО «БЕРАХИМ», г. Обнинск).",
        ]
        for r in refs:
            story.append(P(r, st_ref))

        return story

    def render_pass(out_file: Path, toc_in: dict) -> tuple[int, dict]:
        recorded = {}
        doc = ReportDocTemplate(
            str(out_file),
            recorded,
            pagesize=A4,
            leftMargin=left_m,
            rightMargin=right_m,
            topMargin=top_m,
            bottomMargin=bot_m,
            title="Обоснование выбора L-эрготионеина и метаболического узла SAM для цикла DBTL",
            author="Проект DBTL",
        )
        frame = Frame(left_m, bot_m, usable_w, page_h - top_m - bot_m, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=draw_header_footer)])
        doc.build(build_story(toc_in))
        return doc.page, recorded

    # Pass 1 collects exact page numbers for every bookmark; Pass 2 renders clickable TOC with exact page numbers
    p1_count, toc_exact = render_pass(OUT_PDF, {})
    p2_count, toc_verify = render_pass(OUT_PDF, toc_exact)
    assert p1_count == p2_count and toc_exact == toc_verify, f"Two-pass TOC mismatch: {p1_count} vs {p2_count}"

    # Render PNG pages into .context/ergothioneine-report for visual inspection
    pdf_doc = fitz.open(str(OUT_PDF))
    bad_glyph_pages = [p.number + 1 for p in pdf_doc if chr(65533) in p.get_text()]
    assert not bad_glyph_pages, f"Missing glyphs on pages: {bad_glyph_pages}"
    links_on_p1 = len(pdf_doc[0].get_links())
    for old_png in CTX_DIR.glob("page_*.png"):
        old_png.unlink()
    for idx, page in enumerate(pdf_doc, start=1):
        pix = page.get_pixmap(dpi=160, alpha=False)
        pix.save(str(CTX_DIR / f"page_{idx:02d}.png"))
    pdf_doc.close()

    print(f"SUCCESS: {OUT_PDF} ({p2_count} pages, {links_on_p1} clickable TOC links on page 1, {OUT_PDF.stat().st_size} bytes)")


if __name__ == "__main__":
    build_pdf()
