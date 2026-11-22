# -*- coding: utf-8 -*-
"""
Build standalone publication-grade A4 analytical report PDF on selecting the flavin cascade
Riboflavin (Vitamin B2) -> Riboflavin-5'-monophosphate (5'-FMN) -> FAD in yeast
for the Design-Build-Test-Learn (DBTL) metabolic engineering platform.

Generates:
  - assets/fig_fmn1_screening.png
  - assets/fig_fmn2_pathway.png (composited with RDKit 2D chemical structures)
  - assets/fig_fmn3_chem_vs_bio.png
  - assets/fig_fmn4_dbtl_loop.png
  - fmn-riboflavin-dbtl-report.pdf
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
CTX_DIR = ROOT_DIR / ".context" / "fmn-report"
ASSETS_DIR.mkdir(parents=True, exist_ok=True)
CTX_DIR.mkdir(parents=True, exist_ok=True)

OUT_PDF = ROOT_DIR / "fmn-riboflavin-dbtl-report.pdf"


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
WARN_SOFT = colors.HexColor("#FDF6EC")
WARN_BORDER = colors.HexColor("#B86B1F")

SUB_MAP = {
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
    "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9", "ₐ": "a",
}
SUP_MAP = {"⁺": "+", "⁻": "-", "²": "2", "³": "3", "⁴": "4", "¹": "1"}


def fmt_chem(text: str) -> str:
    out = []
    for ch in text:
        if ch in SUB_MAP:
            out.append(f'<sub size="6.3">{SUB_MAP[ch]}</sub>')
        elif ch in SUP_MAP:
            out.append(f'<super size="6.3">{SUP_MAP[ch]}</super>')
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


def make_fig1_screening(out_path: Path) -> None:
    """Figure 1: Screening of 14 candidate platforms on price vs Russian demand & multi-criteria score."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.35), dpi=300)
    fig.patch.set_facecolor("white")

    # Left panel: Price ($/kg, log) vs Physical Import / Demand Scale in Russia (t/year, log)
    pts = [
        ("Таурин", 2.3, 1800, "#8A969E", (10, -3)),
        ("Витамин B2 (Шаг 1)", 24.0, 400, "#1B5E4A", (12, 11)),
        ("Пуллулан (E1204)", 26.0, 45, "#8A969E", (-76, 15)),
        ("2'-FL (ОГМ)", 38.0, 30, "#6E7E88", (-62, -13)),
        ("ФМН-Na USP (~75%)", 95.0, 18.0, "#2F5573", (-96, -24)),
        ("Ресвератрол 98%", 105.0, 3.6, "#8A969E", (-88, -14)),
        ("Ликопин (10%)", 220.0, 14.0, "#B84A39", (-22, 22)),
        ("(S,S)-SAMe (АФС)", 680.0, 22.0, "#2F5573", (8, 15)),
        ("Биотин (Вит. B7)", 340.0, 5.5, "#8A969E", (14, 9)),
        ("5'-ФМН био (>96%, Шаг 2)", 440.0, 4.2, "#1B5E4A", (-118, -28)),
        ("Кофермент ФАД (≥95%, Шаг 3)", 1700.0, 1.6, "#1B5E4A", (-68, -26)),
        ("L-Эрготионеин 99%", 490.0, 0.08, "#6E7E88", (-96, 7)),
    ]

    # Draw cascade arrow connecting B2 -> 5'-FMN -> FAD
    ax1.annotate(
        "", xy=(440.0, 4.2), xytext=(24.0, 400.0),
        arrowprops=dict(arrowstyle="-|>", color="#1B5E4A", lw=1.4, ls="--", shrinkA=6, shrinkB=6),
        zorder=3,
    )
    ax1.annotate(
        "", xy=(1700.0, 1.6), xytext=(440.0, 4.2),
        arrowprops=dict(arrowstyle="-|>", color="#1B5E4A", lw=1.4, ls="--", shrinkA=6, shrinkB=6),
        zorder=3,
    )

    for label, price, tons, col, offset in pts:
        is_main = "Шаг" in label
        sz = 84 if is_main else 42
        ax1.scatter(price, tons, s=sz, color=col, alpha=0.93, zorder=4, edgecolors="white", linewidth=0.9)
        ax1.annotate(
            label,
            (price, tons),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.3 if is_main else 6.8,
            fontweight="bold" if is_main else "normal",
            color="#1B5E4A" if is_main else "#2E373D",
            arrowprops=dict(arrowstyle="-", color="#B8C2C8", lw=0.55, shrinkA=0, shrinkB=3),
        )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1.1, 8500)
    ax1.set_ylim(0.025, 6000)
    ax1.axhspan(1.0, 4000, color="#EAF2EF", alpha=0.35, zorder=1)
    ax1.grid(True, which="major", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax1.set_xlabel("Оптовая цена субстанции B2B ($/кг, лог. шкала)", fontsize=8.3, color="#2E373D")
    ax1.set_ylabel("Оценка физического спроса в РФ (т/год, лог.)", fontsize=8.3, color="#2E373D")
    ax1.set_title("А. Каскад «тоннажный Шаг 1 → маржинальный Шаг 2/3»", fontsize=9.1, fontweight="bold", color="#181D20", pad=8)
    ax1.tick_params(labelsize=7.8, colors="#4A555E")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Right panel: Multi-criteria score out of 25 (5 criteria x 5 pts)
    cands = [
        ("Таурин / Биотин B7", 9.5, "#9AA6AE"),
        ("Ликопин → Астаксантин", 13.0, "#B84A39"),
        ("п-Кумаровая → Ресвератрол", 13.5, "#9AA6AE"),
        ("Пуллулан (E1204)", 14.0, "#9AA6AE"),
        ("Прямой (S,S)-SAMe", 16.5, "#2F5573"),
        ("2'-Фукозиллактоза (ОГМ)", 17.0, "#2F5573"),
        ("L-Эрготионеин → SAMe", 19.0, "#2F5573"),
        ("B2 → 5'-ФМН → ФАД", 23.5, "#1B5E4A"),
    ]
    y_pos = np.arange(len(cands))
    vals = [c[1] for c in cands]
    cols = [c[2] for c in cands]
    bars = ax2.barh(y_pos, vals, height=0.58, color=cols, edgecolor="none", zorder=3)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([c[0] for c in cands], fontsize=7.7, color="#181D20")
    for idx, tick in enumerate(ax2.get_yticklabels()):
        if idx == len(cands) - 1:
            tick.set_fontweight("bold")
            tick.set_color("#1B5E4A")
    for bar, v in zip(bars, vals):
        ax2.text(v + 0.35, bar.get_y() + bar.get_height() / 2, f"{v:.1f} / 25", va="center", fontsize=7.4, fontweight="bold", color="#2E373D")

    ax2.set_xlim(0, 27.5)
    ax2.axvline(20.0, color="#1B5E4A", linestyle=":", linewidth=1.0, zorder=2)
    ax2.grid(True, axis="x", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax2.set_xlabel("Интегральный балл по 5 критериям проекта DBTL (макс. 25)", fontsize=8.3, color="#2E373D")
    ax2.set_title("Б. Сводный рейтинг продуктовых платформ", fontsize=9.1, fontweight="bold", color="#181D20", pad=8)
    ax2.tick_params(axis="x", labelsize=7.8, colors="#4A555E")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout(pad=1.2)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_fig2_pathway(out_path: Path) -> None:
    """Figure 2: Endogenous yeast metabolic network (PPP + IMP bifurcation + RIB1-5 + FMN1/FAD1 + RDKit structures)."""
    fig = plt.figure(figsize=(11.0, 5.75), dpi=300)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 52.5)
    ax.axis("off")

    outer = FancyBboxPatch(
        (0.6, 0.6), 98.8, 51.3,
        boxstyle="round,pad=0.3,rounding_size=1.1",
        facecolor="#FAFCFD", edgecolor="#CFD8DE", linewidth=1.0, zorder=0
    )
    ax.add_patch(outer)

    ax.text(
        2.0, 49.9,
        "ЭНДОГЕННЫЙ МЕТАБОЛИЧЕСКИЙ КАСКАД ДРОЖЖЕЙ (0 ЧУЖЕРОДНЫХ ГЕНОВ) И КОНФЛИКТ РАЗВИЛКИ ИМФ (ГТФ vs АТФ)",
        fontsize=8.9, fontweight="bold", color="#181D20", va="center"
    )

    # Top-left box: Pentose phosphate pathway
    ppp = FancyBboxPatch(
        (2.0, 35.4), 38.5, 12.6,
        boxstyle="round,pad=0.2,rounding_size=0.8",
        facecolor="#EEF4F8", edgecolor="#8FAAB9", linewidth=1.0, zorder=1
    )
    ax.add_patch(ppp)
    ax.text(3.2, 46.1, "Вход А: Окислительный пентозофосфатный путь (ПФП)", fontsize=7.9, fontweight="bold", color="#2F5573", va="center")
    ax.text(3.2, 43.4, "D-Глюкоза → Глюкозо-6-фосфат → [ZWF1 → SOL3 → GND1]", fontsize=7.1, color="#181D20", va="center")
    ax.text(3.2, 40.8, "→ D-Рибулозо-5-фосфат (Ru5P)  ↔  Рибозо-5-фосфат [RKI1]", fontsize=7.1, color="#181D20", va="center")
    ax.text(3.2, 37.6, "• Поток 1: 2 моль Ru5P → вход в фермент RIB3 (ДБФ-синтаза)\n• Поток 2: Рибозо-5-P + АТФ → [PRS1–5] → PRPP (в синтез пуринов)", fontsize=6.65, color="#2E373D", va="center", linespacing=1.22)

    # Top-right box: Purine trunk and IMP bifurcation (key ML target!)
    imp_box = FancyBboxPatch(
        (43.0, 35.4), 55.0, 12.6,
        boxstyle="round,pad=0.2,rounding_size=0.8",
        facecolor="#FDF6EC", edgecolor="#C78C44", linewidth=1.1, zorder=1
    )
    ax.add_patch(imp_box)
    ax.text(44.3, 46.1, "Вход Б: Ствол пуринов [ADE4] и нелинейная развилка ИМФ (IMP) — ядро задачи ML", fontsize=7.9, fontweight="bold", color="#8A4B08", va="center")
    ax.text(44.3, 43.3, "PRPP + Глутамин → [ADE4_fbr, ADE5,7..ADE1] → Инозинмонофосфат (ИМФ / IMP):", fontsize=7.05, color="#181D20", va="center")
    ax.text(44.3, 40.4, "(1) Гуаниловая ветвь [↑IMD2/3 → GUA1 → GUK1]: ИМФ → ГМФ → 1 ГТФ (субстрат RIB1 на сборку B₂)", fontsize=6.85, color="#1B5E4A", fontweight="bold", va="center")
    ax.text(44.3, 37.5, "(2) Адениловая ветвь [ADE12 → ADE13 → ADK1]: ИМФ → АМФ → 1–2 АТФ (донор фосфата для FMN1 и FAD1)\n    Парадокс: блок ↓ADE12 ускоряет Шаг 1 (B₂), но лишает АТФ реакции Шага 2 (FMN1 / FAD1)!", fontsize=6.65, color="#7A3E06", va="center", linespacing=1.2)

    # Arrow from PPP to IMP box
    a_prpp = FancyArrowPatch((40.6, 38.2), (42.9, 38.2), arrowstyle="-|>,head_length=4,head_width=2.4", color="#5E6A73", linewidth=1.1, zorder=4)
    ax.add_patch(a_prpp)

    def place_card(x, y, w, h, img_key, title, formula, sub, border, bg, title_col):
        card = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.2,rounding_size=0.9",
            facecolor=bg, edgecolor=border, linewidth=1.35, zorder=2
        )
        ax.add_patch(card)
        ax.text(x + w / 2, y + h - 1.75, title, fontsize=8.05, fontweight="bold", color=title_col, ha="center", va="center", zorder=4)
        ax.text(x + w / 2, y + h - 3.6, formula, fontsize=6.8, color="#5E6A73", ha="center", va="center", zorder=4)

        arr = crop_structure_png(STRUCT_DIR / f"{img_key}.png", pad=15)
        ih, iw = arr.shape[:2]
        bot_pad = 4.6
        box_w = w - 2.0
        box_h = h - 4.8 - bot_pad
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
        ax.text(x + w / 2, y + 2.35, sub, fontsize=6.5, color="#2E373D", ha="center", va="center", linespacing=1.2, zorder=4)

    yc, hc = 2.0, 26.6
    place_card(
        2.0, yc, 25.0, hc, "riboflavin",
        "Шаг 1: Рибофлавин (Вит. B2)",
        "C₁₇H₂₀N₄O₆ · р-римость 0,08 г/л",
        "Сборка: RIB1→RIB7→RIB2 + RIB3 → RIB4 → RIB5\nЭкспорт RFE1 · Тест: 445 нм за 10 с без реактивов",
        "#1B5E4A", "#EAF2EF", "#1B5E4A"
    )
    place_card(
        36.8, yc, 25.4, hc, "fmn",
        "Шаг 2 (цель): 5'-ФМН",
        "C₁₇H₂₁N₄O₉P · р-римость 67 г/л",
        "Строго 5'-OH фосфорилирование (рибофлавинкиназа)\nКонтроль холостого гидролиза: ↓PHO8 / ↓фосфатазы",
        "#1B5E4A", "#F3FAF7", "#1B5E4A"
    )
    place_card(
        71.8, yc, 26.2, hc, "fad",
        "Шаг 3 (расширение): ФАД",
        "C₂₇H₃₃N₉O₁₅P₂ · р-римость 50 г/л",
        "Монофункциональная ФАД-синтетаза FAD1 (+АТФ)\nСекреция / высвобождение из флавопротеинов",
        "#2F5573", "#F4F8FB", "#2F5573"
    )

    def rxn_arrow(x0, x1, y, top_t, bot_t):
        arr = FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>,head_length=5,head_width=3", color="#1B5E4A", linewidth=1.6, zorder=5)
        ax.add_patch(arr)
        xm = (x0 + x1) / 2
        ax.text(xm, y + 1.4, top_t, fontsize=7.25, fontweight="bold", color="#1B5E4A", ha="center", va="bottom")
        ax.text(xm, y - 1.4, bot_t, fontsize=6.2, color="#2E373D", ha="center", va="top", linespacing=1.18)

    rxn_arrow(27.4, 36.4, 15.8, "↑FMN1 + АТФ", "−АДФ (строго 5'-O-P)\nобратный срыв:\nфосфатазы PHO8")
    rxn_arrow(62.6, 71.4, 15.8, "↑FAD1 + АТФ", "−PPᵢ (аденилирование)\nкофактор оксидаз\nи FAD-GDH")

    # Orthogonal feeder routing in the open corridor y = 28.8 .. 35.2 (zero card crossings!)
    ax.add_patch(FancyArrowPatch((11.5, 35.3), (11.5, 28.9), arrowstyle="-|>,head_length=4.2,head_width=2.5", color="#2F5573", linewidth=1.2, zorder=4))
    ax.text(12.4, 32.1, "2 Ru5P (RIB3)", fontsize=6.5, fontweight="bold", color="#2F5573", ha="left", va="center")

    ax.plot([49.5, 49.5, 23.2], [35.3, 33.1, 33.1], color="#1B5E4A", linewidth=1.15, linestyle="--", zorder=4)
    ax.add_patch(FancyArrowPatch((23.2, 33.1), (23.2, 28.9), arrowstyle="-|>,head_length=4.2,head_width=2.5", color="#1B5E4A", linewidth=1.15, linestyle="--", zorder=4))
    ax.text(36.5, 33.8, "1 ГТФ из ветви IMD2 (в RIB1)", fontsize=6.45, fontweight="bold", color="#1B5E4A", ha="center", va="center")

    ax.plot([85.0, 85.0, 31.9], [35.3, 30.8, 30.8], color="#B86B1F", linewidth=1.15, linestyle="--", zorder=4)
    ax.add_patch(FancyArrowPatch((31.9, 30.8), (31.9, 19.6), arrowstyle="-|>,head_length=4.2,head_width=2.5", color="#B86B1F", linewidth=1.15, linestyle="--", zorder=4))
    ax.add_patch(FancyArrowPatch((67.0, 30.8), (67.0, 19.6), arrowstyle="-|>,head_length=4.2,head_width=2.5", color="#B86B1F", linewidth=1.15, linestyle="--", zorder=4))
    ax.text(76.0, 31.6, "Пул АТФ из ветви ADE12", fontsize=6.45, fontweight="bold", color="#8A4B08", ha="center", va="center")

    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_fig3_chem_vs_bio(out_path: Path) -> None:
    """Figure 3: Chemical POCl3 semi-synthesis vs Regioselective microbial biosynthesis (purity & pricing)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.15), dpi=300)
    fig.patch.set_facecolor("white")

    # Panel A: Flavin composition (%) in commercial POCl3 USP/Ph.Eur. vs Yeast FMN1 biosynthesis
    routes = ["Хим. полусинтез\n(B2 + 20× POCl3,\nстандарт USP/ЕФ)", "Биосинтез в дрожжах\n(рибофлавинкиназа FMN1,\nпосле водной очистки)"]
    c_5fmn = np.array([75.0, 96.5])
    c_iso = np.array([15.5, 0.5])
    c_diphos = np.array([5.0, 0.8])
    c_freeb2 = np.array([4.5, 2.2])

    x = np.arange(len(routes))
    w = 0.46
    ax1.bar(x, c_5fmn, width=w, color="#1B5E4A", label="Целевой 5'-ФМН (активный кофермент)")
    ax1.bar(x, c_iso, width=w, bottom=c_5fmn, color="#D97B29", label="3'- и 4'-ФМН, цикл. фосфаты (ингибиторы)")
    ax1.bar(x, c_diphos, width=w, bottom=c_5fmn + c_iso, color="#B84A39", label="Рибофлавиндифосфаты (≤6% по ЕФ/USP)")
    ax1.bar(x, c_freeb2, width=w, bottom=c_5fmn + c_iso + c_diphos, color="#9AA6AE", label="Непрореагировавший свободный B2 (≤6%)")

    ax1.text(0, 37.5, "~73–79%\n5'-изомера", ha="center", va="center", color="white", fontsize=8.0, fontweight="bold")
    ax1.text(0, 83.0, "21–27%\nпримесей", ha="center", va="center", color="white", fontsize=7.3, fontweight="bold")
    ax1.text(1, 48.0, "≥95–98%\nрегиоселективный\n5'-ФМН", ha="center", va="center", color="white", fontsize=8.0, fontweight="bold")

    ax1.set_xticks(x)
    ax1.set_xticklabels(routes, fontsize=7.7, color="#181D20")
    ax1.set_ylim(0, 114)
    ax1.set_ylabel("Массовая доля флавиновых фракций (%)", fontsize=8.2, color="#2E373D")
    ax1.set_title("А. Изомерный состав: химия POCl3 против фермента FMN1", fontsize=9.0, fontweight="bold", color="#181D20", pad=8)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False, fontsize=6.9)
    ax1.grid(True, axis="y", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Market value ladder ($/kg, log scale) across flavin grades (2025-2026 verified trade data)
    grades = [
        "Вит. B2 кормовой 80%\n(опт КНР)",
        "Вит. B2 фарм 98%\n(ферментация)",
        "ФМН-Na ~75% USP\n(Индия/КНР, POCl3)",
        "ФМН-Na ~75% Ph.Eur.\n(DSM Франция/Швейц.)",
        "Чистый 5'-ФМН ≥95%\n(биокатализ / реагент)",
        "ФАД-Na2 ≥95% фарм/IVD\n(только биосинтез!)",
    ]
    lows = [15, 38, 63, 203, 450, 1200]
    highs = [22, 46, 120, 384, 1500, 3000]
    mids = [(l + h) / 2 for l, h in zip(lows, highs)]
    y_g = np.arange(len(grades))
    g_cols = ["#9AA6AE", "#7B8B94", "#D97B29", "#2F5573", "#1B5E4A", "#1B5E4A"]

    for i, (y, lo, hi, col) in enumerate(zip(y_g, lows, highs, g_cols)):
        ax2.plot([lo, hi], [y, y], color=col, linewidth=5.5, solid_capstyle="round", zorder=3)
        ax2.scatter([lo, hi], [y, y], color=col, s=26, zorder=4)
        ax2.text(hi * 1.16, y, f"${lo}–{hi}/кг", va="center", fontsize=7.3, fontweight="bold" if i >= 4 else "normal", color="#181D20")

    ax2.set_xscale("log")
    ax2.set_xlim(10, 13000)
    ax2.set_yticks(y_g)
    ax2.set_yticklabels(grades, fontsize=7.3, color="#181D20")
    ax2.set_xlabel("Торгуемый ценовой коридор 2025–2026 гг. ($/кг, лог. шкала)", fontsize=8.2, color="#2E373D")
    ax2.set_title("Б. Ценовая лестница сортов B2, ФМН и ФАД", fontsize=9.0, fontweight="bold", color="#181D20", pad=8)
    ax2.grid(True, axis="x", linestyle="--", linewidth=0.5, color="#E1E6EA", zorder=0)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout(pad=1.15)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)


def make_fig4_dbtl_loop(out_path: Path) -> None:
    """Figure 4: Two-stage Transfer-DBTL architecture with Experimental Controls (Exp 0-6)."""
    fig = plt.figure(figsize=(10.8, 4.45), dpi=300)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 42)
    ax.axis("off")

    outer = FancyBboxPatch(
        (0.6, 0.6), 98.8, 40.8,
        boxstyle="round,pad=0.3,rounding_size=1.0",
        facecolor="#FAFCFD", edgecolor="#CFD8DE", linewidth=1.0, zorder=0
    )
    ax.add_patch(outer)

    ax.text(
        2.0, 39.2,
        "АРХИТЕКТУРА ДВУХЭТАПНОГО КОНВЕЙЕРА TRANSFER-DBTL И СИСТЕМА КОНТРОЛЬНЫХ ОПЫТОВ (ОПЫТЫ 0–6)",
        fontsize=8.8, fontweight="bold", color="#181D20", va="center"
    )

    # Gate 0 box on the left
    gate = FancyBboxPatch(
        (2.0, 4.2), 20.5, 32.2,
        boxstyle="round,pad=0.2,rounding_size=0.8",
        facecolor="#EEF4F8", edgecolor="#7B9BB0", linewidth=1.1, zorder=1
    )
    ax.add_patch(gate)
    ax.text(12.25, 34.4, "Входной контроль\n(Опыты 0–2)", fontsize=7.9, fontweight="bold", color="#2F5573", ha="center", va="center")
    ax.text(
        3.1, 19.0,
        "• Опыт 0 (Аналитика):\n  калибровка 3 стандартов\n  (B₂, 5'-ФМН, ФАД) на\n  445 нм и ВЭЖХ RP-C18;\n  проверка потерь и света.\n\n"
        "• Опыт 1 (Шасси):\n  сравнение K. phaffii /\n  S. cerevisiae (RIB1+3)\n  + контроль ВКПМ B-2694.\n\n"
        "• Опыт 2 (Диагностика):\n  подача внешнего B₂ для\n  отдельной оценки узла\n  B₂ → 5'-ФМН (АТФ, FMN1,\n  фосфатазы PHO8).",
        fontsize=6.45, color="#181D20", va="center", linespacing=1.22
    )

    # Four DBTL boxes in center-right top & bottom
    def dbtl_box(x, y, w, h, title, body, bg, border, tcol):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=0.8", facecolor=bg, edgecolor=border, linewidth=1.2, zorder=2)
        ax.add_patch(b)
        ax.text(x + w / 2, y + h - 2.1, title, fontsize=7.9, fontweight="bold", color=tcol, ha="center", va="center", zorder=3)
        ax.text(x + 1.2, y + (h - 3.2) / 2, body, fontsize=6.45, color="#181D20", va="center", linespacing=1.22, zorder=3)

    dbtl_box(
        26.0, 21.6, 33.5, 14.8,
        "1. DESIGN (Гибрид GEM + Active Learning)",
        "• Ограничения потоков ecYeast8 / iKp1260 (FSEOF)\n• Матрица 6 узлов: [ZWF1/RPE1] × [ADE4] ×\n  [развилка ИМФ: IMD2 vs ADE12] × [RIB1/3] ×\n  [FMN1 + ↓PHO8 + RFE1] × среда (Fe²⁺, фосфат)",
        "#EAF2EF", "#1B5E4A", "#1B5E4A"
    )
    dbtl_box(
        64.0, 21.6, 33.8, 14.8,
        "2. BUILD (Верифицированные штаммы)",
        "• ПЦР собственных генов дрожжей (0 чужеродных)\n• Модульная сборка промоторных кассет (24–36\n  уникальных генотипов на раунд × 3 био-повтора)\n• Паспорт штамма: локус, копийность, пересевы",
        "#F5F7F8", "#8A969E", "#2E373D"
    )
    dbtl_box(
        64.0, 4.2, 33.8, 14.8,
        "3. TEST (Двухуровневая оптика + ВЭЖХ)",
        "• Шаг 1 (B₂): супернатант 445 нм за 10 с без реагентов\n• Шаг 2 (B₂ / 5'-ФМН / ФАД): изократический ВЭЖХ\n  RP-C18 (УФ 445 нм / ФЛД 450→525 нм, 5 мин)\n• Раздельный учёт среды, клеток, биомассы и сахара",
        "#F4F8FB", "#2F5573", "#2F5573"
    )
    dbtl_box(
        26.0, 4.2, 33.5, 14.8,
        "4. LEARN (Transfer-ML и честный контроль)",
        "• GP / Random Forest с разбиением по генотипам\n• Перенос предобученной модели с B₂ на 5'-ФМН\n• Слепое проспективное сравнение с Random и\n  экспертным выбором при равном бюджете опытов",
        "#FDF6EC", "#C78C44", "#8A4B08"
    )

    # Arrows connecting Gate -> Design -> Build -> Test -> Learn -> Design
    arrs = [
        ((22.6, 29.0), (25.8, 29.0), "#2F5573"),
        ((59.7, 29.0), (63.8, 29.0), "#1B5E4A"),
        ((80.9, 21.4), (80.9, 19.2), "#2F5573"),
        ((63.8, 11.6), (59.7, 11.6), "#2F5573"),
        ((42.75, 19.2), (42.75, 21.4), "#C78C44"),
    ]
    for p0, p1, c in arrs:
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>,head_length=4.5,head_width=2.7", color=c, linewidth=1.45, zorder=5))

    ax.text(61.75, 30.2, "Цикл", fontsize=6.3, color="#1B5E4A", ha="center")
    ax.text(46.8, 20.3, "Новые кандидаты", fontsize=6.3, fontweight="bold", color="#8A4B08", ha="left", va="center")

    fig.savefig(out_path, dpi=300)
    plt.close(fig)


class ReportDocTemplate(BaseDocTemplate):
    def __init__(self, filename: str, toc_map: dict, **kwargs):
        super().__init__(filename, **kwargs)
        self.toc_map = toc_map

    def afterFlowable(self, flowable):
        if isinstance(flowable, BookmarkHeading):
            self.toc_map[flowable.toc_key] = self.page


def build_pdf() -> None:
    fig1_path = ASSETS_DIR / "fig_fmn1_screening.png"
    fig2_path = ASSETS_DIR / "fig_fmn2_pathway.png"
    fig3_path = ASSETS_DIR / "fig_fmn3_chem_vs_bio.png"
    fig4_path = ASSETS_DIR / "fig_fmn4_dbtl_loop.png"

    make_fig1_screening(fig1_path)
    make_fig2_pathway(fig2_path)
    make_fig3_chem_vs_bio(fig3_path)
    make_fig4_dbtl_loop(fig4_path)

    page_w, page_h = A4
    left_m = right_m = 18 * mm
    top_m = 19 * mm
    bot_m = 18 * mm
    usable_w = page_w - left_m - right_m

    st_kicker = ParagraphStyle("kicker", fontName="SansBold", fontSize=8.3, leading=11, textColor=ACCENT_TEAL, spaceAfter=3)
    st_title = ParagraphStyle("title", fontName="SansBold", fontSize=15.8, leading=20.5, textColor=INK, spaceAfter=5)
    st_subtitle = ParagraphStyle("subtitle", fontName="Sans", fontSize=9.5, leading=13.5, textColor=SUBINK, spaceAfter=8)
    st_h1 = ParagraphStyle("h1", fontName="SansBold", fontSize=12.6, leading=16.6, textColor=INK, spaceBefore=2, spaceAfter=5.5, keepWithNext=True)
    st_h2 = ParagraphStyle("h2", fontName="SansBold", fontSize=10.3, leading=13.8, textColor=ACCENT_TEAL, spaceBefore=7.5, spaceAfter=3.5, keepWithNext=True)
    st_body = ParagraphStyle("body", fontName="Sans", fontSize=9.05, leading=13.1, textColor=INK, alignment=TA_JUSTIFY, spaceAfter=4.5)
    st_bullet = ParagraphStyle("bullet", fontName="Sans", fontSize=8.9, leading=12.8, textColor=INK, leftIndent=11, firstLineIndent=-7, alignment=TA_JUSTIFY, spaceAfter=3.4)
    st_subbullet = ParagraphStyle("subbullet", fontName="Sans", fontSize=8.65, leading=12.3, textColor=SUBINK, leftIndent=21, firstLineIndent=-7, alignment=TA_LEFT, spaceAfter=2.6)
    st_callout = ParagraphStyle("callout", fontName="Sans", fontSize=8.85, leading=12.7, textColor=INK, alignment=TA_LEFT, spaceAfter=2)
    st_caption = ParagraphStyle("caption", fontName="SansItalic", fontSize=8.0, leading=11.0, textColor=MUTED, spaceBefore=2.5, spaceAfter=6)
    st_th = ParagraphStyle("th", fontName="SansBold", fontSize=7.85, leading=10.4, textColor=INK)
    st_td = ParagraphStyle("td", fontName="Sans", fontSize=7.75, leading=10.4, textColor=INK)
    st_td_bold = ParagraphStyle("td_bold", fontName="SansBold", fontSize=7.75, leading=10.4, textColor=INK)
    st_toc_h = ParagraphStyle("toc_h", fontName="SansBold", fontSize=10.6, leading=14, textColor=INK, spaceBefore=4, spaceAfter=4)
    st_toc_1 = ParagraphStyle("toc_1", fontName="SansBold", fontSize=8.7, leading=12.0, textColor=INK)
    st_toc_2 = ParagraphStyle("toc_2", fontName="Sans", fontSize=8.25, leading=11.2, textColor=SUBINK, leftIndent=8)
    st_toc_p = ParagraphStyle("toc_p", fontName="SansBold", fontSize=8.4, leading=11.2, textColor=ACCENT_TEAL)
    st_ref = ParagraphStyle("ref", fontName="Sans", fontSize=7.8, leading=10.8, textColor=SUBINK, leftIndent=13, firstLineIndent=-13, spaceAfter=2.4)

    toc_entries = [
        ("sec1", 1, "1. Постановка задачи, критерии отбора и архитектура каскада B2 → 5'-ФМН → ФАД"),
        ("sec1_1", 2, "1.1. Целевая функция цикла DBTL и пять обязательных критериев выбора продукта"),
        ("sec1_2", 2, "1.2. Сводный профиль выбранной платформы: разделение ролей Шага 1, Шага 2 и Шага 3"),
        ("sec2", 1, "2. Системный скрининг 16 альтернатив и обоснование выбора флавиновой цепи"),
        ("sec2_1", 2, "2.1. Сравнительная матрица классов молекул-кандидатов"),
        ("sec2_2", 2, "2.2. Почему каскад B2 → 5'-ФМН снимает ограничения других кандидатов"),
        ("sec3", 1, "3. Химический полусинтез (POCl3) против прямой микробной ферментации 5'-ФМН и ФАД"),
        ("sec3_1", 2, "3.1. Как устроена индустрия сегодня: ферментационный B2 и хлорокись фосфора"),
        ("sec3_2", 2, "3.2. Четыре технологических предела химии POCl3 и ниша биосинтеза"),
        ("sec4", 1, "4. Биохимия пути в дрожжах, штаммы БРЦ ВКПМ и референс «Цикла 0»"),
        ("sec4_1", 2, "4.1. Эндогенный путь дрожжей (0 чужеродных генов) против бактериального оперона"),
        ("sec4_2", 2, "4.2. Пять узких мест перехода B2 → 5'-ФМН → ФАД и парадокс развилки ИМФ"),
        ("sec4_3", 2, "4.3. Выбор штамма-шасси и каталог доступных референсов (включая БРЦ ВКПМ)"),
        ("sec5", 1, "5. Аудит мировой литературы по DBTL/ML и математическая постановка Transfer-DBTL"),
        ("sec5_1", 2, "5.1. Что уже сделано в науке к ноябрю 2026 г. и где открытая зона новизны"),
        ("sec5_2", 2, "5.2. Пространство факторов, гибрид GEM + GP и защита от утечки данных в ML"),
        ("sec6", 1, "6. Двухуровневая аналитика стадии Test, контрольные опыты (0–6) и дорожная карта"),
        ("sec6_1", 2, "6.1. Оптический скрининг 445 нм, хроматография RP-C18 и обязательные контроли"),
        ("sec6_2", 2, "6.2. Последовательность экспериментов от диагностических вопросов (Опыты 0–6)"),
        ("sec6_3", 2, "6.3. План 12-недельного пилота (MVP), двухлетний горизонт и критерии GO / NO-GO"),
        ("sec7", 1, "7. Рынки мира и РФ, реестр ГРЛС, Приказ Минсельхоза № 392 и список источников"),
        ("sec7_1", 2, "7.1. Мировой рынок B2, ФМН-Na и ФАД: объёмы и реальные экспортные цены"),
        ("sec7_2", 2, "7.2. Спрос в РФ, структура ГРЛС («Цитофлавин» / «Полисинтез»), господдержка и ТЭО"),
        ("sec7_3", 2, "7.3. Список проверяемых первоисточников, патентов и реестровых записей"),
    ]

    def H1(key: str, title: str) -> BookmarkHeading:
        return BookmarkHeading(title, st_h1, key, 1, title)

    def H2(key: str, title: str) -> BookmarkHeading:
        return BookmarkHeading(title, st_h2, key, 2, title)

    def callout_box(text: str, bg=ACCENT_SOFT, border=ACCENT_TEAL) -> Table:
        tbl = Table([[P(text, st_callout)]], colWidths=[usable_w])
        tbl.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), bg),
                    ("LINEBEFORE", (0, 0), (0, -1), 3.0, border),
                    ("LEFTPADDING", (0, 0), (-1, -1), 9),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 9),
                    ("TOPPADDING", (0, 0), (-1, -1), 6.5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6.5),
                ]
            )
        )
        return tbl

    def styled_table(data, col_widths) -> Table:
        t = Table(data, colWidths=col_widths, repeatRows=1)
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), SOFT_BG),
                    ("LINEABOVE", (0, 0), (-1, 0), 0.85, INK),
                    ("LINEBELOW", (0, 0), (-1, 0), 0.65, INK),
                    ("LINEBELOW", (0, 1), (-1, -1), 0.4, LINE),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4.2),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4.2),
                    ("TOPPADDING", (0, 0), (-1, -1), 3.8),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 3.8),
                ]
            )
        )
        return t

    def draw_header_footer(canvas, doc):
        canvas.saveState()
        if doc.page > 1:
            canvas.setFont("Sans", 7.5)
            canvas.setFillColor(MUTED)
            canvas.drawString(left_m, page_h - 12 * mm, "Проект DBTL · Метаболическая инженерия дрожжей: каскад Рибофлавин (B2) → 5'-ФМН → ФАД")
            canvas.drawRightString(page_w - right_m, page_h - 12 * mm, "Аналитический отчёт · Ноябрь 2026")
            canvas.setStrokeColor(LINE)
            canvas.setLineWidth(0.5)
            canvas.line(left_m, page_h - 13.5 * mm, page_w - right_m, page_h - 13.5 * mm)

        canvas.setStrokeColor(LINE)
        canvas.setLineWidth(0.5)
        canvas.line(left_m, 12.5 * mm, page_w - right_m, 12.5 * mm)
        canvas.setFont("Sans", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(left_m, 8.8 * mm, "Факты первоисточников, расчёты ТЭО и проверяемые гипотезы разделены явно")
        canvas.drawRightString(page_w - right_m, 8.8 * mm, f"Стр. {doc.page}")
        canvas.restoreState()

    def build_story(toc_pages: dict) -> list:
        story = []

        # ==================== COVER + TOC (PAGE 1) ====================
        story.append(P("ПРОЕКТ DBTL · МЕТАБОЛИЧЕСКАЯ ИНЖЕНЕРИЯ ДРОЖЖЕЙ И МАШИННОЕ ОБУЧЕНИЕ · НОЯБРЬ 2026", st_kicker))
        story.append(
            P(
                "Обоснование выбора каскада Рибофлавин (B2) → Рибофлавин-5'-мононуклеотид (5'-ФМН) → ФАД для отработки цикла Design–Build–Test–Learn",
                st_title,
            )
        )
        story.append(
            P(
                "Полная аналитическая записка: сравнение с химическим полусинтезом (POCl₃), эндогенная биохимия дрожжей и конфликт пуриновой развилки ИМФ, референсные штаммы БРЦ ВКПМ, протокол двухуровневой аналитики, контрольные опыты и экономика рынка РФ (Приказ Минсельхоза № 392, п. 32; реестр ГРЛС).",
                st_subtitle,
            )
        )
        story.append(
            callout_box(
                "<b>Главный вывод отчёта:</b> каскад <b>Рибофлавин (B₂) → 5'-ФМН (→ опционально ФАД)</b> в дрожжах решает главную структурную проблему выбора мишени для DBTL. На <b>Шаге 1</b> все 6 генов синтеза B₂ (<i>RIB1–RIB7</i>) уже присутствуют в базовых дрожжах, позволяя без синтеза чужеродной ДНК за 2 недели получить секретирующий жёлтый штамм с 10-секундным безреагентным тестом при <b>445 нм</b> в продукте с российским импортом <b>300–500 т/год</b> (прямая поз. 32 Приказа Минсельхоза РФ № 392). На <b>Шаге 2</b> добавление 1 собственного гена рибофлавинкиназы (<i>FMN1</i>) переводит B₂ в высокомаржинальный региоселективный <b>5'-ФМН (&gt;95–98% чистого 5'-изомера)</b> из отечественного сахара в одном водном биореакторе — без токсичного 20-кратного избытка <i>POCl₃</i>, без 21–27% изомерных примесей химического полусинтеза и с нетривиальной задачей <b>трансферного машинного обучения (Transfer-DBTL)</b> на развилке пуринов <i>ГТФ / АТФ</i>."
            )
        )
        story.append(Spacer(1, 3 * mm))
        story.append(P("Оглавление (кликните по разделу для быстрого перехода)", st_toc_h))

        toc_rows = []
        for key, lvl, title in toc_entries:
            pg = toc_pages.get(key, 1)
            st_item = st_toc_1 if lvl == 1 else st_toc_2
            link_html = f'<a href="#{key}" color="#181D20">{fmt_chem(title)}</a>'
            pg_html = f'<a href="#{key}" color="#1B5E4A">{pg}</a>'
            toc_rows.append([Paragraph(link_html, st_item), Paragraph(pg_html, st_toc_p)])

        toc_tbl = Table(toc_rows, colWidths=[usable_w - 14 * mm, 14 * mm])
        toc_tbl.setStyle(
            TableStyle(
                [
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("ALIGN", (1, 0), (1, -1), "RIGHT"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 2),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                    ("TOPPADDING", (0, 0), (-1, -1), 1.7),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 1.7),
                    ("LINEBELOW", (0, 0), (-1, -1), 0.25, LINE),
                ]
            )
        )
        story.append(toc_tbl)

        # ==================== SECTION 1 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec1", "1. Постановка задачи, критерии отбора и архитектура каскада B2 → 5'-ФМН → ФАД"))

        story.append(H2("sec1_1", "1.1. Целевая функция цикла DBTL и пять обязательных критериев выбора продукта"))
        story.append(
            P(
                "Цель проекта — создать и экспериментально верифицировать в России воспроизводимый вычислительно-лабораторный конвейер <b>Design–Build–Test–Learn (DBTL)</b> для метаболической инженерии эукариотических микроорганизмов (дрожжей). Вместо эмпирического перебора одиночных мутаций система опирается на комбинацию стехиометрической модели метаболизма и алгоритма активного обучения (<i>Active Learning</i>), который по результатам малой серии культур выбирает следующую серию генетических конструкций и условий среды для максимизации объёмной производительности (г/л·ч), выхода на субстрат (г/г сахара) и селективности целевого продукта [1–5].",
                st_body,
            )
        )
        story.append(
            P(
                "Чтобы не застрять на первом этапе из-за недоступной аналитики, отсутствия референсного протокола или выбора коммерчески невостребованного вещества, целевая продуктовая платформа обязана одновременно выполнять пять строгих критериев:",
                st_body,
            )
        )
        crit_bullets = [
            "<b>Критерий 1 — Гарантированный воспроизводимый старт («Цикл 0»):</b> продукт первого шага уже умеют получать в рекомбинантных микроорганизмах по открытым протоколам с минимальным числом генетических правок, чтобы учебно-исследовательская группа за первые недели гарантированно получила измеримый положительный сигнал выше фона.",
            "<b>Критерий 2 — Быстрое и доступное измерение на стадии Test без масс-спектрометрии:</b> на этапе отладки конвейера концентрация продукта Шага 1 должна определяться простым оптическим методом (прямое поглощение или флуоресценция супернатанта в кювете или 96-луночном планшете), а на Шаге 2 — разделяться обычным жидкостным хроматографом с УФ/флуоресцентным детектором (ВЭЖХ-УФ/ФЛД) без дорогостоящего ВЭЖХ-МС/МС.",
            "<b>Критерий 3 — Открытая научная ниша для ML-DBTL и патентная чистота в РФ:</b> выбранная связка не должна быть многократно закрыта опубликованными работами по машинному обучению (как ликопин или бета-каротин [4, 10, 11]) и не должна нарушать действующие патенты РФ и ЕАПО.",
            "<b>Критерий 4 — Прямая химическая преемственность Шага 1 и Шага 2 (Transfer-DBTL):</b> молекула первого шага должна служить буквальным биохимическим предшественником более сложного и дорогого продукта второго шага, чтобы штаммы, аналитика и обученная суррогатная модель переносились вперёд без слома метаболической архитектуры, но при этом содержали нетривиальный нелинейный компромисс потоков для проверки эффективности ML.",
            "<b>Критерий 5 — Подтверждённый спрос в России и прямое попадание в государственные приоритеты:</b> продукты цепочки должны иметь реальных индустриальных потребителей в РФ (кормовая индустрия, фармацевтические заводы, производство медизделий и диагностики) и напрямую входить в отраслевые перечни импортозамещения (Приказ Минсельхоза России от 06.06.2025 № 392 [25], стратегия «Фарма-2030», нацпроект «Биоэкономика»).",
        ]
        for b in crit_bullets:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec1_2", "1.2. Сводный профиль выбранной платформы: разделение ролей Шага 1, Шага 2 и Шага 3"))
        story.append(
            P(
                "Аудит 16 биотехнологических молекул показал, что всем пяти критериям наиболее полно удовлетворяет эндогенный флавиновый каскад дрожжей <b>Рибофлавин (витамин B₂) → Рибофлавин-5'-мононуклеотид (5'-ФМН) → Флавинадениндинуклеотид (ФАД)</b>, в котором каждый этаж выполняет чётко разграниченную технологическую и экономическую функцию:",
                st_body,
            )
        )
        tiers = [
            "<b>Шаг 1 — Рибофлавин (Витамин B₂, C₁₇H₂₀N₄O₆, 376,4 г/моль) как безотказный оптический полигон и сырьевой фундамент:</b> природный жёлтый пигмент и предшественник всех флавиновых коферментов. Все 6 генов его сборки (<i>RIB1–RIB7</i>) изначально присутствуют в геноме любых базовых дрожжей, продукт активно секретируется наружу в минеральную среду и измеряется напрямую за <b>10 секунд при 445 нм без единого реагента</b> (ε₄₄₅ = 12 500 М⁻¹·см⁻¹). В России 100% потребляемого витамина B₂ (<b>300–500 тонн/год</b> для комбикормов, пищевой добавки E101 и фармсубстанций) импортируется при нулевой отечественной ферментации, а позиция <b>№ 32 Приказа Минсельхоза РФ № 392</b> (<i>«Витамин B2 (рибофлавин), его соединения и производные»</i>) открывает доступ к субсидированию капитальных затрат [25]. <i>Честная экономическая граница:</i> сам кормовой B₂ является крупнотоннажным биржевым товаром ($15–22/кг кормовой 80%, $38–46/кг фарм 98% [27]), поэтому на малом реакторе (0,5–2 м³) его не продают как отдельный биржевой товар — он работает идеальной дешёвой тест-системой для обучения модели DBTL и собственным внутренним субстратом для Шага 2.",
            "<b>Шаг 2 (основная прикладная цель) — Региоселективный Рибофлавин-5'-мононуклеотид (5'-ФМН, C₁₇H₂₁N₄O₉P, 456,3 г/моль):</b> прямое ферментативное фосфорилирование первичного 5'-гидроксила рибофлавина собственным дрожжевым геном <b>рибофлавинкиназы (<i>FMN1</i> + 1 АТФ)</b> повышает растворимость в воде более чем в 500 раз (с ~0,08 г/л у B₂ до <b>67–100 г/л</b> у натриевой соли 5'-ФМН [12]). В отличие от применяемого сегодня двухзаводского химического полусинтеза хлорокисью фосфора (<i>POCl₃</i>), который даёт лишь <b>73–79% целевого 5'-изомера</b> и 21–27% изомерных примесей [13, 14], биосинтез из сахара в одном водном биореакторе даёт <b>&gt;95–98% изомерно чистого 5'-ФМН</b> для инъекционных препаратов перечня ЖНВЛП (<b>«Цитофлавин®»</b> НТФФ «ПОЛИСАН», «Рибофлавин-мононуклеотид» [28, 29]), офтальмологического кросслинкинга роговицы («Декстралинк», «Риболинк», FDA-одобренный в окт. 2025 г. <i>Epioxa™</i> [15]) и биокатализа.",
            "<b>Шаг 3 (опциональное расширение 2-го года) — Кофермент ФАД (C₂₇H₃₃N₉O₁₅P₂, 785,6 г/моль):</b> присоединение аденилата вторым собственным геном дрожжей <b>ФАД-синтетазой (<i>FAD1</i> + ещё 1 АТФ)</b> превращает 5'-ФМН в динуклеотид ФАД ($1 200–3 000+/кг за фарм- и IVD-квалификацию ≥95% [26]), для которого промышленного химического синтеза не существует в принципе и который необходим для тест-полосок современных глюкометров на ферменте <i>FAD-GDH</i> и офтальмологии. В дорожной карте проекта ФАД зафиксирован как расширение после достижения целевых показателей на паре <i>B₂ → 5'-ФМН</i>, чтобы не распылять экспериментальный ресурс команды на старте.",
        ]
        for t in tiers:
            story.append(P("• " + t, st_bullet))

        # ==================== SECTION 2 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec2", "2. Системный скрининг 16 альтернатив и обоснование выбора флавиновой цепи"))

        story.append(H2("sec2_1", "2.1. Сравнительная матрица классов молекул-кандидатов"))
        story.append(
            P(
                "Перед фиксацией флавинового каскада был проведён перекрёстный аудит 16 молекул-кандидатов из восьми биохимических классов по открытой научной литературе 2018–2026 гг., патентным базам ФИПС/ЕАПО/WIPO, таможенно-торговой статистике и реестру ГРЛС (Рисунок 1, Таблица 1).",
                st_body,
            )
        )

        img1 = Image(str(fig1_path), width=usable_w, height=usable_w * (4.35 / 10.8))
        story.append(
            KeepTogether(
                [
                    img1,
                    P(
                        "Рисунок 1. Позиционирование продуктовых платформ для цикла DBTL: (А) переход по каскаду от крупнотоннажного импортозависимого витамина B₂ к высокомаржинальным коферментам 5'-ФМН и ФАД; (Б) интегральная оценка по 5 обязательным критериям проекта.",
                        st_caption,
                    ),
                ]
            )
        )

        cand_data = [
            [
                P("Кандидат / связка", st_th),
                P("Старт Цикла 0 и метод Test", st_th),
                P("Новизна ML-DBTL и патенты РФ", st_th),
                P("Связь Шаг 1 → 2 и рынок РФ", st_th),
                P("Итог отбора", st_th),
            ],
            [
                P("<b>Каротиноиды:</b><br/>Ликопин → Астаксантин", st_td_bold),
                P("Вставка 3 чужеродных генов (<i>crtE/B/I</i>); экстракция ацетоном из клеток, 472 нм.", st_td),
                P("Десятки статей ML-DBTL (<i>BioAutomata</i> 2019, <i>Tu</i> 2026 [4, 11]); в промышленности доминирует гриб <i>B. trispora</i>.", st_td),
                P("Прямая цепь, но рынок ликопина мал, нет в Приказе № 392, жёсткая конкуренция с химсинтезом.", st_td),
                P("Отклонён (перегретое поле ML)", st_td),
            ],
            [
                P("<b>Ароматика:</b><br/><i>п</i>-Кумаровая к-та → Ресвератрол", st_td_bold),
                P("Вставка <i>TAL/PAL</i>; УФ 310 нм в среде.", st_td),
                P("Закрыто ML-DBTL на <i>S. cerevisiae</i> (<i>Moreno-Paz</i> 2024, 12,5–22,5 г/л).", st_td),
                P("Узкий рынок БАД в РФ, дешёвый растительный экстракт горца из КНР ($90–120/кг).", st_td),
                P("Отклонён (закрыто ML, малый рынок)", st_td),
            ],
            [
                P("<b>Биотин (Вит. B7)</b> и <b>Таурин</b>", st_td_bold),
                P("Сложный старт: суицидальный <i>BioB</i> (&lt;5 мг/л в дрожжах); у таурина нет УФ-хромофора.", st_td),
                P("Нет ML-DBTL, но крайне высокий риск нулевого сигнала на первых циклах.", st_td),
                P("Биотин есть в п. 34 № 392, но химсинтез таурина стоит $2,2/кг (тысячи т из этиленоксида).", st_td),
                P("Отклонён (риск провала старта / теста)", st_td),
            ],
            [
                P("<b>Прямой синтез (S,S)-SAMe</b> в дрожжах", st_td_bold),
                P("Нативный <i>SAM2</i> + подача L-Met; <b>нет оптического теста</b> (нужен ВЭЖХ на каждую лунку).", st_td),
                P("Сотни статей и плотный патентный ландшафт КНР на <i>S. cerevisiae / K. phaffii</i> (10–15 г/л).", st_td),
                P("Крупный фармрынок РФ (~14–16 млрд ₽ ГЛС), но <b>SAMe химически распадается за часы</b> при pH 6–7 и 30 °C.", st_td),
                P("Высокий риск аналитики и деградации на старте", st_td),
            ],
            [
                P("<b>L-Эрготионеин → (S,S)-SAMe</b>", st_td_bold),
                P("2 чужеродных гена (<i>NcEgt1, CpEgt2</i>); УФ 311 / 343 нм (нужна наработка фермента <i>TdETL</i>).", st_td),
                P("0 работ ML-DBTL в мире; нет патентов РФ на дрожжи (семья DTU WO2020221795 без фазы РФ).", st_td),
                P("Спрос EGT в РФ — единицы–десятки кг/год; <b>конфликт путей</b> (<i>↑HIS1</i> и <i>↑CYS4</i> вредят накоплению SAMe).", st_td),
                P("Рабочий резервный вариант (Платформа № 2)", st_td),
            ],
            [
                P("<b>Рибофлавин (B2) → 5'-ФМН → ФАД</b>", st_td_bold),
                P("<b>0 чужеродных генов</b> (свои <i>RIB1–7, FMN1, FAD1</i>); <b>10 с при 445 нм</b> на Шаге 1 + ВЭЖХ 445 нм на Шаге 2.", st_td),
                P("<b>0 работ ML-DBTL в дрожжах</b> и <b>0 работ по Transfer-DBTL</b> <i>B2 → 5'-ФМН</i> в мире; патенты ВКПМ истекли.", st_td),
                P("Прямой субстрат → продукт; <b>300–500 т/год импорта B2</b> (п. 32 МСХ № 392) + ниша чистого 5'-ФМН без <i>POCl₃</i>.", st_td),
                P("<b>Целевая платформа № 1 (выбрана)</b>", st_td_bold),
            ],
        ]
        story.append(styled_table(cand_data, [31 * mm, 35 * mm, 38 * mm, 44 * mm, 26 * mm]))
        story.append(P("Таблица 1. Сводная матрица отбора целевой продуктовой платформы среди проверенных классов соединений.", st_caption))

        story.append(H2("sec2_2", "2.2. Почему каскад B2 → 5'-ФМН снимает ограничения других кандидатов"))
        adv_bullets = [
            "<b>Термическая и кислотная устойчивость молекул в отличие от лабильных сульфониевых солей:</b> рибофлавин и 5'-ФМН устойчивы в кислых и нейтральных водных средах при температуре ферментации 28–30 °C (единственное простое правило работы в лаборатории — защита колб и пробирок от прямого яркого света фольгой или тёмным стеклом во избежание фотолиза до люмихрома [12, 16]). Это исключает спонтанную рацемизацию и распад проб в очереди на прибор.",
            "<b>Отсутствие метаболической инверсии между Шагом 1 и Шагом 2:</b> в каскаде <i>B₂ → 5'-ФМН</i> весь углеродно-азотный поток, выстроенный на первом шаге через пентозофосфатный путь, пуриновый ствол и ферменты <i>RIB1–RIB5</i>, целиком входит в состав молекулы второго шага. Модели не нужно разрушать предыдущий генотип — она донастраивает баланс адениловой и гуаниловой ветвей на том же штамме.",
            "<b>Чистое аналитическое окно в видимой области спектра (445 нм):</b> ни аминокислоты, ни нуклеотиды (АТФ, ГТФ, НАД⁺), ни сахара, ни соли минеральной среды вообще не поглощают свет при 445 нм и не флуоресцируют при 525 нм. Это даёт минимальный матричный фон как при экспресс-скрининге B₂ в планшете, так и на хроматограмме ВЭЖХ при раздельном измерении B₂, 5'-ФМН и ФАД.",
        ]
        for b in adv_bullets:
            story.append(P("• " + b, st_bullet))

        # ==================== SECTION 3 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec3", "3. Химический полусинтез (POCl3) против прямой микробной ферментации 5'-ФМН и ФАД"))

        story.append(H2("sec3_1", "3.1. Как устроена индустрия сегодня: ферментационный B2 и хлорокись фосфора"))
        story.append(
            P(
                "Ключевой проверочный вопрос для любого биотехнологического проекта: <i>«Не проще и не дешевле ли получать 5'-ФМН обычным химическим синтезом?»</i> Ответ требует строгого разделения двух стадий производства и двух разных квалификаций продукта, принятых в современной химико-фармацевтической промышленности [6, 12–14].",
                st_body,
            )
        )
        chem_status = [
            "<b>Стадия А (получение самого ядра витамина B₂) — 100% победа микробиологической ферментации:</b> исторический полный химический синтез рибофлавина из D-рибозы и 3,4-ксилидина требовал 6–9 многостадийных реакций с токсичными восстановителями и анилиновыми интермедиатами, уступая ферментации на 40–50% по себестоимости и экологической нагрузке [7]. К 2010 г. концерны <i>BASF</i> и <i>DSM</i> полностью закрыли химические производства B₂. Сегодня <b>100% мирового рибофлавина</b> (включая сырьё для синтеза ФМН) получают исключительно микробиологическим синтезом, а Государственная фармакопея РФ XV издания (ФС.2.1.0559) нормативно закрепляет: <i>«Субстанция представляет собой продукт ферментации»</i>.",
            "<b>Стадия Б (перевод B₂ в ФМН-Na) — доминирование неселективного химического полусинтеза с жёстким потолком чистоты:</b> основная масса стандартной коммерческой субстанции <i>Riboflavin Sodium Phosphate</i> (CAS 130-40-5; пищевой краситель E101(ii) и стандартный фарм-грейд USP / Ph. Eur.) сегодня производится путём обработки сухого ферментационного порошка B₂ <b>20-кратным мольным избытком хлорокиси фосфора (POCl₃)</b> или частично гидролизованной хлорфосфорной кислоты в органическом растворителе — пиридине, триалкилфосфате или γ-бутиролактоне — с последующей нейтрализацией NaOH и многократной перекристаллизацией через этаноламиновые/морфолиновые соли (патенты US 2610177, EP 0523582 A1 [13]).",
        ]
        for b in chem_status:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec3_2", "3.2. Четыре технологических предела химии POCl3 и ниша биосинтеза"))
        story.append(
            P(
                "Именно неселективность химического реагента <i>POCl₃</i> объясняет, почему ведущие мировые центры метаболической инженерии в 2020–2026 гг. (<i>Liu et al., 2020</i> [6]; <i>Fedorovych, Dmytruk &amp; Sibirny, 2023</i> [10]; <i>Jiménez-Nava et al., 2026</i> [12]; <i>Hou et al., J. Agric. Food Chem., февр. 2026</i> [11]) активно разрабатывают прямую ферментацию <i>de novo</i> (Рисунок 3, Таблица 2):",
                st_body,
            )
        )

        chem_limits = [
            "<b>1. Неизбежное образование 21–27% изомерных примесей в химическом продукте:</b> рибитильная цепь рибофлавина несёт четыре открытые гидроксильные группы (2'-OH, 3'-OH, 4'-OH и первичную 5'-OH). Хлорокись фосфора атакует их неселективно, образуя сложную смесь из целевого <b>5'-ФМН (лишь 73–79% в готовом фармакопейном порошке USP / Ph. Eur.!)</b>, неактивных региоизомеров <b>3'-ФМН и 4'-ФМН</b>, циклических 4',5'-фосфатов, рибофлавиндифосфатов (до 6%) и остаточного свободного B₂ (до 6%) [13, 14]. Фармакопейные статьи середины XX века (USP ≥73,5%, ЕФ 73,0–79,0%) были вынужденно написаны под этот технический предел химии, так как промышленное хроматографическое разделение региоизомеров 3'-, 4'- и 5'-ФМН экономически нерентабельно. Дрожжевая рибофлавинкиназа <i>FMN1</i> в воде при 30 °C фосфорилирует <b>исключительно первичный 5'-гидроксил со 100%-й региоселективностью</b>, обеспечивая после простой водной очистки чистоту <b>&gt;95–98%</b>.",
            "<b>2. Ингибирование ферментов примесями 3'-/4'-ФМН и спрос высокотехнологичных сегментов:</b> неприродные изомеры 3'- и 4'-ФМН связываются с апо-флавопротеинами как конкурентные ингибиторы, снижая удельную активность биокаталитических каскадов (флавин-зависимых монооксигеназ и оксидаз [17]) и искажая отклик диагностических биосенсоров. В офтальмологии в октябре 2025 г. FDA одобрило первый неинвазивный (<i>epithelium-on</i>) препарат кросслинкинга роговицы при кератоконусе <b>Epioxa™</b> на базе высокоочищенного 5'-ФМН [15], где отсутствие балластных дифосфатов и продуктов деградации критично для воспроизводимой генерации синглетного кислорода под УФ-А.",
            "<b>3. Двухзаводская импортозависимость и фосфатно-хлоридные стоки против одного водного биореактора:</b> даже несмотря на то, что в России химическое фосфорилирование B₂ освоило белгородское <b>ООО «Полисинтез»</b> (дочерний завод АФС компании «ПОЛИСАН» [29]), этот процесс: (а) <b>на 100% зависит от импорта зарубежного ферментационного порошка B₂</b>, которого в РФ нет; (б) требует коррозионностойкого химического цеха под 20-кратный избыток <i>POCl₃</i>, регенерации токсичного пиридина/γ-бутиролактона и утилизации высококонцентрированных солянокислых и фосфатных стоков. Прямая микробная ферментация <i>«отечественный сахар → B₂ → 5'-ФМН»</i> объединяет оба завода в <b>один стандартный водный ферментёр</b> без хлорорганики и одновременно решает проблему низкой растворимости B₂.",
            "<b>4. Полная невозможность химического синтеза следующего продукта цепи — кофермента ФАД:</b> селективно присоединить аденозинмонофосфат пирофосфатной связью в промышленном химическом реакторе невозможно, а попытка ферментативно достроить ФАД из дешёвого химического 75%-го ФМН блокируется примесями 3'-/4'-изомеров и требует дорогого экзогенного АТФ [6, 12]. Поэтому только штамм, генерирующий внутри клетки изомерно чистый 5'-ФМН и собственный пул АТФ, открывает прямой путь к биосинтезу ФАД ($1 200–3 000+/кг).",
        ]
        for b in chem_limits:
            story.append(P("• " + b, st_bullet))

        img3 = Image(str(fig3_path), width=usable_w, height=usable_w * (4.15 / 10.8))
        story.append(
            KeepTogether(
                [
                    Spacer(1, 1 * mm),
                    img3,
                    P(
                        "Рисунок 3. Сравнение химического полусинтеза (POCl₃) и прямой микробной ферментации: (А) распределение флавиновых фракций в готовой субстанции [10, 13, 14]; (Б) реальные торгуемые ценовые коридоры по данным Tridge (май–ноябрь 2026 г.) и каталогов GMP/IVD-поставщиков [26, 27].",
                        st_caption,
                    ),
                ]
            )
        )

        cb_data = [
            [
                P("Параметр сравнения", st_th),
                P("Химический полусинтез<br/>(импортный био-B2 + 20× POCl3)", st_th),
                P("Цельноклеточный биокатализ<br/>(покупной B2 + внешний АТФ)", st_th),
                P("Прямая ферментация <i>de novo</i><br/>(сахар → B2 → 5'-ФМН в дрожжах)", st_th),
            ],
            [
                P("<b>Исходное сырьё и суверенность в РФ</b>", st_td_bold),
                P("100% импорт ферментационного B2 + токсичные POCl3 и пиридин / γ-бутиролактон.", st_td),
                P("Импортный порошок B2 + дорогой кристаллический АТФ (или полифосфат).", st_td),
                P("<b>100% отечественное сырьё</b> (глюкоза, свекловичная меласса, молочная сыворотка, соли).", st_td),
            ],
            [
                P("<b>Содержание 5'-ФМН и примеси</b>", st_td_bold),
                P("<b>73–79%</b> 5'-ФМН + <b>21–27%</b> 3'-, 4'-изомеров, дифосфатов и остаточного B2 [13, 14].", st_td),
                P("90–97% 5'-ФМН, но загрязнение продуктами распада избытка АТФ/АДФ/АМФ [6].", st_td),
                P("<b>&gt;95–98% региоселективного 5'-ФМН</b> (при блоке фосфатаз и отделении от остатка B2).", st_td),
            ],
            [
                P("<b>Переход к синтезу ФАД</b>", st_td_bold),
                P("<b>Невозможен</b> химически; 75%-й сырцовый ФМН отравляет ФАД-синтетазу.", st_td),
                P("Возможен (<i>C. ammoniagenes</i>), но экономически зажат высокой ценой ФМН и АТФ [6, 12].", st_td),
                P("<b>Прямое расширение (+1 ген <i>FAD1</i>)</b> за счёт внутренней регенерации АТФ в клетке [11].", st_td),
            ],
            [
                P("<b>Оптовая цена продукта (2025–2026)</b>", st_td_bold),
                P("$63–120/кг (КНР/Индия USP) и $203–384/кг (DSM Ph.Eur.) за 75%-ю техсмесь [27].", st_td),
                P("Высокая себестоимость предшественников ограничивает масштаб реакторов [12].", st_td),
                P("Целевая ниша: <b>$250–600+/кг</b> (высокочистый 5'-ФМН АФС) и <b>$1 200–3 000/кг</b> (ФАД ≥95%).", st_td),
            ],
        ]
        story.append(styled_table(cb_data, [32 * mm, 47 * mm, 46 * mm, 49 * mm]))
        story.append(P("Таблица 2. Сопоставление трёх промышленных и биотехнологических маршрутов получения 5'-ФМН и ФАД.", st_caption))

        # ==================== SECTION 4 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec4", "4. Биохимия пути в дрожжах, штаммы БРЦ ВКПМ и референс «Цикла 0»"))

        story.append(H2("sec4_1", "4.1. Эндогенный путь дрожжей (0 чужеродных генов) против бактериального оперона"))
        story.append(
            P(
                "Фундаментальное преимущество выбранного каскада состоит в том, что все модельные и промышленные дрожжи (<i>Saccharomyces cerevisiae</i>, <i>Komagataella phaffii / Pichia pastoris</i>, <i>Yarrowia lipolytica</i>, <i>Candida famata</i>) изначально несут в собственных хромосомах <b>все восемь структурных генов</b> синтеза рибофлавина, 5'-ФМН и ФАД [8, 12, 18] (Рисунок 2). Кроме того, эукариотическая архитектура дрожжей даёт два прямых преимущества перед бактериями (<i>B. subtilis</i>, <i>E. coli</i>) именно при целевом синтезе 5'-ФМН [6, 12]:",
                st_body,
            )
        )
        euk_vs_prok = [
            "<b>Разделение ферментов Шага 2 и Шага 3 на два независимых гена (<i>FMN1</i> и <i>FAD1</i>):</b> у большинства бактерий рибофлавинкиназная и ФАД-синтетазная активности слиты в один бифункциональный белок (<i>ribC</i> у <i>B. subtilis</i>, <i>ribF</i> у <i>E. coli</i>), поэтому оверэкспрессия дикого гена неизбежно сбрасывает ФМН дальше в смесь ФМН + ФАД (в работе <i>Hou et al., 2026</i> на <i>E. coli</i> авторам пришлось сначала мутагенезом ломать ФАД-домен белка <i>RibF M-1</i>, а для синтеза ФАД всё равно заимствовать дрожжевые гены <i>FMN1</i> из <i>S. cerevisiae</i> и <i>FAD1</i> из <i>C. glabrata</i>! [11]). В дрожжах монофункциональная киназа <b><i>FMN1</i></b> и синтетаза <b><i>FAD1</i></b> изначально разнесены по разным генам, позволяя чисто остановить каскад на 5'-ФМН.",
            "<b>Отсутствие бактериального ФМН-рибопереключателя (<i>RFN-riboswitch</i>):</b> у <i>B. subtilis</i> накопление внутриклеточного ФМН напрямую связывается с 5'-UTR матричной РНК оперона <i>ribDEAHT</i> и обрубает синтез предшественника B₂ [6, 12]. В дрожжах такого рибопереключателя нет: замена нативных промоторов перед собственными генами <i>RIB</i> на конститутивные промоторы (<i>pGAP</i>, <i>pTEF1</i>) делает синтез B₂ нечувствительным к накоплению ФМН.",
        ]
        for b in euk_vs_prok:
            story.append(P("• " + b, st_bullet))

        img2 = Image(str(fig2_path), width=usable_w, height=usable_w * (5.75 / 11.0))
        story.append(
            KeepTogether(
                [
                    Spacer(1, 1 * mm),
                    img2,
                    P(
                        "Рисунок 2. Полная схема эндогенного биосинтеза B₂, 5'-ФМН и ФАД в дрожжах с 2D-структурами молекул (RDKit) и узлом бифуркации инозинмонофосфата (ИМФ → ГТФ против ИМФ → АТФ).",
                        st_caption,
                    ),
                ]
            )
        )

        story.append(
            P(
                "<b>Экспериментальное подтверждение простоты «Цикла 0» (<i>Marx, Mattanovich &amp; Sauer, Microb. Cell Fact., 2008</i> [8]):</b> в диких дрожжах на глюкозе выход B₂ составляет &lt;0,5 мг/л исключительно из-за транскрипционного ограничения первого фермента пути — ГТФ-циклогидролазы II (<i>RIB1</i>). В базовом лабораторном штамме <i>Komagataella phaffii</i> (<i>X-33 / GS115</i>) замена промотора <b>только перед одним собственным геном <i>PpRIB1</i></b> (амплифицированным обычным ПЦР с геномной ДНК самих дрожжей) на промотор <i>pGAP</i> сразу окрашивает колонии и среду в ярко-жёлтый цвет (<b>3,9–6,9 мг/л</b> в колбах и <b>50 мг/л</b> в супернатанте реактора); комбинация двух нативных генов <b><i>PpRIB1 + PpRIB3</i></b> даёт <b>10,1 мг/л</b> в колбах, а все шесть собственных генов <b><i>RIB1–RIB7</i></b> под <i>pGAP</i> обеспечивают <b>20,2 мг/л в колбах и 175 мг/л в биореакторе</b> ещё до какой-либо перестройки пуринового обмена.",
                st_body,
            )
        )

        story.append(H2("sec4_2", "4.2. Пять узких мест перехода B2 → 5'-ФМН → ФАД и парадокс развилки ИМФ"))
        story.append(
            P(
                "Почему высокое накопление B₂ не превращается в высокий выход 5'-ФМН автоматически и какую именно биологическую задачу решает конвейер DBTL? Анализ первичных работ 2014–2026 гг. выделяет пять конкретных регуляторных узлов:",
                st_body,
            )
        )
        bottlenecks = [
            "<b>1. Нелинейный конфликт ГТФ и АТФ на развилке инозинмонофосфата (<i>IMP</i>):</b> на сборку 1 моль B₂ расходуется 1 моль <b>ГТФ</b> (ветвь <i>IMP → XMP → GMP → GTP</i> через ферменты <i>IMD2/3, GUA1</i>). При инженерии чистого B₂ подавление соседней адениловой ветви (<i>ΔADE12</i>) даёт 5,4-кратный скачок титра B₂ (<i>Ruchala et al., июнь 2025</i> [19]). Однако конверсия <i>B₂ → 5'-ФМН → ФАД</i> требует затраты <b>1–2 моль свободного АТФ</b> на каждый моль флавина плюс АТФ на синтез PRPP (поэтому для наработки 2,31 г/л ФАД в работе <i>Hou et al., февр. 2026</i> адениловую ветвь <i>purA / ADE12</i> пришлось, наоборот, усиливать [11]!). Поиск оптимального соотношения проводимости плеч <i>IMD2/3 : ADE12</i> при сохранном росте клеток — центральная задача для алгоритма оптимизации.",
            "<b>2. Холостой цикл неспецифических клеточных фосфатаз:</b> в стационарной фазе роста вакуолярные и периплазматические фосфатазы дрожжей (<i>PHO8</i>, <i>PHO5</i> и HAD-гидролазы) отщепляют фосфатную группу от 5'-ФМН обратно до свободного B₂ (в классических ферментациях <i>C. famata</i> и <i>E. coli</i> это вызывало падение пика ФМН/ФАД после 40 ч с вторичным выбросом B₂ [6, 11]). Решения в дизайне штаммов и среды: ослабление/делеция <i>PHO8</i> и поддержание репрессирующего пула неорганического фосфата.",
            "<b>3. Мембранная локализация <i>FMN1</i> и белковое связывание кофакторов:</b> нативный белок Fmn1p у дрожжей заякорен на внешней мембране митохондрий и микросом, а Fad1p работает в цитозоле/митохондриях [20]; при этом значительная доля ФАД прочно связывается с клеточными флавопротеинами. Тестирование растворимого цитозольного варианта <i>FMN1</i> (без трансмембранного якоря или ортолога <i>D. hansenii / C. glabrata</i>) снимает пространственное ограничение.",
            "<b>4. Асимметрия мембранного транспорта B₂ и заряженного 5'-ФМН:</b> нейтральный B₂ легко выводится экспортёром <i>RFE1</i> (и пермеазами семейства <i>MCH5</i>), тогда как анионный 5'-ФМН частично задерживается внутри клетки, требуя раздельного учёта внутриклеточной и внеклеточной фракций и проверки транспортных кассет.",
            "<b>5. Чувствительность к ионам железа (Fe²⁺) и отрицательный эпистаз ПФП/пуринов:</b> у флавиногенных дрожжей избыток железа подавляет транскрипционный фактор <i>SEF1</i> [10, 12], а прямое суммирование нокаутов пентозофосфатного и пуринового путей проявляет сильный отрицательный эпистаз (<i>Li et al., июль 2026</i> [21]), обесценивая наивный ручной подход «вставить все сильные промоторы сразу».",
        ]
        for b in bottlenecks:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec4_3", "4.3. Выбор штамма-шасси и каталог доступных референсов (включая БРЦ ВКПМ)"))
        story.append(
            callout_box(
                "<b>Важное практическое ограничение по штаммам <i>Candida famata</i>:</b> рекордные лабораторные штаммы <i>C. famata</i> (<i>AF-4, BRP, BRPI</i> — до 16,4 г/л B₂ и 540 мг/л ФМН [9, 10]) созданы в Институте биологии клетки НАН Украины (Львов) и Жешувском университете (Польша) и физически недоступны для заказа в РФ. Дикие же коллекционные изоляты <i>C. famata</i> (исторически смешивавшиеся с <i>Debaryomyces subglobosus / C. flareri</i> [12]) относятся к CTG-кладе с нестандартным генетическим кодом, бедным набором CRISPR-векторов и известной генетической нестабильностью (из-за которой <i>ADM</i> закрыла завод на <i>C. famata</i> [6]). Поэтому <b>стартовать мокрую программу следует на стандартных доступных шасси <i>K. phaffii</i> (<i>X-33 / GS115</i>) или <i>S. cerevisiae</i> (<i>CEN.PK</i>)</b>, параллельно взяв отечественный штамм <i>B. subtilis</i> из БРЦ ВКПМ за 3 000 ₽ как готовый граммовый контроль аналитики.",
                bg=WARN_SOFT,
                border=WARN_BORDER,
            )
        )
        story.append(Spacer(1, 2 * mm))

        strains_data = [
            [
                P("Группа / коллекция", st_th),
                P("Конкретные штаммы", st_th),
                P("Генотип и ключевые модификации", st_th),
                P("Опубликованный титр", st_th),
                P("Доступность и статус в РФ", st_th),
            ],
            [
                P("<b>Рабочие дрожжевые шасси для старта Цикла 0</b>", st_td_bold),
                P("<i>Komagataella phaffii</i> (<i>X-33, GS115</i>) [8, 22]; <i>S. cerevisiae</i> (<i>CEN.PK113-7D</i>); <i>Y. lipolytica</i> [23]", st_td),
                P("<b>0 чужеродных генов:</b> ПЦР собственных <i>RIB1 + RIB3</i> или <i>RIB1–RIB7</i> под <i>pGAP / pTEF1</i>; затем подключение нативного <i>FMN1</i> и <i>Δpho8</i>.", st_td),
                P("<b>4–20,2 мг/л</b> B₂ в колбах и <b>50–175 мг/л</b> в реакторе уже на <i>pGAP–RIB</i> до инженерии пуринов [8].", st_td),
                P("Штаммы и плазмиды есть в лабораториях РФ и ВКПМ; полные модели <i>iKp1260 / ecYeast8</i>; <b>0 патентных барьеров РФ</b>.", st_td),
            ],
            [
                P("<b>Природные флавиногенные дрожжи (сравнение в Опыте 1)</b>", st_td_bold),
                P("<i>Meyerozyma guilliermondii</i> (ВКПМ / ВКМ); <i>Ogataea polymorpha</i> <b>ВКПМ Y-201</b>; <i>C. famata</i> (<i>AF-4, BRPI</i> — лит. ориентир [9, 10])", st_td),
                P("Природная дерегуляция синтеза B₂ при низком Fe²⁺; у <i>ВКПМ Y-201</i> — высокая нативная активность киназы <i>FMN1</i>.", st_td),
                P("Дикий <i>M. guilliermondii</i>: 8–34 мг/л колбы, 318 мг/л реактор; лит. рекорд <i>C. famata</i>: 16,4 г/л B₂, 0,54 г/л ФМН.", st_td),
                P("Дикие изоляты и <i>ВКПМ Y-201</i> доступны в Москве; требуют проверки паспорта штамма и маркеров.", st_td),
            ],
            [
                P("<b>Отечественные эталоны <i>B. subtilis</i> в фонде БРЦ ВКПМ</b>", st_td_bold),
                P("<b>ВКПМ B-2694</b> (<i>304/pMX45</i>); <b>ВКПМ B-6797</b> (<i>RU 2081906</i>); <b>ВКПМ B-8386 / B-8387</b> (<i>RU 2261273</i>, Миронов, Перумов, Дебабов)", st_td),
                P("Мутации <i>ribC / ribO</i> (снятие блока рибопереключателя, устойчивость к розеофлавину) + плазмида <i>pMX45</i> с <i>rib</i>-опероном.", st_td),
                P("<i>B-2694</i>: <b>3,5–4,5 г/л</b> B₂; <i>B-6797</i>: <b>12,4 г/л</b>; <i>B-8386/87</i>: <b>8,8–9,6 г/л</b> колбы и <b>18–21 г/л</b> реактор (70 ч).", st_td),
                P("<b>Физически лежат в Курчатовском ин-те</b> (<i>B-2694</i> в открытом прайсе за <b>3 000 ₽</b>); все патенты СССР/РФ <b>истекли</b>.", st_td),
            ],
        ]
        story.append(styled_table(strains_data, [30 * mm, 35 * mm, 41 * mm, 34 * mm, 34 * mm]))
        story.append(P("Таблица 3. Реестр штаммов для старта экспериментальной программы, сравнительных контролей и калибровки аналитики.", st_caption))

        # ==================== SECTION 5 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec5", "5. Аудит мировой литературы по DBTL/ML и математическая постановка Transfer-DBTL"))

        story.append(H2("sec5_1", "5.1. Что уже сделано в науке к ноябрю 2026 г. и где открытая зона новизны"))
        story.append(
            P(
                "Для корректного планирования публикации уровня Q1 и патентной стратегии необходимо чётко отделить уже опубликованные мировые результаты от открытой исследовательской ниши (Таблица 4). Само по себе «получение ФМН или ФАД в дрожжах» не является мировой новизной (базовые штаммы <i>C. famata</i> опубликованы в 2009–2023 гг. [10]). Реальная новизна лежит на пересечении <b>замкнутого цикла активного обучения на дрожжах</b> и <b>количественного переноса модели (Transfer Learning) по каскаду <i>B₂ → 5'-ФМН</i></b>:",
                st_body,
            )
        )

        lit_data = [
            [
                P("Год, авторы и объект", st_th),
                P("Что именно сделано и достигнутый результат", st_th),
                P("Отличие от нашей задачи (что остаётся открытым)", st_th),
            ],
            [
                P("<b>2008, Marx et al.</b><br/><i>K. phaffii X-33</i> → B₂ [8]", st_td_bold),
                P("Ручная пошаговая замена промоторов <i>pGAP</i> перед нативными <i>PpRIB1..RIB5</i>: 20,2 мг/л в колбах и 175 мг/л в биореакторе.", st_td),
                P("Классический ручной дизайн (OFAT) без инженерии развилки пуринов, без ML и без перехода к 5'-ФМН.", st_td),
            ],
            [
                P("<b>2009–2024, Dmytruk, Ruchala, Sibirny</b><br/><i>C. famata</i> → B₂, ФМН, ФАД [9, 10]", st_td_bold),
                P("До 16,4 г/л B₂ в 7-л реакторе (2014); до 231 мг/л ФМН в ферментёре и 540 мг/л ФМН в отмытых клетках на сыворотке (2023); 451 мг/л ФАД (2014).", st_td),
                P("Исключительно эмпирическая сборка на труднодоступном CTG-штамме; значительная примесь свободного B₂ и распад ФМН/ФАД после 40 ч; <b>0 работ с ML-DBTL</b>.", st_td),
            ],
            [
                P("<b>2018, Wagner et al.</b><br/><i>Y. lipolytica</i> → B₂ [23]", st_td_bold),
                P("Микрокапельный FACS-сортинг случайных УФ/EMS-мутантов по автофлуоресценции B₂.", st_td),
                P("Слепой мутагенез без известной матрицы «генотип → фенотип» и без обучаемой предсказательной модели.", st_td),
            ],
            [
                P("<b>Февр. 2026, Hou et al.</b><br/><i>E. coli</i> → ФМН и ФАД [11]", st_td_bold),
                P("Мутант <i>RibF M-1</i> + регенерация АТФ (<i>Δndh, ΔappB</i>) → <b>3,46 г/л ФМН</b>; гетерологичные дрожжевые <i>ScFMN1 + CgFAD1</i> + <i>ΔushA + ↑purA</i> → <b>2,31 г/л ФАД</b>.", st_td),
                P("Главный бактериальный ориентир 2026 г., но построен вручную (эндотоксиновое шасси <i>E. coli</i>, разные штаммы под ФМН и ФАД, без цикла ML).", st_td),
            ],
            [
                P("<b>Июль 2026, Li et al.</b><br/><i>B. subtilis</i> → B₂ [21]", st_td_bold),
                P("Стехиометрический расчёт <i>etiBsu1209</i> (OptKnock) + скрининг случайной библиотеки тройных <i>CRISPRi</i> в 96 лунках (<i>↓guaC + ↓ilvD</i>, 20,2 г/л).", st_td),
                P("Нет суррогатной модели Learn по титру; при этом <b>экспериментально доказан сильный отрицательный эпистаз</b> одиночных нокаутов.", st_td),
            ],
        ]
        story.append(styled_table(lit_data, [38 * mm, 70 * mm, 66 * mm]))
        story.append(P("Таблица 4. Систематический аудит релевантных работ по биосинтезу флавинов и смежным методам скрининга (2008–2026 гг.).", st_caption))

        story.append(H2("sec5_2", "5.2. Пространство факторов, гибрид GEM + GP и защита от утечки данных в ML"))
        story.append(
            P(
                "Чтобы вычислительная часть проекта имела строгую доказательную силу для журнала Q1 и не сводилась к формальной аппроксимации шума, архитектура модели строится по четырём методологическим правилам (Рисунок 4):",
                st_body,
            )
        )
        ml_rules = [
            "<b>1. Компактное параметризованное пространство дизайна (Design Space):</b> каждый штамм кодируется вектором дискретных уровней силы промоторов (низкий / средний / сильный) по 6 функциональным модулям: <i>x₁</i> — вход в ПФП (<i>ZWF1 / RPE1</i>); <i>x₂</i> — снятие ретроингибирования пуринового ствола (<i>ADE4_fbr</i>); <i>x₃, x₄</i> — соотношение плеч развилки ИМФ (<i>IMD2/3</i> на ГТФ против <i>ADE12</i> на АТФ); <i>x₅</i> — сборка флавинового ядра (<i>RIB1 + RIB3</i>); <i>x₆</i> — модуль конверсии и сохранения фосфата (<i>FMN1</i> + <i>Δpho8</i> + экспортёр <i>RFE1</i>), дополненным управляемыми параметрами среды (концентрация Fe²⁺, фосфата и режим подачи субстрата).",
            "<b>2. Гибридизация со стехиометрической моделью (Physics-Informed Features):</b> до посадки штаммов для каждого вектора генотипа в ферментно-ограниченной модели дрожжей (<i>ecYeast8</i> для <i>S. cerevisiae</i> или <i>iKp1260</i> для <i>K. phaffii</i>) вычисляются расчётные верхние границы потоков <i>v_GTP</i>, <i>v_ATP</i>, <i>v_Ru5P</i> и теневые цены кофакторов (алгоритмы pFBA / FSEOF). Это даёт суррогатной модели физически осмысленные признаки на малой выборке (<i>N = 24–48</i> уникальных генотипов в раунде).",
            "<b>3. Адекватный класс моделей для малых данных и строгий учёт единиц наблюдения:</b> при объёме серии в 24–36 генотипов × 2 среды × 3 независимых биологических повтора (= 144–216 культур) глубокие нейросети заведомо переобучаются; базовыми алгоритмами выступают <b>гауссовские процессы (Gaussian Process с ядром Матерна для байесовской оптимизации UCB / Expected Improvement)</b>, градиентный бустинг / Random Forest и регуляризованная линейная модель с парными взаимодействиями. При кросс-валидации <b>все повторы одного генотипа и одной партии среды строго группируются (GroupKFold по strain_id и batch_id)</b>, исключая утечку данных.",
            "<b>4. Проспективная проверка гипотезы Transfer-DBTL и контрольные стратегии:</b> главная проверяемая гипотеза для статьи — <i>позволяет ли суррогатная модель, предобученная на быстрых оптических данных Шага 1 (B₂, 445 нм), при переносе на задачу максимизации селективного выхода 5'-ФМН достичь более высокой объёмной производительности за фиксированный бюджет новых сборок (например, 16 новых штаммов), чем (а) модель, обученная только на данных ФМН с нуля, (б) случайный выбор из того же пула и (в) заранее зафиксированное экспертное правило биолога?</i> Все прогнозы замораживаются до постановки проверочного мокрого эксперимента.",
        ]
        for b in ml_rules:
            story.append(P("• " + b, st_bullet))

        img4 = Image(str(fig4_path), width=usable_w, height=usable_w * (4.45 / 10.8))
        story.append(
            KeepTogether(
                [
                    Spacer(1, 1 * mm),
                    img4,
                    P(
                        "Рисунок 4. Схема двухэтапного конвейера Transfer-DBTL с входным аналитическо-диагностическим шлюзом (Опыты 0–2) и проспективным сравнением стратегий отбора.",
                        st_caption,
                    ),
                ]
            )
        )

        # ==================== SECTION 6 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec6", "6. Двухуровневая аналитика стадии Test, контрольные опыты (0–6) и дорожная карта"))

        story.append(H2("sec6_1", "6.1. Оптический скрининг 445 нм, хроматография RP-C18 и обязательные контроли"))
        story.append(
            P(
                "Главная методическая ошибка при работе с флавинами — отождествлять любой жёлтый или флуоресцирующий супернатант с конкретным веществом, поскольку B₂, 5'-ФМН и ФАД несут одинаковое изоаллоксазиновое кольцо [6, 16]. Чтобы совместить высокую пропускную способность и аналитическую строгость без масс-спектрометра, протокол стадии <b>Test</b> разделён на два уровня с фиксированным набором контролей:",
                st_body,
            )
        )
        test_levels = [
            "<b>Работа на прозрачной химически определённой минеральной среде (среда Делфта / Verduyn):</b> в среде отсутствуют дрожжевой экстракт и пептон, содержащие переменные количества природных флавинов, пуринов и окрашенных продуктов реакции Майяра. Это обнуляет оптический фон холостой среды при 445 нм и делает расчёт выхода на грамм сахара (<i>Y_P/S</i>) и материального баланса полностью прозрачным.",
            "<b>Уровень 1 — Экспресс-скрининг супернатанта Шага 1 (B₂) при 445 нм за 10 секунд:</b> для штаммов без оверэкспрессии <i>FMN1</i> (где доля внеклеточных фосфофлавинов пренебрежимо мала, что проверяется в Опыте 0) клетки осаждают центрифугированием (4 000 g, 5 мин), прозрачный супернатант разводят водой или фосфатным буфером pH 7,0 в 5–20 раз и измеряют поглощение при <b>445 нм</b> в стандартной кювете 1 см или 96-луночном планшете с пересчётом по калибровочной прямой стандарта B₂ в той же матрице (при 15 мг/л неразбавленный супернатант даёт A₄₄₅ ≈ 0,50 — в 50 раз выше порога прибора).",
            "<b>Уровень 2 — Раздельное количественное определение B₂, 5'-ФМН и ФАД на обычном ВЭЖХ-УФ/ФЛД (без масс-спектрометра):</b> с момента включения модуля <i>FMN1</i> центральным методом измерения служит стандартная обращённо-фазовая хроматография на колонке <b>RP-C18</b> (элюент: ацетатно-аммонийный или формиатный буфер pH 4,5–5,5 + метанол/ацетонитрил, время анализа <b>5–7 минут на пробу</b>, детекция по поглощению в видимой области <b>445 нм</b> или флуоресценции Ex 450 / Em 525 нм [6, 16]). В отличие от 254 нм, на длине волны 445 нм на хроматограмме видны исключительно три чистых разделённых пика (ФАД → 5'-ФМН → B₂) без помех от клеточных белков и нуклеотидов. <i>Примечание по дифференциальному флуоресцентному тесту Бесси–Лоури</i> (рост свечения в ~10 раз при кислотном гидролизе ФАД → ФМН): метод может использоваться как вспомогательный планшетный фильтр только после валидации на модельных смесях трёх стандартов против ВЭЖХ.",
            "<b>Раздельный учёт двух фракций и обязательные контроли каждой серии:</b> для всех ключевых вариантов анализируются как <b>внеклеточный супернатант</b>, так и <b>внутриклеточный экстракт</b> (горячая буферная/слабокислотная экстракция осадка клеток в темноте для высвобождения белково-связанного и удерживаемого ФМН/ФАД). В каждую серию включаются: (1) индивидуальные калибровки по чистым стандартам B₂, 5'-ФМН и ФАД; (2) холостая среда без клеток; (3) исходный родительский штамм; (4) контроль возврата добавленного стандарта (<i>spike-recovery</i> в реальной матрице); (5) контроль фотостабильности при пробоподготовке.",
        ]
        for b in test_levels:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec6_2", "6.2. Последовательность экспериментов от диагностических вопросов (Опыты 0–6)"))
        story.append(
            P(
                "Вместо немедленного заказа десятков генетических конструкций лабораторная программа выстроена как каскад коротких диагностических опытов, где каждый шаг отвечает на конкретный инженерный вопрос и защищает команду от неверной интерпретации данных:",
                st_body,
            )
        )
        exps = [
            "<b>Опыт 0 — Валидация аналитики и оборудования:</b> построение калибровочных кривых чистых стандартов B₂, 5'-ФМН и ФАД в воде и минеральной среде на спектрофотометре (445 нм) и ВЭЖХ RP-C18; проверка линейности, пределов обнаружения (LOD/LOQ), полноты внутриклеточной экстракции и отсутствия взаимопревращения флавинов на свету.",
            "<b>Опыт 1 — Сравнение доступных стартовых штаммов:</b> параллельное культивирование в идентичных условиях базового дрожжевого референса (<i>K. phaffii X-33 / S. cerevisiae CEN.PK</i> с нативной кассетой <i>RIB1 + RIB3</i>), доступного штамма из ВКПМ (<i>ВКПМ Y-201</i> / <i>M. guilliermondii</i>) и бактериального контроля <i>B. subtilis ВКПМ B-2694</i> с замером роста, расхода сахара и долей B₂/ФМН/ФАД внутри и вне клеток.",
            "<b>Опыт 2 — Диагностическое разделение синтеза B₂ и его фосфорилирования в 5'-ФМН:</b> подача готового экзогенного рибофлавина к клеткам с базовой и усиленной экспрессией <i>FMN1</i>. Этот быстрый опыт без перестройки всего генома однозначно показывает, что именно лимитирует выход 5'-ФМН: нехватка самого предшественника B₂, низкая активность киназы <i>FMN1</i>, дефицит АТФ или обратный распад под действием фосфатаз.",
            "<b>Опыт 3 — Причинная панель метаболических модулей (обучающая серия Шага 1):</b> сборка компактной факторной панели (24–36 верифицированных генотипов) по модулям ПФП (<i>ZWF1/RPE1</i>), пуриновой развилки (<i>ADE4, IMD2 vs ADE12</i>) и ядра <i>RIB1/RIB3</i>; сбор первой размеченной матрицы «генотип × среда → B₂, биомасса, выход на сахар».",
            "<b>Опыт 4 — Проверка робастности на технологически реалистичной среде:</b> тестирование лучших вариантов при варьировании уровня Fe²⁺, фосфата и на доступных отечественных углеродных субстратах (меласса / пермеат молочной сыворотки) для отсева хрупких лабораторных фенотипов.",
            "<b>Опыты 5 и 6 — Честный проспективный раунд DBTL и перенос знаний на 5'-ФМН в биореакторе:</b> заморозка прогнозов GP-модели до опыта, слепое сравнение активного обучения со случайным и экспертным выбором на <i>B₂</i>, а затем — проверка <b>Transfer-DBTL</b> при включении модуля <i>FMN1 + ↓PHO8</i> и масштабирование лидеров в управляемый лабораторный биореактор (1–5 л) с контролем pH, растворённого кислорода (pO₂) и материального баланса очистки.",
        ]
        for b in exps:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec6_3", "6.3. План 12-недельного пилота (MVP), двухлетний горизонт и критерии GO / NO-GO"))
        story.append(
            P(
                "Дорожная карта синхронизирует быстрый 12-недельный стартовый пилот (необходимый для демонстрации первых экспериментальных результатов и открытия основного финансирования) с полной 24-месячной научно-технологической программой под публикацию Q1 и патентование (Таблица 5, Таблица 6). <b>Правило работы с открытым репозиторием:</b> до подачи патентной заявки в публичный GitHub выкладываются только обзорные материалы и архитектура метода, тогда как конкретные новые комбинации промоторов, последовательности конструкций и первичные таблицы титров ведутся в закрытом контуре команды.",
                st_body,
            )
        )

        plan_data = [
            [
                P("Этап и горизонт", st_th),
                P("Ключевые работы команды", st_th),
                P("Проверяемый результат (Deliverable)", st_th),
                P("Стоп-сигнал (NO-GO / пересмотр)", st_th),
            ],
            [
                P("<b>Фаза 1 (MVP старт)</b><br/>Недели 1–4", st_td_bold),
                P("Заказ 3 стандартов и штамма <i>ВКПМ B-2694</i>; валидация 445 нм и ВЭЖХ RP-C18 (Опыт 0); ПЦР-сборка референса <i>pGAP/pTEF1–RIB1+RIB3</i> в дрожжах (Опыт 1).", st_td),
                P("Воспроизводимое разделение 3 флавинов на ВЭЖХ; жёлтый дрожжевой штамм Цикла 0 с выходом B₂ ≥10–20 мг/л в колбах.", st_td),
                P("Нет доступа к рабочему ВЭЖХ-УФ/ФЛД либо ошибка измерения превышает межштаммовую разницу.", st_td),
            ],
            [
                P("<b>Фаза 2 (MVP DBTL-1)</b><br/>Недели 5–12", st_td_bold),
                P("Диагностика конверсии B₂ → ФМН (Опыт 2); сборка первой серии 24–36 генотипов (Опыт 3); обучение GP + GEM и первый проверочный раунд прогнозов.", st_td),
                P("Размеченный датасет ≥150 независимых культур; подтверждённый рост титра/продуктивности B₂ и локализованный лимит 5'-ФМН.", st_td),
                P("Модель активного обучения на новых штаммах не превосходит случайный отбор после проверки генотипов.", st_td),
            ],
            [
                P("<b>Фаза 3 (5'-ФМН)</b><br/>Месяцы 4–12", st_td_bold),
                P("Включение <i>FMN1 + ↓PHO8 + RFE1</i> и балансировка развилки <i>IMD2 / ADE12</i>; проспективный эксперимент Transfer-DBTL (с данными B₂ против обучения с нуля).", st_td),
                P("Дрожжевой продуцент селективного 5'-ФМН (&gt;85–90% доли среди флавинов); доказанное сокращение числа опытов за счёт переноса модели.", st_td),
                P("Неустранимый гидролиз 5'-ФМН до B₂ или отсутствие выигрыша трансферной модели на равном бюджете.", st_td),
            ],
            [
                P("<b>Фаза 4 (Реактор, Q1 и патент)</b><br/>Месяцы 13–24", st_td_bold),
                P("Перенос лидеров в 1–5 л fed-batch биореактор на доступном сырье (Опыты 4, 6); водная очистка образца ≥95% 5'-ФМН; оценка расширения на <i>FAD1</i> (ФАД).", st_td),
                P("Заявка на патент РФ (конструкция/способ), рукопись статьи Q1, лабораторный образец 5'-ФМН с паспортом ВЭЖХ для индустриального партнёра.", st_td),
                P("Потеря продуктивности при переносе из колбы в реактор или на реальное отечественное сырьё.", st_td),
            ],
        ]
        story.append(styled_table(plan_data, [30 * mm, 55 * mm, 50 * mm, 39 * mm]))
        story.append(P("Таблица 5. Интегрированная дорожная карта: 12-недельный результативный пилот (Фазы 1–2) и программа масштабирования до 24 месяцев.", st_caption))

        # ==================== SECTION 7 (NEW PAGE) ====================
        story.append(PageBreak())
        story.append(H1("sec7", "7. Рынки мира и РФ, реестр ГРЛС, Приказ Минсельхоза № 392 и список источников"))

        story.append(H2("sec7_1", "7.1. Мировой рынок B2, ФМН-Na и ФАД: объёмы и реальные экспортные цены"))
        world_mkt = [
            "<b>Витамин B₂ (рибофлавин, код ТН ВЭД / HS 293623):</b> мировой объём производства составляет <b>~12 000–15 000 тонн/год</b> (~$420–500 млн), около 70% которого потребляет кормовая индустрия, а 30% — пищевая промышленность и фармацевтика. Промышленная ферментация в реакторах &gt;100 м³ (<i>Hubei Guangji</i>, <i>DSM-Firmenich</i>, <i>BASF</i>) работает при титрах 26–34 г/л, формируя оптовую цену $15–22/кг на кормовой 80%-й продукт и $38–46/кг на фарм-субстанцию 98% [6, 27].",
            "<b>Рибофлавин-5'-фосфат натрия (ФМН-Na, CAS 130-40-5):</b> по отраслевой оценке <i>Dataintelo</i> мировой рынок ФМН-Na составил <b>$408,6 млн в 2024 г.</b> с прогнозом роста до $671,2 млн к 2033 г. (CAGR 5,7%) [27]. По фактическим таможенным транзакциям платформы <i>Tridge</i> (обновление ноябрь 2026 г., 279 экспортных сделок) средние экспортные цены стандартного 75%-го химического сорта в мае 2026 г. составили: <b>Китай — $62,6–120,3/кг</b>, <b>Индия — $87,9–99,5/кг</b>, <b>Нидерланды — $203,0/кг</b>, <b>Швейцария (фарм-грейд DSM) — $277,7–383,9/кг</b>, тогда как высокоочищенный фармацевтическо-хроматографический сорт без примесных изомеров в малых партиях достигает <b>$1 500–3 800/кг</b> [27].",
            "<b>Флавинадениндинуклеотид (ФАД-Na₂, CAS 84366-81-4):</b> по отчёту <i>Global Info Research</i> (февраль 2026 г.) мировой выпуск ФАД составил <b>~620 тонн ($74,1 млн в 2025 г., прогноз $107 млн к 2032 г.)</b> [26]: от $110–300/кг за технические/кормовые концентраты до <b>$1 200–3 000+/кг</b> за очищенную субстанцию ≥95–99% для офтальмологических капель и диагностических тест-полосок глюкометров (<i>FAD-GDH</i>).",
        ]
        for b in world_mkt:
            story.append(P("• " + b, st_bullet))

        story.append(H2("sec7_2", "7.2. Спрос в РФ, структура ГРЛС («Цитофлавин» / «Полисинтез»), господдержка и ТЭО"))
        ru_mkt = [
            "<b>Фундамент Шага 1 в РФ (300–500 т/год импорта B₂ и п. 32 Приказа № 392):</b> при выпуске в России <b>45,1 млн т комбикормов и 586,6 тыс. т премиксов</b> (2025 г.) отечественное ферментационное производство витамина B₂ равно <b>0%</b>. Пункт <b>32 Приказа Минсельхоза России от 06.06.2025 № 392</b> дословно включает <i>«Витамин B2 (рибофлавин), его соединения и производные»</i> [25] — это даёт юридически безупречное покрытие сразу для всех трёх молекул проекта (B₂, 5'-ФМН и ФАД) по программам субсидирования НИОКР и возмещения до 20% CAPEX.",
            "<b>Верифицированная картина фармрынка РФ по выгрузке ГРЛС Минздрава РФ (ноябрь 2026 г.) [28, 29]:</b> крупнейшим потребителем рибофлавина фосфата натрия в России выступает оригинальный рецептурный препарат перечня ЖНВЛП <b>«Цитофлавин®»</b> (ООО «НТФФ «ПОЛИСАН», Санкт-Петербург; объём продаж >4–5 млн упаковок и >4 млрд руб./год в розничном и госпитальном сегментах, а также дженерик «Флатоцин» компании <i>Solopharm</i>), инъекционный «Рибофлавин-мононуклеотид» («Фармстандарт-УфаВИТА») и медицинские изделия кросслинкинга роговицы («Декстралинк», «Риболинк» Уфимского НИИ глазных болезней). По официальной карточке РУ «Цитофлавина» (ЛП-№(000973)-(РГ-RU)) его субстанцию поставляют французское подразделение <b>DSM Nutritional Products France</b> (П N013263/01), индийская <b>Supriya Lifescience</b> и дочерний химический завод «Полисана» в Белгороде — <b>ООО «Полисинтез»</b> (ЛСР-002295/07), кроме которых в ГРЛС зарегистрированы индийские <b>Harman Finochem</b> (ЛСР-007577/08) и <b>Danashmand Organic</b> (ЛС-002018).",
            "<b>Кто потенциальный индустриальный партнёр в РФ и в чём наше предложение ценности:</b> поскольку белгородский «Полисинтез» (выпускающий 265–285 т/год химических АФС [29]) опирается исключительно на химический синтез из импортного B₂, проект DBTL не конкурирует в лоб с дешёвым китайским 75%-м техпродуктом, а предлагает партнёрам уровня <b>НТФФ «ПОЛИСАН» / ООО «Полисинтез»</b>, <b>Solopharm («Гротекс»)</b> и производителям офтальмологических растворов и диагностических реагентов: (1) <b>отечественную одностадийную биотехнологию изомерно чистого 5'-ФМН (&gt;95–98%)</b> с отечественного сахара без <i>POCl₃</i> и без импорта B₂; (2) выход в полностью импортируемую нишу кофермента <b>ФАД</b>.",
            "<b>Иллюстративный материально-экономический ориентир пилотного масштаба (не обещание до валидации ТЭО):</b> при достижении в управляемом биореакторе рабочего объёма <b>1 м³</b> умеренного целевого титра <b>5 г/л 5'-ФМН</b> (на уровне уже показанных в 2026 г. 3,5 г/л в <i>E. coli</i> [11]), 45 циклах в год и сквозном выходе водной ионообменной/мембранной очистки 70% годовая выработка составит <b>~157,5 кг изомерно чистого 5'-ФМН</b> (эквивалент закрытия годовой потребности десятков тысяч упаковок офтальмологических и инъекционных серий или <b>$50–95 тыс.</b> выручки по цене высокочистой субстанции; при масштабировании в ферментёр <b>10 м³</b> — <b>~1,58 т/год</b>, или <b>$0,5–0,95 млн/год</b>, с возможностью переключения той же линии на <b>ФАД</b> по $1 200–2 500/кг). Точный расчёт CAPEX/OPEX и порога безубыточности выполняется на Фазе 3 по фактическим показателям титра, скорости и извлечения.",
        ]
        for b in ru_mkt:
            story.append(P("• " + b, st_bullet))

        kpi_data = [
            [
                P("Критерий проверки", st_th),
                P("Статус выполнения в платформе «Рибофлавин (B2) → 5'-ФМН (→ ФАД)»", st_th),
            ],
            [
                P("<b>1. Простой старт Цикла 0</b>", st_td_bold),
                P("<b>Выполнен:</b> 0 чужеродных генов; ПЦР своих <i>RIB1+RIB3</i> в <i>K. phaffii / S. cerevisiae</i> даёт 10–20 мг/л жёлтого B₂ + эталон <i>ВКПМ B-2694</i> за 3 000 ₽.", st_td),
            ],
            [
                P("<b>2. Быстрый тест без масс-спектрометра</b>", st_td_bold),
                P("<b>Выполнен:</b> 10 с в супернатанте при 445 нм без реагентов на Шаге 1 + 5–7 мин ВЭЖХ RP-C18 (445 нм) для разделения B₂, 5'-ФМН и ФАД на Шаге 2.", st_td),
            ],
            [
                P("<b>3. Научная новизна Q1 и патент РФ</b>", st_td_bold),
                P("<b>Выполнен:</b> 0 работ с ML-DBTL на флавинах в дрожжах и 0 работ в мире по Transfer-DBTL <i>B₂ → 5'-ФМН</i> на развилке <i>ИМФ (ГТФ/АТФ)</i>; патенты ВКПМ истекли.", st_td),
            ],
            [
                P("<b>4. Превосходство над химсинтезом</b>", st_td_bold),
                P("<b>Чётко разграничено:</b> химия <i>POCl₃</i> даёт дешёвую 75%-ю смесь изомеров из импортного B₂, но не способна дать дешёвый &gt;95% 5'-ФМН из сахара РФ и вообще не применима к ФАД.", st_td),
            ],
            [
                P("<b>5. Спрос в РФ и господдержка</b>", st_td_bold),
                P("<b>Выполнен напрямую:</b> п. 32 Приказа Минсельхоза № 392 покрывает весь каскад; 300–500 т/год импорта B₂ + рынок АФС для «Цитофлавина» (ЖНВЛП), офтальмологии и IVD.", st_td),
            ],
        ]
        story.append(Spacer(1, 1.5 * mm))
        story.append(styled_table(kpi_data, [44 * mm, 130 * mm]))
        story.append(P("Таблица 6. Итоговая верификация выбранной платформы по всем пяти критериям запуска проекта DBTL.", st_caption))

        story.append(H2("sec7_3", "7.3. Список проверяемых первоисточников, патентов и реестровых записей"))
        refs = [
            "[1] Nielsen J., Keasling J. D. Engineering cellular metabolism. <i>Cell</i>. 2016;164(6):1185–1197. doi:10.1016/j.cell.2016.02.004.",
            "[2] Lawson C. E., Martí J. M., Radivojevic T. et al. Machine learning for metabolic engineering: a review. <i>Metabolic Engineering</i>. 2021;63:34–60. doi:10.1016/j.ymben.2020.10.005.",
            "[3] Carbonell P., Jervis A. J., Robinson C. J. et al. An automated Design-Build-Test-Learn pipeline for enhanced microbial production of fine chemicals. <i>Communications Biology</i>. 2018;1:66. doi:10.1038/s42003-018-0076-9.",
            "[4] HamediRad M., Chao R., Weisberg S. et al. Towards a fully automated algorithm driven platform for biosystems design (BioAutomata). <i>Nature Communications</i>. 2019;10:5150. doi:10.1038/s41467-019-13189-z.",
            "[5] Radivojević T., Costello Z., Workman K., Garcia Martin H. A machine learning Automated Recommendation Tool for synthetic biology. <i>Nature Communications</i>. 2020;11:4879. doi:10.1038/s41467-020-18008-4.",
            "[6] Liu S., Hu W., Wang Z., Chen T. Production of riboflavin and related cofactors by biotechnological processes. <i>Microbial Cell Factories</i>. 2020;19:31 (PMC7017516). doi:10.1186/s12934-020-01302-7.",
            "[7] Bretzel W., Schurter W., Ludwig B. et al. Commercial riboflavin production by recombinant <i>Bacillus subtilis</i>: comparison of fermentation and chemical synthesis. <i>J. Ind. Microbiol. Biotechnol.</i> 1999;22(1):19–26; Zhang J.-R. et al. Biotechnological strategies of riboflavin biosynthesis in microbes. <i>Engineering</i>. 2022. doi:10.1016/j.eng.2021.03.018.",
            "[8] Marx H., Mattanovich D., Sauer M. Overexpression of the riboflavin biosynthetic pathway in <i>Pichia pastoris</i> (<i>Komagataella phaffii</i>). <i>Microbial Cell Factories</i>. 2008;7:23. doi:10.1186/1475-2859-7-23.",
            "[9] Dmytruk K. V., Yatsyshyn V. Y., Sybirna N. O., Fedorovych D. V., Sibirny A. A. Metabolic engineering and classic selection of the yeast <i>Candida famata</i> for construction of strains with enhanced riboflavin production (16.4 g/L in 7-L bioreactor). <i>Metabolic Engineering</i>. 2011;13(1):82–88; <i>J. Biotechnol.</i> 2014;172:11–17. doi:10.1016/j.jbiotec.2013.12.005.",
            "[10] Yatsyshyn V. Y., Ishchuk O. P., Voronovsky A. Y., Fedorovych D. V., Sibirny A. A. Production of flavin mononucleotide by metabolically engineered yeast <i>Candida famata</i>. <i>Metabolic Engineering</i>. 2009;11(3):163–167; Fedorovych D. V., Dmytruk K. V., Sibirny A. A. Construction of advanced FMN producers in <i>Candida famata</i> (540 mg/L on cheese whey). <i>Yeast</i>. 2023;40(8):360–366. doi:10.1002/yea.3843; Yatsyshyn V. Y. et al. FAD production in <i>C. famata</i> (451 mg/L). <i>Biochem. Eng. J.</i> 2014 (PMID:24595668).",
            "[11] Hou J., Hu W., Ding M., Liu C., Zou S., Wang Z., Chen T. Metabolic engineering of <i>Escherichia coli</i> for efficient production of redox cofactors flavin mononucleotide (3457.3 mg/L) and flavin adenine dinucleotide (2305.2 mg/L). <i>J. Agric. Food Chem.</i> 2026;74(7):6321–6332. doi:10.1021/acs.jafc.5c07311.",
            "[12] Jiménez-Nava R. A., Chávez-Camarillo G. M., Cristiani-Urbina E. Biosynthesis, regulation, and biotechnological production strategies of riboflavin and its derivatives: a review. <i>Pharmaceuticals</i>. 2026;19(3):389 (PMC13028884). doi:10.3390/ph19030389.",
            "[13] Патенты химического фосфорилирования рибофлавина хлорокисью фосфора и очистки от изомерных примесей: US 2610177 (20-кратный избыток POCl₃) и EP 0523582 A1 (<i>Improved process for preparing riboflavin 5'-phosphate ester monosodium salt</i>); Krauss U. et al. Cofactor trapping to produce highly pure FMN (&gt;95% vs POCl₃ γ-butyrolactone process). <i>Appl. Environ. Microbiol.</i> 2010;76(19):6576–6582. doi:10.1128/AEM.01541-10.",
            "[14] Фармакопейные спецификации <i>Riboflavin Sodium Phosphate</i> (CAS 130-40-5): USP / Ph. Eur. / BP (содержание основного 5'-монофосфата ≥73,5% / 73,0–79,0% на сухое вещество; свободный рибофлавин ≤6,0%; рибофлавиндифосфаты ≤6,0%; прочие региоизомеры 3'- и 4'-ФМН); Государственная фармакопея РФ XV изд. ФС.2.1.0559 «Рибофлавин».",
            "[15] Одобрение FDA препарата <b>Epioxa™</b> (topical ophthalmic riboflavin 5'-phosphate, Oct 2025, UNII 7N464URE7E) для трансэпителиального (<i>epithelium-on</i>) кросслинкинга роговицы при кератоконусе (PubChem CID 643976 / DrugBank, 2025–2026).",
            "[16] Bosch A. M., van Dijk M., Goorden S. M. I. Riboflavin (B2) and FAD/FMN metabolites: chromatographic analysis. <i>Laboratory Guide to the Methods in Biochemical Genetics</i>. Springer, 2024. doi:10.1007/978-3-031-58819-8_23; Bessey O. A., Lowry O. H., Love R. H. Fluorometric measurement of the nucleotides of riboflavin. <i>J. Biol. Chem.</i> 1949;180(2):755–769.",
            "[17] Li Y. et al. Flavin-dependent monooxygenases as versatile biocatalysts in biomanufacturing. <i>Biotechnol. Adv. / PMC13173928</i>. 2026; Santema L. L., Fraaije M. W. Activity assays for flavoprotein oxidases. <i>Appl. Microbiol. Biotechnol.</i> 2025;109:115.",
            "[18] Saccharomyces Genome Database (SGD) Pathway: <i>riboflavin, FMN and FAD biosynthesis (YEAST-RIBOSYN-PWY: RIB1, RIB7, RIB2, RIB3, RIB4, RIB5, FMN1, FAD1)</i>.",
            "[19] Ruchala J., Dmytruk K. V., Sibirny A. A. et al. Blocking AMP branch (<i>ΔADE12</i>) enhances riboflavin synthesis 5.4-fold in engineered yeast; utilization of lignocellulose hydrolysate and cheese whey. <i>Microbial Cell Factories</i> / <i>PubMed 39009031</i> (2022–2025).",
            "[20] Giancaspero T. A., Barile M. et al. Subcellular localization of Fad1p and Fmn1p in <i>Saccharomyces cerevisiae</i>: implications for flavin homeostasis and export. <i>Int. J. Mol. Sci. / PMID 34575116</i> (2021).",
            "[21] Bi X. et al. Multi-omics whole-cell model <i>EM_iBsu1209-ME</i> in <i>B. subtilis</i>. <i>Advanced Science</i>. 2024; Li Y. et al. Genome-scale modeling <i>etiBsu1209</i> and combinatorial <i>CRISPRi</i> (<i>dCpf1</i>) reveal negative epistasis in riboflavin overproduction (20.2 g/L). <i>ENGINEERING Foods</i>. 2026 (Jul 2026); Xu Y. et al. <i>dhfCas13X</i> RNA switches on <i>ribDEAHT</i>. <i>Nucleic Acids Res.</i> 2025.",
            "[22] Baumschabl M., Mattanovich D. et al. Synthetic promoter and CRISPR activation tools in <i>Komagataella phaffii</i>. <i>ACS Synth. Biol.</i> 2020; Pérez-García F. et al. Synthetic <i>ribGCAH</i> operon RBS design in <i>C. glutamicum</i> (8.3 g/L). 2024.",
            "[23] Wagner J. M., Liu L., Yuan S.-F., Venkataraman M. V., Abate A. R., Alper H. S. A comparative analysis of single cell and droplet-based FACS for improving production phenotypes: riboflavin overproduction in <i>Yarrowia lipolytica</i>. <i>Metabolic Engineering</i>. 2018;47:346–356; Elsharawy A. et al. Riboflavin engineering in <i>Y. lipolytica</i>. 2026.",
            "[24] Каталог БРЦ ВКПМ НИЦ «Курчатовский институт» — ГосНИИгенетика (Москва) и истекшие патенты РФ на штаммы <i>B. subtilis</i>: <b>ВКПМ B-2694</b> («ВНИИгенетика-304/pMX45», открытый каталог), <b>ВКПМ B-6797</b> (пат. RU 2081906 C1, 12,4 г/л), <b>ВКПМ B-8386</b> и <b>ВКПМ B-8387</b> (пат. RU 2261273 C2, Миронов А. С., Перумов Д. А., Дебабов В. Г. и др. — до 18–21 г/л), а также дрожжевой штамм <i>Ogataea polymorpha</i> <b>ВКПМ Y-201</b>.",
            "[25] Приказ Минсельхоза России от 06.06.2025 № 392 (рег. Минюста РФ № 83535 от 12.09.2025): <b>пункт 32 «Витамин B2 (рибофлавин), его соединения и производные»</b>; Распоряжение Правительства РФ от 05.10.2026 № 2777-р «Стратегия развития биоэкономики в РФ до 2036 г.».",
            "[26] Global Info Research. <i>Global Flavin Adenine Dinucleotide (FAD) Market 2026–2032</i> (Feb 2026: ~620 t global volume, $74.1M in 2025, $107M forecast); каталоги GMP/IVD-субстанций <i>Glentham Life Sciences</i> (CAS 84366-81-4, ≥95% HPLC) и <i>MP Biomedicals</i> (2026).",
            "[27] Таможенно-ценовая аналитика <i>Tridge</i> по <i>Riboflavin-5'-Phosphate Sodium</i> (обновление 11.11.2026: средние экспортные цены мая 2026 г. — Индия $99,52/кг, КНР $62,58–120,31/кг, Нидерланды $203,02/кг, Швейцария $383,92/кг); котировки опта витаминов <i>Nutri Avenue</i> (B₂ 98% $39–44/кг, ФМН-Na USP $82–85/кг); <i>Dataintelo Riboflavin-5'-Phosphate Sodium Market Report</i> ($408.6M in 2024, CAGR 5.7%).",
            "[28] Государственный реестр лекарственных средств (ГРЛС) Минздрава России, выгрузка ноябрь 2026 г.: регистрационные удостоверения препарата ЖНВЛП <b>«Цитофлавин®»</b> (ООО «НТФФ «ПОЛИСАН», РУ ЛП-№(000973)-(РГ-RU) и ЛС-001767; поставщики АФС рибофлавина мононуклеотида — <i>DSM Nutritional Products France</i> П N013263/01, <i>Supriya Lifescience</i> Индия и <i>ООО «Полисинтез»</i> ЛСР-002295/07), а также реестровые записи АФС <i>Harman Finochem Ltd.</i> (ЛСР-007577/08) и <i>Danashmand Organic Pvt. Ltd.</i> (ЛС-002018).",
            "[29] Данные о производственных мощностях химического синтеза АФС ООО «Полисинтез» (г. Белгород, 100% ДЗО ООО «НТФФ «ПОЛИСАН»): интервью гендиректора А. Г. Пшеничного «Фармвестнику» от 04.08.2026 (265 т АФС в 2025 г., план 285 т в 2026 г., 13 наименований химического синтеза) и аналитика RNC Pharma по продажам «Цитофлавина» и дженерика «Флатоцин» (ООО «Гротекс» / Solopharm).",
            "[30] Научный путеводитель команды проекта DBTL: <i>«Рибофлавин → ФМН → ФАД: как устроен путь, каких результатов мы ищем и что именно проверять в лаборатории»</i> (15 стр., 22.11.2026, локальная копия в репозитории: <i>docs/Riboflavin_FMN_FAD_nauchnyi_putevoditel_2026.pdf</i>).",
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
            title="Обоснование выбора каскада Рибофлавин (B2) -> 5'-ФМН -> ФАД для цикла DBTL",
            author="Проект DBTL",
        )
        frame = Frame(left_m, bot_m, usable_w, page_h - top_m - bot_m, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=draw_header_footer)])
        doc.build(build_story(toc_in))
        return doc.page, recorded

    p1_count, toc_exact = render_pass(OUT_PDF, {})
    p2_count, toc_verify = render_pass(OUT_PDF, toc_exact)
    assert p1_count == p2_count and toc_exact == toc_verify, f"Two-pass TOC mismatch: {p1_count} vs {p2_count}"

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
