"""Render the single-source Russian thermostable PLC DBTL report and its figures."""
from html import escape
from pathlib import Path
import json
import re

import fitz
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image as PILImage, ImageDraw, ImageFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate, Frame, HRFlowable, Image, KeepTogether, PageBreak,
    PageTemplate, Paragraph, Spacer, Table, TableStyle,
)

# Reuse established Cyrillic fonts and typography styles from existing report builders.
from build_fmn_report import MPL_FONT, ROOT_DIR
from build_riboflavin_report import (
    BODY, BULLET, CAPTION, CELL, CELL_HEAD, H1, H2_TOP, H3, INK, MUTED,
    RULE, SMALL, SUBTITLE, TEAL, WIDTH, _build_toc, markup,
)

ROOT = ROOT_DIR
SOURCE = ROOT / 'docs/plc-dbtl-report.md'
OUT = ROOT / 'plc-dbtl-report.pdf'
ASSETS = ROOT / 'assets'
SCRATCH = ROOT / '.context/plc-report'
MARKET_DATA = ROOT / 'docs/plc-market-data.json'

SLATE = '#60849F'
LIGHT_SLATE = '#88A6BE'
MID_TEAL = '#4B8B74'
SOFT_TEAL = '#A7C7BA'
BG_MINT = '#EBF4EF'
BG_SOFT = '#F4F8F6'
BG_BLUE = '#ECF2F7'
BORDER_TEAL = '#B4C9BD'
BORDER_BLUE = '#B8CBD9'

SUB_BULLET = ParagraphStyle(
    'sub_bullet', parent=BODY, leftIndent=21, firstLineIndent=-8, spaceAfter=2.4,
)
CALLOUT_TEXT = ParagraphStyle(
    'callout_text', parent=BODY, fontSize=8.8, leading=12.2, spaceAfter=0,
)


def _setup_mpl():
    plt.rcParams.update({
        'font.family': MPL_FONT,
        'font.size': 11,
        'axes.labelcolor': INK,
        'text.color': INK,
        'axes.edgecolor': '#9AA7AB',
        'xtick.color': INK,
        'ytick.color': INK,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
    })


def _draw_box(ax, x, y, w, h, title, body='', bg=BG_SOFT, edge=BORDER_TEAL,
              title_color=INK, title_fs=9.6, body_fs=8.4, lw=1.15, radius=0.018):
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle=f'round,pad=0.008,rounding_size={radius}',
        facecolor=bg, edgecolor=edge, linewidth=lw,
    )
    ax.add_patch(patch)
    cx = x + w / 2
    if body:
        ax.text(cx, y + h * 0.72, title, ha='center', va='center',
                fontsize=title_fs, fontweight='bold', color=title_color)
        ax.text(cx, y + h * 0.30, body, ha='center', va='center',
                fontsize=body_fs, color=INK, linespacing=1.28)
    else:
        ax.text(cx, y + h / 2, title, ha='center', va='center',
                fontsize=title_fs, fontweight='bold', color=title_color, linespacing=1.25)


def _arrow(ax, start, end, color=TEAL, lw=1.7, rad=0.0):
    ax.annotate(
        '', xy=end, xytext=start,
        arrowprops=dict(
            arrowstyle='-|>', color=color, lw=lw, mutation_scale=13,
            connectionstyle=f'arc3,rad={rad}',
        ),
    )


def fig_benchmarks(data):
    _setup_mpl()
    rows = data['plc_temperature_benchmarks']['variants']
    labels = [r['label'] for r in rows]
    temps = [r['temp_c'] for r in rows]
    times = [r['time_min'] for r in rows]
    doses = [r['dose_mg_kg'] for r in rows]
    palette = {
        'commercial': LIGHT_SLATE,
        'reference': SLATE,
        'benchmark': MID_TEAL,
        'target': TEAL,
    }
    bar_colors = [palette[r['category']] for r in rows]
    y_pos = list(range(len(rows)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0), dpi=180, sharey=True)

    # Left panel: operating temperature
    ax1.axvspan(75, 82, color='#E4F1EA', zorder=0)
    ax1.axvline(75, color=TEAL, ls='--', lw=1.0, alpha=0.65, zorder=1)
    bars1 = ax1.barh(y_pos, temps, color=bar_colors, height=0.56, zorder=2)
    for bar, t in zip(bars1, temps):
        ax1.text(t + 1.6, bar.get_y() + bar.get_height() / 2, f'{t} °C',
                 va='center', fontsize=10, fontweight='bold' if t >= 80 else 'normal', color=INK)
    ax1.set_yticks(y_pos, labels)
    ax1.invert_yaxis()
    ax1.set_xlim(0, 102)
    ax1.set_xlabel('Рабочая температура дегумминга, °C')
    ax1.text(78.5, -0.58, 'Режим МЭЗа\n(75–80 °C)', ha='center', va='bottom',
             fontsize=8.5, color=TEAL, fontweight='bold')
    ax1.grid(axis='x', alpha=0.18, zorder=0)

    # Right panel: reaction time + active protein dose annotation
    bars2 = ax2.barh(y_pos, times, color=bar_colors, height=0.56, zorder=2)
    for bar, tm, dose in zip(bars2, times, doses):
        dose_str = f'{dose:.1f}'.replace('.', ',')
        ax2.text(tm + 2.5, bar.get_y() + bar.get_height() / 2,
                 f'{tm} мин · {dose_str} мг/кг', va='center', fontsize=9.8,
                 fontweight='bold' if tm <= 30 else 'normal', color=INK)
    ax2.set_xlim(0, 185)
    ax2.set_xlabel('Длительность реакции (мин) и доза активного белка')
    ax2.grid(axis='x', alpha=0.18, zorder=0)

    fig.tight_layout(w_pad=1.8)
    fig.subplots_adjust(top=0.88)
    fig.savefig(ASSETS / 'plc-benchmarks.png')
    plt.close(fig)


def fig_reaction_scheme():
    _setup_mpl()
    fig, ax = plt.subplots(figsize=(10.5, 4.7), dpi=180)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Section A banner
    ax.text(0.01, 0.96, 'А. Механизм дегумминга: PLC против PLA1',
            fontsize=10.5, fontweight='bold', color=TEAL, ha='left')

    # Substrate box
    _draw_box(ax, 0.01, 0.56, 0.22, 0.30, 'Фосфолипиды\nсырого масла',
              'PC, PE (глицерин +\n2 жирные кислоты +\nфосфатная головка)',
              bg=BG_BLUE, edge=BORDER_BLUE, title_fs=9.6, body_fs=8.3)

    # Top path: PLC
    _arrow(ax, (0.24, 0.78), (0.32, 0.81), color=TEAL)
    _draw_box(ax, 0.33, 0.72, 0.21, 0.19, 'PLC (EC 3.1.4.3)',
              'Гидролиз связи sn-3\n(80 °C, 30 мин)',
              bg=BG_MINT, edge=TEAL, title_color=TEAL, title_fs=9.5, body_fs=8.2)
    _arrow(ax, (0.55, 0.81), (0.61, 0.81), color=TEAL)
    _draw_box(ax, 0.62, 0.72, 0.37, 0.19, '1,2-DAG остаётся в масле (+1,5–2,0%)',
              'Фосфохолин / фосфоэтаноламин → в воду\nБез роста свободных жирных кислот (СЖК)',
              bg='#DDF0E6', edge=TEAL, title_color=TEAL, title_fs=9.4, body_fs=8.2)

    # Bottom path: PLA1
    _arrow(ax, (0.24, 0.64), (0.32, 0.59), color=SLATE)
    _draw_box(ax, 0.33, 0.49, 0.21, 0.18, 'PLA1 (сравнение)',
              'Гидролиз связи sn-1\n(глубокий P < 10 ppm)',
              bg=BG_SOFT, edge=BORDER_BLUE, title_color=SLATE, title_fs=9.3, body_fs=8.1)
    _arrow(ax, (0.55, 0.58), (0.61, 0.58), color=SLATE)
    _draw_box(ax, 0.62, 0.49, 0.37, 0.18, 'Образуются СЖК + лизофосфолипид',
              'Выше кислотное число и потери масла\nна нейтрализации / дистилляции',
              bg='#F3F6F8', edge=BORDER_BLUE, title_color='#3E5363', title_fs=9.3, body_fs=8.1)

    # Divider line
    ax.plot([0.01, 0.99], [0.43, 0.43], color=RULE, lw=1.0)

    # Section B banner
    ax.text(0.01, 0.38, 'Б. Четыре зоны управляемого тюнинга PLC в K. phaffii (245 а.к., 28,5 кДа)',
            fontsize=10.5, fontweight='bold', color=TEAL, ha='left')

    cards = [
        ('1. Активный центр\n(3 × Zn²⁺)', 'Trp1, His14/69/118/128/142\nСтрогая фиксация ≤ 5 Å\n(запрет мутаций)', '#F5F0EC', '#CBB5A4', '#7A4B2A'),
        ('2. Карман головки\n(селективность)', 'Зона Phe66 (F66Y и др.)\nУсиление гидролиза PE\nи адаптация под сырьё', BG_SOFT, BORDER_TEAL, TEAL),
        ('3. Петли и ядро\n(термостабильность)', 'Пролины в гибких петлях,\nсолевые мостики, консенсус\nЦель: T50 ≥ 80 °C', BG_MINT, TEAL, TEAL),
        ('4. Секреция дрожжей\n(без метанола)', 'Точный стык Kex2–Trp1,\nпромоторы GAP/GTH1/AOX1m,\nсигналы и копийность', BG_BLUE, BORDER_BLUE, SLATE),
    ]
    cw = 0.232
    gap = 0.017
    for idx, (title, body, bg, edge, tc) in enumerate(cards):
        x0 = 0.01 + idx * (cw + gap)
        _draw_box(ax, x0, 0.04, cw, 0.30, title, body,
                  bg=bg, edge=edge, title_color=tc, title_fs=9.0, body_fs=7.9)

    fig.tight_layout()
    fig.savefig(ASSETS / 'plc-reaction-scheme.png')
    plt.close(fig)


def fig_assay_workflow():
    _setup_mpl()
    fig, ax = plt.subplots(figsize=(10.5, 4.3), dpi=180)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Tier 1
    _draw_box(ax, 0.01, 0.70, 0.16, 0.24, 'Уровень 1\n100% вариантов',
              '48–96 лунок\nза 1 день', bg=TEAL, edge=TEAL,
              title_color='white', title_fs=9.4, body_fs=8.3)
    # Re-draw body white on dark teal
    for txt in ax.texts[-1:]:
        txt.set_color('#E5F2EC')

    _draw_box(ax, 0.20, 0.70, 0.22, 0.24, 'Парный термошок',
              'Супернатант → 2 плашки:\n4 °C (контроль) и\n75–80 °C (30 мин)',
              bg=BG_MINT, edge=BORDER_TEAL, title_fs=9.3, body_fs=8.1)
    _arrow(ax, (0.425, 0.82), (0.455, 0.82), color=TEAL)
    _draw_box(ax, 0.46, 0.70, 0.21, 0.24, 'Тест p-NPPC',
              'Планшетный ридер\n405–410 нм, 15–30 мин\n(буфер + 1 мМ ZnCl2)',
              bg=BG_SOFT, edge=BORDER_TEAL, title_fs=9.3, body_fs=8.1)
    _arrow(ax, (0.675, 0.82), (0.705, 0.82), color=TEAL)
    _draw_box(ax, 0.71, 0.70, 0.28, 0.24, 'Разделение двух метрик',
              'R80 = A405(80 °C) / A405(4 °C)\n→ термостабильность (не зав. от копий)\nA405(4 °C) → объёмная активность',
              bg='#DDF0E6', edge=TEAL, title_color=TEAL, title_fs=9.2, body_fs=7.8)

    # Vertical arrow 1 -> 2
    _arrow(ax, (0.50, 0.69), (0.50, 0.61), color=TEAL)

    # Tier 2
    _draw_box(ax, 0.01, 0.36, 0.16, 0.23, 'Уровень 2\n20–30% лидеров',
              '+ контрольные\nклоны', bg=MID_TEAL, edge=MID_TEAL,
              title_color='white', title_fs=9.3, body_fs=8.2)
    ax.texts[-1].set_color('#EAF5F0')

    _draw_box(ax, 0.20, 0.36, 0.37, 0.23, 'SDS-PAGE (28,5 кДа) + BCA',
              'Чистый титр PLC (г/л) и удельная активность (Ед/мг)\nКонтроль протеолиза и гипергликозилирования',
              bg=BG_SOFT, edge=BORDER_TEAL, title_fs=9.3, body_fs=8.1)
    _arrow(ax, (0.575, 0.475), (0.605, 0.475), color=MID_TEAL)
    _draw_box(ax, 0.61, 0.36, 0.38, 0.23, 'Копийность qPCR и Tm (топ-3)',
              'Отделение эффекта промотора от мультикопийности\nОценка температуры плавления Tm (nanoDSF, ЦКП)',
              bg=BG_BLUE, edge=BORDER_BLUE, title_color=SLATE, title_fs=9.3, body_fs=8.1)

    # Vertical arrow 2 -> 3
    _arrow(ax, (0.50, 0.35), (0.50, 0.27), color=SLATE)

    # Tier 3
    _draw_box(ax, 0.01, 0.03, 0.16, 0.23, 'Уровень 3\n6–8 вариантов',
              'Лидеры +\nбаза на цикл', bg=SLATE, edge=SLATE,
              title_color='white', title_fs=9.3, body_fs=8.2)
    ax.texts[-1].set_color('#EBF2F7')

    _draw_box(ax, 0.20, 0.03, 0.37, 0.23, 'Дегумминг сырого масла (10–50 г)',
              'Подсолнечное и соевое масло · 80 °C, 30 мин\n3% воды · доза 3–5 мг активного белка/кг',
              bg=BG_MINT, edge=TEAL, title_color=TEAL, title_fs=9.3, body_fs=8.1)
    _arrow(ax, (0.575, 0.145), (0.605, 0.145), color=SLATE)
    _draw_box(ax, 0.61, 0.03, 0.38, 0.23, 'AOCS Ca 12-55 / ICP-OES + 1,2-DAG',
              'Остаточный фосфор (мг/кг), прирост 1,2-DAG,\nкислотное число (СЖК) и тест на гемолиз',
              bg=BG_BLUE, edge=BORDER_BLUE, title_color=INK, title_fs=9.2, body_fs=8.1)

    fig.tight_layout()
    fig.savefig(ASSETS / 'plc-assay-workflow.png')
    plt.close(fig)


def fig_dbtl_architecture():
    _setup_mpl()
    fig, ax = plt.subplots(figsize=(10.5, 4.7), dpi=180)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Contour 1 container
    c1 = FancyBboxPatch((0.01, 0.54), 0.98, 0.43, boxstyle='round,pad=0.008,rounding_size=0.02',
                        facecolor='#F7FBF9', edgecolor=TEAL, lw=1.2)
    ax.add_patch(c1)
    ax.text(0.03, 0.935, 'КОНТУР 1 · Белковый AI-DBTL: аминокислотная последовательность PLC (цель: 75–80 °C за 30 мин, чистый SEQ ID NO)',
            fontsize=9.4, fontweight='bold', color=TEAL, ha='left', va='center')

    steps1 = [
        ('Design', 'ThermoMPNN + ESM-2\n+ консенсус + фикс. Zn²⁺\nЦикл 1: zero-shot (24–36)'),
        ('Build', 'ПЦР-мутагенез от\nBcPLC-F66Y в K. phaffii\nСтык Kex2–Trp1'),
        ('Test', '96 лунок p-NPPC\nПарный прогрев 80 °C\nТоп-6 → сырое масло'),
        ('Learn', 'Суррогат (GP / Ridge)\nПарето: R80 × активность\nОтбор 4–10 мутаций'),
    ]
    bw, gap = 0.215, 0.028
    for i, (t, b) in enumerate(steps1):
        x0 = 0.03 + i * (bw + gap)
        _draw_box(ax, x0, 0.62, bw, 0.26, t, b,
                  bg=BG_MINT if i in (0, 3) else 'white',
                  edge=TEAL if i in (0, 3) else BORDER_TEAL,
                  title_color=TEAL, title_fs=9.8, body_fs=8.0)
        if i < 3:
            _arrow(ax, (x0 + bw + 0.004, 0.75), (x0 + bw + gap - 0.004, 0.75), color=TEAL)

    # Feedback loop inside Contour 1
    ax.annotate(
        '', xy=(0.135, 0.615), xytext=(0.865, 0.615),
        arrowprops=dict(arrowstyle='-|>', color=TEAL, lw=1.3, ls='--',
                        connectionstyle='arc3,rad=0.07', mutation_scale=12),
    )
    ax.text(0.50, 0.562, 'Замкнутые циклы 1–3 · проспективное сравнение ML vs экспертный и случайный контроль',
            ha='center', va='center', fontsize=8.3, color=TEAL, fontweight='bold')

    # Transition between Contour 1 and Contour 2
    _arrow(ax, (0.50, 0.535), (0.50, 0.455), color=INK, lw=1.9)
    ax.text(0.515, 0.495, 'Топ-2–3 термостабильных варианта PLC', ha='left', va='center',
            fontsize=8.7, fontweight='bold', color=INK)

    # Contour 2 container
    c2 = FancyBboxPatch((0.01, 0.03), 0.98, 0.415, boxstyle='round,pad=0.008,rounding_size=0.02',
                        facecolor='#F5F8FA', edgecolor=SLATE, lw=1.2)
    ax.add_patch(c2)
    ax.text(0.03, 0.41, 'КОНТУР 2 · Штаммовый DBTL: безметанольная секреция в K. phaffii (цель: ≥ 4–5 г/л в 5-л биореакторе, выход UF ≥ 80%)',
            fontsize=9.4, fontweight='bold', color=SLATE, ha='left', va='center')

    steps2 = [
        ('Промоторы', 'P_GAP, P_GTH1,\nдерепресс. P_AOX1\n(без метанола)'),
        ('Секреция', 'Сигналы α, Ost1-α,\n0030-α, SCW11-α +\nHAC1 / PDI1 / Δcog6'),
        ('Копийность', 'Отбор 4–8 копий (qPCR)\nСтабильность 5 пассажей\nбез антибиотика'),
        ('Биореактор 5 л', 'Fed-batch без метанола\nUF 10 кДа (глицерин)\nДегумминг 80 °C, 30 мин'),
    ]
    for i, (t, b) in enumerate(steps2):
        x0 = 0.03 + i * (bw + gap)
        _draw_box(ax, x0, 0.09, bw, 0.27, t, b,
                  bg=BG_BLUE if i < 3 else '#DDF0E6',
                  edge=BORDER_BLUE if i < 3 else TEAL,
                  title_color=SLATE if i < 3 else TEAL, title_fs=9.6, body_fs=8.0)
        if i < 3:
            _arrow(ax, (x0 + bw + 0.004, 0.225), (x0 + bw + gap - 0.004, 0.225), color=SLATE)

    fig.tight_layout()
    fig.savefig(ASSETS / 'plc-dbtl-architecture.png')
    plt.close(fig)


def fig_oil_production(data):
    _setup_mpl()
    seasons = data['russia_vegetable_oils_mln_t']['seasons']
    labels = [s['season'] for s in seasons]
    sun = [s['sunflower'] for s in seasons]
    rape = [s['rapeseed'] for s in seasons]
    soy = [s['soybean'] for s in seasons]
    other = [s['other'] for s in seasons]
    totals = [s['total'] for s in seasons]
    x = list(range(len(seasons)))

    fig, ax = plt.subplots(figsize=(10.0, 3.8), dpi=180)
    w = 0.46
    b1 = ax.bar(x, sun, width=w, color=TEAL, label='Подсолнечное')
    b2 = ax.bar(x, rape, width=w, bottom=sun, color=SLATE, label='Рапсовое')
    bot3 = [a + b for a, b in zip(sun, rape)]
    b3 = ax.bar(x, soy, width=w, bottom=bot3, color=SOFT_TEAL, label='Соевое')
    bot4 = [a + b for a, b in zip(bot3, soy)]
    ax.bar(x, other, width=w, bottom=bot4, color='#D2DCE0', label='Прочие')

    for i in x:
        ax.text(i, sun[i] / 2, f'{sun[i]:.2f}'.replace('.', ','),
                ha='center', va='center', color='white', fontweight='bold', fontsize=10.5)
        ax.text(i, sun[i] + rape[i] / 2, f'{rape[i]:.2f}'.replace('.', ','),
                ha='center', va='center', color='white', fontweight='bold', fontsize=9.8)
        ax.text(i, bot3[i] + soy[i] / 2, f'{soy[i]:.2f}'.replace('.', ','),
                ha='center', va='center', color=INK, fontweight='bold', fontsize=9.5)
        ax.text(i, totals[i] + 0.25, f'{totals[i]:.2f} млн т'.replace('.', ','),
                ha='center', va='bottom', color=INK, fontweight='bold', fontsize=11)

    ax.set_xticks(x, [f'Сезон {lbl}' for lbl in labels])
    ax.set_ylim(0, 12.1)
    ax.set_ylabel('Производство масел в РФ, млн т')
    ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.06), ncol=4)
    ax.grid(axis='y', alpha=0.18)
    fig.tight_layout()
    fig.savefig(ASSETS / 'plc-oil-production-rf.png')
    plt.close(fig)


def fig_reactor_economics(data):
    _setup_mpl()
    sc = data['reactor_1m3_scenarios']
    cycles = sc['annual_hours'] / (sc['cycle_hours'] + sc['turnaround_hours'])
    titers = sc['titers_g_l']
    # Annual active protein (kg) = titer (g/L) * 1000 L / 1000 g/kg * recovery * cycles
    protein_kg = [t * sc['recovery'] * cycles for t in titers]
    # Treated oil (thousand t/yr) = protein_kg * 1000 g/kg / dose_g_per_t / 1000
    oil_kt = [pk / sc['dose_g_per_t_oil'] for pk in protein_kg]

    mez = data['mez_economics_150kt']
    extra_oil_t = mez['annual_oil_t'] * (mez['extra_yield_pct'] / 100.0)
    rev_m = extra_oil_t * mez['oil_price_usd_t'] / 1e6
    cost_mid = (mez['enzyme_cost_low_usd_m'] + mez['enzyme_cost_high_usd_m']) / 2
    cost_half = (mez['enzyme_cost_high_usd_m'] - mez['enzyme_cost_low_usd_m']) / 2
    net_mid = rev_m - cost_mid

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.8), dpi=180)

    # Left panel: 1 m3 bioreactor capacity
    y = list(range(len(titers)))
    bar_colors = [LIGHT_SLATE, MID_TEAL, TEAL]
    bars = ax1.barh(y, oil_kt, color=bar_colors, height=0.52, zorder=2)
    ax1.axvline(150, color='#9C4A32', ls='--', lw=1.2, zorder=3)
    ax1.text(152, 0.05, 'Средний МЭЗ\n150 тыс. т/год', color='#9C4A32',
             fontsize=8.6, va='top', fontweight='bold')
    for bar, val, pk in zip(bars, oil_kt, protein_kg):
        ax1.text(val + 4, bar.get_y() + bar.get_height() / 2,
                 f'{val:.0f} тыс. т ({pk:.0f} кг белка)',
                 va='center', fontsize=9.3, color=INK,
                 fontweight='bold' if val >= 150 else 'normal')
    ax1.set_yticks(y, [f'{t:.1f} г/л'.replace('.', ',') for t in titers])
    ax1.invert_yaxis()
    ax1.set_xlim(0, 275)
    ax1.set_xlabel('Объём обслуживаемого масла, тыс. т/год')
    ax1.set_title('Выпуск 1 м³ реактора (доза 3 г/т, UF 80%)', fontsize=10.2,
                  fontweight='bold', color=INK, pad=10)
    ax1.grid(axis='x', alpha=0.18, zorder=0)

    # Right panel: 150 kt/yr MЭЗ economics
    cats = ['Доп. выручка\n(+1,5% масла)', 'Затраты на PLC\n(диапазон)', 'Чистая выгода\nМЭЗа (CAPEX = 0)']
    vals = [rev_m, cost_mid, net_mid]
    errs = [0.0, cost_half, cost_half]
    cols = [SLATE, '#B9836B', TEAL]
    bars2 = ax2.bar(range(3), vals, color=cols, width=0.50, yerr=errs,
                    capsize=5, error_kw=dict(ecolor=INK, lw=1.2), zorder=2)
    labels_txt = [
        f'+${rev_m:.2f} млн'.replace('.', ','),
        f'${mez["enzyme_cost_low_usd_m"]:.2f}–{mez["enzyme_cost_high_usd_m"]:.2f} млн'.replace('.', ','),
        f'+${rev_m - mez["enzyme_cost_high_usd_m"]:.2f}–{rev_m - mez["enzyme_cost_low_usd_m"]:.2f} млн'.replace('.', ','),
    ]
    for bar, txt, err in zip(bars2, labels_txt, errs):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + err + 0.08,
                 txt, ha='center', va='bottom', fontsize=9.3, fontweight='bold', color=INK)
    ax2.set_xticks(range(3), cats)
    ax2.set_ylim(0, 2.55)
    ax2.set_ylabel('Млн долларов США / год')
    ax2.set_title('Экономика МЭЗа 150 тыс. т/год ($900/т)', fontsize=10.2,
                  fontweight='bold', color=INK, pad=10)
    ax2.grid(axis='y', alpha=0.18, zorder=0)

    fig.tight_layout(w_pad=2.0)
    fig.savefig(ASSETS / 'plc-reactor-economics.png')
    plt.close(fig)


def figures():
    data = json.loads(MARKET_DATA.read_text(encoding='utf-8'))
    fig_benchmarks(data)
    fig_reaction_scheme()
    fig_assay_workflow()
    fig_dbtl_architecture()
    fig_oil_production(data)
    fig_reactor_economics(data)


class PLCReport(BaseDocTemplate):
    def __init__(self, filename):
        super().__init__(
            str(filename), pagesize=A4, leftMargin=18.5 * mm,
            rightMargin=18.5 * mm, topMargin=15.5 * mm, bottomMargin=15.5 * mm,
            title='Фосфолипаза C в дрожжах: от воспроизведения секреции к AI-DBTL для рафинации масел',
            author='DBTL project',
        )
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


def _col_widths(count):
    if count == 2:
        return [WIDTH * 0.30, WIDTH * 0.70]
    if count == 3:
        return [WIDTH * 0.24, WIDTH * 0.38, WIDTH * 0.38]
    if count == 4:
        return [WIDTH * 0.22, WIDTH * 0.16, WIDTH * 0.31, WIDTH * 0.31]
    if count == 5:
        return [WIDTH * 0.20, WIDTH * 0.18, WIDTH * 0.15, WIDTH * 0.24, WIDTH * 0.23]
    return [WIDTH / count] * count


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
        raw_line = lines[index]
        line = raw_line.strip()
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
        elif line.startswith('> '):
            chunks = [line[2:].strip()]
            while index + 1 < len(lines) and lines[index + 1].strip().startswith('> '):
                index += 1
                chunks.append(lines[index].strip()[2:].strip())
            callout_tbl = Table(
                [[Paragraph(markup(' '.join(chunks)), CALLOUT_TEXT)]],
                colWidths=[WIDTH], hAlign='LEFT',
            )
            callout_tbl.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#F1F7F4')),
                ('LINEBEFORE', (0, 0), (0, -1), 2.5, colors.HexColor(TEAL)),
                ('LEFTPADDING', (0, 0), (-1, -1), 8.5),
                ('RIGHTPADDING', (0, 0), (-1, -1), 7.5),
                ('TOPPADDING', (0, 0), (-1, -1), 5.5),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 5.5),
            ]))
            result.append(KeepTogether([Spacer(1, 2.5), callout_tbl, Spacer(1, 5.5)]))
        elif line.startswith('|'):
            rows = []
            while index < len(lines) and lines[index].startswith('|'):
                values = [v.strip() for v in lines[index].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[:\- ]+', v) for v in values):
                    cell_style = CELL_HEAD if not rows else CELL
                    rows.append([Paragraph(markup(v), cell_style) for v in values])
                index += 1
            table = Table(rows, colWidths=_col_widths(len(rows[0])), repeatRows=1, hAlign='LEFT')
            table.setStyle(TableStyle([
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E6EFEA')),
                ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F6F9F7')]),
                ('LEFTPADDING', (0, 0), (-1, -1), 5.5),
                ('RIGHTPADDING', (0, 0), (-1, -1), 5.5),
                ('TOPPADDING', (0, 0), (-1, -1), 3.3),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 3.3),
                ('LINEBELOW', (0, 0), (-1, 0), 0.6, colors.HexColor('#B4C9BD')),
                ('LINEBELOW', (0, -1), (-1, -1), 0.4, colors.HexColor('#DCE5DF')),
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
            image_width = min(WIDTH * 0.96, 138 * mm / ratio)
            img_flow = Image(str(path), width=image_width, height=image_width * ratio)
            lookahead = index + 1
            while lookahead < len(lines) and not lines[lookahead].strip():
                lookahead += 1
            if lookahead < len(lines) and lines[lookahead].strip().startswith('Рисунок '):
                chunks = [lines[lookahead].strip()]
                lookahead += 1
                while lookahead < len(lines) and lines[lookahead].strip() and not lines[lookahead].startswith(('#', '|', '- ', '![', '<!--', '> ')):
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
            is_sub = raw_line.startswith('  ')
            style = SUB_BULLET if is_sub else BULLET
            if line.startswith('- '):
                prefix = '– ' if is_sub else '• '
                result.append(Paragraph(prefix + markup(line[2:]), style))
            else:
                result.append(Paragraph(markup(line), style))
        else:
            chunks = [line]
            while (
                index + 1 < len(lines)
                and lines[index + 1].strip()
                and not lines[index + 1].strip().startswith(('#', '|', '- ', '![', '<!--', '> '))
                and not re.match(r'^\d+\.\s', lines[index + 1].strip())
            ):
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


def validate_report(text):
    """Runnable checks: market/reactor arithmetic, citation coverage, and final PDF geometry."""
    data = json.loads(MARKET_DATA.read_text(encoding='utf-8'))
    for season in data['russia_vegetable_oils_mln_t']['seasons']:
        calc = season['sunflower'] + season['rapeseed'] + season['soybean'] + season['other']
        assert abs(calc - season['total']) < 1e-6, (season['season'], calc, season['total'])
    sc = data['reactor_1m3_scenarios']
    cycles = sc['annual_hours'] / (sc['cycle_hours'] + sc['turnaround_hours'])
    assert abs(cycles - 66.6666667) < 1e-4
    assert round(4.5 * sc['recovery'] * cycles) == 240
    assert round(10.0 * sc['recovery'] * cycles) == 533
    mez = data['mez_economics_150kt']
    extra_oil = mez['annual_oil_t'] * mez['extra_yield_pct'] / 100.0
    assert abs(extra_oil - 2250.0) < 1e-6
    assert abs(extra_oil * mez['oil_price_usd_t'] / 1e6 - 2.025) < 1e-6

    body_part, _, refs_part = text.partition('## 13. Список источников')
    assert refs_part, 'Missing section 13 references'
    citations = set()
    for group in re.findall(r'\[([\d, –-]+)\]', body_part):
        for item in group.split(','):
            bounds = re.split('[–-]', item.strip())
            citations.update(range(int(bounds[0]), int(bounds[-1]) + 1))
    references = set(map(int, re.findall(r'^\[(\d+)\]', refs_part, re.M)))
    assert citations == references == set(range(1, 34)), (citations - references, references - citations)

    for forbidden in ('HVO', 'SAF', 'биотоплив'):
        assert forbidden.lower() not in text.lower(), f'Forbidden term leaked: {forbidden}'
    assert not re.search(r'\bQ1\b', text, re.I), 'Forbidden term leaked: Q1'

    with fitz.open(OUT) as doc:
        toc = doc.get_toc()
        assert len(toc) == 13, toc
        toc_links = doc[0].get_links()
        assert len(toc_links) == 13, len(toc_links)
        destinations = [
            int(link['page']) if link['kind'] == 4 else link['page'] + 1
            for link in toc_links
        ]
        assert destinations == [item[2] for item in toc], (destinations, toc)
        char_counts = []
        for page in doc:
            extracted = page.get_text()
            char_counts.append((page.number + 1, len(extracted.strip())))
            assert '<!--' not in extracted
            assert '\ufffd' not in extracted
            assert '**' not in extracted
            assert '`' not in extracted
            assert not re.search(r'\$[^0-9\s]', extracted), (page.number + 1, extracted[:300])
            for block in page.get_text('blocks'):
                assert block[0] >= 0 and block[1] >= 0, (page.number + 1, block)
                assert block[2] <= page.rect.width + 0.5 and block[3] <= page.rect.height + 0.5, (page.number + 1, block)
        low_pages = [(p, c) for p, c in char_counts if c < 450]
        print(f'PASS: {len(references)} cited references, 13 TOC links on page 1, arithmetic OK, {len(doc)} pages verified')
        print('Page char counts:', char_counts)
        if low_pages:
            print('Low-density pages (<450 chars):', low_pages)


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
    links = sum(len(page.get_links()) for page in doc)
    print(f'Built {OUT.name}: {len(doc)} pages, {links} links')
    print('Review renders:', SCRATCH)


def main():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    figures()
    text = SOURCE.read_text(encoding='utf-8')
    first = PLCReport(SCRATCH / 'first-pass.pdf')
    first.build(story(text, {}))
    final = PLCReport(OUT)
    final.build(story(text, first.pages))
    if final.pages != first.pages:
        raise RuntimeError('TOC page numbers changed between passes')
    validate_report(text)
    render_review()


if __name__ == '__main__':
    main()
