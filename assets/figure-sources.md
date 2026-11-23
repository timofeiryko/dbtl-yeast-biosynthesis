# Происхождение иллюстраций отчёта

## Экспериментальные данные

- `candida-kinetics-2022.png` — рисунок 4 из статьи Petrovska Y., Lyzak O., Ruchala J., Dmytruk K., Sibirny A. *Co-Overexpression of RIB1 and RIB6 Increases Riboflavin Production in the Yeast Candida famata.* Fermentation 2022;8(4):141. DOI: https://doi.org/10.3390/fermentation8040141. Лицензия [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Исходный файл: https://mdpi-res.com/d_attachment/fermentation/fermentation-08-00141/article_deploy/html/images/fermentation-08-00141-g004.png.
- `candida-kinetics-2022-ef.png` — кадрирование нижнего ряда рисунка (панели E–F) без изменения данных, шкал и легенд: динамика биомассы (слева) и рибофлавина (справа) для штаммов BRPI и BRPI/RIB1–RIB6 на среде YNB за 120 ч.

## Графики и схемы

- `riboflavin-market-*.png` (4 графика) — строятся скриптом `scripts/build_riboflavin_report.py` по данным `docs/riboflavin-market-data.json` (WITS/UN Comtrade, FEFAC, VEGA).
- `riboflavin-volume.png` — расчёт объёма культуры на 1 кг выделенного B2 при сценарной степени извлечения 80%.
- `riboflavin-sampling.png` — схема параллельного отбора проб для контроля материального баланса.
- `riboflavin-denovo-pathway.png` — схема биосинтеза рибофлавина de novo в дрожжах с указанием генов *C. famata* и *K. phaffii* и ключевых мишеней DBTL.
- `riboflavin-pathway.png` — схема последовательных ферментативных превращений B2 → FMN → FAD на базе 2D-структур.

## Графики и схемы отчёта по фосфолипазе C (PLC)

- `plc-benchmarks.png`, `plc-oil-production-rf.png`, `plc-reactor-economics.png` — строятся скриптом `scripts/build_plc_report.py` по структурированным данным `docs/plc-market-data.json` (OleoScope / Масложировой союз России, Val et al. 2023–2025, Elena et al. 2016–2017, Marchisio et al. 2022).
- `plc-reaction-scheme.png`, `plc-assay-workflow.png`, `plc-dbtl-architecture.png` — оригинальные инженерные блок-схемы механизма PLC vs PLA1, трёхуровневой аналитической воронки (`p-NPPC` + AOCS Ca 12-55) и двухконтурного AI-DBTL в *K. phaffii*.
