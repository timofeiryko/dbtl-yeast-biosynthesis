# Происхождение иллюстраций текущего отчёта

## Экспериментальная иллюстрация

`candida-kinetics-2022.png` — рисунок 4 из Petrovska Y., Lyzak O., Ruchala J., Dmytruk K., Sibirny A. **Co-Overexpression of RIB1 and RIB6 Increases Riboflavin Production in the Yeast Candida famata.** Fermentation 2022;8(4):141. DOI: https://doi.org/10.3390/fermentation8040141

© 2022 авторы; лицензия [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), подтверждена на странице издателя. Исходное изображение: https://mdpi-res.com/d_attachment/fermentation/fermentation-08-00141/article_deploy/html/images/fermentation-08-00141-g004.png

`candida-kinetics-2022-ef.png` — кадрирование нижнего ряда, панели E–F. Удалены верхние ряды и остатки их подписей; данные, шкалы и легенды не изменены. Слева биомасса, справа B2; BRPI и BRPI/RIB1–RIB6, YNB, 120 ч. Это не рисунок для другой линии BRP и не fed-batch.

## Перестроенные графики

Четыре `riboflavin-market-*.png` строятся скриптом `scripts/build_riboflavin_report.py` по `docs/riboflavin-market-data.json`. В JSON указаны первичные URL, периоды, единицы и границы применимости: WITS/UN Comtrade, FEFAC и VEGA. Изображения сторонних рыночных отчётов не копируются.

`riboflavin-volume.png` — расчёт по явно указанному сценарию извлечения 80%, а не экспериментальные данные.

`riboflavin-sampling.png` — собственная схема контроля материального баланса проб; не изображение из статьи и не валидированный СОП.

`riboflavin-pathway.png` — сборка имеющихся структур B2/FMN/FAD и подписей реакций. Структуры не заменяют полный стехиометрический баланс пути.